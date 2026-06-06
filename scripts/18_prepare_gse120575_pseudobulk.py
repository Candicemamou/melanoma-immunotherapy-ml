from __future__ import annotations

from pathlib import Path
import gzip
import pandas as pd
import numpy as np

print("--- 18 : PRÉPARATION GSE120575 SINGLE-CELL EN PSEUDO-BULK ---")

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"
DATA_PROCESSED.mkdir(exist_ok=True)

EXPR_FILE = DATA_RAW / "8_SINGLE_CELL_GSE120575_expression_TPM_raw.txt.gz"
METADATA_FILE = DATA_PROCESSED / "11_METADATA_GSE120575_clean.csv"

OUT_PSEUDOBULK = DATA_PROCESSED / "18_GSE120575_PSEUDOBULK_TPM_READY.csv"
OUT_SUMMARY = DATA_PROCESSED / "18_GSE120575_PSEUDOBULK_summary.csv"


def normalize_sample_id(value):
    if pd.isna(value):
        return None
    return str(value).strip().replace(" ", "_")


def patient_id_from_cell_id(cell_id: str) -> str | None:
    """Cell IDs look like: A10_P3_M11 -> patient_id = P3"""
    if cell_id is None or pd.isna(cell_id):
        return None
    s = str(cell_id).strip()
    parts = s.split("_")
    if len(parts) < 2:
        return None
    for tok in parts:
        if tok.startswith("P") and tok[1:].isdigit():
            return tok
    if len(parts) >= 2 and parts[1].startswith("P"):
        return parts[1]
    return None


def is_prepost_token(token: str) -> bool:
    if token is None:
        return False
    t = str(token).strip()
    return t.startswith("Pre_") or t.startswith("Post_")


def main():
    metadata = pd.read_csv(METADATA_FILE)
    metadata["sample_id"] = metadata["sample_id"].apply(normalize_sample_id)
    metadata["patient_id"] = metadata["patient_id"].apply(normalize_sample_id)

    metadata_core = metadata[
        ["sample_id", "patient_id", "response", "response_binary", "dataset_role", "notes"]
    ].copy()

    # patient_id -> list of sample_ids
    patient_to_sample_ids: dict[str, list[str]] = {
        pid: sub["sample_id"].tolist() for pid, sub in metadata_core.groupby("patient_id")
    }

    # Parse expression header to get cell ids
    with gzip.open(EXPR_FILE, "rt") as f:
        header_line = f.readline().rstrip("\n")
        header_cols = header_line.split("\t")
        cell_ids = header_cols[1:]

        # Map cell -> sample_id (simple: first sample_id for patient)
        cell_to_sample: dict[str, str] = {}
        for cid in cell_ids:
            pid = patient_id_from_cell_id(cid)
            if pid is None:
                continue
            sample_list = patient_to_sample_ids.get(pid, [])
            if not sample_list:
                continue
            cell_to_sample[cid] = sample_list[0]

        mapped_cells = sum(1 for cid in cell_ids if cid in cell_to_sample)
        print(f"Cells in expression: {len(cell_ids)} | mapped via patient_id: {mapped_cells}")

        if mapped_cells == 0:
            raise RuntimeError(
                "Aucune cellule n'a pu être mappée vers un sample_id via patient_id. Vérifie les conventions des IDs."
            )

        sample_ids = sorted(set(cell_to_sample.values()))

        # Prepare aggregation dicts
        sum_dict: dict[str, dict[str, float]] = {sid: {} for sid in sample_ids}
        cnt_dict: dict[str, dict[str, int]] = {sid: {} for sid in sample_ids}

        # Precompute column->sample_id
        col_samples = [cell_to_sample.get(cid) for cid in cell_ids]

        # IMPORTANT: after the cell header line, the next row(s) are Pre_/Post_ tokens aligned on columns (not genes).
        # Skip until we reach the first real gene row.
        first_gene_fields = None
        while True:
            line = f.readline()
            if not line:
                break
            s = line.rstrip("\n")
            if not s:
                continue
            fields = s.split("\t")
            if len(fields) < 2:
                continue
            gene = fields[0].strip()
            if not gene:
                continue
            if is_prepost_token(gene):
                continue
            first_gene_fields = fields
            break

        if first_gene_fields is None:
            raise RuntimeError("Impossible de trouver une première ligne 'gene' dans le fichier expression.")

        def consume_fields(fields_row: list[str]):
            gene_tok = fields_row[0].strip()
            if not gene_tok or is_prepost_token(gene_tok):
                return
            values_row = fields_row[1:]
            for j, val in enumerate(values_row):
                if j >= len(col_samples):
                    break
                sid = col_samples[j]
                if sid is None:
                    continue
                try:
                    x = float(val)
                except Exception:
                    x = 0.0
                sd = sum_dict[sid]
                cd = cnt_dict[sid]
                sd[gene_tok] = sd.get(gene_tok, 0.0) + x
                cd[gene_tok] = cd.get(gene_tok, 0) + 1

        # consume first gene
        consume_fields(first_gene_fields)

        # consume rest
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            consume_fields(fields)

    # Build dense gene matrix
    all_genes = sorted({g for sid in sample_ids for g in sum_dict[sid].keys()})
    print("Aggregated genes:", len(all_genes))

    if len(all_genes) == 0:
        raise RuntimeError("Aucun gène agrégé. Vérifie le parsing/alignement des lignes 'gene'.")

    df_genes = pd.DataFrame(index=sample_ids, columns=all_genes, dtype=float)
    for sid in sample_ids:
        for gene, ssum in sum_dict[sid].items():
            denom = cnt_dict[sid].get(gene, 1)
            df_genes.at[sid, gene] = ssum / max(denom, 1)

    df_genes = df_genes.reset_index().rename(columns={"index": "sample_id"})

    meta_for_samples = metadata_core.drop_duplicates(subset=["sample_id"]).copy()
    df_ready = meta_for_samples.merge(df_genes, on="sample_id", how="inner")

    clinical_cols = ["sample_id", "patient_id", "response", "response_binary"]
    extra_cols = [c for c in ["dataset_role", "notes"] if c in df_ready.columns]
    ordered = clinical_cols + extra_cols + [c for c in df_ready.columns if c not in clinical_cols + extra_cols]
    df_ready = df_ready[ordered]

    df_ready.to_csv(OUT_PSEUDOBULK, index=False)

    n_genes = df_ready.shape[1] - len(clinical_cols) - len(extra_cols)
    summary = pd.DataFrame(
        [{"dataset": "GSE120575", "n_samples": int(df_ready.shape[0]), "n_genes": int(n_genes)}]
    )
    summary.to_csv(OUT_SUMMARY, index=False)

    print("✅ Pseudo-bulk créé :", OUT_PSEUDOBULK)
    print("✅ Résumé sauvegardé :", OUT_SUMMARY)
    print(summary)


if __name__ == "__main__":
    main()

