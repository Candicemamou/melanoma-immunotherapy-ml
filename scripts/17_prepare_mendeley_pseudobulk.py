from pathlib import Path
import zipfile
import pandas as pd
import numpy as np

print("--- 17 : PRÉPARATION MENDELEY SINGLE-CELL EN PSEUDO-BULK ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"

DATA_PROCESSED.mkdir(exist_ok=True)

MENDELEY_ZIP = DATA_RAW / "8_SINGLE_CELL_MENDELEY_MELANOMA_CITESEQ_raw.zip"

# Fichiers extraits par le script 14
RNA_METADATA_FILE = DATA_PROCESSED / "14_MENDELEY_Melanoma_RNA_metadata.csv"

OUT_COUNTS = DATA_PROCESSED / "17_MENDELEY_RNA_PSEUDOBULK_COUNTS.csv"
OUT_LOGCPM = DATA_PROCESSED / "17_MENDELEY_RNA_PSEUDOBULK_LOGCPM.csv"
OUT_CELL_COUNTS = DATA_PROCESSED / "17_MENDELEY_RNA_PSEUDOBULK_cell_counts.csv"
SUMMARY_FILE = DATA_PROCESSED / "17_MENDELEY_PSEUDOBULK_summary.csv"


# ============================================================
# PARAMÈTRES
# ============================================================

CHUNKSIZE = 1000


# ============================================================
# OUTILS
# ============================================================

def clean_text(value):
    if pd.isna(value):
        return None
    value = str(value).replace('"', "").strip()
    return value if value != "" else None


def standardize_response_binary(value):
    """
    RECIST :
    - Complete Response / Partial Response = 1
    - Progressive Disease = 0
    - Stable Disease = NaN pour l'instant, car biologiquement ambigu
    """
    if pd.isna(value):
        return np.nan

    v = str(value).strip().upper()

    positive = {
        "CR",
        "COMPLETE RESPONSE",
        "COMPLETE_RESPONSE",
        "PR",
        "PARTIAL RESPONSE",
        "PARTIAL_RESPONSE",
        "RESPONDER",
        "RESPONSE",
        "R",
    }

    negative = {
        "PD",
        "PROGRESSIVE DISEASE",
        "PROGRESSIVE_DISEASE",
        "PROGRESSION",
        "NON-RESPONDER",
        "NON_RESPONDER",
        "NONRESPONDER",
        "NR",
    }

    stable = {
        "SD",
        "STABLE DISEASE",
        "STABLE_DISEASE",
    }

    if v in positive:
        return 1

    if v in negative:
        return 0

    if v in stable:
        return np.nan

    if "COMPLETE" in v:
        return 1

    if "PARTIAL" in v:
        return 1

    if "PROGRESS" in v:
        return 0

    if "STABLE" in v:
        return np.nan

    return np.nan


def find_zip_member(zip_obj, pattern):
    """
    Trouve un fichier dans le zip à partir d'un morceau de nom.
    """
    pattern = pattern.lower()

    matches = [
        name for name in zip_obj.namelist()
        if pattern in name.lower()
    ]

    if len(matches) == 0:
        raise FileNotFoundError(f"Aucun fichier contenant '{pattern}' dans le zip.")

    if len(matches) > 1:
        print(f"Plusieurs fichiers trouvés pour {pattern} :")
        for m in matches:
            print("-", m)
        print("On utilise le premier :", matches[0])

    return matches[0]


def read_metadata():
    """
    Lit la metadata RNA extraite par le script 14.
    Si elle n'existe pas, lit directement dans le zip.
    """
    if RNA_METADATA_FILE.exists():
        print("Lecture metadata depuis :", RNA_METADATA_FILE)
        return pd.read_csv(RNA_METADATA_FILE, low_memory=False)

    print("Metadata script 14 absente. Lecture directe dans le zip...")

    with zipfile.ZipFile(MENDELEY_ZIP, "r") as z:
        metadata_member = find_zip_member(z, "Melanoma_RNA_metadata.csv")
        with z.open(metadata_member) as f:
            return pd.read_csv(f, low_memory=False)


def get_counts_header():
    """
    Lit seulement l'en-tête du fichier RNA counts.
    """
    with zipfile.ZipFile(MENDELEY_ZIP, "r") as z:
        counts_member = find_zip_member(z, "Melanoma_RNA_counts.csv")
        with z.open(counts_member) as f:
            header = pd.read_csv(f, nrows=0)
    return counts_member, header.columns.tolist()


def choose_cell_column(metadata, count_cell_cols):
    """
    Essaie de trouver la colonne de metadata qui correspond aux cellules.
    Si aucune colonne ne matche, on utilise l'ordre des lignes.
    """
    count_set = set(map(str, count_cell_cols))

    best_col = None
    best_overlap = 0

    for col in metadata.columns:
        values = metadata[col].dropna().astype(str)
        overlap = len(set(values) & count_set)

        if overlap > best_overlap:
            best_overlap = overlap
            best_col = col

    print("Meilleure colonne cellule trouvée :", best_col)
    print("Overlap avec colonnes counts :", best_overlap)

    if best_col is not None and best_overlap > 0:
        metadata["cell_id_for_counts"] = metadata[best_col].astype(str)
        return metadata, "cell_id_for_counts", "matched_by_column"

    # Fallback : même ordre que les colonnes du fichier counts
    if metadata.shape[0] == len(count_cell_cols):
        print("Aucune colonne cellule claire, mais metadata et counts ont la même longueur.")
        print("On mappe donc les cellules par ordre.")
        metadata = metadata.copy()
        metadata["cell_id_for_counts"] = count_cell_cols
        return metadata, "cell_id_for_counts", "matched_by_order"

    raise ValueError(
        "Impossible d'associer metadata et counts : "
        f"{metadata.shape[0]} lignes metadata vs {len(count_cell_cols)} colonnes counts."
    )


def choose_sample_column(metadata):
    """
    Trouve la colonne sample/patient la plus probable.
    """
    priority_exact = [
        "sample",
        "sample_id",
        "Sample",
        "Sample_ID",
        "patient",
        "patient_id",
        "Patient",
        "Patient_ID",
        "assignment",
        "Assignment",
        "orig.ident",
        "library",
        "Library",
    ]

    for col in priority_exact:
        if col in metadata.columns:
            return col

    # Recherche plus large
    candidates = []
    for col in metadata.columns:
        c = col.lower()
        if any(k in c for k in ["sample", "patient", "assignment", "library", "donor"]):
            if "recist" not in c and "response" not in c:
                candidates.append(col)

    if len(candidates) > 0:
        return candidates[0]

    raise ValueError(
        "Impossible de trouver une colonne sample/patient dans la metadata. "
        f"Colonnes disponibles : {metadata.columns.tolist()}"
    )


def choose_response_column(metadata):
    """
    Trouve la colonne RECIST / réponse.
    """
    for col in metadata.columns:
        if col.lower() == "recist":
            return col

    for col in metadata.columns:
        c = col.lower()
        if "recist" in c or "response" in c:
            return col

    raise ValueError(
        "Impossible de trouver une colonne RECIST/response. "
        f"Colonnes disponibles : {metadata.columns.tolist()}"
    )


def is_valid_sample_id(value):
    """
    Élimine les cellules non assignées ou multiplets visibles dans les barcode files.
    Exemples à exclure : 0/1, 1/0, None, nan.
    """
    value = clean_text(value)

    if value is None:
        return False

    v = value.upper()

    bad_values = {
        "NAN",
        "NA",
        "NONE",
        "UNKNOWN",
        "UNASSIGNED",
        "NEGATIVE",
        "MULTIPLET",
        "DOUBLET",
        "0/1",
        "1/0",
    }

    if v in bad_values:
        return False

    if "/" in v:
        return False

    return True


def build_cell_metadata(metadata, count_cell_cols):
    """
    Construit une table cellule -> sample -> RECIST.
    """
    print("\nColonnes metadata RNA :")
    print(metadata.columns.tolist())

    metadata, cell_col, matching_method = choose_cell_column(metadata, count_cell_cols)

    sample_col = choose_sample_column(metadata)
    response_col = choose_response_column(metadata)

    print("\nColonne cellule utilisée :", cell_col)
    print("Méthode d'association counts/metadata :", matching_method)
    print("Colonne sample utilisée :", sample_col)
    print("Colonne réponse utilisée :", response_col)

    cell_meta = metadata.copy()

    cell_meta["cell_id"] = cell_meta[cell_col].astype(str)
    cell_meta["sample_id"] = cell_meta[sample_col].astype(str).str.strip()
    cell_meta["response"] = cell_meta[response_col]
    cell_meta["response_binary"] = cell_meta["response"].apply(standardize_response_binary)

    before = cell_meta.shape[0]

    cell_meta = cell_meta[cell_meta["sample_id"].apply(is_valid_sample_id)].copy()

    after_sample_filter = cell_meta.shape[0]

    # Garder seulement les cellules présentes dans les counts
    count_set = set(count_cell_cols)
    cell_meta = cell_meta[cell_meta["cell_id"].isin(count_set)].copy()

    after_count_filter = cell_meta.shape[0]

    print("\nCellules metadata avant filtrage :", before)
    print("Cellules après suppression samples invalides :", after_sample_filter)
    print("Cellules présentes dans counts :", after_count_filter)

    print("\nRépartition sample_id :")
    print(cell_meta["sample_id"].value_counts(dropna=False))

    print("\nRépartition response :")
    print(cell_meta["response"].value_counts(dropna=False))

    print("\nRépartition response_binary :")
    print(cell_meta["response_binary"].value_counts(dropna=False))

    # Table sample-level metadata
    sample_meta = (
        cell_meta
        .groupby("sample_id", as_index=False)
        .agg(
            n_cells=("cell_id", "count"),
            response=("response", lambda x: x.dropna().iloc[0] if len(x.dropna()) > 0 else np.nan),
            response_binary=("response_binary", lambda x: x.dropna().iloc[0] if len(x.dropna()) > 0 else np.nan),
        )
    )

    sample_meta["study"] = "MENDELEY_MELANOMA_CITESEQ"
    sample_meta["dataset_role"] = "single_cell_pseudobulk"
    sample_meta["expression_format"] = "pseudobulk RNA counts from single-cell counts"

    # Réordonner
    sample_meta = sample_meta[
        [
            "study",
            "sample_id",
            "response",
            "response_binary",
            "n_cells",
            "expression_format",
            "dataset_role",
        ]
    ]

    return cell_meta, sample_meta


def compute_pseudobulk_counts(counts_member, cell_meta):
    """
    Calcule le pseudo-bulk par somme des counts par sample.
    Le fichier RNA counts est supposé :
    gènes en lignes, cellules en colonnes.
    """
    print("\nCalcul pseudo-bulk RNA counts...")

    cell_to_sample = dict(zip(cell_meta["cell_id"], cell_meta["sample_id"]))

    selected_cells = list(cell_to_sample.keys())

    chunks = []

    with zipfile.ZipFile(MENDELEY_ZIP, "r") as z:
        with z.open(counts_member) as f:
            reader = pd.read_csv(f, chunksize=CHUNKSIZE, low_memory=False)

            for i, chunk in enumerate(reader, start=1):
                print(f"Chunk {i}...")

                gene_col = chunk.columns[0]

                # Colonnes cellules disponibles dans ce chunk
                usable_cells = [col for col in chunk.columns if col in cell_to_sample]

                if len(usable_cells) == 0:
                    raise ValueError("Aucune colonne cellule du counts ne matche avec la metadata.")

                df_values = chunk[[gene_col] + usable_cells].copy()
                df_values[gene_col] = df_values[gene_col].astype(str)

                # Conversion numérique
                values = df_values[usable_cells].apply(pd.to_numeric, errors="coerce").fillna(0)

                # Mapping cellule -> sample
                group_labels = pd.Series(
                    [cell_to_sample[col] for col in usable_cells],
                    index=usable_cells
                )

                # Somme par sample
                summed = values.T.groupby(group_labels).sum().T
                summed.insert(0, "gene", df_values[gene_col].values)

                chunks.append(summed)

    df_gene_sample = pd.concat(chunks, axis=0, ignore_index=True)

    # Si doublons de gènes, on somme
    df_gene_sample = df_gene_sample.groupby("gene", as_index=True).sum(numeric_only=True)

    print("Dimensions pseudo-bulk gènes x samples :", df_gene_sample.shape)

    # Transpose : samples en lignes, gènes en colonnes
    df_sample_gene = df_gene_sample.T.reset_index()
    df_sample_gene = df_sample_gene.rename(columns={"index": "sample_id"})

    print("Dimensions pseudo-bulk samples x gènes :", df_sample_gene.shape)

    return df_sample_gene


def make_logcpm(df_counts_with_meta, metadata_cols):
    """
    Transforme les pseudo-bulk counts en log2(CPM + 1).
    """
    gene_cols = [col for col in df_counts_with_meta.columns if col not in metadata_cols]

    df_log = df_counts_with_meta[metadata_cols].copy()
    counts = df_counts_with_meta[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    library_sizes = counts.sum(axis=1)
    library_sizes = library_sizes.replace(0, np.nan)

    cpm = counts.div(library_sizes, axis=0) * 1_000_000
    logcpm = np.log2(cpm + 1)

    df_log = pd.concat([df_log, logcpm], axis=1)

    return df_log


# ============================================================
# MAIN
# ============================================================

def main():
    if not MENDELEY_ZIP.exists():
        raise FileNotFoundError(f"Zip Mendeley introuvable : {MENDELEY_ZIP}")

    print("Lecture header RNA counts...")
    counts_member, columns = get_counts_header()

    print("Fichier counts utilisé :", counts_member)
    print("Nombre total de colonnes counts :", len(columns))

    gene_col = columns[0]
    count_cell_cols = columns[1:]

    print("Colonne gène :", gene_col)
    print("Nombre de cellules dans counts :", len(count_cell_cols))

    print("\nLecture metadata RNA...")
    metadata = read_metadata()
    print("Dimensions metadata :", metadata.shape)

    cell_meta, sample_meta = build_cell_metadata(metadata, count_cell_cols)

    sample_meta.to_csv(OUT_CELL_COUNTS, index=False)
    print("\n✅ Table sample/cell counts créée :", OUT_CELL_COUNTS)

    df_pseudobulk = compute_pseudobulk_counts(counts_member, cell_meta)

    # Fusion sample metadata + counts
    df_counts_ready = sample_meta.merge(df_pseudobulk, on="sample_id", how="inner")

    print("\nDimensions finale pseudo-bulk counts :", df_counts_ready.shape)

    print("\nRépartition response :")
    print(df_counts_ready["response"].value_counts(dropna=False))

    print("\nRépartition response_binary :")
    print(df_counts_ready["response_binary"].value_counts(dropna=False))

    print("\nNombre de cellules par sample :")
    print(df_counts_ready[["sample_id", "n_cells", "response", "response_binary"]])

    df_counts_ready.to_csv(OUT_COUNTS, index=False)
    print("\n✅ Pseudo-bulk counts créé :", OUT_COUNTS)

    metadata_cols = [
        "study",
        "sample_id",
        "response",
        "response_binary",
        "n_cells",
        "expression_format",
        "dataset_role",
    ]

    df_logcpm = make_logcpm(df_counts_ready, metadata_cols)
    df_logcpm["expression_format"] = "pseudobulk log2(CPM + 1) from single-cell counts"

    df_logcpm.to_csv(OUT_LOGCPM, index=False)
    print("✅ Pseudo-bulk logCPM créé :", OUT_LOGCPM)

    # Résumé
    n_gene_cols = df_counts_ready.shape[1] - len(metadata_cols)

    summary = pd.DataFrame([
        {
            "dataset": "MENDELEY_MELANOMA_CITESEQ",
            "source": "Mendeley",
            "n_cells_used": cell_meta.shape[0],
            "n_samples": df_counts_ready.shape[0],
            "n_genes": n_gene_cols,
            "n_response_0": int((df_counts_ready["response_binary"] == 0).sum()),
            "n_response_1": int((df_counts_ready["response_binary"] == 1).sum()),
            "n_response_missing": int(df_counts_ready["response_binary"].isna().sum()),
            "counts_output": str(OUT_COUNTS),
            "logcpm_output": str(OUT_LOGCPM),
        }
    ])

    summary.to_csv(SUMMARY_FILE, index=False)

    print("\n" + "=" * 100)
    print("RÉSUMÉ GLOBAL SCRIPT 17")
    print("=" * 100)
    print(summary)

    print("\n✅ Résumé sauvegardé :", SUMMARY_FILE)


if __name__ == "__main__":
    main()