from pathlib import Path
import pandas as pd

print("--- 16 : FUSION DES DATASETS BULK CBIOPORTAL ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_PROCESSED = ROOT / "data_processed"

GIDE_FILE = DATA_PROCESSED / "15_CBIO_GIDE_2019_clean.csv"
UCLA_FILE = DATA_PROCESSED / "15_CBIO_UCLA_2016_clean.csv"
LIU_FILE = DATA_PROCESSED / "15_CBIO_LIU_2019_clean.csv"
DFCI_FILE = DATA_PROCESSED / "15_CBIO_DFCI_2019_clean.csv"

OUT_RECIST = DATA_PROCESSED / "16_CBIOPORTAL_BULK_RECIST_ML_READY.csv"
OUT_UCLA = DATA_PROCESSED / "16_CBIOPORTAL_UCLA_2016_ENTREZ_READY.csv"
OUT_DFCI = DATA_PROCESSED / "16_CBIOPORTAL_DFCI_PFS_READY.csv"
SUMMARY_FILE = DATA_PROCESSED / "16_CBIOPORTAL_MERGED_summary.csv"


# ============================================================
# CONFIGURATION
# ============================================================

RECIST_DATASETS = {
    "CBIO_GIDE_2019": GIDE_FILE,
    "CBIO_LIU_2019": LIU_FILE,
}
UCLA_DATASET = {
    "CBIO_UCLA_2016": UCLA_FILE,
}

DFCI_DATASET = {
    "CBIO_DFCI_2019": DFCI_FILE,
}

# Colonnes metadata possibles à conserver.
# Tout ce qui n'est pas dans cette liste sera considéré comme un gène.
METADATA_COLS = [
    "study",
    "sample_id",
    "patient_id",
    "response",
    "response_binary",
    "treatment",
    "expression_format",
    "dataset_role",

    "biopsy_site",
    "biopsy_time",
    "sample_treatment",
    "treatment_response",
    "durable_clinical_benefit",
    "clinical_benefit",
    "progression",
    "pfs_status",
    "pfs_months",
    "os_status",
    "os_months",
    "sex",
    "age_at_diagnosis",
    "m_stage",
    "prior_ici_rx",
    "prior_rx",
    "prior_ctla4",
    "prior_mapk_tx",
    "tmb_nonsynonymous",
    "mutation_load",
]


# ============================================================
# OUTILS
# ============================================================

def read_clean_dataset(path, dataset_name):
    if not path.exists():
        raise FileNotFoundError(f"Fichier introuvable pour {dataset_name} : {path}")

    print("\n" + "=" * 90)
    print(f"LECTURE : {dataset_name}")
    print("=" * 90)

    df = pd.read_csv(path, low_memory=False)

    print("Dimensions :", df.shape)

    if "study" not in df.columns:
        df["study"] = dataset_name

    df["study"] = df["study"].fillna(dataset_name)

    if "response_binary" in df.columns:
        df["response_binary"] = pd.to_numeric(df["response_binary"], errors="coerce")

    print("Répartition response :")
    if "response" in df.columns:
        print(df["response"].value_counts(dropna=False))
    else:
        print("Colonne response absente")

    print("Répartition response_binary :")
    if "response_binary" in df.columns:
        print(df["response_binary"].value_counts(dropna=False))
    else:
        print("Colonne response_binary absente")

    return df


def get_gene_columns(df):
    metadata_present = [col for col in METADATA_COLS if col in df.columns]
    gene_cols = [col for col in df.columns if col not in metadata_present]
    return gene_cols


def align_and_merge_datasets(datasets_dict, output_file, merged_name, drop_missing_response=True):
    """
    Fusionne plusieurs datasets en gardant uniquement les gènes communs.
    """

    print("\n" + "#" * 100)
    print(f"FUSION : {merged_name}")
    print("#" * 100)

    loaded = {}
    gene_sets = {}

    for dataset_name, path in datasets_dict.items():
        df = read_clean_dataset(path, dataset_name)

        if drop_missing_response:
            before = df.shape[0]
            df = df.dropna(subset=["response_binary"]).copy()
            after = df.shape[0]
            print(f"Lignes gardées après suppression response_binary NaN : {after}/{before}")

        loaded[dataset_name] = df
        gene_cols = get_gene_columns(df)
        gene_sets[dataset_name] = set(gene_cols)

        print(f"Nombre de gènes détectés pour {dataset_name} : {len(gene_cols)}")

    # Gènes communs
    common_genes = sorted(set.intersection(*gene_sets.values()))

    print("\nNombre de gènes communs :", len(common_genes))

    if len(common_genes) == 0:
        raise ValueError("Aucun gène commun trouvé entre les datasets.")

    merged_parts = []
    summary_rows = []

    for dataset_name, df in loaded.items():
        metadata_present = [col for col in METADATA_COLS if col in df.columns]

        # Garder metadata + gènes communs
        df_aligned = df[metadata_present + common_genes].copy()

        # S'assurer que les colonnes essentielles sont devant
        preferred_order = [
            "study",
            "sample_id",
            "patient_id",
            "response",
            "response_binary",
            "treatment",
            "expression_format",
            "dataset_role",
        ]

        ordered_metadata = [col for col in preferred_order if col in df_aligned.columns]
        other_metadata = [
            col for col in metadata_present
            if col not in ordered_metadata and col in df_aligned.columns
        ]

        df_aligned = df_aligned[ordered_metadata + other_metadata + common_genes]

        merged_parts.append(df_aligned)

        summary_rows.append({
            "merged_dataset": merged_name,
            "source_dataset": dataset_name,
            "n_samples": df_aligned.shape[0],
            "n_total_columns": df_aligned.shape[1],
            "n_common_genes": len(common_genes),
            "n_response_0": int((df_aligned["response_binary"] == 0).sum()) if "response_binary" in df_aligned.columns else None,
            "n_response_1": int((df_aligned["response_binary"] == 1).sum()) if "response_binary" in df_aligned.columns else None,
            "n_response_missing": int(df_aligned["response_binary"].isna().sum()) if "response_binary" in df_aligned.columns else None,
            "expression_formats": "; ".join(sorted(df_aligned["expression_format"].dropna().astype(str).unique())) if "expression_format" in df_aligned.columns else None,
        })

    df_merged = pd.concat(merged_parts, axis=0, ignore_index=True)

    print("\nDimensions fichier fusionné :", df_merged.shape)

    if "study" in df_merged.columns:
        print("\nRépartition par study :")
        print(df_merged["study"].value_counts(dropna=False))

    if "response" in df_merged.columns:
        print("\nRépartition response :")
        print(df_merged["response"].value_counts(dropna=False))

    if "response_binary" in df_merged.columns:
        print("\nRépartition response_binary :")
        print(df_merged["response_binary"].value_counts(dropna=False))

    if "treatment" in df_merged.columns:
        print("\nRépartition treatment :")
        print(df_merged["treatment"].value_counts(dropna=False))

    df_merged.to_csv(output_file, index=False)

    print("\n✅ Fichier fusionné créé :", output_file)

    return df_merged, summary_rows


# ============================================================
# MAIN
# ============================================================

def main():
    all_summary_rows = []

    # ------------------------------------------------------------
    # 1. Fusion principale RECIST / réponse clinique directe
    # ------------------------------------------------------------

    df_recist, recist_summary = align_and_merge_datasets(
        datasets_dict=RECIST_DATASETS,
        output_file=OUT_RECIST,
        merged_name="CBIOPORTAL_BULK_RECIST_ML_READY",
        drop_missing_response=True
    )

    all_summary_rows.extend(recist_summary)

    # ------------------------------------------------------------
    # 2. DFCI séparé : endpoint progression / censure
    # ------------------------------------------------------------

    df_dfci, dfci_summary = align_and_merge_datasets(
        datasets_dict=DFCI_DATASET,
        output_file=OUT_DFCI,
        merged_name="CBIOPORTAL_DFCI_PFS_READY",
        drop_missing_response=True
    )

    all_summary_rows.extend(dfci_summary)

    # ------------------------------------------------------------
    # 3. UCLA séparé : gènes en Entrez IDs, non fusionné avec Hugo_Symbol
    # ------------------------------------------------------------

    df_ucla, ucla_summary = align_and_merge_datasets(
        datasets_dict=UCLA_DATASET,
        output_file=OUT_UCLA,
        merged_name="CBIOPORTAL_UCLA_2016_ENTREZ_READY",
        drop_missing_response=True
    )

    all_summary_rows.extend(ucla_summary)

    # ------------------------------------------------------------
    # Résumé global
    # ------------------------------------------------------------

    df_summary = pd.DataFrame(all_summary_rows)
    df_summary.to_csv(SUMMARY_FILE, index=False)

    print("\n" + "=" * 100)
    print("RÉSUMÉ GLOBAL SCRIPT 16")
    print("=" * 100)
    print(df_summary)

    print("\n✅ Résumé sauvegardé :", SUMMARY_FILE)

    print("\nFichiers créés :")
    print("-", OUT_RECIST)
    print("-", OUT_UCLA)
    print("-", OUT_DFCI)


if __name__ == "__main__":
    main()