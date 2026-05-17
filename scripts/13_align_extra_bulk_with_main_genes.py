from pathlib import Path
import pandas as pd

print("--- 13 : ALIGNEMENT DES DATASETS BULK SUPPLÉMENTAIRES AVEC LE DATASET PRINCIPAL ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_PROCESSED = ROOT / "data_processed"

MAIN_FILE = DATA_PROCESSED / "3_DATASET_RNA_MULTICOHORTE_ML_READY.csv"

GSE168204_FILE = DATA_PROCESSED / "12_GSE168204_bulk_clean.csv"
GSE244982_FILE = DATA_PROCESSED / "12_GSE244982_bulk_clean.csv"
GSE196434_FILE = DATA_PROCESSED / "12_GSE196434_bulk_clean.csv"

OUT_GSE168204 = DATA_PROCESSED / "13_GSE168204_EXTERNAL_VALIDATION_READY.csv"
OUT_GSE244982 = DATA_PROCESSED / "13_GSE244982_RESISTANCE_READY.csv"
OUT_GSE196434 = DATA_PROCESSED / "13_GSE196434_TARGETED_THERAPY_READY.csv"

SUMMARY_FILE = DATA_PROCESSED / "13_EXTRA_BULK_ALIGNMENT_summary.csv"


# ============================================================
# CONFIGURATION
# ============================================================

MAIN_METADATA_COLS = {
    "patient_id",
    "study",
    "response",
    "survival_event"
}

EXTRA_METADATA_COLS = [
    "study",
    "platform",
    "sample_id",
    "geo_accession",
    "patient_id",
    "treatment",
    "timepoint",
    "response",
    "response_binary",
    "dataset_role",
    "notes",
    "sample_id_merge"
]


# ============================================================
# OUTILS
# ============================================================

def get_main_genes():
    """
    Récupère les gènes du dataset principal ML_READY.
    """
    if not MAIN_FILE.exists():
        raise FileNotFoundError(f"Fichier principal introuvable : {MAIN_FILE}")

    df_main_head = pd.read_csv(MAIN_FILE, nrows=5)

    main_genes = [
        col for col in df_main_head.columns
        if col not in MAIN_METADATA_COLS
    ]

    return main_genes


def align_dataset(input_file, output_file, dataset_name, main_genes):
    """
    Aligne un dataset supplémentaire sur les gènes du dataset principal.
    """

    print("\n" + "=" * 80)
    print(f"ALIGNEMENT : {dataset_name}")
    print("=" * 80)

    if not input_file.exists():
        raise FileNotFoundError(f"Fichier introuvable : {input_file}")

    print("Lecture :", input_file)
    df = pd.read_csv(input_file)

    print("Dimensions avant alignement :", df.shape)

    # Colonnes metadata présentes
    metadata_cols_present = [
        col for col in EXTRA_METADATA_COLS
        if col in df.columns
    ]

    # Gènes présents dans le dataset supplémentaire
    extra_genes = [
        col for col in df.columns
        if col not in metadata_cols_present
    ]

    # Gènes communs avec le dataset principal
    common_genes = [
        gene for gene in main_genes
        if gene in set(extra_genes)
    ]

    print("Nombre de gènes dataset principal :", len(main_genes))
    print("Nombre de gènes dataset supplémentaire :", len(extra_genes))
    print("Nombre de gènes communs :", len(common_genes))

    if len(common_genes) == 0:
        raise ValueError(f"Aucun gène commun trouvé pour {dataset_name}.")

    # Garder metadata + gènes communs dans l'ordre du dataset principal
    df_aligned = df[metadata_cols_present + common_genes].copy()

    # Nettoyage response_binary
    if "response_binary" in df_aligned.columns:
        df_aligned["response_binary"] = pd.to_numeric(
            df_aligned["response_binary"],
            errors="coerce"
        )

    # Sauvegarde
    df_aligned.to_csv(output_file, index=False)

    print("✅ Fichier aligné créé :", output_file)
    print("Dimensions après alignement :", df_aligned.shape)

    if "response_binary" in df_aligned.columns:
        print("\nRépartition response_binary :")
        print(df_aligned["response_binary"].value_counts(dropna=False))

    if "timepoint" in df_aligned.columns:
        print("\nRépartition timepoint :")
        print(df_aligned["timepoint"].value_counts(dropna=False))

    if "treatment" in df_aligned.columns:
        print("\nRépartition treatment :")
        print(df_aligned["treatment"].value_counts(dropna=False))

    print("\nAperçu metadata :")
    preview_cols = [
        col for col in [
            "study",
            "sample_id",
            "patient_id",
            "treatment",
            "timepoint",
            "response",
            "response_binary",
            "dataset_role"
        ]
        if col in df_aligned.columns
    ]
    print(df_aligned[preview_cols].head())

    return {
        "dataset": dataset_name,
        "input_file": str(input_file),
        "output_file": str(output_file),
        "n_samples": df_aligned.shape[0],
        "n_cols_total": df_aligned.shape[1],
        "n_metadata_cols": len(metadata_cols_present),
        "n_common_genes": len(common_genes),
        "n_extra_genes_before_alignment": len(extra_genes),
        "n_main_genes": len(main_genes)
    }


# ============================================================
# MAIN
# ============================================================

def main():
    print("Lecture des gènes du dataset principal...")
    main_genes = get_main_genes()

    print("Nombre de gènes principaux :", len(main_genes))

    summary_rows = []

    summary_rows.append(
        align_dataset(
            input_file=GSE168204_FILE,
            output_file=OUT_GSE168204,
            dataset_name="GSE168204_EXTERNAL_VALIDATION",
            main_genes=main_genes
        )
    )

    summary_rows.append(
        align_dataset(
            input_file=GSE244982_FILE,
            output_file=OUT_GSE244982,
            dataset_name="GSE244982_RESISTANCE",
            main_genes=main_genes
        )
    )

    summary_rows.append(
        align_dataset(
            input_file=GSE196434_FILE,
            output_file=OUT_GSE196434,
            dataset_name="GSE196434_TARGETED_THERAPY",
            main_genes=main_genes
        )
    )

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(SUMMARY_FILE, index=False)

    print("\n" + "=" * 80)
    print("RÉSUMÉ GLOBAL SCRIPT 13")
    print("=" * 80)
    print(df_summary)

    print("\n✅ Résumé sauvegardé :", SUMMARY_FILE)


if __name__ == "__main__":
    main()