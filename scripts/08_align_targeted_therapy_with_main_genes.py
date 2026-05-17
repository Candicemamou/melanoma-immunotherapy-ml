from pathlib import Path
import pandas as pd

print("--- 08 : ALIGNEMENT GSE77940 AVEC LES GÈNES DU DATASET PRINCIPAL ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_PROCESSED = ROOT / "data_processed"

MAIN_ML_FILE = DATA_PROCESSED / "3_DATASET_RNA_MULTICOHORTE_ML_READY.csv"
TARGETED_FILE = DATA_PROCESSED / "5_TARGETED_THERAPY_GSE77940_clean.csv"

OUTPUT_FILE = DATA_PROCESSED / "5_TARGETED_THERAPY_GSE77940_ALIGNED_WITH_MAIN_GENES.csv"
SUMMARY_FILE = DATA_PROCESSED / "5_TARGETED_THERAPY_GSE77940_alignment_summary.txt"


# ============================================================
# MAIN
# ============================================================

def main():
    if not MAIN_ML_FILE.exists():
        raise FileNotFoundError(f"Fichier introuvable : {MAIN_ML_FILE}")

    if not TARGETED_FILE.exists():
        raise FileNotFoundError(f"Fichier introuvable : {TARGETED_FILE}")

    print("Lecture du dataset principal ML ready...")
    df_main = pd.read_csv(MAIN_ML_FILE, nrows=5)

    print("Lecture du dataset GSE77940 nettoyé...")
    df_targeted = pd.read_csv(TARGETED_FILE)

    # Colonnes metadata du dataset principal
    main_metadata_cols = {
        "patient_id",
        "study",
        "response",
        "survival_event"
    }

    # Colonnes metadata du dataset GSE77940
    targeted_metadata_cols = {
        "sample_id",
        "patient_id",
        "treatment_phase",
        "therapy_type",
        "study",
        "treatment_category"
    }

    # Gènes du dataset principal
    main_genes = [
        col for col in df_main.columns
        if col not in main_metadata_cols
    ]

    # Gènes du dataset targeted
    targeted_genes = [
        col for col in df_targeted.columns
        if col not in targeted_metadata_cols
    ]

    main_genes_set = set(main_genes)
    targeted_genes_set = set(targeted_genes)

    common_genes = sorted(list(main_genes_set & targeted_genes_set))

    print("Nombre de gènes dataset principal :", len(main_genes))
    print("Nombre de gènes GSE77940 :", len(targeted_genes))
    print("Nombre de gènes communs :", len(common_genes))

    if len(common_genes) == 0:
        raise ValueError("Aucun gène commun trouvé. Vérifie les noms de colonnes.")

    # Colonnes à garder
    metadata_cols_order = [
        "sample_id",
        "patient_id",
        "study",
        "treatment_category",
        "therapy_type",
        "treatment_phase"
    ]

    metadata_cols_order = [
        col for col in metadata_cols_order
        if col in df_targeted.columns
    ]

    df_aligned = df_targeted[metadata_cols_order + common_genes].copy()

    # Sauvegarde
    df_aligned.to_csv(OUTPUT_FILE, index=False)

    print("\n✅ Fichier aligné créé :")
    print(OUTPUT_FILE)

    print("\nDimensions fichier aligné :", df_aligned.shape)

    print("\nRépartition treatment_phase :")
    print(df_aligned["treatment_phase"].value_counts())

    print("\nRépartition therapy_type :")
    print(df_aligned["therapy_type"].value_counts())

    # Rapport texte
    with open(SUMMARY_FILE, "w", encoding="utf-8") as f:
        f.write("ALIGNEMENT GSE77940 AVEC DATASET PRINCIPAL\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Dataset principal : {MAIN_ML_FILE}\n")
        f.write(f"Dataset ciblé : {TARGETED_FILE}\n\n")
        f.write(f"Nombre de gènes dataset principal : {len(main_genes)}\n")
        f.write(f"Nombre de gènes GSE77940 : {len(targeted_genes)}\n")
        f.write(f"Nombre de gènes communs : {len(common_genes)}\n\n")
        f.write(f"Dimensions fichier aligné : {df_aligned.shape}\n\n")
        f.write("Colonnes metadata conservées :\n")
        for col in metadata_cols_order:
            f.write(f"- {col}\n")

    print("\n✅ Rapport créé :")
    print(SUMMARY_FILE)

    print("\nAperçu :")
    print(df_aligned[metadata_cols_order].head())
    print(df_aligned.shape)


if __name__ == "__main__":
    main()