from pathlib import Path
import re
import pandas as pd

print("--- 07 : NETTOYAGE GSE77940 — THÉRAPIE CIBLÉE PRE/POST ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"

INPUT_FILE = DATA_RAW / "5_TARGETED_THERAPY_GSE77940_pre_post_raw.xlsx"
OUTPUT_FILE = DATA_PROCESSED / "5_TARGETED_THERAPY_GSE77940_clean.csv"
SUMMARY_FILE = DATA_PROCESSED / "5_TARGETED_THERAPY_GSE77940_summary.csv"

DATA_PROCESSED.mkdir(exist_ok=True)


# ============================================================
# FONCTIONS
# ============================================================

def parse_sample_column(col_name):
    """
    Extrait les informations depuis un nom de colonne du type :
    'Patient 1 Pre-RAF/MEKi'
    'Patient 3 Post-RAFi'

    Retourne :
    patient_id, treatment_phase, therapy_type
    """
    col = str(col_name).strip()

    match = re.match(r"Patient\s+(\d+)\s+(Pre|Post)-(.+)", col)

    if match is None:
        return None

    patient_number = match.group(1)
    phase = match.group(2)
    therapy_raw = match.group(3)

    patient_id = f"GSE77940_Patient_{patient_number}"

    if "MEK" in therapy_raw.upper():
        therapy_type = "RAF_MEKi"
    else:
        therapy_type = "RAFi"

    return {
        "sample_id": col,
        "patient_id": patient_id,
        "treatment_phase": phase,
        "therapy_type": therapy_type,
        "study": "GSE77940",
        "treatment_category": "targeted_therapy"
    }


def clean_gene_symbol(value):
    """Nettoie les symboles de gènes."""
    if pd.isna(value):
        return None

    value = str(value).strip()

    if value == "" or value.lower() == "nan":
        return None

    return value


# ============================================================
# MAIN
# ============================================================

def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Fichier introuvable : {INPUT_FILE}")

    print(f"Lecture du fichier : {INPUT_FILE}")

    # Le fichier a une seule feuille utile
    df = pd.read_excel(INPUT_FILE)

    print("Dimensions brutes :", df.shape)
    print("Colonnes :", df.columns.tolist())

    # Vérification des colonnes essentielles
    required_cols = ["Gene Symbol"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Colonne obligatoire absente : {col}")

    # Nettoyer les symboles de gènes
    df["Gene Symbol"] = df["Gene Symbol"].apply(clean_gene_symbol)

    # Supprimer les lignes sans symbole de gène
    before = len(df)
    df = df.dropna(subset=["Gene Symbol"]).copy()
    after = len(df)

    print(f"Lignes sans Gene Symbol supprimées : {before - after}")

    # Identifier les colonnes échantillons
    sample_cols = []
    sample_metadata = []

    for col in df.columns:
        parsed = parse_sample_column(col)
        if parsed is not None:
            sample_cols.append(col)
            sample_metadata.append(parsed)

    print("Nombre d'échantillons trouvés :", len(sample_cols))

    if len(sample_cols) == 0:
        raise ValueError("Aucune colonne échantillon Pre/Post trouvée.")

    # Garder uniquement Gene Symbol + échantillons
    df_expr = df[["Gene Symbol"] + sample_cols].copy()

    # Gérer les doublons de symboles de gènes
    # Plusieurs Ensembl IDs peuvent correspondre au même symbole.
    # On moyenne les valeurs pour obtenir une seule colonne par gène.
    df_expr = df_expr.groupby("Gene Symbol", as_index=True).mean(numeric_only=True)

    print("Nombre de gènes après moyenne des doublons :", df_expr.shape[0])

    # Transposer : échantillons en lignes, gènes en colonnes
    df_t = df_expr.T.reset_index()
    df_t = df_t.rename(columns={"index": "sample_id"})

    # Ajouter les métadonnées échantillons
    df_meta = pd.DataFrame(sample_metadata)

    df_clean = df_meta.merge(df_t, on="sample_id", how="inner")

    print("Dimensions finales :", df_clean.shape)

    # Sauvegarde du fichier propre
    df_clean.to_csv(OUTPUT_FILE, index=False)

    print(f"✅ Fichier nettoyé créé : {OUTPUT_FILE}")

    # Résumé par phase / traitement
    summary = (
        df_clean
        .groupby(["study", "treatment_category", "therapy_type", "treatment_phase"])
        .size()
        .reset_index(name="n_samples")
    )

    summary.to_csv(SUMMARY_FILE, index=False)

    print(f"✅ Résumé créé : {SUMMARY_FILE}")

    print("\nRésumé :")
    print(summary)

    print("\nAperçu du fichier propre :")
    print(df_clean[["sample_id", "patient_id", "treatment_phase", "therapy_type", "study"]].head())

    print("\nNombre de patients uniques :", df_clean["patient_id"].nunique())
    print("Nombre d'échantillons :", df_clean.shape[0])
    print("Nombre de gènes :", df_clean.shape[1] - 5)


if __name__ == "__main__":
    main()