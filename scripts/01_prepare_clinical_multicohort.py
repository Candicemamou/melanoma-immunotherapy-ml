from pathlib import Path
import pandas as pd

print("--- 01 : PRÉPARATION CLINIQUE MULTICOHORTE ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"
ARCHIVE_DATA = ROOT / "archive_data"

DATA_PROCESSED.mkdir(exist_ok=True)

OUTPUT_FILE = DATA_PROCESSED / "1_CLINIQUE_MULTICOHORTE_nettoyee.csv"

# Fichiers attendus
GSE78220_METADATA = DATA_RAW / "2_BULK_RNA_metadata_GSE78220_brut.txt"

# GSE91061 peut être à plusieurs endroits selon ton rangement
GSE91061_CANDIDATES = [
    DATA_RAW / "1_CLINIQUE_GSE91061_nettoyee.csv",
    ARCHIVE_DATA / "1_CLINIQUE_GSE91061_nettoyee.csv",
    DATA_PROCESSED / "1_CLINIQUE_GSE91061_nettoyee.csv",
]


# ============================================================
# FONCTIONS
# ============================================================

def find_existing_file(candidates):
    """Retourne le premier fichier qui existe dans une liste de chemins."""
    for path in candidates:
        if path.exists():
            return path
    return None


def clean_text(value):
    """Nettoyage simple des textes."""
    if pd.isna(value):
        return value
    return str(value).replace('"', "").strip()


def map_response_to_binary(response):
    """
    Convertit la réponse clinique en label binaire.
    Convention utilisée ici :
    0 = non-répondeur / progression
    1 = répondeur ou maladie stable
    """
    if pd.isna(response):
        return None

    r = str(response).strip().upper()

    mapping = {
        "PD": 0,
        "PROGRESSIVE DISEASE": 0,

        "SD": 1,
        "STABLE DISEASE": 1,
        "PR": 1,
        "PARTIAL RESPONSE": 1,
        "CR": 1,
        "COMPLETE RESPONSE": 1,
    }

    return mapping.get(r, None)


def load_gse91061_clinical():
    """
    Charge le fichier clinique déjà nettoyé de GSE91061.
    Ce fichier doit contenir au minimum :
    patient_id, response ou survival_event.
    """
    path = find_existing_file(GSE91061_CANDIDATES)

    if path is None:
        raise FileNotFoundError(
            "Impossible de trouver 1_CLINIQUE_GSE91061_nettoyee.csv.\n"
            "Place ce fichier dans data_raw/ ou archive_data/."
        )

    print(f"Lecture clinique GSE91061 : {path}")
    df = pd.read_csv(path)

    # Nettoyage noms colonnes
    df.columns = df.columns.astype(str).str.strip()

    if "patient_id" not in df.columns:
        raise ValueError("Le fichier GSE91061 doit contenir une colonne patient_id.")

    df["patient_id"] = df["patient_id"].astype(str).str.strip().str.replace('"', "", regex=False)

    # Si response existe, on la garde proprement
    if "response" in df.columns:
        df["response"] = df["response"].apply(clean_text)
    else:
        df["response"] = None

    # Si survival_event existe, on le garde
    if "survival_event" in df.columns:
        df["survival_event"] = pd.to_numeric(df["survival_event"], errors="coerce")
    else:
        # Sinon on essaie de créer depuis response
        df["survival_event"] = df["response"].apply(map_response_to_binary)

    df["response_binary"] = df["survival_event"]
    df["study"] = "GSE91061"

    # On garde uniquement les colonnes essentielles pour le pipeline actuel
    df_final = df[["patient_id", "response", "survival_event", "response_binary", "study"]].copy()

    print(f"GSE91061 : {df_final.shape[0]} lignes")
    return df_final


def load_gse78220_clinical():
    """
    Extrait les métadonnées cliniques principales de GSE78220 depuis le fichier GEO.
    """
    if not GSE78220_METADATA.exists():
        raise FileNotFoundError(f"Fichier introuvable : {GSE78220_METADATA}")

    print(f"Lecture metadata GSE78220 : {GSE78220_METADATA}")

    patient_ids = []
    responses = []
    ages = []
    genders = []
    disease_status = []
    overall_survival_days = []
    vital_status = []
    previous_mapki = []
    biopsy_time = []
    treatment = []

    with open(GSE78220_METADATA, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_clean = line.strip()

            if line_clean.startswith("!Sample_title"):
                values = line_clean.split("\t")[1:]
                patient_ids = [clean_text(v) for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "anti-pd-1 response" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                responses = [clean_text(v).split(":")[-1].strip() for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "age (yrs)" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                ages = [clean_text(v).split(":")[-1].strip() for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "gender" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                genders = [clean_text(v).split(":")[-1].strip() for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "disease status" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                disease_status = [clean_text(v).split(":")[-1].strip() for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "overall survival (days)" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                overall_survival_days = [clean_text(v).split(":")[-1].strip() for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "vital status" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                vital_status = [clean_text(v).split(":")[-1].strip() for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "previous mapki" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                previous_mapki = [clean_text(v).split(":")[-1].strip() for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "biopsy time" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                biopsy_time = [clean_text(v).split(":")[-1].strip() for v in values]

            elif line_clean.startswith("!Sample_characteristics_ch1") and "treatment:" in line_clean.lower():
                values = line_clean.split("\t")[1:]
                treatment = [clean_text(v).split(":")[-1].strip() for v in values]

    if len(patient_ids) == 0:
        raise ValueError("Aucun patient_id trouvé dans le fichier GSE78220.")

    if len(responses) == 0:
        raise ValueError("Aucune réponse anti-PD-1 trouvée dans le fichier GSE78220.")

    if len(patient_ids) != len(responses):
        raise ValueError(
            f"Nombre patient_id différent du nombre response : "
            f"{len(patient_ids)} vs {len(responses)}"
        )

    df = pd.DataFrame({
        "patient_id": patient_ids,
        "response": responses,
    })

    # Colonnes optionnelles : seulement si elles ont la bonne longueur
    optional_columns = {
        "age": ages,
        "gender": genders,
        "disease_status": disease_status,
        "overall_survival_days": overall_survival_days,
        "vital_status": vital_status,
        "previous_mapki": previous_mapki,
        "biopsy_time": biopsy_time,
        "treatment": treatment,
    }

    for col, values in optional_columns.items():
        if len(values) == len(df):
            df[col] = values
        else:
            df[col] = None

    # Nettoyage IDs patients
    df["patient_id"] = (
        df["patient_id"]
        .astype(str)
        .str.replace(".baseline", "", regex=False)
        .str.replace(".OnTx", "", regex=False)
        .str.strip()
    )

    # Création label binaire
    df["survival_event"] = df["response"].apply(map_response_to_binary)
    df["response_binary"] = df["survival_event"]

    df["study"] = "GSE78220"

    print(f"GSE78220 : {df.shape[0]} lignes")
    return df


# ============================================================
# MAIN
# ============================================================

def main():
    df_91061 = load_gse91061_clinical()
    df_78220 = load_gse78220_clinical()

    print("\nFusion clinique GSE91061 + GSE78220...")
    df_final = pd.concat([df_91061, df_78220], ignore_index=True)

    # Nettoyage global
    df_final["patient_id"] = df_final["patient_id"].astype(str).str.strip()
    df_final["study"] = df_final["study"].astype(str).str.strip()

    # Suppression doublons exacts patient + study
    before = len(df_final)
    df_final = df_final.drop_duplicates(subset=["patient_id", "study"], keep="first")
    after = len(df_final)

    print(f"Doublons supprimés : {before - after}")

    # Conversion numérique du label
    df_final["survival_event"] = pd.to_numeric(df_final["survival_event"], errors="coerce")
    df_final["response_binary"] = pd.to_numeric(df_final["response_binary"], errors="coerce")

    # Sauvegarde
    df_final.to_csv(OUTPUT_FILE, index=False)

    print("\n✅ Fichier clinique multicohorte créé :")
    print(OUTPUT_FILE)

    print("\nDimensions finales :", df_final.shape)

    print("\nRépartition par cohorte :")
    print(df_final["study"].value_counts())

    print("\nRépartition survival_event :")
    print(df_final["survival_event"].value_counts(dropna=False))

    print("\nAperçu :")
    print(df_final.head())


if __name__ == "__main__":
    main()