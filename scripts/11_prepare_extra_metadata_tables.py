from pathlib import Path
import re
import pandas as pd

print("--- 11 : STANDARDISATION DES MÉTADONNÉES SUPPLÉMENTAIRES ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_PROCESSED = ROOT / "data_processed"

FILES = {
    "GSE168204_GPL11154": DATA_PROCESSED / "10_METADATA_GSE168204_GPL11154.csv",
    "GSE168204_GPL18573": DATA_PROCESSED / "10_METADATA_GSE168204_GPL18573.csv",
    "GSE196434": DATA_PROCESSED / "10_METADATA_GSE196434.csv",
    "GSE244982": DATA_PROCESSED / "10_METADATA_GSE244982.csv",
    "GSE120575": DATA_PROCESSED / "10_METADATA_GSE120575.csv",
}

OUTPUT_GSE168204 = DATA_PROCESSED / "11_METADATA_GSE168204_clean.csv"
OUTPUT_GSE196434 = DATA_PROCESSED / "11_METADATA_GSE196434_clean.csv"
OUTPUT_GSE244982 = DATA_PROCESSED / "11_METADATA_GSE244982_clean.csv"
OUTPUT_GSE120575 = DATA_PROCESSED / "11_METADATA_GSE120575_clean.csv"

OUTPUT_ALL = DATA_PROCESSED / "11_METADATA_EXTRA_ALL_CLEAN.csv"
SUMMARY_FILE = DATA_PROCESSED / "11_METADATA_EXTRA_summary.csv"


# ============================================================
# OUTILS
# ============================================================

def clean_text(value):
    if pd.isna(value):
        return None
    value = str(value).replace('"', "").strip()
    return value if value != "" else None


def binary_response(value):
    """
    Harmonise les réponses binaires.
    1 = responder / long responder
    0 = non-responder / short responder / resistant
    """
    if pd.isna(value):
        return None

    v = str(value).strip().upper()

    positive = {
        "R",
        "RESPONDER",
        "RESPONSE",
        "LR",
        "LONG RESPONDER",
        "LONG_RESPONSE",
        "LONG-RESPONDER",
    }

    negative = {
        "NR",
        "NON-RESPONDER",
        "NON_RESPONDER",
        "NONRESPONDER",
        "NON RESPONSE",
        "NON-RESPONSE",
        "SR",
        "SHORT RESPONDER",
        "SHORT_RESPONSE",
        "SHORT-RESPONDER",
        "RESISTANT",
    }

    if v in positive:
        return 1

    if v in negative:
        return 0

    return None


def extract_patient_from_description(value):
    """
    Exemple :
    'MGH Patient 200' -> MGH200
    'MGH Patient PDL002' -> MGHPDL002
    """
    value = clean_text(value)
    if value is None:
        return None

    value = value.replace("MGH Patient", "").strip()
    value = value.replace(" ", "")

    if not value.startswith("MGH"):
        value = "MGH" + value

    return value


def extract_patient_from_gse196434_title(title):
    """
    Exemple :
    'EX99 [pt5]' -> pt5
    'EY00 [pt45]' -> pt45
    """
    title = clean_text(title)
    if title is None:
        return None

    match = re.search(r"\[(.*?)\]", title)

    if match:
        return match.group(1).strip()

    return None


def extract_sample_from_gse196434_title(title):
    """
    Exemple :
    'EX99 [pt5]' -> EX99
    """
    title = clean_text(title)
    if title is None:
        return None

    return title.split("[")[0].strip()


def split_gse196434_sample_class(value):
    """
    Exemple :
    'PRE SR'  -> timepoint PRE, response SR
    'POST LR' -> timepoint POST, response LR
    """
    value = clean_text(value)

    if value is None:
        return None, None

    parts = value.split()

    if len(parts) >= 2:
        return parts[0], parts[1]

    return value, None


def parse_gse120575_patient_id(value):
    """
    Exemple :
    Pre_P1  -> patient_id P1, timepoint PRE
    Post_P1 -> patient_id P1, timepoint POST
    """
    value = clean_text(value)

    if value is None:
        return None, None

    if "_" in value:
        timepoint, patient = value.split("_", 1)
        return patient, timepoint.upper()

    return value, None


def standard_columns(df):
    """
    Garantit les colonnes finales dans le même ordre.
    """
    final_cols = [
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
    ]

    for col in final_cols:
        if col not in df.columns:
            df[col] = None

    return df[final_cols].copy()


# ============================================================
# PRÉPARATION GSE168204
# ============================================================

def prepare_gse168204():
    dfs = []

    for key in ["GSE168204_GPL11154", "GSE168204_GPL18573"]:
        path = FILES[key]

        if not path.exists():
            print(f"⚠️ Fichier absent : {path}")
            continue

        df = pd.read_csv(path)

        out = pd.DataFrame(index=df.index)
        out["study"] = "GSE168204"
        out["platform"] = df.get("platform_id")
        out["sample_id"] = df.get("title")
        out["geo_accession"] = df.get("geo_accession")
        out["patient_id"] = df.get("description").apply(extract_patient_from_description)
        out["treatment"] = df.get("characteristics_antibody")
        out["timepoint"] = df.get("characteristics_treatment_state")
        out["response"] = df.get("characteristics_response_immune_checkpoint_blockade_therapy")
        out["response_binary"] = out["response"].apply(binary_response)
        out["dataset_role"] = "extra_immunotherapy_response"
        out["notes"] = df.get("characteristics_timepoint_days")

        dfs.append(out)

    if len(dfs) == 0:
        return pd.DataFrame()

    out = pd.concat(dfs, ignore_index=True)
    out = standard_columns(out)

    out.to_csv(OUTPUT_GSE168204, index=False)
    print(f"✅ GSE168204 metadata clean créé : {OUTPUT_GSE168204}")
    print("Dimensions :", out.shape)
    print(out[["study", "sample_id", "patient_id", "treatment", "timepoint", "response", "response_binary"]].head())

    return out


# ============================================================
# PRÉPARATION GSE196434
# ============================================================

def prepare_gse196434():
    path = FILES["GSE196434"]

    if not path.exists():
        print(f"⚠️ Fichier absent : {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    sample_classes = df["characteristics_sample_class"].apply(split_gse196434_sample_class)

    out = pd.DataFrame(index=df.index)
    out["study"] = "GSE196434"
    out["platform"] = df.get("platform_id")
    out["sample_id"] = df["title"].apply(extract_sample_from_gse196434_title)
    out["geo_accession"] = df.get("geo_accession")
    out["patient_id"] = df["title"].apply(extract_patient_from_gse196434_title)
    out["treatment"] = df.get("characteristics_drug")
    out["timepoint"] = [x[0] for x in sample_classes]
    out["response"] = [x[1] for x in sample_classes]
    out["response_binary"] = out["response"].apply(binary_response)
    out["dataset_role"] = "targeted_therapy_response"
    out["notes"] = df.get("characteristics_pt_id_tm_id_according_to_table_s4,_pubmed_id_35643181")

    out = standard_columns(out)

    out.to_csv(OUTPUT_GSE196434, index=False)
    print(f"\n✅ GSE196434 metadata clean créé : {OUTPUT_GSE196434}")
    print("Dimensions :", out.shape)
    print(out[["study", "sample_id", "patient_id", "treatment", "timepoint", "response", "response_binary"]].head())

    return out


# ============================================================
# PRÉPARATION GSE244982
# ============================================================

def prepare_gse244982():
    path = FILES["GSE244982"]

    if not path.exists():
        print(f"⚠️ Fichier absent : {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    out = pd.DataFrame(index=df.index)
    out["study"] = "GSE244982"
    out["platform"] = df.get("platform_id")
    out["sample_id"] = df.get("title")
    out["geo_accession"] = df.get("geo_accession")
    out["patient_id"] = df.get("title")
    out["treatment"] = df.get("characteristics_treatment")
    out["timepoint"] = "resistant"

    # Ici la cohorte est déjà une cohorte résistante.
    # Ce n'est pas responder/non-responder classique.
    out["response"] = "resistant"
    out["response_binary"] = 0

    out["dataset_role"] = "immunotherapy_resistance"
    out["notes"] = (
        "gender=" + df.get("characteristics_gender").astype(str)
        + "; stage=" + df.get("characteristics_stage").astype(str)
        + "; prior_ctla4i=" + df.get("characteristics_prior.ctla4i").astype(str)
        + "; previous_brafi=" + df.get("characteristics_previous._brafi").astype(str)
    )

    out = standard_columns(out)

    out.to_csv(OUTPUT_GSE244982, index=False)
    print(f"\n✅ GSE244982 metadata clean créé : {OUTPUT_GSE244982}")
    print("Dimensions :", out.shape)
    print(out[["study", "sample_id", "patient_id", "treatment", "timepoint", "response", "response_binary"]].head())

    return out


# ============================================================
# PRÉPARATION GSE120575
# ============================================================

def prepare_gse120575():
    path = FILES["GSE120575"]

    if not path.exists():
        print(f"⚠️ Fichier absent : {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    parsed = df["characteristics_patient_id_pre=baseline;_post=on_treatment"].apply(parse_gse120575_patient_id)

    out = pd.DataFrame(index=df.index)
    out["study"] = "GSE120575"
    out["platform"] = df.get("platform_id")
    out["sample_id"] = df.get("title")
    out["geo_accession"] = df.get("geo_accession")
    out["patient_id"] = [x[0] for x in parsed]
    out["treatment"] = df.get("characteristics_therapy")
    out["timepoint"] = [x[1] for x in parsed]
    out["response"] = df.get("characteristics_response")
    out["response_binary"] = out["response"].apply(binary_response)
    out["dataset_role"] = "single_cell_immunotherapy_response"
    out["notes"] = "single_cell"

    out = standard_columns(out)

    out.to_csv(OUTPUT_GSE120575, index=False)
    print(f"\n✅ GSE120575 metadata clean créé : {OUTPUT_GSE120575}")
    print("Dimensions :", out.shape)
    print(out[["study", "sample_id", "patient_id", "treatment", "timepoint", "response", "response_binary"]].head())

    return out


# ============================================================
# MAIN
# ============================================================

def main():
    df_168204 = prepare_gse168204()
    df_196434 = prepare_gse196434()
    df_244982 = prepare_gse244982()
    df_120575 = prepare_gse120575()

    all_dfs = [
        df for df in [df_168204, df_196434, df_244982, df_120575]
        if not df.empty
    ]

    if len(all_dfs) == 0:
        raise ValueError("Aucun metadata clean généré.")

    df_all = pd.concat(all_dfs, ignore_index=True)
    df_all = standard_columns(df_all)

    df_all.to_csv(OUTPUT_ALL, index=False)

    print(f"\n✅ Metadata global créé : {OUTPUT_ALL}")
    print("Dimensions globales :", df_all.shape)

    summary = (
        df_all
        .groupby(["study", "dataset_role", "response", "response_binary"])
        .size()
        .reset_index(name="n_samples")
    )

    summary.to_csv(SUMMARY_FILE, index=False)

    print(f"\n✅ Résumé créé : {SUMMARY_FILE}")
    print("\nRésumé :")
    print(summary)

    print("\nRépartition par study :")
    print(df_all["study"].value_counts())

    print("\nRépartition response_binary :")
    print(df_all["response_binary"].value_counts(dropna=False))


if __name__ == "__main__":
    main()