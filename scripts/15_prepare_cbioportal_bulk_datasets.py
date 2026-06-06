from pathlib import Path
import tarfile
import pandas as pd
import numpy as np

print("--- 15 : PRÉPARATION DES DATASETS BULK CBIOPORTAL ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"

DATA_PROCESSED.mkdir(exist_ok=True)

ARCHIVES = {
    "CBIO_GIDE_2019": DATA_RAW / "7_CBIOPORTAL_GIDE_2019_raw.tar.gz",
    "CBIO_UCLA_2016": DATA_RAW / "7_CBIOPORTAL_UCLA_2016_raw.tar.gz",
    "CBIO_DFCI_2019": DATA_RAW / "7_CBIOPORTAL_DFCI_2019_raw.tar.gz",
    "CBIO_LIU_2019": DATA_RAW / "7_CBIOPORTAL_LIU_2019_raw.tar.gz",
}

OUTPUTS = {
    "CBIO_GIDE_2019": DATA_PROCESSED / "15_CBIO_GIDE_2019_clean.csv",
    "CBIO_UCLA_2016": DATA_PROCESSED / "15_CBIO_UCLA_2016_clean.csv",
    "CBIO_DFCI_2019": DATA_PROCESSED / "15_CBIO_DFCI_2019_clean.csv",
    "CBIO_LIU_2019": DATA_PROCESSED / "15_CBIO_LIU_2019_clean.csv",
}

SUMMARY_FILE = DATA_PROCESSED / "15_CBIOPORTAL_BULK_summary.csv"


# ============================================================
# CONFIGURATION DES FICHIERS À UTILISER
# ============================================================

# On privilégie les z-scores pour avoir une échelle plus comparable entre études cBioPortal.
# On garde aussi l'information du format d'origine.
DATASET_CONFIG = {
    "CBIO_GIDE_2019": {
        "expression_file": "mel_iatlas_gide_2019/data_mrna_seq_expression_zscores_ref_all_samples.txt",
        "expression_format": "z-score",
        "clinical_sample_file": "mel_iatlas_gide_2019/data_clinical_sample.txt",
        "clinical_patient_file": "mel_iatlas_gide_2019/data_clinical_patient.txt",
        "gene_col_candidates": ["Hugo_Symbol", "Gene Symbol", "GENE_SYMBOL"],
        "sample_id_col": "SAMPLE_ID",
        "patient_id_col": "PATIENT_ID",
        "response_cols": ["RESPONSE", "CLINICAL_BENEFIT"],
        "treatment_cols": ["SAMPLE_TREATMENT"],
    },
    "CBIO_UCLA_2016": {
        "expression_file": "mel_ucla_2016/data_mrna_seq_rpkm_zscores_ref_all_samples.txt",
        "expression_format": "RPKM z-score",
        "clinical_sample_file": "mel_ucla_2016/data_clinical_sample.txt",
        "clinical_patient_file": "mel_ucla_2016/data_clinical_patient.txt",
        "gene_col_candidates": ["Hugo_Symbol", "Gene Symbol", "GENE_SYMBOL"],
        "sample_id_col": "SAMPLE_ID",
        "patient_id_col": "PATIENT_ID",
        "response_cols": ["TREATMENT_RESPONSE", "DURABLE_CLINICAL_BENEFIT"],
        "treatment_cols": ["TREATMENT"],
    },
    "CBIO_DFCI_2019": {
        "expression_file": "mel_dfci_2019/data_mrna_seq_tpm_zscores_ref_all_samples.txt",
        "expression_format": "TPM z-score",
        "clinical_sample_file": "mel_dfci_2019/data_clinical_sample.txt",
        "clinical_patient_file": "mel_dfci_2019/data_clinical_patient.txt",
        "gene_col_candidates": ["Hugo_Symbol", "Gene Symbol", "GENE_SYMBOL"],
        "sample_id_col": "SAMPLE_ID",
        "patient_id_col": "PATIENT_ID",
        "response_cols": ["PFS_STATUS", "OS_STATUS"],
        "treatment_cols": ["IO_THERAPY"],
    },
    "CBIO_LIU_2019": {
        "expression_file": "mel_iatlas_liu_2019/data_mrna_seq_tpm_zscores_ref_all_samples.txt",
        "expression_format": "TPM z-score",
        "clinical_sample_file": "mel_iatlas_liu_2019/data_clinical_sample.txt",
        "clinical_patient_file": "mel_iatlas_liu_2019/data_clinical_patient.txt",
        "gene_col_candidates": ["Hugo_Symbol", "Gene Symbol", "GENE_SYMBOL"],
        "sample_id_col": "SAMPLE_ID",
        "patient_id_col": "PATIENT_ID",
        "response_cols": ["RESPONSE", "CLINICAL_BENEFIT"],
        "treatment_cols": ["SAMPLE_TREATMENT"],
    },
}


# ============================================================
# OUTILS
# ============================================================

def read_table_from_tar(tar, member_name):
    extracted = tar.extractfile(member_name)

    if extracted is None:
        raise ValueError(f"Impossible d'extraire {member_name}")

    df = pd.read_csv(
        extracted,
        sep="\t",
        comment="#",
        low_memory=False
    )

    return df


def clean_text(value):
    if pd.isna(value):
        return None
    value = str(value).replace('"', "").strip()
    return value if value != "" else None


def clean_gene_name(value):
    value = clean_text(value)

    if value is None:
        return None

    return value


def choose_gene_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col

    # fallback : première colonne
    return df.columns[0]


def extract_response(row, response_cols):
    """
    Récupère la première colonne de réponse disponible.
    """
    for col in response_cols:
        if col in row.index:
            value = clean_text(row[col])
            if value is not None and value.lower() not in {"nan", "na", "n/a", "unknown"}:
                return value

    return None


def extract_treatment(row, treatment_cols):
    """
    Récupère la première colonne de traitement disponible.
    """
    for col in treatment_cols:
        if col in row.index:
            value = clean_text(row[col])
            if value is not None and value.lower() not in {"nan", "na", "n/a", "unknown"}:
                return value

    return None


def standardize_response_binary(value):
    """
    Convertit les réponses cBioPortal en binaire quand c'est possible.

    1 = réponse favorable / bénéfice clinique
    0 = réponse défavorable / progression / non bénéfice

    Les règles exactes devront être validées biologiquement avec le prof.
    """
    if pd.isna(value):
        return None

    v = str(value).strip().upper()
    if "PROGRESSION" in v:
        return 0

    if "CENSORED" in v:
        return 1

    positive = {
        "R",
        "RESPONDER",
        "RESPONSE",
        "CR",
        "PR",
        "COMPLETE RESPONSE",
        "PARTIAL RESPONSE",
        "YES",
        "TRUE",
        "BENEFIT",
        "CLINICAL BENEFIT",
        "DURABLE CLINICAL BENEFIT",
        "DCB",
        "CB",
        "NO PROGRESSION",
        "NON-PROGRESSION",
    }

    negative = {
        "NR",
        "NON-RESPONDER",
        "NON_RESPONDER",
        "NONRESPONDER",
        "NO RESPONSE",
        "PROGRESSIVE DISEASE",
        "PD",
        "PROGRESSION",
        "PROGRESSED",
        "NO BENEFIT",
        "NO CLINICAL BENEFIT",
        "NON-DURABLE CLINICAL BENEFIT",
        "NDB",
        "NCB",
        "FALSE",
        "NO",
    }

    # Cas fréquents de cBioPortal : "0:DiseaseFree", "1:Recurred/Progressed", etc.
    if "DISEASEFREE" in v or "DISEASE FREE" in v:
        return 1

    if "RECURRED" in v or "PROGRESSED" in v or "DECEASED" in v or "DEAD" in v:
        return 0

    if v in positive:
        return 1

    if v in negative:
        return 0

    # SD est ambigu : stable disease peut être classé bénéfice clinique ou non selon la définition.
    # On laisse None pour ne pas imposer une règle sans validation.
    if v in {"SD", "STABLE DISEASE"}:
        return None

    return None


def prepare_expression_matrix(df_expr, gene_col):
    """
    cBioPortal expression est souvent :
    gènes en lignes, samples en colonnes.

    On convertit en :
    samples en lignes, gènes en colonnes.
    """

    df_expr[gene_col] = df_expr[gene_col].apply(clean_gene_name)
    df_expr = df_expr.dropna(subset=[gene_col]).copy()

    # Retirer les colonnes techniques autres que le gene symbol
    technical_cols = [
        col for col in df_expr.columns
        if col.lower() in {
            "entrez_gene_id",
            "cytoband",
            "description",
            "gene_id"
        }
    ]

    sample_cols = [
        col for col in df_expr.columns
        if col != gene_col and col not in technical_cols
    ]

    # Conversion numérique
    for col in sample_cols:
        df_expr[col] = pd.to_numeric(df_expr[col], errors="coerce")

    # Moyenne des doublons de gènes
    df_grouped = df_expr[[gene_col] + sample_cols].groupby(gene_col, as_index=True).mean(numeric_only=True)

    # Transpose
    df_t = df_grouped.T.reset_index()
    df_t = df_t.rename(columns={"index": "sample_id"})

    df_t["sample_id"] = df_t["sample_id"].astype(str).str.strip()

    return df_t


def merge_clinical_sample_patient(df_sample, df_patient, patient_id_col):
    """
    Fusionne clinique sample + patient quand possible.
    """

    if patient_id_col not in df_sample.columns:
        raise ValueError(f"Colonne patient absente du fichier sample : {patient_id_col}")

    if patient_id_col not in df_patient.columns:
        print(f"⚠️ Colonne patient absente du fichier patient : {patient_id_col}. On garde seulement clinical_sample.")
        return df_sample.copy()

    df_clin = df_sample.merge(
        df_patient,
        on=patient_id_col,
        how="left",
        suffixes=("", "_patient")
    )

    return df_clin


def prepare_cbio_dataset(dataset_name, archive_path, config):
    print("\n" + "=" * 100)
    print(f"PRÉPARATION CBIOPORTAL : {dataset_name}")
    print("=" * 100)

    if not archive_path.exists():
        raise FileNotFoundError(f"Archive introuvable : {archive_path}")

    with tarfile.open(archive_path, "r:gz") as tar:
        print("Lecture expression...")
        df_expr = read_table_from_tar(tar, config["expression_file"])

        print("Lecture clinical_sample...")
        df_sample = read_table_from_tar(tar, config["clinical_sample_file"])

        print("Lecture clinical_patient...")
        df_patient = read_table_from_tar(tar, config["clinical_patient_file"])

    print("Dimensions expression brute :", df_expr.shape)
    print("Dimensions clinical_sample :", df_sample.shape)
    print("Dimensions clinical_patient :", df_patient.shape)

    gene_col = choose_gene_column(df_expr, config["gene_col_candidates"])
    print("Colonne gène utilisée :", gene_col)
    print("Format expression choisi :", config["expression_format"])

    # Expression : samples x genes
    df_expr_t = prepare_expression_matrix(df_expr, gene_col)

    print("Dimensions expression transposée :", df_expr_t.shape)

    # Clinique sample + patient
    df_clin = merge_clinical_sample_patient(
        df_sample=df_sample,
        df_patient=df_patient,
        patient_id_col=config["patient_id_col"]
    )

    print("Dimensions clinique fusionnée :", df_clin.shape)

    # Standardisation clinique
    sample_id_col = config["sample_id_col"]
    patient_id_col = config["patient_id_col"]

    if sample_id_col not in df_clin.columns:
        raise ValueError(f"Colonne sample absente : {sample_id_col}")

    df_meta = pd.DataFrame(index=df_clin.index)
    df_meta["study"] = dataset_name
    df_meta["sample_id"] = df_clin[sample_id_col].astype(str).str.strip()
    df_meta["patient_id"] = df_clin[patient_id_col].astype(str).str.strip()
    df_meta["response"] = df_clin.apply(lambda row: extract_response(row, config["response_cols"]), axis=1)
    df_meta["response_binary"] = df_meta["response"].apply(standardize_response_binary)
    df_meta["treatment"] = df_clin.apply(lambda row: extract_treatment(row, config["treatment_cols"]), axis=1)
    df_meta["expression_format"] = config["expression_format"]
    df_meta["dataset_role"] = "cbioportal_bulk_immunotherapy"

    # Ajouter quelques colonnes cliniques utiles si elles existent
    optional_cols = [
        "BIOPSY_SITE",
        "BIOPSY_TIME",
        "SAMPLE_TREATMENT",
        "TREATMENT",
        "IO_THERAPY",
        "CLINICAL_BENEFIT",
        "DURABLE_CLINICAL_BENEFIT",
        "RESPONSE",
        "TREATMENT_RESPONSE",
        "PROGRESSION",
        "PFS_STATUS",
        "PFS_MONTHS",
        "OS_STATUS",
        "OS_MONTHS",
        "SEX",
        "AGE_AT_DIAGNOSIS",
        "M_STAGE",
        "PRIOR_ICI_RX",
        "PRIOR_RX",
        "PRIOR_CTLA4",
        "PRIOR_MAPK_TX",
        "TMB_NONSYNONYMOUS",
        "MUTATION_LOAD",
    ]

    for col in optional_cols:
        if col in df_clin.columns and col not in df_meta.columns:
            df_meta[col.lower()] = df_clin[col]

    # Fusion metadata + expression
    df_clean = df_meta.merge(df_expr_t, on="sample_id", how="inner")

    print("Dimensions après fusion metadata + expression :", df_clean.shape)

    missing_from_expr = sorted(set(df_meta["sample_id"]) - set(df_expr_t["sample_id"]))
    extra_expr_samples = sorted(set(df_expr_t["sample_id"]) - set(df_meta["sample_id"]))

    print("Samples metadata absents expression :", len(missing_from_expr))
    if len(missing_from_expr) > 0:
        print(missing_from_expr[:20])

    print("Samples expression absents metadata :", len(extra_expr_samples))
    if len(extra_expr_samples) > 0:
        print(extra_expr_samples[:20])

    # Sauvegarde
    output_path = OUTPUTS[dataset_name]
    df_clean.to_csv(output_path, index=False)

    print("✅ Fichier clean créé :", output_path)

    print("\nRépartition response :")
    print(df_clean["response"].value_counts(dropna=False))

    print("\nRépartition response_binary :")
    print(df_clean["response_binary"].value_counts(dropna=False))

    print("\nRépartition treatment :")
    print(df_clean["treatment"].value_counts(dropna=False))

    print("\nAperçu metadata :")
    preview_cols = [
        "study",
        "sample_id",
        "patient_id",
        "response",
        "response_binary",
        "treatment",
        "expression_format",
    ]
    print(df_clean[preview_cols].head())

    n_metadata_cols = len([col for col in df_clean.columns if col not in df_expr_t.columns or col == "sample_id"])
    n_gene_cols = df_clean.shape[1] - n_metadata_cols

    return {
        "dataset": dataset_name,
        "output_file": str(output_path),
        "expression_format": config["expression_format"],
        "n_expression_genes": df_expr_t.shape[1] - 1,
        "n_expression_samples": df_expr_t.shape[0],
        "n_clinical_samples": df_meta.shape[0],
        "n_clean_samples": df_clean.shape[0],
        "n_clean_columns": df_clean.shape[1],
        "n_gene_cols_estimated": n_gene_cols,
        "n_response_binary_0": int((df_clean["response_binary"] == 0).sum()),
        "n_response_binary_1": int((df_clean["response_binary"] == 1).sum()),
        "n_response_binary_missing": int(df_clean["response_binary"].isna().sum()),
        "n_missing_samples_from_expression": len(missing_from_expr),
        "n_expression_samples_without_metadata": len(extra_expr_samples),
    }


# ============================================================
# MAIN
# ============================================================

def main():
    summary_rows = []

    for dataset_name, archive_path in ARCHIVES.items():
        config = DATASET_CONFIG[dataset_name]

        try:
            summary = prepare_cbio_dataset(
                dataset_name=dataset_name,
                archive_path=archive_path,
                config=config
            )
            summary_rows.append(summary)

        except Exception as e:
            print(f"\n❌ Erreur pour {dataset_name} : {e}")
            summary_rows.append({
                "dataset": dataset_name,
                "error": str(e),
            })

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(SUMMARY_FILE, index=False)

    print("\n" + "=" * 100)
    print("RÉSUMÉ GLOBAL SCRIPT 15")
    print("=" * 100)
    print(df_summary)

    print("\n✅ Résumé sauvegardé :", SUMMARY_FILE)


if __name__ == "__main__":
    main()