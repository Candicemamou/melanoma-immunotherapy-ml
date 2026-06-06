from pathlib import Path
import tarfile
import zipfile
import pandas as pd

print("--- 14 : INSPECTION CBIOPORTAL + MENDELEY SINGLE-CELL ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"

DATA_PROCESSED.mkdir(exist_ok=True)

CBIO_ARCHIVES = {
    "CBIO_GIDE_2019": DATA_RAW / "7_CBIOPORTAL_GIDE_2019_raw.tar.gz",
    "CBIO_UCLA_2016": DATA_RAW / "7_CBIOPORTAL_UCLA_2016_raw.tar.gz",
    "CBIO_DFCI_2019": DATA_RAW / "7_CBIOPORTAL_DFCI_2019_raw.tar.gz",
    "CBIO_LIU_2019": DATA_RAW / "7_CBIOPORTAL_LIU_2019_raw.tar.gz",
}

MENDELEY_ZIP = DATA_RAW / "8_SINGLE_CELL_MENDELEY_MELANOMA_CITESEQ_raw.zip"

REPORT_FILE = DATA_PROCESSED / "14_CBIOPORTAL_MENDELEY_inspection_report.txt"
SUMMARY_FILE = DATA_PROCESSED / "14_CBIOPORTAL_MENDELEY_summary.csv"


# ============================================================
# OUTILS GÉNÉRAUX
# ============================================================

def write_report(text):
    with open(REPORT_FILE, "a", encoding="utf-8") as f:
        f.write(str(text) + "\n")


def print_and_write(text):
    print(text)
    write_report(text)


def safe_read_table_from_tar(tar, member_name, nrows=None):
    """
    Lit un fichier texte/tabulaire à l'intérieur d'une archive tar.gz.
    """
    extracted = tar.extractfile(member_name)

    if extracted is None:
        raise ValueError(f"Impossible d'extraire {member_name}")

    # cBioPortal utilise généralement le tab comme séparateur
    df = pd.read_csv(
        extracted,
        sep="\t",
        comment="#",
        low_memory=False,
        nrows=nrows
    )

    return df


def safe_read_table_from_zip(zip_obj, member_name, nrows=None):
    """
    Lit un fichier texte/csv/tsv à l'intérieur d'un zip.
    """
    with zip_obj.open(member_name) as f:
        if member_name.lower().endswith(".tsv"):
            sep = "\t"
        else:
            sep = ","

        df = pd.read_csv(
            f,
            sep=sep,
            low_memory=False,
            nrows=nrows
        )

    return df


def find_keywords_in_columns(columns):
    """
    Repère les colonnes cliniques intéressantes.
    """
    keywords = [
        "response",
        "recist",
        "radiographic",
        "benefit",
        "progression",
        "pfs",
        "os",
        "survival",
        "vital",
        "death",
        "status",
        "treatment",
        "therapy",
        "immunotherapy",
        "anti",
        "pd1",
        "pd-1",
        "pdl1",
        "pd-l1",
        "ctla",
        "ipilimumab",
        "nivolumab",
        "pembrolizumab",
        "melanoma",
        "stage",
        "sex",
        "gender",
        "age",
        "biopsy",
        "sample",
        "patient",
        "tmb",
        "mutation",
        "braf",
        "mapk",
        "prior",
        "resistance",
        "resistant",
    ]

    matched = []

    for col in columns:
        col_lower = str(col).lower()
        if any(k in col_lower for k in keywords):
            matched.append(col)

    return matched


def summarize_relevant_columns(df, max_cols=25, max_values=12):
    """
    Affiche les valeurs principales des colonnes cliniques pertinentes.
    """
    relevant_cols = find_keywords_in_columns(df.columns)

    print_and_write("\nColonnes potentiellement utiles :")
    print_and_write(relevant_cols)

    for col in relevant_cols[:max_cols]:
        print_and_write(f"\nColonne : {col}")
        try:
            values = df[col].dropna().astype(str).unique().tolist()
            print_and_write(f"Nombre valeurs uniques : {len(values)}")
            print_and_write(values[:max_values])
        except Exception as e:
            print_and_write(f"Impossible de résumer {col} : {e}")

    return relevant_cols


def guess_expression_format(filename):
    """
    Devine le type d'expression à partir du nom du fichier.
    """
    name = filename.lower()

    if "tpm" in name:
        return "TPM"
    if "rpkm" in name:
        return "RPKM"
    if "fpkm" in name:
        return "FPKM"
    if "zscore" in name or "zscores" in name:
        return "z-score"
    if "count" in name or "counts" in name:
        return "counts"
    if "mrna" in name or "expression" in name:
        return "mRNA expression, format exact à confirmer"

    return "unknown"


def is_expression_file(filename):
    """
    Repère les fichiers d'expression moléculaire.
    """
    name = filename.lower()

    if not name.endswith(".txt"):
        return False

    expression_patterns = [
        "data_mrna",
        "expression",
        "rpkm",
        "tpm",
        "fpkm",
        "zscore",
        "zscores",
    ]

    excluded_patterns = [
        "genetic_alteration",
        "case_lists",
        "meta_",
        "clinical",
        "mutation",
        "cna",
        "timeline",
        "signature",
    ]

    return (
        any(p in name for p in expression_patterns)
        and not any(p in name for p in excluded_patterns)
    )


def is_clinical_file(filename):
    name = filename.lower()
    return "clinical" in name and name.endswith(".txt")


# ============================================================
# INSPECTION CBIOPORTAL
# ============================================================

def inspect_cbio_archive(dataset_name, archive_path):
    print_and_write("\n\n" + "#" * 100)
    print_and_write(f"CBIOPORTAL DATASET : {dataset_name}")
    print_and_write("#" * 100)

    if not archive_path.exists():
        print_and_write(f"❌ Archive introuvable : {archive_path}")
        return []

    print_and_write(f"Archive trouvée : {archive_path}")

    rows = []

    with tarfile.open(archive_path, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.isfile()]
        names = [m.name for m in members]

        print_and_write("\nFichiers dans l'archive :")
        for name in names:
            print_and_write(f"- {name}")

        expression_files = [name for name in names if is_expression_file(name)]
        clinical_files = [name for name in names if is_clinical_file(name)]

        print_and_write("\nFichiers expression détectés :")
        print_and_write(expression_files)

        print_and_write("\nFichiers cliniques détectés :")
        print_and_write(clinical_files)

        # ----------------------------
        # Expression files
        # ----------------------------
        for file_name in expression_files:
            print_and_write("\n" + "=" * 80)
            print_and_write(f"EXPRESSION : {file_name}")
            print_and_write("=" * 80)

            expr_format = guess_expression_format(file_name)

            try:
                df_preview = safe_read_table_from_tar(tar, file_name, nrows=5)

                print_and_write(f"Format supposé : {expr_format}")
                print_and_write(f"Dimensions preview : {df_preview.shape[0]} lignes x {df_preview.shape[1]} colonnes")
                print_and_write("Colonnes principales :")
                print_and_write(df_preview.columns[:20].tolist())
                print_and_write("Aperçu :")
                print_and_write(df_preview.head().to_string())

                # Estimation samples/gènes
                first_cols = df_preview.columns[:5].tolist()

                rows.append({
                    "dataset": dataset_name,
                    "source": "cBioPortal",
                    "file": file_name,
                    "file_type": "expression",
                    "expression_format_guess": expr_format,
                    "n_preview_rows": df_preview.shape[0],
                    "n_preview_cols": df_preview.shape[1],
                    "relevant_columns": "; ".join(map(str, first_cols)),
                    "notes": "Inspecter fichier complet si retenu"
                })

            except Exception as e:
                print_and_write(f"❌ Erreur lecture expression {file_name} : {e}")

                rows.append({
                    "dataset": dataset_name,
                    "source": "cBioPortal",
                    "file": file_name,
                    "file_type": "expression",
                    "expression_format_guess": expr_format,
                    "error": str(e)
                })

        # ----------------------------
        # Clinical files
        # ----------------------------
        for file_name in clinical_files:
            print_and_write("\n" + "=" * 80)
            print_and_write(f"CLINIQUE : {file_name}")
            print_and_write("=" * 80)

            try:
                df_clin = safe_read_table_from_tar(tar, file_name)

                print_and_write(f"Dimensions clinique : {df_clin.shape[0]} lignes x {df_clin.shape[1]} colonnes")
                print_and_write("Colonnes :")
                print_and_write(df_clin.columns.tolist())
                print_and_write("Aperçu :")
                print_and_write(df_clin.head().to_string())

                relevant_cols = summarize_relevant_columns(df_clin)

                output_name = f"14_{dataset_name}_{Path(file_name).stem}.csv"
                output_path = DATA_PROCESSED / output_name
                df_clin.to_csv(output_path, index=False)

                print_and_write(f"\n✅ Clinique sauvegardée : {output_path}")

                rows.append({
                    "dataset": dataset_name,
                    "source": "cBioPortal",
                    "file": file_name,
                    "file_type": "clinical",
                    "expression_format_guess": None,
                    "n_rows": df_clin.shape[0],
                    "n_cols": df_clin.shape[1],
                    "relevant_columns": "; ".join(relevant_cols),
                    "output_csv": str(output_path),
                    "notes": "Colonnes RECIST/response à vérifier"
                })

            except Exception as e:
                print_and_write(f"❌ Erreur lecture clinique {file_name} : {e}")

                rows.append({
                    "dataset": dataset_name,
                    "source": "cBioPortal",
                    "file": file_name,
                    "file_type": "clinical",
                    "error": str(e)
                })

    return rows


# ============================================================
# INSPECTION MENDELEY SINGLE-CELL
# ============================================================

def inspect_mendeley_zip(zip_path):
    print_and_write("\n\n" + "#" * 100)
    print_and_write("MENDELEY SINGLE-CELL / CITE-SEQ")
    print_and_write("#" * 100)

    rows = []

    if not zip_path.exists():
        print_and_write(f"❌ Zip introuvable : {zip_path}")
        return rows

    print_and_write(f"Zip trouvé : {zip_path}")

    with zipfile.ZipFile(zip_path, "r") as z:
        names = z.namelist()

        print_and_write("\nFichiers dans le zip :")
        for name in names:
            print_and_write(f"- {name}")

        interesting_files = [
            name for name in names
            if name.lower().endswith((".csv", ".tsv"))
        ]

        print_and_write("\nFichiers tabulaires détectés :")
        print_and_write(interesting_files)

        for file_name in interesting_files:
            print_and_write("\n" + "=" * 80)
            print_and_write(f"MENDELEY FILE : {file_name}")
            print_and_write("=" * 80)

            file_lower = file_name.lower()

            if "rna_counts" in file_lower:
                file_type = "single_cell_rna_counts"
                expr_format = "counts"
                nrows = 5
            elif "adt_counts" in file_lower:
                file_type = "single_cell_adt_counts"
                expr_format = "ADT counts"
                nrows = 5
            elif "metadata" in file_lower:
                file_type = "metadata"
                expr_format = None
                nrows = None
            elif "barcode" in file_lower:
                file_type = "sample_barcode_assignment"
                expr_format = None
                nrows = None
            else:
                file_type = "unknown_tabular"
                expr_format = None
                nrows = 20

            try:
                df = safe_read_table_from_zip(z, file_name, nrows=nrows)

                print_and_write(f"Type détecté : {file_type}")
                print_and_write(f"Format expression supposé : {expr_format}")
                print_and_write(f"Dimensions lues : {df.shape[0]} lignes x {df.shape[1]} colonnes")
                print_and_write("Colonnes principales :")
                print_and_write(df.columns[:30].tolist())
                print_and_write("Aperçu :")
                print_and_write(df.head().to_string())

                relevant_cols = []

                if file_type in ["metadata", "sample_barcode_assignment"]:
                    relevant_cols = summarize_relevant_columns(df)

                    output_name = f"14_MENDELEY_{Path(file_name).stem}.csv"
                    output_path = DATA_PROCESSED / output_name
                    df.to_csv(output_path, index=False)
                    print_and_write(f"\n✅ Metadata/barcode sauvegardé : {output_path}")
                else:
                    # Pour les gros fichiers RNA/ADT, on ne sauvegarde qu'un aperçu
                    output_name = f"14_MENDELEY_{Path(file_name).stem}_preview.csv"
                    output_path = DATA_PROCESSED / output_name
                    df.to_csv(output_path, index=False)
                    print_and_write(f"\n✅ Preview sauvegardé : {output_path}")

                rows.append({
                    "dataset": "MENDELEY_MELANOMA_CITESEQ",
                    "source": "Mendeley",
                    "file": file_name,
                    "file_type": file_type,
                    "expression_format_guess": expr_format,
                    "n_rows_read": df.shape[0],
                    "n_cols_read": df.shape[1],
                    "relevant_columns": "; ".join(relevant_cols),
                    "output_csv": str(output_path),
                    "notes": "RNA/ADT counts à pseudo-bulker par sample/patient"
                })

            except Exception as e:
                print_and_write(f"❌ Erreur lecture {file_name} : {e}")

                rows.append({
                    "dataset": "MENDELEY_MELANOMA_CITESEQ",
                    "source": "Mendeley",
                    "file": file_name,
                    "file_type": file_type,
                    "expression_format_guess": expr_format,
                    "error": str(e)
                })

    return rows


# ============================================================
# MAIN
# ============================================================

def main():
    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write("RAPPORT 14 — INSPECTION CBIOPORTAL + MENDELEY\n")
        f.write("=" * 100 + "\n")

    all_rows = []

    print_and_write("Inspection des archives cBioPortal...")

    for dataset_name, archive_path in CBIO_ARCHIVES.items():
        rows = inspect_cbio_archive(dataset_name, archive_path)
        all_rows.extend(rows)

    print_and_write("\nInspection du dataset Mendeley single-cell...")
    mendeley_rows = inspect_mendeley_zip(MENDELEY_ZIP)
    all_rows.extend(mendeley_rows)

    df_summary = pd.DataFrame(all_rows)
    df_summary.to_csv(SUMMARY_FILE, index=False)

    print_and_write("\n\n" + "=" * 100)
    print_and_write("RÉSUMÉ GLOBAL SCRIPT 14")
    print_and_write("=" * 100)
    print_and_write(df_summary.to_string(index=False))

    print("\n✅ Inspection terminée.")
    print(f"Rapport texte : {REPORT_FILE}")
    print(f"Résumé CSV    : {SUMMARY_FILE}")
    print("Fichiers extraits / previews : data_processed/14_*")


if __name__ == "__main__":
    main()