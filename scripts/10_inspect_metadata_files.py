from pathlib import Path
import gzip
import pandas as pd

print("--- 10 : INSPECTION DES FICHIERS METADATA / SERIES MATRIX ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"

DATA_PROCESSED.mkdir(exist_ok=True)

SERIES_FILES = {
    "GSE145996": DATA_RAW / "6_IMMUNOTHERAPY_GSE145996_series_matrix_raw.txt.gz",
    "GSE168204_GPL11154": DATA_RAW / "6_IMMUNOTHERAPY_GSE168204_GPL11154_series_matrix_raw.txt.gz",
    "GSE168204_GPL18573": DATA_RAW / "6_IMMUNOTHERAPY_GSE168204_GPL18573_series_matrix_raw.txt.gz",
    "GSE244982": DATA_RAW / "6_IMMUNOTHERAPY_GSE244982_series_matrix_raw.txt.gz",
    "GSE196434": DATA_RAW / "5_TARGETED_THERAPY_GSE196434_series_matrix_raw.txt.gz",
    "GSE120575": DATA_RAW / "8_SINGLE_CELL_GSE120575_series_matrix_raw.txt.gz",
}

REPORT_FILE = DATA_PROCESSED / "10_METADATA_inspection_report.txt"
SUMMARY_FILE = DATA_PROCESSED / "10_METADATA_inspection_summary.csv"


# ============================================================
# OUTILS
# ============================================================

def write_report(text):
    with open(REPORT_FILE, "a", encoding="utf-8") as f:
        f.write(str(text) + "\n")


def print_and_write(text):
    print(text)
    write_report(text)


def read_gz_lines(path):
    """Lit un fichier .gz en texte."""
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        return f.readlines()


def clean_value(value):
    if pd.isna(value):
        return None
    value = str(value).replace('"', "").strip()
    return value if value != "" else None


def parse_series_matrix(path):
    """
    Parse les lignes utiles d'un fichier GEO series_matrix.
    On extrait :
    - !Sample_geo_accession
    - !Sample_title
    - !Sample_source_name_ch1
    - toutes les lignes !Sample_characteristics_ch1
    - !Sample_description
    - !Sample_platform_id
    """

    lines = read_gz_lines(path)

    metadata_rows = {}

    useful_prefixes = [
        "!Sample_geo_accession",
        "!Sample_title",
        "!Sample_source_name_ch1",
        "!Sample_description",
        "!Sample_platform_id",
        "!Sample_characteristics_ch1",
    ]

    for line in lines:
        line = line.strip()

        if not line.startswith("!Sample_"):
            continue

        parts = line.split("\t")
        key = parts[0]
        values = [clean_value(v) for v in parts[1:]]

        if any(key.startswith(prefix) for prefix in useful_prefixes):
            # Plusieurs lignes !Sample_characteristics_ch1 existent.
            # On les garde séparées selon leur contenu.
            if key == "!Sample_characteristics_ch1":
                if len(values) > 0 and values[0] is not None and ":" in values[0]:
                    label = values[0].split(":")[0].strip().lower()
                    label = (
                        label.replace(" ", "_")
                        .replace("/", "_")
                        .replace("-", "_")
                        .replace("(", "")
                        .replace(")", "")
                    )
                    clean_key = f"characteristics_{label}"
                else:
                    clean_key = "characteristics_unknown"

                # Éviter écrasement si deux labels identiques
                original_key = clean_key
                i = 2
                while clean_key in metadata_rows:
                    clean_key = f"{original_key}_{i}"
                    i += 1

                cleaned_values = []
                for v in values:
                    if v is None:
                        cleaned_values.append(None)
                    elif ":" in v:
                        cleaned_values.append(v.split(":", 1)[1].strip())
                    else:
                        cleaned_values.append(v)

                metadata_rows[clean_key] = cleaned_values

            else:
                clean_key = key.replace("!Sample_", "").lower()
                metadata_rows[clean_key] = values

    if not metadata_rows:
        return pd.DataFrame()

    # Longueur maximale
    max_len = max(len(v) for v in metadata_rows.values())

    # Normaliser longueurs
    for k, values in metadata_rows.items():
        if len(values) < max_len:
            metadata_rows[k] = values + [None] * (max_len - len(values))

    df = pd.DataFrame(metadata_rows)

    return df


def detect_relevant_columns(df):
    """Repère les colonnes potentiellement utiles."""
    keywords = [
        "response",
        "responder",
        "nonresponder",
        "non_responder",
        "clinical",
        "benefit",
        "recist",
        "progression",
        "progressed",
        "survival",
        "vital",
        "death",
        "alive",
        "treatment",
        "therapy",
        "anti",
        "pd1",
        "pd_1",
        "pdl1",
        "ctla",
        "ipilimumab",
        "nivolumab",
        "pembrolizumab",
        "time",
        "timepoint",
        "baseline",
        "pre",
        "on",
        "post",
        "patient",
        "subject",
        "sample",
        "stage",
        "sex",
        "gender",
        "age",
        "braf",
        "nras",
        "nf1",
        "mapk",
        "mek",
        "raf",
    ]

    matched = []

    for col in df.columns:
        col_lower = str(col).lower()
        if any(k in col_lower for k in keywords):
            matched.append(col)

    return matched


def summarize_column_values(df, cols, max_cols=20, max_values=12):
    """Affiche les valeurs uniques principales de certaines colonnes."""
    for col in cols[:max_cols]:
        print_and_write(f"\nColonne : {col}")
        try:
            values = df[col].dropna().astype(str).unique().tolist()
            print_and_write(f"Nombre valeurs uniques : {len(values)}")
            print_and_write(values[:max_values])
        except Exception as e:
            print_and_write(f"Impossible de résumer la colonne {col} : {e}")


# ============================================================
# MAIN
# ============================================================

def main():
    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write("RAPPORT D'INSPECTION DES METADATA GEO / SERIES MATRIX\n")
        f.write("=" * 80 + "\n")

    summary_rows = []

    for dataset_name, path in SERIES_FILES.items():
        print_and_write("\n\n" + "#" * 80)
        print_and_write(f"DATASET : {dataset_name}")
        print_and_write("#" * 80)

        if not path.exists():
            print_and_write(f"❌ Fichier introuvable : {path}")
            summary_rows.append({
                "dataset": dataset_name,
                "file_exists": False,
                "n_samples": None,
                "n_columns": None,
                "relevant_columns": None,
            })
            continue

        print_and_write(f"Fichier trouvé : {path}")

        try:
            df_meta = parse_series_matrix(path)

            if df_meta.empty:
                print_and_write("❌ Aucun metadata extrait.")
                summary_rows.append({
                    "dataset": dataset_name,
                    "file_exists": True,
                    "n_samples": 0,
                    "n_columns": 0,
                    "relevant_columns": None,
                })
                continue

            print_and_write(f"Dimensions metadata : {df_meta.shape[0]} échantillons x {df_meta.shape[1]} colonnes")

            print_and_write("\nColonnes extraites :")
            print_and_write(df_meta.columns.tolist())

            relevant_cols = detect_relevant_columns(df_meta)

            print_and_write("\nColonnes potentiellement utiles :")
            print_and_write(relevant_cols)

            print_and_write("\nAperçu metadata :")
            print_and_write(df_meta.head().to_string())

            print_and_write("\nRésumé des colonnes utiles :")
            summarize_column_values(df_meta, relevant_cols)

            # Sauvegarde CSV par dataset
            output_csv = DATA_PROCESSED / f"10_METADATA_{dataset_name}.csv"
            df_meta.to_csv(output_csv, index=False)

            print_and_write(f"\n✅ Metadata extrait sauvegardé : {output_csv}")

            summary_rows.append({
                "dataset": dataset_name,
                "file_exists": True,
                "n_samples": df_meta.shape[0],
                "n_columns": df_meta.shape[1],
                "relevant_columns": "; ".join(relevant_cols),
                "output_csv": str(output_csv),
            })

        except Exception as e:
            print_and_write(f"❌ Erreur pendant inspection de {dataset_name} : {e}")
            summary_rows.append({
                "dataset": dataset_name,
                "file_exists": True,
                "error": str(e),
            })

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(SUMMARY_FILE, index=False)

    print_and_write("\n\n" + "=" * 80)
    print_and_write("RÉSUMÉ GLOBAL")
    print_and_write("=" * 80)
    print_and_write(df_summary.to_string(index=False))

    print("\n✅ Inspection metadata terminée.")
    print(f"Rapport texte : {REPORT_FILE}")
    print(f"Résumé CSV    : {SUMMARY_FILE}")
    print("Metadata CSV  : data_processed/10_METADATA_*.csv")


if __name__ == "__main__":
    main()