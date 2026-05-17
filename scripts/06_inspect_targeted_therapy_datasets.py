from pathlib import Path
import pandas as pd

print("--- 06 : INSPECTION DES DATASETS THÉRAPIES CIBLÉES BRAF/MEK ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"

GSE196434_FILE = DATA_RAW / "5_TARGETED_THERAPY_GSE196434_raw.txt.gz"
GSE77940_FILE = DATA_RAW / "5_TARGETED_THERAPY_GSE77940_pre_post_raw.xlsx"

REPORT_FILE = DATA_PROCESSED / "6_TARGETED_THERAPY_inspection_report.txt"


# ============================================================
# OUTILS
# ============================================================

def write_report(text):
    """Ajoute du texte au rapport."""
    with open(REPORT_FILE, "a", encoding="utf-8") as f:
        f.write(text + "\n")


def print_and_write(text):
    """Affiche dans le terminal et écrit dans le rapport."""
    print(text)
    write_report(text)


def inspect_dataframe(df, name):
    """Résumé simple d'un DataFrame."""
    print_and_write("\n" + "=" * 80)
    print_and_write(f"INSPECTION : {name}")
    print_and_write("=" * 80)

    print_and_write(f"Dimensions : {df.shape[0]} lignes x {df.shape[1]} colonnes")

    print_and_write("\nColonnes principales :")
    print_and_write(str(df.columns[:20].tolist()))

    print_and_write("\nTypes des colonnes principales :")
    print_and_write(str(df.dtypes.head(20)))

    print_and_write("\nNombre total de valeurs manquantes :")
    print_and_write(str(df.isna().sum().sum()))

    print_and_write("\nAperçu des 5 premières lignes :")
    print_and_write(str(df.head()))


# ============================================================
# INSPECTION GSE196434
# ============================================================

def inspect_gse196434():
    print_and_write("\n\n" + "#" * 80)
    print_and_write("DATASET GSE196434 — THÉRAPIE CIBLÉE BRAF/MEK")
    print_and_write("#" * 80)

    if not GSE196434_FILE.exists():
        print_and_write(f"❌ Fichier introuvable : {GSE196434_FILE}")
        return

    print_and_write(f"Fichier trouvé : {GSE196434_FILE}")

    try:
        # Le fichier est compressé .txt.gz
        # On laisse pandas détecter le séparateur avec sep=None.
        df = pd.read_csv(
            GSE196434_FILE,
            sep=None,
            engine="python",
            compression="gzip"
        )

        inspect_dataframe(df, "GSE196434 raw")

        # Sauvegarde d'un petit aperçu pour inspection manuelle
        preview_path = DATA_PROCESSED / "6_GSE196434_preview.csv"
        df.head(20).to_csv(preview_path, index=False)

        print_and_write(f"\n✅ Aperçu sauvegardé : {preview_path}")

        # Recherche de colonnes potentiellement utiles
        keywords = [
            "patient", "sample", "response", "resistance", "responder",
            "treatment", "therapy", "braf", "mek", "group", "time"
        ]

        matched_cols = [
            col for col in df.columns
            if any(k in str(col).lower() for k in keywords)
        ]

        print_and_write("\nColonnes potentiellement cliniques / labels :")
        print_and_write(str(matched_cols[:50]))

    except Exception as e:
        print_and_write(f"❌ Erreur pendant la lecture de GSE196434 : {e}")


# ============================================================
# INSPECTION GSE77940
# ============================================================

def inspect_gse77940():
    print_and_write("\n\n" + "#" * 80)
    print_and_write("DATASET GSE77940 — PRE / POST RAF-MEK")
    print_and_write("#" * 80)

    if not GSE77940_FILE.exists():
        print_and_write(f"❌ Fichier introuvable : {GSE77940_FILE}")
        return

    print_and_write(f"Fichier trouvé : {GSE77940_FILE}")

    try:
        # Voir les feuilles du fichier Excel
        excel = pd.ExcelFile(GSE77940_FILE)

        print_and_write("\nFeuilles Excel trouvées :")
        print_and_write(str(excel.sheet_names))

        for sheet in excel.sheet_names:
            print_and_write("\n" + "-" * 80)
            print_and_write(f"Lecture feuille : {sheet}")
            print_and_write("-" * 80)

            df = pd.read_excel(GSE77940_FILE, sheet_name=sheet)

            inspect_dataframe(df, f"GSE77940 - feuille {sheet}")

            # Sauvegarde aperçu de chaque feuille
            safe_sheet_name = sheet.replace(" ", "_").replace("/", "_")
            preview_path = DATA_PROCESSED / f"6_GSE77940_preview_{safe_sheet_name}.csv"
            df.head(20).to_csv(preview_path, index=False)

            print_and_write(f"\n✅ Aperçu sauvegardé : {preview_path}")

            keywords = [
                "patient", "sample", "response", "resistance", "responder",
                "treatment", "therapy", "braf", "mek", "raf", "pre", "post",
                "relapse", "progression", "time", "group"
            ]

            matched_cols = [
                col for col in df.columns
                if any(k in str(col).lower() for k in keywords)
            ]

            print_and_write("\nColonnes potentiellement cliniques / labels :")
            print_and_write(str(matched_cols[:50]))

    except Exception as e:
        print_and_write(f"❌ Erreur pendant la lecture de GSE77940 : {e}")


# ============================================================
# MAIN
# ============================================================

def main():
    # Réinitialiser le rapport à chaque exécution
    DATA_PROCESSED.mkdir(exist_ok=True)

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write("RAPPORT D'INSPECTION — DATASETS THÉRAPIES CIBLÉES BRAF/MEK\n")
        f.write("=" * 80 + "\n")

    inspect_gse196434()
    inspect_gse77940()

    print("\n✅ Inspection terminée.")
    print(f"Rapport créé : {REPORT_FILE}")
    print("\nFichiers preview créés dans data_processed/ si la lecture a fonctionné.")


if __name__ == "__main__":
    main()