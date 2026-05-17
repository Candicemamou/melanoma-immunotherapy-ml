from pathlib import Path
import pandas as pd

print("--- 09 : INSPECTION DES DATASETS IMMUNOTHÉRAPIE SUPPLÉMENTAIRES ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"

MAIN_ML_FILE = DATA_PROCESSED / "3_DATASET_RNA_MULTICOHORTE_ML_READY.csv"

GSE168204_FILE = DATA_RAW / "6_IMMUNOTHERAPY_GSE168204_counts_raw.csv.gz"
GSE244982_FILE = DATA_RAW / "6_IMMUNOTHERAPY_GSE244982_bulkRNAseq_raw.txt.gz"
GSE145996_FILE = DATA_RAW / "6_IMMUNOTHERAPY_GSE145996_FPKM_raw.txt.gz"

REPORT_FILE = DATA_PROCESSED / "9_EXTRA_IMMUNOTHERAPY_inspection_report.txt"


# ============================================================
# OUTILS
# ============================================================

def write_report(text):
    with open(REPORT_FILE, "a", encoding="utf-8") as f:
        f.write(str(text) + "\n")


def print_and_write(text):
    print(text)
    write_report(text)


def detect_separator(file_path):
    """
    Essaie de deviner le séparateur.
    """
    with pd.read_csv(file_path, compression="infer", nrows=1, sep=None, engine="python") as preview:
        return None


def read_table_auto(file_path):
    """
    Lecture robuste d'un fichier texte/csv compressé ou non.
    Essaie plusieurs séparateurs courants.
    """
    separators = [",", "\t", ";"]

    last_error = None

    for sep in separators:
        try:
            df = pd.read_csv(
                file_path,
                sep=sep,
                compression="infer",
                low_memory=False
            )

            # Si tout est lu comme une seule colonne, mauvais séparateur
            if df.shape[1] > 1:
                return df, sep

        except Exception as e:
            last_error = e

    # Dernier essai : détection automatique
    try:
        df = pd.read_csv(
            file_path,
            sep=None,
            engine="python",
            compression="infer"
        )
        return df, "auto"
    except Exception as e:
        raise RuntimeError(f"Impossible de lire {file_path}. Dernière erreur : {last_error}. Auto erreur : {e}")


def clean_gene_name(gene):
    """
    Nettoie les IDs de gènes pour comparaison.
    Exemple :
    ENSG00000141510.15 -> ENSG00000141510
    """
    if pd.isna(gene):
        return None

    gene = str(gene).strip()

    if gene == "" or gene.lower() == "nan":
        return None

    # Retirer version Ensembl éventuelle
    if gene.startswith("ENSG") and "." in gene:
        gene = gene.split(".")[0]

    return gene


def inspect_expression_matrix(df, dataset_name, gene_col_candidates):
    """
    Inspecte une matrice d'expression au format :
    gène en lignes, échantillons en colonnes.
    """

    print_and_write("\n" + "=" * 80)
    print_and_write(f"INSPECTION : {dataset_name}")
    print_and_write("=" * 80)

    print_and_write(f"Dimensions brutes : {df.shape[0]} lignes x {df.shape[1]} colonnes")

    print_and_write("\nColonnes principales :")
    print_and_write(df.columns[:20].tolist())

    # Trouver colonne gène
    gene_col = None
    for candidate in gene_col_candidates:
        if candidate in df.columns:
            gene_col = candidate
            break

    if gene_col is None:
        # fallback : première colonne
        gene_col = df.columns[0]
        print_and_write(f"\n⚠️ Aucune colonne gène standard trouvée. Utilisation de la première colonne : {gene_col}")
    else:
        print_and_write(f"\nColonne gène détectée : {gene_col}")

    # Nettoyer les gènes
    genes = df[gene_col].apply(clean_gene_name)
    n_genes_total = len(genes)
    n_genes_unique = genes.nunique(dropna=True)
    n_genes_missing = genes.isna().sum()

    print_and_write(f"Nombre de lignes gènes : {n_genes_total}")
    print_and_write(f"Nombre de gènes uniques : {n_genes_unique}")
    print_and_write(f"Nombre de gènes manquants : {n_genes_missing}")

    # Colonnes échantillons
    sample_cols = [col for col in df.columns if col != gene_col]

    print_and_write(f"Nombre d'échantillons/colonnes expression : {len(sample_cols)}")
    print_and_write("\nPremiers échantillons :")
    print_and_write(sample_cols[:20])

    # Types de données
    numeric_sample_cols = []
    for col in sample_cols[:20]:
        if pd.api.types.is_numeric_dtype(df[col]):
            numeric_sample_cols.append(col)

    print_and_write(f"\nColonnes numériques parmi les 20 premiers échantillons : {len(numeric_sample_cols)} / {min(20, len(sample_cols))}")

    # Valeurs manquantes globales approximatives
    missing_total = df.isna().sum().sum()
    print_and_write(f"Nombre total de valeurs manquantes : {missing_total}")

    return {
        "dataset": dataset_name,
        "gene_col": gene_col,
        "n_rows": df.shape[0],
        "n_cols": df.shape[1],
        "n_genes_unique": n_genes_unique,
        "n_samples": len(sample_cols),
        "sample_cols": sample_cols,
        "genes_clean": set(genes.dropna().astype(str))
    }


def compare_with_main_genes(dataset_info, main_genes):
    """
    Compare les gènes d'un dataset avec ceux du dataset principal.
    Attention : si un dataset est en Ensembl et l'autre en gene symbols,
    le nombre commun peut être faible.
    """
    genes = dataset_info["genes_clean"]
    common = sorted(list(genes & main_genes))

    print_and_write("\nComparaison avec le dataset principal ML_READY :")
    print_and_write(f"Gènes dans dataset principal : {len(main_genes)}")
    print_and_write(f"Gènes dans {dataset_info['dataset']} : {len(genes)}")
    print_and_write(f"Gènes communs directs : {len(common)}")

    if len(common) > 0:
        print_and_write("Exemples de gènes communs :")
        print_and_write(common[:20])
    else:
        print_and_write("⚠️ Aucun ou très peu de gènes communs directs. Probable différence d'identifiants : Ensembl vs Gene Symbol.")

    return len(common)


# ============================================================
# MAIN
# ============================================================

def main():
    DATA_PROCESSED.mkdir(exist_ok=True)

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write("RAPPORT D'INSPECTION — DATASETS IMMUNOTHÉRAPIE SUPPLÉMENTAIRES\n")
        f.write("=" * 80 + "\n")

    # Charger les gènes du dataset principal
    if not MAIN_ML_FILE.exists():
        raise FileNotFoundError(f"Fichier principal introuvable : {MAIN_ML_FILE}")

    print_and_write("Lecture du dataset principal ML_READY...")
    df_main_head = pd.read_csv(MAIN_ML_FILE, nrows=5)

    main_metadata_cols = {
        "patient_id",
        "study",
        "response",
        "survival_event"
    }

    main_genes = set([
        col for col in df_main_head.columns
        if col not in main_metadata_cols
    ])

    print_and_write(f"Nombre de gènes dans le dataset principal : {len(main_genes)}")

    datasets = [
        {
            "name": "GSE168204",
            "file": GSE168204_FILE,
            "gene_candidates": ["Gene", "gene", "GENE", "gene_id", "Gene_ID", "GENE_ID", "Unnamed: 0"]
        },
        {
            "name": "GSE244982",
            "file": GSE244982_FILE,
            "gene_candidates": ["Gene", "gene", "GENE", "gene_id", "Gene_ID", "GENE_ID", "Unnamed: 0"]
        },
        {
            "name": "GSE145996",
            "file": GSE145996_FILE,
            "gene_candidates": ["Gene", "gene", "GENE", "gene_id", "Gene_ID", "GENE_ID", "Unnamed: 0"]
        },
    ]

    summary_rows = []

    for dataset in datasets:
        name = dataset["name"]
        file_path = dataset["file"]

        print_and_write("\n\n" + "#" * 80)
        print_and_write(f"DATASET : {name}")
        print_and_write("#" * 80)

        if not file_path.exists():
            print_and_write(f"❌ Fichier introuvable : {file_path}")
            continue

        print_and_write(f"Fichier trouvé : {file_path}")

        try:
            df, sep = read_table_auto(file_path)
            print_and_write(f"Séparateur détecté/utilisé : {sep}")

            info = inspect_expression_matrix(
                df=df,
                dataset_name=name,
                gene_col_candidates=dataset["gene_candidates"]
            )

            n_common = compare_with_main_genes(info, main_genes)

            # Sauvegarder un aperçu
            preview_path = DATA_PROCESSED / f"9_{name}_preview.csv"
            df.head(20).to_csv(preview_path, index=False)
            print_and_write(f"\n✅ Aperçu sauvegardé : {preview_path}")

            summary_rows.append({
                "dataset": name,
                "file": str(file_path),
                "separator": sep,
                "gene_col": info["gene_col"],
                "n_rows": info["n_rows"],
                "n_cols": info["n_cols"],
                "n_samples": info["n_samples"],
                "n_unique_genes": info["n_genes_unique"],
                "n_common_genes_with_main_direct": n_common
            })

        except Exception as e:
            print_and_write(f"❌ Erreur pendant l'inspection de {name} : {e}")

            summary_rows.append({
                "dataset": name,
                "file": str(file_path),
                "error": str(e)
            })

    # Sauvegarde résumé
    df_summary = pd.DataFrame(summary_rows)
    summary_file = DATA_PROCESSED / "9_EXTRA_IMMUNOTHERAPY_inspection_summary.csv"
    df_summary.to_csv(summary_file, index=False)

    print_and_write("\n\n" + "=" * 80)
    print_and_write("RÉSUMÉ GLOBAL")
    print_and_write("=" * 80)
    print_and_write(df_summary)

    print("\n✅ Inspection terminée.")
    print(f"Rapport texte : {REPORT_FILE}")
    print(f"Résumé CSV    : {summary_file}")
    print("Aperçus CSV   : data_processed/9_GSE..._preview.csv")


if __name__ == "__main__":
    main()