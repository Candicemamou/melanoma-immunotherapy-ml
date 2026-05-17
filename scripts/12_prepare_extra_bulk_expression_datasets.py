from pathlib import Path
import pandas as pd
import numpy as np

print("--- 12 : PRÉPARATION DES DATASETS BULK SUPPLÉMENTAIRES ---")

# ============================================================
# CHEMINS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data_raw"
DATA_PROCESSED = ROOT / "data_processed"

DATA_PROCESSED.mkdir(exist_ok=True)

# Expressions brutes
GSE168204_EXPR = DATA_RAW / "6_IMMUNOTHERAPY_GSE168204_counts_raw.csv.gz"
GSE244982_EXPR = DATA_RAW / "6_IMMUNOTHERAPY_GSE244982_bulkRNAseq_raw.txt.gz"
GSE196434_EXPR = DATA_RAW / "5_TARGETED_THERAPY_GSE196434_raw.txt.gz"

# Métadonnées standardisées générées par le script 11
GSE168204_META = DATA_PROCESSED / "11_METADATA_GSE168204_clean.csv"
GSE244982_META = DATA_PROCESSED / "11_METADATA_GSE244982_clean.csv"
GSE196434_META = DATA_PROCESSED / "11_METADATA_GSE196434_clean.csv"

# Mapping STRING pour essayer de convertir Ensembl -> gene symbol si nécessaire
STRING_MAPPING = DATA_RAW / "4_GRAPH_STRING_mapping_brut.txt.gz"

# Sorties
OUT_GSE168204 = DATA_PROCESSED / "12_GSE168204_bulk_clean.csv"
OUT_GSE244982 = DATA_PROCESSED / "12_GSE244982_bulk_clean.csv"
OUT_GSE196434 = DATA_PROCESSED / "12_GSE196434_bulk_clean.csv"

SUMMARY_FILE = DATA_PROCESSED / "12_EXTRA_BULK_summary.csv"


# ============================================================
# OUTILS
# ============================================================

def clean_text(value):
    if pd.isna(value):
        return None
    value = str(value).replace('"', "").strip()
    return value if value != "" else None


def clean_gene_name(value):
    """
    Nettoie les noms de gènes.
    Si c'est un Ensembl ID avec version :
    ENSG00000141510.15 -> ENSG00000141510
    """
    value = clean_text(value)

    if value is None:
        return None

    if value.startswith("ENSG") and "." in value:
        value = value.split(".")[0]

    return value


def read_expression_table(path):
    """
    Lecture robuste csv / txt / gz.
    Essaie plusieurs séparateurs.
    """
    if not path.exists():
        raise FileNotFoundError(f"Fichier introuvable : {path}")

    separators = [",", "\t", ";"]

    for sep in separators:
        try:
            df = pd.read_csv(
                path,
                sep=sep,
                compression="infer",
                low_memory=False
            )

            if df.shape[1] > 1:
                return df, sep

        except Exception:
            pass

    # Dernier essai en auto
    df = pd.read_csv(
        path,
        sep=None,
        engine="python",
        compression="infer"
    )
    return df, "auto"


def guess_gene_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col

    # fallback : première colonne
    return df.columns[0]


def load_ensembl_to_symbol_mapping():
    """
    Essaie de créer un mapping Ensembl gene ID -> gene symbol à partir du fichier STRING.
    Si ça échoue, retourne un dictionnaire vide.
    """
    if not STRING_MAPPING.exists():
        print("⚠️ Mapping STRING absent. Conversion Ensembl -> Symbol impossible.")
        return {}

    try:
        mapping = pd.read_csv(STRING_MAPPING, sep="\t", compression="infer", low_memory=False)

        if "#string_protein_id" not in mapping.columns or "alias" not in mapping.columns or "source" not in mapping.columns:
            print("⚠️ Colonnes attendues absentes dans le mapping STRING.")
            return {}

        mapping["alias"] = mapping["alias"].astype(str).str.strip()
        mapping["source"] = mapping["source"].astype(str).str.strip()

        # Candidats Ensembl gene IDs
        ens = mapping[mapping["alias"].str.startswith("ENSG", na=False)].copy()
        ens["ensembl_id"] = ens["alias"].str.split(".").str[0]

        # Candidats gene symbols : on exclut les alias Ensembl et Entrez numériques
        sym = mapping[
            (~mapping["alias"].str.startswith("ENSG", na=False)) &
            (~mapping["alias"].str.fullmatch(r"\d+", na=False))
        ].copy()

        # On privilégie les sources HGNC quand elles existent
        sym_hgnc = sym[sym["source"].str.contains("HGNC", case=False, na=False)].copy()

        if len(sym_hgnc) > 0:
            sym = sym_hgnc

        ens = ens[["#string_protein_id", "ensembl_id"]].drop_duplicates()
        sym = sym[["#string_protein_id", "alias"]].rename(columns={"alias": "gene_symbol"}).drop_duplicates()

        merged = ens.merge(sym, on="#string_protein_id", how="inner")
        merged = merged.dropna(subset=["ensembl_id", "gene_symbol"])
        merged = merged.drop_duplicates(subset=["ensembl_id"], keep="first")

        mapping_dict = dict(zip(merged["ensembl_id"], merged["gene_symbol"]))

        print(f"Mapping Ensembl -> Symbol chargé : {len(mapping_dict)} correspondances")
        return mapping_dict

    except Exception as e:
        print(f"⚠️ Impossible de charger le mapping Ensembl -> Symbol : {e}")
        return {}


def prepare_expression_matrix(
    expr_path,
    meta_path,
    output_path,
    dataset_name,
    gene_col_candidates,
    convert_ensembl=False
):
    """
    Convertit une matrice expression gènes x échantillons en :
    échantillons x gènes + metadata.
    """

    print("\n" + "=" * 80)
    print(f"PRÉPARATION : {dataset_name}")
    print("=" * 80)

    if not expr_path.exists():
        raise FileNotFoundError(f"Expression introuvable : {expr_path}")

    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata introuvable : {meta_path}")

    print("Lecture expression...")
    df_expr, sep = read_expression_table(expr_path)

    print("Lecture metadata...")
    df_meta = pd.read_csv(meta_path)

    print("Séparateur utilisé expression :", sep)
    print("Dimensions expression brute :", df_expr.shape)
    print("Dimensions metadata :", df_meta.shape)

    gene_col = guess_gene_column(df_expr, gene_col_candidates)
    print("Colonne gène utilisée :", gene_col)

    # Nettoyer gènes
    df_expr[gene_col] = df_expr[gene_col].apply(clean_gene_name)
    df_expr = df_expr.dropna(subset=[gene_col]).copy()

    # Conversion Ensembl -> gene symbol si demandé
    if convert_ensembl:
        mapping_dict = load_ensembl_to_symbol_mapping()

        if len(mapping_dict) > 0:
            before_unique = df_expr[gene_col].nunique()

            df_expr[gene_col] = df_expr[gene_col].map(
                lambda g: mapping_dict.get(g, g)
            )

            after_unique = df_expr[gene_col].nunique()

            print("Conversion Ensembl -> Symbol appliquée.")
            print("Gènes uniques avant :", before_unique)
            print("Gènes uniques après :", after_unique)
        else:
            print("⚠️ Conversion Ensembl demandée mais mapping vide. IDs conservés.")

    # Colonnes échantillons
    sample_cols = [col for col in df_expr.columns if col != gene_col]

    print("Nombre de colonnes échantillons expression :", len(sample_cols))
    print("Premiers échantillons expression :", sample_cols[:10])

    # Conversion numérique
    for col in sample_cols:
        df_expr[col] = pd.to_numeric(df_expr[col], errors="coerce")

    # Moyenne des doublons de gènes
    df_expr_grouped = df_expr.groupby(gene_col, as_index=True).mean(numeric_only=True)

    print("Nombre de gènes après moyenne doublons :", df_expr_grouped.shape[0])

    # Transposition : samples en lignes
    df_t = df_expr_grouped.T.reset_index()
    df_t = df_t.rename(columns={"index": "sample_id"})

    def normalize_sample_id_for_merge(value):
        """
        Normalise les sample_id pour corriger les différences simples :
        - tirets vs underscores
        - suffixe _bam
        - suffixe .bam
        """
        if pd.isna(value):
            return None

        value = str(value).strip()

        value = value.replace(".bam", "")
        value = value.replace("_bam", "")
        value = value.replace("-", "_")

        return value


    df_t["sample_id_original"] = df_t["sample_id"].astype(str).str.strip()
    df_meta["sample_id_original"] = df_meta["sample_id"].astype(str).str.strip()

    df_t["sample_id_merge"] = df_t["sample_id_original"].apply(normalize_sample_id_for_merge)
    df_meta["sample_id_merge"] = df_meta["sample_id_original"].apply(normalize_sample_id_for_merge)

        # Fusion expression + metadata
    df_clean = df_meta.merge(df_t, on="sample_id_merge", how="inner")

    # On garde le sample_id original côté metadata comme identifiant principal
    if "sample_id_original_x" in df_clean.columns:
        df_clean["sample_id"] = df_clean["sample_id_original_x"]

    # Supprimer colonnes techniques de merge inutiles
    cols_to_drop = [
        "sample_id_original_x",
        "sample_id_original_y",
        "sample_id_y"
    ]

    cols_to_drop = [col for col in cols_to_drop if col in df_clean.columns]
    df_clean = df_clean.drop(columns=cols_to_drop)

    print("Dimensions après fusion expression + metadata :", df_clean.shape)

    missing_samples = sorted(set(df_meta["sample_id_merge"]) - set(df_t["sample_id_merge"]))
    extra_expr_samples = sorted(set(df_t["sample_id_merge"]) - set(df_meta["sample_id_merge"]))

    print("Échantillons metadata absents de l'expression :", len(missing_samples))
    if len(missing_samples) > 0:
        print(missing_samples[:20])

    print("Échantillons expression absents de la metadata :", len(extra_expr_samples))
    if len(extra_expr_samples) > 0:
        print(extra_expr_samples[:20])

    # Sécurité : supprimer lignes sans response_binary si la colonne existe
    if "response_binary" in df_clean.columns:
        df_clean["response_binary"] = pd.to_numeric(df_clean["response_binary"], errors="coerce")

    # Sauvegarde
    df_clean.to_csv(output_path, index=False)

    print(f"✅ Fichier clean créé : {output_path}")

    n_gene_cols = df_clean.shape[1] - len(df_meta.columns)

    print("Nombre d'échantillons clean :", df_clean.shape[0])
    print("Nombre de gènes clean :", n_gene_cols)

    if "response_binary" in df_clean.columns:
        print("Répartition response_binary :")
        print(df_clean["response_binary"].value_counts(dropna=False))

    if "timepoint" in df_clean.columns:
        print("Répartition timepoint :")
        print(df_clean["timepoint"].value_counts(dropna=False))

    if "treatment" in df_clean.columns:
        print("Répartition treatment :")
        print(df_clean["treatment"].value_counts(dropna=False))

    return {
        "dataset": dataset_name,
        "output_file": str(output_path),
        "n_expression_rows_raw": df_expr.shape[0],
        "n_expression_cols_raw": df_expr.shape[1],
        "n_metadata_samples": df_meta.shape[0],
        "n_clean_samples": df_clean.shape[0],
        "n_clean_genes": n_gene_cols,
        "n_missing_samples_from_expression": len(missing_samples),
        "n_expression_samples_without_metadata": len(extra_expr_samples),
    }


# ============================================================
# MAIN
# ============================================================

def main():
    summary_rows = []

    # GSE168204 : expression avec gene symbols, samples MGH...
    summary_rows.append(
        prepare_expression_matrix(
            expr_path=GSE168204_EXPR,
            meta_path=GSE168204_META,
            output_path=OUT_GSE168204,
            dataset_name="GSE168204",
            gene_col_candidates=["Gene", "gene", "GENE", "Gene Name", "Unnamed: 0"],
            convert_ensembl=False
        )
    )

    # GSE244982 : expression avec gene symbols, samples Pat...
    summary_rows.append(
        prepare_expression_matrix(
            expr_path=GSE244982_EXPR,
            meta_path=GSE244982_META,
            output_path=OUT_GSE244982,
            dataset_name="GSE244982",
            gene_col_candidates=["Gene", "gene", "GENE", "Gene Name", "Unnamed: 0"],
            convert_ensembl=False
        )
    )

    # GSE196434 : expression souvent en Ensembl IDs, on tente conversion
    summary_rows.append(
        prepare_expression_matrix(
            expr_path=GSE196434_EXPR,
            meta_path=GSE196434_META,
            output_path=OUT_GSE196434,
            dataset_name="GSE196434",
            gene_col_candidates=["GENE_ID", "Gene", "gene", "GENE", "Gene Name", "Unnamed: 0"],
            convert_ensembl=True
        )
    )

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(SUMMARY_FILE, index=False)

    print("\n" + "=" * 80)
    print("RÉSUMÉ GLOBAL SCRIPT 12")
    print("=" * 80)
    print(df_summary)

    print(f"\n✅ Résumé sauvegardé : {SUMMARY_FILE}")


if __name__ == "__main__":
    main()