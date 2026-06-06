from pathlib import Path
import pandas as pd

print("--- 19 : FUSION PSEUDO-BULK SINGLE-CELL ML-READY (BÉTON) ---")

# -----------------------------
# CHEMINS
# -----------------------------
ROOT = Path(".")
DATA_PROCESSED = ROOT / "data_processed"

MENDELEY_FILE = DATA_PROCESSED / "17_MENDELEY_RNA_PSEUDOBULK_LOGCPM.csv"
GSE_EXPR_FILE = DATA_PROCESSED / "18_GSE120575_PSEUDOBULK_TPM_READY.csv"

OUT_FUSION = DATA_PROCESSED / "19_SINGLE_CELL_PSEUDOBULK_ML_READY_FULL.csv"
OUT_SUMMARY = DATA_PROCESSED / "19_SINGLE_CELL_PSEUDOBULK_summary.csv"

# -----------------------------
# FONCTIONS
# -----------------------------
def normalize_sample_id(value):
    if pd.isna(value):
        return None
    return str(value).strip().replace(" ", "_")

# -----------------------------
# SCRIPT PRINCIPAL
# -----------------------------
def main():
    # 1. Charger Mendeley
    df_mendeley = pd.read_csv(MENDELEY_FILE)
    df_mendeley['sample_id'] = df_mendeley['sample_id'].apply(normalize_sample_id)
    df_mendeley['dataset_source'] = 'MENDELEY'
    print("MENDELEY :", df_mendeley.shape)

    # 2. Charger GSE pseudo-bulk (Déjà au bon format échantillon x gènes !)
    df_gse_ready = pd.read_csv(GSE_EXPR_FILE)
    df_gse_ready['sample_id'] = df_gse_ready['sample_id'].apply(normalize_sample_id)
    df_gse_ready['dataset_source'] = 'GSE120575'
    print("GSE120575 :", df_gse_ready.shape)

    # 3. Identifier proprement les gènes (en excluant les colonnes de métadonnées)
    metadata_cols_mendeley = ['study', 'sample_id', 'patient_id', 'response', 'response_binary', 'n_cells', 'expression_format', 'dataset_role', 'dataset_source']
    gene_cols_mendeley = [c for c in df_mendeley.columns if c not in metadata_cols_mendeley]
    
    metadata_cols_gse = ['sample_id', 'patient_id', 'response', 'response_binary', 'dataset_role', 'notes', 'dataset_source']
    gene_cols_gse = [c for c in df_gse_ready.columns if c not in metadata_cols_gse]
    
    # Trouver l'intersection des gènes
    common_genes = sorted(list(set(gene_cols_mendeley) & set(gene_cols_gse)))
    print("Nombre de gènes communs :", len(common_genes))

    # 4. Sélection des métadonnées communes + gènes communs
    metadata_to_keep = ['sample_id', 'response', 'response_binary', 'dataset_source']

    df_mendeley_sel = df_mendeley[metadata_to_keep + common_genes]
    df_gse_sel = df_gse_ready[metadata_to_keep + common_genes]

    # 5. Fusion finale par concaténation verticale
    df_fusion = pd.concat([df_mendeley_sel, df_gse_sel], axis=0, ignore_index=True)
    df_fusion.to_csv(OUT_FUSION, index=False)
    print("✅ Fusion sauvegardée :", OUT_FUSION)

    # 6. Résumé final
    summary = pd.DataFrame([
        {'dataset_source': 'MENDELEY', 'n_samples': df_mendeley_sel.shape[0], 'n_genes': len(common_genes)},
        {'dataset_source': 'GSE120575', 'n_samples': df_gse_sel.shape[0], 'n_genes': len(common_genes)},
        {'dataset_source': 'TOTAL_FUSION', 'n_samples': df_fusion.shape[0], 'n_genes': len(common_genes)}
    ])
    summary.to_csv(OUT_SUMMARY, index=False)
    print("✅ Résumé sauvegardé :", OUT_SUMMARY)
    print(summary)

if __name__ == "__main__":
    main()