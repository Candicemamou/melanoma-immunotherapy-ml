import pandas as pd

print("--- CONVERSION GSE91061 : ENTREZ ID -> GENE SYMBOL ---")

# =========================
# FICHIERS
# =========================
FICHIER_RNA = "2_BULK_RNA_expression_GSE91061_nettoyee.csv"
FICHIER_MAPPING = "data/4_GRAPH_STRING_mapping_brut.txt.gz"
FICHIER_SORTIE = "2_BULK_RNA_expression_GSE91061_symboles.csv"

# =========================
# 1. CHARGEMENT RNA GSE91061
# =========================
print("Lecture du fichier RNA GSE91061...")
df_rna = pd.read_csv(FICHIER_RNA)

print("Dimensions RNA initiales :", df_rna.shape)

# Nettoyage patient_id
df_rna["patient_id"] = df_rna["patient_id"].astype(str).str.strip()

# Les colonnes de gènes sont toutes les colonnes sauf patient_id
gene_cols = [col for col in df_rna.columns if col != "patient_id"]

# On force les noms de colonnes en string
df_rna.columns = df_rna.columns.astype(str)
gene_cols = [str(col) for col in gene_cols]

# =========================
# 2. CHARGEMENT MAPPING STRING
# =========================
print("Lecture du fichier de mapping STRING...")
mapping = pd.read_csv(FICHIER_MAPPING, sep="\t")

print("Colonnes mapping :", mapping.columns.tolist())

# =========================
# 3. EXTRACTION ENTREZ -> PROTEIN STRING
# =========================
mapping_entrez = mapping[mapping["source"] == "Ensembl_HGNC_entrez_id"].copy()
mapping_entrez = mapping_entrez[["alias", "#string_protein_id"]]
mapping_entrez = mapping_entrez.rename(columns={"alias": "entrez_id"})

mapping_entrez["entrez_id"] = mapping_entrez["entrez_id"].astype(str)

# =========================
# 4. EXTRACTION PROTEIN STRING -> GENE SYMBOL
# =========================
mapping_symbol = mapping[mapping["source"] == "Ensembl_HGNC"].copy()
mapping_symbol = mapping_symbol[["#string_protein_id", "alias"]]
mapping_symbol = mapping_symbol.rename(columns={"alias": "gene_symbol"})

# =========================
# 5. FUSION DES DEUX MAPPINGS
# =========================
mapping_final = mapping_entrez.merge(
    mapping_symbol,
    on="#string_protein_id",
    how="inner"
)

mapping_final = mapping_final[["entrez_id", "gene_symbol"]].drop_duplicates()

print("Nombre de correspondances Entrez -> Symbol :", len(mapping_final))

# Dictionnaire : "1234" -> "GENE"
dict_mapping = dict(zip(mapping_final["entrez_id"], mapping_final["gene_symbol"]))

# =========================
# 6. RENOMMAGE DES COLONNES RNA
# =========================
rename_dict = {}

for col in gene_cols:
    if col in dict_mapping:
        rename_dict[col] = dict_mapping[col]

print("Nombre de gènes renommés :", len(rename_dict))
print("Nombre de gènes non retrouvés :", len(gene_cols) - len(rename_dict))

df_rna_renamed = df_rna.rename(columns=rename_dict)

# =========================
# 7. SUPPRESSION DES COLONNES NON CONVERTIES
# =========================
colonnes_converties = ["patient_id"] + list(rename_dict.values())
df_rna_renamed = df_rna_renamed[colonnes_converties]

# =========================
# 8. GÉRER LES DOUBLONS DE SYMBOLS
# =========================
# Certains Entrez IDs peuvent pointer vers le même symbole.
# On moyenne les colonnes dupliquées.
print("Suppression / moyenne des doublons de noms de gènes...")

patient_ids = df_rna_renamed["patient_id"]
genes_only = df_rna_renamed.drop(columns=["patient_id"])

genes_only = genes_only.groupby(level=0, axis=1).mean()

df_final = pd.concat([patient_ids, genes_only], axis=1)

# =========================
# 9. SAUVEGARDE
# =========================
df_final.to_csv(FICHIER_SORTIE, index=False)

print("✅ Conversion terminée.")
print("Fichier créé :", FICHIER_SORTIE)
print("Dimensions finales :", df_final.shape)
print(df_final.head())