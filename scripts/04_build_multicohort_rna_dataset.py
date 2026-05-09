import pandas as pd
import numpy as np

print("--- FUSION RNA MULTICOHORTE CORRIGÉE ---")

# =========================
# FICHIERS
# =========================
FICHIER_RNA_91061 = "2_BULK_RNA_expression_GSE91061_symboles.csv"
FICHIER_RNA_78220 = "2_BULK_RNA_expression_GSE78220_nettoyee.csv"
FICHIER_CLINIQUE = "1_CLINIQUE_MULTICOHORTE_nettoyee.csv"
FICHIER_SCORE = "data/2_BULK_RNA_score_immunitaire_cytolytique.txt"

FICHIER_SORTIE = "3_DATASET_RNA_MULTICOHORTE_FINAL.csv"

# =========================
# 1. CHARGEMENT DES FICHIERS
# =========================
print("Lecture RNA GSE91061...")
df_91061 = pd.read_csv(FICHIER_RNA_91061)

print("Lecture RNA GSE78220...")
df_78220 = pd.read_csv(FICHIER_RNA_78220)

print("Lecture clinique multicohorte...")
df_clin = pd.read_csv(FICHIER_CLINIQUE)

print("Lecture score cytolytique...")
df_score = pd.read_csv(FICHIER_SCORE, sep="\t")
df_score.columns = ["patient_id", "cytolytic_score"]

print("Dimensions GSE91061 :", df_91061.shape)
print("Dimensions GSE78220 :", df_78220.shape)
print("Dimensions clinique :", df_clin.shape)

# =========================
# 2. NETTOYAGE DES IDS
# =========================
for df in [df_91061, df_78220, df_clin, df_score]:
    df["patient_id"] = df["patient_id"].astype(str).str.strip().str.replace('"', '', regex=False)

# =========================
# 3. AJOUT DE LA COHORTE
# =========================
df_91061["study"] = "GSE91061"
df_78220["study"] = "GSE78220"

# =========================
# 4. TROUVER LES GÈNES COMMUNS
# =========================
genes_91061 = set(df_91061.columns) - {"patient_id", "study"}
genes_78220 = set(df_78220.columns) - {"patient_id", "study"}

genes_communs = sorted(list(genes_91061 & genes_78220))

print("Nombre de gènes communs :", len(genes_communs))

if len(genes_communs) == 0:
    raise ValueError("Aucun gène commun trouvé. Vérifie les noms de colonnes.")

# =========================
# 5. GARDER LES MÊMES COLONNES DANS LES DEUX COHORTES
# =========================
colonnes_finales = ["patient_id", "study"] + genes_communs

df_91061_common = df_91061[colonnes_finales].copy()
df_78220_common = df_78220[colonnes_finales].copy()

# =========================
# 6. FUSION RNA VERTICALE
# =========================
df_rna_all = pd.concat([df_91061_common, df_78220_common], ignore_index=True)

print("Dimensions RNA fusionné :", df_rna_all.shape)

# =========================
# 7. AJOUT DES DONNÉES CLINIQUES
# =========================
df_final = pd.merge(
    df_rna_all,
    df_clin,
    on=["patient_id", "study"],
    how="inner"
)

print("Après fusion clinique :", df_final.shape)

# =========================
# 8. AJOUT SCORE CYTOLYTIQUE GSE91061
# =========================
df_final = pd.merge(
    df_final,
    df_score,
    on="patient_id",
    how="left"
)

# =========================
# 9. NETTOYAGE DE LA CIBLE
# =========================
# On enlève les lignes sans réponse claire
df_final = df_final.dropna(subset=["survival_event"]).copy()

# Conversion sécurité
df_final["survival_event"] = df_final["survival_event"].astype(int)

# =========================
# 10. SAUVEGARDE
# =========================
df_final.to_csv(FICHIER_SORTIE, index=False)

print("✅ Fusion terminée.")
print("Fichier créé :", FICHIER_SORTIE)
print("Dimensions finales :", df_final.shape)

print("\nDistribution des cohortes :")
print(df_final["study"].value_counts())

print("\nDistribution de la cible survival_event :")
print(df_final["survival_event"].value_counts())

print("\nColonnes principales :")
print(df_final[["patient_id", "study", "response", "survival_event", "cytolytic_score"]].head())