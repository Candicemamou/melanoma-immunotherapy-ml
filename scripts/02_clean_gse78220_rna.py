import pandas as pd

print("--- NETTOYAGE RNA GSE78220 ---")

# 1. Lire le fichier Excel
df = pd.read_excel(
    "data/2_BULK_RNA_expression_GSE78220_brut.xlsx",
    sheet_name="FPKM"
)

print("Dimensions brutes :", df.shape)

# 2. Mettre les gènes en index
df = df.set_index("Gene")

# 3. Transposer : patients en lignes, gènes en colonnes
df = df.T

# 4. Remettre patient_id comme colonne
df.reset_index(inplace=True)
df.rename(columns={"index": "patient_id"}, inplace=True)

# 5. Nettoyer les IDs patients : Pt1.baseline -> Pt1
df["patient_id"] = (
    df["patient_id"]
    .astype(str)
    .str.replace(".baseline", "", regex=False)
    .str.replace(".OnTx", "", regex=False)
    .str.strip()
)

# 6. Nettoyer les noms de colonnes
df.columns = df.columns.astype(str).str.strip()

# 7. Sauvegarder
df.to_csv("2_BULK_RNA_expression_GSE78220_nettoyee.csv", index=False)

print("✅ Fichier créé : 2_BULK_RNA_expression_GSE78220_nettoyee.csv")
print("Dimensions finales :", df.shape)
print(df.head())