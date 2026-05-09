


import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder

# =============================================================================
# CONFIGURATION
# =============================================================================

DATA_DIR   = "."        # dossier contenant tes fichiers nettoyés
OUTPUT_DIR = "ml_ready" # dossier de sortie

os.makedirs(OUTPUT_DIR, exist_ok=True)

PATH_CLINIQUE = os.path.join(DATA_DIR, "1_CLINIQUE_donnees_principales_nettoyees.csv")
PATH_RNA      = os.path.join(DATA_DIR, "2_BULK_RNA_expression_nettoyee.csv")

LABEL_COL     = "survival_event"   # colonne cible dans le fichier clinique
RANDOM_STATE  = 42
TEST_SIZE     = 0.20
VAL_SIZE      = 0.15               # fraction du train set mise de côté pour val

META_COLS_RNA = ["patient_id", "treatment_status", "raw_sample_id"]


# =============================================================================
# UTILITAIRES
# =============================================================================

def log2_norm(df_genes: pd.DataFrame) -> pd.DataFrame:
    """Log2(x + 1) sur toutes les colonnes (suppose que df ne contient que des gènes)."""
    return np.log2(df_genes + 1)


def filter_by_variance(df_genes: pd.DataFrame, top_n: int = 5000) -> pd.DataFrame:
    """
    Garde les top_n gènes les plus variables (variance inter-échantillons).
    Réduit la dimensionnalité avant PCA sans perte d'information utile.
    """
    variances = df_genes.var(axis=0)
    top_genes = variances.nlargest(top_n).index
    print(f"  Sélection top {top_n} gènes les plus variables (sur {df_genes.shape[1]} total)")
    return df_genes[top_genes]


def apply_pca(df_genes: pd.DataFrame, n_components: int = 50,
              prefix: str = "PC") -> pd.DataFrame:
    """
    PCA → réduit les gènes en composantes principales.
    Retourne un DataFrame avec les composantes comme colonnes.
    """
    pca = PCA(n_components=n_components, random_state=RANDOM_STATE)
    coords = pca.fit_transform(df_genes.values)
    explained = pca.explained_variance_ratio_.cumsum()
    print(f"  PCA {n_components} composantes → variance expliquée cumulée : "
          f"{explained[-1]*100:.1f}%")
    cols = [f"{prefix}{i+1}" for i in range(n_components)]
    return pd.DataFrame(coords, index=df_genes.index, columns=cols), pca


def train_val_test_split(X: pd.DataFrame, y: pd.Series):
    """
    Découpe stratifiée en train / val / test.
    test  = TEST_SIZE  du total
    val   = VAL_SIZE   du reste (train+val)
    """
    X_tv, X_test, y_tv, y_test = train_test_split(
        X, y, test_size=TEST_SIZE, stratify=y, random_state=RANDOM_STATE
    )
    X_train, X_val, y_train, y_val = train_test_split(
        X_tv, y_tv, test_size=VAL_SIZE/(1-TEST_SIZE),
        stratify=y_tv, random_state=RANDOM_STATE
    )
    print(f"  Train : {X_train.shape[0]} | Val : {X_val.shape[0]} | Test : {X_test.shape[0]}")
    return X_train, X_val, X_test, y_train, y_val, y_test


def save_split(X_train, X_val, X_test, y_train, y_val, y_test, prefix: str):
    """Sauvegarde les 6 fichiers train/val/test X et y."""
    for split_name, X, y in [("train", X_train, y_train),
                              ("val",   X_val,   y_val),
                              ("test",  X_test,  y_test)]:
        X.to_csv(os.path.join(OUTPUT_DIR, f"{prefix}_X_{split_name}.csv"))
        y.to_csv(os.path.join(OUTPUT_DIR, f"{prefix}_y_{split_name}.csv"))
    print(f"  ✅ Fichiers {prefix}_X/y_train/val/test.csv sauvegardés dans '{OUTPUT_DIR}/'")


# =============================================================================
# ÉTAPE 1 — Chargement et normalisation RNA
# =============================================================================

def etape1_charger_normaliser_rna():
    print("\n" + "="*60)
    print("ÉTAPE 1 — Chargement et normalisation log2 des données RNA")
    print("="*60)

    df_rna = pd.read_csv(PATH_RNA)
    print(f"  Dimensions brutes : {df_rna.shape}  ({df_rna.shape[0]} échantillons × colonnes)")

    # Séparer métadonnées et expression
    meta = df_rna[META_COLS_RNA].copy()
    genes = df_rna.drop(columns=META_COLS_RNA)

    print(f"  Colonnes gènes : {genes.shape[1]}")
    print(f"  Patients uniques : {meta['patient_id'].nunique()}")
    print(f"  Distribution Pre/On : {meta['treatment_status'].value_counts().to_dict()}")

    # Vérifier valeurs négatives ou NaN résiduels
    n_neg = (genes < 0).sum().sum()
    n_nan = genes.isna().sum().sum()
    print(f"  Valeurs négatives : {n_neg} | NaN résiduels : {n_nan}")

    if n_neg > 0:
        print("  ⚠️  Valeurs négatives remplacées par 0 avant log2")
        genes = genes.clip(lower=0)

    # Log2-normalisation
    genes_log = log2_norm(genes)

    # Remettre les métadonnées
    df_rna_norm = pd.concat([meta.reset_index(drop=True),
                             genes_log.reset_index(drop=True)], axis=1)

    out = os.path.join(OUTPUT_DIR, "RNA_log2_normalise.csv")
    df_rna_norm.to_csv(out, index=False)
    print(f"  ✅ Sauvegardé : {out}")

    return df_rna_norm, meta, genes_log


# =============================================================================
# ÉTAPE 2 — Réduction dimensionnelle RNA
# =============================================================================

def etape2_reduction_dimensionnelle(genes_log: pd.DataFrame, meta: pd.DataFrame):
    print("\n" + "="*60)
    print("ÉTAPE 2 — Réduction dimensionnelle RNA (Variance → PCA)")
    print("="*60)

    # Sélection des gènes les plus variables
    genes_top = filter_by_variance(genes_log, top_n=5000)

    # PCA → 50 composantes
    genes_pca, pca_model = apply_pca(genes_top, n_components=50, prefix="PC")
    genes_pca.index = meta["raw_sample_id"].values

    # Rajouter les métadonnées
    df_pca = pd.concat([meta.reset_index(drop=True), genes_pca.reset_index(drop=True)], axis=1)

    out = os.path.join(OUTPUT_DIR, "RNA_PCA50.csv")
    df_pca.to_csv(out, index=False)
    print(f"  ✅ Sauvegardé : {out}")

    return genes_pca, genes_top, pca_model


# =============================================================================
# ÉTAPE 3 — Feature engineering : delta Pre → On
# =============================================================================

def etape3_delta_pre_on(genes_log: pd.DataFrame, meta: pd.DataFrame):
    """
    Pour chaque patient ayant une mesure Pre ET On :
    calcule delta = On - Pre pour chaque gène.
    Cela capture le changement d'expression dû au traitement.
    """
    print("\n" + "="*60)
    print("ÉTAPE 3 — Feature engineering : Δ expression (On - Pre)")
    print("="*60)

    df = pd.concat([meta.reset_index(drop=True), genes_log.reset_index(drop=True)], axis=1)

    pre = df[df["treatment_status"] == "Pre"].set_index("patient_id")
    on  = df[df["treatment_status"] == "On"].set_index("patient_id")

    gene_cols = [c for c in pre.columns if c not in META_COLS_RNA]

    # Patients avec les deux mesures
    patients_both = pre.index.intersection(on.index)
    print(f"  Patients avec Pre ET On : {len(patients_both)} / {meta['patient_id'].nunique()}")

    if len(patients_both) == 0:
        print("  ⚠️  Aucun patient avec les deux mesures — delta non calculable.")
        return None

    delta = on.loc[patients_both, gene_cols] - pre.loc[patients_both, gene_cols]
    delta.index.name = "patient_id"
    delta.columns = ["delta_" + c for c in delta.columns]

    # Aussi garder les valeurs Pre seules (baseline) pour les patients avec Pre uniquement
    df_pre_only = pre.loc[~pre.index.isin(patients_both), gene_cols].copy()
    df_pre_only.columns = ["pre_" + c for c in df_pre_only.columns]

    out_delta = os.path.join(OUTPUT_DIR, "RNA_delta_On_minus_Pre.csv")
    delta.to_csv(out_delta)
    print(f"  ✅ Delta sauvegardé : {out_delta}  {delta.shape}")

    # Top gènes les plus modifiés en moyenne
    top_delta_genes = delta.abs().mean().nlargest(20)
    print(f"\n  Top 10 gènes les plus modifiés Pre→On :")
    for gene, val in top_delta_genes.head(10).items():
        print(f"    {gene.replace('delta_',''):>12}  Δ moyen absolu = {val:.3f}")

    return delta


# =============================================================================
# ÉTAPE 4 — Préparation dataset ML Clinique (TCGA)
# =============================================================================

def etape4_dataset_ml_clinique():
    print("\n" + "="*60)
    print("ÉTAPE 4 — Préparation dataset ML — Clinique TCGA")
    print("="*60)

    df_clin = pd.read_csv(PATH_CLINIQUE)
    print(f"  Dimensions : {df_clin.shape}")
    print(f"  Label '{LABEL_COL}' : {df_clin[LABEL_COL].value_counts().to_dict()}")

    # Vérification : les données sont déjà normalisées (StandardScaler)
    means = df_clin.drop(columns=[LABEL_COL]).mean()
    stds  = df_clin.drop(columns=[LABEL_COL]).std()
    already_scaled = (means.abs() < 0.1).all() and (stds.between(0.9, 1.1)).all()
    if already_scaled:
        print("  ✅ Données déjà normalisées (StandardScaler détecté) — aucune action")
    else:
        print("  ⚠️  Données non normalisées — application StandardScaler")
        feat_cols = [c for c in df_clin.columns if c != LABEL_COL]
        scaler = StandardScaler()
        df_clin[feat_cols] = scaler.fit_transform(df_clin[feat_cols])

    # Vérification NaN résiduels
    n_nan = df_clin.isna().sum().sum()
    if n_nan > 0:
        print(f"  ⚠️  {n_nan} NaN résiduels → imputation par médiane")
        df_clin.fillna(df_clin.median(numeric_only=True), inplace=True)

    X = df_clin.drop(columns=[LABEL_COL])
    y = df_clin[LABEL_COL]

    print(f"\n  Features : {X.shape[1]}  |  Échantillons : {X.shape[0]}")
    print(f"  Classe 0 (survie) : {(y==0).sum()}  |  Classe 1 (décès) : {(y==1).sum()}")

    X_train, X_val, X_test, y_train, y_val, y_test = train_val_test_split(X, y)
    save_split(X_train, X_val, X_test, y_train, y_val, y_test, prefix="CLINIQUE")

    return X, y


# =============================================================================
# ÉTAPE 5 — Préparation dataset ML RNA (GSE91061)
# =============================================================================

def etape5_dataset_ml_rna(genes_pca: pd.DataFrame, meta: pd.DataFrame,
                           delta: pd.DataFrame = None):
    """
    Deux variantes du dataset RNA pour le ML :
      A) PCA 50 composantes, par échantillon (Pre + On séparés)
      B) Delta Pre→On (un vecteur par patient)

    Le label ici est treatment_status (Pre=0 / On=1) —
    utile pour des analyses de changement d'expression.
    ⚠️  Pour un vrai label clinique (réponse immuno), il faudra
    fusionner avec les métadonnées GSE91061 (responder/non-responder).
    """
    print("\n" + "="*60)
    print("ÉTAPE 5 — Préparation dataset ML — RNA GSE91061")
    print("="*60)

    # --- Variante A : PCA par échantillon, label = Pre/On ---
    print("\n  [5A] PCA par échantillon — label = Pre/On (0/1)")
    le = LabelEncoder()
    y_preOn = pd.Series(
        le.fit_transform(meta["treatment_status"].values),
        name="treatment_status_encoded"
    )
    print(f"  Encodage : {dict(zip(le.classes_, le.transform(le.classes_)))}")
    print(f"  Échantillons : {genes_pca.shape[0]}  |  Features PCA : {genes_pca.shape[1]}")

    X_train, X_val, X_test, y_train, y_val, y_test = train_val_test_split(genes_pca, y_preOn)
    save_split(X_train, X_val, X_test, y_train, y_val, y_test, prefix="RNA_PCA_PreOn")

    # --- Variante B : Delta par patient ---
    if delta is not None:
        print("\n  [5B] Delta Pre→On par patient")
        print(f"  Patients : {delta.shape[0]}  |  Features delta : {delta.shape[1]}")

        # Réduction variance sur le delta
        top_delta = filter_by_variance(delta, top_n=min(2000, delta.shape[1]))
        delta_pca, _ = apply_pca(top_delta, n_components=min(30, delta.shape[0]-1), prefix="DPC")

        # Note : sans label de réponse au traitement ici, on sauvegarde X seulement
        # Le label réel (responder) sera ajouté à l'étape 6
        out = os.path.join(OUTPUT_DIR, "RNA_delta_PCA30.csv")
        delta_pca.to_csv(out)
        print(f"  ✅ Sauvegardé : {out}")

    return genes_pca, y_preOn


# =============================================================================
# ÉTAPE 6 — Dataset combiné Clinique + RNA
# =============================================================================

def etape6_dataset_combine(X_clin: pd.DataFrame, y_clin: pd.Series,
                            genes_pca: pd.DataFrame, meta: pd.DataFrame):
    """
    Le fichier clinique (TCGA, 436 patients) et le fichier RNA (GSE91061, 65 patients)
    ne partagent PAS les mêmes patients — les IDs sont incompatibles.

    Stratégie recommandée pour le ML :
      → Entraîner sur TCGA (clinique) → valider sur GSE91061 (RNA)
      → Ou entraîner deux modèles séparés et comparer

    Ce bloc prépare un dataset RNA avec normalisation StandardScaler
    indépendante, prêt pour une validation croisée inter-cohortes.
    """
    print("\n" + "="*60)
    print("ÉTAPE 6 — Dataset combiné / inter-cohortes")
    print("="*60)

    print("  ℹ️  Les cohortes TCGA et GSE91061 sont indépendantes.")
    print("  → Stratégie : entraîner sur TCGA, évaluer la transférabilité sur GSE91061.")
    print()

    # Normaliser les PCA RNA avec un nouveau scaler (indépendant du clinique)
    scaler_rna = StandardScaler()
    pca_scaled = pd.DataFrame(
        scaler_rna.fit_transform(genes_pca.values),
        index=genes_pca.index,
        columns=genes_pca.columns
    )

    # Résumé des datasets disponibles
    print("  Récapitulatif datasets ML prêts :")
    print(f"    CLINIQUE (TCGA)    : {X_clin.shape[0]} patients × {X_clin.shape[1]} features")
    print(f"    RNA PCA (GSE91061) : {pca_scaled.shape[0]} échantillons × {pca_scaled.shape[1]} features")
    print()
    print("  Fichiers de sortie générés dans le dossier 'ml_ready/' :")
    print("    CLINIQUE_X_train.csv / _val.csv / _test.csv")
    print("    CLINIQUE_y_train.csv / _val.csv / _test.csv")
    print("    RNA_PCA_PreOn_X_train.csv / _val.csv / _test.csv")
    print("    RNA_PCA_PreOn_y_train.csv / _val.csv / _test.csv")
    print("    RNA_delta_PCA30.csv")
    print("    RNA_log2_normalise.csv")
    print("    RNA_PCA50.csv")

    out = os.path.join(OUTPUT_DIR, "RNA_PCA50_scaled.csv")
    pca_scaled_full = pd.concat(
        [meta[["patient_id", "treatment_status"]].reset_index(drop=True),
         pca_scaled.reset_index(drop=True)], axis=1
    )
    pca_scaled_full.to_csv(out, index=False)
    print(f"\n  ✅ RNA normalisé final sauvegardé : {out}")


# =============================================================================
# ÉTAPE 7 — Résumé et vérifications finales
# =============================================================================

def etape7_verification_finale():
    print("\n" + "="*60)
    print("ÉTAPE 7 — Vérification des fichiers générés")
    print("="*60)

    fichiers = os.listdir(OUTPUT_DIR)
    total_size = 0
    for f in sorted(fichiers):
        path = os.path.join(OUTPUT_DIR, f)
        size_kb = os.path.getsize(path) / 1024
        total_size += size_kb
        # Vérifier qu'on peut le relire
        try:
            df_check = pd.read_csv(path, nrows=2)
            print(f"  ✅ {f:<45} {size_kb:>8.0f} KB  {df_check.shape[1]} colonnes")
        except Exception as e:
            print(f"  ❌ {f} — Erreur lecture : {e}")

    print(f"\n  Total : {len(fichiers)} fichiers  |  {total_size/1024:.1f} MB")
    print()
    print("  PROCHAINES ÉTAPES RECOMMANDÉES :")
    print("  ┌─────────────────────────────────────────────────────┐")
    print("  │ 1. Baseline ML sur clinique TCGA                    │")
    print("  │    → LogisticRegression, RandomForest, XGBoost      │")
    print("  │    → Métriques : AUC-ROC, F1, accuracy              │")
    print("  │                                                     │")
    print("  │ 2. Modèle sur RNA GSE91061                          │")
    print("  │    → Même pipeline sur RNA_PCA50_scaled.csv         │")
    print("  │    → Comparer avec baseline clinique                │")
    print("  │                                                     │")
    print("  │ 3. Deep Learning (CNN / GNN)                        │")
    print("  │    → Input : RNA_log2_normalise.csv (21894 gènes)   │")
    print("  │    → Ou delta On-Pre pour capturer la dynamique     │")
    print("  │                                                     │")
    print("  │ 4. Si disponible : ajouter labels responder/non-    │")
    print("  │    responder depuis métadonnées GSE91061            │")
    print("  └─────────────────────────────────────────────────────┘")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("PIPELINE ML — MÉLANOME — DÉMARRAGE")
    print(f"Lecture depuis : '{os.path.abspath(DATA_DIR)}'")
    print(f"Sortie vers    : '{os.path.abspath(OUTPUT_DIR)}'")

    # Étape 1 — Normalisation RNA
    df_rna_norm, meta, genes_log = etape1_charger_normaliser_rna()

    # Étape 2 — Réduction dimensionnelle
    genes_pca, genes_top, pca_model = etape2_reduction_dimensionnelle(genes_log, meta)

    # Étape 3 — Delta Pre→On
    delta = etape3_delta_pre_on(genes_log, meta)

    # Étape 4 — Dataset ML clinique
    X_clin, y_clin = etape4_dataset_ml_clinique()

    # Étape 5 — Dataset ML RNA
    genes_pca_ml, y_preOn = etape5_dataset_ml_rna(genes_pca, meta, delta)

    # Étape 6 — Combinaison inter-cohortes
    etape6_dataset_combine(X_clin, y_clin, genes_pca_ml, meta)

    # Étape 7 — Vérification finale
    etape7_verification_finale()

    print("\n✅ PIPELINE TERMINÉ — Tous les fichiers ML sont prêts dans 'ml_ready/'")