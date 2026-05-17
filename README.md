<p align="center">
  <img src="assets/logoMelanoma.png" alt="Logo Melanoma" width="300"/>
</p>

# Melanoma Immunotherapy Data Pipeline

## 1. Résumé du projet

Ce projet porte sur la préparation de données transcriptomiques de patients atteints de mélanome afin de permettre, dans une étape ultérieure, le développement de modèles de machine learning ou de deep learning pour étudier la réponse aux traitements.

Le mélanome est un cancer de la peau qui peut devenir agressif lorsqu’il atteint un stade avancé ou métastatique. Plusieurs familles de traitements existent :

- **immunothérapie** : traitement qui aide le système immunitaire à reconnaître et attaquer les cellules tumorales ;
- **anti-PD-1 / anti-PD-L1 / anti-CTLA4** : types d’immunothérapies appelées immune checkpoint inhibitors ;
- **thérapies ciblées BRAFi / MEKi** : traitements qui ciblent des voies moléculaires spécifiques, notamment la voie BRAF/MEK ;
- **résistance thérapeutique** : situation où la tumeur ne répond pas, ou cesse de répondre, au traitement.

L’objectif de cette partie du projet n’est pas encore de construire le modèle final.  
L’objectif est de préparer des données propres, reproductibles, documentées et exploitables par les membres de l’équipe qui s’occuperont ensuite :

- du machine learning ;
- de l’analyse des résultats ;
- de l’interprétation biologique ;
- de la présentation de l’avancement au mentor.

---

## 2. Ce qui a été fait dans cette partie

Cette partie du projet correspond à la préparation des données.

Les tâches réalisées sont :

1. récupération de plusieurs cohortes publiques de mélanome ;
2. nettoyage des matrices d’expression RNA ;
3. harmonisation des noms de gènes ;
4. extraction des métadonnées cliniques ;
5. standardisation des variables importantes ;
6. fusion de deux cohortes principales ;
7. création d’un dataset principal prêt pour le machine learning ;
8. préparation de datasets externes pour validation ou analyse complémentaire ;
9. séparation claire entre :
   - dataset principal immunothérapie ;
   - cohorte externe de validation ;
   - cohorte de résistance ;
   - cohortes de thérapies ciblées ;
   - données single-cell à part.

La modélisation finale n’est pas réalisée dans cette partie.  
Le script `05_train_baseline_models.py` existe, mais la partie ML/deep learning doit être gérée par un autre membre de l’équipe.

---

## 3. Notions biologiques de base

### 3.1 Qu’est-ce qu’un gène ?

Un gène est une portion d’ADN qui contient une instruction biologique.  
Dans les données de ce projet, chaque gène est représenté par une colonne.

Exemples de gènes importants en cancer/immunologie :

```text
BRAF
TP53
CD8A
CD8B
GZMA
PRF1
PDCD1
CTLA4
LAG3
IFNG


### 3.2 Qu’est-ce que le RNA-seq ?

Le RNA-seq est une technique qui mesure l’activité des gènes dans un échantillon biologique.

Dans notre cas, les échantillons sont généralement des biopsies tumorales de patients atteints de mélanome.

Une matrice RNA-seq ressemble à ceci :

| sample_id | BRAF | TP53 | CD8A | GZMA | response |
|---|---:|---:|---:|---:|---|
| Patient_1 | 5.2 | 10.1 | 8.4 | 4.3 | responder |
| Patient_2 | 8.9 | 2.1 | 1.5 | 0.7 | non-responder |

Chaque ligne représente un échantillon ou un patient.  
Chaque colonne de gène représente l’expression d’un gène.  
La valeur numérique indique à quel point le gène est exprimé.

Dans ce projet, les matrices RNA ont souvent été transformées pour obtenir le format suivant :

```text
échantillons en lignes
gènes en colonnes
métadonnées cliniques ajoutées à gauche
```

C’est le format le plus pratique pour les modèles de machine learning.

---

### 3.3 Qu’est-ce qu’une cohorte ?

Une cohorte est un groupe de patients provenant d’une même étude scientifique.

Dans ce projet, les cohortes viennent principalement de GEO, une base de données publique qui stocke des données biologiques issues d’études scientifiques.

Les cohortes utilisées ont des noms comme :

```text
GSE91061
GSE78220
GSE168204
GSE244982
GSE145996
GSE196434
GSE77940
GSE120575
```

Chaque cohorte peut avoir :

- un format de fichier différent ;
- des noms de patients différents ;
- des noms de gènes différents ;
- des labels cliniques différents ;
- des traitements différents ;
- des définitions différentes de la réponse ;
- des plateformes techniques différentes.

Le travail principal a donc été de rendre ces données comparables.

---

### 3.4 Qu’est-ce qu’un répondeur ?

Un patient est dit “répondeur” si son état clinique s’améliore ou reste contrôlé après traitement.

Les études utilisent plusieurs notations :

| Terme | Signification |
|---|---|
| `R` | Responder |
| `NR` | Non-responder |
| `CR` | Complete response |
| `PR` | Partial response |
| `SD` | Stable disease |
| `PD` | Progressive disease |
| `LR` | Long responder |
| `SR` | Short responder |
| `resistant` | Résistant au traitement |

Dans les fichiers nettoyés, ces réponses ont souvent été transformées en variable binaire :

```text
1 = réponse favorable / répondeur / long responder
0 = réponse défavorable / non-répondeur / short responder / résistant
```

Cette variable est appelée :

```text
response_binary
```

Dans le dataset principal, la cible s’appelle :

```text
survival_event (1 survie 0 décès)
```

---

## 4. Organisation du projet

Arborescence générale :

```text
Python/
├── data_raw/
│   └── fichiers bruts téléchargés depuis GEO ou autres sources
│
├── data_processed/
│   └── fichiers générés par les scripts : données propres, fusionnées ou alignées
│
├── scripts/
│   └── scripts Python du pipeline
│
├── archive_data/
│   └── anciens fichiers ou fichiers conservés mais non utilisés directement
│
├── archive_old/
│   └── anciennes versions de scripts ou fichiers obsolètes
│
├── venv/
│   └── environnement virtuel Python local
│
├── .gitignore
├── README.md
└── scripts/requirements.txt
```

Les dossiers suivants ne sont normalement pas envoyés sur GitHub :

```text
data_raw/
data_processed/
archive_data/
archive_old/
venv/
```

Ils sont ignorés par `.gitignore`, car ils contiennent soit des données lourdes, soit des fichiers générés automatiquement, soit l’environnement virtuel local.

Le dépôt GitHub contient surtout :

```text
scripts/
README.md
.gitignore
scripts/requirements.txt
```

---

## 5. Installation et environnement Python

### 5.1 Créer un environnement virtuel

Depuis la racine du projet :

```powershell
python -m venv venv
```

Cela crée un dossier :

```text
venv/
```

---

### 5.2 Activer l’environnement virtuel sous Windows PowerShell

Si l’exécution des scripts PowerShell est bloquée, lancer :

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned
```

Puis activer l’environnement :

```powershell
.\venv\Scripts\Activate.ps1
```

Quand l’environnement est activé, PowerShell affiche :

```text
(venv)
```

---

### 5.3 Installer les dépendances

Depuis la racine du projet :

```powershell
pip install -r scripts/requirements.txt
```

Les bibliothèques principales utilisées sont :

```text
pandas
numpy
openpyxl
scikit-learn
```

Selon les scripts, d’autres bibliothèques peuvent être ajoutées plus tard.

---

## 6. Données brutes utilisées

Les fichiers bruts sont placés dans :

```text
data_raw/
```

Ces fichiers ne sont pas envoyés sur GitHub, car ils peuvent être volumineux.

Le dossier `data_raw/` contient notamment :

```text
1_CLINIQUE_donnees_principales.tsv

2_BULK_RNA_expression_GSE78220_brut.xlsx
2_BULK_RNA_metadata_GSE78220_brut.txt
2_BULK_RNA_metadata_GSE91061_brut.txt
2_BULK_RNA_score_immunitaire_cytolytique.txt

4_GRAPH_STRING_mapping_brut.txt.gz
4_GRAPH_STRING_interactions_brut.txt.gz

5_TARGETED_THERAPY_GSE77940_pre_post_raw.xlsx
5_TARGETED_THERAPY_GSE196434_raw.txt.gz
5_TARGETED_THERAPY_GSE196434_series_matrix_raw.txt.gz

6_IMMUNOTHERAPY_GSE145996_FPKM_raw.txt.gz
6_IMMUNOTHERAPY_GSE145996_series_matrix_raw.txt.gz
6_IMMUNOTHERAPY_GSE168204_counts_raw.csv.gz
6_IMMUNOTHERAPY_GSE168204_GPL11154_series_matrix_raw.txt.gz
6_IMMUNOTHERAPY_GSE168204_GPL18573_series_matrix_raw.txt.gz
6_IMMUNOTHERAPY_GSE244982_bulkRNAseq_raw.txt.gz
6_IMMUNOTHERAPY_GSE244982_series_matrix_raw.txt.gz

8_SINGLE_CELL_GSE120575_expression_TPM_raw.txt.gz
8_SINGLE_CELL_GSE120575_patient_mapping_raw.txt.gz
8_SINGLE_CELL_GSE120575_series_matrix_raw.txt.gz
```

---

## 7. Sources de données

### 7.1 GSE91061

Type de données :

```text
bulk RNA-seq
mélanome
immunothérapie
réponse clinique
```

Rôle :

```text
cohorte principale du dataset ML
```

Utilisation :

- expression RNA ;
- labels de réponse ;
- conversion des identifiants de gènes ;
- fusion avec GSE78220.

Fichiers utilisés :

```text
data_raw/2_BULK_RNA_metadata_GSE91061_brut.txt
data_processed/2_BULK_RNA_expression_GSE91061_symboles.csv
```

Remarque :

Le fichier d’expression propre de GSE91061 existait déjà dans le projet.  
Le script `03_convert_gse91061_entrez_to_symbol.py` convertit les identifiants numériques en symboles de gènes.

---

### 7.2 GSE78220

Type de données :

```text
bulk RNA-seq
mélanome
anti-PD-1
réponse clinique
```

Rôle :

```text
cohorte principale du dataset ML
```

Fichiers bruts :

```text
data_raw/2_BULK_RNA_expression_GSE78220_brut.xlsx
data_raw/2_BULK_RNA_metadata_GSE78220_brut.txt
```

Fichier généré :

```text
data_processed/2_BULK_RNA_expression_GSE78220_nettoyee.csv
```

Informations cliniques extraites :

```text
patient_id
response
age
gender
disease_status
overall_survival_days
vital_status
previous_mapki
biopsy_time
treatment
```

---

### 7.3 GSE168204

Type de données :

```text
bulk RNA-seq
mélanome métastatique
immunothérapie
anti-PD1 / anti-PDL1 / anti-PD1+anti-CTLA4
réponse R / NR
timepoint PRE / ON
```

Rôle :

```text
validation externe potentielle pour le modèle immunothérapie
```

Fichiers bruts :

```text
data_raw/6_IMMUNOTHERAPY_GSE168204_counts_raw.csv.gz
data_raw/6_IMMUNOTHERAPY_GSE168204_GPL11154_series_matrix_raw.txt.gz
data_raw/6_IMMUNOTHERAPY_GSE168204_GPL18573_series_matrix_raw.txt.gz
```

Fichiers générés :

```text
data_processed/10_METADATA_GSE168204_GPL11154.csv
data_processed/10_METADATA_GSE168204_GPL18573.csv
data_processed/11_METADATA_GSE168204_clean.csv
data_processed/12_GSE168204_bulk_clean.csv
data_processed/13_GSE168204_EXTERNAL_VALIDATION_READY.csv
```

Résultat final :

```text
27 échantillons
17 322 gènes communs avec le dataset principal
18 non-répondeurs
9 répondeurs
```

Utilisation recommandée :

```text
validation externe indépendante du modèle principal
```

---

### 7.4 GSE244982

Type de données :

```text
bulk RNA-seq
mélanome
résistance à l’immunothérapie
anti-PD1 resistant
anti-CTLA4 resistant
```

Rôle :

```text
analyse complémentaire de résistance
```

Fichiers bruts :

```text
data_raw/6_IMMUNOTHERAPY_GSE244982_bulkRNAseq_raw.txt.gz
data_raw/6_IMMUNOTHERAPY_GSE244982_series_matrix_raw.txt.gz
```

Fichiers générés :

```text
data_processed/10_METADATA_GSE244982.csv
data_processed/11_METADATA_GSE244982_clean.csv
data_processed/12_GSE244982_bulk_clean.csv
data_processed/13_GSE244982_RESISTANCE_READY.csv
```

Résultat final :

```text
41 échantillons
17 216 gènes communs avec le dataset principal
21 anti-CTLA4 resistant
20 anti-PD1 resistant
```

Utilisation recommandée :

```text
analyse des profils transcriptomiques de résistance
```

Important :

Ce dataset ne doit pas être utilisé comme simple cohorte responder/non-responder.  
Tous les échantillons sont des cas résistants.

---

### 7.5 GSE145996

Type de données :

```text
bulk RNA-seq
mélanome
immunothérapie potentielle
```

Fichiers bruts :

```text
data_raw/6_IMMUNOTHERAPY_GSE145996_FPKM_raw.txt.gz
data_raw/6_IMMUNOTHERAPY_GSE145996_series_matrix_raw.txt.gz
```

Résultat de l’inspection :

```text
expression RNA disponible
métadonnées partielles
pas de label clair de réponse dans les fichiers récupérés
```

Informations récupérées :

```text
sample type
age at diagnosis
gender
stage at diagnosis
```

Statut :

```text
non intégré au dataset final
gardé en réserve
```

Raison :

Le fichier metadata récupéré ne contient pas suffisamment d’information pour créer une cible claire de type responder/non-responder.

---

### 7.6 GSE77940

Type de données :

```text
bulk RNA-seq
mélanome
thérapie ciblée RAF / MEK
échantillons Pre et Post
```

Rôle :

```text
analyse complémentaire pré/post thérapie ciblée
```

Fichier brut :

```text
data_raw/5_TARGETED_THERAPY_GSE77940_pre_post_raw.xlsx
```

Fichiers générés :

```text
data_processed/5_TARGETED_THERAPY_GSE77940_clean.csv
data_processed/5_TARGETED_THERAPY_GSE77940_summary.csv
data_processed/5_TARGETED_THERAPY_GSE77940_ALIGNED_WITH_MAIN_GENES.csv
data_processed/5_TARGETED_THERAPY_GSE77940_alignment_summary.txt
```

Résultat :

```text
12 échantillons
6 patients
Pre / Post
RAFi ou RAF_MEKi
16 564 gènes communs avec le dataset principal
```

Utilisation recommandée :

```text
analyse exploratoire pré/post traitement ciblé
```

Limite :

```text
seulement 6 patients
pas adapté à un modèle ML robuste
```

---

### 7.7 GSE196434

Type de données :

```text
bulk RNA-seq
mélanome
thérapie ciblée BRAFi / BRAFi+MEKi
long responders / short responders
PRE / POST
```

Rôle :

```text
analyse complémentaire sur les thérapies ciblées
```

Fichiers bruts :

```text
data_raw/5_TARGETED_THERAPY_GSE196434_raw.txt.gz
data_raw/5_TARGETED_THERAPY_GSE196434_series_matrix_raw.txt.gz
```

Fichiers générés :

```text
data_processed/10_METADATA_GSE196434.csv
data_processed/11_METADATA_GSE196434_clean.csv
data_processed/12_GSE196434_bulk_clean.csv
data_processed/13_GSE196434_TARGETED_THERAPY_READY.csv
```

Résultat :

```text
91 échantillons
48 long responders
43 short responders
69 PRE
22 POST
55 BRAFi
36 BRAFi+MEKi
```

Important :

Le fichier brut utilise principalement des identifiants Ensembl.  
Une conversion vers symboles de gènes a été tentée avec STRING, mais elle n’est pas parfaite.

Après alignement avec le dataset principal :

```text
seulement 622 gènes communs directs
```

Donc, pour une analyse large de GSE196434, il vaut mieux utiliser :

```text
data_processed/12_GSE196434_bulk_clean.csv
```

plutôt que :

```text
data_processed/13_GSE196434_TARGETED_THERAPY_READY.csv
```

---

### 7.8 GSE120575

Type de données :

```text
single-cell RNA-seq
mélanome
immunothérapie
responder / non-responder
```

Rôle :

```text
piste complémentaire single-cell
non fusionnée avec le bulk RNA-seq principal
```

Fichiers bruts :

```text
data_raw/8_SINGLE_CELL_GSE120575_expression_TPM_raw.txt.gz
data_raw/8_SINGLE_CELL_GSE120575_patient_mapping_raw.txt.gz
data_raw/8_SINGLE_CELL_GSE120575_series_matrix_raw.txt.gz
```

Fichier généré :

```text
data_processed/11_METADATA_GSE120575_clean.csv
```

Résultat metadata :

```text
48 échantillons
17 responders
31 non-responders
anti-CTLA4
anti-PD1
anti-CTLA4+PD1
```

Important :

Le single-cell RNA-seq est différent du bulk RNA-seq.  
Il ne faut donc pas mélanger directement GSE120575 avec les datasets bulk sans traitement spécifique.

---

### 7.9 STRING

Type de données :

```text
mapping de gènes
interactions gène/protéine
```

Fichiers :

```text
data_raw/4_GRAPH_STRING_mapping_brut.txt.gz
data_raw/4_GRAPH_STRING_interactions_brut.txt.gz
```

Utilisation actuelle :

```text
conversion Entrez ID vers gene symbol pour GSE91061
tentative de conversion Ensembl vers gene symbol pour GSE196434
```

Utilisation future possible :

```text
construction d’un graphe de gènes
GNN
analyse réseau
```

---

## 8. Variables importantes

### 8.1 Variables d’identification

| Variable | Description |
|---|---|
| `study` | Cohorte d’origine |
| `sample_id` | Identifiant de l’échantillon |
| `patient_id` | Identifiant patient |
| `geo_accession` | Identifiant GEO |
| `platform` | Plateforme GEO |

---

### 8.2 Variables cliniques

| Variable | Description |
|---|---|
| `response` | Réponse clinique brute |
| `response_binary` | Réponse harmonisée : 1 = favorable, 0 = défavorable |
| `survival_event` | Cible binaire du dataset principal |
| `treatment` | Type de traitement |
| `timepoint` | Moment du prélèvement |
| `dataset_role` | Rôle du dataset |
| `notes` | Informations complémentaires |

---

### 8.3 Variables transcriptomiques

Toutes les autres colonnes sont généralement des gènes.

Exemples :

```text
TP53
BRAF
CD8A
CD8B
GZMA
PRF1
PDCD1
CTLA4
IFNG
```

---

## 9. Fichiers finaux obtenus

### 9.1 Dataset principal ML-ready

```text
data_processed/3_DATASET_RNA_MULTICOHORTE_ML_READY.csv
```

Contenu :

```text
GSE91061 + GSE78220
133 échantillons
17 355 gènes
patient_id
study
response
survival_event
```

Distribution :

```text
survival_event = 1 : 72
survival_event = 0 : 61
```

Utilisation :

```text
dataset principal pour le machine learning
```

---

### 9.2 Dataset principal complet

```text
data_processed/3_DATASET_RNA_MULTICOHORTE_FINAL.csv
```

Contenu :

```text
version complète du dataset principal
colonnes cliniques supplémentaires
cytolytic_score quand disponible
```

Utilisation :

```text
analyse descriptive
vérification
archive propre
```

---

### 9.3 Validation externe

```text
data_processed/13_GSE168204_EXTERNAL_VALIDATION_READY.csv
```

Contenu :

```text
27 échantillons
17 322 gènes communs avec le dataset principal
response_binary
treatment
timepoint
```

Utilisation :

```text
validation externe potentielle du modèle
```

---

### 9.4 Résistance immunothérapie

```text
data_processed/13_GSE244982_RESISTANCE_READY.csv
```

Contenu :

```text
41 échantillons
17 216 gènes communs avec le dataset principal
anti-PD1 resistant
anti-CTLA4 resistant
```

Utilisation :

```text
analyse complémentaire de résistance
```

---

### 9.5 Thérapie ciblée GSE196434

```text
data_processed/12_GSE196434_bulk_clean.csv
```

Contenu :

```text
91 échantillons
59 729 gènes
LR / SR
PRE / POST
BRAFi / BRAFi+MEKi
```

Utilisation :

```text
analyse complémentaire des thérapies ciblées
```

---

### 9.6 Thérapie ciblée GSE77940

```text
data_processed/5_TARGETED_THERAPY_GSE77940_ALIGNED_WITH_MAIN_GENES.csv
```

Contenu :

```text
12 échantillons
6 patients
16 564 gènes communs
Pre / Post
RAFi / RAF_MEKi
```

Utilisation :

```text
analyse exploratoire pré/post
```

---

### 9.7 Métadonnées standardisées

```text
data_processed/11_METADATA_EXTRA_ALL_CLEAN.csv
```

Contenu :

```text
207 échantillons supplémentaires
GSE168204
GSE196434
GSE244982
GSE120575
```

Distribution :

```text
response_binary = 0 : 133
response_binary = 1 : 74
```

---

## 10. Scripts et ordre d’exécution

Tous les scripts doivent être lancés depuis la racine du projet.

Exemple :

```powershell
python scripts/01_prepare_clinical_multicohort.py
```

---

### 10.1 Pipeline principal

#### 1. `01_prepare_clinical_multicohort.py`

Objectif :

```text
préparer les données cliniques principales pour GSE91061 et GSE78220
```

Entrées :

```text
data_raw/2_BULK_RNA_metadata_GSE78220_brut.txt
archive_data/1_CLINIQUE_GSE91061_nettoyee.csv
```

Sortie :

```text
data_processed/1_CLINIQUE_MULTICOHORTE_nettoyee.csv
```

Ce script extrait :

```text
patient_id
response
survival_event
response_binary
study
age
gender
disease_status
overall_survival_days
vital_status
previous_mapki
biopsy_time
treatment
```

---

#### 2. `02_clean_gse78220_rna.py`

Objectif :

```text
nettoyer la matrice RNA de GSE78220
```

Entrée :

```text
data_raw/2_BULK_RNA_expression_GSE78220_brut.xlsx
```

Sortie :

```text
data_processed/2_BULK_RNA_expression_GSE78220_nettoyee.csv
```

Ce script :

```text
lit le fichier Excel
identifie les gènes
transpose la matrice
met les patients en lignes
met les gènes en colonnes
nettoie les identifiants patients
```

---

#### 3. `03_convert_gse91061_entrez_to_symbol.py`

Objectif :

```text
convertir les identifiants Entrez de GSE91061 en symboles de gènes
```

Entrées :

```text
data_processed/2_BULK_RNA_expression_GSE91061_nettoyee.csv
data_raw/4_GRAPH_STRING_mapping_brut.txt.gz
```

Sortie :

```text
data_processed/2_BULK_RNA_expression_GSE91061_symboles.csv
```

Ce script :

```text
lit la matrice GSE91061
utilise STRING comme mapping
convertit les identifiants numériques
gère les doublons de gènes
```

---

#### 4. `04_build_multicohort_rna_dataset.py`

Objectif :

```text
fusionner GSE91061 et GSE78220
```

Entrées :

```text
data_processed/2_BULK_RNA_expression_GSE91061_symboles.csv
data_processed/2_BULK_RNA_expression_GSE78220_nettoyee.csv
data_processed/1_CLINIQUE_MULTICOHORTE_nettoyee.csv
data_raw/2_BULK_RNA_score_immunitaire_cytolytique.txt
```

Sorties :

```text
data_processed/3_DATASET_RNA_MULTICOHORTE_FINAL.csv
data_processed/3_DATASET_RNA_MULTICOHORTE_ML_READY.csv
```

Ce script :

```text
identifie les gènes communs
fusionne les deux cohortes
ajoute les métadonnées cliniques
ajoute le cytolytic_score
retire les échantillons sans label
crée une version complète
crée une version ML-ready
```

Résultat final du pipeline principal :

```text
133 échantillons
17 355 gènes
72 classe 1
61 classe 0
```

---

### 10.2 Script de modélisation

#### `05_train_baseline_models.py`

Objectif :

```text
script prévu pour la partie machine learning
```

Statut :

```text
non traité dans cette partie
réservé au membre responsable du ML/deep learning
```

---

### 10.3 Pipeline thérapies ciblées

#### `06_inspect_targeted_therapy_datasets.py`

Objectif :

```text
inspecter GSE196434 et GSE77940 avant nettoyage
```

Entrées :

```text
data_raw/5_TARGETED_THERAPY_GSE196434_raw.txt.gz
data_raw/5_TARGETED_THERAPY_GSE77940_pre_post_raw.xlsx
```

Sorties :

```text
data_processed/6_TARGETED_THERAPY_inspection_report.txt
data_processed/6_GSE196434_preview.csv
data_processed/6_GSE77940_preview_*.csv
```

---

#### `07_clean_gse77940_targeted_therapy.py`

Objectif :

```text
nettoyer GSE77940
```

Entrée :

```text
data_raw/5_TARGETED_THERAPY_GSE77940_pre_post_raw.xlsx
```

Sorties :

```text
data_processed/5_TARGETED_THERAPY_GSE77940_clean.csv
data_processed/5_TARGETED_THERAPY_GSE77940_summary.csv
```

Résultat :

```text
12 échantillons
6 patients
48 643 gènes
Pre / Post
RAFi / RAF_MEKi
```

---

#### `08_align_targeted_therapy_with_main_genes.py`

Objectif :

```text
aligner GSE77940 sur les gènes du dataset principal
```

Entrées :

```text
data_processed/3_DATASET_RNA_MULTICOHORTE_ML_READY.csv
data_processed/5_TARGETED_THERAPY_GSE77940_clean.csv
```

Sorties :

```text
data_processed/5_TARGETED_THERAPY_GSE77940_ALIGNED_WITH_MAIN_GENES.csv
data_processed/5_TARGETED_THERAPY_GSE77940_alignment_summary.txt
```

Résultat :

```text
12 échantillons
16 564 gènes communs
```

---

### 10.4 Pipeline datasets supplémentaires

#### `09_inspect_extra_immunotherapy_datasets.py`

Objectif :

```text
inspecter GSE168204, GSE244982 et GSE145996
```

Entrées :

```text
data_raw/6_IMMUNOTHERAPY_GSE168204_counts_raw.csv.gz
data_raw/6_IMMUNOTHERAPY_GSE244982_bulkRNAseq_raw.txt.gz
data_raw/6_IMMUNOTHERAPY_GSE145996_FPKM_raw.txt.gz
```

Sorties :

```text
data_processed/9_EXTRA_IMMUNOTHERAPY_inspection_report.txt
data_processed/9_EXTRA_IMMUNOTHERAPY_inspection_summary.csv
data_processed/9_GSE168204_preview.csv
data_processed/9_GSE244982_preview.csv
data_processed/9_GSE145996_preview.csv
```

Résultat important :

```text
GSE168204 : 17 322 gènes communs directs
GSE244982 : 17 216 gènes communs directs
GSE145996 : 15 137 gènes communs directs
```

---

#### `10_inspect_metadata_files.py`

Objectif :

```text
extraire les métadonnées depuis les fichiers series_matrix GEO
```

Entrées :

```text
data_raw/6_IMMUNOTHERAPY_GSE145996_series_matrix_raw.txt.gz
data_raw/6_IMMUNOTHERAPY_GSE168204_GPL11154_series_matrix_raw.txt.gz
data_raw/6_IMMUNOTHERAPY_GSE168204_GPL18573_series_matrix_raw.txt.gz
data_raw/6_IMMUNOTHERAPY_GSE244982_series_matrix_raw.txt.gz
data_raw/5_TARGETED_THERAPY_GSE196434_series_matrix_raw.txt.gz
data_raw/8_SINGLE_CELL_GSE120575_series_matrix_raw.txt.gz
```

Sorties :

```text
data_processed/10_METADATA_inspection_report.txt
data_processed/10_METADATA_inspection_summary.csv
data_processed/10_METADATA_GSE145996.csv
data_processed/10_METADATA_GSE168204_GPL11154.csv
data_processed/10_METADATA_GSE168204_GPL18573.csv
data_processed/10_METADATA_GSE244982.csv
data_processed/10_METADATA_GSE196434.csv
data_processed/10_METADATA_GSE120575.csv
```

Ce script repère :

```text
response
treatment
timepoint
patient_id
sample_id
platform
age
gender
stage
```

---

#### `11_prepare_extra_metadata_tables.py`

Objectif :

```text
standardiser les métadonnées des datasets supplémentaires
```

Entrées :

```text
data_processed/10_METADATA_GSE168204_GPL11154.csv
data_processed/10_METADATA_GSE168204_GPL18573.csv
data_processed/10_METADATA_GSE196434.csv
data_processed/10_METADATA_GSE244982.csv
data_processed/10_METADATA_GSE120575.csv
```

Sorties :

```text
data_processed/11_METADATA_GSE168204_clean.csv
data_processed/11_METADATA_GSE196434_clean.csv
data_processed/11_METADATA_GSE244982_clean.csv
data_processed/11_METADATA_GSE120575_clean.csv
data_processed/11_METADATA_EXTRA_ALL_CLEAN.csv
data_processed/11_METADATA_EXTRA_summary.csv
```

Colonnes standardisées :

```text
study
platform
sample_id
geo_accession
patient_id
treatment
timepoint
response
response_binary
dataset_role
notes
```

Résultat global :

```text
207 échantillons supplémentaires
133 classe 0
74 classe 1
```

---

#### `12_prepare_extra_bulk_expression_datasets.py`

Objectif :

```text
associer expression RNA et métadonnées pour GSE168204, GSE244982 et GSE196434
```

Entrées :

```text
data_raw/6_IMMUNOTHERAPY_GSE168204_counts_raw.csv.gz
data_raw/6_IMMUNOTHERAPY_GSE244982_bulkRNAseq_raw.txt.gz
data_raw/5_TARGETED_THERAPY_GSE196434_raw.txt.gz

data_processed/11_METADATA_GSE168204_clean.csv
data_processed/11_METADATA_GSE244982_clean.csv
data_processed/11_METADATA_GSE196434_clean.csv
```

Sorties :

```text
data_processed/12_GSE168204_bulk_clean.csv
data_processed/12_GSE244982_bulk_clean.csv
data_processed/12_GSE196434_bulk_clean.csv
data_processed/12_EXTRA_BULK_summary.csv
```

Résultats :

```text
GSE168204 : 27 échantillons, 75 253 gènes
GSE244982 : 41 échantillons, 18 759 gènes
GSE196434 : 91 échantillons, 59 729 gènes
```

---

#### `13_align_extra_bulk_with_main_genes.py`

Objectif :

```text
aligner les datasets supplémentaires avec les gènes du dataset principal
```

Entrées :

```text
data_processed/3_DATASET_RNA_MULTICOHORTE_ML_READY.csv
data_processed/12_GSE168204_bulk_clean.csv
data_processed/12_GSE244982_bulk_clean.csv
data_processed/12_GSE196434_bulk_clean.csv
```

Sorties :

```text
data_processed/13_GSE168204_EXTERNAL_VALIDATION_READY.csv
data_processed/13_GSE244982_RESISTANCE_READY.csv
data_processed/13_GSE196434_TARGETED_THERAPY_READY.csv
data_processed/13_EXTRA_BULK_ALIGNMENT_summary.csv
```

Résultats :

```text
GSE168204 : 27 échantillons, 17 322 gènes communs
GSE244982 : 41 échantillons, 17 216 gènes communs
GSE196434 : 91 échantillons, 622 gènes communs directs
```

---

## 11. Ordre complet d’exécution

### 11.1 Recréer le dataset principal

```powershell
python scripts/01_prepare_clinical_multicohort.py
python scripts/02_clean_gse78220_rna.py
python scripts/03_convert_gse91061_entrez_to_symbol.py
python scripts/04_build_multicohort_rna_dataset.py
```

---

### 11.2 Préparer les datasets de thérapies ciblées

```powershell
python scripts/06_inspect_targeted_therapy_datasets.py
python scripts/07_clean_gse77940_targeted_therapy.py
python scripts/08_align_targeted_therapy_with_main_genes.py
```

---

### 11.3 Préparer les datasets supplémentaires

```powershell
python scripts/09_inspect_extra_immunotherapy_datasets.py
python scripts/10_inspect_metadata_files.py
python scripts/11_prepare_extra_metadata_tables.py
python scripts/12_prepare_extra_bulk_expression_datasets.py
python scripts/13_align_extra_bulk_with_main_genes.py
```

---

## 12. Que doit utiliser chaque membre de l’équipe ?

### 12.1 Pour la personne qui fait le machine learning

Dataset principal :

```text
data_processed/3_DATASET_RNA_MULTICOHORTE_ML_READY.csv
```

Validation externe possible :

```text
data_processed/13_GSE168204_EXTERNAL_VALIDATION_READY.csv
```

Ne pas utiliser directement pour le modèle principal immunothérapie :

```text
GSE120575
GSE196434
GSE77940
```

car ces datasets correspondent à d’autres types de données ou d’autres traitements.

---

### 12.2 Pour les personnes qui analysent les résultats

Fichiers utiles :

```text
data_processed/3_DATASET_RNA_MULTICOHORTE_FINAL.csv
data_processed/3_DATASET_RNA_MULTICOHORTE_ML_READY.csv
data_processed/13_GSE168204_EXTERNAL_VALIDATION_READY.csv
data_processed/13_GSE244982_RESISTANCE_READY.csv
data_processed/12_GSE196434_bulk_clean.csv
data_processed/5_TARGETED_THERAPY_GSE77940_ALIGNED_WITH_MAIN_GENES.csv
```

Points à analyser :

```text
distribution des classes
différences entre cohortes
gènes communs
traitements représentés
réponse vs non-réponse
résistance
limites des datasets
```

---

### 12.3 Pour les personnes qui présentent au mentor

Message principal à transmettre :

```text
Nous avons construit un pipeline de préparation de données transcriptomiques de mélanome.
Le dataset principal combine deux cohortes immunothérapie : GSE91061 et GSE78220.
Nous avons harmonisé les noms de gènes et les labels cliniques.
Nous avons généré un dataset ML-ready de 133 échantillons et 17 355 gènes.
Nous avons aussi préparé GSE168204 comme validation externe potentielle.
Nous avons préparé GSE244982 pour l’analyse de résistance.
Nous avons préparé GSE196434 et GSE77940 pour une analyse complémentaire des thérapies ciblées.
Les données single-cell GSE120575 sont disponibles comme piste complémentaire mais non fusionnées au bulk RNA.
```

---

## 13. Limites actuelles

### 13.1 Taille des données

Le dataset principal contient 133 échantillons.  
C’est suffisant pour des modèles classiques, mais limité pour du deep learning complexe.

---

### 13.2 Hétérogénéité des cohortes

Les datasets viennent d’études différentes.  
Les protocoles expérimentaux, plateformes et méthodes de normalisation peuvent varier.

Cela peut créer des effets de batch.

---

### 13.3 GSE196434

La conversion des identifiants Ensembl vers gene symbols est incomplète.  
Cela limite l’alignement direct avec le dataset principal.

---

### 13.4 GSE145996

Le fichier expression est disponible, mais les métadonnées récupérées ne contiennent pas encore de label clair de réponse.  
Le dataset n’est donc pas utilisé dans la version finale actuelle.

---

### 13.5 GSE120575

Ce dataset est single-cell.  
Il ne doit pas être fusionné directement avec les datasets bulk RNA-seq.

---

## 14. État final du projet à ce stade

À la fin de cette partie, on dispose de :

```text
1 dataset principal ML-ready
1 cohorte externe de validation immunothérapie
1 cohorte de résistance immunothérapie
2 datasets de thérapies ciblées
1 dataset single-cell en réserve
des métadonnées standardisées
un pipeline reproductible en scripts Python
```

La partie data preparation est donc prête pour la suite.

La suite du projet peut maintenant porter sur :

```text
analyse exploratoire
feature selection
modélisation classique
validation externe
interprétation biologique
présentation au mentor
deep learning si justifié
```

---

## 15. Commandes Git utiles

Après modification d’un script ou du README :

```powershell
git status
git add scripts/nom_du_script.py
git commit -m "Message clair du changement"
git push
```

Pour ajouter le README :

```powershell
git add README.md
git commit -m "Add complete project documentation"
git push
```

Les fichiers de données ne doivent normalement pas être ajoutés à GitHub.
