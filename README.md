# Projet 6
## Description 
Application Streamlit permettant:
- d'explorer un dataset scRNA-seq `.h5ad`,
- de charger une signature de gènes `.txt`,
- de calculer un score de signature par cellule,
- de visualiser les résultats (UMAP, distribution, violin plot).

## Instalation 
### 1. Cloner le projet 
```bash
git clone <url_du_projet>
cd projet-6
```
### 2. Créer un environnement virtuel
```bash
python -m venv venv
source venv/bin/activate
```
### 3. Installer des dépendances
```bash
pip install -r requirements.txt
```
### 4. Lancer l'application
```bash
streamlit run app.py
```
## Données 
Placer:
- le fichier `.h5ad`
- le fichier `.txt` de signature 
dans le dossier data/.
L'application va automatiquement charger ces fichiers. 
Attention, les fichiers de données ne sont pas versionnés dans le dépôt pour éviter de surcharger Git. 



## Visualisation et analyse en cellule unique de tumoroïdes de neuroblastome

---

## Equipe

- **Etudiants** :
  - ALIOUAT Yasmine - yasmine.aliouat@etu.univ-lyon1.fr
  - CAUCHOIS Fanny - fanny.cauchois@etu.univ-lyon1.fr
  - MERAL Duzguncan - duzguncan.meral@etu.univ-lyon1.fr
  - NAVARRO Janice - janice.navarro@etu.univ-lyon1.fr
- **Maîtres d'Ouvrage** :
  - GANDRILLON Olivier (CNRS) - olivier.gandrillon@ens-lyon.fr
- **Responsable Pédagogique** :
  - MARY Arnaud - arnaud.mary@univ-lyon1.fr
