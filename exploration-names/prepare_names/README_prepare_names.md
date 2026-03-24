# Prepare_Names : 
## Description :
Ce script permet de préparer un objet AnnData pour faciliter la recherche de gènes en ajoutant des colonnes utiles comme :  
- `base_name` (nom sans version)  
- `has_versions` (présence d’un suffixe .X)  
- annotations BioMart (optionnel)  


## Fonctionnement :
Le script repose sur deux fichiers:  
- `run_prepare_names.py`: Le script principal.
- `normalise_names.py` : Contenant les fonctions de préparation des noms de gènes.  

Cette préparation se fait en étapes:
1. Lecture du fichier .h5ad.  
2. Préparation des noms de gènes :
   - création de `base_name` (nom sans version)
   - création de `has_versions`  (booléen indiquant si présence de versions)
3. Enrichissement avec les annotations BioMart (optionnel).    
4. Affichage des colonnes ajoutées.    
5. Sauvegarde d’un nouveau fichier (optionnel).


## Prérequis :
Les bibliothéques `scanpy` et `pandas` sont nécessaires pour le fonctionnement de ce script .

Si nécessaire, installer avec:

```bash
pip install scanpy pandas
```

## Utilisation :

```bash
python run_prepare_names.py chemin/vers/fichier/AnnData [--biomart] [--organism  <nom>] [--save] [--output chemin/sortie.h5ad]
```

**Options (facultatives) :**  
- `--biomart` : Pour ajouter les annotations BioMart.
- `--organism  <nom>` : Pour spécifier l'organisme utilisé pour BioMart (par défaut : hsapiens).
- `--save` : Pour sauvegarder le nouveau fichier .h5ad.
- `--output` : Pour spécifier le chemin et le nom du fichier de sortie. (uniquement si `--save` est utilisée)

**Exemple :**
```bash
python run_prepare_names.py data.h5ad --biomart --organism hsapiens --save --output data_enrechit.h5ad
```

### Comportement:

**Sans `--save` :**  
Les modifications sont faites uniquement en mémoire et aucun fichier n’est créé.  

**Avec `--save` :**  
Un nouveau fichier `.h5ad` est créé, son nom par défaut est : `nom_fichier_de_départ.prep_versions.h5ad`  
On peut aussi personnalisé le chemin et le nom du nouveau fichier avec `--output`.

### Colonnes ajoutées:
**Sans `--biomart` :**  
Les seules colonnes ajoutées sont:
- `base_name` : nom du gène sans suffixe (ex: AL627309)
- `has_versions` : booléen indiquant la présence de version (.1, .2, etc.)

**Avec `--biomart` :**   
En plus des colonnes `base_name` et `has_versions`, l'application de biomart ajoute :  

- `gene_symbol`
- `alias_symbol`
- `hgnc_symbol`
- `hgnc_id`
- `NCBI_symbol`
- `refseq_mrna`
- `uniprot_swissprot`
- `uniprot_sptrembl`

> [!NOTE]  
> Le fichier original n’est jamais modifié sauf si `--save` est utilisé.  
> BioMart a besoin d'une connection internet pour fonctionner.  
> La colonne `gene_ids` est requise pour utiliser BioMart.  Si elle est absente, elle sera automatiquement remplacée par `adata.var_names`.