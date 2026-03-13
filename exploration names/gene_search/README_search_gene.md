# Description :
Ce script permet de rechercher un gène dans un fichier AnnData(.h5ad) à partir de différents identifiants biologiques.


## Fonctionnement :
La recherche est effectuée en plusieurs étapes :

**1. Recherche exacte**

Le script cherche d'abord une correspondance exacte dans les différentes colonnes d'annotation du gène.

**2. Recherche partielle**

Si aucune correspondance exacte n'est trouvée, une recherche par correspondance partielle est effectuée.

**3. Suggestions**

Si aucun résultat n'est trouvé, le script propose des noms de gènes proches.

---

***Gestion des versions de gènes***

Certains gènes peuvent exister avec plusieurs versions dans le dataset.

Exemple :
```
AL390719.1
AL390719.3
```

Si l'utilisateur tape `AL390719` le script retourne toutes les versions disponibles et l'utilisateur peut choisir celle à utiliser.


## Prérequis :
Les bibliothéques `scanpy`, `pandas` et `difflib` sont nécessaires pour le fonctionnement de ce script .

Si necessaie, installer avec:

```bash
pip install scanpy pandas difflib
```

## Utilisation :

```bash
python search_gene.py chemin/vers/fichier/AnnData <nom_du_gene> --max <n> --choose <version> --interactive --show-partial
```

- `<nom_du_gene>` : Nom du gène recherché.
- `--max <n>` : Nombre maximum de résultats affichés (par défaut : 30)
- `--choose <version>` : Pour sélectionner directement un résultat.
- `--interactive` : Pour demander à l'utilisateur de choisir.
- `--show-partial` : Pour afficher les résultats partiels même si un résultat exact existe.

## Structure:

Le script utilise les fonctions définies dans le script `gene_search_utils.py` qui permettent:
- la recherche exacte
- la recherche partielle
- la suggestions de gènes proches
- l'affichage les résultats

**Exemple de sortie**

Recherche:

```
python3 search_gene.py /home/lili/Documents/M1/S2/Projet6/adata_3583_new.h5ad TTL7A
```
Résultat:

```
Résultat partiel :

1. METTL7A | gene_ids=ENSG00000185432 | gene_symbol=TMT1A | hgnc_symbol=TMT1A | alias_symbol=DKFZP586A0522|METTL7A | refseq_mrna=NM_014033 | uniprot_swissprot=Q9H8H3 | base_name=METTL7A

Gène sélectionné : METTL7A
```