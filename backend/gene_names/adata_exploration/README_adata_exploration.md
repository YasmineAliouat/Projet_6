# Adata_Exploration :
## Description :
Ce script permet d’explorer les noms de gènes présents dans un fichier AnnData (.h5ad) afin de comprendre leur structure et les annotations disponibles.

Il est utiler pour :

- Vérifier le format des noms de gènes.
- Détecter la présence de versions (.1, .2, etc.)
- Identifier les colonnes d’annotation disponibles.
- Evaluer la qualité des métadonnées (BioMart)


## Affichage :
Ce script affiche :
- Les n premiers noms des gènes (`var_names`)
- Les n premiers identifiants des gènes (`gene_ids`)
- Si les noms des gènes sont unique
- Le nombre de gènes avec suffixe de versions (.x).
- Les colonnes d’annotation disponibles. (`gene_symbol`, `hgnc_symbol`, etc.)
- Le nombre de gènes annotés pour chaque type d’identifiant.

## Prérequis :
La bibliothèque scanpy est nécessaire pour le fonctionnement de ce script .

Si nécessaire, installer avec:

```bash
pip install scanpy
```

## Utilisation :

```bash
python adata_exploration.py chemin/vers/fichier/AnnData [-n <nbr>]
```

- `-n <nbr>` : Pour spécifier le nombre de gènes à afficher. (10 par défaut)

**Exemple :**
```bash
python adata_exploration.py data/adata_3583_new.h5ad -n 2
```

**Résultat :**
```bash
Premies 2 var_names:
Index(['AL627309.1', 'AL627309.3'], dtype='object')
 Premies 2 gene_ids:
AL627309.1    ENSG00000238009
AL627309.3    ENSG00000239945
Name: gene_ids, dtype: object
var_names uniques?: True
Nombre de noms avec plusieurs versions (.x) 8554
                   gene_ids gene_symbol alias_symbol hgnc_symbol hgnc_id NCBI_symbol refseq_mrna uniprot_swissprot uniprot_sptrembl  base_name  has_versions
AL627309.1  ENSG00000238009         NaN          NaN         NaN     NaN         NaN         NaN               NaN              NaN  AL6273091          True
AL627309.3  ENSG00000239945                                                                                    NaN              NaN  AL6273093          True
Nombre de gènes avec gene_symbol: 29499 / 30034
Nombre de gènes avec alias_symbol: 29499 / 30034
Nombre de gènes avec hgnc_symbol: 29499 / 30034
Nombre de gènes avec hgnc_id: 29499 / 30034
Nombre de gènes avec NCBI_symbol: 29499 / 30034
Nombre de gènes avec refseq_mrna: 29499 / 30034
Nombre de gènes avec uniprot_swissprot: 15853 / 30034

```
