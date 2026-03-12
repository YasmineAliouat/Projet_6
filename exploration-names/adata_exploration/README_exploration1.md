# Description :
Ce script permet d'afficher les noms des gènes dans un fichier AnnData(.h5ad).


## Affichage :
Ce script affiche :
- Les n premiers noms des gènes (`var_names`)
- Les n premiers identifiants des gènes (`gene_ids`)
- Si les noms des gènes sont unique
- Le nombre de noms de gènes contenant plusieurs versions (.x)

## Prérequis :
La bibliothéque scanpy est nécessaire pour le fonctionnement de ce script .

Si necessaie, installer avec:

```bash
pip install scanpy
```

## Utilisation :

```bash
python exploration1.py chemin/vers/fichier/AnnData -n <nbr>
```

- `<nbr>` : Nombre de gènes à afficher. (10 par défaut)
