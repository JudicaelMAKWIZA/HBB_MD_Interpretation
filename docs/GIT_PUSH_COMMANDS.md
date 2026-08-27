# Publication du projet sur GitHub

Dépôt distant : https://github.com/JudicaelMAKWIZA/HBB_MD_Interpretation

Depuis WSL, même après avoir fermé l'ancien terminal :

```bash
cd ~/HBB_MD_Interpretation

git status
git remote -v
```

Vérifier que les fichiers lourds GROMACS sont ignorés :

```bash
git status --short | grep -E '\.(xtc|trr|cpt|edr|tpr|log)$' || true
```

Préparer le commit :

```bash
git add README.md .gitignore data mdp scripts workflow tests figures presentation_assets results docs

git status
```

Créer le commit :

```bash
git commit -m "Finalize HBB WT vs HbS molecular dynamics practical project"
```

Publier :

```bash
git push origin main
```

Vérification finale :

```bash
git status
git log --oneline -3
```

Le `git status` final doit indiquer un arbre de travail propre.
