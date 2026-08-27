# Publication du projet sur GitHub

Dépôt distant : <https://github.com/JudicaelMAKWIZA/HBB_MD_Interpretation>

Ce ZIP contient déjà le projet nettoyé, le README mis à jour, les figures, le diaporama PDF final et le dossier LaTeX/Overleaf. Il contient aussi les commits préparés localement.

## Option recommandée : pousser depuis le ZIP final

Dans un nouveau terminal WSL :

```bash
WINUSER=$(cmd.exe /c echo %USERNAME% | tr -d '\r')
mkdir -p ~/HBB_publish
unzip -o "/mnt/c/Users/$WINUSER/Downloads/HBB_MD_Interpretation_final_v2_ready_to_push.zip" -d ~/HBB_publish
cd ~/HBB_publish/HBB_MD_Interpretation
```

Vérifie l'état :

```bash
git status -sb
git remote -v
git log --oneline -5
```

Tu dois voir une branche `main` en avance sur `origin/main`, par exemple `[ahead 2]`.

Puis pousse directement :

```bash
git push origin main
```

## Option si tu travailles dans ton dossier local initial

Si tu veux plutôt utiliser ton dossier local `~/HBB_MD_Interpretation`, remplace-le par le contenu du ZIP final ou copie les fichiers modifiés, puis fais :

```bash
cd ~/HBB_MD_Interpretation

git status -sb
git add README.md .gitignore data mdp scripts workflow tests figures presentation_assets results docs
git commit -m "Revise final Beamer presentation and README"
git push origin main
```

## Vérification finale

```bash
git status
git log --oneline -5
```

Le `git status` final doit indiquer que la branche `main` est à jour avec `origin/main` et que l'arbre de travail est propre.
