# HBB MD Interpretation

## Mini-projet de Bioinformatique — Master 2

### Sujet
Amélioration de l'interprétation de variants génétiques par l'intégration de la structure 3D de protéine et de données de dynamique moléculaire.

### Cas d'étude
Le projet étudie le variant pathogène du gène **HBB** responsable de l'hémoglobine S (HbS), en comparant une structure de référence de l'hémoglobine normale à une structure portant la substitution Glu→Val.

### Objectif général
Intégrer :
- l'annotation génétique et clinique du variant ;
- l'analyse structurale tridimensionnelle ;
- la dynamique moléculaire ;
- des métriques quantitatives issues des simulations ;

afin d'améliorer l'interprétation fonctionnelle du variant étudié.

### Méthodologie générale
1. Identification et annotation du variant.
2. Récupération des structures 3D WT et mutante.
3. Préparation des systèmes avec GROMACS.
4. Minimisation énergétique.
5. Équilibration NVT et NPT.
6. Dynamique moléculaire.
7. Analyse RMSD, RMSF et rayon de giration.
8. Comparaison WT vs mutant.
9. Interprétation biologique.
10. Construction d'un pipeline reproductible avec Nextflow.

### Outils
- Ubuntu 24.04 / WSL2
- Git et GitHub
- Python 3
- GROMACS
- Nextflow
- Java
- RCSB Protein Data Bank
- NCBI ClinVar

### Auteurs
Groupe de Master 2 — Informatique, Data Science et Intelligence Artificielle  
Université de Kinshasa

### Statut
Projet en cours.
