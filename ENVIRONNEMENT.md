# ENVIRONNEMENT — EthnoHack 2026, construction D11 (épinglé le 02/11/2026)

- Python 3.12 (venv projet : `~/workspace/projets-revenu/ethnohack/.venv-ethnohack/`)
- **RDKit 2026.09.1** — même version que le projet Enveda (venv de référence) et que les fondations D11 éprouvées devant jury (Certifications D11, fronts C/E). Aucun écart de version à documenter pour ce chantier.
- scikit-learn 1.9.1 · pandas 3.1.0 · pyarrow 26.0.1 · numpy (dépendance RDKit)
- Vérification exécutée le 02/11 : imports RDKit/sklearn/pandas/pyarrow OK, versions lues depuis l'environnement lui-même.
- Reconstruction après purge générale du 03/11 : venv rebâti aux mêmes épingles (pyarrow résolu en 26.0.1, patch sans effet sur nos usages), imports revérifiés le 03/11 (EXIT=0).
- Données : COCONUT Lite 10-2026 (cache VM vérifié par D11 le 02/11 : 199 096 787 o, MD5 `c86adbd2beaa0c8d1841ce6c544f1306`) ; ChEMBL Pf `data/chembl_pf_activities.csv` (MD5 `ca6baae652fabb620887c9bd145321cc`).
- Note machine (02/11 ~08h05) : charge VM élevée (load ~14, RAM disponible ~1 Go) — la passe descripteurs complète (738 794 molécules) attend la porte du Gardien (charge < 12 et RAM dispo ≥ 1,5 Gio) ; validation sur échantillon 5 000 lignes d'abord.
