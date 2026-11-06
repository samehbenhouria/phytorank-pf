# PhytoRank-Pf — Rapport de méthodes et résultats (siège D11, 2026-11-06)

Priorisation computationnelle de produits naturels contre *Plasmodium
falciparum* (EthnoHack 2026, Track 2). **Portée exacte : des hypothèses de
priorisation pour une revue experte — aucune affirmation d'activité
biologique, d'efficacité ou d'usage médical.** Aucun candidat n'a été testé
expérimentalement dans ce projet.

## 1. Données

| Jeu | Source | Contenu |
|---|---|---|
| COCONUT Lite 10-2026 | COCONUT 2.0 (téléchargement officiel, MD5 vérifié `c86adbd2…`) | 738 794 molécules, 40 colonnes ; organismes sources absents du Lite |
| Activités *P. falciparum* | ChEMBL, cible CHEMBL364, IC50/EC50 | 31 652 molécules uniques, pChEMBL médian par molécule (médiane 6,19) |
| Organismes du top-N | COCONUT Full 10-2026 (même millésime), jointure locale par le fil projet | voir §7 |

Limite des données d'entraînement, dite d'emblée : CHEMBL364 agrège des
essais sur organisme entier — souches et protocoles mélangés. Le modèle
apprend « ressemble aux molécules rapportées actives dans ces essais »,
rien de plus.

## 2. Protocole (gelé avant entraînement)

- Étiquette : actif si pChEMBL médian ≥ seuil ; seuil de référence **6**,
  sensibilité rapportée à 5 et 7.
- Validation : **split par échafaudage de Bemis–Murcko** (70/15/15,
  glouton déterministe, graine 42) avec **assertion de disjonction**
  exécutée dans le code — aucun échafaudage partagé entre plis. Aucun
  split aléatoire n'est présenté comme validation (mesures du domaine
  sur le banc BBBP, Certifications D11 front E : en passant du split
  aléatoire au split par échafaudage, le ROC-AUC chute de 0,116 pour
  ECFP4+forêt aléatoire et de 0,044 pour un modèle à descripteurs).
- Représentations comparées : ECFP4 (Morgan r=2, 2048 bits) seul ;
  descripteurs RDKit seuls (masse, logP, TPSA, HBD, HBA, liaisons
  rotatives, fraction Csp3, cycles aromatiques, QED) ; **combiné**.
- Modèle : forêt aléatoire 300 arbres, `class_weight=balanced_subsample`.
- Calibration : Platt (régression logistique) ajustée sur le pli de
  validation seulement.
- Métriques : ROC-AUC, PR-AUC **toujours publiée avec sa prévalence**,
  précision@50/@100, enrichissement à 1 %, score de Brier avant/après
  calibration, table de fiabilité.

## 3. Descripteurs COCONUT recalculés — validation

Recalcul RDKit 2026.09.1 des descripteurs et du QED des 738 794
molécules (passe en blocs fermés un à un, assemblage publié seulement
après relecture de contrôle : 738 794 lignes, pied de fichier valide) :
**738 790 valides, 4 illisibles**. QED recalculé vs QED fourni par
COCONUT : Pearson **0,99993**, écart absolu médian **0,0027**, 99,99 %
des écarts ≤ 0,01 (2 molécules > 0,1, max 0,16) — les valeurs fournies
sont cohérentes au niveau de l'arrondi ; masses : écart médian 0,002.

## 4. Résultats du modèle (pli test, split par échafaudage, n = 4 748)

Prévalence du pli test au seuil 6 : **0,607** (le hasard vaut donc 0,607
en précision et en PR-AUC).

| Représentation | ROC-AUC | PR-AUC | P@50 | P@100 | Brier |
|---|---|---|---|---|---|
| ECFP4 seul | 0,8700 | 0,9076 | 0,98 | 0,99 | 0,1470 |
| Descripteurs seuls | 0,7444 | 0,8204 | 0,98 | 0,98 | 0,1989 |
| **Combiné** | **0,8710** | **0,9093** | **0,98** | **0,98** | 0,1477 → **0,1420** après Platt |

Sensibilité (combiné) : seuil 5 → ROC 0,9084 (prévalence 0,897) ;
seuil 7 → ROC 0,8779 (prévalence 0,262).

Lecture honnête :
- Le combiné ne bat l'ECFP seul que d'un souffle (+0,001) sur ce jeu —
  contrairement à d'autres bancs du domaine où les descripteurs aident
  nettement plus. Nous retenons le combiné (protocole) sans en exagérer
  l'apport.
- L'enrichissement à 1 % est **saturé** (1,61 pour les trois
  représentations, plafond théorique 1/0,607 = 1,65) : non discriminant
  à cette prévalence, il n'est pas utilisé comme argument.
- La table de fiabilité du modèle calibré est fournie dans
  `resultats/metriques_modele.json` (écarts par bac ≤ 0,064 ; max aux
  bacs 4 et 6, n = 353 et 343 — correction du 05/11 après relecture du
  fil projet, la première rédaction disait ~0,03).

## 5. Criblage et score de priorisation

Score transparent, formule écrite dans le code **avant** exécution :

```
score = 0,55 × P calibrée + 0,20 × QED + 0,15 × (1 − violations Lipinski/4)
        + 0,10 × nouveauté d'échafaudage
```

738 790 molécules scorées. Chaque candidat du top 500 est justifié par :
probabilité calibrée, QED, violations de Lipinski, nouveauté, et plus
proche voisin actif ChEMBL (Tanimoto ECFP4, avec son pChEMBL médian).
Le classement brut est saturé de stéréoisomères (classement v1 : top 50
brut = 5 connectivités uniques ; classement v2 courant : top 50 brut =
22 connectivités uniques, top 500 = 214) : la liste de présentation est
**dédoublonnée par connectivité** (InChIKey-14) — 100 candidats
distincts (rangs d'origine ≤ 237 en v2).

## 6. Incident de la v1 et correctif v2 (correction complète)

La première passe de criblage (v1, 05/11 03h59) calculait la nouveauté
en comparant la chaîne `murcko_framework` **fournie par COCONUT** à
l'ensemble des échafaudages **calculés par RDKit** côté ChEMBL. Les deux
générateurs écrivent le même squelette sous des chaînes différentes
(prouvé au rang 1 : candidat CNP0144043.8 — chaîne COCONUT
`n1ccc(c2ccccc12)CC3N4CCC(CC4)C3` vs chaîne RDKit
`c1ccc2c(CC3CC4CCN3CC4)ccnc2c1`, cette dernière identique à
l'échafaudage du voisin actif CHEMBL2092739, Tanimoto exact 1,0). La
comparaison entre générateurs déclarait donc « nouveau » ce qui ne
l'était pas : **nouveauté du top 50 affichée dans les sorties v1 :
1,00 — réelle : 0,00** ; top 500 : 1,0 → 0,13.

Correctif (v2, même jour) : échafaudage du candidat **recalculé par
RDKit** dans le criblage, comparaison à générateur unique. Contrôle
indépendant du fil projet : concordance 1,0000 de la colonne nouveauté
sur les 251 identifiants communs avec le fichier d'audit du correctif.
Conséquences : recouvrement des top 50 v1/v2 = **0/50** ; la famille
quinine (CNP0144043), tête de la v1, ressort au **rang 61** de la v2,
sans aucun bonus de nouveauté (son terme de nouveauté vaut 0), portée
par sa probabilité calibrée (0,9699) et son QED (0,9352). **Nouveauté du top 50 v2 = 1,0, cette fois au
sens corrigé** : le nouveau haut de classement est constitué
d'échafaudages réellement absents de ChEMBL. Toute affirmation de
nouveauté fondée sur la v1 est retirée ; les fichiers v1 restent au
disque comme trace d'audit.

## 7. Candidats et organismes sources

Liste de présentation : `resultats/top100_dedoublonne_v2_organismes.csv`
(100 connectivités uniques, rangs d'origine ≤ 237). Organismes joints
depuis le COCONUT Full du même millésime (voie de repli du fil projet,
la voie ciblée par page web ayant échoué — site en application
JavaScript sans données accessibles) : **15/100 candidats avec organisme
documenté** (ex. rang 5 : CNP0228647.0, ascidiathiazone B — *Aplidium*
sp.) ; top 500 v2 : 99/500 documentés ; top 50 v2 : 5/50. Un candidat
sans organisme documenté est présenté tel quel : **« organisme source
non documenté dans COCONUT »** — jamais d'organisme inféré.

## 8. Limites générales

- Aucune validation expérimentale ; les probabilités sont des
  probabilités *du modèle*, calibrées sur des données d'essais hétérogènes.
- Graine unique (42) pour le split et le modèle ; pas d'intervalle
  d'incertitude inter-graines.
- Les stéréoisomères sont scorés individuellement mais présentés
  dédoublonnés : la stéréochimie peut changer l'activité réelle, le
  modèle ne la distingue que par l'empreinte.
- COCONUT contient des structures rapportées sans organisme (85 % du
  top 100 dédoublonné) : la dimension « produit naturel sourcé » de ces
  candidats repose sur leur présence dans COCONUT, pas sur une fiche
  d'organisme.

## 9. Reproductibilité

Environnement épinglé (`ENVIRONNEMENT.md`) : Python 3.12, RDKit
2026.09.1, scikit-learn 1.9.1, pandas 3.1.0, pyarrow 26.0.1. Code :
`code/calculer_descripteurs_v2.py`, `code/entrainer_modele.py`,
`code/criblage.py` (v2). Sorties : `resultats/metriques_modele.json`,
`resultats/criblage_scores_v2.parquet`,
`resultats/criblage_top500_v2.csv`, `resultats/resume_descripteurs.json`,
`resultats/split_scaffolds.json`, `resultats/chembl_features.parquet`.
Le modèle de criblage est un ré-ajustement déterministe (même graine,
mêmes plis) depuis les features gelées — aucun artefact opaque.
Divulgation IA : ce projet a été construit avec l'assistance d'IA
(pipeline d'agents sous direction humaine) ; l'équipe humaine a revu les
méthodes, vérifié les résultats et répond de cette soumission.
