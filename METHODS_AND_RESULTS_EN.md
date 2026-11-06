# PhytoRank-Pf — Methods and Results Report (D11 seat, 2026-11-06)

Computational prioritization of natural products against *Plasmodium
falciparum* (EthnoHack 2026, Track 2). **Exact scope: prioritization
hypotheses for expert review — no claim of biological activity, efficacy,
or medical use.** No candidate was experimentally tested in this project.

*(English translation of the original French report,
`RAPPORT_METHODES_RESULTATS.md`; the French text governs.)*

## 1. Data

| Dataset | Source | Content |
|---|---|---|
| COCONUT Lite 10-2026 | COCONUT 2.0 (official download, MD5 verified `c86adbd2…`) | 738,794 molecules, 40 columns; source organisms absent from the Lite release |
| *P. falciparum* activities | ChEMBL, target CHEMBL364, IC50/EC50 | 31,652 unique molecules, median pChEMBL per molecule (median 6.19) |
| Source organisms of the top-N | COCONUT Full 10-2026 (same release), local join by the project thread | see §7 |

Training-data limitation, stated up front: CHEMBL364 aggregates
whole-organism assays — mixed strains and protocols. The model learns
"resembles molecules reported active in those assays", nothing more.

## 2. Protocol (frozen before training)

- Label: active if median pChEMBL ≥ threshold; reference threshold **6**,
  sensitivity reported at 5 and 7.
- Validation: **Bemis–Murcko scaffold split** (70/15/15, greedy
  deterministic, seed 42) with a **disjointness assertion** executed in the
  code — no scaffold shared between folds. No random split is presented as
  validation (domain measurements on the BBBP benchmark, Certifications
  D11 front E: moving from a random split to a scaffold split, ROC-AUC
  drops by 0.116 for ECFP4 + random forest and by 0.044 for a
  descriptor-based model).
- Representations compared: ECFP4 (Morgan r=2, 2048 bits) alone; RDKit
  descriptors alone (mass, logP, TPSA, HBD, HBA, rotatable bonds, fraction
  Csp3, aromatic rings, QED); **combined**.
- Model: random forest, 300 trees, `class_weight=balanced_subsample`.
- Calibration: Platt (logistic regression) fitted on the validation fold
  only.
- Metrics: ROC-AUC, PR-AUC **always published with its prevalence**,
  precision@50/@100, enrichment at 1%, Brier score before/after
  calibration, reliability table.

## 3. Recomputed COCONUT descriptors — validation

RDKit 2026.09.1 recomputation of the descriptors and QED of all 738,794
molecules (pass in blocks closed one by one; assembly published only after
a control re-read: 738,794 rows, valid file footer): **738,790 valid, 4
unreadable**. Recomputed QED vs COCONUT-supplied QED: Pearson **0.99993**,
median absolute gap **0.0027**, 99.99% of gaps ≤ 0.01 (2 molecules > 0.1,
max 0.16) — the supplied values are consistent at rounding level; masses:
median gap 0.002.

## 4. Model results (test fold, scaffold split, n = 4,748)

Test-fold prevalence at threshold 6: **0.607** (chance therefore equals
0.607 in precision and in PR-AUC).

| Representation | ROC-AUC | PR-AUC | P@50 | P@100 | Brier |
|---|---|---|---|---|---|
| ECFP4 alone | 0.8700 | 0.9076 | 0.98 | 0.99 | 0.1470 |
| Descriptors alone | 0.7444 | 0.8204 | 0.98 | 0.98 | 0.1989 |
| **Combined** | **0.8710** | **0.9093** | **0.98** | **0.98** | 0.1477 → **0.1420** after Platt |

Sensitivity (combined): threshold 5 → ROC 0.9084 (prevalence 0.897);
threshold 7 → ROC 0.8779 (prevalence 0.262).

Honest reading:
- The combined model beats ECFP alone by a hair (+0.001) on this dataset —
  unlike other domain benchmarks where descriptors help much more. We keep
  the combined model (protocol) without exaggerating its contribution.
- Enrichment at 1% is **saturated** (1.61 for all three representations,
  theoretical ceiling 1/0.607 = 1.65): non-discriminative at this
  prevalence, it is not used as an argument.
- The calibrated model's reliability table is provided in
  `results/metriques_modele.json` (per-bin gaps ≤ 0.064; maximum at bins 4
  and 6, n = 353 and 343 — corrected note, 2026-11-05: as flagged by the
  project thread, the first draft said ~0.03).

## 5. Screening and prioritization score

Transparent score, formula written in the code **before** execution:

```
score = 0.55 × calibrated P + 0.20 × QED + 0.15 × (1 − Lipinski violations/4)
        + 0.10 × scaffold novelty
```

738,790 molecules scored. Each top-500 candidate is justified by:
calibrated probability, QED, Lipinski violations, novelty, and nearest
active ChEMBL neighbour (ECFP4 Tanimoto, with its median pChEMBL). The raw
ranking is saturated with stereoisomers (v1 ranking: raw top 50 = 5 unique
connectivities; current v2 ranking: raw top 50 = 22 unique connectivities,
top 500 = 214): the presentation list is **deduplicated by connectivity**
(InChIKey-14) — 100 distinct candidates (original ranks ≤ 237 in v2).

## 6. The v1 incident and the v2 fix (complete correction)

The first screening pass (v1, 2026-11-05 03:59) computed novelty by
comparing the `murcko_framework` string **supplied by COCONUT** against
the set of scaffolds **computed by RDKit** on the ChEMBL side. The two
generators write the same skeleton under different strings (proven at
rank 1: candidate CNP0144043.8 — COCONUT string
`n1ccc(c2ccccc12)CC3N4CCC(CC4)C3` vs RDKit string
`c1ccc2c(CC3CC4CCN3CC4)ccnc2c1`, the latter identical to the scaffold of
the active neighbour CHEMBL2092739, exact Tanimoto 1.0). The
cross-generator comparison therefore declared "novel" what was not:
**top-50 novelty as displayed in the v1 outputs: 1.00 — real: 0.00**;
top 500: 1.0 → 0.13.

Fix (v2, same day): the candidate's scaffold **recomputed by RDKit** in
the screening, single-generator comparison. Independent check by the
project thread: 1.0000 agreement of the novelty column over the 251
shared identifiers with the correction's audit file. Consequences:
top-50 overlap v1/v2 = **0/50**; the quinine family (CNP0144043), the
v1's top hit, comes out at **rank 61** in v2 — it stays in the elite
through its activity, QED and Lipinski terms, with a novelty bonus of
exactly zero, which is presentable and verified. **Top-50 v2
novelty = 1.0, this time in the corrected sense**: the new top of the
ranking is made of scaffolds genuinely absent from ChEMBL. Every novelty
claim based on v1 is withdrawn; the v1 files remain on disk as an audit
trail.

## 7. Candidates and source organisms

Presentation list: `results/top100_dedoublonne_v2_organismes.csv` (100
unique connectivities, original ranks ≤ 237). Organisms joined from the
COCONUT Full of the same release (fallback route of the project thread —
the targeted per-page web route failed: JavaScript application site with
no accessible data): **15/100 candidates with a documented organism**
(e.g. rank 5: CNP0228647.0, ascidiathiazone B — *Aplidium* sp.); v2 top
500: 99/500 documented; v2 top 50: 5/50. A candidate without a documented
organism is presented as such: **"source organism not documented in
COCONUT"** — an organism is never inferred.

## 8. General limits

- No experimental validation; the probabilities are *model* probabilities,
  calibrated on heterogeneous assay data.
- Single seed (42) for the split and the model; no cross-seed uncertainty
  interval.
- Stereoisomers are scored individually but presented deduplicated:
  stereochemistry can change real activity; the model distinguishes it
  only through the fingerprint.
- COCONUT contains structures reported without an organism (85% of the
  deduplicated top 100): for those candidates, the "sourced natural
  product" dimension rests on their presence in COCONUT, not on an
  organism record.

## 9. Reproducibility

Pinned environment (`ENVIRONNEMENT.md`): Python 3.12, RDKit 2026.09.1,
scikit-learn 1.9.1, pandas 3.1.0, pyarrow 26.0.1. Code:
`code/calculer_descripteurs_v2.py`, `code/entrainer_modele.py`,
`code/criblage.py` (v2). Outputs: `results/metriques_modele.json`,
`results/criblage_scores_v2.parquet`,
`results/criblage_top500_v2.csv`, `results/resume_descripteurs.json`,
`results/split_scaffolds.json`, `results/chembl_features.parquet`. The
screening model is a deterministic refit (same seed, same folds) from the
frozen features — no opaque artefact.
AI disclosure: this project was built with AI assistance (an agent
pipeline under human direction); the human team reviewed the methods,
verified the results, and is responsible for this submission.
