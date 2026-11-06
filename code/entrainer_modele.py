#!/usr/bin/env python3
"""PhytoRank-Pf — jeu d'entraînement figé + baselines + modèle (D11).

Protocole gelé avant entraînement (CHANTIER_ETHNOHACK D11) :
- source : data/chembl_pf_activities.csv (CHEMBL364, IC50/EC50, pChEMBL médian)
- étiquette : actif si pChEMBL médian >= seuil ; seuil de référence 6,
  sensibilité rapportée à 5 et 7
- split par échafaudage de Bemis–Murcko (assertion de disjonction),
  jamais de split aléatoire présenté comme validation
- représentations comparées : (a) ECFP4 seul, (b) descripteurs seuls,
  (c) ECFP4 + descripteurs (candidat obligatoire, mesures front E D11)
- métriques : ROC-AUC, PR-AUC publiée AVEC sa prévalence, précision@k,
  enrichissement à 1 %, Brier + calibration Platt ajustée sur le pli val
Sorties : resultats/chembl_features.parquet, resultats/split_scaffolds.json,
resultats/metriques_modele.json — écrits en tmp puis renommés.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, QED, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, brier_score_loss,
                             precision_score, roc_auc_score)

RDLogger.DisableLog("rdApp.*")
RACINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
GRAINE = 42

DESC = ["mw", "logp", "tpsa", "hbd", "hba", "rotbonds", "frac_csp3",
        "arom_rings", "qed"]


def descripteurs(mol):
    mw = Descriptors.MolWt(mol)
    return [float(mw), float(Descriptors.MolLogP(mol)),
            float(rdMolDescriptors.CalcTPSA(mol)), float(Lipinski.NumHDonors(mol)),
            float(Lipinski.NumHAcceptors(mol)), float(Lipinski.NumRotatableBonds(mol)),
            float(Lipinski.FractionCSP3(mol)), float(Lipinski.NumAromaticRings(mol)),
            float(QED.qed(mol))]


def main():
    src = pd.read_csv(os.path.join(RACINE, "data/chembl_pf_activities.csv"))
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    lignes, fps = [], []
    invalides = 0
    for r in src.itertuples(index=False):
        mol = Chem.MolFromSmiles(r.canonical_smiles)
        if mol is None:
            invalides += 1
            continue
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        lignes.append([r.molecule_chembl_id, r.canonical_smiles, scaf,
                       float(r.median_pchembl_value)] + descripteurs(mol))
        fps.append(gen.GetFingerprintAsNumPy(mol))
    df = pd.DataFrame(lignes, columns=["id", "smiles", "scaffold",
                                       "pchembl_median"] + DESC)
    X_fp = np.array(fps, dtype=np.int8)
    print(f"molécules valides {len(df)} | invalides {invalides} | "
          f"échafaudages uniques {df.scaffold.nunique()}", flush=True)

    # --- split par échafaudage, glouton déterministe 70/15/15 ---
    groupes = df.groupby("scaffold").size().sort_values(ascending=False)
    cibles = {"train": 0.70 * len(df), "val": 0.15 * len(df), "test": 0.15 * len(df)}
    comptes = {"train": 0, "val": 0, "test": 0}
    pli = {}
    for scaf, taille in groupes.items():
        choix = max(cibles, key=lambda k: (cibles[k] - comptes[k]) / cibles[k])
        pli[scaf] = choix
        comptes[choix] += taille
    df["pli"] = df.scaffold.map(pli)
    ens = {k: set(df.loc[df.pli == k, "scaffold"]) for k in cibles}
    assert not (ens["train"] & ens["val"] or ens["train"] & ens["test"]
                or ens["val"] & ens["test"]), "FUITE : échafaudage partagé entre plis"
    split_info = {k: int((df.pli == k).sum()) for k in cibles}
    print("split :", split_info, flush=True)

    os.makedirs(os.path.join(RACINE, "resultats"), exist_ok=True)
    feat = df.copy()
    feat["ecfp4"] = [",".join(map(str, row.nonzero()[0])) for row in X_fp]
    df = feat  # correctif 04/11 : evaluer/X_de lisent la colonne ecfp4, portée par feat
    tmp = os.path.join(RACINE, "resultats/chembl_features.parquet.tmp")
    feat.to_parquet(tmp)
    os.replace(tmp, os.path.join(RACINE, "resultats/chembl_features.parquet"))
    with open(os.path.join(RACINE, "resultats/split_scaffolds.json"), "w") as f:
        json.dump({"graine": GRAINE, "tailles": split_info,
                   "assertion_disjonction": "vérifiée dans le code",
                   "assignation": {s: p for s, p in sorted(pli.items())}}, f)

    def X_de(df_, rep):
        fp = np.zeros((len(df_), 2048), dtype=np.float32)
        for i, s in enumerate(df_.ecfp4):
            if s:
                fp[i, [int(x) for x in s.split(",")]] = 1.0
        if rep == "ecfp":
            return fp
        if rep == "desc":
            return df_[DESC].to_numpy(dtype=np.float32)
        return np.hstack([fp, df_[DESC].to_numpy(dtype=np.float32)])

    def ef_1pct(y, p):
        n = max(1, int(0.01 * len(y)))
        top = np.argsort(-p)[:n]
        return float(y.iloc[top].mean() / y.mean()) if y.mean() > 0 else 0.0

    def evaluer(seuil, rep, calibrer=False):
        y = (df.pchembl_median >= seuil).astype(int)
        tr, va, te = df.pli == "train", df.pli == "val", df.pli == "test"
        clf = RandomForestClassifier(n_estimators=300, class_weight="balanced_subsample",
                                     n_jobs=2, random_state=GRAINE)
        clf.fit(X_de(df[tr], rep), y[tr])
        pv = clf.predict_proba(X_de(df[va], rep))[:, 1]
        pt = clf.predict_proba(X_de(df[te], rep))[:, 1]
        res = {"seuil": seuil, "representation": rep,
               "prevalence_test": float(y[te].mean()),
               "roc_auc_test": float(roc_auc_score(y[te], pt)),
               "pr_auc_test": float(average_precision_score(y[te], pt)),
               "precision_at_50": float(y[te].iloc[np.argsort(-pt)[:50]].mean()),
               "precision_at_100": float(y[te].iloc[np.argsort(-pt)[:100]].mean()),
               "ef_1pct": ef_1pct(y[te].reset_index(drop=True),
                                  pd.Series(pt).reset_index(drop=True)),
               "brier_brut": float(brier_score_loss(y[te], pt))}
        if calibrer:
            platt = LogisticRegression(C=1.0, max_iter=1000)
            platt.fit(pv.reshape(-1, 1), y[va])
            pc = platt.predict_proba(pt.reshape(-1, 1))[:, 1]
            res["brier_platt"] = float(brier_score_loss(y[te], pc))
            bacs = np.clip((pc * 10).astype(int), 0, 9)
            rel = []
            for b in range(10):
                m = bacs == b
                if m.sum() >= 10:
                    rel.append({"bac": b, "n": int(m.sum()),
                                "proba_moy": float(pc[m].mean()),
                                "freq_obs": float(y[te].to_numpy()[m].mean())})
            res["fiabilite_platt"] = rel
        return res

    resultats = []
    for rep in ["ecfp", "desc", "combine"]:
        r = evaluer(6, rep, calibrer=(rep == "combine"))
        resultats.append(r)
        print(f"seuil 6 | {rep:8s} | ROC {r['roc_auc_test']:.4f} | "
              f"PR {r['pr_auc_test']:.4f} (prév {r['prevalence_test']:.3f}) | "
              f"P@100 {r['precision_at_100']:.3f} | EF1% {r['ef_1pct']:.2f}", flush=True)
    for seuil in [5, 7]:
        r = evaluer(seuil, "combine")
        resultats.append(r)
        print(f"seuil {seuil} | combine  | ROC {r['roc_auc_test']:.4f} | "
              f"PR {r['pr_auc_test']:.4f} (prév {r['prevalence_test']:.3f})", flush=True)

    sortie = {"protocole": "scaffold split Bemis-Murcko, assertion de disjonction, "
                           "étiquette = pChEMBL médian >= seuil, RF 300 arbres, graine 42",
              "n_molecules": int(len(df)), "n_invalides": int(invalides),
              "split": split_info, "resultats": resultats}
    tmp = os.path.join(RACINE, "resultats/metriques_modele.json.tmp")
    with open(tmp, "w") as f:
        json.dump(sortie, f, indent=2)
    os.replace(tmp, os.path.join(RACINE, "resultats/metriques_modele.json"))
    print("ÉCRIT : resultats/metriques_modele.json", flush=True)


if __name__ == "__main__":
    main()
