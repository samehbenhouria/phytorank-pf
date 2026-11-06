#!/usr/bin/env python3
"""PhytoRank-Pf — criblage COCONUT complet (D11, 05/11/2026).

Pipeline figée AVANT tout résultat (même protocole que l'entraînement) :
- modèle : RF 300 arbres, graine 42, représentation combinée ECFP4 +
  descripteurs, étiquette pChEMBL médian >= 6, entraîné sur le pli train
  du split par échafaudage ; calibration Platt ajustée sur le pli val.
  (Re-ajustement déterministe depuis resultats/chembl_features.parquet —
  aucun artefact de modèle à stocker, rejeu exact.)
- score de priorisation TRANSPARENT, formule écrite ici avant exécution :
    score = 0,55 * P_calibrée + 0,20 * QED + 0,15 * (1 - violations_Lipinski/4)
            + 0,10 * nouveauté (1 si l'échafaudage de Murcko est absent
              des molécules ChEMBL d'entraînement, 0 sinon)
- justification du top 200 : plus proche voisin actif ChEMBL (pChEMBL >= 6)
  par Tanimoto sur ECFP4, avec son pChEMBL médian.

Sorties (tmp puis renommage) : resultats/criblage_scores.parquet (toutes
les molécules valides), resultats/criblage_top500.csv,
resultats/criblage_resume.json.
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

RDLogger.DisableLog("rdApp.*")
RACINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
GRAINE = 42
DESC = ["mw", "logp", "tpsa", "hbd", "hba", "rotbonds", "frac_csp3",
        "arom_rings", "qed"]


def matrice_fp(liste_sparse, n):
    X = np.zeros((n, 2048), dtype=np.float32)
    for i, s in enumerate(liste_sparse):
        if isinstance(s, str) and s:
            X[i, [int(x) for x in s.split(",")]] = 1.0
    return X


_gen = None


def _fp(smiles):
    """v2 (correctif 05/11) : retourne (empreinte, Murcko RDKit).

    Le Murcko du candidat est recalculé par RDKit, JAMAIS lu dans la
    colonne murcko_framework de COCONUT : les conventions d'écriture des
    deux sources diffèrent (prouvé au rang 1) et la comparaison de chaînes
    entre générateurs produisait une nouveauté fictive (top 50 publié à
    1,00, corrigé à 0,00)."""
    global _gen
    if _gen is None:
        _gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return (_gen.GetFingerprintAsNumPy(mol),
            MurckoScaffold.MurckoScaffoldSmiles(mol=mol))


def main():
    feat = pd.read_parquet(os.path.join(RACINE, "resultats/chembl_features.parquet"))
    with open(os.path.join(RACINE, "resultats/split_scaffolds.json")) as f:
        pli = json.load(f)["assignation"]
    feat["pli"] = feat.scaffold.map(pli)
    y = (feat.pchembl_median >= 6).astype(int)
    X = np.hstack([matrice_fp(feat.ecfp4, len(feat)),
                   feat[DESC].to_numpy(dtype=np.float32)])
    tr, va = feat.pli == "train", feat.pli == "val"
    clf = RandomForestClassifier(n_estimators=300, class_weight="balanced_subsample",
                                 n_jobs=2, random_state=GRAINE)
    clf.fit(X[tr], y[tr])
    platt = LogisticRegression(C=1.0, max_iter=1000)
    platt.fit(clf.predict_proba(X[va])[:, 1].reshape(-1, 1), y[va])
    print("modèle ré-ajusté (train) + Platt (val)", flush=True)

    scafs_connus = set(feat.scaffold)
    actifs = feat[y == 1]
    fp_actifs = matrice_fp(actifs.ecfp4, len(actifs)) > 0
    ids_actifs = actifs.id.to_numpy()
    pchembl_actifs = actifs.pchembl_median.to_numpy()
    print(f"voisins de référence : {len(actifs)} actifs ChEMBL", flush=True)

    chemin = os.path.join(RACINE, "resultats/descripteurs_coconut.parquet")
    pf = pq.ParquetFile(chemin)
    colonnes = ["identifier", "inchi_key", "smiles", "formule", "murcko",
                "np_pathway", "np_superclass", "np_class", "valide",
                "qed_rdkit", "lipinski_violations"] + DESC[:-1]
    resultats = []
    with Pool(2) as pool:
        for lot in pf.iter_batches(batch_size=40000, columns=colonnes):
            d = lot.to_pandas()
            d = d[d.valide == 1].reset_index(drop=True)
            if len(d) == 0:
                continue
            fps = pool.map(_fp, d.smiles, chunksize=1000)
            ok = np.array([f is not None for f in fps])
            d = d[ok].reset_index(drop=True)
            Xf = np.array([f[0] for f in fps if f is not None], dtype=np.float32)
            murckos_rdkit = [f[1] for f in fps if f is not None]
            Xc = np.hstack([Xf, d[DESC[:-1] + ["qed_rdkit"]].to_numpy(dtype=np.float32)])
            p_brut = clf.predict_proba(Xc)[:, 1]
            p_cal = platt.predict_proba(p_brut.reshape(-1, 1))[:, 1]
            nouv = np.array([0.0 if m in scafs_connus else 1.0
                             for m in murckos_rdkit])
            score = (0.55 * p_cal + 0.20 * d.qed_rdkit.to_numpy()
                     + 0.15 * (1 - d.lipinski_violations.to_numpy() / 4.0)
                     + 0.10 * nouv)
            sortie = pd.DataFrame({
                "identifier": d.identifier, "inchi_key": d.inchi_key,
                "smiles": d.smiles, "formule": d.formule, "murcko": d.murcko,
                "np_pathway": d.np_pathway, "np_superclass": d.np_superclass,
                "np_class": d.np_class, "qed": d.qed_rdkit,
                "lipinski_violations": d.lipinski_violations,
                "p_calibree": p_cal, "nouveaute": nouv, "score": score})
            sortie["_fp"] = [",".join(map(str, row.nonzero()[0])) for row in Xf]
            resultats.append(sortie)
            print(f"lot traité — cumul {sum(len(r) for r in resultats)} molécules",
                  flush=True)
    tout = pd.concat(resultats, ignore_index=True)
    print(f"criblage terminé : {len(tout)} molécules scorées", flush=True)

    tmp = os.path.join(RACINE, "resultats/criblage_scores_v2.parquet.tmp")
    tout.drop(columns=["_fp"]).to_parquet(tmp)
    os.replace(tmp, os.path.join(RACINE, "resultats/criblage_scores_v2.parquet"))

    top = tout.sort_values("score", ascending=False).head(500).copy()
    # justification : plus proche voisin actif pour le top 200
    lignes = []
    for rang, (_, r) in enumerate(top.iterrows(), start=1):
        vo_id, vo_tan, vo_pc = "", np.nan, np.nan
        if rang <= 200:
            fp = np.zeros(2048, dtype=bool)
            if r._fp:
                fp[[int(x) for x in r._fp.split(",")]] = True
            inter = fp_actifs[:, fp].sum(axis=1) if fp.any() else np.zeros(len(actifs))
            union = fp_actifs.sum(axis=1) + fp.sum() - inter
            tan = np.where(union > 0, inter / np.maximum(union, 1), 0.0)
            j = int(np.argmax(tan))
            vo_id, vo_tan, vo_pc = ids_actifs[j], float(tan[j]), float(pchembl_actifs[j])
        lignes.append({"rang": rang, "identifier": r.identifier,
                       "inchi_key": r.inchi_key, "smiles": r.smiles,
                       "formule": r.formule, "score": round(float(r.score), 4),
                       "p_calibree": round(float(r.p_calibree), 4),
                       "qed": round(float(r.qed), 4),
                       "lipinski_violations": int(r.lipinski_violations),
                       "nouveaute_scaffold": int(r.nouveaute),
                       "murcko": r.murcko, "np_pathway": r.np_pathway,
                       "np_superclass": r.np_superclass, "np_class": r.np_class,
                       "voisin_actif_chembl": vo_id,
                       "tanimoto_voisin": round(float(vo_tan), 3) if vo_tan == vo_tan else "",
                       "pchembl_voisin": vo_pc if vo_pc == vo_pc else ""})
    top_df = pd.DataFrame(lignes)
    tmp = os.path.join(RACINE, "resultats/criblage_top500_v2.csv.tmp")
    top_df.to_csv(tmp, index=False)
    os.replace(tmp, os.path.join(RACINE, "resultats/criblage_top500_v2.csv"))

    res = {"version": "v2 — correctif nouveauté du 05/11 (Murcko RDKit des deux côtés) ; "
                      "remplace la v1 dont la nouveauté comparait des chaînes de deux générateurs",
           "n_scorées": int(len(tout)),
           "formule_score": "0,55*P_cal + 0,20*QED + 0,15*(1-viol/4) + 0,10*nouveauté",
           "score_max": float(tout.score.max()),
           "score_top50_min": float(top.score.iloc[49]) if len(top) >= 50 else None,
           "part_nouveaute_top50": float(top.nouveaute.head(50).mean())}
    tmp = os.path.join(RACINE, "resultats/criblage_resume_v2.json.tmp")
    with open(tmp, "w") as f:
        json.dump(res, f, indent=2)
    os.replace(tmp, os.path.join(RACINE, "resultats/criblage_resume_v2.json"))
    print(json.dumps(res, indent=2), flush=True)


if __name__ == "__main__":
    main()
