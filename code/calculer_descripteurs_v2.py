#!/usr/bin/env python3
"""PhytoRank-Pf — descripteurs v2, ÉCRITURE DURCIE (D11, 03/11/2026).

Leçon de la passe du 02/11 : un fichier unique écrit en flux a été tronqué
par une interruption machine — illisible, et pris pour complet. Ici :

- calcul par blocs de 20 000, chaque bloc = un fichier parquet PARTIE
  fermé proprement, écrit sous nom temporaire puis renommé (atomique) ;
- progression journalisée (descripteurs_progres.jsonl) après chaque bloc ;
- reprise : un bloc dont la partie existe ET se relit avec le bon nombre
  de lignes est sauté, jamais recalculé ;
- assemblage séparé : le fichier final n'est renommé à son nom définitif
  qu'après relecture de contrôle (pied valide + nombre de lignes égal au
  nombre de lignes de données du CSV source).

Usage :
  calculer_descripteurs_v2.py calculer ENTREE.csv REPERTOIRE_SORTIE [NB_PROC]
  calculer_descripteurs_v2.py assembler ENTREE.csv REPERTOIRE_SORTIE FINAL.parquet RESUME.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calculer_descripteurs import COLS, SCHEMA, traiter  # noqa: E402

import pandas as pd  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
from multiprocessing import Pool  # noqa: E402

TAILLE_BLOC = 20000


def compter_lignes_donnees(chemin):
    n = 0
    with open(chemin, "rb") as f:
        for _ in f:
            n += 1
    return n - 1  # en-tête


def partie_valide(chemin, lignes_attendues):
    try:
        return pq.ParquetFile(chemin).metadata.num_rows == lignes_attendues
    except Exception:
        return False


def calculer(entree, rep, nb_proc):
    os.makedirs(rep, exist_ok=True)
    journal = os.path.join(os.path.dirname(rep), "descripteurs_progres.jsonl")
    attendu_total = compter_lignes_donnees(entree)
    print(f"lignes de données attendues : {attendu_total}", flush=True)
    with Pool(nb_proc) as pool:
        for i, bloc in enumerate(pd.read_csv(
                entree, usecols=COLS, chunksize=TAILLE_BLOC,
                dtype=str, keep_default_na=False)):
            nom = os.path.join(rep, f"bloc_{i:03d}.parquet")
            if partie_valide(nom, len(bloc)):
                print(f"bloc {i} déjà valide ({len(bloc)} lignes) — sauté", flush=True)
                continue
            qed_f = pd.to_numeric(bloc["qed_drug_likeliness"], errors="coerce")
            mw_f = pd.to_numeric(bloc["molecular_weight"], errors="coerce")
            lignes = list(zip(bloc["identifier"], bloc["canonical_smiles"],
                              bloc["standard_inchi_key"], bloc["molecular_formula"],
                              bloc["murcko_framework"], bloc["np_classifier_pathway"],
                              bloc["np_classifier_superclass"], bloc["np_classifier_class"],
                              qed_f, mw_f))
            resultats = pool.map(traiter, lignes, chunksize=500)
            table = pa.Table.from_pylist(
                [dict(zip(SCHEMA.names, r)) for r in resultats], schema=SCHEMA)
            tmp = nom + ".tmp"
            pq.write_table(table, tmp)
            os.replace(tmp, nom)
            with open(journal, "a") as j:
                j.write(json.dumps({"bloc": i, "lignes": len(bloc)}) + "\n")
            print(f"bloc {i} écrit et fermé ({len(bloc)} lignes)", flush=True)
    print("calcul terminé — lancer le mode assembler pour le fichier final", flush=True)


def assembler(entree, rep, final, resume):
    import numpy as np
    attendu = compter_lignes_donnees(entree)
    parties = sorted(f for f in os.listdir(rep)
                     if f.startswith("bloc_") and f.endswith(".parquet"))
    if not parties:
        sys.exit("aucune partie à assembler")
    tmp = final + ".tmp"
    writer = pq.ParquetWriter(tmp, SCHEMA)
    total = 0
    for p in parties:
        t = pq.read_table(os.path.join(rep, p))
        writer.write_table(t)
        total += t.num_rows
    writer.close()
    # Relecture de contrôle AVANT tout renommage définitif.
    controle = pq.ParquetFile(tmp).metadata.num_rows
    if controle != attendu or total != attendu:
        sys.exit(f"CONTRÔLE ÉCHOUÉ : assemblé {total}, relu {controle}, attendu {attendu} "
                 f"— fichier final NON publié (reste {tmp})")
    df = pq.read_table(tmp, columns=["valide", "qed_fourni", "qed_rdkit",
                                     "mw_fourni", "mw"]).to_pandas()
    v = df[df.valide == 1].dropna(subset=["qed_fourni", "qed_rdkit"])
    d = (v.qed_fourni - v.qed_rdkit).abs()
    w = df[df.valide == 1].dropna(subset=["mw_fourni", "mw"])
    res = {
        "n_total": int(attendu),
        "n_valides": int(df.valide.sum()),
        "n_invalides": int((df.valide == 0).sum()),
        "controle_lignes": int(controle),
        "qed_comparaison": {
            "n": int(len(v)),
            "pearson": float(v.qed_fourni.corr(v.qed_rdkit)),
            "spearman": float(v.qed_fourni.corr(v.qed_rdkit, method="spearman")),
            "ecart_abs_median": float(d.median()),
            "ecart_abs_moyen": float(d.mean()),
            "ecart_abs_p95": float(d.quantile(0.95)),
            "ecart_abs_max": float(d.max()),
            "frac_ecart_le_0_01": float((d <= 0.01).mean()),
            "frac_ecart_gt_0_1": float((d > 0.1).mean()),
        },
        "mw_comparaison": {
            "n": int(len(w)),
            "ecart_abs_median": float((w.mw_fourni - w.mw).abs().median()),
            "ecart_abs_max": float((w.mw_fourni - w.mw).abs().max()),
        },
    }
    with open(resume, "w") as f:
        json.dump(res, f, indent=2)
    os.replace(tmp, final)
    print(json.dumps(res, indent=2))
    print(f"FICHIER FINAL PUBLIÉ : {final} ({controle} lignes vérifiées)", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "calculer":
        calculer(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 2)
    elif mode == "assembler":
        assembler(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        sys.exit("mode inconnu : calculer | assembler")
