#!/usr/bin/env python3
"""Fetch ChEMBL activities for Plasmodium falciparum (CHEMBL364).

Filters (server-side, re-checked client-side):
  target_chembl_id=CHEMBL364, standard_type in {IC50, EC50},
  standard_units=nM, standard_value not null, pchembl_value not null.
Activity records carry canonical_smiles directly, so no molecule
lookup is needed. Pages are fetched concurrently (totals known upfront).
Output CSV deduplicated by molecule_chembl_id: median pChEMBL, best
(lowest) standard_value, n_activities, assay/type info, one SMILES.
"""
import csv, json, statistics, sys, time, urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

BASE = "https://www.ebi.ac.uk/chembl/api/data"
OUT = sys.argv[1] if len(sys.argv) > 1 else "chembl_pf_activities.csv"
LIMIT = 1000

def get(url, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(3 * (i + 1))

def page_url(stype, offset):
    return (f"{BASE}/activity.json?target_chembl_id=CHEMBL364"
            f"&standard_type={stype}&standard_units=nM"
            f"&pchembl_value__isnull=false&limit={LIMIT}&offset={offset}")

def fetch_activities():
    # Sequential pagination: concurrent requests (4 workers) were tried on
    # 2026-10-31 and repeatedly hit RemoteDisconnected from the API/proxy;
    # the sequential run is slower (~18 min) but reliable.
    acts = []
    for stype in ("IC50", "EC50"):
        offset = 0
        kept = 0
        while True:
            d = get(page_url(stype, offset))
            batch = d["activities"]
            total = d["page_meta"]["total_count"]
            for a in batch:
                if a.get("standard_value") is None or a.get("pchembl_value") is None:
                    continue
                if a.get("standard_units") != "nM" or not a.get("canonical_smiles"):
                    continue
                acts.append(a)
                kept += 1
            offset += len(batch)
            print(f"{stype}: {offset}/{total} kept={kept}", flush=True)
            if offset >= total or not batch:
                break
    return acts

def main():
    acts = fetch_activities()
    print(f"raw filtered activities: {len(acts)}")
    by_mol = defaultdict(list)
    smiles = {}
    for a in acts:
        by_mol[a["molecule_chembl_id"]].append(a)
        smiles.setdefault(a["molecule_chembl_id"], a["canonical_smiles"])
    print(f"unique molecules: {len(by_mol)}")
    rows = []
    for mid, lst in by_mol.items():
        pvals = [float(a["pchembl_value"]) for a in lst]
        svals = [float(a["standard_value"]) for a in lst]
        types = sorted({a["standard_type"] for a in lst})
        assays = sorted({a.get("assay_chembl_id") for a in lst if a.get("assay_chembl_id")})
        rows.append({
            "molecule_chembl_id": mid,
            "canonical_smiles": smiles[mid],
            "n_activities": len(lst),
            "standard_types": ";".join(types),
            "median_pchembl_value": round(statistics.median(pvals), 4),
            "max_pchembl_value": max(pvals),
            "min_standard_value_nM": min(svals),
            "median_standard_value_nM": statistics.median(svals),
            "n_assays": len(assays),
            "assay_chembl_ids": ";".join(assays[:10]),
        })
    rows.sort(key=lambda r: -r["median_pchembl_value"])
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {OUT}")

if __name__ == "__main__":
    main()
