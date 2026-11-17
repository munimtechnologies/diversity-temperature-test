# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Run the pre-registered test in PREREG.md end to end.

    python job_diversity.py --data-dir /data --out-dir /out

Downloads (or reuses) PBDB marine invertebrate genus occurrences and the PhanDA
GMST percentiles (Judd et al. 2024, CC BY 4.0), computes coverage-standardised
genus richness and three-timer origination on PhanDA's time slices, and runs
the registered tests.
Made with AI agents working in MV2's name.
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

import diversity
import job_thermo as jt
import thermo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261117)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    prov = {"started_utc": jt.now(), "pbdb": {}, "phanda": {}, "groups": jt.GROUPS,
            "query_common": jt.COMMON, "n_perm": args.n_perm, "seed": args.seed,
            "quorum": diversity.QUORUM, "quorums_alt": list(diversity.QUORUMS_ALT)}
    res = {"prereg": "PREREG.md (diversity-temperature-test v0.1.0)", "status": "ok"}

    try:
        temps = jt.get_phanda(args.data_dir, prov["phanda"])
    except Exception as e:  # noqa: BLE001
        prov["phanda"]["error"] = repr(e)[:500]
        temps = {}
    got_slices = set(temps)
    prov["phanda"]["n_slices"] = len(got_slices)
    prov["phanda"]["missing_expected"] = sorted(set(thermo.PHANDA_SLICES) - got_slices)
    prov["phanda"]["unexpected"] = sorted(got_slices - set(thermo.PHANDA_SLICES))
    bins = thermo.build_bins()
    needed = {b["phanda"] for b in bins if b["phanda"]}
    if not needed <= got_slices:
        res["status"] = "phanda_incomplete: missing " + ",".join(sorted(needed - got_slices))
        return _write(args.out_dir, res, prov, 2)

    paths = jt.get_pbdb(args.data_dir, prov["pbdb"])
    got = {g for g, _ in paths}
    if not jt.REQUIRED <= got:
        res["status"] = "pbdb_incomplete: missing " + ",".join(sorted(jt.REQUIRED - got))
        return _write(args.out_dir, res, prov, 2)

    assigner = thermo.BinAssigner()
    acc = diversity.CountAccumulator(assigner.n)
    for _, p in paths:
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            acc.add_rows(csv.DictReader(f), assigner)
    prov["parse"] = acc.stats

    counts, n_coll = diversity.merge_counts(acc.counts, bins, acc.n_coll)
    sets = [set(c) for c in counts]
    lam = diversity.three_timer_origination(sets)
    quorums = (diversity.QUORUM,) + tuple(diversity.QUORUMS_ALT)
    cr = {q: [diversity.coverage_richness(list(c.values()), q) for c in counts]
          for q in quorums}
    rich = {q: [d["richness"] for d in cr[q]] for q in quorums}
    raw = [len(c) for c in counts]
    rows = diversity.level_rows(bins, rich, raw, n_coll, lam, temps)

    with open(os.path.join(args.out_dir, "bins.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["bin", "stages", "phanda_slice", "base_ma", "top_ma", "n_collections",
                    "n_genus_collection_records", "genera_observed", "coverage_estimate",
                    "m_at_0.7", "R_0.5", "R_0.7", "R_0.8", "origination_3t_base",
                    "GMST_16", "GMST_50", "GMST_84"])
        for j, b in enumerate(bins):
            t = temps.get(b["phanda"], ("", "", "")) if b["phanda"] else ("", "", "")
            c7 = cr[diversity.QUORUM][j]
            w.writerow([b["name"], "+".join(b["stages"]), b["phanda"] or "", b["base"], b["top"],
                        n_coll[j], c7["n"], c7["s_obs"], c7["coverage"], c7["m"],
                        rich[0.5][j], rich[diversity.QUORUM][j], rich[0.8][j], lam[j],
                        t[0], t[1], t[2]])

    d = diversity.diff_rows(rows, "R")
    with open(os.path.join(args.out_dir, "differences.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["boundary_top_of", "dlogR_0.7", "dT", "big_five"])
        for r in d:
            w.writerow([r["bin"], r["dlogR"], r["dT"], r["big_five"]])

    res["analysis"] = diversity.analyze(rows, args.n_perm, rng)
    res["decision"] = res["analysis"]["decision"]
    prov["finished_utc"] = jt.now()
    return _write(args.out_dir, res, prov, 0)


def _write(out_dir, res, prov, code):
    def conv(o):
        return o.item() if hasattr(o, "item") else str(o)
    with open(os.path.join(out_dir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=conv)
    with open(os.path.join(out_dir, "provenance.json"), "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=1, default=str)
    a = res.get("analysis", {})
    print(json.dumps({"status": res.get("status"), "decision": res.get("decision"),
                      "sanity": a.get("sanity"), "D1": a.get("D1_levels"),
                      "D2": a.get("D2_differences")}, default=conv), flush=True)
    return code


if __name__ == "__main__":
    sys.exit(main())
