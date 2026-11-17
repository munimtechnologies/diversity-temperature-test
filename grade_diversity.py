# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""grade_diversity.py: the mechanical reading of diversity-temperature-test.

Written and published 2026-11-17 (about 18:30 UTC), BEFORE the analysis job
(diversity-temperature-v1) had produced any result and before any PBDB
richness or origination rate was paired with any PhanDA temperature.

It turns the job's results.json into the graded prior table of PREREG.md
("Hypotheses, tests and honest priors") and the registered decision label,
row by row. No new statistic is computed; it only reads numbers the job
wrote. Reading choices fixed now, in writing (PRE-RESULT-NOTE-1.md):

  * D1+/D1-/D5a use the one-sided circular-shift p (p_shift / p_shift_lower);
    D2+/D2-/D5b use the one-sided permutation p (p_perm / p_perm_lower), as in
    the PREREG table. A row holds only if rho is finite, has the stated sign,
    and the p value in that direction is < 0.05.
  * D2b ("significant in either direction") holds if the D2 test without the
    Big Five boundaries has rho > 0 with p_perm < 0.05 OR rho < 0 with
    p_perm_lower < 0.05. It is NOT read as a two-sided test at 0.10.
  * D3 holds only if the duration-partial rho and the D1 rho are both finite,
    non-zero and of the same sign.
  * D4 holds only if D1's rho is finite and non-zero and at least 90% of the
    PhanDA draws give that sign.
  * D6 holds only if, at quorums 0.5 and 0.8, the levels rho has D1's sign and
    the differences rho has D2's sign (all non-zero).
  * D7a holds if rho_raw < rho_R (both finite); D7b if rho_raw < 0. As
    registered, these two rows have no significance requirement.

The script re-derives the label from its own row grades and compares it, and
the four signals D1+/D1-/D2+/D2-, with what the job wrote. Any disagreement
gives status CHECK_CODE. Optional: --bins bins.csv checks the job's PhanDA
medians against the GMST_50 values read by hand on 2026-11-17 (tolerance
0.011 C); a mismatch gives status CHECK_DATA. Neither status changes the
label; both stop the result from being reported until a person has looked.

    python grade_diversity.py --results results.json [--bins bins.csv] --out-dir .

Outputs: graded.json and GRADED.md. Made with AI agents working in MV2's name.
"""
import argparse
import csv
import json
import math
import os

ALPHA = 0.05
QUORUMS_ALT = ("0.5", "0.8")

# (id, prior that it holds, short text) -- copied from PREREG.md.
PRIORS = [
    ("S", 0.80, "Sanity: R defined for >= 60 bins AND R(Induan) < R(Changhsingian)"),
    ("D1+", 0.15, "Mayhew 2012: richness level tracks temperature level (rho > 0, p_shift < 0.05)"),
    ("D1-", 0.30, "Reverse: warm = fewer genera (rho < 0, lower p_shift < 0.05)"),
    ("D2+", 0.10, "Richness changes track temperature changes (rho > 0, p_perm < 0.05)"),
    ("D2-", 0.30, "Warming steps coincide with richness losses (rho < 0, lower p_perm < 0.05)"),
    ("D2b", 0.20, "D2 without the 6 Big Five boundaries is significant in either direction"),
    ("D3", 0.75, "D1's sign survives removing the rank effect of bin duration"),
    ("D4", 0.50, "D1 robust to PhanDA uncertainty (>= 90% of draws give D1's sign)"),
    ("D5a", 0.20, "Origination rises with temperature level (rho > 0, p_shift < 0.05)"),
    ("D5b", 0.25, "Origination rises after warming into the bin (rho > 0, p_perm < 0.05)"),
    ("D6", 0.70, "D1 and D2 keep their signs at quorums 0.5 and 0.8"),
    ("D7a", 0.50, "Raw richness correlates less positively with T than standardised richness"),
    ("D7b", 0.45, "Raw richness level is negatively related to T"),
]
LABELS = ("DATA_LIMITED", "SUPPORTED", "PARTLY_SUPPORTED", "REVERSED", "MIXED",
          "NOT_SUPPORTED")
READING = {
    "SUPPORTED": "On the best global temperature history, marine diversity is higher when "
                 "Earth is warmer: Mayhew et al. (2012) confirmed on independent temperature data.",
    "PARTLY_SUPPORTED": "One of the two registered positive signals (levels or changes) holds "
                        "and no negative one does.",
    "REVERSED": "On PhanDA, warm intervals and/or warming steps go with fewer genera: the "
                "opposite of the 2012 result, with sampling controlled.",
    "MIXED": "A positive and a negative signal both hold: the link depends on timescale "
             "(compare Spiridonov & Lovejoy 2022).",
    "NOT_SUPPORTED": "No registered signal holds: consistent with Balembois et al. (2025), "
                     "global temperature barely sets long-term global richness.",
    "DATA_LIMITED": "The sanity gate failed; nothing is interpreted.",
}

# PhanDA GMST_50 (deg C), read by hand from PhanDA_GMSTandCO2_percentiles.csv on
# 2026-11-17, rounded to 2 decimals (same table as extinction-temperature-test).
HAND_GMST50 = {
    "Tremadocian": 32.10, "Floian": 31.53, "Dapingian": 30.90, "Darriwilian": 31.00,
    "Sandbian": 27.31, "Katian": 27.85, "Hirnantian": 20.80, "Aeronian/Rhuddanian": 19.48,
    "Telychian": 24.36, "Sheinwoodian": 19.90, "Homerian": 20.97, "Gorstian": 25.60,
    "Ludfordian": 34.05, "Pridoli": 33.25, "Pragian/Lochkovian": 28.81, "Emsian": 25.29,
    "Eifelian": 22.49, "Givetian": 25.67, "Frasnian": 27.91, "Famennian": 23.76,
    "Tournaisian": 21.75, "Visean": 12.83, "Serpukhovian": 12.02, "Bashkirian": 14.08,
    "Moscovian": 18.09, "Kasimovian": 14.46, "Gzhelian": 16.55, "Asselian": 20.25,
    "Sakmarian": 20.80, "Artinskian": 20.76, "Kungurian": 18.70, "Roadian": 13.50,
    "Capitanian/Wordian": 18.68, "Wuchiapingian": 22.19, "Changhsingian": 18.02,
    "Induan": 28.50, "Anisian/Olenekian": 25.82, "Ladinian": 19.18, "Carnian": 23.43,
    "Norian": 25.02, "Rhaetian": 22.32, "Sinemurian/Hettangian": 25.06,
    "Pliensbachian": 24.10, "Toarcian": 26.37, "Aalenian": 23.20, "Bajocian": 21.39,
    "Bathonian": 17.77, "Callovian": 20.11, "Oxfordian": 23.88, "Kimmeridgian": 24.47,
    "Tithonian": 27.27, "Berriasian": 27.00, "Valanginian": 28.25, "Hauterivian": 27.27,
    "Barremian": 30.52, "Aptian": 32.00, "Albian": 26.01, "Cenomanian": 34.14,
    "Turonian": 35.56, "Coniacian": 32.26, "Santonian": 31.42, "Campanian": 28.13,
    "Maastrichtian": 24.51, "Selandian/Danian": 28.85, "Thanetian": 25.60,
    "Ypresian": 32.28, "Lutetian": 26.99, "Bartonian": 28.03, "Priabonian": 26.08,
    "Rupelian": 23.40, "Chattian": 23.44, "Aquitanian": 22.07, "Burdigalian": 24.06,
    "Langhian": 23.78, "Serravallian": 23.60, "Tortonian": 20.62, "Messinian": 16.94,
    "Zanclean": 16.27, "Piacenzian": 15.63,
}
HAND_TOL = 0.011


def _f(x):
    """float or NaN."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _finite(x):
    v = _f(x)
    return not (math.isnan(v) or math.isinf(v))


def _sign(x):
    v = _f(x)
    if not _finite(v) or v == 0:
        return 0
    return 1 if v > 0 else -1


def _sig(d, direction, p_key):
    """True if d's rho has `direction` (+1/-1) and the one-sided p in that
    direction is < ALPHA. p_key is 'p_perm' or 'p_shift'; the lower-tail key
    is p_key + '_lower'. False if anything is missing or NaN."""
    if not isinstance(d, dict):
        return False
    rho = _f(d.get("rho"))
    p = _f(d.get(p_key if direction > 0 else p_key + "_lower"))
    return bool(_finite(rho) and _finite(p) and _sign(rho) == direction and p < ALPHA)


def grade_rows(a):
    """{id: True/False/None} from results.json's "analysis" block. None means
    the quantity was not computed (analysis stopped early)."""
    g = {k: None for k, _, _ in PRIORS}
    if not isinstance(a, dict):
        return g
    san = a.get("sanity")
    if isinstance(san, dict) and "passed" in san:
        g["S"] = bool(san["passed"])
    if "D1_levels" not in a:      # analyze() stopped at MIN_RUN
        return g
    d1, d2 = a.get("D1_levels") or {}, a.get("D2_differences") or {}
    g["D1+"] = _sig(d1, +1, "p_shift")
    g["D1-"] = _sig(d1, -1, "p_shift")
    g["D2+"] = _sig(d2, +1, "p_perm")
    g["D2-"] = _sig(d2, -1, "p_perm")
    d2b = a.get("D2b_differences_without_big_five") or {}
    g["D2b"] = bool(_sig(d2b, +1, "p_perm") or _sig(d2b, -1, "p_perm"))
    s1 = _sign(d1.get("rho"))
    s2 = _sign(d2.get("rho"))
    s3 = _sign((a.get("D3_duration_partial") or {}).get("rho"))
    g["D3"] = bool(s1 != 0 and s3 == s1)
    d4 = a.get("D4_uncertainty") or {}
    fr = _f(d4.get("frac_same_sign_as_D1"))
    g["D4"] = bool(s1 != 0 and _finite(fr) and fr >= 0.9)
    g["D5a"] = _sig(a.get("D5a_origination_vs_level"), +1, "p_shift")
    g["D5b"] = _sig(a.get("D5b_origination_vs_change_in"), +1, "p_perm")
    d6 = a.get("D6_quorum") or {}
    ok6 = s1 != 0 and s2 != 0
    for q in QUORUMS_ALT:
        blk = d6.get("R_" + q) or {}
        ok6 = ok6 and _sign((blk.get("levels") or {}).get("rho")) == s1 \
            and _sign((blk.get("differences") or {}).get("rho")) == s2
    g["D6"] = bool(ok6)
    raw = _f((a.get("D7_raw_levels") or {}).get("rho"))
    r1 = _f(d1.get("rho"))
    g["D7a"] = bool(_finite(raw) and _finite(r1) and raw < r1)
    g["D7b"] = bool(_finite(raw) and raw < 0)
    return g


def label_from(g):
    """PREREG decision rule from the graded rows."""
    if not g.get("S"):
        return "DATA_LIMITED"
    pos = bool(g.get("D1+")) or bool(g.get("D2+"))
    neg = bool(g.get("D1-")) or bool(g.get("D2-"))
    if g.get("D1+") and g.get("D2+"):
        return "SUPPORTED"
    if pos and not neg:
        return "PARTLY_SUPPORTED"
    if neg and not pos:
        return "REVERSED"
    if pos and neg:
        return "MIXED"
    return "NOT_SUPPORTED"


def brier(g, priors=PRIORS):
    """Mean squared error of the registered priors over graded rows, and the
    same for a know-nothing 0.5 forecaster."""
    sq, n = 0.0, 0
    for k, p, _ in priors:
        if g.get(k) is None:
            continue
        o = 1.0 if g[k] else 0.0
        sq += (p - o) ** 2
        n += 1
    if n == 0:
        return {"n": 0, "brier": float("nan"), "brier_coinflip": float("nan")}
    return {"n": n, "brier": sq / n, "brier_coinflip": 0.25}


def check_bins(rows):
    """Compare bins.csv rows (list of dicts) with HAND_GMST50."""
    seen, bad = {}, []
    for r in rows:
        name = (r.get("bin") or "").strip()
        val = (r.get("GMST_50") or "").strip()
        if not name:
            continue
        if val == "" and name not in HAND_GMST50:
            continue            # Cambrian and Quaternary bins carry no temperature
        seen[name] = _f(val)
    for name, v in sorted(seen.items()):
        if name not in HAND_GMST50:
            bad.append({"bin": name, "job": v, "hand": None})
        elif not (_finite(v) and abs(v - HAND_GMST50[name]) <= HAND_TOL):
            bad.append({"bin": name, "job": v, "hand": HAND_GMST50[name]})
    missing = sorted(set(HAND_GMST50) - set(seen))
    return {"n_checked": len(seen), "mismatches": bad, "missing_bins": missing,
            "passed": not bad and not missing}


def grade(results, bin_rows=None):
    """Everything graded.json holds."""
    out = {"source": "grade_diversity.py (published before results)",
           "job_status": results.get("status") if isinstance(results, dict) else None,
           "label": None, "rows": []}
    if not isinstance(results, dict) or results.get("status") != "ok" \
            or "analysis" not in results:
        out["status"] = "JOB_FAILED"
        return out
    a = results["analysis"]
    g = grade_rows(a)
    out["label"] = label_from(g)
    out["reading"] = READING[out["label"]]
    out["label_from_job"] = results.get("decision", a.get("decision"))
    out["label_agrees_with_job"] = out["label"] == out["label_from_job"]
    sig = a.get("signals") or {}
    mine = {"D1_pos": g["D1+"], "D1_neg": g["D1-"], "D2_pos": g["D2+"], "D2_neg": g["D2-"]}
    out["signals_disagreeing_with_job"] = sorted(
        k for k, v in mine.items() if k in sig and v is not None and bool(sig[k]) != v)
    out["qualifiers_from_job"] = a.get("qualifiers")
    for k, p, text in PRIORS:
        out["rows"].append({"id": k, "prior": p, "prediction": text, "holds": g[k]})
    out["brier"] = brier(g)
    status = "GRADED"
    if not g.get("S"):
        status = "GRADED_SANITY_FAILED"
    if not out["label_agrees_with_job"] or out["signals_disagreeing_with_job"]:
        status = "CHECK_CODE"
    if bin_rows is not None:
        out["phanda_check"] = check_bins(bin_rows)
        if not out["phanda_check"]["passed"]:
            status = "CHECK_DATA"
    out["status"] = status
    out["n_levels"] = a.get("n_levels")
    return out


def _fmt(x, nd=3):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    v = _f(x)
    if not _finite(v):
        return "n/a"
    return f"{v:.{nd}f}"


def _test_line(name, d):
    d = d if isinstance(d, dict) else {}
    parts = [f"rho {_fmt(d.get('rho'))}", f"n {d.get('n', 'n/a')}"]
    for k in ("p_perm", "p_perm_lower", "p_shift", "p_shift_lower"):
        if k in d:
            parts.append(f"{k} {_fmt(d.get(k), 4)}")
    return f"- {name}: " + ", ".join(parts)


def render_md(gr, results):
    a = results.get("analysis", {}) if isinstance(results, dict) else {}
    L = ["# Graded result: diversity-temperature-test", "",
         "Produced mechanically by `grade_diversity.py`, which was published before the",
         "result existed. Made with AI agents working in MV2's name.", "",
         f"- Status: **{gr.get('status')}**",
         f"- Decision label (PREREG rule): **{gr.get('label')}**"]
    if gr.get("status") == "JOB_FAILED":
        L += ["", f"The job did not finish its analysis (status `{gr.get('job_status')}`).",
              "Nothing is graded."]
        return "\n".join(L) + "\n"
    L += [f"- Registered reading: {gr.get('reading')}",
          f"- Label written by the job: {gr.get('label_from_job')} "
          f"(agrees: {_fmt(gr.get('label_agrees_with_job'))})",
          f"- Signals disagreeing with the job: "
          f"{', '.join(gr.get('signals_disagreeing_with_job') or []) or 'none'}",
          f"- Bins with a standardised richness: {gr.get('n_levels')}", ""]
    if gr.get("status") == "GRADED_SANITY_FAILED":
        L += ["The sanity gate failed, so the label is DATA_LIMITED and the rows below are",
              "reported but not interpreted.", ""]
    L += ["| # | Prediction | Prior | Holds? |", "|---|---|---|---|"]
    for r in gr["rows"]:
        h = "n/a" if r["holds"] is None else ("YES" if r["holds"] else "no")
        L.append(f"| {r['id']} | {r['prediction']} | {int(round(100 * r['prior']))}% | {h} |")
    b = gr["brier"]
    L += ["", f"Brier score of the registered priors: {_fmt(b['brier'])} over {b['n']} rows "
          f"(a 50% forecaster scores {_fmt(b['brier_coinflip'])}).", "", "## Key numbers", ""]
    for key, name in (("D1_levels", "D1 levels (detrended ln R vs T)"),
                      ("D2_differences", "D2 differences (d ln R vs dT)"),
                      ("D2b_differences_without_big_five", "D2b without Big Five"),
                      ("D3_duration_partial", "D3 duration partial"),
                      ("D5a_origination_vs_level", "D5a origination vs T"),
                      ("D5b_origination_vs_change_in", "D5b origination vs dT into bin"),
                      ("D7_raw_levels", "D7 raw richness levels"),
                      ("D7_collections_levels", "D7 collections (effort) levels")):
        L.append(_test_line(name, a.get(key)))
    d4 = a.get("D4_uncertainty") or {}
    L.append(f"- D4: fraction of PhanDA draws with D1's sign: "
             f"{_fmt(d4.get('frac_same_sign_as_D1'))} (rho 5-95%: {_fmt(d4.get('rho_p05'))} "
             f"to {_fmt(d4.get('rho_p95'))})")
    d6 = a.get("D6_quorum") or {}
    for q in QUORUMS_ALT:
        blk = d6.get("R_" + q) or {}
        L.append(f"- D6 quorum {q}: levels rho {_fmt((blk.get('levels') or {}).get('rho'))}, "
                 f"differences rho {_fmt((blk.get('differences') or {}).get('rho'))}")
    L += ["", "## Qualifiers written by the job (never change the label)", ""]
    for k, v in (gr.get("qualifiers_from_job") or {}).items():
        L.append(f"- {k}: {_fmt(v)}")
    if "phanda_check" in gr:
        pc = gr["phanda_check"]
        L += ["", "## PhanDA hand check", "",
              f"{pc['n_checked']} bins compared with the hand-read GMST_50 values; "
              f"mismatches: {len(pc['mismatches'])}; missing: {len(pc['missing_bins'])}."]
        for m in pc["mismatches"]:
            L.append(f"- {m['bin']}: job {_fmt(m['job'], 2)}, hand {_fmt(m['hand'], 2)}")
        for m in pc["missing_bins"]:
            L.append(f"- {m}: missing from bins.csv")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--bins", default=None)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args(argv)
    with open(args.results, encoding="utf-8") as f:
        results = json.load(f)
    rows = None
    if args.bins and os.path.exists(args.bins):
        with open(args.bins, newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
    gr = grade(results, rows)
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "graded.json"), "w", encoding="utf-8") as f:
        json.dump(gr, f, indent=1, default=str)
    with open(os.path.join(args.out_dir, "GRADED.md"), "w", encoding="utf-8") as f:
        f.write(render_md(gr, results))
    print(json.dumps({"status": gr["status"], "label": gr["label"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
