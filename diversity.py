# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Pre-registered test (PREREG.md): does marine biodiversity track global
temperature over the Phanerozoic? Mayhew, Bell, Benton & McGowan (2012, PNAS
109:15141, "Biodiversity tracks temperature over time") found that, once
sampling is controlled, marine invertebrate richness rises with temperature.
This module tests that on PBDB genera standardised by coverage-based
rarefaction (Chao & Jost 2012) against PhanDA global mean surface temperature
(Judd et al. 2024, Science 385:eadk3705).

Pure functions only (no network); job_diversity.py does the downloads.
Made with AI agents working in MV2's name.
"""
import math
from collections import Counter

import numpy as np

import thermo

QUORUM = 0.7              # primary coverage quorum
QUORUMS_ALT = (0.5, 0.8)  # robustness quorums
MIN_BINS = 60             # sanity: bins with a standardised richness
MIN_RUN = 12              # below this a test is not computed
ALPHA = 0.05
NAN = float("nan")


# ------------------------------------------------------------- counting
class CountAccumulator:
    """Per-stage genus frequencies. A genus counts once per collection, so its
    frequency is the number of collections (occurrence lists) it appears in."""

    def __init__(self, n):
        self.counts = [Counter() for _ in range(n)]
        self._pairs = [set() for _ in range(n)]
        self._colls = [set() for _ in range(n)]
        self.stats = {"rows": 0, "assigned": 0, "no_genus": 0, "no_bin": 0,
                      "same_genus_same_collection": 0}

    def add_rows(self, rows, assigner):
        for row in rows:
            self.stats["rows"] += 1
            key = thermo.genus_key(row)
            if key is None:
                self.stats["no_genus"] += 1
                continue
            b = assigner(thermo.to_float(row.get("max_ma")), thermo.to_float(row.get("min_ma")))
            if b is None:
                self.stats["no_bin"] += 1
                continue
            self.stats["assigned"] += 1
            coll = (row.get("collection_no") or "").strip()
            if coll:
                self._colls[b].add(coll)
                pair = (key, coll)
                if pair in self._pairs[b]:
                    self.stats["same_genus_same_collection"] += 1
                    continue
                self._pairs[b].add(pair)
            self.counts[b][key] += 1

    @property
    def n_coll(self):
        return [len(c) for c in self._colls]


def merge_counts(stage_counts, bins, stage_ncoll):
    """Sum stage-level genus frequencies (and collection counts) into bins."""
    idx = {s[0]: i for i, s in enumerate(thermo.STAGES)}
    counts, ncoll = [], []
    for b in bins:
        c = Counter()
        for s in b["stages"]:
            c.update(stage_counts[idx[s]])
        counts.append(c)
        ncoll.append(sum(stage_ncoll[idx[s]] for s in b["stages"]))
    return counts, ncoll


# ------------------------------------------- coverage-based rarefaction
def freq_counts(counts):
    """{k: number of genera seen exactly k times}."""
    return Counter(int(v) for v in counts if v > 0)


def _lcr(a, b, m):
    """log(C(a, m) / C(b, m)) for a, b >= m >= 0."""
    return (math.lgamma(a + 1) - math.lgamma(a - m + 1)
            - math.lgamma(b + 1) + math.lgamma(b - m + 1))


def sample_coverage(fk, n):
    """Chao & Jost (2012) estimate of the coverage of the full sample."""
    if n <= 0:
        return NAN
    f1, f2 = fk.get(1, 0), fk.get(2, 0)
    if f1 == 0:
        return 1.0
    if f2 > 0:
        a = (n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
    else:
        a = (n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2)
    return 1.0 - f1 / n * a


def rarefied(fk, n, m):
    """(expected coverage, expected richness) of a random subsample of m of the
    n occurrences, without replacement (Chao & Jost 2012, abundance form).
    Coverage is valid for 1 <= m <= n - 1."""
    loss, s = 0.0, 0.0
    for k, f in fk.items():
        if n - k >= m:
            r_s = math.exp(_lcr(n - k, n, m))        # C(n-k, m) / C(n, m)
            r_c = math.exp(_lcr(n - k, n - 1, m))    # C(n-k, m) / C(n-1, m)
        else:
            r_s = r_c = 0.0
        s += f * (1.0 - r_s)
        loss += f * (k / n) * r_c
    return 1.0 - loss, s


def coverage_richness(counts, quorum):
    """Coverage-standardised richness: expected richness of the smallest
    subsample whose expected coverage reaches `quorum`. Interpolation only:
    NaN if the bin's estimated coverage is below the quorum."""
    counts = [int(c) for c in counts if c > 0]
    n = sum(counts)
    fk = freq_counts(counts)
    out = {"n": n, "s_obs": len(counts), "coverage": sample_coverage(fk, n),
           "m": None, "richness": NAN}
    if n < 2 or not (out["coverage"] >= quorum):
        return out
    if rarefied(fk, n, n - 1)[0] < quorum:
        out["m"], out["richness"] = n, float(len(counts))
        return out
    lo, hi = 1, n - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if rarefied(fk, n, mid)[0] >= quorum:
            hi = mid
        else:
            lo = mid + 1
    out["m"] = lo
    out["richness"] = rarefied(fk, n, lo)[1]
    return out


# --------------------------------------------------------------- rates
def three_timer_origination(sets):
    """Alroy (2008) three-timer origination rate at the BASE of bin i:
    lambda_i = ln(2t_{i+1} / 3t_i) + ln(Ps_{i-1}), with 2t_{i+1} = |S_i & S_{i+1}|,
    3t_i = |S_{i-1} & S_i & S_{i+1}| and Ps_j = 3t_j / (3t_j + pt_j); NaN where
    undefined. Mirror image of thermo.three_timer."""
    n = len(sets)
    ps = [NAN] * n
    three = [0] * n
    for j in range(1, n - 1):
        a, b, c = sets[j - 1], sets[j], sets[j + 1]
        three[j] = len(a & b & c)
        part = len((a & c) - b)
        if three[j] + part > 0:
            ps[j] = three[j] / (three[j] + part)
    lam = [NAN] * n
    for i in range(2, n - 1):
        two = len(sets[i] & sets[i + 1])
        p = ps[i - 1]
        if two > 0 and three[i] > 0 and p == p and p > 0:
            lam[i] = math.log(two / three[i]) + math.log(p)
    return lam


# ---------------------------------------------------------------- table
def level_rows(bins, rich, raw, n_coll, lam, temps):
    """One row per PhanDA bin. rich = {quorum: [richness per bin]}."""
    rows = []
    for j, b in enumerate(bins):
        if b["phanda"] is None or b["phanda"] not in temps:
            continue
        t = temps[b["phanda"]]
        prev = bins[j - 1] if j > 0 else None
        if prev is not None and prev["phanda"] is not None and prev["phanda"] in temps:
            dt_in = t[1] - temps[prev["phanda"]][1]
        else:
            dt_in = NAN
        r = {"j": j, "bin": b["name"], "mid_ma": 0.5 * (b["base"] + b["top"]),
             "dur": b["base"] - b["top"], "T16": t[0], "T": t[1], "T84": t[2],
             "dT_in": dt_in, "R": rich[QUORUM][j], "raw": float(raw[j]),
             "n_coll": float(n_coll[j]), "lam": lam[j]}
        for q in QUORUMS_ALT:
            r[f"R_{q}"] = rich[q][j]
        rows.append(r)
    return rows


def _pos(v):
    return v is not None and v == v and v > 0


def _arr(rows, key):
    return np.array([r[key] for r in rows], float)


def detrend(y, x):
    """Residuals of y after a least-squares straight line in x."""
    y, x = np.asarray(y, float), np.asarray(x, float)
    A = np.column_stack([np.ones(len(x)), x])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ coef


def diff_rows(rows, key="R"):
    """First differences between adjacent bins: the boundary at the top of
    bin a; dlogR = ln R(next) - ln R(a), dT = T(next) - T(a)."""
    out = []
    for a, b in zip(rows[:-1], rows[1:]):
        if b["j"] != a["j"] + 1 or not (_pos(a[key]) and _pos(b[key])):
            continue
        out.append({"bin": a["bin"], "dlogR": math.log(b[key]) - math.log(a[key]),
                    "dT": b["T"] - a["T"], "big_five": a["bin"] in thermo.BIG_FIVE_TOPS})
    return out


# ----------------------------------------------------------- statistics
def _empty(n):
    return {"n": int(n), "rho": NAN, "p_perm": 1.0, "p_perm_lower": 1.0,
            "p_shift": 1.0, "p_shift_lower": 1.0}


def signed_test(x, y, n_perm, rng, shift=True):
    """Spearman test with one-sided p values in both directions."""
    x = np.asarray(x, float)
    if len(x) < MIN_RUN:
        return _empty(len(x))
    up = thermo.spearman_test(x, y, n_perm, rng, shift)
    lo = thermo.spearman_test(-x, y, n_perm, rng, shift)
    out = dict(up)
    out["p_perm_lower"] = lo["p_perm"]
    if shift:
        out["p_shift_lower"] = lo["p_shift"]
    return out


def level_test(rows, key, n_perm, rng):
    """Detrended ln(key) against detrended T (both linear in bin mid-age)."""
    lv = [r for r in rows if _pos(r[key])]
    if len(lv) < MIN_RUN:
        return _empty(len(lv))
    age = _arr(lv, "mid_ma")
    y = detrend(np.log(_arr(lv, key)), age)
    x = detrend(_arr(lv, "T"), age)
    return signed_test(x, y, n_perm, rng)


def diff_test(rows, key, n_perm, rng, drop_big_five=False):
    d = diff_rows(rows, key)
    if drop_big_five:
        d = [r for r in d if not r["big_five"]]
    if len(d) < MIN_RUN:
        return _empty(len(d))
    return signed_test(_arr(d, "dT"), _arr(d, "dlogR"), n_perm, rng)


def _sig_pos(t, key):
    return bool(t["rho"] == t["rho"] and t["rho"] > 0 and t.get(key, 1.0) < ALPHA)


def _sig_neg(t, key):
    return bool(t["rho"] == t["rho"] and t["rho"] < 0 and t.get(key + "_lower", 1.0) < ALPHA)


def _sign(v):
    return 0 if v != v or v == 0 else (1 if v > 0 else -1)


def analyze(rows, n_perm, rng, n_draws=2000):
    """All registered tests (PREREG.md) on the per-bin table."""
    lv = [r for r in rows if _pos(r["R"])]
    byname = {r["bin"]: r for r in lv}
    ind, chg = byname.get("Induan"), byname.get("Changhsingian")
    pt_drop = bool(ind is not None and chg is not None and ind["R"] < chg["R"])
    res = {"n_levels": len(lv), "bins_with_richness": [r["bin"] for r in lv]}
    res["sanity"] = {"n_ok": len(lv) >= MIN_BINS, "pt_drop": pt_drop,
                     "passed": bool(len(lv) >= MIN_BINS and pt_drop)}
    if len(lv) < MIN_RUN:
        res["decision"] = "DATA_LIMITED"
        return res

    d1 = level_test(rows, "R", n_perm, rng)
    d2 = diff_test(rows, "R", n_perm, rng)
    res["D1_levels"] = d1
    res["D2_differences"] = d2
    res["D2b_differences_without_big_five"] = diff_test(rows, "R", n_perm, rng,
                                                        drop_big_five=True)

    age = _arr(lv, "mid_ma")
    y = detrend(np.log(_arr(lv, "R")), age)
    x = detrend(_arr(lv, "T"), age)
    res["D3_duration_partial"] = thermo.partial_spearman_test(
        x, y, np.log(_arr(lv, "dur")), n_perm, rng)

    # D4: PhanDA uncertainty (independent two-piece normal draws per slice).
    z = rng.standard_normal((n_draws, len(lv)))
    Td = thermo.split_normal(_arr(lv, "T16"), _arr(lv, "T"), _arr(lv, "T84"), z)
    A = np.column_stack([np.ones(len(lv)), age])
    coef, *_ = np.linalg.lstsq(A, Td.T, rcond=None)
    resid = Td.T - A @ coef
    ry = thermo.rankdata(y)
    rhos = np.array([thermo.corr(thermo.rankdata(resid[:, k]), ry) for k in range(n_draws)])
    sgn = 1.0 if d1["rho"] >= 0 else -1.0
    res["D4_uncertainty"] = {"n_draws": int(n_draws),
                             "rho_p05": float(np.quantile(rhos, 0.05)),
                             "rho_median": float(np.median(rhos)),
                             "rho_p95": float(np.quantile(rhos, 0.95)),
                             "frac_rho_pos": float(np.mean(rhos > 0)),
                             "frac_same_sign_as_D1": float(np.mean(sgn * rhos > 0))}

    # D5: origination against temperature level and against change into the bin.
    lr = [r for r in rows if np.isfinite(r["lam"])]
    res["D5a_origination_vs_level"] = signed_test(_arr(lr, "T"), _arr(lr, "lam"), n_perm, rng)
    lr2 = [r for r in lr if np.isfinite(r["dT_in"])]
    res["D5b_origination_vs_change_in"] = signed_test(_arr(lr2, "dT_in"), _arr(lr2, "lam"),
                                                      n_perm, rng)

    # D6: other quorums.
    d6 = {}
    for q in QUORUMS_ALT:
        k = f"R_{q}"
        d6[k] = {"levels": level_test(rows, k, n_perm, rng),
                 "differences": diff_test(rows, k, n_perm, rng)}
    d6["same_signs"] = bool(all(
        _sign(d6[f"R_{q}"]["levels"]["rho"]) == _sign(d1["rho"]) != 0
        and _sign(d6[f"R_{q}"]["differences"]["rho"]) == _sign(d2["rho"]) != 0
        for q in QUORUMS_ALT))
    res["D6_quorum"] = d6

    # D7: raw (unstandardised) richness and sampling effort against temperature.
    raw = level_test(rows, "raw", n_perm, rng)
    res["D7_raw_levels"] = raw
    res["D7_raw_below_standardised"] = bool(raw["rho"] == raw["rho"] and d1["rho"] == d1["rho"]
                                            and raw["rho"] < d1["rho"])
    res["D7_collections_levels"] = level_test(rows, "n_coll", n_perm, rng)

    pos1, neg1 = _sig_pos(d1, "p_shift"), _sig_neg(d1, "p_shift")
    pos2, neg2 = _sig_pos(d2, "p_perm"), _sig_neg(d2, "p_perm")
    res["signals"] = {"D1_pos": pos1, "D1_neg": neg1, "D2_pos": pos2, "D2_neg": neg2}
    if not res["sanity"]["passed"]:
        res["decision"] = "DATA_LIMITED"
    elif pos1 and pos2:
        res["decision"] = "SUPPORTED"
    elif (pos1 or pos2) and not (neg1 or neg2):
        res["decision"] = "PARTLY_SUPPORTED"
    elif (neg1 or neg2) and not (pos1 or pos2):
        res["decision"] = "REVERSED"
    elif (pos1 or pos2) and (neg1 or neg2):
        res["decision"] = "MIXED"
    else:
        res["decision"] = "NOT_SUPPORTED"
    res["qualifiers"] = {
        "D2_autocorrelation_robust": bool(_sig_pos(d2, "p_shift") or _sig_neg(d2, "p_shift")),
        "D2_without_big_five_same_sign": bool(
            _sign(res["D2b_differences_without_big_five"]["rho"]) == _sign(d2["rho"]) != 0),
        "D1_duration_same_sign": bool(_sign(res["D3_duration_partial"]["rho"]) == _sign(d1["rho"])),
        "D1_uncertainty_robust": bool(res["D4_uncertainty"]["frac_same_sign_as_D1"] >= 0.9),
        "quorum_robust": d6["same_signs"]}
    return res
