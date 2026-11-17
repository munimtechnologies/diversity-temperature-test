# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Pre-registered test (PREREG.md) of Song et al. (2021, Nature Communications
12:4694, "Thresholds of temperature change for mass extinctions") on an
independent temperature history: PhanDA (Judd et al. 2024, Science 385:eadk3705).

Question: does the size of global temperature change between consecutive
stages predict marine genus extinction at the boundary between them, beyond
the Big Five mass extinctions, and do only the Big Five exceed 5.2 C?

Pure functions only (no network); job_thermo.py does the downloads.
Made with AI agents working in MV2's name.

(Copied from extinction-temperature-test v0.1.0 without changes to any function
kept, so that diversity-temperature-test uses the same stage table, bins, PhanDA
parser and statistics; the extinction-only analysis functions were left out.
diversity.py imports it.)
"""
import bisect
import csv
import io
import math

import numpy as np

# International stages, GTS2020 ages as served by the Paleobiology Database
# (same table as extinction-cycle-test/cycle.py). (name, base_ma, top_ma),
# oldest first. All Quaternary ages are one bin.
STAGES = [
    ("Fortunian", 538.8, 529.0), ("Stage 2", 529.0, 521.0),
    ("Stage 3", 521.0, 514.5), ("Stage 4", 514.5, 506.5),
    ("Wuliuan", 506.5, 504.5), ("Drumian", 504.5, 500.5),
    ("Guzhangian", 500.5, 497.0), ("Paibian", 497.0, 494.2),
    ("Jiangshanian", 494.2, 491.0), ("Stage 10", 491.0, 486.85),
    ("Tremadocian", 486.85, 477.1), ("Floian", 477.1, 471.3),
    ("Dapingian", 471.3, 469.4), ("Darriwilian", 469.4, 458.2),
    ("Sandbian", 458.2, 452.8), ("Katian", 452.8, 445.2),
    ("Hirnantian", 445.2, 443.1), ("Rhuddanian", 443.1, 440.5),
    ("Aeronian", 440.5, 438.6), ("Telychian", 438.6, 432.9),
    ("Sheinwoodian", 432.9, 430.6), ("Homerian", 430.6, 426.7),
    ("Gorstian", 426.7, 425.0), ("Ludfordian", 425.0, 422.7),
    ("Pridoli", 422.7, 419.62), ("Lochkovian", 419.62, 413.02),
    ("Pragian", 413.02, 410.62), ("Emsian", 410.62, 393.47),
    ("Eifelian", 393.47, 387.95), ("Givetian", 387.95, 382.31),
    ("Frasnian", 382.31, 372.15), ("Famennian", 372.15, 358.86),
    ("Tournaisian", 358.86, 346.7), ("Visean", 346.7, 330.3),
    ("Serpukhovian", 330.3, 323.4), ("Bashkirian", 323.4, 315.2),
    ("Moscovian", 315.2, 307.0), ("Kasimovian", 307.0, 303.7),
    ("Gzhelian", 303.7, 298.9), ("Asselian", 298.9, 293.52),
    ("Sakmarian", 293.52, 290.1), ("Artinskian", 290.1, 283.3),
    ("Kungurian", 283.3, 274.4), ("Roadian", 274.4, 266.9),
    ("Wordian", 266.9, 264.28), ("Capitanian", 264.28, 259.857),
    ("Wuchiapingian", 259.857, 254.14), ("Changhsingian", 254.14, 251.902),
    ("Induan", 251.902, 250.8), ("Olenekian", 250.8, 247.0),
    ("Anisian", 247.0, 241.464), ("Ladinian", 241.464, 237.0),
    ("Carnian", 237.0, 227.3), ("Norian", 227.3, 205.7),
    ("Rhaetian", 205.7, 201.4), ("Hettangian", 201.4, 199.5),
    ("Sinemurian", 199.5, 192.9), ("Pliensbachian", 192.9, 184.2),
    ("Toarcian", 184.2, 174.7), ("Aalenian", 174.7, 170.9),
    ("Bajocian", 170.9, 168.2), ("Bathonian", 168.2, 165.3),
    ("Callovian", 165.3, 161.5), ("Oxfordian", 161.5, 154.8),
    ("Kimmeridgian", 154.8, 149.2), ("Tithonian", 149.2, 143.1),
    ("Berriasian", 143.1, 137.05), ("Valanginian", 137.05, 132.6),
    ("Hauterivian", 132.6, 125.77), ("Barremian", 125.77, 121.4),
    ("Aptian", 121.4, 113.2), ("Albian", 113.2, 100.5),
    ("Cenomanian", 100.5, 93.9), ("Turonian", 93.9, 89.8),
    ("Coniacian", 89.8, 85.7), ("Santonian", 85.7, 83.6),
    ("Campanian", 83.6, 72.2), ("Maastrichtian", 72.2, 66.0),
    ("Danian", 66.0, 61.66), ("Selandian", 61.66, 59.24),
    ("Thanetian", 59.24, 56.0), ("Ypresian", 56.0, 48.07),
    ("Lutetian", 48.07, 41.03), ("Bartonian", 41.03, 37.71),
    ("Priabonian", 37.71, 33.9), ("Rupelian", 33.9, 27.3),
    ("Chattian", 27.3, 23.04), ("Aquitanian", 23.04, 20.45),
    ("Burdigalian", 20.45, 15.98), ("Langhian", 15.98, 13.82),
    ("Serravallian", 13.82, 11.63), ("Tortonian", 11.63, 7.246),
    ("Messinian", 7.246, 5.333), ("Zanclean", 5.333, 3.6),
    ("Piacenzian", 3.6, 2.58), ("Quaternary", 2.58, 0.0),
]

# All 85 PhanDA time slices (Stage column of PhanDA_GMSTandCO2_percentiles.csv),
# youngest first, as listed in the file on 2026-11-17.
PHANDA_SLICES = (
    "Holocene", "Upper Pleistocene", "Chibanian", "Calabrian", "Gelasian",
    "Piacenzian", "Zanclean", "Messinian", "Tortonian", "Serravallian", "Langhian",
    "Burdigalian", "Aquitanian", "Chattian", "Rupelian", "Priabonian", "Bartonian",
    "Lutetian", "Ypresian", "PETM", "Thanetian", "Selandian/Danian", "Maastrichtian",
    "Campanian", "Santonian", "Coniacian", "Turonian", "Cenomanian", "Albian", "Aptian",
    "Barremian", "Hauterivian", "Valanginian", "Berriasian", "Tithonian", "Kimmeridgian",
    "Oxfordian", "Callovian", "Bathonian", "Bajocian", "Aalenian", "Toarcian",
    "Pliensbachian", "Sinemurian/Hettangian", "Rhaetian", "Norian", "Carnian",
    "Ladinian", "Anisian/Olenekian", "Induan", "Changhsingian", "Wuchiapingian",
    "Capitanian/Wordian", "Roadian", "Kungurian", "Artinskian", "Sakmarian", "Asselian",
    "Gzhelian", "Kasimovian", "Moscovian", "Bashkirian", "Serpukhovian", "Visean",
    "Tournaisian", "Famennian", "Frasnian", "Givetian", "Eifelian", "Emsian",
    "Pragian/Lochkovian", "Pridoli", "Ludfordian", "Gorstian", "Homerian",
    "Sheinwoodian", "Telychian", "Aeronian/Rhuddanian", "Hirnantian", "Katian",
    "Sandbian", "Darriwilian", "Dapingian", "Floian", "Tremadocian",
)
# Slices merged by PhanDA, mapped to the PBDB stages they contain (oldest first).
MERGED = {
    "Aeronian/Rhuddanian": ("Rhuddanian", "Aeronian"),
    "Pragian/Lochkovian": ("Lochkovian", "Pragian"),
    "Capitanian/Wordian": ("Wordian", "Capitanian"),
    "Anisian/Olenekian": ("Olenekian", "Anisian"),
    "Sinemurian/Hettangian": ("Hettangian", "Sinemurian"),
    "Selandian/Danian": ("Danian", "Selandian"),
}
# Not used as temperature bins: the five Quaternary slices (PBDB is binned as
# one Quaternary bin) and the 0.3-Myr PETM event slice (registered choice).
DROPPED_SLICES = ("Holocene", "Upper Pleistocene", "Chibanian", "Calabrian",
                  "Gelasian", "PETM")

# Boundaries (bin TOPS) of the Big Five; the end-Ordovician has two pulses.
BIG_FIVE_TOPS = ("Katian", "Hirnantian", "Frasnian", "Changhsingian", "Rhaetian",
                 "Maastrichtian")
BIG_FIVE_EVENTS = (("Katian", "Hirnantian"), ("Frasnian",), ("Changhsingian",),
                   ("Rhaetian",), ("Maastrichtian",))
THRESHOLD_C = 5.2      # Song et al. 2021 magnitude threshold
MIN_N = 60             # sanity: boundaries with a defined rate
MIN_RUN = 12           # below this nothing is computed
MIN_GROUP = 10         # warming / cooling subsets


# ------------------------------------------------------------------ bins
def build_bins():
    """Composite bins, oldest first: the 10 Cambrian stages (no temperature),
    79 PhanDA temperature bins (Tremadocian..Piacenzian), then Quaternary
    (no temperature). Each bin: name, stages, phanda, base, top."""
    idx = {s[0]: i for i, s in enumerate(STAGES)}
    used = [p for p in reversed(PHANDA_SLICES) if p not in DROPPED_SLICES]
    bins = []
    for s in STAGES[:10]:
        bins.append({"name": s[0], "stages": (s[0],), "phanda": None})
    for p in used:
        bins.append({"name": p, "stages": MERGED.get(p, (p,)), "phanda": p})
    bins.append({"name": "Quaternary", "stages": ("Quaternary",), "phanda": None})
    for b in bins:
        b["base"] = max(STAGES[idx[s]][1] for s in b["stages"])
        b["top"] = min(STAGES[idx[s]][2] for s in b["stages"])
    return bins


class BinAssigner:
    """Assign an occurrence age range to one stage: the stage containing the
    range midpoint; the range must lie inside it within `tol` Myr, else None."""

    def __init__(self, stages=STAGES, tol=0.1):
        self.stages = stages
        self.n = len(stages)
        self.tol = tol
        self.bases = [s[1] for s in stages]
        self.tops = [s[2] for s in stages]
        self._tops_asc = self.tops[::-1]

    def __call__(self, max_ma, min_ma):
        if max_ma is None or min_ma is None:
            return None
        if max_ma < min_ma:
            max_ma, min_ma = min_ma, max_ma
        m = 0.5 * (max_ma + min_ma)
        if m < 0 or m >= self.bases[0]:
            return None
        k = bisect.bisect_right(self._tops_asc, m) - 1
        if k < 0:
            return None
        i = self.n - 1 - k
        if max_ma <= self.bases[i] + self.tol and min_ma >= self.tops[i] - self.tol:
            return i
        return None


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def genus_key(row):
    """'phylum|Genus' for a PBDB occurrence row (pbdb vocabulary), or None."""
    g = (row.get("genus") or "").strip()
    if not g or g.upper().startswith("NO_"):
        rank = (row.get("accepted_rank") or "").strip()
        name = (row.get("accepted_name") or "").strip()
        g = name.split()[0] if (name and rank in (
            "genus", "subgenus", "species", "subspecies")) else ""
    if not g:
        return None
    g = g.split()[0]
    ph = (row.get("phylum") or "").strip()
    return ph + "|" + g


class Accumulator:
    def __init__(self, n):
        self.sets = [set() for _ in range(n)]
        self.n_occ = [0] * n
        self._colls = [set() for _ in range(n)]
        self.stats = {"rows": 0, "assigned": 0, "no_genus": 0, "no_bin": 0}

    def add_rows(self, rows, assigner):
        for row in rows:
            self.stats["rows"] += 1
            key = genus_key(row)
            if key is None:
                self.stats["no_genus"] += 1
                continue
            b = assigner(to_float(row.get("max_ma")), to_float(row.get("min_ma")))
            if b is None:
                self.stats["no_bin"] += 1
                continue
            self.stats["assigned"] += 1
            self.sets[b].add(key)
            self.n_occ[b] += 1
            self._colls[b].add(row.get("collection_no"))

    @property
    def n_coll(self):
        return [len(c) for c in self._colls]


def merge_sets(stage_sets, bins, stage_counts=None):
    """Union the stage-level genus sets (and sum counts) into composite bins."""
    idx = {s[0]: i for i, s in enumerate(STAGES)}
    sets = [set().union(*[stage_sets[idx[s]] for s in b["stages"]]) for b in bins]
    if stage_counts is None:
        return sets
    counts = [sum(stage_counts[idx[s]] for s in b["stages"]) for b in bins]
    return sets, counts


# ------------------------------------------------------------------ rates
def three_timer(sets):
    """Alroy (2008) three-timer extinction rate at the TOP of each bin:
    mu_i = ln(2t_i / 3t_i) + ln(Ps_{i+1}); NaN where undefined."""
    n = len(sets)
    nan = float("nan")
    ps = [nan] * n
    mu = [nan] * n
    two = [0] * n
    three = [0] * n
    part = [0] * n
    for j in range(1, n - 1):
        a, b, c = sets[j - 1], sets[j], sets[j + 1]
        ab = a & b
        two[j] = len(ab)
        three[j] = len(ab & c)
        part[j] = len((a & c) - b)
        if three[j] + part[j] > 0:
            ps[j] = three[j] / (three[j] + part[j])
    for i in range(1, n - 1):
        p = ps[i + 1]
        if two[i] > 0 and three[i] > 0 and p == p and p > 0:
            mu[i] = math.log(two[i] / three[i]) + math.log(p)
    return mu


def gap_filler(sets):
    """Alroy (2014) gap-filler extinction rate at the TOP of bin i:
    ln((2t + pt) / (3t + pt + gf)), 2t = |S(i-1) & S(i)|, 3t = |S(i-1) & S(i) & S(i+1)|,
    pt = |(S(i-1) & S(i+1)) - S(i)|, gf = |(S(i-1) & S(i+2)) - S(i+1)|."""
    n = len(sets)
    mu = [float("nan")] * n
    for i in range(1, n - 2):
        a, b, c, d = sets[i - 1], sets[i], sets[i + 1], sets[i + 2]
        t2 = len(a & b)
        t3 = len(a & b & c)
        pt = len((a & c) - b)
        gf = len((a & d) - c)
        num, den = t2 + pt, t3 + pt + gf
        if num > 0 and den > 0:
            mu[i] = math.log(num / den)
    return mu


# ------------------------------------------------------------ temperature
def parse_phanda(text):
    """{slice: (GMST_16, GMST_50, GMST_84)} from the PhanDA percentiles CSV."""
    out = {}
    if text.startswith(chr(0xFEFF)):
        text = text[1:]
    rd = csv.DictReader(io.StringIO(text.lstrip("﻿")))
    for row in rd:
        row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
        name = row.get("Stage", "")
        vals = [to_float(row.get(c)) for c in ("GMST_16", "GMST_50", "GMST_84")]
        if name and all(v is not None for v in vals):
            out[name] = tuple(vals)
    return out


def split_normal(t16, t50, t84, z):
    """Value at standard-normal deviate z of a two-piece normal with the given
    16th/50th/84th percentiles."""
    t16, t50, t84, z = (np.asarray(v, float) for v in (t16, t50, t84, z))
    return np.where(z >= 0, t50 + z * (t84 - t50), t50 + z * (t50 - t16))


def boundary_rows(bins, mu3t, mugf, temps):
    """One row per boundary whose two neighbouring bins both have a PhanDA
    temperature. The boundary is the TOP of bin j; dT = T(j+1) - T(j)."""
    rows = []
    for j in range(len(bins) - 1):
        a, b = bins[j], bins[j + 1]
        if a["phanda"] is None or b["phanda"] is None:
            continue
        if a["phanda"] not in temps or b["phanda"] not in temps:
            continue
        ta, tb = temps[a["phanda"]], temps[b["phanda"]]
        d = tb[1] - ta[1]
        rows.append({"bin": a["name"], "top_ma": a["top"], "mu3t": mu3t[j],
                     "mugf": mugf[j], "T16_bin": ta[0], "T_bin": ta[1], "T84_bin": ta[2],
                     "T16_next": tb[0], "T_next": tb[1], "T84_next": tb[2],
                     "dT": d, "abs_dT": abs(d), "dur_bin": a["base"] - a["top"],
                     "dur_next": b["base"] - b["top"],
                     "big_five": a["name"] in BIG_FIVE_TOPS})
    return rows


# ------------------------------------------------------------- statistics
def rankdata(x):
    """Ranks 1..n, ties get their average rank."""
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    xs = x[order]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return r


def _unit(v):
    v = np.asarray(v, float) - np.mean(v)
    nv = float(np.linalg.norm(v))
    return v / nv if nv > 0 else None


def corr(a, b):
    ua, ub = _unit(a), _unit(b)
    if ua is None or ub is None:
        return float("nan")
    return float(ua @ ub)


def spearman(x, y):
    return corr(rankdata(x), rankdata(y))


def perm_p(a, b, obs, n_perm, rng):
    """One-sided permutation p for corr(a, b) >= obs."""
    ua, ub = _unit(a), _unit(b)
    if ua is None or ub is None or obs != obs:
        return 1.0
    exceed, done = 0, 0
    while done < n_perm:
        m = min(5000, n_perm - done)
        P = rng.permuted(np.tile(ub, (m, 1)), axis=1)
        exceed += int(np.sum(P @ ua >= obs - 1e-12))
        done += m
    return (1 + exceed) / (n_perm + 1)


def shift_p(a, b, obs, min_shift=3):
    """One-sided p from all circular shifts of b (|shift| >= min_shift), which
    keep each series' own autocorrelation. Returns (p, number of shifts)."""
    ua, ub = _unit(a), _unit(b)
    n = len(a)
    shifts = list(range(min_shift, n - min_shift + 1))
    if ua is None or ub is None or obs != obs or not shifts:
        return 1.0, len(shifts)
    vals = np.array([float(np.roll(ub, s) @ ua) for s in shifts])
    return (1 + int(np.sum(vals >= obs - 1e-12))) / (len(shifts) + 1), len(shifts)


def spearman_test(x, y, n_perm, rng, shift=True):
    x, y = np.asarray(x, float), np.asarray(y, float)
    out = {"n": int(len(x))}
    if len(x) < 4:
        out.update({"rho": float("nan"), "p_perm": 1.0})
        return out
    rx, ry = rankdata(x), rankdata(y)
    rho = corr(rx, ry)
    out["rho"] = rho
    out["p_perm"] = perm_p(rx, ry, rho, n_perm, rng)
    if shift:
        out["p_shift"], out["n_shifts"] = shift_p(rx, ry, rho)
    return out


def residualize(r, Z):
    Z = np.asarray(Z, float)
    if Z.ndim == 1:
        Z = Z[:, None]
    A = np.column_stack([np.ones(len(r))] + [rankdata(Z[:, k]) for k in range(Z.shape[1])])
    coef, *_ = np.linalg.lstsq(A, r, rcond=None)
    return r - A @ coef


def partial_spearman_test(x, y, Z, n_perm, rng):
    """Spearman correlation of x and y after removing the (rank) linear effect
    of the columns of Z from both; permutation p on the residuals."""
    ex = residualize(rankdata(x), Z)
    ey = residualize(rankdata(y), Z)
    rho = corr(ex, ey)
    return {"n": int(len(ex)), "rho": rho, "p_perm": perm_p(ex, ey, rho, n_perm, rng)}
