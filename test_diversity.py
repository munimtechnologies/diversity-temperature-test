# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for diversity.py. Run: python -m unittest -v test_diversity"""
import math
import unittest
from collections import Counter

import numpy as np

import diversity
import thermo


class Coverage(unittest.TestCase):
    def test_enumeration(self):
        # Occurrences A, A, b, c. The 6 pairs: AA (1 genus), 4 x A+single (2), bc (2).
        fk = diversity.freq_counts([2, 1, 1])
        cov, s = diversity.rarefied(fk, 4, 2)
        self.assertAlmostEqual(s, 11.0 / 6.0, places=10)
        self.assertAlmostEqual(cov, 1.0 / 3.0, places=10)

    def test_one_draw_is_one_genus(self):
        for counts in ([2, 1, 1], [5, 3, 3, 1], [40, 1, 1, 1, 7]):
            fk = diversity.freq_counts(counts)
            self.assertAlmostEqual(diversity.rarefied(fk, sum(counts), 1)[1], 1.0, places=10)

    def test_full_sample(self):
        counts = [5, 3, 3, 1]
        fk = diversity.freq_counts(counts)
        self.assertAlmostEqual(diversity.rarefied(fk, 12, 12)[1], 4.0, places=10)

    def test_sample_coverage(self):
        fk = diversity.freq_counts([2, 1, 1])
        self.assertAlmostEqual(diversity.sample_coverage(fk, 4), 0.625, places=12)
        self.assertEqual(diversity.sample_coverage(diversity.freq_counts([3, 2]), 5), 1.0)

    def test_quorum_search(self):
        # 4 genera x 50: coverage(m=4) = 0.6804 < 0.7 <= coverage(m=5) = 0.7607.
        out = diversity.coverage_richness([50, 50, 50, 50], 0.7)
        self.assertEqual(out["m"], 5)
        self.assertEqual(out["coverage"], 1.0)
        self.assertTrue(3.0 < out["richness"] < 3.1)

    def test_low_coverage_is_nan(self):
        out = diversity.coverage_richness([1] * 10, 0.7)
        self.assertTrue(math.isnan(out["richness"]))
        self.assertTrue(math.isnan(diversity.coverage_richness([], 0.7)["richness"]))

    def test_monotone(self):
        rng = np.random.default_rng(7)
        counts = list(rng.integers(1, 30, size=60))
        fk = diversity.freq_counts(counts)
        n = int(sum(counts))
        vals = [diversity.rarefied(fk, n, m) for m in range(1, n, 37)]
        covs, ss = [v[0] for v in vals], [v[1] for v in vals]
        self.assertTrue(all(b >= a - 1e-12 for a, b in zip(covs, covs[1:])))
        self.assertTrue(all(b >= a - 1e-12 for a, b in zip(ss, ss[1:])))
        self.assertTrue(ss[-1] <= 60 + 1e-9)

    def test_higher_quorum_more_genera(self):
        rng = np.random.default_rng(8)
        counts = list(rng.integers(1, 12, size=80)) + [200, 150]
        lo = diversity.coverage_richness(counts, 0.5)["richness"]
        hi = diversity.coverage_richness(counts, 0.8)["richness"]
        self.assertTrue(lo < hi <= 82)


class Counting(unittest.TestCase):
    def test_accumulator_dedup(self):
        a = thermo.BinAssigner()
        rows = [
            {"genus": "Claraia", "phylum": "Mollusca", "max_ma": "251.5", "min_ma": "251.0",
             "collection_no": "1"},
            {"genus": "Claraia", "phylum": "Mollusca", "max_ma": "251.5", "min_ma": "251.0",
             "collection_no": "1"},
            {"genus": "Claraia", "phylum": "Mollusca", "max_ma": "251.5", "min_ma": "251.0",
             "collection_no": "2"},
            {"genus": "Lingula", "phylum": "Brachiopoda", "max_ma": "251.5", "min_ma": "251.0",
             "collection_no": "2"},
            {"genus": "", "phylum": "Mollusca", "max_ma": "251.5", "min_ma": "251.0",
             "collection_no": "3"},
            {"genus": "Lingula", "phylum": "Brachiopoda", "max_ma": "300", "min_ma": "200",
             "collection_no": "4"},
        ]
        acc = diversity.CountAccumulator(a.n)
        acc.add_rows(rows, a)
        i = [s[0] for s in thermo.STAGES].index("Induan")
        self.assertEqual(acc.counts[i], Counter({"Mollusca|Claraia": 2, "Brachiopoda|Lingula": 1}))
        self.assertEqual(acc.n_coll[i], 2)
        self.assertEqual(acc.stats["same_genus_same_collection"], 1)
        self.assertEqual(acc.stats["no_genus"], 1)
        self.assertEqual(acc.stats["no_bin"], 1)

    def test_merge_counts(self):
        bins = thermo.build_bins()
        names = [s[0] for s in thermo.STAGES]
        sc = [Counter() for _ in names]
        sc[names.index("Danian")]["x"] = 2
        sc[names.index("Selandian")]["x"] = 3
        sc[names.index("Selandian")]["y"] = 1
        nc = [0] * len(names)
        nc[names.index("Danian")] = 4
        nc[names.index("Selandian")] = 5
        counts, ncoll = diversity.merge_counts(sc, bins, nc)
        j = [b["name"] for b in bins].index("Selandian/Danian")
        self.assertEqual(counts[j], Counter({"x": 5, "y": 1}))
        self.assertEqual(ncoll[j], 9)

    def test_origination_mirror(self):
        # Time-mirror of the three-timer extinction test in test_thermo.py.
        S = [{"a", "b"}, {"a", "b"}, {"a", "b"}, {"a", "b", "c"}, {"a", "b", "c"}]
        lam = diversity.three_timer_origination(S)
        self.assertAlmostEqual(lam[3], math.log(1.5), places=12)
        self.assertAlmostEqual(lam[2], 0.0, places=12)
        self.assertTrue(math.isnan(lam[0]) and math.isnan(lam[1]) and math.isnan(lam[4]))

    def test_detrend(self):
        x = np.arange(20.0)
        self.assertTrue(np.allclose(diversity.detrend(2.0 + 3.0 * x, x), 0.0))
        r = diversity.detrend(np.sin(x) + 0.5 * x, x)
        self.assertAlmostEqual(float(np.mean(r)), 0.0, places=10)
        self.assertAlmostEqual(float(r @ (x - x.mean())), 0.0, places=8)


def make_rows(seed, slope, pt_drop=True, noise=0.01):
    rng = np.random.default_rng(seed)
    bins = thermo.build_bins()
    temps, k = {}, 20.0
    for b in bins:
        if b["phanda"]:
            k += rng.choice([-1.0, 1.0]) * rng.uniform(0.5, 6.0)
            temps[b["phanda"]] = (k - 1.0, k, k + 1.0)
    base = []
    for b in bins:
        if not b["phanda"]:
            base.append(float("nan"))
        elif slope == 0.0:
            base.append(math.exp(rng.standard_normal()))
        else:
            base.append(math.exp(slope * temps[b["phanda"]][1] + noise * rng.standard_normal()))
    names = [b["name"] for b in bins]
    i_ind, i_chg = names.index("Induan"), names.index("Changhsingian")
    base[i_ind] = base[i_chg] * (0.5 if pt_drop else 2.0)
    rich = {q: list(base) for q in (diversity.QUORUM,) + tuple(diversity.QUORUMS_ALT)}
    raw = [int(rng.integers(50, 500)) for _ in bins]
    n_coll = [int(rng.integers(20, 300)) for _ in bins]
    lam = [float(rng.normal()) for _ in bins]
    return diversity.level_rows(bins, rich, raw, n_coll, lam, temps)


class Analyze(unittest.TestCase):
    def test_rows(self):
        rows = make_rows(1, 0.1)
        self.assertEqual(len(rows), 79)
        self.assertEqual(rows[0]["bin"], "Tremadocian")
        self.assertTrue(math.isnan(rows[0]["dT_in"]))
        self.assertTrue(np.isfinite(rows[1]["dT_in"]))
        d = diversity.diff_rows(rows, "R")
        self.assertEqual(len(d), 78)
        self.assertEqual(sum(r["big_five"] for r in d), 6)

    def test_planted_positive(self):
        res = diversity.analyze(make_rows(3, 0.1), 300, np.random.default_rng(4), n_draws=20)
        self.assertTrue(res["sanity"]["passed"])
        self.assertGreater(res["D1_levels"]["rho"], 0.8)
        self.assertGreater(res["D2_differences"]["rho"], 0.8)
        self.assertEqual(res["D2_differences"]["n"], 78)
        self.assertEqual(res["D2b_differences_without_big_five"]["n"], 72)
        self.assertEqual(res["decision"], "SUPPORTED")
        self.assertTrue(res["D6_quorum"]["same_signs"])

    def test_planted_negative(self):
        res = diversity.analyze(make_rows(5, -0.1), 300, np.random.default_rng(6), n_draws=20)
        self.assertLess(res["D1_levels"]["rho"], -0.8)
        self.assertLess(res["D2_differences"]["rho"], -0.8)
        self.assertEqual(res["decision"], "REVERSED")

    def test_null(self):
        res = diversity.analyze(make_rows(9, 0.0), 300, np.random.default_rng(10), n_draws=20)
        for key in ("D1_levels", "D2_differences", "D5a_origination_vs_level",
                    "D5b_origination_vs_change_in", "D7_raw_levels"):
            self.assertTrue(-1.0 <= res[key]["rho"] <= 1.0)
            self.assertTrue(0.0 < res[key]["p_perm"] <= 1.0)
            self.assertTrue(0.0 < res[key]["p_perm_lower"] <= 1.0)
        self.assertTrue(0.0 <= res["D4_uncertainty"]["frac_rho_pos"] <= 1.0)
        self.assertIn(res["decision"], ("SUPPORTED", "PARTLY_SUPPORTED", "REVERSED", "MIXED",
                                        "NOT_SUPPORTED"))

    def test_sanity_gate(self):
        res = diversity.analyze(make_rows(3, 0.1, pt_drop=False), 100,
                                np.random.default_rng(1), n_draws=10)
        self.assertFalse(res["sanity"]["passed"])
        self.assertEqual(res["decision"], "DATA_LIMITED")

    def test_too_few(self):
        rows = make_rows(3, 0.1)
        for r in rows:
            r["R"] = float("nan")
        res = diversity.analyze(rows, 100, np.random.default_rng(1))
        self.assertEqual(res["decision"], "DATA_LIMITED")
        self.assertEqual(res["n_levels"], 0)


if __name__ == "__main__":
    unittest.main()
