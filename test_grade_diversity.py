# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for grade_diversity.py (standard library only).
Made with AI agents working in MV2's name."""
import copy
import json
import os
import tempfile
import unittest

import grade_diversity as gd

NAN = float("nan")


def _t(rho, pp=0.5, ppl=0.5, ps=None, psl=None, n=79):
    d = {"n": n, "rho": rho, "p_perm": pp, "p_perm_lower": ppl}
    if ps is not None:
        d["p_shift"] = ps
        d["p_shift_lower"] = psl
    return d


def make_results(d1=None, d2=None, passed=True, decision=None, signals=None):
    """A results.json shaped like job_diversity.py's output. Default: nothing
    significant, D1 and D2 positive but weak."""
    a = {
        "n_levels": 75, "bins_with_richness": [],
        "sanity": {"n_ok": passed, "pt_drop": True, "passed": passed},
        "D1_levels": d1 or _t(0.10, 0.2, 0.8, 0.3, 0.7),
        "D2_differences": d2 or _t(0.05, 0.3, 0.7, 0.4, 0.6, n=74),
        "D2b_differences_without_big_five": _t(0.04, 0.35, 0.65, 0.4, 0.6, n=68),
        "D3_duration_partial": {"n": 75, "rho": 0.08, "p_perm": 0.25},
        "D4_uncertainty": {"n_draws": 2000, "rho_p05": -0.05, "rho_median": 0.09,
                           "rho_p95": 0.2, "frac_rho_pos": 0.85,
                           "frac_same_sign_as_D1": 0.85},
        "D5a_origination_vs_level": _t(0.02, 0.4, 0.6, 0.45, 0.55, n=70),
        "D5b_origination_vs_change_in": _t(0.03, 0.4, 0.6, 0.45, 0.55, n=69),
        "D6_quorum": {"R_0.5": {"levels": _t(0.12), "differences": _t(0.06)},
                      "R_0.8": {"levels": _t(0.09), "differences": _t(0.04)},
                      "same_signs": True},
        "D7_raw_levels": _t(-0.2, 0.9, 0.1, 0.9, 0.1),
        "D7_raw_below_standardised": True,
        "D7_collections_levels": _t(-0.3, 0.95, 0.05, 0.9, 0.1),
        "qualifiers": {"quorum_robust": True},
    }
    g = gd.grade_rows(a)
    a["signals"] = signals if signals is not None else {
        "D1_pos": g["D1+"], "D1_neg": g["D1-"], "D2_pos": g["D2+"], "D2_neg": g["D2-"]}
    lab = decision if decision is not None else gd.label_from(g)
    a["decision"] = lab
    return {"status": "ok", "decision": lab, "analysis": a}


class TestRows(unittest.TestCase):
    def test_default_not_supported(self):
        gr = gd.grade(make_results())
        self.assertEqual(gr["label"], "NOT_SUPPORTED")
        self.assertEqual(gr["status"], "GRADED")
        rows = {r["id"]: r["holds"] for r in gr["rows"]}
        self.assertFalse(rows["D1+"])
        self.assertFalse(rows["D1-"])
        self.assertTrue(rows["D3"])
        self.assertFalse(rows["D4"])     # 0.85 < 0.9
        self.assertTrue(rows["D6"])
        self.assertTrue(rows["D7a"])     # -0.2 < 0.10
        self.assertTrue(rows["D7b"])

    def test_supported(self):
        r = make_results(d1=_t(0.4, 0.001, 0.999, 0.02, 0.98),
                         d2=_t(0.3, 0.004, 0.996, 0.1, 0.9))
        gr = gd.grade(r)
        self.assertEqual(gr["label"], "SUPPORTED")

    def test_d1_uses_shift_not_perm(self):
        # permutation p tiny, shift p not: D1+ must NOT hold.
        r = make_results(d1=_t(0.4, 0.0001, 0.9999, 0.2, 0.8))
        rows = {x["id"]: x["holds"] for x in gd.grade(r)["rows"]}
        self.assertFalse(rows["D1+"])
        self.assertEqual(gd.grade(r)["label"], "NOT_SUPPORTED")

    def test_reversed_and_mixed(self):
        rev = make_results(d1=_t(-0.4, 0.999, 0.001, 0.98, 0.02),
                           d2=_t(-0.05, 0.6, 0.4))
        self.assertEqual(gd.grade(rev)["label"], "REVERSED")
        mix = make_results(d1=_t(0.4, 0.001, 0.999, 0.02, 0.98),
                           d2=_t(-0.3, 0.996, 0.004))
        self.assertEqual(gd.grade(mix)["label"], "MIXED")
        part = make_results(d2=_t(0.3, 0.004, 0.996))
        self.assertEqual(gd.grade(part)["label"], "PARTLY_SUPPORTED")

    def test_d2b_either_direction(self):
        r = make_results()
        r["analysis"]["D2b_differences_without_big_five"] = _t(-0.3, 0.99, 0.01)
        rows = {x["id"]: x["holds"] for x in gd.grade(r)["rows"]}
        self.assertTrue(rows["D2b"])
        # Not a two-sided test at 0.10: one-sided 0.07 does not count.
        r["analysis"]["D2b_differences_without_big_five"] = _t(-0.3, 0.93, 0.07)
        rows = {x["id"]: x["holds"] for x in gd.grade(r)["rows"]}
        self.assertFalse(rows["D2b"])

    def test_nan_rho_never_holds(self):
        r = make_results(d1=_t(NAN, 1.0, 1.0, 1.0, 1.0))
        rows = {x["id"]: x["holds"] for x in gd.grade(r)["rows"]}
        for k in ("D1+", "D1-", "D3", "D4", "D6", "D7a"):
            self.assertFalse(rows[k], k)

    def test_d6_sign_flip(self):
        r = make_results()
        r["analysis"]["D6_quorum"]["R_0.8"]["differences"]["rho"] = -0.01
        rows = {x["id"]: x["holds"] for x in gd.grade(r)["rows"]}
        self.assertFalse(rows["D6"])

    def test_d4_threshold(self):
        r = make_results()
        r["analysis"]["D4_uncertainty"]["frac_same_sign_as_D1"] = 0.95
        rows = {x["id"]: x["holds"] for x in gd.grade(r)["rows"]}
        self.assertTrue(rows["D4"])


class TestStatus(unittest.TestCase):
    def test_sanity_failed(self):
        gr = gd.grade(make_results(passed=False))
        self.assertEqual(gr["label"], "DATA_LIMITED")
        self.assertEqual(gr["status"], "GRADED_SANITY_FAILED")

    def test_early_stop(self):
        r = {"status": "ok", "decision": "DATA_LIMITED",
             "analysis": {"n_levels": 5, "sanity": {"n_ok": False, "pt_drop": False,
                                                    "passed": False},
                          "decision": "DATA_LIMITED"}}
        gr = gd.grade(r)
        self.assertEqual(gr["label"], "DATA_LIMITED")
        self.assertIsNone({x["id"]: x["holds"] for x in gr["rows"]}["D1+"])

    def test_job_failed(self):
        gr = gd.grade({"status": "pbdb_incomplete: missing x"})
        self.assertEqual(gr["status"], "JOB_FAILED")
        self.assertIn("Nothing is graded", gd.render_md(gr, {}))

    def test_check_code_on_label(self):
        gr = gd.grade(make_results(decision="SUPPORTED"))
        self.assertEqual(gr["status"], "CHECK_CODE")
        self.assertEqual(gr["label"], "NOT_SUPPORTED")

    def test_check_code_on_signal(self):
        sig = {"D1_pos": True, "D1_neg": False, "D2_pos": False, "D2_neg": False}
        gr = gd.grade(make_results(signals=sig, decision="NOT_SUPPORTED"))
        self.assertEqual(gr["status"], "CHECK_CODE")
        self.assertEqual(gr["signals_disagreeing_with_job"], ["D1_pos"])


class TestBins(unittest.TestCase):
    def _rows(self):
        rows = [{"bin": "Fortunian", "GMST_50": ""}]
        rows += [{"bin": k, "GMST_50": f"{v:.4f}"} for k, v in gd.HAND_GMST50.items()]
        rows.append({"bin": "Quaternary", "GMST_50": ""})
        return rows

    def test_pass(self):
        pc = gd.check_bins(self._rows())
        self.assertTrue(pc["passed"])
        self.assertEqual(pc["n_checked"], 79)

    def test_mismatch_and_missing(self):
        rows = self._rows()
        rows = [r for r in rows if r["bin"] != "Induan"]
        for r in rows:
            if r["bin"] == "Visean":
                r["GMST_50"] = "13.50"
        pc = gd.check_bins(rows)
        self.assertFalse(pc["passed"])
        self.assertEqual(pc["missing_bins"], ["Induan"])
        self.assertEqual([m["bin"] for m in pc["mismatches"]], ["Visean"])
        gr = gd.grade(make_results(), rows)
        self.assertEqual(gr["status"], "CHECK_DATA")

    def test_hand_table_size(self):
        self.assertEqual(len(gd.HAND_GMST50), 79)


class TestMainAndBrier(unittest.TestCase):
    def test_brier(self):
        b = gd.brier({"S": True, "D1+": False}, gd.PRIORS)
        self.assertEqual(b["n"], 2)
        self.assertAlmostEqual(b["brier"], ((0.8 - 1) ** 2 + 0.15 ** 2) / 2, places=9)

    def test_main_writes_files(self):
        r = make_results(d1=_t(0.4, 0.001, 0.999, 0.02, 0.98))
        r["analysis"]["D3_duration_partial"]["rho"] = NAN
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "results.json")
            with open(p, "w", encoding="utf-8") as f:
                f.write(json.dumps(r))        # NaN written as NaN, as the job does
            self.assertEqual(gd.main(["--results", p, "--out-dir", d]), 0)
            with open(os.path.join(d, "graded.json"), encoding="utf-8") as f:
                gr = json.load(f)
            self.assertEqual(gr["label"], "PARTLY_SUPPORTED")
            with open(os.path.join(d, "GRADED.md"), encoding="utf-8") as f:
                md = f.read()
            self.assertIn("PARTLY_SUPPORTED", md)
            self.assertIn("| D1+ |", md)

    def test_input_not_mutated(self):
        r = make_results()
        before = copy.deepcopy(r)
        gd.grade(r)
        self.assertEqual(json.dumps(before, sort_keys=True), json.dumps(r, sort_keys=True))


if __name__ == "__main__":
    unittest.main()
