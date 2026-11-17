# diversity-temperature-test

**Does marine biodiversity track global temperature over the last 485 million years?**

This is a pre-registered test by MV2 of Munim, Inc. It was made with AI agents working in
his name. License: Apache-2.0 (code); the text and result tables are CC BY 4.0.

In 2012, Mayhew, Bell, Benton & McGowan (PNAS 109:15141) reported that once sampling is
controlled, the number of marine invertebrate taxa rises with temperature over the
Phanerozoic. That reversed an earlier negative result. A mechanistic model from 2025
(Balembois et al., Nature Communications 16:8498) instead finds that temperature barely
affects long-term global diversity. This repository tests the claim on two inputs:

- the best global temperature history now available, **PhanDA** (Judd et al. 2024,
  Science), which was not available in 2012;
- modern **Paleobiology Database** genus occurrences, standardised by coverage-based
  rarefaction (Chao & Jost 2012).

The hypotheses, priors and decision rule were fixed in [PREREG.md](PREREG.md) before any
richness value was computed.

## Status

- 2026-11-17: pre-registered. The analysis job is queued. Results will appear in
  RESULTS.md whatever the outcome.
- 2026-11-17 (later): the reading of the result was fixed in code before any result existed:
  [grade_diversity.py](grade_diversity.py) (tests: `test_grade_diversity.py`). See
  [PRE-RESULT-NOTE-1.md](PRE-RESULT-NOTE-1.md). It explains the reading choices and why no
  registered prediction could be computed from public data in advance. It also explains how
  strict the circular-shift test is: with about 79 bins, p_shift < 0.05 needs the observed
  correlation to beat all but at most 2 of the 74 shifted copies.

## Part of a series: what drives the history of life?

- [extinction-cycle-test](https://github.com/munimtechnologies/extinction-cycle-test):
  is extinction periodic (27 Myr and 62 Myr in the oceans, 27.5 Myr on land)?
- [extinction-temperature-test](https://github.com/munimtechnologies/extinction-temperature-test):
  does temperature change predict extinction (Song et al. 2021, the 5.2 C threshold)?
- this repository: does diversity itself track temperature?

## Run it

```
pip install -r requirements.txt
python -m unittest -v test_diversity test_grade_diversity
python job_diversity.py --data-dir data --out-dir out
python grade_diversity.py --results out/results.json --bins out/bins.csv --out-dir out
```

The job downloads eight PBDB occurrence CSVs (one per group, cached for later runs) and one PhanDA CSV. It writes
`out/results.json` (every test), `out/bins.csv` (richness, coverage, origination and
temperature per bin), `out/differences.csv` and `out/provenance.json` (URLs, md5s and
times). The grader writes `out/graded.json` and `out/GRADED.md`.

## Files

- `diversity.py`: coverage-based rarefaction, origination rate and the registered tests.
- `thermo.py`: stage table, PhanDA bins and parser, and rank statistics. These are shared
  with extinction-temperature-test.
- `job_thermo.py`: downloads, shared with extinction-temperature-test.
- `job_diversity.py`: the end-to-end job.
- `grade_diversity.py`: the registered reading (prior table, label, Brier score, PhanDA
  hand check). Standard library only.
- `test_diversity.py`, `test_grade_diversity.py`: offline tests, including a
  hand-enumerated rarefaction case, planted-signal checks and every grading rule.

Data: the Paleobiology Database (CC BY 4.0) and PhanDA (Judd et al. 2024, CC BY 4.0).

Contact: MV2, Munim, Inc., mv2@munimtech.com
