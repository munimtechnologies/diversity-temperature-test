# Pre-result note 1 (2026-11-17, about 18:30 UTC): how the result will be read

Written by MV2 of Munim, Inc. (made with AI agents working in his name) **before the analysis
job `diversity-temperature-v1` had produced any output**. No PBDB richness, coverage or
origination value from this project has been computed or seen. Nothing in `PREREG.md` is
changed by this note. It does three things: it fixes the reading of the result in code
(`grade_diversity.py`), it says which registered predictions could already be computed from
public inputs (none), and it states a property of the registered shift test that readers
should know before they see a p value.

## 1. The reading is now code

`grade_diversity.py` (tests: `test_grade_diversity.py`) turns the job's `results.json` into
the graded prior table and the decision label. It computes no new statistic. Where the
PREREG table left a choice open, the script fixes it here:

- **D1+/D1-** and **D5a** use the one-sided circular-shift p in the stated direction
  (`p_shift` for rho > 0, `p_shift_lower` for rho < 0); **D2+/D2-** and **D5b** use the
  one-sided permutation p. A row holds only if rho is finite, has the stated sign, and that
  p is below 0.05.
- **D2b** ("significant in either direction") holds if the Big-Five-free differences test
  has rho > 0 with `p_perm` < 0.05 or rho < 0 with `p_perm_lower` < 0.05. It is *not* read as
  a two-sided test at 0.10.
- **D3**, **D4**, **D6** require D1's (and for D6 also D2's) rho to be finite and non-zero;
  a zero or undefined rho makes them false, never true.
- **D7a/D7b** compare raw and standardised rho values only, with no significance
  requirement, exactly as registered.
- The script re-derives the label and the four signals (D1+, D1-, D2+, D2-) and compares
  them with what the job wrote. Any disagreement gives status `CHECK_CODE`. With
  `--bins bins.csv` it also checks the job's 79 PhanDA medians against the values read by
  hand on 2026-11-17 (the table in extinction-temperature-test's PRE-RESULT-NOTE-1.md);
  any mismatch or missing bin gives `CHECK_DATA`. Neither status changes the label; both stop
  the result from being reported until a person has looked.

## 2. Nothing in the table can be computed from public inputs alone

For the sister test (extinction-temperature-test), two registered rows depended only on
PhanDA and were computed and disclosed before the result. Here, every row (S, D1 to D7)
needs coverage-standardised richness, raw richness or three-timer origination computed by
this project's pipeline from PBDB. None of those numbers exists yet. So no row has been
answered in advance.

One fact from public PhanDA bears on how to read D2, so it is stated now. The largest
temperature step between adjacent bins is Changhsingian to Induan (GMST_50 18.02 to
28.50 C, +10.48 C). The sanity gate requires richness to fall across that same boundary.
So, if the gate passes, the single largest warming step comes with a richness loss. That
pushes D2 towards a negative rho. In a rank correlation over about 74 boundaries, one point
moves rho by only a few hundredths, and D2b repeats the test without the six Big Five
boundaries. D2b is the row to read if D2- holds.

## 3. What the circular-shift test can and cannot show

D1 and D5a get their p values from all circular shifts of 3 or more steps (`thermo.shift_p`):
with n bins there are n - 5 shifts, and p = (1 + number of shifts at least as extreme) /
(n - 4). So:

- the smallest possible p_shift is 1/(n - 4), about 0.013 for n = 79;
- for n >= 65 bins, p_shift < 0.05 needs the observed rho to be beaten or tied by **at most
  2** of the shifts; for 60 to 64 bins, by at most 1.

This is a deliberately strict null. PhanDA has long, slow swings: Ordovician to Devonian
warmth, the late Palaeozoic ice age, Mesozoic hothouse, Cenozoic cooling. A richness curve
with any slow swing of its own will match some shifted copy of the temperature curve well.
So a non-significant p_shift together with a large rho means "not distinguishable from two
slowly wandering series". It does not mean "no relationship". The permutation p for D1 is
reported next to it but does not decide D1, as registered.

## What was seen before writing this note

The PhanDA medians of the 79 bins (as in the sister note), the project's own code, and the
PREREG. No PBDB-derived number.
