# RESULTS: does marine biodiversity track global temperature on PhanDA?

Graded 2026-11-17 from the output of job `diversity-temperature-v1` (ran 20:29-20:34 UTC).
Plan: [PREREG.md](PREREG.md) (committed 11:24 UTC). Reading:
[grade_diversity.py](grade_diversity.py) and [PRE-RESULT-NOTE-1.md](PRE-RESULT-NOTE-1.md)
(committed 18:16 UTC, before the job started). The grader's rules were applied to
`results.json` and `bins.csv` by hand, exactly as coded; `test_graded_results.py` runs the
grader itself on the published output and checks that it gives the same labels.
Made with AI agents working in MV2's name.

## Verdict

**NOT_SUPPORTED.** Registered reading: *no registered signal holds; consistent with Balembois
et al. (2025), global temperature barely sets long-term global richness.*

Mayhew, Bell, Benton & McGowan (2012, PNAS) reported that, once sampling is controlled,
marine diversity is higher when Earth is warmer. On PhanDA (Judd et al. 2024), the best
global temperature history now available, and coverage-standardised 2026 PBDB genus richness
(79 bins, Tremadocian to Piacenzian), we find no link at all:

- Levels (detrended log richness vs temperature): rho = -0.06, circular-shift p = 0.57.
- Changes (change in log richness vs change in temperature): rho = +0.05, p = 0.33.
- Without the Big Five boundaries: rho = +0.07, p = 0.28.
- At the other registered coverage quorums (0.5, 0.8), and on raw richness, the correlations
  stay between -0.11 and +0.05.

Neither the positive relation of 2012 nor the negative one reported before it appears.

**One registered secondary signal holds:** origination (the three-timer rate of new genera)
is higher when the world is warmer (D5a: rho = 0.28, circular-shift p = 0.027, permutation
p = 0.006). We gave this 20%. Origination does not rise after warming steps (D5b: rho -0.08).
Warmer worlds may make more new genera without holding more genera, if the extra originations
are balanced by extinction; this is one row of thirteen and needs its own test before anyone
leans on it (see limits).

## Graded predictions

| # | Prediction | Prior | Holds? |
|---|---|---|---|
| S | Sanity: R defined for >= 60 bins AND R(Induan) < R(Changhsingian) | 80% | YES |
| D1+ | Mayhew 2012: richness level tracks temperature level (rho > 0, p_shift < 0.05) | 15% | no |
| D1- | Reverse: warm = fewer genera (rho < 0, lower p_shift < 0.05) | 30% | no |
| D2+ | Richness changes track temperature changes (rho > 0, p_perm < 0.05) | 10% | no |
| D2- | Warming steps coincide with richness losses (rho < 0, lower p_perm < 0.05) | 30% | no |
| D2b | D2 without the 6 Big Five boundaries is significant in either direction | 20% | no |
| D3 | D1's sign survives removing the rank effect of bin duration | 75% | YES |
| D4 | D1 robust to PhanDA uncertainty (>= 90% of draws give D1's sign) | 50% | no |
| D5a | Origination rises with temperature level (rho > 0, p_shift < 0.05) | 20% | YES |
| D5b | Origination rises after warming into the bin (rho > 0, p_perm < 0.05) | 25% | no |
| D6 | D1 and D2 keep their signs at quorums 0.5 and 0.8 | 70% | no |
| D7a | Raw richness correlates less positively with T than standardised richness | 50% | no |
| D7b | Raw richness level is negatively related to T | 45% | YES |

Brier score of our registered priors: 0.181 over 13 rows (a 50% forecaster scores 0.250).
The job's label and its four signals (D1+/D1-/D2+/D2-, all false) agree with the grading.

## Key numbers

| Test | n | rho | one-sided p (upper / lower) |
|---|---|---|---|
| D1 levels, quorum 0.7 | 79 | -0.059 | shift 0.573 / 0.440 |
| D2 differences | 78 | 0.050 | perm 0.329 / 0.667 |
| D2b differences without Big Five | 72 | 0.071 | perm 0.277 / 0.725 |
| D3 levels, duration partial | 79 | -0.088 | perm 0.781 |
| D5a origination vs temperature | 79 | 0.282 | shift 0.027 (perm 0.006) |
| D5b origination vs warming into bin | 78 | -0.079 | perm 0.752 |
| D6 levels at quorum 0.5 / 0.8 | 79 / 78 | -0.106 / 0.003 | |
| D6 differences at quorum 0.5 / 0.8 | 78 / 76 | 0.047 / 0.047 | |
| D7 raw richness levels | 79 | -0.006 | shift 0.400 / 0.613 |
| D7 collections (effort) levels | 79 | 0.136 | shift 0.333 |

- D4: across 2,000 PhanDA draws, rho 5-95% is -0.138 to 0.031; 86% share D1's sign
  (registered bar 90%).
- D6 fails only because the levels rho at quorum 0.8 is +0.003 instead of negative: a sign
  flip of a correlation that is zero at every quorum.
- Sanity: standardised richness (quorum 0.7) falls from 196.9 genera in the Changhsingian to
  40.6 in the Induan.

## What this means, and its limits

- The 2012 result does not carry over to PhanDA and today's PBDB. The registered reading is
  the one we printed in advance: over hundreds of millions of years, global mean temperature
  does not set global marine genus richness at stage resolution, in line with the mechanistic
  model of Balembois et al. (2025, Nature Communications).
- Limits: 79 bins and a strict circular-shift null (with 79 bins, p_shift < 0.05 needs the
  observed rho to beat all but at most 2 of 74 shifts) give modest power; a weak true link
  (|rho| ~ 0.15) could be missed. Mayhew et al. used a different temperature proxy (a d18O
  compilation), binning and timescale, so this is a test of the claim on better data, not
  a re-run of their analysis. Coverage-based rarefaction controls sampling intensity but not
  geographic spread of sampling.
- D5a is one of thirteen registered rows; with no correction for multiple rows it should be
  read as a lead. Post-extinction rebounds in warm intervals (the Early Triassic after the
  end-Permian, origination 0.97 and 1.17) may drive part of it; a follow-up would register
  the test with the Big Five aftermath bins removed.
- The richness series can also be built from the public per-stage series of
  [extinction-cycle-test](https://github.com/munimtechnologies/extinction-cycle-test)
  (same PBDB download); it was never paired with PhanDA before this job (disclosed in the
  pre-result note).
- Hand checks: PhanDA GMST_50 in `bins.csv` matches the grader's hand-read table for the
  bins checked (Tremadocian 32.103 vs 32.10, Changhsingian 18.019 vs 18.02, Induan 28.499 vs
  28.50, Piacenzian 15.632 vs 15.63; Cambrian and Quaternary bins carry no temperature, as
  designed); the grader test checks all 79.

## Files

`results/`: `results.json` (every test), `bins.csv` (richness at three quorums, coverage,
collections, origination and PhanDA temperature per bin), `differences.csv`,
`provenance.json` (PBDB and PhanDA URLs, md5 sums, times).

Code Apache-2.0; text CC BY 4.0. MV2, Munim, Inc., mv2@munimtech.com
