# PREREG: does marine biodiversity track global temperature on an independent temperature history?

Registered 2026-11-17 by MV2 of Munim, Inc., before the analysis job ran and before any
PBDB richness or origination rate was paired with any PhanDA temperature. Code:
`diversity.py`, `job_diversity.py` (with `thermo.py` and `job_thermo.py` shared from
extinction-temperature-test), tests `test_diversity.py`. Made with AI agents working in
MV2's name.

Disclosure of what was seen before registering: while preparing the sister test
(extinction-temperature-test, registered earlier the same day) the project read the
PhanDA percentile CSV's header, its 85 slice names and ages, and the GMST values of 20
rows (Holocene to Rupelian; Sandbian to Tremadocian). No richness, coverage or origination
value from this project's data has been computed or looked at. The priors below draw on
general knowledge of the published PhanDA curve (Judd et al. 2024, Fig. 2) and of
published sampling-standardised diversity curves (e.g. Alroy et al. 2008, Alroy 2010;
Fan et al. 2020).

## The claim under test

Mayhew, Bell, Benton & McGowan (2012), *Biodiversity tracks temperature over time*,
PNAS 109:15141-15145 (doi:10.1073/pnas.1200844109):

- With sampling effort controlled (sample-standardised marine invertebrate richness),
  "taxonomic richness increases, not decreases, with temperature" over the Phanerozoic.
- This reverses the negative association found with unstandardised data by Mayhew, Jenkins
  & Benton (2008, Proc. R. Soc. B 275:47-53), which the 2012 paper attributes to
  sampling bias.
- Increasing temperatures also predict origination and extinction rates.

Their temperatures were tropical sea-surface series from oxygen isotopes and a modelled
seawater temperature (time series of N = 51 bins in their figures), detrended before
correlation. The question matters beyond palaeontology: it is often cited when asking
whether a warmer world is, in the long run, a more or less diverse one.

## Why this test is new

- **PhanDA** (Judd et al. 2024, Science 385:eadk3705; data CC BY 4.0,
  github.com/EJJudd/PhanDA) is a global mean surface temperature (GMST) history for
  485 Myr, built by assimilating the PhanSST proxy database into HadCM3L climate model
  priors. It is much better constrained than the curves available in 2012, and it
  shows a far larger range (about 11-36 C).
- Searches on 2026-11-17 found no study that correlates PBDB global richness with PhanDA.
  Closest work:
  - Balembois et al. (2025, Nature Communications 16:8498) use a mechanistic model and
    conclude that temperature "barely affects" long-term global biodiversity. They used
    HadCM3 climate runs, including two of Judd et al.'s series as a sensitivity check,
    not the PhanDA GMST product, and they did not correlate fossil richness with
    temperature directly.
  - Benson et al. (2025, Palaeontology 68:e70006) study spatial (latitudinal and
    temperature) gradients in grid cells, not the global time series.
  - Spiridonov & Lovejoy (2022, Nature, doi:10.1038/s41586-022-04867-y) find that
    climate (from d18O) drives genus origination and extinction below about 40 Myr and
    that life drives it above.
  - Mathes et al. (2021, PNAS, doi:10.1073/pnas.2105769118) relate origination to
    multi-scale cooling.
  - Shang et al. (2025, Palaeontology, doi:10.1111/pala.12734) correlate
    "accelerations" (second derivatives) of richness and temperature. That is a different
    test, and we could not confirm which temperature series it uses.
- So this is the first test of Mayhew et al.'s level and change claims on PhanDA, with
  coverage-based standardisation of modern PBDB data.

## What is measured

- **Bins.** The same 79 PhanDA bins as extinction-temperature-test (Tremadocian to
  Piacenzian; the five Quaternary slices and the PETM slice are dropped). PBDB occurrences
  are binned into GTS2020 stages, and stages are merged where PhanDA merged them.
- **Data.** PBDB marine genus occurrences of Porifera, Cnidaria, Bryozoa, Brachiopoda,
  Mollusca, Echinodermata, Trilobita and Graptolithina (envtype=marine, genus certain,
  regular preservation). These are the same queries and cache as the sister tests.
- **Standardised richness R.** Within each bin, a genus's frequency is the number of
  collections it occurs in. R is the expected number of genera in a random subsample of
  these genus-in-collection records whose expected coverage reaches the quorum: 0.7
  primary, 0.5 and 0.8 as robustness checks. It uses Chao & Jost's (2012) interpolation
  formulas, the analytic counterpart of shareholder quorum subsampling. Bins whose
  estimated coverage is below the quorum get no R (no extrapolation).
- **Temperature T.** PhanDA GMST median (50th percentile) of the bin. Uncertainty comes
  from the 16th and 84th percentiles.
- **Origination.** Alroy (2008) three-timer origination rate at the base of each bin,
  computed on the eight-group genus lists.
- **Raw richness.** Genera observed in the bin, unstandardised. **Effort** is the number
  of collections.

## Hypotheses, tests and honest priors

All correlations are Spearman, with seed 20261117 and 20,000 permutations.
- "Levels" are ln R and T, each detrended by a straight line in bin mid-age. Their
  significance comes from all circular shifts of 3 or more steps (p_shift), which keeps
  each series' autocorrelation. Plain permutation p is also reported.
- "Differences" are d ln R and dT between adjacent bins (78 boundaries). Their
  significance is from permutation (p_perm), with p_shift reported too.
- Every test is reported in both directions (one-sided p for rho > 0 and for rho < 0).

| # | Prediction | Test | Prior that it holds |
|---|---|---|---|
| S | Sanity: R defined for >= 60 bins AND R(Induan) < R(Changhsingian) | gate | 80% |
| D1+ | Mayhew 2012: richness level tracks temperature level | levels rho > 0, p_shift < 0.05 | 15% |
| D1- | Reverse (Mayhew 2008 direction): warm = fewer genera | levels rho < 0, lower p_shift < 0.05 | 30% |
| D2+ | Richness changes track temperature changes | differences rho > 0, p_perm < 0.05 | 10% |
| D2- | Warming steps coincide with richness losses | differences rho < 0, lower p_perm < 0.05 | 30% |
| D2b | D2 with the 6 Big Five boundaries removed is significant in either direction | p < 0.05 | 20% |
| D3 | D1's sign survives removing the rank effect of bin duration | same sign | 75% |
| D4 | D1 is robust to PhanDA uncertainty (2,000 independent two-piece normal draws per slice) | >= 90% of draws give D1's sign | 50% |
| D5a | Mayhew 2012: origination rises with temperature level | rho > 0, p_shift < 0.05 | 20% |
| D5b | Origination rises after warming into the bin | rho > 0, p_perm < 0.05 | 25% |
| D6 | D1 and D2 keep their signs at quorums 0.5 and 0.8 | same signs | 70% |
| D7a | Mayhew 2012's sampling argument: raw richness correlates less positively with T than standardised richness | rho_raw < rho_R | 50% |
| D7b | Raw richness level is negatively related to T | rho_raw < 0 | 45% |

## Decision rule (fixed now)

- `DATA_LIMITED` if the sanity gate fails.
- `SUPPORTED` if D1+ and D2+ both hold.
- `PARTLY_SUPPORTED` if exactly one of D1+ or D2+ holds and neither D1- nor D2- holds.
- `REVERSED` if D1- or D2- holds and neither D1+ nor D2+ holds.
- `MIXED` if a positive and a negative signal both hold (for example, levels positive but
  changes negative).
- `NOT_SUPPORTED` otherwise.

Qualifiers are reported alongside but never change the label: D2 autocorrelation-robust,
D2 same sign without the Big Five, D1 duration-robust, D1 uncertainty-robust and
quorum-robust. Overall expectation: `NOT_SUPPORTED` or `REVERSED` about three times in four
and `SUPPORTED` about one time in fifteen. Our reasons: the Cenozoic combines cooling with
high standardised richness; the end-Permian combines extreme heat with a richness crash;
and warming pulses coincide with several extinctions (the sister test examines that link
directly).

## How the result will be read

- `SUPPORTED`: on the best global temperature history, marine diversity is higher when
  Earth is warmer, which confirms Mayhew et al. on fully independent temperature data.
- `REVERSED`: on PhanDA, warm intervals and warming steps go with fewer genera. That is
  the opposite of the 2012 result and the direction of the 2008 one, but this time with
  sampling controlled.
- `MIXED` (for example, positive levels but negative changes) would mean the link depends
  on timescale, as Spiridonov & Lovejoy (2022) argue.
- `NOT_SUPPORTED` agrees with Balembois et al. (2025): global temperature barely sets
  long-term global richness.

Every number is reported in RESULTS.md, whatever the outcome.

## Known limits (stated before the result)

- Bins are stages of unequal length (about 1-20 Myr). Coverage standardisation
  removes much, but not all, of the duration effect; D3 checks this.
- PhanDA GMST is not the same as the tropical sea-surface temperatures Mayhew et al.
  used. The claim is tested in its global-mean form.
- PhanDA draws are treated as independent between slices. Its ensemble is correlated in
  time, so D4 overstates the noise (a conservative test).
- PBDB is queried live, so its content changes over time. File md5s are saved in
  provenance.json. Global richness is still shaped by where fossils were collected
  (Europe and North America dominate); coverage standardisation does not fix spatial
  sampling.
- Linear detrending is a registered choice. Other detrendings are not run, to avoid
  forking paths.
