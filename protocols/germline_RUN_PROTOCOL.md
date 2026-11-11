# Germline first-hit simulation: run protocol

Purpose (Reviewer comment 11; revision plan §1): test whether an **inherited first hit** produces what a uniform increase in mutation input does not, namely a relative risk that **falls with age** and a MEDA shift **beyond** the proportional-hazards (PH) expectation (negative residual). This is the BRCA-like pattern.

## Why the 6-gene multistage model

- In the 179-gene ESCC model a second hit in one specific gene is too rare (biallelic TP53 in ≤0.1% of runs), so a "same gene, second allele" design shows nothing.
- esiCancer has no homeostatic control, so any growth-promoting event given to *all* cells makes the whole tissue grow. A strong driver given as a germline event (e.g. TP53, CDKN2A deletion) would cross the 10× threshold within a few years. That is not tumorigenesis.
- In the 6-gene model (Fig 3C), gene A alone changes net growth by only +3×10⁻⁵ per generation (about +16% tissue size over 92 y), and every later hit in B–F interacts with A (division effect ×3, death effect ×1/3). A germline A (given to 100% of cells at generation 0) therefore means **one stage fewer**, the textbook multistage set-up (Armitage–Doll; Knudson).

## Files (unchanged from Fig 3C; confirmed by the authors 2026-10-16)

- esiTable: `esitables/esiTable_six_genes.csv` (genome size 750; A–F, ADD +0.000015 division / −0.000015 death, dominance 1; aging records present but not scheduled)
- Interaction table: `esitables/interaction_six_genes.csv` (all ordered pairs, DIVMOD 3, DEMOD 0.333, MAXDIVMOD 1000)

## Settings common to all runs (as `esicancer/default_config.png`, without aging)

| Setting | Value |
|---|---|
| Prob of division / death | 0.0100 / 0.0100 |
| Max divisions | 1000 |
| Number of cells (initial) | 1000 |
| Stop: generations / cells | 4800 / 10000 |
| Seeds | 1–1000, increment 1 |
| Tumor environment, asymmetric mitosis | off |
| Manual events | **none** except the germline event below. **No aging events** (Fig 3C had none: no extinctions) |
| Output files | for entire simulations; tick "Output sequence" (needed for driver-order analysis) |

## The runs

Events per division: the codes in the Fig 3C file names `6_genes_mut_002/003/004/006` mean **2, 3, 4 and 6 events per division**. This is inferred, not stated: the same zero-padded convention in `revision/constant_events_escc_aging/mut_rate_const_001…010` gives seed-for-seed identical tumours (before 60 y) to the 1, 3, 5 and 10 events/division runs, and the share of cells carrying gene A scales accordingly (27/45/72/87%). Non-carrier runs at 002, 003 and 004 already exist (Fig 3C). Rerun them anyway in the same session and software version, so both arms are generated identically.

**Number of runs, from a power check** (`power_check.py`, `output/germline_power_results.txt`). On synthetic data with the theoretical multistage effect, 150 replicates per cell:

| Non-carrier occurrence | Runs per arm | Power (trend CI excludes 1) | False positives, PH control |
|---|---|---|---|
| 5% ("002") | 1,000 | 23% | 7% |
| 5% | 3,000 | 67% | 7% |
| 20% ("003") | 1,000 | 68% | 4% |
| 20% | 3,000 | **99%** | 4% |

The estimators are unbiased: multistage trend 0.84–0.86 vs theory 0.85; PH trend 1.00; PH residual ≈ 0. Hence:

| Priority | File name (Simulation name) | Events per division | Seeds | Manual event |
|---|---|---|---|---|
| **1 (primary)** | `germ_003_noncarrier` | **3** | **1–3000** | none |
| **1 (primary)** | `germ_003_carrier` | **3** | **1–3000** | generation **0**, percentage **100**, `A - mutation` |
| 2 | `germ_004_noncarrier` | **4** | 1–1000 (3000 if time allows) | none |
| 2 | `germ_004_carrier` | **4** | same | gen 0, 100%, `A - mutation` |
| 3 (optional) | `germ_002_noncarrier` / `germ_002_carrier` | **2** | ≥ 5000 to be useful | none / as above |

Save all outputs in `data/simulations/germline/`. Run the **carrier at 003 with seeds 1–20 first** and do the checks below before launching the full set.

Caveat: power was computed for the classic multistage effect size. If esiCancer's germline effect is weaker, power will be lower. The first 003 results will show this, and the run count can be raised then.

**Checks after the first carrier run** (before running the rest):
1. In `germ_003_carrier_automatics_genesMutationResults.csv`, column `a monoAl` (or `a biAl`) should equal the final population in (nearly) every surviving run. That confirms the manual event reached all cells.
2. Tissue size in carrier runs without tumour should be only slightly larger than in non-carriers (expected about 1.2×). If carriers cross the threshold within the first years, the germline effect is acting as tissue-wide growth, and the design must be revisited.
3. Generation 0 is accepted by the GUI (confirmed), so the germline event is present from the start.

**Why 002–004.** At low baseline occurrence (002: 3.6%; 003: 21%), depletion of the at-risk pool is small, so any MEDA shift beyond the PH expectation is easy to see. 004 (61%) shows the high-occurrence regime.

## What the analysis will report (`germline_analysis.py`)

For each input level, carriers vs non-carriers:
- cumulative occurrence (Wilson CI) and MEDA (bootstrap CI);
- observed MEDA shift, PH-expected shift, and residual (bootstrap CI): the §1 decomposition;
- hazard ratio at age 60 and its change per decade (Poisson model of the cause-specific hazard), so the result can be placed on Fig 1G / Fig 3 next to smoking, BRCA1/2 and the mutation-input scenarios;
- age-specific hazard ratio and cumulative-incidence plots.

## The comparison that decides the result

The 6-gene model already shows some age dependence when only mutation input is raised. Non-carriers at "003" vs "004" (existing Fig 3C data): HR trend 0.79 per decade (0.63–1.00), residual −2.1 y (−3.7 to −0.6). So the germline hit must be judged **against a mutation-input increase in the same model**, not only against the flat line. The analysis does this automatically: whenever non-carrier runs exist at two levels, it adds the contrast "more mutation input" (non-carriers 003 → 004) next to "germline first hit" (non-carriers 003 → carriers 003).

The non-carrier occurrences differ (003: 21%, 004: 61%). If the carriers at 003 land near 61%, the two contrasts produce the **same rise in incidence by two different mechanisms**. That is the cleanest possible comparison, because then only the mechanism differs. If carriers land far from 61%, I will add a non-carrier run at an intermediate or higher input to match.

## Prediction, stated before running

Multistage theory: if non-carriers need *k* rate-limiting steps and carriers *k−1*, the hazards scale as t^(k−1) and t^(k−2). Their ratio then falls as 1/t, which is about **0.86 per decade around age 60**, independent of *k*. So the predictions are:
- **HR change per decade < 1** (around 0.85 if the model is close to classic multistage);
- **residual < 0** (carriers' MEDA lower than the PH expectation);
- **both more extreme than for the matched mutation-input increase** (see above);
- carriers sit **below** the flat line on panel G, the BRCA side.

The analysis was validated on synthetic data with a known answer (`germline_analysis.py --synthetic`). Multistage carriers (4 vs 5 stages): residual −5.9 y (−10.7 to −2.2), trend 0.72 (0.53–0.98). PH control: residual compatible with 0.

If the result is null (trend ≈ 1), that is also informative and reportable. In esiCancer's 6-gene model, a germline stage would then act as a uniform multiplier, and the BRCA-like shift would need another mechanism (e.g. an effect on cell turnover or on the early, non-expanding phase).

## Optional second design (only if the first is positive)

To test robustness outside the toy model, use the ESCC table with one added record `GERM` (dominance 0: no effect when monoallelic; strong biallelic effect, e.g. TP53's values; PROBEVENT raised enough that a second hit occurs in a realistic share of carrier runs). Carriers get `GERM - mutation` at 100% at generation 0. This needs a decision on PROBEVENT, so discuss before running.

---

## Results, first batch (2026-10-19; seeds 1–1000 per group)

Files: `germ_003_carrier_*`, `germ_003_non_carrier_*` (in `revision_analysis/germline/`). Comparator `germ_004_non_carrier_*` is a symlink to the Fig 3C run `6_genes_mut_004`. Scripts: `germline_analysis.py --dir .` (→ `output/germline_results_003.txt`, `figures/germline_germline_003.png`) and `germline_vs_input.py` (→ `output/germline_vs_input_results.txt`).

**Checks (all passed)**
- The new non-carrier run is identical, seed by seed, to the Fig 3C run at 3 events/division (1,000/1,000). This confirms the same software and settings, and justifies reusing the Fig 3C 004 run as the comparator.
- Gene A is present in 100% of final cells in all 1,000 carrier runs.
- No early tissue-wide growth: earliest carrier tumour at 38 y, none before 20 y, no extinctions.
- **Report this:** carriers' surviving tissue roughly doubles by 92 y (median 3,310 cells vs 1,606 in non-carriers), because A interacts with the B–F hits that accumulate in every cell. Part of the carrier effect is therefore field-wide expansion, not only clonal.

**Main results** (reference: non-carriers at 3 events/division; occurrence 21%, MEDA 84.8 y)

| | Germline first hit | More mutation input (3 → 4) |
|---|---|---|
| Occurrence | 66.3% | 61.4% |
| MEDA | 80.0 y | 81.2 y |
| Observed MEDA shift | −4.9 y | −3.6 y |
| PH-expected shift | −1.8 y | −1.5 y |
| Residual (95% CI) | **−3.1 (−4.5 to −1.5)** | −2.1 (−3.6 to −0.5) |
| HR at 60 y | 9.5 | 6.7 |
| HR change per decade (95% CI) | **0.69 (0.55–0.86)** | 0.79 (0.63–1.00) |

**Direct comparison (joint model / joint bootstrap; the shared baseline is accounted for)**
- HR trend, germline relative to input: **0.87 (0.76–0.99), p = 0.03**. The germline hit concentrates the extra risk at younger ages more than a matched increase in mutation supply.
- Residual difference, germline − input: −1.0 y (−2.4 to +0.5), bootstrap P(≥0) = 0.08. Same direction, not yet resolved.
- At similar occurrence, carriers' MEDA is 1.3 y lower (−2.6 to +0.2).

**Against the predictions stated before running:** trend < 1 ✔; residual < 0 ✔; more extreme than matched input ✔ for the trend (p = 0.03), same direction but not significant for the residual. Magnitude: 0.69 per decade lies between BRCA2 (0.72) and BRCA1 (0.54) in BCAC. Also, mutation input alone is not exactly proportional in this 6-gene model (0.79), unlike the 179-gene model (0.97–1.03). Report this.

## Next batch (pre-specified: completes the 3,000 seeds per arm in the original protocol)

These runs are not a second look for significance. The protocol fixed 3,000 seeds per arm before any results, and the final analysis will be the pooled 3,000. Same settings as before, seeds **1001–3000**:

| Simulation name | Events/division | Seeds | Manual event |
|---|---|---|---|
| `germ_003_carrier_s1001` | 3 | 1001–3000 | gen 0, 100%, `A - mutation` |
| `germ_003_non_carrier_s1001` | 3 | 1001–3000 | none |
| `germ_004_non_carrier_s1001` | 4 | 1001–3000 | none |

Save them in `revision_analysis/germline/`. Both scripts automatically pool files of the same group (`germ_<level>_<arm>_<anything>_automatics_…`) and stop if any seed appears twice. Expected precision with 3,000 per arm: CIs about √3 narrower (residual difference ±0.85 y instead of ±1.45 y).

---

## Final results (2026-10-19; 3,000 seeds per arm, pre-specified)

The two 3-events/division arms were re-run with seeds 1–3000 (one file each; the first batch files were replaced). Checks: no duplicate seeds; seeds 1–1000 reproduce the first batch exactly (663 and 210 tumours) and the Fig 3C run; gene A is in every final cell of all 3,000 carrier runs; earliest carrier tumour 38.1 y; no extinctions. The first-batch outputs are kept as `output/*_batch1.txt`.

**Primary result: carriers vs non-carriers (3,000 vs 3,000)**

| | Non-carriers | Carriers |
|---|---|---|
| Occurrence by 92 y | 21.5% | 67.1% |
| MEDA | 84.3 y | 79.6 y |

| Decomposition | Value |
|---|---|
| Observed MEDA shift | −4.7 y |
| PH-expected shift (depletion only) | −2.4 y |
| **Residual (age-dependent part)** | **−2.3 y (95% CI −3.5 to −1.3)** |
| HR at 60 y | 8.4 |
| **HR change per decade** | **0.74 (0.65–0.84)** |

Both predictions stated before running are confirmed: the hazard ratio falls with age, and MEDA shifts beyond the PH expectation. The estimates moved towards smaller effects from the 1,000-seed batch (trend 0.69 → 0.74; residual −3.1 → −2.3), as expected when samples grow, and both intervals now exclude the null clearly. The trend of 0.74 per decade lies close to BRCA2 (0.72) and short of BRCA1 (0.54) in BCAC.

**Secondary result: germline vs matched mutation input** (comparator: 4 events/division, **still 1,000 seeds**, Fig 3C run; the `germ_004_non_carrier_s1001` batch was not run)

| | Value |
|---|---|
| HR trend: germline 0.74 vs input 0.83; ratio | **0.89 (0.80–0.99), p = 0.03** |
| MEDA residual: germline −2.3 vs input −1.25; difference | −1.0 y (−2.2 to +0.2); bootstrap P(≥0) = 0.04, two-sided ≈ 0.09 |
| MEDA at similar occurrence (67% vs 61%) | carriers 1.6 y earlier (−2.6 to −0.4) |

The last line is not a clean test on its own, because carriers' occurrence is somewhat higher, and higher occurrence alone lowers MEDA through depletion. The residual difference is the proper comparison.

**What can be claimed.** *In the 6-gene multistage model, a germline first hit concentrates the extra risk at younger ages, so the hazard ratio declines by about 26% per decade and MEDA shifts 2.3 y beyond what the rise in occurrence alone predicts. A mutation-input increase that raises occurrence similarly gives a significantly flatter hazard-ratio trend (p = 0.03), and a smaller MEDA residual whose difference from the germline residual is in the same direction but not statistically resolved.* Also report: in this model, more mutation input is itself not perfectly proportional (trend 0.83), unlike the 179-gene model (0.97–1.03), and carriers' tissue roughly doubles by 92 y.

**To resolve the secondary comparison** (optional, about 2,000 runs): `germ_004_non_carrier_s1001`, 4 events/division, seeds 1001–3000, no manual event. With 3,000 comparator runs, the CI for the residual difference narrows by about a third. If you skip it, report the secondary result as above, without claiming a MEDA difference between mechanisms.
