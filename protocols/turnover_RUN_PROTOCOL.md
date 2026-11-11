# Phase 2 controls: run protocol

Updated 2026-10-19; Part B redesigned 2026-10-28 for integer events per division. The no-mutator control (old Experiment 1) was **dropped**: the authors report that the esiTable `mutation` column (mutation-rate modifier, used by TP53 and ATM) has no effect in the esiCancer version used. A short verification run (Part A) documents this for the reviewer (Comment 8). Part B is the matched-input turnover experiment (Comment 9).

Save every output in `revision_analysis/phase2_controls/runs/`, with the names given below; the analysis script finds runs by name.

---

## Part A: verification that the mutation-rate field is inert (2 short runs, about 5 minutes)

**Outcome (2026-11-10).** Not run: the esiCancer source settles it. The file loader reads the mutation-rate modifier only from columns named `mut` and `mutfunc`, while the esiTables use `mutation` and `mutationFUNC`, so the field is never read (see `../esicancer/CONFIGURATIONS.md`).

**Why.** The reviewer asked whether TP53/ATM change the mutation rate. Answering "the field does nothing" needs evidence. The 6-gene germline set-up is the most sensitive place to test it, because gene A is in every cell from generation 0 and genes are hit often (genome size 750). If the field worked, +10 events per division would multiply carriers' mutation rate by 4.3 (3 → 13 events) and the runs would diverge almost immediately.

**Run** exactly as `germ_003_carrier` (`germline/RUN_PROTOCOL.md`: 3 events per division, manual event at generation 0, 100% of cells, `A - mutation`, no aging, interaction table `esitables/interaction_six_genes.csv`), except:

| Simulation name | esiTable | Seeds |
|---|---|---|
| `verify_mutator10` | `phase2_controls/esiTable_6genes_A_mutator10_TEST.csv` (= `esitables/esiTable_six_genes.csv` with gene A `mutation` ADD **10** instead of 0; one cell changed) | 1–50 |

**Reading.** `python phase2_analysis.py verify` compares it with seeds 1–50 of `germline/germ_003_carrier`.
- **Identical** generation and population in all 50 seeds → the field is inert. Report this in the response and set the TP53/ATM `mutation` values to 0 in the published Table S2, so the table describes what was actually simulated.
- **Different** → the field works after all. Then TP53/ATM did raise the mutation rate in the Fig 3 runs, and the no-mutator control has to come back. Tell me before going on.

---

## Part B: matched-input turnover (Comment 9)

**Question.** With a fixed number of events per division, higher turnover (CTOR) also means more mutations per unit time. Does turnover still change occurrence and MEDA once mutation input per unit time is held equal?

**Design choices**
- **179-gene ESCC table, unchanged:** `esitables/esiTable_escc.csv` (identical to Table S2), as in Fig 3A.
- **Constant events, no aging:** isolates turnover from the aging-driven extinction that confounded Fig 3E/G/H, and needs no modified code.
- **CTOR = probability of division = probability of death** (as in Fig 3E/G/H; the divisions per cell per generation in their outputs are 0.01, 0.02 and 0.2).
- **Max divisions = 10000 in the new runs.** At CTOR 0.2 a lineage divides about 960 times in 92 y, so the default limit of 1000 would start to bind only at high turnover. At CTOR 0.01 (about 48 divisions per lineage) it never binds, so the existing Fig 3A runs (limit 1000) serve as the CTOR 0.01 arms.
- **Matching rule:** mutation input per unit time = CTOR × events per division ÷ genome factor. Matched at the two Fig 3A levels, 3 and 10 events per division at CTOR 0.01.
- **Integer events only (updated 2026-10-28).** esiCancer accepts only integer events per division, so the input is lowered with the **genome size** in the esiTable instead. Each event lands on a random genome position, so a gene's chance of being hit per event is (its weight) ÷ (genome size). Doubling the genome size halves every gene's mutation probability, which is the same as halving the events per division. Only the genome-size cell differs from `esiTable_escc.csv` (16,647,721):
  - `phase2_controls/esiTable_escc_genome_x2.csv`: genome size 33,295,442 (×2), for CTOR 0.02;
  - `phase2_controls/esiTable_escc_genome_x20.csv`: genome size 332,954,420 (×20), for CTOR 0.2.

**Common settings:** initial cells 1000; stop at 4800 generations or 10000 cells; seeds 1–1000; no manual events; no interaction table (as Fig 3A); tumour environment and asymmetric mitosis off.

**Matched arms (priority 1).** Same events per division as Fig 3A; the larger genome compensates for the extra divisions.

| Matched input (events per cell per generation) | CTOR 0.01 (existing) | CTOR 0.02 (new, genome ×2) | CTOR 0.2 (new, genome ×20) |
|---|---|---|---|
| 0.03 | 3 ev/div: Fig 3A `constant_3` | 3 ev/div → `ctor002_e3_g2` | 3 ev/div → `ctor02_e3_g20` |
| 0.10 | 10 ev/div: Fig 3A `constant_10` | 10 ev/div → `ctor002_e10_g2` | 10 ev/div → `ctor02_e10_g20` |

**Unmatched arms (priority 2, recommended).** Same events per division **and the original table** (`esiTable_escc.csv`), so input per unit time rises with turnover. This is the confounded design of the submitted Fig 3, but without aging. Each unmatched arm differs from its matched arm only in genome size, so together they split the total turnover effect into a mutation-supply part and a population-dynamics part.

| CTOR 0.02 | CTOR 0.2 |
|---|---|
| 3 ev/div → `ctor002_e3` | 3 → `ctor02_e3` |
| 10 → `ctor002_e10` | 10 → `ctor02_e10` |

### Pilot first: does genome size scale the mutation input?

All at CTOR 0.01 (the Fig 3A settings), seeds 1–50, no aging:
1. `pilot_e10`: 10 events per division, `esiTable_escc.csv`.
2. `pilot_e10_g2`: 10 events per division, `esiTable_escc_genome_x2.csv`.
3. `pilot_e5`: 5 events per division, `esiTable_escc.csv`.

`python phase2_analysis.py pilot` compares the mutated gene-events per cell. If genome ×2 and 5 events both give about half of 10 events, the matched design above stands. If genome ×2 gives about the same as 10 events, genome size does not set the hit probability: use **Plan B**.

**Plan B (only if the pilot fails):** match at 0.2 and 0.4 events per cell per generation with the original table: CTOR 0.01 × 20 / 40, CTOR 0.02 × 10 / 20, CTOR 0.2 × 1 / 2 (all new runs, max divisions 10000). The cost is high occurrence at CTOR 0.01, where depletion dominates MEDA.

### Expected result (stated before running)

If the turnover effect in the submitted Fig 3 was mainly mutation supply, the matched arms will show similar occurrence and MEDA across CTOR, and the unmatched arms will carry the difference. If population dynamics matters on its own (as our 2024 Carcinogenesis paper argues), the matched arms will still differ: faster clonal selection at high turnover gives higher occurrence and lower MEDA. Also watch for:
- **Drift extinction:** with no density control, high turnover makes whole-tissue extinction by drift much more frequent (in Fig 3H, 34–46% of runs). It is analysed as a competing event.
- **Realised input:** the analysis reports total events / total years at risk, overall and in tumour-free runs. In the existing CTOR 0.01 runs, tumour-free runs get 96–99% of the expected input, which validates the matching rule.

---

## Run list and order

| Order | Runs | Size |
|---|---|---|
| 1 | ~~Part A: `verify_mutator10`~~ (skipped 2026-10-28: stated from earlier tests) | – |
| 2 | Pilot: `pilot_e10`, `pilot_e10_g2`, `pilot_e5` | 3 × 50 seeds |
| 3 | Part B matched: `ctor002_e3_g2`, `ctor002_e10_g2`, `ctor02_e3_g20`, `ctor02_e10_g20` | 4 × 1,000 |
| 4 | Part B unmatched: `ctor002_e3`, `ctor002_e10`, `ctor02_e3`, `ctor02_e10` | 4 × 1,000 |

## Analysis

In this folder: `python phase2_analysis.py verify | pilot | exp2`.
