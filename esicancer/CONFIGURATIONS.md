# esiCancer configurations

Software: esiCancer (Minussi et al., Cancer Res 2019; www.ufrgs.br/labsinal/esiCancer/), Linux build (Qt; source files `cell.cpp`, `cellsystem.cpp`, `mutationtable.cpp`).
The increasing schedules of events per division (Figure 3B, E, G, H; Figure 2C) used a modified build in which `cellsystem.cpp` sets the events per division from the generation; the change is in `increasing_events_per_division.patch` (see below).

Settings shown in the GUI are in `default_config.png`. Unless stated otherwise, every condition used:

| Setting | Value |
|---|---|
| Initial number of cells | 1,000 |
| Probability of division / death | 0.01 / 0.01 (CTOR 0.01) |
| Maximum divisions per cell | 1,000 |
| Stop | 4,800 generations (92 years, 1 generation = 1 week) or 10,000 cells (esiTumor) |
| Seeds | 1–1000, increment 1 (one run per seed) |
| Tumour environment, asymmetric mitosis | off |
| Output | for entire simulations (`*_automatics_ancestralResults.csv`, `*_eventsMutationResults.csv`, `*_genesMutationResults.csv`) |

esiTables are in `../esitables/`. The esiTable of each series was checked against the gene list written in its `*_genesMutationResults.csv`: every ESCC series has the 177 genes (179 gene–event records plus 4 aging records) of `esiTable_escc.csv`, and the six-gene series have genes A–F.

## Conditions in Figure 3 (folders in `../data/simulations/`)

| Panel | Folder | Files (prefix) | esiTable | Interaction table | Events per division | Aging events | CTOR |
|---|---|---|---|---|---|---|---|
| A | `revision/constant_events_escc_no_aging` | `constant_1`, `constant_3`, `constant_5`, `constant_10` | `esiTable_escc.csv` | none | 1, 3, 5, 10 (constant) | none | 0.01 |
| B | `revision/linear_events_escc_no_aging` | `linear_5`, `linear_10`, `linear_14`, `linear_20` | `esiTable_escc.csv` | none | increasing: 1 + floor(age / step), steps of 20, 10, 7 and 5 years (modified esiCancer) | none | 0.01 |
| C | `revision/6_genes_interaction_mut_rate` | `6_genes_mut_002`, `_003`, `_004`, `_006` | `esiTable_six_genes.csv` | `interaction_six_genes.csv` | 2, 3, 4, 6 (constant) | none | 0.01 |
| D | `revision/constant_events_escc_aging` | `mut_rate_const_001`, `_003`, `_005`, `_010` | `esiTable_escc.csv` | none | 1, 3, 5, 10 (constant) | 60, 70, 80, 90 y | 0.01 |
| E | `revision/mut_rate_PA_mut` | `mut_rate_PA5_mut_2`, `PA10_mut`, `PA14_mut`, `PA20_mut` | `esiTable_escc.csv` | none | increasing, as B (modified esiCancer) | 60, 70, 80, 90 y | 0.01 |
| G | `revision/ctor_002` | `ctor_002_5`, `_10`, `_14`, `_20` | `esiTable_escc.csv` | none | increasing, as B (modified esiCancer) | 60, 70, 80, 90 y | 0.02 (division = death) |
| H | `revision/ctor_02` | `ctor_02_5`, `_10`, `_14`, `_20` | `esiTable_escc.csv` | none | increasing, as B (modified esiCancer) | 60, 70, 80, 90 y | 0.2 (division = death) |

**Increasing schedule.** Events per division = 1 + floor(generation / step), with the step in generations (522 for 10 years in the archived source). The level in the file names is nominal (100 / step): at 92 years the schedules reach 5, 10, 14 and 19 events per division. The value typed in the interface is not used by this build.

**CTOR.** Division and death probabilities were equal. In the outputs, the divisions per cell per generation are about 0.01, 0.02 and 0.2 for the three CTOR levels.

**Aging events** (records `aging_1`…`aging_4` of the esiTable, applied as manual events to 100% of cells): probability of division ×0.9 and of death ×1.1 at each event. The drops in net growth at 60 and 70 years are visible in the single runs of Figure 2B (`revision/single_runs_aging/`); after 80 years too few cells remain to see the later events in individual runs.

**Figure 2B** (`revision/single_runs_aging/`): 14 single runs (seeds 1–14), `esiTable_escc.csv`, 20 events per division, aging events at 60, 70, 80 and 90 years, CTOR 0.01. The file `single_runs_aging_population.csv` holds the population of every run at every generation (extracted from the esiCancer ancestral output); the parameter files of the runs are included.

## Rules of esiCancer relevant to the reviewer's questions (checked in the source code)

- **Dominance.** A first event on one allele adds effect × DOMINANCE; the same event on the second allele adds effect × (1 − DOMINANCE), so both alleles together give the full effect (`MutationEvent::computeNewProlRate`, `computeNewDeathRate`).
- **ADD / MULT.** ADD adds the effect to the cell's probability; MULT multiplies the probability by the effect.
- **Bounds.** Probabilities of division and death are capped at 1, and when their sum exceeds 1 both are rescaled so that it equals 1 (`CellSystem::addNewCell`, `Cell::reproduce`). There is no lower clip: a negative probability acts as 0 (the cell never divides, or never dies), and later effects are added to the stored value.
- **Interactions.** For each partner event already present in the cell, the effect of a new event on division is multiplied by DIVMOD and its effect on death by DEMOD (`computeProlSinergyModifier`, `computeDeathSinergyModifier`). In the six-gene table DIVMOD = 3 and DEMOD = 1/3, so the division effect is multiplied by 3 and the death effect divided by 3 for each partner already mutated.
- **Mutation-rate modifier.** esiCancer reads it only from columns named `mutfunc` and `mut` (`MutationTable`, file loader). The esiTables name them `mutationFUNC` and `mutation`, so the +1 of TP53 and ATM was never read and had no effect in any simulation. In the modified build the events per division are also set by the schedule.

## Inherited first hit (Figure 3J, K) – `../data/simulations/germline`

Protocol: `../protocols/germline_RUN_PROTOCOL.md`. Six-gene model as panel C (`esiTable_six_genes.csv`, `interaction_six_genes.csv`), 3 events per division, no aging.

| Arm | File prefix | Manual event | Seeds |
|---|---|---|---|
| Carriers | `germ_003_carrier` | `A - mutation` in 100% of cells at generation 0 | 1–3000 |
| Non-carriers | `germ_003_non_carrier` | none | 1–3000 |
| Comparator, 4 events per division | `germ_004_non_carrier` (= panel C `6_genes_mut_004`) | none | 1–1000 |

## Matched-input turnover (Comment 9; Figure S4) – `../data/simulations/phase2_turnover`

Protocol: `../protocols/turnover_RUN_PROTOCOL.md`. `esiTable_escc.csv` with the genome size multiplied by 2 (`esiTable_escc_genome_x2.csv`, CTOR 0.02) or 20 (`esiTable_escc_genome_x20.csv`, CTOR 0.2); maximum divisions 10,000 (it never binds: lineages divided fewer than about 500 times before extinction or tumor); 1,000 runs per condition (seeds 1–1000). Files: `ctor002_e3_g2`, `ctor002_e10_g2`, `ctor_02_e3_g20`, `ctor_02_e10_g20`. The genome-scaled esiTables were used instead of fractional events, which esiCancer does not accept.
