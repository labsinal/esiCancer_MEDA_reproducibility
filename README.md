# Cancer incidence and median age at detection – reproducibility package

Code, input data, simulation outputs and configurations for Kraemer-Mattos et al., *The relation between cancer median age of detection and incidence is differently impacted by mutation rate* (revised version). Every number, table and figure of the revision can be regenerated from this folder with one command.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
bash run_all.sh
```

Runtime about 10 minutes on a laptop. Tables and text results go to `output/`, figures to `figures/`. Bootstrap analyses use fixed random seeds, so results are identical between runs.

## Folder structure

| Folder | Content |
|---|---|
| `code/` | All analysis and figure scripts; `paths.py` defines the folder layout |
| `data/epidemiology/` | Published inputs: NMSC (GBD 2023, NCI UV), lung cancer (Halpern et al. 1993, Table 4), BRCA (BCAC 2021, Tables S11, S12) |
| `data/simulations/` | esiCancer outputs, one folder per simulated series (one row per run in `*_ancestralResults.csv`) |
| `esitables/` | esiTables and interaction table used in the simulations |
| `esicancer/` | Simulation settings and model rules (`CONFIGURATIONS.md`, `default_config.png`) and the modification used for the increasing schedules (`increasing_events_per_division.patch`) |
| `protocols/` | Run protocols written before the new simulations (inherited first hit, matched-input turnover) |
| `output/`, `figures/` | Created by `run_all.sh` |

## Manuscript item → script → output

| Manuscript item | Script | Output |
|---|---|---|
| Fig 1A, B (NMSC vs UV) | `code/nmsc_analysis.ipynb` | `output/nmsc/state_level_results.csv`, `regression_summary.csv`, `high_vs_low_uv_ratios.csv`, `provenance.json` |
| Fig 1C (lung cancer) | `code/lung_halpern_reanalysis.py` | `output/lung_halpern_results.txt` |
| Fig 1E, change in relative risk per decade | `code/age_trend_rr.py` | `output/fig1_RR_trend_summary.csv`, `fig1_age_specific_RR.csv` |
| Fig 1 (all panels) | `code/fig1_revised_numbers.py`, `code/fig1_revised.py` | `output/fig1_revised_numbers.json`, `figures/Figure1_revised.png/.pdf` |
| Fig 3A–I, Table S3 (occurrence, MEDA, PH benchmark) | `code/ph_benchmark.py` | `output/conditions_occurrence_meda.csv`, `output/pairs_ph_benchmark.csv` (T = 60, 75, 92 y) |
| Hazard ratio by age (Results), MEDA decomposition | `code/sim_age_trend_rr.py`, `code/meda_decomposition.py` | `output/fig3_HR_trend_summary.csv`, `fig3_age_specific_HR.csv`, `meda_decomposition.csv` |
| Fig 3J, K (inherited first hit) | `code/germline_analysis.py`, `code/germline_vs_input.py` | `output/germline/germline_summary.csv`, `germline_vs_input.csv` |
| Validation on synthetic data, power (Fig S2) | `code/germline_analysis.py --synthetic`, `code/power_check.py`, `code/figS2_synthetic.py` | `output/germline/germline_SYNTHETIC_validation.csv`, `germline_power_check.csv` |
| Fig 3 (all panels) | `code/fig3_revised.py` | `figures/Figure3_revised.png/.pdf` |
| Matched-input turnover (Comment 9; Fig S4) | `code/phase2_analysis.py exp2`, `code/figS4_turnover.py` | `output/turnover_results.txt`, `output/turnover/`, `figures/FigureS4_turnover.png/.pdf` |
| Fig 2A, Fig S1 (schematics) | `code/schematics.py` | `figures/Figure2A_schematic.png/.pdf`, `figures/FigureS1_guide.png/.pdf` |
| Fig 2B, C | `code/fig2_BC.py` | `figures/Figure2_BC.png/.pdf` |
| Sensitivity to the end of follow-up (Fig S3) | `code/figS3_endpoint.py` | `figures/FigureS3_endpoint.png/.pdf` |
| Supplementary Tables S1–S4 | `code/supp_tables.py` | `tables/TableS1_parameters.xlsx` … `TableS4_endpoints.xlsx` |

## Definitions used in the analysis

- **Run fate:** tumour (population ≥ 10,000 cells, tenfold the initial population), extinction (population 0) or neither (censored at T = 92 years). Age = generation / 52.143.
- **Cumulative occurrence:** proportion of runs with a tumour by T; 95% Wilson interval.
- **MEDA:** median detection age among runs with a tumour by T, median(τ | τ ≤ T); 95% bootstrap interval (2,000 resamples of runs).
- **Proportional-hazards (PH) benchmark:** for a pair of conditions with occurrences p₀ and p₁, a = ln(1 − p₁)/ln(1 − p₀); the PH-predicted MEDA is the age at which the baseline cumulative occurrence reaches q* = 1 − (1 − p₁/2)^(1/a). Residual = observed − predicted MEDA, with a bootstrap 95% interval over the runs of both conditions.
- **Hazard ratio by age:** Poisson regression of tumour counts in 10-year bands with log person-time at risk as offset; change per decade from the age × condition term.

## Simulation output files

`*_automatics_ancestralResults.csv` (semicolon-separated), one row per run: `Seed`, `Generation` (end of the run), `NumberOfDivisions`, `Population` (≥ 10,000 = tumour, 0 = extinction), followed by per-ancestor counts. `*_genesMutationResults.csv` and `*_eventsMutationResults.csv`: mutated cells per gene and per event at the end of each run. The per-cell sequence files (`*_sequenceEachCell.csv`, several GB per condition) are not included because of their size; they are available from the authors on request.

## esiCancer version

Simulations used the esiCancer Linux build (Minussi et al., 2019). The increasing schedules of events per division used a modified build; its source is at https://github.com/labsinal/esiCancer_dynamic_events, and the change to `cellsystem.cpp` is in `esicancer/increasing_events_per_division.patch`, and `esicancer/CONFIGURATIONS.md` describes the settings of every condition and the model rules checked in the source code.

## Sources of the epidemiological data

- NMSC: Global Burden of Disease 2023 results, Institute for Health Metrics and Evaluation (GHDx results tool); cause "Non-melanoma skin cancer", incidence, number and rate by age and sex and age-standardised rate, 2023, 48 contiguous US states and DC (`data/epidemiology/nmsc/IHME_*.csv`). UV: National Cancer Institute, county ultraviolet exposure (`uv-county.xlsx`), unweighted mean of county values per state.
- Lung cancer: Halpern MT, Gillespie BW, Warner KE. J Natl Cancer Inst 1993;85:457–64, Table 4 (deaths and person-years by age band, current and never smokers), checked against the printed table.
- BRCA: Breast Cancer Association Consortium. N Engl J Med 2021;384:428–39, Supplementary Appendix Table S12 (population-based studies: carrier cases and odds ratios by age group) and Table S11 (case-only odds ratio per year of age).

## License

Code: MIT License (see `LICENSE`). Simulation outputs and tables: CC BY 4.0. Third-party epidemiological data in `data/epidemiology/` remain under the terms of their sources (IHME GBD terms of use; NCI; Halpern et al., 1993; Breast Cancer Association Consortium, 2021).
