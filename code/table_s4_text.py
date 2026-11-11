"""Text of Supplementary Table S4 (endpoints, denominators and calculations), as in the manuscript."""

COLS = ['',
 'NMSC incidence, high vs low UV exposure',
 'Lung cancer mortality, current vs never smokers',
 'Breast cancer, BRCA1/2 PTV carriers vs non-carriers',
 'esiCancer simulations']

ROWS = [('Figures', '1A, 1B, 1E', '1C, 1E', '1D, 1E', '2, 3 (3J, K: inherited first hit); S2–S5'),
 ('Data source',
  'Global Burden of Disease Study 2023 results (IHME): cause "Non-melanoma skin cancer", measure incidence, '
  'metrics number and rate by age and sex, and age-standardised rate (both sexes), year 2023. NCI '
  'county-level UV exposure dataset (column 2020–2024).',
  'Halpern MT, Gillespie BW, Warner KE (1993) J Natl Cancer Inst 85:457–464, Appendix Table 4 (American '
  'Cancer Society Cancer Prevention Study II): published aggregate person-years and lung cancer deaths by '
  'age band and smoking status, checked against the printed table. NLST data not used.',
  'Breast Cancer Association Consortium (2021) N Engl J Med 384:428–439, Supplementary Appendix: Table S12 '
  '(carrier cases and controls and odds ratios by age group) and Table S11 (case-only odds ratio per year of '
  'age). Population-based studies only. Published aggregate values; no individual-level data.',
  'Simulated with esiCancer (ref. 17); time-varying event schedules (Figure 3B, E, G, H) with a modified '
  'version of esiCancer (esicancer/increasing_events_per_division.patch). 1,000 runs (seeds 1–1000) per condition in Figure 3A–H; 3,000 '
  'runs per arm in the inherited first-hit experiment. Configuration of every condition in Table S3.'),
 ('Study design',
  'Ecological cross-sectional comparison across states: linear regression of each state-level outcome on '
  'mean state UV exposure.',
  'Prospective volunteer cohort, enrolled 1982, biennial follow-up through 1988 (6 years), 98.2% traced; '
  'about 900,000 of 1.2 million participants analysed.',
  'Case–control study (BCAC): 34 genes sequenced in 60,466 breast cancer cases and 53,461 controls from 44 '
  'studies; values used here come from the population-based studies (cases and controls not selected for '
  'family history).',
  'In silico stochastic evolutionary simulation (esiCancer); each run is one independent simulated '
  'individual.'),
 ('Unit of analysis',
  'One US mainland location: 48 states plus DC (n = 49); Alaska and Hawaii excluded.',
  'Age band × smoking status (sexes pooled in the source table).',
  'Age group at diagnosis × carrier status (carriers of protein-truncating germline variants in BRCA1 or '
  'BRCA2). Non-carriers enter only through the odds ratios.',
  'One simulation run (seed).'),
 ('Event definition',
  'Incident non-melanoma skin cancer (GBD cause "Non-melanoma skin cancer"), both sexes combined.',
  'Death from cancer of the trachea, bronchus or lung (ICD-9). Mortality endpoint, not diagnosis.',
  'Breast cancer case status (cases vs unaffected controls); age at diagnosis in four groups.',
  'Tumour (esiTumor): cell population reaches ≥10× its initial size (≥10,000 cells from 1,000). Other fates: '
  'extinction of the whole cell population (competing event) or neither by the end of the simulation '
  '(censored).'),
 ('Risk measure',
  'Age-standardised incidence rate, cases per 100,000 population (GBD 2023, both sexes; y-axis of Figure '
  '1A). Exposure: UV, mean of county values 2020–2024 (Wh/m²), unweighted per state.',
  'Age-specific lung cancer death rate (deaths / person-years) and rate ratio current ÷ never: age-adjusted '
  '21.1 (95% CI 19.0–23.6), Poisson model with log person-years offset; age-specific ratios 12.6–22.7.',
  'Odds ratio of carrier cases in each age group vs all controls (BCAC Table S12), interpreted as the '
  'incidence rate ratio vs non-carriers of the same age: BRCA1 46.3 (<40 y), 14.2, 6.52, 2.33 (≥60 y); BRCA2 '
  '18.7, 7.85, 5.14, 2.81.',
  'Cumulative occurrence by T = 92 y: proportion of runs reaching the tumour threshold, with Wilson 95% CI.'),
 ('Population at risk',
  'Resident population of each location, both sexes, all ages.',
  'Aged 40–80; excluded for incomplete data or age outside 40–80. Smoking status recorded at baseline only.',
  'Women in the population-based BCAC studies (cases with breast cancer and unaffected controls).',
  '1,000 initial cells per run; probability of division = probability of death = CTOR (0.01 by default; 0.02 '
  'and 0.2 in Figure 3G, H); maximum 1,000 divisions per cell; events per division as given for each '
  'condition (Table S3).'),
 ('Denominator',
  'State population (modelled GBD estimate). For the age-specific data, population implied as Number ÷ Rate '
  '× 100,000.',
  'Person-years of exposure (published).',
  'None (case–control design). Odds ratios compare the frequency of carriers among cases and among controls '
  'and do not depend on the sampling fractions.',
  'Total number of runs per condition (1,000; 3,000 per arm in the inherited first-hit experiment).'),
 ('Observation period',
  'GBD estimates for 2023. UV exposure: mean of county values for 2020–2024.',
  'Enrolment 1982; biennial follow-up through 1988 (6-year follow-up window).',
  'Not applicable (case–control sampling; no follow-up).',
  '4,800 generations ≈ 92 years (1 generation = 1 week; not clinically calibrated). Endpoint sensitivity: '
  'follow-up ending at T = 60 and 75 years (Figure S3).'),
 ('Age range / bins',
  '19 GBD age classes: <1, 2–4, 5–14, 15–19, 20–24 … 70–74, 75–84, 85–89, 90–94, 95+. Open-ended 95+ class '
  'bounded at 100 y (sensitivity: 105, 110 y).',
  '40–43, 44–48, 49–53, 54–58, 59–63, 64–68, 69–73, 74–80 (top band closed at 81 y for interpolation).',
  '<40, 40–49, 50–59, ≥60 (four groups; the first and last are open-ended).',
  '0–92 simulated years; 10-year bands for the hazard-ratio model.'),
 ('Standardisation',
  'Incidence is age-standardised (GBD standard population). Median age of cases (Figure 1B) is not '
  "standardised and depends on each state's age structure; the rate-profile median is free of age structure.",
  'Rate ratio age-adjusted by Poisson regression. Distribution of ages at death compared at a common age '
  "structure: each group's age-specific death rates applied to the pooled person-years of both groups.",
  'None by us; the BCAC odds ratios are adjusted for country/strata.',
  'Not applicable (all runs start from identical populations).'),
 ('Age endpoint',
  'Age at diagnosis. Primary: median age of diagnosed cases (mean 71.5 y; range 70.5–73.6). Second endpoint: '
  'median of the cumulative age-specific rate profile (mean 84.8 y; range 83.6–85.6).',
  'Age at death (lung cancer mortality).',
  'Age at diagnosis of carrier cases (grouped).',
  'Simulated detection age: generation at which the threshold is crossed ÷ 52.143 (weeks per year).'),
 ('Median age calculation',
  'Grouped-data median using the median class and class-specific widths, uniform within class. Cases: GBD '
  'Number summed over sexes. Profile: both-sex rate × class width.',
  'Median age at death at the common age structure, linear interpolation within age bands: never smokers '
  '67.5 y, current smokers 67.6 y; difference +0.1 y (95% CI −1.2 to +1.3; parametric Poisson bootstrap). '
  'Deaths as observed in the cohort (not standardised): 68.5 vs 64.9 y.',
  'Median of carrier cases only, linear interpolation within the age group containing the median (40–49 y '
  'for both genes), uniform within group: BRCA1 44.5 y (95% CI 43.5–45.5), BRCA2 49.4 y (48.4–50.6); '
  'multinomial bootstrap. Unchanged (<0.1 y) for open-group bounds of 15–25 y and 80–100 y. No non-carrier '
  'median (not published by BCAC).',
  'MEDA: median simulated detection age among runs with a tumour before T, median(τ | τ ≤ T), a conditional '
  'median; bootstrap 95% CI (2,000 resamples of runs).'),
 ('Statistical test',
  'Linear regression on UV (n = 49), 95% slope CI and Spearman. Incidence: R² = 0.209, p = 0.00096. Median '
  'age of cases: R² = 0.004, p = 0.68 (slope +0.06 y per 1,000 Wh/m², 95% CI −0.25 to +0.38). Median of rate '
  'profile: R² = 0.096, p = 0.031 (p = 0.053–0.088 if the 95+ bound is 105–110 y). High (>5,000, n = 8) vs '
  'low (<4,000, n = 11) UV: ratio of incidence 1.164 (1.06–1.27); of median age of cases 1.0015 (0.99–1.01); '
  'of rate profile median 1.0056 (1.001–1.010).',
  'Poisson GLM: deaths ~ smoking + age band, offset log(person-years). Smoking × age interaction, '
  'likelihood-ratio test: χ² = 3.86, 7 df, p = 0.80 (tests whether the rate ratio varies with age; not a '
  'test of equal median age).',
  'No tests by us: case–control counts are not incidence. Age dependence from the BCAC case-only odds ratio '
  'per year of age (BCAC Table S11): BRCA1 0.941 (0.94–0.95), BRCA2 0.968 (0.96–0.97).',
  'Wilson CI (occurrence); bootstrap CIs (MEDA and between-condition differences); proportional-hazards '
  'benchmark: residual = observed − PH-predicted MEDA, bootstrap CI; Poisson model of the hazard ratio by '
  '10-year age band (change per decade).'),
 ('Change in relative risk per decade of age',
  '1.04 (descriptive): weighted fit of the high- vs low-UV incidence rate ratios by age, which rise from '
  '1.07 at 30–34 y to 1.36 at ≥95 y; modelled GBD values, no CI.',
  '1.10 (95% CI 0.98–1.25): Poisson model with a smoking × age (continuous) term.',
  'BRCA1 0.54 (0.54–0.60); BRCA2 0.72 (0.66–0.74): BCAC case-only odds ratio per year raised to the 10th '
  'power.',
  'Hazard ratio vs the lowest input of each series: 0.97–1.14 per decade for changes in mutation input '
  '(multistage series imprecise, 0.77–1.24); 0.74 (0.65–0.84) for the inherited first hit vs non-carriers.'),
 ('Comparison with the proportional-hazards (PH) expectation',
  'Not computed (see change in relative risk per decade).',
  'Not computed (see change in relative risk per decade).',
  'Not computed (requires absolute incidence; see change in relative risk per decade).',
  'Residual = observed − PH-predicted MEDA. All 21 comparisons vs the lowest input of each series include 0 '
  '(e.g. 1→10 events per division: −0.9 y, 95% CI −5.3 to +3.0). Multistage, adjacent inputs: −2.1 y (−3.6 '
  'to −0.5) and −4.3 y (−5.8 to −2.2). Inherited first hit: −2.3 y (−3.5 to −1.3).'),
 ('Data and code',
  'IHME_data_USA_states.csv; IHME-GBD_2023_DATA-05dbcd0f-1.csv; uv-county.xlsx; NMSC notebook '
  '(nmsc_analysis.ipynb); provenance.json.',
  'lung_halpern_table4.csv; lung_halpern_reanalysis.py.',
  'Values transcribed from BCAC Tables S11–S12 (bcac2021_tableS11_case_only.csv, '
  'bcac2021_tableS12_population_based.csv; fig1_revised_numbers.py).',
  'Per-run outputs (*_ancestralResults.csv); ph_benchmark.py, meda_decomposition.py, sim_age_trend_rr.py, '
  'germline_analysis.py, germline_vs_input.py (code/ of this package).')]
