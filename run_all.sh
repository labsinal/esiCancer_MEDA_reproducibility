#!/usr/bin/env bash
# Reproduce every table and figure of the revision from the data in this repository.
# Usage: bash run_all.sh   (from the repository root; Python environment as in requirements.txt)
set -euo pipefail
cd "$(dirname "$0")/code"
mkdir -p ../output/nmsc ../output/germline ../figures

echo "1/9 NMSC (Figure 1A, B): GBD 2023 incidence vs NCI UV exposure"
NMSC_DATA_DIR=../data/epidemiology/nmsc NMSC_OUT_DIR=../output/nmsc \
  jupyter nbconvert --to notebook --execute nmsc_analysis.ipynb --output-dir ../output/nmsc --output nmsc_analysis_executed.ipynb

echo "2/9 Lung cancer (Figure 1C): Halpern et al. 1993, Table 4"
python lung_halpern_reanalysis.py > ../output/lung_halpern_results.txt

echo "3/9 Relative risk by age, epidemiological data (Figure 1E)"
python age_trend_rr.py > ../output/fig1_RR_trend_results.txt

echo "4/9 Figure 1 numbers and figure"
python fig1_revised_numbers.py > ../output/fig1_revised_numbers.txt
python fig1_revised.py

echo "5/9 Simulations: occurrence, MEDA and proportional-hazards benchmark (Figure 3A-I, Table S3)"
python ph_benchmark.py > ../output/ph_benchmark_results.txt

echo "6/9 Simulations: hazard ratio by age and MEDA decomposition"
python sim_age_trend_rr.py > ../output/fig3_HR_trend_results.txt
python meda_decomposition.py > ../output/meda_decomposition_results.txt

echo "7/9 Inherited first hit (Figure 3J, K) and validation on artificial data (Figure S2)"
python germline_analysis.py > ../output/germline/germline_results.txt
python germline_analysis.py --synthetic > ../output/germline/germline_synthetic_validation.txt
python germline_vs_input.py > ../output/germline/germline_vs_input_results.txt
python power_check.py > ../output/germline/germline_power_results.txt

echo "8/9 Figures 2, 3, S1, S2 and S3"
python schematics.py
python fig2_BC.py
python fig3_revised.py
python figS2_synthetic.py
python figS3_endpoint.py

echo "9/9 Matched-input turnover, Figure S4 (only if the runs are in data/simulations/phase2_turnover)"
if ls ../data/simulations/phase2_turnover/ctor*_automatics_ancestralResults.csv >/dev/null 2>&1; then
  python phase2_analysis.py exp2 > ../output/turnover_results.txt
  python figS4_turnover.py
else
  echo "   no turnover runs yet: skipped"
fi
echo "Supplementary Tables S1-S4"
python supp_tables.py
echo "done: tables in output/, figures in figures/"
