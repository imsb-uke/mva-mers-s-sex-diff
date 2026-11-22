#!/usr/bin/env bash
# Run the revision analysis from the input tables in ../data.
# Step 14 reads the WGCNA networks written by ../original_analysis; run that first.
#   PYTHON   Python interpreter (default: python)
#   RSCRIPT  Rscript (default: Rscript)
#   THREADS  cores for DESeq2, dream and the WGCNA reconstructions (default 16)
set -euo pipefail
cd "$(dirname "$0")/scripts"
PY=${PYTHON:-python}
RS=${RSCRIPT:-Rscript}
mkdir -p ../logs

run() { echo "== $1"; "${@:2}" > "../logs/${1%.*}.log" 2>&1; }
run 01_prepare_data.py             $PY 01_prepare_data.py
run 02_vst.R                       $RS 02_vst.R
run 03_population.py               $PY 03_population.py
run 04_gene_level_deseq2.R         $RS 04_gene_level_deseq2.R
run 05_gene_level_mixed_model.R    $RS 05_gene_level_mixed_model.R
run 06_signature_scores.py         $PY 06_signature_scores.py
run 07_signature_mixed_model.R     $RS 07_signature_mixed_model.R
run 08_gene_level_summary.py       $PY 08_gene_level_summary.py
run 09_variance_partitioning.py    $PY 09_variance_partitioning.py
run 10_variance_summary_and_cells.py $PY 10_variance_summary_and_cells.py
run 11_irf7_isg15.py               $PY 11_irf7_isg15.py
run 12_irf7_isg15_baseline.R       $RS 12_irf7_isg15_baseline.R
run 13_deg_counts_and_permutation.R $RS 13_deg_counts_and_permutation.R
run 14_wgcna_robustness.R          $RS 14_wgcna_robustness.R
run 15_wgcna_foxp3_summary.py      $PY 15_wgcna_foxp3_summary.py
run 16_antibody.R                  $RS 16_antibody.R
run 17_figures.py                  $PY 17_figures.py
echo "done"
