#!/usr/bin/env bash
# Runs the whole-blood RNA-seq analysis in order. Run from original_analysis/.
# PYTHON and RSCRIPT can be set to the interpreters of the environments in ../environment/.
# Steps 06 (GO terms) and 07 query the Enrichr web service and need internet access.
set -euo pipefail
cd "$(dirname "$0")"
PYTHON=${PYTHON:-python}
RSCRIPT=${RSCRIPT:-Rscript}

$PYTHON 01_sample_set.py
$RSCRIPT 02_vst.R
$PYTHON 03_pca.py
$RSCRIPT 04_deseq2_per_sex.R
$PYTHON 05_deg_tables_and_counts.py
$PYTHON 06_de_gene_clustering.py
$PYTHON 07_enrichr_per_comparison.py
$PYTHON 08_wgcna_input.py
$RSCRIPT 09_wgcna.R
$PYTHON 10_wgcna_modules.py
$RSCRIPT 11_wgcna_module_enrichment.R
$PYTHON 12_figure6_tables.py
