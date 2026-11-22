# Original analysis

Whole-blood RNA-seq analyses as reported in the paper. Run `bash run_all.sh` from this folder;
outputs are written to `results/` and `figures/`.

| Step | Script | Content |
|---|---|---|
| 01 | `01_sample_set.py` | Sample sets: analysis (138), PCA (174), WGCNA (160) |
| 02 | `02_vst.R` | Variance-stabilising transformation (DESeq2) |
| 03 | `03_pca.py` | PCA on the 2,000 most variable genes |
| 04 | `04_deseq2_per_sex.R` | DESeq2 D1 vs D0 per arm, vaccination and sex, design `~ Timepoint` |
| 05 | `05_deg_tables_and_counts.py` | DEG tables (adjusted p ≤ 0.05, fold change ≥ 1.5) and counts |
| 06 | `06_de_gene_clustering.py` | Clustering of the DEGs, heatmap, Enrichr GO terms per cluster |
| 07 | `07_enrichr_per_comparison.py` | Enrichr GO terms of each SI comparison, shared and sex-specific terms |
| 08 | `08_wgcna_input.py` | WGCNA input (10,000 genes, 80 samples per arm) |
| 09 | `09_wgcna.R` | WGCNA networks and module-trait correlations |
| 10 | `10_wgcna_modules.py` | Module codes, gene lists, module-trait heatmaps |
| 11 | `11_wgcna_module_enrichment.R` | GO enrichment of the modules (clusterProfiler) |
| 12 | `12_figure6_tables.py` | Tables for Figure 6 (TPM, per-participant fold changes) |

Steps 06 and 07 query the Enrichr web service and need internet access.
