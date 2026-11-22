# Revision analysis

Direct comparisons of females and males for the transcriptional response to vaccination and for
antibody titres. Run `THREADS=16 bash run_all.sh` from this folder after `original_analysis/`.
Summary tables are written to `tables/`; intermediate results, per-participant values and figures
to `results/` and `figures/`.

| Step | Script | Content |
|---|---|---|
| 01 | `01_prepare_data.py` | Sample sets, counts, cell-type fractions |
| 02 | `02_vst.R` | Variance-stabilising transformation |
| 03 | `03_population.py` | Participants and samples per arm and sex |
| 04 | `04_gene_level_deseq2.R` | Sex difference in the D1 vs D0 change per gene, DESeq2 with participant as blocking factor |
| 05 | `05_gene_level_mixed_model.R` | Same comparison with participant as random effect (dream) |
| 06 | `06_signature_scores.py` | Vaccine-response and ISG scores, sex difference in their change |
| 07 | `07_signature_mixed_model.R` | Scores, mixed model |
| 08 | `08_gene_level_summary.py` | Summary of the gene-level results |
| 09 | `09_variance_partitioning.py` | Partitioning of expression variance |
| 10 | `10_variance_summary_and_cells.py` | Variance summary, change in estimated cell proportions |
| 11 | `11_irf7_isg15.py` | IRF7 and ISG15 fold change by sex |
| 12 | `12_irf7_isg15_baseline.R` | IRF7 and ISG15 adjusted for D0, levels at each visit |
| 13 | `13_deg_counts_and_permutation.R` | Per-sex DEG counts at equal group size, sex-label permutation |
| 14 | `14_wgcna_robustness.R` | WGCNA module-trait FDR, networks rebuilt without each participant |
| 15 | `15_wgcna_foxp3_summary.py` | WGCNA summary, FOXP3 |
| 16 | `16_antibody.R` | Antibody titres by sex: regression at V3D28, longitudinal mixed model |
| 17 | `17_figures.py` | Figures |

Main model (step 04), per arm on all three vaccinations:
`counts ~ participant + vaccination + female:vaccination:D1 + male:vaccination:D1`; the sex
difference is the contrast `female:Tk:D1 - male:Tk:D1`.
