#!/usr/bin/env python
"""Summary of the WGCNA robustness analysis, and FOXP3.

WGCNA, per arm (from 14_wgcna_robustness.R):
  module-trait cells of the matrix shown in the paper: number tested, p < 0.05, BH-FDR < 0.05
  cells with p < 0.05 that stay at p < 0.05 when any single participant is left out
  module eigengene mixed models: modules with a D1 vs D0 change and with a sex difference in the
    change at FDR < 0.05
  the vaccine-response modules of the paper (SI M19, LI M23): size and median Jaccard index with
    the best-matching module when the network is rebuilt without one participant
FOXP3: D1 vs D0 change per sex and the sex difference, DESeq2 basic model (04_gene_level_deseq2.R).

Inputs:  results/wgcna_robustness/, results/gene_level_deseq2/gene_level_results.csv.gz
Outputs: tables/wgcna_summary.csv, tables/foxp3_fold_change.csv
"""
import os

import pandas as pd

from paths import RES, TAB

W = os.path.join(RES, "wgcna_robustness")
mt = pd.read_csv(os.path.join(W, "module_trait_with_FDR.csv"))
mm = pd.read_csv(os.path.join(W, "module_mixed_model_sex_x_day.csv"))
lo = pd.read_csv(os.path.join(W, "loso_stability_nominal_cells.csv"))
jac = pd.read_csv(os.path.join(W, "loso_module_reconstruction_jaccard.csv"))
HIGHLIGHT = {"SI": "M19", "LI": "M23"}

rows = []
for arm in ["SI", "LI"]:
    t, m, l, j = (x[x.arm == arm] for x in (mt, mm, lo, jac))
    h = j[j.module == HIGHLIGHT[arm]].iloc[0]
    rows += [
        (arm, "Module-trait cells tested (matrix shown in the paper)", len(t)),
        (arm, "Cells with nominal p < 0.05", int((t.p < 0.05).sum())),
        (arm, "Cells with FDR < 0.05", int((t.p_BH_matrix < 0.05).sum())),
        (arm, "Nominal cells with p < 0.05 in all leave-one-participant-out networks",
         int((l.LOSO_frac_p05 == 1).sum())),
        (arm, "Modules (mixed model)", len(m)),
        (arm, "Modules with D1 vs D0 change, FDR < 0.05", int((m.padj_day < 0.05).sum())),
        (arm, "Modules with sex difference in D1 vs D0 change, FDR < 0.05", int((m.padj_sex_x_day < 0.05).sum())),
        (arm, f"{HIGHLIGHT[arm]}: genes", int(h.n_genes)),
        (arm, f"{HIGHLIGHT[arm]}: median Jaccard index when one participant is left out", round(h.jaccard_LOSO_median, 2)),
    ]
out = pd.DataFrame(rows, columns=["arm", "measure", "value"], dtype=object)
out.to_csv(os.path.join(TAB, "wgcna_summary.csv"), index=False)
print(out.to_string(index=False))

gl = pd.read_csv(os.path.join(RES, "gene_level_deseq2", "gene_level_results.csv.gz"))
f = gl[(gl.gene == "FOXP3") & (gl.model == "basic")]
f = f[["arm", "vaccination", "estimate", "log2FC", "CI_low", "CI_high", "pvalue", "padj"]]
f.to_csv(os.path.join(TAB, "foxp3_fold_change.csv"), index=False)
print(f.round(3).to_string(index=False))
