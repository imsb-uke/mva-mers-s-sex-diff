#!/usr/bin/env python
"""Summary of the gene-level sex-difference tests per arm and vaccination.

From the DESeq2 models (04_gene_level_deseq2.R) and the mixed model (05_gene_level_mixed_model.R):
  responsive genes: D1 vs D0 change at FDR < 0.05 in females or in males (basic DESeq2 model)
  genes with a sex difference in the response at FDR < 0.05, per model version
  share of the responsive genes with a sex difference at p < 0.05, DESeq2 and mixed model
  Storey estimate of the share of responsive genes with a non-zero sex difference
    (1 - share of p > 0.5 / 0.5)
  smallest female vs male difference in fold change detectable with 80% power, from the median
    standard error of the sex difference over the responsive genes:
      fold difference = 2 ^ ((z[1 - alpha/2] + z[0.80]) x SE)
    for alpha = 0.05 (one gene selected in advance) and alpha = 0.001 (approximately the per-gene
    p value that passes BH-FDR 0.05 over about 15,500 genes when a few hundred genes pass)
  genes with a female vs male difference at D0 in the mixed model

Inputs:  results/gene_level_deseq2/, results/gene_level_mixed/
Output:  tables/gene_level_summary.csv
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

from paths import RES, TAB

OD = os.path.join(RES, "gene_level_deseq2")
ARMS, VACS = ["SI", "LI"], ["T1", "T2", "T3"]
MODELS = {"basic": "Unadjusted", "age": "Adjusted for age",
          "age_cells": "Adjusted for age and cell proportions"}

gl = pd.read_csv(os.path.join(OD, "gene_level_results.csv.gz"))
summ = pd.read_csv(os.path.join(OD, "gene_level_summary.csv"))

rows = []
for arm in ARMS:
    for vac in VACS:
        b = gl[(gl.arm == arm) & (gl.vaccination == vac) & (gl.model == "basic")]
        w = b.pivot(index="gene", columns="estimate", values="padj")
        resp = w.index[(w.female_D1_vs_D0 < 0.05) | (w.male_D1_vs_D0 < 0.05)]
        sd = b[(b.estimate == "sex_difference_F_minus_M") & b.gene.isin(resp)]
        se = sd.SE.median()
        row = {"arm": arm, "vaccination": vac, "n_responsive_genes": len(resp),
               "detectable_fold_difference_alpha0.05": 2 ** ((stats.norm.ppf(0.975) + stats.norm.ppf(0.8)) * se),
               "detectable_fold_difference_alpha0.001": 2 ** ((stats.norm.ppf(0.9995) + stats.norm.ppf(0.8)) * se),
               "pct_responsive_genes_nominal_p05": 100 * (sd.pvalue < 0.05).mean(),
               "storey_share_with_sex_difference": 1 - min(1, (sd.pvalue > 0.5).mean() / 0.5)}
        for mk in MODELS:
            s = summ[(summ.arm == arm) & (summ.vaccination == vac) & (summ.model == mk) &
                     (summ.estimate == "sex_difference_F_minus_M")].iloc[0]
            row[f"n_sex_difference_FDR05_{mk}"] = s.n_FDR05
            for e in ["female_D1_vs_D0", "male_D1_vs_D0"]:
                row[f"n_{e}_FDR05_{mk}"] = summ[(summ.arm == arm) & (summ.vaccination == vac) &
                                                (summ.model == mk) & (summ.estimate == e)].n_FDR05.iloc[0]
        rows.append(row)
tab = pd.DataFrame(rows)
# mixed model with the participant as random effect (05_gene_level_mixed_model.R): genes with a
# female vs male difference at D0 and in the response, FDR < 0.05
mx = pd.read_csv(os.path.join(RES, "gene_level_mixed", "mixed_model_summary.csv"))
mres = pd.read_csv(os.path.join(RES, "gene_level_mixed", "mixed_model_results.csv.gz"))
pct = []
for arm, vac in zip(tab.arm, tab.vaccination):
    b = gl[(gl.arm == arm) & (gl.vaccination == vac) & (gl.model == "basic")]
    w = b.pivot(index="gene", columns="estimate", values="padj")
    resp = w.index[(w.female_D1_vs_D0 < 0.05) | (w.male_D1_vs_D0 < 0.05)]
    x = mres[(mres.arm == arm) & (mres.vaccination == vac) & (mres.contrast == "response") & mres.gene.isin(resp)]
    pct.append(100 * (x.pvalue < 0.05).mean())
tab["pct_responsive_genes_p05_mixed_model"] = pct      # responsive genes (as for DESeq2)
for c in ["D0", "response"]:
    x = mx[mx.contrast == c].set_index(["arm", "vaccination"])
    tab[f"n_{c}_sex_difference_FDR05_mixed_model"] = [x.loc[(a, v), "n_FDR05"] for a, v in zip(tab.arm, tab.vaccination)]
    tab[f"pct_{c}_sex_difference_p05_mixed_model"] = [x.loc[(a, v), "pct_p05"] for a, v in zip(tab.arm, tab.vaccination)]
tab.to_csv(os.path.join(TAB, "gene_level_summary.csv"), index=False)
print(tab.round(2).T.to_string())
