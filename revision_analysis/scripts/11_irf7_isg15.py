#!/usr/bin/env python
"""IRF7 and ISG15: expression per participant at D0 and D1, and the sex difference in the change.

Expression: log2(TPM + 1) at D0 and D1 of each vaccination for every participant (the values shown
as expression levels in the paper). Sex difference in the D1 vs D0 log2 fold change: from the
DESeq2 model with the participant as blocking factor (04_gene_level_deseq2.R, basic model), with
95% CI, Wald p value and Holm correction over the 12 comparisons (2 genes x 2 arms x 3 vaccinations).

Inputs:  data/tpm.csv.gz, results/prepared/analysis_samples_vaccination.csv,
         results/gene_level_deseq2/gene_level_results.csv.gz
Outputs: tables/irf7_isg15_fold_change.csv (estimates per sex and sex difference),
         results/participant_level/irf7_isg15_per_participant.csv (log2(TPM + 1) at D0 and D1)
"""
import os

import numpy as np
import pandas as pd

from paths import DATA, PART, PREP, RES, TAB

GENES, ARMS, VACS = ["IRF7", "ISG15"], ["SI", "LI"], ["T1", "T2", "T3"]

meta = pd.read_csv(os.path.join(PREP, "analysis_samples_vaccination.csv"), index_col=0)
tpm = pd.read_csv(os.path.join(DATA, "tpm.csv.gz"), index_col=0,
                  usecols=lambda c: c.startswith("Unnamed") or c in GENES).T[meta.index]
gl = pd.read_csv(os.path.join(RES, "gene_level_deseq2", "gene_level_results.csv.gz"))
gl = gl[gl.gene.isin(GENES) & (gl.model == "basic")]

rows = []
for g in GENES:
    for arm in ARMS:
        for vac in VACS:
            m = meta[(meta.Interval == arm) & (meta.Vaccination == vac)]
            for s, d in m.groupby("Subject"):
                d0 = np.log2(tpm.loc[g, d[d.Day == "D0"].index[0]] + 1)
                d1 = np.log2(tpm.loc[g, d[d.Day == "D1"].index[0]] + 1)
                rows.append({"gene": g, "arm": arm, "vaccination": vac, "Subject": s,
                             "sex": d.Sex.iloc[0], "D0": d0, "D1": d1, "log2FC": d1 - d0})
pp = pd.DataFrame(rows)
st = gl.pivot_table(index=["gene", "arm", "vaccination"], columns="estimate",
                    values=["log2FC", "CI_low", "CI_high", "pvalue"])
st.columns = [f"{v}_{e}" for v, e in st.columns]
st = st.reset_index()
# Holm over the 12 sex-difference tests (2 genes x 2 arms x 3 vaccinations)
pv = st["pvalue_sex_difference_F_minus_M"].values
o = np.argsort(pv)
st.loc[st.index[o], "p_holm_12_sex_difference_F_minus_M"] = np.minimum(1, np.maximum.accumulate((len(pv) - np.arange(len(pv))) * pv[o]))
st.to_csv(os.path.join(TAB, "irf7_isg15_fold_change.csv"), index=False)
pp.to_csv(os.path.join(PART, "irf7_isg15_per_participant.csv"), index=False)
