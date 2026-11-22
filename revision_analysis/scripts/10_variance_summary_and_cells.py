#!/usr/bin/env python
"""Variance partitioning in eight categories, and change in estimated cell proportions.

a  Share of expression variance per factor, for all expressed genes and for the first principal
   component of the 2,000 most variable genes, from 09_variance_partitioning.py, grouped into:
   D1 vs D0 (vaccination), individual differences, sex, age, arm, vaccination number and its
   interactions with arm and day, sex x response (sex x D1 and sex x vaccination x D1), unexplained.
b  Change from D0 to D1 in estimated neutrophil, monocyte and lymphocyte proportions (CIBERSORTx,
   percentage points) per participant and vaccination.

Inputs:  results/variance_partitioning/, results/prepared/analysis_samples_vaccination.csv,
         cell_fractions_cibersortx.csv
Outputs: tables/variance_partitioning.csv, results/participant_level/cell_proportion_change.csv
"""
import os

import pandas as pd

from paths import PART, PREP, RES, TAB

VD = os.path.join(RES, "variance_partitioning")
meta = pd.read_csv(os.path.join(PREP, "analysis_samples_vaccination.csv"), index_col=0)
summ = pd.read_csv(os.path.join(VD, "variance_decomposition_summary.csv"))
pcs = pd.read_csv(os.path.join(VD, "pca_per_pc_decomposition.csv"), index_col=0)
cells = pd.read_csv(os.path.join(PREP, "cell_fractions_cibersortx.csv"), index_col=0).loc[meta.index]

GROUPS = {  # plain category -> components of the decomposition
    "D1 vs D0 (vaccination)": ["within: Day 1 vs day 0 (vaccination)"],
    "Individual": ["between: individual (residual)"],
    "Sex": ["between: Sex"],
    "Age": ["between: Age"],
    "Arm (SI, LI)": ["between: Arm"],
    "Vaccination number\nand arm interactions": [
        "within: Vaccination number (visit baseline)", "within: Arm x vaccination number",
        "within: Vaccination number x day", "within: Arm x day", "within: Arm x vaccination number x day"],
    "Sex x response": ["within: Sex x day", "within: Sex x vaccination number x day"],
    "Unexplained": ["within: residual"],
}
allg = summ[summ.gene_set == "all expressed genes"].set_index("component").aggregate_pct
rows = []
for g, comps in GROUPS.items():
    rows.append({"factor": g.replace("\n", " "), "all_genes_pct": allg[comps].sum(),
                 "PC1_pct": 100 * pcs.loc["PC1", comps].sum()})
var = pd.DataFrame(rows)
var.to_csv(os.path.join(TAB, "variance_partitioning.csv"), index=False)
print(var.round(1).to_string(index=False), "\nPC1 explains", round(pcs.loc["PC1", "pct_variance_of_PC"], 1), "% of HVG variance")

# cell proportions
c = cells.copy()
c["Lymphocytes"] = c[["T_CD4", "T_CD8_gd", "NK", "B"]].sum(1)
c = c.rename(columns={"Neutrophil": "Neutrophils", "Monocyte": "Monocytes"})
CT = ["Neutrophils", "Monocytes", "Lymphocytes"]
d = c[CT].join(meta[["Subject", "Sex", "Interval", "Vaccination", "Day"]])
w = d.pivot_table(index=["Subject", "Sex", "Interval", "Vaccination"], columns="Day", values=CT)
ch = pd.DataFrame({ct: 100 * (w[(ct, "D1")] - w[(ct, "D0")]) for ct in CT}).reset_index()
ch.to_csv(os.path.join(PART, "cell_proportion_change.csv"), index=False)
