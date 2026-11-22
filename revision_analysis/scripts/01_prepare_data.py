#!/usr/bin/env python
"""Sample sets, count matrix and cell-type fractions used by all later steps.

Cohort: high-dose (HD) participants of the short interval (SI) and long interval (LI) arms. Each
vaccinated participant received four injections, three vaccinations and one placebo: SI was
vaccinated at visits V1, V2 and V4 (placebo at V3), LI at V1, V3 and V4 (placebo at V2). Blood was
sampled before (D0) and one day after (D1) each injection.

Samples are excluded as in the analysis of the paper (original_analysis/):
  placebo visit of the arm         SI V3, LI V2 (Comment "Internal_placebo")
  incomplete visit series          participants 141, 158 and 165 (Comment "Incomplete")
  protocol deviation               placebo participant 133
  placebo arm                      not part of the vaccine-response analysis
The three vaccinations are labelled T1, T2 and T3 in both arms (SI: V1, V2, V4; LI: V1, V3, V4).

Analysis set: 23 participants (SI 8 female and 4 male, LI 3 female and 8 male), 6 samples each
(D0 and D1 of T1 to T3), 138 samples.
WGCNA set: the participants with all eight visits including the placebo visit (10 per arm,
160 samples), the samples of the networks shown in the paper.

Cell-type fractions: CIBERSORTx (LM22) estimates from the RNA-seq data, rescaled to sum to 1 and
collapsed to neutrophils, monocytes (monocytes, macrophages, dendritic cells), CD4 T cells,
CD8 and gamma-delta T cells, NK cells, B cells and other (eosinophils, mast cells).

Inputs (data/): sample_metadata.csv, counts_raw.csv.gz, gene_annotation.csv, cibersortx_fractions.csv
Outputs (results/prepared/):
  counts_raw_HD.csv.gz              genes x samples, integer counts, all 244 HD-cohort samples
  sample_metadata_HD.csv            one row per sample, with T1 to T3 labels and inclusion flags
  analysis_samples_vaccination.csv  the 138 samples of the analysis set
  wgcna_samples.csv                 the 160 samples of the WGCNA set
  cell_fractions_cibersortx.csv     lineage fractions of the 138 analysis samples
  gene_annotation.csv               Ensembl gene id, gene symbol and the unique symbol used as row name
"""
import os

import numpy as np
import pandas as pd
from anndata.utils import make_index_unique

from paths import DATA, PREP, WGCNA_ORIG

VACC_VISITS = {"SI": {"V1": "T1", "V2": "T2", "V4": "T3"},
               "LI": {"V1": "T1", "V3": "T2", "V4": "T3"}}

# ---------------------------------------------------------------- metadata
meta = pd.read_csv(os.path.join(DATA, "sample_metadata.csv"), index_col=0)
meta.index.name = "sample"
meta["Comment"] = meta["Comment"].astype(str)
assert meta.shape[0] == 244 and set(meta["Dose"]) == {"HD", "P"}

meta["visit"] = meta["Timepoint"].str[:2]
meta["Day"] = "D" + meta["Timepoint"].str[-1]
meta["Vaccination"] = [VACC_VISITS.get(i, {}).get(v, np.nan)
                       for i, v in zip(meta["Interval"], meta["visit"])]

placebo_visit = (((meta.Interval == "SI") & meta.Timepoint.str.contains("V3")) |
                 ((meta.Interval == "LI") & meta.Timepoint.str.contains("V2")))
# participants 141, 158 and 165 are labelled "Incomplete" at every visit, including the placebo visit
chk = meta.Comment != "Incomplete"
assert (placebo_visit[chk] == (meta.Comment[chk] == "Internal_placebo")).all()
assert set(meta.loc[meta.Comment == "Incomplete", "Subject"]) == {141, 158, 165}
assert set(meta.loc[meta.Comment.str.startswith("PD"), "Subject"]) == {133}

meta["exclusion_reason"] = ""
meta.loc[meta.Comment == "Internal_placebo", "exclusion_reason"] = "placebo visit of the arm"
meta.loc[meta.Comment == "Incomplete", "exclusion_reason"] = "incomplete visit series"
meta.loc[meta.Comment.str.startswith("PD"), "exclusion_reason"] = "protocol deviation"
meta.loc[meta.Interval == "Placebo", "exclusion_reason"] = (
    meta.loc[meta.Interval == "Placebo", "exclusion_reason"].replace("", "placebo arm"))

meta["in_vaccination_DE"] = (meta.Comment == "nan") & meta.Interval.isin(["SI", "LI"])
assert meta["Vaccination"][meta.in_vaccination_DE].notna().all()

w = meta[meta.Comment.isin(["Internal_placebo", "nan"]) & meta.Interval.isin(["SI", "LI"])]
n_tp = w.groupby("Subject").size()
meta["in_WGCNA"] = meta.index.isin(w.index) & meta.Subject.isin(n_tp[n_tp == 8].index)
wg_meta = [os.path.join(WGCNA_ORIG, f"metadata_{a}.csv") for a in ["SI", "LI"]]
if all(os.path.exists(f) for f in wg_meta):
    wg = pd.concat([pd.read_csv(f, index_col=0) for f in wg_meta])
    assert set(wg.index) == set(meta.index[meta.in_WGCNA]), "WGCNA sample set differs from the paper"

# the analysis set equals the HD samples of the CIBERSORTx table
cbx = pd.read_csv(os.path.join(DATA, "cibersortx_fractions.csv"), index_col=0)
de = set(meta.index[meta.in_vaccination_DE])
assert de == set(cbx.index[cbx.Dose == "HD"])

# ---------------------------------------------------------------- counts
counts = pd.read_csv(os.path.join(DATA, "counts_raw.csv.gz"), index_col=0).T   # genes x samples
assert set(counts.columns) == set(meta.index)
counts = counts[meta.index]

ann = pd.read_csv(os.path.join(DATA, "gene_annotation.csv"))[["gene_id", "gene_name"]]
uniq = list(make_index_unique(pd.Index(ann["gene_name"].astype(str))))
assert len(uniq) == counts.shape[0]
assert all(a == b for a, b in zip(uniq, counts.index)), "gene order differs from gene_annotation.csv"
ann["gene"] = uniq
ann.to_csv(os.path.join(PREP, "gene_annotation.csv"), index=False)

counts.index.name = "gene"
counts.astype(int).to_csv(os.path.join(PREP, "counts_raw_HD.csv.gz"))

# ---------------------------------------------------------------- tables
cols = ["Subject", "Sex", "Age", "Age_group_2", "Dose", "Interval", "Timepoint", "visit",
        "Day", "Vaccination", "Comment", "exclusion_reason", "in_vaccination_DE", "in_WGCNA"]
meta[cols].to_csv(os.path.join(PREP, "sample_metadata_HD.csv"))
meta.loc[meta.in_vaccination_DE, cols].to_csv(os.path.join(PREP, "analysis_samples_vaccination.csv"))
meta.loc[meta.in_WGCNA, cols].to_csv(os.path.join(PREP, "wgcna_samples.csv"))

ct_cols = ["NK", "CD4 T Naive", "CD8 T cells", "DC", "CD4 T memory (activated)", "Eosinophil",
           "TfH", "Mac", "Mono", "Neutrophil", "B cells", "Mast", "gdT", "Tregs",
           "CD4 T memory (resting)"]
fr = cbx.loc[sorted(de), ct_cols]
fr = fr.div(fr.sum(axis=1), axis=0)
lin = pd.DataFrame({
    "Neutrophil": fr["Neutrophil"],
    "Monocyte": fr["Mono"] + fr["Mac"] + fr["DC"],
    "T_CD4": fr[["CD4 T Naive", "CD4 T memory (activated)", "CD4 T memory (resting)", "TfH", "Tregs"]].sum(1),
    "T_CD8_gd": fr["CD8 T cells"] + fr["gdT"],
    "NK": fr["NK"],
    "B": fr["B cells"],
    "Other": fr["Eosinophil"] + fr["Mast"],
})
lin.index.name = "sample"
lin.to_csv(os.path.join(PREP, "cell_fractions_cibersortx.csv"))

# ---------------------------------------------------------------- summary
d = meta[meta.in_vaccination_DE]
print(d.groupby(["Interval", "Sex"]).Subject.nunique().unstack())
print(pd.crosstab([d.Interval, d.Vaccination, d.Day], d.Sex))
print(f"WGCNA set: {meta.in_WGCNA.sum()} samples, "
      f"{meta[meta.in_WGCNA].groupby('Interval').Subject.nunique().to_dict()} participants")
