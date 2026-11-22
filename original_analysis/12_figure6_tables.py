"""12 Tables behind Figure 6 (TPM per participant and per-participant fold changes).

The Figure 6 heatmap of log2 fold changes is drawn in GraphPad Prism from the per-sex DESeq2
tables (04_deseq2_per_sex.R, exported by 05_deg_tables_and_counts.py); the gene panels are
drawn from the tables written here.

1. ``TPM_DEGs_<ARM>_with_metadata.csv``: TPM values of the DE genes of the arm (the gene set
   clustered in 06_de_gene_clustering.py: SI 2,426 genes, LI 1,857 genes) for the arm's samples
   of the 138-sample analysis set (SI 72, LI 66), followed by the sample metadata. Samples are
   ordered by visit, females before males.
2. ``<ARM>_log2fc_per_subject.xlsx``: for each participant and injection (sheets V1 to V4,
   including the arm's placebo-injection visit) the per-participant fold change
   log2(TPM_D1 + 1) - log2(TPM_D0 + 1) and the two TPM values, for the 10,000 genes of the
   WGCNA input (08_wgcna_input.py). Participants are all SI or LI participants without the
   "Incomplete" label (SI 12, LI 11). Where the D0 or D1 sample of a visit is missing, the row
   is NA: SI participants 106 and 145 at V3, LI participant 154 at V2, and SI participant 150 at
   V3 (two libraries labelled V3D0, none labelled V3D1).

TPM values are read from ``../data/tpm.csv.gz`` (quantification pipeline output, 6 significant
digits).

Outputs:
  results/figure6/TPM_DEGs_SI_with_metadata.csv, TPM_DEGs_LI_with_metadata.csv
  results/figure6/SI_log2fc_per_subject.xlsx, LI_log2fc_per_subject.xlsx
"""

import numpy as np
import pandas as pd

from common import RESULTS, read_metadata, read_sample_sets, read_tpm, setup

setup()
outdir = RESULTS / "figure6"
outdir.mkdir(parents=True, exist_ok=True)

meta = read_metadata()
sets = read_sample_sets()
tpm = read_tpm().loc[meta.index]
obs_cols = ["Subject", "Timepoint", "Project", "Age", "Sex", "Dose", "Interval", "Comment",
            "Age_group_4", "Age_group_2"]

# --- 1. TPM of the clustered DE genes -------------------------------------------------
for arm in ["SI", "LI"]:
    genes = pd.read_csv(RESULTS / "de_gene_clustering" / f"{arm}_gene_clusters.csv")["gene"].tolist()
    order = pd.read_csv(RESULTS / "de_gene_clustering" / f"{arm}_scaled_expression_DEGs.csv",
                        index_col=0, usecols=[0]).index
    assert sets.loc[order, "analysis"].all() and (sets.loc[order, "Interval"] == arm).all()
    df = pd.concat([tpm.loc[order, genes], meta.loc[order, obs_cols]], axis=1)
    df.to_csv(outdir / f"TPM_DEGs_{arm}_with_metadata.csv")
    print(f"TPM_DEGs_{arm}: {df.shape[0]} samples, {len(genes)} genes")

# --- 2. per-participant log2 fold changes -----------------------------------------------
wgcna_genes = pd.read_csv(RESULTS / "wgcna" / "counts_SI.csv", index_col=0, nrows=1).columns
assert (pd.read_csv(RESULTS / "wgcna" / "counts_LI.csv", index_col=0, nrows=1).columns == wgcna_genes).all()

pool = meta[meta.Comment.isin(["Internal_placebo", "nan"]) & meta.Interval.isin(["SI", "LI"])].copy()
pool["Vaccine_No"] = pool.Timepoint.str.extract(r"V(\d+)")[0]
pool["Day"] = pool.Timepoint.str.extract(r"D(\d+)")[0]
X = tpm.loc[pool.index, wgcna_genes]

for arm in ["SI", "LI"]:
    sub = pool[pool.Interval == arm]
    donors = sub.Subject.unique()
    sheets, missing = {}, {}
    for vaccine in sorted(sub.Vaccine_No.unique()):
        parts = []
        for donor in donors:
            d = sub[sub.Subject == donor]
            d0 = d[(d.Vaccine_No == vaccine) & (d.Day == "0")].index
            d1 = d[(d.Vaccine_No == vaccine) & (d.Day == "1")].index
            part = pd.DataFrame({"Subject": donor, "Interval": d.Interval.iloc[0],
                                 "Sex": d.Sex.iloc[0], "Gene": wgcna_genes})
            if len(d0) > 0 and len(d1) > 0:
                day0 = X.loc[d0[0]].values
                day1 = X.loc[d1[0]].values
                part["log2FC"] = np.log2(day1 + 1) - np.log2(day0 + 1)
                part["TPM_D0"] = day0
                part["TPM_D1"] = day1
            else:
                missing.setdefault(donor, []).append(vaccine)
                part["log2FC"] = np.nan
                part["TPM_D0"] = np.nan
                part["TPM_D1"] = np.nan
            parts.append(part)
        sheets[f"V{vaccine}"] = pd.concat(parts, ignore_index=True)
    with pd.ExcelWriter(outdir / f"{arm}_log2fc_per_subject.xlsx", engine="openpyxl") as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)
    print(f"{arm}_log2fc_per_subject: {len(donors)} participants; missing D0/D1 pairs: {missing}")
