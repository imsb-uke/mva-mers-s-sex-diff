#!/usr/bin/env python
"""Partitioning of whole-blood expression variance into between- and within-participant sources.

Every participant contributes the same set of six samples (D0 and D1 of three vaccinations), so the
total variance of each gene splits exactly into
  between-participant variance (participant means): arm, sex, age (sequential sums of squares),
    remaining individual differences
  within-participant variance (participant-centred data): vaccination number (D0 level of each
    vaccination), arm x vaccination number, D1 vs D0, vaccination number x D1, arm x D1,
    arm x vaccination number x D1, sex x D1, sex x vaccination number x D1, residual
Within-participant terms are estimated from participant-centred data, so stable differences
between individuals cannot absorb the vaccination effect. Sex terms are entered last; the unique
(drop-one) sum of squares of every term is reported as well.

Summaries per gene set (all expressed genes, and the 2,000 highly variable genes used for the PCA):
  aggregate share  = sum over genes of SS(term) / sum over genes of SS(total), the share of the total
                     variance of the expression matrix
  median per-gene share
The same partitioning is applied to the first 10 principal components of the 2,000 genes
(principal variance component analysis, PCs weighted by their variance). Sensitivity: the
within-participant partitioning with participant-centred CIBERSORTx fractions entered before or
after the vaccination terms.

Inputs:  results/prepared/analysis_samples_vaccination.csv, vst_vaccination_samples.csv.gz,
         cell_fractions_cibersortx.csv
Outputs: results/variance_partitioning/ (variance_decomposition_summary.csv,
         pca_per_pc_decomposition.csv, hvg2000.csv, per-gene shares)
"""
import os

import numpy as np
import pandas as pd

from paths import PREP, res_dir

OD = res_dir("variance_partitioning")

meta = pd.read_csv(os.path.join(PREP, "analysis_samples_vaccination.csv"), index_col=0)
vst = pd.read_csv(os.path.join(PREP, "vst_vaccination_samples.csv.gz"), index_col=0)[meta.index]
cells = pd.read_csv(os.path.join(PREP, "cell_fractions_cibersortx.csv"), index_col=0).loc[meta.index]
assert meta.groupby("Subject").size().eq(6).all(), "the decomposition requires 6 samples per participant"

donor = meta.Subject.values
D1 = (meta.Day == "D1").astype(float).values
T2 = (meta.Vaccination == "T2").astype(float).values
T3 = (meta.Vaccination == "T3").astype(float).values
FEM = (meta.Sex == "Female").astype(float).values
SI = (meta.Interval == "SI").astype(float).values
AGE = (meta.Age - meta.Age.mean()).values.astype(float)


def donor_center(x):
    x = pd.DataFrame(x)
    return (x - x.groupby(donor).transform("mean")).values


WITHIN_TERMS = [  # (label, columns) - entered in this order
    ("Vaccination number (visit baseline)", np.c_[T2, T3]),
    ("Arm x vaccination number", np.c_[SI * T2, SI * T3]),
    ("Day 1 vs day 0 (vaccination)", np.c_[D1]),
    ("Vaccination number x day", np.c_[T2 * D1, T3 * D1]),
    ("Arm x day", np.c_[SI * D1]),
    ("Arm x vaccination number x day", np.c_[SI * T2 * D1, SI * T3 * D1]),
    ("Sex x day", np.c_[FEM * D1]),
    ("Sex x vaccination number x day", np.c_[FEM * T2 * D1, FEM * T3 * D1]),
]
BETWEEN_TERMS = [("Arm", np.c_[SI]), ("Sex", np.c_[FEM]), ("Age", np.c_[AGE])]


def proj_ss(Y, X):
    """Sum of squares of the projection of every column of Y onto span(X). Y: samples x genes."""
    if X.shape[1] == 0:
        return np.zeros(Y.shape[1])
    Q, R = np.linalg.qr(X)
    Q = Q[:, np.abs(np.diag(R)) > 1e-10]
    return ((Q.T @ Y) ** 2).sum(axis=0)


def sequential(Y, terms, centre):
    """Sequential and unique SS for a list of (label, X); X is centred with `centre`."""
    Xs = [(lab, centre(X)) for lab, X in terms]
    seq, prev, cols = {}, np.zeros(Y.shape[1]), np.zeros((Y.shape[0], 0))
    for lab, X in Xs:
        cols = np.c_[cols, X]
        cur = proj_ss(Y, cols)
        seq[lab] = cur - prev
        prev = cur
    full = prev
    uniq = {lab: full - proj_ss(Y, np.column_stack([X for l2, X in Xs if l2 != lab]))
            for lab, _ in Xs}
    return seq, uniq, full


def decompose(Y, extra_within=None, extra_first=True):
    """Y: samples x features (columns centred inside). Returns DataFrame features x components
    (SS), sequential, plus unique SS for within/between terms."""
    Y = Y - Y.mean(axis=0)
    total = (Y ** 2).sum(axis=0)
    Yw = donor_center(Y)
    within_tot = (Yw ** 2).sum(axis=0)
    Yb = Y - Yw                                   # participant means, repeated per sample
    terms = list(WITHIN_TERMS)
    if extra_within is not None:
        terms = ([("Cell-type fractions", extra_within)] + terms) if extra_first else \
                (terms + [("Cell-type fractions", extra_within)])
    wseq, wuniq, wfull = sequential(Yw, terms, donor_center)
    bseq, buniq, bfull = sequential(Yb, BETWEEN_TERMS, lambda X: X - X.mean(axis=0))
    out = {f"between: {k}": v for k, v in bseq.items()}
    out["between: individual (residual)"] = (Yb ** 2).sum(axis=0) - bfull
    out.update({f"within: {k}": v for k, v in wseq.items()})
    out["within: residual"] = within_tot - wfull
    seq = pd.DataFrame(out)
    seq["total"] = total
    uniq = pd.DataFrame({**{f"between: {k}": v for k, v in buniq.items()},
                         **{f"within: {k}": v for k, v in wuniq.items()}})
    assert np.allclose(seq.drop(columns="total").sum(axis=1), total)
    return seq, uniq


def summarise(seq, uniq, label):
    comps = [c for c in seq.columns if c != "total"]
    frac = seq[comps].div(seq["total"], axis=0)
    rows = []
    for c in comps:
        rows.append({"gene_set": label, "component": c,
                     "aggregate_pct": 100 * seq[c].sum() / seq["total"].sum(),
                     "median_gene_pct": 100 * frac[c].median(),
                     "q75_gene_pct": 100 * frac[c].quantile(0.75),
                     "unique_aggregate_pct": (100 * uniq[c].sum() / seq["total"].sum()
                                              if c in uniq else np.nan),
                     "n_genes_gt10pct": int((frac[c] > 0.10).sum())})
    within = [c for c in comps if c.startswith("within")]
    wtot = seq[within].sum(axis=1)
    for c in within:
        rows.append({"gene_set": label, "component": c.replace("within:", "share of within-participant:"),
                     "aggregate_pct": 100 * seq[c].sum() / wtot.sum(),
                     "median_gene_pct": 100 * (seq[c] / wtot).median()})
    return pd.DataFrame(rows), frac


# ------------------------------------------------------------------ gene sets
import scanpy as sc
ad = sc.AnnData(vst.T.values.copy())
ad.var_names = vst.index
sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat")
hvg = vst.index[ad.var.highly_variable.values]
pd.Series(hvg).to_csv(os.path.join(OD, "hvg2000.csv"), index=False, header=["gene"])

summ, per_gene = [], {}
for label, genes in [("all expressed genes", vst.index), ("2,000 HVGs (PCA genes)", hvg)]:
    Y = vst.loc[genes].T.values
    seq, uniq = decompose(Y)
    s, frac = summarise(seq, uniq, label)
    summ.append(s)
    per_gene[label] = frac.set_index(pd.Index(genes))

# PCA / PVCA on HVGs
Y = vst.loc[hvg].T.values
Yc = Y - Y.mean(axis=0)
U, S, Vt = np.linalg.svd(Yc, full_matrices=False)
pcs = U * S
var_exp = S ** 2 / (S ** 2).sum()
k = 10
seq_pc, uniq_pc = decompose(pcs[:, :k])
frac_pc = seq_pc.drop(columns="total").div(seq_pc["total"], axis=0)
frac_pc.index = [f"PC{i+1}" for i in range(k)]
frac_pc.insert(0, "pct_variance_of_PC", 100 * var_exp[:k])
frac_pc.to_csv(os.path.join(OD, "pca_per_pc_decomposition.csv"))
w = var_exp[:k] / var_exp[:k].sum()
pv = (frac_pc.drop(columns="pct_variance_of_PC").mul(w, axis=0)).sum()
summ.append(pd.DataFrame({"gene_set": f"PVCA (top {k} PCs of HVGs, {100*var_exp[:k].sum():.0f}% var)",
                          "component": pv.index, "aggregate_pct": 100 * pv.values}))

# cell-type sensitivity (within-participant)
cf = cells[["Neutrophil", "Monocyte", "T_CD4", "T_CD8_gd", "NK", "B"]].values
for label, genes in [("all expressed genes", vst.index), ("2,000 HVGs (PCA genes)", hvg)]:
    Y = vst.loc[genes].T.values
    for first in (True, False):
        seq, uniq = decompose(Y, cf, extra_first=first)
        s, _ = summarise(seq, uniq, f"{label} + cell fractions ({'first' if first else 'last'})")
        summ.append(s)

S_ = pd.concat(summ, ignore_index=True)
S_.to_csv(os.path.join(OD, "variance_decomposition_summary.csv"), index=False)
for label, f in per_gene.items():
    f.to_csv(os.path.join(OD, f"per_gene_fractions_{'hvg' if 'HVG' in label else 'all'}.csv.gz"))
pd.set_option("display.width", 200)
print(S_[S_.gene_set.isin(["all expressed genes", "2,000 HVGs (PCA genes)"]) |
         S_.gene_set.str.startswith("PVCA")].round(2).to_string(index=False))
print(frac_pc.round(3).to_string())
print("done")
