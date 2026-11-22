#!/usr/bin/env python
"""Sex difference in the response to vaccination at the participant level (signature scores).

One response value per participant and vaccination: the change from D0 to D1 (VST, log2 scale)
averaged over a gene set.
  Vaccine-response signature: mean of the 200 most induced genes minus mean of the 200 most
    repressed genes of that arm and vaccination. Genes are ranked by the paired D1 vs D0 t
    statistic over all participants of the arm and vaccination, without using sex.
  ISG signature: mean of 30 interferon-stimulated genes (fixed list below, includes IRF7 and ISG15).
Female minus male difference in the change per arm and vaccination, linear regression:
  Unadjusted                              change ~ sex
  Adjusted for age                        change ~ sex + age
  Adjusted for D0 level                   change ~ sex + D0 score
  Adjusted for age and cell proportions   change ~ sex + age + D1 - D0 change in the neutrophil
                                          and monocyte proportions
The unadjusted comparison is also tested by exact permutation of the sex labels (all assignments of
the female labels within the arm: 495 in SI, 165 in LI).
The score of every sample (D0 and D1) is written for the mixed model in 07_signature_mixed_model.R.

Inputs:  results/prepared/analysis_samples_vaccination.csv, vst_vaccination_samples.csv.gz,
         cell_fractions_cibersortx.csv
Outputs: results/participant_level/participant_signature_scores.csv (one row per participant and
         vaccination),
         tables/signature_sex_difference.csv, results/signature_scores/sample_scores.csv
"""
import itertools
import os

import numpy as np
import pandas as pd
from scipy import stats

from paths import PART, PREP, TAB, res_dir

ISGS = ["IFI27", "IFI44L", "IFI44", "RSAD2", "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "MX1",
        "MX2", "OAS1", "OAS2", "OAS3", "OASL", "HERC5", "SIGLEC1", "IRF7", "STAT1", "EPSTI1",
        "LY6E", "XAF1", "CMPK2", "USP18", "SPATS2L", "IFITM3", "BST2", "LAMP3", "SERPING1", "GBP1"]
N_SIG = 200
ARMS, VACS = ["SI", "LI"], ["T1", "T2", "T3"]
MODELS = {"Unadjusted": ["female"], "Adjusted for age": ["female", "age"],
          "Adjusted for D0 level": ["female", "D0"],
          "Adjusted for age and cell proportions": ["female", "age", "d_neut", "d_mono"]}
OD = res_dir("signature_scores")


def ols_female(d, covars):
    """Least squares fit of score on an intercept + covars; estimate, 95% CI and p for 'female'."""
    X = np.column_stack([np.ones(len(d))] + [d[c].values.astype(float) for c in covars])
    y = d.score.values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    df = len(y) - X.shape[1]
    s2 = ((y - X @ beta) ** 2).sum() / df
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
    q = stats.t.ppf(0.975, df)
    return beta[1], beta[1] - q * se, beta[1] + q * se, 2 * stats.t.sf(abs(beta[1] / se), df)

meta = pd.read_csv(os.path.join(PREP, "analysis_samples_vaccination.csv"), index_col=0)
vst = pd.read_csv(os.path.join(PREP, "vst_vaccination_samples.csv.gz"), index_col=0)[meta.index]
cells = pd.read_csv(os.path.join(PREP, "cell_fractions_cibersortx.csv"), index_col=0).loc[meta.index]

# ---------------------------------------------------------------- scores
rows, long = [], []
for arm in ARMS:
    for vac in VACS:
        m = meta[(meta.Interval == arm) & (meta.Vaccination == vac)]
        d0 = m[m.Day == "D0"].reset_index().set_index("Subject")["sample"]
        d1 = m[m.Day == "D1"].reset_index().set_index("Subject")["sample"]
        subj = d0.index.intersection(d1.index)
        r = pd.DataFrame(vst[d1[subj]].values - vst[d0[subj]].values, index=vst.index, columns=subj)
        t = r.mean(axis=1) / (r.std(axis=1, ddof=1) / np.sqrt(r.shape[1]))
        up, down = t.nlargest(N_SIG).index, t.nsmallest(N_SIG).index
        isg = [g for g in ISGS if g in r.index]
        dc = cells.loc[d1[subj]].values - cells.loc[d0[subj]].values
        dc = pd.DataFrame(dc, index=subj, columns=cells.columns)
        lev = {"vaccine_response_signature": vst.loc[up].mean() - vst.loc[down].mean(),
               "ISG_signature": vst.loc[isg].mean()}          # score of every sample
        for s in subj:
            info = m[m.Subject == s].iloc[0]
            row = {"Subject": s, "arm": arm, "vaccination": vac, "sex": info.Sex, "age": info.Age,
                   "d_neut": dc.loc[s, "Neutrophil"], "d_mono": dc.loc[s, "Monocyte"]}
            for k, v in lev.items():
                row[k] = v[d1[s]] - v[d0[s]]
                row[f"{k}_D0"] = v[d0[s]]
                for day, smp in (("D0", d0[s]), ("D1", d1[s])):
                    long.append({"Subject": s, "arm": arm, "vaccination": vac, "sex": info.Sex,
                                 "age": info.Age, "day": day, "score": k, "value": v[smp]})
            rows.append(row)
pd.DataFrame(long).to_csv(os.path.join(OD, "sample_scores.csv"), index=False)
sc = pd.DataFrame(rows)
sc["female"] = (sc.sex == "Female").astype(int)
sc.to_csv(os.path.join(PART, "participant_signature_scores.csv"), index=False)


def perm_p(y, f):
    """Exact two-sided permutation p for the female minus male difference in means."""
    n, k = len(y), int(f.sum())
    obs = y[f == 1].mean() - y[f == 0].mean()
    diffs = []
    for idx in itertools.combinations(range(n), k):
        g = np.zeros(n, bool); g[list(idx)] = True
        diffs.append(y[g].mean() - y[~g].mean())
    return np.mean(np.abs(diffs) >= abs(obs) - 1e-12)


res = []
for score in ["vaccine_response_signature", "ISG_signature"]:
    for arm in ARMS:
        for vac in VACS:
            d = sc[(sc.arm == arm) & (sc.vaccination == vac)].rename(columns={score: "score", f"{score}_D0": "D0"})
            for lab, f in MODELS.items():
                est, lo, hi, p = ols_female(d, f)
                res.append({"score": score, "arm": arm, "vaccination": vac, "model": lab,
                            "n_female": int(d.female.sum()), "n_male": int((1 - d.female).sum()),
                            "mean_female": d.score[d.female == 1].mean(),
                            "mean_male": d.score[d.female == 0].mean(),
                            "difference_F_minus_M": est, "CI_low": lo, "CI_high": hi,
                            "p_regression": p,
                            "p_exact_permutation": perm_p(d.score.values, d.female.values)
                            if lab == "Unadjusted" else np.nan})
res = pd.DataFrame(res)
res.to_csv(os.path.join(TAB, "signature_sex_difference.csv"), index=False)
print(res.round(3).to_string(index=False))

