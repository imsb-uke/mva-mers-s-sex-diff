"""01 Sample sets.

The sequenced high-dose (HD) cohort comprises 244 whole-blood samples from the short-interval
(SI) arm, the long-interval (LI) arm and the placebo arm. The ``Comment`` column of
``sample_metadata.csv`` marks the samples that are excluded:

* ``Internal_placebo``: D0 and D1 of the visit at which a vaccinated participant received the
  placebo injection (SI V3, LI V2).
* ``Incomplete``: participants 141, 158 and 165, whose sample series is incomplete.
* ``PD - No third vaccination``: placebo-arm participant 133 (protocol deviation).

The script derives the first two labels from visit and participant, checks them against the
metadata, and writes the membership of every sample in the sample sets used downstream:

* ``analysis`` (138 samples): SI and LI, D0 and D1 of the three vaccinations, participants
  without exclusion (SI 12 participants x 6 = 72, LI 11 x 6 = 66). Used for the per-sex
  differential expression, the DEG clustering and the Figure 6 tables.
* ``pca`` (174 samples): all samples without exclusion label, i.e. the analysis set plus the
  placebo arm (without participant 133).
* ``wgcna`` (160 samples): SI and LI participants with eight samples (D0 and D1 of all four
  injections, the placebo-injection visit included), 10 participants per arm. Participants
  106, 145 (SI) and 154 (LI) lack one sample and are not included. Participant 150 (SI) has two
  libraries labelled V3D0 and none labelled V3D1; the selection counts samples per participant,
  so participant 150 is included with both V3D0 libraries.

Outputs: results/sample_sets/sample_sets.csv (one row per sample, one column per set) and
results/sample_sets/sample_counts.csv (samples per set, arm, visit and sex).
"""

import pandas as pd

from common import RESULTS, PLACEBO_VISIT, read_metadata, setup

setup()
outdir = RESULTS / "sample_sets"
outdir.mkdir(parents=True, exist_ok=True)

meta = read_metadata()
assert meta.shape[0] == 244

# --- derive the exclusion labels from visit and participant ---------------------------
derived = pd.Series("nan", index=meta.index)
derived[meta.Subject == 133] = "PD - No third vaccination"
for arm, visit in PLACEBO_VISIT.items():
    derived[(meta.Interval == arm) & meta.Timepoint.str.contains(visit)] = "Internal_placebo"
derived[meta.Subject.isin([141, 158, 165])] = "Incomplete"
assert (derived == meta.Comment).all(), "exclusion labels differ from sample_metadata.csv"

# --- sample sets ---------------------------------------------------------------------
sets = meta[["Subject", "Timepoint", "Sex", "Age", "Age_group_2", "Interval", "Comment"]].copy()

# analysis: no exclusion label, vaccinated arms only
sets["analysis"] = (meta.Comment == "nan") & (meta.Interval != "Placebo")

# pca: no exclusion label (placebo arm kept)
sets["pca"] = meta.Comment == "nan"

# wgcna: participants with 8 samples among the
# vaccinated arms, internal-placebo visits included
pool = meta[meta.Comment.isin(["Internal_placebo", "nan"]) & meta.Interval.isin(["SI", "LI"])]
n_per_subject = pool.groupby("Subject").size()
complete = n_per_subject[n_per_subject == 8].index
sets["wgcna"] = meta.Subject.isin(complete)
print("Participants without 8 samples (not in WGCNA):",
      sorted(n_per_subject[n_per_subject != 8].index.tolist()))

assert sets["analysis"].sum() == 138
assert sets["pca"].sum() == 174
assert sets["wgcna"].sum() == 160

sets.to_csv(outdir / "sample_sets.csv")

# --- sample numbers -----------------------------------------------------------------
rows = []
for name in ["analysis", "pca", "wgcna"]:
    sub = sets[sets[name]]
    tab = sub.groupby(["Interval", "Timepoint", "Sex"]).size().rename("n_samples").reset_index()
    tab.insert(0, "sample_set", name)
    rows.append(tab)
    print(f"\n{name}: {sub.shape[0]} samples, "
          f"{sub.groupby('Interval').Subject.nunique().to_dict()} participants")
    print(pd.crosstab([sub.Interval, sub.Timepoint], sub.Sex))
pd.concat(rows).to_csv(outdir / "sample_counts.csv", index=False)
