#!/usr/bin/env python
"""Participants and samples of the transcriptomic analysis.

Per arm and sex: high-dose participants of the mITT cohort, participants with RNA-seq data,
participants excluded for an incomplete visit series, analysed participants and samples, and age.
Age and sex of the analysed participants are checked against the mITT table.

Inputs:  data/participants_mITT.xlsx, results/prepared/sample_metadata_HD.csv,
         results/prepared/analysis_samples_vaccination.csv
Outputs: tables/population_transcriptome.csv, tables/samples_per_visit.csv
"""
import os

import pandas as pd

from paths import DATA, PREP, TAB

meta = pd.read_csv(os.path.join(PREP, "analysis_samples_vaccination.csv"), index_col=0)
allm = pd.read_csv(os.path.join(PREP, "sample_metadata_HD.csv"), index_col=0)
mitt = pd.read_excel(os.path.join(DATA, "participants_mITT.xlsx"))

subj = meta.drop_duplicates("Subject").set_index("Subject")[["Sex", "Age", "Interval"]]
chk = subj.join(mitt.set_index("Subject")[["Sex", "Age", "Dose", "Interval"]], rsuffix="_mITT")
assert (chk.Sex == chk.Sex_mITT).all() and (chk.Age == chk.Age_mITT).all()
assert (chk.Dose == "HD").all() and (chk.Interval == chk.Interval_mITT).all()
print("age and sex of the", len(subj), "analysed participants agree with the mITT table")

rows = []
for arm in ["SI", "LI"]:
    for sex in ["Female", "Male"]:
        m = mitt[(mitt.Dose == "HD") & (mitt.Interval == arm) & (mitt.Sex == sex)]
        s = allm[(allm.Interval == arm) & (allm.Sex == sex)]
        a = subj[(subj.Interval == arm) & (subj.Sex == sex)]
        rows.append({"arm": arm, "sex": sex, "mITT_high_dose": len(m),
                     "RNAseq_participants": s.Subject.nunique(),
                     "excluded_incomplete_series": s[s.Comment == "Incomplete"].Subject.nunique(),
                     "analysed_participants": len(a),
                     "analysed_samples": int(((meta.Interval == arm) & (meta.Sex == sex)).sum()),
                     "age_median": a.Age.median(), "age_range": f"{a.Age.min()} to {a.Age.max()}"})
pop = pd.DataFrame(rows)
pop.to_csv(os.path.join(TAB, "population_transcriptome.csv"), index=False)
print(pop.to_string(index=False))

vis = (meta.groupby(["Interval", "Vaccination", "Timepoint", "Sex"]).size().unstack(fill_value=0)
       .reset_index().rename(columns={"Interval": "arm"}))
vis.to_csv(os.path.join(TAB, "samples_per_visit.csv"), index=False)
