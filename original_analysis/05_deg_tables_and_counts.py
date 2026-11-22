"""05 DEG tables and DEG counts per sex.

Input: the 12 per-sex DESeq2 tables from 04_deseq2_per_sex.R (design ~ Timepoint).

A gene is called differentially expressed (DEG) if padj <= 0.05 and |log2FC| >= log2(1.5).
Genes whose padj is NA (removed by DESeq2 independent filtering or Cook's outlier detection)
are not DEGs.

Outputs:
  results/deg_tables/complete_tables/<ARM>/<ARM>_<Sex>_<Vx>D1_vs_<Vx>D0.xlsx  full tables
  results/deg_tables/filtered_tables/<ARM>/<ARM>_<Sex>_<Vx>D1_vs_<Vx>D0.xlsx  DEGs only
  results/deg_tables/DE_counts.csv                                          DEGs per table
  figures/deg_counts/deg_counts_<ARM>.pdf|png|svg                           bar plots
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from common import (ARMS, FIGURES, RESULTS, deg_table_names, read_deg_tables, setup,
                    significant)

setup()
outdir = RESULTS / "deg_tables"
figdir = FIGURES / "deg_counts"
figdir.mkdir(parents=True, exist_ok=True)

names = deg_table_names()
degs = read_deg_tables(names)

rows = []
for name in names:
    arm, sex, visit = name.split("_")[:3]
    full = degs[name].drop(columns="abs_log2FoldChange")
    filt = significant(degs[name]).drop(columns="abs_log2FoldChange")
    for kind, df in [("complete_tables", full), ("filtered_tables", filt)]:
        d = outdir / kind / arm
        d.mkdir(parents=True, exist_ok=True)
        df.to_excel(d / name.replace(".csv", ".xlsx"), index=True, engine="openpyxl")
    rows.append({"table": name, "n_degs": filt.shape[0],
                 "n_up": int((filt.log2FoldChange > 0).sum()),
                 "n_down": int((filt.log2FoldChange < 0).sum()),
                 "treatment": arm, "timepoint": visit, "sex": sex})

df_counts = pd.DataFrame(rows).set_index("table").sort_values(by=["treatment", "timepoint", "sex"])
df_counts.to_csv(outdir / "DE_counts.csv")
print(df_counts)

sns.set(style="ticks", font_scale=1)
for arm in ARMS:
    plt.figure(figsize=(3, 3))
    sns.barplot(data=df_counts[df_counts.treatment == arm], x="timepoint", y="n_degs",
                hue="sex", hue_order=["Male", "Female"])
    plt.title(arm)
    plt.legend(frameon=False)
    plt.ylabel("Number of degs\n(Adj. P<=0.05, FC>=1.5)", fontsize=10)
    plt.xlabel("")
    for ext, kw in [("pdf", {}), ("png", {"dpi": 500}), ("svg", {})]:
        plt.savefig(figdir / f"deg_counts_{arm}.{ext}", bbox_inches="tight", **kw)
    plt.close()
