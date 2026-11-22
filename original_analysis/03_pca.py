"""03 Principal component analysis of VST expression.

Samples: the 174 samples without exclusion label (``pca`` set of 01_sample_set.py: SI, LI and
the placebo arm). Values: VST expression from 02_vst.R (computed on all 244 samples).

For each group (all 174 samples, placebo arm, SI, LI) separately:

1. the 2,000 most variable genes are selected with ``scanpy.pp.highly_variable_genes``
   (``n_top_genes=2000``, default ``flavor="seurat"``) on the VST values of that group;
2. ``scanpy.tl.pca`` with 30 components is computed on these genes (centred, not scaled);
3. PC1 and PC2 are plotted, coloured by sex, visit, day, arm (all samples only), age group
   and participant. In the participant panel, participants 109, 132, 154 and 156, which are
   set apart from the other samples in the PCA, are coloured individually and all others are
   grey ("Other").

Outputs:
  figures/pca/pca_<group>.pdf|png
  results/pca/pca_<group>_coordinates.csv   (30 PCs per sample)
  results/pca/pca_<group>_variance_ratio.csv
  results/pca/pca_<group>_hvg.txt           (the 2,000 genes)
"""

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

from common import FIGURES, RESULTS, read_sample_sets, read_vst, setup

setup()
figdir = FIGURES / "pca"
outdir = RESULTS / "pca"
figdir.mkdir(parents=True, exist_ok=True)
outdir.mkdir(parents=True, exist_ok=True)

N_TOP_GENES = 2000
N_COMPS = 30
HIGHLIGHT_SUBJECTS = [132, 156, 109, 154]

sets = read_sample_sets()
obs = sets[sets["pca"]].copy()
vst = read_vst().loc[obs.index]

adata = AnnData(vst.astype("float64"), obs=obs)
adata.obs["Subject_info"] = "Other"
for subject in HIGHLIGHT_SUBJECTS:
    adata.obs.loc[adata.obs.Subject == subject, "Subject_info"] = str(subject)
adata.obs["Day"] = [t[-1] for t in adata.obs["Timepoint"]]
for col in ["Sex", "Timepoint", "Day", "Interval", "Age_group_2", "Subject_info"]:
    adata.obs[col] = adata.obs[col].astype(str).astype("category")

groups = ["all", "Placebo", "SI", "LI"]
for group in groups:
    if group == "all":
        sub = adata.copy()
        cols = ["Sex", "Timepoint", "Day", "Interval", "Age_group_2", "Subject_info"]
        sc.set_figure_params(dpi=100, figsize=(13, 2.7))
    else:
        sub = adata[adata.obs.Interval == group].copy()
        cols = ["Sex", "Timepoint", "Day", "Age_group_2", "Subject_info"]
        sc.set_figure_params(dpi=100, figsize=(11, 2.7))
    sns.set(style="ticks", font_scale=0.9)
    for col in cols:
        sub.obs[col] = sub.obs[col].cat.remove_unused_categories()

    sc.pp.highly_variable_genes(sub, n_top_genes=N_TOP_GENES)
    sc.tl.pca(sub, n_comps=N_COMPS)

    pcs = pd.DataFrame(sub.obsm["X_pca"], index=sub.obs_names,
                       columns=[f"PC{i + 1}" for i in range(N_COMPS)])
    pd.concat([sub.obs[["Subject", "Timepoint", "Sex", "Interval"]], pcs], axis=1).to_csv(
        outdir / f"pca_{group}_coordinates.csv")
    pd.Series(sub.uns["pca"]["variance_ratio"], index=pcs.columns, name="variance_ratio").to_csv(
        outdir / f"pca_{group}_variance_ratio.csv")
    (outdir / f"pca_{group}_hvg.txt").write_text(
        "\n".join(sub.var_names[sub.var["highly_variable"]]) + "\n")

    plt.figure(figsize=(17, 3))
    for i, col in enumerate(cols, start=1):
        ax = plt.subplot(1, 6, i)
        sub.uns["Subject_info_colors"] = sns.color_palette(
            "colorblind", len(sub.obs["Subject_info"].cat.categories)).as_hex()
        sub.uns["Timepoint_colors"] = sns.color_palette(
            "colorblind", len(sub.obs["Timepoint"].cat.categories)).as_hex()
        sc.pl.pca(sub, color=[col], show=False, s=100, ax=ax, annotate_var_explained=True)
        ax.spines[["right", "top"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_edgecolor("lightgray")
        ax.set_xlabel(ax.get_xlabel(), color="gray")
        ax.set_ylabel(ax.get_ylabel(), color="gray")
        if col != cols[0]:
            ax.set_ylabel("")
        ax.set_title(col.capitalize())
    if group == "all":
        plt.suptitle(f"{N_TOP_GENES} most variable genes")
    else:
        plt.suptitle(f"{group}: {N_TOP_GENES} most variable genes")
    plt.tight_layout()
    plt.savefig(figdir / f"pca_{group}.pdf", bbox_inches="tight")
    plt.savefig(figdir / f"pca_{group}.png", bbox_inches="tight", dpi=500)
    plt.close()
    var = sub.uns["pca"]["variance_ratio"]
    print(f"{group}: {sub.n_obs} samples, PC1 {var[0]:.1%}, PC2 {var[1]:.1%}")
