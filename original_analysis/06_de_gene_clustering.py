"""06 Clustering of DE genes and GO enrichment per cluster (Enrichr), SI and LI.

Per arm:

1. Gene set: union of the genes with padj <= 0.05 and |log2FC| >= log2(1.5) in the per-sex
   D1 vs D0 tables of the arm (04_deseq2_per_sex.R), both sexes combined, without genes that
   have no symbol (Ensembl IDs) and without the 5.8S rRNA gene.
   SI: union of the DEGs of the V1 and V4 comparisons (2,426 genes). ``SI_INCLUDE_V2 = True``
   adds the DEGs of the V2 comparisons (5 further genes, 2,431).
   LI: union of the DEGs of the V1, V3 and V4 comparisons (1,857 genes).
2. Expression: VST values (02_vst.R) of the arm's samples of the 138-sample analysis set
   (SI 72, LI 66), scaled per gene to zero mean and unit variance across these samples
   (``scanpy.pp.scale``, max_value=10; all genes are scaled, then the DE genes are selected).
3. Clustering: hierarchical clustering of the genes (rows; samples are not clustered) with
   ``scipy.cluster.hierarchy.linkage(method="ward", metric="euclidean")`` and
   ``fcluster(..., 3, criterion="maxclust")``, i.e. three gene clusters. Genes are sorted
   alphabetically before clustering, which fixes the cluster numbering.
4. Heatmap: samples ordered by visit, females before males within each visit; colour scale
   -4 to 4; participant IDs as column labels.
5. Enrichment: each cluster's gene list is submitted to Enrichr through GSEApy
   (``gseapy.enrichr``, libraries GO_Biological_Process_2023, GO_Molecular_Function_2023 and
   GO_Cellular_Component_2023, organism human, Enrichr default background). Terms with
   adjusted p <= 0.05 are kept; the first ten rows of this table (GSEApy returns the libraries
   in the order above, each sorted by p value) are printed next to the cluster in the heatmap.
   This needs GSEApy and internet access; with ``RUN_ENRICHR = False``, or if Enrichr cannot be
   reached, the clustering and the heatmap are produced without terms.

Outputs (per arm):
  results/de_gene_clustering/<ARM>_gene_clusters.csv           gene, cluster
  results/de_gene_clustering/<ARM>_scaled_expression_DEGs.csv  matrix used for clustering
  results/de_gene_clustering/<ARM>_enrichr_cluster<k>.csv      Enrichr results (all terms)
  figures/de_gene_clustering/<ARM>_heatmap.svg|pdf|png
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage

from common import (FIGURES, RESULTS, VACCINATION_VISITS, deg_table_names, read_deg_tables,
                    read_sample_sets, read_vst, setup, significant)

SI_INCLUDE_V2 = False
RUN_ENRICHR = True
N_CLUSTERS = 3
ENRICHR_LIBRARIES = ["GO_Biological_Process_2023",
                     "GO_Molecular_Function_2023",
                     "GO_Cellular_Component_2023"]

setup()
outdir = RESULTS / "de_gene_clustering"
figdir = FIGURES / "de_gene_clustering"
outdir.mkdir(parents=True, exist_ok=True)
figdir.mkdir(parents=True, exist_ok=True)

sets = read_sample_sets()
vst = read_vst()

T_COLORS = ["#00CED1", "#40E0D0", "#FF1493", "#FF69B4", "#8A2BE2", "#9370DB"]
SEX_COLORS = {"Male": "#000080", "Female": "#FFBF00"}
CLUSTER_COLORS = ["#32CD32", "#8B4513", "#800080"]


def de_gene_union(arm):
    names = deg_table_names(arm)
    if arm == "SI" and not SI_INCLUDE_V2:
        names = [n for n in names if "V2" not in n]
    degs = read_deg_tables(names)
    genes = set()
    for name in names:
        genes |= set(significant(degs[name]).index)
    genes = [g for g in genes if "ENSG" not in g and g != "X5_8S_rRNA"]
    # R column names back to gene symbols ("-" had become ".")
    return sorted(g.replace(".", "-") for g in genes), names


def run_enrichr(gene_list):
    import gseapy as gp
    gs = gp.enrichr(gene_list=gene_list, gene_sets=ENRICHR_LIBRARIES,
                    organism="human", outdir=None)
    return gs.results.copy()


for arm in ["SI", "LI"]:
    genes, used_tables = de_gene_union(arm)
    print(f"{arm}: {len(genes)} DE genes from {len(used_tables)} tables")

    obs = sets[sets["analysis"] & (sets.Interval == arm)]
    adata = AnnData(vst.loc[obs.index].astype("float64"), obs=obs.copy())
    sc.pp.scale(adata, max_value=10)
    adata = adata[:, genes]

    # sample order: visit, then females before males (metadata order within groups)
    idxs = []
    for tp in sorted(adata.obs.Timepoint.unique()):
        sub = adata.obs[adata.obs.Timepoint == tp]
        idxs += sub[sub.Sex == "Female"].index.tolist() + sub[sub.Sex == "Male"].index.tolist()
    adata = adata[idxs]
    df_to_plot = adata.to_df().T  # genes x samples

    # --- hierarchical clustering of genes ---------------------------------------
    linkage_matrix = linkage(df_to_plot, method="ward", metric="euclidean")
    clusters = fcluster(linkage_matrix, N_CLUSTERS, criterion="maxclust")
    gene_clusters = pd.DataFrame({"gene": df_to_plot.index, "cluster": clusters})
    gene_clusters.to_csv(outdir / f"{arm}_gene_clusters.csv", index=False)
    pd.concat([df_to_plot.T, adata.obs[["Subject", "Timepoint", "Sex"]]], axis=1).to_csv(
        outdir / f"{arm}_scaled_expression_DEGs.csv")
    print(gene_clusters.cluster.value_counts().sort_index().to_dict())

    # --- enrichment per cluster -------------------------------------------------
    pathways = {}
    for k in sorted(gene_clusters.cluster.unique()):
        pathways[k] = []
        if not RUN_ENRICHR:
            continue
        try:
            res = run_enrichr(gene_clusters.loc[gene_clusters.cluster == k, "gene"].tolist())
        except Exception as e:  # no internet access or service unavailable
            print(f"  Enrichr not available ({type(e).__name__}); heatmap without terms")
            RUN_ENRICHR = False
            continue
        res.to_csv(outdir / f"{arm}_enrichr_cluster{k}.csv", index=False)
        sig = res[res["Adjusted P-value"] <= 0.05]
        pathways[k] = [t.split(" (")[0] for t in sig["Term"].tolist()[:10]]

    # --- heatmap -------------------------------------------------------------------
    timepoints = sorted(adata.obs.Timepoint.unique())
    t_colordict = dict(zip(timepoints, T_COLORS))
    colors_df = pd.DataFrame({"Sex": [SEX_COLORS[s] for s in adata.obs.Sex],
                              "Timepoint": [t_colordict[t] for t in adata.obs.Timepoint]},
                             index=df_to_plot.columns)
    order = leaves_list(linkage_matrix)
    # clusters from top to bottom of the heatmap, coloured in that order
    cl_top_to_bottom = list(dict.fromkeys(clusters[order]))
    cluster_colordict = dict(zip(cl_top_to_bottom, CLUSTER_COLORS))
    cluster_colors = [cluster_colordict[c] for c in clusters]
    xlabels = adata.obs.loc[df_to_plot.columns, "Subject"].astype(str).tolist()

    sns.set(style="white", font_scale=1.2)
    kws = dict(cbar_kws=dict(ticks=[-4, 0, 4], orientation="horizontal"))
    g = sns.clustermap(df_to_plot, row_linkage=linkage_matrix, center=0, vmax=4, vmin=-4,
                       cmap="RdBu_r", figsize=(15, 15), col_colors=colors_df, col_cluster=False,
                       row_colors=cluster_colors, rasterized=True, xticklabels=True, **kws)
    g.ax_heatmap.set_xticklabels(xlabels, fontsize=10, ha="center", rotation=90)
    g.ax_heatmap.set_yticks([], [])
    g.ax_heatmap.set_ylabel(f"{df_to_plot.shape[0]} genes")
    g.ax_heatmap.set_xlabel("")

    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_position([x0, 0.85, g.ax_row_dendrogram.get_position().width, 0.01])
    g.ax_cbar.set_xlim(-4, 4)
    g.ax_cbar.set_title("Scaled gene expression")
    g.ax_cbar.tick_params(axis="x", length=10)
    for spine in g.ax_cbar.spines:
        g.ax_cbar.spines[spine].set_color("crimson")
        g.ax_cbar.spines[spine].set_linewidth(2)

    fig = plt.gcf()
    leg1 = fig.legend([Patch(facecolor=c) for c in SEX_COLORS.values()], SEX_COLORS,
                      title="Sex", bbox_to_anchor=(1, 1), loc="upper right", frameon=False)
    leg2 = fig.legend([Patch(facecolor=c) for c in t_colordict.values()], t_colordict,
                      title="Timepoint", bbox_to_anchor=(0.8, 1), loc="upper right",
                      frameon=False)
    leg3 = fig.legend([Patch(facecolor=cluster_colordict[k]) for k in sorted(cluster_colordict)],
                      [str(k) for k in sorted(cluster_colordict)], title="Gene cluster",
                      bbox_to_anchor=(-0.1, 0.5), loc="upper right", frameon=False)
    fig.suptitle(arm)

    # GO terms next to each cluster block (heatmap rows are in dendrogram order)
    ordered_clusters = clusters[order]
    ax = g.ax_heatmap
    x_text = df_to_plot.shape[1] + 12
    ax.text(x_text, 0, "GO Biological Process", fontweight="bold",
            va="bottom", clip_on=False)
    for k in cl_top_to_bottom:
        rows = np.where(ordered_clusters == k)[0]
        ax.text(x_text, rows.mean(), "\n".join(pathways[k]), fontsize=13, va="center",
                clip_on=False)

    for ext in ["svg", "pdf"]:
        plt.savefig(figdir / f"{arm}_heatmap.{ext}", bbox_inches="tight", dpi=1000)
    plt.savefig(figdir / f"{arm}_heatmap.png", bbox_inches="tight", dpi=150)
    plt.close("all")
