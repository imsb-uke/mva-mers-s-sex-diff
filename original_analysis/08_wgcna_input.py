"""08 WGCNA input matrices.

1. Genes: the 10,000 most variable genes, selected with ``scanpy.pp.highly_variable_genes``
   (``n_top_genes=10000``, default ``flavor="seurat"``) on the VST values (02_vst.R) of all
   244 samples of the high-dose cohort. The same 10,000 genes are used for both arms.
   The gene list of the networks in the paper is provided as ``../data/wgcna_genes_10000.txt``
   and used when present; recomputing the selection with another DESeq2 version can exchange
   a few genes with tied dispersions at the cut-off. The script reports how many recomputed
   genes differ from the list.
2. Samples: the ``wgcna`` set of 01_sample_set.py, i.e. SI and LI participants with eight
   samples (D0 and D1 of all four injections, including the arm's placebo-injection visit);
   80 samples (10 participants) per arm.
3. Trait label: ``Timepoint`` is replaced by "<visit>_<sex>" (e.g. ``V2D1_Female``); 09_wgcna.R
   correlates module eigengenes with the indicator of each of these 16 groups.

Outputs:
  results/wgcna/counts_SI.csv, counts_LI.csv      VST values, samples x 10,000 genes
  results/wgcna/metadata_SI.csv, metadata_LI.csv  sample metadata with visit_sex label
"""

import pandas as pd
import scanpy as sc
from anndata import AnnData

from common import DATA, RESULTS, read_metadata, read_sample_sets, read_vst, setup

N_TOP_GENES = 10000

setup()
outdir = RESULTS / "wgcna"
outdir.mkdir(parents=True, exist_ok=True)

meta = read_metadata()
sets = read_sample_sets()
vst = read_vst()
assert (vst.index == meta.index).all()

adata = AnnData(vst.astype("float64"), obs=meta.copy())
sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)
hvg = adata.var_names[adata.var["highly_variable"]]

gene_file = DATA / "wgcna_genes_10000.txt"
if gene_file.exists():
    paper_genes = set(gene_file.read_text().split())
    print(f"Recomputed vs provided gene set: {len(set(hvg) - paper_genes)} genes differ")
    genes = adata.var_names[adata.var_names.isin(paper_genes)]  # keeps the matrix gene order
    assert len(genes) == len(paper_genes)
else:
    genes = hvg
adata = adata[:, genes]

to_save = adata[sets.loc[adata.obs_names, "wgcna"].values].copy()
to_save.obs["Timepoint"] = to_save.obs["Timepoint"].astype(str) + "_" + to_save.obs["Sex"].astype(str)

for arm in to_save.obs.Interval.unique():
    tmp = to_save[to_save.obs.Interval == arm]
    tmp.to_df().to_csv(outdir / f"counts_{arm}.csv")
    tmp.obs.to_csv(outdir / f"metadata_{arm}.csv")
    print(f"{arm}: {tmp.n_obs} samples, {tmp.obs.Subject.nunique()} participants, "
          f"{tmp.n_vars} genes")
