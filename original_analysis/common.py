"""Shared paths, constants and helpers for the Python steps.

All scripts are run from the ``original_analysis/`` folder, for example
``python 03_pca.py``. Inputs are read from ``../data``; outputs go to
``results/`` and ``figures/``.
"""

import os
import random
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"

# Thresholds for calling a gene differentially expressed (D1 vs D0)
PADJ_CUTOFF = 0.05
FC_CUTOFF = 1.5
LOG2FC_CUTOFF = np.log2(FC_CUTOFF)

# Vaccination visits per arm. Each participant received four injections; the
# visit not listed here was the arm's placebo injection (SI V3, LI V2).
VACCINATION_VISITS = {"SI": ["V1", "V2", "V4"], "LI": ["V1", "V3", "V4"]}
PLACEBO_VISIT = {"SI": "V3", "LI": "V2"}
ARMS = ["SI", "LI"]
SEXES = ["Male", "Female"]

SEED = 42


def setup(seed=SEED):
    """Fix random seeds and set the plotting style used for the figures."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set(style="ticks", font_scale=1)
    # The figures were drawn in Arial; fall back to a metric-compatible or
    # generic sans-serif font when Arial is not installed.
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    for d in (RESULTS, FIGURES):
        d.mkdir(parents=True, exist_ok=True)


def read_metadata():
    """Sample metadata (244 high-dose cohort samples incl. the placebo arm)."""
    meta = pd.read_csv(DATA / "sample_metadata.csv", index_col=0)
    # Samples without exclusion label are compared against the string "nan"
    # (as written in the metadata file).
    meta["Comment"] = meta["Comment"].astype(str)
    return meta


def read_sample_sets():
    """Sample-set membership written by 01_sample_set.py."""
    return pd.read_csv(RESULTS / "sample_sets" / "sample_sets.csv", index_col=0)


def read_vst():
    """VST expression of all 244 samples (samples x genes), from 02_vst.R."""
    return pd.read_csv(RESULTS / "vst" / "vst_counts.csv.gz", index_col=0)


def read_tpm():
    """TPM matrix of all 244 samples (samples x genes)."""
    return pd.read_csv(DATA / "tpm.csv.gz", index_col=0)


def deg_dir():
    return RESULTS / "deseq2_per_sex"


def deg_table_names(arm=None):
    """File names of the 12 per-sex D1 vs D0 tables used in the paper."""
    names = []
    for a in ARMS if arm is None else [arm]:
        for sex in ["Female", "Male"]:
            for v in VACCINATION_VISITS[a]:
                names.append(f"{a}_{sex}_{v}D1_vs_{v}D0.csv")
    return names


def read_deg_tables(names):
    """Read DESeq2 result tables, indexed by gene, with |log2FC| added.

    Gene identifiers in these tables are R column names of the count matrix
    (``read.csv`` default ``check.names = TRUE``): "-" became "." and names
    starting with a digit got an "X" prefix. Downstream steps convert "." back
    to "-" where they look genes up in the expression matrices.
    """
    degs = {}
    for name in names:
        df = pd.read_csv(deg_dir() / name, index_col=0).set_index("gene_id")
        assert df.index.is_unique, f"duplicated genes in {name}"
        df["abs_log2FoldChange"] = df["log2FoldChange"].abs()
        degs[name] = df
    return degs


def significant(df):
    """Rows with padj <= 0.05 and |log2FC| >= log2(1.5)."""
    return df[(df.padj <= PADJ_CUTOFF) & (df.abs_log2FoldChange >= LOG2FC_CUTOFF)]
