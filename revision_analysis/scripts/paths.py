"""Folder layout shared by the Python scripts of the revision analysis."""
import os

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(SCRIPTS)                       # revision_analysis/
REPO = os.path.dirname(ANALYSIS)
DATA = os.path.join(REPO, "data")                         # input tables (see data/README.md)
RES = os.path.join(ANALYSIS, "results")                   # intermediate results, one folder per step
PREP = os.path.join(RES, "prepared")                      # sample sets, counts, VST, cell fractions
TAB = os.path.join(ANALYSIS, "tables")                    # result tables (summary statistics)
PART = os.path.join(RES, "participant_level")             # per-participant values, used by the figures
FIG = os.path.join(ANALYSIS, "figures")                   # figures
WGCNA_ORIG = os.path.join(REPO, "original_analysis", "results", "wgcna")   # networks of the paper

for _d in (PREP, TAB, PART, FIG):
    os.makedirs(_d, exist_ok=True)


def res_dir(name):
    """Create and return results/<name>."""
    d = os.path.join(RES, name)
    os.makedirs(d, exist_ok=True)
    return d
