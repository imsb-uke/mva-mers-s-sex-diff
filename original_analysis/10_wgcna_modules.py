"""10 WGCNA module codes, module gene tables and module-trait heatmaps.

Per arm, from the outputs of 09_wgcna.R:

1. ``<ARM>_module_genes.csv``: one column per module colour (columns in order of first
   appearance along the gene list), genes listed below, shorter columns padded with empty cells.
2. Module codes: the module colours are sorted alphabetically and numbered M1, M2, ...
   ("grey", the unassigned genes, receives a code like any other module: SI M7, LI M10). SI has
   20 modules (M1-M20) and LI 24 (M1-M24). The mapping, the genes per code and the module-trait
   tables with codes are written to ``codenames/``.
3. Module-trait heatmap with all 16 visit_sex traits (Pearson r, colour scale -0.6 to 0.6,
   star for p < 0.05, unadjusted).
4. Module-trait heatmap with 12 traits as shown in the paper: the arm's placebo-injection visit
   is dropped (SI V3, LI V2) and the remaining visits are renumbered so that V1, V2 and V3
   denote the first, second and third vaccination in both arms (SI: V4 -> V3; LI: V3 -> V2,
   V4 -> V3). Only the columns change; the correlations and p values are those of the network
   fitted on all 80 samples of the arm, placebo-visit samples included.

Outputs:
  results/wgcna/wgcna_results_<ARM>/<ARM>_module_genes.csv
  results/wgcna/wgcna_results_<ARM>/codenames/mapping_module_codenames.json
  results/wgcna/wgcna_results_<ARM>/codenames/<ARM>_module_genes_codenames.json|csv
  results/wgcna/wgcna_results_<ARM>/codenames/<ARM>_module_assignments_codenames.csv
  results/wgcna/wgcna_results_<ARM>/codenames/<ARM>_moduleTraitCor_codenames.csv
  results/wgcna/wgcna_results_<ARM>/codenames/<ARM>_moduleTraitPvalue_codenames.csv
  results/wgcna/wgcna_results_<ARM>/codenames/<ARM>_moduleTrait{Cor,Pvalue}_relabelled.csv
  figures/wgcna/<ARM>_module_trait_heatmap_codenames.pdf        (16 traits)
  figures/wgcna/<ARM>_module_trait_heatmap_relabelled.pdf|png   (12 traits)
"""

import json

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from common import FIGURES, RESULTS, setup

setup()
figdir = FIGURES / "wgcna"
figdir.mkdir(parents=True, exist_ok=True)

RELABEL = {
    "SI": {"drop": ["V3D0", "V3D1"], "rename": {"V4D0": "V3D0", "V4D1": "V3D1"}},
    "LI": {"drop": ["V2D0", "V2D1"],
           "rename": {"V3D0": "V2D0", "V3D1": "V2D1", "V4D0": "V3D0", "V4D1": "V3D1"}},
}
KEEP_ORDER = ["V1D0", "V1D1", "V2D0", "V2D1", "V3D0", "V3D1"]


def json_to_csv(module_data, path):
    """Modules as columns, genes as rows, padded with empty strings (json_to_csv.py)."""
    max_genes = max(len(g) for g in module_data.values())
    df = pd.DataFrame({m: g + [""] * (max_genes - len(g)) for m, g in module_data.items()})
    df.to_csv(path, index=False)


def convert_trait_table(path_in, path_out, mapping):
    """Module eigengene rows (ME<colour>) -> module code; drop 'timepoint' from columns."""
    df = pd.read_csv(path_in, index_col=0)
    df.index = [mapping[idx[2:]] for idx in df.index]
    df.columns = [c.replace("timepoint", "") for c in df.columns]
    df.to_csv(path_out)
    return df


def relabel_columns(df, arm):
    spec = RELABEL[arm]
    visit = df.columns.str.split("_").str[0]
    df = df.loc[:, ~visit.isin(spec["drop"])].copy()
    df.columns = [f'{spec["rename"].get(c.split("_")[0], c.split("_")[0])}_{c.split("_")[1]}'
                  for c in df.columns]
    order = [f"{v}_{s}" for v in KEEP_ORDER for s in ["Female", "Male"]]
    assert sorted(df.columns) == sorted(order)
    return df[order]


def heatmap(df_r, df_p, title, path, figsize, cbar_kws):
    short = [c.replace("Female", "F").replace("Male", "M") for c in df_r.columns]
    annots = df_p.map(lambda p: "*" if p < 0.05 else "")
    sns.set(style="ticks", font_scale=1.2)
    plt.figure(figsize=figsize)
    ax = sns.heatmap(df_r.set_axis(short, axis=1), cmap="RdBu_r", center=0, vmin=-0.6,
                     vmax=0.6, square=False, annot=annots.values, fmt="s",
                     annot_kws={"color": "white"}, cbar_kws=cbar_kws)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.title(title)
    for p in path:
        plt.savefig(p, bbox_inches="tight", dpi=300)
    plt.close()


for arm in ["SI", "LI"]:
    rdir = RESULTS / "wgcna" / f"wgcna_results_{arm}"
    cdir = rdir / "codenames"
    cdir.mkdir(parents=True, exist_ok=True)

    with open(rdir / f"{arm}_module_genes.json") as f:
        module_data = json.load(f)
    json_to_csv(module_data, rdir / f"{arm}_module_genes.csv")

    mapping = {name: f"M{i + 1}" for i, name in enumerate(sorted(module_data.keys()))}
    with open(cdir / "mapping_module_codenames.json", "w") as f:
        json.dump(mapping, f, indent=2)
    coded = {mapping[c]: g for c, g in module_data.items()}
    with open(cdir / f"{arm}_module_genes_codenames.json", "w") as f:
        json.dump(coded, f, indent=2)
    json_to_csv(coded, cdir / f"{arm}_module_genes_codenames.csv")

    assign = pd.read_csv(rdir / f"{arm}_module_assignments.csv", index_col=0)
    assign["module_code"] = assign["module"].map(mapping)
    assign.to_csv(cdir / f"{arm}_module_assignments_codenames.csv", index=False)

    df_r = convert_trait_table(rdir / f"{arm}_moduleTraitCor.csv",
                               cdir / f"{arm}_moduleTraitCor_codenames.csv", mapping)
    df_p = convert_trait_table(rdir / f"{arm}_moduleTraitPvalue.csv",
                               cdir / f"{arm}_moduleTraitPvalue_codenames.csv", mapping)
    print(f"{arm}: {len(mapping)} modules, sizes "
          f"{assign.module_code.value_counts().reindex(list(mapping.values())).to_dict()}")

    heatmap(df_r, df_p,
            f"{arm}\nCorrelation of WGCNA modules with traits\nStars indicate p-value < 0.05",
            [figdir / f"{arm}_module_trait_heatmap_codenames.pdf"], (7, 7),
            {"label": "Pearson's r"})

    order = sorted(df_r.index, key=lambda m: int(m[1:]))
    r12 = relabel_columns(df_r, arm).loc[order]
    p12 = relabel_columns(df_p, arm).loc[order]
    r12.to_csv(cdir / f"{arm}_moduleTraitCor_relabelled.csv")
    p12.to_csv(cdir / f"{arm}_moduleTraitPvalue_relabelled.csv")
    heatmap(r12, p12, arm,
            [figdir / f"{arm}_module_trait_heatmap_relabelled.pdf",
             figdir / f"{arm}_module_trait_heatmap_relabelled.png"], (6.4, 7.5),
            {"label": "Pearson's r", "shrink": 0.6})
