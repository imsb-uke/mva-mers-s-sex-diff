"""07 GO enrichment of the up-regulated DEGs of each SI comparison, and overlap between sexes.

For each SI per-sex table (V1, V2, V4; 04_deseq2_per_sex.R), the up-regulated DEGs
(padj <= 0.05, log2FC >= log2(1.5)) are submitted to Enrichr through GSEApy
(GO_Biological_Process_2023, GO_Molecular_Function_2023, GO_Cellular_Component_2023;
organism human; Enrichr default background). Terms with adjusted p <= 0.05 are kept. Per
vaccination, the significant terms of females and males are compared: number of terms shared
by both sexes and number found only in females or only in males. The analysis covers the SI
arm.

Gene identifiers are submitted as they appear in the DESeq2 tables (R column names).

Requires GSEApy and internet access (Enrichr web service).

Outputs:
  results/enrichr_per_comparison/<table>_up_enrichr.csv      all terms
  results/enrichr_per_comparison/SI_term_overlap_by_sex.csv  shared / female-only / male-only
"""

import pandas as pd

from common import RESULTS, VACCINATION_VISITS, deg_table_names, read_deg_tables, setup, significant

LIBRARIES = ["GO_Biological_Process_2023",
             "GO_Molecular_Function_2023",
             "GO_Cellular_Component_2023"]

setup()
import gseapy as gp  # noqa: E402

outdir = RESULTS / "enrichr_per_comparison"
outdir.mkdir(parents=True, exist_ok=True)

names = deg_table_names("SI")
degs = read_deg_tables(names)

enr = {}
for name in names:
    sub = significant(degs[name])
    geneset_up = sub[sub.log2FoldChange > 0].index.tolist()
    gs = gp.enrichr(gene_list=geneset_up, gene_sets=LIBRARIES, organism="human", outdir=None)
    res = gs.results.copy()
    res.to_csv(outdir / name.replace(".csv", "_up_enrichr.csv"), index=False)
    enr[name] = res[res["Adjusted P-value"] <= 0.05]
    print(f"{name}: {len(geneset_up)} up-regulated DEGs, {enr[name].shape[0]} terms")

rows = []
for v in VACCINATION_VISITS["SI"]:
    terms_f = enr[f"SI_Female_{v}D1_vs_{v}D0.csv"].Term.tolist()
    terms_m = enr[f"SI_Male_{v}D1_vs_{v}D0.csv"].Term.tolist()
    rows.append({"visit": v,
                 "shared": len(set(terms_f) & set(terms_m)),
                 "female_only": len([t for t in terms_f if t not in terms_m]),
                 "male_only": len([t for t in terms_m if t not in terms_f])})
overlap = pd.DataFrame(rows)
overlap.to_csv(outdir / "SI_term_overlap_by_sex.csv", index=False)
print(overlap)
