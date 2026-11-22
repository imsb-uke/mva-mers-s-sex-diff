# Input data

Both analysis folders start from the tables in this folder. The files contain participant-level
data and are not part of the repository.

| File | Content |
|---|---|
| `counts_raw.csv.gz` | Gene-level read counts of the 244 whole-blood RNA-seq samples of the high-dose (HD) cohort and the placebo arm. Rows are samples, columns are gene symbols (62,754 genes, made unique by appending `-1`, `-2`, ... to repeated symbols). Quantification: nf-core/rnaseq 3.12.0, STAR 2.7.9a and Salmon 1.10.1 against GRCh38 (Ensembl 110, GENCODE annotation). |
| `tpm.csv.gz` | Transcripts per million of the same samples and genes, same layout as the counts. |
| `sample_metadata.csv` | One row per sample: participant (`Subject`), study visit and day (`Timepoint`, for example `V2D1` = one day after the injection at visit 2), age, sex, dose (`HD`, `P` for placebo), interval (`SI` short, `LI` long, `Placebo`), age group, and `Comment` marking the placebo visit of each arm (`Internal_placebo`), participants with an incomplete visit series (`Incomplete`) and a protocol deviation. |
| `gene_annotation.csv` | Ensembl gene id and gene symbol in the row order of the count matrix. |
| `cibersortx_fractions.csv` | CIBERSORTx (LM22 signature) estimates of blood-cell fractions for 338 RNA-seq samples of the high-dose, low-dose and placebo participants, with the sample metadata. The revision analysis uses the 138 high-dose analysis samples. |
| `antibody_titres.xlsx` | Sheet `Antibody data`: MERS-CoV S1 binding antibodies (`ELISA_S1`) and neutralising antibodies (`VSV_PRNT50`) per participant (`Subject ID`) and study day (`T`, days 0 to 252), in IU/mL (`Unit` = `IU/mL`; for PRNT50 the raw titre is also given, `Unit` = `titer`). Values at the lower limit of the assay are 4.6 IU/mL (ELISA) and 16 IU/mL (PRNT50). |
| `wgcna_genes_10000.txt` | The 10,000 genes of the WGCNA networks shown in the paper (the most variable genes of the VST values of all 244 samples), one symbol per line. |
| `participants_mITT.xlsx` | The 125 participants of the modified intention-to-treat cohort: age, sex, dose (`HD`, `LD`, `P`) and interval (`SI`, `LI`). |

Samples are named by their sequencing identifier (for example `221000000120.A04`); participants by
their study number.
