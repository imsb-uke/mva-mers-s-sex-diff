# Sex differences in immune responses to the MVA-MERS-S vaccine

Code for the whole-blood transcriptome and antibody analyses of the paper.

| Folder | Content |
|---|---|
| `original_analysis/` | Transcriptome analyses as reported in the paper (PCA, per-sex differential expression, clustering, WGCNA) |
| `revision_analysis/` | Direct comparisons of females and males (transcriptome and antibody titres) |
| `environment/` | Conda environments (R and Python) |
| `data/` | Input tables, not included (see `data/README.md`) |

## Run

```
mamba env create -f environment/r_environment.yml
mamba env create -f environment/python_environment.yml
export PYTHON="conda run --no-capture-output -n mva_mers_sex_py python"
export RSCRIPT="conda run --no-capture-output -n mva_mers_sex_r Rscript"
cd original_analysis && bash run_all.sh && cd ..
cd revision_analysis && THREADS=16 bash run_all.sh
```

`revision_analysis/` uses the WGCNA networks of `original_analysis/`, so run that first. Each
script starts with a description of its inputs, methods and outputs. Outputs are written to
`results/`, `figures/` and `logs/` of each analysis folder. License: MIT.
