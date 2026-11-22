# Folder layout shared by the R scripts of the revision analysis. Source this file first.
script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
SCRIPTS <- script_dir()
ANALYSIS <- dirname(SCRIPTS)                                   # revision_analysis/
REPO <- dirname(ANALYSIS)
DATA <- file.path(REPO, "data")                                # input tables (see data/README.md)
RES <- file.path(ANALYSIS, "results")                          # intermediate results, one folder per step
PREP <- file.path(RES, "prepared")                             # sample sets, counts, VST, cell fractions
TAB <- file.path(ANALYSIS, "tables")                           # result tables (summary statistics)
PART <- file.path(RES, "participant_level")                    # per-participant values, used by the figures
FIG <- file.path(ANALYSIS, "figures")
WGCNA_ORIG <- file.path(REPO, "original_analysis", "results", "wgcna")   # networks of the paper
for (d in c(PREP, TAB, PART, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

res_dir <- function(name) {
  d <- file.path(RES, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
