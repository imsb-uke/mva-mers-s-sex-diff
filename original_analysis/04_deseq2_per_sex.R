# 04 Per-sex differential expression, D1 vs D0 of each vaccination
#
# For every arm (SI, LI), vaccination visit (SI: V1, V2, V4; LI: V1, V3, V4) and sex, DESeq2 is
# fitted on the raw counts of the D0 and D1 samples of that visit (138-sample analysis set of
# 01_sample_set.py), and the contrast Timepoint VxD1 vs VxD0 is extracted: Wald test,
# Benjamini-Hochberg adjusted p values, DESeq2 defaults for independent filtering and Cook's
# cutoff, no log2 fold change shrinkage, no pre-filtering of genes. Each table comes from its
# own fit with its own size factors and dispersions.
#
# Design ~ Timepoint (no covariates; D0 and D1 samples are not paired by participant).
#
# Gene identifiers: the count matrix is read with read.csv defaults (check.names = TRUE), so the
# gene_id column holds R syntactic names of the gene symbols ("-" written as ".", an "X" prefix
# for symbols starting with a digit). The Python steps map "." back to "-".
#
# Input : ../data/counts_raw.csv.gz, ../data/sample_metadata.csv
# Output: results/deseq2_per_sex/<ARM>_<Sex>_<Vx>D1_vs_<Vx>D0.csv  (12 tables)
#
# Run from original_analysis/:  Rscript 04_deseq2_per_sex.R

suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(parallel)
})

N_CORES <- 6

out_dir <- "results/deseq2_per_sex"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df_counts <- read.csv(gzfile("../data/counts_raw.csv.gz"), row.names = 1)
df_meta <- read.csv("../data/sample_metadata.csv", row.names = 1, header = TRUE)
stopifnot(identical(rownames(df_counts), rownames(df_meta)))

# Samples without exclusion label have Comment "nan" in the metadata file. Keep these
# samples from the two vaccinated arms.
no_comment <- is.na(df_meta$Comment) | df_meta$Comment %in% c("", "nan")
df_meta <- df_meta[no_comment & df_meta$Interval %in% c("SI", "LI"), ]
df_counts <- df_counts[rownames(df_meta), ]
stopifnot(nrow(df_meta) == 138)

visits <- list(SI = c("1", "2", "4"), LI = c("1", "3", "4"))
sexes <- c("Male", "Female")

jobs <- list()
for (condition in names(visits)) {
  for (Timepoint in visits[[condition]]) {
    for (sex in sexes) {
      jobs[[length(jobs) + 1]] <- list(condition = condition, Timepoint = Timepoint, sex = sex)
    }
  }
}

fit_and_save <- function(inmeta, indata, design, pair_1, pair_0, file) {
  dds <- DESeqDataSetFromMatrix(countData = indata, colData = inmeta, design = design)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("Timepoint", pair_1, pair_0))
  res <- as.data.frame(res)
  res <- tibble::rownames_to_column(res, "gene_id")
  res <- dplyr::arrange(res, padj)
  write.csv(res, file)
}

run_job <- function(job) {
  pair_0 <- paste0("V", job$Timepoint, "D0")
  pair_1 <- paste0("V", job$Timepoint, "D1")
  inmeta <- df_meta[df_meta$Interval == job$condition & df_meta$Sex == job$sex &
                      df_meta$Timepoint %in% c(pair_0, pair_1), ]
  indata <- t(df_counts[rownames(inmeta), ])
  name <- paste0(job$condition, "_", job$sex, "_", pair_1, "_vs_", pair_0, ".csv")

  fit_and_save(inmeta, indata, ~ Timepoint, pair_1, pair_0, file.path(out_dir, name))
  sprintf("%s: %d samples (%s)", name, nrow(inmeta),
          paste(names(table(inmeta$Timepoint)), table(inmeta$Timepoint), collapse = ", "))
}

msgs <- mclapply(jobs, run_job, mc.cores = N_CORES, mc.preschedule = FALSE)
for (m in msgs) message(m)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
