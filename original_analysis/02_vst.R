# 02 Variance stabilising transformation
#
# DESeq2::varianceStabilizingTransformation (blind = TRUE) on the raw (integer) Salmon gene
# counts of all 244 samples of the high-dose cohort, placebo arm and excluded samples included,
# all genes. Every later step that works on VST values (PCA, DEG clustering, WGCNA) takes the
# samples it needs from this matrix.
#
# Input : ../data/counts_raw.csv.gz  (samples x gene symbols)
# Output: results/vst/vst_counts.csv.gz  (samples x genes)
#
# Run from original_analysis/:  Rscript 02_vst.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

stopifnot(file.exists("../data/counts_raw.csv.gz"))
out_dir <- "results/vst"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- fread("../data/counts_raw.csv.gz", header = TRUE)
samples <- counts[[1]]
counts <- as.matrix(counts[, -1])
rownames(counts) <- samples

vsd <- varianceStabilizingTransformation(t(counts), blind = TRUE)

out <- data.table(sample = samples, t(vsd))
# fwrite of this data.table build cannot compress; write plain CSV and gzip it
out_file <- file.path(out_dir, "vst_counts.csv")
fwrite(out, out_file)
unlink(paste0(out_file, ".gz"))
stopifnot(system2("gzip", out_file) == 0)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
