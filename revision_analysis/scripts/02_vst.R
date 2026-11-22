# Variance-stabilising transformation (DESeq2, blind) of the 138 analysis samples.
# Genes with at least 10 reads in at least 3 samples. Used for the signature scores, the
# variance partitioning and the per-gene linear models.
#
# Input:  results/prepared/counts_raw_HD.csv.gz, analysis_samples_vaccination.csv
# Output: results/prepared/vst_vaccination_samples.csv.gz (genes x samples)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "paths.R"))
suppressPackageStartupMessages({ library(DESeq2); library(data.table) })

counts <- fread(file.path(PREP, "counts_raw_HD.csv.gz"), data.table = FALSE)
rownames(counts) <- counts$gene; counts <- as.matrix(counts[, -1])
meta <- read.csv(file.path(PREP, "analysis_samples_vaccination.csv"), row.names = 1)
counts <- counts[, rownames(meta)]
counts <- counts[rowSums(counts >= 10) >= 3, ]
vst <- assay(vst(DESeqDataSetFromMatrix(counts, meta, ~ 1), blind = TRUE))
write.csv(vst, gzfile(file.path(PREP, "vst_vaccination_samples.csv.gz")))
message(sprintf("VST: %d genes x %d samples", nrow(vst), ncol(vst)))
