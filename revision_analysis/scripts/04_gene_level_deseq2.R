# Gene-level test of a sex difference in the response to vaccination (DESeq2, one model per arm).
#
# For each gene: is the change from D0 to D1 after vaccination different between females and males?
#
# Model, fitted separately in SI and LI on all three vaccinations (user-supplied design matrix):
#   counts ~ participant                 blocking factor: D1 is compared with D0 within each participant
#          + vaccination [T2, T3]        D0 level of the second and third vaccination vs the first
#          + female : vaccination : D1   D1 vs D0 change in females, per vaccination
#          + male   : vaccination : D1   D1 vs D0 change in males, per vaccination
# Sex difference = female change minus male change (log2 fold change), per vaccination, tested
# with a Wald test on the contrast female_Tk - male_Tk. Sex itself is constant within a participant
# and is absorbed by the participant terms.
#
# Three versions of the model:
#   basic          as above
#   age            + age x D1 (the D1 change may depend on age; age itself is absorbed by the
#                  participant terms)
#   age_cells      + age x D1 + estimated neutrophil and monocyte proportion of each sample
#                  (CIBERSORTx). Lymphocytes are largely the complement of neutrophils (r = -0.9).
#                  Because the proportions shift at D1 as part of the response, this version
#                  removes most of the D1 response and is read as the response not explained by
#                  changes in blood-cell composition.
#
# Genes: at least 10 reads in at least half of the female and half of the male samples of the arm
# (15,507 genes in SI, 15,681 in LI). BH correction over genes within each estimate.
#
# Inputs:  results/prepared/counts_raw_HD.csv.gz, analysis_samples_vaccination.csv,
#          cell_fractions_cibersortx.csv
# Outputs: results/gene_level_deseq2/gene_level_results.csv.gz (all genes, all estimates),
#          gene_level_summary.csv (genes at FDR < 0.05), design_and_session.txt
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "paths.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(data.table)
})
register(MulticoreParam(as.integer(Sys.getenv("THREADS", "16"))))
OD <- res_dir("gene_level_deseq2")

cnt <- fread(file.path(PREP, "counts_raw_HD.csv.gz"), data.table = FALSE)
rownames(cnt) <- cnt$gene; cnt <- as.matrix(cnt[, -1])
meta <- read.csv(file.path(PREP, "analysis_samples_vaccination.csv"), row.names = 1)
cells <- read.csv(file.path(PREP, "cell_fractions_cibersortx.csv"), row.names = 1)[rownames(meta), ]

design_matrix <- function(m, model) {
  subj <- model.matrix(~ 0 + factor(Subject), m)
  colnames(subj) <- paste0("participant_", sort(unique(m$Subject)))
  D1 <- as.numeric(m$Day == "D1"); fem <- as.numeric(m$Sex == "Female")
  cols <- list(base_T2 = as.numeric(m$Vaccination == "T2"),
               base_T3 = as.numeric(m$Vaccination == "T3"))
  for (v in c("T1", "T2", "T3")) {
    cols[[paste0("female_", v)]] <- fem * (m$Vaccination == v) * D1
    cols[[paste0("male_", v)]] <- (1 - fem) * (m$Vaccination == v) * D1
  }
  mm <- cbind(subj, do.call(cbind, cols))
  if (model != "basic") mm <- cbind(mm, age_x_D1 = as.numeric(scale(m$Age)) * D1)
  if (model == "age_cells")
    mm <- cbind(mm, neutrophils = m$Neutrophil - mean(m$Neutrophil),
                monocytes = m$Monocyte - mean(m$Monocyte))
  stopifnot(qr(mm)$rank == ncol(mm))
  mm
}

tidy <- function(r, arm, model, vac, what) {
  d <- as.data.frame(r)
  data.frame(gene = rownames(d), arm = arm, model = model, vaccination = vac, estimate = what,
             baseMean = d$baseMean, log2FC = d$log2FoldChange, SE = d$lfcSE,
             CI_low = d$log2FoldChange - 1.96 * d$lfcSE, CI_high = d$log2FoldChange + 1.96 * d$lfcSE,
             pvalue = d$pvalue, padj = d$padj)
}

out <- list(); summ <- list()
for (arm in c("SI", "LI")) {
  m <- cbind(meta[meta$Interval == arm, ], cells[meta$Interval == arm, c("Neutrophil", "Monocyte")])
  cm <- cnt[, rownames(m)]
  keep <- Reduce(`&`, lapply(c("Female", "Male"), function(s)
    rowMeans(cm[, m$Sex == s] >= 10) >= 0.5))
  cm <- cm[keep, ]
  message(arm, ": ", nrow(cm), " genes, ", ncol(cm), " samples")
  for (model in c("basic", "age", "age_cells")) {
    mm <- design_matrix(m, model)
    dds <- DESeqDataSetFromMatrix(cm, m, ~ 1)
    dds <- DESeq(dds, full = mm, betaPrior = FALSE, parallel = TRUE, quiet = TRUE)
    rn <- resultsNames(dds)
    for (v in c("T1", "T2", "T3")) {
      fe <- paste0("female_", v); ma <- paste0("male_", v)
      w <- setNames(numeric(length(rn)), rn); w[fe] <- 1; w[ma] <- -1
      res <- rbind(tidy(results(dds, name = fe), arm, model, v, "female_D1_vs_D0"),
                   tidy(results(dds, name = ma), arm, model, v, "male_D1_vs_D0"),
                   tidy(results(dds, contrast = w), arm, model, v, "sex_difference_F_minus_M"))
      out[[paste(arm, model, v)]] <- res
      for (e in unique(res$estimate)) {
        d <- res[res$estimate == e, ]
        summ[[paste(arm, model, v, e)]] <- data.frame(
          arm = arm, model = model, vaccination = v, estimate = e,
          n_female = length(unique(m$Subject[m$Sex == "Female"])),
          n_male = length(unique(m$Subject[m$Sex == "Male"])), n_genes_tested = nrow(d),
          n_FDR05 = sum(d$padj < 0.05, na.rm = TRUE),
          n_FDR05_up = sum(d$padj < 0.05 & d$log2FC > 0, na.rm = TRUE),
          n_FDR05_down = sum(d$padj < 0.05 & d$log2FC < 0, na.rm = TRUE),
          n_FDR10 = sum(d$padj < 0.10, na.rm = TRUE),
          median_SE = median(d$SE, na.rm = TRUE))
      }
    }
    message("  ", model, " done")
  }
}
res <- do.call(rbind, out)
write.csv(res, gzfile(file.path(OD, "gene_level_results.csv.gz")), row.names = FALSE)
S <- do.call(rbind, summ); rownames(S) <- NULL
write.csv(S, file.path(OD, "gene_level_summary.csv"), row.names = FALSE)
print(S[S$estimate == "sex_difference_F_minus_M", ])
writeLines(c("Design columns (SI, age + cells):",
             paste(colnames(design_matrix(cbind(meta[meta$Interval == "SI", ],
                   cells[meta$Interval == "SI", c("Neutrophil", "Monocyte")]), "age_cells")), collapse = ", "),
             "", capture.output(sessionInfo())), file.path(OD, "design_and_session.txt"))
