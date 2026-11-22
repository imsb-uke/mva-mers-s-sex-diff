# Gene-level linear mixed model with the participant as a random effect (variancePartition::dream).
#
# Compares females and males before vaccination (D0) and in the D1 vs D0 response, with the repeated
# samples of each participant modelled by a random intercept.
#
# Model per arm (all three vaccinations, 6 samples per participant):
#   log2 expression (voom with dream precision weights)
#     ~ 0 + sex:vaccination:day      one mean per sex, vaccination and day
#       + age
#       + (1 | participant)
# Contrasts per vaccination:
#   D0         female D0 - male D0
#   response   (female D1 - female D0) - (male D1 - male D0)
# The random intercept allows the D0 comparison between the sexes, which the DESeq2 model with a
# fixed term per participant (04_gene_level_deseq2.R) cannot make. For the response both models use
# the within-participant change. Genes: same filter as 04_gene_level_deseq2.R.
#
# Inputs:  results/prepared/counts_raw_HD.csv.gz, analysis_samples_vaccination.csv
# Outputs: results/gene_level_mixed/mixed_model_results.csv.gz, mixed_model_summary.csv,
#          design_and_session.txt
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "paths.R"))
suppressPackageStartupMessages({
  library(variancePartition)
  library(edgeR)
  library(BiocParallel)
  library(data.table)
})
param <- MulticoreParam(as.integer(Sys.getenv("THREADS", "16")))
OD <- res_dir("gene_level_mixed")

cnt <- fread(file.path(PREP, "counts_raw_HD.csv.gz"), data.table = FALSE)
rownames(cnt) <- cnt$gene; cnt <- as.matrix(cnt[, -1])
meta <- read.csv(file.path(PREP, "analysis_samples_vaccination.csv"), row.names = 1)

out <- list(); summ <- list()
for (arm in c("SI", "LI")) {
  m <- meta[meta$Interval == arm, ]
  m$grp <- factor(paste(m$Sex, m$Vaccination, m$Day, sep = "_"))
  m$age <- as.numeric(scale(m$Age))
  m$Subject <- factor(m$Subject)
  cm <- cnt[, rownames(m)]
  keep <- Reduce(`&`, lapply(c("Female", "Male"), function(s) rowMeans(cm[, m$Sex == s] >= 10) >= 0.5))
  dge <- calcNormFactors(DGEList(cm[keep, ]))
  form <- ~ 0 + grp + age + (1 | Subject)
  vobj <- voomWithDreamWeights(dge, form, m, BPPARAM = param)
  L <- list()
  for (v in c("T1", "T2", "T3")) {
    g <- function(s, d) paste0("grp", s, "_", v, "_", d)
    L[[paste0("D0_", v)]] <- paste(g("Female", "D0"), "-", g("Male", "D0"))
    L[[paste0("response_", v)]] <- sprintf("(%s - %s) - (%s - %s)", g("Female", "D1"), g("Female", "D0"),
                                           g("Male", "D1"), g("Male", "D0"))
  }
  Lm <- do.call(makeContrastsDream, c(list(formula = form, data = m), L))
  fit <- eBayes(dream(vobj, form, m, L = Lm, BPPARAM = param))
  for (k in names(L)) {
    tt <- topTable(fit, coef = k, number = Inf, sort.by = "none", confint = TRUE)
    d <- data.frame(gene = rownames(tt), arm = arm, contrast = sub("_T[123]$", "", k),
                    vaccination = sub("^.*_", "", k), log2FC = tt$logFC,
                    CI_low = tt$CI.L, CI_high = tt$CI.R,
                    pvalue = tt$P.Value, padj = tt$adj.P.Val)
    out[[paste(arm, k)]] <- d
    summ[[paste(arm, k)]] <- data.frame(arm = arm, contrast = d$contrast[1], vaccination = d$vaccination[1],
      n_genes_tested = nrow(d), n_FDR05 = sum(d$padj < 0.05), n_FDR10 = sum(d$padj < 0.10),
      pct_p05 = 100 * mean(d$pvalue < 0.05))
  }
  message(arm, " done")
}
write.csv(do.call(rbind, out), gzfile(file.path(OD, "mixed_model_results.csv.gz")), row.names = FALSE)
S <- do.call(rbind, summ); rownames(S) <- NULL
write.csv(S, file.path(OD, "mixed_model_summary.csv"), row.names = FALSE)
print(S)
writeLines(c("dream formula: ~ 0 + grp + age + (1 | Subject); grp = sex_vaccination_day", "",
             capture.output(sessionInfo())), file.path(OD, "design_and_session.txt"))
