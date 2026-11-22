# 09 Weighted gene co-expression network analysis (WGCNA), per arm
#
# Input per arm (08_wgcna_input.py): VST values of the 10,000 most variable genes for the
# 80 samples of the arm (10 participants x D0/D1 of four injections) and the visit_sex label.
#
# Steps, per arm separately:
#   1. Sample clustering (hclust, average linkage on Euclidean distances) for outlier
#      inspection (plot only; no sample is removed).
#   2. Soft-thresholding power: pickSoftThreshold over powers 1-20 with default settings
#      (unsigned network, scale-free fit R^2 cut-off 0.85); the estimated power
#      (sft$powerEstimate) is used, 6 if none is found. Both arms give power 3.
#   3. Network and modules: blockwiseModules(power, TOMType = "unsigned", minModuleSize = 30,
#      reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
#      pamRespectsDendro = FALSE), other arguments at their defaults (maxBlockSize = 5000, so
#      the 10,000 genes are split into blocks by projective k-means pre-clustering;
#      randomSeed = 54321; deepSplit = 2; corType = "pearson"). Module labels are converted to
#      colours with labels2colors; "grey" holds the unassigned genes. SI gives 20 modules and
#      LI 24, grey included.
#   4. Module-trait relationships: module eigengenes (moduleEigengenes) are correlated
#      (Pearson) with the 0/1 indicator of each of the 16 visit_sex groups
#      (model.matrix(~ timepoint - 1)); p values from corPvalueStudent with n = 80. The p values
#      are not adjusted for multiple testing.
#   5. Export: module-trait correlations and p values, module assignment of every gene, and the
#      genes per module (JSON; converted to CSV and module codes M1... in 10_wgcna_modules.py).
#
# Gene identifiers: the counts are read with read.csv defaults (check.names = TRUE), so in the
# module tables gene symbols are R syntactic names ("-" written as ".").
#
# Outputs:
#   results/wgcna/wgcna_results_<ARM>/<ARM>_soft_threshold_analysis.csv
#   results/wgcna/wgcna_results_<ARM>/<ARM>_moduleTraitCor.csv, <ARM>_moduleTraitPvalue.csv
#   results/wgcna/wgcna_results_<ARM>/<ARM>_module_assignments.csv
#   results/wgcna/wgcna_results_<ARM>/<ARM>_module_genes.json
#   figures/wgcna/<ARM>_sample_dendrogram.pdf, <ARM>_scale_independence.pdf,
#                 <ARM>_module_dendrogram.pdf
# The topological overlap matrices are written to a temporary folder and deleted.
#
# Run from original_analysis/:  Rscript 09_wgcna.R

suppressPackageStartupMessages({
  library(WGCNA)
  library(jsonlite)
})

run_wgcna <- function(counts_file, metadata_file, output_dir, fig_dir, condition) {
  cor <- WGCNA::cor
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  tom_dir <- file.path(tempdir(), paste0("tom_", condition))
  dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)

  counts <- read.csv(counts_file, row.names = 1)
  metadata <- read.csv(metadata_file, row.names = 1)

  common_samples <- intersect(rownames(counts), rownames(metadata))
  datExpr <- counts[common_samples, ]
  metadata <- metadata[common_samples, ]

  # 1. sample clustering
  sampleTree <- hclust(dist(datExpr), method = "average")
  pdf(file = file.path(fig_dir, paste0(condition, "_sample_dendrogram.pdf")), width = 12, height = 9)
  par(cex = 0.6)
  par(mar = c(0, 4, 2, 0))
  plot(sampleTree, main = paste("Sample clustering to detect outliers -", condition),
       sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()

  # 2. soft-thresholding power
  powers <- c(1:20)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

  pdf(file = file.path(fig_dir, paste0(condition, "_scale_independence.pdf")), width = 9, height = 5)
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
       type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
  dev.off()

  write.csv(sft$fitIndices, file = file.path(output_dir, paste0(condition, "_soft_threshold_analysis.csv")))

  power <- sft$powerEstimate
  if (is.na(power)) power <- 6
  message(condition, ": soft-thresholding power ", power)

  # 3. network construction and module detection
  net <- blockwiseModules(datExpr,
                          power = power,
                          TOMType = "unsigned",
                          minModuleSize = 30,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          numericLabels = TRUE,
                          pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = file.path(tom_dir, paste0(condition, "_TOM")),
                          verbose = 3)
  unlink(tom_dir, recursive = TRUE)

  moduleColors <- labels2colors(net$colors)
  pdf(file = file.path(fig_dir, paste0(condition, "_module_dendrogram.pdf")), width = 12, height = 9)
  # with several blocks, the dendrogram of the first block is shown
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  # 4. module-trait relationships
  timepoint <- metadata$Timepoint
  MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  moduleTraitCor <- WGCNA::cor(MEs, model.matrix(~ timepoint - 1), use = "p")
  moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(datExpr))

  # 5. export
  write.csv(moduleTraitCor, file = file.path(output_dir, paste0(condition, "_moduleTraitCor.csv")))
  write.csv(moduleTraitPvalue, file = file.path(output_dir, paste0(condition, "_moduleTraitPvalue.csv")))
  module_assignments <- data.frame(gene = colnames(counts), module = moduleColors)
  write.csv(module_assignments, file = file.path(output_dir, paste0(condition, "_module_assignments.csv")))
  modules <- unique(moduleColors)
  module_genes <- lapply(modules, function(m) colnames(counts)[moduleColors == m])
  names(module_genes) <- modules
  write_json(module_genes, file.path(output_dir, paste0(condition, "_module_genes.json")))

  message(condition, ": ", length(modules), " modules (incl. grey)")
  invisible(list(moduleColors = moduleColors, power = power))
}

for (condition in c("SI", "LI")) {
  run_wgcna(counts_file = file.path("results/wgcna", paste0("counts_", condition, ".csv")),
            metadata_file = file.path("results/wgcna", paste0("metadata_", condition, ".csv")),
            output_dir = file.path("results/wgcna", paste0("wgcna_results_", condition)),
            fig_dir = "figures/wgcna",
            condition = condition)
}
writeLines(capture.output(sessionInfo()), "results/wgcna/sessionInfo.txt")
