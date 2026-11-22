# 11 GO enrichment of the WGCNA modules
#
# For every module (codes M1... from 10_wgcna_modules.py, grey included):
#   * gene symbols are mapped to Entrez IDs with org.Hs.eg.db (mapIds, multiVals = "first";
#     unmapped symbols are dropped). Symbols that contain "-" are stored as R names with "."
#     (see 09_wgcna.R) and do not map;
#   * clusterProfiler::enrichGO(ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
#     qvalueCutoff = 0.2), universe = all annotated genes (default);
#   * if terms remain, redundant terms are collapsed with simplify(cutoff = 0.7,
#     by = "p.adjust"); if terms remain after that, the full enrichGO table (before
#     simplification) is written and the top 10 simplified terms are shown in a dot plot
#     (font size 12), one page per module. Modules without significant terms get a page
#     stating this.
#
# Requires clusterProfiler (version 4.2.0 was used for the paper) and org.Hs.eg.db.
#
# Outputs:
#   results/wgcna/wgcna_results_<ARM>/codenames/enrichment_tables/<ARM>_<M>_enrichment.csv
#   figures/wgcna/<ARM>_module_enrichment_codenames.pdf
#
# Run from original_analysis/:  Rscript 11_wgcna_module_enrichment.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(jsonlite)
})

create_empty_plot <- function(module) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("No significant enrichment found for ", module), size = 6) +
    theme_void() +
    ggtitle(paste("Module:", module)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
}

dir.create("figures/wgcna", recursive = TRUE, showWarnings = FALSE)

for (condition in c("SI", "LI")) {
  result_dir <- file.path("results/wgcna", paste0("wgcna_results_", condition), "codenames")
  table_dir <- file.path(result_dir, "enrichment_tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  module_genes <- fromJSON(file.path(result_dir, paste0(condition, "_module_genes_codenames.json")))
  pdf(file.path("figures/wgcna", paste0(condition, "_module_enrichment_codenames.pdf")),
      width = 7, height = 5)

  for (module in names(module_genes)) {
    message("Processing ", condition, " ", module)
    entrez_ids <- mapIds(org.Hs.eg.db, keys = module_genes[[module]], column = "ENTREZID",
                         keytype = "SYMBOL", multiVals = "first")
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]

    ego <- tryCatch(
      enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2),
      error = function(e) { message("enrichGO failed for ", module, ": ", e$message); NULL })

    if (is.null(ego) || nrow(ego@result) == 0) {
      print(create_empty_plot(module))
      next
    }

    ego_simp <- simplify(ego, cutoff = 0.7, by = "p.adjust")
    if (nrow(ego_simp@result) == 0) {
      print(create_empty_plot(module))
      next
    }
    write.csv(ego@result, file.path(table_dir, paste0(condition, "_", module, "_enrichment.csv")))
    p <- tryCatch(
      dotplot(ego_simp, showCategory = min(10, nrow(ego_simp@result)), font.size = 12) +
        ggtitle(module) +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              axis.text.y = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 12)),
      error = function(e) create_empty_plot(module))
    print(p)
  }
  dev.off()
}
