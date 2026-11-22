# WGCNA: multiple-testing correction of the module-trait correlations, repeated-measures module
# tests, and robustness of the networks to individual participants.
#
# Uses the WGCNA input and module assignments of the networks shown in the paper
# (original_analysis/results/wgcna; SI modules M1 to M20, LI M1 to M24; 10 participants x 8 visits
# per arm; unsigned network, minimum module size 30, merge height 0.25). Module codes M1 to Mk
# number the module colours in alphabetical order.
#
#  1. Recompute the eigengenes and module-trait correlations and compare them with the saved ones.
#  2. BH-FDR over the module x trait matrix shown in the paper (12 traits: T1 to T3 x D0/D1 x
#     female/male, placebo visit not shown), per arm.
#  3. Mixed model per module on the vaccination visits (T1 to T3, 60 samples per arm):
#       eigengene ~ vaccination * day + sex + sex:day + (1 | participant)
#     day (response) and sex:day (sex difference in the response), BH-FDR across modules.
#  4. Leave one participant out:
#     a) module-trait correlations recomputed without each participant;
#     b) the network rebuilt without each participant (same parameters and soft-threshold power);
#        each module is matched to its best module in the rebuilt network by Jaccard index.
#
# Outputs: results/wgcna_robustness/ (module_trait_with_FDR.csv,
#          module_mixed_model_sex_x_day.csv, loso_*.csv, check_recomputed_vs_saved.csv)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "paths.R"))
suppressPackageStartupMessages({
  library(WGCNA)
  library(lmerTest)
  library(parallel)
})
set.seed(42)
cor <- WGCNA::cor
allowWGCNAThreads(4)
od <- res_dir("wgcna_robustness")
WD <- WGCNA_ORIG

VACC <- list(SI = c(V1 = "T1", V2 = "T2", V4 = "T3"), LI = c(V1 = "T1", V3 = "T2", V4 = "T3"))
PUB_TRAITS <- list(  # 12 trait columns shown in the paper -> sample time point
  SI = c(V1D0 = "V1D0", V1D1 = "V1D1", V2D0 = "V2D0", V2D1 = "V2D1", V3D0 = "V4D0", V3D1 = "V4D1"),
  LI = c(V1D0 = "V1D0", V1D1 = "V1D1", V2D0 = "V3D0", V2D1 = "V3D1", V3D0 = "V4D0", V3D1 = "V4D1"))

net_params <- function(datExpr, power) blockwiseModules(
  datExpr, power = power, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0,
  mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 0)

jaccard_match <- function(orig, new) {
  # orig/new: named vectors gene -> module; returns best Jaccard for each original module
  mods <- setdiff(unique(orig), "grey")
  sapply(mods, function(m) {
    g <- names(orig)[orig == m]
    best <- 0
    for (n in setdiff(unique(new), c("grey", 0))) {
      h <- names(new)[new == n]
      best <- max(best, length(intersect(g, h)) / length(union(g, h)))
    }
    best
  })
}

all_fdr <- list(); all_lmm <- list(); all_loso_tr <- list(); all_loso_mod <- list(); checks <- list()

for (arm in c("SI", "LI")) {
  # read.csv with check.names = TRUE, as in the network construction
  datExpr <- read.csv(file.path(WD, sprintf("counts_%s.csv", arm)), row.names = 1)
  md <- read.csv(file.path(WD, sprintf("metadata_%s.csv", arm)), row.names = 1)
  datExpr <- datExpr[rownames(md), ]
  ma <- read.csv(file.path(WD, sprintf("wgcna_results_%s/%s_module_assignments.csv", arm, arm)),
                 row.names = 1)
  colors <- setNames(ma$module, ma$gene)
  stopifnot(identical(names(colors), colnames(datExpr)))
  cmap <- setNames(paste0("M", seq_along(sort(unique(colors)))), sort(unique(colors)))

  # 1. eigengenes and module-trait correlations
  MEs <- moduleEigengenes(datExpr, colors)$eigengenes
  colnames(MEs) <- cmap[sub("^ME", "", colnames(MEs))]
  tp <- md$Timepoint                                     # e.g. "V1D0_Female"
  traits <- model.matrix(~ tp - 1); colnames(traits) <- sub("^tp", "", colnames(traits))
  mtc <- cor(MEs, traits, use = "p")
  saved <- read.csv(file.path(WD, sprintf("wgcna_results_%s/%s_moduleTraitCor.csv", arm, arm)),
                    row.names = 1)
  rownames(saved) <- cmap[sub("^ME", "", rownames(saved))]
  colnames(saved) <- sub("^timepoint", "", colnames(saved))
  checks[[arm]] <- data.frame(arm = arm, max_abs_diff_cor_vs_saved =
                                max(abs(mtc[rownames(saved), colnames(saved)] - as.matrix(saved))))

  # 2. BH-FDR over the module-trait matrix shown in the paper (grey = unassigned genes, excluded)
  pub_cols <- unlist(lapply(names(PUB_TRAITS[[arm]]), function(k)
    paste0(PUB_TRAITS[[arm]][[k]], c("_Female", "_Male"))))
  pub_lab <- unlist(lapply(names(PUB_TRAITS[[arm]]), function(k) paste0(k, c("_Female", "_Male"))))
  mods <- setdiff(colnames(MEs), cmap[["grey"]])
  r <- mtc[mods, pub_cols]; p <- corPvalueStudent(r, nrow(datExpr))
  fdr <- matrix(p.adjust(p, "BH"), nrow(p), dimnames = dimnames(p))
  long <- data.frame(arm = arm, module = rep(mods, times = length(pub_cols)),
                     trait = rep(pub_lab, each = length(mods)),
                     original_trait = rep(pub_cols, each = length(mods)),
                     r = as.vector(r), p = as.vector(p), p_BH_matrix = as.vector(fdr))
  all_fdr[[arm]] <- long

  # 3. mixed models on vaccination visits
  vis <- substr(md$Timepoint, 1, 2)
  keep <- vis %in% names(VACC[[arm]])
  d <- data.frame(Subject = factor(md$Subject[keep]), Sex = factor(md$Sex[keep], c("Male", "Female")),
                  Day = factor(substr(md$Timepoint[keep], 3, 4), c("D0", "D1")),
                  Vaccination = factor(VACC[[arm]][vis[keep]]))
  lmm <- do.call(rbind, lapply(mods, function(m) {
    d$ME <- MEs[keep, m]
    f <- lmer(ME ~ Vaccination * Day + Sex + Sex:Day + (1 | Subject), data = d, REML = TRUE)
    an <- anova(f)
    cf <- summary(f)$coefficients
    ix <- grep("^(DayD1:SexFemale|SexFemale:DayD1)$", rownames(cf))
    data.frame(arm = arm, module = m,
               day_effect_T1 = cf["DayD1", "Estimate"], p_day = an["Day", "Pr(>F)"],
               sex_x_day = cf[ix, "Estimate"],
               sex_x_day_SE = cf[ix, "Std. Error"],
               p_sex_x_day = cf[ix, "Pr(>|t|)"],
               p_sex_baseline = an["Sex", "Pr(>F)"])
  }))
  lmm$padj_day <- p.adjust(lmm$p_day, "BH")
  lmm$padj_sex_x_day <- p.adjust(lmm$p_sex_x_day, "BH")
  lmm$padj_sex_baseline <- p.adjust(lmm$p_sex_baseline, "BH")
  lmm$n_female <- length(unique(d$Subject[d$Sex == "Female"]))
  lmm$n_male <- length(unique(d$Subject[d$Sex == "Male"]))
  all_lmm[[arm]] <- lmm

  # 4a. leave-one-participant-out module-trait correlations (eigengenes of the full network)
  subjects <- unique(md$Subject)
  lt <- do.call(rbind, lapply(subjects, function(s) {
    k <- md$Subject != s
    rr <- cor(MEs[k, mods], traits[k, pub_cols], use = "p")
    pp <- corPvalueStudent(rr, sum(k))
    data.frame(arm = arm, left_out = s, module = rep(mods, times = length(pub_cols)),
               trait = rep(pub_lab, each = length(mods)), r = as.vector(rr), p = as.vector(pp))
  }))
  all_loso_tr[[arm]] <- lt

  # 4b. network rebuilt without each participant
  sft <- read.csv(file.path(WD, sprintf("wgcna_results_%s/%s_soft_threshold_analysis.csv", arm, arm)))
  power <- sft$Power[which(sft$SFT.R.sq >= 0.85)[1]]
  if (is.na(power)) power <- 6
  full <- net_params(as.matrix(datExpr), power)
  full_col <- setNames(labels2colors(full$colors), colnames(datExpr))
  full_j <- jaccard_match(colors, full_col)
  cache <- file.path(od, sprintf("cache_loso_networks_%s.rds", arm))
  if (file.exists(cache)) {
    lm <- readRDS(cache)
  } else {
    lm <- mclapply(subjects, function(s) {
      k <- md$Subject != s
      n <- net_params(as.matrix(datExpr[k, ]), power)
      jaccard_match(colors, setNames(labels2colors(n$colors), colnames(datExpr)))
    }, mc.cores = as.integer(Sys.getenv("THREADS", "10")))
    saveRDS(lm, cache)
  }
  jm <- do.call(cbind, lm); colnames(jm) <- paste0("LOSO_", subjects)
  write.csv(data.frame(module = cmap[rownames(jm)], jm, check.names = FALSE),
            file.path(od, sprintf("loso_jaccard_per_subject_%s.csv", arm)), row.names = FALSE)
  mod_df <- data.frame(arm = arm, module = cmap[rownames(jm)], colour = rownames(jm),
                       n_genes = as.integer(table(colors)[rownames(jm)]),
                       power = power, jaccard_full_rerun = full_j[rownames(jm)],
                       jaccard_LOSO_median = apply(jm, 1, median),
                       jaccard_LOSO_min = apply(jm, 1, min),
                       n_LOSO_jaccard_ge_0.5 = rowSums(jm >= 0.5), n_LOSO = ncol(jm))
  all_loso_mod[[arm]] <- mod_df
  message(arm, " done")
}

fdr <- do.call(rbind, all_fdr)
write.csv(fdr, file.path(od, "module_trait_with_FDR.csv"), row.names = FALSE)
lmm <- do.call(rbind, all_lmm)
write.csv(lmm, file.path(od, "module_mixed_model_sex_x_day.csv"), row.names = FALSE)
lt <- do.call(rbind, all_loso_tr)
write.csv(lt, gzfile(file.path(od, "loso_module_trait.csv.gz")), row.names = FALSE)
# stability of the module-trait cells with p < 0.05
sig <- fdr[fdr$p < 0.05, ]
stab <- do.call(rbind, lapply(seq_len(nrow(sig)), function(i) {
  x <- lt[lt$arm == sig$arm[i] & lt$module == sig$module[i] & lt$trait == sig$trait[i], ]
  data.frame(sig[i, c("arm", "module", "trait", "r", "p", "p_BH_matrix")],
             LOSO_frac_p05 = mean(x$p < 0.05), LOSO_r_min = min(x$r), LOSO_r_max = max(x$r),
             LOSO_same_sign = all(sign(x$r) == sign(sig$r[i])))
}))
write.csv(stab, file.path(od, "loso_stability_nominal_cells.csv"), row.names = FALSE)
mod <- do.call(rbind, all_loso_mod)
write.csv(mod, file.path(od, "loso_module_reconstruction_jaccard.csv"), row.names = FALSE)
write.csv(do.call(rbind, checks), file.path(od, "check_recomputed_vs_saved.csv"), row.names = FALSE)

cat("\nchecks\n"); print(do.call(rbind, checks))
cat("\nnominal p<0.05 cells:", sum(fdr$p < 0.05), " BH<0.05:", sum(fdr$p_BH_matrix < 0.05), "\n")
print(table(fdr$arm, fdr$p_BH_matrix < 0.05))
print(stab)
print(lmm[order(lmm$p_sex_x_day), ][1:10, ])
print(mod)
