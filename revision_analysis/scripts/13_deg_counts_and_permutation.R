# Per-sex DEG counts at equal group size, per-gene linear models of the response, and a
# permutation null for the number of genes with a sex difference.
#
# Data: VST expression (results/prepared/vst_vaccination_samples.csv.gz); response = D1 - D0 per
# participant and vaccination. Genes: same filter as 04_gene_level_deseq2.R (at least 10 reads in
# half of the female and half of the male samples of the arm).
#
# A  Genes changed D1 vs D0 within each sex: one-sample moderated t test (limma) of the response,
#    BH-FDR < 0.05. All participants, and every subset of the larger sex reduced to the size of the
#    smaller sex (SI: 4 of 8 females, 70 subsets; LI: 3 of 8 males, 56 subsets).
# B  Per-gene sex difference in the response, ordinary least squares (the D0 covariate is the D0
#    level of each gene):
#      unadjusted     response ~ sex
#      D0             response ~ sex + D0 level of the gene
#      D0 + age       response ~ sex + D0 level of the gene + age
# C  Permutation null. For each model in B: number of genes at p < 0.05, p < 0.01 and BH-FDR < 0.05,
#    0.10, 0.20 and 0.30, observed and after every assignment of the sex labels within the arm (SI 495,
#    LI 165). Permuting labels keeps the correlation between genes, so the null gives the number of
#    genes expected without a sex difference.
#    Direction: over the genes that change after vaccination in all participants together (one-sample
#    t test, BH-FDR < 0.05, computed without sex), the mean female minus male difference in the
#    response, signed by the direction of the response (positive = larger response in females).
#
# Outputs: tables/deg_counts_per_sex_equal_n.csv, tables/sex_difference_permutation_null.csv,
#          results/deg_counts_permutation/ (DEG counts of every subset, per-gene results of B,
#          permutation values of the unadjusted model)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "paths.R"))
suppressPackageStartupMessages({ library(limma); library(data.table) })
OD <- res_dir("deg_counts_permutation")

vst <- fread(file.path(PREP, "vst_vaccination_samples.csv.gz"), data.table = FALSE)
rownames(vst) <- vst[[1]]; vst <- as.matrix(vst[, -1])
cnt <- fread(file.path(PREP, "counts_raw_HD.csv.gz"), data.table = FALSE)
rownames(cnt) <- cnt$gene; cnt <- as.matrix(cnt[, -1])
meta <- read.csv(file.path(PREP, "analysis_samples_vaccination.csv"), row.names = 1)

CUT <- list(p05 = c("p", 0.05), p01 = c("p", 0.01), fdr05 = c("q", 0.05), fdr10 = c("q", 0.10),
            fdr20 = c("q", 0.20), fdr30 = c("q", 0.30))
count_cuts <- function(p) {
  q <- p.adjust(p, "BH")
  sapply(CUT, function(k) sum((if (k[1] == "p") p else q) < as.numeric(k[2])))
}

# OLS sex coefficient per gene: y (genes x n), s (n), Z common covariates (n x k), G optional
# gene-specific covariate (genes x n). Frisch-Waugh: residualise on Z, then on G per gene.
ols_sex <- function(y, s, Z, G = NULL) {
  P <- diag(length(s)) - Z %*% solve(crossprod(Z), t(Z))
  yr <- y %*% P; sr <- as.numeric(P %*% s); df <- length(s) - ncol(Z) - 1
  if (is.null(G)) {
    b <- as.numeric(yr %*% sr) / sum(sr^2)
    res <- yr - outer(b, sr); se <- sqrt(rowSums(res^2) / df / sum(sr^2))
  } else {
    gr <- G %*% P; df <- df - 1
    gg <- rowSums(gr^2)
    yg <- yr - gr * (rowSums(yr * gr) / gg)
    sg <- matrix(sr, nrow(gr), length(sr), byrow = TRUE) - gr * (as.numeric(gr %*% sr) / gg)
    ss <- rowSums(sg^2); b <- rowSums(sg * yg) / ss
    res <- yg - sg * b; se <- sqrt(rowSums(res^2) / df / ss)
  }
  t <- b / se
  list(b = b, se = se, p = 2 * pt(-abs(t), df))
}

A <- list(); C <- list(); Bg <- list(); D <- list()
for (arm in c("SI", "LI")) {
  m <- meta[meta$Interval == arm, ]
  cm <- cnt[, rownames(m)]
  keep <- Reduce(`&`, lapply(c("Female", "Male"), function(s) rowMeans(cm[, m$Sex == s] >= 10) >= 0.5))
  genes <- intersect(rownames(cm)[keep], rownames(vst))
  for (v in c("T1", "T2", "T3")) {
    mv <- m[m$Vaccination == v, ]
    subj <- sort(unique(mv$Subject))
    s0 <- sapply(subj, function(i) rownames(mv)[mv$Subject == i & mv$Day == "D0"])
    s1 <- sapply(subj, function(i) rownames(mv)[mv$Subject == i & mv$Day == "D1"])
    d0 <- vst[genes, s0]; y <- vst[genes, s1] - d0
    sex <- mv$Sex[match(subj, mv$Subject)]; age <- mv$Age[match(subj, mv$Subject)]
    fem <- as.numeric(sex == "Female"); nf <- sum(fem); nm <- sum(1 - fem)

    # A: per-sex DEG counts, full and equal n
    ndeg <- function(cols) {
      fit <- eBayes(lmFit(y[, cols, drop = FALSE], matrix(1, length(cols), 1)))
      sum(p.adjust(fit$p.value[, 1], "BH") < 0.05)
    }
    nmin <- min(nf, nm)
    for (sx in c("Female", "Male")) {
      idx <- which(sex == sx)
      A[[length(A) + 1]] <- data.frame(arm = arm, vaccination = v, sex = sx, n = length(idx),
                                       sample = "full", subset = 0, n_DEG = ndeg(idx))
      if (length(idx) > nmin) {
        cb <- combn(idx, nmin)
        A[[length(A) + 1]] <- data.frame(arm = arm, vaccination = v, sex = sx, n = nmin,
                                         sample = "equal n", subset = seq_len(ncol(cb)),
                                         n_DEG = apply(cb, 2, ndeg))
      }
    }

    # B and C
    Z1 <- matrix(1, length(subj), 1); Za <- cbind(1, scale(age))
    models <- list(unadjusted = function(s) ols_sex(y, s, Z1),
                   D0 = function(s) ols_sex(y, s, Z1, d0),
                   `D0 + age` = function(s) ols_sex(y, s, Za, d0))
    resp <- p.adjust(eBayes(lmFit(y, Z1))$p.value[, 1], "BH") < 0.05
    rdir <- sign(rowMeans(y))
    perms <- combn(length(subj), nf)
    for (mn in names(models)) {
      r <- models[[mn]](fem)
      Bg[[paste(arm, v, mn)]] <- data.frame(gene = genes, arm = arm, vaccination = v, model = mn,
        difference_F_minus_M = r$b, SE = r$se, pvalue = r$p, padj = p.adjust(r$p, "BH"), responsive = resp)
      obs <- c(count_cuts(r$p), direction = mean((r$b * rdir)[resp]))
      null <- t(apply(perms, 2, function(ix) {
        s <- numeric(length(subj)); s[ix] <- 1
        rr <- models[[mn]](s)
        c(count_cuts(rr$p), direction = mean((rr$b * rdir)[resp]))
      }))
      for (k in names(obs)) {
        C[[length(C) + 1]] <- data.frame(arm = arm, vaccination = v, model = mn, measure = k,
          observed = obs[[k]], null_median = median(null[, k]), null_q95 = quantile(null[, k], 0.95),
          null_max = max(null[, k]), n_permutations = nrow(null),
          p_permutation = mean(null[, k] >= obs[[k]]), n_genes = length(genes), n_responsive = sum(resp))
      }
      if (mn == "unadjusted") D[[paste(arm, v)]] <- data.frame(arm = arm, vaccination = v, null)
    }
    message(arm, " ", v, " done")
  }
}
A <- do.call(rbind, A); C <- do.call(rbind, C); rownames(C) <- NULL
write.csv(A, file.path(OD, "deg_counts_equal_n_all_subsets.csv"), row.names = FALSE)
As <- aggregate(n_DEG ~ arm + vaccination + sex + n + sample, A, function(x) c(median = median(x), min = min(x), max = max(x)))
As <- do.call(data.frame, As); names(As) <- sub("n_DEG\\.", "n_DEG_", names(As))
write.csv(As, file.path(TAB, "deg_counts_per_sex_equal_n.csv"), row.names = FALSE)
write.csv(C, file.path(TAB, "sex_difference_permutation_null.csv"), row.names = FALSE)
write.csv(do.call(rbind, Bg), gzfile(file.path(OD, "gene_level_D0_adjusted.csv.gz")), row.names = FALSE)
write.csv(do.call(rbind, D), file.path(OD, "permutation_null_unadjusted.csv"), row.names = FALSE)
print(As)
print(C[C$measure %in% c("p05", "fdr05", "fdr10", "fdr20", "direction"), c(1:5, 6, 7, 10)], digits = 3)
