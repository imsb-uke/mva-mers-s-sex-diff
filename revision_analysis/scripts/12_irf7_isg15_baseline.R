# IRF7 and ISG15: sex difference at D1 for the same D0 level (ANCOVA), and expression levels.
#
# Data: log2(TPM + 1) per participant (results/participant_level/irf7_isg15_per_participant.csv,
# from 11_irf7_isg15.py)
# and VST values (results/prepared/vst_vaccination_samples.csv.gz) as a check.
#
# D0-adjusted comparison, per gene, arm and vaccination (one row per participant, least squares):
#   D0           D1 ~ sex + D0
#   D0 + age     D1 ~ sex + D0 + age
# The sex coefficient is the female minus male difference in log2 expression at D1 for the same D0
# level (identical to the sex coefficient of change ~ sex + D0).
# p values: t test of the coefficient, and exact permutation over all assignments of the sex labels
# within the arm (SI 495, LI 165; two-sided, |t|). Holm over the 12 tests per model.
#
# Expression levels at each visit (T1 to T3, D0 and D1), female minus male, mixed model per gene and
# arm on all six visits (lme4/lmerTest, REML, Satterthwaite df):
#   expression ~ 0 + sex:vaccination:day (+ age) + (1 | participant)
# p values: t test of each contrast, and exact permutation of the sex labels (|t|). Holm over 24.
#
# Outputs: tables/irf7_isg15_d0_adjusted.csv, tables/irf7_isg15_levels.csv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "paths.R"))
suppressPackageStartupMessages({ library(data.table); library(lmerTest) })

pp <- read.csv(file.path(PART, "irf7_isg15_per_participant.csv"))
meta <- read.csv(file.path(PREP, "analysis_samples_vaccination.csv"), row.names = 1)
vst <- fread(file.path(PREP, "vst_vaccination_samples.csv.gz"), data.table = FALSE)
rownames(vst) <- vst[[1]]; vst <- as.matrix(vst[, -1])
age <- tapply(meta$Age, meta$Subject, `[`, 1)

# VST values per participant in the same layout
vv <- do.call(rbind, lapply(seq_len(nrow(pp)), function(i) {
  r <- pp[i, ]; m <- meta[meta$Subject == r$Subject & meta$Vaccination == r$vaccination, ]
  data.frame(D0 = vst[r$gene, rownames(m)[m$Day == "D0"]], D1 = vst[r$gene, rownames(m)[m$Day == "D1"]])
}))
dat <- list(`log2(TPM+1)` = pp, VST = transform(pp, D0 = vv$D0, D1 = vv$D1))

fit_sex <- function(d1, fem, X) {
  f <- lm.fit(cbind(1, fem, X), d1); df <- length(d1) - f$rank
  V <- sum(f$residuals^2) / df * chol2inv(qr.R(f$qr))
  c(b = f$coefficients[["fem"]], se = sqrt(V[2, 2]), df = df)
}

out <- list()
for (scale in names(dat)) for (mn in c("D0", "D0 + age")) {
  d <- dat[[scale]]
  for (g in unique(d$gene)) for (arm in c("SI", "LI")) for (v in c("T1", "T2", "T3")) {
    x <- d[d$gene == g & d$arm == arm & d$vaccination == v, ]
    fem <- as.numeric(x$sex == "Female")
    X <- cbind(D0 = x$D0); if (mn == "D0 + age") X <- cbind(X, age = age[as.character(x$Subject)])
    r <- fit_sex(x$D1, fem, X); t0 <- r[["b"]] / r[["se"]]
    perm <- combn(length(fem), sum(fem))
    tn <- apply(perm, 2, function(ix) { s <- numeric(length(fem)); s[ix] <- 1; q <- fit_sex(x$D1, s, X); q[["b"]] / q[["se"]] })
    tq <- qt(0.975, r[["df"]])
    out[[length(out) + 1]] <- data.frame(scale = scale, model = mn, gene = g, arm = arm, vaccination = v,
      n_female = sum(fem), n_male = sum(1 - fem),
      D1_difference_F_minus_M = r[["b"]], CI_low = r[["b"]] - tq * r[["se"]], CI_high = r[["b"]] + tq * r[["se"]],
      pvalue = 2 * pt(-abs(t0), r[["df"]]), p_permutation = mean(abs(tn) >= abs(t0) - 1e-12),
      n_permutations = ncol(perm))
  }
}
out <- do.call(rbind, out)
out$p_holm_12 <- ave(out$pvalue, out$scale, out$model, FUN = function(p) p.adjust(p, "holm"))
out$p_permutation_holm_12 <- ave(out$p_permutation, out$scale, out$model, FUN = function(p) p.adjust(p, "holm"))
write.csv(out, file.path(TAB, "irf7_isg15_d0_adjusted.csv"), row.names = FALSE)
print(out[, c(1:5, 8, 11, 12, 14, 15)], digits = 3)

# Levels at each visit, mixed model
lev <- list()
for (scale in names(dat)) for (mn in c("unadjusted", "age")) {
  d <- dat[[scale]]
  long <- rbind(transform(d, day = "D0", y = D0), transform(d, day = "D1", y = D1))
  long$age_c <- as.numeric(scale(age[as.character(long$Subject)]))
  for (g in unique(long$gene)) for (arm in c("SI", "LI")) {
    x <- long[long$gene == g & long$arm == arm, ]
    subj <- sort(unique(x$Subject)); sx <- x$sex[match(subj, x$Subject)]
    tv <- function(fem_subj) {
      x$sex2 <- ifelse(x$Subject %in% fem_subj, "F", "M")
      x$grp <- factor(paste(x$sex2, x$vaccination, x$day, sep = "_"))
      f <- suppressMessages(lmer(if (mn == "age") y ~ 0 + grp + age_c + (1 | Subject) else y ~ 0 + grp + (1 | Subject),
                                 data = x, REML = TRUE))
      cf <- fixef(f)
      sapply(c("T1_D0", "T1_D1", "T2_D0", "T2_D1", "T3_D0", "T3_D1"), function(k) {
        L <- setNames(numeric(length(cf)), names(cf)); L[paste0("grpF_", k)] <- 1; L[paste0("grpM_", k)] <- -1
        ct <- contest1D(f, L, confint = TRUE)
        c(ct$Estimate, ct$`t value`, ct$`Pr(>|t|)`, ct$lower, ct$upper)
      })
    }
    obs <- tv(subj[sx == "Female"])
    perm <- combn(subj, sum(sx == "Female"))
    tn <- apply(perm, 2, function(ix) tv(ix)[2, ])
    for (j in seq_len(ncol(obs))) {
      k <- strsplit(colnames(obs)[j], "_")[[1]]
      lev[[length(lev) + 1]] <- data.frame(scale = scale, model = mn, gene = g, arm = arm, vaccination = k[1], day = k[2],
        n_female = sum(sx == "Female"), n_male = sum(sx == "Male"), difference_F_minus_M = obs[1, j],
        CI_low = obs[4, j], CI_high = obs[5, j], pvalue = obs[3, j],
        p_permutation = mean(abs(tn[j, ]) >= abs(obs[2, j]) - 1e-12), n_permutations = ncol(perm))
    }
  }
}
lev <- do.call(rbind, lev)
lev$p_holm_24 <- ave(lev$pvalue, lev$scale, lev$model, FUN = function(p) p.adjust(p, "holm"))
write.csv(lev, file.path(TAB, "irf7_isg15_levels.csv"), row.names = FALSE)
print(lev[lev$scale == "log2(TPM+1)", c(2:6, 9, 12, 13, 15)], digits = 3)
