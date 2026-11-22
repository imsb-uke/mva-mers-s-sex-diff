# Antibody titres by sex: comparison at each time point with adjustment for arm and age, per-arm
# analysis, and a longitudinal mixed model. Binding (ELISA S1) and neutralising (VSV-PRNT50)
# antibodies.
#
# Data: data/antibody_titres.xlsx (sheet "Antibody data", 125 mITT participants, 10 time points;
#       rows with unit IU/mL are used) and data/participants_mITT.xlsx (sex, age, dose, interval).
# Outcome: log10 titre (IU/mL). Values at the lower limit of the assay (ELISA 4.6 IU/mL, PRNT50
# 16 IU/mL) are used as given, with a censored (Tobit) model as sensitivity analysis.
# Arms: HD SI, HD LI, LD SI, LD LI (vaccinated). Placebo participants appear in the population
# table only.
# Time points: study days 0, 14, 28, 42, 56, 70, 84, 168, 224, 252. Vaccinations: SI days 0, 28 and
# 224 (placebo injection day 56), LI days 0, 56 and 224 (placebo injection day 28). Day 252 is 28 days
# after the third vaccination (V3D28 in the paper, which counts vaccinations; the trial data label it
# V4D28 because they count the placebo visit).
#
# A  Population: participants per arm and sex, age, samples per time point.
# B  Geometric mean titre (GMT) per arm, sex and time point, with 95% CI.
# C  Each time point separately, vaccinated arms, male / female geometric mean ratio (GMR):
#      unadjusted           log10 titre ~ sex
#      arm                  log10 titre ~ sex + arm
#      arm + age            log10 titre ~ sex + arm + age
#      arm + age + D0       log10 titre ~ sex + arm + age + log10 day-0 titre
#      Tobit                as arm + age, left-censored at the lower limit (survival::survreg)
#      per arm              log10 titre ~ sex, each arm separately
#    plus the sex x arm interaction (F test) and a Wilcoxon test of the unadjusted comparison.
#    Holm over the 10 time points per assay, model and arm. The D0 model is not fitted at day 0;
#    models that cannot be fitted (all values at the lower limit) are left empty.
# D  Longitudinal model, vaccinated arms, all time points (lme4/lmerTest, REML, Satterthwaite df):
#      log10 titre ~ arm * time + sex * time + age + (1 | participant)
#    time as a factor (10 levels). Sex x time: F test (Satterthwaite) and likelihood-ratio test (ML).
#    Sex x arm x time: likelihood-ratio test of adding sex:arm and sex:arm:time.
#    Male / female GMR at each time point (emmeans), Holm over the 10 time points.
#    Sensitivity: time points after day 0 only, with the log10 day-0 titre as covariate.
#
# Outputs: tables/antibody_population.csv, antibody_gmt.csv, antibody_sex_difference_by_day.csv,
#          antibody_sex_difference_mixed_model.csv, antibody_model_tests.csv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "paths.R"))
suppressPackageStartupMessages({ library(readxl); library(lmerTest); library(emmeans); library(survival) })

ab <- as.data.frame(read_excel(file.path(DATA, "antibody_titres.xlsx"), sheet = "Antibody data"))
mitt <- as.data.frame(read_excel(file.path(DATA, "participants_mITT.xlsx")))
ab <- ab[ab$Unit == "IU/mL", ]
ab$Subject <- ab$`Subject ID`; ab$day <- ab$T
ab <- merge(ab, mitt[, c("Subject", "Sex", "Age", "Dose", "Interval")], by = "Subject")
stopifnot(!anyNA(ab$Sex))
mitt$arm <- ifelse(mitt$Dose == "P", "Placebo", paste(mitt$Dose, mitt$Interval))
ab$arm <- ifelse(ab$Dose == "P", "Placebo", paste(ab$Dose, ab$Interval))
ab$assay <- c(ELISA_S1 = "Binding (ELISA S1)", VSV_PRNT50 = "Neutralising (PRNT50)")[ab$Assay]
LLOQ <- c(`Binding (ELISA S1)` = 4.6, `Neutralising (PRNT50)` = 16)
VISIT <- c(`0` = "V1D0", `14` = "V1D14", `28` = "V2D0", `42` = "V2D14", `56` = "V3D0", `70` = "V3D14",
           `84` = "V3D28", `168` = "V3D112", `224` = "V4D0", `252` = "V4D28")
ARMS <- c("HD SI", "HD LI", "LD SI", "LD LI")
ab$visit <- factor(VISIT[as.character(ab$day)], levels = VISIT)
ab$y <- log10(ab$Value)
ab$Sex <- factor(ab$Sex, levels = c("Female", "Male"))
ab$age_c <- ab$Age - mean(mitt$Age)
d0 <- ab[ab$day == 0, c("Subject", "assay", "y")]; names(d0)[3] <- "y0"
ab <- merge(ab, d0, by = c("Subject", "assay"), all.x = TRUE)
vac <- ab[ab$arm %in% ARMS, ]; vac$arm <- factor(vac$arm, levels = ARMS)

# A: population
pop <- do.call(rbind, lapply(c(ARMS, "Placebo"), function(a) do.call(rbind, lapply(c("Female", "Male"), function(s) {
  m <- mitt[mitt$arm == a & mitt$Sex == s, ]
  x <- ab[ab$arm == a & ab$Sex == s, ]
  data.frame(arm = a, sex = s, mITT = nrow(m), age_median = median(m$Age), age_range = paste(range(m$Age), collapse = " to "),
             binding_day252 = sum(x$assay == names(LLOQ)[1] & x$day == 252),
             neutralising_day252 = sum(x$assay == names(LLOQ)[2] & x$day == 252),
             binding_day84 = sum(x$assay == names(LLOQ)[1] & x$day == 84),
             neutralising_day84 = sum(x$assay == names(LLOQ)[2] & x$day == 84),
             binding_samples = sum(x$assay == names(LLOQ)[1]), neutralising_samples = sum(x$assay == names(LLOQ)[2]))
}))))
write.csv(pop, file.path(TAB, "antibody_population.csv"), row.names = FALSE)

# B: GMT per arm, sex, time point
gm <- aggregate(y ~ arm + Sex + assay + visit + day, ab, function(v) c(n = length(v), m = mean(v), se = sd(v) / sqrt(length(v))))
gm <- do.call(data.frame, gm)
tq <- qt(0.975, pmax(gm$y.n - 1, 1))
gmt <- data.frame(gm[, c("arm", "Sex", "assay", "visit", "day")], n = gm$y.n, GMT = 10^gm$y.m,
                  CI_low = 10^(gm$y.m - tq * gm$y.se), CI_high = 10^(gm$y.m + tq * gm$y.se))
gmt$pct_above_LLOQ <- mapply(function(a, s, as, dd) 100 * mean(ab$Value[ab$arm == a & ab$Sex == s & ab$assay == as & ab$day == dd] > LLOQ[[as]]),
                             gmt$arm, as.character(gmt$Sex), gmt$assay, gmt$day)
write.csv(gmt[order(gmt$assay, gmt$arm, gmt$Sex, gmt$day), ], file.path(TAB, "antibody_gmt.csv"), row.names = FALSE)

# C: each time point separately
gmr <- function(est, lo, hi, p, ...) data.frame(..., GMR_male_vs_female = 10^est, CI_low = 10^lo, CI_high = 10^hi, pvalue = p)
lm_sex <- function(f, x, ...) {
  nn <- list(n_female = sum(x$Sex == "Female"), n_male = sum(x$Sex == "Male"))
  r <- tryCatch({
    fit <- lm(f, x); ci <- confint(fit)["SexMale", ]; s <- summary(fit)$coefficients["SexMale", ]
    c(s[1], ci[1], ci[2], s[4])
  }, error = function(e) rep(NA, 4), warning = function(w) rep(NA, 4))
  do.call(gmr, c(list(r[1], r[2], r[3], r[4], ...), nn))
}
C <- list(); E <- list()
DAYS <- setNames(paste("day", names(VISIT)), names(VISIT)); DAYS["252"] <- "day 252 (V3D28)"
for (as in names(LLOQ)) for (dd in names(DAYS)) {
  x <- vac[vac$assay == as & vac$day == as.numeric(dd), ]
  C[[length(C) + 1]] <- lm_sex(y ~ Sex, x, assay = as, timepoint = DAYS[[dd]], model = "unadjusted", arm = "all vaccinated")
  C[[length(C) + 1]] <- lm_sex(y ~ Sex + arm, x, assay = as, timepoint = DAYS[[dd]], model = "arm", arm = "all vaccinated")
  C[[length(C) + 1]] <- lm_sex(y ~ Sex + arm + age_c, x, assay = as, timepoint = DAYS[[dd]], model = "arm + age", arm = "all vaccinated")
  if (dd != "0")
    C[[length(C) + 1]] <- lm_sex(y ~ Sex + arm + age_c + y0, x[!is.na(x$y0), ], assay = as, timepoint = DAYS[[dd]], model = "arm + age + D0", arm = "all vaccinated")
  s <- tryCatch(summary(survreg(Surv(y, y > log10(LLOQ[[as]]) + 1e-9, type = "left") ~ Sex + arm + age_c,
                                data = x, dist = "gaussian"))$table["SexMale", ],
                error = function(e) rep(NA, 4), warning = function(w) rep(NA, 4))
  if (mean(x$y > log10(LLOQ[[as]]) + 1e-9) < 0.2) s <- rep(NA, 4)  # < 20% above the limit: not estimable
  C[[length(C) + 1]] <- gmr(s[1], s[1] - 1.96 * s[2], s[1] + 1.96 * s[2], s[4], assay = as, timepoint = DAYS[[dd]], model = "arm + age, Tobit",
                            arm = "all vaccinated", n_female = sum(x$Sex == "Female"), n_male = sum(x$Sex == "Male"))
  for (a in ARMS) C[[length(C) + 1]] <- lm_sex(y ~ Sex, x[x$arm == a, ], assay = as, timepoint = DAYS[[dd]], model = "per arm", arm = a)
  E[[length(E) + 1]] <- data.frame(assay = as, analysis = DAYS[[dd]], test = "sex x arm interaction (F test, + age)",
                                   pvalue = tryCatch(anova(lm(y ~ Sex + arm + age_c, x), lm(y ~ Sex * arm + age_c, x))$`Pr(>F)`[2],
                                                     error = function(e) NA))
  E[[length(E) + 1]] <- data.frame(assay = as, analysis = DAYS[[dd]], test = "Wilcoxon rank-sum, unadjusted",
                                   pvalue = wilcox.test(y ~ Sex, x, exact = FALSE)$p.value)
}
C <- do.call(rbind, C); rownames(C) <- NULL
C$day <- as.numeric(sub("day (\\d+).*", "\\1", C$timepoint))
C$p_holm_10_days <- ave(C$pvalue, C$assay, C$model, C$arm, FUN = function(p) p.adjust(p, "holm"))
write.csv(C, file.path(TAB, "antibody_sex_difference_by_day.csv"), row.names = FALSE)

# D: longitudinal
D <- list()
for (as in names(LLOQ)) {
  x <- vac[vac$assay == as, ]
  f1 <- lmer(y ~ arm * visit + Sex * visit + age_c + (1 | Subject), x, REML = TRUE)
  a1 <- anova(f1)
  E[[length(E) + 1]] <- data.frame(assay = as, analysis = "longitudinal", test = "sex x time (F test, Satterthwaite)", pvalue = a1["visit:Sex", "Pr(>F)"])
  E[[length(E) + 1]] <- data.frame(assay = as, analysis = "longitudinal", test = "sex, main effect (F test, Satterthwaite)", pvalue = a1["Sex", "Pr(>F)"])
  m1 <- update(f1, REML = FALSE); m0 <- update(m1, . ~ . - visit:Sex)
  E[[length(E) + 1]] <- data.frame(assay = as, analysis = "longitudinal", test = "sex x time (likelihood ratio)", pvalue = anova(m0, m1)$`Pr(>Chisq)`[2])
  m2 <- update(m1, . ~ . + Sex:arm + Sex:arm:visit)
  E[[length(E) + 1]] <- data.frame(assay = as, analysis = "longitudinal", test = "sex x arm and sex x arm x time (likelihood ratio)", pvalue = anova(m1, m2)$`Pr(>Chisq)`[2])
  for (mn in c("all time points", "after D0, D0 titre as covariate")) {
    f <- if (mn == "all time points") f1 else
      lmer(y ~ arm * visit + Sex * visit + age_c + y0 + (1 | Subject), droplevels(x[x$day > 0 & !is.na(x$y0), ]), REML = TRUE)
    em <- summary(contrast(emmeans(f, ~ Sex | visit), list(male_vs_female = c(-1, 1))), infer = TRUE)
    D[[length(D) + 1]] <- data.frame(assay = as, model = mn, visit = em$visit, day = as.numeric(names(VISIT)[match(em$visit, VISIT)]),
      GMR_male_vs_female = 10^em$estimate, CI_low = 10^em$lower.CL, CI_high = 10^em$upper.CL, pvalue = em$p.value,
      p_holm = p.adjust(em$p.value, "holm"))
  }
}
D <- do.call(rbind, D); E <- do.call(rbind, E)
write.csv(D, file.path(TAB, "antibody_sex_difference_mixed_model.csv"), row.names = FALSE)
write.csv(E, file.path(TAB, "antibody_model_tests.csv"), row.names = FALSE)
print(pop); print(C, digits = 3); print(D, digits = 3); print(E, digits = 3)
