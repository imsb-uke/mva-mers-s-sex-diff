# Sex difference in the signature scores, linear mixed model with the participant as a random effect.
#
# Data: score of every sample (D0 and D1 of each vaccination) from 06_signature_scores.py.
# Model per arm and score (lme4/lmerTest, REML, Satterthwaite degrees of freedom):
#   score ~ 0 + sex:vaccination:day + age + (1 | participant)
# Contrasts:
#   D0 sex difference per vaccination         female D0 - male D0
#   response sex difference per vaccination   (female D1 - D0) - (male D1 - D0)
#   response sex difference, mean over T1 to T3
# With complete data the response contrast equals the female minus male difference in the paired
# change. The random intercept also allows the D0 comparison and pools the three vaccinations.
#
# Input:  results/signature_scores/sample_scores.csv
# Output: tables/signature_sex_difference_mixed_model.csv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "paths.R"))
suppressPackageStartupMessages({ library(lmerTest) })
d <- read.csv(file.path(RES, "signature_scores", "sample_scores.csv"))
d$grp <- factor(paste(d$sex, d$vaccination, d$day, sep = "_"))
d$age_c <- as.numeric(scale(d$age))

out <- list()
for (sc in unique(d$score)) for (arm in c("SI", "LI")) {
  x <- d[d$score == sc & d$arm == arm, ]
  fit <- lmer(value ~ 0 + grp + age_c + (1 | Subject), data = x, REML = TRUE)
  cn <- names(fixef(fit))
  vec <- function(w) { v <- setNames(numeric(length(cn)), cn); v[paste0("grp", names(w))] <- w; v }
  C <- list()
  for (v in c("T1", "T2", "T3")) {
    C[[paste("D0", v)]] <- vec(setNames(c(1, -1), paste0(c("Female", "Male"), "_", v, "_D0")))
    C[[paste("response", v)]] <- vec(setNames(c(1, -1, -1, 1),
      c(paste0("Female_", v, "_D1"), paste0("Female_", v, "_D0"), paste0("Male_", v, "_D1"), paste0("Male_", v, "_D0"))))
  }
  C[["response mean"]] <- (C[["response T1"]] + C[["response T2"]] + C[["response T3"]]) / 3
  for (k in names(C)) {
    ct <- contest(fit, matrix(C[[k]], nrow = 1), joint = FALSE, confint = TRUE)
    out[[paste(sc, arm, k)]] <- data.frame(score = sc, arm = arm, contrast = sub(" .*", "", k),
      vaccination = sub(".* ", "", k), difference_F_minus_M = ct$Estimate, CI_low = ct$lower,
      CI_high = ct$upper, df = ct$df, pvalue = ct$`Pr(>|t|)`,
      participant_variance = as.data.frame(VarCorr(fit))$vcov[1],
      residual_variance = sigma(fit)^2)
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL
write.csv(res, file.path(TAB, "signature_sex_difference_mixed_model.csv"), row.names = FALSE)
print(res[, 1:9], digits = 3)
