### A bit of data exploration
# ---- 1. What are the variables? ----
cat("Variable overview:\n")
str(sipp_clean[, c("age_5", "nvstat", "rsncode", "kwage", "educ")])

# Quick summary of key variables
cat("\nSummary of key variables:\n")
summary(sipp_clean[, c("age_5", "nvstat", "rsncode", "kwage", "educ")])

# the key variables are 
# rsncode = Z (draft lottery eligibility: 1 = eligible, 0 = not)
# nvstat  = D (veteran status: 1 = veteran)
# kwage   = Y (hourly wage in dollars, before log)


# ---- 2. Treatment-instrument cross-table ----
# This is the most important table — it shows you who compliers, 
# always-takers, never-takers are (in aggregate)
cat("\n2x2 cross-table: Z (rows) x D (columns)\n")
tab <- with(sipp_clean, table(Z = rsncode, D = nvstat))
print(tab)
print(prop.table(tab, margin = 1))  # row proportions



# First-stage: P(D=1|Z=1) vs P(D=1|Z=0)
p_d_z1 <- mean(sipp_clean$nvstat[sipp_clean$rsncode == 1])
p_d_z0 <- mean(sipp_clean$nvstat[sipp_clean$rsncode == 0])
cat(sprintf("\nP(D=1|Z=1) = %.3f  [took up treatment when eligible]\n", p_d_z1))
cat(sprintf("P(D=1|Z=0) = %.3f  [took up treatment when NOT eligible]\n", p_d_z0))
cat(sprintf("First stage (compliance rate) = %.3f\n", p_d_z1 - p_d_z0))




# ---- 3. Who are the always-takers, never-takers, compliers? ----
# Under monotonicity (no defiers):
#   Always-takers:  D=1 regardless → P(D=1|Z=0)
#   Never-takers:   D=0 regardless → P(D=0|Z=1) = 1 - P(D=1|Z=1)
#   Compliers:      switch when encouraged → P(D=1|Z=1) - P(D=1|Z=0)

cat("\nLatent subgroup shares (unconditional, under monotonicity):\n")
cat(sprintf("  Always-takers  (D=1 always):    %.1f%%\n", 100 * p_d_z0))
cat(sprintf("  Never-takers   (D=0 always):    %.1f%%\n", 100 * (1 - p_d_z1)))
cat(sprintf("  Compliers      (respond to Z):  %.1f%%\n", 100 * (p_d_z1 - p_d_z0)))
cat(sprintf("  --> LATE is estimated on ~%.0f people (in expectation)\n",
            (p_d_z1 - p_d_z0) * nrow(sipp_clean)))



# ---- 4. Does compliance vary by age? ----
# This is why controlling for age matters — each cohort had different draft ceilings
comp_by_age <- aggregate(
  cbind(Z = rsncode, D = nvstat) ~ age_5,
  data = sipp_clean,
  FUN = mean
)
comp_by_age$first_stage <- NA
for (a in unique(sipp_clean$age_5)) {
  sub <- sipp_clean[sipp_clean$age_5 == a, ]
  if (nrow(sub) > 10)
    comp_by_age$first_stage[comp_by_age$age_5 == a] <-
      mean(sub$nvstat[sub$rsncode == 1]) - mean(sub$nvstat[sub$rsncode == 0])
}
cat("\nFirst stage by age group:\n")
print(comp_by_age[order(comp_by_age$age_5), ], row.names = FALSE)




# ---- 5. Now look at the propensity score (instrument ps) ----
X5_kappa <- model.matrix(~ age + I(age^2) + I(age^3), data = sipp_clean)
p_ml5    <- logit_mle(sipp_clean$rsncode, X5_kappa)

cat("\nInstrument propensity score p(X) = P(Z=1|X) summary:\n")
print(summary(p_ml5))
hist(p_ml5, breaks = 30, main = "Distribution of p(X)", xlab = "p(X)")
abline(v = mean(sipp_clean$rsncode), col = "red", lty = 2)

# Note: with only age as X, p(X) is a smooth curve over age
# The variation in p(X) reflects the age-specific draft ceilings
plot(sipp_clean$age_5, p_ml5, 
     xlab = "Age", ylab = "p(X) = P(Z=1|X)",
     main = "Instrument propensity score by age",
     pch = 16, col = rgb(0,0,1,0.3))




# ---- 6. Compute kappa weights and look at them ----
kw <- kappa_weights(sipp_clean$rsncode, sipp_clean$nvstat, p_ml5)

cat("\nKappa weight summaries:\n")
cat("kappa (overall):\n");  print(summary(kw$kappa))
cat("kappa1 (treated/complier):\n"); print(summary(kw$kappa1))
cat("kappa0 (control/complier):\n"); print(summary(kw$kappa0))

# Who gets a high kappa weight?
# High kappa1 = treated AND had lottery number below ceiling for their cohort
# High kappa0 = untreated AND had lottery number above ceiling for their cohort
cat(sprintf("\nShare with kappa > 0:  %.1f%%\n", 100 * mean(kw$kappa > 0)))
cat(sprintf("Share with kappa1 > 0: %.1f%%\n", 100 * mean(kw$kappa1 > 0)))
cat(sprintf("Share with kappa0 > 0: %.1f%%\n", 100 * mean(kw$kappa0 > 0)))

# The mean of kappa estimates P(complier) — check against first stage
cat(sprintf("\nmean(kappa) = %.4f  [estimate of P(complier)]\n", mean(kw$kappa)))
cat(sprintf("First stage = %.4f  [also estimates P(complier)]\n", p_d_z1 - p_d_z0))



##### ----- Maximum Likelihood estimation not perfectly aligned 
# Let's derive mean(kappa) analytically from the formula
# mean(kappa) = mean(1) - mean(D*(1-Z)/(1-p)) - mean((1-D)*Z/p)
#             = 1 - E[D(1-Z)/(1-p)] - E[(1-D)Z/p]

# Compute each piece
term1 <- mean(sipp_clean$nvstat * (1 - sipp_clean$rsncode) / (1 - p_ml5))
term2 <- mean((1 - sipp_clean$nvstat) * sipp_clean$rsncode / p_ml5)

cat("1 - term1 - term2 =", 1 - term1 - term2, "\n")
cat("= mean(kappa)     =", mean(kw$kappa), "\n")

# Now compare to first stage
cat("\nFirst stage (raw):", p_d_z1 - p_d_z0, "\n")

# The difference comes from: first stage is unweighted by p(X),
# mean(kappa) uses estimated p(X). If p(X) is badly calibrated,
# they diverge. Check calibration:
cat("\nCalibration check:\n")
cat("mean(Z/p):       ", mean(sipp_clean$rsncode / p_ml5), "\n")
cat("mean((1-Z)/(1-p)):", mean((1-sipp_clean$rsncode)/(1-p_ml5)), "\n")
# Both should equal 1 if ps is well calibrated (logit with intercept guarantees this
# only marginally, not conditionally)



### Ki for the 4 populations exampls
# Show that Z=1,D=1 and Z=0,D=0 always get kappa=1
cat("kappa for Z=1,D=1:", 
    summary(kw$kappa[sipp_clean$rsncode==1 & sipp_clean$nvstat==1]), "\n")
cat("kappa for Z=0,D=0:", 
    summary(kw$kappa[sipp_clean$rsncode==0 & sipp_clean$nvstat==0]), "\n")
cat("kappa for Z=1,D=0:", 
    summary(kw$kappa[sipp_clean$rsncode==1 & sipp_clean$nvstat==0]), "\n")
cat("kappa for Z=0,D=1:", 
    summary(kw$kappa[sipp_clean$rsncode==0 & sipp_clean$nvstat==1]), "\n")

cat("Distribution of p(X):\n")
print(quantile(p_ml5, c(0.01, 0.05, 0.10, 0.90, 0.95, 0.99)))

# Extreme low p(X): young cohorts, rarely eligible
# These people get huge weights if they are never-takers (Z=1,D=0)
low_p  <- sipp_clean[p_ml5 < 0.25, c("age_5","rsncode","nvstat")]
high_p <- sipp_clean[p_ml5 > 0.75, c("age_5","rsncode","nvstat")]

cat("\nAge distribution for low p(X) < 0.25:\n")
print(table(low_p$age_5))
cat("\nAge distribution for high p(X) > 0.75:\n")  
print(table(high_p$age_5))