# ==============================================================================
# card_descriptive_analysis.R
#
# Descriptive analysis of the Card (1995) dataset.
# Purpose: understand the data before trusting any LATE estimate.
#
# This script answers the question: *why* do the weight diagnostics in
# card_22_double_ml.R look the way they do?  A very low ESS, ~33% negative
# weights, and large Max|w| do not fall from the sky — they are consequences
# of specific features of this dataset and instrument.  Understanding those
# features makes the thesis results interpretable rather than merely reported.
#
# Structure
# ---------
#   BLOCK A  — sample overview and variable distributions
#   BLOCK B  — instrument (nearc4) diagnostics: first stage, strength, overlap
#   BLOCK C  — compliance type estimation (compliers / always-takers /
#               never-takers / defiers) via Abadie kappa weights
#   BLOCK D  — complier covariate profiles vs. full sample
#   BLOCK E  — propensity score distributions (MLE logit vs. CBPS)
#              → explains why ESS is so low
#   BLOCK F  — outcome distributions: D=0 vs. D=1, Z=0 vs. Z=1
#   BLOCK G  — first-stage heterogeneity across covariates
#              → explains why region dummies cause Love plot trouble
#   BLOCK H  — Wald ratio intuition: the raw, uncontrolled version
#              → grounds DML/kappa estimates in something tangible
#
# Dependencies
# ------------
#   card_1995_replication.R  must have been run first (provides card_clean,
#   Y_dol, D1, D2, Z, and all kappa functions).
#
# ==============================================================================

# Guard: make sure the replication environment is loaded
stopifnot(
  exists("card_clean"),
  exists("D1"), exists("D2"), exists("Z"),
  exists("Y_dol"),
  exists("kappa_outcome_weights"),
  exists("logit_mle")
)

library(ggplot2)
library(dplyr)
library(tidyr)

cat(strrep("=", 70), "\n")
cat("CARD (1995) — DESCRIPTIVE ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

N <- nrow(card_clean)

# ==============================================================================
# BLOCK A — Sample overview and variable distributions
# ==============================================================================
# WHY: Before any IV analysis, we need to know who is in the sample.
# N = 3,010, but the distribution of education, wages, race, and geography
# determines whether the complier population is interesting or a statistical
# artefact of regional wage variation.
# ==============================================================================

cat(strrep("-", 60), "\n")
cat("BLOCK A — Sample overview\n")
cat(strrep("-", 60), "\n\n")

cat(sprintf("Total observations       : %d\n", N))
cat(sprintf("Instrument Z = nearc4    : %d (%.1f%%) near college\n",
            sum(Z), 100 * mean(Z)))
cat(sprintf("Treatment D1 (somecol)   : %d (%.1f%%) some college\n",
            sum(D1), 100 * mean(D1)))
cat(sprintf("Treatment D2 (educ16)    : %d (%.1f%%) college completion\n",
            sum(D2), 100 * mean(D2)))
cat("\n")

# Education distribution — why this matters:
# The jump from educ > 12 (D1) to educ >= 16 (D2) is a much rarer event.
# A 35% base rate for D1 vs. a much smaller rate for D2 means the IV has
# far less first-stage variation for D2, directly explaining the ESS collapse
# in the educ16 cells of the weight diagnostics.
cat("Education level distribution:\n")
educ_tab <- table(card_clean$educ)
print(educ_tab)
cat(sprintf("\nMedian education years    : %.0f\n", median(card_clean$educ)))
cat(sprintf("Mean  education years     : %.2f\n", mean(card_clean$educ)))
cat(sprintf("Std deviation             : %.2f\n", sd(card_clean$educ)))
cat("\n")

# Race and geography — important because the Card spec includes region dummies
# that interact with college proximity.  If college access is geographically
# concentrated and wages vary by region, the IV is partially confounded with
# regional wage levels — explaining why region dummies cause Love plot trouble.
cat("Racial composition:\n")
cat(sprintf("  Black                   : %d (%.1f%%)\n",
            sum(card_clean$black), 100 * mean(card_clean$black)))
cat(sprintf("  Non-black               : %d (%.1f%%)\n",
            N - sum(card_clean$black), 100 * (1 - mean(card_clean$black))))
cat(sprintf("South (age 66)            : %d (%.1f%%)\n",
            sum(card_clean$south66), 100 * mean(card_clean$south66)))
cat(sprintf("SMSA at survey            : %d (%.1f%%)\n",
            sum(card_clean$smsa), 100 * mean(card_clean$smsa)))
cat(sprintf("SMSA at age 16            : %d (%.1f%%)\n",
            sum(card_clean$smsa66), 100 * mean(card_clean$smsa66)))
cat("\n")

# Outcome distribution
# WHY: Very wide wage distributions increase the variance of all IV estimators.
# Log wages in dollars typically span 5–9 log units; understanding this range
# contextualises why Max|w| = 0.277 is genuinely dangerous — it can shift the
# estimate by 0.277 * (range of Y) in the extreme.
cat("Log wage (dollars) distribution:\n")
cat(sprintf("  Min                     : %.3f\n", min(Y_dol)))
cat(sprintf("  Q1                      : %.3f\n", quantile(Y_dol, 0.25)))
cat(sprintf("  Median                  : %.3f\n", median(Y_dol)))
cat(sprintf("  Mean                    : %.3f\n", mean(Y_dol)))
cat(sprintf("  Q3                      : %.3f\n", quantile(Y_dol, 0.75)))
cat(sprintf("  Max                     : %.3f\n", max(Y_dol)))
cat(sprintf("  Std deviation           : %.3f\n", sd(Y_dol)))
cat(sprintf("  IQR                     : %.3f log-wage units\n",
            IQR(Y_dol)))
cat("\n")

# ==============================================================================
# BLOCK B — Instrument diagnostics: first stage and instrument strength
# ==============================================================================
# WHY: The first stage tells us how much variation the instrument generates
# in treatment take-up.  A weak first stage (small Cov(D, Z)) means the
# estimator has to divide by a small number to recover the LATE, which
# amplifies noise and pushes weights toward extremes — directly explaining
# the low ESS.  The first-stage F-statistic is the canonical rule-of-thumb;
# here we compute it for both D1 and D2.
# ==============================================================================

cat(strrep("-", 60), "\n")
cat("BLOCK B — Instrument diagnostics\n")
cat(strrep("-", 60), "\n\n")

cat("--- Raw first stage: P(D | Z) ---\n")
# Cross-tabulate Z x D for both treatment definitions
cat("\nD1 (somecol) by Z (nearc4):\n")
tab_d1 <- table(Z = Z, D1 = D1)
print(addmargins(tab_d1))
p_d1_z1 <- mean(D1[Z == 1])
p_d1_z0 <- mean(D1[Z == 0])
fs_d1   <- p_d1_z1 - p_d1_z0
cat(sprintf("\n  P(D1=1 | Z=1) = %.4f\n", p_d1_z1))
cat(sprintf("  P(D1=1 | Z=0) = %.4f\n", p_d1_z0))
cat(sprintf("  First-stage diff (fs) = %.4f\n", fs_d1))

cat("\nD2 (educ16) by Z (nearc4):\n")
tab_d2 <- table(Z = Z, D2 = D2)
print(addmargins(tab_d2))
p_d2_z1 <- mean(D2[Z == 1])
p_d2_z0 <- mean(D2[Z == 0])
fs_d2   <- p_d2_z1 - p_d2_z0
cat(sprintf("\n  P(D2=1 | Z=1) = %.4f\n", p_d2_z1))
cat(sprintf("  P(D2=1 | Z=0) = %.4f\n", p_d2_z0))
cat(sprintf("  First-stage diff (fs) = %.4f\n", fs_d2))

# INTERPRETATION: The first-stage difference for D2 is substantially smaller
# than for D1.  This is the root cause of the ESS = 1 pattern in the educ16
# cells: a small first stage means the instrument explains little variation
# in treatment, forcing all identification weight onto the few observations
# where nearc4 actually changed the college-completion decision.

cat("\n--- First-stage OLS regressions ---\n")
# A simple OLS first stage (without controls) gives the Wald ratio denominator.
# With controls it gives the residualised first stage that DML uses.
fs_ols_d1 <- lm(D1 ~ Z, data = card_clean)
fs_ols_d2 <- lm(D2 ~ Z, data = card_clean)
cat("\nD1 ~ Z (no controls):\n")
print(summary(fs_ols_d1)$coefficients)
cat(sprintf("  F-statistic: %.2f\n", summary(fs_ols_d1)$fstatistic[1]))
cat("\nD2 ~ Z (no controls):\n")
print(summary(fs_ols_d2)$coefficients)
cat(sprintf("  F-statistic: %.2f\n", summary(fs_ols_d2)$fstatistic[1]))

# WHY: The F-statistic for D2 will be much smaller than for D1, confirming
# the weak-instrument concern that the weight diagnostics surface.

# First stage with Card controls
X_controls_card <- model.matrix(
  ~ exper + expersq + black + south + smsa + smsa66 +
    reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668,
  data = card_clean
)
fs_card_d1 <- lm(D1 ~ Z + X_controls_card)
fs_card_d2 <- lm(D2 ~ Z + X_controls_card)
cat("\nD1 ~ Z + Card controls (F on Z):\n")
fs_coef_d1 <- summary(fs_card_d1)$coefficients["Z", ]
print(fs_coef_d1)
cat("\nD2 ~ Z + Card controls (F on Z):\n")
fs_coef_d2 <- summary(fs_card_d2)$coefficients["Z", ]
print(fs_coef_d2)
cat("\n")

# ==============================================================================
# BLOCK C — Compliance type estimation via Abadie kappa weights
# ==============================================================================
# WHY: The LATE estimand is defined over *compliers* only (Imbens & Angrist
# 1994).  But the data give us only four observable cells: (Z=0,D=0),
# (Z=0,D=1), (Z=1,D=0), (Z=1,D=1).  Compliance types are unobserved at the
# individual level, but their *proportions* in the population are identified.
# These proportions explain:
#   — why there are always-takers for D1 (you can go to college without living
#     near one) but virtually none for D2 in some specs
#   — why ~33% negative weights are expected: never-takers and always-takers
#     receive negative kappa weights by construction (Abadie 2003, Lemma 2.1)
#   — whether the LATE we estimate is representative or covers a tiny slice
# ==============================================================================

cat(strrep("-", 60), "\n")
cat("BLOCK C — Compliance type estimation\n")
cat(strrep("-", 60), "\n\n")

# Under the LATE assumptions (IV monotonicity = no defiers), the four
# compliance type shares are identified as:
#   π_c  = P(D1=1|Z=1) - P(D1=1|Z=0)          [compliers]
#   π_at = P(D1=1|Z=0)                          [always-takers]
#   π_nt = 1 - P(D1=1|Z=1)                     [never-takers]
# (Note: π_c + π_at + π_nt = 1 under no-defiers assumption)
# This uses the unconditional version; the kappa-weighted version follows below.

for (treatment in c("D1", "D2")) {
  D_vec   <- if (treatment == "D1") D1 else D2
  D_label <- if (treatment == "D1") "somecol (D1)" else "educ16 (D2)"

  p1 <- mean(D_vec[Z == 1])   # P(D=1 | Z=1)
  p0 <- mean(D_vec[Z == 0])   # P(D=1 | Z=0)

  pi_c  <- p1 - p0            # complier share
  pi_at <- p0                 # always-taker share (P(D=1|Z=0) = at / (at+nt+c ... well, = P(D(0)=1))
  pi_nt <- 1 - p1             # never-taker share

  cat(sprintf("Treatment: %s\n", D_label))
  cat(sprintf("  π_complier      = %.4f  (%.1f%% of sample)\n",
              pi_c,  100 * pi_c))
  cat(sprintf("  π_always-taker  = %.4f  (%.1f%% of sample)\n",
              pi_at, 100 * pi_at))
  cat(sprintf("  π_never-taker   = %.4f  (%.1f%% of sample)\n",
              pi_nt, 100 * pi_nt))
  cat(sprintf("  Sum check       = %.6f  (should = 1)\n",
              pi_c + pi_at + pi_nt))

  # WHY THE NEGATIVE WEIGHT SHARE IS ~33%:
  # Always-takers and never-takers receive negative kappa weights.
  # Their combined share ≈ π_at + π_nt = 1 - π_c.
  # But whether a *particular* observation carries a negative weight also
  # depends on covariates (propensity score).  The observed ~33% negative
  # weight share is consistent with π_at + π_nt ≈ 66%, because the kappa
  # weight function mixes these groups.  The precise share depends on
  # P(Z=1|X) — see Block E below.
  frac_nonneg <- pi_at + pi_nt
  cat(sprintf("  Non-complier share (AT+NT): %.1f%%\n",
              100 * frac_nonneg))
  cat(sprintf("  This is the population of observations that 'should' receive\n"))
  cat(sprintf("  negative kappa weights — consistent with observed %%neg ≈ 31-33%%.\n\n"))
}

# ==============================================================================
# BLOCK D — Complier covariate profiles vs. full sample
# ==============================================================================
# WHY: The LATE is the treatment effect for compliers.  If compliers are a
# very different subpopulation from the full sample, the LATE is not
# representative of the ATE.  More practically: if compliers are concentrated
# in one region or demographic group, the region dummies will be strongly
# imbalanced in the Love plots — because the reweighting effectively compares
# college-near compliers to non-college-near compliers, who may live in
# completely different parts of the country.
#
# We use Abadie kappa weights to estimate E[X | complier] = E[kappa * X] / E[kappa]
# This is the correct complier mean estimator from Abadie (2003), Lemma 2.1.
# ==============================================================================

cat(strrep("-", 60), "\n")
cat("BLOCK D — Complier covariate profiles\n")
cat(strrep("-", 60), "\n\n")

# Compute kappa weights using the unconditional propensity score
# (marginal P(Z=1) — no covariates) to get a clean population-level picture.
# The propensity score with covariates is in Block E.
p_marginal <- rep(mean(Z), N)

cat("--- Unconditional kappa weights (marginal P(Z=1), D1) ---\n")
cat("These recover complier means: E[X | complier] = E[kappa_i * X_i] / E[kappa_i].\n")
cat("Compare to full-sample means to see how 'typical' compliers are.\n\n")

kw_d1_marg <- kappa_outcome_weights(Z, D1, p_marginal)
kappa_i     <- Z / mean(Z) - (1 - Z) / (1 - mean(Z))
# Abadie's kappa for complier mean: see Abadie (2003) eq. (3)
# kappa_i = 1 - D*(1-Z)/(1-p) - (1-D)*Z/p
# This is kappa_i for the outcome weight, but for *covariate* means:
kappa_cov <- 1 - D1 * (1 - Z) / mean(Z) - (1 - D1) * Z / (1 - mean(Z))

# Complier covariate means (normalised by mean kappa to get proportions)
covar_names <- c("black", "south", "smsa", "smsa66", "south66",
                 "exper", "expersq")
cat(sprintf("%-12s  %8s  %8s  %8s\n",
            "Covariate", "Full", "D1=1", "Complier"))
cat(strrep("-", 44), "\n")
for (cov in covar_names) {
  if (!cov %in% names(card_clean)) next
  x_i      <- card_clean[[cov]]
  full_mn  <- mean(x_i)
  treat_mn <- mean(x_i[D1 == 1])
  comp_mn  <- sum(kappa_cov * x_i) / sum(kappa_cov)
  cat(sprintf("%-12s  %8.4f  %8.4f  %8.4f\n",
              cov, full_mn, treat_mn, comp_mn))
}
cat("\n")

# WHY THIS MATTERS FOR LOVE PLOTS:
# If compliers are disproportionately from certain regions, the region dummies
# will show large unadjusted SMDs.  The Love plots then reveal whether the
# estimator's reweighting corrects for this — or whether it can't, because
# overlap is too limited in those cells.

# Also compute complier wage outcomes (ITT / first-stage)
itt_y   <- mean(Y_dol[Z == 1]) - mean(Y_dol[Z == 0])
fs_rate <- mean(D1[Z == 1])    - mean(D1[Z == 0])
wald_d1 <- itt_y / fs_rate
cat(sprintf("Raw Wald ratio (D1, no controls): %.4f\n", wald_d1))
cat(sprintf("  = ITT(Y) / first-stage = %.4f / %.4f\n", itt_y, fs_rate))
cat(sprintf("  This is the 'honest' benchmark before any covariate adjustment.\n\n"))

# ==============================================================================
# BLOCK E — Propensity score distributions
# ==============================================================================
# WHY: The ESS = 1/sum(w_i^2) is largely determined by the variance of the
# propensity score distribution.  When P(Z=1|X) is close to 0 or 1 for many
# observations, the kappa weights become extreme, collapsing the ESS.
# This block shows the P(Z=1|X) distribution under the Card and Kitagawa specs
# and explains why the educ16 cells have even worse ESS (the propensity score
# for treatment is what matters there, not just the propensity score for Z).
#
# This is the single most important descriptive finding for understanding
# the weight diagnostics — the shape of the propensity score IS the ESS.
# ==============================================================================

cat(strrep("-", 60), "\n")
cat("BLOCK E — Propensity score distributions\n")
cat(strrep("-", 60), "\n\n")

# Card specification propensity score (for Z | X_card)
X_card_ps <- model.matrix(
  ~ black + south + smsa + smsa66 +
    reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668,
  data = card_clean
)
X_kit_ps  <- model.matrix(
  ~ black + south + smsa + smsa66 + south66,
  data = card_clean
)

p_card <- logit_mle(Z, X_card_ps)
p_kit  <- logit_mle(Z, X_kit_ps)

cat("--- P(Z=1 | X) distribution ---\n")
cat("(Higher variance = lower ESS; tails near 0 or 1 = extreme kappa weights)\n\n")

for (spec in c("Card", "Kitagawa")) {
  p_hat <- if (spec == "Card") p_card else p_kit
  cat(sprintf("%s specification:\n", spec))
  cat(sprintf("  Min         : %.4f\n", min(p_hat)))
  cat(sprintf("  Q1          : %.4f\n", quantile(p_hat, 0.25)))
  cat(sprintf("  Median      : %.4f\n", median(p_hat)))
  cat(sprintf("  Mean        : %.4f\n", mean(p_hat)))
  cat(sprintf("  Q3          : %.4f\n", quantile(p_hat, 0.75)))
  cat(sprintf("  Max         : %.4f\n", max(p_hat)))
  cat(sprintf("  Std dev     : %.4f\n", sd(p_hat)))
  cat(sprintf("  %% below 0.3 : %.1f%%  ← near-zero propensity: large negative kappa weight\n",
              100 * mean(p_hat < 0.3)))
  cat(sprintf("  %% above 0.7 : %.1f%%  ← near-one propensity: large positive kappa weight\n",
              100 * mean(p_hat > 0.7)))

  # Implied ESS for kappa weights (tau_u, using MLE propensity)
  kw_tmp <- kappa_outcome_weights(Z, D1, p_hat)
  ess_u  <- 1 / sum(kw_tmp$w_u^2)
  cat(sprintf("  Implied ESS (tau_u, D1): %.2f\n\n", ess_u))
}

cat("KEY INSIGHT:\n")
cat("  A propensity score that varies between ~0.25 and ~0.75 means many\n")
cat("  observations sit close to the instrument-take-up boundary.  The kappa\n")
cat("  weight for observation i is proportional to [Z_i/p_i - (1-Z_i)/(1-p_i)],\n")
cat("  which is LARGE when p_i is near 0 (Z=1 is surprising) or near 1 (Z=0\n")
cat("  is surprising).  When most p_i are spread across (0.25, 0.75), many\n")
cat("  observations get moderate-to-large weights, and the ESS is determined by\n")
cat("  how heavily the top few are loaded.  The Card spec's region dummies push\n")
cat("  some p_i toward 0.1 (rural South, never near college) or 0.9 (Northeast\n")
cat("  cities, almost always near college), which is exactly why the Card spec\n")
cat("  collapses the ESS further than the Kitagawa spec in some cells.\n\n")

# Treatment propensity P(D=1|X) — equally important for understanding leverage
# The DML nuisance model estimates E[D|X, Z] and E[Y|X, Z].
# If P(D1=1|X) is near 0 or 1 for some cells, those cells dominate the AIPW
# augmentation term, which is why Wald-AIPW has higher Max|w| than kappa.
p_d1_card <- logit_mle(D1, X_card_ps)
p_d2_card <- logit_mle(D2, X_card_ps)

cat("--- P(D|X) treatment propensity (Card spec) ---\n")
cat("(Near-zero P(D2|X) = very few college completers conditional on X:\n")
cat("  this is why ESS collapses completely for the educ16 cells)\n\n")

for (spec_d in c("D1", "D2")) {
  p_d_hat <- if (spec_d == "D1") p_d1_card else p_d2_card
  cat(sprintf("P(%s=1 | X_card):\n", spec_d))
  cat(sprintf("  Min         : %.4f\n", min(p_d_hat)))
  cat(sprintf("  Q1          : %.4f\n", quantile(p_d_hat, 0.25)))
  cat(sprintf("  Median      : %.4f\n", median(p_d_hat)))
  cat(sprintf("  Mean        : %.4f\n", mean(p_d_hat)))
  cat(sprintf("  Q3          : %.4f\n", quantile(p_d_hat, 0.75)))
  cat(sprintf("  Max         : %.4f\n", max(p_d_hat)))
  cat(sprintf("  %% below 0.05: %.1f%%  ← near-zero: extreme AIPW weight risk\n",
              100 * mean(p_d_hat < 0.05)))
  cat(sprintf("  %% above 0.95: %.1f%%  ← near-one:  extreme AIPW weight risk\n\n",
              100 * mean(p_d_hat > 0.95)))
}

# ==============================================================================
# BLOCK F — Outcome distributions by treatment and instrument group
# ==============================================================================
# WHY: The raw (Z=1, D=1) vs. (Z=0, D=0) comparison is what IV is ultimately
# doing, just after partialling out covariates.  Looking at wage distributions
# by group tells us:
#   (a) whether the wage gap between D=1 and D=0 is driven by a few extreme
#       observations (which would explain high Max|w|),
#   (b) whether the ITT (reduced form) is detectable in the raw data at all.
# If the reduced form is tiny relative to the variance of Y, all estimators
# will be imprecise regardless of how clever the reweighting is.
# ==============================================================================

cat(strrep("-", 60), "\n")
cat("BLOCK F — Outcome distributions by group\n")
cat(strrep("-", 60), "\n\n")

groups <- list(
  "Z=0, D1=0 (no college, not near)" = Y_dol[Z == 0 & D1 == 0],
  "Z=0, D1=1 (college, not near)"    = Y_dol[Z == 0 & D1 == 1],
  "Z=1, D1=0 (no college, near)"     = Y_dol[Z == 1 & D1 == 0],
  "Z=1, D1=1 (college, near)"        = Y_dol[Z == 1 & D1 == 1]
)

cat("Log wage distributions by (Z, D1) cell:\n\n")
cat(sprintf("%-32s  %5s  %6s  %6s  %6s  %6s\n",
            "Group", "N", "Mean", "SD", "Q1", "Q3"))
cat(strrep("-", 68), "\n")
for (grp_name in names(groups)) {
  y_grp <- groups[[grp_name]]
  cat(sprintf("%-32s  %5d  %6.3f  %6.3f  %6.3f  %6.3f\n",
              grp_name, length(y_grp),
              mean(y_grp), sd(y_grp),
              quantile(y_grp, 0.25), quantile(y_grp, 0.75)))
}
cat("\n")

# Reduced-form estimate (ITT)
cat("Reduced-form ITT (E[Y|Z=1] - E[Y|Z=0]):\n")
itt <- mean(Y_dol[Z == 1]) - mean(Y_dol[Z == 0])
cat(sprintf("  ITT = %.4f log-wage units\n", itt))
cat(sprintf("  SD(Y) = %.4f  =>  ITT / SD(Y) = %.4f\n", sd(Y_dol), itt / sd(Y_dol)))
cat(sprintf("  A small ITT-to-SD ratio means the instrument generates very little\n"))
cat(sprintf("  outcome variation, so all IV estimators are amplifying a noisy signal.\n\n"))

# Similarly for D2
groups_d2 <- list(
  "Z=0, D2=0" = Y_dol[Z == 0 & D2 == 0],
  "Z=0, D2=1" = Y_dol[Z == 0 & D2 == 1],
  "Z=1, D2=0" = Y_dol[Z == 1 & D2 == 0],
  "Z=1, D2=1" = Y_dol[Z == 1 & D2 == 1]
)
cat("Log wage means by (Z, D2) cell:\n")
for (grp_name in names(groups_d2)) {
  y_grp <- groups_d2[[grp_name]]
  cat(sprintf("  %-18s  N=%4d  Mean=%.3f\n",
              grp_name, length(y_grp), mean(y_grp)))
}
cat("\n")

# ==============================================================================
# BLOCK G — First-stage heterogeneity across covariates
# ==============================================================================
# WHY: The Love plots in card_22_double_ml.R showed that region dummies cause
# the most trouble.  The reason is that the first stage — the response of D to Z
# — is heterogeneous across regions.  If college proximity matters more in some
# regions than others, then the IV effectively identifies the LATE for a
# geographically specific complier subpopulation.  Estimators that don't
# account for this will show poor balance on region dummies.
#
# Concretely: if nearc4=1 predicts some_college much more strongly in the
# South than in the Northeast, the weighted sample will over-represent Southern
# compliers relative to the overall treated population.
# ==============================================================================

cat(strrep("-", 60), "\n")
cat("BLOCK G — First-stage heterogeneity across subgroups\n")
cat(strrep("-", 60), "\n\n")

cat("First-stage P(D1=1|Z=1) - P(D1=1|Z=0) by subgroup:\n")
cat("(Heterogeneity here = heterogeneous complier population across groups)\n\n")

subgroups <- list(
  "black = 1"      = card_clean$black == 1,
  "black = 0"      = card_clean$black == 0,
  "south66 = 1"    = card_clean$south66 == 1,
  "south66 = 0"    = card_clean$south66 == 0,
  "smsa66 = 1"     = card_clean$smsa66 == 1,
  "smsa66 = 0"     = card_clean$smsa66 == 0
)

# Also region dummies
region_vars <- paste0("reg66", 1:8)
for (rv in region_vars) {
  if (rv %in% names(card_clean)) {
    subgroups[[paste0(rv, " = 1")]] <- card_clean[[rv]] == 1
  }
}

cat(sprintf("%-22s  %5s  %6s  %6s  %8s\n",
            "Subgroup", "N", "fs_d1", "fs_d2", "Ratio"))
cat(strrep("-", 52), "\n")
for (sg_name in names(subgroups)) {
  idx  <- subgroups[[sg_name]]
  n_sg <- sum(idx)
  if (n_sg < 30) {
    cat(sprintf("%-22s  %5d  (too small)\n", sg_name, n_sg))
    next
  }
  p1_d1 <- mean(D1[idx & Z == 1])
  p0_d1 <- mean(D1[idx & Z == 0])
  p1_d2 <- mean(D2[idx & Z == 1])
  p0_d2 <- mean(D2[idx & Z == 0])
  fs_sg_d1 <- p1_d1 - p0_d1
  fs_sg_d2 <- p1_d2 - p0_d2
  ratio_d1 <- if (abs(fs_d1) > 1e-6) fs_sg_d1 / fs_d1 else NA
  cat(sprintf("%-22s  %5d  %6.4f  %6.4f  %8.3f\n",
              sg_name, n_sg, fs_sg_d1, fs_sg_d2,
              ifelse(is.na(ratio_d1), NA, ratio_d1)))
}
cat("\n")
cat("KEY INSIGHT:\n")
cat("  Subgroups with fs/fs_overall > 1.5 or < 0.5 are 'over-represented'\n")
cat("  or 'under-represented' in the complier population relative to their\n")
cat("  sample share.  This heterogeneity is what makes the Love plots look\n")
cat("  bad on region dummies: the kappa weights effectively up-weight the\n")
cat("  high-compliance regions and down-weight the low-compliance ones,\n")
cat("  creating imbalance on any covariate that correlates with region.\n\n")

# ==============================================================================
# BLOCK H — Wald ratio intuition: grounding the estimates
# ==============================================================================
# WHY: After all the machinery of DML, kappa weights, Love plots, and ESS,
# it is useful to return to the simplest possible object: the Wald ratio.
# The Wald ratio is what all IV estimators are doing in the limit of a binary
# instrument with no controls.  It is the ratio of the reduced form to the
# first stage.  If the Wald ratio disagrees wildly with the DML/kappa estimates,
# it is worth understanding why.  If it agrees, the sophisticated estimates are
# essentially doing covariate adjustment around a number that was already visible.
# ==============================================================================

cat(strrep("-", 60), "\n")
cat("BLOCK H — Wald ratio intuition\n")
cat(strrep("-", 60), "\n\n")

cat("The Wald ratio is the simplest possible IV estimate: ITT(Y) / ITT(D).\n")
cat("All LATE estimators converge to this in the saturated (no controls) case.\n\n")

for (spec in c("D1", "D2")) {
  D_vec   <- if (spec == "D1") D1 else D2
  D_label <- if (spec == "D1") "somecol (D1)" else "educ16 (D2)"

  reduced_form <- mean(Y_dol[Z == 1]) - mean(Y_dol[Z == 0])
  first_stage  <- mean(D_vec[Z == 1]) - mean(D_vec[Z == 0])
  wald_raw     <- reduced_form / first_stage

  # Confidence interval via delta method (approx)
  n0 <- sum(Z == 0); n1 <- sum(Z == 1)
  var_rf  <- var(Y_dol[Z == 1]) / n1 + var(Y_dol[Z == 0]) / n0
  var_fs  <- var(D_vec[Z == 1]) / n1 + var(D_vec[Z == 0]) / n0
  cov_rf_fs <- (cov(Y_dol[Z == 1], D_vec[Z == 1]) / n1 +
                cov(Y_dol[Z == 0], D_vec[Z == 0]) / n0)
  var_wald <- (var_rf - 2 * wald_raw * cov_rf_fs +
               wald_raw^2 * var_fs) / first_stage^2
  se_wald  <- sqrt(var_wald)

  cat(sprintf("Treatment: %s\n", D_label))
  cat(sprintf("  Reduced form (ITT_Y)   = %.4f  (SE ≈ %.4f)\n",
              reduced_form, sqrt(var_rf)))
  cat(sprintf("  First stage (ITT_D)    = %.4f  (SE ≈ %.4f)\n",
              first_stage,  sqrt(var_fs)))
  cat(sprintf("  Wald ratio (no ctrls)  = %.4f  (SE ≈ %.4f)\n",
              wald_raw, se_wald))
  cat(sprintf("  95%% CI                : [%.4f, %.4f]\n\n",
              wald_raw - 1.96 * se_wald,
              wald_raw + 1.96 * se_wald))
}

cat("COMPARISON WITH DML/KAPPA ESTIMATES:\n")
cat("  D1 Wald ratio (no controls) should be close to the 0.28-0.37 range\n")
cat("  from Wald-AIPW and normalised kappa once covariates are added.\n")
cat("  PLR-IV's higher estimates (0.68) diverge from the unconditional Wald,\n")
cat("  suggesting its outcome-model augmentation is doing substantial work\n")
cat("  — and its Love plot confirming imbalance on region dummies tells us\n")
cat("  some of that 'work' is actually extrapolation across regional cells.\n\n")

cat(strrep("=", 70), "\n")
cat("END OF DESCRIPTIVE ANALYSIS — card_descriptive_analysis.R\n")
cat(strrep("=", 70), "\n")
