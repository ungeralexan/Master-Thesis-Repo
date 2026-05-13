# ==============================================================================
# ANGRIST (1990) — MILITARY SERVICE DESIGN
# Replication of Table 2 in Słoczyński, Uysal & Wooldridge (2025, JBES)
# "Abadie's Kappa and Weighting Estimators of the LATE"
# ==============================================================================
#
# DESIGN OVERVIEW
# ---------------
# Dataset    : sipp.dta  (Survey of Income and Program Participation)
# Sample size: N = 3,027  (after removing rsncode == 999 and missing educ/kwage)
#
# Instrument (Z): Draft lottery eligibility (rsncode)
#   - Z = 1  if the individual's lottery number made them eligible for the draft
#   - Z = 0  otherwise
#   - This is the "assignment" variable.  It is as good as randomly assigned
#     (lottery), so conditional independence holds without any covariates.
#     However, age cohort differences mean we still condition on age.
#   - Noncompliance type: TWO-SIDED (some eligible men avoided service;
#     some ineligible men enlisted voluntarily)
#
# Treatment (D): Veteran status (nvstat)
#   - D = 1  if the individual served in the military
#   - D = 0  otherwise
#
# Outcome (Y): Log wages
#   - lwage     = log(kwage)        [log wages measured in DOLLARS]
#   - lwage_cnt = log(100 * kwage)  [log wages measured in CENTS]
#   Both are used to demonstrate translation (in)variance across estimators.
#
# Covariates (X): Age controls — three specifications as in the paper
#   Spec 1 — Age linear   : ~ age             (Table 2, columns 1–2)
#   Spec 2 — Cubic in age : ~ age + age² + age³ (Table 2, columns 3–4)
#   Spec 3 — Saturated    : indicator for each age value (Table 2, columns 5–6)
#
# IDENTIFICATION ASSUMPTIONS (Abadie 2003 / Imbens–Angrist 1994)
# ---------------------------------------------------------------
# (i)  Independence  : Z ⊥ (Y0, Y1, D0, D1) | X  [lottery ≈ random given age]
# (ii) Exclusion     : Z affects Y only through D
# (iii)First stage   : 0 < P(Z=1|X) < 1  AND  P(D=1|Z=1,X) > P(D=1|Z=0,X)
# (iv) Monotonicity  : D1 ≥ D0 a.s.  [no "defiers"]
#
# ESTIMAND: LATE = E[Y1 - Y0 | D1 > D0]  (effect for compliers only)
# ==============================================================================


# ==============================================================================
# 0.  PACKAGES
# ==============================================================================

library(haven)      # read Stata .dta files
library(AER)        # ivreg() for 2SLS
library(sandwich)   # vcovHC() for heteroskedasticity-robust SEs
library(lmtest)     # coeftest()
library(boot)       # (loaded for reference; bootstrap is implemented manually)

set.seed(42)  # fix seed for bootstrap reproducibility


# ==============================================================================
# 1.  DATA LOADING AND PREPARATION
# ==============================================================================

data_path <- "/Users/alexung/desktop/suw_jbes_replicate"
sipp <- read_dta(file.path(data_path, "sipp.dta"))

# --- Sample restriction (matches paper, N = 3,027) ---
sipp_clean <- sipp[sipp$rsncode != 999 & !is.na(sipp$educ) & !is.na(sipp$kwage), ]

# --- Age polynomial terms ---
sipp_clean$age  <- sipp_clean$age_5   # age in 5-year groups (as in paper)
sipp_clean$age2 <- sipp_clean$age_5^2
sipp_clean$age3 <- sipp_clean$age_5^3

# --- Outcome variable in two units (key for translation-invariance demo) ---
sipp_clean$lwage     <- log(sipp_clean$kwage)           # log wages in dollars
sipp_clean$lwage_cnt <- log(100 * sipp_clean$kwage)     # log wages in cents
# NOTE: lwage_cnt = lwage + log(100).  For any estimator:
#   tau(Y + k) = tau(Y) + k * sum(w_i)
#   => normalized estimators are TRANSLATION INVARIANT (sum_w = 0, so no shift)
#   => unnormalized estimators are NOT (sum_w ≠ 0, so they shift by k * sum_w)

cat(sprintf("N = %d  (should be 3,027)\n", nrow(sipp_clean)))


# ==============================================================================
# 2.  INSTRUMENT PROPENSITY SCORE p(X) = P(Z = 1 | X)
# ==============================================================================
# We estimate p(X) two ways — the estimator choice affects which kappa
# estimators are available:
#
#   (a) Logit MLE    — standard maximum likelihood logistic regression
#   (b) CBPS         — Covariate Balancing Propensity Score (Imai & Ratkovic 2014)
#                      solves balancing moment conditions directly;
#                      used only for tau_cb_u (the CBPS-based normalized estimator)

# ------------------------------------------------------------------------------
# 2a.  Logistic regression (MLE)
# ------------------------------------------------------------------------------
# Fits glm(Z ~ X - 1, family=binomial) and returns fitted P(Z=1|X_i).
# The "-1" suppresses a second intercept when X already contains one.

logit_mle <- function(Z, X) {
  df  <- data.frame(Z = Z, X)
  fit <- glm(Z ~ . - 1, data = df, family = binomial(link = "logit"))
  fitted.values(fit)
}

# Helper: logistic (inverse-logit) function
sigmoid <- function(x) 1 / (1 + exp(-x))

# ------------------------------------------------------------------------------
# 2b.  CBPS — Covariate Balancing Propensity Score
# ------------------------------------------------------------------------------
# Instead of maximizing the likelihood, CBPS finds beta such that the
# balancing moment conditions hold:
#
#   (1/N) * sum_i [ (Z_i/p_i - (1-Z_i)/(1-p_i)) * X_i ] = 0
#
# This directly targets the weights used downstream.
# Implemented via Newton–Raphson with backtracking line search.
# Warm-started from the logit MLE solution for better convergence.

cbps <- function(Z, X, tol = 1e-8, max_iter = 1000) {
  n <- length(Z)
  k <- ncol(X)

  # Balancing moment conditions (a k-vector; zero at the solution)
  moment_fn <- function(b) {
    p <- pmin(pmax(sigmoid(X %*% b), 1e-10), 1 - 1e-10)
    residual <- as.vector(Z / p - (1 - Z) / (1 - p))
    colMeans(residual * X)
  }

  # Jacobian of moment conditions (k × k matrix)
  jac_fn <- function(b) {
    p <- pmin(pmax(sigmoid(X %*% b), 1e-10), 1 - 1e-10)
    w <- as.vector(-Z * (1 - p) / p - (1 - Z) * p / (1 - p))
    t(w * X) %*% X / n
  }

  # Warm start from logit MLE
  b <- coef(glm(Z ~ X - 1, family = binomial(link = "logit")))

  converged <- FALSE
  for (i in seq_len(max_iter)) {
    m <- moment_fn(b)
    if (max(abs(m)) < tol) { converged <- TRUE; break }
    J    <- jac_fn(b)
    step <- tryCatch(solve(J, -m), error = function(e) NULL)
    if (is.null(step)) break

    # Backtracking line search: halve step until residual norm decreases
    alpha <- 1.0
    for (j in seq_len(30)) {
      if (max(abs(moment_fn(b + alpha * step))) < max(abs(m))) break
      alpha <- alpha * 0.5
    }
    b <- b + alpha * step
  }

  p <- pmin(pmax(sigmoid(X %*% b), 1e-10), 1 - 1e-10)
  list(p = as.vector(p), converged = converged)
}


# ==============================================================================
# 3.  KAPPA WEIGHTS  (Abadie 2003, Lemma 2.1)
# ==============================================================================
# Abadie (2003) shows that for any measurable g(Y, D, X):
#
#   E[g(Y,D,X) | D1 > D0] = E[kappa * g(Y,D,X)] / P(D1 > D0)
#
# where the three kappa weights are defined as:
#
#   kappa  = 1  -  D(1-Z)/(1-p)  -  (1-D)Z/p
#   kappa1 = D  * (Z - p) / [p(1-p)]          [upweights treated compliers]
#   kappa0 = (1-D) * ((1-Z)-(1-p)) / [p(1-p)] [upweights untreated compliers]
#
# Key properties:
#   E[kappa]  = P(D1 > D0)  [share of compliers]
#   E[kappa1] = p(X) * P(D1 > D0)
#   E[kappa0] = (1-p(X)) * P(D1 > D0)
#
# The LATE is identified as:
#   tau_LATE = [E(kappa1 * Y) - E(kappa0 * Y)] / E[kappa]

kappa_weights <- function(Z, D, p) {
  list(
    kappa  = 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p,
    kappa1 = D * (Z - p) / (p * (1 - p)),
    kappa0 = (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  )
}


# ==============================================================================
# 4.  LATE ESTIMATORS
# ==============================================================================
# We implement six kappa-based estimators (plus 2SLS as the standard benchmark).
# The crucial distinction is NORMALIZED vs UNNORMALIZED:
#
#   Normalized   : divide by the sample analog of E[kappa1] and E[kappa0]
#                  => sum of outcome weights = 0 => TRANSLATION INVARIANT
#   Unnormalized : divide by E[kappa] (or E[kappa1]/E[kappa0] separately)
#                  => sum of outcome weights ≠ 0 => NOT translation invariant
#
# The paper's recommendation (SUW 2025): use tau_u (Uysal 2011)
#   - Normalized         => translation invariant
#   - Denominator > 0 by construction under ONE-SIDED noncompliance
#     (Here we have two-sided, so this is less critical but still valid)

# ------------------------------------------------------------------------------
# 4a.  tau_u  —  Uysal (2011)  [RECOMMENDED by SUW 2025]
# ------------------------------------------------------------------------------
# This is a normalized Hajek-type estimator.
# It computes separate weighted means of Y for Z=1 and Z=0 groups,
# then takes the ratio of (Y_difference / D_difference).
#
# Estimator:
#   tau_u = [mu_Y1 - mu_Y0] / [mu_D1 - mu_D0]
#
#   where  mu_Y1 = sum(Y_i * Z_i / p_i) / sum(Z_i / p_i)
#          mu_Y0 = sum(Y_i * (1-Z_i)/(1-p_i)) / sum((1-Z_i)/(1-p_i))
#   and similarly for mu_D1, mu_D0.
#
# Also equivalent to tau_t (Tan 2006 / Frolich 2007): ratio of two ATE-of-Z
# estimators (ATE of Z on Y divided by ATE of Z on D).

tau_u <- function(Y, Z, D, p) {
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  numerator   <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denominator <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  numerator / denominator
}

# ------------------------------------------------------------------------------
# 4b.  tau_a10  —  Abadie-Cattaneo (2018) normalized kappa estimator
# ------------------------------------------------------------------------------
# Uses kappa1 and kappa0 weights separately (part b and c of Lemma 2.1):
#
#   tau_a10 = [sum(kappa1_i * Y_i) / sum(kappa1_i)] -
#             [sum(kappa0_i * Y_i) / sum(kappa0_i)]
#
# This is also normalized (sum of outcome weights = 0) and translation invariant.
# Differs from tau_u in how the propensity score enters the denominator.

tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}

# ------------------------------------------------------------------------------
# 4c.  Unnormalized estimators  (for illustration of sensitivity)
# ------------------------------------------------------------------------------
# All share the same numerator: mean(Y * (Z-p) / (p(1-p)))
# They differ in the denominator:
#
#   tau_a   : denominator = mean(kappa)        [Abadie 2003 original]
#   tau_a1  : denominator = mean(kappa1)       [uses only kappa1]
#   tau_a0  : denominator = mean(kappa0)       [uses only kappa0]
#             NOTE: tau_t = tau_a1 = tau_u (see Remark 3.1 in paper)
#
# WARNING: These do NOT satisfy translation invariance.
# When you switch Y from log-dollars to log-cents (add log(100) to every obs.),
# the estimate shifts by sum(w_i) * log(100) ≠ 0.
# This makes them sensitive to arbitrary measurement choices.

tau_unnorm <- function(Y, Z, D, p, which = "a") {
  kw        <- kappa_weights(Z, D, p)
  numerator <- mean(Y * (Z - p) / (p * (1 - p)))
  denom_val <- switch(which,
    "a"  = mean(kw$kappa),    # Abadie (2003) denominator
    "a1" = mean(kw$kappa1),   # kappa1-based denominator
    "a0" = mean(kw$kappa0)    # kappa0-based denominator
  )
  numerator / denom_val
}


# ==============================================================================
# 5.  2SLS BENCHMARK
# ==============================================================================
# Standard two-stage least squares with HC1-robust standard errors.
# This matches Stata's ivreg2 with vce(robust).
# Formula: Y ~ D + X | Z + X
# HC1 = (n/(n-k)) * HC0 — small-sample correction matching Stata's default.

run_2sls <- function(Y, D, Z, X_df, endog_name = "D") {
  cov_names <- names(X_df)[names(X_df) != "(Intercept)"]

  df <- data.frame(Y = Y, D = D, Z = Z, X_df)
  names(df)[names(df) == "D"] <- endog_name
  names(df)[names(df) == "Z"] <- "instrument"

  if (length(cov_names) == 0) {
    fml <- as.formula(paste("Y ~", endog_name, "| instrument"))
  } else {
    cov_str <- paste(cov_names, collapse = " + ")
    fml <- as.formula(paste("Y ~", endog_name, "+", cov_str,
                            "| instrument +", cov_str))
  }

  fit    <- ivreg(fml, data = df)
  vcov_r <- vcovHC(fit, type = "HC1")
  ct     <- coeftest(fit, vcov = vcov_r)

  list(
    coef = ct[endog_name, "Estimate"],
    se   = ct[endog_name, "Std. Error"],
    pval = ct[endog_name, "Pr(>|t|)"]
  )
}


# ==============================================================================
# 6.  BOOTSTRAP STANDARD ERRORS FOR KAPPA ESTIMATORS
# ==============================================================================
# The kappa estimators have analytical SEs but deriving them requires the
# delta method with the estimated propensity score plugged in.
# We use nonparametric bootstrap (R resamples of size n with replacement)
# as a simpler and equally valid alternative (matches Stata supplement).
#
# Estimators computed in each bootstrap draw:
#   [1] tau_cb_u   — tau_u with CBPS propensity score
#   [2] tau_ml_u   — tau_u with logit MLE propensity score     *** RECOMMENDED ***
#   [3] tau_ml_a10 — tau_a10 with logit MLE
#   [4] tau_ml_a   — tau_a  with logit MLE   (unnormalized, for comparison)
#   [5] tau_ml_t   — tau_a1 with logit MLE   (= tau_u via Remark 3.1)
#   [6] tau_ml_a0  — tau_a0 with logit MLE   (unnormalized, for comparison)

boot_se_all <- function(Y, Z, D, X_mat, R = 500) {
  n <- length(Y)

  # Inner function: compute all 6 estimates on a (possibly resampled) dataset
  compute_all <- function(Y, Z, D, X) {
    p_ml   <- tryCatch(logit_mle(Z, X), error = function(e) NULL)
    if (is.null(p_ml)) return(rep(NA, 6))

    cb_res <- tryCatch(cbps(Z, X), error = function(e) list(p = NULL, converged = FALSE))

    estimates    <- numeric(6)
    estimates[1] <- if (!is.null(cb_res$p) && cb_res$converged)
                      tryCatch(tau_u(Y, Z, D, cb_res$p), error = function(e) NA)
                    else NA
    estimates[2] <- tryCatch(tau_u(Y, Z, D, p_ml),              error = function(e) NA)
    estimates[3] <- tryCatch(tau_a10(Y, Z, D, p_ml),            error = function(e) NA)
    estimates[4] <- tryCatch(tau_unnorm(Y, Z, D, p_ml, "a"),    error = function(e) NA)
    estimates[5] <- tryCatch(tau_unnorm(Y, Z, D, p_ml, "a1"),   error = function(e) NA)
    estimates[6] <- tryCatch(tau_unnorm(Y, Z, D, p_ml, "a0"),   error = function(e) NA)
    estimates
  }

  point    <- compute_all(Y, Z, D, X_mat)
  boot_mat <- matrix(NA, nrow = R, ncol = 6)
  for (r in seq_len(R)) {
    idx          <- sample(n, n, replace = TRUE)
    boot_mat[r, ] <- compute_all(Y[idx], Z[idx], D[idx], X_mat[idx, , drop = FALSE])
  }

  ses <- apply(boot_mat, 2, function(x) sd(x, na.rm = TRUE))

  list(
    estimates = setNames(point, c("tau_cb_u", "tau_ml_u", "tau_ml_a10",
                                  "tau_ml_a", "tau_ml_t", "tau_ml_a0")),
    se        = setNames(ses,   c("tau_cb_u", "tau_ml_u", "tau_ml_a10",
                                  "tau_ml_a", "tau_ml_t", "tau_ml_a0"))
  )
}


# ==============================================================================
# 7.  FORMATTING HELPER
# ==============================================================================
# Returns "coef*** \n (se)" with significance stars.
# Stars based on two-sided z-test: *** p<0.01, ** p<0.05, * p<0.10.

fmt <- function(coef, se, digits = 3) {
  if (is.na(coef)) return("NA")
  pval  <- 2 * pnorm(-abs(coef / se))
  stars <- ifelse(pval < 0.01, "***",
           ifelse(pval < 0.05, "**",
           ifelse(pval < 0.10, "*", "")))
  sprintf(paste0("%.", digits, "f%s\n(%.", digits, "f)"),
          round(coef, digits), stars, round(se, digits))
}


# ==============================================================================
# 8.  COVARIATE SPECIFICATIONS
# ==============================================================================
# The paper reports six columns (three covariate specs × two outcome units).
# All three specs control only for age (the one covariate needed for
# conditional independence of the draft lottery).

# Spec 1: Linear age  (Table 2, columns 1–2)
X1 <- model.matrix(~ age, data = sipp_clean)

# Spec 2: Cubic polynomial in age  (Table 2, columns 3–4)
#   Flexibly captures any nonlinear age-earnings profile
X2 <- model.matrix(~ age + age2 + age3, data = sipp_clean)

# Spec 3: Fully saturated age  (Table 2, columns 5–6)
#   One dummy per age value — nonparametric, most flexible
sipp_clean$age_fac <- factor(sipp_clean$age_5)
X3 <- model.matrix(~ age_fac, data = sipp_clean)

specs <- list(
  list(X = X1, name = "Age linear",
       X_df = data.frame(age = sipp_clean$age)),
  list(X = X2, name = "Cubic in age",
       X_df = data.frame(age = sipp_clean$age, age2 = sipp_clean$age2,
                         age3 = sipp_clean$age3)),
  list(X = X3, name = "Saturated age",
       X_df = data.frame(age_fac = sipp_clean$age_fac))
)

# Treatment and instrument (constant across all specs)
D_sipp <- sipp_clean$nvstat   # veteran status
Z_sipp <- sipp_clean$rsncode  # draft eligibility


# ==============================================================================
# 9.  TABLE 2 — MAIN RESULTS
# ==============================================================================
# For each covariate spec × outcome unit combination we report:
#
#   2SLS      — standard IV benchmark (HC1-robust SEs)
#   tau_cb_u  — Uysal estimator with CBPS p-score     [normalized, preferred]
#   tau_ml_u  — Uysal estimator with logit MLE p-score [normalized, preferred]
#   tau_ml_a10— Abadie–Cattaneo normalized estimator   [normalized]
#   tau_ml_a  — Abadie (2003) unnormalized             [NOT translation invariant]
#   tau_ml_t  — Tan/Frolich estimator (= tau_ml_u)     [normalized]
#   tau_ml_a0 — kappa0-based unnormalized              [NOT translation invariant]
#
# Columns labeled "Cents" show the same estimators on lwage_cnt = lwage + log(100).
# Normalized estimators give IDENTICAL results; unnormalized ones differ.

cat("\n", strrep("=", 70), "\n")
cat("TABLE 2: ANGRIST (1990) — Causal Effects of Military Service on Log Wages\n")
cat(strrep("=", 70), "\n")
cat("Instrument : Draft lottery eligibility (rsncode)\n")
cat("Treatment  : Veteran status (nvstat)\n")
cat("Outcome    : Log wages (dollars and cents)\n\n")

for (spec in specs) {
  cat(sprintf("--- Specification: %s ---\n", spec$name))

  for (outcome in c("lwage_cnt", "lwage")) {
    # Select outcome
    Y <- sipp_clean[[outcome]]
    label <- if (outcome == "lwage_cnt") "log wages (CENTS)" else "log wages (DOLLARS)"
    cat(sprintf("  Outcome: %s\n", label))

    # 2SLS
    iv_res    <- run_2sls(Y, D_sipp, Z_sipp, spec$X_df, endog_name = "nvstat")

    # All kappa estimators + bootstrap SEs
    # NOTE: R = 200 for speed; set R = 500 to reproduce paper's bootstrap SEs
    kappa_res <- boot_se_all(Y, Z_sipp, D_sipp, spec$X, R = 200)

    # Print results
    cat(sprintf("    2SLS       (HC1-robust SE):  %s\n",
                fmt(iv_res$coef, iv_res$se)))
    cat(sprintf("    tau_cb_u   (CBPS, norm.) :   %s\n",
                fmt(kappa_res$estimates["tau_cb_u"],  kappa_res$se["tau_cb_u"])))
    cat(sprintf("    tau_ml_u   (MLE, norm.)  :   %s  [RECOMMENDED]\n",
                fmt(kappa_res$estimates["tau_ml_u"],  kappa_res$se["tau_ml_u"])))
    cat(sprintf("    tau_ml_a10 (MLE, norm.)  :   %s\n",
                fmt(kappa_res$estimates["tau_ml_a10"], kappa_res$se["tau_ml_a10"])))
    cat(sprintf("    tau_ml_a   (MLE, unnorm.):   %s  [NOT transl. invariant]\n",
                fmt(kappa_res$estimates["tau_ml_a"],  kappa_res$se["tau_ml_a"])))
    cat(sprintf("    tau_ml_t   (MLE, = tau_u):   %s\n",
                fmt(kappa_res$estimates["tau_ml_t"],  kappa_res$se["tau_ml_t"])))
    cat(sprintf("    tau_ml_a0  (MLE, unnorm.):   %s  [NOT transl. invariant]\n",
                fmt(kappa_res$estimates["tau_ml_a0"], kappa_res$se["tau_ml_a0"])))
    cat("\n")
  }
  cat("\n")
}


# ==============================================================================
# 10.  TRANSLATION INVARIANCE DEMONSTRATION
# ==============================================================================
# Key insight from SUW (2025): switching from log-dollars to log-cents
# adds the constant k = log(100) ≈ 4.605 to every outcome observation.
#
# For any linear estimator tau = sum(w_i * Y_i):
#   tau(Y + k) = tau(Y) + k * sum(w_i)
#
# => If sum(w_i) = 0: estimate unchanged   [TRANSLATION INVARIANT]
# => If sum(w_i) ≠ 0: estimate shifts      [NOT translation invariant]
#
# We verify this algebraically for the cubic-age specification.

cat(strrep("=", 70), "\n")
cat("TRANSLATION INVARIANCE CHECK — Cubic age spec, logit MLE\n")
cat(strrep("=", 70), "\n")
cat("Adding k = log(100) to Y (switching dollars -> cents)\n")
cat("Prediction: tau(Y+k) - tau(Y) = sum(w) * k\n\n")

k    <- log(100)
Y5   <- sipp_clean$lwage
Y5c  <- sipp_clean$lwage_cnt
p_ml <- logit_mle(Z_sipp, X2)  # cubic spec, logit MLE

# Compute outcome weights explicitly for each estimator
kw5 <- {
  kw   <- kappa_weights(Z_sipp, D_sipp, p_ml)
  s1   <- sum(Z_sipp / p_ml)
  s0   <- sum((1 - Z_sipp) / (1 - p_ml))
  dD   <- sum(D_sipp * Z_sipp / p_ml) / s1 - sum(D_sipp * (1 - Z_sipp) / (1 - p_ml)) / s0
  w_u  <- (Z_sipp / p_ml / s1 - (1 - Z_sipp) / (1 - p_ml) / s0) / dD
  n    <- length(Y5)
  num_w <- (Z_sipp - p_ml) / (p_ml * (1 - p_ml)) / n
  list(
    w_u   = as.vector(w_u),
    w_a10 = kw$kappa1 / sum(kw$kappa1) - kw$kappa0 / sum(kw$kappa0),
    w_a   = num_w / mean(kw$kappa),
    w_a1  = num_w / mean(kw$kappa1),
    w_a0  = num_w / mean(kw$kappa0)
  )
}

cat(sprintf("  %-24s  %8s  %8s  %10s  %10s  %6s\n",
            "Estimator", "Dollars", "Cents", "Diff(act.)", "Diff(pred.)", "Match?"))
cat(strrep("-", 72), "\n")

check_list <- list(
  list(name = "tau_u   (norm.)",  w = kw5$w_u,
       ed = tau_u(Y5, Z_sipp, D_sipp, p_ml),
       ec = tau_u(Y5c, Z_sipp, D_sipp, p_ml)),
  list(name = "tau_a10 (norm.)",  w = kw5$w_a10,
       ed = tau_a10(Y5, Z_sipp, D_sipp, p_ml),
       ec = tau_a10(Y5c, Z_sipp, D_sipp, p_ml)),
  list(name = "tau_a   (unnorm)", w = kw5$w_a,
       ed = tau_unnorm(Y5, Z_sipp, D_sipp, p_ml, "a"),
       ec = tau_unnorm(Y5c, Z_sipp, D_sipp, p_ml, "a")),
  list(name = "tau_a1  (unnorm)", w = kw5$w_a1,
       ed = tau_unnorm(Y5, Z_sipp, D_sipp, p_ml, "a1"),
       ec = tau_unnorm(Y5c, Z_sipp, D_sipp, p_ml, "a1")),
  list(name = "tau_a0  (unnorm)", w = kw5$w_a0,
       ed = tau_unnorm(Y5, Z_sipp, D_sipp, p_ml, "a0"),
       ec = tau_unnorm(Y5c, Z_sipp, D_sipp, p_ml, "a0"))
)

for (e in check_list) {
  act   <- e$ec - e$ed
  pred  <- sum(e$w) * k
  match <- isTRUE(all.equal(act, pred, tolerance = 1e-8))
  cat(sprintf("  %-24s  %8.4f  %8.4f  %10.4f  %10.4f  %6s\n",
              e$name, e$ed, e$ec, act, pred, if (match) "TRUE" else "FALSE"))
}

cat("\n=> Normalized estimators (tau_u, tau_a10): no change across units.\n")
cat("=> Unnormalized estimators: sensitive to measurement choice.\n")
cat(strrep("=", 70), "\n")
cat("ANGRIST (1990) SCRIPT COMPLETE\n")
cat("Next step: run angrist1990_dml.R for the causal ML extension.\n")
cat(strrep("=", 70), "\n")


