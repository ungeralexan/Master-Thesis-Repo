# ==============================================================================
# REPLICATION: CARD (1995) — CAUSAL EFFECTS OF COLLEGE EDUCATION ON LOG WAGES
# SUW Table 3 reproduction
#
# Reference:
#   Słoczyński, Uysal & Wooldridge (2025), "Abadie's Kappa and Weighting
#   Estimators of the Local Average Treatment Effect," Journal of Business &
#   Economic Statistics, 43(1), 164–177.
#   Table 3: Causal effects of college education on log wages.
#
# Data source:
#   Card (1995) NLSYM subsample. N = 3,010 men with valid wage and education
#   information, interviewed in 1976.
#
# Design:
#   Instrument:   Z = nearc4, indicator for growing up near a 4-year college
#   Treatments:   D1 = somecol = 1{educ > 12}   ("some college attendance")
#                 D2 = educ16  = 1{educ >= 16}   ("college completion")
#   Outcomes:     lwage_cnt = log(wage)           (wage in cents before log)
#                 lwage     = log(wage / 100)     (wage in dollars before log)
#   NOTE: lwage_cnt = lwage + log(100), i.e. the two outcomes differ by a
#         constant shift of log(100). This is SUW's translation-invariance test.
#
# Covariate specifications (matching SUW Table 3 notes):
#   Card spec    : exper, expersq, black, south, smsa, smsa66,
#                  region dummies reg661–reg668  (as in Card 1995)
#   Kitagawa spec: black, south, smsa, smsa66, south66  (as in Kitagawa 2015)
#
# Estimators (matching SUW Table 3, Panels A–C):
#   Panel A:  2SLS with HC1-robust standard errors
#   Panel B:  Normalized kappa estimators
#             tau_cb_u   — Uysal (2011) normalized, CBPS propensity score
#             tau_ml_u   — Uysal (2011) normalized, logit MLE propensity score
#             tau_ml_a10 — Abadie-Cattaneo normalized (kappa1/kappa0 separate)
#   Panel C:  Unnormalized kappa estimators
#             tau_ml_a   — Abadie (2003) unnormalized, kappa denominator
#             tau_ml_t   — Tan (2006) / Frölich (2007), kappa1 denominator
#             tau_ml_a0  — unnormalized, kappa0 denominator
#
# All kappa standard errors are from the M-estimation sandwich formula in the
# SUW (2025) online appendix. Identical functions to vietnam_14_05.R are used
# to ensure consistency across applications.
#
# Consistency note:
#   The functions below are IDENTICAL in structure and naming to those used in
#   the Vietnam (Angrist 1990) replication. Variable names (Y, Z, D, p, X) are
#   kept generic to make the parallel clear. Any professor reading both files
#   can verify that the same estimation pipeline is applied.
# ==============================================================================


# ==============================================================================
# 0. PACKAGES
# ==============================================================================
# haven     : read Stata .dta files (Card dataset is distributed as card.dta)
# AER       : ivreg() for the 2SLS benchmark
# sandwich  : HC1 robust variance for 2SLS
# lmtest    : coeftest() for extracting robust 2SLS standard errors

library(haven)
library(AER)
library(sandwich)
library(lmtest)


# ==============================================================================
# 1. DATA LOADING AND PREPARATION
# ==============================================================================
# The replication uses the card.dta file distributed with the SUW replication
# package. Set data_path to the folder containing card.dta.
#
# Key variable definitions (Card 1995 NLSYM):
#   wage    : hourly wage in cents (as distributed in the original file)
#   educ    : completed years of schooling
#   nearc4  : 1 if a 4-year college was in the county of residence in 1966
#   exper   : potential labor market experience (age - educ - 6)
#   expersq : exper^2
#   black   : 1 if Black
#   south   : 1 if lived in South in 1976
#   smsa    : 1 if lived in SMSA in 1976
#   smsa66  : 1 if lived in SMSA in 1966
#   south66 : 1 if lived in South in 1966
#   reg661–reg668 : region indicators for 1966 (9 Census divisions minus one)
# ==============================================================================

data_path <- "/Users/alexung/Desktop/master_thesis/code/"          # <--- USER: set this path
card <- read_dta(file.path(data_path, "card.dta"))

# --- Sample restriction ---
# SUW state N = 3,010. Card (1995) restricts to men interviewed in 1976 with
# valid wage and education. The 'iqscore' condition (IQ score available) is NOT
# imposed because it reduces the sample. The SUW code uses the full 3,010 sample.
card_clean <- card[!is.na(card$wage) & !is.na(card$educ), ]
cat(sprintf("Sample size: N = %d\n", nrow(card_clean)))
# Expected: N = 3010

# --- Outcome variable in two units ---
# lwage_cnt = log(wage)       — wage measured in CENTS before taking logs
# lwage     = log(wage / 100) — wage measured in DOLLARS before taking logs
# Relationship: lwage_cnt = lwage + log(100)
# This is the key for the translation-invariance demonstration in Panel C.

card_clean$lwage_cnt <- log(card_clean$wage)           # cents (original units)
card_clean$lwage     <- log(card_clean$wage / 100)     # dollars

cat(sprintf("log(100) = %.6f\n", log(100)))
cat(sprintf("Mean diff (should equal log(100)): %.6f\n",
            mean(card_clean$lwage_cnt - card_clean$lwage)))

# --- Binary treatment indicators (SUW binarizations) ---
# Treatment 1: "some college attendance" — at least one year beyond high school
#   D1 = somecol = 1{educ > 12}
# Treatment 2: "college completion" — completed at least 4 years of college
#   D2 = educ16 = 1{educ >= 16}

card_clean$somecol <- as.integer(card_clean$educ > 12)
card_clean$educ16  <- as.integer(card_clean$educ >= 16)

cat(sprintf("somecol: mean = %.4f (share with > 12 yrs schooling)\n",
            mean(card_clean$somecol)))
cat(sprintf("educ16:  mean = %.4f (share with >= 16 yrs schooling)\n",
            mean(card_clean$educ16)))

# --- Instrument ---
# Z = nearc4: indicator for presence of a 4-year college in local labor market
# in 1966. This is Card's (1995) proximity-to-college instrument.

card_clean$Z_inst <- as.integer(card_clean$nearc4)
cat(sprintf("nearc4 (Z): mean = %.4f (first stage >  0 required)\n",
            mean(card_clean$Z_inst)))

# --- Experience squared (may not exist in some versions of the dataset) ---
if (!"expersq" %in% names(card_clean)) {
  card_clean$expersq <- card_clean$exper^2
}

# --- Region dummies: reg661–reg668 ---
# Card (1995) uses 9 Census division indicators, with one omitted as reference.
# In the NLSYM dataset these are often named reg661, reg662, ..., reg668 (8 dummies,
# leaving region 9 as reference). Check column names and adapt if needed.
region_vars <- grep("^reg66", names(card_clean), value = TRUE)
cat("Region variables found:", paste(region_vars, collapse = ", "), "\n")


# ==============================================================================
# 2. COVARIATE DESIGN MATRICES
# ==============================================================================
# Two specifications (Card and Kitagawa) correspond to the two blocks of columns
# in SUW Table 3. Each specification produces two design matrices:
#   X_kappa : WITH intercept — used for propensity score estimation (logit, CBPS)
#   X_df    : WITHOUT intercept as a data.frame — passed to run_2sls()
#
# Card specification (columns 1–2 and 5–6 in SUW Table 3):
#   Controls: exper, expersq, black, south, smsa, smsa66, reg661–reg668
#   This matches the covariate set in Card (1995) Table 3.
#
# Kitagawa specification (columns 3–4 and 7–8 in SUW Table 3):
#   Controls: black, south, smsa, smsa66, south66
#   This is the parsimonious specification from Kitagawa (2015).
# ==============================================================================

# --- Card specification: design matrices ---

# For the kappa propensity score (logit / CBPS): intercept + covariates
X_card_kappa <- model.matrix(
  ~ exper + expersq + black + south + smsa + smsa66 +
    reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668,
  data = card_clean
)
# NOTE: model.matrix() automatically includes an intercept column "(Intercept)".
# This is the standard design matrix passed to logit_mle() and cbps() below.

# For 2SLS: same covariates as a data.frame, WITHOUT the intercept column
X_card_df <- as.data.frame(X_card_kappa[, -1, drop = FALSE])
# The intercept column is removed because ivreg() adds its own intercept.

# --- Kitagawa specification: design matrices ---

X_kit_kappa <- model.matrix(
  ~ black + south + smsa + smsa66 + south66,
  data = card_clean
)

X_kit_df <- as.data.frame(X_kit_kappa[, -1, drop = FALSE])

# --- Convenience aliases for the two outcome vectors ---
Y_cnt <- card_clean$lwage_cnt   # log wages, cents before log
Y_dol <- card_clean$lwage       # log wages, dollars before log
Z     <- card_clean$Z_inst      # instrument: nearc4
D1    <- card_clean$somecol     # treatment 1: some college
D2    <- card_clean$educ16      # treatment 2: college completion


# ==============================================================================
# 3. TOOLKIT FUNCTIONS
# ==============================================================================
# These functions are IDENTICAL to those used in the Vietnam replication
# (vietnam_14_05.R). The naming convention and argument order are preserved
# so that the two applications can be compared function-by-function.
#
# Functions sourced here:
#   §3a  safe_logit()              — numerically stable sigmoid
#   §3b  prep_design_for_mest()    — standardise non-intercept columns
#   §3c  logit_mle()               — logit propensity score (glm)
#   §3d  fit_logit_alpha()         — logit for M-estimation (returns alpha + X_used)
#   §3e  fit_cbps_alpha()          — CBPS propensity score (Newton with line search)
#   §3f  get_cbps_p()              — convenience wrapper returning p vector
#   §3g  kappa_weights()           — Abadie (2003) kappa, kappa1, kappa0
#   §3h  tau_u()                   — Uysal normalized estimator
#   §3i  tau_a10()                 — Abadie-Cattaneo normalized estimator
#   §3j  tau_unnorm()              — unnormalized kappa estimators (a, a1, a0)
#   §3k  kappa_outcome_weights()   — per-observation outcome weights omega_i
#   §3l  num_jacobian()            — central-difference numerical Jacobian
#   §3m  matrix_inverse_safe()     — safe matrix inversion with fallbacks
#   §3n  sandwich_se_mest()        — sandwich SE from stacked moment function
#   §3o  alpha_moment_matrix()     — propensity-score moment contribution
#   §3p  kappa_analytic_se_one()   — SE for a single kappa estimator
#   §3q  safe_kappa_se()           — error-safe wrapper for kappa SE
#   §3r  kappa_analytic_se_all()   — point estimates + SEs for all estimators
#   §3s  run_2sls()                — 2SLS with HC1-robust SE
#   §3t  weight_diag()             — ESS, sum(w), % negative, max|w|
#   §3u  check_weight_identity()   — algebraic check sum(w*Y) == tau_hat
#   §3v  fmt()                     — coefficient + stars + SE formatting
# ==============================================================================


# ==============================================================================
# §3a  NUMERICALLY STABLE LOGISTIC FUNCTION
# ==============================================================================
# Clips the linear predictor to [-35, 35] before applying sigma(eta) = 1/(1+e^{-eta})
# to avoid exp() overflow/underflow. Used in all propensity score computations.
#
# Context (Card application):
#   With many region dummies and binary covariates the propensity score can be
#   close to 0 or 1 for some observations, making numerical stability important.

safe_logit <- function(eta) {
  eta <- pmin(pmax(eta, -35), 35)
  1 / (1 + exp(-eta))
}


# ==============================================================================
# §3b  DESIGN MATRIX STANDARDISATION
# ==============================================================================
# Centres and scales non-intercept columns for better numerical conditioning.
# With polynomial or dummy-heavy designs the Jacobian of the M-estimation system
# can be ill-conditioned without standardisation. Because the intercept is kept
# fixed at 1 this is a pure reparameterisation: fitted propensity scores are
# identical before and after standardisation.
#
# Context (Card application):
#   The Card spec includes exper, expersq, and 8 region dummies — a mix of
#   continuous and binary variables that benefits from standardisation.

prep_design_for_mest <- function(X) {
  X <- as.matrix(X)
  is_intercept <- apply(X, 2, function(v) all(abs(v - 1) < 1e-12))
  X_new <- X
  for (j in seq_len(ncol(X))) {
    if (!is_intercept[j]) {
      mu  <- mean(X[, j])
      sdj <- sd(X[, j])
      if (is.finite(sdj) && sdj > 1e-12)
        X_new[, j] <- (X[, j] - mu) / sdj
    }
  }
  X_new
}


# ==============================================================================
# §3c  LOGIT MLE PROPENSITY SCORE  [point estimates only]
# ==============================================================================
# Standard glm() logit. The "-1" suppresses a duplicate intercept when X already
# contains a column of 1s. Returns fitted P(Z = 1 | X_i) for all observations.
#
# Context (Card application):
#   Z = nearc4. The propensity score P(nearc4 = 1 | X_i) is estimated once per
#   covariate specification (Card or Kitagawa). The same p_ml vector is then
#   used by tau_ml_u, tau_ml_a10, and all three unnormalized estimators.

logit_mle <- function(Z, X) {
  df  <- data.frame(Z = Z, X)
  fit <- glm(Z ~ . - 1, data = df, family = binomial(link = "logit"))
  fitted.values(fit)
}


# ==============================================================================
# §3d  LOGIT MLE FOR M-ESTIMATION  [returns alpha + standardised X]
# ==============================================================================
# Used internally by kappa_analytic_se_one(). Standardises X via
# prep_design_for_mest() for numerical stability of the Jacobian computation.
# The returned X_used and alpha are on the standardised scale; fitted p values
# match those from logit_mle() when evaluated at the unstandardised X.

fit_logit_alpha <- function(Z, X) {
  Z     <- as.numeric(Z)
  X     <- prep_design_for_mest(X)
  fit   <- glm.fit(x = X, y = Z, family = binomial(link = "logit"))
  alpha <- as.numeric(coef(fit));  alpha[is.na(alpha)] <- 0
  p     <- as.vector(pmin(pmax(safe_logit(X %*% alpha), 1e-8), 1 - 1e-8))
  list(alpha = alpha, p = p, X_used = X)
}


# ==============================================================================
# §3e  CBPS PROPENSITY SCORE  [Newton with line search]
# ==============================================================================
# Solves the covariate-balancing moment condition
#
#   E[ (Z - p(X)) / {p(X)(1-p(X))} * X ] = 0
#
# via Newton steps with backtracking line search. Initialised at the logit MLE.
# Returns the best iterate (lowest max-moment norm) if tolerance is not reached.
#
# Context (Card application):
#   CBPS is used exclusively for the tau_cb_u estimator (Panel B, row 1 of
#   SUW Table 3). All other kappa estimators use the logit MLE propensity score.
#   The CBPS moment condition enforces that the instrument is balanced across
#   covariate cells, which can yield more stable estimates than pure MLE.

fit_cbps_alpha <- function(Z, X, tol = 1e-9, max_iter = 5000) {
  Z <- as.numeric(Z)
  X <- prep_design_for_mest(X)
  n <- length(Z)
  k <- ncol(X)

  b <- tryCatch(
    as.numeric(glm.fit(x = X, y = Z, family = binomial())$coefficients),
    error = function(e) rep(0, k)
  )
  b[is.na(b)] <- 0

  moment_fn <- function(b) {
    p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
    colMeans(as.vector((Z - p) / (p * (1 - p))) * X)
  }

  jac_fn <- function(b) {
    p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
    w <- as.vector(-Z * (1 - p) / p - (1 - Z) * p / (1 - p))
    crossprod(X, w * X) / n
  }

  best_b    <- b
  best_norm <- max(abs(moment_fn(b)))
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    m <- moment_fn(b);  m_norm <- max(abs(m))
    if (m_norm < best_norm) { best_norm <- m_norm;  best_b <- b }
    if (m_norm < tol)       { converged <- TRUE;    best_b <- b;  break }

    J    <- jac_fn(b)
    step <- tryCatch(qr.solve(J, -m), error = function(e) NULL)

    if (is.null(step) || any(!is.finite(step))) {
      for (ridge in c(1e-10, 1e-8, 1e-6, 1e-4)) {
        step <- tryCatch(solve(J + ridge * diag(k), -m), error = function(e) NULL)
        if (!is.null(step) && all(is.finite(step))) break
      }
    }
    if (is.null(step) || any(!is.finite(step))) break

    alpha_step <- 1
    for (j in seq_len(50)) {
      b_new <- b + alpha_step * step
      if (is.finite(max(abs(moment_fn(b_new)))) &&
          max(abs(moment_fn(b_new))) < m_norm) { b <- b_new;  break }
      alpha_step <- alpha_step * 0.5
    }
  }

  b <- best_b
  p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
  list(alpha = as.numeric(b), p = p, X_used = X,
       converged = converged, max_moment = best_norm)
}


# §3f  Convenience wrapper: returns only the propensity score vector p
cbps <- function(Z, X, tol = 1e-9, max_iter = 5000, verbose = FALSE) {
  fit_cbps_alpha(Z, X, tol = tol, max_iter = max_iter)
}

get_cbps_p <- function(Z, X) {
  out <- cbps(Z, X)
  if (is.list(out) && !is.null(out$p)) return(as.vector(out$p))
  as.vector(out)
}


# ==============================================================================
# §3g  KAPPA WEIGHTS  (Abadie 2003, Lemma 2.1)
# ==============================================================================
# The three kappa weights identify complier moments:
#
#   kappa  = 1 - D(1-Z)/(1-p) - (1-D)Z/p
#   kappa1 = D(Z - p) / [p(1-p)]        — identifies complier treated moments
#   kappa0 = (1-D)((1-Z)-(1-p)) / [p(1-p)] — identifies complier control moments
#
# In population:
#   E[kappa]  = E[kappa1] = E[kappa0] = P(D(1) > D(0)) = complier share
#
# Context (Card application):
#   Z = nearc4 (college proximity), D = somecol or educ16, p = P(Z=1|X).
#   With two-sided noncompliance (both compliers and always-takers can exist in
#   the Card design), all three kappa weights are relevant.

kappa_weights <- function(Z, D, p) {
  list(
    kappa  = 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p,
    kappa1 = D * (Z - p) / (p * (1 - p)),
    kappa0 = (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  )
}


# ==============================================================================
# §3h  tau_u — UYSAL (2011) NORMALIZED ESTIMATOR  [translation invariant]
# ==============================================================================
# Computes IPW-weighted means of Y and D separately for Z=1 and Z=0, then takes
# the ratio of the two differences:
#
#   tau_u = [mu_Y1 - mu_Y0] / [mu_D1 - mu_D0],
#
# where mu_Y1 = sum(Z*Y/p) / sum(Z/p),  mu_Y0 = sum((1-Z)*Y/(1-p)) / sum((1-Z)/(1-p)).
#
# Translation invariance: adding a constant k to Y shifts both mu_Y1 and mu_Y0
# by k, so the difference (and hence tau_u) is unchanged. This is confirmed by
# the cents/dollars comparison in Panel B of SUW Table 3 (columns 1 vs 2 etc.).
#
# Context (Card application):
#   tau_u is the recommended estimator in SUW (2025). For the Card design it
#   appears in row 1 of Panel B (with CBPS) and row 2 of Panel B (with MLE).

tau_u <- function(Y, Z, D, p) {
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  numerator   <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denominator <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  numerator / denominator
}


# ==============================================================================
# §3i  tau_a10 — ABADIE-CATTANEO NORMALIZED ESTIMATOR  [translation invariant]
# ==============================================================================
# Separately normalizes the kappa1 and kappa0 weighted outcome means:
#
#   tau_a10 = sum(kappa1 * Y) / sum(kappa1)   <-- complier E[Y(1)|complier]
#           - sum(kappa0 * Y) / sum(kappa0)   <-- complier E[Y(0)|complier]
#
# Also translation invariant because the constant k cancels in each ratio.
#
# Context (Card application):
#   tau_ml_a10 appears in row 3 of Panel B in SUW Table 3.
#   Unlike tau_u, it identifies the two complier potential outcome means
#   separately — kappa1 targets treated compliers, kappa0 targets untreated.

tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}


# ==============================================================================
# §3j  tau_unnorm — UNNORMALIZED KAPPA ESTIMATORS  [NOT translation invariant]
# ==============================================================================
# Common numerator:  Delta = mean[ Y * (Z - p) / {p(1-p)} ]
# Three denominator choices (SUW notation):
#   "a"  : Gamma = mean(kappa)    — Abadie (2003) original
#   "a1" : Gamma = mean(kappa1)   — Tan (2006) / Frölich (2007); tau_t in SUW
#   "a0" : Gamma = mean(kappa0)
#
# NOT translation invariant: adding k to Y shifts Delta by k * mean((Z-p)/(p(1-p))),
# which does not generally cancel with the denominator. This causes the large
# discrepancy between cents and dollars in Panel C of SUW Table 3.
#
# Context (Card application):
#   tau_ml_a, tau_ml_t, tau_ml_a0 in Panel C of Table 3.
#   The Card specification shows particularly large cents/dollars differences
#   (e.g. -0.319 vs 0.170 for tau_ml_a, columns 1–2), which demonstrates the
#   severity of the translation-invariance failure for the Card spec.
#   The Kitagawa spec is more stable because it has fewer covariates.

tau_unnorm <- function(Y, Z, D, p, which = "a") {
  kw        <- kappa_weights(Z, D, p)
  numerator <- mean(Y * (Z - p) / (p * (1 - p)))
  denom_val <- switch(which,
    "a"  = mean(kw$kappa),    # Abadie (2003) denominator
    "a1" = mean(kw$kappa1),   # kappa1-based; tau_t in SUW notation
    "a0" = mean(kw$kappa0)    # kappa0-based
  )
  numerator / denom_val
}


# ==============================================================================
# §3k  OUTCOME-WEIGHT CONSTRUCTORS  (omega_i representation)
# ==============================================================================
# Every kappa estimator can be written as tau_hat = sum_i omega_i * Y_i.
# The weight vectors below make this explicit. Key properties:
#   sum(w) == 0  <=> translation invariant (holds for w_u and w_a10)
#   ESS = 1 / sum(w^2)  — effective sample size
#   % negative weights  — diagnostic for potential instability
#
# Context (Card application):
#   Outcome weights are computed for the unnormalized check (Section 7 below)
#   and for the weight diagnostics table.

kappa_outcome_weights <- function(Z, D, p) {
  n  <- length(Z)
  kw <- kappa_weights(Z, D, p)

  # tau_u weights
  s1  <- sum(Z / p)
  s0  <- sum((1 - Z) / (1 - p))
  dD  <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  w_u <- (Z / p / s1 - (1 - Z) / (1 - p) / s0) / dD

  # tau_a10 weights
  w_a10 <- kw$kappa1 / sum(kw$kappa1) - kw$kappa0 / sum(kw$kappa0)

  # common numerator weight for unnormalized estimators
  num_w <- (Z - p) / (p * (1 - p)) / n

  list(
    w_u   = as.vector(w_u),
    w_a10 = as.vector(w_a10),
    w_a   = as.vector(num_w / mean(kw$kappa)),
    w_a1  = as.vector(num_w / mean(kw$kappa1)),
    w_a0  = as.vector(num_w / mean(kw$kappa0))
  )
}


# ==============================================================================
# §3l  NUMERICAL JACOBIAN  (central differences)
# ==============================================================================
# Used inside sandwich_se_mest() to avoid hand-coding the Jacobian of every
# stacked moment system. Adaptive step size: h = eps * (|theta_j| + 1).
#
# Context (Card application):
#   With 15 propensity-score parameters (Card spec) the Jacobian matrix is
#   large; numerical differentiation avoids hard-coded expressions.

num_jacobian <- function(f, theta, eps = 1e-6) {
  theta <- as.numeric(theta)
  f0 <- f(theta)
  m  <- length(f0)
  k  <- length(theta)
  J  <- matrix(NA_real_, nrow = m, ncol = k)

  for (j in seq_len(k)) {
    h  <- eps * (abs(theta[j]) + 1)
    tp <- theta;  tp[j] <- tp[j] + h
    tm <- theta;  tm[j] <- tm[j] - h
    J[, j] <- (f(tp) - f(tm)) / (2 * h)
  }
  J
}


# ==============================================================================
# §3m  SAFE MATRIX INVERSE  (fallback chain)
# ==============================================================================
# Tries solve() -> qr.solve() -> ridge-regularised solve() -> SVD pseudo-inverse.
# Returns the first finite inverse found.
#
# Context (Card application):
#   The A matrix in the sandwich formula can be near-singular when the
#   propensity score is nearly constant or when region dummies are collinear.

matrix_inverse_safe <- function(A, tol = 1e-10) {
  inv <- tryCatch(solve(A), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)

  inv <- tryCatch(qr.solve(A), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)

  for (ridge in c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
    inv <- tryCatch(solve(A + ridge * diag(ncol(A))), error = function(e) NULL)
    if (!is.null(inv) && all(is.finite(inv))) return(inv)
  }

  # Moore-Penrose via SVD
  sv    <- svd(A)
  d     <- sv$d
  d_inv <- ifelse(d > tol * max(d), 1 / d, 0)
  sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
}


# ==============================================================================
# §3n  SANDWICH SE FROM STACKED MOMENT FUNCTION  (M-estimation)
# ==============================================================================
# Computes the SE of theta[tau_index] from the sandwich formula
#
#   Avar(sqrt(n) * theta_hat) = A^{-1} V (A^{-1})'
#
# where:
#   A = E[d/dtheta' psi(O_i, theta)] — Jacobian of mean moment (estimated numerically)
#   V = Var(psi(O_i, theta_0))       — variance of moment contributions
#
# The standard error for theta[tau_index] is sqrt(vcov[tau_index, tau_index]).
#
# Context (Card application):
#   All kappa standard errors in SUW Table 3 (Panels B and C) use this formula.
#   The routine is called by kappa_analytic_se_one() for each estimator type.

sandwich_se_mest <- function(moment_matrix_fn, theta_hat, tau_index) {
  theta_hat    <- as.numeric(theta_hat)
  psi_hat      <- moment_matrix_fn(theta_hat)
  n            <- nrow(psi_hat)

  A            <- num_jacobian(function(th) colMeans(moment_matrix_fn(th)), theta_hat)
  psi_centered <- scale(psi_hat, center = TRUE, scale = FALSE)
  V            <- crossprod(psi_centered) / n

  A_inv        <- matrix_inverse_safe(A)
  vcov_theta   <- A_inv %*% V %*% t(A_inv) / n

  se2 <- vcov_theta[tau_index, tau_index]
  if (!is.finite(se2)) return(NA_real_)
  sqrt(abs(se2))
}


# ==============================================================================
# §3o  PROPENSITY-SCORE MOMENT MATRIX
# ==============================================================================
# Returns the n x k matrix of per-observation propensity score moment
# contributions psi_alpha_i, which are stacked with the LATE estimator moments.
#
#   logit MLE: psi_alpha_i = (Z_i - p_i) * X_i   [score of logit log-likelihood]
#   CBPS:      psi_alpha_i = (Z_i - p_i) / {p_i(1-p_i)} * X_i
#
# Context (Card application):
#   alpha has dimension equal to ncol(X): 15 for Card spec, 6 for Kitagawa.
#   These moments account for the uncertainty in the estimated propensity score.

alpha_moment_matrix <- function(Z, p, X_used, method) {
  p <- as.vector(p);  Z <- as.vector(Z)
  if (method == "ml") return(as.vector(Z - p) * X_used)
  if (method == "cb") return(as.vector((Z - p) / (p * (1 - p))) * X_used)
  stop("method must be 'ml' or 'cb'")
}


# ==============================================================================
# §3p  ANALYTICAL M-ESTIMATION SE FOR ONE KAPPA ESTIMATOR
# ==============================================================================
# Stacks the propensity-score moments with the LATE estimator moments and calls
# sandwich_se_mest(). The tau_index is the last component of theta_hat.
#
# Supported estimators:
#   "u"   — tau_u  (Uysal normalized)
#     theta = (alpha, mu1, mu0, m1, m0, tau), dim = k + 5
#   "a10" — tau_a10
#     theta = (alpha, Delta1, Gamma1, Delta0, Gamma0, tau), dim = k + 5
#   "a"   — tau_a  (Abadie unnormalized, kappa denominator)
#   "a1"  — tau_a1 = tau_t (Tan/Frölich)
#   "a0"  — tau_a0
#     theta = (alpha, Delta, Gamma, tau), dim = k + 3
#
# Context (Card application):
#   For the Card specification k = 15 (intercept + 14 covariates).
#   For the Kitagawa specification k = 6 (intercept + 5 covariates).
#   The moment systems are identical in structure to the Vietnam replication.

kappa_analytic_se_one <- function(Y, Z, D, X, estimator, method = "ml") {
  Y <- as.numeric(Y);  Z <- as.numeric(Z)
  D <- as.numeric(D);  X <- as.matrix(X)
  n <- length(Y)

  fit       <- if (method == "ml") fit_logit_alpha(Z, X) else fit_cbps_alpha(Z, X)
  alpha_hat <- fit$alpha
  X_used    <- fit$X_used
  k         <- length(alpha_hat)

  am <- function(p) alpha_moment_matrix(Z, p, X_used, method)

  # ---------------------------------------------------------------------------
  # tau_u: theta = (alpha [k], mu1, mu0, m1, m0, tau) — length k+5
  # Moment system:
  #   psi_alpha : propensity score score (logit or CBPS)
  #   psi_mu1   : Z*(Y - mu1)/p                   — IPW mean of Y for Z=1
  #   psi_mu0   : (1-Z)*(Y - mu0)/(1-p)           — IPW mean of Y for Z=0
  #   psi_m1    : Z*(D - m1)/p                    — IPW mean of D for Z=1
  #   psi_m0    : (1-Z)*(D - m0)/(1-p)            — IPW mean of D for Z=0
  #   psi_tau   : (mu1-mu0)/(m1-m0) - tau         — tau_u definition
  # ---------------------------------------------------------------------------
  if (estimator == "u") {
    p   <- fit$p
    s1  <- sum(Z / p);  s0 <- sum((1 - Z) / (1 - p))
    mu1 <- sum(Z * Y / p) / s1;            mu0 <- sum((1 - Z) * Y / (1 - p)) / s0
    m1  <- sum(Z * D / p) / s1;            m0  <- sum((1 - Z) * D / (1 - p)) / s0
    tau <- (mu1 - mu0) / (m1 - m0)

    theta_hat <- c(alpha_hat, mu1, mu0, m1, m0, tau)

    moment_fn <- function(theta) {
      a   <- theta[seq_len(k)]
      mu1 <- theta[k+1]; mu0 <- theta[k+2]
      m1  <- theta[k+3]; m0  <- theta[k+4]
      tau <- theta[k+5]
      p   <- as.vector(pmin(pmax(safe_logit(X_used %*% a), 1e-8), 1-1e-8))
      cbind(am(p),
            psi_mu1 = Z * (Y - mu1) / p,
            psi_mu0 = (1 - Z) * (Y - mu0) / (1 - p),
            psi_m1  = Z * (D - m1) / p,
            psi_m0  = (1 - Z) * (D - m0) / (1 - p),
            psi_tau = (mu1 - mu0) / (m1 - m0) - tau)
    }
    return(sandwich_se_mest(moment_fn, theta_hat, length(theta_hat)))
  }

  # ---------------------------------------------------------------------------
  # tau_a, tau_a1, tau_a0: theta = (alpha [k], Delta, Gamma, tau) — length k+3
  # Moment system:
  #   psi_alpha  : propensity score score
  #   psi_Delta  : Z*Y/p - (1-Z)*Y/(1-p) - Delta   — weighted outcome mean diff
  #   psi_Gamma  : kappa_i (or kappa1_i or kappa0_i) - Gamma
  #   psi_tau    : Delta/Gamma - tau
  # ---------------------------------------------------------------------------
  if (estimator %in% c("a", "a1", "a0")) {
    p     <- fit$p
    kw    <- kappa_weights(Z, D, p)
    Delta <- mean(Y * (Z - p) / (p * (1 - p)))
    Gamma <- switch(estimator,
                    "a"  = mean(kw$kappa),
                    "a1" = mean(kw$kappa1),
                    "a0" = mean(kw$kappa0))
    tau       <- Delta / Gamma
    theta_hat <- c(alpha_hat, Delta, Gamma, tau)

    moment_fn <- function(theta) {
      a     <- theta[seq_len(k)]
      Delta <- theta[k+1]; Gamma <- theta[k+2]; tau <- theta[k+3]
      p     <- as.vector(pmin(pmax(safe_logit(X_used %*% a), 1e-8), 1-1e-8))
      psi_Delta <- Z * Y / p - (1 - Z) * Y / (1 - p) - Delta
      psi_Gamma <- switch(estimator,
        "a"  = 1 - (1-Z)*D/(1-p) - Z*(1-D)/p - Gamma,
        "a1" = Z*D/p - (1-Z)*D/(1-p) - Gamma,
        "a0" = Z*(D-1)/p - (1-Z)*(D-1)/(1-p) - Gamma)
      cbind(am(p),
            psi_Delta = psi_Delta,
            psi_Gamma = psi_Gamma,
            psi_tau   = Delta / Gamma - tau)
    }
    return(sandwich_se_mest(moment_fn, theta_hat, length(theta_hat)))
  }

  # ---------------------------------------------------------------------------
  # tau_a10: theta = (alpha [k], Delta1, Gamma1, Delta0, Gamma0, tau) — length k+5
  # Moment system:
  #   psi_alpha  : propensity score score
  #   psi_Delta1 : kappa1_i * Y_i - Delta1      — kappa1-weighted outcome mean
  #   psi_Gamma1 : Z*D/p - (1-Z)*D/(1-p) - Gamma1
  #   psi_Delta0 : kappa0_i * Y_i - Delta0      — kappa0-weighted outcome mean
  #   psi_Gamma0 : Z*(D-1)/p - (1-Z)*(D-1)/(1-p) - Gamma0
  #   psi_tau    : Delta1/Gamma1 - Delta0/Gamma0 - tau
  # ---------------------------------------------------------------------------
  if (estimator == "a10") {
    p   <- fit$p
    kw  <- kappa_weights(Z, D, p)
    Delta1 <- mean(kw$kappa1 * Y);  Gamma1 <- mean(kw$kappa1)
    Delta0 <- mean(kw$kappa0 * Y);  Gamma0 <- mean(kw$kappa0)
    tau    <- Delta1/Gamma1 - Delta0/Gamma0
    theta_hat <- c(alpha_hat, Delta1, Gamma1, Delta0, Gamma0, tau)

    moment_fn <- function(theta) {
      a      <- theta[seq_len(k)]
      Delta1 <- theta[k+1]; Gamma1 <- theta[k+2]
      Delta0 <- theta[k+3]; Gamma0 <- theta[k+4]; tau <- theta[k+5]
      p      <- as.vector(pmin(pmax(safe_logit(X_used %*% a), 1e-8), 1-1e-8))
      k1_i   <- D * (Z - p) / (p * (1 - p))
      k0_i   <- (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
      cbind(am(p),
            psi_Delta1 = k1_i * Y - Delta1,
            psi_Gamma1 = Z*D/p - (1-Z)*D/(1-p) - Gamma1,
            psi_Delta0 = k0_i * Y - Delta0,
            psi_Gamma0 = Z*(D-1)/p - (1-Z)*(D-1)/(1-p) - Gamma0,
            psi_tau    = Delta1/Gamma1 - Delta0/Gamma0 - tau)
    }
    return(sandwich_se_mest(moment_fn, theta_hat, length(theta_hat)))
  }

  stop("Unknown estimator: must be one of 'u', 'a10', 'a', 'a1', 'a0'")
}


# ==============================================================================
# §3q  ERROR-SAFE WRAPPER FOR KAPPA SE
# ==============================================================================
# Returns NA_real_ (with a warning) if the sandwich system is singular, rather
# than stopping the script mid-table.

safe_kappa_se <- function(Y, Z, D, X_mat, estimator, method) {
  tryCatch(
    kappa_analytic_se_one(Y, Z, D, X_mat, estimator = estimator, method = method),
    error = function(e) {
      warning(sprintf("SE failed: estimator=%s, method=%s. %s",
                      estimator, method, e$message))
      NA_real_
    }
  )
}


# ==============================================================================
# §3r  FULL KAPPA TABLE: ESTIMATES + ANALYTICAL SEs
# ==============================================================================
# Convenience wrapper that computes all six kappa point estimates and their
# corresponding M-estimation standard errors in one call.
#
# Returns:
#   list(estimates = named numeric(6), se = named numeric(6))
#
# Context (Card application):
#   Called once per (treatment, outcome-unit, spec) combination.
#   For SUW Table 3 this means 2 x 2 x 2 = 8 calls (two treatments, two
#   outcome units, two specs), but the cents/dollars comparison uses the same
#   propensity score — only Y changes.

kappa_analytic_se_all <- function(Y, Z, D, X_mat) {
  Y <- as.numeric(Y); Z <- as.numeric(Z)
  D <- as.numeric(D); X_mat <- as.matrix(X_mat)

  p_ml <- logit_mle(Z, X_mat)
  p_cb <- get_cbps_p(Z, X_mat)

  estimates <- c(
    tau_cb_u   = tau_u(Y, Z, D, p_cb),
    tau_ml_u   = tau_u(Y, Z, D, p_ml),
    tau_ml_a10 = tau_a10(Y, Z, D, p_ml),
    tau_ml_a   = tau_unnorm(Y, Z, D, p_ml, "a"),
    tau_ml_t   = tau_unnorm(Y, Z, D, p_ml, "a1"),   # tau_t = tau_a,1 in SUW
    tau_ml_a0  = tau_unnorm(Y, Z, D, p_ml, "a0")
  )

  se <- c(
    tau_cb_u   = safe_kappa_se(Y, Z, D, X_mat, "u",   "cb"),
    tau_ml_u   = safe_kappa_se(Y, Z, D, X_mat, "u",   "ml"),
    tau_ml_a10 = safe_kappa_se(Y, Z, D, X_mat, "a10", "ml"),
    tau_ml_a   = safe_kappa_se(Y, Z, D, X_mat, "a",   "ml"),
    tau_ml_t   = safe_kappa_se(Y, Z, D, X_mat, "a1",  "ml"),
    tau_ml_a0  = safe_kappa_se(Y, Z, D, X_mat, "a0",  "ml")
  )

  list(estimates = estimates, se = se)
}


# ==============================================================================
# §3s  2SLS BENCHMARK WITH HC1-ROBUST SEs
# ==============================================================================
# Standard ivreg() from the AER package with sandwich HC1 robust SEs.
# HC1 = (n/(n-k)) * HC0 — matches Stata's ivreg2 vce(robust) default.
# Formula: Y ~ D + X_controls | Z + X_controls
#
# Context (Card application):
#   Z = nearc4, D = somecol or educ16, X_controls = Card or Kitagawa covariates.
#   2SLS estimates in SUW Table 3 (Panel A) are the benchmark. Note that 2SLS
#   is translation invariant by construction (the constant from the outcome unit
#   shift cancels in the numerator), so columns 1 and 2 should be identical, etc.

run_2sls <- function(Y, D, Z, X_df, endog_name = "D") {
  cov_names <- names(X_df)[names(X_df) != "(Intercept)"]
  df <- data.frame(Y = Y, D = D, Z = Z, X_df)
  names(df)[names(df) == "D"] <- endog_name
  names(df)[names(df) == "Z"] <- "instrument"

  fml <- if (length(cov_names) == 0) {
    as.formula(paste("Y ~", endog_name, "| instrument"))
  } else {
    cs <- paste(cov_names, collapse = " + ")
    as.formula(paste("Y ~", endog_name, "+", cs, "| instrument +", cs))
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
# §3t  WEIGHT DIAGNOSTICS
# ==============================================================================
# Key diagnostics for a weight vector omega:
#   Sum_w     : sum(omega_i) — should be ~0 for translation-invariant estimators
#   ESS       : effective sample size = 1 / sum(omega_i^2)
#   Pct_neg   : % observations with omega_i < 0
#   Max_abs_w : max |omega_i| — outlier detection

weight_diag <- function(w, name) {
  data.frame(
    Estimator = name,
    Sum_w     = round(sum(w), 8),
    ESS       = round(1 / sum(w^2), 0),
    Pct_neg   = round(mean(w < 0) * 100, 1),
    Max_abs_w = round(max(abs(w)), 6),
    stringsAsFactors = FALSE
  )
}

# Algebraic identity check: sum(w * Y) == tau_hat  (should always hold)
check_weight_identity <- function(w, Y, estimate, tol = 1e-8) {
  isTRUE(all.equal(sum(w * Y), estimate, tolerance = tol))
}


# ==============================================================================
# §3v  OUTPUT FORMATTER
# ==============================================================================
# Formats a coefficient with significance stars and standard error in
# parentheses. Stars: *** p<0.01, ** p<0.05, * p<0.10 (two-sided z-test).

fmt <- function(coef, se, digits = 3) {
  if (is.na(coef)) return("NA")
  if (is.na(se) || !is.finite(se) || se <= 0)
    return(sprintf(paste0("%.", digits, "f\n(NA)"), round(coef, digits)))
  pval  <- 2 * pnorm(-abs(coef / se))
  stars <- ifelse(pval < 0.01, "***", ifelse(pval < 0.05, "**",
           ifelse(pval < 0.10, "*", "")))
  sprintf(paste0("%.", digits, "f%s\n(%.", digits, "f)"),
          round(coef, digits), stars, round(se, digits))
}


# ==============================================================================
# 4. ESTIMATION: ALL EIGHT COLUMNS OF SUW TABLE 3
# ==============================================================================
# Table 3 has 8 columns structured as:
#   Cols 1-4 : treatment = somecol (some college, educ > 12)
#   Cols 5-8 : treatment = educ16  (college completion, educ >= 16)
#   Odd cols : outcome = lwage_cnt (cents)
#   Even cols: outcome = lwage     (dollars)
#   Cols 1,2,5,6 : Card specification
#   Cols 3,4,7,8 : Kitagawa specification
#
# We compute estimates for 4 distinct (treatment, spec) combinations, each with
# two outcome units. This gives 4 x 2 = 8 column groups.
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("CARD (1995) REPLICATION — SUW TABLE 3\n")
cat("==============================================================================\n")

# ------------------------------------------------------------------------------
# Helper: run one complete column (2SLS + all kappa estimators)
# ------------------------------------------------------------------------------
# Arguments:
#   Y         : outcome vector (lwage_cnt or lwage)
#   Z         : instrument vector (nearc4)
#   D         : treatment vector (somecol or educ16)
#   X_kappa   : design matrix WITH intercept (for logit/CBPS)
#   X_df      : design matrix WITHOUT intercept as data.frame (for 2SLS)
#   label     : column label for printing
# Returns:
#   list(tsls, kappa_res) for further use

run_column <- function(Y, Z, D, X_kappa, X_df, label) {
  cat(sprintf("  Computing column: %s\n", label))

  # 2SLS
  tsls      <- run_2sls(Y, D, Z, X_df, endog_name = "D")

  # Kappa estimators (point estimates + M-estimation SEs)
  kappa_res <- kappa_analytic_se_all(Y, Z, D, X_kappa)

  list(tsls = tsls, kappa = kappa_res)
}


# --- Compute propensity scores once per specification ---
# (The same p_ml / p_cb is reused for cents and dollars within a spec —
#  the propensity score does not depend on Y, only on Z and X.)

cat("\nEstimating propensity scores...\n")
p_ml_card <- logit_mle(Z, X_card_kappa)
p_cb_card <- get_cbps_p(Z, X_card_kappa)
p_ml_kit  <- logit_mle(Z, X_kit_kappa)
p_cb_kit  <- get_cbps_p(Z, X_kit_kappa)
cat("  Done.\n")


# ==============================================================================
# 4a. TREATMENT 1: somecol = 1{educ > 12}  (some college attendance)
#     Columns 1–4 of SUW Table 3
# ==============================================================================

cat("\n--- Treatment: somecol (educ > 12) ---\n")

# Col 1: somecol, cents, Card spec
col1 <- run_column(Y_cnt, Z, D1, X_card_kappa, X_card_df, "col1: somecol, cents, Card")

# Col 2: somecol, dollars, Card spec
col2 <- run_column(Y_dol, Z, D1, X_card_kappa, X_card_df, "col2: somecol, dollars, Card")

# Col 3: somecol, cents, Kitagawa spec
col3 <- run_column(Y_cnt, Z, D1, X_kit_kappa, X_kit_df, "col3: somecol, cents, Kitagawa")

# Col 4: somecol, dollars, Kitagawa spec
col4 <- run_column(Y_dol, Z, D1, X_kit_kappa, X_kit_df, "col4: somecol, dollars, Kitagawa")


# ==============================================================================
# 4b. TREATMENT 2: educ16 = 1{educ >= 16}  (college completion)
#     Columns 5–8 of SUW Table 3
# ==============================================================================

cat("\n--- Treatment: educ16 (educ >= 16) ---\n")

# Col 5: educ16, cents, Card spec
col5 <- run_column(Y_cnt, Z, D2, X_card_kappa, X_card_df, "col5: educ16, cents, Card")

# Col 6: educ16, dollars, Card spec
col6 <- run_column(Y_dol, Z, D2, X_card_kappa, X_card_df, "col6: educ16, dollars, Card")

# Col 7: educ16, cents, Kitagawa spec
col7 <- run_column(Y_cnt, Z, D2, X_kit_kappa, X_kit_df, "col7: educ16, cents, Kitagawa")

# Col 8: educ16, dollars, Kitagawa spec
col8 <- run_column(Y_dol, Z, D2, X_kit_kappa, X_kit_df, "col8: educ16, dollars, Kitagawa")


# ==============================================================================
# 5. DISPLAY: REPLICATION OF SUW TABLE 3
# ==============================================================================

cols <- list(col1, col2, col3, col4, col5, col6, col7, col8)

print_table3 <- function(cols) {

  cat("\n")
  cat("==============================================================================\n")
  cat("Table 3. Causal effects of college education on log wages.\n")
  cat("Replication of Słoczyński, Uysal & Wooldridge (2025), Table 3.\n")
  cat("==============================================================================\n")
  cat(sprintf("%-18s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
              "", "(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)"))
  cat(sprintf("%-18s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
              "Treatment:", "somecol", "somecol", "somecol", "somecol",
              "educ16", "educ16", "educ16", "educ16"))
  cat(sprintf("%-18s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
              "Outcome:", "cents", "dollars", "cents", "dollars",
              "cents", "dollars", "cents", "dollars"))
  cat(sprintf("%-18s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
              "Spec:", "Card", "Card", "Kitagawa", "Kitagawa",
              "Card", "Card", "Kitagawa", "Kitagawa"))
  cat(strrep("-", 90), "\n")

  # Panel A: 2SLS
  cat("\nPanel A: 2SLS\n")
  tsls_coefs <- sapply(cols, function(x) x$tsls$coef)
  tsls_ses   <- sapply(cols, function(x) x$tsls$se)
  cat(sprintf("%-18s", "2SLS"))
  for (i in 1:8) cat(sprintf("  %8.3f", tsls_coefs[i]))
  cat("\n")
  cat(sprintf("%-18s", ""))
  for (i in 1:8) cat(sprintf("  %8.3f", tsls_ses[i]))
  cat("  [SE]\n")

  # Panel B: Normalized estimators
  cat("\nPanel B: Normalized kappa estimators\n")

  estimators_b <- c("tau_cb_u", "tau_ml_u", "tau_ml_a10")
  labels_b     <- c("tau_cb_u", "tau_ml_u", "tau_ml_a10")

  for (est in estimators_b) {
    ests <- sapply(cols, function(x) x$kappa$estimates[est])
    ses  <- sapply(cols, function(x) x$kappa$se[est])
    cat(sprintf("%-18s", est))
    for (i in 1:8) cat(sprintf("  %8.3f", ests[i]))
    cat("\n")
    cat(sprintf("%-18s", ""))
    for (i in 1:8) cat(sprintf("  %8.3f", ses[i]))
    cat("  [SE]\n")
  }

  # Panel C: Unnormalized estimators
  cat("\nPanel C: Unnormalized kappa estimators\n")

  estimators_c <- c("tau_ml_a", "tau_ml_t", "tau_ml_a0")

  for (est in estimators_c) {
    ests <- sapply(cols, function(x) x$kappa$estimates[est])
    ses  <- sapply(cols, function(x) x$kappa$se[est])
    cat(sprintf("%-18s", est))
    for (i in 1:8) cat(sprintf("  %8.3f", ests[i]))
    cat("\n")
    cat(sprintf("%-18s", ""))
    for (i in 1:8) cat(sprintf("  %8.3f", ses[i]))
    cat("  [SE]\n")
  }

  cat(strrep("-", 90), "\n")
  cat("Notes:\n")
  cat("  2SLS: HC1-robust SEs (vcovHC type='HC1'), matches Stata ivreg2 vce(robust).\n")
  cat("  Kappa: M-estimation sandwich SEs from SUW (2025) online appendix formulas.\n")
  cat("  Propensity score: logit MLE (tau_ml_*) or CBPS (tau_cb_u).\n")
  cat("  Card spec: exper, expersq, black, south, smsa, smsa66, reg661-reg668.\n")
  cat("  Kitagawa spec: black, south, smsa, smsa66, south66.\n")
  cat("  Significance: *** p<0.01, ** p<0.05, * p<0.10 (two-sided z-test).\n")
}

print_table3(cols)


# ==============================================================================
# 6. NUMERICAL VERIFICATION AGAINST SUW TABLE 3
# ==============================================================================
# The SUW Table 3 values are hard-coded for verification. Deviations > 0.001
# indicate a problem with the data, variable construction, or estimation.

cat("\n\n")
cat("==============================================================================\n")
cat("NUMERICAL VERIFICATION — SUW Table 3 target values\n")
cat("==============================================================================\n")

# Target values from SUW Table 3 (printed in the paper)
# Format: list(panel, estimator, col, target_coef, target_se)
targets <- list(
  # Panel A: 2SLS
  list("A", "2SLS",       1, 0.661, 0.294),
  list("A", "2SLS",       2, 0.661, 0.294),
  list("A", "2SLS",       3, 0.575, 0.308),
  list("A", "2SLS",       4, 0.575, 0.308),
  list("A", "2SLS",       5, 1.392, 0.798),
  list("A", "2SLS",       6, 1.392, 0.798),
  list("A", "2SLS",       7, 0.991, 0.610),
  list("A", "2SLS",       8, 0.991, 0.610),
  # Panel B: tau_cb_u (row 1)
  list("B", "tau_cb_u",   1, 0.376, 0.223),
  list("B", "tau_cb_u",   2, 0.376, 0.223),
  list("B", "tau_cb_u",   3, 0.331, 0.236),
  list("B", "tau_cb_u",   4, 0.331, 0.236),
  list("B", "tau_cb_u",   5, 0.853, 0.549),
  list("B", "tau_cb_u",   6, 0.853, 0.549),
  list("B", "tau_cb_u",   7, 0.588, 0.433),
  list("B", "tau_cb_u",   8, 0.588, 0.433),
  # Panel B: tau_ml_u (row 2)
  list("B", "tau_ml_u",   1, 0.331, 0.202),
  list("B", "tau_ml_u",   2, 0.331, 0.202),
  list("B", "tau_ml_u",   3, 0.356, 0.244),
  list("B", "tau_ml_u",   4, 0.356, 0.244),
  list("B", "tau_ml_u",   5, 0.619, 0.387),
  list("B", "tau_ml_u",   6, 0.619, 0.387),
  list("B", "tau_ml_u",   7, 0.628, 0.448),
  list("B", "tau_ml_u",   8, 0.628, 0.448),
  # Panel B: tau_ml_a10 (row 3)
  list("B", "tau_ml_a10", 1, 0.346, 0.200),
  list("B", "tau_ml_a10", 2, 0.346, 0.200),
  list("B", "tau_ml_a10", 3, 0.293, 0.252),
  list("B", "tau_ml_a10", 4, 0.293, 0.252),
  list("B", "tau_ml_a10", 5, 0.586, 0.356),
  list("B", "tau_ml_a10", 6, 0.586, 0.356),
  list("B", "tau_ml_a10", 7, 0.836, 0.821),
  list("B", "tau_ml_a10", 8, 0.836, 0.821),
  # Panel C: tau_ml_a (row 1)
  list("C", "tau_ml_a",   1, -0.319, 1.182),
  list("C", "tau_ml_a",   2,  0.170, 0.370),
  list("C", "tau_ml_a",   3,  2.248, 0.971),
  list("C", "tau_ml_a",   4,  0.842, 0.362),
  list("C", "tau_ml_a",   5, -0.594, 2.184),
  list("C", "tau_ml_a",   6,  0.315, 0.696),
  list("C", "tau_ml_a",   7,  4.317, 2.485),
  list("C", "tau_ml_a",   8,  1.617, 0.891),
  # Panel C: tau_ml_t (row 2)
  list("C", "tau_ml_t",   1, -0.321, 1.201),
  list("C", "tau_ml_t",   2,  0.171, 0.367),
  list("C", "tau_ml_t",   3,  2.053, 0.813),
  list("C", "tau_ml_t",   4,  0.769, 0.308),
  list("C", "tau_ml_t",   5, -0.601, 2.251),
  list("C", "tau_ml_t",   6,  0.319, 0.687),
  list("C", "tau_ml_t",   7,  3.651, 1.780),
  list("C", "tau_ml_t",   8,  1.367, 0.648),
  # Panel C: tau_ml_a0 (row 3)
  list("C", "tau_ml_a0",  1, -0.290, 1.036),
  list("C", "tau_ml_a0",  2,  0.154, 0.354),
  list("C", "tau_ml_a0",  3,  2.846, 1.592),
  list("C", "tau_ml_a0",  4,  1.066, 0.574),
  list("C", "tau_ml_a0",  5, -0.501, 1.728),
  list("C", "tau_ml_a0",  6,  0.266, 0.639),
  list("C", "tau_ml_a0",  7,  7.241, 7.246),
  list("C", "tau_ml_a0",  8,  2.712, 2.577)
)

# Check each target
n_pass <- 0; n_fail <- 0
for (t in targets) {
  panel <- t[[1]]; est <- t[[2]]; col_idx <- t[[3]]
  tgt_coef <- t[[4]]; tgt_se <- t[[5]]

  col_data <- cols[[col_idx]]
  if (panel == "A") {
    got_coef <- col_data$tsls$coef
    got_se   <- col_data$tsls$se
  } else {
    got_coef <- col_data$kappa$estimates[est]
    got_se   <- col_data$kappa$se[est]
  }

  coef_ok <- abs(got_coef - tgt_coef) < 0.002
  se_ok   <- is.na(got_se) || abs(got_se - tgt_se) < 0.005

  status <- if (coef_ok && se_ok) { n_pass <- n_pass + 1; "PASS" } else { n_fail <- n_fail + 1; "FAIL" }
  cat(sprintf("  [%s]  Panel %s, %-12s col %d:  coef %.3f (tgt %.3f)  SE %.3f (tgt %.3f)\n",
              status, panel, est, col_idx, got_coef, tgt_coef,
              ifelse(is.na(got_se), NA, got_se), tgt_se))
}
cat(sprintf("\nTotal: %d PASS, %d FAIL out of %d checks.\n",
            n_pass, n_fail, length(targets)))


# ==============================================================================
# 7. UNNORMALIZED ESTIMATOR CHECK (translation-invariance demonstration)
# ==============================================================================
# For the normalized estimators (tau_cb_u, tau_ml_u, tau_ml_a10):
#   estimate(cents) == estimate(dollars)   [up to numerical tolerance]
#
# For the unnormalized estimators (tau_ml_a, tau_ml_t, tau_ml_a0):
#   estimate(cents) != estimate(dollars)   [demonstrates the failure]
#
# The difference is: estimate(cents) - estimate(dollars) = sum(omega_i) * log(100)
# For translation-invariant estimators: sum(omega_i) == 0, so the diff is 0.
# For unnormalized: sum(omega_i) != 0, so the diff equals sum(omega_i) * log(100).
#
# This exactly replicates the demonstration in SUW (2025) Section 4.2.
# The unnormalized check is performed for both covariate specifications and both
# treatment variables (4 design cells).

k_shift <- log(100)   # constant by which lwage_cnt exceeds lwage

cat("\n\n")
cat("==============================================================================\n")
cat("TRANSLATION-INVARIANCE CHECK\n")
cat("Demonstrates: normalized estimators are cents/dollars invariant;\n")
cat("unnormalized are NOT (SUW 2025, Section 4.2).\n")
cat("==============================================================================\n")

unnorm_check <- function(D, X_kappa, spec_label) {
  p_ml <- logit_mle(Z, X_kappa)
  p_cb <- get_cbps_p(Z, X_kappa)
  kw_ml <- kappa_outcome_weights(Z, D, p_ml)
  kw_cb <- kappa_outcome_weights(Z, D, p_cb)

  cat(sprintf("\n=== %s ===\n", spec_label))
  cat(sprintf("%-18s  %8s  %8s  %10s  %10s  %10s\n",
              "Estimator", "Cents", "Dollars", "Diff(act)",
              "Diff(pred)", "Match?"))
  cat(strrep("-", 72), "\n")

  # Helper: print one row
  pr <- function(name, w, est_cnt, est_dol) {
    act  <- est_cnt - est_dol
    pred <- sum(w) * k_shift
    ok   <- isTRUE(all.equal(act, pred, tolerance = 1e-6))
    cat(sprintf("%-18s  %8.4f  %8.4f  %10.4f  %10.4f  %10s\n",
                name, est_cnt, est_dol, act, pred,
                if (ok) "TRUE" else "FALSE"))
  }

  # Normalized (sum(w) ~ 0, diff ~ 0)
  pr("tau_cb_u",    kw_cb$w_u,
     tau_u(Y_cnt, Z, D, p_cb),    tau_u(Y_dol, Z, D, p_cb))
  pr("tau_ml_u",    kw_ml$w_u,
     tau_u(Y_cnt, Z, D, p_ml),    tau_u(Y_dol, Z, D, p_ml))
  pr("tau_ml_a10",  kw_ml$w_a10,
     tau_a10(Y_cnt, Z, D, p_ml),  tau_a10(Y_dol, Z, D, p_ml))

  # Unnormalized (sum(w) != 0, diff = sum(w) * log(100) != 0)
  pr("tau_ml_a",    kw_ml$w_a,
     tau_unnorm(Y_cnt, Z, D, p_ml, "a"),
     tau_unnorm(Y_dol, Z, D, p_ml, "a"))
  pr("tau_ml_t",    kw_ml$w_a1,
     tau_unnorm(Y_cnt, Z, D, p_ml, "a1"),
     tau_unnorm(Y_dol, Z, D, p_ml, "a1"))
  pr("tau_ml_a0",   kw_ml$w_a0,
     tau_unnorm(Y_cnt, Z, D, p_ml, "a0"),
     tau_unnorm(Y_dol, Z, D, p_ml, "a0"))
}

unnorm_check(D1, X_card_kappa, "somecol — Card specification")
unnorm_check(D1, X_kit_kappa,  "somecol — Kitagawa specification")
unnorm_check(D2, X_card_kappa, "educ16  — Card specification")
unnorm_check(D2, X_kit_kappa,  "educ16  — Kitagawa specification")


# ==============================================================================
# 8. WEIGHT DIAGNOSTICS TABLE
# ==============================================================================
# Reports sum(omega), ESS, % negative, and max|omega| for all kappa estimators
# across all 4 (treatment, spec) cells. This is the input for the thesis weight
# diagnostics table (analogous to Table 4.4 in thesis_gliederung_updated.pdf).

cat("\n\n")
cat("==============================================================================\n")
cat("WEIGHT DIAGNOSTICS TABLE\n")
cat("==============================================================================\n")

diag_block <- function(D, X_kappa, label) {
  p_ml  <- logit_mle(Z, X_kappa)
  p_cb  <- get_cbps_p(Z, X_kappa)
  kw_ml <- kappa_outcome_weights(Z, D, p_ml)
  kw_cb <- kappa_outcome_weights(Z, D, p_cb)

  cat(sprintf("\n--- %s ---\n", label))
  tab <- rbind(
    weight_diag(kw_cb$w_u,   "tau_cb_u"),
    weight_diag(kw_ml$w_u,   "tau_ml_u"),
    weight_diag(kw_ml$w_a10, "tau_ml_a10"),
    weight_diag(kw_ml$w_a,   "tau_ml_a"),
    weight_diag(kw_ml$w_a1,  "tau_ml_t"),
    weight_diag(kw_ml$w_a0,  "tau_ml_a0")
  )
  print(tab, row.names = FALSE)
}

diag_block(D1, X_card_kappa, "somecol — Card specification")
diag_block(D1, X_kit_kappa,  "somecol — Kitagawa specification")
diag_block(D2, X_card_kappa, "educ16  — Card specification")
diag_block(D2, X_kit_kappa,  "educ16  — Kitagawa specification")


# ==============================================================================
# 9. ALGEBRAIC WEIGHT IDENTITY CHECK
# ==============================================================================
# For every estimator and every (treatment, outcome, spec) combination, verify
# that sum(omega_i * Y_i) == tau_hat. This is a pure algebraic identity and
# must hold to machine precision. If any check fails, the weight construction
# is wrong.

cat("\n\n")
cat("==============================================================================\n")
cat("ALGEBRAIC WEIGHT IDENTITY CHECK: sum(omega * Y) == tau_hat\n")
cat("==============================================================================\n")

identity_check_block <- function(D, X_kappa, Y, label) {
  p_ml  <- logit_mle(Z, X_kappa)
  p_cb  <- get_cbps_p(Z, X_kappa)
  kw_ml <- kappa_outcome_weights(Z, D, p_ml)
  kw_cb <- kappa_outcome_weights(Z, D, p_cb)

  checks <- list(
    list("tau_cb_u",   kw_cb$w_u,   tau_u(Y, Z, D, p_cb)),
    list("tau_ml_u",   kw_ml$w_u,   tau_u(Y, Z, D, p_ml)),
    list("tau_ml_a10", kw_ml$w_a10, tau_a10(Y, Z, D, p_ml)),
    list("tau_ml_a",   kw_ml$w_a,   tau_unnorm(Y, Z, D, p_ml, "a")),
    list("tau_ml_t",   kw_ml$w_a1,  tau_unnorm(Y, Z, D, p_ml, "a1")),
    list("tau_ml_a0",  kw_ml$w_a0,  tau_unnorm(Y, Z, D, p_ml, "a0"))
  )

  cat(sprintf("\n  %s\n", label))
  for (ch in checks) {
    ok <- check_weight_identity(ch[[2]], Y, ch[[3]])
    cat(sprintf("    %-14s : %s  [sum(w*Y)=%.8f, tau=%.8f]\n",
                ch[[1]], if (ok) "PASS" else "FAIL",
                sum(ch[[2]] * Y), ch[[3]]))
  }
}

identity_check_block(D1, X_card_kappa, Y_cnt, "somecol, cents, Card")
identity_check_block(D1, X_card_kappa, Y_dol, "somecol, dollars, Card")
identity_check_block(D1, X_kit_kappa,  Y_cnt, "somecol, cents, Kitagawa")
identity_check_block(D1, X_kit_kappa,  Y_dol, "somecol, dollars, Kitagawa")
identity_check_block(D2, X_card_kappa, Y_cnt, "educ16, cents, Card")
identity_check_block(D2, X_card_kappa, Y_dol, "educ16, dollars, Card")
identity_check_block(D2, X_kit_kappa,  Y_cnt, "educ16, cents, Kitagawa")
identity_check_block(D2, X_kit_kappa,  Y_dol, "educ16, dollars, Kitagawa")


cat("\n\n")
cat("==============================================================================\n")
cat("END OF CARD (1995) REPLICATION\n")
cat("==============================================================================\n")
