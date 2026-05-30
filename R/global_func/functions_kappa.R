# ==============================================================================
# functions_kappa.R
# ==============================================================================
# Shared function library for Abadie-kappa LATE estimators.
# Used identically in all three empirical applications:
#   - Angrist (1990)  — Vietnam draft lottery
#   - Card (1995)     — College proximity
#   - Angrist & Evans (1998) — Childbearing and labor supply
#
# Source this file at the top of every replication script:
#   source("functions_kappa.R")
#
# Variable conventions across all applications:
#   Y   : outcome variable
#   Z   : instrument (binary)
#   D   : treatment (binary)
#   p   : propensity score P(Z = 1 | X_i)
#   X   : covariate design matrix (including intercept column)
#   X_df: covariate data frame (no intercept) for run_2sls()
#
# Section map:
#   §a  safe_logit()                — numerically stable logistic function
#   §b  prep_design_for_mest()      — standardise X for M-estimation
#   §c  logit_mle()                 — logit MLE propensity score (point est.)
#   §d  fit_logit_alpha()           — logit MLE for M-estimation sandwich
#   §e  fit_cbps_alpha()            — CBPS propensity score (Newton, line search)
#   §f  cbps() / get_cbps_p()       — CBPS convenience wrappers
#   §g  kappa_weights()             — kappa, kappa1, kappa0 (Abadie 2003)
#   §h  tau_u()                     — Uysal (2011) normalized estimator
#   §i  tau_a10()                   — Abadie-Cattaneo normalized estimator
#   §j  tau_unnorm()                — unnormalized kappa estimators (a/a1/a0)
#   §k  kappa_outcome_weights()     — omega_i weight vectors for all estimators
#   §l  num_jacobian()              — numerical Jacobian (central differences)
#   §m  matrix_inverse_safe()       — robust matrix inverse (fallback chain)
#   §n  sandwich_se_mest()          — sandwich SE from stacked moment function
#   §o  alpha_moment_matrix()       — propensity-score moment matrix (MLE/CBPS)
#   §p  kappa_analytic_se_one()     — M-estimation SE for one estimator
#   §q  safe_kappa_se()             — error-safe wrapper for kappa SE
#   §r  kappa_analytic_se_all()     — all six estimates + analytical SEs
#   §s  run_2sls()                  — 2SLS benchmark with HC1-robust SEs
#   §t  weight_diag()               — weight diagnostics (ESS, % neg, max|w|)
#   §u  check_weight_identity()     — algebraic check sum(w*Y) == tau_hat
#   §v  fmt()                       — coefficient formatter with stars
#
# CHANGELOG vs earlier per-application versions
# -----------------------------------------------
# Vietnam (vietnam_14_05.R) → functions_kappa.R:
#   • Added: get_cbps_p()          — convenience wrapper (was inline in Card)
#   • Added: kappa_outcome_weights() — weight constructor (new in Card)
#   • Added: alpha_moment_matrix() — propensity moment matrix (new in Card)
#   • Added: weight_diag()         — weight diagnostics (new in Card)
#   • Added: check_weight_identity() — algebraic check (new in Card)
#   • Changed: num_jacobian() — variable names harmonised:
#       Vietnam used th_plus/th_minus; Card used tp/tm.
#       Global version uses tp/tm (Card convention) — no behavioural change.
#   • All other functions are byte-for-byte identical to card_21_05.R.
#
# Card (card_21_05.R) → functions_kappa.R:
#   • No changes. All functions copied verbatim.
#   • get_cbps_p(), weight_diag(), check_weight_identity(),
#     kappa_outcome_weights(), alpha_moment_matrix() are present in Card
#     and now promoted to the global file unchanged.
#
# Angrist & Evans (ae98_replication.R):
#   • This application already sources functions_kappa.R unchanged.
#   • No application-specific function variants exist.
# ==============================================================================


# ==============================================================================
# §a  NUMERICALLY STABLE LOGISTIC FUNCTION
# ==============================================================================
# Clips eta to [-35, 35] to prevent exp() overflow or underflow.
# Used in all propensity score computations.

safe_logit <- function(eta) {
  eta <- pmin(pmax(eta, -35), 35)
  1 / (1 + exp(-eta))
}


# ==============================================================================
# §b  DESIGN MATRIX STANDARDISATION
# ==============================================================================
# Centres and scales non-intercept columns for better numerical conditioning
# of the Jacobian in M-estimation. The intercept column (all 1s) is left
# unchanged. Fitted propensity scores are identical before and after
# standardisation (pure reparameterisation).

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
# §c  LOGIT MLE PROPENSITY SCORE  [point estimates only]
# ==============================================================================
# Standard glm() logit. The "-1" suppresses a duplicate intercept when X
# already contains a column of 1s. Returns fitted P(Z = 1 | X_i).
# Used whenever point estimates only are needed (not M-estimation SEs).

logit_mle <- function(Z, X) {
  df  <- data.frame(Z = Z, X)
  fit <- glm(Z ~ . - 1, data = df, family = binomial(link = "logit"))
  fitted.values(fit)
}


# ==============================================================================
# §d  LOGIT MLE FOR M-ESTIMATION  [returns alpha + standardised X]
# ==============================================================================
# Used internally by kappa_analytic_se_one(). Standardises X via
# prep_design_for_mest() for numerical stability of the Jacobian.
# Returns alpha (on standardised scale), fitted p, and X_used.

fit_logit_alpha <- function(Z, X) {
  Z     <- as.numeric(Z)
  X     <- prep_design_for_mest(X)
  fit   <- glm.fit(x = X, y = Z, family = binomial(link = "logit"))
  alpha <- as.numeric(coef(fit));  alpha[is.na(alpha)] <- 0
  p     <- as.vector(pmin(pmax(safe_logit(X %*% alpha), 1e-8), 1 - 1e-8))
  list(alpha = alpha, p = p, X_used = X)
}


# ==============================================================================
# §e  CBPS PROPENSITY SCORE  [Newton with backtracking line search]
# ==============================================================================
# Solves the covariate-balancing moment condition
#
#   E[ (Z - p(X)) / {p(X)(1-p(X))} * X ] = 0
#
# via Newton steps with backtracking. Initialised at the logit MLE.
# Returns the best iterate (lowest max-moment norm) if tol is not reached.
# Used exclusively for tau_cb_u (CBPS + Uysal normalisation).

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


# ==============================================================================
# §f  CBPS CONVENIENCE WRAPPERS
# ==============================================================================
# cbps()      : calls fit_cbps_alpha(), returns the full list.
# get_cbps_p(): calls cbps() and extracts only the p vector.
#   Used wherever a plain propensity score vector is needed.
#
# CHANGELOG: get_cbps_p() was absent in vietnam_14_05.R (CBPS p was extracted
# inline). Promoted here from card_21_05.R — no behavioural change.

cbps <- function(Z, X, tol = 1e-9, max_iter = 5000, verbose = FALSE) {
  fit_cbps_alpha(Z, X, tol = tol, max_iter = max_iter)
}

get_cbps_p <- function(Z, X) {
  out <- cbps(Z, X)
  if (is.list(out) && !is.null(out$p)) return(as.vector(out$p))
  as.vector(out)
}


# ==============================================================================
# §g  KAPPA WEIGHTS  (Abadie 2003, Lemma 2.1)
# ==============================================================================
# The three kappa weights identify complier moments:
#
#   kappa  = 1 - D(1-Z)/(1-p) - (1-D)Z/p
#   kappa1 = D(Z - p) / [p(1-p)]
#   kappa0 = (1-D)((1-Z)-(1-p)) / [p(1-p)]
#
# In population: E[kappa] = E[kappa1] = E[kappa0] = P(complier).

kappa_weights <- function(Z, D, p) {
  list(
    kappa  = 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p,
    kappa1 = D * (Z - p) / (p * (1 - p)),
    kappa0 = (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  )
}


# ==============================================================================
# §h  tau_u — UYSAL (2011) NORMALIZED ESTIMATOR  [translation invariant]
# ==============================================================================
# Computes separately normalised IPW means for Z=1 and Z=0, then takes
# the ratio of differences:
#
#   tau_u = [mu_Y1 - mu_Y0] / [mu_D1 - mu_D0]
#
# Translation invariant: adding a constant c to Y shifts both mu_Y1 and
# mu_Y0 by c, so the numerator (and hence tau_u) is unchanged.

tau_u <- function(Y, Z, D, p) {
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  numerator   <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denominator <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  numerator / denominator
}


# ==============================================================================
# §i  tau_a10 — ABADIE-CATTANEO NORMALIZED ESTIMATOR  [translation invariant]
# ==============================================================================
# Separately normalises kappa1 and kappa0 weighted outcome means:
#
#   tau_a10 = sum(kappa1 * Y) / sum(kappa1)
#           - sum(kappa0 * Y) / sum(kappa0)
#
# Translation invariant because the constant c cancels in each ratio.

tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}


# ==============================================================================
# §j  tau_unnorm — UNNORMALIZED KAPPA ESTIMATORS  [NOT translation invariant]
# ==============================================================================
# Common numerator: Delta = mean[ Y * (Z - p) / {p(1-p)} ]
# Three denominator choices (SUW notation):
#   "a"  : Gamma = mean(kappa)    — Abadie (2003) original
#   "a1" : Gamma = mean(kappa1)   — Tan (2006) / Frölich (2007); tau_t in SUW
#   "a0" : Gamma = mean(kappa0)
#
# NOT translation invariant: adding c to Y shifts Delta but not Gamma.

tau_unnorm <- function(Y, Z, D, p, which = "a") {
  kw        <- kappa_weights(Z, D, p)
  numerator <- mean(Y * (Z - p) / (p * (1 - p)))
  denom_val <- switch(which,
                      "a"  = mean(kw$kappa),
                      "a1" = mean(kw$kappa1),
                      "a0" = mean(kw$kappa0)
  )
  numerator / denom_val
}


# ==============================================================================
# §k  OUTCOME-WEIGHT CONSTRUCTORS  (omega_i representation)
# ==============================================================================
# Every kappa estimator can be written as tau_hat = sum_i omega_i * Y_i.
# Key properties:
#   sum(w) == 0  ⟺  translation invariant  (holds for w_u and w_a10)
#   ESS = 1 / sum(w^2)
#   % negative weights = diagnostic for instability
#
# CHANGELOG: absent in vietnam_14_05.R; added in card_21_05.R; copied verbatim.

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
# §l  NUMERICAL JACOBIAN  (central differences)
# ==============================================================================
# Used inside sandwich_se_mest() to avoid hand-coding the Jacobian of every
# stacked moment system. Adaptive step size: h = eps * (|theta_j| + 1).
#
# CHANGELOG: vietnam_14_05.R used variable names th_plus/th_minus.
#            card_21_05.R used tp/tm. Global version uses tp/tm (Card
#            convention). Behaviour is identical.

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
# §m  SAFE MATRIX INVERSE  (fallback chain)
# ==============================================================================
# Tries solve() → qr.solve() → ridge-regularised solve() → SVD pseudo-inverse.
# Returns the first finite result found.

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
# §n  SANDWICH SE FROM STACKED MOMENT FUNCTION  (M-estimation)
# ==============================================================================
# Computes the SE of theta[tau_index] from the sandwich formula
#
#   Avar(sqrt(n) * theta_hat) = A^{-1} V (A^{-1})'
#
# where:
#   A = E[d/dtheta' psi(O_i, theta)]  — Jacobian of mean moment (numerical)
#   V = Var(psi(O_i, theta_0))        — variance of moment contributions
#
# The SE for tau_index is sqrt(vcov[tau_index, tau_index]).

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
# §o  PROPENSITY-SCORE MOMENT MATRIX
# ==============================================================================
# Returns the n x k matrix of per-observation propensity-score moment
# contributions psi_alpha_i, stacked with the LATE estimator moments.
#
#   logit MLE: psi_alpha_i = (Z_i - p_i) * X_i   [score of logit log-lik]
#   CBPS:      psi_alpha_i = (Z_i - p_i) / {p_i(1-p_i)} * X_i
#
# CHANGELOG: absent in vietnam_14_05.R (was inline inside kappa_analytic_se_one);
#            extracted into its own function in card_21_05.R; copied verbatim.

alpha_moment_matrix <- function(Z, p, X_used, method) {
  p <- as.vector(p);  Z <- as.vector(Z)
  if (method == "ml") return(as.vector(Z - p) * X_used)
  if (method == "cb") return(as.vector((Z - p) / (p * (1 - p))) * X_used)
  stop("method must be 'ml' or 'cb'")
}


# ==============================================================================
# §p  ANALYTICAL M-ESTIMATION SE FOR ONE KAPPA ESTIMATOR
# ==============================================================================
# Stacks the propensity-score moments with the LATE estimator moments and calls
# sandwich_se_mest(). The tau_index is the last component of theta_hat.
#
# Supported estimators:
#   "u"   — tau_u  (Uysal normalized)
#     theta = (alpha, mu1, mu0, m1, m0, tau), dim = k+5
#   "a10" — tau_a10
#     theta = (alpha, Delta1, Gamma1, Delta0, Gamma0, tau), dim = k+5
#   "a"   — tau_a  (Abadie unnormalized, kappa denominator)
#   "a1"  — tau_a1 = tau_t (Tan/Frölich)
#   "a0"  — tau_a0
#     theta = (alpha, Delta, Gamma, tau), dim = k+3

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
  # tau_a10: theta = (alpha [k], Delta1, Gamma1, Delta0, Gamma0, tau) — k+5
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
# §q  ERROR-SAFE WRAPPER FOR KAPPA SE
# ==============================================================================
# Returns NA_real_ (with a warning) if the sandwich system is singular,
# rather than stopping the script mid-table.

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
# §r  FULL KAPPA TABLE: ESTIMATES + ANALYTICAL SEs
# ==============================================================================
# Convenience wrapper that computes all six kappa point estimates and their
# corresponding M-estimation standard errors in one call.
#
# Returns:
#   list(estimates = named numeric(6), se = named numeric(6))
#
# Estimator names (matching SUW Table notation):
#   tau_cb_u   : CBPS + Uysal normalisation      [Panel B, row 1]
#   tau_ml_u   : MLE  + Uysal normalisation      [Panel B, row 2]
#   tau_ml_a10 : MLE  + Abadie-Cattaneo          [Panel B, row 3]
#   tau_ml_a   : MLE  + unnorm, kappa denom      [Panel C, row 4]
#   tau_ml_t   : MLE  + unnorm, kappa1 denom     [Panel C, row 5]
#   tau_ml_a0  : MLE  + unnorm, kappa0 denom     [Panel C, row 6]

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
# §s  2SLS BENCHMARK WITH HC1-ROBUST SEs
# ==============================================================================
# Standard ivreg() from the AER package with sandwich HC1 robust SEs.
# HC1 = (n/(n-k)) * HC0 — matches Stata's ivreg2 vce(robust) default.
# Formula: Y ~ D + X_controls | Z + X_controls
#
# Arguments:
#   Y          : outcome vector
#   D          : treatment vector
#   Z          : instrument vector
#   X_df       : covariate data.frame (NO intercept column)
#   endog_name : name for D in the formula (default "D")

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
# §t  WEIGHT DIAGNOSTICS
# ==============================================================================
# Key diagnostics for a weight vector omega:
#   Sum_w     : sum(omega_i) — should be ~0 for translation-invariant estimators
#   ESS       : effective sample size = 1 / sum(omega_i^2)
#   Pct_neg   : % observations with omega_i < 0
#   Max_abs_w : max |omega_i| — outlier detection
#
# CHANGELOG: absent in vietnam_14_05.R; added in card_21_05.R; copied verbatim.

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


# ==============================================================================
# §u  ALGEBRAIC IDENTITY CHECK
# ==============================================================================
# Verifies: sum(w * Y) == tau_hat  (should hold to machine precision).
# Returns TRUE/FALSE. Use after computing outcome weights to confirm that
# the weight vector correctly represents the estimator.
#
# CHANGELOG: absent in vietnam_14_05.R; added in card_21_05.R; copied verbatim.

check_weight_identity <- function(w, Y, estimate, tol = 1e-8) {
  isTRUE(all.equal(sum(w * Y), estimate, tolerance = tol))
}


# ==============================================================================
# §v  OUTPUT FORMATTER
# ==============================================================================
# Formats a coefficient with significance stars and standard error in
# parentheses. Stars: *** p<0.01, ** p<0.05, * p<0.10 (two-sided z-test).
# Returns a character string for cat() or table output.

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
