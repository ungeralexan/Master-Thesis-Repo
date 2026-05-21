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


