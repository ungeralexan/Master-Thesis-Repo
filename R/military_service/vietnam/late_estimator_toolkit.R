# ==============================================================================
# LATE ESTIMATOR TOOLKIT
# Extracted from: vietnam_14_05.html, vietmam_14_05_double_ml.html,
#                 vietnam_14_05_double_ml_extended.html
#
# Purpose:
#   Reusable building blocks for LATE estimation under the kappa-weighting
#   framework (Abadie 2003; Słoczyński, Uysal & Wooldridge 2025) and the
#   DML/DoubleML extension (Knaus 2024 / OutcomeWeights package).
#
#   Everything here is design-agnostic: variable names (Y, D, Z, X, p) are
#   generic. When porting to a new design (e.g. Imbens-Angrist 401k), just
#   supply the correct vectors and matrices.
#
# Contents:
#   §1  Numerical / matrix utilities
#   §2  Propensity-score estimators  (logit MLE + CBPS)
#   §3  Kappa weights
#   §4  LATE point estimators        (tau_u, tau_a10, tau_unnorm)
#   §5  Outcome-weight constructors  (per-observation weights)
#   §6  Sandwich / M-estimation SEs
#   §7  2SLS benchmark
#   §8  Weight diagnostics
#   §9  DML helpers                  (OutcomeWeights / DoubleML wrappers)
#   §10 Output helpers               (fmt, print_inv_table, translation-inv check)
# ==============================================================================


# ==============================================================================
# §1  NUMERICAL / MATRIX UTILITIES
# ==============================================================================

#' Numerically stable logistic function
#'
#' Clips the linear predictor to [-35, 35] before applying sigma to avoid
#' exp() overflow/underflow. Use whenever you need sigma(X %*% alpha).
#'
#' @param eta  numeric vector of linear predictors
#' @return     vector of probabilities in (0, 1)
safe_logit <- function(eta) {
  eta <- pmin(pmax(eta, -35), 35)
  1 / (1 + exp(-eta))
}


#' Centre and scale non-intercept columns of a design matrix
#'
#' Improves numerical conditioning for logit/CBPS M-estimation, especially
#' with polynomial terms (age, age^2, age^3). Because the intercept is kept
#' fixed at 1 this is a pure reparameterisation: fitted propensity scores in
#' the population are unchanged.
#'
#' @param X  numeric matrix (may contain an intercept column of all 1s)
#' @return   matrix with non-intercept columns standardised
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


#' Numerical Jacobian (central differences)
#'
#' Used in sandwich_se_mest() to avoid hand-coding the Jacobian of the
#' stacked moment vector. Central differences with adaptive step size.
#'
#' @param f      function R^k -> R^m
#' @param theta  parameter vector (length k)
#' @param eps    base step-size
#' @return       m x k Jacobian matrix
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


#' Safe matrix inverse with fallback chain
#'
#' Tries solve() -> qr.solve() -> ridge-regularised solve() -> SVD
#' pseudo-inverse. Returns the first finite inverse found.
#'
#' @param A    square numeric matrix
#' @param tol  relative tolerance for SVD truncation
#' @return     inverse (or pseudo-inverse) of A
matrix_inverse_safe <- function(A, tol = 1e-10) {
  inv <- tryCatch(solve(A),      error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)

  inv <- tryCatch(qr.solve(A),   error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)

  for (ridge in c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
    inv <- tryCatch(solve(A + ridge * diag(ncol(A))), error = function(e) NULL)
    if (!is.null(inv) && all(is.finite(inv))) return(inv)
  }

  # Moore-Penrose via SVD
  sv     <- svd(A)
  d      <- sv$d
  d_inv  <- ifelse(d > tol * max(d), 1 / d, 0)
  sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
}


# ==============================================================================
# §2  PROPENSITY-SCORE ESTIMATORS
# ==============================================================================
# Both estimators return p = P(Z = 1 | X_i) for every observation.
# Swap the generic name "Z" for your instrument variable.
# The design matrix X should INCLUDE an intercept column (model.matrix style).

#' Logit MLE propensity score
#'
#' Standard glm(Z ~ X - 1, binomial). The "-1" suppresses a duplicate
#' intercept when X already contains a column of 1s.
#'
#' @param Z  binary instrument vector (0/1)
#' @param X  design matrix (with intercept column)
#' @return   fitted P(Z = 1 | X_i) for each observation
logit_mle <- function(Z, X) {
  df  <- data.frame(Z = Z, X)
  fit <- glm(Z ~ . - 1, data = df, family = binomial(link = "logit"))
  fitted.values(fit)
}


#' Fit logit MLE and return coefficients + fitted p
#'
#' Used internally by kappa_analytic_se_one(). Standardises X for numerical
#' stability; see prep_design_for_mest().
#'
#' @param Z  binary instrument vector
#' @param X  design matrix
#' @return   list(alpha = coefficient vector, p = fitted probabilities,
#'                X_used = standardised design matrix)
fit_logit_alpha <- function(Z, X) {
  Z     <- as.numeric(Z)
  X     <- prep_design_for_mest(X)
  fit   <- glm.fit(x = X, y = Z, family = binomial(link = "logit"))
  alpha <- as.numeric(coef(fit));  alpha[is.na(alpha)] <- 0
  p     <- as.vector(pmin(pmax(safe_logit(X %*% alpha), 1e-8), 1 - 1e-8))
  list(alpha = alpha, p = p, X_used = X)
}


#' Covariate Balancing Propensity Score (CBPS)
#'
#' Solves the moment condition
#'   E[(Z - p(X)) / {p(X)(1-p(X))} * X] = 0
#' via Newton steps with line search.  Initialised at the logit MLE.
#' Returns the best iterate (lowest max-moment norm) if the tolerance
#' is not reached within max_iter steps.
#'
#' @param Z         binary instrument vector
#' @param X         design matrix (with intercept)
#' @param tol       convergence tolerance on max |moment|
#' @param max_iter  maximum Newton iterations
#' @return  list(alpha, p, X_used, converged, max_moment)
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


#' Convenience wrapper — drop-in replacement for old cbps(Z, X)$p
cbps <- function(Z, X, tol = 1e-9, max_iter = 5000, verbose = FALSE) {
  fit_cbps_alpha(Z, X, tol = tol, max_iter = max_iter)
}

#' Safely extract P(Z=1|X) from cbps() output regardless of return type
#'
#' @param Z  instrument vector
#' @param X  design matrix (with intercept)
get_cbps_p <- function(Z, X) {
  out <- cbps(Z, X)
  if (is.list(out) && !is.null(out$p)) return(as.vector(out$p))
  as.vector(out)
}


# ==============================================================================
# §3  KAPPA WEIGHTS
# ==============================================================================
# Abadie (2003): complier moments can be recovered by weighting with
#
#   kappa  = 1 - D(1-Z)/(1-p) - (1-D)Z/p
#   kappa1 = D(Z-p) / [p(1-p)]
#   kappa0 = (1-D)((1-Z)-(1-p)) / [p(1-p)]
#
# In population E[kappa] = E[kappa1] = E[kappa0] = P(D(1) > D(0)).
# Variable names are generic: Z = instrument, D = treatment, p = P(Z=1|X).

#' Compute the three Abadie kappa weight vectors
#'
#' @param Z  binary instrument  (0/1)
#' @param D  binary treatment   (0/1)
#' @param p  estimated P(Z=1|X_i)
#' @return   list(kappa, kappa1, kappa0) — each a length-n numeric vector
kappa_weights <- function(Z, D, p) {
  list(
    kappa  = 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p,
    kappa1 = D * (Z - p) / (p * (1 - p)),
    kappa0 = (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  )
}


# ==============================================================================
# §4  LATE POINT ESTIMATORS
# ==============================================================================
# All estimators target  LATE = E[Y(1)-Y(0) | D(1)>D(0)].
# The difference between them is normalization of the denominator.

#' tau_u — Uysal normalized estimator  [RECOMMENDED; translation invariant]
#'
#' Computes normalized IPW means for each instrument arm and takes their ratio:
#'   tau_u = [E_hat(Y|Z=1) - E_hat(Y|Z=0)] / [E_hat(D|Z=1) - E_hat(D|Z=0)]
#' where E_hat(·|Z=z) = weighted average with weights proportional to Z/p or (1-Z)/(1-p).
#'
#' @param Y  outcome vector
#' @param Z  binary instrument
#' @param D  binary treatment
#' @param p  estimated instrument propensity score
#' @return   scalar LATE estimate
tau_u <- function(Y, Z, D, p) {
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  numerator   <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denominator <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  numerator / denominator
}


#' tau_a10 — normalized Abadie-Cattaneo estimator  [translation invariant]
#'
#' Separately normalises the kappa1 and kappa0 weighted outcome means:
#'   tau_a10 = sum(kappa1 * Y)/sum(kappa1) - sum(kappa0 * Y)/sum(kappa0)
#'
#' @inheritParams tau_u
tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}


#' tau_unnorm — unnormalized kappa estimators  [NOT translation invariant]
#'
#' The common numerator is mean[Y*(Z-p)/(p(1-p))].
#' Three denominator choices (SUW notation):
#'   "a"  : mean(kappa)   — Abadie (2003) original
#'   "a1" : mean(kappa1)  — Tan/Frolich (also tau_t in SUW)
#'   "a0" : mean(kappa0)
#'
#' These are fragile when outcome units change (e.g. dollars vs cents).
#'
#' @inheritParams tau_u
#' @param which  one of "a", "a1", "a0"
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


#' Compute all six kappa point estimates at once
#'
#' Convenience wrapper: returns a named numeric vector.
#'
#' @param Y         outcome vector
#' @param Z         instrument vector
#' @param D         treatment vector
#' @param X_kappa   design matrix WITH intercept (for logit_mle / cbps)
#' @return          named numeric vector of six estimates
kappa_estimates <- function(Y, Z, D, X_kappa) {
  p_ml <- logit_mle(Z, X_kappa)
  p_cb <- get_cbps_p(Z, X_kappa)
  c(
    tau_cb_u   = tau_u(Y, Z, D, p_cb),
    tau_ml_u   = tau_u(Y, Z, D, p_ml),
    tau_ml_a10 = tau_a10(Y, Z, D, p_ml),
    tau_ml_a   = tau_unnorm(Y, Z, D, p_ml, "a"),
    tau_ml_a1  = tau_unnorm(Y, Z, D, p_ml, "a1"),
    tau_ml_a0  = tau_unnorm(Y, Z, D, p_ml, "a0")
  )
}


# ==============================================================================
# §5  OUTCOME-WEIGHT CONSTRUCTORS
# ==============================================================================
# For any weighted-outcome representation  tau_hat = sum_i omega_i * Y_i,
# the per-observation weights omega_i are listed below.
# Note: sum(omega) = 0 <=> translation invariance.

#' Per-observation outcome weights for all five kappa estimators
#'
#' Returns a list of length-n weight vectors.  Use these to:
#'   - verify sum(w * Y) == estimate  (algebraic identity)
#'   - check translation invariance:  sum(w) should be ~0 for tau_u and tau_a10
#'   - compute ESS = 1/sum(w^2)
#'   - create Love plots via cobalt
#'
#' @param Z  binary instrument
#' @param D  binary treatment
#' @param p  estimated instrument propensity score
#' @return   list(w_u, w_a10, w_a, w_a1, w_a0)
kappa_outcome_weights <- function(Z, D, p) {
  n  <- length(Z)
  kw <- kappa_weights(Z, D, p)

  # tau_u weights  (Uysal normalized)
  s1  <- sum(Z / p)
  s0  <- sum((1 - Z) / (1 - p))
  dD  <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  w_u <- (Z / p / s1 - (1 - Z) / (1 - p) / s0) / dD

  # tau_a10 weights  (Abadie-Cattaneo normalized)
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


#' Bundle all kappa outcome weights for a given design matrix
#'
#' Computes both the MLE and CBPS propensity scores and returns all weights
#' as a named list — convenient for diagnostic tables and Love plots.
#'
#' @param Z        instrument vector
#' @param D        treatment vector
#' @param X_kappa  design matrix WITH intercept
#' @return         named list of length-n weight vectors
kappa_weights_bundle <- function(Z, D, X_kappa) {
  p_ml  <- logit_mle(Z, X_kappa)
  p_cb  <- get_cbps_p(Z, X_kappa)
  kw_ml <- kappa_outcome_weights(Z, D, p_ml)
  kw_cb <- kappa_outcome_weights(Z, D, p_cb)
  list(
    tau_cb_u   = kw_cb$w_u,
    tau_ml_u   = kw_ml$w_u,
    tau_ml_a10 = kw_ml$w_a10,
    tau_ml_a   = kw_ml$w_a,
    tau_ml_a1  = kw_ml$w_a1,
    tau_ml_a0  = kw_ml$w_a0
  )
}


# ==============================================================================
# §6  SANDWICH / M-ESTIMATION STANDARD ERRORS
# ==============================================================================
# The SUW paper treats each kappa estimator as an M-estimator and derives
# analytical SEs from the sandwich formula:
#
#   Avar(sqrt(n) * (theta_hat - theta_0)) = A^{-1} V (A^{-1})'
#
# where A = E[d/dtheta psi(O_i, theta)] and V = Var(psi(O_i, theta_0)).
# Numerical Jacobians (num_jacobian) are used throughout.

#' Helper: alpha-equation moments (score of the propensity model)
#'
#' Returns an n x k matrix of per-observation propensity score moments.
#' These are stacked with the estimator moments in kappa_analytic_se_one().
#'
#' @param Z       instrument vector
#' @param p       fitted propensity scores
#' @param X_used  standardised design matrix
#' @param method  "ml" (logit score) or "cb" (CBPS score)
alpha_moment_matrix <- function(Z, p, X_used, method) {
  p <- as.vector(p);  Z <- as.vector(Z)
  if (method == "ml") {
    return(as.vector(Z - p) * X_used)         # logit ML score: (Z-p)*X
  }
  if (method == "cb") {
    return(as.vector((Z - p) / (p * (1 - p))) * X_used)  # CBPS moment
  }
  stop("method must be 'ml' or 'cb'")
}


#' Sandwich SE from a stacked moment function (M-estimation)
#'
#' @param moment_matrix_fn  function(theta) -> n x dim(theta) matrix of
#'                           per-observation moment contributions psi_i(theta)
#' @param theta_hat          point estimate of the full parameter vector
#' @param tau_index          integer — which component of theta is tau (the LATE)
#' @return                   scalar standard error for theta[tau_index]
sandwich_se_mest <- function(moment_matrix_fn, theta_hat, tau_index) {
  theta_hat   <- as.numeric(theta_hat)
  psi_hat     <- moment_matrix_fn(theta_hat)
  n           <- nrow(psi_hat)

  A           <- num_jacobian(function(th) colMeans(moment_matrix_fn(th)), theta_hat)
  psi_centered <- scale(psi_hat, center = TRUE, scale = FALSE)
  V           <- crossprod(psi_centered) / n

  A_inv       <- matrix_inverse_safe(A)
  vcov_theta  <- A_inv %*% V %*% t(A_inv) / n

  se2 <- vcov_theta[tau_index, tau_index]
  if (!is.finite(se2)) return(NA_real_)
  sqrt(abs(se2))
}


#' Analytical M-estimation SE for one kappa estimator
#'
#' Stacks the propensity-score moments with the LATE estimator moments,
#' then calls sandwich_se_mest(). Supported estimators:
#'   "u"   — tau_u  (Uysal normalized)
#'   "a10" — tau_a10
#'   "a"   — tau_a  (Abadie unnormalized, kappa denominator)
#'   "a1"  — tau_a1 (Tan/Frolich)
#'   "a0"  — tau_a0
#'
#' @param Y          outcome vector
#' @param Z          instrument vector
#' @param D          treatment vector
#' @param X          design matrix WITH intercept
#' @param estimator  character, one of "u", "a10", "a", "a1", "a0"
#' @param method     "ml" or "cb" (propensity score fitting method)
#' @return           scalar SE
kappa_analytic_se_one <- function(Y, Z, D, X, estimator, method = "ml") {
  Y <- as.numeric(Y);  Z <- as.numeric(Z)
  D <- as.numeric(D);  X <- as.matrix(X)
  n <- length(Y)

  fit <- if (method == "ml") fit_logit_alpha(Z, X) else fit_cbps_alpha(Z, X)
  alpha_hat <- fit$alpha
  X_used    <- fit$X_used
  k         <- length(alpha_hat)

  am <- function(p) alpha_moment_matrix(Z, p, X_used, method)

  # ---------------------------------------------------------------------------
  # tau_u: theta = (alpha, mu1, mu0, m1, m0, tau)
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
  # tau_a, tau_a1, tau_a0: theta = (alpha, Delta, Gamma, tau)
  # ---------------------------------------------------------------------------
  if (estimator %in% c("a", "a1", "a0")) {
    p   <- fit$p
    kw  <- kappa_weights(Z, D, p)
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
      cbind(am(p), psi_Delta = psi_Delta,
            psi_Gamma = psi_Gamma,
            psi_tau   = Delta / Gamma - tau)
    }
    return(sandwich_se_mest(moment_fn, theta_hat, length(theta_hat)))
  }

  # ---------------------------------------------------------------------------
  # tau_a10: theta = (alpha, Delta1, Gamma1, Delta0, Gamma0, tau)
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


#' Error-safe wrapper around kappa_analytic_se_one
#'
#' Returns NA_real_ (with a warning) if the matrix system is singular rather
#' than stopping the script mid-table.
safe_kappa_se <- function(Y, Z, D, X_mat, estimator, method) {
  tryCatch(
    kappa_analytic_se_one(Y, Z, D, X_mat, estimator = estimator, method = method),
    error = function(e) {
      warning(sprintf("SE failed: estimator=%s, method=%s. %s", estimator, method, e$message))
      NA_real_
    }
  )
}


#' Compute all six kappa point estimates AND their analytical SEs
#'
#' @param Y      outcome vector
#' @param Z      instrument vector
#' @param D      treatment vector
#' @param X_mat  design matrix WITH intercept
#' @return       list(estimates = named numeric(6), se = named numeric(6))
kappa_analytic_se_all <- function(Y, Z, D, X_mat) {
  Y <- as.numeric(Y); Z <- as.numeric(Z)
  D <- as.numeric(D); X_mat <- as.matrix(X_mat)

  p_ml <- logit_mle(Z, X_mat)
  p_cb <- cbps(Z, X_mat)$p

  estimates <- c(
    tau_cb_u   = tau_u(Y, Z, D, p_cb),
    tau_ml_u   = tau_u(Y, Z, D, p_ml),
    tau_ml_a10 = tau_a10(Y, Z, D, p_ml),
    tau_ml_a   = tau_unnorm(Y, Z, D, p_ml, "a"),
    tau_ml_t   = tau_unnorm(Y, Z, D, p_ml, "a1"),
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
# §7  2SLS BENCHMARK
# ==============================================================================

#' Run 2SLS with HC1-robust standard errors
#'
#' Uses AER::ivreg + sandwich::vcovHC(type="HC1").
#' HC1 = (n/(n-k)) * HC0 — matches Stata's ivreg2 vce(robust) default.
#'
#' @param Y          outcome vector
#' @param D          treatment vector
#' @param Z          instrument vector
#' @param X_df       data.frame of control variables (WITHOUT intercept column)
#' @param endog_name character name used for D in the formula
#' @return           list(coef, se, pval) for the endogenous variable
run_2sls <- function(Y, D, Z, X_df, endog_name = "D") {
  # Requires: library(AER); library(sandwich); library(lmtest)
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
# §8  WEIGHT DIAGNOSTICS
# ==============================================================================

#' Summary diagnostics for a weight vector
#'
#' Key diagnostics:
#'   Sum_w     : should be ~0 for translation-invariant estimators
#'   ESS       : effective sample size = 1 / sum(w_i^2)
#'   Pct_neg   : % of observations with negative weight
#'   Max_abs_w : largest absolute weight (outlier detection)
#'
#' @param w     numeric weight vector
#' @param name  label for the Estimator column
#' @return      one-row data.frame
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


#' Algebraic identity check: sum(w * Y) == estimate
#'
#' @param w        weight vector
#' @param Y        outcome vector
#' @param estimate scalar point estimate
#' @param tol      numerical tolerance
#' @return         logical
check_weight_identity <- function(w, Y, estimate, tol = 1e-8) {
  isTRUE(all.equal(sum(w * Y), estimate, tolerance = tol))
}


# ==============================================================================
# §9  DML HELPERS  (OutcomeWeights / DoubleML wrappers)
# ==============================================================================
# These functions assume OutcomeWeights (dev version from GitHub) and DoubleML
# are loaded. They wrap the package calls with the checks and extraction
# patterns established in the Vietnam replication.

#' Strip non-numeric attributes from model.matrix output
#'
#' grf::validate_X() rejects matrices with factor-level attributes.
#' This wrapper produces a clean numeric matrix safe for all ML learners.
#' NOTE: requires ncol >= 2 (single-column matrices fail grf in R 4.5.x).
#'
#' @param formula  one-sided formula, e.g. ~ 0 + age + I(age^2)
#' @param data     data.frame
#' @return         plain numeric matrix with column names
make_X <- function(formula, data) {
  mm <- model.matrix(formula, data = data)
  m  <- matrix(as.numeric(mm), nrow = nrow(mm), ncol = ncol(mm))
  colnames(m) <- colnames(mm)
  m
}


#' Extract a point estimate from OutcomeWeights summary() output
#'
#' summary(dml_with_smoother(...)) can return either "Estimate" or the first
#' column as the estimate column depending on version. This handles both.
#'
#' @param res       matrix returned by summary()
#' @param row_name  character, e.g. "PLR-IV" or "Wald-AIPW"
get_estimate <- function(res, row_name) {
  if ("Estimate" %in% colnames(res)) return(as.numeric(res[row_name, "Estimate"]))
  as.numeric(res[row_name, 1])
}


#' Extract scalar coefficient from a DoubleML object
get_dml_coef <- function(obj) {
  as.numeric(obj$coef)
}


#' Verify omega'Y == tau_hat for DoubleML outcome weights
#'
#' @param omega_obj  object from get_outcome_weights(iivm, dml_data)
#' @param Y          outcome vector
#' @param coef       scalar point estimate from iivm$coef
#' @param tol        tolerance
#' @return           logical
check_doubleml_identity <- function(omega_obj, Y, coef, tol = 1e-8) {
  isTRUE(all.equal(as.numeric(omega_obj$omega %*% Y),
                   as.numeric(coef), tolerance = tol))
}


#' Check that expected rows are present in an OutcomeWeights omega object
#'
#' @param omega_obj  object from get_outcome_weights()
#' @param rows       character vector of required row names
check_omega_rows <- function(omega_obj, rows = c("PLR-IV", "Wald-AIPW")) {
  missing <- setdiff(rows, rownames(omega_obj$omega))
  if (length(missing) > 0)
    stop("Missing expected outcome-weight rows: ", paste(missing, collapse = ", "))
}


# ==============================================================================
# §10  OUTPUT HELPERS
# ==============================================================================

#' Format coefficient with significance stars and SE in parentheses
#'
#' Stars use two-sided z-test: *** p<0.01, ** p<0.05, * p<0.10.
#'
#' @param coef    point estimate
#' @param se      standard error
#' @param digits  decimal places
#' @return        character string, e.g. "0.153***\n(0.042)"
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


#' Build one row of a translation-invariance table
#'
#' @param name  estimator label
#' @param w     weight vector
#' @param Y_dol outcome in original units
#' @param Y_cnt outcome in shifted units (Y_dol + k)
#' @param k     shift constant, e.g. log(100)
#' @return      one-row data.frame
translation_row <- function(name, w, Y_dol, Y_cnt, k) {
  ed <- sum(w * Y_dol)
  ec <- sum(w * Y_cnt)
  data.frame(
    Estimator   = name,
    Dollars     = round(ed, 4),
    Cents       = round(ec, 4),
    Diff_actual = round(ec - ed, 6),
    Diff_pred   = round(sum(w) * k, 6),
    Match       = isTRUE(all.equal(ec - ed, sum(w) * k, tolerance = 1e-8)),
    stringsAsFactors = FALSE
  )
}


#' Print a full translation-invariance comparison table
#'
#' Lists DML weights (w_plriv, w_waldaipw) next to all kappa estimators
#' for a given covariate specification.
#'
#' @param spec_label  character label for the table header
#' @param w_plriv     PLR-IV weight vector from OutcomeWeights
#' @param w_waldaipw  Wald-AIPW weight vector from OutcomeWeights
#' @param Z           instrument vector
#' @param D           treatment vector
#' @param X_kappa     kappa design matrix WITH intercept
#' @param Y_dol       outcome in base units
#' @param Y_cnt       outcome in shifted units
#' @param k           shift constant
print_inv_table <- function(spec_label, w_plriv, w_waldaipw,
                            Z, D, X_kappa, Y_dol, Y_cnt, k) {
  cat(sprintf("--- %s ---\n", spec_label))
  cat(sprintf("  %-22s  %9s  %9s  %10s  %10s  %7s\n",
              "Estimator", "Dollars", "Cents", "Diff(act)", "Diff(pred)", "Match?"))
  cat(strrep("-", 78), "\n")

  p_ml  <- logit_mle(Z, X_kappa)
  p_cb  <- get_cbps_p(Z, X_kappa)
  kw_ml <- kappa_outcome_weights(Z, D, p_ml)
  kw_cb <- kappa_outcome_weights(Z, D, p_cb)

  rows <- list(
    list(name = "PLR-IV",      w = w_plriv,
         ed = sum(w_plriv * Y_dol),    ec = sum(w_plriv * Y_cnt)),
    list(name = "Wald-AIPW",   w = w_waldaipw,
         ed = sum(w_waldaipw * Y_dol), ec = sum(w_waldaipw * Y_cnt)),
    list(name = "tau_cb_u",    w = kw_cb$w_u,
         ed = tau_u(Y_dol, Z, D, p_cb),   ec = tau_u(Y_cnt, Z, D, p_cb)),
    list(name = "tau_ml_u",    w = kw_ml$w_u,
         ed = tau_u(Y_dol, Z, D, p_ml),   ec = tau_u(Y_cnt, Z, D, p_ml)),
    list(name = "tau_ml_a10",  w = kw_ml$w_a10,
         ed = tau_a10(Y_dol, Z, D, p_ml), ec = tau_a10(Y_cnt, Z, D, p_ml)),
    list(name = "tau_ml_a",    w = kw_ml$w_a,
         ed = tau_unnorm(Y_dol, Z, D, p_ml, "a"),
         ec = tau_unnorm(Y_cnt, Z, D, p_ml, "a")),
    list(name = "tau_ml_a1/t", w = kw_ml$w_a1,
         ed = tau_unnorm(Y_dol, Z, D, p_ml, "a1"),
         ec = tau_unnorm(Y_cnt, Z, D, p_ml, "a1")),
    list(name = "tau_ml_a0",   w = kw_ml$w_a0,
         ed = tau_unnorm(Y_dol, Z, D, p_ml, "a0"),
         ec = tau_unnorm(Y_cnt, Z, D, p_ml, "a0"))
  )

  for (r in rows) {
    act  <- r$ec - r$ed
    pred <- sum(r$w) * k
    ok   <- isTRUE(all.equal(act, pred, tolerance = 1e-8))
    cat(sprintf("  %-22s  %9.4f  %9.4f  %10.4f  %10.4f  %7s\n",
                r$name, r$ed, r$ec, act, pred, if (ok) "TRUE" else "FALSE"))
  }
  cat("\n")
}


#' Build a Love-plot for covariate balance (cobalt::love.plot wrapper)
#'
#' Signed outcome weights are mapped to cobalt's treated-vs-control convention
#' by multiplying by (2D - 1). The covariate matrix X must be supplied.
#'
#' @param title_str  plot title
#' @param w_vec      outcome weight vector
#' @param D          binary treatment vector
#' @param X          covariate matrix (no intercept; passed to cobalt)
#' @return           ggplot2 object from love.plot()
make_love <- function(title_str, w_vec, D, X) {
  # Requires: library(cobalt); library(viridis)
  cobalt::love.plot(
    D ~ X,
    weights    = w_vec * (2 * D - 1),
    position   = "bottom",
    title      = title_str,
    thresholds = c(m = 0.1),
    var.order  = "unadjusted",
    binary     = "std",
    abs        = TRUE,
    line       = TRUE,
    colors     = viridis::viridis(2),
    shapes     = c("circle", "triangle")
  )
}


# ==============================================================================
# END OF TOOLKIT
# ==============================================================================
# PORTING CHECKLIST for a new design (e.g. Imbens 401k):
#
#  1. Replace Z with your instrument, D with treatment, Y with outcome.
#  2. Build X matrices via model.matrix(); include intercept for logit/CBPS,
#     omit intercept for grf/OutcomeWeights (use make_X()).
#  3. Call kappa_analytic_se_all(Y, Z, D, X_kappa) for the full kappa table.
#  4. For DML: use dml_with_smoother() or DoubleMLIIVM$new(), then
#     get_outcome_weights() to extract weights.
#  5. Use weight_diag() to compare Sum_w, ESS, Pct_neg across estimators.
#  6. Use translation_row() / print_inv_table() to verify/demonstrate
#     translation invariance.
#  7. Use make_love() for covariate balance Love plots.
# ==============================================================================
