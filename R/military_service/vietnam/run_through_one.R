# ==============================================================================
# ANGRIST (1990) — MILITARY SERVICE DESIGN
# Replication of Table 2 in Słoczyński, Uysal, and Wooldridge (2025, JBES)
# "Abadie's Kappa and Weighting Estimators of the Local Average Treatment Effect"
# ==============================================================================
#
# Object of interest:
#   LATE = E[Y(1) - Y(0) | D(1) > D(0)]
#
# Data:
#   Survey of Income and Program Participation, as used in Angrist (1990).
#   The final sample excludes rsncode == 999 and observations with missing
#   education or wage information. The target sample size is N = 3,027.
#
# Instrument:
#   Z = 1{draft eligible}, based on the Vietnam draft lottery number.
#
# Treatment:
#   D = 1{veteran}, indicating military service.
#
# Outcome:
#   Y = log wage. The paper reports results when wages are measured in:
#     (i) cents before taking logs, and
#     (ii) dollars before taking logs.
#   Since log(100 * wage) = log(wage) + log(100), this comparison tests
#   translation invariance of the estimators.
#
# Covariate specifications:
#   1. Linear age
#   2. Cubic polynomial in age
#   3. Saturated age indicators
#
# Identification assumptions:
#   conditional independence of Z given X, exclusion, overlap/first stage,
#   and monotonicity. These are the standard LATE assumptions.
# ==============================================================================
# ==============================================================================
# 0. PACKAGES
# ==============================================================================
# haven     : read Stata .dta files
# AER       : ivreg() for the 2SLS benchmark
# sandwich  : HC1 robust variance for 2SLS
# lmtest    : coeftest() for extracting robust 2SLS results

library(haven)
library(AER)
library(sandwich)
library(lmtest)

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
# 2. INSTRUMENT PROPENSITY SCORE p(X) = P(Z = 1 | X)
# ==============================================================================
# The kappa estimators require the conditional probability of instrument
# assignment. The paper estimates p(X) in two ways:
#
#   1. Logit MLE:
#      standard logistic maximum likelihood.
#
#   2. CBPS:
#      covariate balancing moment conditions, corresponding to the paper's
#      moment-based propensity score estimator.
#
# The reported Table 2 rows use:
#   tau_cb_u for the CBPS version of tau_u,
#   tau_ml_* for all logit-MLE-based estimators.
# ==============================================================================

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





# ==============================================================================
# 3. KAPPA WEIGHTS
# ==============================================================================
# Abadie (2003) shows that complier moments can be recovered using weights
# that depend on the instrument propensity score p(X) = P(Z = 1 | X).
#
# The three weights used by SUW are:
#
#   kappa  = 1 - D(1-Z)/(1-p) - (1-D)Z/p
#   kappa1 = D(Z-p) / [p(1-p)]
#   kappa0 = (1-D)((1-Z)-(1-p)) / [p(1-p)]
#
# In population:
#   E[kappa] = E[kappa1] = E[kappa0] = P(D(1) > D(0)).
#
# In finite samples, their sample means need not be equal. This finite-sample
# difference is one reason why normalized and unnormalized estimators can
# behave differently.
# ==============================================================================

kappa_weights <- function(Z, D, p) {
  list(
    kappa  = 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p,
    kappa1 = D * (Z - p) / (p * (1 - p)),
    kappa0 = (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  )
}


# ==============================================================================
# 4. LATE ESTIMATORS
# ==============================================================================
# The paper reports six kappa-based rows:
#
#   tau_cb_u   : Uysal normalized estimator using CBPS p-score
#   tau_ml_u   : Uysal normalized estimator using logit MLE p-score
#   tau_ml_a10 : normalized Abadie-Cattaneo estimator using kappa1 and kappa0
#   tau_ml_a   : unnormalized Abadie estimator
#   tau_ml_t   : Tan/Frolich estimator; in SUW notation tau_t = tau_a,1
#   tau_ml_a0  : unnormalized kappa0-denominator estimator
#
# The key distinction is normalization:
#
#   Normalized estimators:
#     tau_u and tau_a10. These are translation invariant, so adding a constant
#     to Y does not change the estimate.
#
#   Unnormalized estimators:
#     tau_a, tau_t = tau_a1, and tau_a0. These are generally not translation
#     invariant. Therefore, when Y is a log outcome, estimates can change when
#     wages are measured in cents rather than dollars before taking logs.
#
# SUW recommend tau_u for wider use because it is normalized and has favorable
# finite-sample properties, especially under one-sided noncompliance.
# ==============================================================================

# ------------------------------------------------------------------------------
# 4a. tau_u — Uysal normalized estimator
# ------------------------------------------------------------------------------
# This estimator forms inverse-propensity weighted means of Y and D by instrument
# status and takes the ratio:
#
#   tau_u = [mu_Y1 - mu_Y0] / [mu_D1 - mu_D0],
#
# where mu_Y1 and mu_Y0 are normalized IPW means for Z = 1 and Z = 0,
# and mu_D1 and mu_D0 are defined analogously for D.
#
# This is translation invariant because the additive constant cancels between
# the two normalized outcome means.

tau_u <- function(Y, Z, D, p) {
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  numerator   <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denominator <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  numerator / denominator
}

# ------------------------------------------------------------------------------
# 4b. tau_a10 — normalized Abadie-Cattaneo estimator
# ------------------------------------------------------------------------------
# This estimator separately normalizes the kappa1 and kappa0 weighted outcome
# means:
#
#   tau_a10 = sum(kappa1 * Y) / sum(kappa1)
#           - sum(kappa0 * Y) / sum(kappa0).
#
# Like tau_u, it is translation invariant.

tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}

# ------------------------------------------------------------------------------
# 4c. Unnormalized estimators: tau_a, tau_t = tau_a1, tau_a0
# ------------------------------------------------------------------------------
# These estimators use the common numerator
#
#   mean[ Y * (Z - p(X)) / {p(X)(1-p(X))} ],
#
# but divide by different sample analogues of the complier share:
#
#   tau_a  : denominator mean(kappa)
#   tau_t  : denominator mean(kappa1); this is tau_a1 in the paper
#   tau_a0 : denominator mean(kappa0)
#
# These estimators are generally not translation invariant. This is exactly what
# Table 2 illustrates when wages are measured in cents versus dollars.

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
# 6. ANALYTICAL STANDARD ERRORS FOR KAPPA ESTIMATORS
# ==============================================================================
# The paper does not use bootstrap standard errors for the kappa estimators.
# It treats each estimator as an M-estimator and computes standard errors from
# the sandwich variance formula in the online appendix.
#
# For a stacked moment vector psi(O_i, theta), the estimator solves:
#
#   N^{-1} sum_i psi(O_i, theta_hat) = 0.
#
# The asymptotic variance is based on:
#
#   A^{-1} V A^{-1}',
#
# where A is the Jacobian of the mean moment vector and V is the variance of the
# moment functions.
#
# The code below stacks the appropriate propensity-score moments and estimator
# moments for each kappa estimator. Numerical derivatives are used for the
# Jacobian, with stable matrix-inversion fallbacks for the cubic-age design.
# ==============================================================================

safe_logit <- function(eta) {
  eta <- pmin(pmax(eta, -35), 35)
  1 / (1 + exp(-eta))
}

# Centering and scaling non-intercept columns improves numerical conditioning
# for the cubic-age specification. Because the model includes an intercept, this
# is only a reparameterization of the same logit/CBPS problem and does not change
# the fitted propensity scores in population.

prep_design_for_mest <- function(X) {
  X <- as.matrix(X)
  
  # Detect intercept column: all values equal to 1
  is_intercept <- apply(X, 2, function(v) all(abs(v - 1) < 1e-12))
  
  X_new <- X
  
  for (j in seq_len(ncol(X))) {
    if (!is_intercept[j]) {
      mu <- mean(X[, j])
      sdj <- sd(X[, j])
      
      if (is.finite(sdj) && sdj > 1e-12) {
        X_new[, j] <- (X[, j] - mu) / sdj
      }
    }
  }
  
  X_new
}

fit_logit_alpha <- function(Z, X) {
  Z <- as.numeric(Z)
  X <- prep_design_for_mest(X)
  
  fit <- glm.fit(
    x = X,
    y = Z,
    family = binomial(link = "logit")
  )
  
  alpha <- as.numeric(coef(fit))
  alpha[is.na(alpha)] <- 0
  
  p <- as.vector(pmin(pmax(safe_logit(X %*% alpha), 1e-8), 1 - 1e-8))
  
  list(alpha = alpha, p = p, X_used = X)
}


fit_cbps_alpha <- function(Z, X, tol = 1e-9, max_iter = 5000) {
  Z <- as.numeric(Z)
  X <- prep_design_for_mest(X)
  n <- length(Z)
  k <- ncol(X)
  
  # Start from ML logit on the same scaled X
  start <- tryCatch(
    glm.fit(x = X, y = Z, family = binomial(link = "logit"))$coefficients,
    error = function(e) rep(0, k)
  )
  start[is.na(start)] <- 0
  b <- as.numeric(start)
  
  moment_fn <- function(b) {
    p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
    w <- as.vector((Z - p) / (p * (1 - p)))
    colMeans(w * X)
  }
  
  jac_fn <- function(b) {
    p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
    
    # derivative of (Z-p)/(p(1-p)) wrt alpha for logit p
    # equals: - [Z(1-p)/p + (1-Z)p/(1-p)] * X
    w <- as.vector(-Z * (1 - p) / p - (1 - Z) * p / (1 - p))
    crossprod(X, w * X) / n
  }
  
  best_b <- b
  best_norm <- max(abs(moment_fn(b)))
  converged <- FALSE
  
  for (iter in seq_len(max_iter)) {
    m <- moment_fn(b)
    m_norm <- max(abs(m))
    
    if (m_norm < best_norm) {
      best_norm <- m_norm
      best_b <- b
    }
    
    if (m_norm < tol) {
      converged <- TRUE
      best_b <- b
      break
    }
    
    J <- jac_fn(b)
    
    step <- tryCatch(
      qr.solve(J, -m),
      error = function(e) NULL
    )
    
    if (is.null(step) || any(!is.finite(step))) {
      for (ridge in c(1e-10, 1e-8, 1e-6, 1e-4)) {
        step <- tryCatch(
          solve(J + ridge * diag(k), -m),
          error = function(e) NULL
        )
        if (!is.null(step) && all(is.finite(step))) break
      }
    }
    
    if (is.null(step) || any(!is.finite(step))) break
    
    alpha_step <- 1
    accepted <- FALSE
    
    for (j in seq_len(50)) {
      b_new <- b + alpha_step * step
      new_norm <- max(abs(moment_fn(b_new)))
      
      if (is.finite(new_norm) && new_norm < m_norm) {
        b <- b_new
        accepted <- TRUE
        break
      }
      
      alpha_step <- alpha_step * 0.5
    }
    
    if (!accepted) break
  }
  
  b <- best_b
  p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
  
  list(
    alpha = as.numeric(b),
    p = as.vector(p),
    X_used = X,
    converged = converged,
    max_moment = best_norm
  )
}


# Keep this wrapper so old downstream code using cbps(...)$p still works.
cbps <- function(Z, X, tol = 1e-9, max_iter = 5000, verbose = FALSE) {
  fit_cbps_alpha(Z, X, tol = tol, max_iter = max_iter)
}


num_jacobian <- function(f, theta, eps = 1e-6) {
  theta <- as.numeric(theta)
  f0 <- f(theta)
  m <- length(f0)
  k <- length(theta)
  J <- matrix(NA_real_, nrow = m, ncol = k)
  
  for (j in seq_len(k)) {
    h <- eps * (abs(theta[j]) + 1)
    th_plus <- theta
    th_minus <- theta
    th_plus[j] <- th_plus[j] + h
    th_minus[j] <- th_minus[j] - h
    
    J[, j] <- (f(th_plus) - f(th_minus)) / (2 * h)
  }
  
  J
}


matrix_inverse_safe <- function(A, tol = 1e-10) {
  # 1. Try exact inverse
  inv <- tryCatch(
    solve(A),
    error = function(e) NULL
  )
  if (!is.null(inv) && all(is.finite(inv))) return(inv)
  
  # 2. Try QR inverse
  inv <- tryCatch(
    qr.solve(A),
    error = function(e) NULL
  )
  if (!is.null(inv) && all(is.finite(inv))) return(inv)
  
  # 3. Try ridge-regularized inverse
  ridge_grid <- c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
  
  for (ridge in ridge_grid) {
    inv <- tryCatch(
      solve(A + ridge * diag(ncol(A))),
      error = function(e) NULL
    )
    
    if (!is.null(inv) && all(is.finite(inv))) return(inv)
  }
  
  # 4. Last fallback: Moore-Penrose pseudo-inverse via SVD
  sv <- svd(A)
  d <- sv$d
  
  cutoff <- tol * max(d)
  d_inv <- ifelse(d > cutoff, 1 / d, 0)
  
  inv <- sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
  inv
}


sandwich_se_mest <- function(moment_matrix_fn, theta_hat, tau_index) {
  theta_hat <- as.numeric(theta_hat)
  
  psi_hat <- moment_matrix_fn(theta_hat)
  n <- nrow(psi_hat)
  
  gbar_fn <- function(th) {
    colMeans(moment_matrix_fn(th))
  }
  
  A <- num_jacobian(gbar_fn, theta_hat)
  
  # Robust empirical variance of moment functions
  psi_centered <- scale(psi_hat, center = TRUE, scale = FALSE)
  V <- crossprod(psi_centered) / n
  
  A_inv <- matrix_inverse_safe(A)
  
  vcov_theta <- A_inv %*% V %*% t(A_inv) / n
  
  se2 <- vcov_theta[tau_index, tau_index]
  
  if (!is.finite(se2)) return(NA_real_)
  
  sqrt(abs(se2))
}




kappa_analytic_se_one <- function(Y, Z, D, X, estimator, method = "ml") {
  Y <- as.numeric(Y)
  Z <- as.numeric(Z)
  D <- as.numeric(D)
  X <- as.matrix(X)
  
  n <- length(Y)
  
  if (method == "ml") {
    fit <- fit_logit_alpha(Z, X)
    alpha_hat <- fit$alpha
    X_used <- fit$X_used
  } else if (method == "cb") {
    fit <- fit_cbps_alpha(Z, X)
    alpha_hat <- fit$alpha
    X_used <- fit$X_used
  } else {
    stop("method must be 'ml' or 'cb'")
  }
  
  k <- length(alpha_hat)
  
  alpha_moment <- function(Z, p, X_used, method) {
    p <- as.vector(p)
    Z <- as.vector(Z)
    
    if (method == "ml") {
      # logit ML score: (Z - p) X
      w <- as.vector(Z - p)
      return(w * X_used)
    }
    
    if (method == "cb") {
      # CBPS moment: (Z - p) / [p(1-p)] X
      w <- as.vector((Z - p) / (p * (1 - p)))
      return(w * X_used)
    }
    
    stop("Unknown method")
  }
  
  # ---------------------------------------------------------------------------
  # tau_u: normalized Uysal estimator
  # theta = (alpha, mu1, mu0, m1, m0, tau)
  # ---------------------------------------------------------------------------
  if (estimator == "u") {
    p <- fit$p
    
    mu1 <- sum(Z * Y / p) / sum(Z / p)
    mu0 <- sum((1 - Z) * Y / (1 - p)) / sum((1 - Z) / (1 - p))
    m1  <- sum(Z * D / p) / sum(Z / p)
    m0  <- sum((1 - Z) * D / (1 - p)) / sum((1 - Z) / (1 - p))
    tau <- (mu1 - mu0) / (m1 - m0)
    
    theta_hat <- c(alpha_hat, mu1, mu0, m1, m0, tau)
    tau_index <- length(theta_hat)
    
    moment_fn <- function(theta) {
      alpha <- theta[seq_len(k)]
      mu1 <- theta[k + 1]
      mu0 <- theta[k + 2]
      m1  <- theta[k + 3]
      m0  <- theta[k + 4]
      tau <- theta[k + 5]
      
      p <- as.vector(pmin(pmax(safe_logit(X_used %*% alpha), 1e-8), 1 - 1e-8))
      
      cbind(
        alpha_moment(Z, p, X_used, method),
        psi_mu1 = Z * (Y - mu1) / p,
        psi_mu0 = (1 - Z) * (Y - mu0) / (1 - p),
        psi_m1  = Z * (D - m1) / p,
        psi_m0  = (1 - Z) * (D - m0) / (1 - p),
        psi_tau = (mu1 - mu0) / (m1 - m0) - tau
      )
    }
    
    return(sandwich_se_mest(moment_fn, theta_hat, tau_index))
  }
  
  # ---------------------------------------------------------------------------
  # tau_a, tau_a1, tau_a0:
  # theta = (alpha, Delta, Gamma_variant, tau)
  # ---------------------------------------------------------------------------
  if (estimator %in% c("a", "a1", "a0")) {
    p <- fit$p
    kw <- kappa_weights(Z, D, p)
    
    Delta <- mean(Y * (Z - p) / (p * (1 - p)))
    
    Gamma <- switch(
      estimator,
      "a"  = mean(kw$kappa),
      "a1" = mean(kw$kappa1),
      "a0" = mean(kw$kappa0)
    )
    
    tau <- Delta / Gamma
    
    theta_hat <- c(alpha_hat, Delta, Gamma, tau)
    tau_index <- length(theta_hat)
    
    moment_fn <- function(theta) {
      alpha <- theta[seq_len(k)]
      Delta <- theta[k + 1]
      Gamma <- theta[k + 2]
      tau   <- theta[k + 3]
      
      p <- as.vector(pmin(pmax(safe_logit(X_used %*% alpha), 1e-8), 1 - 1e-8))
      
      psi_Delta <- Z * Y / p - (1 - Z) * Y / (1 - p) - Delta
      
      psi_Gamma <- switch(
        estimator,
        "a" = 1 - (1 - Z) * D / (1 - p) - Z * (1 - D) / p - Gamma,
        "a1" = Z * D / p - (1 - Z) * D / (1 - p) - Gamma,
        "a0" = Z * (D - 1) / p - (1 - Z) * (D - 1) / (1 - p) - Gamma
      )
      
      cbind(
        alpha_moment(Z, p, X_used, method),
        psi_Delta = psi_Delta,
        psi_Gamma = psi_Gamma,
        psi_tau   = Delta / Gamma - tau
      )
    }
    
    return(sandwich_se_mest(moment_fn, theta_hat, tau_index))
  }
  
  # ---------------------------------------------------------------------------
  # tau_a10:
  # theta = (alpha, Delta1, Gamma1, Delta0, Gamma0, tau)
  # ---------------------------------------------------------------------------
  if (estimator == "a10") {
    p <- fit$p
    kw <- kappa_weights(Z, D, p)
    
    Delta1 <- mean(kw$kappa1 * Y)
    Gamma1 <- mean(kw$kappa1)
    Delta0 <- mean(kw$kappa0 * Y)
    Gamma0 <- mean(kw$kappa0)
    tau <- Delta1 / Gamma1 - Delta0 / Gamma0
    
    theta_hat <- c(alpha_hat, Delta1, Gamma1, Delta0, Gamma0, tau)
    tau_index <- length(theta_hat)
    
    moment_fn <- function(theta) {
      alpha  <- theta[seq_len(k)]
      Delta1 <- theta[k + 1]
      Gamma1 <- theta[k + 2]
      Delta0 <- theta[k + 3]
      Gamma0 <- theta[k + 4]
      tau    <- theta[k + 5]
      
      p <- as.vector(pmin(pmax(safe_logit(X_used %*% alpha), 1e-8), 1 - 1e-8))
      
      kappa1_i <- D * (Z - p) / (p * (1 - p))
      kappa0_i <- (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
      
      psi_Delta1 <- kappa1_i * Y - Delta1
      psi_Gamma1 <- Z * D / p - (1 - Z) * D / (1 - p) - Gamma1
      
      psi_Delta0 <- kappa0_i * Y - Delta0
      psi_Gamma0 <- Z * (D - 1) / p - (1 - Z) * (D - 1) / (1 - p) - Gamma0
      
      cbind(
        alpha_moment(Z, p, X_used, method),
        psi_Delta1 = psi_Delta1,
        psi_Gamma1 = psi_Gamma1,
        psi_Delta0 = psi_Delta0,
        psi_Gamma0 = psi_Gamma0,
        psi_tau    = Delta1 / Gamma1 - Delta0 / Gamma0 - tau
      )
    }
    
    return(sandwich_se_mest(moment_fn, theta_hat, tau_index))
  }
  
  stop("Unknown estimator")
}

safe_kappa_se <- function(Y, Z, D, X_mat, estimator, method) {
  tryCatch(
    kappa_analytic_se_one(Y, Z, D, X_mat, estimator = estimator, method = method),
    error = function(e) {
      warning(sprintf(
        "SE failed for estimator=%s, method=%s. Returning NA. Original error: %s",
        estimator, method, e$message
      ))
      NA_real_
    }
  )
}


kappa_analytic_se_all <- function(Y, Z, D, X_mat) {
  Y <- as.numeric(Y)
  Z <- as.numeric(Z)
  D <- as.numeric(D)
  X_mat <- as.matrix(X_mat)
  
  # Propensity scores for point estimates
  p_ml <- logit_mle(Z, X_mat)
  p_cb <- cbps(Z, X_mat)$p
  
  # Point estimates
  estimates <- c(
    tau_cb_u   = tau_u(Y, Z, D, p_cb),
    tau_ml_u   = tau_u(Y, Z, D, p_ml),
    tau_ml_a10 = tau_a10(Y, Z, D, p_ml),
    tau_ml_a   = tau_unnorm(Y, Z, D, p_ml, "a"),
    tau_ml_t   = tau_unnorm(Y, Z, D, p_ml, "a1"),
    tau_ml_a0  = tau_unnorm(Y, Z, D, p_ml, "a0")
  )
  
  # Analytical M-estimation SEs
  # safe_kappa_se() prevents the table from stopping if one matrix is singular.
  se <- c(
    tau_cb_u   = safe_kappa_se(Y, Z, D, X_mat, estimator = "u",   method = "cb"),
    tau_ml_u   = safe_kappa_se(Y, Z, D, X_mat, estimator = "u",   method = "ml"),
    tau_ml_a10 = safe_kappa_se(Y, Z, D, X_mat, estimator = "a10", method = "ml"),
    tau_ml_a   = safe_kappa_se(Y, Z, D, X_mat, estimator = "a",   method = "ml"),
    tau_ml_t   = safe_kappa_se(Y, Z, D, X_mat, estimator = "a1",  method = "ml"),
    tau_ml_a0  = safe_kappa_se(Y, Z, D, X_mat, estimator = "a0",  method = "ml")
  )
  
  list(
    estimates = estimates,
    se = se
  )
}
# ==============================================================================
# 7.  FORMATTING HELPER
# ==============================================================================
# Returns "coef*** \n (se)" with significance stars.
# Stars based on two-sided z-test: *** p<0.01, ** p<0.05, * p<0.10.

fmt <- function(coef, se, digits = 3) {
  if (is.na(coef)) return("NA")
  if (is.na(se) || !is.finite(se) || se <= 0) {
    return(sprintf(paste0("%.", digits, "f\n(NA)"), round(coef, digits)))
  }
  
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
#   tau_ml_t  — Tan/Frolich estimator, equal to tau_ml_a1 [unnormalized]
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
    
    # All kappa point estimates and analytical M-estimation SEs
    # Kappa SEs follow the online appendix; 2SLS SEs are HC1 robust.
    X_mat <- spec$X
    Z <- Z_sipp
    D <- D_sipp
    kappa_res <- kappa_analytic_se_all(Y, Z, D, X_mat)
    
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
    cat(sprintf("    tau_ml_t   (MLE, = tau_a1):  %s\n",
                fmt(kappa_res$estimates["tau_ml_t"],  kappa_res$se["tau_ml_t"])))
    cat(sprintf("    tau_ml_a0  (MLE, unnorm.):   %s  [NOT transl. invariant]\n",
                fmt(kappa_res$estimates["tau_ml_a0"], kappa_res$se["tau_ml_a0"])))
    cat("\n")
  }
  cat("\n")
}


# ==============================================================================
# 10. TRANSLATION INVARIANCE CHECK
# ==============================================================================
# Switching from log-dollars to log-cents adds the constant log(100) to every
# outcome observation.
#
# For any estimator that can be written as a weighted outcome sum,
#
#   tau_hat(Y + k) = tau_hat(Y) + k * sum_i omega_i.
#
# Therefore:
#   sum_i omega_i = 0  -> translation invariant
#   sum_i omega_i != 0 -> not translation invariant
#
# The paper proves that tau_u and tau_a10 are translation invariant, while
# tau_a, tau_t = tau_a1, and tau_a0 generally are not. The calculations below
# verify this algebraically for the cubic-age specification.
# ==============================================================================

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
