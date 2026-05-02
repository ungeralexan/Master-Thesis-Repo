# ==============================================================================
# 00_functions.R
# Core functions for kappa weighting estimators of the LATE
# Based on: Słoczyński, Uysal & Wooldridge (2025, JBES)
#           Abadie (2003, JoE)
# ==============================================================================

# ---- 0. Data path ------------------------------------------------------------
# Set this to your local data directory
DATA_PATH <- "/Users/alexung/desktop/suw_jbes_replicate"

# ---- 1. Packages -------------------------------------------------------------
library(haven)      # read_dta
library(AER)        # ivreg (2SLS)
library(sandwich)   # vcovHC (robust SE)
library(lmtest)     # coeftest
library(boot)       # bootstrap
set.seed(42)

# ---- 2. Propensity score estimation ------------------------------------------

#' Estimate instrument propensity score p(X) = P(Z=1|X) via logit (ML)
#'
#' @param Z  Binary instrument vector (0/1)
#' @param X  Design matrix including intercept
#' @return   Fitted probabilities p̂(Xᵢ) ∈ (0,1)
fit_propensity_ml <- function(Z, X) {
  fit <- glm.fit(X, Z, family = binomial(link = "logit"))
  pmax(pmin(fitted(fit), 1 - 1e-6), 1e-6)  # trim to (0,1)
}

#' Estimate instrument propensity score via covariate balancing (CBPS)
#' Solves: (1/N) Σ [Z/p(X) - (1-Z)/(1-p(X))] X = 0
#' Following Imai & Ratkovic (2014), as described in SUW (2025) eq. (9)
#'
#' @param Z  Binary instrument vector
#' @param X  Design matrix including intercept
#' @return   Fitted probabilities p̂(Xᵢ)
fit_propensity_cb <- function(Z, X) {
  # Initial values from logit
  alpha_init <- glm.fit(X, Z, family = binomial())$coefficients
  sigmoid    <- function(u) 1 / (1 + exp(-u))

  # Moment condition: E[Z/p(X) * X] = E[(1-Z)/(1-p(X)) * X]
  # Equivalent to: E[(Z - p(X)) / (p(X)(1-p(X))) * X] = 0
  obj <- function(alpha) {
    p   <- sigmoid(as.vector(X %*% alpha))
    p   <- pmax(pmin(p, 1 - 1e-6), 1e-6)
    res <- colMeans((Z - p) / (p * (1 - p)) * X)
    sum(res^2)
  }

  opt   <- optim(alpha_init, obj, method = "BFGS",
                 control = list(maxit = 500, reltol = 1e-10))
  p_hat <- sigmoid(as.vector(X %*% opt$par))
  pmax(pmin(p_hat, 1 - 1e-6), 1e-6)
}

# ---- 3. Kappa weights --------------------------------------------------------

#' Compute the three kappa weights from Abadie (2003) / SUW (2025) Lemma 2.1
#'
#' @param Z  Binary instrument
#' @param D  Binary treatment
#' @param p  Propensity scores p(Xᵢ)
#' @return   List with κ, κ₁, κ₀ vectors
kappa_weights <- function(Z, D, p) {
  kappa1 <- D       * (Z - p)       / (p * (1 - p))
  kappa0 <- (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  kappa  <- 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p
  list(kappa = kappa, kappa1 = kappa1, kappa0 = kappa0)
}

# ---- 4. Point estimators -----------------------------------------------------

#' τ̂ᵤ — Uysal (2011) normalized estimator [RECOMMENDED]
#' Translation invariant + scale equivariant + positive denominator (one-sided)
#' SUW (2025) equation (3)
tau_u <- function(Y, Z, D, p) {
  s1   <- sum(Z / p);       s0   <- sum((1 - Z) / (1 - p))
  num  <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denom <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  num / denom
}

#' τ̂ₐ,₁₀ — Abadie-Cattaneo (2018) normalized estimator
#' Translation invariant + scale equivariant, but may have near-zero denom
tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}

#' τ̂ₐ, τ̂ₐ,₁, τ̂ₐ,₀ — unnormalized estimators [NOT RECOMMENDED]
#' Fail translation invariance and scale equivariance
#' SUW (2025) equations (4)–(6)
tau_unnorm <- function(Y, Z, D, p, type = c("a", "a1", "a0")) {
  type  <- match.arg(type)
  kw    <- kappa_weights(Z, D, p)
  n     <- length(Y)
  num   <- mean(Y * (Z - p) / (p * (1 - p)))

  denom <- switch(type,
    "a"  = mean(kw$kappa),
    "a1" = mean(kw$kappa1),   # = tau_t (Tan 2006 / Frölich 2007)
    "a0" = mean(kw$kappa0)
  )
  num / denom
}

# ---- 5. All estimators together ----------------------------------------------

#' Compute all five kappa estimators for a given (Y, Z, D, X)
#'
#' @param Y  Outcome vector
#' @param Z  Binary instrument
#' @param D  Binary treatment
#' @param X  Design matrix (with intercept)
#' @param method  "ml" (logit MLE) or "cb" (covariate balancing)
#' @return  Named vector of point estimates
all_kappa_estimators <- function(Y, Z, D, X, method = c("both", "ml", "cb")) {
  method <- match.arg(method)
  results <- c()

  if (method %in% c("both", "ml")) {
    p_ml <- fit_propensity_ml(Z, X)
    results <- c(results,
      tau_ml_u    = tau_u(Y, Z, D, p_ml),
      tau_ml_a10  = tau_a10(Y, Z, D, p_ml),
      tau_ml_a    = tau_unnorm(Y, Z, D, p_ml, "a"),
      tau_ml_t    = tau_unnorm(Y, Z, D, p_ml, "a1"),
      tau_ml_a0   = tau_unnorm(Y, Z, D, p_ml, "a0")
    )
  }

  if (method %in% c("both", "cb")) {
    p_cb <- fit_propensity_cb(Z, X)
    results <- c(results,
      tau_cb_u    = tau_u(Y, Z, D, p_cb),
      tau_cb_a10  = tau_a10(Y, Z, D, p_cb)
    )
  }

  results
}

# ---- 6. Bootstrap standard errors -------------------------------------------

#' Bootstrap SE for all kappa estimators
#'
#' @param Y, Z, D, X  Data vectors/matrix
#' @param R  Number of bootstrap replications (200 for speed, 500 for final)
#' @return  List with estimates and bootstrap SEs
boot_se_all <- function(Y, Z, D, X, R = 200) {
  ests <- all_kappa_estimators(Y, Z, D, X, method = "both")

  boot_fn <- function(data, idx) {
    Yi <- data$Y[idx]; Zi <- data$Z[idx]
    Di <- data$D[idx]; Xi <- data$X[idx, , drop = FALSE]
    all_kappa_estimators(Yi, Zi, Di, Xi, method = "both")
  }

  dat  <- list(Y = Y, Z = Z, D = D, X = X)
  res  <- boot::boot(dat, boot_fn, R = R)
  se   <- apply(res$t, 2, sd, na.rm = TRUE)
  names(se) <- names(ests)

  list(estimates = ests, se = se)
}

# ---- 7. 2SLS helper ----------------------------------------------------------

#' Run 2SLS and return coefficient + robust SE for the endogenous variable
#'
#' @param Y  Outcome
#' @param D  Endogenous treatment (binary)
#' @param Z  Instrument (binary)
#' @param X_df  Data frame of covariates (no intercept)
#' @param endog_name  Name of D in formula (default "D")
run_2sls <- function(Y, D, Z, X_df, endog_name = "D") {
  df         <- data.frame(Y = Y, D = D, Z = Z, X_df)
  cov_names  <- names(X_df)
  rhs_endo   <- paste(endog_name, paste(cov_names, collapse = " + "), sep = " + ")
  rhs_inst   <- paste("Z", paste(cov_names, collapse = " + "), sep = " + ")
  fml        <- as.formula(paste("Y ~", rhs_endo, "| ", rhs_inst))
  fit        <- AER::ivreg(fml, data = df)
  se         <- sqrt(diag(sandwich::vcovHC(fit, type = "HC1")))[endog_name]
  list(coef = coef(fit)[endog_name], se = se)
}

# ---- 8. Formatting helpers ---------------------------------------------------

#' Format coefficient with SE and significance stars: "0.234**\n(0.123)"
fmt <- function(coef, se, digits = 3) {
  if (is.na(coef) || is.na(se)) return("NA")
  pval  <- 2 * pnorm(-abs(coef / se))
  stars <- ifelse(pval < 0.01, "***", ifelse(pval < 0.05, "**",
           ifelse(pval < 0.10, "*", "")))
  sprintf(paste0("%.", digits, "f%s\n(%.", digits, "f)"),
          round(coef, digits), stars, round(se, digits))
}

cat("00_functions.R loaded successfully.\n")
