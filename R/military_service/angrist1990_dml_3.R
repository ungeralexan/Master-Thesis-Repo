# ==============================================================================
# ANGRIST (1990) — MILITARY SERVICE: DML EXTENSION
# Builds on angrist1990_military_kappa.R — run that first.
# ==============================================================================
#
# WHY ONLY TWO SPECS (CUBIC AND SATURATED)?
# ------------------------------------------
# The linear age spec (single column X) consistently fails grf's internal
# validate_X() in R 4.5.x: a 1-column matrix triggers an edge case in the
# random forest implementation where the splitting algorithm has no useful
# variation to work with.  Beyond the technical issue, a single linear
# age term is also the least interesting spec for a nonparametric method
# like DML — the kernel smoother adds no value over 2SLS when X has only
# one column.  The cubic and saturated specs are the substantively
# relevant ones for DML and match Table 2 columns 3-6 of SUW (2025).
#
# STRUCTURE (run BLOCKS in order; B and C are independent of each other):
#   BLOCK A — Setup, helpers, X matrices
#   BLOCK B — Spec 2: Cubic in age   (main spec)
#   BLOCK C — Spec 3: Saturated age  (most flexible)
#   BLOCK D — Translation invariance demo (both specs)
#   BLOCK E — Combined weight diagnostics (cubic spec)
#   BLOCK F — Love plots (cubic spec)
#   BLOCK G — Summary point-estimates table
# ==============================================================================


# ==============================================================================
# BLOCK A — SETUP
# Run first. Requires angrist1990_military_kappa.R to be in memory.
# ==============================================================================

library(OutcomeWeights)  # dml_with_smoother(), get_outcome_weights(), plot()
library(cobalt)          # love.plot()
library(viridis)         # colour palette
library(gridExtra)       # grid.arrange()
library(ggplot2)         # required by OutcomeWeights

stopifnot(
  exists("sipp_clean"), exists("D_sipp"), exists("Z_sipp"),
  exists("kappa_weights"), exists("tau_u"), exists("tau_a10"),
  exists("tau_unnorm"), exists("logit_mle")
)

D     <- D_sipp
Z     <- Z_sipp
Y_dol <- sipp_clean$lwage        # log wages in dollars
Y_cnt <- sipp_clean$lwage_cnt    # log wages in cents = Y_dol + log(100)
k     <- log(100)                 # translation constant

cat(sprintf("BLOCK A ready. N=%d | mean(D)=%.3f | mean(Z)=%.3f\n\n",
            length(D), mean(D), mean(Z)))

# ------------------------------------------------------------------------------
# make_X(): strips all attributes from model.matrix to produce a clean
# numeric matrix that passes grf's validate_X() checks.
# Requires ncol >= 2 — single-column matrices fail grf in R 4.5.x.
# ------------------------------------------------------------------------------
make_X <- function(formula, data) {
  mm  <- model.matrix(formula, data = data)
  m   <- matrix(as.numeric(mm), nrow = nrow(mm), ncol = ncol(mm))
  colnames(m) <- colnames(mm)
  m
}

sipp_clean$age_fac <- factor(sipp_clean$age_5)

# Spec 2: Cubic in age — 3 columns, works fine
X_dml_cub <- make_X(~ 0 + age + I(age^2) + I(age^3), sipp_clean)

# Spec 3: Saturated age — 10 columns (one per age level minus reference)
X_dml_sat_full <- make_X(~ 0 + age_fac, sipp_clean)
X_dml_sat      <- X_dml_sat_full[, -1, drop = FALSE]

cat("DML X matrices:\n")
cat(sprintf("  X_dml_cub: %s | %dx%d | numeric=%s\n",
            class(X_dml_cub)[1], nrow(X_dml_cub), ncol(X_dml_cub),
            is.numeric(X_dml_cub)))
cat(sprintf("  X_dml_sat: %s | %dx%d | numeric=%s\n",
            class(X_dml_sat)[1], nrow(X_dml_sat), ncol(X_dml_sat),
            is.numeric(X_dml_sat)))

# Kappa X matrices — WITH intercept (for logit_mle)
X_kappa_lin <- model.matrix(~ age,               data = sipp_clean)
X_kappa_cub <- model.matrix(~ age + age2 + age3, data = sipp_clean)
X_kappa_sat <- model.matrix(~ age_fac,           data = sipp_clean)

# ------------------------------------------------------------------------------
# Helper: per-observation kappa outcome weights
# ------------------------------------------------------------------------------
kappa_outcome_weights <- function(Y, Z, D, p) {
  n   <- length(Y)
  kw  <- kappa_weights(Z, D, p)
  s1  <- sum(Z / p);  s0 <- sum((1 - Z) / (1 - p))
  dD  <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  w_u <- (Z / p / s1 - (1 - Z) / (1 - p) / s0) / dD
  w_a10 <- kw$kappa1 / sum(kw$kappa1) - kw$kappa0 / sum(kw$kappa0)
  num_w <- (Z - p) / (p * (1 - p)) / n
  list(w_u   = as.vector(w_u),
       w_a10 = as.vector(w_a10),
       w_a   = as.vector(num_w / mean(kw$kappa)),
       w_a1  = as.vector(num_w / mean(kw$kappa1)),
       w_a0  = as.vector(num_w / mean(kw$kappa0)))
}

# ------------------------------------------------------------------------------
# Helper: weight diagnostics row
# Sum_w = 0 => translation invariant | ESS = 1/sum(w^2)
# ------------------------------------------------------------------------------
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

cat("\nBLOCK A complete. Ready for BLOCKS B and C.\n\n")


# ==============================================================================
# BLOCK B — SPEC 2: CUBIC IN AGE  (main spec)
# X = [age, age^2, age^3]   (Table 2, columns 3-4)
# ==============================================================================
#
# Three age terms give the random forest enough variation to adapt.
# This is the main spec used for the weight diagnostics and Love plots.

set.seed(42)
cat(strrep("=", 70), "\n")
cat("BLOCK B — DML: Cubic in age specification\n")
cat(strrep("=", 70), "\n\n")

dml_cub <- dml_with_smoother(Y_dol, D, X_dml_cub, Z,
                              n_cf_folds = 5, tune.parameters = "all")
res_cub <- summary(dml_cub)

cat("Point estimates (IV-relevant rows only):\n")
print(res_cub[c("PLR-IV", "Wald-AIPW"), , drop = FALSE])

# omega$omega is (4 x N): rows = PLR / PLR-IV / AIPW-ATE / Wald-AIPW
omega_cub      <- get_outcome_weights(dml_cub)
w_plriv_cub    <- as.vector(omega_cub$omega["PLR-IV",    ])
w_waldaipw_cub <- as.vector(omega_cub$omega["Wald-AIPW", ])

# Algebraic identity check: omega'Y = tau_hat exactly (Knaus 2024)
cat("\nAlgebraic check (omega'Y = tau_hat, must be TRUE):\n")
cat("  PLR-IV:   ",
    isTRUE(all.equal(sum(w_plriv_cub    * Y_dol), res_cub["PLR-IV",    1])), "\n")
cat("  Wald-AIPW:",
    isTRUE(all.equal(sum(w_waldaipw_cub * Y_dol), res_cub["Wald-AIPW", 1])), "\n\n")

cat("Weight diagnostics (Sum_w = 0 => translation invariant):\n")
print(rbind(
  weight_diag(w_plriv_cub,    "PLR-IV     (DML, cubic)"),
  weight_diag(w_waldaipw_cub, "Wald-AIPW  (DML, cubic)")
), row.names = FALSE)

cat("\nKnaus weight distribution plot (all 4 DML estimators):\n")
plot(dml_cub)

cat("\nBLOCK B done. Objects: dml_cub, w_plriv_cub, w_waldaipw_cub\n\n")


# ==============================================================================
# BLOCK C — SPEC 3: SATURATED AGE  (most flexible)
# X = one dummy per age value   (Table 2, columns 5-6)
# ==============================================================================
#
# Fully nonparametric in age. Equivalent to within-age-cell Wald for kappa.
# The random forest has many binary features — splits cleanly.

set.seed(42)
cat(strrep("=", 70), "\n")
cat("BLOCK C — DML: Saturated age specification\n")
cat(strrep("=", 70), "\n\n")

dml_sat <- dml_with_smoother(Y_dol, D, X_dml_sat, Z,
                              n_cf_folds = 5, tune.parameters = "all")
res_sat <- summary(dml_sat)

cat("Point estimates (IV-relevant rows only):\n")
print(res_sat[c("PLR-IV", "Wald-AIPW"), , drop = FALSE])

omega_sat      <- get_outcome_weights(dml_sat)
w_plriv_sat    <- as.vector(omega_sat$omega["PLR-IV",    ])
w_waldaipw_sat <- as.vector(omega_sat$omega["Wald-AIPW", ])

cat("\nAlgebraic check:\n")
cat("  PLR-IV:   ",
    isTRUE(all.equal(sum(w_plriv_sat    * Y_dol), res_sat["PLR-IV",    1])), "\n")
cat("  Wald-AIPW:",
    isTRUE(all.equal(sum(w_waldaipw_sat * Y_dol), res_sat["Wald-AIPW", 1])), "\n\n")

cat("Weight diagnostics:\n")
print(rbind(
  weight_diag(w_plriv_sat,    "PLR-IV     (DML, saturated)"),
  weight_diag(w_waldaipw_sat, "Wald-AIPW  (DML, saturated)")
), row.names = FALSE)

cat("\nKnaus weight distribution plot:\n")
plot(dml_sat)

cat("\nBLOCK C done. Objects: dml_sat, w_plriv_sat, w_waldaipw_sat\n\n")


# ==============================================================================
# BLOCK D — TRANSLATION INVARIANCE & SCALE EQUIVARIANCE
# Requires: BLOCKS B and C.
# ==============================================================================
#
# For tau = sum(w_i * Y_i):
#   tau(Y + k) = tau(Y) + k * sum(w_i)   where k = log(100)
#
# sum(w_i) = 0  =>  no shift  [TRANSLATION INVARIANT]
# sum(w_i) != 0 =>  shifts    [NOT invariant]
#
# DML PLR-IV and Wald-AIPW have sum(w) = 0 by construction of the
# residual-on-residual estimator — they are invariant like normalized kappa.
#
# For DML we do NOT re-run the model: apply same omega weights to Y_cnt.
# This demonstrates that invariance is a property of the WEIGHTS, not
# of the estimation procedure.

cat(strrep("=", 70), "\n")
cat("BLOCK D — TRANSLATION INVARIANCE & SCALE EQUIVARIANCE\n")
cat(sprintf("k = log(100) = %.6f\n\n", k))

print_inv_table <- function(spec_label, w_plriv, w_waldaipw, X_kappa) {
  cat(sprintf("--- %s ---\n", spec_label))
  cat(sprintf("  %-28s  %8s  %8s  %10s  %10s  %7s\n",
              "Estimator", "Dollars", "Cents", "Diff(act)", "Diff(pred)", "Match?"))
  cat(strrep("-", 80), "\n")
  p_ml <- logit_mle(Z, X_kappa)
  kw   <- kappa_outcome_weights(Y_dol, Z, D, p_ml)
  rows <- list(
    list(name = "PLR-IV     (DML)",
         w = w_plriv,
         ed = sum(w_plriv    * Y_dol), ec = sum(w_plriv    * Y_cnt)),
    list(name = "Wald-AIPW  (DML)",
         w = w_waldaipw,
         ed = sum(w_waldaipw * Y_dol), ec = sum(w_waldaipw * Y_cnt)),
    list(name = "tau_u      (norm.)",
         w = kw$w_u,
         ed = tau_u(Y_dol,Z,D,p_ml),    ec = tau_u(Y_cnt,Z,D,p_ml)),
    list(name = "tau_a10    (norm.)",
         w = kw$w_a10,
         ed = tau_a10(Y_dol,Z,D,p_ml),  ec = tau_a10(Y_cnt,Z,D,p_ml)),
    list(name = "tau_a      (unnorm.)",
         w = kw$w_a,
         ed = tau_unnorm(Y_dol,Z,D,p_ml,"a"),
         ec = tau_unnorm(Y_cnt,Z,D,p_ml,"a")),
    list(name = "tau_a1     (unnorm.)",
         w = kw$w_a1,
         ed = tau_unnorm(Y_dol,Z,D,p_ml,"a1"),
         ec = tau_unnorm(Y_cnt,Z,D,p_ml,"a1")),
    list(name = "tau_a0     (unnorm.)",
         w = kw$w_a0,
         ed = tau_unnorm(Y_dol,Z,D,p_ml,"a0"),
         ec = tau_unnorm(Y_cnt,Z,D,p_ml,"a0"))
  )
  for (r in rows) {
    act  <- r$ec - r$ed
    pred <- sum(r$w) * k
    ok   <- isTRUE(all.equal(act, pred, tolerance = 1e-8))
    cat(sprintf("  %-28s  %8.4f  %8.4f  %10.4f  %10.4f  %7s\n",
                r$name, r$ed, r$ec, act, pred, if (ok) "TRUE" else "FALSE"))
  }
  cat("\n")
}

print_inv_table("Cubic in age",   w_plriv_cub, w_waldaipw_cub, X_kappa_cub)
print_inv_table("Saturated age",  w_plriv_sat, w_waldaipw_sat, X_kappa_sat)

cat("DML PLR-IV and Wald-AIPW: sum(w) ~= 0  =>  INVARIANT\n")
cat("Kappa tau_a / tau_a1 / tau_a0: sum(w) != 0  =>  NOT invariant\n\n")


# ==============================================================================
# BLOCK E — COMBINED WEIGHT DIAGNOSTICS (cubic spec)
# Requires: BLOCK B.
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("BLOCK E — COMBINED WEIGHT DIAGNOSTICS (cubic spec)\n")
cat(strrep("=", 70), "\n\n")

p_ml_cub <- logit_mle(Z, X_kappa_cub)
kw_cub   <- kappa_outcome_weights(Y_dol, Z, D, p_ml_cub)

cat("Sum_w = 0 => translation invariant | ESS = effective sample size\n\n")
print(rbind(
  weight_diag(w_plriv_cub,    "PLR-IV     (DML)"),
  weight_diag(w_waldaipw_cub, "Wald-AIPW  (DML)"),
  weight_diag(kw_cub$w_u,    "tau_u      (kappa, norm.)"),
  weight_diag(kw_cub$w_a10,  "tau_a10    (kappa, norm.)"),
  weight_diag(kw_cub$w_a,    "tau_a      (kappa, unnorm.)"),
  weight_diag(kw_cub$w_a1,   "tau_a1     (kappa, unnorm.)"),
  weight_diag(kw_cub$w_a0,   "tau_a0     (kappa, unnorm.)")
), row.names = FALSE)
cat("\n")


# ==============================================================================
# BLOCK F — LOVE PLOTS (cubic spec)
# Requires: BLOCK B.
# Knaus convention: multiply outcome weights by (2*D-1) for cobalt sign
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("BLOCK F — LOVE PLOTS (cubic spec)\n")
cat(strrep("=", 70), "\n\n")

make_love <- function(title_str, w_vec) {
  love.plot(
    D ~ X_dml_cub,
    weights    = w_vec * (2 * D - 1),
    position   = "bottom",
    title      = title_str,
    thresholds = c(m = 0.1),
    var.order  = "unadjusted",
    binary     = "std",
    abs        = TRUE,
    line       = TRUE,
    colors     = viridis(2),
    shapes     = c("circle", "triangle")
  )
}

lp1 <- make_love("Angrist (1990): PLR-IV",    w_plriv_cub)
lp2 <- make_love("Angrist (1990): Wald-AIPW", w_waldaipw_cub)
lp3 <- make_love("Angrist (1990): tau_u",     kw_cub$w_u)
lp4 <- make_love("Angrist (1990): tau_a10",   kw_cub$w_a10)
grid.arrange(lp1, lp2, lp3, lp4, nrow = 2)


# ==============================================================================
# BLOCK G — SUMMARY TABLE: DML VS KAPPA (cubic + saturated)
# Requires: BLOCKS B and C.
# ==============================================================================
#
# Note: the linear age kappa estimates are included for completeness
# since they run fine from the kappa script.  DML linear age is skipped
# for the reason documented at the top of this script.

cat(strrep("=", 70), "\n")
cat("BLOCK G — POINT ESTIMATES: DML vs KAPPA, ALL AVAILABLE SPECS\n")
cat("(log wages in dollars)\n")
cat(strrep("=", 70), "\n\n")

p_lin <- logit_mle(Z, X_kappa_lin)
p_cub <- logit_mle(Z, X_kappa_cub)
p_sat <- logit_mle(Z, X_kappa_sat)

cat(sprintf("  %-28s  %10s  %10s  %10s\n",
            "Estimator", "Lin.(kappa)", "Cubic", "Saturated"))
cat(strrep("-", 62), "\n")

# DML: only cubic and saturated available
cat(sprintf("  %-28s  %10s  %10.4f  %10.4f\n", "PLR-IV (DML)",
            "n/a",
            res_cub["PLR-IV",    1],
            res_sat["PLR-IV",    1]))
cat(sprintf("  %-28s  %10s  %10.4f  %10.4f\n", "Wald-AIPW (DML)",
            "n/a",
            res_cub["Wald-AIPW", 1],
            res_sat["Wald-AIPW", 1]))
cat(strrep("-", 62), "\n")

# Kappa: all three specs
krow <- function(lbl, fn) {
  cat(sprintf("  %-28s  %10.4f  %10.4f  %10.4f\n",
              lbl, fn(p_lin), fn(p_cub), fn(p_sat)))
}
krow("tau_u      (norm.)",   function(p) tau_u(Y_dol, Z, D, p))
krow("tau_a10    (norm.)",   function(p) tau_a10(Y_dol, Z, D, p))
krow("tau_a      (unnorm.)", function(p) tau_unnorm(Y_dol, Z, D, p, "a"))
krow("tau_a1     (unnorm.)", function(p) tau_unnorm(Y_dol, Z, D, p, "a1"))
krow("tau_a0     (unnorm.)", function(p) tau_unnorm(Y_dol, Z, D, p, "a0"))

cat(strrep("=", 70), "\n")
cat("DML extension complete.\n")
cat(strrep("=", 70), "\n")
