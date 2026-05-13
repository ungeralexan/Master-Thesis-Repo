# ==============================================================================
# ANGRIST (1990) — MILITARY SERVICE: DoubleML EXTENSION
# Builds on angrist1990_military_kappa.R and angrist1990_dml.R
# Run those first — this script assumes their objects are in memory.
# ==============================================================================
#
# CRITICAL: INSTALL THE GITHUB DEV VERSION OF OutcomeWeights
# -----------------------------------------------------------
# get_outcome_weights() support for DoubleML objects is NOT in the CRAN
# version of OutcomeWeights. It is only in the development version on GitHub.
# Your supervisor explicitly mentioned this ("die Development Version vom
# OutcomeWeights package auf GitHub ist jetzt auch kompatibel mit DoubleML").
#
# You MUST run the installation in BLOCK A before anything else works.
# The CRAN version gives: "no applicable method for 'get_outcome_weights'
# applied to an object of class c('DoubleMLIIVM', 'DoubleML', 'R6')"
#
# STRUCTURE:
#   BLOCK A — Install dev version + setup
#   BLOCK B — DoubleMLIIVM: linear + logistic (parametric baseline)
#   BLOCK C — DoubleMLIIVM: ranger random forest
#   BLOCK D — DoubleMLIIVM: XGBoost
#   BLOCK E — Comparison table: all estimators
#   BLOCK F — Love plots: learner comparison
# ==============================================================================


# ==============================================================================
# BLOCK A — SETUP
# MUST be run first. Installs dev version of OutcomeWeights from GitHub.
# ==============================================================================

# Step 1: install remotes if needed (for GitHub installation)
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

# Step 2: install the GitHub DEV version of OutcomeWeights
# This overwrites the CRAN version and adds DoubleML compatibility
# Only needs to be done once — comment out after first successful install
remotes::install_github("MCKnaus/OutcomeWeights", force = FALSE)

# Step 3: other dependencies
if (!require("DoubleML",  quietly = TRUE)) install.packages("DoubleML")
if (!require("mlr3verse", quietly = TRUE)) install.packages("mlr3verse")
if (!require("xgboost",   quietly = TRUE)) install.packages("xgboost")
if (!require("ranger",    quietly = TRUE)) install.packages("ranger")

# Step 4: load libraries (OutcomeWeights must be loaded AFTER dev install)
library(OutcomeWeights)   # dev version — now has DoubleML support
library(DoubleML)
library(mlr3verse)
library(cobalt)
library(viridis)
library(gridExtra)

# Verify dev version is loaded (should show > 0.1.0 or similar)
cat("OutcomeWeights version:", as.character(packageVersion("OutcomeWeights")), "\n")
cat("Check get_outcome_weights methods available:\n")
print(methods("get_outcome_weights"))
cat("\n")

# Check prerequisites from kappa script
stopifnot(
  exists("sipp_clean"), exists("D_sipp"), exists("Z_sipp"),
  exists("kappa_weights"), exists("tau_u"), exists("tau_a10"),
  exists("tau_unnorm"), exists("logit_mle"),
  exists("kappa_outcome_weights"), exists("weight_diag")
)

D     <- D_sipp
Z     <- Z_sipp
Y_dol <- sipp_clean$lwage
Y_cnt <- sipp_clean$lwage_cnt
k     <- log(100)

cat(sprintf("N=%d | mean(D)=%.3f | mean(Z)=%.3f\n\n",
            length(D), mean(D), mean(Z)))

# ------------------------------------------------------------------------------
# DoubleML data object (cubic covariate spec)
# d and z must be numeric 0/1 (not haven_labelled or factor)
# ------------------------------------------------------------------------------
df_dml <- data.frame(
  y    = as.numeric(sipp_clean$lwage),
  d    = as.numeric(sipp_clean$nvstat),
  z    = as.numeric(sipp_clean$rsncode),
  age  = as.numeric(sipp_clean$age),
  age2 = as.numeric(sipp_clean$age)^2,
  age3 = as.numeric(sipp_clean$age)^3
)

dml_data <- DoubleMLData$new(
  data   = df_dml,
  y_col  = "y",
  d_cols = "d",
  z_cols = "z",
  x_cols = c("age", "age2", "age3")
)

cat("DoubleML data object:\n")
print(dml_data)
cat("\n")

# Kappa matrices and weights for comparison
X_kappa_cub <- model.matrix(~ age + age2 + age3, data = sipp_clean)
X_dml_cub   <- model.matrix(~ 0 + age + I(age^2) + I(age^3), data = sipp_clean)
p_ml_cub    <- logit_mle(Z, X_kappa_cub)
kw_cub      <- kappa_outcome_weights(Y_dol, Z, D, p_ml_cub)

cat("BLOCK A complete. Ready for BLOCKS B, C, D.\n\n")


# ==============================================================================
# BLOCK B — LINEAR + LOGISTIC (parametric baseline)
# ==============================================================================
#
# ml_g (outcome regression E[Y|X,D,Z]):  regr.lm       — regression
# ml_m (instrument propensity P(Z=1|X)): classif.log_reg — classification
# ml_r (treatment propensity P(D=1|X,Z)):classif.log_reg — classification
#
# This is the parametric baseline. With linear/logistic nuisance functions,
# DoubleML essentially recovers a 2SLS-type estimate.

set.seed(42)
cat(strrep("=", 70), "\n")
cat("BLOCK B — DoubleMLIIVM: linear + logistic (parametric baseline)\n")
cat(strrep("=", 70), "\n\n")

iivm_lm <- DoubleMLIIVM$new(
  data    = dml_data,
  ml_g    = lrn("regr.lm"),
  ml_m    = lrn("classif.log_reg"),
  ml_r    = lrn("classif.log_reg"),
  n_folds = 5,
  score   = "LATE"
)
iivm_lm$fit(store_models = TRUE, store_predictions = TRUE)

cat("Point estimate:\n")
print(iivm_lm$summary())

omega_lm   <- get_outcome_weights(object = iivm_lm, dml_data = dml_data)
w_waipw_lm <- as.vector(omega_lm$omega)

ok_lm <- isTRUE(all.equal(
  as.numeric(omega_lm$omega %*% Y_dol),
  as.numeric(iivm_lm$all_coef)
))
cat(sprintf("Algebraic check (omega'Y = tau_hat): %s\n\n", ok_lm))
print(weight_diag(w_waipw_lm, "Wald-AIPW (DoubleML, linear+logistic)"), row.names = FALSE)
cat("\nBLOCK B done. Objects: iivm_lm, w_waipw_lm\n\n")


# ==============================================================================
# BLOCK C — RANGER RANDOM FOREST
# ==============================================================================
#
# ml_g: regr.ranger    (keep.inbag = TRUE — required for outcome weights)
# ml_m: classif.ranger (keep.inbag = TRUE)
# ml_r: classif.ranger (keep.inbag = TRUE)
#
# keep.inbag = TRUE is mandatory: get_outcome_weights() needs the in-bag
# sample indices to reconstruct the smoother matrix from the forest.

set.seed(42)
cat(strrep("=", 70), "\n")
cat("BLOCK C — DoubleMLIIVM: ranger random forest\n")
cat(strrep("=", 70), "\n\n")

iivm_rf <- DoubleMLIIVM$new(
  data    = dml_data,
  ml_g    = lrn("regr.ranger",    keep.inbag = TRUE),
  ml_m    = lrn("classif.ranger", keep.inbag = TRUE),
  ml_r    = lrn("classif.ranger", keep.inbag = TRUE),
  n_folds = 5,
  score   = "LATE"
)
iivm_rf$fit(store_models = TRUE, store_predictions = TRUE)

cat("Point estimate:\n")
print(iivm_rf$summary())

omega_rf   <- get_outcome_weights(object = iivm_rf, dml_data = dml_data)
w_waipw_rf <- as.vector(omega_rf$omega)

ok_rf <- isTRUE(all.equal(
  as.numeric(omega_rf$omega %*% Y_dol),
  as.numeric(iivm_rf$all_coef)
))
cat(sprintf("Algebraic check (omega'Y = tau_hat): %s\n\n", ok_rf))
print(weight_diag(w_waipw_rf, "Wald-AIPW (DoubleML, ranger)"), row.names = FALSE)
cat("\nBLOCK C done. Objects: iivm_rf, w_waipw_rf\n\n")


# ==============================================================================
# BLOCK D — XGBOOST
# ==============================================================================
#
# ml_g: regr.xgboost    with alpha=0, subsample=1, max_delta_step=0, base_score=0
# ml_m: classif.xgboost with alpha=0, subsample=1, max_delta_step=0
# ml_r: classif.xgboost with alpha=0, subsample=1, max_delta_step=0
#
# The four hyperparameters on ml_g are required for omega'Y = tau_hat to hold
# (documented in DML_smoothers_2.rmd by Knaus & Rakov).
# base_score only applies to regr.xgboost (regression learner), not classif.

set.seed(42)
cat(strrep("=", 70), "\n")
cat("BLOCK D — DoubleMLIIVM: XGBoost\n")
cat(strrep("=", 70), "\n\n")

iivm_xgb <- DoubleMLIIVM$new(
  data    = dml_data,
  ml_g    = lrn("regr.xgboost",
                alpha = 0, subsample = 1, max_delta_step = 0, base_score = 0),
  ml_m    = lrn("classif.xgboost",
                alpha = 0, subsample = 1, max_delta_step = 0),
  ml_r    = lrn("classif.xgboost",
                alpha = 0, subsample = 1, max_delta_step = 0),
  n_folds = 5,
  score   = "LATE"
)
iivm_xgb$fit(store_models = TRUE, store_predictions = TRUE)

cat("Point estimate:\n")
print(iivm_xgb$summary())

omega_xgb   <- get_outcome_weights(object = iivm_xgb, dml_data = dml_data)
w_waipw_xgb <- as.vector(omega_xgb$omega)

ok_xgb <- isTRUE(all.equal(
  as.numeric(omega_xgb$omega %*% Y_dol),
  as.numeric(iivm_xgb$all_coef)
))
cat(sprintf("Algebraic check (omega'Y = tau_hat): %s\n\n", ok_xgb))
print(weight_diag(w_waipw_xgb, "Wald-AIPW (DoubleML, XGBoost)"), row.names = FALSE)
cat("\nBLOCK D done. Objects: iivm_xgb, w_waipw_xgb\n\n")


# ==============================================================================
# BLOCK E — FULL COMPARISON TABLE (cubic age spec)
# Requires: BLOCKS B, C, D + angrist1990_dml.R BLOCK B (w_waldaipw_cub)
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("BLOCK E — COMPARISON: all estimators, cubic age, log wages (dollars)\n")
cat(strrep("=", 70), "\n\n")

# Point estimates
cat(sprintf("  %-44s  %8s\n", "Estimator", "Estimate"))
cat(strrep("-", 56), "\n")
cat(sprintf("  %-44s  %8.4f\n", "tau_u      (kappa, MLE, normalized)",
            tau_u(Y_dol, Z, D, p_ml_cub)))
cat(sprintf("  %-44s  %8.4f\n", "tau_a10    (kappa, MLE, normalized)",
            tau_a10(Y_dol, Z, D, p_ml_cub)))
if (exists("res_cub"))
  cat(sprintf("  %-44s  %8.4f\n", "Wald-AIPW  (grf, dml_with_smoother)",
              res_cub["Wald-AIPW", 1]))
cat(sprintf("  %-44s  %8.4f\n", "Wald-AIPW  (DoubleML, linear+logistic)",
            as.numeric(iivm_lm$coef)))
cat(sprintf("  %-44s  %8.4f\n", "Wald-AIPW  (DoubleML, ranger)",
            as.numeric(iivm_rf$coef)))
cat(sprintf("  %-44s  %8.4f\n", "Wald-AIPW  (DoubleML, XGBoost)",
            as.numeric(iivm_xgb$coef)))
cat("\n")

# Weight diagnostics
cat("Weight diagnostics (Sum_w = 0 => translation invariant):\n\n")
diag_rows <- rbind(
  weight_diag(kw_cub$w_u,    "tau_u          (kappa, norm.)"),
  weight_diag(kw_cub$w_a10,  "tau_a10        (kappa, norm.)"),
  weight_diag(kw_cub$w_a,    "tau_a          (kappa, unnorm.)"),
  weight_diag(kw_cub$w_a0,   "tau_a0         (kappa, unnorm.)")
)
if (exists("w_waldaipw_cub"))
  diag_rows <- rbind(diag_rows,
                     weight_diag(w_waldaipw_cub, "Wald-AIPW      (grf)"))
diag_rows <- rbind(diag_rows,
                   weight_diag(w_waipw_lm,  "Wald-AIPW      (DoubleML, linear+logistic)"),
                   weight_diag(w_waipw_rf,  "Wald-AIPW      (DoubleML, ranger)"),
                   weight_diag(w_waipw_xgb, "Wald-AIPW      (DoubleML, XGBoost)")
)
print(diag_rows, row.names = FALSE)

cat("\nKEY QUESTIONS FOR THESIS SECTION 4.1 / 5.2:\n")
cat("  1. Convergence of point estimates across learners? (robustness)\n")
cat("  2. ESS comparison: kappa vs all Wald-AIPW variants\n")
cat("  3. % negative weights: does ML learner choice affect extrapolation?\n")
cat("  4. Sum_w ~= 0 for all Wald-AIPW? (translation invariance holds)\n\n")


# ==============================================================================
# BLOCK F — LOVE PLOTS: LEARNER COMPARISON
# Requires: BLOCKS B, C, D.
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("BLOCK F — LOVE PLOTS: learner comparison (cubic spec)\n")
cat(strrep("=", 70), "\n\n")

make_love_dml <- function(title_str, w_vec) {
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

lp_lm  <- make_love_dml("Wald-AIPW (linear + logistic)", w_waipw_lm)
lp_rf  <- make_love_dml("Wald-AIPW (ranger)",             w_waipw_rf)
lp_xgb <- make_love_dml("Wald-AIPW (XGBoost)",            w_waipw_xgb)
lp_ku  <- make_love_dml("tau_u (kappa, MLE)",              kw_cub$w_u)

grid.arrange(lp_lm, lp_rf, lp_xgb, lp_ku, nrow = 2)

cat(strrep("=", 70), "\n")
cat("DoubleML extension complete.\n")
cat(strrep("=", 70), "\n")