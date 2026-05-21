# before double ml
library(haven)
library(AER)
library(sandwich)
library(lmtest)

data_path <- "/Users/alexung/Desktop/master_thesis/code/"          # <--- USER: set this path
card <- read_dta(file.path(data_path, "card.dta"))

# check whether the sample size sums up to 3010
card_clean <- card[!is.na(card$wage) & !is.na(card$educ), ]
cat(sprintf("Sample size: N = %d\n", nrow(card_clean)))


# define log wage in cents and in dollars 
card_clean$lwage_cnt <- log(card_clean$wage)           # cents (original units)
card_clean$lwage     <- log(card_clean$wage / 100)     # dollars

cat(sprintf("log(100) = %.6f\n", log(100)))
cat(sprintf("Mean diff (should equal log(100)): %.6f\n",
            mean(card_clean$lwage_cnt - card_clean$lwage)))
card_clean$somecol <- as.integer(card_clean$educ > 12)
card_clean$educ16  <- as.integer(card_clean$educ >= 16)
card_clean$Z_inst <- as.integer(card_clean$nearc4)

# --- Experience squared (may not exist in some versions of the dataset) ---
if (!"expersq" %in% names(card_clean)) {
  card_clean$expersq <- card_clean$exper^2
}

# --- Region dummies: reg661–reg668 ---
# Card (1995) uses 9 Census division indicators, with one omitted as reference.
# In the NLSYM dataset these are often named reg661, reg662, ..., reg668 (8 dummies,
# leaving region 9 as reference). Check column names and adapt if needed.
region_vars <- grep("^reg66", names(card_clean), value = TRUE)



# --- Card specification: design matrices ---

# For the kappa propensity score (logit / CBPS): intercept + covariates
X_card_kappa <- model.matrix(
  ~ exper + expersq + black + south + smsa + smsa66 +
    reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668,
  data = card_clean
)


# For 2SLS: same covariates as a data.frame, WITHOUT the intercept column
X_card_df <- as.data.frame(X_card_kappa[, -1, drop = FALSE])
# The intercept column is removed because ivreg() adds its own intercept.

# --- Kitagawa specification: design matrices ---

X_kit_kappa <- model.matrix(
  ~ black + south + smsa + smsa66 + south66,
  data = card_clean
)

X_kit_df <- as.data.frame(X_kit_kappa[, -1, drop = FALSE])

# convenience short hand notations 
Y_cnt <- card_clean$lwage_cnt   # log wages, cents before log
Y_dol <- card_clean$lwage       # log wages, dollars before log
Z     <- card_clean$Z_inst      # instrument: nearc4
D1    <- card_clean$somecol     # treatment 1: some college
D2    <- card_clean$educ16      # treatment 2: college completion


# for the next step alreading reading in all functions I need 
source("/Users/alexung/Desktop/master_thesis/showing/functions_card.R")


# this is to run every column
run_column <- function(Y, Z, D, X_kappa, X_df, label) {
  cat(sprintf("  Computing column: %s\n", label))
  
  # 2SLS
  tsls      <- run_2sls(Y, D, Z, X_df, endog_name = "D")
  
  # Kappa estimators (point estimates + M-estimation SEs)
  kappa_res <- kappa_analytic_se_all(Y, Z, D, X_kappa)
  
  list(tsls = tsls, kappa = kappa_res)
}



# compuing propensity score and maximum likelihood
p_ml_card <- logit_mle(Z, X_card_kappa)
p_cb_card <- get_cbps_p(Z, X_card_kappa)
p_ml_kit  <- logit_mle(Z, X_kit_kappa)
p_cb_kit  <- get_cbps_p(Z, X_kit_kappa)




# Col 1: somecol, cents, Card spec
col1 <- run_column(Y_cnt, Z, D1, X_card_kappa, X_card_df, "col1: somecol, cents, Card")

# Col 2: somecol, dollars, Card spec
col2 <- run_column(Y_dol, Z, D1, X_card_kappa, X_card_df, "col2: somecol, dollars, Card")

# Col 3: somecol, cents, Kitagawa spec
col3 <- run_column(Y_cnt, Z, D1, X_kit_kappa, X_kit_df, "col3: somecol, cents, Kitagawa")

# Col 4: somecol, dollars, Kitagawa spec
col4 <- run_column(Y_dol, Z, D1, X_kit_kappa, X_kit_df, "col4: somecol, dollars, Kitagawa")





# Col 5: educ16, cents, Card spec
col5 <- run_column(Y_cnt, Z, D2, X_card_kappa, X_card_df, "col5: educ16, cents, Card")

# Col 6: educ16, dollars, Card spec
col6 <- run_column(Y_dol, Z, D2, X_card_kappa, X_card_df, "col6: educ16, dollars, Card")

# Col 7: educ16, cents, Kitagawa spec
col7 <- run_column(Y_cnt, Z, D2, X_kit_kappa, X_kit_df, "col7: educ16, cents, Kitagawa")

# Col 8: educ16, dollars, Kitagawa spec
col8 <- run_column(Y_dol, Z, D2, X_kit_kappa, X_kit_df, "col8: educ16, dollars, Kitagawa")


cols <- list(col1, col2, col3, col4, col5, col6, col7, col8)




library(OutcomeWeights)  # dml_with_smoother(), get_outcome_weights()
library(cobalt)          # love.plot()
library(viridis)         # colour palette for Love plots
library(gridExtra)       # grid.arrange()
library(ggplot2)         # required by OutcomeWeights


k_shift <- log(100)  # Y_cnt = Y_dol + k_shift (used in Block D)



# ==============================================================================
# BLOCK A: X matrices for DML
# ==============================================================================

make_X <- function(formula, data) {
  # Strips all attributes from model.matrix output to produce a plain numeric
  # matrix safe for grf::validate_X(). This is identical to the Vietnam DML
  # implementation (vietmam_14_05_double_ml.R, Block A).
  mm <- model.matrix(formula, data = data)
  m  <- matrix(as.numeric(mm), nrow = nrow(mm), ncol = ncol(mm))
  colnames(m) <- colnames(mm)
  m
}

# Card spec DML matrix (14 covariates, no intercept)
X_dml_card <- make_X(
  ~ 0 + exper + expersq + black + south + smsa + smsa66 +
    reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668,
  data = card_clean
)

# Kitagawa spec DML matrix (5 covariates, no intercept)
X_dml_kit <- make_X(
  ~ 0 + black + south + smsa + smsa66 + south66,
  data = card_clean
)



get_estimate <- function(res, row_name) {
  # Extract point estimate from OutcomeWeights summary() output.
  # Handles both "Estimate" column name and positional access.
  if ("Estimate" %in% colnames(res)) return(as.numeric(res[row_name, "Estimate"]))
  as.numeric(res[row_name, 1])
}

check_omega_rows <- function(omega_obj, rows = c("PLR-IV", "Wald-AIPW")) {
  # Verify that the expected outcome-weight rows exist in the omega object.
  missing_rows <- setdiff(rows, rownames(omega_obj$omega))
  if (length(missing_rows) > 0)
    stop("Missing expected outcome-weight rows: ", paste(missing_rows, collapse = ", "))
}






set.seed(42)



dml_card_d1 <- dml_with_smoother(Y_dol, D1, X_dml_card, Z,
                                 n_cf_folds = 5, tune.parameters = "all")
res_card_d1 <- summary(dml_card_d1)

cat("Full summary:\n")
print(res_card_d1)

cat("\nIV-relevant rows (PLR-IV, Wald-AIPW):\n")
print(res_card_d1[c("PLR-IV", "Wald-AIPW"), , drop = FALSE])

# Extract outcome weights
# omega$omega is a (rows x N) matrix; rows = PLR / PLR-IV / AIPW-ATE / Wald-AIPW
omega_card_d1       <- get_outcome_weights(dml_card_d1)
check_omega_rows(omega_card_d1)
w_plriv_card_d1     <- as.vector(omega_card_d1$omega["PLR-IV",    ])
w_waldaipw_card_d1  <- as.vector(omega_card_d1$omega["Wald-AIPW", ])

# Algebraic identity: omega'Y = tau_hat must hold exactly (Knaus 2024)
cat("\nAlgebraic check (omega'Y = tau_hat, must be TRUE):\n")
cat("  PLR-IV:   ",
    check_weight_identity(w_plriv_card_d1,    Y_dol, get_estimate(res_card_d1, "PLR-IV")), "\n")
cat("  Wald-AIPW:",
    check_weight_identity(w_waldaipw_card_d1, Y_dol, get_estimate(res_card_d1, "Wald-AIPW")), "\n\n")

cat("Weight diagnostics (Sum_w = 0 => translation invariant):\n")
print(rbind(
  weight_diag(w_plriv_card_d1,    "PLR-IV     (DML, Card, somecol)"),
  weight_diag(w_waldaipw_card_d1, "Wald-AIPW  (DML, Card, somecol)")
), row.names = FALSE)

cat("\nBLOCK B1 done.\n\n")


# --- B2: educ16, Card spec ---
set.seed(42)
cat(strrep("=", 70), "\n")
cat("BLOCK B2 — DML: Card specification, educ16 (D2 = educ >= 16)\n")
cat(strrep("=", 70), "\n\n")

# Treatment: D2 = educ16 = 1{educ >= 16}
# Outcome:   Y_dol = log(wage/100)   (dollars)
# Instrument: Z = nearc4
# Covariates: X_dml_card (Card specification, no intercept)
# Note: The same DML machinery is applied; only D changes.

dml_card_d2 <- dml_with_smoother(Y_dol, D2, X_dml_card, Z,
                                 n_cf_folds = 5, tune.parameters = "all")
res_card_d2 <- summary(dml_card_d2)

cat("Full summary:\n")
print(res_card_d2)

cat("\nIV-relevant rows (PLR-IV, Wald-AIPW):\n")
print(res_card_d2[c("PLR-IV", "Wald-AIPW"), , drop = FALSE])

omega_card_d2       <- get_outcome_weights(dml_card_d2)
check_omega_rows(omega_card_d2)
w_plriv_card_d2     <- as.vector(omega_card_d2$omega["PLR-IV",    ])
w_waldaipw_card_d2  <- as.vector(omega_card_d2$omega["Wald-AIPW", ])

cat("\nAlgebraic check (omega'Y = tau_hat, must be TRUE):\n")
cat("  PLR-IV:   ",
    check_weight_identity(w_plriv_card_d2,    Y_dol, get_estimate(res_card_d2, "PLR-IV")), "\n")
cat("  Wald-AIPW:",
    check_weight_identity(w_waldaipw_card_d2, Y_dol, get_estimate(res_card_d2, "Wald-AIPW")), "\n\n")

cat("Weight diagnostics:\n")
print(rbind(
  weight_diag(w_plriv_card_d2,    "PLR-IV     (DML, Card, educ16)"),
  weight_diag(w_waldaipw_card_d2, "Wald-AIPW  (DML, Card, educ16)")
), row.names = FALSE)

cat("\nBLOCK B2 done.\n\n")



# ==============================================================================
# BLOCK C — DML: KITAGAWA SPECIFICATION
# X = [black, south, smsa, smsa66, south66]  (5 covariates, no intercept)
# Corresponds to columns 3-4 and 7-8 of SUW Table 3.
# ==============================================================================

set.seed(42)
cat(strrep("=", 70), "\n")
cat("BLOCK C1 — DML: Kitagawa specification, somecol (D1 = educ > 12)\n")
cat(strrep("=", 70), "\n\n")


dml_kit_d1 <- dml_with_smoother(Y_dol, D1, X_dml_kit, Z,
                                n_cf_folds = 5, tune.parameters = "all")
res_kit_d1 <- summary(dml_kit_d1)

cat("Full summary:\n")
print(res_kit_d1)

cat("\nIV-relevant rows (PLR-IV, Wald-AIPW):\n")
print(res_kit_d1[c("PLR-IV", "Wald-AIPW"), , drop = FALSE])

omega_kit_d1       <- get_outcome_weights(dml_kit_d1)
check_omega_rows(omega_kit_d1)
w_plriv_kit_d1     <- as.vector(omega_kit_d1$omega["PLR-IV",    ])
w_waldaipw_kit_d1  <- as.vector(omega_kit_d1$omega["Wald-AIPW", ])

cat("\nAlgebraic check (omega'Y = tau_hat, must be TRUE):\n")
cat("  PLR-IV:   ",
    check_weight_identity(w_plriv_kit_d1,    Y_dol, get_estimate(res_kit_d1, "PLR-IV")), "\n")
cat("  Wald-AIPW:",
    check_weight_identity(w_waldaipw_kit_d1, Y_dol, get_estimate(res_kit_d1, "Wald-AIPW")), "\n\n")

cat("Weight diagnostics (Sum_w = 0 => translation invariant):\n")
print(rbind(
  weight_diag(w_plriv_kit_d1,    "PLR-IV     (DML, Kitagawa, somecol)"),
  weight_diag(w_waldaipw_kit_d1, "Wald-AIPW  (DML, Kitagawa, somecol)")
), row.names = FALSE)

cat("\nBLOCK C1 done.\n\n")


# --- C2: educ16, Kitagawa spec ---
set.seed(42)
cat(strrep("=", 70), "\n")
cat("BLOCK C2 — DML: Kitagawa specification, educ16 (D2 = educ >= 16)\n")
cat(strrep("=", 70), "\n\n")

dml_kit_d2 <- dml_with_smoother(Y_dol, D2, X_dml_kit, Z,
                                n_cf_folds = 5, tune.parameters = "all")
res_kit_d2 <- summary(dml_kit_d2)

cat("Full summary:\n")
print(res_kit_d2)

cat("\nIV-relevant rows (PLR-IV, Wald-AIPW):\n")
print(res_kit_d2[c("PLR-IV", "Wald-AIPW"), , drop = FALSE])

omega_kit_d2       <- get_outcome_weights(dml_kit_d2)
check_omega_rows(omega_kit_d2)
w_plriv_kit_d2     <- as.vector(omega_kit_d2$omega["PLR-IV",    ])
w_waldaipw_kit_d2  <- as.vector(omega_kit_d2$omega["Wald-AIPW", ])

cat("\nAlgebraic check (omega'Y = tau_hat, must be TRUE):\n")
cat("  PLR-IV:   ",
    check_weight_identity(w_plriv_kit_d2,    Y_dol, get_estimate(res_kit_d2, "PLR-IV")), "\n")
cat("  Wald-AIPW:",
    check_weight_identity(w_waldaipw_kit_d2, Y_dol, get_estimate(res_kit_d2, "Wald-AIPW")), "\n\n")

cat("Weight diagnostics:\n")
print(rbind(
  weight_diag(w_plriv_kit_d2,    "PLR-IV     (DML, Kitagawa, educ16)"),
  weight_diag(w_waldaipw_kit_d2, "Wald-AIPW  (DML, Kitagawa, educ16)")
), row.names = FALSE)

cat("\nBLOCK C2 done.\n\n")



# ==============================================================================
# BLOCK D — TRANSLATION INVARIANCE: DML VS KAPPA (all 4 DML fits)
# Requires: BLOCKS B and C.
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("BLOCK D — TRANSLATION INVARIANCE: DML VS KAPPA\n")
cat(strrep("=", 70), "\n")
cat(sprintf("k = log(100) = %.6f\n\n", k_shift))

# Helper: one row of the translation-invariance table
# (Adapted from kappa replication; now accepts any weight vector generically.)
print_inv_table_dml <- function(spec_label, D,
                                w_plriv, w_waldaipw,
                                X_kappa, Y_dol, Y_cnt) {
  cat(sprintf("\n--- %s ---\n", spec_label))
  cat(sprintf("  %-26s  %9s  %9s  %10s  %10s  %7s\n",
              "Estimator", "Dollars", "Cents", "Diff(act)", "Diff(pred)", "Match?"))
  cat(strrep("-", 80), "\n")
  
  p_ml  <- logit_mle(Z, X_kappa)
  p_cb  <- get_cbps_p(Z, X_kappa)
  kw_ml <- kappa_outcome_weights(Z, D, p_ml)
  kw_cb <- kappa_outcome_weights(Z, D, p_cb)
  
  rows <- list(
    list(name = "PLR-IV (DML)",
         w = w_plriv,
         ed = sum(w_plriv * Y_dol),
         ec = sum(w_plriv * Y_cnt)),
    list(name = "Wald-AIPW (DML)",
         w = w_waldaipw,
         ed = sum(w_waldaipw * Y_dol),
         ec = sum(w_waldaipw * Y_cnt)),
    list(name = "tau_cb_u",
         w = kw_cb$w_u,
         ed = tau_u(Y_dol, Z, D, p_cb),
         ec = tau_u(Y_cnt, Z, D, p_cb)),
    list(name = "tau_ml_u",
         w = kw_ml$w_u,
         ed = tau_u(Y_dol, Z, D, p_ml),
         ec = tau_u(Y_cnt, Z, D, p_ml)),
    list(name = "tau_ml_a10",
         w = kw_ml$w_a10,
         ed = tau_a10(Y_dol, Z, D, p_ml),
         ec = tau_a10(Y_cnt, Z, D, p_ml)),
    list(name = "tau_ml_a",
         w = kw_ml$w_a,
         ed = tau_unnorm(Y_dol, Z, D, p_ml, "a"),
         ec = tau_unnorm(Y_cnt, Z, D, p_ml, "a")),
    list(name = "tau_ml_a1/t",
         w = kw_ml$w_a1,
         ed = tau_unnorm(Y_dol, Z, D, p_ml, "a1"),
         ec = tau_unnorm(Y_cnt, Z, D, p_ml, "a1")),
    list(name = "tau_ml_a0",
         w = kw_ml$w_a0,
         ed = tau_unnorm(Y_dol, Z, D, p_ml, "a0"),
         ec = tau_unnorm(Y_cnt, Z, D, p_ml, "a0"))
  )
  
  for (r in rows) {
    act  <- r$ec - r$ed
    pred <- sum(r$w) * k_shift
    ok   <- isTRUE(all.equal(act, pred, tolerance = 1e-6))
    cat(sprintf("  %-26s  %9.4f  %9.4f  %10.4f  %10.4f  %7s\n",
                r$name, r$ed, r$ec, act, pred,
                if (ok) "TRUE" else "FALSE"))
  }
}


# ==============================================================================
# BLOCK E — WEIGHT DIAGNOSTICS: DML VS KAPPA (all 4 DML fits)
# Requires: BLOCKS B, C, and the kappa_weights_bundle() function.
# ==============================================================================


cat(strrep("=", 70), "\n")
cat("BLOCK E — WEIGHT DIAGNOSTICS: DML VS KAPPA\n")
cat(strrep("=", 70), "\n\n")
cat("Sum_w = 0 => translation invariant | ESS = 1 / sum(w_i^2)\n\n")


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

# --- E1: somecol, Card spec ---
cat("--- Card specification, somecol (D1) ---\n")
kw_card_d1 <- kappa_weights_bundle(Z, D1, X_card_kappa)

diag_card_d1 <- rbind(
  weight_diag(w_plriv_card_d1,       "PLR-IV        (DML, Card, somecol)"),
  weight_diag(w_waldaipw_card_d1,    "Wald-AIPW     (DML, Card, somecol)"),
  weight_diag(kw_card_d1$tau_cb_u,   "tau_cb_u      (kappa, Card, somecol)"),
  weight_diag(kw_card_d1$tau_ml_u,   "tau_ml_u      (kappa, Card, somecol)"),
  weight_diag(kw_card_d1$tau_ml_a10, "tau_ml_a10    (kappa, Card, somecol)"),
  weight_diag(kw_card_d1$tau_ml_a,   "tau_ml_a      (kappa, Card, somecol)"),
  weight_diag(kw_card_d1$tau_ml_a1,  "tau_ml_a1/t   (kappa, Card, somecol)"),
  weight_diag(kw_card_d1$tau_ml_a0,  "tau_ml_a0     (kappa, Card, somecol)")
)


# --- E2: educ16, Card spec ---
cat("\n--- Card specification, educ16 (D2) ---\n")
kw_card_d2 <- kappa_weights_bundle(Z, D2, X_card_kappa)

diag_card_d2 <- rbind(
  weight_diag(w_plriv_card_d2,       "PLR-IV        (DML, Card, educ16)"),
  weight_diag(w_waldaipw_card_d2,    "Wald-AIPW     (DML, Card, educ16)"),
  weight_diag(kw_card_d2$tau_cb_u,   "tau_cb_u      (kappa, Card, educ16)"),
  weight_diag(kw_card_d2$tau_ml_u,   "tau_ml_u      (kappa, Card, educ16)"),
  weight_diag(kw_card_d2$tau_ml_a10, "tau_ml_a10    (kappa, Card, educ16)"),
  weight_diag(kw_card_d2$tau_ml_a,   "tau_ml_a      (kappa, Card, educ16)"),
  weight_diag(kw_card_d2$tau_ml_a1,  "tau_ml_a1/t   (kappa, Card, educ16)"),
  weight_diag(kw_card_d2$tau_ml_a0,  "tau_ml_a0     (kappa, Card, educ16)")
)


# --- E3: somecol, Kitagawa spec ---
cat("\n--- Kitagawa specification, somecol (D1) ---\n")
kw_kit_d1 <- kappa_weights_bundle(Z, D1, X_kit_kappa)

diag_kit_d1 <- rbind(
  weight_diag(w_plriv_kit_d1,       "PLR-IV        (DML, Kitagawa, somecol)"),
  weight_diag(w_waldaipw_kit_d1,    "Wald-AIPW     (DML, Kitagawa, somecol)"),
  weight_diag(kw_kit_d1$tau_cb_u,   "tau_cb_u      (kappa, Kitagawa, somecol)"),
  weight_diag(kw_kit_d1$tau_ml_u,   "tau_ml_u      (kappa, Kitagawa, somecol)"),
  weight_diag(kw_kit_d1$tau_ml_a10, "tau_ml_a10    (kappa, Kitagawa, somecol)"),
  weight_diag(kw_kit_d1$tau_ml_a,   "tau_ml_a      (kappa, Kitagawa, somecol)"),
  weight_diag(kw_kit_d1$tau_ml_a1,  "tau_ml_a1/t   (kappa, Kitagawa, somecol)"),
  weight_diag(kw_kit_d1$tau_ml_a0,  "tau_ml_a0     (kappa, Kitagawa, somecol)")
)


# --- E4: educ16, Kitagawa spec ---
cat("\n--- Kitagawa specification, educ16 (D2) ---\n")
kw_kit_d2 <- kappa_weights_bundle(Z, D2, X_kit_kappa)

diag_kit_d2 <- rbind(
  weight_diag(w_plriv_kit_d2,       "PLR-IV        (DML, Kitagawa, educ16)"),
  weight_diag(w_waldaipw_kit_d2,    "Wald-AIPW     (DML, Kitagawa, educ16)"),
  weight_diag(kw_kit_d2$tau_cb_u,   "tau_cb_u      (kappa, Kitagawa, educ16)"),
  weight_diag(kw_kit_d2$tau_ml_u,   "tau_ml_u      (kappa, Kitagawa, educ16)"),
  weight_diag(kw_kit_d2$tau_ml_a10, "tau_ml_a10    (kappa, Kitagawa, educ16)"),
  weight_diag(kw_kit_d2$tau_ml_a,   "tau_ml_a      (kappa, Kitagawa, educ16)"),
  weight_diag(kw_kit_d2$tau_ml_a1,  "tau_ml_a1/t   (kappa, Kitagawa, educ16)"),
  weight_diag(kw_kit_d2$tau_ml_a0,  "tau_ml_a0     (kappa, Kitagawa, educ16)")
)



# ==============================================================================
# BLOCK F — LOVE PLOTS: DML VS KAPPA
# Requires: BLOCKS B and E.
# ==============================================================================





kw_card_d1 <- kappa_weights_bundle(Z, D1, X_card_kappa)
kw_card_d2 <- kappa_weights_bundle(Z, D2, X_card_kappa)
kw_kit_d1  <- kappa_weights_bundle(Z, D1, X_kit_kappa)
kw_kit_d2  <- kappa_weights_bundle(Z, D2, X_kit_kappa)

# --- F1: Card spec, somecol — all 8 estimators ---
cat("F1: Love plots — Card specification, somecol (D1)\n\n")

# Love plot helper for Card spec (somecol, D1)
# Uses X_dml_card for the covariate reference matrix.
make_love_card_d1 <- function(title_str, w_vec) {
  # Multiply by (2D1 - 1) to convert signed outcome weights to the
  # treated-vs-control convention expected by cobalt::love.plot().
  # This is identical to the Vietnam DML Block F approach.
  cobalt::love.plot(
    D1 ~ X_dml_card,
    weights    = w_vec * (2 * D1 - 1),
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

lp_f1_1 <- make_love_card_d1("PLR-IV ",        w_plriv_card_d1)
lp_f1_2 <- make_love_card_d1("Wald-AIPW ",     w_waldaipw_card_d1)
lp_f1_3 <- make_love_card_d1("tau_cb_u ",      kw_card_d1$tau_cb_u)
lp_f1_4 <- make_love_card_d1("tau_ml_u ",      kw_card_d1$tau_ml_u)
lp_f1_5 <- make_love_card_d1("tau_ml_a10 ",    kw_card_d1$tau_ml_a10)
lp_f1_6 <- make_love_card_d1("tau_ml_a ",      kw_card_d1$tau_ml_a)
lp_f1_7 <- make_love_card_d1("tau_ml_a1/t ",   kw_card_d1$tau_ml_a1)
lp_f1_8 <- make_love_card_d1("tau_ml_a0 ",     kw_card_d1$tau_ml_a0)




# --- F2: Card vs Kitagawa — DML only, somecol ---

cat("\nF2: Love plots — DML only, Card vs Kitagawa, somecol (D1)\n\n")

make_love_kit_d1 <- function(title_str, w_vec) {
  cobalt::love.plot(
    D1 ~ X_dml_kit,
    weights    = w_vec * (2 * D1 - 1),
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

lp_f2_1 <- make_love_card_d1("PLR-IV Card",      w_plriv_card_d1)
lp_f2_2 <- make_love_card_d1("Wald-AIPW Card ",   w_waldaipw_card_d1)
lp_f2_3 <- make_love_kit_d1( "PLR-IV Kitagawa ",  w_plriv_kit_d1)
lp_f2_4 <- make_love_kit_d1( "Wald-AIPW Kitagawa", w_waldaipw_kit_d1)




# --- F3: Card spec, educ16 — all 8 estimators ---
cat("\nF3: Love plots — Card specification, educ16 (D2)\n\n")

make_love_card_d2 <- function(title_str, w_vec) {
  cobalt::love.plot(
    D2 ~ X_dml_card,
    weights    = w_vec * (2 * D2 - 1),
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

lp_f3_1 <- make_love_card_d2("PLR-IV",        w_plriv_card_d2)
lp_f3_2 <- make_love_card_d2("Wald-AIPW",     w_waldaipw_card_d2)
lp_f3_3 <- make_love_card_d2("tau_cb_u",      kw_card_d2$tau_cb_u)
lp_f3_4 <- make_love_card_d2("tau_ml_u",      kw_card_d2$tau_ml_u)
lp_f3_5 <- make_love_card_d2("tau_ml_a10",    kw_card_d2$tau_ml_a10)
lp_f3_6 <- make_love_card_d2("tau_ml_a",      kw_card_d2$tau_ml_a)
lp_f3_7 <- make_love_card_d2("tau_ml_a1/t",   kw_card_d2$tau_ml_a1)
lp_f3_8 <- make_love_card_d2("tau_ml_a0",     kw_card_d2$tau_ml_a0)



# ==============================================================================
# BLOCK G — SUMMARY TABLE: DML VS KAPPA POINT ESTIMATES
# Requires: BLOCKS B and C.
# ==============================================================================



# Helper: extract PLR-IV and Wald-AIPW from a dml summary matrix
dml_row <- function(res, estimator) {
  coef_val <- get_estimate(res, estimator)
  se_val   <- if ("SE" %in% colnames(res)) as.numeric(res[estimator, "SE"]) else NA_real_
  c(coef = coef_val, se = se_val)
}

# Recompute kappa SEs for dollars outcome (fresh call for self-containedness)
cat("Computing kappa point estimates + SEs (dollars, all 4 cells)...\n")
kappa_card_d1 <- kappa_analytic_se_all(Y_dol, Z, D1, X_card_kappa)
kappa_card_d2 <- kappa_analytic_se_all(Y_dol, Z, D2, X_card_kappa)
kappa_kit_d1  <- kappa_analytic_se_all(Y_dol, Z, D1, X_kit_kappa)
kappa_kit_d2  <- kappa_analytic_se_all(Y_dol, Z, D2, X_kit_kappa)
cat("Done.\n\n")

# Print function for one cell column
print_cell <- function(coef, se) {
  if (is.na(se) || !is.finite(se) || se <= 0)
    return(sprintf("%7.3f (  NA  )", coef))
  pval  <- 2 * pnorm(-abs(coef / se))
  stars <- ifelse(pval < 0.01, "***", ifelse(pval < 0.05, "**",
                                             ifelse(pval < 0.10, "*",   "")))
  sprintf("%7.3f%3s (%5.3f)", coef, stars, se)
}



# DML rows
plriv_card_d1   <- dml_row(res_card_d1, "PLR-IV")
plriv_card_d2   <- dml_row(res_card_d2, "PLR-IV")
plriv_kit_d1    <- dml_row(res_kit_d1,  "PLR-IV")
plriv_kit_d2    <- dml_row(res_kit_d2,  "PLR-IV")
waipw_card_d1   <- dml_row(res_card_d1, "Wald-AIPW")
waipw_card_d2   <- dml_row(res_card_d2, "Wald-AIPW")
waipw_kit_d1    <- dml_row(res_kit_d1,  "Wald-AIPW")
waipw_kit_d2    <- dml_row(res_kit_d2,  "Wald-AIPW")


# Kappa rows
kappa_estimators <- c("tau_cb_u", "tau_ml_u", "tau_ml_a10",
                      "tau_ml_a", "tau_ml_t", "tau_ml_a0")
kappa_labels     <- c("tau_cb_u (kappa)",   "tau_ml_u (kappa)",
                      "tau_ml_a10 (kappa)", "tau_ml_a (kappa)",
                      "tau_ml_t (kappa)",   "tau_ml_a0 (kappa)")









