# ==============================================================================
# MANUAL LOVE PLOTS — Knaus-style, no cobalt package
# Works for both SECTION 5 (Angrist 1990) and SECTION 6 (Card 1995)
#
# A Love plot shows, for each covariate X_k:
#   - Unadjusted SMD: standardized mean difference between D=1 and D=0 groups
#     (raw, ignoring any weights)
#   - Adjusted SMD:   weighted SMD using outcome weights omega_i
#
# SMD_k = |mean_treated(X_k) - mean_control(X_k)| / SD_pooled(X_k)
#
# For outcome weights, the "weighted mean" of X_k in treated group is:
#   mean_treated_w(X_k) = sum(omega_i * X_ik * I(D_i=1)) / sum(omega_i * I(D_i=1))
# But omega_i are signed (negative for non-compliers), so we follow Knaus:
#   multiply omega_i by (2*D_i - 1) to get "participation weights" that are
#   always positive in the treated group — then use those for cobalt-style SMD.
#
# The threshold line at |SMD| = 0.1 is the conventional "good balance" cutoff.
# ==============================================================================

library(ggplot2)
library(viridis)   # color-blind-friendly palette


# ------------------------------------------------------------------------------
# Core function: compute weighted SMD for one covariate
# omega_signed: outcome weights (signed, as returned by get_outcome_weights)
# D:            treatment indicator (0/1)
# x:            covariate vector (numeric)
# ------------------------------------------------------------------------------
weighted_smd <- function(omega_signed, D, x) {
  # Convert signed outcome weights to "participation weights" (Knaus convention)
  # w_i = omega_i * (2*D_i - 1)
  # This flips the sign for control units, making all weights positive for
  # the group they belong to — compatible with cobalt-style SMD computation
  w <- omega_signed * (2 * D - 1)

  # Weighted means in each group
  w_treat <- w[D == 1];  x_treat <- x[D == 1]
  w_ctrl  <- w[D == 0];  x_ctrl  <- x[D == 0]

  mean_t <- sum(w_treat * x_treat) / sum(w_treat)
  mean_c <- sum(w_ctrl  * x_ctrl)  / sum(w_ctrl)

  # Pooled SD: use the UNWEIGHTED standard deviation (standard cobalt convention)
  # so the denominator is the same for unadjusted and adjusted SMD
  sd_pool <- sd(x)

  if (sd_pool < 1e-10) return(NA_real_)   # constant covariate — skip
  abs(mean_t - mean_c) / sd_pool
}


# Unadjusted SMD (no weights — same for all estimators, just raw group diff)
unadjusted_smd <- function(D, x) {
  sd_pool <- sd(x)
  if (sd_pool < 1e-10) return(NA_real_)
  abs(mean(x[D == 1]) - mean(x[D == 0])) / sd_pool
}


# ------------------------------------------------------------------------------
# Main Love plot function
#
# X_mat       : covariate matrix (N x K), NO intercept column
# var_names   : character vector of length K with display names for covariates
# D           : treatment indicator (length N)
# weight_list : named list of outcome weight vectors, e.g.:
#               list("tau_u" = w_u, "PLR-IV" = w_plriv, ...)
# title       : plot title string
# threshold   : |SMD| threshold line (default 0.1)
# ------------------------------------------------------------------------------
love_plot_manual <- function(X_mat, var_names, D, weight_list,
                             title = "Love Plot",
                             threshold = 0.1) {

  K <- ncol(X_mat)
  stopifnot(length(var_names) == K)

  # --- 1. Unadjusted SMDs (one per covariate, same regardless of estimator) ---
  smd_unadj <- sapply(seq_len(K), function(k) unadjusted_smd(D, X_mat[, k]))

  # --- 2. Adjusted SMDs for each estimator ---
  estimator_names <- names(weight_list)
  n_est <- length(weight_list)

  # Build long-format data frame: one row per (estimator, covariate) pair
  rows <- vector("list", n_est * K)
  idx  <- 1L
  for (e in seq_len(n_est)) {
    omega <- weight_list[[e]]
    for (k in seq_len(K)) {
      smd_adj <- weighted_smd(omega, D, X_mat[, k])
      rows[[idx]] <- data.frame(
        Covariate  = var_names[k],
        Estimator  = estimator_names[e],
        SMD        = smd_adj,
        Type       = "Adjusted",
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  df_adj <- do.call(rbind, rows)

  # Add unadjusted row (one per covariate, labelled "Unadjusted")
  df_unadj <- data.frame(
    Covariate = var_names,
    Estimator = "Unadjusted",
    SMD       = smd_unadj,
    Type      = "Unadjusted",
    stringsAsFactors = FALSE
  )

  df_all <- rbind(df_unadj, df_adj)

  # Order covariates by unadjusted SMD (largest at top — Knaus convention)
  covariate_order <- var_names[order(smd_unadj, decreasing = FALSE)]
  df_all$Covariate <- factor(df_all$Covariate, levels = covariate_order)

  # --- 3. Colors and shapes ---
  all_estimators <- c("Unadjusted", estimator_names)
  n_colors       <- length(all_estimators)

  # Unadjusted always gets a neutral dark grey; adjusted get viridis palette
  pal <- c("grey30", viridis(n_est, option = "D", begin = 0.1, end = 0.85))
  names(pal) <- all_estimators

  # Shapes: circle for unadjusted, then triangle/diamond/square/cross for adjusted
  shape_vals <- c(16, 17, 18, 15, 8, 3)[seq_len(n_colors)]
  names(shape_vals) <- all_estimators

  # Line types: dotted for unadjusted, solid for adjusted
  lty_vals <- c("dotted", rep("solid", n_est))
  names(lty_vals) <- all_estimators

  # Size: unadjusted slightly smaller
  size_vals <- c(2.5, rep(3.2, n_est))
  names(size_vals) <- all_estimators

  # Make "Unadjusted" always last in legend (it's the baseline)
  df_all$Estimator <- factor(df_all$Estimator,
                             levels = c(estimator_names, "Unadjusted"))

  # --- 4. Build ggplot ---
  p <- ggplot(df_all, aes(x = SMD, y = Covariate,
                          colour = Estimator, shape = Estimator,
                          linetype = Estimator, size = Estimator)) +

    # Threshold reference line
    geom_vline(xintercept = threshold,
               linetype = "dashed", colour = "firebrick", linewidth = 0.6) +

    # Zero reference line
    geom_vline(xintercept = 0,
               linetype = "solid", colour = "grey70", linewidth = 0.4) +

    # Connect dots across estimators for each covariate (one line per covariate)
    geom_line(aes(group = interaction(Covariate, Type)),
              linewidth = 0.5, alpha = 0.5) +

    # Dots
    geom_point(alpha = 0.9) +

    scale_colour_manual(values = pal)   +
    scale_shape_manual(values  = shape_vals) +
    scale_linetype_manual(values = lty_vals) +
    scale_size_manual(values   = size_vals) +

    # x-axis starts at 0
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +

    # Vertical dashed annotation label
    annotate("text", x = threshold + 0.003, y = 0.6,
             label = paste0("threshold = ", threshold),
             angle = 90, hjust = 0, size = 3, colour = "firebrick") +

    labs(
      title    = title,
      x        = "Absolute Standardized Mean Difference",
      y        = NULL,
      colour   = "Estimator", shape = "Estimator",
      linetype = "Estimator", size  = "Estimator"
    ) +

    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      legend.position = "bottom",
      legend.title    = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92"),
      panel.grid.major.x = element_line(colour = "grey88")
    )

  p
}


# ==============================================================================
# SECTION 5: ANGRIST (1990) — Love plots
# Run AFTER Section 5 code has executed (dml_5, kw5, w_plriv5, w_waldaipw5 exist)
# ==============================================================================

# Covariate matrix for balance check: same as DML input (no intercept)
# age, age^2, age^3
X5_bal      <- as.matrix(X5)   # already defined in Section 5 as data.frame
var_names_5 <- c("Age", "Age²", "Age³")

# Weight list: only IV-relevant DML estimators + recommended kappa estimators
weights_5 <- list(
  "tau_u (kappa)"    = kw5$w_u,
  "tau_a10 (kappa)"  = kw5$w_a10,
  "PLR-IV (DML)"     = w_plriv5,
  "Wald-AIPW (DML)"  = w_waldaipw5
)

love5 <- love_plot_manual(
  X_mat        = X5_bal,
  var_names    = var_names_5,
  D            = D5,
  weight_list  = weights_5,
  title        = "Angrist (1990): Covariate Balance — Outcome Weights",
  threshold    = 0.1
)

print(love5)


# ==============================================================================
# SECTION 6: CARD (1995) — Love plots
# Run AFTER Section 6 code has executed (dml_6, kw6, w_plriv6, w_waldaipw6 exist)
# ==============================================================================

# Covariate matrix: Card (1995) full spec, no intercept
X6_bal <- as.matrix(X6)

# Cleaner display names (strip "reg66x" repetition)
var_names_6 <- c(
  "Black", "South", "SMSA", "SMSA 1966",
  "Region NE", "Region MW", "Region S", "Region W",
  "Region NE66", "Region MW66", "Region S66", "Region W66",
  "Experience", "Experience²"
)

weights_6 <- list(
  "tau_u (kappa)"    = kw6$w_u,
  "tau_a10 (kappa)"  = kw6$w_a10,
  "PLR-IV (DML)"     = w_plriv6,
  "Wald-AIPW (DML)"  = w_waldaipw6
)

love6 <- love_plot_manual(
  X_mat        = X6_bal,
  var_names    = var_names_6,
  D            = D6,
  weight_list  = weights_6,
  title        = "Card (1995): Covariate Balance — Outcome Weights",
  threshold    = 0.1
)

print(love6)




####################################################
### sign flip in outcome weights check
###################################################

cat("sum(w_waldaipw6):", sum(w_waldaipw6), "\n")
cat("sum(w_plriv6):",    sum(w_plriv6), "\n")
cat("sum(w_waldaipw6 * (2*D6-1)):", sum(w_waldaipw6 * (2*D6-1)), "\n")
cat("sum(w_plriv6 * (2*D6-1)):",    sum(w_plriv6 * (2*D6-1)), "\n")
cat("sum(w_waldaipw6[D6==1]):", sum(w_waldaipw6[D6==1]), "\n")
cat("sum(w_waldaipw6[D6==0]):", sum(w_waldaipw6[D6==0]), "\n")
cat("sum(w_plriv6[D6==1]):", sum(w_plriv6[D6==1]), "\n")
cat("sum(w_plriv6[D6==0]):", sum(w_plriv6[D6==0]), "\n")
