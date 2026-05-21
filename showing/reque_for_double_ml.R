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

