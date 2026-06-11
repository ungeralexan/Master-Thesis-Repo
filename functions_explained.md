# Function Documentation — `functions_all.R`
For
### Thesis: *Comparing Kappa Weighting and Causal Machine Learning Estimators via Weight-Based Diagnostics*
 
> This document will show the central functions used in my Master Thesis and explain what it entails


 
---
 

## 1 `safe_logit()`
 
**Section:** §a | **Part:** Part 1 — Kappa Estimators
---
 
### The code
 
```r
safe_logit <- function(eta) {
  eta <- pmin(pmax(eta, -35), 35)
  1 / (1 + exp(-eta))
}
```

#### What it does 
`safe_logit()` takes a real number (or a vector of real numbers) and maps it to a probability between 0 and 1. It is the standard **logistic (sigmoid) function**. It adds one protection as it clips extreme values before computing, so the calculation never crashes due to floating-point overflow.
It is the **numerical primitive** shared by every propensity score routine in the file.
the input is a numeric scalar or a vector and the output are probability values between 0 and 1, numeric same size as the input eta.





## 2 `prep_design_for_mest()`

---
 
### The code
 
```r
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
```

#### What it does 
It takes the covariate matrix X and standardizes non-intercept column to have mean 0 and standard deviation 1.
The purpose is that when it comes to estimation we beed to compute the Jacobian and invert so we can get the standard errors, If teh Covariates have scales that are different that could make the Jacobian matrix ill conditioned. 
The inputs are a numeric matrix X or data frame. Must contain an intercept column (all 1s) plus covariate columns.
The output is a numeric matrix Same dimensions as `X`. All other columns transformed to mean 0, standard deviation 1.



## 3 `logit_mle()`

---
### The Code
```r
logit_mle <- function(Z, X) {
  df  <- data.frame(Z = Z, X)
  fit <- glm(Z ~ . - 1, data = df, family = binomial(link = "logit"))
  fitted.values(fit)
}

```
### What it does 
The input is Z which is a vector full of Z and 1 per obs, X which is the design covariate matrix. It already includes an intercept column (which is full of 1)
You receive a vector of probabilities, where each one stands for one observation.
They are between 0 and 1 as propensity scores. 
Important we surpress automatic the intercept. X should already have ones. The glm of R would already add another intercept thats why we surpress the automatic one. 
In Słoczyński, Uysal & Wooldridge (2025) mention that if X includes a constant the tau estimators all collapse to the same number, even though in final sample maybe not hold. 


## 4 `fit_logit_alpha()`

---
### The code
```r
fit_logit_alpha <- function(Z, X) {
  Z     <- as.numeric(Z)
  X     <- prep_design_for_mest(X)
  fit   <- glm.fit(x = X, y = Z, family = binomial(link = "logit"))
  alpha <- as.numeric(coef(fit));  alpha[is.na(alpha)] <- 0
  p     <- as.vector(pmin(pmax(safe_logit(X %*% alpha), 1e-8), 1 - 1e-8))
  list(alpha = alpha, p = p, X_used = X)
}
```

### What it does:
Gives coefficients, p and the standardizes X matrix. Will be needed for sandwhich SE formula.
It again takes the instrument vector, the covariate Matrix, inclusive intercept.
It outputs the estimated coeffs, the propensity score which is clipped and the standardised version of the X that was used to get the alphas. 
Important it returns X_used as the coefs (alpha) are on the standardized covariate sclae. If I want later compute X*alpha to reconstruct the linear predictor it is necessary to used the standardized not the origanl X. 
The clip  prevents p=0 or p=1 exactly, which would cause division by zero in the kappa formulas. 
We also prevent perfect collinearity and replace NA with zero which would mean that this varible has no effect. 


## 5 `fit_cbps_alpha()`

---
```r
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
```


### What it does:
This is a handcoded DBPS solver using Newtons method optimization loop.
We wanna find alha such that the wieghted covariate emans are balanced between Z=1 and Z=0.
Input are again instrument and and covariate matrix X a tolerance of 1e-9 and a maximum of iterations of 5000 same as in the SUW paper defined.
As output we have the CBPS coeffs on a standardized sclae, the CBPS propensity scores that are clipped. The X_used which is standardized. Moreover, it shows whether the algorithm converged and actually reached tolerance. Moreover the maximum absolut value of the moment conditions at the final iterate. Showing how close the exact balance is I got. 

The Newton Loop entails:
Initialising from the MLE solution as it is a already decnt guess. Will converge faster then starting from 0.
We first evaluate the moment condition and check how far are we from balance
If the norm is msaller than the best seen so far we save it as best_b
If the norm is below the tolerance we say convergence and stop
Then we compute the Jacobian (the matrix of derivatives of the moment conditions with respect to the coefficients) and solve J * step.
Then a backtracking line search is activated in which we try the full step then half, then quater up to 50 halvings until the momet norm actually decreases. 
I also implemented a ridge fallback as if the Jacobian is singular (collinearity), I add a small positive number along the diagonal (ridge) which makes the matrix invertible again