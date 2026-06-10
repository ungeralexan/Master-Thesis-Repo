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