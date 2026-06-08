# Knaus OutcomeWeights: Function-by-Function Walkthrough
# and Vietnam Implementation Consistency Check

---

## Part 1 — What each function does, referenced to the paper

This walks through every function Knaus uses in the 401k notebook, in the
order he uses them, with the exact paper section that motivates each one.

---

### `dml_with_smoother(Y, D, X, Z, n_cf_folds, tune.parameters)`

**What it does.**
This is the entry point for all DML estimation in the OutcomeWeights package.
It runs four estimators simultaneously: PLR, PLR-IV, AIPW-ATE, and Wald-AIPW,
all via cross-fitting with honest generalised random forests (grf).

The critical thing that makes it different from a standard DML implementation
is the name: it computes and *stores* the full N×N **smoother matrices** for
every nuisance function. A smoother matrix S satisfies Ŷ = S·Y — it expresses
each fitted value as a linear combination of all observed outcomes. Standard
DML packages discard the smoother after fitting (they only need Ŷ, not S).
Knaus keeps S because Proposition 1 of the paper shows that any PIVE whose
pseudo-outcome is a linear transformation of Y admits outcome weights exactly
equal to ω = (Z̃'D̃)⁻¹ · T'·ω̃, where T is the smoother matrix stack.

**Paper reference.** Section 2 (PIVE framework) and Section 3 (DML as a
special case of PIVE). The claim that outcome weights exist *because* the
pseudo-outcome is a linear map of Y is Proposition 1. The n_cf_folds argument
controls the cross-fitting structure described in Section 3.2.

**Memory note.** Knaus explicitly warns in the notebook: "the many required
N×N smoother matrices make this notebook relatively memory intensive (32GB
RAM should suffice)." This is why he calls rm(dml_2f); gc() after extracting
the weights — he is freeing the stored smoother matrices.

**In Knaus's code:**
```r
dml_2f = dml_with_smoother(Y, D, X, Z, n_cf_folds = 2)
dml_5f = dml_with_smoother(Y, D, X, Z, n_cf_folds = 5)
```

---

### `summary(dml_object)`

**What it does.**
Returns a data frame of point estimates and standard errors for all four
estimators (PLR, PLR-IV, AIPW-ATE, Wald-AIPW). The row names are the
estimator labels used throughout; you index into this object with
`res["Wald-AIPW", ]` to pull a specific estimate later.

**Paper reference.** Table 2 of the paper reports these estimates. The
summary method is the bridge between the fitted object and the numbers.

**In Knaus's code:**
```r
results_dml_2f = summary(dml_2f)
results_dml_5f = summary(dml_5f)
```

---

### `plot(dml_object)`

**What it does.**
Produces a coefficient plot showing point estimates and confidence intervals
for all four estimators. Used for visual inspection only, not for the paper's
figures.

**Paper reference.** No direct paper section — this is a convenience
diagnostic, not used in the main analysis.

---

### `get_outcome_weights(dml_object)`

**What it does.**
This is the core function of the package. It takes the fitted dml object
(which contains the stored smoother matrices) and computes the per-observation
outcome weights ω_i for each of the four estimators. The returned object
`omega` has a `$omega` slot which is a (4 × N) matrix: 4 rows (one per
estimator) and N columns (one per observation). Row names are
c("PLR", "PLR-IV", "AIPW", "Wald-AIPW").

**What it computes internally.** For a PIVE with smoother T, the outcome
weights are ω = (Z̃'D̃)⁻¹ · T'z̃ where z̃ is the partialled-out instrument.
For Wald-AIPW specifically this simplifies to the expression in Proposition 3
of the paper (the IV smoother weights). The function handles all four
estimator types and the cross-fitting fold structure automatically.

**Paper reference.** Proposition 1 (general PIVE outcome weights), Proposition
3 (Wald-AIPW specific form), and Section 5.1 (the 401k application). The
sentence "we use the get_outcome_weights() method to extract the outcome
weights as described in the paper" in the notebook is the direct bridge.

**In Knaus's code:**
```r
omega_dml_2f = get_outcome_weights(dml_2f)
omega_dml_5f = get_outcome_weights(dml_5f)
```

---

### `omega$omega %*% Y` — the algebraic identity check

**What it does.**
Multiplies the (4 × N) weight matrix by the (N × 1) outcome vector Y.
The result is a (4 × 1) vector of estimated effects — one per estimator.
Knaus then checks `all.equal()` against `results_dml[,1]` (the point
estimates from summary()). If the check passes, the outcome weights are exact,
not approximate. This is the finite-sample verification of Proposition 1.

**Paper reference.** This is the direct numerical verification of the key
claim of the paper: "the weights multiplied by the outcome vector reproduce
the conventionally generated point estimates." Proposition 1 states this
algebraically; `all.equal()` confirms it numerically. Any deviation from TRUE
means the smoother reconstruction is approximate, not exact.

**In Knaus's code:**
```r
all.equal(as.numeric(omega_dml_2f$omega %*% Y),
          as.numeric(results_dml_2f[,1]))
```

---

### `love.plot()` with `weights = omega$omega[index,] * (2*D-1)`

**What it does.**
The `cobalt` package's love.plot() function plots Absolute Standardized Mean
Differences (ASMD) for each covariate, comparing the unadjusted sample to
the weight-adjusted sample.

The critical transformation is `* (2*D-1)`. Outcome weights ω_i are *signed*
— treated units typically carry positive weights and control units typically
carry negative weights (since the estimator is a contrast). The cobalt package
expects *unsigned* importance weights that are non-negative. Multiplying by
(2D_i - 1) flips the sign for control observations (D_i=0 → 2*0-1 = -1,
so the negative weight becomes positive), making the weight vector compatible
with cobalt's balance convention. This is documented in the notebook:
"we need to flip the sign of the untreated outcome weights to make them
compatible with the package framework."

**Paper reference.** Section 4 (covariate balance diagnostics) and Figure 2
of the paper. The Love plots are how the paper assesses whether the implied
outcome weights achieve covariate balance.

**In Knaus's code:**
```r
love.plot(D ~ X,
  weights = list(
    "2-fold" = omega_dml_2f$omega[index, ] * (2*D-1),
    "5-fold" = omega_dml_5f$omega[index, ] * (2*D-1)
  ),
  position = "bottom", title = title,
  thresholds = c(m = 0.1), var.order = "unadjusted",
  binary = "std", abs = TRUE, line = TRUE,
  colors = viridis(3), shapes = c("circle","triangle","diamond")
)
```

`index` selects which estimator row: 1=PLR, 2=PLR-IV, 3=AIPW, 4=Wald-AIPW.

---

### `summary(omega_object)` — weight descriptives

**What it does.**
Reports descriptive statistics of the weight vectors: sum of weights,
ESS (effective sample size = 1/Σω²), share of negative weights, and
maximum absolute weight. The notebook uses this to illustrate that PLR(-IV)
and Wald-AIPW are "only scale-normalized" (their weight sums are close to
but not exactly zero), contrasting with AIPW-ATE which is fully normalised.

**Paper reference.** Tables 4–5 of the paper (normalization classification)
and Section 3.4 (the distinction between scale normalization and full
normalization). This is where Knaus shows that different estimators achieve
different normalization properties.

---

## Part 2 — Vietnam implementation: complete consistency check

I now check every function usage in your Vietnam DML file against Knaus's
exact pattern. Verdict for each: ✓ = correct, ✗ = deviation, ⚠ = minor issue.

---

### `dml_with_smoother()` — **✓ Correct**

**Knaus:** `dml_with_smoother(Y, D, X, Z, n_cf_folds = 5)`
**Vietnam:**
```r
dml_cub <- dml_with_smoother(Y_dol, D, X_dml_cub, Z,
                              n_cf_folds = 5, tune.parameters = "all")
dml_sat <- dml_with_smoother(Y_dol, D, X_dml_sat, Z,
                              n_cf_folds = 5, tune.parameters = "all")
```

✓ Argument order Y, D, X, Z is correct.
✓ n_cf_folds = 5 matches the main spec.
✓ You added `tune.parameters = "all"` — Knaus uses default (no tuning in the
  401k notebook) but tuning is valid and documented. Knaus notes in the 401k
  notebook: "default honest random forest (tuning only increases running time
  without changing the insights)". Your choice to tune is defensible and more
  rigorous for a thesis.
✓ X_dml_cub and X_dml_sat are correctly built via make_X() as numeric
  matrices without intercept — this is required by grf's validate_X().

One thing to confirm: Knaus does NOT drop a reference category from his X
in the 401k notebook (his X has no dummy variables). You correctly drop one
column from the saturated spec: `X_dml_sat <- X_dml_sat_full[, -1, drop=FALSE]`.
This is necessary for grf (it tolerates near-collinearity better than exact
collinearity). ✓

---

### `summary(dml_object)` — **✓ Correct**

**Knaus:** `results_dml_2f = summary(dml_2f)`
**Vietnam:**
```r
res_cub <- summary(dml_cub)
res_sat <- summary(dml_sat)
```

✓ Identical pattern.
✓ You index correctly: `res_cub[c("PLR-IV", "Wald-AIPW"), , drop=FALSE]`
  — the `drop=FALSE` prevents R from collapsing to a vector when you select
  a single row, which is good defensive practice.

---

### `get_outcome_weights()` — **✓ Correct, with one note**

**Knaus:** `omega_dml_2f = get_outcome_weights(dml_2f)`
**Vietnam:**
```r
omega_cub      <- get_outcome_weights(dml_cub)
check_omega_rows(omega_cub)                       # ← your addition
w_plriv_cub    <- as.vector(omega_cub$omega["PLR-IV",    ])
w_waldaipw_cub <- as.vector(omega_cub$omega["Wald-AIPW", ])
```

✓ `get_outcome_weights(dml_cub)` — correct, no extra arguments needed.
✓ You correctly extract individual weight vectors by row name from
  `omega_cub$omega["PLR-IV", ]` and wrap in `as.vector()` to strip the
  matrix dimension.
✓ Your `check_omega_rows()` helper is a good defensive addition — it stops
  the script immediately if the expected row names are absent, catching any
  version-related changes in the package output structure.

**One note:** Knaus indexes by integer position (`omega_dml_2f$omega[2,]`
for PLR-IV in the love.plot code). You index by name ("PLR-IV"). Indexing
by name is strictly safer and better practice — if the package ever reorders
rows, integer indexing silently fails while name indexing errors loudly. ✓

---

### Algebraic identity check — **⚠ Functionally correct, minor difference**

**Knaus:**
```r
all.equal(as.numeric(omega_dml_2f$omega %*% Y),
          as.numeric(results_dml_2f[,1]))
```
He checks ALL four estimators at once via the full (4×N) %*% (N×1) product.

**Vietnam:**
```r
check_weight_identity(w_plriv_cub,    Y_dol, get_estimate(res_cub, "PLR-IV"))
check_weight_identity(w_waldaipw_cub, Y_dol, get_estimate(res_cub, "Wald-AIPW"))
```
You check each extracted weight vector individually via your custom helper:
```r
check_weight_identity <- function(w, Y, estimate, tol = 1e-8) {
  isTRUE(all.equal(sum(w * Y), estimate, tolerance = tol))
}
```

✓ The check is correct: `sum(w * Y)` is the same as `omega[row,] %*% Y` for
  a single row. The result is identical.
⚠ You are only checking PLR-IV and Wald-AIPW (the IV estimators). Knaus
  checks all four including PLR and AIPW-ATE. Since your thesis focuses on
  IV estimators this is fine, but be aware you have not verified the ATE
  estimators. Not a problem unless you use those rows anywhere.
✓ Your tolerance of 1e-8 matches what `all.equal()` uses by default.

---

### `weight_diag()` — **⚠ Local redefinition — potential inconsistency**

**Knaus:** Knaus uses `summary(omega_dml_2f)` (the summary method of the
omega object) for weight descriptives. He does not have a standalone
`weight_diag()` function.

**Vietnam (Block A, local definition):**
```r
weight_diag <- function(w, name) {
  data.frame(
    Estimator = name,
    Sum_w     = round(sum(w), 8),
    ESS       = round(1 / sum(w^2), 0),
    Pct_neg   = round(mean(w < 0) * 100, 1),
    Max_abs_w = round(max(abs(w)), 6)
  )
}
```

**Vietnam (functions_kappa.R, global version):**
The global file also has `weight_diag()` in §t. 

**⚠ CRITICAL ISSUE:** You define `weight_diag()` both locally in Block A of
the DML file AND in the global `functions_kappa.R`. Since you source
`functions_kappa.R` first and then define it again locally, the local
definition silently **overrides** the global one for the rest of the DML
script. If the two definitions differ at all, this will cause inconsistencies
between the kappa diagnostics (which may use the global version) and the
DML diagnostics (which use the local version).

**Check needed:** Verify that both definitions are byte-identical. If they
are, no problem — the redefinition is harmless redundancy. If they differ in
rounding precision, column names, or ESS calculation, your kappa and DML
weight tables will not be directly comparable.

**Recommended fix:** Remove the local definition from Block A of the DML
script. Source `functions_kappa.R` at the top and rely on the global version.
This eliminates the override risk entirely.

---

### `love.plot()` with `* (2*D-1)` — **✓ Correct**

**Knaus:**
```r
weights = omega_dml_2f$omega[index, ] * (2*D-1)
```

**Vietnam:**
```r
make_love <- function(title_str, w_vec) {
  love.plot(D ~ X_dml_cub,
    weights    = w_vec * (2 * D - 1), ...)
}
w_plriv_cub |> make_love("PLR-IV")
```

✓ The `* (2*D-1)` sign flip is correctly applied.
✓ You use the same cobalt arguments: `abs=TRUE`, `binary="std"`,
  `var.order="unadjusted"`, `thresholds=c(m=0.1)`, `line=TRUE`.
✓ You balance on `X_dml_cub` (age, age², age³) — the same covariate matrix
  passed to grf. Knaus balances on `X` (the full covariate set). In your case
  X_dml_cub IS the full covariate set for Vietnam, so this is correct.

**One note:** Knaus overlays 2-fold and 5-fold on the same plot using a named
list in `weights`. You plot each estimator separately. This is a presentational
choice, not a correctness issue. Your approach is cleaner for comparing across
estimators rather than across fold counts. ✓

---

### `kappa_outcome_weights()` — **⚠ Defined locally in DML file, conflicts with global**

This is the most important finding.

**In `functions_kappa.R` (global, §k):** The function is defined there.

**In Block A of your DML file:** You define it again locally:
```r
kappa_outcome_weights <- function(Z, D, p) {
  n  <- length(Z)
  kw <- kappa_weights(Z, D, p)
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  dD <- sum(D * Z / p) / s1 -
        sum(D * (1 - Z) / (1 - p)) / s0
  w_u <- (Z / p / s1 - (1 - Z) / (1 - p) / s0) / dD
  w_a10 <- kw$kappa1 / sum(kw$kappa1) -
           kw$kappa0 / sum(kw$kappa0)
  num_w <- (Z - p) / (p * (1 - p)) / n
  list(w_u=..., w_a10=..., w_a=..., w_a1=..., w_a0=...)
}
```

**⚠ CRITICAL:** Same override issue as `weight_diag()`. The local definition
takes precedence. You need to check whether this local version is identical to
the global §k version. If they differ by even one line — different variable
names, different normalisation in one of the w_ components, different handling
of n — your kappa weights in the DML comparison will be computed from a
different formula than the kappa weights in the standalone kappa script.
This would silently corrupt the cross-table comparisons in Block E.

**Recommended fix:** Delete both the local `kappa_outcome_weights()` and
`weight_diag()` definitions from Block A of the DML file. Source only
`functions_kappa.R`. The global file already has both.

---

## Summary: complete verdict

| Function | Status | Action needed |
|---|---|---|
| `dml_with_smoother()` | ✓ Correct | None |
| `summary(dml)` | ✓ Correct | None |
| `get_outcome_weights()` | ✓ Correct | None |
| `omega %*% Y` identity check | ✓ Correct | None (IV-only is fine for thesis) |
| `love.plot()` with `*(2*D-1)` | ✓ Correct | None |
| `weight_diag()` | ⚠ Override risk | Remove local definition from DML file |
| `kappa_outcome_weights()` | ⚠ Override risk | Remove local definition from DML file |
| `check_weight_identity()` | ✓ Correct | None (local but only in DML file) |
| `get_estimate()` / `check_omega_rows()` | ✓ Good additions | Could promote to global file |

**The two items that need action** are both the same class of issue: you have
local redefinitions of functions that already exist in `functions_kappa.R`.
The fix in both cases is identical — remove the local definition from Block A
of the DML file and rely on the global version. This takes two deletions and
makes the entire codebase consistent across all three applications.

No errors in the actual estimation logic were found. The DML output you
produced (identity checks TRUE for grf, the XGBoost caveat, the weight
diagnostics table) is consistent with correctly applied Knaus methodology.
