# Why the Linear Learner Fails on the Saturated Design

## The error, line by line

```
Error in solve.default(crossprod(X), tol = 2.225074e-308) :
  Lapack routine dgesv: system is exactly singular: U[12,12] = 0
```

This single line contains four pieces of information. Read it from the inside out.

- **`crossprod(X)`** is R shorthand for the matrix product $X'X$ (X-transpose times X), the so-called *Gram matrix* of the covariate design $X$. For the saturated specification, $X$ is the matrix of age-cell dummies, so $X'X$ is a square matrix whose dimension equals the number of age cells.
- **`solve.default(...)`** is R's general-purpose matrix inverter. `solve(A)` computes $A^{-1}$; `solve(A, b)` solves $Ax = b$. Here it is being asked to invert (or solve a system involving) $X'X$.
- **`Lapack routine dgesv`** is the underlying numerical linear-algebra routine (from the LAPACK library) that actually performs the solve via LU decomposition. `dgesv` = **d**ouble-precision **ge**neral **s**ol**v**e.
- **`system is exactly singular: U[12,12] = 0`** is the failure. During LU decomposition, the matrix is factored into a lower-triangular $L$ and an upper-triangular $U$. A matrix is invertible if and only if every diagonal entry of $U$ is non-zero. LAPACK reports that the **12th diagonal pivot is exactly zero** — meaning column 12 of the design carries no information that is not already contained in the other columns. The matrix has no inverse.

In plain terms: **R was asked to invert a matrix that cannot be inverted, because one of its columns is a perfect linear combination of the others.**

---

## What "singular" means and why `solve()` needs a non-singular matrix

A matrix $X'X$ is *singular* when its columns are **linearly dependent** — at least one column can be written as a weighted sum of the others. Geometrically, the columns do not span the full space; they collapse onto a lower-dimensional subspace. Algebraically, the determinant is zero, and a zero-determinant matrix has no inverse, exactly as the number 0 has no reciprocal.

The linear smoother that `get_outcome_weights()` reconstructs for an OLS-type learner is

$$
S \;=\; X\,(X'X)^{-1}X',
$$

the **hat matrix** of linear regression. The whole construction hinges on the term $(X'X)^{-1}$. If $X'X$ is singular, that inverse does not exist, `solve()` aborts, and there is no smoother matrix — hence no outcome weights. This is why the crash occurs *inside* `get_outcome_weights()` and not during model fitting.

---

## Why it appears specifically in the *saturated* specification

The three specifications differ only in how age enters the design:

| Specification | Design columns | Approximate column count |
|---|---|---|
| Linear | intercept, age | 2 |
| Cubic | intercept, age, age², age³ | 4 |
| **Saturated** | **one dummy per age cell** | **as many as there are cells** |

The saturated design is, by construction, the **most flexible** — it places a separate parameter on every distinct age value. That flexibility is exactly what causes the problem. Two mechanisms combine:

**1. The dummy-variable trap.** A full set of category dummies plus an intercept is always perfectly collinear: the dummies sum to the intercept column (every observation falls into exactly one cell, so the row-sum of the dummies is 1 for every row, which *is* the intercept). Unless one cell is dropped as the reference, $X'X$ is automatically rank-deficient. The pivot failure at column 12 is the LU decomposition hitting the first dummy that is redundant given the intercept and the earlier dummies.

**2. Cross-fitting makes it worse.** DoubleML splits the sample into 5 folds and re-fits the nuisance models on each training partition. Within a single fold, a sparse age cell may contain **very few or even zero observations**. A dummy column that is all zeros (or constant) on a fold's training data contributes nothing — it is a zero column, the cleanest possible source of singularity. The full sample might support all cells, but a 4/5 subsample need not. So even if you dropped the reference category, a fold-specific empty cell can still produce an exactly singular $X'X$ for that fold.

The cubic and linear designs never hit this because their columns (age, age², age³) are smooth, continuous, and always linearly independent for any subsample containing at least four distinct ages — which every fold does. There is simply nothing to collapse.

---

## Why ranger and XGBoost survive but the linear learner does not

This is the crucial diagnostic point, and it is not a coincidence:

- **The linear learner (`regr.lm`)** builds its smoother by an explicit matrix inversion $(X'X)^{-1}$. That operation requires full column rank. The saturated design denies it that rank within at least one fold, so the inversion fails.
- **Tree-based learners (`ranger`, `xgboost`)** build their smoother matrices from **in-bag sample membership** — which training observations land in which leaf — not from any linear projection. They never form $X'X$ and never invert anything. Collinear or empty dummy columns are completely harmless to them; a tree simply does not split on a column that carries no information. This is why your `ranger` block returned `TRUE` and `xgboost` returned `FALSE` (the affine-smoother caveat) but **both still produced weight vectors**, while only the **linear** learner crashed.

Note also that `regr.lm` *fit* without error and printed a valid point estimate (0.2468). That is because `lm()` internally uses a **rank-revealing QR decomposition** that silently drops aliased columns (it sets their coefficients to `NA` and proceeds). The point estimate is therefore well-defined. It is only the *separate, explicit* smoother reconstruction in `get_outcome_weights()` — which calls `solve(crossprod(X))` directly, without QR pivoting — that has no such fallback and fails hard.

---

## What this means substantively (and why it is a finding, not a bug)

The failure carries a genuine methodological message that fits the *design dominates learner* theme:

> In the most nonparametric covariate specification, it is the **parametric** DML learner whose outcome-weight representation breaks down, while the **flexible** learners handle it without difficulty.

The parametric linear smoother is simply not well-defined on a rank-deficient age-cell design once cross-fitting thins out the sparse cells. This is not a coding mistake — it is an intrinsic property of trying to fit a saturated linear projection on a finite sample split into folds. The flexible learners, which the literature often treats as the "riskier" or "less transparent" choice, are precisely the ones that remain valid here.

---

## What to report

The honest and clean choice is to **omit the linear+logistic learner from the saturated comparison** and state why. Concretely:

1. Report the linear+logistic **point estimate** (0.2468) from `$summary()` — it is valid, since `lm()` handled the collinearity via pivoting.
2. Report **no outcome-weight diagnostics** for the linear learner in the saturated spec, with a one-line note that the OLS smoother is rank-deficient on the age-cell design within folds.
3. Keep ranger (identity check `TRUE`, $\sum_i\omega_i = 0$, ESS = 5) and XGBoost (identity `FALSE`, $\sum_i\omega_i = -1.1\times 10^{-5}$, the affine-smoother caveat) as the two DoubleML learners with valid weight representations in the saturated case.

A suggested footnote for the thesis:

> In the saturated specification the linear+logistic learner is omitted from
> the outcome-weight comparison. Its point estimate is well-defined
> ($\hat{\tau} = 0.247$), but the explicit OLS smoother matrix
> $S = X(X'X)^{-1}X'$ required by \texttt{get\_outcome\_weights()} is not:
> the age-cell design is rank-deficient within cross-fitting folds, so
> $X'X$ is singular and cannot be inverted. Tree-based learners are
> unaffected, as their smoother matrices derive from in-bag sample
> membership rather than a linear projection.

---

## The code fix (so the script does not halt)

Wrap the weight extraction in `tryCatch()` so a singular fold is handled gracefully rather than stopping the whole script:

```r
w_waipw_lm_sat <- tryCatch({
  omega_lm_sat <- get_outcome_weights(object = iivm_lm_sat, dml_data = dml_data_sat)
  ok_lm_sat <- check_doubleml_identity(omega_lm_sat, Y_dol, iivm_lm_sat$coef)
  cat(sprintf("Algebraic check (omega'Y = tau_hat): %s\n\n", ok_lm_sat))
  as.vector(omega_lm_sat$omega)
}, error = function(e) {
  cat("NOTE: OLS smoother is singular on the saturated age-cell design\n")
  cat("      (rank-deficient X'X within folds). Outcome weights undefined.\n")
  cat("      Point estimate remains valid:",
      round(as.numeric(iivm_lm_sat$coef), 4), "\n\n")
  NULL
})
```

Then in Block E-SAT, guard the linear row so the comparison table simply
skips it instead of erroring on `object 'w_waipw_lm_sat' not found`:

```r
if (!is.null(w_waipw_lm_sat)) {
  diag_rows_sat <- rbind(
    diag_rows_sat,
    weight_diag(w_waipw_lm_sat, "Wald-AIPW (DoubleML linear/logit)")
  )
}
```
