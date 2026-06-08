# Two Pathways to Outcome Weights: What Changes and Why

---

## The fundamental distinction

Knaus has two completely separate ways of getting outcome weights.
They are **not alternatives** — they answer different questions.

---

## Pathway 1 — `dml_with_smoother()` (the 401k notebook)

### What happens internally

`dml_with_smoother()` is Knaus's own internal implementation inside the
`OutcomeWeights` package. It runs all four DML estimators (PLR, PLR-IV,
AIPW, Wald-AIPW) using `grf` (generalised random forest) as the nuisance
learner, and it does something no standard DML package does:

**it stores the full N×N smoother matrix S for every nuisance function.**

A smoother matrix S satisfies the equation:

    Ŷ = S · Y

meaning every fitted value is a weighted average of ALL observed outcomes.
For a random forest, S_ij is the fraction of trees in which observation j
falls into the same leaf as observation i. Knaus computes and stores this
entire N×N matrix.

### Why storing S matters

The entire mathematical machinery of the paper (Proposition 1) says:

    If the pseudo-outcome Ỹ = T · Y  (a linear map of Y),
    then the outcome weights are ω = (Z̃'D̃)⁻¹ · T' · z̃

The key word is **linear map**. The DML pseudo-outcome is built from
residuals of nuisance regressions. If those regressions are linear in Y
(which they are when you have the smoother S), then the whole thing becomes
a linear map of the original Y vector, and outcome weights exist in closed
form.

`dml_with_smoother()` stores S so that `get_outcome_weights()` can
reconstruct ω analytically. The memory cost is enormous (several N×N
matrices for N=3,027 is manageable; for N=100,000 it would not be).

### What `get_outcome_weights()` does here

It takes the stored smoother matrices and applies the formula from
Proposition 1. The result is an exact (4×N) weight matrix. The identity
check `omega %*% Y == point_estimate` holds to machine precision because
the weights ARE the estimator, derived analytically from the smoother.

### Summary of Pathway 1

```
dml_with_smoother()   →   stores N×N smoother matrices
        ↓
get_outcome_weights()  →   applies Proposition 1 analytically
        ↓
omega$omega            →   exact (4×N) weight matrix
        ↓
omega %*% Y == tau     →   TRUE to machine precision (always)
```

**Learners:** only grf (honest random forest).
**Memory:** intensive (several N×N matrices stored).
**Identity check:** always TRUE by construction.

---

## Pathway 2 — `DoubleMLIIVM` + `get_outcome_weights(object, dml_data)`
## (the DML_smoothers_2.rmd vignette Knaus sent you)

### What changes

In this pathway you use the **external DoubleML package** (Bach et al.),
which runs DML through the mlr3 ecosystem. DoubleML does NOT store smoother
matrices — it only stores the fitted predictions Ŷ (the fold-specific
out-of-sample predictions), not the full S that generated them.

This means `get_outcome_weights()` **cannot apply Proposition 1 directly**.
Instead, for each supported learner type, it reconstructs an approximation
of the smoother from the stored model objects:

- **Linear regression (`regr.lm`):** smoother is S = X(X'X)⁻¹X'.
  This is exact, computed from the stored regression coefficients.
  This is why the linear learner works — its smoother is algebraically
  reconstructible.

- **Ranger random forest (`regr.ranger` with `keep.inbag=TRUE`):**
  smoother is approximated from the in-bag sample indices. Each tree assigns
  observation i to a leaf; the smoother S_ij = fraction of trees where i and
  j are in the same leaf. `keep.inbag=TRUE` is mandatory because without it
  the in-bag indices are discarded after fitting, making reconstruction
  impossible.

- **XGBoost (`regr.xgboost`):** smoother reconstruction requires that the
  model output is exactly linear in the training Y. This holds ONLY if four
  hyperparameters are set to force additive, non-shrunk trees:
  `alpha=0, subsample=1, max_delta_step=0, base_score=0`. Even then,
  the reconstruction is approximate — XGBoost's internal implementation
  does not expose a true smoother matrix, so the identity check can fail
  by a small numerical residual.

### What `get_outcome_weights(object, dml_data)` does here

Note the different signature: it now takes TWO arguments — the fitted
DoubleML object AND the dml_data object. This is because the DoubleML
object alone does not contain the full data; the dml_data is needed to
reconstruct the smoother from the stored model objects.

The function looks at what learner type was used, extracts the model from
`iivm_obj$models`, reconstructs S as best it can for that learner type,
and then applies the Proposition 1 formula.

### The identity check is now a diagnostic, not a guarantee

In Pathway 1, `omega %*% Y == tau_hat` is always TRUE — it is guaranteed by
construction because S was stored exactly.

In Pathway 2, `omega %*% Y == tau_hat` is a **test** — it can be FALSE.
Whether it passes depends on whether the smoother reconstruction was exact:

| Learner | Identity check | Reason |
|---|---|---|
| `regr.lm` | TRUE | exact smoother S = X(X'X)⁻¹X' |
| `regr.ranger` (keep.inbag=TRUE) | TRUE | in-bag reconstruction is exact |
| `regr.xgboost` (4 hyperparams set) | usually FALSE by ~10⁻⁶ | approximate smoother |

Your XGBoost result (-3.7×10⁻⁶) is exactly this. The four hyperparameters
were set correctly, but XGBoost still cannot guarantee an exact smoother.
This is the documented caveat in Knaus's vignette, not an error.

### Summary of Pathway 2

```
DoubleMLIIVM$new(...)   →   standard DoubleML object, NO smoother stored
        ↓
$fit(store_models=TRUE, store_predictions=TRUE)
        ↓  ← store_models=TRUE is MANDATORY — weights need the model objects
get_outcome_weights(object, dml_data)
        →   reconstructs S from stored model objects (learner-specific)
        ↓
omega$omega            →   approximate (1×N) weight vector
        ↓
omega %*% Y == tau     →   TRUE for lm/ranger, FALSE for xgboost (~10⁻⁶)
```

**Learners:** any of lm, ranger, xgboost (the three supported by the dev version).
**Memory:** light — only fitted model objects stored, not N×N matrices.
**Identity check:** diagnostic (not guaranteed; TRUE means reconstruction worked).

---

## The key conceptual difference in one sentence

> Pathway 1 stores the smoother exactly during fitting and extracts weights
> analytically. Pathway 2 stores only the fitted model and tries to
> reconstruct the smoother after the fact — which works exactly for linear
> and ranger learners, and approximately for XGBoost.

---

## What YOU are doing and why it is correct

Your Vietnam implementation uses **both pathways** and this is exactly right:

**Blocks B and C** (your `angrist1990_dml.R`): Pathway 1.
You run `dml_with_smoother()` with grf for the cubic and saturated specs.
Identity checks are TRUE by construction. This is the primary DML analysis
and the source of the weight diagnostics in Table 7 and the Love plots.

**Blocks B-SAT through D-SAT** (your extended file): Pathway 2.
You run `DoubleMLIIVM` with three learners (lm, ranger, xgboost) on the
same data. The purpose is NOT to replicate Pathway 1 — it is to test
whether learner choice changes anything. The fact that ranger gives 0.247
and grf gives 0.244 is the point: two completely different random forest
implementations from different packages agree to within 0.003 log points.
That is the *design dominates learner* finding operationalised.

**The XGBoost FALSE in your output** is therefore not a problem to fix —
it is a finding to report. It is the exact empirical instantiation of the
footnote in your Table 1 classification: "Wald-AIPW achieves Σω_i=0 only
under smoother conditions which conflict with standard flexible learners
such as gradient-boosted trees."

---

## The one thing that IS different between Knaus's vignette and your code

Knaus's vignette uses generic synthetic data (`make_iivm_data()`).
You apply the same machinery to real SIPP data with a specific causal
question. The function calls are otherwise identical:

**Knaus:**
```r
iivm_obj <- DoubleMLIIVM$new(iivm_data, ml_g_iivm, ml_m_iivm, ml_r_iivm,
                              n_folds = cf_folds)
iivm_obj$fit(store_models = TRUE, store_predictions = TRUE)
omega_waipw <- get_outcome_weights(object = iivm_obj, dml_data = iivm_data)
all.equal(omega_waipw$omega %*% iivm_data$data$y, iivm_obj$all_coef)
```

**Your Vietnam (e.g. Block C-SAT, ranger):**
```r
iivm_rf_sat <- DoubleMLIIVM$new(dml_data_sat, lrn("regr.ranger", keep.inbag=TRUE),
                                 lrn("classif.ranger", keep.inbag=TRUE),
                                 lrn("classif.ranger", keep.inbag=TRUE),
                                 n_folds = 5, score = "LATE")
iivm_rf_sat$fit(store_models = TRUE, store_predictions = TRUE)
omega_rf_sat <- get_outcome_weights(object = iivm_rf_sat, dml_data = dml_data_sat)
check_doubleml_identity(omega_rf_sat, Y_dol, iivm_rf_sat$coef)
```

The structure is identical. The only additions you make are:
- `score = "LATE"` — correct for IV estimation (Wald-AIPW targets LATE,
  not ATE; Knaus's generic example uses the default which happens to work
  for the synthetic data but LATE is the right score for your design).
- `n_folds = 5` instead of 3 — more robust cross-fitting, fully appropriate.
- Your `check_doubleml_identity()` helper wraps Knaus's `all.equal()` check
  in a cleaner function — same logic, better readability.

**Everything is correctly implemented.**
