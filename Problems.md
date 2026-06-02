# Problems and the ideas to solve them 
Typical here I adress problems that arise during the analysis and try to solve them


## Problem status

### Top Tier 
- [x] **Smoother matrix memory ceiling — Angrist & Evans (1998) DML not feasible**
  `dml_with_smoother()` builds a full N×N smoother matrix to extract outcome weights.
  With N ≈ 254,654 (LFP sample) this requires ~500 GB in float64 — structurally
  impossible regardless of machine memory. Confirmed with the error:
  `Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()`
  on both the LFP outcome (`workedm`, N ≈ 254k) and the income subsample.
  **Resolution:** DML outcome-weight diagnostics are not reported for the AE design.
  This is documented as a design limitation in the thesis (Section 5.2). The near-random
  samesex instrument means "design dominates learner" anyway — the DML comparison
  would add little beyond what the kappa diagnostics already show.
  Subsampling (~5k obs) was considered but rejected: point estimates would no longer
  match the kappa replication table, muddying the comparison.
  **To discuss with professor:** whether OutcomeWeights has a sparse/approximate
  smoother path, or whether a cluster could be used.
- [ ] **Tuning vs. no-tuning and fold choice (5-fold vs. 2-fold)**
  Find a way to implement and compare `tune.parameters = "all"` vs. default and
  `n_cf_folds = 5` vs. `n_cf_folds = 2` in a clean side-by-side design.

- [ ] **XGBoost base_score = 0 distortion**
  In XGBoost ≥ 1.6 the `base_score` parameter was changed: it no longer produces a
  zero-mean base learner when set to 0. This biases initial predictions toward zero and
  distorts all subsequent tree residuals, violating the affine smoother requirement for
  `omega'Y = tau_hat` (Knaus & Rakov). XGBoost < 1.6 initialised every observation's
  predicted value to 0, making the first tree's residual equal to the raw outcome.
  Right now cannot be fixed within the current design. An alternative approach exists
  but must be checked for correctness before use.

- [ ] **PLR-IV vs. kappa: not always the same target estimator**
  PLR-IV and the kappa weights do not always estimate exactly the same target
  (LATE vs. something else depending on compliance structure). Needs to be clarified
  and stated explicitly in the thesis.


### What maybe you have to change in expression
- [ ] Col 1 → Col 3: this is NOT a pure translation. It's Y' = −Y + 1, i.e. a sign flip (scale by −1) plus a translation (+1). A translation-invariant estimator will correctly handle the +1 part, but the scale-by-−1 part means the estimate should flip sign: τ̂(−Y + 1) = −τ̂(Y). Combined: the display column should be −τ̂(blk_lfp_neg).
You're testing:
What happens when you code the outcome as the inverse?
The conceptual content of the displayed column is: "the negative effect of morekids on working", which is the same causal quantity as the positive effect on not-working. For translation-invariant estimators, multiplying displayed values by −1 recovers column 1 exactly.
For unnormalized estimators, this fails in two ways:

1. The +1 translation in 1 − workedm = −workedm + 1 contaminates the estimate.
2. The sign flip then doesn't recover the original, because the translation contamination is asymmetric.