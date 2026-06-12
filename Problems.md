# Problems and the ideas to solve them 
Typical here I adress problems that arise during the analysis and try to solve them


## Problem status

### Top Tier 

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





