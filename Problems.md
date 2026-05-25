# Problems and the ideas to solve them 
Typical here I adress problems that arise during the analysis and try to solve them


## Problem status

### Top Tier 
- [ ] the third design I want to implement does not work out yet as the sample is too big and I think the smoother matrix explodes 
- [ ] find a way how to implement the tuning vs not tuning design yet and the 5 or 2 fold design
- [ ] base 0 parameter might be dissroting the estimates 
In XGBoost ≥ 1.6 (the version in current mlr3/CRAN) the base_score parameter was changed: it now defaults to 0.5 for regression but the internal centering is done differently, and explicitly setting base_score=0 does NOT produce a zero-mean base learner anymore  instead it biases the initial predictions severely toward zero, which distorts all subsequent tree residuals.
XGBoost < 1.6 it initialised every observation's predicted value to 0, making the first tree's residual equal to the raw outcome. This made the smoother matrix affine (Knaus & Rakov's requirement for omega'Y = tau_hat).
Right now cannot be fixed if I work with my design, there is an alternative but it must be exploited and checked for correctness
-[ ] I figured for instance the PLR-IV and the Kappa weighst do not estimate the exact same target estimator or atleast not always. 
