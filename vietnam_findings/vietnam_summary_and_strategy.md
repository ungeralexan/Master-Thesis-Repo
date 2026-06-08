# Vietnam Section: Findings Summary + Strategic Plan

---

## Part 1 — Key findings as bullet points

Organised by whether they confirm, extend, or qualify the SUW and Knaus papers.

### Findings that CONFIRM Słoczyński, Uysal & Wooldridge (2025)

- **Normalised kappa estimators are translation invariant; unnormalised
  ones are not.** Across the cents/dollars recoding, τ̂ᵤᶜᵇ, τ̂ᵤᵐˡ, and τ̂ₐ,₁₀
  are unchanged, while τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ shift. This reproduces the central
  SUW result exactly (their Table 2, their Proposition 3.2).

- **The linear-age unnormalised estimates change sign under recoding.**
  τ̂ₐ moves from −0.429 (cents) to +0.015 (dollars) — a gap of 0.444 log
  points that reverses the economic conclusion. SUW report precisely this
  sign reversal in the linear-age Angrist application.

- **The saturated specification collapses all estimators to a common value
  (≈0.241).** SUW note this as a demonstration that fully flexible
  propensity-score specifications can remove finite-sample pathologies.
  Your weight diagnostics confirm the mechanism: Σωᵢ → 0 even for the
  unnormalised estimators in the saturated design.

### Findings that CONFIRM Knaus (2024)

- **Every estimator admits an exact outcome-weight representation.** The
  identity τ̂ = ΣωᵢYᵢ holds to machine precision for all kappa estimators
  and for the grf-based DML estimators (PLR-IV, Wald-AIPW), verified via
  check_weight_identity(). This is the core Knaus (2024) claim, confirmed
  on real IV data.

- **Σωᵢ = 0 is the exact finite-sample criterion for translation
  invariance.** Your cents/dollars exercise operationalises Knaus's
  equation (3): the empirical estimate shift equals exactly k·Σωᵢ. This
  links the SUW property to the Knaus weight framework — the bridge that
  is the conceptual core of your thesis.

- **Wald-AIPW with gradient-boosted trees fails the affine-smoother
  condition.** Your XGBoost identity check returns FALSE (residual
  ≈ −3.7×10⁻⁶ cubic, −1.1×10⁻⁵ saturated). This is the empirical
  instantiation of the caveat Knaus documents: XGBoost does not produce an
  exact smoother matrix, so Σωᵢ = 0 is not guaranteed even with the four
  required hyperparameters set.

### Findings that are YOUR OWN CONTRIBUTION (extend both papers)

- **"Design dominates learner."** ESS (≈5–6) and the negative-weight share
  (exactly 54.4%) are virtually identical across ALL estimators — DML and
  kappa, normalised and unnormalised, grf and ranger and XGBoost. Neither
  SUW nor Knaus report this uniformity; it is specific to your cross-estimator,
  cross-learner comparison in a low-dimensional near-random design.

- **Learner choice is empirically irrelevant in this design.** The four
  Wald-AIPW variants (grf 0.229, linear+logit 0.218, ranger 0.246,
  XGBoost 0.246) span only ≈0.03 log points, all consistent with the
  normalised kappa range of 0.20–0.21. Two independent random-forest
  implementations (grf and ranger) from different packages agree to within
  0.003.

- **Covariate balance is achieved equally by all estimators.** Unadjusted
  ASMD ≈0.5 on all three age moments; all eight estimators reduce it below
  0.1. The estimators that FAIL translation invariance achieve the same
  balance as those that pass it — balance is set by the design, not the
  normalisation choice.

- **In the saturated design, the PARAMETRIC learner is the one that
  breaks.** The linear+logistic Wald-AIPW cannot produce outcome weights
  in the saturated spec (singular X'X within folds), while the flexible
  learners (ranger, XGBoost) handle it without issue. A nice inversion of
  the usual intuition that flexible learners are the "risky" choice.

### Findings that QUALIFY / nuance the papers

- **The saturated collapse is design-specific, NOT a general robustness
  property.** Your text correctly warns against over-reading the saturated
  result: unnormalised estimators are not "safe" in general — they are safe
  here only because the age-cell propensity score mechanically forces
  Σωᵢ → 0. This is a more careful framing than a casual reading of SUW
  might suggest.

- **Translation invariance (Σωᵢ=0) is only ONE of three Knaus conditions.**
  Your cents/dollars exercise tests only the first. The other two
  (ΣωᵢDᵢ=+1, Σωᵢ(1−Dᵢ)=−1) are theoretical guarantees for the normalised
  estimators, verified numerically but not the focus of the empirical
  exercise. Worth stating explicitly so the scope of the empirical claim
  is clear.

---

## Part 2 — Suggested conclusion paragraph for the Vietnam section

```latex
\subsection{Summary}
\label{subsec:vietnam_summary}

The Vietnam application establishes the benchmark behaviour of the unified
outcome-weight diagnostics in a clean, low-dimensional, near-random-assignment
design. Three results stand out. First, the cents/dollars recoding reproduces
the central finding of \citet{sloczynski2025abadie}: normalised kappa
estimators are invariant to additive outcome shifts, while unnormalised
estimators are not, with the linear-age specification reversing sign entirely.
The outcome-weight representation explains this mechanically --- the estimate
shift equals exactly $k\sum_i\omega_i$, so translation invariance is
equivalent to $\sum_i\omega_i = 0$ \citep{knaus2024}. Second, point estimates,
effective sample sizes, and negative-weight shares are virtually identical
across all kappa and DML estimators and across grf, ranger, and XGBoost
learners. In this design the instrument structure, not the choice of estimator
or learner, determines the weight properties --- a pattern we term
\emph{design dominates learner}. Third, the diagnostics confirm rather than
challenge the underlying theory: covariate balance is achieved equally by all
estimators, and the only estimator-specific failure --- the inexact
outcome-weight identity for gradient-boosted Wald-AIPW --- is the expected
finite-sample manifestation of the affine-smoother caveat in
\citet{knaus2024}. The Vietnam design is therefore a benchmark precisely
because the diagnostics are almost too clean. The Card application introduces
a richer covariate-adjustment problem, an observational instrument, and
two-sided noncompliance; it is there that the same diagnostics begin to
differentiate across estimators in a substantively meaningful way.
```

---

## Part 3 — Strategic recommendation: how to proceed

You asked whether to (A) globalise + bulletproof the Vietnam R files now, or
(B) start the Card framework first and globalise later. Here is my honest
assessment given your end-of-August deadline.

### Recommendation: **Start Card now (Option B), but build the global file AS you do it.**

The reasoning:

**1. The global file is cheapest to build once, against two use cases.**
If you globalise now against Vietnam only, you will almost certainly discover
when you start Card that some function needs a slightly different signature
(e.g. Card has continuous covariates, two-sided noncompliance, a different
instrument structure). You would then refactor the global file a second time.
Building it while writing Card means you design each function against BOTH
applications from the start — one refactor, not two.

**2. Card is the section with genuine intellectual risk; Vietnam is done.**
Your Vietnam findings are solid and the numbers will barely move under
globalisation (you correctly note the tables are easy to swap). The substance
is settled. Card is where new things can go wrong: richer covariates mean the
propensity score actually matters, ESS may diverge, Kitagawa-style ESS ≈ 1
reliability flags may appear, and the "design dominates learner" story may
break (which would itself be a finding). You want maximum time on the part
that is uncertain, not the part that is finished.

**3. Globalisation is mechanical; Card writing is creative.** Mechanical work
(refactoring, swapping table values) is low-risk and can be done in the final
weeks even under time pressure. Creative work (interpreting Card results,
writing the comparison, handling whatever surprises Card throws up) needs
slack and a fresh mind. Do the creative, uncertain work first while you have
runway.

### Concrete sequencing for the weeks ahead

**Phase 1 (now → mid-July): Card descriptives + design + first results.**
Write the Card data/design section in the same structure as Vietnam (you
already have the template). Run the kappa replication and the DML extension.
AS you write the Card code, extract every shared function into the global
file — `kappa_outcome_weights()`, `weight_diag()`, `kappa_estimates()`,
`kappa_weights_bundle()`, the DML helpers. Card forces you to write them
generically.

**Phase 2 (mid-July → early August): Angrist–Evans + globalise Vietnam.**
Run the third application. By now the global file is stable (tested against
Card), so retrofitting Vietnam to use it is a quick swap. Replace the Vietnam
table values with the globalised output (they should be identical or nearly
so). Remove the local function redefinitions in the Vietnam DML file.

**Phase 3 (August): cross-application synthesis + polish.**
Write the cross-application comparison (your Table 8). This is where the
three applications pay off: Vietnam (clean), Card (richer), Angrist–Evans
(near one-sided noncompliance). Final proofreading, figure polish, appendix
tables, abstract.

### One caveat on the global file

Do the **two deletions I flagged earlier** (local `weight_diag()` and
`kappa_outcome_weights()` in the Vietnam DML file) BEFORE you start Card,
regardless of which option you choose. It takes five minutes and prevents
the override bug from propagating into Card if you copy the Vietnam DML file
as your Card template. This is the one piece of globalisation worth doing
immediately.

### Bottom line

You are advanced enough that the Vietnam section is effectively done pending
mechanical globalisation. The deadline risk is entirely in Card and
Angrist–Evans, not in Vietnam polish. Spend your best hours on the uncertain
work: start Card now, build the global file as a byproduct of writing Card
generically, and retrofit Vietnam at the end when the global file is already
battle-tested. Polishing a finished section first would be optimising the
wrong thing.
