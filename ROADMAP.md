# Thesis Roadmap

**Last updated:** May 2026  
**Phase:** Replication complete → Writing active

---

## Reading status

### ✅ Tier 0 — Core papers (DONE)

| Paper | Status | Key takeaways |
|---|---|---|
| Imbens & Angrist (1994) | ✅ Read | 3 conditions, 3 theorems, LATE as Wald ratio |
| Angrist, Imbens & Rubin (1996) | ✅ Read | Four compliance types, AIR notation, defier exclusion |
| Abadie (2003) | ✅ Read | Kappa theorem, κ₀/κ₁/κ, Lemma 2.1, Prop 5.1 |
| Słoczyński, Uysal & Wooldridge (2025) | ✅ Read | All sections. Five estimators, Prop 3.2 (TI + SE), Prop 3.3–3.4 (one-sided noncompliance), three applications |
| Knaus (2024) | ✅ Read | PIVE framework (Sec 2), concrete outcome weights (Sec 3), normalization properties (Sec 4), **Appendix A.4** (kappa in PIVE — directly feeds Chapter 3) |

### 🔲 Tier 1 — Foundational (NOT-DONE)

| Paper | Status | Key takeaways |
|---|---|---|
| Chernozhukov et al. (2018) | ✅ Read | Section 2.3: Wald-AIPW formula, cross-fitting, PLR-IV |

### ✅ Tier 2 — Empirical application papers (DONE)

| Paper | Status | Key facts extracted |
|---|---|---|
| Angrist (1990) | ✅ Read | Z = draft lottery, D = veteran, Y = log wages. No always-takers → one-sided noncompliance. N = 3,027 from SIPP 1984 |
| Card (1993/1995) | ✅ Read | Z = proximity to 4-year college. D = some college (educ > 12) or completion (educ ≥ 16). Two covariate specs: Card full controls vs Kitagawa parsimonious |
| Angrist & Evans (1998) | ✅ Read | Z = same-sex siblings (+ twin second birth as robustness). D = third child. Y = LFP and log income. Near one-sided noncompliance — no always-takers |

### 🔲 Tier 3 — Cite-and-move-on (do not read deeply)

| Paper | Purpose | Status |
|---|---|---|
| Uysal (2011) | Origin of τ̂ᵤ. Formula fully in SUW 2025 | 🔲 Not needed |
| Frölich (2007) / Tan (2006) | Origin of τ̂ₜ = τ̂ₐ,₁. One paragraph in lit review | 🔲 Not needed |
| Heiler (2022) | CBPS for LATE. Skim Sec 2 only for CBPS definition | 🔲 Skim when writing Sec 2.6 |
| Blandhol et al. (2022) | "When is TSLS LATE?" Abstract + intro — motivation | 🔲 Not needed |
| Angrist & Pischke (2009) | General IV textbook reference. Chapter 4 if needed | 🔲 Not needed |

---

## Central Research Question

To what extent do machine learning-based and classical kappa weighting estimators of the LATE implicitly target different subpopulations, and what do their outcome weights reveal about covariate balance and estimator reliability across empirical applications?

### Sub-questions

**RQ1.** What are the structural properties of outcome weights for kappa-based LATE estimators (ESS, negative weight share, sum-to-zero as translation invariance criterion), and how do they differ from those of DML-based estimators independently of any specific dataset?

**RQ2.** When applied to empirical datasets, how do outcome weights of DML-based estimators (Wald-AIPW with cross-fitted random forests) and kappa estimators compare in terms of covariate balance (Love plots / SMDs), effective sample size, and negative weight patterns?

**RQ3.** Can outcome-weight diagnostics guide practitioners toward more robust estimator choices, and does the outcome-weights lens explain divergences between classical kappa, normalized kappa, and DML-based IV estimators?

---

## Thesis structure

### Chapter 1 — Introduction (3–4 pages)

**Section 1.1 — Motivation**
Hook: the same dataset yields wildly different LATE estimates depending on which estimator is used — not because of different assumptions, but because of how the outcome variable is coded. Unnormalized kappa estimators violate translation invariance: adding a constant to every outcome changes the treatment effect estimate. This is the practical failure mode that motivates the thesis.

**Section 1.2 — Background**
Standard practice uses 2SLS with additive covariates, which Blandhol et al. (2022) and Słoczyński (2021) show does not generally recover the LATE under heterogeneous effects. Kappa weighting estimators (Abadie 2003; SUW 2025) are a flexible alternative.

**Section 1.3 — Research gap**
Knaus (2024) introduces the PIVE framework and derives outcome weights ωᵢ such that τ̂ = Σᵢ ωᵢYᵢ for DML/GRF estimators, enabling covariate balance diagnostics via Love plots. Appendix A.4 of Knaus (2024) sketches the same derivation for kappa estimators but does not apply it empirically. This thesis fills that gap.

**Section 1.4 — Contribution**
1. Derive closed-form outcome weights for τ̂ᵤ and τ̂ₐ,₁₀ in the Knaus PIVE framework; show analytically why Σωᵢ = 0 iff translation invariant.
2. Clarify the distinction between Abadie's kappa weights (identification objects) and outcome weights in the PIVE sense (ωᵢ such that τ̂ = ΣωᵢYᵢ).
3. Apply Love plots and ESS diagnostics to kappa estimators for the first time, using the same pipeline as Knaus (2024).
4. Compare kappa estimators (τ̂ᵤᵐˡ, τ̂ᵤᶜᵇ, τ̂ₐ,₁₀) with DML Wald-AIPW across three empirical applications.
5. Discuss implications for the OutcomeWeights package: show how kappa outcome weights can be computed in the same format as `get_outcome_weights()`, enabling unified Love-plot diagnostics.

**Section 1.5 — Road map**
Brief chapter-by-chapter overview.

---

### Chapter 2 — Econometric Framework (8–10 pages)

**Section 2.1 — IV, LATE, and compliers**
- Potential outcomes notation: Yᵢ(0), Yᵢ(1), Dᵢ(0), Dᵢ(1)
- Four compliance types (AIR 1996): always-takers, never-takers, compliers, defiers
- IV assumptions (i)–(iv): conditional independence, exclusion restriction, first stage / overlap, monotonicity
- LATE definition: τᴸᴬᵀᴱ = E[Y₁ − Y₀ | D₁ > D₀]
- Why 2SLS may not recover LATE under heterogeneous effects (one-sentence reference to Blandhol et al. 2022)

**Section 2.2 — Abadie's kappa theorem**
- Lemma 2.1 (Abadie 2003, restated in SUW 2025 notation): the three weights κ, κ₁, κ₀ and their cell-by-cell values (Table 1 of SUW 2025)
- Parts (a), (b), (c) of the kappa theorem: any complier moment is identified
- Remark 2.2: E(κ) = E(κ₁) = E(κ₀) = P(D₁ > D₀) in population; why they diverge in finite samples

**Section 2.3 — Kappa-based LATE estimators**
- The five estimators: τ̂ᵤ (Uysal 2011), τ̂ₐ,₁₀ (Abadie & Cattaneo 2018), unnormalized τ̂ₐ, τ̂ₜ (= τ̂ₐ,₁, Frölich/Tan), τ̂ₐ,₀
- Normalized vs. unnormalized: what the distinction means mechanically

**Section 2.4 — Why normalization matters**
- Definition TI (translation invariance): τ̂(Y, W) = τ̂(Y+k, W) for all k
- Proposition 3.2 (SUW 2025): τ̂ᵤ and τ̂ₐ,₁₀ pass; τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ fail
- Definition SE (scale equivariance): brief statement, linked to log-unit sensitivity
- Concrete example: cents vs. dollars failure from Table 2 of SUW 2025

**Section 2.5 — One-sided noncompliance**
- Definition and examples (401k, draft lottery, twin births)
- Table 1 cell-by-cell signs of κ₁ and κ₀ under no-always-takers / no-never-takers
- Proposition 3.3 (SUW 2025): positive denominators guaranteed under one-sided noncompliance
- Proposition 3.4: τ̂ᵤ denominator positive in both one-sided cases; τ̂ₐ,₁₀ fails one case

**Section 2.6 — Estimation of the instrument propensity score**
- ML logit (τ̂ᵤᵐˡ) vs. covariate balancing CBPS (τ̂ᵤᶜᵇ)
- Proposition 3.5 (SUW 2025): with CBPS all normalized estimators coincide
- Why CBPS pushes weights away from extremes (Heiler 2022 argument, one paragraph)

**Section 2.7 — Double Machine Learning and Wald-AIPW**
- DML framework (Chernozhukov et al. 2018): PLR-IV model, nuisance parameters Ê[Y|Z,X] and Ê[D|Z,X], cross-fitting (K-fold), why it matters for valid inference
- The Wald-AIPW estimator: DML analogue of the Wald ratio, augmented with outcome and treatment regressions for efficiency and double robustness
- This is the ML benchmark against which kappa estimators are compared in Chapter 4
- Brief note on the OutcomeWeights R package (Knaus 2024): `dml_with_smoother()`, `get_outcome_weights()`, and the GitHub dev version's new DoubleML compatibility

---

### Chapter 3 — Connecting the Frameworks (5–7 pages)

*This is the thesis's unique theoretical contribution.*

**Section 3.1 — Kappa weights vs. outcome weights: clarifying the distinction**
- Kappa weights κᵢ: identification weights from Abadie (2003). They turn population expectations into complier expectations. They are not the same as outcome weights.
- Outcome weights ωᵢ (Knaus 2024): the scalar weights such that τ̂ = ΣᵢωᵢYᵢ exactly. Derived from kappa weights but a different object.
- The PIVE framework (Definition 1 of Knaus 2024): estimators solving Eₙ[(Ỹᵢ − τ̂D̃ᵢ)Z̃ᵢ] = 0. The two-step to outcome weights: (i) identify pseudo-instrument Z̃ and transformation matrix T; (ii) ω' = (Z̃'D̃)⁻¹Z̃'T.

**Section 3.2 — Analytical derivation of outcome weights for τ̂ᵤ and τ̂ₐ,₁₀**
- Express τ̂ᵤ (Equation 3 of SUW 2025) as Σᵢωᵢᵘ Yᵢ. Closed form:
  ωᵢᵘ = (1/D̂) · [Zᵢ/(Ŝ₁ p(Xᵢ)) − (1−Zᵢ)/(Ŝ₀(1−p(Xᵢ)))]
  where Ŝ₁ = (1/N)Σⱼ Zⱼ/p(Xⱼ), Ŝ₀ = (1/N)Σⱼ (1−Zⱼ)/(1−p(Xⱼ)), D̂ = estimated complier share
- Similarly for τ̂ₐ,₁₀: ωᵢᵃ'¹⁰ = κᵢ₁/Σⱼκⱼ₁ − κᵢ₀/Σⱼκⱼ₀
- Algebraic proof: Σᵢωᵢᵘ = 0 and Σᵢωᵢᵃ'¹⁰ = 0 ⟺ translation invariant
- Contrast: for unnormalized τ̂ₐ, Σᵢωᵢᵃ ≠ 0 in general (finite sample)
- Place in PIVE framework following Knaus (2024) Appendix A.4: identify Z̃, D̃, T for each kappa estimator

**Section 3.3 — Weight diagnostics: comparison across estimators**
- Sum-to-zero check (Σωᵢ) as empirical translation invariance diagnostic
- ESS = 1/Σωᵢ² comparison: τ̂ᵤ vs. τ̂ₐ,₁₀ vs. Wald-AIPW
- Negative weight share: always-takers and never-takers receive negative κ weights by construction (Table 1 of SUW 2025); does Wald-AIPW also assign negative ωᵢ, and to whom?
- Maximum absolute weight: which observations are most leveraged?
- Summary table of theoretical properties (to appear before any empirical application):

| Estimator | Σωᵢ = 0? | ESS | Neg. weights | Near-zero denom |
|---|---|---|---|---|
| τ̂ᵤ | ✓ exact | high | yes (AT + NT) | safe (one-sided) |
| τ̂ₐ,₁₀ | ✓ exact | high | yes | risk (one-sided) |
| τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ | ✗ finite sample | moderate | yes | risk |
| Wald-AIPW (DML) | ✓ approx. | moderate | yes | safe |

- Connection to Knaus (2024) Table 5 normalization classification: kappa normalized estimators are fully-normalized (Σωᵢ=0, Σωᵢ Dᵢ=1) by construction; Wald-AIPW is only scale-normalized in standard grf implementation unless C5b holds

---

### Chapter 4 — Empirical Applications (12–15 pages)

*For each application the structure is: (1) data and instrument; (2) replication of SUW 2025 table; (3) translation invariance check (cents/dollars/thousands); (4) outcome weight diagnostics (Σωᵢ, ESS, % negative, max weight); (5) Love plots for kappa estimators and Wald-AIPW; (6) interpretation.*

**Section 4.1 — Military service and wages (Angrist 1990)**
- Z = draft lottery eligibility, D = veteran status, Y = log wages. N = 3,027 (SIPP 1984)
- Covariate specs: linear age, cubic age, saturated age (three specifications)
- Replication: Table 2 of SUW 2025 — normalized estimators stable across cents/dollars; unnormalized flip sign
- New analysis: Love plots and ESS for the three age specifications. Does the saturated spec (where unnormalized = normalized) produce better covariate balance? How does Wald-AIPW compare?
- Note on one-sided noncompliance: no always-takers → Proposition 3.3 applies for τ̂ᵤ

**Section 4.2 — College education and wages (Card 1995)**
- Z = proximity to 4-year college, D = some college (educ > 12) and completion (educ ≥ 16), Y = log wages
- Two covariate specs: Card (1995) full controls; Kitagawa (2015) parsimonious
- Replication: Table 3 of SUW 2025 — large divergence of unnormalized estimates between specs; normalized more consistent
- New analysis: Love plots for both specs and both treatment definitions. Do weight diagnostics explain why estimates diverge between Card and Kitagawa? Does τ̂ᵤᶜᵇ outperform τ̂ᵤᵐˡ in covariate balance?

**Section 4.3 — Childbearing and labor supply (Angrist & Evans 1998)**
- Z = same-sex siblings, D = third child, Y = LFP and log income
- Near one-sided noncompliance (no always-takers) → Proposition 3.3 in practice
- Replication: Table 4 of SUW 2025 — most dramatic translation invariance failure; income estimates flip sign across cents/dollars/thousands
- New analysis: demonstrate Σκᵢ₁ > 0 by construction under no-always-takers; Wald-AIPW comparison; Love plots for both outcomes (LFP and log income)

**Section 4.4 — Cross-application comparison of outcome weight diagnostics**
*This section is the thesis's empirical synthesis — it does not appear in SUW 2025.*
- Tabulate ESS, % negative weights, Σωᵢ, and max weight across all three applications and all estimators
- Compare Love plots: which estimator achieves |SMD| ≤ 0.1 most reliably across covariates and applications?
- Do τ̂ᵤ and Wald-AIPW target the same subpopulation, or do their weight distributions look structurally different?
- Where do weight diagnostics reveal problems that point estimates alone do not?

---

### Chapter 5 — Discussion (3–4 pages)

**Section 5.1 — What outcome weights add**
- Point estimates alone do not reveal why estimators differ
- Outcome weights show which observations drive the estimate and whether the estimator targets the intended subpopulation
- Love plots make IV balance properties visible for the first time for kappa estimators

**Section 5.2 — DML vs. kappa: do they target the same compliers?**
- Weight distribution comparison: τ̂ᵤ vs. Wald-AIPW
- ESS comparison: kappa estimators typically have higher ESS (more observations contribute) vs. Wald-AIPW which concentrates weight more
- Sensitivity to covariate specification: Wald-AIPW adapts flexibly (RF); kappa estimators depend on propensity score specification (ML logit vs. CBPS)
- Practical guidance: when to prefer τ̂ᵤᶜᵇ vs. Wald-AIPW vs. 2SLS

**Section 5.3 — Implications for the OutcomeWeights package**
- The thesis shows that kappa outcome weights can be computed in the same format as `get_outcome_weights()` returns, enabling unified Love-plot diagnostics
- Concretely: the `kappa_outcome_weights()` function developed here returns a weight vector ωᵢ with the same structure as the `omega` matrix rows in OutcomeWeights — it can therefore be passed directly to `cobalt::love.plot()` using the same wrapper
- With Knaus's GitHub dev version now compatible with DoubleML, a natural next step would be a PR or companion vignette adding kappa estimators to the package workflow — this is flagged as an extension, not part of the thesis itself
- The translation invariance check (Σωᵢ = 0) and normalization classification from Knaus (2024) Table 5 are useful diagnostics that could be added as a `check_normalization()` utility

**Section 5.4 — Limitations**
- Bootstrap inference is computationally expensive
- Analytical standard errors not implemented for all estimators
- Love plots are descriptive, not inferential
- Weak overlap creates extreme weights in all estimators
- DML estimates depend on tuning choices; RF hyperparameter tuning affects balance (cf. Knaus 2024 Figure 3)

---

### Chapter 6 — Conclusion (1–2 pages)

- Summary of findings: which estimators are translation invariant, which achieve covariate balance, what the outcome-weights lens adds
- Recommendation: τ̂ᵤᶜᵇ preferred for robustness; Wald-AIPW for flexibility when large N and rich X; 2SLS defensible only with saturated covariate specification
- Extensions: doubly robust kappa estimators; formal tests of covariate balance; heterogeneous effects (translation invariance of causal forest estimates)

---

## Code task tracker

### ✅ Done

- [x] Read Angrist (1990): 2SLS + all five kappa estimators, cents vs. dollars — replicated Table 2 of SUW 2025
- [x] Replicate Card (1995): two treatments (somecol, educ16), two covariate specs (Card, Kitagawa), two outcomes — replicated Table 3 of SUW 2025
- [x] Replicate Angrist & Evans (1998): LFP + log income, all unit transformations — replicated Table 4 of SUW 2025
- [x] Bootstrap standard errors for all kappa estimators (R = 500 with set.seed(42))
- [x] `kappa_weights()` function: κ, κ₁, κ₀ from Abadie (2003) Lemma 2.1
- [x] `logit_mle()` propensity score estimation
- [x] `cbps()` covariate balancing propensity score (Newton step with backtracking)
- [x] `tau_u()`, `tau_a10()`, `tau_unnorm()`: all five kappa estimators
- [x] `kappa_outcome_weights()`: closed-form ωᵢ for all five kappa estimators with `stopifnot` verification that τ̂ = Σωᵢ Yᵢ
- [x] `weight_diag_table()`: Σωᵢ, ESS, % negative, max absolute weight
- [x] DML via `OutcomeWeights::dml_with_smoother()` with 5-fold cross-fitting for Angrist (1990) and Card (1995) applications
- [x] `get_outcome_weights()` extraction and verification (ω'Y = point estimates, TRUE)
- [x] Love plots via `cobalt::love.plot()` for DML estimators (Angrist, Card)
- [x] Translation invariance check: Y → Y + constant, verified normalized estimators stable and unnormalized estimators not
- [x] Comparison tables: kappa weight diagnostics vs. DML Wald-AIPW side by side

### 🟡 In progress

- [ ] Angrist & Evans (1998): DML Wald-AIPW with `dml_with_smoother()` — convergence issues under near one-sided noncompliance; verify AIPW-ATE NaN is expected
- [ ] Love plots for kappa estimators in Angrist & Evans (1998) application
- [ ] Modularize code: split the single Rmd into one script per application + shared functions file
- [ ] Check compatibility with OutcomeWeights GitHub dev version (`remotes::install_github("mknaus/OutcomeWeights")`) for DoubleML integration

### 🔲 TODO

- [ ] **Chapter 3 (theory):** write the analytical proof that Σωᵢᵘ = 0; write the PIVE representation for τ̂ᵤ following Knaus (2024) Appendix A.4
- [ ] **Section 4.4:** cross-application comparison table and discussion — tabulate ESS, % negative, Σωᵢ across all applications and estimators
- [ ] **DoubleML integration:** install GitHub dev version of OutcomeWeights; run DoubleML Wald-AIPW on all three applications; extract and compare outcome weights against grf-based Wald-AIPW
- [ ] **Love plots for kappa in Card (1995):** add cobalt love.plot calls using `kappa_outcome_weights()` output — currently only DML love plots exist for Card
- [ ] **Tuning sensitivity check (Section 4, inspired by Knaus 2024 Fig. 3):** compare `dml_with_smoother()` with default vs. `tune.parameters = "all"` — already coded for Angrist (1990) but not yet written up
- [ ] **Translation invariance for Wald-AIPW:** add the Y → Y+c check to DML estimators as a complement to the kappa translation invariance results
- [ ] **Section 5.3 code:** write a minimal `kappa_to_outcome_weights_format()` wrapper that returns output in the same `$omega` matrix structure as `get_outcome_weights()` — this unifies the pipeline and is the basis for the OutcomeWeights package discussion
- [ ] **Write-up:** Chapter 2 draft (framework sections — most reading is done)
- [ ] **Write-up:** Chapter 3 draft (theoretical derivations — `kappa_outcome_weights()` already implements this, needs to be written up formally)
- [ ] **Write-up:** Chapter 4 draft (empirical sections — replication done, narrative needed)

---

## Who introduced what — quick reference

| Concept | Source |
|---|---|
| Potential outcomes + LATE definition | Imbens & Angrist (1994) |
| Four subpopulations (compliers etc.) | Angrist, Imbens & Rubin (1996) |
| Kappa theorem (κ identifies any complier moment) | Abadie (2003) |
| τ̂ₜ estimator (ratio of two IPW estimators) | Frölich (2007) / Tan (2006) |
| τ̂ᵤ estimator (normalized, recommended) | Uysal (2011) |
| Translation invariance criterion + all five estimators compared | Słoczyński, Uysal & Wooldridge (2025) |
| One-sided noncompliance propositions (3.3, 3.4) | Słoczyński, Uysal & Wooldridge (2025) |
| CBPS for LATE (Proposition 3.5) | Słoczyński et al. (2025), after Heiler (2022) |
| 2SLS doesn't recover LATE | Blandhol et al. (2022) |
| Wald-AIPW / DML framework | Chernozhukov et al. (2018) |
| Outcome weights framework for ML estimators (PIVE) | Knaus (2024) |
| Kappa estimators in PIVE framework | Knaus (2024) Appendix A.4 |
| Normalization properties and Table 5 classification | Knaus (2024) Sections 4.2–4.3 |
| OutcomeWeights R package | Knaus (2024) / GitHub dev version (DoubleML compatible) |

---

## Page budget

| Chapter | Target pages |
|---|---|
| 1. Introduction | 3–4 |
| 2. Econometric Framework | 8–10 |
| 3. Connecting the Frameworks | 5–7 |
| 4. Empirical Applications | 12–15 |
| 5. Discussion | 3–4 |
| 6. Conclusion | 1–2 |
| References | 2–3 |
| **Total** | **34–45** |
