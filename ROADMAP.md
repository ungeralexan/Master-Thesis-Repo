# Thesis Roadmap

**Last updated:** May 23, 2026  
**Phase:** Vietnam + Card fully coded → Writing active | Professor meeting recommended

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

### ✅ Tier 1 — Foundational (DONE)

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

## Thesis structure (Gliederung)

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
4. Compare kappa estimators (τ̂ᵤᵐˡ, τ̂ᵤᶜᵇ, τ̂ₐ,₁₀) with DML Wald-AIPW across three empirical applications, using multiple ML learners for the nuisance parameters.
5. **NEW — Design dominates learner finding:** document empirically that in low-dimensional near-random-assignment designs, ESS and weight structure are identical across all estimators and ML learners — normalization, not covariate adjustment, drives estimate divergence.
6. Discuss implications for the OutcomeWeights package: show how kappa outcome weights can be computed in the same format as `get_outcome_weights()`, enabling unified Love-plot diagnostics.

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
- Propensity score estimation: MLE logit (τ̂ᵤᵐˡ) vs. CBPS (τ̂ᵤᶜᵇ); Proposition 3.5: with CBPS all normalized estimators coincide

**Section 2.4 — Why normalization matters**
- Definition TI (translation invariance): τ̂(Y, W) = τ̂(Y+k, W) for all k
- Proposition 3.2 (SUW 2025): τ̂ᵤ and τ̂ₐ,₁₀ pass; τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ fail
- Definition SE (scale equivariance): brief statement, linked to log-unit sensitivity
- Concrete example: cents vs. dollars failure from Table 2 of SUW 2025
- One-sided noncompliance (Propositions 3.3, 3.4 of SUW 2025): near-zero denominators; why τ̂ᵤ is safe but τ̂ₐ,₁₀ is not in the κ₁-only case

**Section 2.5 — DML and Wald-AIPW**
- DML framework (Chernozhukov et al. 2018): PLR-IV and Wald-AIPW estimator
- Two nuisance parameters: E[Y|Z, X] and E[D|Z, X], estimated via K-fold cross-fitting
- Brief description of cross-fitting: why it matters for valid inference

**Section 2.6 — PIVE framework and outcome weights (Knaus 2024)**
- Any estimator fitting the pseudo-IV structure: τ̂ = Σᵢ ωᵢYᵢ
- The two-step: (i) form pseudo-instrument Z̃ and transformation matrix T; (ii) ω' = (Z̃'D̃)⁻¹Z̃'T
- "Fully normalized" in Knaus (2024): Σ_{D=1} ωᵢ = +1, Σ_{D=0} ωᵢ = −1 (Table 5)
- ML learners for DML: Wald-AIPW via grf (OutcomeWeights::dml_with_smoother) and DoubleMLIIVM (ranger, XGBoost, linear+logistic)

**Section 2.7 — Covariate balance diagnostics**
- Standardized Mean Difference (SMD): |X̄ₜᵣₑₐₜₑ_ₖ − X̄_cₒₙₜᵣₒₗ_ₖ| / SD(Xₖ), computed with outcome weights ωᵢ
- Love plots: one dot per covariate, unadjusted vs. weighted SMD; threshold at |SMD| ≤ 0.1
- Effective Sample Size (ESS): ESS = 1 / Σᵢ ωᵢ², measures how many observations effectively contribute
- Negative weight share: % of observations with ωᵢ < 0
- Connection to translation invariance: Σᵢ ωᵢ = 0 ⟺ translation invariant (the sum-to-zero condition)

---

### Chapter 3 — Connecting the Frameworks (5–7 pages)

*This is the thesis's unique theoretical contribution. Knaus (2024) Appendix A.4 derives the PIVE representation for τ̂ₐ, τ̂ₐ,₀, τ̂ₐ,₁₀ but does not explicitly derive ωᵢ for τ̂ᵤ, nor does it apply Love plots or ESS diagnostics to kappa estimators empirically. This chapter bridges that gap.*

**Section 3.1 — Kappa weights vs. outcome weights: clarifying the distinction**
- Kappa weights κᵢ: identification weights from Abadie (2003). Population objects that turn expectations into complier expectations. Not directly τ̂ = ΣωᵢYᵢ.
- Outcome weights ωᵢ (Knaus 2024): sample-level weights such that τ̂ = Σᵢ ωᵢYᵢ. Derived from κ but not the same object.
- This distinction matters practically: Love plots use ωᵢ, not κᵢ.

**Section 3.2 — Analytical derivation of outcome weights for τ̂ᵤ and τ̂ₐ,₁₀**
- Express τ̂ᵤ as weighted sum Σᵢ ωᵢᵘ Yᵢ; closed-form ωᵢᵘ in terms of p(Xᵢ), Zᵢ, and sample normalizations Ŝ₁, Ŝ₀, D̂
- Similarly for τ̂ₐ,₁₀: weights ωᵢᵃ'¹⁰ = κᵢ₁/Σκⱼ₁ − κᵢ₀/Σκⱼ₀
- Algebraic proof: Σᵢ ωᵢᵘ = 0 and Σᵢ ωᵢᵃ'¹⁰ = 0 (translation invariant); for unnormalized τ̂ₐ, Σᵢ ωᵢᵃ ≠ 0 in general
- Place in PIVE framework following Knaus (2024) Appendix A.4: identify Z̃, D̃, T for each kappa estimator
- `kappa_outcome_weights()` function implements this; verified τ̂ = Σωᵢ Yᵢ up to machine precision for all five estimators

**Section 3.3 — Theoretical weight properties: a summary table**

| Estimator | Σωᵢ = 0? | ESS | Negative weights | Near-zero denom. risk |
|---|---|---|---|---|
| τ̂ᵤᶜᵇ (CBPS, normalized) | ✓ exact | high | yes (AT+NT) | safe (one-sided) |
| τ̂ᵤᵐˡ (MLE, normalized) | ✓ exact | high | yes (AT+NT) | safe (one-sided) |
| τ̂ₐ,₁₀ (normalized) | ✓ exact | high | yes | risk (one-sided) |
| τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ (unnormalized) | ✗ finite sample | moderate | yes | risk |
| Wald-AIPW (DML, grf) | ✓ approx. | moderate | yes | safe |
| Wald-AIPW (DoubleML) | ✓ approx. | moderate | yes | safe |

---

### Chapter 4 — Empirical Application: Military Service and Wages (12–15 pages)

*For each estimator group and specification, the structure is: (1) data and setup; (2) point estimates and replication; (3) translation invariance check; (4) outcome weight diagnostics (Σωᵢ, ESS, % negative, max |ω|); (5) Love plots; (6) interpretation.*

**Section 4.1 — Data and design (Angrist 1990)**
- Setting: SIPP 1984, N = 3,027 white men born 1950–1953
- Instrument Z: Vietnam-era draft lottery eligibility (rsncode, binary)
- Treatment D: veteran status (nvstat, binary)
- Outcome Y: log wages (dollars and cents — both coded to demonstrate translation invariance)
- Key feature: no always-takers (one-sided noncompliance) → Proposition 3.3 of SUW 2025 applies; τ̂ₐ,₁₀ has near-zero denominator risk
- Only covariate: age (the one variable needed for conditional independence of the lottery)

**Section 4.2 — Kappa estimators: replication and covariate specifications**

Three covariate specifications are compared:
- **Spec 1 — Linear age:** logit p-score on age (one-dimensional linear control). Replicates SUW 2025 Table 2, columns 1–2. *Note: dropped from DML comparison — a single continuous predictor is uninformative for flexible ML smoothers and showed unstable behavior with grf; kappa estimates still reported.*
- **Spec 2 — Cubic age:** logit p-score on age + age² + age³. Replicates SUW 2025 Table 2, columns 3–4. **Main specification for DML comparison.**
- **Spec 3 — Saturated age:** full set of age dummies (one per age value). Replicates SUW 2025 Table 2, columns 5–6. Also used in DML comparison.

*Key finding:* Normalized estimators (τ̂ᵤ, τ̂ₐ,₁₀) are stable across cents/dollars in all specs. Unnormalized estimators (τ̂ₐ, τ̂ₜ, τ̂ₐ,₀) flip sign dramatically: ~+0.5 in cents vs. ~+0.3 in dollars for cubic spec, and large negative values in linear spec. Saturated spec is special: all estimators — including unnormalized — agree, because the fully nonparametric propensity score forces Σωᵢ ≈ 0 even for unnormalized estimators.

*Estimates (cubic spec, dollars):*
- 2SLS: 0.243; τ̂ᵤᶜᵇ: 0.210; τ̂ᵤᵐˡ: 0.202; τ̂ₐ,₁₀: 0.204
- All normalized estimators cluster around 0.20–0.24; unnormalized range from 0.30–0.32 (dollars) to 0.52–0.54 (cents)

**Section 4.3 — DML Wald-AIPW: grf (OutcomeWeights) vs. DoubleMLIIVM**

Two DML frameworks are compared on cubic and saturated specs:

*Framework 1 — grf via OutcomeWeights::dml_with_smoother() (5-fold cross-fitting):*
- PLR-IV (cubic): 0.243; Wald-AIPW (cubic): 0.229
- PLR-IV (saturated): 0.242; Wald-AIPW (saturated): 0.244
- Translation invariance: verified exactly (Σωᵢ = 0 for both)

*Framework 2 — DoubleMLIIVM with three learners (cubic spec):*
- Linear + logistic (parametric baseline): Wald-AIPW ≈ 0.218
- Ranger (random forest): Wald-AIPW ≈ 0.246
- XGBoost: Wald-AIPW ≈ 0.246
- *Key result:* point estimates converge (0.218–0.246) across all learners; ESS and % negative weights are identical (ESS = 5, 54.4% negative) — learner choice is immaterial in this low-dimensional near-random-assignment design
- *XGBoost flag:* algebraic identity Σωᵢ Yᵢ = τ̂ fails marginally (Sum_w = −3.73e−06), consistent with Knaus (2024) Table 6 — XGBoost violates Condition 3 (non-affine smoother); documented for completeness, does not affect point estimate meaningfully

**Section 4.4 — Weight diagnostics and comparison**

Summary table across all estimators (cubic spec, dollar outcome):

| Estimator | Estimate | Σωᵢ | ESS | % neg. | max|ω| |
|---|---|---|---|---|---|
| τ̂ᵤᶜᵇ (kappa, CBPS) | 0.210 | 0 | 5 | 54.4 | 0.025 |
| τ̂ᵤᵐˡ (kappa, MLE) | 0.202 | 0 | 5 | 54.4 | 0.028 |
| τ̂ₐ,₁₀ (kappa, MLE) | 0.204 | 0 | 5 | 54.4 | 0.029 |
| τ̂ₐ (kappa, unnorm.) | 0.314 | 0.048 | 5 | 54.4 | 0.029 |
| Wald-AIPW (grf, cubic) | 0.229 | 0 | 5 | 54.4 | 0.018 |
| Wald-AIPW (DoubleML, ranger) | 0.246 | 0 | 5 | 54.4 | 0.022 |
| Wald-AIPW (DoubleML, XGBoost) | 0.246 | ~0 | 5 | 54.4 | 0.022 |

*Interpretation:* ESS and % negative weights are **identical** across all estimators and all ML learners. In this near-random-assignment, low-dimensional setting, the instrument almost perfectly randomizes treatment — all estimators converge to the same weighting structure. The **design dominates learner choice**: ESS uniformity is not a failure of the diagnostics but a correct reflection of the data-generating process. Estimate divergence is driven entirely by normalization, not by covariate adjustment. Love plots confirm near-perfect balance for all estimators.

**Section 4.5 — Love plots**
- Cubic and saturated specs: Love plots for kappa estimators and DML Wald-AIPW side-by-side
- All estimators achieve near-perfect balance (age is the only covariate; draft lottery is near-random)
- Diagnostic value: confirms the weight pipeline is correctly implemented before applying to richer designs

**Section 4.6 — Status**
- ✅ Point estimates: full replication (kappa all specs, DML cubic + saturated, DoubleML three learners)
- ✅ Translation invariance: verified for kappa and DML (cents vs. dollars)
- ✅ Weight diagnostics table: Σωᵢ, ESS, % negative, max|ω| for all estimators
- ✅ Love plots: kappa + DML (cubic + saturated) + DoubleML learner comparison
- ✅ DoubleML learner comparison (ranger, XGBoost, linear+logistic): all coded and results documented
- ✅ **New insight documented:** ESS uniformity across learners = design dominates learner, not a bug
- 🟡 Write-up of tuning sensitivity: compare `dml_with_smoother()` default vs. `tune.parameters="all"` — coded, not yet written
- 🟡 XGBoost non-affine issue: flag and cite Knaus (2024) Table 6 formally in text

---

### Chapter 5 — Empirical Applications: Card (1995) and Angrist & Evans (1998) (6–8 pages)

**Section 5.1 — Card (1995): College education and wages**
- Z = proximity to four-year college; D = some college (educ>12) or completion (educ≥16); Y = log wages
- Two covariate specs: Card (1995) full controls vs. Kitagawa (2015) parsimonious
- Key finding: large divergence in unnormalized estimates across specs; normalized more consistent. Kitagawa ESS = 1 (extreme weight concentration) vs. Card spec ESS higher — covariate spec matters here unlike Vietnam
- DML comparison done: DoubleML with three learners (linear+logistic, ranger, XGBoost)
- Love plots done: Card spec vs. Kitagawa spec, both treatments (somecol D1, educ16 D2), 8-panel grids
- **Status:** ✅ All coding complete (kappa + DML + love plots + weight diagnostics); 🟡 write-up pending

*Key results Card (1995):*
- Kitagawa spec: ESS ≈ 1, very high % negative weights — extreme weight concentration, raises reliability concerns
- Card spec: ESS somewhat higher, better-spread weights
- Love plots: Kitagawa spec shows more imbalance than Card spec → covariate richness matters for complier balance
- **Contrast with Vietnam:** here, unlike Vietnam, learner and spec choice *do* affect weight diagnostics — richer design reveals diagnostic value of the outcome-weights lens

**Section 5.2 — Angrist & Evans (1998): Childbearing and labor supply**
- Z = same-sex siblings; D = third child; Y = LFP and log income
- Near one-sided noncompliance (no always-takers) → Proposition 3.3 applies
- Key finding to replicate: most dramatic translation invariance failure; income estimates flip sign across units
- DML comparison: Wald-AIPW only (AIPW-ATE gives NaN under near one-sided noncompliance — documented as expected)
- **Status:** ✅ kappa done + love plots done; 🟡 DML Wald-AIPW coding in progress

---

### Chapter 6 — Discussion (3–4 pages)

**Section 6.1 — Cross-application summary**
- Comparative table: ESS, % negative, Σωᵢ, Love plot quality across all three applications and all estimators
- Key contrast: Vietnam (near-random, low-dim) → ESS uniform, design dominates; Card (observational, rich X) → ESS diverges, diagnostics informative; Angrist & Evans (one-sided noncompliance) → τ̂ₐ,₁₀ near-zero denominator risk confirmed
- Does τ̂ᵤᶜᵇ (CBPS) consistently outperform τ̂ᵤᵐˡ (MLE)?

**Section 6.2 — The outcome weights lens: what it adds**
- What do Love plots reveal beyond point estimates alone?
- When do weight diagnostics signal problems the estimates don't? (Answer: Card Kitagawa ESS ≈ 1 is a red flag even if point estimate looks reasonable)
- DML vs. kappa: are Wald-AIPW and τ̂ᵤ targeting the same complier subpopulation?

**Section 6.3 — DML learner comparison: implications**
- Point estimates from ranger, XGBoost, and linear+logistic converge for Vietnam (0.218–0.246) and Card
- **Design dominates learner:** in low-dimensional near-random-assignment settings, learner flexibility doesn't add value — ESS and % negative weights are identical across learners; this is the expected outcome, not a methodological failure
- XGBoost's non-affine smoother issue: recommend always verifying Σωᵢ Yᵢ = τ̂ numerically
- Where learner matters: richer observational designs (Card) may show more divergence — worth flagging as open question

**Section 6.4 — Practical guidance for practitioners**
- Decision framework: when to use τ̂ᵤ vs. τ̂ᵤᶜᵇ vs. Wald-AIPW; when is 2SLS still defensible (saturated covariate spec)
- Always use normalized estimators; run the sum-to-zero check (Σωᵢ = 0)
- Check ESS: ESS ≈ 1 is a reliability red flag even if the point estimate looks plausible
- For DML: use ranger or grf; XGBoost needs additional algebraic verification
- Package note: `kappa_to_outcome_weights_format()` wrapper enables unified Love-plot pipeline

**Section 6.5 — Limitations**
- Bootstrap inference is computationally expensive
- Analytical standard errors not implemented for all estimators
- Love plots are descriptive, not inferential
- Weak overlap creates extreme weights in all estimators
- DML estimates depend on tuning choices; RF hyperparameter tuning affects balance (cf. Knaus 2024 Figure 3)

---

### Extension — Heterogeneous Treatment Effects via Instrumental Forest

**Condition:** only pursued if ATE results motivate it (professor's gate).
Motivation exists if Wald-AIPW and κ-based estimators diverge meaningfully
in point estimates or weight diagnostics across covariate subgroups.

**Current assessment:** Vietnam ATE results show tight convergence — weak motivation. Card (1995) has richer covariates and Kitagawa vs. Card spec divergence in ESS — stronger motivation. Make this decision after professor meeting.

**Status:** 🔲 Conditional on Card (1995) ATE write-up and professor feedback

---

### Chapter 7 — Conclusion (1–2 pages)

- Summary of findings: which estimators are translation invariant, which achieve covariate balance, what the outcome-weights lens adds
- **New:** summarize the "design dominates learner" finding as a cross-application empirical result
- Recommendation: τ̂ᵤᶜᵇ preferred for robustness; Wald-AIPW for flexibility when large N and rich X; 2SLS defensible only with saturated covariate specification; check ESS before trusting any point estimate
- Extensions: doubly robust kappa estimators; formal tests of covariate balance; heterogeneous effects; `kappa_to_outcome_weights_format()` as package contribution

---

## Code task tracker

### ✅ Done

- [x] Read Angrist (1990): 2SLS + all five kappa estimators, cents vs. dollars — replicated Table 2 of SUW 2025
- [x] Replicate Card (1995): two treatments (somecol, educ16), two covariate specs (Card, Kitagawa), two outcomes — replicated Table 3 of SUW 2025
- [x] Replicate Angrist & Evans (1998): LFP + log income, all unit transformations — replicated Table 4 of SUW 2025
- [x] `kappa_weights()` function: κ, κ₁, κ₀ from Abadie (2003) Lemma 2.1
- [x] `logit_mle()` propensity score estimation
- [x] `cbps()` covariate balancing propensity score (Newton step with backtracking)
- [x] `tau_u()`, `tau_a10()`, `tau_unnorm()`: all five kappa estimators
- [x] `kappa_outcome_weights()`: closed-form ωᵢ for all five kappa estimators with `stopifnot` verification that τ̂ = Σωᵢ Yᵢ
- [x] Chapter 3 (theory): write the analytical proof that Σωᵢᵘ = 0; write the PIVE representation for τ̂ᵤ following Knaus (2024) Appendix A.4
- [x] `weight_diag_table()`: Σωᵢ, ESS, % negative, max absolute weight
- [x] Analytical M-estimation standard errors for all kappa estimators (following SUW 2025 online appendix)
- [x] DML via `OutcomeWeights::dml_with_smoother()` (grf) with 5-fold cross-fitting for Angrist (1990) cubic + saturated specs
- [x] `get_outcome_weights()` extraction and verification (ω'Y = point estimates, TRUE) for grf objects
- [x] Love plots via `cobalt::love.plot()` for grf DML estimators — Angrist (1990) cubic + saturated specs
- [x] Love plots for kappa estimators in Angrist (1990)
- [x] Translation invariance check: Y → Y + constant, verified normalized estimators stable and unnormalized not — for all three covariate specs and both kappa and DML estimators; saturated spec special case documented
- [x] Comparison tables: kappa weight diagnostics vs. DML Wald-AIPW side by side (Σωᵢ, ESS, % negative, max|ω|) — Vietnam
- [x] **DoubleML learner comparison — Angrist (1990), cubic spec:** ranger, XGBoost, linear+logistic; weight diagnostics + love plots; algebraic check; comprehensive summary table
- [x] **DoubleML learner comparison — Card (1995):** same three-learner pipeline; all four spec × treatment combinations (Card/Kitagawa × somecol/educ16)
- [x] Love plots for Card (1995): 8-panel grids for all four spec × treatment combinations
- [x] Weight diagnostics table for Card (1995): Kitagawa ESS ≈ 1 documented; Card spec comparison
- [x] Love plots for kappa estimators in Angrist & Evans (1998)
- [x] Translation invariance check for Card (1995) and Vietnam DML
- [x] Created reference file showing all functions used

### 🟡 In progress

- [ ] Angrist & Evans (1998): DML Wald-AIPW with `dml_with_smoother()` — AIPW-ATE NaN confirmed expected under near one-sided noncompliance; only Wald-AIPW reported; coding in progress

### 🔲 TODO — Code

- [ ] **Angrist & Evans (1998):** complete DML Wald-AIPW coding; weight diagnostics table; love plots for DML estimators

### 🔲 TODO — Writing (priority order)

1. [ ] **NEW: Write the "design dominates learner" paragraph** — explain *why* ESS = 5 and 54.4% negative weights across all Vietnam estimators is the *expected* result for near-random-assignment low-dimensional designs, not a diagnostic failure. This reframes the finding for the reader before the Chapter 4 results are presented. ~1 paragraph, fits in Section 4.3 or 4.4 interpretation.

2. [ ] **Section 5.1 (Card) write-up** — point estimates, translation invariance, weight diagnostics, love plots. Key narrative: contrast with Vietnam (ESS uniform there, ESS diverges here; Kitagawa ESS ≈ 1 is a reliability flag). This is the most important remaining writing task.

3. [ ] **DML learner write-up (Section 4.3 / 6.3)** — document convergence of estimates across learners; identical ESS and % negative weights; XGBoost non-affine smoother issue; cite Knaus (2024) Table 6 + Figure 1

4. [ ] **Chapter 2 draft** — framework sections; most reading is done; should be relatively fast

5. [ ] **Chapter 3 draft** — theoretical derivations; `kappa_outcome_weights()` already implements this; needs formal write-up

6. [ ] **Chapter 4 draft** — Angrist (1990) in full; fully coded; write-up should be straightforward

7. [ ] **Tuning sensitivity write-up (Section 4.2)** — compare `dml_with_smoother()` default vs. `tune.parameters="all"` — already coded, not yet written

8. [ ] **Section 5.2 (Angrist & Evans) write-up** — after DML coding complete

9. [ ] **Cross-application comparison table (Section 6.1)** — ESS, % negative, Σωᵢ across all applications; needs Angrist & Evans DML to be complete

10. [ ] **Chapter 6 Discussion draft** — builds on all empirical chapters; do last

### 🔲 TODO — Before professor meeting

- [ ] Prepare 1-page summary of key empirical findings across Vietnam and Card: ESS table, Love plot highlights, "design dominates learner" finding stated clearly
- [ ] Decide with professor: is Angrist & Evans DML worth finishing, or can it be described as "kappa only, DML pending" in Chapter 5?
- [ ] Decide with professor: pursue heterogeneous effects extension (Card instrumental forest) or stay with ATE?
- [ ] Clarify thesis deadline and remaining time budget → prioritize writing vs. remaining coding

### 💡 Later / extension ideas

- [ ] **Package contribution — `kappa_to_outcome_weights_format()` wrapper:** makes kappa outcome weights directly passable to `cobalt::love.plot()` and OutcomeWeights-compatible functions. Could be proposed to Knaus as a PR or companion vignette
- [ ] **Package contribution — `check_normalization()` utility:** takes any ω vector and classifies it as fully-normalized / scale-normalized / unnormalized following Knaus (2024) Table 4
- [ ] **Package contribution — kappa vignette:** worked example showing how to compute kappa outcome weights and pass them alongside DML weights into the same Love plot pipeline

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
| **Design dominates learner (ESS uniformity in low-dim near-random designs)** | **This thesis — empirical finding, Vietnam + Card** |

---

## Page budget

| Chapter | Target pages |
|---|---|
| 1. Introduction | 3–4 |
| 2. Econometric Framework | 8–10 |
| 3. Connecting the Frameworks | 5–7 |
| 4. Empirical Application: Angrist (1990) | 7–9 |
| 5. Empirical Applications: Card + Angrist & Evans | 6–8 |
| 6. Discussion | 3–4 |
| 7. Conclusion | 1–2 |
| References | 2–3 |
| **Total** | **35–47** |
