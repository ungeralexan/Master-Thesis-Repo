## Central Research Question

To what extent do machine learning-based and classical kappa weighting estimators
of the LATE implicitly target different subpopulations, and what do their outcome
weights reveal about covariate balance and estimator reliability across empirical
applications?

### Sub-questions

**RQ1.** What are the structural properties of outcome weights for kappa-based LATE
estimators (ESS, negative weight share, sum-to-zero as translation invariance
criterion), and how do they differ from those of DML-based estimators independently
of any specific dataset?

**RQ2.** When applied to empirical datasets, how do outcome weights of DML-based
estimators (Wald-AIPW with cross-fitted random forests) and kappa estimators compare
in terms of covariate balance (Love plots / SMDs), effective sample size, and
negative weight patterns?

**RQ3.** Can outcome-weight diagnostics guide practitioners toward more robust
estimator choices, and does the outcome-weights lens explain divergences between
classical kappa, normalized kappa, and DML-based IV estimators?

---

## Thesis structure (Gliederung)

### Chapter 1 — Introduction (3–4 pages) 🔲

**Section 1.1 — Motivation**
Hook: the same dataset yields wildly different LATE estimates depending on which
estimator is used — not because of different assumptions, but because of how the
outcome variable is coded. Unnormalized kappa estimators violate translation
invariance: adding a constant to every outcome changes the treatment effect estimate.
This is the practical failure mode that motivates the thesis.

**Section 1.2 — Background**
Standard practice uses 2SLS with additive covariates, which Blandhol et al. (2022)
and Słoczyński (2021) show does not generally recover the LATE under heterogeneous
effects. Kappa weighting estimators (Abadie 2003; SUW 2025) are a flexible
alternative.

**Section 1.3 — Research gap**
Knaus (2024) introduces the PIVE framework and derives outcome weights ωᵢ such that
τ̂ = Σᵢ ωᵢYᵢ for DML/GRF estimators, enabling covariate balance diagnostics via
Love plots. Appendix A.4 of Knaus (2024) sketches the same derivation for kappa
estimators but does not apply it empirically. This thesis fills that gap.

**Section 1.4 — Contribution**
1. Derive closed-form outcome weights for τ̂ᵤ in the Knaus PIVE framework; show
   analytically why Σωᵢ = 0 iff translation invariant, check this also for all the applications in the code, noticing that there are genuinly to methods how to proceed there. 
2. Clarify the distinction between Abadie's kappa weights (identification objects)
   and outcome weights in the PIVE sense (ωᵢ such that τ̂ = ΣωᵢYᵢ).
3. Apply Love plots and ESS diagnostics to kappa estimators for the first time,
   using the same pipeline as Knaus (2024).
4. Compare kappa estimators (τ̂ᵤᵐˡ, τ̂ᵤᶜᵇ, τ̂ₐ,₁₀) with DML Wald-AIPW across
   three empirical applications, using multiple ML learners for the nuisance
   parameters.
5. Apply the DML estimators using Knaus outcome weights package. And make a full diagnostics section meaning really look sometimes at the smoother matrix and maybe look also at the descritptives or so.  
6. Implement the Wald AIPW estimation with double machine learning using XGBOOST and linear regression and the ranger function of the package and try there to set up a tuning section which tunes the parameters for the xgboost in order to get clean results for the DMl estimator section.
7. New task from Knaus whihc he summarised as 2. Eine sehr spannende Erweiterung wäre noch Wald-AIPW und PLR-IV mit logit pscore und OLS outcome regression zu implementieren, quasi als Zwischenstufe zwischen den parametrischen Methoden im Wooldridge Paper und den DML Methoden. Das ist nicht im OutcomeWeights Paket abgedeckt und wäre gerade deswegen ideal, um auch ein bisschen zu Coden. Am einfachsten wäre es vermutlich, das über das DoubleML package zu machen. Ich hänge dir mal an, wie man es ohne Instrument innerhalb von DoubeML implementieren würde. Ich denke es würde dann darauf hinauslaufen diese Funktion https://github.com/MCKnaus/OutcomeWeights/blob/0c94f940b04c14d0247b46842af37752e306b79e/R/DoubleML.R#L191 kompatibel zu machen mit lrn("regr.lm"). Wenn du das sauber hinbekommst (ich würde dir bei Problemen auch helfen) hättest du deutlich mehr erreicht, als wenn du die instrumental_forest Funktionen einfach anwendest und es würde sich exzellent in die übergeordnete Fragestellung einfügen. Was meinst du?



**Section 1.5 — Road map**
Brief chapter-by-chapter overview. The three applications are chosen to span the
structural variation theory predicts should matter: near-random assignment with
low-dimensional covariates (Angrist 1990), a richer observational design with
multiple covariate specifications (Card 1995), and near one-sided noncompliance
where unnormalized estimators face their sharpest failure mode (Angrist & Evans
1998).

---

### Chapter 2 — Econometric Framework (8–10 pages) 🔲

**Section 2.1 — IV, LATE, and compliers**
- Potential outcomes notation: Yᵢ(0), Yᵢ(1), Dᵢ(0), Dᵢ(1)
- Four compliance types (AIR 1996): always-takers, never-takers, compliers, defiers
- IV assumptions (i)–(iv): conditional independence, exclusion restriction,
  first stage / overlap, monotonicity
- LATE definition: τᴸᴬᵀᴱ = E[Y₁ − Y₀ | D₁ > D₀]
- Why 2SLS may not recover LATE under heterogeneous effects (one-sentence reference
  to Blandhol et al. 2022)

**Section 2.2 — Abadie's kappa theorem**
- Lemma 2.1 (Abadie 2003, restated in SUW 2025 notation): the three weights
  κ, κ₁, κ₀ and their cell-by-cell values (Table 1 of SUW 2025)
- Parts (a), (b), (c) of the kappa theorem: any complier moment is identified
- Remark 2.2: E(κ) = E(κ₁) = E(κ₀) = P(D₁ > D₀) in population; why they
  diverge in finite samples

**Section 2.3 — Kappa-based LATE estimators**
- The five estimators: τ̂ᵤ (Uysal 2011), τ̂ₐ,₁₀ (Abadie & Cattaneo 2018),
  unnormalized τ̂ₐ, τ̂ₜ (= τ̂ₐ,₁, Frölich/Tan), τ̂ₐ,₀
- Normalized vs. unnormalized: what the distinction means mechanically
- Propensity score estimation: MLE logit (τ̂ᵤᵐˡ) vs. CBPS (τ̂ᵤᶜᵇ);
  Proposition 3.5: with CBPS all normalized estimators coincide

**Section 2.4 — Why normalization matters**
- Definition TI (translation invariance): τ̂(Y, W) = τ̂(Y+k, W) for all k
- Proposition 3.2 (SUW 2025): τ̂ᵤ and τ̂ₐ,₁₀ pass; τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ fail
- Definition SE (scale equivariance): brief statement, linked to log-unit sensitivity
- Concrete example: cents vs. dollars failure from Table 2 of SUW 2025
- One-sided noncompliance (Propositions 3.3, 3.4 of SUW 2025): near-zero
  denominators; why τ̂ᵤ is safe but τ̂ₐ,₁₀ is not in the κ₁-only case

**Section 2.5 — DML and Wald-AIPW**
- DML framework (Chernozhukov et al. 2018): PLR-IV and Wald-AIPW estimator
- The difference between LATE and the constant structural treatment theta and when
  they are the same
- Two nuisance parameters: E[Y|Z, X] and E[D|Z, X], estimated via K-fold
  cross-fitting
- Brief description of cross-fitting: why it matters for valid inference

**Section 2.6 — PIVE framework and outcome weights (Knaus 2024)**
- Any estimator fitting the pseudo-IV structure: τ̂ = Σᵢ ωᵢYᵢ
- The two-step: (i) form pseudo-instrument Z̃ and transformation matrix T;
  (ii) ω' = (Z̃'D̃)⁻¹Z̃'T
- "Fully normalized" in Knaus (2024): Σ_{D=1} ωᵢ = +1, Σ_{D=0} ωᵢ = −1
  (Table 5)
- ML learners for DML: Wald-AIPW via grf (OutcomeWeights::dml_with_smoother)
  and DoubleMLIIVM (ranger, XGBoost, linear+logistic)

**Section 2.7 — Covariate balance diagnostics**
- Standardized Mean Difference (SMD): |X̄ₜᵣₑₐₜₑ_ₖ − X̄_cₒₙₜᵣₒₗ_ₖ| / SD(Xₖ),
  computed with outcome weights ωᵢ
- Love plots: one dot per covariate, unadjusted vs. weighted SMD; threshold at
  |SMD| ≤ 0.1
- Effective Sample Size (ESS): ESS = 1 / Σᵢ ωᵢ², measures how many observations
  effectively contribute
- Negative weight share: % of observations with ωᵢ < 0
- Connection to translation invariance: Σᵢ ωᵢ = 0 ⟺ translation invariant
  (the sum-to-zero condition)

---

### Chapter 3 — Connecting the Frameworks (5–7 pages) ✅ COMPLETE

*This is the thesis's unique theoretical contribution. Knaus (2024) Appendix A.4
derives the PIVE representation for τ̂ₐ, τ̂ₐ,₀, τ̂ₐ,₁₀ but does not explicitly
derive ωᵢ for τ̂ᵤ, nor does it apply Love plots or ESS diagnostics to kappa
estimators empirically. This chapter bridges that gap.*

*Status: first draft complete — flagged for one reread pass before submission.*

**Section 3.1 — Kappa weights vs. outcome weights** ✅
- Distinction between κᵢ (identification objects, Abadie 2003) and ωᵢ
  (sample-level PIVE weights, Knaus 2024) clearly established
- Subsection 3.1.1 (Two notions of normalization): SUW normalization (estimator
  construction) vs. Knaus normalization (final weight properties);
  translation-invariance identity derived algebraically as eq:ti_identity —
  the sum-to-zero condition shown as necessary and sufficient in finite samples
- Love plots use ωᵢ, not κᵢ — stated explicitly

**Section 3.2 — PIVE framework** ✅
- Proposition 1 (Knaus 2024) stated formally with label prop:pive_weights
- Diagonal T matrix structure for kappa estimators established — no smoother
  condition needed; existence of ωᵢ is purely algebraic
- Sets up the derivation in Section 3.3

**Section 3.3 — Outcome weights for τ̂ᵤ** ✅ (main derivation)
- Full derivation: normalized IPW contrast → diagonal T^u → closed-form ωᵢᵘ
  (eq:u_omega_scalar, boxed)
- Three normalization conditions verified algebraically:
  Σωᵢ = 0 (via equal-mass property), ΣωᵢDᵢ = 1, Σωᵢ(1−Dᵢ) = −1
- Remark (rem:hajek_contrast): why τ̂ₐ,₁ fails where τ̂ᵤ succeeds — the Hájek
  normalization is the single algebraic step that determines translation
  invariance; connects to Appendix E derivation
- NOTE: τ̂ₐ,₁₀ handled via the Remark and classification table, not via a
  standalone derivation section.
- % TODO (compile): \ref{sssec:derive_a10} on line 863 — label does not exist,
  produces ??; fix to \ref{rem:hajek_contrast} or remove

**Section 3.4 — Normalization classification and summary** ✅
- Subsection 3.4.1: how SUW normalization and Knaus normalization align for the
  kappa family — co-occurrence is algebraic, not a general theorem
- Subsection 3.4.2: Table 1 (tab:normalization) — all six estimators classified
  by SUW norm., Σωᵢ=0, ΣωᵢDᵢ, Σωᵢ(1−Dᵢ), Knaus class
  % TODO: apply resizebox fix — table currently overflows right margin
- Three observations from the table written in prose

**Section 3.5 — From Derivation to Diagnostics: Computational Implementation** ✅
- kappa_outcome_weights(Z, D, p) described: returns all five ωᵢ vectors in
  closed form, no numerical optimisation
- check_weight_identity() and weight_diag() described as companion functions
- Pipeline: propensity score → weights → verify identity → Love plots and
  diagnostics tables
- Forward pointer to Chapter 4 written
- % TODO (compile): \ref{eq:ti_algebra} → \eqref{eq:ti_identity} on line 892
- % TODO (compile): remove forced \\ line break before \texttt on line 876

**Appendix E — Outcome Weights of τ̂ₐ,₁** ✅
- Full derivation: PIVE representation, closed-form ωᵢᵃ'¹ (boxed)
- Three normalization properties proved: Σωᵢ ≠ 0, ΣωᵢDᵢ = 1, Σωᵢ(1−Dᵢ) ≠ −1
- Classification: untreated-unnormalized in Knaus sense
- Cross-referenced from Table 1 footnote b

**Open compile issues — fix before next render:**
1. \ref{sssec:derive_a10} line 863 → label missing, produces ??
2. \ref{eq:ti_algebra} line 892 → change to \eqref{eq:ti_identity}
3. Forced \\ line 876 before \texttt → remove
4. Double \appendix call (lines 1112 + 1171) → merge into one block
5. \bibliography placed after second \appendix → move before first \appendix
6. tab:normalization overflows right margin → apply resizebox fix

---

### Chapter 4 — Empirical Application: Angrist (1990) Vietnam Draft Lottery
(7–9 pages) 🟡 IN PROGRESS

*Purpose in the thesis: establishes the diagnostic pipeline on a near-ideal
experiment. The design is so clean that the outcome-weights framework should show
complete uniformity across estimators and learners — ESS, % negative weights, and
Love plots are near-identical for all. This is the "design dominates learner"
application. It also serves as the replication benchmark: all estimates should
match SUW (2025) Table 2.*

**Section 4.1 — Data and design diagnostics** 🔲 WRITE FIRST
- Dataset: SIPP 1984, white men born 1950–1953, N = 3,027
- Instrument Z: Vietnam draft lottery eligibility (Angrist 1990)
- Treatment D: veteran status
- Outcome Y: log earnings (dollars and cents — both reported)
- Covariates: age (cubic or saturated year-of-birth dummies)
- Design diagnostics table — report explicitly:
  - Pr(Z=1), Pr(D=1), Pr(D=1|Z=1), Pr(D=1|Z=0)
  - First stage: Pr(D=1|Z=1) − Pr(D=1|Z=0) + first-stage F-statistic
  - Estimated complier share under IV assumptions
  - Compliance structure: one-sided noncompliance (no always-takers since
    draft-ineligible cannot be induced into service by lottery)
  - Instrument propensity: near Pr(Z=1) ≈ 0.5 (close to random)
- Opening sentence for this section:
  "Before comparing estimators, I summarize the empirical features of the design
  that are directly relevant for interpreting outcome weights: the assignment
  structure of the instrument, the strength of the first stage, the compliance
  pattern, and the covariate adjustment problem."
- NOTE: first-stage F is reported here as a design fact; its connection to
  ESS = 5 is made in Section 4.4, not here.

**Section 4.2 — Point estimates and translation-invariance replication** 🔲
- Replicate SUW (2025) Table 2: all five kappa estimators, cubic and saturated
  specs, dollars and cents
- Show normalized estimators (τ̂ᵤ, τ̂ₐ,₁₀) stable under outcome recoding;
  unnormalized (τ̂ₐ, τ̂ₐ,₁, τ̂ₐ,₀) change when Y → Y×100 (cents)


**Section 4.3 — DML learner comparison** 🔲
- Three learners (ranger, XGBoost, linear+logistic): estimates converge at
  0.218–0.246
- XGBoost flag: algebraic identity Σωᵢ Yᵢ = τ̂ fails marginally
  (Sum_w = −3.73e−06), consistent with Knaus (2024) Table 6 — XGBoost violates
  Condition 3 (non-affine smoother); documented for completeness, does not affect
  point estimate meaningfully
- "Design dominates learner" paragraph goes here: explain why ESS uniformity
  across learners is the expected result for this design, not a diagnostic failure

**Section 4.4 — Outcome weight diagnostics** 🔲
- Diagnostics table (cubic spec, dollar outcome):

| Estimator              | Estimate | Σωᵢ    | ESS | % neg. | max\|ω\| |
|------------------------|----------|--------|-----|--------|----------|
| τ̂ᵤᶜᵇ (kappa, CBPS)   | 0.210    | 0      | 5   | 54.4   | 0.025    |
| τ̂ᵤᵐˡ (kappa, MLE)    | 0.202    | 0      | 5   | 54.4   | 0.028    |
| τ̂ₐ,₁₀ (kappa)        | 0.204    | 0      | 5   | 54.4   | 0.029    |
| τ̂ₐ (unnorm.)         | 0.314    | 0.048  | 5   | 54.4   | 0.029    |
| Wald-AIPW (grf)       | 0.229    | 0      | 5   | 54.4   | 0.018    |
| Wald-AIPW (ranger)    | 0.246    | 0      | 5   | 54.4   | 0.022    |
| Wald-AIPW (XGBoost)   | 0.246    | ~0     | 5   | 54.4   | 0.022    |

- Interpretation narrative: ESS = 5 across all estimators reflects the ~16%
  complier share documented in Section 4.1, not a failure — connect back to
  first-stage F explicitly. % negative weights ≈ 54% expected under near-random
  assignment (compliers are a minority). Normalized estimators: Σωᵢ = 0 exact;
  τ̂ₐ Σωᵢ = 0.048 ≠ 0 confirms translation-invariance failure.
- check_weight_identity() passes for all estimators except XGBoost marginal case

**Section 4.5 — Love plots and covariate balance** 🔲
- Cubic and saturated specs: Love plots for kappa estimators and DML Wald-AIPW
  side-by-side
- All estimators achieve near-perfect balance (age is the only covariate;
  draft lottery is near-random)
- Diagnostic value: confirms the weight pipeline is correctly implemented before
  applying to richer designs — this is the baseline

**Open / 🟡 pending:**
- 🟡 Write-up of tuning sensitivity: dml_with_smoother() default vs.
  tune.parameters="all" — coded, not yet written
- 🟡 XGBoost non-affine issue: cite Knaus (2024) Table 6 formally in text

---

### Chapter 5 — Empirical Applications: Card (1995) and Angrist & Evans (1998)
(6–8 pages) 🔲

*Purpose in the thesis: Chapter 5 stress-tests the diagnostic framework. Card (1995)
shows what happens when the design is observational and covariate richness matters —
ESS diverges across specs, Kitagawa ESS ≈ 1 is a reliability flag. Angrist & Evans
(1998) provides the most demanding compliance test — near one-sided noncompliance,
most dramatic translation-invariance failure, sign flips across outcome units.*

**Section 5.1 — Card (1995): College education and wages** 🔲

*5.1.1 — Data and design diagnostics*
- Dataset: NLS Young Men cohort (Card 1995)
- Instrument Z: proximity to four-year college
- Treatment D: two definitions — some college (educ > 12) and completion (educ ≥ 16)
- Outcome Y: log wages
- Two covariate specifications:
  - Card (1995) full controls: age, race, family background, region, urban status
  - Kitagawa (2015) parsimonious: age + race only
- Design diagnostics — report explicitly:
  - Pr(Z=1), Pr(D=1), first stage for each treatment definition
  - First-stage F for Card vs. Kitagawa spec (does covariate adjustment change it?)
  - Estimated complier share: both treatment definitions
  - Compliance structure: two-sided (proximity to college creates always-takers
    and compliers; no one-sided noncompliance)
  - Instrument propensity distribution: overlap of p̂(X) across specs —
    Kitagawa parsimonious spec has less propensity score variation, contributes
    to ESS ≈ 1

*5.1.2 — Point estimates and translation-invariance replication*
- Replicate SUW (2025) Table 3: four spec × treatment combinations
  (Card/Kitagawa × somecol/educ16), kappa estimators
- Show normalization divergence: unnormalized estimates vary more across specs
  than normalized

*5.1.3 — Outcome weight diagnostics*
- Kitagawa spec: ESS ≈ 1 — extreme weight concentration; connect back to
  propensity score overlap noted in 5.1.1
- Card spec: ESS higher, better-spread weights
- Key contrast with Vietnam: here ESS varies by spec and learner — the design
  does not dominate; the outcome-weights lens is genuinely informative

*5.1.4 — Love plots and covariate balance*
- 8-panel grids: Card vs. Kitagawa spec, both treatment definitions
- Kitagawa spec shows more imbalance → covariate richness matters for complier
  balance in this observational design

**Status:** ✅ All coding complete (kappa + DML + love plots + weight diagnostics);
🔲 write-up pending

---

**Section 5.2 — Angrist & Evans (1998): Childbearing and labor supply** 🔲

*5.2.1 — Data and design diagnostics*
- Dataset: US Census 1980 (Angrist & Evans 1998)
- Instrument Z: same-sex siblings (also: twin second birth as robustness)
- Treatment D: third child (more than two children)
- Outcome Y: female labor force participation (LFP) and log income
- Covariates: age, race, education, family background controls
- Design diagnostics — report explicitly:
  - Pr(Z=1), Pr(D=1|Z=1), Pr(D=1|Z=0)
  - First stage and estimated complier share
  - Compliance structure: near one-sided noncompliance — no always-takers
    (parents with opposite-sex children are not induced into having a third child
    by the same-sex instrument, but same-sex parents may be; verify empirically)
  - → Proposition 3.3 (SUW 2025) applies: τ̂ₐ,₁₀ faces near-zero denominator
    risk; τ̂ᵤ safe
  - Instrument propensity: near 0.5 (same-sex probability ≈ 50%)

*5.2.2 — Point estimates and translation-invariance replication*
- Replicate SUW (2025) Table 4: LFP and log income outcomes
- Most dramatic TI failure: log income estimates flip sign when moving from
  dollars to log-dollars or alternative unit scalings
- τ̂ₐ,₁ and τ̂ₐ,₀ show sign reversal; τ̂ᵤ stable

*5.2.3 — Outcome weight diagnostics*
- DML comparison: Wald-AIPW only (AIPW-ATE gives NaN under near one-sided
  noncompliance — documented as expected behavior, not an error)
- Weight diagnostics for kappa estimators: check τ̂ₐ,₁₀ denominator behavior
  under near one-sided noncompliance

*5.2.4 — Love plots*
- Love plots for kappa estimators: near one-sided noncompliance affects weight
  structure vs. Vietnam and Card

**Status:** ✅ kappa done + love plots done;
🟡 DML Wald-AIPW coding in progress

---

### Chapter 6 — Discussion (3–4 pages) 🔲

**Section 6.1 — Cross-application summary**
- Comparative design diagnostics table across all three applications:

| Quantity            | Vietnam (1990)    | Card (1995)      | A&E (1998)         |
|---------------------|-------------------|------------------|--------------------|
| N                   | 3,027             | ~3,000           | ~250,000           |
| Instrument          | Draft lottery     | College proximity| Same-sex siblings  |
| Pr(Z=1)             | ≈ 0.5             | ...              | ≈ 0.5              |
| First stage         | ...               | ...              | ...                |
| Compliance          | One-sided         | Two-sided        | Near one-sided     |
| Covariate dim.      | Low               | Medium/high      | Medium             |
| ESS (normalized)    | 5 (uniform)       | Varies by spec   | ...                |
| Design dominates?   | Yes               | No               | Partial            |

- Key contrast narrative: Vietnam (near-random, low-dim) → ESS uniform, design
  dominates; Card (observational, rich X) → ESS diverges, diagnostics informative;
  AE (one-sided noncompliance) → Proposition 3.3 risk confirmed empirically

**Section 6.2 — The outcome weights lens: what it adds**
- What do Love plots reveal beyond point estimates alone?
- When do weight diagnostics signal problems the estimates don't?
  (Card Kitagawa ESS ≈ 1 is a red flag even if point estimate looks reasonable)
- DML vs. kappa: are Wald-AIPW and τ̂ᵤ targeting the same complier subpopulation?

**Section 6.3 — DML learner comparison: implications**
- Point estimates converge for Vietnam (0.218–0.246) and Card
- Design dominates learner: in low-dimensional near-random-assignment settings,
  learner flexibility doesn't add value — ESS and % negative weights identical
  across learners; expected outcome, not methodological failure
- XGBoost non-affine smoother issue: always verify Σωᵢ Yᵢ = τ̂ numerically
- Where learner matters: richer observational designs (Card) may show divergence
  — open question

**Section 6.4 — Practical guidance for practitioners**
- Decision framework: when to use τ̂ᵤ vs. τ̂ᵤᶜᵇ vs. Wald-AIPW; 2SLS defensible
  only with saturated covariate specification
- Always use normalized estimators; run sum-to-zero check (Σωᵢ = 0)
- Check ESS: ESS ≈ 1 is a reliability red flag even if point estimate looks
  plausible; connect to first-stage F and complier share
- For DML: use ranger or grf; XGBoost needs additional algebraic verification
- Package note: kappa_to_outcome_weights_format() wrapper enables unified
  Love-plot pipeline

**Section 6.5 — Limitations**
- Bootstrap inference computationally expensive; analytical SEs not implemented
  for all estimators
- Love plots descriptive, not inferential
- Weak overlap creates extreme weights in all estimators
- DML estimates depend on tuning choices; RF hyperparameter tuning affects
  balance (cf. Knaus 2024 Figure 3)
- AE DML analysis incomplete (memory constraints at full N)

---

### Extension — Heterogeneous Treatment Effects via Instrumental Forest 🔲

**Condition:** only pursued if ATE results motivate it (professor's gate).

**Current assessment:** Vietnam ATE convergence → weak motivation. Card (1995)
richer covariates and Kitagawa vs. Card ESS divergence → stronger motivation.
Decision after professor meeting.

**Status:** 🔲 Conditional on Card write-up and professor feedback

---

### Chapter 7 — Conclusion (1–2 pages) 🔲

- Summary of findings: which estimators are translation invariant, which achieve
  covariate balance, what the outcome-weights lens adds
- Summarize "design dominates learner" as a cross-application empirical result
- Recommendation: τ̂ᵤᶜᵇ preferred for robustness; Wald-AIPW for flexibility when
  large N and rich X; 2SLS defensible only with saturated covariate specification;
  check ESS before trusting any point estimate
- Extensions: doubly robust kappa estimators; formal tests of covariate balance;
  heterogeneous effects; kappa_to_outcome_weights_format() as package contribution

---

## What to write next — recommended order

Given where you are (Ch. 3 complete, all coding done, Ch. 4 tables exist as
skeletons):

**Start here: Section 4.1 Data and design diagnostics (Vietnam)**
This is the right entry point because:
- It is purely factual prose — no interpretation, no theory, just design facts
- The numbers are all in the roadmap already (N, Pr(Z), first stage, compliance)
- Writing it forces you to produce the design diagnostics table that you will
  reuse as a template for 5.1.1 and 5.2.1
- It is short (≈ 1 page) and completable in one sitting
- Once it exists, Section 4.4 (weight diagnostics) can reference it by name

**Then: "Design dominates learner" paragraph (Section 4.3)**
Writing plan correctly identifies this as the single highest-value short piece.
Write it as a standalone 200-word block immediately after 4.1. The numbers are
in the diagnostics table above. Do not wait for the rest of Ch. 4 to exist.

**Then: Section 4.2 and 4.4 (point estimates + weight diagnostics)**
The table skeletons exist. Fill the numbers, write the interpretation. Use
% REREAD comments liberally — do not polish, just produce.

**Then: Chapter 2 (framework)**
Once Ch. 4.1–4.4 draft exists, Ch. 2 writes itself as setup. All definitions
are in the roadmap reading notes. Keep the discipline: Ch. 2 contains no results.

**Do not touch yet: Ch. 1, Ch. 6, Ch. 7**
These are genuinely last. Ch. 1 needs all findings to exist; Ch. 6 needs all
empirical chapters; Ch. 7 follows Ch. 6.

**Perfectionist protocol:**
Add % REREAD at the top of every section you write in first-draft mode.
Do not go back to polish until all chapters have a first draft.
Ch. 3 is the only chapter that should not have a % REREAD tag — it is done.
