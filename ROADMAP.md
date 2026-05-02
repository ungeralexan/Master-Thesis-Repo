# Thesis Roadmap

**Last updated:** May 2026  
**Phase:** Reading complete → Replication active

---

## Reading status

### ✅ Tier 0 — Core papers (DONE)

| Paper | Status | Notes |
|---|---|---|
| Imbens & Angrist (1994) | ✅ Read | 3 conditions, 3 theorems, LATE as Wald ratio — see `output/notes/imbens_angrist_1994_notes.pdf` |
| Abadie (2003) | ✅ Read | Kappa theorem, κ₀/κ₁/κ, Lemma 2.1, Prop 5.1 — see `output/notes/abadie_2003_notes.pdf` |
| Słoczyński, Uysal & Wooldridge (2025) | ✅ Read | All sections. Five estimators, Prop 3.2 (TI + SE), Sec 3.4 (near-zero denom), three applications |

### 🟡 Tier 1 — Foundational (in progress)

| Paper | Status | What to read |
|---|---|---|
| Angrist, Imbens & Rubin (1996) | ✅ Read | Sections 1–2 only: four compliance types, AIR notation |
| Chernozhukov et al. (2018) | 🔲 TODO | Section 2.3 only: Wald-AIPW estimator formula |
| Knaus (2024) | 🟡 Started | Read fully. Key: Sec 2 (PIVE), Sec 3.2.1 (instrumental forest), **Appendix A.4** (kappa in PIVE — your Section 3) |

### 🔲 Tier 2 — Empirical application papers

| Paper | Status | What to extract |
|---|---|---|
| Angrist (1990) | 🔲 TODO | Intro + data: Z = draft lottery, D = veteran, Y = wages. No always-takers? |
| Card (1993/1995) | 🔲 TODO | Intro + instrument: Z = college proximity. Two treatment specs (somecol, educ16) |
| Angrist & Evans (1998) | 🔲 TODO | Intro + Table 1: two instruments (twins, same-sex). AE98 = near one-sided noncompliance |

### 🔲 Tier 3 — Cite-and-move-on (don't read deeply)

| Paper | Purpose |
|---|---|
| Uysal (2011) | Origin of τ̂ᵤ. Formula fully in SUW (2025) |
| Frölich (2007) | Origin of τ̂ₜ = τ̂ₐ,₁. One paragraph in lit review |
| Heiler (2022) | Covariate balancing for LATE. Skim Sec 2 for CBPS definition |
| Blandhol et al. (2022) | "When is TSLS LATE?" Abstract + intro only — motivation for kappa |
| Angrist & Pischke (2009) | General IV textbook reference. Chapter 4 if IV intuition needed |

---

## Thesis outline (draft)

### Chapter 1 — Introduction (2–3 pages)
- IV estimation is ubiquitous but 2SLS has limitations under heterogeneous effects
- Kappa weighting as a flexible alternative (Abadie 2003, SUW 2025)
- Contribution: replication + extension via outcome weights + Love plots
- Road map of the thesis

### Chapter 2 — Theoretical Framework (6–8 pages)

**2.1 The LATE setup**
- Potential outcomes notation: Y(0), Y(1), D(0), D(1)
- Four compliance types (AIR 1996)
- Assumption IV (i)–(iv): conditional independence, exclusion, first stage, monotonicity
- LATE definition: E[Y₁ − Y₀ | D₁ > D₀]

**2.2 The kappa theorem**
- Lemma 2.1 (Abadie 2003 restated in SUW notation)
- The three weights κ, κ₁, κ₀ and their cell-by-cell values (Table 1 of SUW)
- Parts (a), (b), (c): any complier moment identified
- Remark 2.2: E(κ) = E(κ₁) = E(κ₀) in population; why this matters for finite samples

**2.3 The five estimators**
- Equations (1) and (2) as distinct starting points
- τ̂ₐ, τ̂ₐ,₁ (= τ̂ₜ), τ̂ₐ,₀: unnormalized, three denominator choices
- τ̂ₐ,₁₀: normalized, Abadie-Cattaneo (2018)
- τ̂ᵤ: normalized, Uysal (2011) — recommended estimator

**2.4 Why normalization matters**
- Definition TI (translation invariance) and Definition SE (scale equivariance)
- Proposition 3.2: τ̂ᵤ and τ̂ₐ,₁₀ pass; τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ fail
- Mechanical explanation: finite-sample non-cancellation of instrument residuals

**2.5 One-sided noncompliance**
- Definition and examples (401k, randomized trials, twin births)
- Table 1 cell-by-cell signs of κ₁ and κ₀
- Proposition 3.3: positive denominators guaranteed
- Proposition 3.4: τ̂ᵤ denominator positive in both one-sided cases
- Why τ̂ₐ,₁₀ fails here but τ̂ᵤ doesn't

**2.6 Estimation of the propensity score**
- ML logit (τ̂ᵤᵐˡ) vs covariate balancing CBPS (τ̂ᵤᶜᵇ)
- Proposition 3.5: with CBPS, all normalized estimators are identical
- Why CBPS regularizes weights away from extremes (Heiler 2022 argument)

### Chapter 3 — Outcome Weights Framework (4–5 pages)
- Knaus (2024) PIVE framework: τ̂ = Σᵢ ωᵢYᵢ
- Deriving outcome weights for τ̂ᵤ, τ̂ₐ,₁₀, τ̂ₐ (following Knaus Appendix A.4)
- Normalization properties of ωᵢ: do treated weights sum to +1, untreated to −1?
- Love plots as covariate balance diagnostics using ωᵢ

### Chapter 4 — Empirical Applications (8–10 pages)

**4.1 Angrist (1990) — Military service**
- Setup: Z = draft lottery, D = veteran, Y = log wages
- Covariate specs: linear age, cubic age, saturated age
- Table 2 replication + discussion of cents vs dollars failure

**4.2 Card (1995) — College education**
- Setup: Z = college proximity, D = some college / completion, Y = log wages
- Two covariate specs: Card, Kitagawa
- Table 3 replication + discussion of Card vs Kitagawa divergence

**4.3 Angrist & Evans (1998) — Childbearing**
- Setup: Z = same-sex, D = third child, Y = LFP + log income
- Table 4 replication + most dramatic illustration of unnormalized failure
- One-sided noncompliance structure and Proposition 3.3 in practice

**4.4 Outcome weights and Love plots** *(extension beyond SUW 2025)*
- For each application: compute ωᵢ for τ̂ᵤ, τ̂ₐ,₁₀, τ̂ᵤᶜᵇ
- Love plots: standardized mean differences before/after kappa reweighting
- Covariate balance comparison: normalized vs unnormalized estimators
- Connection to Knaus (2024) Appendix A.4

### Chapter 5 — Conclusion (1–2 pages)
- Summary of findings
- Policy implications: which estimator to use and when
- Limitations and future work (doubly robust estimators, weak overlap)

---

## Code task tracker

### ✅ Done
- [x] Replicate Angrist (1990): 2SLS + all five kappa estimators, cents vs dollars
- [x] Replicate Card (1995): two treatments, two specs, two outcomes
- [x] Replicate Angrist & Evans (1998): LFP + log income, all transformations
- [x] Bootstrap standard errors for all kappa estimators
- [x] Outcome weights functions (kappa_outcome_weights)
- [x] Basic weight diagnostics (sum of weights, ESS, % negative, max weight)

### 🟡 In progress
- [ ] Love plots for all three applications (cobalt package)
- [ ] Clean modular R scripts (currently in one Rmd)

### 🔲 TODO
- [ ] Replicate SUW simulation study (Designs A–D from Heiler 2022)
- [ ] Analytical standard errors (M-estimator sandwich) as alternative to bootstrap
- [ ] Final replication tables in LaTeX format
- [ ] Thesis LaTeX draft (`tex/thesis.tex`)
- [ ] Optional: 401(k) application via Knaus OutcomeWeights package

---

## Who introduced what — quick reference

| Concept | Source |
|---|---|
| Potential outcomes + LATE definition | Imbens & Angrist (1994) |
| Four subpopulations (compliers etc.) | Angrist, Imbens & Rubin (1996) |
| Kappa theorem (κ identifies any complier moment) | Abadie (2003) |
| τ̂ₜ estimator (ratio of two IPW estimators) | Frölich (2007) / Tan (2006) |
| τ̂ᵤ estimator (normalized, recommended) | Uysal (2011) |
| Translation invariance criterion + all five compared | Słoczyński, Uysal & Wooldridge (2025) |
| CBPS for LATE | Heiler (2022) |
| 2SLS doesn't recover LATE | Blandhol et al. (2022) |
| Wald-AIPW / DML framework | Chernozhukov et al. (2018) |
| Outcome weights framework for ML estimators | Knaus (2024) |
| Kappa estimators in PIVE framework | Knaus (2024) Appendix A.4 |
