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
