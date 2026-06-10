# Supervisor Meeting — Wednesday
## Meeting Prep & Agenda

---

##  Topics to Discuss

---

### 1. Low ESS in Vietnam and Card — Why and Does It Matter?

**The finding:** ESS is extremely low across both settings — sometimes 3–5, occasionally near 1 (especially Kitagawa spec in Card). This holds for both the kappa estimators and the Wald-AIPW via Double ML.

**Explanation:**
- **Vietnam:** The design is near-random-assignment and low-dimensional (only age controls). The instrument propensity score is tightly clustered around 0.456 across all specs. When p̂(Xᵢ) is nearly constant, the kappa weights collapse onto a small number of compliers.
- **Card:** The Kitagawa spec (ESS ≈ 1) reflects a weak-instrument problem compounded by near-overlap violations. When the first stage is weak, the denominator (E[D|Z=1,X] − E[D|Z=0,X]) is close to zero for many observations, producing extreme weights. This is *qualitatively different* from Vietnam — here the spec choice actually matters.

**Key question:**
> Does the low ESS in Vietnam invalidate or undermine the thesis narrative, or is it precisely the *expected* result for a near-random low-dimensional design and should be framed as a confirmatory diagnostic rather than a warning sign?

**Sub-question (Card-specific):**
> In the Card setting, p̂(Xᵢ) can be high even for units with Z=0, which directly inflates the kappa denominators. Is this the right framing for why the Kitagawa spec collapses to ESS ≈ 1, and how should it be communicated in the write-up?

---

### 2. XGBoost Produces Anomalous Results — What to Do?

**The finding:** The algebraic check `ω'Y = τ̂` fails for XGBoost in both Vietnam and Card. The issue is traced to `base_score = 0` in XGBoost ≥ 1.6, which changed behavior and violates the affine smoother condition required by Knaus & Rakov.

**Your current position:** Results are reported but flagged. Recommendation in the thesis will be to use ranger or grf for standard use; XGBoost requires additional algebraic verification.

**Question for Michael:**
> Is it sufficient to document the XGBoost non-affine smoother issue and flag it as a limitation, or does the non-passing algebraic check make the XGBoost results unpresentable/unreliable entirely? Should those rows be dropped from the learner comparison table or kept with a caveat?

---

### 3. Smoother Matrix Memory Problem — Angrist & Evans (1998)

**The problem:** `dml_with_smoother()` requires building a full N×N smoother matrix. With N ≈ 254,000 (LFP sample), this is ~500 GB in float64. Already hit `vector memory limit of 16 Gb` on both the LFP and income subsamples. DML outcome-weight diagnostics cannot be reported for Angrist & Evans as a result.

**Current resolution documented in thesis:** Design limitation. 

**Questions for Michael:**
> 1. Is it acceptable to simply exclude DML diagnostics for Angrist & Evans and document this transparently
> 2. Subsampling (~5k obs) was considered and rejected because point estimates would no longer match the kappa replication table, does he agree with that reasoning?

---

### 4. Heterogeneous Treatment Effects — Does It Fit?

**Background:** Discussed previously as a potential extension (instrumental forest / CATE). The ROADMAP flags it as conditional on ATE results motivating it.

**Current assessment:**
- Vietnam: estimates converge tightly across learners, weak motivation
- Card: richer covariates + Kitagawa vs. Card spec divergence in ESS 
- Angrist & Evans: near one-sided noncompliance makes CATE translation invariance especially interesting, but it's the most demanding setting

**The honest tension:** It feels like it could only fit naturally in the Angrist & Evans section but that section is the least complete, and it's unclear whether it fits the storyline or just adds bulk.

**Questions for Michael:**
> 1. Is heterogeneous effects via instrumental forest expected or is it genuinely optional?
> 2. If pursued, which application is the most defensible entry point Card (richer X, ESS divergence) or Angrist & Evans (compliance structure motivation)?
> 3. Or is it better to mention it as a clean extension in Section 6.5 / Conclusion and not pursue the code at all?

---

### 5. Career / What's Next — Causal ML in Practice


**What I actually want to know:**
- Is this kind of empirical causal ML work more of an academic niche (→ PhD path) or is there meaningful demand in industry too?
- Which departments / sectors are actually doing this right now tech, policy evaluation, pharma, finance?
- Given how fast the AI/data landscape is moving, does he see causal inference becoming more or less central?

**How to phrase it:**
>  I'd genuinely love your read on where this is heading. Not just for me personally, but because I'm trying to understand the landscape while I still have the chance to ask someone who's shaping it.


---

## 🔲 Things Still to Prepare Before Wednesday

- [ ] **Clean up the notebooks** — the current HTML files display functions inline rather than using `functions_all.R`. At minimum, make sure the two notebooks you show him (Vietnam + Card) render clearly and the key tables/plots are visible without scrolling through raw code.
- [ ] **Prepare a 2-sentence narrative for each key result** so you can walk through the tables confidently without reading off the screen.
- [ ] **Double-check the XGBoost algebraic check output** — know the exact numbers so you can present the failure clearly.
- [ ] **Review Problems.md** before the meeting — make sure you can articulate the smoother memory issue concisely (one minute, not five).
- [ ] If time allows: write the "design dominates learner" framing paragraph (ROADMAP TODO #1) so you can show him you already know how to interpret the Vietnam ESS result.

---

## 💡 Additional Things Worth Raising (if there's time)

- **PLR-IV vs. kappa target estimator ambiguity** — your Problems.md flags that PLR-IV and kappa don't always estimate exactly the same target. Is this worth a footnote or does it need more formal treatment?
- **Translation invariance check Col 1 → Col 3** —the sign-flip issue (Y' = −Y + 1 isn't a pure translation). You've flagged this in Problems.md. Make sure the expression in the thesis is correct before he reads it.

