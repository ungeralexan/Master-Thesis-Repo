## Phase 1 — Two short pieces that unlock everything else (this week)
Start here: Chapter 3, Section 3.2 — the outcome weight derivations.
You already wrote this in the vignette. Convert it from Rmd to LaTeX prose. This is ~2 pages, you have the math, and it is the unique theoretical contribution of the thesis. Writing it first gives you the clearest possible statement of what Chapter 3 is doing — which then makes Chapter 2 easy to frame ("this chapter provides the building blocks for Section 3.2").
Then: the "design dominates learner" paragraph (the single highest-value piece of writing you have, currently just a bullet on the roadmap). Write it as a self-contained 200-word block — what the finding is, why it is expected, what it means for interpretation. This will slot into Section 4.4 but once it is written you can reference it everywhere.
These two pieces together take one focused day and unlock the rest.


## Phase 2 — Chapter 4, Vietnam (next 3–4 days)
This is your most complete chapter: fully coded, three covariate specs, DML comparison, Love plots, TI check, weight diagnostics, learner comparison — everything is there. Write it as:
4.1 Data and design (~1 page) — SIPP 1984, N, instrument design, one-sided noncompliance. Already in your roadmap bullet by bullet, just needs prose.
4.2 Point estimates and replication (~1.5 pages) — reproduce your Table 2 replication results, cents vs. dollars comparison, normalized vs. unnormalized divergence. The numbers are in child_repl_30 equivalents for Vietnam. Write what the table shows.
4.3 DML comparison and learner invariance (~1.5 pages) — PLR-IV, Wald-AIPW estimates, the three-learner comparison. The "design dominates learner" paragraph goes here.
4.4 Outcome weight diagnostics and Love plots (~2 pages) — the weight diagnostics table, the Love plot grid, what ESS = 5 means in this context, why % negative weights are expected.
Writing Chapter 4 first is the right move because it is the most straightforward chapter and finishing it gives you momentum and a template for Chapter 5.


## Phase 3 — Chapter 2, the framework (after Chapter 4 exists)
Once you have Chapter 4, you know exactly what tools the reader needs to have understood before reaching it. Write Chapter 2 as a technical preparation for Chapter 4, not as a free-standing literature review. Sections 2.1–2.4 (IV, kappa weights, estimators, translation invariance) are close to mechanical given your reading notes — these are 1 page each. Sections 2.5–2.7 (DML, PIVE, covariate balance) follow the same pattern.
The key discipline here: Chapter 2 should contain no results, no interpretation, only definitions and facts. If you catch yourself writing "which is why Vietnam shows..." — stop, that belongs in Chapter 4.

## Phase 4 — Chapter 3, the theoretical bridge (after Chapter 2)
Section 3.1 (kappa weights vs. outcome weights distinction) is one page. Section 3.2 you already have from Phase 1. Section 3.3 is a summary table — also already drafted in the vignette. This chapter practically writes itself once Chapter 2 gives you the notation.


## Phase 5 — Chapter 5, Card + AE (after Chapter 4)
Section 5.1 (Card) follows the exact same structure as Chapter 4 — data, estimates, TI check, weight diagnostics, Love plots — but the narrative is different: ESS diverges here across specs, Kitagawa ESS ≈ 1 is a reliability flag, contrast with Vietnam. Write this by copying Chapter 4's structure and replacing the findings.
Section 5.2 (AE) is shorter: kappa only, DML infeasible at full N, one-sided noncompliance makes tau_ml_a0 explode, the most dramatic TI failure. The memory limitation becomes a one-paragraph limitation statement, not a gap.


## Phase 6 — Introduction and Discussion (last)
Chapter 6 Discussion writes itself from the cross-application comparison table (Section 6.1) and the three narrative threads you already have: design dominates learner, normalization matters, outcome weights add diagnostic value. Write it after all empirical chapters exist.
Chapter 1 Introduction is genuinely last. At that point you know exactly what you found, what the gap was, and what you contributed. The hook (cents vs. dollars failure), background (2SLS doesn't recover LATE), gap (Knaus Appendix A.4 not applied empirically), and contribution bullet points are already in your Gliederung — you just need to convert them to flowing prose.