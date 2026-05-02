# Kappa Weighting Estimators of the LATE — Thesis Repository

**Author:** Alexander Unger  
**Institution:** [Your University]  
**Supervisor:** [Supervisor Name]  
**Status:** 🟡 In progress — reading phase complete, replication phase active

---

## What this thesis is about

This thesis studies **kappa weighting estimators of the Local Average Treatment Effect (LATE)**, building on Słoczyński, Uysal & Wooldridge (2025) and Knaus (2024).

The central question: when estimating the LATE with covariates using instrumental variables, which estimator should you use — and why?

The thesis makes three contributions:

1. **Replicates** the three empirical applications in Słoczyński et al. (2025) in R, comparing 2SLS to five kappa-weighted estimators across normalized and unnormalized variants
2. **Extends** the analysis by placing the kappa estimators into the outcome-weights framework of Knaus (2024), computing Love plots for covariate balance diagnostics
3. **Documents** the finite-sample failure modes (translation invariance, scale equivariance, near-zero denominators) with concrete numerical evidence

---

## Repository structure

```
thesis-kappa-late/
│
├── README.md                    ← this file
├── ROADMAP.md                   ← reading list + thesis outline + current status
│
├── R/                           ← all R scripts, numbered by execution order
│   ├── 00_functions.R           ← kappa weight functions, estimators, bootstrap SE
│   ├── 01_angrist1990.R         ← Angrist (1990): military service / draft lottery
│   ├── 02_card1995.R            ← Card (1995): college proximity / education
│   ├── 03_angrist_evans1998.R   ← Angrist & Evans (1998): childbearing / labor supply
│   ├── 04_outcome_weights.R     ← Knaus (2024) extension: Love plots, weight diagnostics
│   └── applied_kappa.Rmd        ← master Rmd notebook (renders all results)
│
├── data/                        ← raw data files (NOT committed — see below)
│   ├── sipp.dta                 ← Angrist (1990) — SIPP 1984
│   ├── card.dta                 ← Card (1995) — NLSYM
│   └── ae98.dta                 ← Angrist & Evans (1998) — 1980 Census
│
├── output/
│   ├── tables/                  ← replication tables (CSV + LaTeX)
│   ├── figures/                 ← Love plots, weight distribution plots
│   └── notes/                   ← compiled reading notes (PDF)
│       ├── imbens_angrist_1994_notes.pdf
│       ├── abadie_2003_notes.pdf
│       └── suw_2025_notes.pdf   ← (to add)
│
├── tex/                         ← thesis LaTeX source (to be added)
│   └── thesis.tex
│
└── docs/                        ← rendered HTML outputs for sharing
    └── applied_kappa.html       ← current working notebook
```

---

## Data

Raw `.dta` files are **not committed** (see `.gitignore`) to keep the repo lightweight and respect data terms. Place the three Stata files from the SUW (2025) replication package into `data/` after cloning:

```
data/sipp.dta    # N = 3,027  — Angrist (1990) SIPP subsample
data/card.dta    # N = 3,010  — Card (1995) NLSYM subsample
data/ae98.dta    # N = 394,840 — AE (1998) 1980 Census subsample
```

Data source: supplementary materials of Słoczyński, Uysal & Wooldridge (2025), JBES.

---

## How to run

```r
# 1. Install required packages
install.packages(c("haven", "AER", "sandwich", "lmtest", "boot", 
                   "OutcomeWeights", "ggplot2", "cobalt"))

# 2. Set your data path in R/00_functions.R (line ~10)

# 3. Render the full notebook
rmarkdown::render("R/applied_kappa.Rmd")

# 4. Or run scripts individually in order
source("R/00_functions.R")
source("R/01_angrist1990.R")
source("R/02_card1995.R")
source("R/03_angrist_evans1998.R")
source("R/04_outcome_weights.R")
```

---

## Key papers

| Paper | Role |
|---|---|
| Imbens & Angrist (1994), *Econometrica* | LATE identification: compliers, monotonicity, Wald ratio |
| Angrist, Imbens & Rubin (1996), *JASA* | Four compliance types; exclusion restriction |
| Abadie (2003), *JoE* | Kappa theorem: any complier moment identified with covariates |
| **Słoczyński, Uysal & Wooldridge (2025), *JBES*** | **Core paper: five estimators, normalization, translation invariance** |
| **Knaus (2024), *JoE*** | **Extension: outcome weights framework, Love plots for IV** |

---

## Current status

See `ROADMAP.md` for the full reading list, thesis outline, and task tracker.
