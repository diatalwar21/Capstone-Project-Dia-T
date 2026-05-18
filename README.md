# California Paid Family Leave and Female Labor Force Participation

Replication package for the ECON 620 capstone paper:  
**"The Effect of California's Paid Family Leave on Female Labor Force Participation: A Synthetic Difference-in-Differences Approach"**  
Dia Talwar, University of San Francisco, 2026

Overleaf (paper source): https://www.overleaf.com/read/gctxvjjgtktf#a3e322

---

## Research Summary

California's 2004 Paid Family Leave (PFL) program was the first state-level PFL law in the United States. Using Synthetic Difference-in-Differences (Arkhangelsky et al., *American Economic Review* 2021), I estimate the program's effect on the female labor force participation (LFP) rate among women aged 25–54, using CPS Annual Social and Economic Supplement microdata from 1995–2015.

**Main result:** PFL increased the female LFP rate by **+1.19 percentage points** (95% CI: [+0.23, +2.15]; jackknife p = 0.019) over the 2007–2015 post-treatment window, after excluding the 2004–2006 adjustment period during which program take-up was low.

---

## Repository Structure

```
.
├── code/
│   ├── 01_data_construction.py      # Clean CPS microdata → state×year panel
│   ├── 02_main_estimation.R         # SDiD via synthdid (Arkhangelsky et al.)
│   ├── 03_robustness_donors.py      # Donor pool sensitivity (leave-one-out, subsets)
│   ├── 04_robustness_windows.py     # Post-window sensitivity (4 specifications)
│   ├── 05_placebo_inference.py      # Per-year jackknife CIs + in-space placebo
│   └── 06_figures.py                # All paper figures + summary statistics table
├── data/
│   ├── raw/                         # CPS microdata (not tracked — see data/raw/README_raw_data.md)
│   └── derived/
│       ├── state_year_panel_deduped_1995_2015.csv   # Main estimation panel
│       └── state_year_predictors.csv                # Panel for R synthdid spec
├── outputs/
│   ├── figures/                     # All pres_*.png/.pdf figures
│   └── tables/                      # All sdid_*.csv and .tex output tables
└── paper/
    └── capstone_paper.tex           # LaTeX source (figures auto-resolved via \graphicspath)
```

---

## Data

**Source:** Current Population Survey (CPS) Annual Social and Economic Supplement (ASEC), 1995–2015, accessed via IPUMS CPS.

**Sample:** Women aged 25–54. Aggregate to state × year cells; outcome is the female LFP rate (%).

**Analytical sample:** California + 43 donor states. Excluded jurisdictions:
- New York, Rhode Island, Hawaii — existing Temporary Disability Insurance programs structurally confound the counterfactual
- New Jersey (PFL 2008), Washington (PFL 2007) — enacted own PFL laws during the post-treatment window
- Alaska — structural outlier in labor market trends
- District of Columbia — city-state, not comparable to state units

**Derived files (committed):**
- `data/derived/state_year_panel_deduped_1995_2015.csv` — 1,071 state×year cells (51 jurisdictions × 21 years), outcome + demographic covariates, individual-level observations deduplicated per `01_data_construction.py`
- `data/derived/state_year_predictors.csv` — subset used by the R `synthdid` spec (`02_main_estimation.R`)

**Raw CPS microdata:** Not committed (files exceed GitHub's 100 MB limit). See `data/raw/README_raw_data.md` for download instructions.

---

## Replication Instructions

### Step 0 — Clone and install dependencies

```bash
git clone https://github.com/diatalwar21/Capstone-Project-Dia-T.git
cd Capstone-Project-Dia-T
pip install -r requirements.txt
```

For R:
```r
install.packages(c("synthdid", "dplyr"))
```

### Step 1 — (Optional) Rebuild the panel from raw CPS data

Only needed if you want to re-run data construction from scratch. Follow the download instructions in `data/raw/README_raw_data.md` first, then:

```bash
python code/01_data_construction.py
```

The derived CSVs in `data/derived/` are already committed, so this step can be skipped.

### Step 2 — Main estimation (R)

```bash
Rscript code/02_main_estimation.R
```

Produces SDiD point estimates and unit/time weights via the `synthdid` package.

### Step 3 — Robustness checks (Python)

```bash
python code/03_robustness_donors.py     # donor pool sensitivity
python code/04_robustness_windows.py    # post-window sensitivity → outputs/tables/window_sensitivity_table.csv
```

### Step 4 — Placebo inference (Python)

```bash
python code/05_placebo_inference.py
```

Produces per-year jackknife 95% CIs and in-space permutation null bounds. Outputs to `outputs/tables/` and `outputs/figures/`.

### Step 5 — Figures and summary statistics table (Python)

```bash
python code/06_figures.py
```

Produces all `pres_*.png/.pdf` figures and `summary_stats_1995_2015.csv` in `outputs/`.

### Step 6 — Compile the paper

From the `paper/` directory, compile with any LaTeX distribution (pdflatex or XeLaTeX). The `\graphicspath{{../outputs/figures/}}` declaration in the preamble means figures are resolved automatically — no path editing needed.

```bash
cd paper && pdflatex capstone_paper.tex && pdflatex capstone_paper.tex
```

---

## Key Results

| Specification | ATT (pp) | SE | 95% CI | Jack. p |
|---|---|---|---|---|
| Baseline SDiD (2004–2015) | +0.66 | 0.50 | [−0.31, +1.63] | 0.191 |
| **Robust SDiD (2007–2015)** | **+1.19** | **0.49** | **[+0.23, +2.15]** | **0.019** |

Donor pool: 43 states. Pre-treatment RMSPE: 0.726 pp (robust spec). Permutation p (in-space placebo): 0.326.

---

## Identification Strategy

**Synthetic DiD** (Arkhangelsky et al. 2021) combines:
- **Unit weights (ω):** which donor states best replicate California's pre-treatment trend in female LFP
- **Time weights (λ):** which pre-treatment years best approximate the post-treatment period

The ATT is a doubly-weighted 2×2 DiD using both weight sets. Standard errors use a leave-one-unit-out jackknife. The in-space placebo test assigns the treatment to each of the 43 donor states in turn and checks whether California's ATT is unusual relative to the null distribution.

---

## Citation

If you use this code or data, please cite:

> Talwar, Dia. "The Effect of California's Paid Family Leave on Female Labor Force Participation: A Synthetic Difference-in-Differences Approach." ECON 620 Capstone, University of San Francisco, 2026.

And the estimator:

> Arkhangelsky, Dmitry, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. "Synthetic Difference-in-Differences." *American Economic Review* 111, no. 12 (2021): 4088–4118.

---

## License

Code: MIT. Data: subject to IPUMS CPS terms of use (non-commercial, non-redistribution of raw microdata).


