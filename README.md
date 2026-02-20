# Capstone-Project-Dia-T
 
## California Paid Family Leave and Female Labor Force Participation Analysis

# California Paid Family Leave and Female Labor Force Participation Analysis

## Data Source and Sample Construction

### Primary Data Source
This analysis uses monthly microdata from the **Current Population Survey (CPS)**, accessed through **IPUMS-CPS**. The CPS is a monthly household survey conducted by the U.S. Census Bureau and Bureau of Labor Statistics and represents the primary source of U.S. labor force statistics.

### Sample Period
**January 2000 – December 2010 (132 months)**

- Original observations: 17,894,047  
- Final analytical sample: 3,899,679 women aged 25–54
- 
## Sample Restrictions

Three filters are applied:

1. **Sex**: Female respondents only (`SEX == 2`)
2. **Age**: 25–54 (prime working age; excludes education transitions and early retirement)
3. **Geographic Coverage**: All 50 states + District of Columbia

After restrictions:
- 21.8% of original data retained  
- 3,899,679 observations  
- Aggregated into **561 state-year cells** (51 jurisdictions × 11 years)

## Variable Definitions

### Primary Outcome Variable

**Female Labor Force Participation Rate (FLFP)**  
Proportion of women aged 25–54 either employed or actively seeking employment.

Constructed as:

```python
in_lf = (LABFORCE == 2)
'''
### Labor Force Coding

Where:

- `LABFORCE == 1` → Not in labor force  
- `LABFORCE == 2` → In labor force  

All aggregates use CPS sampling weights (`WTFINL`).

---

## Predictor Variables (Synthetic Control Matching)

### Education

```python
college = (EDUC >= 111)
Predictor Variables (Synthetic Control Matching)

Education

college = (EDUC >= 111)

(EDUC 111–125 = Bachelor's degree or higher)

Marital Status

married = MARST in [1,2]

Age

Continuous (25–54)

Mean age calculated per state-year

Lagged FLFP

Annual averages (2000–2003)

Monthly pre-treatment values

All state-year aggregates are weighted using WTFINL.

Data Quality and Limitations
Strengths

Balanced Panel: All 51 jurisdictions observed 2000–2010

Large Sample Sizes:

Smallest states: ~3,500 observations/year

California: ~28,577 observations/year

National Representativeness via CPS weights

Key Limitations
1. Treatment Timing Ambiguity

California PFL began July 1, 2004 (mid-year).

Current coding:

treated = (state == "California") & (YEAR >= 2004)

This treats 2004 as partially treated.

Sensitivity Analyses:

Exclude 2004 (treat 2005 as first full treatment year)

Monthly-level analysis (future extension)

Observed Anomaly:
FLFP drops from 72.9% (2003) to 70.5% (2004)

Possible explanations:

Early housing slowdown

CPS compositional shifts

Policy endogeneity

Measurement error

2. Great Recession Confounding (2007–2009)

Post-treatment window includes recession period.

Mitigation strategy:

Primary focus: 2004–2007

Robustness checks including/excluding 2008–2010

SDID time weights

Dynamic treatment effects (year-by-year)

3. California’s Structurally Low FLFP

California average FLFP: 72.1%

Other states average: 77.5%

Despite higher education levels

Likely drivers:

High housing costs

Tech-sector spousal income

Demographic composition

Implication:
Donor pool should weight low-FLFP states (e.g., WV, AR, MS, LA, NM)
High-FLFP states should receive near-zero weights.

Methodology
Synthetic Difference-in-Differences (SDID)

Primary method: SDID
Also implemented:

Standard Difference-in-Differences (DID)

Synthetic Control (SC)

Reasons:

Parallel trends assumption implausible

Double robustness (unit + time weights)

Strong simulation performance in literature

Aggregation Decision

Data aggregated to state-year level rather than state-month:

Reasons:

Predictor availability (annual)

Reduced sampling noise

Computational feasibility (561 vs 6,732 cells)

Trade-off:
Reduced precision around July 2004 implementation.

Monthly extension remains future work.


