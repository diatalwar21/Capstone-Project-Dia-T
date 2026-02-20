# Capstone-Project-Dia-T
 
## California Paid Family Leave and Female Labor Force Participation Analysis

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

in_lf = (LABFORCE == 2)

### Labor Force Coding

Where:

- `LABFORCE == 1` → Not in labor force  
- `LABFORCE == 2` → In labor force  

All aggregates use CPS sampling weights (`WTFINL`).

---

## Predictor Variables (Synthetic Control Matching)

### Education

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



