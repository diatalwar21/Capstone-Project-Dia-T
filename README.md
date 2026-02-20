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

## Predictor variables for synthetic control matching:
#### Education: 
College degree indicator constructed as college = (EDUC >= 111), where EDUC values 111-125 correspond to Bachelor's degree or higher in IPUMS coding
### Marital status: 
Married indicator defined as married = MARST in [1,2], capturing women with spouse present or absent
### Age: 
Continuous variable (25-54) and calculated mean age per state-year
### Labor force participation (lagged): 
Annual averages for 2000-2003 and monthly values for pre-treatment period
All state-year aggregates are calculated using CPS person-level weights (WTFINL) to ensure population-representative estimates.

## Data Quality and Limitations
### Strengths
Complete temporal coverage: All 51 jurisdictions have data for all 11 years (2000-2010), yielding a perfectly balanced panel with no missing state-year cells.
Adequate sample sizes: Even the smallest states (Mississippi, Montana) average 3,500+ observations per year, sufficient for reliable state-year FLFP estimates. California, with 28,577 observations per year on average, provides exceptionally precise estimates.
National representativeness: CPS sampling weights account for demographic oversampling and non-response, making state-year aggregates population-representative.

