# Capstone-Project-Dia-T
 
## California Paid Family Leave and Female Labor Force Participation Analysis
### Data Source and Sample Construction
#### Primary Data Source : My analysis uses monthly microdata from the Current Population Survey (CPS), accessed through the Integrated Public Use Microdata Series (IPUMS-CPS). The CPS is a monthly household survey conducted by the U.S. Census Bureau and Bureau of Labor Statistics, representing the primary source for U.S. labor force statistics.
Sample Period: January 2000 through December 2010 (132 months)
 Original observations: 17,894,047
 Final analytical sample: 3,899,679 women aged 25-54
Sample Restrictions
I apply three filters to construct the analytical sample:
Sex: Female respondents only (SEX == 2)
Age: Restrict to ages 25-54 to focus on prime working-age women and avoid confounding from education transitions (under 25) and early retirement decisions (over 54)
Geographic coverage: All 50 states plus District of Columbia
These restrictions reduce the sample to 21.8% of the original data, yielding 3,899,679 individual-level observations aggregated into 561 state-year cells (51 jurisdictions × 11 years).
Variable Definitions
Primary outcome variable:
Female Labor Force Participation Rate (FLFP): Proportion of women aged 25-54 who are either employed or actively seeking employment, calculated using CPS sampling weights (WTFINL). Constructed as: in_lf = (LABFORCE == 2), where LABFORCE codes 1=Not in labor force, 2=In labor force.
