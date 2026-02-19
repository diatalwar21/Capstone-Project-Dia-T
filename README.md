# Capstone-Project-Dia-T
California's Paid Family Leave and Female Labor Force Participation 

# Project Proposal: California's Paid Family Leave and Female Labor Force Participation 

## Research Question: 
What is the causal impact of California's Paid Family Leave program on labor force participation rates among women aged 25-54?
## My Motivation: 
California's 2004 Paid Family Leave (PFL) program was the first comprehensive state-level paid leave policy in the United States, creating a critical natural experiment for understanding whether government-mandated leave can help women balance career and caregiving.
### Policy relevance: 
Understanding whether PFL increases female labor force participation (FLFP) is essential for evidence-based policymaking. If effective, it provides a proven model for reducing the "motherhood penalty"; if not, it suggests paid leave alone may be insufficient without complementary policies.
### Economic relevance: 
Women's labor supply affects GDP growth, household income inequality, and social insurance sustainability. 
### Academic relevance: 
Prior research on California's PFL yields mixed findings due to methodological challenges in constructing valid counterfactuals for California's unique economy and demographics. My project uses synthetic diff-in-diff methods to address these limitations.
## Expected Contribution
I will apply the new synthetic difference-in-differences (SDID) method (Arkhangelsky et al. 2021) to construct a data-driven counterfactual for California. SDID combines synthetic control unit weighting with difference-in-differences regression and time weights, providing double robustness not available in standard approaches. This directly addresses the challenge that California's distinctive characteristics make finding valid comparison states difficult.

## Primary Data Source: IPUMS Current Population Survey (CPS) monthly microdata, January 2000 through December 2010.
### Outcome Variable: 
Female labor force participation rate (FLFP) for women aged 25-54 in state s at month t, calculated using CPS sampling weights. This age range avoids confounding from education decisions (under 25) and retirement (over 54).
Treatment Variable: Binary indicator equal to 1 for California observations after July 2004, 0 otherwise.
Predictor Variables for Synthetic Control Matching:
Lagged outcomes (key): FLFP rates at multiple pre-treatment points (annual means for 2000, 2001, 2002, 2003; monthly values for 6-12 months immediately pre-treatment)
Education: Mean years of schooling for women / share with college degree
Demographics: Average female age, age distribution, marriage rate among women
Economic indicators (potential): Industry composition, unemployment rate, median wages
Donor Pool: Approximately 40 control states after exclusions:
Exclude: States that implemented PFL during study period (New Jersey 2009, Washington 2007)
Exclude: States with Temporary Disability Insurance that may contaminate controls (New York, Rhode Island, Hawaii)
Consider excluding: Outliers with fundamentally different labor markets (Alaska)
Feasibility: CPS microdata is publicly available through IPUMS and widely used for state-level labor market analysis. I have it downloaded.
## Initial Empirical Strategy
I will implement a synthetic difference-in-differences (SDID) approach to estimate the causal effect of California's PFL on FLFP.
### Anticipated challenges:
The Great Recession (2008-2009)
Gradual implementation effects
Pre-treatment match quality
