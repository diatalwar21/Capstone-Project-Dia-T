# Capstone-Project-Dia-T
California's Paid Family Leave and Female Labor Force Participation 

\documentclass[12pt]{article}
\usepackage{setspace}
\usepackage{geometry}
\usepackage{amsmath}
\usepackage{booktabs}
\usepackage{enumitem}
\usepackage{hyperref}

\geometry{margin=1in}
\onehalfspacing

\title{California Paid Family Leave and Female Labor Force Participation Analysis}
\date{}

\begin{document}

\maketitle

\section{Data Source and Sample Construction}

\textbf{Primary Data Source.} 
This analysis uses monthly microdata from the \textit{Current Population Survey (CPS)}, accessed through the \textit{Integrated Public Use Microdata Series (IPUMS-CPS)}. The CPS is a monthly household survey conducted by the U.S. Census Bureau and the Bureau of Labor Statistics and represents the primary source of U.S. labor force statistics.

\textbf{Sample Period.} January 2000 through December 2010 (132 months).

\begin{itemize}
    \item Original observations: 17,894,047
    \item Final analytical sample: 3,899,679 women aged 25--54
\end{itemize}

\subsection*{Sample Restrictions}

Three filters are applied to construct the analytical sample:

\begin{enumerate}
    \item \textbf{Sex:} Female respondents only ($\texttt{SEX} = 2$)
    \item \textbf{Age:} Restrict to ages 25--54 to focus on prime working-age women and avoid confounding from education transitions (under 25) and early retirement decisions (over 54)
    \item \textbf{Geographic Coverage:} All 50 states plus the District of Columbia
\end{enumerate}

These restrictions reduce the sample to 21.8\% of the original data, yielding 3,899,679 individual-level observations aggregated into 561 state-year cells (51 jurisdictions $\times$ 11 years).

\section{Variable Definitions}

\subsection*{Primary Outcome Variable}

\textbf{Female Labor Force Participation Rate (FLFP):}  
Proportion of women aged 25--54 who are either employed or actively seeking employment, calculated using CPS sampling weights (\texttt{WTFINL}).

Constructed as:
\[
\texttt{in\_lf} = (\texttt{LABFORCE} = 2)
\]
where \texttt{LABFORCE} codes:
\begin{itemize}
    \item 1 = Not in labor force
    \item 2 = In labor force
\end{itemize}

\subsection*{Predictor Variables for Synthetic Control Matching}

\begin{itemize}
    \item \textbf{Education:} College degree indicator constructed as
    \[
    \texttt{college} = (\texttt{EDUC} \geq 111)
    \]
    where \texttt{EDUC} values 111--125 correspond to Bachelor's degree or higher in IPUMS coding.

    \item \textbf{Marital Status:} Married indicator defined as
    \[
    \texttt{married} = \texttt{MARST} \in \{1,2\}
    \]
    capturing women with spouse present or absent.

    \item \textbf{Age:} Continuous variable (25--54); mean age calculated per state-year.

    \item \textbf{Lagged Labor Force Participation:} Annual averages for 2000--2003 and monthly values for the pre-treatment period.
\end{itemize}

All state-year aggregates are calculated using CPS person-level weights (\texttt{WTFINL}) to ensure population-representative estimates.

\section{Data Quality and Limitations}

\subsection*{Strengths}

\begin{itemize}
    \item \textbf{Complete Temporal Coverage:} All 51 jurisdictions have data for all 11 years (2000--2010), yielding a perfectly balanced panel with no missing state-year cells.
    
    \item \textbf{Adequate Sample Sizes:} Even the smallest states average 3,500+ observations per year. California averages 28,577 observations per year, providing highly precise FLFP estimates.
    
    \item \textbf{National Representativeness:} CPS sampling weights account for demographic oversampling and non-response, making aggregates population-representative.
\end{itemize}

\subsection*{Key Limitations}

\subsubsection*{Treatment Timing Ambiguity}

California's Paid Family Leave (PFL) became effective July 1, 2004, creating mid-year implementation. Although CPS includes a \texttt{MONTH} variable, the current analysis uses annual aggregates for consistency with state-year predictors. This introduces measurement error in treatment timing.

\textbf{Treatment Coding Decision:}
\[
\texttt{treated} = (\texttt{state = California} \ \& \ \texttt{YEAR} \geq 2004)
\]

Thus, all of 2004 is treated as partially treated.

\textbf{Sensitivity Analyses:}
\begin{enumerate}
    \item Exclude 2004 entirely (treat 2005 as first full treatment year)
    \item Conduct monthly-level analysis if computationally feasible
\end{enumerate}

\textbf{Observed Anomaly:} California's FLFP declines from 72.9\% (2003) to 70.5\% (2004), opposite the expected direction if PFL increases labor force attachment. Potential explanations include:

\begin{itemize}
    \item Early recession effects in California’s housing market
    \item CPS compositional changes
    \item Policy endogeneity
    \item Measurement error
\end{itemize}

This anomaly necessitates careful pre-treatment fit diagnostics and dynamic treatment effect estimation.

\subsubsection*{Great Recession Confounding}

The post-treatment period (2004--2010) includes the Great Recession (2007--2009), which disproportionately affected California due to the housing market collapse.

National female FLFP declined from 77.5\% (2000) to 77.0\% (2010), with sharper drops during 2008--2009.

\textbf{Mitigation Strategy:}

\begin{itemize}
    \item Primary estimates focus on 2004--2007 (pre-recession)
    \item Robustness checks include/exclude 2008--2010
    \item SDID time weights down-weight poorly predicted periods
    \item Dynamic treatment effects estimated year-by-year
\end{itemize}

\subsubsection*{California's Structural Low FLFP}

California’s average FLFP (72.1\%) is 5.4 percentage points lower than other states (77.5\%), despite higher educational attainment (31.9\% vs.\ 30.9\% college-educated).

Likely explanations include high housing costs, dual-income substitution patterns, tech-sector spousal earnings, and demographic composition.

\textbf{Implication for SDID:} Donor states should resemble California’s low FLFP profile (e.g., West Virginia, Arkansas, Mississippi, Louisiana, New Mexico). High-FLFP states are poor counterfactuals and should receive near-zero weights. This will be verified through unit-weight inspection and leave-one-out robustness checks.

\section{Methodological Decisions and Justifications}

\subsection*{Synthetic Difference-in-Differences (SDID)}

The primary method is Synthetic Difference-in-Differences (SDID), supplemented by standard Difference-in-Differences (DID) and Synthetic Control (SC).

Justifications:

\begin{itemize}
    \item \textbf{Parallel Trends Implausible:} California’s unique economic structure makes parallel trends assumptions tenuous.
    \item \textbf{Double Robustness:} SDID reduces bias if either:
    \begin{enumerate}
        \item Unit weights balance characteristics, or
        \item Time weights predict post-treatment outcomes
    \end{enumerate}
    \item \textbf{Empirical Performance:} Prior literature shows SDID outperforms DID and SC in simulations and policy applications.
\end{itemize}

\subsection*{Aggregation to State-Year Level}

Individual CPS data are aggregated to state-year level rather than state-month for three reasons:

\begin{itemize}
    \item Predictor variables available annually
    \item Reduced sampling noise
    \item Computational tractability (561 vs.\ 6,732 cells)
\end{itemize}

\textbf{Trade-Off:} Annual aggregation reduces precision regarding July 2004 treatment timing. Monthly analysis remains a future extension.

\end{document}
