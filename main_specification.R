# =============================================================================
# MAIN SPECIFICATION + ROBUSTNESS CHECK
# Effect of California Paid Family Leave (2004) on Female Labor Force Participation
# Estimator: Arkhangelsky et al. (AER 2021) via synthdid package
#
# ROBUSTNESS: Exclude 2004-2005 post-treatment window
# Motivation: The 2005 FLFP dip in CA is a confound (housing bubble + immigration
# composition), not the PFL treatment effect. Treating 2004-2005 as additional
# pre-periods and estimating the ATT only over 2006-2010 isolates the sustained
# medium-term PFL effect from the transitory 2005 shock.
# =============================================================================

# ── Block 1: Packages ─────────────────────────────────────────────────────────
# synthdid: official package by Arkhangelsky et al. — SDiD, SC, DiD estimators
# dplyr:    data cleaning only

library(synthdid)
library(dplyr)

# ── Block 2: Load data ────────────────────────────────────────────────────────
# Drop NJ, DC, and Alaska — they passed similar paid leave laws near 2004
# and would contaminate the counterfactual if left in the donor pool.

df <- read.csv("state_year_predictors.csv") %>%
  filter(!state %in% c("New Jersey", "District of Columbia", "Alaska"))

# ── Block 3: Main specification (post = 2004-2010) ───────────────────────────
# panel.matrices() converts the long panel into Y (outcome matrix), N0 (number
# of donor states), and T0 (number of pre-treatment periods).
# The 'treated' column in the data is already: 1 for CA from 2004 onward, 0 else.

setup_main <- panel.matrices(
  df,
  unit      = "state",
  time      = "YEAR",
  outcome   = "flfp_pct",
  treatment = "treated"
)

tau_did  <- did_estimate(setup_main$Y, setup_main$N0, setup_main$T0)
tau_sc   <- sc_estimate(setup_main$Y,  setup_main$N0, setup_main$T0)
tau_sdid <- synthdid_estimate(setup_main$Y, setup_main$N0, setup_main$T0)

cat("MAIN SPEC — pre: 2000-2003 | post: 2004-2010\n")
cat("  DiD  ATT:", round(as.numeric(tau_did),  4), "pp\n")
cat("  SC   ATT:", round(as.numeric(tau_sc),   4), "pp\n")
cat("  SDiD ATT:", round(as.numeric(tau_sdid), 4), "pp  <- main\n")

# vcov() returns a 1×1 variance matrix; [1,1] extracts the scalar variance.
# We use the "placebo" method — in-space permutation inference — because with
# only T0=4 pre-periods the jackknife over time is unreliable.
se_did  <- sqrt(vcov(tau_did,  method = "placebo")[1,1])
se_sc   <- sqrt(vcov(tau_sc,   method = "placebo")[1,1])
se_sdid <- sqrt(vcov(tau_sdid, method = "placebo")[1,1])

ci_lo  <- as.numeric(tau_sdid) - 1.96 * se_sdid
ci_hi  <- as.numeric(tau_sdid) + 1.96 * se_sdid
t_stat <- as.numeric(tau_sdid) / se_sdid
pval   <- 2 * pt(abs(t_stat), df = setup_main$N0 - 1, lower.tail = FALSE)

cat("  SE:", round(se_sdid, 4), "| 95% CI: [", round(ci_lo,4), ",", round(ci_hi,4), "] | p:", round(pval,4), "\n")

# ── Block 4: ROBUSTNESS CHECK — exclude 2004-2005 ────────────────────────────
# Strategy: recode treated = 1 for California only from 2006 onward.
# 2004 and 2005 become additional pre-treatment periods (T0 grows from 4 to 6).
# SDiD then estimates the ATT only over 2006-2010, after the housing/immigration
# confound has dissipated. If the robustness ATT > main ATT, it confirms the
# 2005 dip was masking a real positive PFL effect.

df_robust <- df %>%
  mutate(treated_robust = ifelse(state == "California" & YEAR >= 2006, 1, 0))

setup_robust <- panel.matrices(
  df_robust,
  unit      = "state",
  time      = "YEAR",
  outcome   = "flfp_pct",
  treatment = "treated_robust"
)

tau_did_r  <- did_estimate(setup_robust$Y, setup_robust$N0, setup_robust$T0)
tau_sc_r   <- sc_estimate(setup_robust$Y,  setup_robust$N0, setup_robust$T0)
tau_sdid_r <- synthdid_estimate(setup_robust$Y, setup_robust$N0, setup_robust$T0)

cat("\nROBUSTNESS — pre: 2000-2005 | post: 2006-2010 (excl. 2004-2005 confound)\n")
cat("  DiD  ATT:", round(as.numeric(tau_did_r),  4), "pp\n")
cat("  SC   ATT:", round(as.numeric(tau_sc_r),   4), "pp\n")
cat("  SDiD ATT:", round(as.numeric(tau_sdid_r), 4), "pp  <- robustness\n")

se_sdid_r <- sqrt(vcov(tau_sdid_r, method = "placebo")[1,1])
ci_lo_r   <- as.numeric(tau_sdid_r) - 1.96 * se_sdid_r
ci_hi_r   <- as.numeric(tau_sdid_r) + 1.96 * se_sdid_r
t_stat_r  <- as.numeric(tau_sdid_r) / se_sdid_r
pval_r    <- 2 * pt(abs(t_stat_r), df = setup_robust$N0 - 1, lower.tail = FALSE)

cat("  SE:", round(se_sdid_r, 4), "| 95% CI: [", round(ci_lo_r,4), ",", round(ci_hi_r,4), "] | p:", round(pval_r,4), "\n")

# ── Block 5: Side-by-side comparison table ────────────────────────────────────
# This table is what you'd put in your paper as Table 1.
# It shows how much the estimate changes once the confound is removed.

results <- data.frame(
  Specification = c("DiD",
                    "Synthetic Control",
                    "SDiD — Main (2004-2010)",
                    "SDiD — Robust (2006-2010, excl. confound)"),
  ATT_pp        = round(c(as.numeric(tau_did),
                           as.numeric(tau_sc),
                           as.numeric(tau_sdid),
                           as.numeric(tau_sdid_r)), 4),
  SE            = round(c(se_did, se_sc, se_sdid, se_sdid_r), 4)
)
results$CI_lo <- round(results$ATT_pp - 1.96 * results$SE, 4)
results$CI_hi <- round(results$ATT_pp + 1.96 * results$SE, 4)

cat("\n--- FULL COMPARISON TABLE ---\n")
print(results, row.names = FALSE)

# ── Block 6: Unit and time weights for both specs ─────────────────────────────
# We extract and print the weights so you can see:
# (a) which states form synthetic CA in each spec
# (b) which pre-periods get the most weight (λ)
# For the robustness spec, 2004-2005 are now pre-periods so λ covers 2000-2005.

extract_weights <- function(tau_est, setup, label) {
  w     <- attr(tau_est, "weights")
  omega <- w$omega
  lam   <- w$lambda
  donor_names <- rownames(setup$Y)[1:setup$N0]
  pre_years   <- colnames(setup$Y)[1:setup$T0]
  omega_df <- data.frame(State = donor_names, omega = round(omega, 4))
  omega_df <- omega_df[omega_df$omega > 0.005, ]
  omega_df <- omega_df[order(-omega_df$omega), ]
  lam_df   <- data.frame(Year = pre_years, lambda = round(lam, 4))
  cat("\nUNIT WEIGHTS (ω) —", label, "\n"); print(omega_df, row.names = FALSE)
  cat("TIME WEIGHTS (λ) —", label, "\n");  print(lam_df,   row.names = FALSE)
  list(omega = omega_df, lambda = lam_df)
}

w_main   <- extract_weights(tau_sdid,   setup_main,   "Main spec (2004-2010)")
w_robust <- extract_weights(tau_sdid_r, setup_robust, "Robustness (2006-2010)")

# ── Block 7: Plots ────────────────────────────────────────────────────────────
# synthdid_plot() shows actual CA vs the synthetic counterfactual over time.
# The shaded post-treatment region is where the ATT is measured.
# We make one plot per specification so you can see the confound visually.

png("sdid_main_spec_trend.png", width = 900, height = 500, res = 120)
synthdid_plot(tau_sdid,
              treated.name = "California",
              control.name = "Synthetic California")
dev.off()
cat("\nSaved: sdid_main_spec_trend.png\n")

png("sdid_robust_trend.png", width = 900, height = 500, res = 120)
synthdid_plot(tau_sdid_r,
              treated.name = "California",
              control.name = "Synthetic California (Robust)")
dev.off()
cat("Saved: sdid_robust_trend.png\n")

png("sdid_main_spec_weights.png", width = 700, height = 500, res = 120)
synthdid_units_plot(tau_sdid, negligible.threshold = 0.001)
dev.off()
cat("Saved: sdid_main_spec_weights.png\n")

png("sdid_robust_weights.png", width = 700, height = 500, res = 120)
synthdid_units_plot(tau_sdid_r, negligible.threshold = 0.001)
dev.off()
cat("Saved: sdid_robust_weights.png\n")

# ── Block 8: Save CSVs ────────────────────────────────────────────────────────
write.csv(results,          "sdid_comparison_table.csv",     row.names = FALSE)
write.csv(w_main$omega,     "sdid_main_spec_omega.csv",      row.names = FALSE)
write.csv(w_main$lambda,    "sdid_main_spec_lambda.csv",     row.names = FALSE)
write.csv(w_robust$omega,   "sdid_robust_omega.csv",         row.names = FALSE)
write.csv(w_robust$lambda,  "sdid_robust_lambda.csv",        row.names = FALSE)
cat("Saved: sdid_comparison_table.csv + weight CSVs\n")

# ── Final summary ─────────────────────────────────────────────────────────────
cat("\n================================================================\n")
cat("  RESULTS SUMMARY\n")
cat("================================================================\n")
cat("  MAIN SPEC   SDiD ATT:", round(as.numeric(tau_sdid),   4), "pp",
    "| SE:", round(se_sdid,   4), "| p:", round(pval,   4), "\n")
cat("  ROBUSTNESS  SDiD ATT:", round(as.numeric(tau_sdid_r), 4), "pp",
    "| SE:", round(se_sdid_r, 4), "| p:", round(pval_r, 4), "\n")
cat("----------------------------------------------------------------\n")
delta <- as.numeric(tau_sdid_r) - as.numeric(tau_sdid)
cat("  Difference (robust - main):", round(delta, 4), "pp\n")
if (delta > 0) {
  cat("  -> Robustness ATT is larger. The 2004-2005 confound was\n")
  cat("     masking a positive PFL effect, as hypothesized.\n")
} else {
  cat("  -> Robustness ATT is not larger. Investigate further.\n")
}
cat("================================================================\n")
