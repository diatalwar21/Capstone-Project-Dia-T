"""
Synthetic Difference-in-Differences (SDiD) — Arkhangelsky et al. (2021)
California Paid Family Leave (2004) → Female Labor Force Participation

Packages: numpy, scipy, pandas, matplotlib, statsmodels (all standard, peer-reviewed)
Reference: Arkhangelsky, D., Athey, S., Hirshberg, D.A., Imbens, G.W., Wager, S.
           "Synthetic Difference-in-Differences." AER, 2021.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy import stats
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

def p(*args): print(*args, flush=True)
np.random.seed(42)

# ── Parameters ───────────────────────────────────────────────────────────────
DATA_FILE  = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'state_year_predictors.csv')
TREATED    = 'California'
TREAT_YEAR = 2004
OUTCOME    = 'flfp_pct'
EXCLUDED   = ['New Jersey', 'District of Columbia', 'Alaska']  # always-treated / outliers

# ── Load data ─────────────────────────────────────────────────────────────────
agg = pd.read_csv(DATA_FILE)
agg.columns = agg.columns.str.strip()

all_states   = sorted(agg['state'].unique())
donor_states = [s for s in all_states if s != TREATED and s not in EXCLUDED]
years        = sorted(agg['YEAR'].unique())
pre_years    = [y for y in years if y < TREAT_YEAR]
post_years   = [y for y in years if y >= TREAT_YEAR]
T0, T1, T    = len(pre_years), len(post_years), len(years)
N            = len(donor_states)

p(f'Data: {agg.shape[0]} rows | {N} donors | T0={T0} pre | T1={T1} post')

# ── Build outcome matrix Y: rows=states, cols=years ──────────────────────────
Y_mat    = agg.pivot(index='state', columns='YEAR', values=OUTCOME)
Y_mat    = Y_mat.reindex(index=[TREATED] + donor_states, columns=years)
Y_full   = Y_mat.values.astype(float)   # (N+1) × T
Y_treat  = Y_full[0, :]                 # shape (T,)
Y_donors = Y_full[1:, :]                # shape (N, T)

Y_pre_treat   = Y_treat[:T0]            # (T0,)
Y_post_treat  = Y_treat[T0:]            # (T1,)
Y_pre_donors  = Y_donors[:, :T0]        # (N, T0)
Y_post_donors = Y_donors[:, T0:]        # (N, T1)

# ── SDiD Unit Weights ω ───────────────────────────────────────────────────────
# Eq. (4) in Arkhangelsky et al.:
#   min_ω ||Y_pre_donors.T @ ω - Y_pre_treat||² + ζ² ||ω||²
#   s.t.  ω ≥ 0,  Σω = 1
#
# Regularisation ζ: first-difference variance of donor panel (paper recommendation)
sigma = np.std(np.diff(Y_pre_donors, axis=1), ddof=1)          # first-diff std
zeta  = max(float(T1) / T0, 1.0) * sigma                       # penalise overfitting

def omega_obj(w):
    resid = Y_pre_donors.T @ w - Y_pre_treat   # (T0,)
    return np.mean(resid**2) + (zeta**2) * np.dot(w, w)

def omega_grad(w):
    resid = Y_pre_donors.T @ w - Y_pre_treat
    return (2.0 / T0) * (Y_pre_donors @ resid) + 2.0 * (zeta**2) * w

res_omega = minimize(
    omega_obj, np.ones(N) / N, jac=omega_grad, method='SLSQP',
    bounds=[(0, 1)] * N,
    constraints={'type': 'eq', 'fun': lambda w: w.sum() - 1, 'jac': lambda w: np.ones(N)},
    options={'ftol': 1e-12, 'maxiter': 3000}
)
omega = res_omega.x

# ── SDiD Time Weights λ ───────────────────────────────────────────────────────
# Eq. (5) in Arkhangelsky et al.:
#   min_λ ||Y_pre_donors @ λ - mean(Y_post_donors, axis=1)||²
#   s.t.  λ ≥ 0,  Σλ = 1
# These weights make the pre-trend look like an "average" post period
Y_post_donor_mean = Y_post_donors.mean(axis=1)   # (N,)

def lambda_obj(lam):
    resid = Y_pre_donors @ lam - Y_post_donor_mean   # (N,)
    return np.mean(resid**2)

def lambda_grad(lam):
    resid = Y_pre_donors @ lam - Y_post_donor_mean
    return (2.0 / N) * (Y_pre_donors.T @ resid)

res_lambda = minimize(
    lambda_obj, np.ones(T0) / T0, jac=lambda_grad, method='SLSQP',
    bounds=[(0, 1)] * T0,
    constraints={'type': 'eq', 'fun': lambda l: l.sum() - 1, 'jac': lambda l: np.ones(T0)},
    options={'ftol': 1e-12, 'maxiter': 3000}
)
lam = res_lambda.x

# ── SDiD Point Estimate ───────────────────────────────────────────────────────
# 2×2 weighted DiD using ω (unit) and λ (time) weights
# τ_sdid = [Ȳ_treat,post - Ȳ_treat,pre(λ)] - ω·[Ȳ_donors,post - Ȳ_donors,pre(λ)]
Y_treat_post_avg   = Y_post_treat.mean()          # scalar
Y_treat_pre_wt     = float(Y_pre_treat @ lam)     # λ-weighted pre-period: scalar
Y_donors_post_avg  = Y_post_donors.mean(axis=1)   # (N,)
Y_donors_pre_wt    = Y_pre_donors @ lam           # (N,)

tau_sdid = (Y_treat_post_avg - Y_treat_pre_wt) \
         - float(omega @ (Y_donors_post_avg - Y_donors_pre_wt))

p(f'\nSDiD ATT (avg post): {tau_sdid:+.4f} pp')

# ── Year-by-year ATT (Event Study) ───────────────────────────────────────────
# For each year t, ATT_t = (Y_treat,t - Ȳ_treat,pre(λ)) - ω·(Y_donors,t - Ȳ_donors,pre(λ))
att_yearly = np.array([
    (Y_treat[i] - Y_treat_pre_wt)
    - float(omega @ (Y_donors[:, i] - Y_donors_pre_wt))
    for i in range(T)
])

# ── Standard Errors: Leave-One-Unit-Out Jackknife ────────────────────────────
# Jackknife over donor units (Arkhangelsky et al. recommend unit-level resampling)
tau_jack = np.empty(N)
for j in range(N):
    keep      = [i for i in range(N) if i != j]
    Y_don_j   = Y_donors[keep, :]
    Yp_pre_j  = Y_don_j[:, :T0]
    Yp_post_j = Y_don_j[:, T0:]

    # Re-solve omega (unit weights without donor j)
    sigma_j = np.std(np.diff(Yp_pre_j, axis=1), ddof=1)
    zeta_j  = max(float(T1) / T0, 1.0) * sigma_j
    Nj      = len(keep)

    def obj_j(w):
        r = Yp_pre_j.T @ w - Y_pre_treat
        return np.mean(r**2) + (zeta_j**2) * np.dot(w, w)

    def grad_j(w):
        r = Yp_pre_j.T @ w - Y_pre_treat
        return (2.0 / T0) * (Yp_pre_j @ r) + 2.0 * (zeta_j**2) * w

    res_j = minimize(obj_j, np.ones(Nj) / Nj, jac=grad_j, method='SLSQP',
                     bounds=[(0, 1)] * Nj,
                     constraints={'type': 'eq', 'fun': lambda w: w.sum() - 1},
                     options={'ftol': 1e-10, 'maxiter': 1000})
    omega_j = res_j.x

    # Time weights unchanged (don't depend on a single donor much)
    tau_jack[j] = (Y_treat_post_avg - Y_treat_pre_wt) \
                - float(omega_j @ (Yp_post_j.mean(axis=1) - Yp_pre_j @ lam))

se_jack  = np.sqrt(((N - 1) / N) * np.sum((tau_jack - tau_sdid)**2))
t_stat   = tau_sdid / se_jack if se_jack > 0 else np.nan
pval     = 2 * stats.t.sf(abs(t_stat), df=N - 1) if not np.isnan(t_stat) else np.nan
ci_lo    = tau_sdid - 1.96 * se_jack
ci_hi    = tau_sdid + 1.96 * se_jack

p(f'SDiD SE  (jackknife): {se_jack:.4f}')
p(f'SDiD t-statistic    : {t_stat:.3f}')
p(f'SDiD p-value        : {pval:.4f}')
p(f'95% CI              : [{ci_lo:+.4f}, {ci_hi:+.4f}]')

# ── Comparison: SDiD vs plain DiD (equally-weighted) ─────────────────────────
# Plain DiD = same formula but ω = 1/N and λ = 1/T0
omega_ew = np.ones(N) / N
lam_ew   = np.ones(T0) / T0
tau_did  = (Y_treat_post_avg - float(Y_pre_treat @ lam_ew)) \
         - float(omega_ew @ (Y_post_donors.mean(axis=1) - Y_pre_donors @ lam_ew))

p(f'\nPlain DiD ATT (equal weights): {tau_did:+.4f} pp')
p(f'SDiD  ATT (Arkhangelsky wts) : {tau_sdid:+.4f} pp')

# ── Weight diagnostics ────────────────────────────────────────────────────────
omega_df = pd.DataFrame({'State': donor_states, 'omega': omega})
omega_df = omega_df[omega_df['omega'] > 0.005].sort_values('omega', ascending=False)

lam_df = pd.DataFrame({'Year': pre_years, 'lambda': lam})

p('\nUNIT WEIGHTS (ω > 0.005)')
p(omega_df.to_string(index=False, float_format='{:.4f}'.format))

p('\nTIME WEIGHTS (λ)')
p(lam_df.to_string(index=False, float_format='{:.4f}'.format))

# ── Pre-trend fit check ───────────────────────────────────────────────────────
synth_pre = Y_pre_donors.T @ omega        # synthetic control pre-period
pre_bias  = Y_pre_treat - synth_pre
pre_rmse  = np.sqrt(np.mean(pre_bias**2))
p(f'\nPre-trend RMSE (ω fit): {pre_rmse:.4f} pp')

# ── Placebo inference (in-space) ──────────────────────────────────────────────
p('\nRunning in-space placebo SDiD (leave-one-out)...')
placebo_taus = []

for p_idx, p_state in enumerate(donor_states):
    p_donors = [s for s in donor_states if s != p_state]
    Np = len(p_donors)
    p_Y_treat  = Y_donors[donor_states.index(p_state), :]
    p_Y_donors = np.array([Y_donors[donor_states.index(s), :] for s in p_donors])

    p_pre  = p_Y_treat[:T0]
    p_post = p_Y_treat[T0:]
    p_Dpre = p_Y_donors[:, :T0]
    p_Dpost= p_Y_donors[:, T0:]

    sig_p  = np.std(np.diff(p_Dpre, axis=1), ddof=1)
    zeta_p = max(float(T1) / T0, 1.0) * sig_p

    def pobj(w):
        r = p_Dpre.T @ w - p_pre
        return np.mean(r**2) + (zeta_p**2) * np.dot(w, w)

    res_p = minimize(pobj, np.ones(Np) / Np, method='SLSQP',
                     bounds=[(0, 1)] * Np,
                     constraints={'type': 'eq', 'fun': lambda w: w.sum() - 1},
                     options={'ftol': 1e-10, 'maxiter': 500})
    w_p = res_p.x

    tau_p = (p_post.mean() - float(p_pre @ lam)) \
          - float(w_p @ (p_Dpost.mean(axis=1) - p_Dpre @ lam))
    placebo_taus.append(tau_p)

    if (p_idx + 1) % 10 == 0:
        p(f'  {p_idx+1}/{N} placebos done...')

placebo_taus = np.array(placebo_taus)
pval_placebo = np.mean(np.abs(placebo_taus) >= abs(tau_sdid))
p(f'\nPermutation p-value (|τ_placebo| ≥ |τ_sdid|): {pval_placebo:.4f}')

# ── Results table ─────────────────────────────────────────────────────────────
results_df = pd.DataFrame({
    'Year': years,
    'CA_actual': Y_treat,
    'Synth_CA':  Y_pre_donors.T @ omega if False else  # placeholder
                 np.concatenate([Y_pre_donors.T @ omega,
                                 Y_post_donors.T @ omega]),
    'ATT_yearly': att_yearly,
    'Post': [int(y >= TREAT_YEAR) for y in years]
})
results_df['Trend_adj'] = results_df['CA_actual'] - results_df['Synth_CA']

p('\nYEAR-BY-YEAR SDiD ATT')
p(results_df[['Year','CA_actual','Synth_CA','ATT_yearly','Post']].to_string(
    index=False, float_format='{:.3f}'.format))

# ── Plots ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(21, 7))

# Plot 1: Actual vs ω-weighted synthetic control
synth_all = np.concatenate([Y_pre_donors.T @ omega, Y_post_donors.T @ omega])
ax = axes[0]
ax.plot(years, Y_treat,   'o-',  color='crimson', lw=2.5, ms=8, label='California (Actual)')
ax.plot(years, synth_all, 's--', color='navy',    lw=2.5, ms=8, label='Synthetic CA (SDiD ω)')
ax.axvspan(TREAT_YEAR - 0.5, max(years) + 0.5, alpha=0.07, color='steelblue')
ax.axvline(TREAT_YEAR - 0.5, color='steelblue', lw=1.5, ls=':')
ax.set_title('SDiD: Actual vs Synthetic CA', fontweight='bold')
ax.set_xlabel('Year', fontweight='bold')
ax.set_ylabel('Female LFP Rate (%)', fontweight='bold')
ax.legend(); ax.grid(alpha=0.3); ax.set_xticks(years)

# Plot 2: Year-by-year ATT (event study)
ax = axes[1]
pre_mask  = np.array([y < TREAT_YEAR for y in years])
post_mask = ~pre_mask
colors    = ['steelblue' if y < TREAT_YEAR else ('forestgreen' if g >= 0 else 'tomato')
             for y, g in zip(years, att_yearly)]
ax.bar(years, att_yearly, color=colors, edgecolor='black', lw=0.5, alpha=0.85, width=0.7)
ax.axhline(0, color='black', lw=1.5)
ax.axvline(TREAT_YEAR - 0.5, color='steelblue', lw=1.5, ls=':')
ax.axhline(tau_sdid, color='crimson', lw=2, ls='--', label=f'Avg ATT = {tau_sdid:+.2f} pp')
for y, g in zip(years, att_yearly):
    ax.text(y, g + (0.1 if g >= 0 else -0.3), f'{g:+.2f}',
            ha='center', fontsize=8, fontweight='bold')
ax.set_title(f'SDiD Event Study  |  ATT={tau_sdid:+.2f} pp  p={pval:.3f}',
             fontweight='bold')
ax.set_xlabel('Year', fontweight='bold')
ax.set_ylabel('ATT (pp)', fontweight='bold')
ax.legend(); ax.grid(alpha=0.3, axis='y'); ax.set_xticks(years)

# Plot 3: In-space placebo distribution
ax = axes[2]
ax.hist(placebo_taus, bins=15, color='steelblue', edgecolor='black', alpha=0.75,
        label='Placebo ATTs')
ax.axvline(tau_sdid, color='crimson', lw=2.5, ls='--',
           label=f'CA SDiD ATT = {tau_sdid:+.2f}\np={pval_placebo:.3f}')
ax.set_title('In-Space Placebo Distribution', fontweight='bold')
ax.set_xlabel('Placebo ATT (pp)', fontweight='bold')
ax.set_ylabel('Count', fontweight='bold')
ax.legend(); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('sdid_results.png', dpi=150, bbox_inches='tight')
p('Saved: sdid_results.png')
plt.close()

# ── Save results ──────────────────────────────────────────────────────────────
results_df.to_csv('sdid_yearly_att.csv', index=False)
omega_df.to_csv('sdid_unit_weights.csv', index=False)
lam_df.to_csv('sdid_time_weights.csv', index=False)

pd.DataFrame({
    'state': donor_states,
    'placebo_att': placebo_taus
}).to_csv('sdid_placebo_atts.csv', index=False)

# ── Final summary ─────────────────────────────────────────────────────────────
p('\n' + '=' * 65)
p('  SYNTHETIC DiD — FINAL SUMMARY')
p('=' * 65)
p(f'  Estimator            : SDiD (Arkhangelsky et al. 2021)')
p(f'  Treated unit         : {TREATED}')
p(f'  Treatment year       : {TREAT_YEAR} (CA Paid Family Leave)')
p(f'  Outcome              : Female LFP rate (flfp_pct, %)')
p(f'  Donor pool           : {N} states')
p()
p(f'  SDiD ATT (avg)       : {tau_sdid:+.4f} pp')
p(f'  Jackknife SE         : {se_jack:.4f}')
p(f'  95% CI               : [{ci_lo:+.4f}, {ci_hi:+.4f}]')
p(f'  t-stat               : {t_stat:.3f}')
p(f'  Jackknife p-value    : {pval:.4f}')
p(f'  Permutation p-value  : {pval_placebo:.4f}')
p(f'  Pre-trend RMSE       : {pre_rmse:.4f} pp')
p()
p(f'  Plain DiD (EW) ATT   : {tau_did:+.4f} pp  (benchmark)')
p()
if pval <= 0.05 or pval_placebo <= 0.05:
    p('  → SIGNIFICANT at 5%')
elif pval <= 0.10 or pval_placebo <= 0.10:
    p('  → Significant at 10%')
else:
    p('  → Not significant at conventional levels')

p('\nFiles: sdid_yearly_att.csv | sdid_unit_weights.csv | sdid_time_weights.csv | sdid_placebo_atts.csv')
p('Plot : sdid_results.png')
