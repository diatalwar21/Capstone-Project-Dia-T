"""
Per-year confidence intervals for the SDiD event study
California Paid Family Leave (2004) -> Female LFP, Women 25-54

Method 1: Leave-one-unit-out jackknife SE per year
  tau_jack_yr[j, t]: ATT in year t from the fold that drops donor j
  SE_t = sqrt((N-1)/N * sum_j (tau_jack_yr[j,t] - mean_j tau_jack_yr[:,t])^2)
  CI  = att_t +/- 1.96 * SE_t

Method 2: In-space placebo null-distribution bounds per year
  tau_plac_yr[p, t]: ATT in year t when donor p is treated as pseudo-CA
  Bounds = [2.5th percentile, 97.5th percentile] of tau_plac_yr[:, t]
  Main-sample lam used for placebo (identical convention to run_sdid_robust_restricted.py).

Spec: 43-state donor pool, post window 2007-2015 (T1=9, T0=9).
Zeta conventions identical to run_sdid_robust_restricted.py.

Outputs:
  outputs/tables/sdid_jackknife_per_year_ci.csv
  outputs/tables/sdid_permutation_per_year_ci.csv
  outputs/figures/pres_event_study_jackknife_ci.png/.pdf   (300 DPI)
  outputs/figures/pres_event_study_permutation_ci.png/.pdf (300 DPI)
  outputs/figures/pres_event_study_line.png/.pdf           (300 DPI)
"""

import os
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

DIR       = os.path.dirname(os.path.abspath(__file__))
ROOT      = os.path.abspath(os.path.join(DIR, ".."))
FIG_DIR   = os.path.join(ROOT, 'outputs', 'figures')
TABLE_DIR = os.path.join(ROOT, 'outputs', 'tables')
os.makedirs(FIG_DIR,   exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

def p(*args): print(*args, flush=True)
np.random.seed(42)
t_start = time.time()

# ── Parameters ────────────────────────────────────────────────────────────────
DATA_FILE  = os.path.join(ROOT, 'data', 'derived', 'state_year_panel_deduped_1995_2015.csv')
TREATED    = 'California'
TREAT_YEAR = 2004
POST_START = 2007
OUTCOME    = 'flfp_pct'

EXCLUDED = [
    'New York',             # TDI state
    'Rhode Island',         # TDI state
    'Hawaii',               # TDI state
    'New Jersey',           # own PFL law effective 2008
    'Washington',           # own PFL law effective 2007
    'Alaska',               # structural outlier
    'District of Columbia', # city-state
]

# ── Load data ─────────────────────────────────────────────────────────────────
agg = pd.read_csv(DATA_FILE)
agg.columns = agg.columns.str.strip()

all_states   = sorted(agg['state'].unique())
donor_states = [s for s in all_states if s != TREATED and s not in EXCLUDED]
years        = sorted(agg['YEAR'].unique())
N            = len(donor_states)

years_arr = np.array(years)
pre_idx   = np.where(years_arr < TREAT_YEAR)[0]    # 1995-2003
post_idx  = np.where(years_arr >= POST_START)[0]   # 2007-2015
pre_years  = [years[i] for i in pre_idx]
post_years = [years[i] for i in post_idx]
T0 = len(pre_idx)
T1 = len(post_idx)

p(f'{N} donors | T0={T0} pre ({pre_years[0]}-{pre_years[-1]}) | '
  f'T1={T1} post ({post_years[0]}-{post_years[-1]})')

# ── Build outcome matrix ──────────────────────────────────────────────────────
Y_mat    = agg.pivot(index='state', columns='YEAR', values=OUTCOME)
Y_mat    = Y_mat.reindex(index=[TREATED] + donor_states, columns=years)
Y_full   = Y_mat.values.astype(float)
Y_treat  = Y_full[0, :]
Y_donors = Y_full[1:, :]

Y_pre_treat   = Y_treat[pre_idx]
Y_pre_donors  = Y_donors[:, pre_idx]
Y_post_donors = Y_donors[:, post_idx]

# ── Main estimation: unit weights omega ───────────────────────────────────────
# zeta convention: (N * T1)^0.25 * sigma  [Arkhangelsky et al. main spec]
sigma = np.std(np.diff(Y_pre_donors, axis=1), ddof=1)
zeta  = (N * T1) ** 0.25 * sigma
p(f'zeta = {zeta:.4f}  [N={N}, T1={T1}]')

def omega_obj(w):
    resid = Y_pre_donors.T @ w - Y_pre_treat
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

# ── Main estimation: time weights lambda ──────────────────────────────────────
Y_post_donor_mean = Y_post_donors.mean(axis=1)

def lambda_obj(lam):
    resid = Y_pre_donors @ lam - Y_post_donor_mean
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

# ── Main per-year ATTs ────────────────────────────────────────────────────────
Y_treat_pre_wt  = float(Y_pre_treat @ lam)
Y_donors_pre_wt = Y_pre_donors @ lam      # (N,) — each donor's lambda-weighted pre average

att_yearly = np.array([
    (Y_treat[post_idx[t]] - Y_treat_pre_wt)
    - float(omega @ (Y_donors[:, post_idx[t]] - Y_donors_pre_wt))
    for t in range(T1)
])
tau_sdid = att_yearly.mean()
p(f'Main SDiD ATT (avg 2007-2015): {tau_sdid:+.4f} pp')

# Cross-check against the saved CSV from run_sdid_robust_restricted.py
ref_csv  = pd.read_csv(os.path.join(TABLE_DIR, 'sdid_robust_restricted_yearly_att.csv'))
ref_atts = ref_csv['ATT_yearly'].values
max_diff = np.max(np.abs(att_yearly - ref_atts))
p(f'Max deviation from saved yearly ATTs: {max_diff:.6f} pp  '
  f'({"OK" if max_diff < 1e-3 else "MISMATCH -- check optimizer convergence"})')

# ── Jackknife: per-year ATT matrix (N x T1) ──────────────────────────────────
# Each fold j drops donor j, re-optimizes omega_j and lambda_j, then records
# the per-year ATT for all T1 years. Same zeta_j = (Nj * T1)^0.25 * sigma_j
# as run_sdid_robust_restricted.py.
p(f'\nRunning jackknife ({N} leave-one-out folds)...')
tau_jack_yr = np.empty((N, T1))

for j in range(N):
    keep      = [i for i in range(N) if i != j]
    Y_don_j   = Y_donors[keep, :]
    Yp_pre_j  = Y_don_j[:, pre_idx]
    Yp_post_j = Y_don_j[:, post_idx]
    Nj        = len(keep)

    sigma_j = np.std(np.diff(Yp_pre_j, axis=1), ddof=1)
    zeta_j  = (Nj * T1) ** 0.25 * sigma_j

    def obj_j(w):
        r = Yp_pre_j.T @ w - Y_pre_treat
        return np.mean(r**2) + (zeta_j**2) * np.dot(w, w)

    def grad_j(w):
        r = Yp_pre_j.T @ w - Y_pre_treat
        return (2.0 / T0) * (Yp_pre_j @ r) + 2.0 * (zeta_j**2) * w

    res_j = minimize(obj_j, np.ones(Nj) / Nj, jac=grad_j, method='SLSQP',
                     bounds=[(0, 1)] * Nj,
                     constraints={'type': 'eq', 'fun': lambda w: w.sum() - 1,
                                  'jac': lambda w: np.ones(Nj)},
                     options={'ftol': 1e-10, 'maxiter': 1000})
    omega_j = res_j.x

    Yp_post_mean_j = Yp_post_j.mean(axis=1)

    def lobj_j(lam_):
        r = Yp_pre_j @ lam_ - Yp_post_mean_j
        return np.mean(r**2)

    def lgrad_j(lam_):
        r = Yp_pre_j @ lam_ - Yp_post_mean_j
        return (2.0 / Nj) * (Yp_pre_j.T @ r)

    res_lj = minimize(lobj_j, np.ones(T0) / T0, jac=lgrad_j, method='SLSQP',
                      bounds=[(0, 1)] * T0,
                      constraints={'type': 'eq', 'fun': lambda l: l.sum() - 1,
                                   'jac': lambda l: np.ones(T0)},
                      options={'ftol': 1e-10, 'maxiter': 1000})
    lam_j = res_lj.x

    Y_treat_pre_wt_j  = float(Y_pre_treat @ lam_j)
    Y_donors_pre_wt_j = Yp_pre_j @ lam_j    # (Nj,)

    for t in range(T1):
        tau_jack_yr[j, t] = (
            (Y_treat[post_idx[t]] - Y_treat_pre_wt_j)
            - float(omega_j @ (Yp_post_j[:, t] - Y_donors_pre_wt_j))
        )

    if (j + 1) % 10 == 0:
        p(f'  {j+1}/{N} jackknife folds done...')

p('Jackknife complete.')

# Per-year SE: jackknife finite-sample formula
# SE_t = sqrt((N-1)/N * sum_j (tau_jack_yr[j,t] - mean_j tau_jack_yr[:,t])^2)
tau_jack_yr_mean = tau_jack_yr.mean(axis=0)  # (T1,) — cross-fold mean per year
se_jack_yr = np.sqrt(
    ((N - 1) / N)
    * np.sum((tau_jack_yr - tau_jack_yr_mean[np.newaxis, :]) ** 2, axis=0)
)
ci_lo_jack = att_yearly - 1.96 * se_jack_yr
ci_hi_jack = att_yearly + 1.96 * se_jack_yr

# ── Placebo: per-year ATT matrix (N x T1) ────────────────────────────────────
# Each donor p is treated as pseudo-CA. omega_p is re-estimated; lam is reused
# (same convention as run_sdid_robust_restricted.py lines 278-280).
# Per-year: tau_plac_yr[p, t] = (p_post[t] - p_pre @ lam) - w_p @ (p_Dpost[:,t] - p_Dpre @ lam)
p(f'\nRunning in-space placebo ({N} folds)...')
tau_plac_yr = np.empty((N, T1))

for pi, p_state in enumerate(donor_states):
    p_donors   = [s for s in donor_states if s != p_state]
    Np         = len(p_donors)
    sidx       = donor_states.index(p_state)
    p_Y_treat  = Y_donors[sidx, :]
    p_Y_donors = np.array([Y_donors[donor_states.index(s), :] for s in p_donors])

    p_pre   = p_Y_treat[pre_idx]
    p_post  = p_Y_treat[post_idx]
    p_Dpre  = p_Y_donors[:, pre_idx]
    p_Dpost = p_Y_donors[:, post_idx]

    # Placebo zeta: max(T1/T0, 1) * sigma_p  [same as run_sdid_robust_restricted.py]
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

    # Pre-period weighted baselines using main-sample lam
    p_pre_wt     = float(p_pre @ lam)
    p_don_pre_wt = p_Dpre @ lam    # (Np,)

    for t in range(T1):
        tau_plac_yr[pi, t] = (
            (p_post[t] - p_pre_wt)
            - float(w_p @ (p_Dpost[:, t] - p_don_pre_wt))
        )

    if (pi + 1) % 10 == 0:
        p(f'  {pi+1}/{N} placebo folds done...')

p('Placebo complete.')

# Per-year null distribution bounds: 2.5th and 97.5th percentiles across 43 states
ci_lo_perm = np.percentile(tau_plac_yr, 2.5,  axis=0)
ci_hi_perm = np.percentile(tau_plac_yr, 97.5, axis=0)

# ── Save CSVs ─────────────────────────────────────────────────────────────────
post_years_arr = np.array(post_years, dtype=int)

jack_df = pd.DataFrame({
    'year':    post_years_arr,
    'att':     att_yearly,
    'se':      se_jack_yr,
    'ci_low':  ci_lo_jack,
    'ci_high': ci_hi_jack,
})
jack_df.to_csv(os.path.join(TABLE_DIR, 'sdid_jackknife_per_year_ci.csv'), index=False)
p('\nSaved: sdid_jackknife_per_year_ci.csv  ->  outputs/tables/')

perm_df = pd.DataFrame({
    'year':    post_years_arr,
    'att':     att_yearly,
    'ci_low':  ci_lo_perm,
    'ci_high': ci_hi_perm,
})
perm_df.to_csv(os.path.join(TABLE_DIR, 'sdid_permutation_per_year_ci.csv'), index=False)
p('Saved: sdid_permutation_per_year_ci.csv  ->  outputs/tables/')

# ── Print per-year tables ─────────────────────────────────────────────────────
p('\nPer-year ATTs -- jackknife 95% CIs:')
p(f"  {'Year':>4}  {'ATT':>6}  {'SE':>6}  {'CI Low':>7}  {'CI High':>8}  {'CI excl 0?':>10}")
jack_excl = 0
for t in range(T1):
    excl = (ci_lo_jack[t] > 0) or (ci_hi_jack[t] < 0)
    jack_excl += int(excl)
    p(f"  {post_years_arr[t]:>4d}  {att_yearly[t]:>+6.3f}  {se_jack_yr[t]:>6.4f}  "
      f"{ci_lo_jack[t]:>+7.3f}  {ci_hi_jack[t]:>+8.3f}  "
      f"{'YES' if excl else 'no':>10}")
p(f'\n  {jack_excl} of {T1} years: CI excludes zero (jackknife)')

p('\nPer-year ATTs -- permutation null bounds [2.5th, 97.5th pctile across 43 donors]:')
p(f"  {'Year':>4}  {'ATT':>6}  {'Null lo':>8}  {'Null hi':>8}  {'ATT outside null?':>17}")
perm_excl = 0
for t in range(T1):
    excl = (att_yearly[t] > ci_hi_perm[t]) or (att_yearly[t] < ci_lo_perm[t])
    perm_excl += int(excl)
    p(f"  {post_years_arr[t]:>4d}  {att_yearly[t]:>+6.3f}  "
      f"{ci_lo_perm[t]:>+8.3f}  {ci_hi_perm[t]:>+8.3f}  "
      f"{'YES' if excl else 'no':>17}")
p(f'\n  {perm_excl} of {T1} years: CA ATT outside null distribution 95% bounds')

# ── Color palette and typography (matches make_presentation_figures.py) ───────
C_CA    = '#1B4F72'   # deep navy
C_POS   = '#2874A6'   # steel blue
C_REF   = '#17202A'   # near-black
C_GRID  = '#E8E8E8'   # light gray
C_NULL  = '#717D7E'   # medium gray -- null distribution bounds

plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.size':         13,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.linewidth':    0.8,
    'xtick.major.size':  4,
    'ytick.major.size':  4,
    'xtick.labelsize':   11,
    'ytick.labelsize':   11,
})
FT  = 16
FL  = 13
FA  = 10
DPI = 300

def save_fig(fig, name):
    for ext in ('png', 'pdf'):
        kw = dict(dpi=DPI) if ext == 'png' else {}
        fig.savefig(os.path.join(FIG_DIR, f'pres_{name}.{ext}'),
                    bbox_inches='tight', **kw)
    plt.close(fig)
    p(f'  Saved: pres_{name}.png + .pdf  ->  outputs/figures/')

tau_avg = att_yearly.mean()   # same as tau_sdid

# ── Figure 1: Jackknife per-year CIs ─────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 6))

ax.axhline(0,       color=C_REF, lw=1.5, zorder=2)
ax.axhline(tau_avg, color=C_CA,  lw=2.5, ls='--', zorder=4)

ax.bar(post_years_arr, att_yearly, color=C_POS,
       edgecolor=C_REF, lw=0.5, alpha=0.90, width=0.65, zorder=3)

# Per-year jackknife error bars (deep navy, capsize to distinguish from bars)
err_lo = att_yearly - ci_lo_jack
err_hi = ci_hi_jack - att_yearly
ax.errorbar(post_years_arr, att_yearly,
            yerr=[err_lo, err_hi],
            fmt='none', color=C_CA, lw=2.0, capsize=6, capthick=2.0, zorder=5)

# ATT value labels above the CI cap so they never overlap with error bars
for t in range(T1):
    y_label = ci_hi_jack[t] + 0.10
    ax.text(post_years_arr[t], y_label, f'{att_yearly[t]:+.2f}',
            ha='center', va='bottom',
            fontsize=10, fontweight='bold', color=C_REF, zorder=6)

ax.set_xticks(post_years_arr)
ax.set_xticklabels([str(y) for y in post_years_arr], fontsize=12)
ax.tick_params(axis='both', labelsize=12)
ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('ATT (percentage points)', fontsize=FL, labelpad=8)
ax.set_title(
    'SDiD Event Study — Year-by-Year Effects with Jackknife 95% CIs\n'
    'CA Paid Family Leave on Female LFP Rate, Women 25–54',
    fontsize=FT, fontweight='bold', pad=12
)

legend_handles = [
    mpatches.Patch(facecolor=C_POS, alpha=0.9, edgecolor=C_REF, lw=0.5,
                   label='ATT per year (2007–15)'),
    plt.Line2D([0], [0], color=C_CA, lw=2.5, ls='--',
               label=f'Avg ATT (2007–15) = +{tau_avg:.2f} pp'),
    plt.Line2D([0], [0], color=C_CA, lw=2.0,
               label='95% CI (leave-one-out jackknife)'),
]
ax.legend(handles=legend_handles, fontsize=11.5, loc='upper left',
          framealpha=0.9, edgecolor='#cccccc')

note = (f'Leave-one-out jackknife, N = {N} donor folds\n'
        f'{jack_excl} of {T1} years: CI excludes zero')
ax.text(0.97, 0.04, note,
        transform=ax.transAxes, fontsize=FA,
        ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                  edgecolor='#cccccc', alpha=0.90))

ax.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)
ax.set_axisbelow(True)

ymax_jack = max(3.5, float(np.max(ci_hi_jack)) + 1.0)
ax.set_ylim(-0.5, ymax_jack)

fig.tight_layout()
save_fig(fig, 'event_study_jackknife_ci')

# ── Figure 2: Permutation null-distribution bounds per year ───────────────────
fig, ax = plt.subplots(figsize=(13, 6))

ax.axhline(0,       color=C_REF, lw=1.5, zorder=2)
ax.axhline(tau_avg, color=C_CA,  lw=2.5, ls='--', zorder=4)

ax.bar(post_years_arr, att_yearly, color=C_POS,
       edgecolor=C_REF, lw=0.5, alpha=0.90, width=0.65, zorder=3)

# Per-year permutation bounds: gray whiskers showing the null distribution range
err_lo_p = att_yearly - ci_lo_perm
err_hi_p = ci_hi_perm - att_yearly
ax.errorbar(post_years_arr, att_yearly,
            yerr=[err_lo_p, err_hi_p],
            fmt='none', color=C_NULL, lw=2.0, capsize=6, capthick=2.0, zorder=5)

# ATT value labels above the null bound cap
for t in range(T1):
    y_label = max(ci_hi_perm[t], att_yearly[t]) + 0.10
    ax.text(post_years_arr[t], y_label, f'{att_yearly[t]:+.2f}',
            ha='center', va='bottom',
            fontsize=10, fontweight='bold', color=C_REF, zorder=6)

ax.set_xticks(post_years_arr)
ax.set_xticklabels([str(y) for y in post_years_arr], fontsize=12)
ax.tick_params(axis='both', labelsize=12)
ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('ATT (percentage points)', fontsize=FL, labelpad=8)
ax.set_title(
    'SDiD Event Study — Year-by-Year Effects with Permutation Null Bounds\n'
    'CA Paid Family Leave on Female LFP Rate, Women 25–54',
    fontsize=FT, fontweight='bold', pad=12
)

legend_handles_perm = [
    mpatches.Patch(facecolor=C_POS, alpha=0.9, edgecolor=C_REF, lw=0.5,
                   label='ATT per year (2007–15)'),
    plt.Line2D([0], [0], color=C_CA, lw=2.5, ls='--',
               label=f'Avg ATT (2007–15) = +{tau_avg:.2f} pp'),
    plt.Line2D([0], [0], color=C_NULL, lw=2.0,
               label='95% null bounds (in-space placebo)'),
]
ax.legend(handles=legend_handles_perm, fontsize=11.5, loc='upper left',
          framealpha=0.9, edgecolor='#cccccc')

note_p = (f'In-space placebo, N = {N} donor states\n'
          f'Whiskers: [2.5th, 97.5th pctile] of null ATTs\n'
          f'{perm_excl} of {T1} years: CA ATT outside null bounds')
ax.text(0.97, 0.04, note_p,
        transform=ax.transAxes, fontsize=FA,
        ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                  edgecolor='#cccccc', alpha=0.90))

ax.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)
ax.set_axisbelow(True)

ymax_perm = max(3.5, float(np.max(ci_hi_perm)) + 1.0,
                float(np.max(att_yearly)) + 1.5)
ax.set_ylim(-0.5, ymax_perm)

fig.tight_layout()
save_fig(fig, 'event_study_permutation_ci')

# ── Figure 3: Line-graph event study with jackknife CI band ──────────────────
# Reads the saved CSV rather than reusing in-memory arrays, so this block is
# self-contained and can also be run in isolation.
_jack_csv = pd.read_csv(os.path.join(TABLE_DIR, 'sdid_jackknife_per_year_ci.csv'))
_yr  = _jack_csv['year'].values.astype(int)
_att = _jack_csv['att'].values
_clo = _jack_csv['ci_low'].values
_chi = _jack_csv['ci_high'].values
_avg = _att.mean()

fig, ax = plt.subplots(figsize=(13, 6))

# Shaded CI band drawn first so the line plots on top
ax.fill_between(_yr, _clo, _chi, color=C_CA, alpha=0.15, zorder=1,
                label='95% CI (leave-one-out jackknife)')

# Horizontal reference lines
ax.axhline(0,    color='#566573', lw=1.0, zorder=2)
ax.axhline(_avg, color=C_CA, lw=1.5, ls='--', zorder=3,
           label=f'Avg ATT (2007–15) = +{_avg:.2f} pp')

# ATT line with circular markers
ax.plot(_yr, _att, color=C_CA, lw=2.5, marker='o', markersize=6,
        markerfacecolor=C_CA, markeredgecolor=C_CA, zorder=4,
        label='ATT per year (2007–15)')

# Point labels offset above each marker
for t in range(len(_yr)):
    offset = 0.13 + (_chi[t] - _att[t]) * 0.05   # nudge above CI cap
    ax.text(_yr[t], _att[t] + offset, f'{_att[t]:+.2f}',
            ha='center', va='bottom',
            fontsize=10, fontweight='bold', color=C_REF, zorder=5)

ax.set_xticks(_yr)
ax.set_xticklabels([str(y) for y in _yr], fontsize=11)
ax.tick_params(axis='both', labelsize=11)
ax.set_xlabel('Year', fontsize=13, labelpad=8)
ax.set_ylabel('ATT (percentage points)', fontsize=13, labelpad=8)
ax.set_title(
    'SDiD Event Study — Year-by-Year ATT with 95% Jackknife CIs\n'
    'CA Paid Family Leave on Female LFP Rate, Women 25–54',
    fontsize=14, fontweight='bold', pad=12
)
ax.set_ylim(-1.0, 3.5)

# Legend: reorder so ATT line is first
_handles, _labels = ax.get_legend_handles_labels()
# current order: fill_between, axhline, plot — reorder to plot/axhline/fill
_h = [_handles[2], _handles[1], _handles[0]]
_l = [_labels[2],  _labels[1],  _labels[0]]
ax.legend(_h, _l, fontsize=11, loc='upper left', framealpha=0.9,
          edgecolor='#cccccc')

# Annotation box lower right
_jack_excl = int(np.sum((_clo > 0) | (_chi < 0)))
ax.text(0.97, 0.04,
        f'Leave-one-out jackknife, N = {N} donor folds\n'
        f'{_jack_excl} of {T1} years: CI excludes zero',
        transform=ax.transAxes, fontsize=9,
        ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                  edgecolor='#cccccc', alpha=0.90))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.yaxis.grid(False)
ax.xaxis.grid(False)

fig.tight_layout()
save_fig(fig, 'event_study_line')

# ── Confirm pres_event_study_clean.png was not touched ───────────────────────
clean_path = os.path.join(FIG_DIR, 'pres_event_study_clean.png')
p(f'\npres_event_study_clean.png: '
  f'{"present (untouched)" if os.path.exists(clean_path) else "NOT FOUND"}')

elapsed = time.time() - t_start
p(f'Run time: {elapsed:.1f} seconds')
