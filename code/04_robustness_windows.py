"""
Window Sensitivity Battery: Synthetic DiD
California Paid Family Leave (2004), Female Labor Force Participation

Four post-window specifications on the 1995-2015 panel with the
43-state donor pool from run_sdid_robust_restricted.py.
Section 5.3 robustness table.

Donor exclusions: NY (TDI), RI (TDI), HI (TDI), NJ (PFL 2008),
                  WA (PFL 2007), AK (structural outlier), DC (city-state)

Zeta conventions preserved exactly from run_sdid_robust_restricted.py:
  Main estimation : zeta   = (N  * T1)^0.25 * sigma
  Jackknife fold  : zeta_j = (Nj * T1)^0.25 * sigma_j
  Placebo fold    : zeta_p = max(T1 / T0, 1.0) * sig_p
Each spec propagates its own T1 into all three formulas.

Outputs (new files only, no overwrites):
  outputs/tables/window_sensitivity_table.csv
  outputs/tables/window_sensitivity_table.tex
"""

import os
import time
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

DIR  = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(DIR, ".."))

def p(*args): print(*args, flush=True)
np.random.seed(42)
t_start = time.time()

# ── Parameters (identical to run_sdid_robust_restricted.py) ───────────────────
DATA_FILE = os.path.join(ROOT, 'data', 'derived', 'state_year_panel_deduped_1995_2015.csv')
TREATED    = 'California'
TREAT_YEAR = 2004
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

OUT_DIR  = os.path.join(ROOT, 'outputs', 'tables')
os.makedirs(OUT_DIR, exist_ok=True)

CSV_PATH = os.path.join(OUT_DIR, 'window_sensitivity_table.csv')
TEX_PATH = os.path.join(OUT_DIR, 'window_sensitivity_table.tex')
for _path in (CSV_PATH, TEX_PATH):
    if os.path.exists(_path):
        raise FileExistsError(
            f'Output already exists: {_path}  -- move or delete it first.')

# ── Load data ─────────────────────────────────────────────────────────────────
agg = pd.read_csv(DATA_FILE)
agg.columns = agg.columns.str.strip()

all_states   = sorted(agg['state'].unique())
donor_states = [s for s in all_states if s != TREATED and s not in EXCLUDED]
years        = sorted(agg['YEAR'].unique())   # [1995, ..., 2015], length 21
N            = len(donor_states)

years_arr = np.array(years)
pre_idx   = np.where(years_arr < TREAT_YEAR)[0]   # [0..8], 1995 to 2003
T0        = len(pre_idx)                           # 9

# ── Build full outcome matrix (all 21 years) ──────────────────────────────────
Y_mat    = agg.pivot(index='state', columns='YEAR', values=OUTCOME)
Y_mat    = Y_mat.reindex(index=[TREATED] + donor_states, columns=years)
Y_full   = Y_mat.values.astype(float)
Y_treat  = Y_full[0, :]
Y_donors = Y_full[1:, :]

Y_pre_treat  = Y_treat[pre_idx]
Y_pre_donors = Y_donors[:, pre_idx]

# ── Define the four window specifications ─────────────────────────────────────
# post_idx uses np.where with boolean conditions so non-contiguous windows
# (specs 2 and 3) are handled correctly with index arrays, not slices.
SPECS = [
    {
        'name':        'Full window',
        'window_desc': '2004 to 2015',
        'post_idx':    np.where(years_arr >= 2004)[0],
    },
    {
        'name':        'Drop 2005 only',
        'window_desc': '2004, 2006 to 2015',
        'post_idx':    np.where((years_arr >= 2004) & (years_arr != 2005))[0],
    },
    {
        'name':        'Drop 2006 only',
        'window_desc': '2004 to 2005, 2007 to 2015',
        'post_idx':    np.where((years_arr >= 2004) & (years_arr != 2006))[0],
    },
    {
        'name':        'Drop both (preferred)',
        'window_desc': '2007 to 2015',
        'post_idx':    np.where(years_arr >= 2007)[0],
    },
]

# ── Pre-run verification ───────────────────────────────────────────────────────
p('=' * 65)
p('  PRE-RUN VERIFICATION')
p('=' * 65)
p(f'  Data file   : {DATA_FILE}')
p(f'  Treated     : {TREATED}')
p(f'  Donor pool  : {N} states after exclusions')
p(f'  Excluded    : {", ".join(EXCLUDED)}')
p(f'  Pre-period  : {years[pre_idx[0]]} to {years[pre_idx[-1]]}  (T0={T0})')
p()
p('  Post-window index arrays:')
for spec in SPECS:
    idx = spec['post_idx']
    T1  = len(idx)
    yrs = [years[i] for i in idx]
    p(f'    {spec["name"]:30s}  T1={T1}')
    p(f'      post_idx = {list(idx)}')
    p(f'      years    = {yrs}')
p('=' * 65)


# ── SDiD estimation function ───────────────────────────────────────────────────
# Structure mirrors run_sdid_robust_restricted.py exactly.
# Zeta conventions are preserved and T1 is propagated from the caller.
def run_spec(spec_post_idx, T1,
             Y_treat, Y_donors, Y_pre_treat, Y_pre_donors,
             pre_idx, T0, N, donor_states):

    Y_post_treat  = Y_treat[spec_post_idx]
    Y_post_donors = Y_donors[:, spec_post_idx]

    # ── omega (unit weights) ──────────────────────────────────────────────────
    # zeta = (N * T1)^0.25 * sigma
    sigma = np.std(np.diff(Y_pre_donors, axis=1), ddof=1)
    zeta  = (N * T1) ** 0.25 * sigma

    def omega_obj(w):
        resid = Y_pre_donors.T @ w - Y_pre_treat
        return np.mean(resid**2) + (zeta**2) * np.dot(w, w)

    def omega_grad(w):
        resid = Y_pre_donors.T @ w - Y_pre_treat
        return (2.0 / T0) * (Y_pre_donors @ resid) + 2.0 * (zeta**2) * w

    res_omega = minimize(
        omega_obj, np.ones(N) / N, jac=omega_grad, method='SLSQP',
        bounds=[(0, 1)] * N,
        constraints={'type': 'eq', 'fun': lambda w: w.sum() - 1,
                     'jac': lambda w: np.ones(N)},
        options={'ftol': 1e-12, 'maxiter': 3000}
    )
    omega = res_omega.x

    # ── lambda (time weights) ─────────────────────────────────────────────────
    # Target is the post-period donor mean over spec_post_idx columns.
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
        constraints={'type': 'eq', 'fun': lambda l: l.sum() - 1,
                     'jac': lambda l: np.ones(T0)},
        options={'ftol': 1e-12, 'maxiter': 3000}
    )
    lam = res_lambda.x

    # ── SDiD point estimate ───────────────────────────────────────────────────
    Y_treat_post_avg  = Y_post_treat.mean()
    Y_treat_pre_wt    = float(Y_pre_treat @ lam)
    Y_donors_post_avg = Y_post_donors.mean(axis=1)
    Y_donors_pre_wt   = Y_pre_donors @ lam

    tau_sdid = (Y_treat_post_avg - Y_treat_pre_wt) \
             - float(omega @ (Y_donors_post_avg - Y_donors_pre_wt))

    # ── Leave-one-unit-out jackknife SE ───────────────────────────────────────
    # zeta_j = (Nj * T1)^0.25 * sigma_j  -- T1 propagated from this spec
    # Index arrays (pre_idx, spec_post_idx) used throughout, no slices.
    tau_jack = np.empty(N)
    for j in range(N):
        keep      = [i for i in range(N) if i != j]
        Y_don_j   = Y_donors[keep, :]
        Yp_pre_j  = Y_don_j[:, pre_idx]
        Yp_post_j = Y_don_j[:, spec_post_idx]

        Nj      = len(keep)
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

        Y_treat_pre_wt_j = float(Y_pre_treat @ lam_j)
        tau_jack[j] = (Y_treat_post_avg - Y_treat_pre_wt_j) \
                    - float(omega_j @ (Yp_post_j.mean(axis=1) - Yp_pre_j @ lam_j))

    se_jack = np.sqrt(((N - 1) / N) * np.sum((tau_jack - tau_sdid)**2))
    t_stat  = tau_sdid / se_jack if se_jack > 0 else np.nan
    pval    = 2 * stats.t.sf(abs(t_stat), df=N - 1) if not np.isnan(t_stat) else np.nan
    ci_lo   = tau_sdid - 1.96 * se_jack
    ci_hi   = tau_sdid + 1.96 * se_jack

    # ── Pre-trend RMSE ────────────────────────────────────────────────────────
    synth_pre = Y_pre_donors.T @ omega
    pre_rmse  = np.sqrt(np.mean((Y_pre_treat - synth_pre)**2))

    # ── In-space placebo (same post window as main spec) ──────────────────────
    # Uses the main-sample lam (this spec's time weights) for placebo tau,
    # matching the convention in run_sdid_robust_restricted.py line 279.
    # zeta_p = max(T1 / T0, 1.0) * sig_p  -- T1 propagated from this spec
    placebo_taus = []
    for plac_i, p_state in enumerate(donor_states):
        p_donors   = [s for s in donor_states if s != p_state]
        Np         = len(p_donors)
        p_Y_treat  = Y_donors[donor_states.index(p_state), :]
        p_Y_donors = np.array([Y_donors[donor_states.index(s), :] for s in p_donors])

        p_pre   = p_Y_treat[pre_idx]
        p_post  = p_Y_treat[spec_post_idx]
        p_Dpre  = p_Y_donors[:, pre_idx]
        p_Dpost = p_Y_donors[:, spec_post_idx]

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

        if (plac_i + 1) % 10 == 0:
            p(f'    placebo {plac_i + 1}/{N}...')

    placebo_taus = np.array(placebo_taus)
    n_extreme    = int(np.sum(np.abs(placebo_taus) >= abs(tau_sdid)))
    pval_placebo = n_extreme / len(placebo_taus)

    return {
        'tau':          tau_sdid,
        'se':           se_jack,
        'pval':         pval,
        'ci_lo':        ci_lo,
        'ci_hi':        ci_hi,
        'pval_placebo': pval_placebo,
        'n_extreme':    n_extreme,
        'pre_rmse':     pre_rmse,
    }


# ── Run all four specifications ────────────────────────────────────────────────
results_rows = []

for spec_num, spec in enumerate(SPECS, start=1):
    post_idx_s = spec['post_idx']
    T1_s       = len(post_idx_s)
    p()
    p(f'Running spec {spec_num} of {len(SPECS)}: {spec["name"]}')
    p(f'  Post window : {spec["window_desc"]}  (T1={T1_s})')
    p(f'  Jackknife   : {N} folds  |  Placebo: {N} folds')

    t0_s = time.time()
    res  = run_spec(
        post_idx_s, T1_s,
        Y_treat, Y_donors, Y_pre_treat, Y_pre_donors,
        pre_idx, T0, N, donor_states
    )
    elapsed_s = time.time() - t0_s

    p(f'  ATT         = {res["tau"]:+.4f} pp')
    p(f'  SE          = {res["se"]:.4f}   p (jack.) = {res["pval"]:.4f}')
    p(f'  95% CI      = [{res["ci_lo"]:+.4f}, {res["ci_hi"]:+.4f}]')
    p(f'  p (perm.)   = {res["pval_placebo"]:.4f}  ({res["n_extreme"]} of {N} donors)')
    p(f'  Pre-RMSE    = {res["pre_rmse"]:.4f} pp')
    p(f'  Done in {elapsed_s:.1f}s')

    results_rows.append({
        'spec_name':        spec['name'],
        'post_window':      spec['window_desc'],
        'T1':               T1_s,
        'ATT_pp':           round(res['tau'],          4),
        'jackknife_SE':     round(res['se'],            4),
        'jackknife_pval':   round(res['pval'],          4),
        'CI_95_lo':         round(res['ci_lo'],         4),
        'CI_95_hi':         round(res['ci_hi'],         4),
        'permutation_pval': round(res['pval_placebo'],  4),
        'pre_trend_RMSE':   round(res['pre_rmse'],      4),
    })


# ── ATT gradient diagnostic ────────────────────────────────────────────────────
p()
p('=' * 65)
p('  ATT GRADIENT  (full window to preferred window)')
p('=' * 65)
for row in results_rows:
    p(f'  {row["spec_name"]:32s}  ATT = {row["ATT_pp"]:+.4f} pp'
      f'   p = {row["jackknife_pval"]:.4f}')
p('=' * 65)


# ── Spec 4 replication check against run_sdid_robust_restricted.py ────────────
# Reference ATT: 1.1913 is the restricted-window SDiD ATT baked into
# make_presentation_figures.py (tau_ca = 1.1913, line 453).
# Note: make_presentation_figures.py also prints tau_restr = +1.1845, but that
# is the mean of event-study yearly ATTs from the full-window run (a different
# estimator), not the re-optimized restricted-window SDiD ATT.
REF_ATT     = 1.1913
spec4_att   = results_rows[3]['ATT_pp']
discrepancy = abs(spec4_att - REF_ATT)
p()
if discrepancy < 1e-3:
    p(f'  Spec 4 replication check: {spec4_att:+.4f} vs reference {REF_ATT:+.4f}   OK')
else:
    p(f'  Spec 4 replication check: {spec4_att:+.4f} vs reference {REF_ATT:+.4f}   ** MISMATCH **')
    p(f'  Discrepancy = {discrepancy:.6f} pp  (threshold 0.001 pp)')


# ── Save CSV ───────────────────────────────────────────────────────────────────
df = pd.DataFrame(results_rows)
df.to_csv(CSV_PATH, index=False)
p(f'\nSaved: {CSV_PATH}')


# ── Save LaTeX tabular fragment ────────────────────────────────────────────────
# Produces a bare tabular environment suitable for \input{} in Overleaf.
# No \begin{table}, \caption, or \label wrapper.
def fmt_ci(lo, hi):
    return f'[{lo:+.2f}, {hi:+.2f}]'

tex_lines = []
tex_lines.append(r'\begin{tabular}{llccccccc}')
tex_lines.append(r'\hline\hline')
tex_lines.append(
    r'Specification & Post Window & $T_1$ & ATT (pp) & SE & '
    r'$p$ (jack.) & 95\% CI & $p$ (perm.) & Pre-RMSE \\'
)
tex_lines.append(r'\hline')

for row in results_rows:
    ci   = fmt_ci(row['CI_95_lo'], row['CI_95_hi'])
    line = (
        f'{row["spec_name"]} & '
        f'{row["post_window"]} & '
        f'{row["T1"]} & '
        f'${row["ATT_pp"]:+.2f}$ & '
        f'${row["jackknife_SE"]:.2f}$ & '
        f'${row["jackknife_pval"]:.3f}$ & '
        f'${ci}$ & '
        f'${row["permutation_pval"]:.3f}$ & '
        f'${row["pre_trend_RMSE"]:.3f}$ \\\\'
    )
    tex_lines.append(line)

tex_lines.append(r'\hline\hline')
tex_lines.append(r'\end{tabular}')

with open(TEX_PATH, 'w') as fh:
    fh.write('\n'.join(tex_lines) + '\n')

p(f'Saved: {TEX_PATH}')

p(f'\nTotal runtime: {time.time() - t_start:.1f}s')
