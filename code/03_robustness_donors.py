"""
Donor-Pool Robustness Checks — Synthetic DiD
CA Paid Family Leave (2004) → Female LFP, Restricted Post Window 2007-2015

Baseline: 43-state pool, excl. TDI/own-PFL/DC/AK.
Headline: ATT = +1.19 pp, SE = 0.49, jackknife p = 0.019.

Four robustness checks:
  Check 1 — Exclude top Appalachian donors (WV, LA, KY) → 40 states
  Check 2 — Western census region only (AZ, CO, ID, MT, NV, NM, OR, UT, WY) → 9 states
  Check 3 — Demographically matched pool (AZ, CO, FL, IL, NV, NM, OR, TX) → 8 states
  Check 4 — Leave-one-out across top-5 donors (WV, LA, KY, AZ, OK) → 5 × 42-state specs

In-space placebo is NOT run for alternative specs (small pools → unstable distributions).

Outputs:
  sdid_donor_robustness.csv
  pres_robustness_table.png + pres_robustness_table.pdf
"""

import os, time, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import minimize
from scipy import stats
warnings.filterwarnings('ignore')

def p(*a): print(*a, flush=True)
t_start = time.time()
np.random.seed(42)

# ── Configuration ─────────────────────────────────────────────────────────────
DIR       = os.path.dirname(os.path.abspath(__file__))
ROOT      = os.path.abspath(os.path.join(DIR, ".."))
DATA_FILE = os.path.join(ROOT, 'data/derived/state_year_panel_deduped_1995_2015.csv')
TREATED   = 'California'
TREAT_YEAR = 2004
POST_START = 2007
OUTCOME   = 'flfp_pct'

# Expected baseline numbers for sanity check
BASELINE_ATT  = 1.19
BASELINE_SE   = 0.49
BASELINE_PVAL = 0.019
TOL_ATT, TOL_SE, TOL_PVAL = 0.02, 0.02, 0.005

# Presentation colors — match make_presentation_figures.py exactly
C_CA  = '#1B4F72'
C_REF = '#17202A'
FT, FL = 16, 13

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 11,
    'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.labelsize': 11, 'ytick.labelsize': 11,
})

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_FILE)
df.columns = df.columns.str.strip()

BASELINE_EXCL = [
    'New York', 'Rhode Island', 'Hawaii',
    'New Jersey', 'Washington',
    'Alaska', 'District of Columbia',
]

all_states      = sorted(df['state'].unique())
baseline_donors = [s for s in all_states if s != TREATED and s not in BASELINE_EXCL]
years           = sorted(df['YEAR'].unique())
years_arr       = np.array(years)
pre_idx         = np.where(years_arr < TREAT_YEAR)[0]   # 1995-2003
post_idx        = np.where(years_arr >= POST_START)[0]  # 2007-2015
T0, T1          = len(pre_idx), len(post_idx)

Y_mat = df.pivot(index='state', columns='YEAR', values=OUTCOME)
Y_mat = Y_mat.reindex(index=[TREATED] + baseline_donors, columns=years)
Y_full = Y_mat.values.astype(float)

Y_treat_full  = Y_full[0, :]       # CA × 21 years
Y_donors_full = Y_full[1:, :]      # 43 baseline donors × 21 years

p(f'Panel: {len(all_states)} states, {len(years)} years ({years[0]}-{years[-1]})')
p(f'Baseline donor pool: {len(baseline_donors)} states')
p(f'Pre-period:  {years_arr[pre_idx[0]]}-{years_arr[pre_idx[-1]]}  (T0={T0})')
p(f'Post-period: {years_arr[post_idx[0]]}-{years_arr[post_idx[-1]]}  (T1={T1})')

# ── Closure-safe helper factories ─────────────────────────────────────────────
# Using factory functions avoids Python closure-in-loop bugs.
# Each factory captures current argument values by copying them into the closure.

def _omega_fns(Yp, Yt, z, t0):
    def obj(w):  return np.mean((Yp.T @ w - Yt)**2) + (z**2) * np.dot(w, w)
    def grad(w): r = Yp.T @ w - Yt; return (2/t0)*(Yp @ r) + 2*(z**2)*w
    return obj, grad

def _lam_fns(Yp, tgt, n):
    def obj(lam):  return np.mean((Yp @ lam - tgt)**2)
    def grad(lam): r = Yp @ lam - tgt; return (2/n)*(Yp.T @ r)
    return obj, grad

# ── Core SDiD estimation function ─────────────────────────────────────────────
def run_spec(donor_list, label=''):
    """
    Run SDiD with restricted post window 2007-2015 for a given donor pool.

    Reuses the pre-loaded Y_treat_full and Y_donors_full (baseline 43-state
    matrix). Extracts the relevant rows for the requested donor_list.

    Returns a dict with: N, ATT, SE, CI_lo, CI_hi, pval, RMSPE, omega (Series),
    converged (bool), df (jackknife degrees of freedom).
    """
    N     = len(donor_list)
    d_idx = [baseline_donors.index(s) for s in donor_list]

    Y_don    = Y_donors_full[d_idx, :]
    Y_pre_t  = Y_treat_full[pre_idx]          # (T0,) — CA pre-period
    Y_post_t = Y_treat_full[post_idx]         # (T1,) — CA post-period
    Y_pre_d  = Y_don[:, pre_idx]              # (N × T0)
    Y_post_d = Y_don[:, post_idx]             # (N × T1)

    # Regularization — same formula as run_sdid_robust_restricted.py
    sigma = np.std(np.diff(Y_pre_d, axis=1), ddof=1)
    zeta  = (N * T1) ** 0.25 * sigma

    # ── Omega (unit weights) ──────────────────────────────────────────────────
    o_fn, o_gn = _omega_fns(Y_pre_d, Y_pre_t, zeta, T0)
    res_w = minimize(o_fn, np.ones(N)/N, jac=o_gn, method='SLSQP',
                     bounds=[(0, 1)]*N,
                     constraints={'type': 'eq', 'fun': lambda w: w.sum()-1},
                     options={'ftol': 1e-12, 'maxiter': 3000})
    omega = res_w.x

    # ── Lambda (time weights) ─────────────────────────────────────────────────
    target = Y_post_d.mean(axis=1)           # mean over 2007-2015
    l_fn, l_gn = _lam_fns(Y_pre_d, target, N)
    res_l = minimize(l_fn, np.ones(T0)/T0, jac=l_gn, method='SLSQP',
                     bounds=[(0, 1)]*T0,
                     constraints={'type': 'eq', 'fun': lambda l: l.sum()-1},
                     options={'ftol': 1e-12, 'maxiter': 3000})
    lam = res_l.x

    # ── Point estimate ────────────────────────────────────────────────────────
    tau = (Y_post_t.mean() - float(Y_pre_t @ lam)) \
        - float(omega @ (Y_post_d.mean(axis=1) - Y_pre_d @ lam))

    # ── Pre-treatment RMSPE ───────────────────────────────────────────────────
    rmspe = np.sqrt(np.mean((Y_pre_t - Y_pre_d.T @ omega)**2))

    # ── Jackknife SE (leave-one-unit-out) ────────────────────────────────────
    # Each fold drops one donor, re-estimates omega and lambda, recomputes ATT.
    # The same restricted post window (post_idx) is used in every fold.
    tau_j = np.empty(N)
    for j in range(N):
        keep   = [i for i in range(N) if i != j]
        Yp     = Y_pre_d[keep, :]         # (N-1 × T0)
        Yk     = Y_post_d[keep, :]        # (N-1 × T1)
        Nj     = len(keep)
        sig_j  = np.std(np.diff(Yp, axis=1), ddof=1)
        zeta_j = (Nj * T1) ** 0.25 * sig_j
        tgt_j  = Yk.mean(axis=1)

        oo, og = _omega_fns(Yp, Y_pre_t, zeta_j, T0)
        rj = minimize(oo, np.ones(Nj)/Nj, jac=og, method='SLSQP',
                      bounds=[(0, 1)]*Nj,
                      constraints={'type': 'eq', 'fun': lambda w: w.sum()-1},
                      options={'ftol': 1e-10, 'maxiter': 1000})
        om_j = rj.x

        lo, lg = _lam_fns(Yp, tgt_j, Nj)
        rl = minimize(lo, np.ones(T0)/T0, jac=lg, method='SLSQP',
                      bounds=[(0, 1)]*T0,
                      constraints={'type': 'eq', 'fun': lambda l: l.sum()-1},
                      options={'ftol': 1e-10, 'maxiter': 1000})
        lm_j = rl.x

        tau_j[j] = (Y_post_t.mean() - float(Y_pre_t @ lm_j)) \
                 - float(om_j @ (Yk.mean(axis=1) - Yp @ lm_j))

    se    = np.sqrt(((N-1)/N) * np.sum((tau_j - tau)**2))
    tstat = tau / se if se > 0 else np.nan
    pval  = 2 * stats.t.sf(abs(tstat), df=N-1) if not np.isnan(tstat) else np.nan

    return {
        'N': N, 'ATT': tau, 'SE': se,
        'CI_lo': tau - 1.96*se, 'CI_hi': tau + 1.96*se,
        'pval': pval, 'RMSPE': rmspe,
        'omega': pd.Series(omega, index=donor_list).sort_values(ascending=False),
        'converged': res_w.success and res_l.success,
        'df': N-1,
    }

# ── STEP 0: Baseline sanity check ────────────────────────────────────────────
p('\n' + '='*65)
p('  STEP 0 — BASELINE SANITY CHECK')
p('='*65)
p(f'  Running baseline spec ({len(baseline_donors)} donors)...')

base = run_spec(baseline_donors, 'baseline')
p(f'  ATT   : {base["ATT"]:+.4f} pp  (expected ~+{BASELINE_ATT})')
p(f'  SE    : {base["SE"]:.4f}   (expected ~{BASELINE_SE})')
p(f'  p-val : {base["pval"]:.4f}   (expected ~{BASELINE_PVAL})')

drift = []
if abs(base['ATT'] - BASELINE_ATT) > TOL_ATT:
    drift.append(f'ATT: got {base["ATT"]:+.4f}, expected ~+{BASELINE_ATT} (tol {TOL_ATT})')
if abs(base['SE'] - BASELINE_SE) > TOL_SE:
    drift.append(f'SE: got {base["SE"]:.4f}, expected ~{BASELINE_SE} (tol {TOL_SE})')
if abs(base['pval'] - BASELINE_PVAL) > TOL_PVAL:
    drift.append(f'p-val: got {base["pval"]:.4f}, expected ~{BASELINE_PVAL} (tol {TOL_PVAL})')

if drift:
    p('\n  *** BASELINE DRIFT DETECTED — STOPPING ***')
    for d in drift: p(f'    {d}')
    p('  Reconcile run_sdid_robust_restricted.py before proceeding.')
    raise SystemExit(1)

p('  Baseline reproduced within tolerance. Proceeding with robustness checks.')

# ── CHECK 1: Exclude top Appalachian donors ───────────────────────────────────
p('\n' + '='*65)
p('  CHECK 1 — Exclude top Appalachian donors (WV, LA, KY) → 40 states')
p('  Rationale: these three carry the heaviest omega weights and represent')
p('  a "low FLFP from deindustrialization" mechanism different from CA.')
p('='*65)

APPALACHIAN = ['West Virginia', 'Louisiana', 'Kentucky']
pool1 = [s for s in baseline_donors if s not in APPALACHIAN]
p(f'  Donor pool: {len(pool1)} states')
r1 = run_spec(pool1)
p(f'  ATT={r1["ATT"]:+.4f}  SE={r1["SE"]:.4f}  p={r1["pval"]:.4f}  '
  f'CI=[{r1["CI_lo"]:+.4f}, {r1["CI_hi"]:+.4f}]  RMSPE={r1["RMSPE"]:.4f}')
p(f'  Converged: {r1["converged"]}')
p('  Top-5 omega (new weights after removing WV/LA/KY):')
p(r1['omega'].head(5).to_string(float_format='{:.4f}'.format))

# ── CHECK 2: Western census region only ──────────────────────────────────────
p('\n' + '='*65)
p('  CHECK 2 — Western census region only → 9 states')
p('  Rationale: CA economy is structurally western; western peers may be')
p('  better counterfactuals than Appalachian/Southern states.')
p('  WA excluded (own PFL 2007). HI, AK excluded per baseline.')
p('='*65)

WESTERN = ['Arizona', 'Colorado', 'Idaho', 'Montana', 'Nevada',
           'New Mexico', 'Oregon', 'Utah', 'Wyoming']
missing_w = [s for s in WESTERN if s not in baseline_donors]
if missing_w:
    p(f'  WARNING: {missing_w} not in baseline pool — removing')
    WESTERN = [s for s in WESTERN if s in baseline_donors]
pool2 = WESTERN
p(f'  Donor pool: {len(pool2)} states: {pool2}')
p(f'  Note: N={len(pool2)} → df={len(pool2)-1} for t-test; CI uses ±1.96 SE for table comparability.')

r2 = run_spec(pool2)
p(f'  ATT={r2["ATT"]:+.4f}  SE={r2["SE"]:.4f}  p={r2["pval"]:.4f}  '
  f'CI=[{r2["CI_lo"]:+.4f}, {r2["CI_hi"]:+.4f}]  RMSPE={r2["RMSPE"]:.4f}')
p(f'  Converged: {r2["converged"]}')
p('  Top-5 omega:')
p(r2['omega'].head(5).to_string(float_format='{:.4f}'.format))
p('  In-space placebo NOT run (small pool → unstable permutation distribution).')

# ── CHECK 3: Demographically matched pool ────────────────────────────────────
p('\n' + '='*65)
p('  CHECK 3 — Demographically matched pool → 8 states')
p('  Selection: 2010 Census Hispanic share > 15% or foreign-born share > 10%.')
p('  States: AZ, CO, FL, IL, NV, NM, OR, TX.')
p('='*65)

DEMOGRAPHIC = ['Arizona', 'Colorado', 'Florida', 'Illinois',
               'Nevada', 'New Mexico', 'Oregon', 'Texas']
missing_d = [s for s in DEMOGRAPHIC if s not in baseline_donors]
if missing_d:
    p(f'  WARNING: {missing_d} not in baseline pool — removing')
    DEMOGRAPHIC = [s for s in DEMOGRAPHIC if s in baseline_donors]
pool3 = DEMOGRAPHIC
p(f'  Donor pool: {len(pool3)} states: {pool3}')
p(f'  Note: N={len(pool3)} → df={len(pool3)-1}; CI uses ±1.96 SE for table comparability.')

r3 = run_spec(pool3)
p(f'  ATT={r3["ATT"]:+.4f}  SE={r3["SE"]:.4f}  p={r3["pval"]:.4f}  '
  f'CI=[{r3["CI_lo"]:+.4f}, {r3["CI_hi"]:+.4f}]  RMSPE={r3["RMSPE"]:.4f}')
p(f'  Converged: {r3["converged"]}')
p('  Top-5 omega:')
p(r3['omega'].head(min(5, len(pool3)), ).to_string(float_format='{:.4f}'.format))
p('  In-space placebo NOT run (small pool → unstable permutation distribution).')

# ── CHECK 4: Leave-one-out across top-5 donors ───────────────────────────────
p('\n' + '='*65)
p('  CHECK 4 — Leave-one-out across top-5 donors')
p('  Top-5 from baseline: West Virginia, Louisiana, Kentucky, Arizona, Oklahoma')
p('  Each run uses 42 donors (43 baseline minus one state).')
p('='*65)

TOP5 = ['West Virginia', 'Louisiana', 'Kentucky', 'Arizona', 'Oklahoma']
loo_results = {}
for state in TOP5:
    pool_loo = [s for s in baseline_donors if s != state]
    p(f'\n  LOO: drop {state} → {len(pool_loo)} donors')
    r = run_spec(pool_loo, f'LOO_{state}')
    loo_results[state] = r
    p(f'    ATT={r["ATT"]:+.4f}  SE={r["SE"]:.4f}  p={r["pval"]:.4f}  '
      f'CI=[{r["CI_lo"]:+.4f}, {r["CI_hi"]:+.4f}]  RMSPE={r["RMSPE"]:.4f}')
    p(f'    Converged: {r["converged"]}')
    p(f'    New top-3 omega: '
      + ', '.join(f'{s} ({w:.4f})' for s, w in r['omega'].head(3).items()))

# ── Compile results ───────────────────────────────────────────────────────────
p('\n' + '='*65)
p('  COMPILING RESULTS TABLE')
p('='*65)

specs_ordered = [
    ('Baseline (43 states)',            base),
    ('No Appalachian: drop WV, LA, KY', r1),
    ('Western region only (9 states)',  r2),
    ('Demographic match (8 states)',    r3),
]
loo_labels = {
    'West Virginia': 'LOO: drop West Virginia',
    'Louisiana':     'LOO: drop Louisiana',
    'Kentucky':      'LOO: drop Kentucky',
    'Arizona':       'LOO: drop Arizona',
    'Oklahoma':      'LOO: drop Oklahoma',
}
for state in TOP5:
    if state in loo_results:
        specs_ordered.append((loo_labels[state], loo_results[state]))

csv_rows = []
for label, res in specs_ordered:
    csv_rows.append({
        'Specification': label,
        'N_donors':      res['N'],
        'ATT_pp':        round(res['ATT'], 4),
        'SE':            round(res['SE'],  4),
        'CI_lo':         round(res['CI_lo'], 4),
        'CI_hi':         round(res['CI_hi'], 4),
        'Jackknife_p':   round(res['pval'],  4),
        'RMSPE':         round(res['RMSPE'], 4),
    })
    p(f'  {label}:  ATT={res["ATT"]:+.3f}  SE={res["SE"]:.3f}  p={res["pval"]:.3f}  '
      f'RMSPE={res["RMSPE"]:.3f}')

results_df = pd.DataFrame(csv_rows)
csv_path = os.path.join(ROOT, 'outputs/tables/sdid_donor_robustness.csv')
results_df.to_csv(csv_path, index=False)
p(f'\n  Saved: sdid_donor_robustness.csv')

# ── NARRATIVE SUMMARY ─────────────────────────────────────────────────────────
THRESHOLD_ROBUST = 0.3
THRESHOLD_MOD    = 0.5
BASELINE_ATT_ACT = base['ATT']   # use actual computed baseline (not hardcoded)

p('\n' + '='*65)
p('  ROBUSTNESS NARRATIVE SUMMARY')
p('='*65)

alt_atts  = []
n_robust  = 0
n_mod     = 0
n_sens    = 0

for label, res in specs_ordered[1:]:   # skip baseline
    att   = res['ATT']
    delta = abs(att - BASELINE_ATT_ACT)
    alt_atts.append(att)

    if delta <= THRESHOLD_ROBUST:
        verdict = f'Robust. ATT = {att:+.2f} pp, consistent with baseline.'
        n_robust += 1
    elif delta <= THRESHOLD_MOD:
        verdict = f'Moderately sensitive. ATT = {att:+.2f} pp, broadly consistent but worth noting.'
        n_mod += 1
    else:
        verdict = f'Sensitive. ATT = {att:+.2f} pp, substantively different from baseline. Investigate.'
        n_sens += 1

    p(f'\n  {label}:')
    p(f'    {verdict}')
    if res['N'] <= 10:
        p(f'    [Small pool: N={res["N"]}, df={res["df"]}. SE and p-value have wide uncertainty.]')

n_total = len(alt_atts)
p()
p(f'  {'─'*55}')
p(f'  OVERALL ASSESSMENT  ({n_robust}/{n_total} within ±{THRESHOLD_ROBUST} pp of baseline)')
p(f'  {'─'*55}')

if n_robust == n_total:
    p('  Strong robustness. Headline result does not depend on donor pool composition.')
elif n_robust >= 6:
    p(f'  Generally robust. {n_total - n_robust} specification(s) show modest sensitivity.'
      f' Note in defense narrative.')
else:
    p('  Sensitive to donor pool. Rewrite narrative to acknowledge this.')

if alt_atts:
    att_min = min(alt_atts)
    att_max = max(alt_atts)
    if n_robust == n_total:
        robustness_word = 'robust'
    elif n_robust >= 6:
        robustness_word = 'moderately robust'
    else:
        robustness_word = 'sensitive'

    p()
    p(f'  COMMIT-READY SUMMARY STATEMENT:')
    p(f'  "Headline ATT of +{BASELINE_ATT_ACT:.2f} pp is {robustness_word} to donor pool specification.')
    p(f'  Alternative specifications produce estimates ranging from {att_min:+.2f} to {att_max:+.2f} pp,')
    p(f'  with {n_robust} of {n_total} falling within ±{THRESHOLD_ROBUST} pp of baseline."')

# ── FIGURE: Comparison table ──────────────────────────────────────────────────
p('\nBuilding comparison table figure...')

def fmt(val, fmt_str, na='—'):
    return fmt_str.format(val) if not (isinstance(val, float) and np.isnan(val)) else na

display_rows = []
for label, res in specs_ordered:
    att  = res['ATT']
    se   = res['SE']
    cilo = res['CI_lo']
    cihi = res['CI_hi']
    pv   = res['pval']
    rm   = res['RMSPE']

    display_rows.append((
        label,
        str(res['N']),
        fmt(att,  '{:+.2f}'),
        fmt(se,   '({:.2f})'),
        f'[{fmt(cilo, "{:+.2f}")}, {fmt(cihi, "{:+.2f}")}]' if not np.isnan(cilo) else '—',
        fmt(pv,   '{:.3f}'),
        fmt(rm,   '{:.3f}'),
    ))

# specs_ordered index: [0]=baseline, [1]=no-App, [2]=western, [3]=demo, [4-8]=LOO
# New display order: baseline → LOO (strongest evidence) → alt pools
# Subheader sentinels are 2-tuples (None, label_text); data rows are 7-tuples.
_SH = None   # sentinel tag for subheader rows
display_rows_final = (
    [display_rows[0]] +
    [(_SH, 'Leave-one-out (top 5 donors)')] +
    display_rows[4:9] +
    [(_SH, 'Alternative donor pools')] +
    display_rows[1:4]
)
# 1 baseline + 1 subhdr + 5 LOO + 1 subhdr + 3 alt = 11 rows
p(f'  [robustness table] total display rows: {len(display_rows_final)}  (expected 11)')

col_labels = ['Specification', 'N', 'ATT (pp)', '(SE)', '95% CI', 'Jack. p', 'RMSPE']
col_x      = [0.01, 0.33, 0.41, 0.49, 0.58, 0.76, 0.88]
col_ha     = ['left', 'center', 'center', 'center', 'center', 'center', 'center']

n_total = len(display_rows_final)              # 11
row_y   = np.linspace(0.90, 0.06, n_total + 1) # [0]=header, [1..11]=rows
step    = (0.90 - 0.06) / n_total

fig, ax = plt.subplots(figsize=(15, 8.0))
ax.axis('off')

# Header row
for j, (lbl, ha) in enumerate(zip(col_labels, col_ha)):
    ax.text(col_x[j], row_y[0], lbl, transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='center', ha=ha, color=C_REF)

ax.plot([0.01, 0.99], [row_y[0] - 0.042]*2,
        color=C_REF, lw=1.2, transform=ax.transAxes)

zebra_i = 0   # counts only data rows for alternating background
for i, row in enumerate(display_rows_final):
    y           = row_y[i + 1]
    is_baseline = (i == 0)
    is_subhdr   = (row[0] is _SH)

    if is_subhdr:
        # Dashed group separator sits between this row and the one above it
        sep_y = y + step * 0.5
        ax.plot([0.01, 0.99], [sep_y]*2,
                color='#BBBBBB', lw=0.8, ls='--', transform=ax.transAxes)
        # Italic gray group label — no background patch
        ax.text(col_x[0], y, row[1], transform=ax.transAxes,
                fontsize=10, va='center', ha='left',
                color='#666666', style='italic')
        continue   # subheader rows don't count toward zebra striping

    bg = '#EBF5FB' if is_baseline else ('#F4F6F7' if zebra_i % 2 == 0 else 'white')
    zebra_i += 1
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.01, y - 0.030), 0.98, 0.058,
        boxstyle='square,pad=0', linewidth=0,
        facecolor=bg, transform=ax.transAxes, zorder=0))

    for j, (val, ha) in enumerate(zip(row, col_ha)):
        weight = 'bold' if is_baseline else 'normal'
        color  = C_CA   if is_baseline else C_REF
        ax.text(col_x[j], y, val, transform=ax.transAxes,
                fontsize=11, va='center', ha=ha,
                color=color, fontweight=weight)

# Bottom rule
ax.plot([0.01, 0.99], [row_y[-1] - 0.036]*2,
        color=C_REF, lw=0.8, transform=ax.transAxes)

# Note
ax.text(0.5, -0.01,
        'All specifications use the restricted post-treatment window 2007\u20132015. '
        'Jackknife p-values reflect leave-one-unit-out inference. '
        'In-space placebo p-values are not reported for alternative specifications '
        'due to small donor pools. '
        'CI uses \u00b11.96\u00d7SE across all rows for direct comparability.',
        transform=ax.transAxes, fontsize=9, ha='center', va='top',
        color='#555555', style='italic', wrap=True)

ax.set_title(
    'Donor-Pool Robustness Checks \u2014 Synthetic DiD\n'
    'CA Paid Family Leave (2004), Effect on Female LFP Rate, Women 25\u201354',
    fontsize=FT, fontweight='bold', pad=14, color=C_REF)

fig.tight_layout()
for ext in ('png', 'pdf'):
    kw = dict(dpi=300) if ext == 'png' else {}
    fig.savefig(os.path.join(ROOT, f'outputs/figures/pres_robustness_table.{ext}'),
                bbox_inches='tight', **kw)
    p(f'  Saved: pres_robustness_table.{ext}')
plt.close(fig)

elapsed = time.time() - t_start
p(f'\nTotal runtime: {elapsed:.1f} seconds')
