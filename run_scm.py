import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # non-interactive backend (no display needed)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import minimize
import statsmodels.formula.api as smf
import warnings, sys
warnings.filterwarnings('ignore')

def p(*args): print(*args, flush=True)

np.random.seed(42)
p('Libraries loaded.')

# ── Parameters ──────────────────────────────────────────────────
agg        = pd.read_csv('state_year_predictors.csv')
TREATED    = 'California'
TREAT_YEAR = 2004
OUTCOME    = 'flfp_pct'
REF_YEAR   = 2003

years      = sorted(agg['YEAR'].unique())
pre_years  = [y for y in years if y < TREAT_YEAR]
post_years = [y for y in years if y >= TREAT_YEAR]
T0, T1, T  = len(pre_years), len(post_years), len(years)

p(f'Data loaded: {agg.shape[0]} rows | pre={pre_years} | post={post_years}')

# ── Donor pool ───────────────────────────────────────────────────
EXCLUDED     = ['New Jersey', 'District of Columbia', 'Alaska']
all_states   = sorted(agg['state'].unique())
donor_states = [s for s in all_states if s != TREATED and s not in EXCLUDED]
J = len(donor_states)
p(f'Donor pool: {J} states  (excluded: {EXCLUDED})')

# ── Predictor matrices ───────────────────────────────────────────
def build_predictors(state_name, data, pre_yrs):
    s_pre = data[(data['state'] == state_name) & (data['YEAR'].isin(pre_yrs))].sort_values('YEAR')
    return np.concatenate([
        s_pre['flfp_pct'].values,
        [s_pre['college_pct'].mean(), s_pre['married_pct'].mean(), s_pre['mean_age'].mean()]
    ])

predictor_names = [f'FLFP_{y}' for y in pre_years] + ['Avg_College', 'Avg_Married', 'Avg_Age']
K  = len(predictor_names)

X1 = build_predictors(TREATED, agg, pre_years)
X0 = np.column_stack([build_predictors(s, agg, pre_years) for s in donor_states])

X_scale        = X0.std(axis=1, ddof=1)
X_scale[X_scale < 1e-8] = 1.0
X1_s, X0_s     = X1 / X_scale, X0 / X_scale[:, None]

Y_matrix  = agg.pivot(index='state', columns='YEAR', values=OUTCOME)
Y_matrix  = Y_matrix.reindex(index=[TREATED] + donor_states, columns=years)
Y_treated = Y_matrix.loc[TREATED, :].values
Y_donors  = Y_matrix.loc[donor_states, :].values
Y1_pre    = Y_treated[:T0]
Y0_pre    = Y_donors[:, :T0]

p(f'Predictor matrix built: K={K}, J={J}')

# ── Optimization functions ───────────────────────────────────────
def get_W(V_diag, X1_s, X0_s):
    V  = np.diag(V_diag)
    J_ = X0_s.shape[1]
    def obj(w):  return (X1_s - X0_s @ w) @ V @ (X1_s - X0_s @ w)
    def grad(w): return -2 * X0_s.T @ V @ (X1_s - X0_s @ w)
    res = minimize(obj, np.ones(J_)/J_, jac=grad, method='SLSQP',
                   bounds=[(0,1)]*J_,
                   constraints={'type':'eq','fun':lambda w: w.sum()-1,'jac':lambda w:np.ones(J_)},
                   options={'ftol':1e-12,'maxiter':2000})
    return res.x

def outer_loss(V_raw, X1_s, X0_s, Y1_pre, Y0_pre):
    V_diag = np.abs(V_raw) / (np.abs(V_raw).sum() + 1e-12)
    W      = get_W(V_diag, X1_s, X0_s)
    return np.mean((Y1_pre - Y0_pre.T @ W)**2)

# ── Main optimization (30 starts) ───────────────────────────────
p('Optimizing V weights (30 multi-start runs)...')
best_loss, best_V_raw = np.inf, None
for i in range(30):
    v0  = np.random.dirichlet(np.ones(K))
    res = minimize(outer_loss, v0, args=(X1_s, X0_s, Y1_pre, Y0_pre),
                   method='Nelder-Mead',
                   options={'xatol':1e-10,'fatol':1e-10,'maxiter':10000,'maxfev':20000})
    if res.fun < best_loss:
        best_loss, best_V_raw = res.fun, res.x

V_opt = np.abs(best_V_raw) / np.abs(best_V_raw).sum()
W_opt = get_W(V_opt, X1_s, X0_s)
p(f'Optimization done. Best MSPE={best_loss:.6f}')

# ── Weights & balance ────────────────────────────────────────────
v_df = pd.DataFrame({'Predictor': predictor_names, 'V_weight': V_opt}).sort_values('V_weight', ascending=False)
p('\nPREDICTOR IMPORTANCE WEIGHTS (V)')
p(v_df.to_string(index=False))

w_df = pd.DataFrame({'State': donor_states, 'Weight': W_opt})
w_df = w_df[w_df['Weight'] > 0.001].sort_values('Weight', ascending=False)
p('\nDONOR STATE WEIGHTS (Synthetic CA composition)')
p(w_df.to_string(index=False))

synth_X_orig = X0 @ W_opt
balance = pd.DataFrame({
    'Predictor'  : predictor_names,
    'California' : X1,
    'Synthetic'  : synth_X_orig,
    'Donor Mean' : X0.mean(axis=1),
    'CA-Synth'   : X1 - synth_X_orig
})
p('\nPREDICTOR BALANCE TABLE')
p(balance.to_string(index=False, float_format='{:.3f}'.format))

# ── Synthetic CA & gap ───────────────────────────────────────────
synth_CA  = Y_donors.T @ W_opt
gap       = Y_treated - synth_CA
pre_mask  = np.array([y < TREAT_YEAR for y in years])
post_mask = ~pre_mask

pre_rmspe    = np.sqrt(np.mean(gap[pre_mask]**2))
post_rmspe   = np.sqrt(np.mean(gap[post_mask]**2))
mspe_ratio_ca = (post_rmspe**2) / max(pre_rmspe**2, 1e-8)
att_avg      = gap[post_mask].mean()

results_df = pd.DataFrame({
    'Year': years, 'CA_actual': Y_treated, 'Synthetic_CA': synth_CA,
    'Gap': gap, 'Post': post_mask.astype(int)
})
p('\nYEAR-BY-YEAR RESULTS')
p(results_df.to_string(index=False, float_format='{:.3f}'.format))
p(f'\nPre-RMSPE : {pre_rmspe:.4f} pp')
p(f'Post-RMSPE: {post_rmspe:.4f} pp')
p(f'MSPE ratio: {mspe_ratio_ca:.2f}')
p(f'Avg ATT   : {att_avg:+.4f} pp')

# ── Plot 1: Actual vs Synthetic ──────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
ax = axes[0]
ax.plot(years, Y_treated, 'o-', color='red',   lw=2.5, ms=8, label='California (Actual)')
ax.plot(years, synth_CA,  's--', color='black', lw=2.5, ms=8, label='Synthetic California')
ax.axvspan(TREAT_YEAR-0.5, max(years)+0.5, alpha=0.07, color='blue')
ax.axvline(TREAT_YEAR-0.5, color='blue', lw=1.5, linestyle=':')
ax.set_xlabel('Year', fontweight='bold'); ax.set_ylabel('Female LFP Rate (%)', fontweight='bold')
ax.set_title('Synthetic Control: CA vs Synthetic CA', fontweight='bold')
ax.legend(); ax.grid(alpha=0.3); ax.set_xticks(years)

ax = axes[1]
bar_colors = ['steelblue' if y < TREAT_YEAR else ('green' if g >= 0 else 'tomato')
              for y, g in zip(years, gap)]
ax.bar(years, gap, color=bar_colors, edgecolor='black', lw=0.5, alpha=0.85, width=0.7)
ax.axhline(0, color='black', lw=1.5)
ax.axvline(TREAT_YEAR-0.5, color='blue', lw=1.5, linestyle=':')
for y, g in zip(years, gap):
    ax.text(y, g+(0.1 if g>=0 else -0.3), f'{g:+.2f}', ha='center', fontsize=8.5, fontweight='bold')
ax.set_xlabel('Year', fontweight='bold'); ax.set_ylabel('Gap: Actual − Synthetic (pp)', fontweight='bold')
ax.set_title(f'Treatment Effect Gap  |  Avg ATT: {att_avg:+.2f} pp', fontweight='bold')
ax.grid(alpha=0.3, axis='y'); ax.set_xticks(years)
plt.tight_layout()
plt.savefig('scm_main_results.png', dpi=150, bbox_inches='tight')
p('Saved: scm_main_results.png')
plt.close()

# ── DiD regression ───────────────────────────────────────────────
panel_df = pd.DataFrame({
    'year': years*2,
    'flfp': np.concatenate([Y_treated, synth_CA]),
    'ca'  : [1]*T + [0]*T,
    'post': [1 if y >= TREAT_YEAR else 0 for y in years]*2
})
panel_df['ca_post'] = panel_df['ca'] * panel_df['post']

did_model = smf.ols('flfp ~ ca + post + ca_post', data=panel_df).fit(cov_type='HC3')
att_did = did_model.params['ca_post']
se_did  = did_model.bse['ca_post']
pval_did = did_model.pvalues['ca_post']

p('\nDIFFERENCE-IN-DIFFERENCES RESULTS')
p(did_model.summary().tables[1])
p(f'\n*** DiD ATT: {att_did:+.4f} pp  |  SE={se_did:.4f}  |  p={pval_did:.4f} ***')

# ── Event study ──────────────────────────────────────────────────
ref_gap  = results_df.loc[results_df['Year']==REF_YEAR, 'Gap'].values[0]
es_coefs = results_df['Gap'].values - ref_gap
p('\nEVENT STUDY COEFFICIENTS (gap_t - gap_2003)')
for y, c, post in zip(years, es_coefs, post_mask):
    tag = 'POST' if post else 'pre '
    ref = ' <- ref' if y == REF_YEAR else ''
    p(f'  {y} [{tag}]: {c:+.4f} pp{ref}')

# ── Placebo tests (5 starts per state — fast) ─────────────────────
p('\nRunning placebo tests (47 states, 5 starts each)...')
placebo_gaps, placebo_prermspe, placebo_mspe_ratio = {}, {}, {}

for idx, p_state in enumerate(donor_states):
    p_donors = [s for s in all_states if s != p_state and s not in EXCLUDED]
    p_Y_tr   = Y_matrix.loc[p_state, :].values
    p_Y_don  = Y_matrix.loc[p_donors, :].values
    p_Y1_pre = p_Y_tr[:T0]
    p_Y0_pre = p_Y_don[:, :T0]

    p_X1 = build_predictors(p_state, agg, pre_years)
    p_X0 = np.column_stack([build_predictors(s, agg, pre_years) for s in p_donors])
    p_sc = p_X0.std(axis=1, ddof=1); p_sc[p_sc < 1e-8] = 1.0
    p_X1_s, p_X0_s = p_X1/p_sc, p_X0/p_sc[:,None]

    best_m, best_v = np.inf, None
    for _ in range(5):
        v0  = np.random.dirichlet(np.ones(K))
        res = minimize(outer_loss, v0, args=(p_X1_s, p_X0_s, p_Y1_pre, p_Y0_pre),
                       method='Nelder-Mead',
                       options={'xatol':1e-7,'fatol':1e-7,'maxiter':3000})
        if res.fun < best_m:
            best_m, best_v = res.fun, res.x

    p_V    = np.abs(best_v) / np.abs(best_v).sum()
    p_W    = get_W(p_V, p_X1_s, p_X0_s)
    p_synth = p_Y_don.T @ p_W
    p_gap   = p_Y_tr - p_synth

    pr = np.sqrt(np.mean(p_gap[pre_mask]**2))
    po = np.sqrt(np.mean(p_gap[post_mask]**2))
    placebo_gaps[p_state]       = p_gap
    placebo_prermspe[p_state]   = pr
    placebo_mspe_ratio[p_state] = (po**2) / max(pr**2, 1e-8)

    if (idx+1) % 10 == 0:
        p(f'  {idx+1}/{J} placebos done...')

# Add California
placebo_gaps[TREATED]        = gap
placebo_prermspe[TREATED]    = pre_rmspe
placebo_mspe_ratio[TREATED]  = mspe_ratio_ca
p('Placebo tests complete.')

# ── Inference ────────────────────────────────────────────────────
cutoff      = 2.0 * pre_rmspe
good_states = [s for s, r in placebo_prermspe.items() if r <= cutoff or s == TREATED]
donor_ratios = [placebo_mspe_ratio[s] for s in good_states if s != TREATED]
pvalue       = np.mean([r >= mspe_ratio_ca for r in donor_ratios])

p(f'\nINFERENCE')
p(f'  Pre-RMSPE cutoff (2× CA): {cutoff:.4f} pp')
p(f'  Placebos included: {len(good_states)-1}')
p(f'  CA MSPE ratio: {mspe_ratio_ca:.2f}')
p(f'  Permutation p-value: {pvalue:.4f}')

# Placebo event study bands
p_es = {}
ref_idx = years.index(REF_YEAR)
for s in good_states:
    if s != TREATED:
        p_es[s] = placebo_gaps[s] - placebo_gaps[s][ref_idx]
p_es_array = np.array(list(p_es.values()))
ci_lo = np.percentile(p_es_array, 5,  axis=0)
ci_hi = np.percentile(p_es_array, 95, axis=0)

# ── Plot 2: Placebo + MSPE + Event Study ──────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(21, 7))

ax = axes[0]
for s in good_states:
    if s != TREATED:
        ax.plot(years, placebo_gaps[s], color='lightgray', lw=0.8, alpha=0.6)
ax.plot(years, gap, 'o-', color='red', lw=2.5, ms=7, zorder=5, label='California')
ax.axhline(0, color='black', lw=1.5)
ax.axvline(TREAT_YEAR-0.5, color='blue', lw=1.5, linestyle=':')
ax.axvspan(TREAT_YEAR-0.5, max(years)+0.5, alpha=0.07, color='blue')
ax.set_title('In-Space Placebo Test', fontweight='bold')
ax.set_xlabel('Year', fontweight='bold'); ax.set_ylabel('Gap (pp)', fontweight='bold')
ax.legend(); ax.grid(alpha=0.3); ax.set_xticks(years)

ax = axes[1]
all_r  = sorted([(placebo_mspe_ratio[s], s) for s in good_states])
ratios = [r for r,_ in all_r]
colors = ['red' if s == TREATED else 'steelblue' for _,s in all_r]
ax.barh(range(len(ratios)), ratios, color=colors, edgecolor='black', lw=0.4, alpha=0.85)
ax.axvline(mspe_ratio_ca, color='red', lw=2, linestyle='--',
           label=f'CA: {mspe_ratio_ca:.1f}\np={pvalue:.3f}')
ax.set_title('MSPE Ratio Test', fontweight='bold')
ax.set_xlabel('Post/Pre MSPE Ratio', fontweight='bold')
ax.legend(); ax.grid(alpha=0.3, axis='x')

ax = axes[2]
ax.fill_between(years, ci_lo, ci_hi, alpha=0.2, color='gray', label='90% placebo band')
ax.axhline(0, color='black', lw=1.5)
ax.axvline(TREAT_YEAR-0.5, color='blue', lw=1.5, linestyle=':')
ax.axvspan(TREAT_YEAR-0.5, max(years)+0.5, alpha=0.06, color='blue')
pre_idx_  = [i for i,y in enumerate(years) if y < TREAT_YEAR and y != REF_YEAR]
post_idx_ = [i for i,y in enumerate(years) if y >= TREAT_YEAR]
ax.plot([years[i] for i in pre_idx_],  es_coefs[pre_idx_],  'o', color='steelblue', ms=9, zorder=5, label='Pre-treatment')
ax.plot(years[ref_idx], es_coefs[ref_idx], 's', color='black', ms=10, zorder=5, label=f'Ref ({REF_YEAR})')
ax.plot([years[i] for i in post_idx_], es_coefs[post_idx_], 'D', color='red',      ms=9, zorder=5, label='Post-treatment')
ax.plot(years, es_coefs, '-', color='black', lw=1.5, alpha=0.5)
for y, c in zip(years, es_coefs):
    ax.text(y, c+(0.12 if c>=0 else -0.32), f'{c:+.2f}', ha='center', fontsize=8, fontweight='bold')
ax.set_title(f'Event Study  |  p={pvalue:.3f}', fontweight='bold')
ax.set_xlabel('Year', fontweight='bold'); ax.set_ylabel(f'Effect rel. to {REF_YEAR} (pp)', fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_xticks(years)

plt.tight_layout()
plt.savefig('scm_placebo_eventstudy.png', dpi=150, bbox_inches='tight')
p('Saved: scm_placebo_eventstudy.png')
plt.close()

# ── Save CSVs ────────────────────────────────────────────────────
results_df.to_csv('scm_results.csv', index=False)
w_df.to_csv('scm_weights.csv', index=False)
mspe_df = pd.DataFrame([{'State': s, 'pre_rmspe': placebo_prermspe[s],
                          'mspe_ratio': placebo_mspe_ratio[s]} for s in placebo_gaps]
                       ).sort_values('mspe_ratio', ascending=False)
mspe_df.to_csv('scm_placebo_mspe.csv', index=False)

# ── Final summary ────────────────────────────────────────────────
p('\n' + '='*65)
p('  FINAL SUMMARY')
p('='*65)
p(f'  Pre-treatment RMSPE  : {pre_rmspe:.4f} pp')
p(f'  Avg. ATT (2004-2010) : {att_avg:+.4f} pp')
p(f'  DiD ATT (OLS)        : {att_did:+.4f} pp  (SE={se_did:.4f}, p={pval_did:.4f})')
p(f'  MSPE ratio (CA)      : {mspe_ratio_ca:.2f}')
p(f'  Permutation p-value  : {pvalue:.4f}')
if pvalue <= 0.05:  p('  → SIGNIFICANT at 5%')
elif pvalue <= 0.10: p('  → Significant at 10%')
else:               p('  → Not significant at conventional levels')
p(f'\n  Top donor states:')
for _, row in w_df.head(5).iterrows():
    p(f"    {row['State']:22s}: {row['Weight']*100:.1f}%")
p('\nFiles saved: scm_results.csv | scm_weights.csv | scm_placebo_mspe.csv')
p('Plots saved: scm_main_results.png | scm_placebo_eventstudy.png')
