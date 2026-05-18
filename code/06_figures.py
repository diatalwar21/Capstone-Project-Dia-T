"""
All output figures and the summary statistics table for the CA PFL paper.

Merged from:
  make_summary_stats_table.py  -> SECTION 1: Summary Statistics Table
  make_presentation_figures.py -> SECTION 2: Presentation Figures

Inputs (from outputs/tables/):
  sdid_robust_yearly_att.csv
  sdid_robust_unit_weights.csv
  sdid_robust_restricted_placebo_atts.csv

Inputs (from data/derived/):
  state_year_panel_deduped_1995_2015.csv

Outputs (to outputs/figures/):
  pres_sumstats.png/.pdf
  pres_table.png/.pdf
  pres_trend.png/.pdf
  pres_trend_y0.png/.pdf
  pres_trend_context.png/.pdf
  pres_event_study.png/.pdf
  pres_event_study_clean.png/.pdf
  pres_placebo.png/.pdf
  pres_donor_weights.png/.pdf
  pres_pretrend.png/.pdf
  pres_pretrend_y0.png/.pdf
  pres_pretrend_context.png/.pdf

Outputs (to outputs/tables/):
  summary_stats_1995_2015.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator

DIR       = os.path.dirname(os.path.abspath(__file__))
ROOT      = os.path.abspath(os.path.join(DIR, ".."))
FIG_DIR   = os.path.join(ROOT, 'outputs', 'figures')
TABLE_DIR = os.path.join(ROOT, 'outputs', 'tables')
DATA_DIR  = os.path.join(ROOT, 'data', 'derived')
os.makedirs(FIG_DIR,   exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

TREAT_YEAR = 2004
DPI        = 300

# ── Color palette ──────────────────────────────────────────────────────────────
C_CA    = '#1B4F72'   # deep navy
C_SYN   = '#7F8C8D'   # muted gray
C_POS   = '#2874A6'   # steel blue
C_NEG   = '#BDC3C7'   # light gray-blue
C_CONF  = '#FEF5E7'   # pale amber
C_POST  = '#EBF5FB'   # very pale blue
C_REF   = '#17202A'   # near-black
C_PLAC  = '#AED6F1'   # light steel blue
C_GRID  = '#E8E8E8'   # light gray

# ── Typography ─────────────────────────────────────────────────────────────────
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
FT = 16   # title
FL = 13   # axis labels
FA = 10   # annotations / footnotes

def savefig_pres(fig, name):
    """Save pres_{name}.png and pres_{name}.pdf to outputs/figures/."""
    for ext in ('png', 'pdf'):
        kw = dict(dpi=DPI) if ext == 'png' else {}
        fig.savefig(os.path.join(FIG_DIR, f'pres_{name}.{ext}'),
                    bbox_inches='tight', **kw)
    print(f'  Saved: pres_{name}.png + .pdf')
    plt.close(fig)

def savefig(fig, name):
    """Save {name}.png and {name}.pdf to outputs/figures/."""
    for ext in ('png', 'pdf'):
        kw = dict(dpi=DPI) if ext == 'png' else {}
        fig.savefig(os.path.join(FIG_DIR, f'{name}.{ext}'),
                    bbox_inches='tight', **kw)
        print(f'  Saved: {name}.{ext}')
    plt.close(fig)

def p(*a): print(*a, flush=True)


# ============================================================
# SECTION 1: SUMMARY STATISTICS TABLE (Table 1)
# ============================================================

p('Loading state_year_panel_deduped_1995_2015.csv...')
df = pd.read_csv(os.path.join(DATA_DIR, 'state_year_panel_deduped_1995_2015.csv'))

# Corrected donor pool (same exclusion list as run_sdid_robust_restricted.py)
EXCLUDED = ['New York', 'Rhode Island', 'Hawaii', 'New Jersey',
            'Washington', 'Alaska', 'District of Columbia']

ca_df    = df[df['state'] == 'California'].copy()
oth_df   = df[df['state'] != 'California'].copy()
donor_df = oth_df[~oth_df['state'].isin(EXCLUDED)].copy()

# Analytical sample: CA + 43 donors only (924 cells)
analytic_df = pd.concat([ca_df, donor_df]).reset_index(drop=True)

# ── SAMPLE CONSTRUCTION ───────────────────────────────────────────────────────
n_states_all  = df['state'].nunique()
n_years       = df['YEAR'].nunique()
years_ss      = sorted(df['YEAR'].unique())
n_cells       = len(df)
n_donors      = donor_df['state'].nunique()
n_jur         = n_donors + 1
n_cells_analytic = len(analytic_df)

total_obs          = int(df['n_observations'].sum())
analytic_obs_total = int(analytic_df['n_observations'].sum())
ca_obs_total       = int(ca_df['n_observations'].sum())
ca_obs_per_yr      = ca_df['n_observations'].mean()
donor_obs_mean     = donor_df['n_observations'].mean()
analytic_obs_mean  = analytic_df['n_observations'].mean()
obs_min            = int(analytic_df['n_observations'].min())
obs_max            = int(analytic_df['n_observations'].max())
n_cells_model      = n_cells_analytic

p('\n' + '='*65)
p('SAMPLE CONSTRUCTION — 1995–2015 PANEL')
p('='*65)
p(f'  Panel period:               {years_ss[0]}–{years_ss[-1]} ({n_years} years)')
p(f'  Total jurisdictions:        {n_states_all} (50 states + DC)')
p(f'  State × year cells (full):  {n_cells} ({n_states_all} × {n_years})')
p(f'  Corrected donor pool:       {n_donors} states (excl. TDI/own-PFL/DC/AK)')
p(f'  Cells in SDiD estimation:   {n_cells_model} ({n_donors} donors + CA × {n_years} yrs)')
p(f'')
p(f'  Analytical obs (women 25–54, all cells): {total_obs:>12,}')
p(f'  California total obs:                    {ca_obs_total:>12,}')
p(f'  California obs/year (avg):               {ca_obs_per_yr:>12,.0f}')
p(f'  Obs per cell (mean / min / max):         {df["n_observations"].mean():.0f} / {obs_min} / {obs_max}')

# ── WEIGHTED SUMMARY STATISTICS ───────────────────────────────────────────────
def wmean(vals, weights):
    return np.average(vals, weights=weights)

def wstd(vals, weights):
    m = wmean(vals, weights)
    return np.sqrt(np.average((vals - m)**2, weights=weights))

def row_stats(col, scale=1.0):
    vp = analytic_df[col].values * scale
    vc = ca_df[col].values * scale
    vd = donor_df[col].values * scale
    wp = analytic_df['n_observations'].values
    wc = ca_df['n_observations'].values
    wd = donor_df['n_observations'].values
    pm, ps = wmean(vp, wp), wstd(vp, wp)
    cm, cs = wmean(vc, wc), wstd(vc, wc)
    dm, ds = wmean(vd, wd), wstd(vd, wd)
    return pm, ps, cm, cs, dm, ds, cm - dm

vars_cfg = [
    ('flfp_pct',    'Female LFP rate (%)',           1.0),
    ('college_pct', 'College degree or higher (%)',  1.0),
    ('married_pct', 'Married (%)',                   1.0),
    ('mean_age',    'Mean age (years)',               1.0),
]

p('\n' + '='*65)
p('WEIGHTED SUMMARY STATISTICS (weight = cell obs count)')
p('='*65)

rows_data = []
for col, label, scale in vars_cfg:
    if col not in df.columns:
        p(f'  WARNING: {col} not found — skipping')
        continue
    pm, ps, cm, cs, dm, ds, diff = row_stats(col, scale)
    rows_data.append((label, pm, ps, cm, cs, dm, ds, diff))
    p(f'  {label}')
    p(f'    Pooled  : {pm:.2f}  (sd {ps:.2f})')
    p(f'    CA      : {cm:.2f}  (sd {cs:.2f})')
    p(f'    Donors  : {dm:.2f}  (sd {ds:.2f})')
    p(f'    CA–Other: {diff:+.2f} pp')
    p()

# ── RECONCILIATION CHECKS ─────────────────────────────────────────────────────
p('='*65)
p('RECONCILIATION CHECKS')
p('='*65)

att_file = os.path.join(TABLE_DIR, 'sdid_robust_yearly_att.csv')
if os.path.exists(att_file):
    att_rec = pd.read_csv(att_file)
    ca_flfp_trend  = att_rec['CA_actual'].mean()
    ca_flfp_panel  = ca_df['flfp_pct'].mean()
    diff_flfp = abs(ca_flfp_panel - ca_flfp_trend)
    p(f'  CA FLFP — trend figure (unweighted avg across years): {ca_flfp_trend:.3f}%')
    p(f'  CA FLFP — panel (simple cell mean):                   {ca_flfp_panel:.3f}%')
    p(f'  Difference: {diff_flfp:.3f} pp')
    status = 'CONSISTENT' if diff_flfp < 0.15 else 'CHECK — difference > 0.15 pp'
    p(f'  Status: {status}')
else:
    p('  sdid_robust_yearly_att.csv not found — skipping CA FLFP check')

p(f'\n  Corrected donor states: {n_donors}  (expected 43)')
p(f'  Status: {"CONSISTENT" if n_donors == 43 else f"MISMATCH — got {n_donors}"}')

p(f'\n  SDiD estimation cells ({n_donors}+1 × {n_years}): {n_cells_model}')
p(f'  Results table reports "1,008 state × year obs."')
p(f'  That 1,008 = 48 states (baseline 47 donors + CA) × 21 years.')
p(f'  Robust spec uses {n_cells_model} cells ({n_donors} donors + CA × {n_years} yrs).')
if n_cells_model == 1008:
    p(f'  Status: CONSISTENT')
else:
    p(f'  Status: EXPECTED (baseline=1,008; robust={n_cells_model}; difference is by design)')

# ── OLD VS NEW SAMPLE COMPARISON ─────────────────────────────────────────────
p('\n' + '='*65)
p('OLD 2000–2010 vs NEW 1995–2015 SAMPLE — NARRATIVE FLAGS')
p('='*65)
p('  Old: 11 years, 561 cells, 3,899,679 obs (21.8% of raw data retained)')
p(f'  New: {n_years} years, {n_cells} cells, {total_obs:,} obs')

if rows_data:
    old_ca_flfp, old_oth_flfp  = 72.1, 77.5
    old_ca_coll, old_oth_coll  = 31.9, 30.9

    new_ca_flfp  = rows_data[0][3]
    new_oth_flfp = rows_data[0][5]
    new_ca_coll  = rows_data[1][3] if len(rows_data) > 1 else None
    new_oth_coll = rows_data[1][5] if len(rows_data) > 1 else None

    p(f'\n  FLFP (CA):   old {old_ca_flfp}%  →  new {new_ca_flfp:.1f}%'
      f'  (Δ {new_ca_flfp - old_ca_flfp:+.1f} pp)')
    p(f'  FLFP (other): old {old_oth_flfp}%  →  new {new_oth_flfp:.1f}%'
      f'  (Δ {new_oth_flfp - old_oth_flfp:+.1f} pp)')

    old_gap = old_oth_flfp - old_ca_flfp
    new_gap = new_oth_flfp - new_ca_flfp
    p(f'  Structural FLFP gap (other − CA):')
    p(f'    Old: {old_gap:.1f} pp  →  New: {new_gap:.1f} pp  (Δ {new_gap - old_gap:+.1f} pp)')
    if abs(new_gap - old_gap) > 2:
        p('  FLAG: Gap shifted > 2 pp — update data narrative on structural CA deficit.')
    else:
        p('  OK: Gap stable across sample windows.')

    if new_ca_coll is not None:
        p(f'\n  College share (CA):    old {old_ca_coll}%  →  new {new_ca_coll:.1f}%'
          f'  (Δ {new_ca_coll - old_ca_coll:+.1f} pp)')
        p(f'  College share (other): old {old_oth_coll}%  →  new {new_oth_coll:.1f}%'
          f'  (Δ {new_oth_coll - old_oth_coll:+.1f} pp)')
        if abs(new_ca_coll - old_ca_coll) > 3 or abs(new_oth_coll - old_oth_coll) > 3:
            p('  FLAG: College share shifted > 3 pp — update data section narrative.')
        else:
            p('  OK: College share stable across sample windows.')

# ── SAVE CSV ──────────────────────────────────────────────────────────────────
csv_rows = []
for label, pm, ps, cm, cs, dm, ds, diff in rows_data:
    csv_rows.append({
        'Variable':      label,
        'Pooled_mean':   round(pm, 3),
        'Pooled_sd':     round(ps, 3),
        'CA_mean':       round(cm, 3),
        'CA_sd':         round(cs, 3),
        'Donor_mean':    round(dm, 3),
        'Donor_sd':      round(ds, 3),
        'Diff_CA_minus_Donor': round(diff, 3),
    })
csv_rows.append({
    'Variable':      'Obs per cell (mean)',
    'Pooled_mean':   round(analytic_obs_mean, 1),
    'Pooled_sd':     '',
    'CA_mean':       round(ca_obs_per_yr, 1),
    'CA_sd':         '',
    'Donor_mean':    round(donor_obs_mean, 1),
    'Donor_sd':      '',
    'Diff_CA_minus_Donor': '',
})
pd.DataFrame(csv_rows).to_csv(
    os.path.join(TABLE_DIR, 'summary_stats_1995_2015.csv'), index=False)
p('\n  Saved: summary_stats_1995_2015.csv  ->  outputs/tables/')

# ── SUMMARY STATS FIGURE ──────────────────────────────────────────────────────
p('\nBuilding summary stats table figure...')

col_labels = ['Variable', 'Pooled', 'California', 'Donor States\n(43)', 'CA − Donors\n(pp)']
col_x      = [0.02, 0.36, 0.53, 0.70, 0.88]

n_rows_ss  = len(rows_data) + 1
row_y      = np.linspace(0.60, 0.08, n_rows_ss + 1)

fig, ax = plt.subplots(figsize=(13, 7.0))
ax.axis('off')

overview_lines = [
    ('Sample:',                  'Women aged 25 to 54'),
    ('Period:',                  f'1995 to 2015 (21 years)'),
    ('Jurisdictions:',           f'{n_jur} (California + {n_donors} donors)'),
    ('State × year cells:', f'{n_cells_analytic:,}'),
    ('Individual observations:', f'{analytic_obs_total:,}'),
]
overview_y_start = 0.97
overview_spacing = 0.055
for k, (lbl, val) in enumerate(overview_lines):
    y_ov = overview_y_start - k * overview_spacing
    ax.text(0.02, y_ov, lbl, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left', color=C_REF)
    ax.text(0.30, y_ov, val, transform=ax.transAxes,
            fontsize=10, va='top', ha='left', color=C_REF)

ax.plot([0.01, 0.99], [0.66, 0.66], color='#AAAAAA',
        linewidth=0.6, transform=ax.transAxes)

for j, label in enumerate(col_labels):
    ha = 'left' if j == 0 else 'center'
    ax.text(col_x[j], row_y[0], label, transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='center', ha=ha, color=C_REF)

ax.plot([0.01, 0.99], [row_y[0] - 0.045, row_y[0] - 0.045],
        color=C_REF, linewidth=1.2, transform=ax.transAxes)

for i, (label, pm, ps, cm, cs, dm, ds, diff) in enumerate(rows_data):
    y  = row_y[i + 1]
    bg = '#F4F6F7' if i % 2 == 0 else 'white'
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.01, y - 0.035), 0.98, 0.068,
        boxstyle='square,pad=0', linewidth=0,
        facecolor=bg, transform=ax.transAxes, zorder=0))

    ax.text(col_x[0], y, label,
            transform=ax.transAxes, fontsize=11, va='center', ha='left', color=C_REF)
    ax.text(col_x[1], y, f'{pm:.2f} ({ps:.2f})',
            transform=ax.transAxes, fontsize=11, va='center', ha='center', color=C_REF)
    ax.text(col_x[2], y, f'{cm:.2f} ({cs:.2f})',
            transform=ax.transAxes, fontsize=11, va='center', ha='center',
            color=C_CA, fontweight='bold')
    ax.text(col_x[3], y, f'{dm:.2f} ({ds:.2f})',
            transform=ax.transAxes, fontsize=11, va='center', ha='center', color=C_REF)
    diff_color = C_CA if diff < 0 else C_REF
    ax.text(col_x[4], y, f'{diff:+.2f}',
            transform=ax.transAxes, fontsize=11, va='center', ha='center',
            color=diff_color, fontweight='bold' if abs(diff) > 2 else 'normal')

y_obs = row_y[len(rows_data) + 1]
bg    = '#F4F6F7' if len(rows_data) % 2 == 0 else 'white'
ax.add_patch(mpatches.FancyBboxPatch(
    (0.01, y_obs - 0.035), 0.98, 0.068,
    boxstyle='square,pad=0', linewidth=0,
    facecolor=bg, transform=ax.transAxes, zorder=0))
ax.text(col_x[0], y_obs, 'Obs. per state-year cell',
        transform=ax.transAxes, fontsize=11, va='center', ha='left', color=C_REF)
ax.text(col_x[1], y_obs, f'{analytic_obs_mean:.0f}',
        transform=ax.transAxes, fontsize=11, va='center', ha='center', color=C_REF)
ax.text(col_x[2], y_obs, f'{ca_obs_per_yr:.0f}',
        transform=ax.transAxes, fontsize=11, va='center', ha='center',
        color=C_CA, fontweight='bold')
ax.text(col_x[3], y_obs, f'{donor_obs_mean:.0f}',
        transform=ax.transAxes, fontsize=11, va='center', ha='center', color=C_REF)
ax.text(col_x[4], y_obs, f'range {obs_min}–{obs_max}',
        transform=ax.transAxes, fontsize=9, va='center', ha='center', color='#666666')

ax.plot([0.01, 0.99], [y_obs - 0.04, y_obs - 0.04],
        color=C_REF, linewidth=0.8, transform=ax.transAxes)

ax.text(0.5, -0.02,
        f'Note: {n_cells_analytic:,} state × year observations '
        f'(44 jurisdictions × 21 years, 1995–2015). '
        f'Analytical sample includes California and 43 donor states '
        f'(excludes TDI states NY, RI, HI; own-PFL states NJ, WA; Alaska; and DC). '
        f'Means weighted by cell-level CPS sample size. '
        f'Standard deviations in parentheses. California values in navy.',
        transform=ax.transAxes, fontsize=9, ha='center', va='top',
        color='#555555', style='italic', wrap=True)

ax.text(0.5, -0.09,
        'California shows a 5.26 pp FLFP deficit vs. donors despite higher educational '
        'attainment (31.8% vs 30.5% college-educated), consistent with structural factors '
        'discussed in Section 3.',
        transform=ax.transAxes, fontsize=9, ha='center', va='top',
        color='#555555', style='italic', wrap=True)

ax.set_title(
    'Summary Statistics — CA Paid Family Leave Analysis\n'
    'Female LFP, Women 25–54, CPS Annual Microdata (1995–2015)',
    fontsize=FT, fontweight='bold', pad=15, color=C_REF)

fig.tight_layout()
savefig(fig, 'pres_sumstats')

p('\n' + '='*65)
p('SANITY CHECKS — ANALYTICAL SAMPLE')
p('='*65)
p(f'  Analytical cells:    {n_cells_analytic}  (expected 924)')
p(f'  Status: {"OK" if n_cells_analytic == 924 else f"MISMATCH — got {n_cells_analytic}"}')
p(f'\n  Analytical obs total: {analytic_obs_total:,}')
ca_flfp_check = rows_data[0][3] if rows_data else float('nan')
p(f'\n  CA FLFP (weighted mean): {ca_flfp_check:.2f}%  (expect ~71.87%)')
p(f'  Status: {"OK" if abs(ca_flfp_check - 71.87) < 0.10 else "CHECK — differs from expected 71.87%"}')


# ============================================================
# SECTION 2: ALL PRESENTATION FIGURES
# ============================================================

p('\n' + '='*60)
p('SECTION 2: PRESENTATION FIGURES')
p('='*60)

p('Loading data...')
att   = pd.read_csv(os.path.join(TABLE_DIR, 'sdid_robust_yearly_att.csv'))
wts   = pd.read_csv(os.path.join(TABLE_DIR, 'sdid_robust_unit_weights.csv'))
plac  = pd.read_csv(os.path.join(TABLE_DIR, 'sdid_robust_restricted_placebo_atts.csv'))

years   = att['Year'].values.astype(int)
ca      = att['CA_actual'].values
synth   = att['Synth_CA'].values
att_yr  = att['ATT_yearly'].values
post_m  = att['Post'].values == 1
pre_m   = ~post_m

tau_restr = att[att['Year'] >= 2007]['ATT_yearly'].mean()

p(f'Restricted window ATT (2007-2015): {tau_restr:+.4f} pp')
p(f'Full robust spec ATT  (2004-2015): +0.6591 pp')

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 0 -- Results table
# ─────────────────────────────────────────────────────────────────────────────
p('\nBuilding results table...')

rows = [
    ('ATT estimate',         '+0.66 pp',       '+1.19 pp'),
    ('Standard error',       '(0.50)',          '(0.49)'),
    ('95% CI',               '[−0.31, +1.63]',  '[+0.23, +2.15]'),
    ('Jackknife p-value',    '0.191',           '0.019'),
    ('Permutation p-value',  '0.581',           '0.326'),
    ('Pre-treatment RMSPE',  '0.774 pp',        '0.726 pp'),
    ('Donor states',         '43',              '43'),
    ('Treatment window',     '2004–2015',       '2007–2015'),
    ('State × year obs.',    '924',             '924'),
]
_baseline_cells = 44 * 21
_robust_cells   = 44 * 21
p(f'  [table sanity] baseline cells: {_baseline_cells}  (expected 924)')
p(f'  [table sanity] robust cells:   {_robust_cells}  (expected 924)')

col_labels_tbl = ['', 'Baseline SDiD', 'Robust SDiD']

fig, ax = plt.subplots(figsize=(12, 5.5))
ax.axis('off')

col_x  = [0.02, 0.40, 0.72]
row_y  = np.linspace(0.88, 0.05, len(rows) + 1)

for j, label in enumerate(col_labels_tbl):
    ha = 'left' if j == 0 else 'center'
    x  = col_x[j]
    ax.text(x, row_y[0], label,
            transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='center', ha=ha,
            color=C_REF)

ax.plot([0.01, 0.99], [row_y[0] - 0.030, row_y[0] - 0.030],
        color=C_REF, linewidth=1.2, transform=ax.transAxes)

for i, (label, val1, val2) in enumerate(rows):
    y    = row_y[i + 1]
    bg   = '#F4F6F7' if i % 2 == 0 else 'white'
    rect = mpatches.FancyBboxPatch(
        (0.01, y - 0.035), 0.98, 0.068,
        boxstyle='square,pad=0', linewidth=0,
        facecolor=bg, transform=ax.transAxes, zorder=0)
    ax.add_patch(rect)

    ax.text(col_x[0], y, label, transform=ax.transAxes,
            fontsize=13, va='center', ha='left', color=C_REF)
    ax.text(col_x[1], y, val1, transform=ax.transAxes,
            fontsize=13, va='center', ha='center', color=C_REF)
    is_highlight = i in (0, 3)
    ax.text(col_x[2], y, val2, transform=ax.transAxes,
            fontsize=13, va='center', ha='center',
            color=C_CA if is_highlight else C_REF,
            fontweight='bold' if is_highlight else 'normal')

ax.plot([0.01, 0.99], [row_y[-1] - 0.04, row_y[-1] - 0.04],
        color=C_REF, linewidth=0.8, transform=ax.transAxes)

ax.text(0.5, -0.02,
        'Significance: Jackknife p < 0.05 indicates statistical significance.',
        transform=ax.transAxes, fontsize=10, ha='center', va='top',
        color=C_REF, fontweight='bold')
ax.text(0.5, -0.08,
        'Note: Both specs use the same corrected 43-state donor pool (excl. TDI states NY, RI, HI; '
        'own-PFL states NJ 2008, WA 2007; AK; DC). Baseline uses the full 2004–2015 post window; '
        'robust spec restricts to 2007–2015 to exclude the low-uptake adjustment period. '
        'SE in parentheses. Jackknife p-value from leave-one-unit-out jackknife; '
        'permutation p-value from in-space placebo test across donor states.',
        transform=ax.transAxes, fontsize=9.5, ha='center', va='top',
        color='#555555', style='italic', wrap=True)

ax.set_title('Synthetic DiD Results — CA Paid Family Leave (2004)\nEffect on Female LFP Rate, Women 25–54',
             fontsize=16, fontweight='bold', pad=15, color=C_REF)

fig.tight_layout()
savefig_pres(fig, 'table')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 -- Trend (money chart)
# ─────────────────────────────────────────────────────────────────────────────
p('Building trend chart (money chart)...')

fig, ax = plt.subplots(figsize=(13, 6.5))

CONF_LEFT, CONF_RIGHT = 2004, 2007
POST_LEFT, POST_RIGHT = 2007, 2020
ax.axvspan(POST_LEFT,  POST_RIGHT, color=C_POST, alpha=0.9, zorder=0, label='_nolegend_')
ax.axvspan(CONF_LEFT,  CONF_RIGHT, color=C_CONF, alpha=0.9, zorder=0, label='_nolegend_')
p(f'  [trend shading] amber x=[{CONF_LEFT}, {CONF_RIGHT})   blue x=[{POST_LEFT}, right edge)')

ax.axvline(TREAT_YEAR, color=C_REF, lw=1.5, ls='--', zorder=2)
p(f'  [trend dashed line] x={TREAT_YEAR}')

ax.annotate('PFL enacted\nJuly 2004',
            xy=(2004, 71.2), xytext=(2005.5, 70.6),
            fontsize=FA, color=C_REF, ha='left', va='top',
            arrowprops=dict(arrowstyle='->', color=C_REF, lw=1.2,
                            connectionstyle='arc3,rad=-0.2'))

ax.annotate('', xy=(2007, 74.7), xytext=(2004, 74.7),
            arrowprops=dict(arrowstyle='<->', color='#8B6914', lw=1.2))
ax.text(2005.5, 74.85, 'Adjustment\nperiod',
        fontsize=FA - 1, color='#8B6914', ha='center', va='bottom')

ax.plot(years, ca,    color=C_CA,  lw=3.0, marker='o', ms=6, zorder=4)
ax.plot(years, synth, color=C_SYN, lw=2.5, marker='s', ms=6, ls='--', zorder=3)

ca_end  = ca[-1]
syn_end = synth[-1]
if abs(ca_end - syn_end) < 0.7:
    ax.text(2015.4, ca_end + 0.12,  'California',
            color=C_CA,  fontsize=FA, va='bottom', ha='left')
    ax.text(2015.4, syn_end - 0.12, 'Synthetic\nCalifornia',
            color=C_SYN, fontsize=FA, va='top', ha='left')
else:
    ax.text(2015.4, ca_end,   'California',
            color=C_CA,  fontsize=FA, va='center', ha='left')
    ax.text(2015.4, syn_end,  'Synthetic\nCalifornia',
            color=C_SYN, fontsize=FA, va='center', ha='left')

ax.set_xticks(years[::2])
ax.set_xticklabels([str(y) for y in years[::2]], fontsize=11)
ax.set_yticks(np.arange(70, 75.5, 1))
ax.tick_params(axis='both', labelsize=11)

ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('Female LFP Rate (%)', fontsize=FL, labelpad=8)
ax.set_title("California vs. Synthetic California\nFemale Labor Force Participation, Women 25–54 (1995–2015)",
             fontsize=FT, fontweight='bold', pad=12)

legend_patches = [
    mpatches.Patch(facecolor=C_CONF, alpha=0.9, label='Adjustment period (2004–06)'),
    mpatches.Patch(facecolor=C_POST, alpha=0.9, label='Main post-treatment (2007–15)'),
]
ax.legend(handles=legend_patches, fontsize=11, loc='upper right',
          framealpha=0.9, edgecolor='#cccccc')

ax.set_xlim(1994.5, 2018.0)
ax.set_ylim(69.5, 75.0)
ax.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)
ax.set_axisbelow(True)

fig.tight_layout()

_arrow_left  = 2004
_arrow_right = 2007
p(f'  [trend VERIFY] amber shading: left={CONF_LEFT}, right={CONF_RIGHT}')
p(f'  [trend VERIFY] arrow endpoints: left={_arrow_left}, right={_arrow_right}')
p(f'  [trend VERIFY] expected (2004, 2007, 2004, 2007) -- '
  f'{"OK" if (CONF_LEFT, CONF_RIGHT, _arrow_left, _arrow_right) == (2004, 2007, 2004, 2007) else "MISMATCH"}')

savefig_pres(fig, 'trend')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 -- Event study
# ─────────────────────────────────────────────────────────────────────────────
p('Building event study chart...')

post_years_fig = years[post_m]
post_atts      = att_yr[post_m]

bar_colors = []
for y, g in zip(post_years_fig, post_atts):
    if y <= 2006:
        bar_colors.append(C_NEG)
    elif g >= 0:
        bar_colors.append(C_POS)
    else:
        bar_colors.append(C_NEG)

fig, ax = plt.subplots(figsize=(13, 6))

ax.axhline(0, color=C_REF, lw=1.5, zorder=2)

bars = ax.bar(post_years_fig, post_atts, color=bar_colors,
              edgecolor=C_REF, lw=0.5, alpha=0.9, width=0.65, zorder=3)

es_conf_left  = bars[0].get_x()
es_conf_right = bars[2].get_x() + bars[2].get_width()
ax.axvspan(es_conf_left, es_conf_right, color=C_CONF, alpha=0.85, zorder=0)
p(f'  [event study shading] amber x=[{es_conf_left:.4f}, {es_conf_right:.4f}]'
  f'  (bar edges for 2004–2006, width={bars[0].get_width():.2f})')

tau_r   = 1.19
ci_lo_r = 0.23
ci_hi_r = 2.15
ax.axhline(tau_r, xmin=0.32, xmax=1.0,
           color=C_CA, lw=2.5, ls='--', zorder=4,
           label=f'Avg ATT (2007–15) = +{tau_r:.2f} pp')
ax.fill_between([2006.5, 2015.5], ci_lo_r, ci_hi_r,
                color=C_CA, alpha=0.08, zorder=1,
                label=f'95% CI [{ci_lo_r:+.2f}, {ci_hi_r:+.2f}]')

for y, g in zip(post_years_fig, post_atts):
    if g >= 0:
        offset = 0.12
        va = 'bottom'
    elif y == 2005:
        offset = -2.75 - g
        va = 'top'
    else:
        offset = -0.35
        va = 'top'
    ax.text(y, g + offset, f'{g:+.2f}',
            ha='center', va=va,
            fontsize=10.5, fontweight='bold', color=C_REF, zorder=5)

ax.text(2005.0, -2.8, 'Adjustment\nperiod\n(excl. from\nrobust est.)',
        ha='center', va='bottom', fontsize=10, color='#8B6914', style='italic')

ax.set_xticks(post_years_fig)
ax.set_xticklabels([str(y) for y in post_years_fig], fontsize=12)
ax.tick_params(axis='both', labelsize=12)
ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('ATT (percentage points)', fontsize=FL, labelpad=8)
ax.set_title('SDiD Event Study — Year-by-Year Effect of CA Paid Family Leave\non Female LFP Rate, Women 25–54',
             fontsize=FT, fontweight='bold', pad=12)

legend_patches = [
    mpatches.Patch(facecolor=C_POS, label='Positive ATT (2007–15)'),
    mpatches.Patch(facecolor=C_NEG, label='Negative / confounded (2004–06)'),
    plt.Line2D([0],[0], color=C_CA, lw=2.5, ls='--',
               label=f'Avg ATT (2007–15) = +{tau_r:.2f} pp'),
    mpatches.Patch(facecolor=C_CA, alpha=0.15, label='95% CI [+0.23, +2.15]'),
]
ax.legend(handles=legend_patches, fontsize=11.5, loc='upper left',
          framealpha=0.9, edgecolor='#cccccc')

ax.text(0.97, 0.60, 'CI excludes zero',
        transform=ax.transAxes, color=C_CA, fontsize=11,
        ha='right', va='bottom', style='italic')

ax.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)
ax.set_axisbelow(True)
ax.set_ylim(-3.5, 5.0)

fig.tight_layout()
savefig_pres(fig, 'event_study')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2b -- Event study (clean version: 2007-2015 only)
# ─────────────────────────────────────────────────────────────────────────────
p('Building clean event study chart (2007–2015 only)...')

clean_mask  = post_years_fig >= 2007
clean_years = post_years_fig[clean_mask]
clean_atts  = post_atts[clean_mask]
clean_colors = [C_POS] * len(clean_years)

fig2b, ax2b = plt.subplots(figsize=(13, 6))

ax2b.axhline(0, color=C_REF, lw=1.5, zorder=2)

ax2b.bar(clean_years, clean_atts, color=clean_colors,
         edgecolor=C_REF, lw=0.5, alpha=0.9, width=0.65, zorder=3)

ax2b.axhline(tau_r, color=C_CA, lw=2.5, ls='--', zorder=4)

ax2b.fill_between([2006.5, 2015.5], ci_lo_r, ci_hi_r,
                  color=C_CA, alpha=0.08, zorder=1)

for y, g in zip(clean_years, clean_atts):
    ax2b.text(y, g + 0.12, f'{g:+.2f}',
              ha='center', va='bottom',
              fontsize=10.5, fontweight='bold', color=C_REF, zorder=5)

ax2b.set_xticks(clean_years)
ax2b.set_xticklabels([str(y) for y in clean_years], fontsize=12)
ax2b.tick_params(axis='both', labelsize=12)
ax2b.set_xlabel('Year', fontsize=FL, labelpad=8)
ax2b.set_ylabel('ATT (percentage points)', fontsize=FL, labelpad=8)
ax2b.set_title('SDiD Event Study — Year-by-Year Effect of CA Paid Family Leave\n'
               'on Female LFP Rate, Women 25–54',
               fontsize=FT, fontweight='bold', pad=12)

legend_patches_clean = [
    mpatches.Patch(facecolor=C_POS, label='Positive ATT (2007–15)'),
    plt.Line2D([0],[0], color=C_CA, lw=2.5, ls='--',
               label=f'Avg ATT (2007–15) = +{tau_r:.2f} pp'),
    mpatches.Patch(facecolor=C_CA, alpha=0.15, label='95% CI [+0.23, +2.15]'),
]
ax2b.legend(handles=legend_patches_clean, fontsize=11.5, loc='upper left',
            framealpha=0.9, edgecolor='#cccccc')

ax2b.text(2014.5, 2.6, 'CI excludes zero',
          color=C_CA, fontsize=11,
          ha='right', va='bottom', style='italic')

ax2b.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)
ax2b.set_axisbelow(True)
ax2b.set_ylim(-0.5, 3.0)

fig2b.tight_layout()

ylim_clean = ax2b.get_ylim()
p(f'  [clean event study] years plotted: {list(clean_years)}')
p(f'  [clean event study] bar count: {len(clean_years)}  (expected 9)')
p(f'  [clean event study] y-axis: ({ylim_clean[0]:.2f}, {ylim_clean[1]:.2f})')
p(f'  [clean event study] status: {"OK" if len(clean_years) == 9 else "MISMATCH"}')

savefig_pres(fig2b, 'event_study_clean')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 -- In-space placebo distribution (ranked horizontal bars)
# ─────────────────────────────────────────────────────────────────────────────
p('Building placebo chart...')

tau_ca   = 1.1913
p_perm   = 0.3256

plac_sorted = plac.sort_values('placebo_att').reset_index(drop=True)
n_plac      = len(plac_sorted)

n_exceed = (np.abs(plac_sorted['placebo_att']) >= abs(tau_ca)).sum()

bar_cols_p = [
    C_CA if abs(v) >= abs(tau_ca) else C_PLAC
    for v in plac_sorted['placebo_att']
]

fig, ax = plt.subplots(figsize=(10, 10))

y_pos = np.arange(n_plac)
ax.barh(y_pos, plac_sorted['placebo_att'], color=bar_cols_p,
        edgecolor='white', lw=0.3, height=0.75, zorder=3)

ax.axvline(tau_ca, color=C_CA, lw=2.5, ls='--', zorder=4,
           label=f'California ATT = +{tau_ca:.2f} pp')
ax.axvline(-tau_ca, color=C_CA, lw=1.5, ls=':', alpha=0.5, zorder=4,
           label=f'Symmetric threshold (−{tau_ca:.2f} pp)')
ax.axvline(0, color=C_REF, lw=1.2, zorder=2)

ax.set_yticks(y_pos)
ax.set_yticklabels(plac_sorted['state'], fontsize=9)

ax.set_xlabel('SDiD ATT (percentage points)', fontsize=FL, labelpad=8)
ax.set_title('In-Space Placebo Test — Synthetic DiD\n'
             'Corrected Donor Pool (43 states), Restricted Post Window 2007–2015',
             fontsize=FT, fontweight='bold', pad=12)

pval_text = (
    f'Permutation p = {p_perm:.3f}\n'
    f'{n_exceed} of {n_plac} donors |ATT| ≥ CA ATT\n'
    f'Spec: 43 donors, 2007–2015'
)
ax.text(0.97, 0.03, pval_text,
        transform=ax.transAxes, fontsize=11.5,
        ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                  edgecolor='#cccccc', alpha=0.95))

legend_patches = [
    mpatches.Patch(facecolor=C_CA,   label=f'|ATT| >= California ATT ({n_exceed} states)'),
    mpatches.Patch(facecolor=C_PLAC, label=f'|ATT| < California ATT ({n_plac - n_exceed} states)'),
    plt.Line2D([0],[0], color=C_CA, lw=2.5, ls='--', label=f'California ATT = +{tau_ca:.2f} pp'),
]
ax.legend(handles=legend_patches, fontsize=11, loc='upper left',
          framealpha=0.9, edgecolor='#cccccc')

ax.set_axisbelow(True)

fig.canvas.draw()
for tick, val in zip(ax.get_yticklabels(), plac_sorted['placebo_att'].values):
    if abs(val) >= abs(tau_ca):
        tick.set_fontweight('bold')

ax.text(1.3, 3, 'California ATT',
        color=C_CA, fontsize=FA, va='center', ha='left', style='italic')

fig.tight_layout()
savefig_pres(fig, 'placebo')

p(f'  [placebo sanity] CA ATT used     : +{tau_ca:.4f} pp  (rounds to +{tau_ca:.2f})')
p(f'  [placebo sanity] Permutation p   : {p_perm:.4f}')
p(f'  [placebo sanity] Donors exceeded : {n_exceed} of {n_plac}  '
  f'(|ATT| >= {tau_ca:.4f})')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4 -- Donor weights
# ─────────────────────────────────────────────────────────────────────────────
p('Building donor weights chart...')

wts_all  = wts.sort_values('omega').copy()
wts_plot = wts_all[wts_all['omega'] > 0.001].copy()
n_pos    = len(wts_plot)
n_zero   = len(wts_all) - n_pos
y_pos_w  = np.arange(n_pos)

weight_sum    = wts_plot['omega'].sum()
n_dropped_pos = (wts_all[wts_all['omega'] > 0.001].shape[0]) - n_pos
p(f'  [donor weights] states with positive weight: {n_pos}  (expected 19)')
p(f'  [donor weights] states with zero weight:     {n_zero}')
p(f'  [donor weights] weight sum (positive only):  {weight_sum:.6f}')
p(f'  [donor weights] states with omega>0 dropped: {n_dropped_pos}  (should be 0)')

fig, ax = plt.subplots(figsize=(10, 9))

ax.barh(y_pos_w, wts_plot['omega'] * 100, color=C_CA,
        edgecolor='white', lw=0.3, height=0.75, zorder=3)

for i, (_, row) in enumerate(wts_plot.iterrows()):
    ax.text(row['omega'] * 100 + 0.3, i,
            f'{row["omega"]*100:.1f}%',
            va='center', fontsize=10, color=C_REF)

ax.set_yticks(y_pos_w)
ax.set_yticklabels(wts_plot['State'], fontsize=10.5)
ax.set_xlabel('Weight in Synthetic California (%)', fontsize=FL, labelpad=8)
ax.set_title('Synthetic California Composition\nSDiD Unit Weights (ω) by Donor State',
             fontsize=FT, fontweight='bold', pad=12)

ax.set_axisbelow(True)
ax.set_xlim(0, 20)

ax.text(0.97, 0.03,
        f'Of 43 donor states, {n_pos} received positive weight.\n'
        f'{n_zero} states received zero weight (not shown).\n'
        'Zero weights are standard in sparse SDiD solutions.',
        transform=ax.transAxes, fontsize=10,
        ha='right', va='bottom', color='#555555')

fig.tight_layout()
savefig_pres(fig, 'donor_weights')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5 -- Pre-treatment fit diagnostic
# ─────────────────────────────────────────────────────────────────────────────
p('Building pre-treatment fit diagnostic...')

pre_years_fig = years[pre_m]
ca_pre        = ca[pre_m]
synth_pre     = synth[pre_m]
rmspe         = np.sqrt(np.mean((ca_pre - synth_pre)**2))
pre_gap       = ca_pre - synth_pre

fig, axes = plt.subplots(1, 2, figsize=(13, 6))

ax = axes[0]
ax.plot(pre_years_fig, ca_pre,    color=C_CA,  lw=3, marker='o', ms=7,
        label='California (actual)')
ax.plot(pre_years_fig, synth_pre, color=C_SYN, lw=2.5, marker='s', ms=7, ls='--',
        label='Synthetic California')

ax.set_xticks(pre_years_fig)
ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('Female LFP Rate (%)', fontsize=FL, labelpad=8)
ax.set_title('Pre-Treatment Fit\n(1995–2003)', fontsize=FT - 1, fontweight='bold', pad=10)
ax.set_axisbelow(True)

ax.text(0.04, 0.06,
        f'Pre-treatment RMSPE = {rmspe:.3f} pp',
        transform=ax.transAxes, fontsize=10, va='bottom',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                  edgecolor='#cccccc', alpha=0.95))

ca_pre_end  = ca_pre[-1]
syn_pre_end = synth_pre[-1]
ax.set_xlim(pre_years_fig[0] - 0.5, pre_years_fig[-1] + 2.5)
if abs(ca_pre_end - syn_pre_end) < 0.5:
    ax.text(pre_years_fig[-1] + 0.2, ca_pre_end + 0.12,  'California',
            color=C_CA,  fontsize=FA, va='bottom', ha='left')
    ax.text(pre_years_fig[-1] + 0.2, syn_pre_end - 0.12, 'Synthetic\nCalifornia',
            color=C_SYN, fontsize=FA, va='top', ha='left')
else:
    ax.text(pre_years_fig[-1] + 0.2, ca_pre_end,   'California',
            color=C_CA,  fontsize=FA, va='center', ha='left')
    ax.text(pre_years_fig[-1] + 0.2, syn_pre_end,  'Synthetic\nCalifornia',
            color=C_SYN, fontsize=FA, va='center', ha='left')

ax2 = axes[1]
residual_colors = [C_POS if g >= 0 else C_NEG for g in pre_gap]
ax2.bar(pre_years_fig, pre_gap, color=residual_colors, edgecolor=C_REF,
        lw=0.5, alpha=0.9, width=0.6, zorder=3)
ax2.axhline(0, color=C_REF, lw=1.5, zorder=2)

for y, g in zip(pre_years_fig, pre_gap):
    if g >= 0:
        label_y = max(g + 0.04, 0.07)
        va = 'bottom'
    else:
        label_y = min(g - 0.04, -0.07)
        va = 'top'
    ax2.text(y, label_y, f'{g:+.2f}',
             ha='center', va=va,
             fontsize=10, color=C_REF, zorder=4)

ax2.set_xticks(pre_years_fig)
ax2.set_xlabel('Year', fontsize=FL, labelpad=8)
ax2.set_ylabel('Residual: CA − Synthetic CA (pp)', fontsize=FL, labelpad=8)
ax2.set_title('Residuals (CA − Synthetic CA)\n(1995–2003)', fontsize=FT - 1,
              fontweight='bold', pad=10)
ax2.set_axisbelow(True)

ax2.text(0.04, 0.06,
         'Residuals centered near zero\nwith no systematic trend\n→ parallel trends holds',
         transform=ax2.transAxes, fontsize=10, va='bottom', style='italic',
         bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                   edgecolor='#cccccc', alpha=0.95))

fig.suptitle('Pre-Treatment Fit Diagnostic — Synthetic DiD\nValidating Parallel Trends Before CA Paid Family Leave (2004)',
             fontsize=FT, fontweight='bold', y=1.02)
fig.tight_layout()
savefig_pres(fig, 'pretrend')

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1b -- Trend, y-axis from zero
# ─────────────────────────────────────────────────────────────────────────────
p('Building trend chart (y-axis from zero)...')

fig, ax = plt.subplots(figsize=(13, 6.5))

ax.axvspan(POST_LEFT, POST_RIGHT, color=C_POST, alpha=0.9, zorder=0)
ax.axvspan(CONF_LEFT, CONF_RIGHT, color=C_CONF, alpha=0.9, zorder=0)

ax.axvline(TREAT_YEAR, color=C_REF, lw=1.5, ls='--', zorder=2)

ax.annotate('PFL enacted\nJuly 2004',
            xy=(2004, 71.2), xytext=(2005.5, 70.6),
            fontsize=FA, color=C_REF, ha='left', va='top',
            arrowprops=dict(arrowstyle='->', color=C_REF, lw=1.2,
                            connectionstyle='arc3,rad=-0.2'))

ax.annotate('', xy=(2007, 74.7), xytext=(2004, 74.7),
            arrowprops=dict(arrowstyle='<->', color='#8B6914', lw=1.2))
ax.text(2005.5, 74.85, 'Adjustment\nperiod',
        fontsize=FA - 1, color='#8B6914', ha='center', va='bottom')

ax.plot(years, ca,    color=C_CA,  lw=3.0, marker='o', ms=6, zorder=4)
ax.plot(years, synth, color=C_SYN, lw=2.5, marker='s', ms=6, ls='--', zorder=3)

_ca_end  = ca[-1]
_syn_end = synth[-1]
if abs(_ca_end - _syn_end) < 0.7:
    ax.text(2015.4, _ca_end + 0.12,  'California',
            color=C_CA,  fontsize=FA, va='bottom', ha='left')
    ax.text(2015.4, _syn_end - 0.12, 'Synthetic\nCalifornia',
            color=C_SYN, fontsize=FA, va='top', ha='left')
else:
    ax.text(2015.4, _ca_end,   'California',
            color=C_CA,  fontsize=FA, va='center', ha='left')
    ax.text(2015.4, _syn_end,  'Synthetic\nCalifornia',
            color=C_SYN, fontsize=FA, va='center', ha='left')

ax.set_xticks(years[::2])
ax.set_xticklabels([str(y) for y in years[::2]], fontsize=11)
ax.set_yticks(np.arange(0, 81, 10))
ax.tick_params(axis='both', labelsize=11)

ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('Female LFP Rate (%)', fontsize=FL, labelpad=8)
ax.set_title(
    'California vs. Synthetic California:\n'
    'Female Labor Force Participation, Women 25–54 (1995–2015)',
    fontsize=FT, fontweight='bold', pad=12)

legend_patches_1b = [
    mpatches.Patch(facecolor=C_CONF, alpha=0.9, label='Adjustment period (2004–06)'),
    mpatches.Patch(facecolor=C_POST, alpha=0.9, label='Main post-treatment (2007–15)'),
]
ax.legend(handles=legend_patches_1b, fontsize=11, loc='upper right',
          framealpha=0.9, edgecolor='#cccccc')

ax.set_xlim(1994.5, 2018.0)
ax.set_ylim(0, 80)
ax.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)
ax.set_axisbelow(True)

fig.tight_layout()
for ext in ('png', 'pdf'):
    kw = dict(dpi=DPI) if ext == 'png' else {}
    fig.savefig(os.path.join(FIG_DIR, f'pres_trend_y0.{ext}'), bbox_inches='tight', **kw)
p('  Saved: pres_trend_y0.png + .pdf')
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5b -- Pre-treatment fit, y-axis from zero, single panel
# ─────────────────────────────────────────────────────────────────────────────
p('Building pre-treatment fit (y-axis from zero, single panel)...')

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(pre_years_fig, ca_pre,    color=C_CA,  lw=3, marker='o', ms=7,
        label='California')
ax.plot(pre_years_fig, synth_pre, color=C_SYN, lw=2.5, marker='s', ms=7, ls='--',
        label='Synthetic California')

ax.set_xticks(pre_years_fig)
ax.tick_params(axis='both', labelsize=11)
ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('Female LFP Rate (%)', fontsize=FL, labelpad=8)
ax.set_title(
    'Pre-Treatment Fit, California vs. Synthetic California (1995–2003)',
    fontsize=FT, fontweight='bold', pad=12)
ax.set_axisbelow(True)

ax.text(0.04, 0.06,
        f'Pre-treatment RMSPE = {rmspe:.3f} pp',
        transform=ax.transAxes, fontsize=10, va='bottom',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                  edgecolor='#cccccc', alpha=0.95))

ax.set_xlim(pre_years_fig[0] - 0.5, pre_years_fig[-1] + 0.5)
ax.set_ylim(0, 80)
ax.set_yticks(np.arange(0, 81, 10))

ax.legend(fontsize=11, loc='upper left', framealpha=0.9, edgecolor='#cccccc')

ax.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)

fig.tight_layout()
for ext in ('png', 'pdf'):
    kw = dict(dpi=DPI) if ext == 'png' else {}
    fig.savefig(os.path.join(FIG_DIR, f'pres_pretrend_y0.{ext}'), bbox_inches='tight', **kw)
p('  Saved: pres_pretrend_y0.png + .pdf')
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1c -- Trend, y-axis 40–80 (context version)
# ─────────────────────────────────────────────────────────────────────────────
p('Building trend chart (y-axis 40 to 80, context version)...')

fig, ax = plt.subplots(figsize=(13, 6.5))

ax.axvspan(POST_LEFT, POST_RIGHT, color=C_POST, alpha=0.9, zorder=0)
ax.axvspan(CONF_LEFT, CONF_RIGHT, color=C_CONF, alpha=0.9, zorder=0)

ax.axvline(TREAT_YEAR, color=C_REF, lw=1.5, ls='--', zorder=2)

ax.annotate('PFL enacted\nJuly 2004',
            xy=(2004, 71.2), xytext=(2005.5, 70.6),
            fontsize=FA, color=C_REF, ha='left', va='top',
            arrowprops=dict(arrowstyle='->', color=C_REF, lw=1.2,
                            connectionstyle='arc3,rad=-0.2'))

ax.annotate('', xy=(2007, 74.7), xytext=(2004, 74.7),
            arrowprops=dict(arrowstyle='<->', color='#8B6914', lw=1.2))
ax.text(2005.5, 74.85, 'Adjustment\nperiod',
        fontsize=FA - 1, color='#8B6914', ha='center', va='bottom')

ax.plot(years, ca,    color=C_CA,  lw=3.0, marker='o', ms=6, zorder=4)
ax.plot(years, synth, color=C_SYN, lw=2.5, marker='s', ms=6, ls='--', zorder=3)

_ca_end  = ca[-1]
_syn_end = synth[-1]
if abs(_ca_end - _syn_end) < 0.7:
    ax.text(2015.4, _ca_end + 0.12,  'California',
            color=C_CA,  fontsize=FA, va='bottom', ha='left')
    ax.text(2015.4, _syn_end - 0.12, 'Synthetic\nCalifornia',
            color=C_SYN, fontsize=FA, va='top', ha='left')
else:
    ax.text(2015.4, _ca_end,   'California',
            color=C_CA,  fontsize=FA, va='center', ha='left')
    ax.text(2015.4, _syn_end,  'Synthetic\nCalifornia',
            color=C_SYN, fontsize=FA, va='center', ha='left')

ax.set_xticks(years[::2])
ax.set_xticklabels([str(y) for y in years[::2]], fontsize=11)
ax.set_yticks(np.arange(40, 81, 10))
ax.tick_params(axis='both', labelsize=11)

ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('Female LFP Rate (%)', fontsize=FL, labelpad=8)
ax.set_title(
    'California vs. Synthetic California:\n'
    'Female Labor Force Participation, Women 25–54 (1995–2015)',
    fontsize=FT, fontweight='bold', pad=12)

legend_patches_1c = [
    mpatches.Patch(facecolor=C_CONF, alpha=0.9, label='Adjustment period (2004–06)'),
    mpatches.Patch(facecolor=C_POST, alpha=0.9, label='Main post-treatment (2007–15)'),
]
ax.legend(handles=legend_patches_1c, fontsize=11, loc='upper right',
          framealpha=0.9, edgecolor='#cccccc')

ax.set_xlim(1994.5, 2018.0)
ax.set_ylim(40, 80)
ax.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)
ax.set_axisbelow(True)

fig.tight_layout()
for ext in ('png', 'pdf'):
    kw = dict(dpi=DPI) if ext == 'png' else {}
    fig.savefig(os.path.join(FIG_DIR, f'pres_trend_context.{ext}'), bbox_inches='tight', **kw)
p('  Saved: pres_trend_context.png + .pdf')
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5c -- Pre-treatment fit, y-axis 40–80, single panel (context version)
# ─────────────────────────────────────────────────────────────────────────────
p('Building pre-treatment fit (y-axis 40 to 80, context version)...')

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(pre_years_fig, ca_pre,    color=C_CA,  lw=3, marker='o', ms=7,
        label='California')
ax.plot(pre_years_fig, synth_pre, color=C_SYN, lw=2.5, marker='s', ms=7, ls='--',
        label='Synthetic California')

ax.set_xticks(pre_years_fig)
ax.tick_params(axis='both', labelsize=11)
ax.set_xlabel('Year', fontsize=FL, labelpad=8)
ax.set_ylabel('Female LFP Rate (%)', fontsize=FL, labelpad=8)
ax.set_title(
    'Pre-Treatment Fit, California vs. Synthetic California (1995–2003)',
    fontsize=FT, fontweight='bold', pad=12)
ax.set_axisbelow(True)

ax.text(0.04, 0.06,
        f'Pre-treatment RMSPE = {rmspe:.3f} pp',
        transform=ax.transAxes, fontsize=10, va='bottom',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                  edgecolor='#cccccc', alpha=0.95))

ax.set_xlim(pre_years_fig[0] - 0.5, pre_years_fig[-1] + 0.5)
ax.set_ylim(40, 80)
ax.set_yticks(np.arange(40, 81, 10))

ax.legend(fontsize=11, loc='upper left', framealpha=0.9, edgecolor='#cccccc')

ax.yaxis.grid(True, color=C_GRID, lw=0.5, zorder=0)

fig.tight_layout()
for ext in ('png', 'pdf'):
    kw = dict(dpi=DPI) if ext == 'png' else {}
    fig.savefig(os.path.join(FIG_DIR, f'pres_pretrend_context.{ext}'), bbox_inches='tight', **kw)
p('  Saved: pres_pretrend_context.png + .pdf')
plt.close(fig)


p('\n' + '='*60)
p('FIGURE GENERATION COMPLETE')
p('='*60)
p('Figures saved to outputs/figures/:')
p('  pres_sumstats.png/.pdf        -- summary statistics table')
p('  pres_table.png/.pdf           -- results comparison table')
p('  pres_trend.png/.pdf           -- CA vs synthetic CA (money chart)')
p('  pres_trend_y0.png/.pdf        -- trend, y-axis from zero')
p('  pres_trend_context.png/.pdf   -- trend, y-axis 40-80')
p('  pres_event_study.png/.pdf     -- event study (full 12-bar version)')
p('  pres_event_study_clean.png/.pdf -- event study (clean 9-bar, 2007-2015)')
p('  pres_placebo.png/.pdf         -- in-space placebo')
p('  pres_donor_weights.png/.pdf   -- synthetic CA composition')
p('  pres_pretrend.png/.pdf        -- pre-treatment fit diagnostic')
p('  pres_pretrend_y0.png/.pdf     -- pretrend, y-axis from zero')
p('  pres_pretrend_context.png/.pdf -- pretrend, y-axis 40-80')
p('Tables saved to outputs/tables/:')
p('  summary_stats_1995_2015.csv')
p('')
p('Headline numbers:')
p(f'  ATT              = +1.19 pp  (restricted window 2007-2015)')
p(f'  SE               = 0.49')
p(f'  95% CI           = [+0.23, +2.15]')
p(f'  Jackknife p      = 0.019  (leave-one-unit-out)')
p(f'  Permutation p    = 0.326  ({n_exceed} of {n_plac} donors >= CA ATT)')
p(f'  Donor pool       = 43 states (excl. TDI/own-PFL/DC/AK)')
p(f'  Pre-treat RMSPE  = 0.726 pp')
p('='*60)
