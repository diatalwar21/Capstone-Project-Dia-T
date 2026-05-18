"""
CPS Deduplication + UHRSWORKT Quarantine
=========================================
Corrects two issues from the first cleaning pass:

1. March rows (ASECFLAG=2) were incorrectly dropped in the original script.
   ASECFLAG=2 means "Basic Monthly record" — not an ASEC supplement record.
   This pass keeps them and instead drops ASECFLAG==1 (pure ASEC supplement),
   of which there are none in this extract.

2. Deduplicates on CPSID + PERNUM + YEAR + MONTH using the raw file
   (CPSID/PERNUM were not saved in the first individual file).

3. Drops mean_uhrs_worked from the state-year panel — implausible values,
   flagged pending codebook verification.

Outputs (new files, originals untouched)
-----------------------------------------
cps_deduped_individual_1995_2015.csv
state_year_panel_deduped_1995_2015.csv
"""

import os
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def p(*args): print(*args, flush=True)

HERE = os.path.dirname(os.path.abspath(__file__))

STATE_MAP = {
    1:"Alabama", 2:"Alaska", 4:"Arizona", 5:"Arkansas", 6:"California",
    8:"Colorado", 9:"Connecticut", 10:"Delaware", 11:"District of Columbia",
    12:"Florida", 13:"Georgia", 15:"Hawaii", 16:"Idaho", 17:"Illinois",
    18:"Indiana", 19:"Iowa", 20:"Kansas", 21:"Kentucky", 22:"Louisiana",
    23:"Maine", 24:"Maryland", 25:"Massachusetts", 26:"Michigan", 27:"Minnesota",
    28:"Mississippi", 29:"Missouri", 30:"Montana", 31:"Nebraska", 32:"Nevada",
    33:"New Hampshire", 34:"New Jersey", 35:"New Mexico", 36:"New York",
    37:"North Carolina", 38:"North Dakota", 39:"Ohio", 40:"Oklahoma",
    41:"Oregon", 42:"Pennsylvania", 44:"Rhode Island", 45:"South Carolina",
    46:"South Dakota", 47:"Tennessee", 48:"Texas", 49:"Utah", 50:"Vermont",
    51:"Virginia", 53:"Washington", 54:"West Virginia", 55:"Wisconsin", 56:"Wyoming"
}

# ── Load raw data ─────────────────────────────────────────────────────────────
p("=" * 65)
p("  STEP 1 — LOAD RAW DATA")
p("=" * 65)

RAW = os.path.join(ROOT, "data/raw/ipums_cps_1995_2015.csv")
df = pd.read_csv(RAW)
n_raw = len(df)
p(f"Raw rows loaded: {n_raw:,}")

# ── ASEC identification ───────────────────────────────────────────────────────
p("\n" + "=" * 65)
p("  STEP 2 — ASEC IDENTIFICATION")
p("=" * 65)

p(f"\nASECFLAG distribution:")
p(df["ASECFLAG"].value_counts(dropna=False).to_string())

# ASECFLAG=1 → pure ASEC supplement records (drop these)
n_asec1 = (df["ASECFLAG"] == 1).sum()
p(f"\nASECFLAG==1 (pure ASEC supplement): {n_asec1:,}")
if n_asec1 > 0:
    df = df[df["ASECFLAG"] != 1].copy()
    p(f"  → Dropped. Rows remaining: {len(df):,}")
else:
    p("  → None found. No rows dropped.")

# ASECFLAG=2 → Basic Monthly records (keep — these are valid March observations)
n_asec2 = (df["ASECFLAG"] == 2).sum()
p(f"\nASECFLAG==2 (Basic Monthly March records): {n_asec2:,}")
p("  → Keeping. These are valid Basic Monthly observations.")
p("  NOTE: The prior pipeline incorrectly dropped these rows.")
p(f"  This restores {n_asec2:,} March observations (~{100*n_asec2/n_raw:.1f}% of raw).")

# Cross-check: ASECFLAG=2 should match MONTH==3 exactly
n_march = (df["MONTH"] == 3).sum()
p(f"\nSanity check — MONTH==3 count: {n_march:,} | ASECFLAG==2 count: {n_asec2:,}")
if n_march == n_asec2:
    p("  ✓ Perfect match — ASECFLAG=2 is exactly the March sample.")
else:
    p("  [FLAG] Counts differ — investigate before proceeding.")

# March sample size vs other months
p("\nObs per month (raw, should be roughly equal):")
p(df.groupby("MONTH").size().to_string())

# ── Deduplication ─────────────────────────────────────────────────────────────
p("\n" + "=" * 65)
p("  STEP 3 — DEDUPLICATION")
p("=" * 65)

n_before_dedup = len(df)
p(f"Rows before deduplication: {n_before_dedup:,}")

# Primary key: CPSID + PERNUM + YEAR + MONTH
# CPSID=0 means the person could not be linked across months — these are
# deduplicated by SERIAL + PERNUM + YEAR + MONTH as a fallback.
has_cpsid = (df["CPSID"] != 0).sum()
no_cpsid  = (df["CPSID"] == 0).sum()
p(f"\nCPSID != 0 (linkable):  {has_cpsid:,} ({100*has_cpsid/n_before_dedup:.1f}%)")
p(f"CPSID == 0 (unlinkable): {no_cpsid:,}  ({100*no_cpsid/n_before_dedup:.1f}%)")

# Deduplicate: keep first occurrence
# For CPSID=0 rows, SERIAL+PERNUM+YEAR+MONTH is the fallback key
df["_dedup_key"] = np.where(
    df["CPSID"] != 0,
    df["CPSID"].astype(str) + "_" + df["PERNUM"].astype(str) + "_" +
    df["YEAR"].astype(str) + "_" + df["MONTH"].astype(str),
    "serial_" + df["SERIAL"].astype(str) + "_" + df["PERNUM"].astype(str) + "_" +
    df["YEAR"].astype(str) + "_" + df["MONTH"].astype(str),
)

n_dupes = df.duplicated("_dedup_key").sum()
df = df[~df.duplicated("_dedup_key")].copy()
df.drop(columns=["_dedup_key"], inplace=True)

n_after_dedup = len(df)
p(f"\nDuplicate rows removed   : {n_dupes:,} ({100*n_dupes/n_before_dedup:.3f}%)")
p(f"Rows after deduplication : {n_after_dedup:,}")

# Per-month check after dedup (should still be roughly equal)
p("\nObs per month after deduplication:")
p(df.groupby("MONTH").size().to_string())

# ── Sample restriction: women 25–54 ──────────────────────────────────────────
p("\n" + "=" * 65)
p("  STEP 4 — SAMPLE RESTRICTION")
p("=" * 65)

n_before_filter = len(df)
df = df[(df["SEX"] == 2) & df["AGE"].between(25, 54)].copy()
n_after_filter = len(df)

p(f"Before filter (all persons) : {n_before_filter:,}")
p(f"After  filter (women 25–54) : {n_after_filter:,}")
p(f"Kept: {100*n_after_filter/n_before_filter:.1f}%")

# ── Variable recodes (identical to prior script) ──────────────────────────────
p("\n" + "=" * 65)
p("  STEP 5 — VARIABLE RECODES")
p("=" * 65)

df["state"]    = df["STATEFIP"].map(STATE_MAP)
df["in_lf"]   = (df["LABFORCE"] == 2).astype(int)
df["college"]  = (df["EDUC"] >= 111).astype(int)
df["married"]  = df["MARST"].isin([1, 2]).astype(int)

# UHRSWORKT — quarantined. Recode 999→NA but do NOT aggregate.
df["uhrs_worked"] = df["UHRSWORKT"].where(df["UHRSWORKT"] != 999, other=np.nan)
p("UHRSWORKT → uhrs_worked: 999 recoded to NaN. Variable is QUARANTINED.")
p("  It will appear in the individual file for audit but is excluded from the panel.")

# YNGCH: 99 → NA
df["yngch_clean"] = df["YNGCH"].where(df["YNGCH"] != 99, other=np.nan)
df["has_child"]   = (df["YNGCH"] != 99).astype(int)

# CLASSWKR
df["private_sector"] = df["CLASSWKR"].isin([22, 23]).astype(int)
df["govt_sector"]    = df["CLASSWKR"].isin([25, 26, 27]).astype(int)
df["self_employed"]  = df["CLASSWKR"].isin([13, 14]).astype(int)

# RELATE
df["is_head"] = (df["RELATE"] == 101).astype(int)

# OCC/IND era label
def occ_era(year):
    if year < 2003:  return "1950-basis"
    if year < 2011:  return "2002-basis"
    return "2010-basis"
df["occ_era"] = df["YEAR"].map(occ_era)

p(f"\nin_lf    : {df['in_lf'].mean()*100:.2f}% in labor force")
p(f"college  : {df['college'].mean()*100:.2f}% college degree")
p(f"married  : {df['married'].mean()*100:.2f}% married")
p(f"has_child: {df['has_child'].mean()*100:.2f}% have at least one child")

# MARST==9 flag
n_marst9 = (df["MARST"] == 9).sum()
p(f"\n[FLAG] MARST==9 (unknown): {n_marst9:,} rows ({100*n_marst9/len(df):.3f}%) → treated as not married")

# ── Before/after summary for individual file ──────────────────────────────────
p("\n" + "=" * 65)
p("  BEFORE / AFTER SUMMARY — INDIVIDUAL FILE")
p("=" * 65)
p(f"  Original cleaned file (first pass)  : 6,620,137 rows")
p(f"  This file (women 25–54, deduped)    : {len(df):,} rows")
p(f"  Difference                          : {len(df) - 6_620_137:+,} rows")
p(f"  Reason: March rows restored (~{n_asec2//21:,}/yr × 21 yrs × 21.6% filter rate)")

# ── Save individual file ──────────────────────────────────────────────────────
ind_cols = [
    "YEAR", "MONTH", "CPSID", "PERNUM", "STATEFIP", "state",
    "AGE", "RACE", "EDUC", "LABFORCE", "in_lf",
    "MARST", "married", "college",
    "UHRSWORKT", "uhrs_worked",  # quarantined — present for audit only
    "NCHILD", "YNGCH", "yngch_clean", "has_child",
    "FAMSIZE", "RELATE", "is_head",
    "CLASSWKR", "private_sector", "govt_sector", "self_employed",
    "IND", "OCC", "occ_era", "HOURWAGE", "WTFINL"
]
ind_out = os.path.join(ROOT, "data/derived/cps_deduped_individual_1995_2015.csv")
df[ind_cols].to_csv(ind_out, index=False)
p(f"\nSaved: {ind_out}")
p(f"  {len(df):,} rows × {len(ind_cols)} columns")

# ── Aggregate to state×year panel ─────────────────────────────────────────────
p("\n" + "=" * 65)
p("  STEP 6 — AGGREGATE TO STATE×YEAR PANEL")
p("=" * 65)

def weighted_mean(x, w):
    mask = x.notna()
    if mask.sum() == 0:
        return np.nan
    return np.average(x[mask], weights=w[mask])

agg_rows = []
for (year, state), grp in df.groupby(["YEAR", "state"]):
    w = grp["WTFINL"]
    row = {
        "YEAR"           : year,
        "state"          : state,
        "flfp"           : weighted_mean(grp["in_lf"],         w),
        "college_share"  : weighted_mean(grp["college"],       w),
        "married_share"  : weighted_mean(grp["married"],       w),
        "mean_age"       : weighted_mean(grp["AGE"],           w),
        "n_observations" : len(grp),
        "sum_weights"    : w.sum(),
        "nchild_mean"    : weighted_mean(grp["NCHILD"],        w),
        "has_child_pct"  : weighted_mean(grp["has_child"],     w) * 100,
        "yngch_mean"     : weighted_mean(grp["yngch_clean"],   w),
        "famsize_mean"   : weighted_mean(grp["FAMSIZE"],       w),
        "pct_head"       : weighted_mean(grp["is_head"],       w) * 100,
        "pct_private"    : weighted_mean(grp["private_sector"],w) * 100,
        "pct_govt"       : weighted_mean(grp["govt_sector"],   w) * 100,
        "pct_selfempl"   : weighted_mean(grp["self_employed"], w) * 100,
        "occ_era"        : grp["occ_era"].iloc[0],
    }
    agg_rows.append(row)

agg = pd.DataFrame(agg_rows)
agg["flfp_pct"]    = agg["flfp"]          * 100
agg["college_pct"] = agg["college_share"] * 100
agg["married_pct"] = agg["married_share"] * 100
agg["treated"]     = ((agg["state"] == "California") & (agg["YEAR"] >= 2004)).astype(int)

p(f"Panel: {agg.shape[0]} rows × {agg.shape[1]} columns")
p(f"States: {agg['state'].nunique()} | Years: {agg['YEAR'].nunique()} ({agg['YEAR'].min()}–{agg['YEAR'].max()})")

# Thin cell check
thin = agg[agg["n_observations"] < 200]
if len(thin):
    p(f"\n[FLAG] {len(thin)} thin cells (<200 obs):")
    p(thin[["YEAR","state","n_observations"]].sort_values("n_observations").to_string(index=False))
else:
    p("✓ All state×year cells have ≥200 observations")

# ── Before/after for panel ────────────────────────────────────────────────────
p("\n" + "=" * 65)
p("  BEFORE / AFTER SUMMARY — PANEL FILE")
p("=" * 65)
p("  Original panel (first pass, 1,071 rows × 22 cols)")
p(f"  This panel (deduped,        {agg.shape[0]} rows × {agg.shape[1]-1} cols)") # -1: no mean_uhrs_worked
p("  Column removed: mean_uhrs_worked (QUARANTINED — implausible values,")
p("                  pending IPUMS Basic Monthly codebook verification)")

# ── Summary statistics ─────────────────────────────────────────────────────────
p("\n" + "=" * 65)
p("  SUMMARY STATISTICS — FINAL PANEL")
p("=" * 65)

stat_vars = [
    ("flfp_pct",       "Female LFP rate (%)"),
    ("college_pct",    "College degree (%)"),
    ("married_pct",    "Married (%)"),
    ("mean_age",       "Mean age (years)"),
    ("nchild_mean",    "Mean # children"),
    ("has_child_pct",  "Has any child (%)"),
    ("yngch_mean",     "Age youngest child (mothers)"),
    ("famsize_mean",   "Mean family size"),
    ("pct_private",    "Private sector (%)"),
    ("pct_govt",       "Govt sector (%)"),
    ("pct_selfempl",   "Self-employed (%)"),
]

print(f"\n{'Variable':<35} {'Mean':>7} {'SD':>7} {'Min':>7} {'Max':>7} {'N':>6}")
print("-" * 75)
for col, label in stat_vars:
    s = agg[col]
    print(f"{label:<35} {s.mean():7.2f} {s.std():7.2f} {s.min():7.2f} {s.max():7.2f} {s.notna().sum():6d}")

# CA vs others
p("\nCalifornia vs. other states (panel means):")
ca   = agg[agg["state"] == "California"]
rest = agg[agg["state"] != "California"]
print(f"\n{'Variable':<25} {'CA':>8} {'Others':>8} {'Diff':>8}")
print("-" * 50)
for col, label in stat_vars[:6]:
    c = ca[col].mean()
    r = rest[col].mean()
    print(f"{label:<25} {c:8.2f} {r:8.2f} {c-r:+8.2f}")

# Trend check: CA FLFP over time
p("\nCalifornia FLFP (%) by year:")
ca_trend = ca.set_index("YEAR")["flfp_pct"].sort_index()
p(ca_trend.to_string())

# ── Save panel ────────────────────────────────────────────────────────────────
panel_cols = [
    "YEAR", "state",
    "flfp", "college_share", "married_share", "mean_age",
    "n_observations", "sum_weights",
    "flfp_pct", "college_pct", "married_pct",
    # mean_uhrs_worked EXCLUDED (quarantined)
    "nchild_mean", "has_child_pct", "yngch_mean",
    "famsize_mean", "pct_head", "pct_private", "pct_govt", "pct_selfempl",
    "occ_era", "treated"
]
panel_out = os.path.join(ROOT, "data/derived/state_year_panel_deduped_1995_2015.csv")
agg[panel_cols].to_csv(panel_out, index=False)
p(f"\nSaved: {panel_out}")
p(f"  {agg.shape[0]} rows × {len(panel_cols)} columns")
p(f"  Columns: {panel_cols}")
p("\nDone.")
