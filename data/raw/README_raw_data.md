# Raw CPS Data

The raw CPS microdata files are **not committed** to this repository because they exceed GitHub's 100 MB file size limit and are subject to IPUMS CPS terms of use, which prohibit redistribution of the raw microdata.

The derived panel files in `data/derived/` are committed and are sufficient to replicate all estimation results without re-downloading the raw data.

---

## If You Want to Re-Run Data Construction from Scratch

### 1. Create an IPUMS CPS account

Go to https://cps.ipums.org and register for a free account.

### 2. Build the extract

Create a new data extract with the following specifications:

**Sample:** CPS ASEC (Annual Social and Economic Supplement)  
**Years:** 1995–2015 (all ASEC years in this range)  
**Unit of analysis:** Person

**Variables to include:**

| Variable | Description |
|---|---|
| `YEAR` | Survey year |
| `ASECFLAG` | ASEC sample flag (keep all values including 2) |
| `STATEFIP` | State FIPS code |
| `SEX` | Sex |
| `AGE` | Age |
| `EMPSTAT` | Employment status |
| `LABFORCE` | Labor force status |
| `EDUC` | Educational attainment |
| `MARST` | Marital status |
| `WTFINL` | Final person weight |
| `PERNUM` | Person number within household |
| `CPSIDP` | CPS longitudinal person identifier (for deduplication) |

**Select cases:** No case selection in IPUMS — the `01_data_construction.py` script filters to women aged 25–54.

### 3. Download and place files

Download the extract as a `.csv.gz` file. Place it in this directory (`data/raw/`) with the filename:

```
data/raw/cps_asec_1995_2015.csv.gz
```

Or update the `DATA_FILE` path at the top of `code/01_data_construction.py` to point to your downloaded file.

### 4. Run data construction

```bash
python code/01_data_construction.py
```

This produces:
- `data/derived/state_year_panel_deduped_1995_2015.csv` — main estimation panel
- `data/derived/state_year_predictors.csv` — panel for the R synthdid spec

### Notes on the ASEC sample flag

IPUMS CPS includes two overlapping ASEC samples in some years: the main ASEC sample (`ASECFLAG = 1`) and a larger oversample (`ASECFLAG = 2`). The construction script **retains both** (`ASECFLAG` is not filtered), consistent with the approach used in the published analysis. The `CPSIDP`-based deduplication step removes the small number of individuals who appear in both samples in the same year.

### Notes on deduplication

The CPS ASEC sample contains person-level records, and some individuals appear in the rotation group overlap. The script deduplicates using `CPSIDP` (the longitudinal person identifier), keeping the first occurrence per person per year. See `code/01_data_construction.py` for the full deduplication logic.
