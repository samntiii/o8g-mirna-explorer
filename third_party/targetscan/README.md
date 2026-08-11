# TargetScan offline site finder (Bartel lab)

Vendored from the community mirror of TargetScan 7.0 prediction code:
https://github.com/nsoranzo/targetscan

Original distribution: http://www.targetscan.org/ (Data Download pages).

## Files

- `targetscan_70.pl` — seed-site finder (TargetScanS algorithm)
- `README_70.txt` — input/output formats

Context++ / PCT scripts are **not** vendored here; this explorer’s
**TargetScan de novo** mode uses the site finder (Python `TargetScanner`
backend by default; optional `O8G_TS_DENOVO_BACKEND=perl`).

## o8G encoding

Oxidized seed Gs are encoded as **T** for Watson–Crick search so sites
require **A** opposite o8G (see `o8g_ts_denovo.wc_seed_for_oxidation`).
