# DB schema migration — precision v2

## Summary

Precision filters (conservation, multiplicity / `SITE_WEIGHT` score, context-style
features) are applied **identically** to unmodified and oxidized target tables;
gained/lost/shared are always set differences **after** filtering
(`o8g_precision.partition_after_filter`).

## TargetScan release (Consensus)

- **TargetScanHuman 8.0** (`conservation.TARGETSCAN_RELEASE = TargetScanHuman_8.0`)
- Files (from https://www.targetscan.org/vert_80/vert_80_data_download/):
  - `Conserved_Family_Info.txt`
  - `Conserved_Site_Context_Scores.txt`
- We do **not** recompute phastCons; conservation ≡ membership in these tables
  (human taxid 9606; strong seed matches 8mer / 7mer-m8).

## Schema

### v1 (original)

`states(gene_blob, rank_blob)` — zlib int32 gene_idx + int8 best rank for rank≥3.

### v2 (precision)

Same as v1 plus optional columns (nullable for backward-compatible reads):

| Column | Type | Meaning |
|---|---|---|
| `n8_blob` | BLOB | zlib int8 `n_8mer` per gene (parallel to gene_blob) |
| `n7m8_blob` | BLOB | zlib int8 `n_7mer_m8` |
| `cons_blob` | BLOB | zlib int8 `is_conserved` (TargetScan family/mature) |

`meta(key,value)` with at least:

- `schema_version=2`
- `targetscan_release=TargetScanHuman_8.0`

`o8g_db.TargetDB` detects columns via `PRAGMA table_info` and degrades gracefully
when blobs are absent (approximates n8/n7 from best rank; conservation can be
joined at query time from `conservation.ConservedIndex`).

## How to migrate / rebuild

**Option A — migrate existing DB** (re-scans UTRs):

```bash
python scripts/migrate_precision_schema.py
# optional: --limit-seeds 20 for a smoke migrate
```

**Option B — full rebuild** (`precompute_db.WRITE_PRECISION_BLOBS = True`):

```bash
python precompute_db.py
python scripts/build_reverse_index.py   # if gene lists changed
```

Reverse index (`o8g_reverse.db`) is unchanged in structure; precision modes filter
at query time from `o8g_targets.db`.

## Modes (`o8g_precision.PrecisionMode`)

| Mode | Rule | Default use |
|---|---|---|
| Sensitive | rank ≥ 3 | UI discovery |
| Stringent | rank == 4 | Paper / Fig.4 |
| Consensus | rank ≥ 3 + TS conserved on **unmodified**; ox uses rank gate | Main-text claims |

`HighConf` was removed: without multiplicity blobs it collapsed to Stringent
(`score ≥ 1.0` ≈ one 8mer). `PrecisionConfig.from_mode("HighConf")` still
aliases to Stringent. SITE_WEIGHT / context scores remain available on live
scanner rows for ranking/export, not as a separate UI tier.

## Feature flags

- `O8G_USE_CONSERVATION=0|1`

## Tests

```bash
python scripts/test_precision_regression.py
python paper/scripts/benchmark_precision_modes.py
```
