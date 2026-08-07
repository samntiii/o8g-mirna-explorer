# o8G-miRNA Retargeting Explorer

Predict how **8-oxoguanine (o8G) oxidation of a microRNA seed** rewires its mRNA
target list, compare oxidation states by pathway enrichment, and cross-check
predictions against external unmodified-miRNA catalogs.

**Public site:** [oxomir.samnti.com](https://oxomir.samnti.com)
(alias [o8g.samnti.com](https://o8g.samnti.com)) — Cloudflare Tunnel from the
lab Mac; see [`DEPLOY.md`](DEPLOY.md). Apex `samnti.com` stays on a separate
Foundry tunnel.

> **Biology.** A seed guanine normally pairs **C**; oxidized to 8-oxoG it pairs
> **A** (Hoogsteen). Each of the *k* guanines in a seed can independently be
> normal or o8G, giving **2^k** target repertoires per miRNA. Details in
> `VALIDATION_REPORT.md` and the in-app Methods caption.

## Quick start

```bash
cd ~/o8g_mirna_explorer
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./start.sh
# open http://localhost:8501
```

**Required next to the code**

| Asset | Needed for |
|---|---|
| `o8g_targets.db` | All target browsing |
| `genesets/` (or `genesets.tar.gz`) | Pathway enrichment |
| `gene_aliases.parquet` + `o8g_reverse.db` | Gene → miRNA/oxomiR |
| `utr3_human.parquet` | Optional live UTR rescans / ViennaRNA windows |
| `paper/data/Conserved_Site_Context_Scores.txt` | TargetScan context++ columns |
| ViennaRNA (`RNAduplex` / `RNAup` or Python `ViennaRNA`) | Opt-in energetics |

## App surfaces

Single Streamlit app (`app.py`). Views are a **radio switch** (not tabs) so only
the active section runs — Streamlit re-executes the whole script on every
widget change.

| View | What it does |
|---|---|
| **Single state — targets** | One miRNA oxidation state; BE ranking; optional ViennaRNA / TargetScan energetics |
| **Compare two states** | Lost / gained / shared genes + pathway volcano |
| **All states** | Heatmap / enrichment across every o8G combination |
| **Gene → miRNA/oxomiR** | Reverse lookup (symbol / Ensembl / Entrez); optional duplex / RNAup / context++ on the first *N* rows |
| **External DB comparison** | Explorer vs TargetScan / miRDB / DIANA / miRmap / ENCORI / miRTarBase (UpSet + master list) |

**Sidebar.** Precision mode (**Sensitive** / **Stringent** / **Consensus**),
pathway library, and an optional min site-quality floor.

Energetics checkboxes default **off** (RNAduplex / RNAup / TargetScan scans are
expensive). Turn them on only when you need the columns.

## Binding metrics

| Column | Source | Notes |
|---|---|---|
| `binding_efficiency` (BE) | Engine site score + multiplicity + context analog + conservation | Always available from the DB |
| `dG_RNAduplex` / `dG_hybrid` | ViennaRNA `duplexfold` / `RNAduplex` | Opt-in |
| `dG_RNAup` / `dG_open` / `ddG` | RNAup (hybridization + opening) | Opt-in; slower |
| `contextpp_TargetScan` | TargetScan 8 weighted context++ | Unmodified mature only |

**o8G folding caveat.** ViennaRNA has no native 8-oxoG parameters. For oxidized
states, oxidized seed Gs are temporarily replaced with **U** so Watson–Crick
U·A stands in for Hoogsteen o8G·A. That is a **pairing-code proxy for ranking**,
not a chemical force field — do not treat those kcal/mol as measured o8G ΔG.
Site-type / gained–lost calls use the biological **G·C → o8G·A** rule and do
**not** depend on the U proxy.

## Package layout

| File | Role |
|---|---|
| `app.py` | Streamlit explorer |
| `o8g_engine.py` | Seed extraction, 2^k state enumeration, motifs |
| `o8g_scanner.py` | k-mer index over 3′UTRs |
| `o8g_db.py` | Read-only SQLite accessor (+ reverse gene query) |
| `o8g_precision.py` | Sensitive / Stringent / Consensus filters |
| `o8g_binding.py` | Binding efficiency (BE) |
| `o8g_thermo.py` | ViennaRNA + TargetScan context++ for target tables |
| `o8g_pubthermo.py` | Publication-safe duplex / RNAup / context++ for Gene→miRNA |
| `o8g_refsets.py` | External DB loaders (TargetScan, miRDB, DIANA, miRmap, ENCORI, miRTarBase) |
| `o8g_enrich.py` | Hypergeometric ORA + BH |
| `o8g_plots.py` | Volcano, heatmap, diverging bar, dot plot, UpSet |
| `o8g_genes.py` | Local gene ID resolver (NCBI/HGNC map) |
| `conservation.py` | TargetScan conserved-family helpers for Consensus mode |
| `precompute_db.py` | Rebuild `o8g_targets.db` |
| `start.sh` | Local Streamlit |
| `start_public.sh` | Streamlit + Cloudflare Tunnel (+ wait loop for LaunchAgent) |
| `deploy/com.samnti.o8g-explorer.plist` | macOS LaunchAgent (KeepAlive) |
| `scripts/fetch_mirbase.py` | Human mature miRNAs → `hsa_mature.parquet` |
| `scripts/fetch_utr3.py` | Ensembl BioMart 3′UTRs → `utr3_human.parquet` |
| `scripts/fetch_gene_aliases.py` | NCBI + HGNC → `gene_aliases.parquet` |
| `scripts/build_reverse_index.py` | Invert targets DB → `o8g_reverse.db` |
| `scripts/validate_db.py` | miR-1 / miR-124 sanity checks |
| `scripts/benchmark_nullmodel.py` | Gold recovery null models (local `paper/`) |

## Rebuild

```bash
python scripts/fetch_mirbase.py
python scripts/fetch_utr3.py       # ~40 min if not cached
python precompute_db.py            # ~13 min after the UTR index
python scripts/validate_db.py
```

### Refresh gene aliases + reverse index (monthly / after DB rebuild)

No runtime APIs — download maps once, then query locally:

```bash
python scripts/fetch_gene_aliases.py   # NCBI gene_info + gene2ensembl + HGNC
python scripts/build_reverse_index.py  # invert o8g_targets.db → o8g_reverse.db
```

Smoke check (HDAC4 targeted by unmodified miR-1 / seed `GGAATGT`, **absent** at
`o8G@7`):

```python
from o8g_genes import GeneResolver
from o8g_db import TargetDB

hit = GeneResolver().resolve("HDAC4").hits[0]
df = TargetDB().states_targeting_gene(hit.gene_idx)
mir1 = df[df["seed"] == "GGAATGT"]
assert (mir1["state_label"] == "none").any()
assert not (mir1["state_label"] == "o8G@7").any()
```

## Public deploy

```bash
./start_public.sh                 # foreground wait (LaunchAgent-friendly)
DETACH=1 ./start_public.sh        # fire-and-forget
```

Recommended persistence on the lab Mac:

```bash
cp deploy/com.samnti.o8g-explorer.plist ~/Library/LaunchAgents/
launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/com.samnti.o8g-explorer.plist
```

Unload: `launchctl bootout gui/$(id -u)/com.samnti.o8g-explorer`. Full tunnel
setup: [`DEPLOY.md`](DEPLOY.md).

## Programmatic use

```python
from o8g_engine import enumerate_from_mature
from o8g_db import TargetDB

db = TargetDB("o8g_targets.db")
info = db.mirna_info("hsa-miR-1-3p")
for state in db.states_for_seed(info["seed"]).itertuples():
    genes = db.target_symbols(info["seed"], state.label, min_rank=3)

# reverse: which oxomiR states target HDAC4?
from o8g_genes import GeneResolver
g = GeneResolver().resolve("HDAC4").hits[0]
hits = db.states_targeting_gene(g.gene_idx)
```

## Caveats

- Seed-match prediction produces large target lists; treat pathway output as
  **direction of retargeting**, not validated targets.
- Reverse lookup uses the same strong-site (7mer-m8 + 8mer) precomputed DB —
  not live TargetScan/miRDB consensus.
- External databases catalog **unmodified** mature miRNAs only; only Explorer
  changes with o8G state.
- Oxidized-state ViennaRNA ΔG uses a **G→U proxy** (see above). Prefer site type,
  BE, and gold/null benchmarks for claims; use duplex numbers exploratorily.
