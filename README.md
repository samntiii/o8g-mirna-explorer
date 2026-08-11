# o8G-miRNA Retargeting Explorer

Predict how **8-oxoguanine (o8G) oxidation of a microRNA seed** rewires its mRNA
target list, compare oxidation states by pathway enrichment, and cross-check
predictions against external unmodified-miRNA catalogs.

**Live demo:** [oxomir.samnti.com](https://oxomir.samnti.com)
([o8g.samnti.com](https://o8g.samnti.com))

> **Biology.** A seed guanine normally pairs **C**; oxidized to 8-oxoG it pairs
> **A** (Hoogsteen). Each of the *k* guanines in a seed can independently be
> normal or o8G, giving **2^k** target repertoires per miRNA. See
> `VALIDATION_REPORT.md` and the in-app Methods caption for details.

## Quick start

```bash
git clone https://github.com/samntiii/o8g-mirna-explorer.git
cd o8g-mirna-explorer
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./start.sh
# open http://localhost:8501
```

### Required data (not in git)

Large SQLite / parquet assets are gitignored. Place them next to the code
(or rebuild — see below). Views that need a missing file report which path
is required.

| Asset | Needed for |
|---|---|
| `o8g_targets.db` | All target browsing |
| `genesets/` (or `genesets.tar.gz`) | Pathway enrichment |
| `gene_aliases.parquet` + `o8g_reverse.db` | Gene → miRNA/oxomiR |
| `mirdb_ref.parquet` | Loss-of-function (miRDB-anchored attrition); build with `python scripts/build_mirdb_ref.py` |
| `utr3_human.parquet` | Optional live UTR rescans / ViennaRNA windows |
| `paper/data/Conserved_Family_Info.txt` | Consensus precision mode (TargetScanHuman 8.0) |
| `paper/data/Conserved_Site_Context_Scores.txt` | TargetScan context++ columns |
| ViennaRNA (`RNAduplex` / `RNAup` or Python `ViennaRNA`) | Opt-in energetics |

## App views

Single Streamlit app (`app.py`). Navigation is a **radio switch** (not tabs) so
only the active section runs on each rerun.

| View | Description | Optional data |
|---|---|---|
| **Single state — targets** | One oxidation state; site-type ranking; optional energetics | `utr3_human.parquet`, ViennaRNA |
| **Compare two states** | Venn + lost/gained/shared + pathway volcano with within-pool control + pathway gene lists | `genesets/` |
| **All states** | Lost/gained/shared per ox state (+ optional miRDB attrition); pathway heatmap/dotplot | `genesets/`, `mirdb_ref.parquet` |
| **Overlap (Venn/UpSet)** | Multi-state set overlap after precision filtering | — |
| **Transcription factors** | Enriched TF regulons (TRANSFAC/JASPAR, JASPAR 2025, Genome Browser PWMs, ChEA/ENCODE/TRRUST) among lost/gained | `genesets/*TF*.gmt`, `scripts/fetch_tf_genesets.py` |
| **Antagomir design** | Oxo-selective oligo discrimination (G→U proxy); explains weak full-length folds | ViennaRNA optional |
| **RNA-seq / DEG upload** | Exploratory UP↔lost / DOWN↔gained concordance for a miRNA panel (CSV/TSV/XLSX with symbol, log2FC, padj) | session upload only |
| **Gene → miRNA/oxomiR** | Reverse lookup (symbol / Ensembl / Entrez) | `gene_aliases.parquet`, `o8g_reverse.db` |
| **External DB comparison** | Explorer vs TargetScan / miRDB / DIANA / miRmap / ENCORI / miRTarBase | local ref files |

**Precision modes** (sidebar): **Sensitive** / **Stringent** / **TargetScan** / **Consensus**.

- **TargetScan** anchors the unmodified baseline to TargetScanHuman 8.0 *predicted* strong sites (`Predicted_Targets_Info`) — no conserved-family file required.
- **Consensus** requires TargetScan conserved-family tables. If they are missing, the
app **refuses Consensus** (persistent banner) and falls back to Sensitive — it
does **not** silently intersect with an empty conserved set.

**OBOE-style prior** (sidebar + All states): ranks which seed Gs / ox states to inspect using a local GC-rich motif prior inspired by OBOE; live OBOE web inference is optional and often offline.

Energetics checkboxes default **off** (RNAduplex / RNAup / TargetScan scans are
expensive).

## Binding metrics

| Column | Source | Notes |
|---|---|---|
| `site_type` / `site_rank` | Engine (8mer > 7mer-m8 > 7mer-A1 > 6mer) | Default UI ranking |
| `dG_RNAduplex` / `dG_hybrid` | ViennaRNA `duplexfold` / `RNAduplex` | Opt-in |
| `dG_RNAup` / `dG_open` / `ddG` | RNAup (hybridization + opening) | Opt-in; slower |
| `contextpp_TargetScan` | TargetScan 8 weighted context++ | Unmodified mature only |

**o8G folding caveat.** ViennaRNA has no native 8-oxoG parameters. For oxidized
states, oxidized seed Gs are temporarily replaced with **U** so Watson–Crick
U·A stands in for Hoogsteen o8G·A. That is a **pairing-code proxy for ranking**,
not a chemical force field — do not treat those kcal/mol as measured o8G ΔG.
Site-type / gained–lost calls use the biological **G·C → o8G·A** rule and do
**not** depend on the U proxy.

## Repository layout

| Path | Role |
|---|---|
| `app.py` | Streamlit explorer |
| `o8g_engine.py` | Seed extraction, 2^k state enumeration, motifs |
| `o8g_scanner.py` | k-mer index over 3′UTRs |
| `o8g_db.py` | Read-only SQLite accessor (+ reverse gene query) |
| `o8g_precision.py` | Sensitive / Stringent / Consensus filters |
| `o8g_sections.py` | Section context + dispatch for optional views |
| `o8g_venn.py`, `o8g_tf.py`, `o8g_lof.py`, `o8g_stats.py`, `o8g_anti.py`, `o8g_energy.py` | Optional section views |
| `o8g_thermo.py` | ViennaRNA + TargetScan context++ for target tables |
| `o8g_pubthermo.py` | Duplex / RNAup / context++ for Gene → miRNA |
| `o8g_refsets.py` | External catalog loaders |
| `o8g_enrich.py` | Hypergeometric ORA + BH (+ within-pool control) |
| `o8g_plots.py` | Volcano, heatmap, diverging bar, dot plot, UpSet |
| `o8g_genes.py` | Local gene ID resolver (NCBI/HGNC map) |
| `conservation.py` | TargetScan conserved-family helpers for Consensus |
| `precompute_db.py` | Rebuild `o8g_targets.db` |
| `start.sh` | Local Streamlit launcher |
| `scripts/` | Data fetch, reverse index, validation, benchmarks |

## Rebuild

```bash
python scripts/fetch_mirbase.py
python scripts/fetch_utr3.py       # slow if not cached
python precompute_db.py
python scripts/validate_db.py
```

### Gene aliases + reverse index

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

## Programmatic use

```python
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
- Oxidized-state ViennaRNA ΔG uses a **G→U proxy** (see above). Prefer site type
  and gold/null benchmarks for claims; use duplex numbers exploratorily.

## License

See repository metadata / `LICENSE` if present.
