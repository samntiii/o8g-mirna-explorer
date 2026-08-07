# o8G-miRNA Retargeting Explorer

Predict how **8-oxoguanine (o8G) oxidation of a microRNA seed** rewires its mRNA
target list, and compare target lists by pathway enrichment.

Public concept site: **[oxomir.samnti.com](https://oxomir.samnti.com)**
(Cloudflare Tunnel from the lab Mac — see [`DEPLOY.md`](DEPLOY.md)).
Apex `samnti.com` is reserved for an existing Foundry tunnel.

> **Biology.** A seed guanine normally pairs C; oxidized to 8-oxoG it pairs A
> (Hoogsteen). Each of the *k* guanines in a seed can independently be normal
> or o8G, giving 2^k target repertoires per miRNA. Details in
> `VALIDATION_REPORT.md` and the in-app Methods caption.

## Quick start

```bash
cd ~/o8g_mirna_explorer
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./start.sh
# open http://localhost:8501
```

Required next to the code: `o8g_targets.db`, `genesets/` (or `genesets.tar.gz`).
For the **Gene → miRNA/oxomiR** tab also need `gene_aliases.parquet` +
`o8g_reverse.db` (see below). `utr3_human.parquet` enables optional on-the-fly
rescans.

## App surfaces

Single Streamlit app (`app.py`) with four tabs:

- **Single state — targets** — browse one miRNA oxidation state.
- **Compare two states** — lost / gained / shared genes + pathway volcano.
- **All states** — heatmap / enrichment across every o8G combination.
- **Gene → miRNA/oxomiR** — offline reverse lookup: enter a gene (symbol /
  Ensembl / Entrez) and list every strong-site miRNA/oxomiR state that targets
  it, with retargeting vs unmodified.

## Package layout

| File | Role |
|---|---|
| `app.py` | Streamlit explorer (4 tabs) |
| `o8g_engine.py` | Seed extraction, 2^k state enumeration, motifs |
| `o8g_scanner.py` | k-mer index over 3′UTRs |
| `o8g_enrich.py` | Hypergeometric ORA + BH |
| `o8g_plots.py` | Volcano, heatmap, diverging bar, dot plot |
| `o8g_db.py` | Read-only SQLite accessor (+ reverse gene query) |
| `o8g_genes.py` | Local gene ID resolver (NCBI/HGNC map) |
| `precompute_db.py` | Rebuild `o8g_targets.db` |
| `scripts/fetch_mirbase.py` | Human mature miRNAs → `hsa_mature.parquet` |
| `scripts/fetch_utr3.py` | Ensembl BioMart 3′UTRs → `utr3_human.parquet` |
| `scripts/fetch_gene_aliases.py` | NCBI + HGNC → `gene_aliases.parquet` |
| `scripts/build_reverse_index.py` | Invert targets DB → `o8g_reverse.db` |
| `scripts/validate_db.py` | miR-1 / miR-124 sanity checks |

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

Smoke check (HDAC4 should be targeted by unmodified miR-1 / seed `GGAATGT`,
and **absent** at `o8G@7`):

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

Seed-match prediction produces large target lists; interpret pathway output as
**direction of retargeting**, not validated targets. Confirm candidate sites
experimentally. Reverse lookup uses the same strong-site (7mer-m8 + 8mer)
precomputed DB — not live TargetScan/miRDB consensus.
