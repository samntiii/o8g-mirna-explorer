# o8G-miRNA Explorer — running on a new machine

## Path A — Run the finished app

```bash
cd ~/o8g_mirna_explorer
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./start.sh
```

Open http://localhost:8501. For **samnti.com** via Cloudflare Tunnel, follow
[`DEPLOY.md`](DEPLOY.md) Option A (`./start_public.sh`).

Needs `o8g_targets.db` + `genesets/` next to the code.

## Path B — Rebuild the database

Only if you change UTRs, site rules, or miRBase snapshot:

```bash
source .venv/bin/activate
python scripts/fetch_mirbase.py
python scripts/fetch_utr3.py          # Ensembl BioMart; ~40 min first time
python precompute_db.py               # ~13 min after the 6-mer index
python scripts/validate_db.py         # miR-1 HDAC4 lost at o8G@7; miR-124 recovery
```

Or paste `REBUILD_PROMPT.md` into a coding agent to reconstruct from scratch.
