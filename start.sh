#!/usr/bin/env bash
# Local launch: o8G-miRNA Explorer on http://localhost:8501
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -d genesets ] && [ -f genesets.tar.gz ]; then
  echo "Unpacking gene-set libraries…"
  tar -xzf genesets.tar.gz
fi

if [ ! -f o8g_targets.db ]; then
  echo "ERROR: o8g_targets.db is missing." >&2
  echo "Fetch inputs and rebuild:" >&2
  echo "  source .venv/bin/activate" >&2
  echo "  python scripts/fetch_mirbase.py" >&2
  echo "  python scripts/fetch_utr3.py          # slow (~40 min) if parquet absent" >&2
  echo "  python precompute_db.py" >&2
  exit 1
fi

if [ -x .venv/bin/streamlit ]; then
  PY=.venv/bin/python
  ST=.venv/bin/streamlit
elif command -v streamlit >/dev/null 2>&1; then
  PY=python3
  ST=streamlit
else
  echo "ERROR: streamlit not found. Create .venv and pip install -r requirements.txt" >&2
  exit 1
fi

echo "Starting o8G-miRNA Explorer → http://localhost:8501"
exec "$ST" run app.py \
  --server.port=8501 \
  --server.address=0.0.0.0 \
  --server.headless=true
