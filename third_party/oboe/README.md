# OBOE (Xia et al.) — vendored Figshare package

Source: [Figshare 29634239](https://doi.org/10.6084/m9.figshare.29634239)
Paper: Sequence determinant and functional relevance of 8-oxoguanine (o8G) RNA
modification… (OBOE web tool at http://www.rnamd.org/o8GPredictor/).

## Contents

The deposit includes **training / inference scripts and labeled CSVs**, not
published weight files. Layout mirrors the zip:

- `rnabert/` — preferred backbone used by OBOE (`multimolecule/rnabert`)
- `bert/`, `biobert/`, `dnabert/`, `bert2ome/` — alternate fine-tunes from the paper

## How this explorer uses it

1. Fine-tune RNABERT on `rnabert/data/train_0.9.csv` (+ valid/test):
   ```bash
   python scripts/train_oboe_rnabert.py
   ```
2. Checkpoint → `models/oboe_rnabert/checkpoint` (held-out test AUC ~0.98 in our run).
3. Runtime ranking via `o8g_oboe_model.py` / `o8g_oboe.py`: for each seed G,
   score a 51-nt N-padded window → P(o8G), then rank oxidation states.

Cite Xia et al. when reporting OBOE-prioritized LoF oxomiR candidates.
