#!/usr/bin/env python3
"""Build mirdb_ref.parquet for the Loss-of-function view.

Reads paper/data/miRDB_v6.0_prediction_result.txt.gz + refseq_to_symbol.tsv
and writes a compact human-only table next to the app:

  mirna, symbol, score

Default keeps score >= 50 so the LoF slider (50–100) can filter without
re-scanning the gzip. Re-run after updating miRDB dumps.
"""
from __future__ import annotations

import gzip
import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper" / "data"
PRED = DATA / "miRDB_v6.0_prediction_result.txt.gz"
MAP = DATA / "refseq_to_symbol.tsv"
OUT = ROOT / "mirdb_ref.parquet"
SCORE_FLOOR = 50.0


def main() -> int:
    if not PRED.exists() or not MAP.exists():
        print(
            f"Missing miRDB inputs.\n  {PRED}\n  {MAP}\n"
            "Download miRDB v6.0 predictions and a RefSeq→symbol map into paper/data/.",
            file=sys.stderr,
        )
        return 1

    mp = dict(
        pd.read_csv(MAP, sep="\t")[["refseq", "symbol"]].itertuples(index=False)
    )
    # strip version suffixes in map keys
    mp = {str(k).split(".")[0]: str(v) for k, v in mp.items()}

    # Optional: restrict to matures present in o8g_targets.db
    keep: set[str] | None = None
    db = ROOT / "o8g_targets.db"
    if db.exists():
        con = sqlite3.connect(db)
        keep = {r[0] for r in con.execute("SELECT mirna FROM mirnas")}
        con.close()
        print(f"Restricting to {len(keep):,} matures in o8g_targets.db")

    rows: list[tuple[str, str, float]] = []
    n_in = n_hsa = n_kept = 0
    with gzip.open(PRED, "rt") as f:
        for line in f:
            n_in += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            mir, ref, score_s = parts[0], parts[1], parts[2]
            if not mir.startswith("hsa-"):
                continue
            n_hsa += 1
            if keep is not None and mir not in keep:
                continue
            try:
                score = float(score_s)
            except ValueError:
                continue
            if score < SCORE_FLOOR:
                continue
            sym = mp.get(ref.split(".")[0])
            if not sym:
                continue
            rows.append((mir, sym.upper(), score))
            n_kept += 1
            if n_kept % 500_000 == 0:
                print(f"  … {n_kept:,} rows", flush=True)

    df = pd.DataFrame(rows, columns=["mirna", "symbol", "score"])
    # collapse duplicate mirna-symbol to max score
    df = df.groupby(["mirna", "symbol"], as_index=False)["score"].max()
    df.to_parquet(OUT, index=False)
    print(
        f"Wrote {OUT}  rows={len(df):,}  miRNAs={df['mirna'].nunique():,}  "
        f"(scanned {n_in:,} lines, {n_hsa:,} hsa)"
    )
    # smoke: miR-1
    n1 = len(df[(df.mirna == "hsa-miR-1-3p") & (df.score >= 80)])
    print(f"Smoke hsa-miR-1-3p score>=80: {n1} genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
