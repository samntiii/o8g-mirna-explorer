#!/usr/bin/env python3
"""
Migrate o8g_targets.db to precision schema v2.

Adds optional parallel blobs on ``states``:
  n8_blob, n7m8_blob  — per-gene site counts (int8)
  cons_blob           — TargetScan-conserved flag for unmodified family (int8)

Also creates ``meta(key,value)`` with schema_version=2 and targetscan_release.

Usage:
  python scripts/migrate_precision_schema.py [--db o8g_targets.db] [--limit-seeds N]

Re-scan uses utr3_human.parquet. Conservation from paper/data TargetScan 8.0 tables.
Backward-compatible: readers tolerate missing blobs (o8g_db.py).
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
import time
import zlib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from conservation import (  # noqa: E402
    TARGETSCAN_RELEASE,
    build_seed_family_map,
    get_conserved_index,
)
from o8g_engine import SeedState, enumerate_states  # noqa: E402
from o8g_scanner import TargetScanner  # noqa: E402

STRONG = 3


def ensure_columns(con: sqlite3.Connection) -> None:
    cols = {r[1] for r in con.execute("PRAGMA table_info(states)").fetchall()}
    for col, typ in (
        ("n8_blob", "BLOB"),
        ("n7m8_blob", "BLOB"),
        ("cons_blob", "BLOB"),
    ):
        if col not in cols:
            con.execute(f"ALTER TABLE states ADD COLUMN {col} {typ}")
    con.execute(
        "CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT)"
    )
    con.commit()


def migrate(db_path: Path, utr_path: Path, limit_seeds: int | None = None) -> None:
    print(f"Loading UTRs from {utr_path}…", flush=True)
    sc = TargetScanner.from_parquet(str(utr_path))
    cons = get_conserved_index(ROOT / "paper" / "data")
    cons.ensure_loaded()
    fam_map = build_seed_family_map(db_path)

    con = sqlite3.connect(str(db_path))
    ensure_columns(con)
    mir = {
        r[0]: r[1]
        for r in con.execute("SELECT seed, mirna FROM mirnas").fetchall()
    }
    # pick a representative mature per seed
    seed_mirna: dict[str, str] = {}
    for seed, m in con.execute("SELECT seed, mirna FROM mirnas").fetchall():
        seed_mirna.setdefault(seed, m)

    seeds = [r[0] for r in con.execute("SELECT DISTINCT seed FROM states").fetchall()]
    if limit_seeds:
        seeds = seeds[:limit_seeds]
    print(f"Migrating {len(seeds)} seeds…", flush=True)
    t0 = time.time()
    n_states = 0
    for si, seed in enumerate(seeds):
        mirna = seed_mirna.get(seed, "")
        fam = fam_map.get(seed)
        conserved = cons.conserved_symbols_for_mirna(mirna, fam) if mirna else set()
        for st in enumerate_states(seed):
            idx, rank, n8, n7 = sc.scan_state_arrays_full(st)
            strong = rank >= STRONG
            s_idx, s_rank = idx[strong], rank[strong]
            s_n8, s_n7 = n8[strong], n7[strong]
            # conservation: gene symbol in TargetScan conserved set (family/mirna)
            # Meaningful for unmodified motifs; for oxidized states still stored
            # (usually False) so both states have the column for identical schemas.
            cons_arr = np.array(
                [1 if sc.symbols[int(i)] in conserved else 0 for i in s_idx],
                dtype=np.int8,
            )
            con.execute(
                "UPDATE states SET n8_blob=?, n7m8_blob=?, cons_blob=?, "
                "gene_blob=?, rank_blob=? WHERE seed=? AND label=?",
                (
                    zlib.compress(s_n8.tobytes()),
                    zlib.compress(s_n7.tobytes()),
                    zlib.compress(cons_arr.tobytes()),
                    zlib.compress(s_idx.astype(np.int32).tobytes()),
                    zlib.compress(s_rank.astype(np.int8).tobytes()),
                    seed,
                    st.label,
                ),
            )
            n_states += 1
        if si % 50 == 0:
            con.commit()
            print(f"  {si}/{len(seeds)} seeds, {n_states} states, {time.time()-t0:.0f}s", flush=True)
    con.execute(
        "INSERT OR REPLACE INTO meta(key,value) VALUES (?,?)",
        ("schema_version", "2"),
    )
    con.execute(
        "INSERT OR REPLACE INTO meta(key,value) VALUES (?,?)",
        ("targetscan_release", TARGETSCAN_RELEASE),
    )
    con.execute(
        "INSERT OR REPLACE INTO meta(key,value) VALUES (?,?)",
        ("precision_migrated_at", time.strftime("%Y-%m-%dT%H:%M:%S")),
    )
    con.commit()
    con.close()
    print(f"done: {n_states} states → schema v2 ({time.time()-t0:.0f}s)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path, default=ROOT / "o8g_targets.db")
    ap.add_argument("--utr", type=Path, default=ROOT / "utr3_human.parquet")
    ap.add_argument("--limit-seeds", type=int, default=None)
    args = ap.parse_args()
    if not args.utr.exists():
        sys.exit(f"UTR parquet missing: {args.utr}")
    migrate(args.db, args.utr, args.limit_seeds)


if __name__ == "__main__":
    main()
