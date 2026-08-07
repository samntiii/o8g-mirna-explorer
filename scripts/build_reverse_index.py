#!/usr/bin/env python3
"""Invert o8g_targets.db state→genes into gene→(seed,state,site) reverse index.

Writes o8g_reverse.db with:
  gene_targets(gene_idx, seed, label, site_rank, oxidized_positions)
  seed_mirnas(seed, mirna)   -- all mature names sharing a seed
"""
from __future__ import annotations

import argparse
import sqlite3
import time
import zlib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def build(db_path: Path, out_path: Path):
    src = sqlite3.connect(str(db_path))
    if out_path.exists():
        out_path.unlink()
    dst = sqlite3.connect(str(out_path))
    dst.execute("PRAGMA journal_mode=WAL")
    dst.executescript(
        """
        CREATE TABLE gene_targets(
            gene_idx INTEGER NOT NULL,
            seed TEXT NOT NULL,
            label TEXT NOT NULL,
            site_rank INTEGER NOT NULL,
            oxidized_positions TEXT NOT NULL
        );
        CREATE TABLE seed_mirnas(
            seed TEXT NOT NULL,
            mirna TEXT NOT NULL,
            PRIMARY KEY(seed, mirna)
        );
        """
    )

    mir = src.execute("SELECT seed, mirna FROM mirnas").fetchall()
    dst.executemany("INSERT OR IGNORE INTO seed_mirnas(seed, mirna) VALUES (?,?)", mir)
    dst.commit()
    print(f"seed_mirnas: {len(mir)} rows", flush=True)

    cur = src.execute(
        "SELECT seed, label, oxidized_positions, gene_blob, rank_blob FROM states"
    )
    batch = []
    n_states = 0
    n_pairs = 0
    t0 = time.time()
    while True:
        rows = cur.fetchmany(200)
        if not rows:
            break
        for seed, label, ox, gblob, rblob in rows:
            gidx = np.frombuffer(zlib.decompress(gblob), dtype=np.int32)
            rnk = np.frombuffer(zlib.decompress(rblob), dtype=np.int8)
            ox = ox or ""
            for gi, rk in zip(gidx.tolist(), rnk.tolist()):
                batch.append((int(gi), seed, label, int(rk), ox))
            n_states += 1
            n_pairs += len(gidx)
        if batch:
            dst.executemany(
                "INSERT INTO gene_targets(gene_idx, seed, label, site_rank, oxidized_positions) "
                "VALUES (?,?,?,?,?)",
                batch,
            )
            dst.commit()
            batch.clear()
        if n_states % 1000 == 0:
            print(
                f"  {n_states} states, {n_pairs:,} pairs, {time.time()-t0:.0f}s",
                flush=True,
            )

    print("creating indexes…", flush=True)
    dst.execute("CREATE INDEX idx_gt_gene ON gene_targets(gene_idx)")
    dst.execute("CREATE INDEX idx_gt_seed ON gene_targets(seed)")
    dst.execute("CREATE INDEX idx_sm_seed ON seed_mirnas(seed)")
    dst.commit()
    src.close()
    dst.close()
    print(f"done → {out_path}: {n_states} states, {n_pairs:,} gene-state pairs, {time.time()-t0:.0f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path, default=ROOT / "o8g_targets.db")
    ap.add_argument("--out", type=Path, default=ROOT / "o8g_reverse.db")
    args = ap.parse_args()
    build(args.db, args.out)


if __name__ == "__main__":
    main()
