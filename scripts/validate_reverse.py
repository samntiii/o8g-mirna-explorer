#!/usr/bin/env python3
"""Validate reverse gene lookup is a faithful invert of o8g_targets.db (no hallucination)."""
from __future__ import annotations

import random
import sqlite3
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from o8g_db import RANK_SITE, TargetDB  # noqa: E402
from o8g_genes import GeneResolver, detect_id_type  # noqa: E402

random.seed(42)
PASS: list[str] = []
FAIL: list[str] = []
WARN: list[str] = []


def ok(m: str) -> None:
    PASS.append(m)
    print(f"  PASS  {m}", flush=True)


def bad(m: str) -> None:
    FAIL.append(m)
    print(f"  FAIL  {m}", flush=True)


def warn(m: str) -> None:
    WARN.append(m)
    print(f"  WARN  {m}", flush=True)


def main() -> int:
    t0 = time.time()
    db = TargetDB(str(ROOT / "o8g_targets.db"))
    rev = sqlite3.connect(str(ROOT / "o8g_reverse.db"))
    src = sqlite3.connect(str(ROOT / "o8g_targets.db"))
    for c in (rev, src):
        c.execute("PRAGMA temp_store=MEMORY")
    resolver = GeneResolver()

    print("\n=== 1. INDEX ORIENTATION ===", flush=True)
    n_states = src.execute("SELECT COUNT(*) FROM states").fetchone()[0]
    n_genes = src.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
    n_mir = src.execute("SELECT COUNT(*) FROM mirnas").fetchone()[0]
    n_sm = rev.execute("SELECT COUNT(*) FROM seed_mirnas").fetchone()[0]
    print(
        f"  {n_mir} mirnas, {n_genes} genes, {n_states} states, seed_mirnas={n_sm}",
        flush=True,
    )

    sm = pd.read_sql("SELECT seed, mirna FROM seed_mirnas", rev)
    mir = pd.read_sql("SELECT seed, mirna FROM mirnas", src)
    diff = len(
        sm.merge(mir, on=["seed", "mirna"], how="outer", indicator=True).query(
            "_merge!='both'"
        )
    )
    (ok if diff == 0 and len(sm) == len(mir) else bad)(
        f"seed_mirnas exact vs mirnas (diff={diff})"
    )

    bad_seeds = rev.execute(
        "SELECT COUNT(*) FROM seed_mirnas WHERE length(seed)!=7 OR seed GLOB '*[^ACGT]*'"
    ).fetchone()[0]
    bad_names = rev.execute(
        "SELECT COUNT(*) FROM seed_mirnas WHERE mirna NOT LIKE 'hsa-%'"
    ).fetchone()[0]
    (ok if bad_seeds == 0 and bad_names == 0 else bad)(
        f"column orientation (bad_seeds={bad_seeds}, bad_names={bad_names})"
    )

    print("\n=== 2. FORWARD → REVERSE (40 states × ≤25 genes) ===", flush=True)
    ids = [r[0] for r in src.execute("SELECT state_id FROM states").fetchall()]
    miss = rank_bad = checked = 0
    for sid in random.sample(ids, 40):
        seed, label, gblob, rblob = src.execute(
            "SELECT seed, label, gene_blob, rank_blob FROM states WHERE state_id=?",
            [sid],
        ).fetchone()
        gidx = np.frombuffer(zlib.decompress(gblob), dtype=np.int32)
        rnk = np.frombuffer(zlib.decompress(rblob), dtype=np.int8)
        pick = random.sample(range(len(gidx)), min(25, len(gidx)))
        for i in pick:
            checked += 1
            row = rev.execute(
                "SELECT site_rank FROM gene_targets WHERE gene_idx=? AND seed=? AND label=?",
                [int(gidx[i]), seed, label],
            ).fetchone()
            if row is None:
                miss += 1
            elif row[0] != int(rnk[i]):
                rank_bad += 1
    (ok if miss == 0 and rank_bad == 0 else bad)(
        f"forward→reverse exact on {checked} pairs (miss={miss}, rank_bad={rank_bad})"
    )

    print("\n=== 3. REVERSE → FORWARD (≤200 rows) ===", flush=True)
    g_sample = random.sample(range(n_genes), 40)
    rows: list[tuple] = []
    for gi in g_sample:
        rows.extend(
            rev.execute(
                "SELECT gene_idx, seed, label, site_rank FROM gene_targets "
                "WHERE gene_idx=? LIMIT 5",
                [gi],
            ).fetchall()
        )
    rows = rows[:200]
    rmiss = rrank = 0
    for gi, seed, label, rk in rows:
        r = src.execute(
            "SELECT gene_blob, rank_blob FROM states WHERE seed=? AND label=?",
            [seed, label],
        ).fetchone()
        if not r:
            rmiss += 1
            continue
        gidx = np.frombuffer(zlib.decompress(r[0]), dtype=np.int32)
        rnk = np.frombuffer(zlib.decompress(r[1]), dtype=np.int8)
        pos = np.where(gidx == gi)[0]
        if len(pos) == 0:
            rmiss += 1
        elif int(rnk[pos[0]]) != int(rk):
            rrank += 1
    (ok if rmiss == 0 and rrank == 0 else bad)(
        f"reverse→forward on {len(rows)} rows (miss={rmiss}, rank_bad={rrank})"
    )

    print("\n=== 4. PER-STATE COUNTS (100 states) ===", flush=True)
    count_bad = 0
    for sid in random.sample(ids, 100):
        seed, label, gblob = src.execute(
            "SELECT seed, label, gene_blob FROM states WHERE state_id=?", [sid]
        ).fetchone()
        n = len(np.frombuffer(zlib.decompress(gblob), dtype=np.int32))
        n2 = rev.execute(
            "SELECT COUNT(*) FROM gene_targets WHERE seed=? AND label=?",
            [seed, label],
        ).fetchone()[0]
        if n != n2:
            count_bad += 1
    (ok if count_bad == 0 else bad)(
        f"gene count match on 100 states (mismatches={count_bad})"
    )

    print("\n=== 5. API / HDAC4 / miR-1 BIOLOGY ===", flush=True)
    hit = resolver.resolve("HDAC4").hits[0]
    api = db.states_targeting_gene(hit.gene_idx)
    raw_n = rev.execute(
        "SELECT COUNT(*) FROM gene_targets WHERE gene_idx=?", [hit.gene_idx]
    ).fetchone()[0]
    (ok if len(api) == raw_n else bad)(f"API rows == SQL for HDAC4 ({len(api)})")

    bad_mir = sum(
        1
        for s, m in zip(api["seed"], api["mirna"])
        if not m
        or rev.execute(
            "SELECT 1 FROM seed_mirnas WHERE seed=? AND mirna=?", [s, m]
        ).fetchone()
        is None
    )
    (ok if bad_mir == 0 else bad)(f"mirna belongs to seed ({bad_mir} bad)")
    (ok if (api["site_type"] == api["site_rank"].map(RANK_SITE)).all() else bad)(
        "site_type from rank"
    )

    logic_bad = 0
    for seed, sub in api.groupby("seed"):
        has_none = (sub["state_label"] == "none").any()
        for _, row in sub.iterrows():
            lab = row["state_label"]
            expect = (
                "unmodified"
                if lab == "none"
                else ("also in unmodified" if has_none else "gained on oxidation")
            )
            if row["vs_unmodified"] != expect:
                logic_bad += 1
    (ok if logic_bad == 0 else bad)(
        f"vs_unmodified deterministic ({logic_bad} bad / {len(api)})"
    )

    unmod = set(db.target_symbols("GGAATGT", "none", min_rank=3))
    ox7 = set(db.target_symbols("GGAATGT", "o8G@7", min_rank=3))
    mir1 = api[api["seed"] == "GGAATGT"]
    (ok if "HDAC4" in unmod else bad)("forward: HDAC4 in unmodified miR-1")
    (ok if "HDAC4" not in ox7 else bad)("forward: HDAC4 NOT in o8G@7")
    (ok if (mir1["state_label"] == "none").any() else bad)(
        "reverse: HDAC4 has unmodified miR-1"
    )
    (ok if not (mir1["state_label"] == "o8G@7").any() else bad)(
        "reverse: HDAC4 lacks o8G@7"
    )

    xf = 0
    for _, row in mir1.iterrows():
        if "HDAC4" not in set(
            db.target_symbols("GGAATGT", row["state_label"], min_rank=3)
        ):
            xf += 1
    (ok if xf == 0 else bad)(
        f"all reverse miR-1 HDAC4 states confirmed forward ({xf} miss)"
    )

    st = db.states_for_seed("GGAATGT")
    none_m = st.loc[st.label == "none", "motif_7mer_m8"].iloc[0]
    ox7_m = st.loc[st.label == "o8G@7", "motif_7mer_m8"].iloc[0]
    (ok if none_m != ox7_m else bad)(f"o8G@7 changes motif {none_m} → {ox7_m}")
    (ok if len(st) == 8 else bad)(f"miR-1 has 2^3=8 states (got {len(st)})")

    print("\n=== 6. GENE ID RESOLUTION ===", flush=True)
    for q, ens, entrez in [
        ("HDAC4", "ENSG00000068024", 9759),
        ("ENSG00000068024", "ENSG00000068024", 9759),
        ("9759", "ENSG00000068024", 9759),
        ("hdac4", "ENSG00000068024", 9759),
    ]:
        h = resolver.resolve(q).hits[0]
        (
            ok
            if h.ensembl == ens and h.entrez == entrez and h.symbol == "HDAC4"
            else bad
        )(f"resolve({q!r})")
    for q, e in [
        ("ENSG00000068024", "Ensembl"),
        ("9759", "Entrez"),
        ("HDAC4", "Symbol"),
    ]:
        (ok if detect_id_type(q) == e else bad)(f"detect({q!r})={detect_id_type(q)}")
    row = src.execute(
        "SELECT gene_idx,symbol FROM genes WHERE gene_id='ENSG00000068024'"
    ).fetchone()
    (ok if row[0] == hit.gene_idx else bad)(
        f"gene_idx {hit.gene_idx} matches genes table"
    )

    alias_idxs = set(pd.read_parquet(ROOT / "gene_aliases.parquet")["gene_idx"].tolist())
    actual = set(pd.read_sql("SELECT gene_idx FROM genes", src)["gene_idx"])
    (ok if not (alias_idxs - actual) else bad)(
        f"aliases ⊆ DB genes (orphans={len(alias_idxs - actual)})"
    )
    miss_cov = actual - alias_idxs
    (ok if not miss_cov else warn)(f"DB genes without aliases: {len(miss_cov)}")

    print("\n=== 7. ARTIFACT CHECKS ===", flush=True)
    dup_n = 0
    for gi in random.sample(range(n_genes), 200):
        d = rev.execute(
            "SELECT COUNT(*) FROM ("
            "SELECT seed,label,COUNT(*) c FROM gene_targets "
            "WHERE gene_idx=? GROUP BY seed,label HAVING c>1)",
            [gi],
        ).fetchone()[0]
        dup_n += d
    (ok if dup_n == 0 else bad)(
        "no duplicate (seed,label) per gene in 200-gene sample"
    )

    mn, mx = rev.execute(
        "SELECT MIN(gene_idx), MAX(gene_idx) FROM gene_targets"
    ).fetchone()
    (ok if mn >= 0 and mx < n_genes else bad)(
        f"gene_idx range [{mn},{mx}] < {n_genes}"
    )

    ranks = [
        r[0]
        for r in rev.execute(
            "SELECT site_rank FROM gene_targets WHERE gene_idx=? LIMIT 1000",
            [hit.gene_idx],
        )
    ]
    (ok if min(ranks) >= 3 and max(ranks) <= 4 else warn)(
        f"HDAC4 site_ranks in {{{min(ranks)}..{max(ranks)}}} (strong only)"
    )

    weak = 0
    for sid in random.sample(ids, 30):
        rblob = src.execute(
            "SELECT rank_blob FROM states WHERE state_id=?", [sid]
        ).fetchone()[0]
        if (np.frombuffer(zlib.decompress(rblob), dtype=np.int8) < 3).any():
            weak += 1
    (ok if weak == 0 else warn)(f"state blobs strong-only ({weak}/30 had rank<3)")

    g_test = 500
    api2 = db.states_targeting_gene(g_test)
    sym = db.symbols[g_test]
    cf = 0
    for _, row in api2.sample(min(40, len(api2)), random_state=1).iterrows():
        if sym not in set(
            db.target_symbols(row["seed"], row["state_label"], min_rank=3)
        ):
            cf += 1
    (ok if cf == 0 else bad)(
        f"{sym} (idx {g_test}): {min(40, len(api2))} reverse hits ⊆ forward"
    )

    gained = api[api["vs_unmodified"] == "gained on oxidation"].head(20)
    g_bad = 0
    for _, row in gained.iterrows():
        if "HDAC4" in set(db.target_symbols(row["seed"], "none", min_rank=3)):
            g_bad += 1
    (ok if g_bad == 0 else bad)(
        f"'gained' rows truly absent from unmodified forward ({g_bad}/≤20 false)"
    )

    print(f"\n=== SUMMARY ({time.time() - t0:.1f}s) ===", flush=True)
    print(f"PASS={len(PASS)}  WARN={len(WARN)}  FAIL={len(FAIL)}", flush=True)
    for m in FAIL:
        print(f"  FAIL: {m}")
    for m in WARN:
        print(f"  WARN: {m}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
