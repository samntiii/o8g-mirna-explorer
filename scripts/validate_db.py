#!/usr/bin/env python3
"""Sanity-check o8g_targets.db against miR-1 / miR-124 canonical targets."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from o8g_db import TargetDB
from o8g_engine import extract_seed, g_positions

CANON = {
    "hsa-miR-1-3p": ["HDAC4", "TWF1", "GJA1", "KCNJ2", "BDNF", "CDK6", "PURB", "SRSF9", "PTBP1"],
    "hsa-miR-124-3p": ["CDK6", "SP1", "PTBP1", "ROCK1", "VAMP3", "RAB27A", "SLC16A1", "IQGAP1", "CEBPA", "STAT3", "SNAI2"],
}


def main():
    db = TargetDB(str(ROOT / "o8g_targets.db"))
    mir = db.mirnas()
    print(f"miRNAs: {len(mir)}   unique seeds: {mir['seed'].nunique()}")
    n_states = db._con.execute("SELECT COUNT(*) FROM states").fetchone()[0]
    n_genes = db._con.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
    print(f"states: {n_states}   genes: {n_genes}")

    info = db.mirna_info("hsa-miR-1-3p")
    assert info, "hsa-miR-1-3p missing"
    assert info["seed"] == "GGAATGT", info["seed"]
    assert g_positions(info["seed"]) == [2, 3, 7]
    none = set(db.target_symbols("GGAATGT", "none", min_rank=3))
    ox7 = set(db.target_symbols("GGAATGT", "o8G@7", min_rank=3))
    print(f"miR-1 none strong={len(none)}  o8G@7 strong={len(ox7)}")
    for g in CANON["hsa-miR-1-3p"]:
        print(f"  {g:8s} none={'yes' if g in none else 'NO ':3s}  o8G@7={'yes' if g in ox7 else 'no'}")
    assert "HDAC4" in none, "HDAC4 should be a strong-site target of unmodified miR-1"
    assert "HDAC4" not in ox7, "HDAC4 should be lost at o8G@7"
    assert "TWF1" in none and "GJA1" in none

    info124 = db.mirna_info("hsa-miR-124-3p")
    assert info124, "hsa-miR-124-3p missing"
    assert info124["seed"] == "AAGGCAC", info124["seed"]
    none124 = set(db.target_symbols("AAGGCAC", "none", min_rank=3))
    print(f"miR-124 none strong={len(none124)}")
    missing = [g for g in CANON["hsa-miR-124-3p"] if g not in none124]
    for g in CANON["hsa-miR-124-3p"]:
        print(f"  {g:8s} none={'yes' if g in none124 else 'NO'}")
    if missing:
        print("WARNING missing miR-124 strong-site canonicals:", missing)
    else:
        print("miR-124 canonical strong-site recovery: all present")

    print("validation OK")
    db.close()


if __name__ == "__main__":
    main()
