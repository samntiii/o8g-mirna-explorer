#!/usr/bin/env python3
"""Regression: precision filters must not erase miR-1 o8G@7 retargeting signal."""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from o8g_db import TargetDB
from o8g_precision import PrecisionMode, assert_retargeting_signal


def main() -> int:
    db = TargetDB(str(ROOT / "o8g_targets.db"))
    seed = "GGAATGT"
    mirna = "hsa-miR-1-3p"
    info = db.mirna_info(mirna)
    assert info and info["seed"] == seed

    scanner = None
    if os.environ.get("O8G_TEST_LIVE_SCANNER", "0") == "1":
        utr = ROOT / "utr3_human.parquet"
        if utr.exists():
            from o8g_scanner import TargetScanner

            print("Loading UTR scanner…", flush=True)
            scanner = TargetScanner.from_parquet(str(utr))

    print("Precision regression — miR-1-3p none vs o8G@7\n", flush=True)
    for mode in PrecisionMode:
        parts = db.retarget_partition(
            seed,
            "o8G@7",
            mode,
            scanner=scanner,
            mature_dna=info["seq_dna"],
            mirna=mirna,
        )
        print(
            f"  {mode.value:10s}  gained={len(parts['gained']):5d}  "
            f"lost={len(parts['lost']):5d}  shared={len(parts['shared']):5d}",
            flush=True,
        )
        mins = {
            PrecisionMode.SENSITIVE: (500, 500, 500),
            PrecisionMode.STRINGENT: (50, 50, 50),
            PrecisionMode.CONSENSUS: (20, 20, 20),
        }[mode]
        assert_retargeting_signal(
            parts,
            min_gained=mins[0],
            min_lost=mins[1],
            min_shared=mins[2],
            label=f"miR-1 o8G@7 [{mode.value}]",
        )
    print("\nPASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
