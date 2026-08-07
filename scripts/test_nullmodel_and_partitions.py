#!/usr/bin/env python3
"""
Sanity tests for upgraded null-model pipeline.

(a) decoy=real-seed reproduces observed gold recovery counts
(b) pooled observed equals sum of per-miRNA observed
(c) both nulls run without altering miR-1 o8G@7 partitions
(d) Wilson CI is within [0,1] and contains the point estimate
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from o8g_db import TargetDB
from o8g_precision import PrecisionMode
from o8g_scanner import TargetScanner
from benchmark_nullmodel import (  # noqa: E402
    load_gold_by_mirna,
    real_seed_gold_recovery,
    observed_partition,
    wilson_ci,
)


def test_mir1_partitions_unchanged():
    db = TargetDB(str(ROOT / "o8g_targets.db"))
    expected = {
        PrecisionMode.SENSITIVE: (2843, 2251, 1066),
        PrecisionMode.STRINGENT: (1484, 1146, 209),
        PrecisionMode.CONSENSUS: (2843, 1257, 650),
    }
    for mode, (g, l, s) in expected.items():
        parts = observed_partition(db, "hsa-miR-1-3p", "GGAATGT", "o8G@7", mode)
        assert len(parts["gained"]) == g, (mode, len(parts["gained"]), g)
        assert len(parts["lost"]) == l, (mode, len(parts["lost"]), l)
        assert len(parts["shared"]) == s, (mode, len(parts["shared"]), s)
        print(f"  OK partitions {mode.value}: gained={g} lost={l} shared={s}")


def test_real_seed_decoy_matches_gold_recovery():
    db = TargetDB(str(ROOT / "o8g_targets.db"))
    scanner = TargetScanner.from_parquet(str(ROOT / "utr3_human.parquet"))
    detail = Path(ROOT / "paper/benchmarks/gold_recovery_by_mode.csv")
    assert detail.exists(), "run benchmark_precision_modes.py first"
    rec = pd.read_csv(detail)
    for mirna in load_gold_by_mirna():
        for mode in PrecisionMode:
            g_scan, l_scan = real_seed_gold_recovery(scanner, db, mirna, mode)
            sub = rec[(rec["mirna"] == mirna) & (rec["mode"] == mode.value)]
            g_db = int(
                sub.loc[
                    (sub["effect"] == "gained_on_oxidation")
                    & (sub["recovered"] == True),  # noqa: E712
                    "gene",
                ].shape[0]
            )
            l_db = int(
                sub.loc[
                    (sub["effect"] == "lost_on_oxidation")
                    & (sub["recovered"] == True),  # noqa: E712
                    "gene",
                ].shape[0]
            )
            assert g_scan == g_db, (mirna, mode, "gained", g_scan, g_db)
            assert l_scan == l_db, (mirna, mode, "lost", l_scan, l_db)
            print(
                f"  OK real-seed decoy {mirna} {mode.value}: gained={g_scan} lost={l_scan}"
            )


def test_pooled_equals_sum_per_mirna():
    per = ROOT / "paper/benchmarks/nullmodel_recovery.csv"
    pool = ROOT / "paper/benchmarks/nullmodel_pooled.csv"
    assert per.exists() and pool.exists(), "run benchmark_nullmodel.py first"
    per_df = pd.read_csv(per)
    pool_df = pd.read_csv(pool)
    for mode in ["Sensitive", "Stringent", "Consensus"]:
        for effect in ["gained", "lost"]:
            s = int(
                per_df[
                    (per_df["mode"] == mode) & (per_df["effect_type"] == effect)
                ]["observed"].sum()
            )
            p = int(
                pool_df[
                    (pool_df["mode"] == mode) & (pool_df["effect_type"] == effect)
                ]["observed"].iloc[0]
            )
            assert s == p, (mode, effect, s, p)
            print(f"  OK pooled {mode} {effect}: {p} == sum per-miRNA")


def test_wilson_basics():
    lo0, hi0 = wilson_ci(0, 3)
    assert 0.0 <= lo0 <= hi0 <= 1.0 and lo0 == 0.0
    lo1, hi1 = wilson_ci(3, 3)
    assert 0.0 <= lo1 <= hi1 <= 1.0
    assert lo1 < 1.0  # Wilson lower bound < 1 even for 3/3
    lo, hi = wilson_ci(11, 28)
    assert lo < 11 / 28 < hi
    print(f"  OK Wilson examples: 0/3→[{lo0:.3f},{hi0:.3f}]; 11/28→[{lo:.3f},{hi:.3f}]")


def test_pooled_has_both_nulls():
    pool = pd.read_csv(ROOT / "paper/benchmarks/nullmodel_pooled.csv")
    for col in ("decoy_p", "labelperm_p", "agreement_flag", "wilson_ci_low", "ci_low"):
        assert col in pool.columns, col
    print(f"  OK pooled columns present ({len(pool)} rows)")


def main() -> int:
    print("1) miR-1 o8G@7 partitions unchanged")
    test_mir1_partitions_unchanged()
    print("2) real seed as decoy reproduces gold recovery")
    test_real_seed_decoy_matches_gold_recovery()
    print("3) Wilson CI sanity")
    test_wilson_basics()
    print("4) pooled == sum per-miRNA (requires prior null run)")
    test_pooled_equals_sum_per_mirna()
    print("5) pooled has both null columns")
    test_pooled_has_both_nulls()
    print("\nPASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
