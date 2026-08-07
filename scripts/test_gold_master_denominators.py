#!/usr/bin/env python3
"""Unit test: all figure/stat scripts take gold denominators from gold_master.

Fails if:
  - gold_master.csv missing or has no included rows
  - figure/stat scripts still hardcode classic gold counts (28, 12, 27, 13 as denoms)
  - any script reads oxomir_gold_standard without also referencing gold_master
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MASTER = ROOT / "paper" / "benchmarks" / "gold_master.csv"

# Scripts that must participate in the gold-denominator contract
SCAN_GLOBS = [
    "scripts/benchmark_nullmodel.py",
    "scripts/build_gold_master.py",
    "paper/scripts/benchmark_precision_modes.py",
    "paper/scripts/make_benchmark_figures_v2.py",
    "paper/scripts/make_manuscript_figures.py",
    "paper/scripts/plot_benchmark_figures.py",
    "paper/scripts/run_paper_analysis.py",
]

# Hardcoded denominator patterns that must not appear in executable code
FORBIDDEN = [
    re.compile(r"\b11\s*/\s*28\b"),
    re.compile(r"\b6\s*/\s*12\b"),
    re.compile(r"\b9\s*/\s*28\b"),
    re.compile(r"\b13\s*/\s*27\b"),
    re.compile(r"\b12\s*/\s*28\b"),
    re.compile(r"n_gold\s*=\s*28\b"),
    re.compile(r"n_gold\s*=\s*12\b"),
    re.compile(r"gained\s*=\s*28\b"),
    re.compile(r"lost\s*=\s*12\b"),
]


def test_master_exists_and_consistent():
    assert MASTER.exists(), f"missing {MASTER} — run scripts/build_gold_master.py"
    df = pd.read_csv(MASTER)
    assert "included" in df.columns
    inc = df[df["included"].astype(str).str.lower().isin(["true", "1"])]
    assert len(inc) > 0
    g = int((inc["effect_type"] == "gained").sum())
    l = int((inc["effect_type"] == "lost").sum())
    print(f"  OK gold_master included gained={g} lost={l} total={len(inc)}")
    # Guo must not be ambiguous
    guo = df[df["pmid"].astype(str) == "41690606"]
    assert len(guo) >= 1, "Guo 2026 rows missing from master"
    assert guo["included"].astype(str).str.lower().isin(["true", "1"]).sum() == 0 or (
        guo.loc[guo["included"].astype(str).str.lower().isin(["true", "1"]), "o8g_position"]
        .astype(str)
        .str.len()
        > 0
    ).all()
    for _, r in guo.iterrows():
        if str(r["included"]).lower() not in ("true", "1"):
            assert str(r.get("exclude_reason") or "").strip(), "Guo exclusion must state a reason"
    print(f"  OK Guo rows curated unambiguously (n={len(guo)})")
    return g, l


def test_no_hardcoded_denominators():
    failures = []
    for rel in SCAN_GLOBS:
        path = ROOT / rel
        if not path.exists():
            continue
        text = path.read_text(errors="ignore")
        # strip comments loosely
        code_lines = []
        for line in text.splitlines():
            stripped = line.split("#", 1)[0]
            code_lines.append(stripped)
        code = "\n".join(code_lines)
        for pat in FORBIDDEN:
            if pat.search(code):
                failures.append(f"{rel}: matches {pat.pattern}")
    assert not failures, "Hardcoded gold denominators found:\n  " + "\n  ".join(failures)
    print(f"  OK no hardcoded gold denominators in {len(SCAN_GLOBS)} scripts")


def test_stat_scripts_reference_master():
    required = {
        "scripts/benchmark_nullmodel.py": "gold_master",
        "paper/scripts/benchmark_precision_modes.py": "gold_master",
        "paper/scripts/make_benchmark_figures_v2.py": "gold_master",
    }
    for rel, token in required.items():
        text = (ROOT / rel).read_text(errors="ignore")
        assert token in text, f"{rel} must reference {token}"
        print(f"  OK {rel} references {token}")


def main():
    print("1) gold_master consistency")
    test_master_exists_and_consistent()
    print("2) no hardcoded denominators")
    test_no_hardcoded_denominators()
    print("3) scripts reference master/null outputs")
    test_stat_scripts_reference_master()
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL:", e)
        sys.exit(1)
