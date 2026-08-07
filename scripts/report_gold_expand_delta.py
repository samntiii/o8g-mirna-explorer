#!/usr/bin/env python3
"""Before/after table: original vs expanded gold_master pooled null results."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "paper" / "benchmarks"
BEFORE = BENCH / "nullmodel_pooled_before_gold_expand.csv"
AFTER = BENCH / "nullmodel_pooled.csv"
OUT = BENCH / "gold_expand_delta.csv"
OUT_MD = BENCH / "gold_expand_delta.md"


def main():
    if not BEFORE.exists():
        raise SystemExit(f"Missing before snapshot: {BEFORE}")
    if not AFTER.exists():
        raise SystemExit(f"Missing after pooled null: {AFTER} — run benchmark_nullmodel.py")

    b = pd.read_csv(BEFORE).assign(set="original")
    a = pd.read_csv(AFTER).assign(set="expanded")
    keys = ["mode", "effect_type"]
    cols = [
        "n_gold",
        "observed",
        "recovery_rate",
        "decoy_p",
        "labelperm_p",
        "agreement_flag",
        "empirical_p",
    ]
    m = b[keys + cols].merge(a[keys + cols], on=keys, suffixes=("_before", "_after"))
    m["crossed_p05"] = (m["decoy_p_before"] >= 0.05) & (m["decoy_p_after"] < 0.05)
    m["delta_obs"] = m["observed_after"] - m["observed_before"]
    m["delta_n_gold"] = m["n_gold_after"] - m["n_gold_before"]
    m.to_csv(OUT, index=False)

    lines = [
        "# Gold expansion — pooled null before/after",
        "",
        f"Before snapshot: `{BEFORE.name}` (Seok+Eom scored set, 28 gained / 12 lost).",
        f"After: `{AFTER.name}` (gold_master included=True).",
        "",
        "| mode | effect | n_before | obs_before | decoy_p_before | n_after | obs_after | decoy_p_after | crossed p<0.05? |",
        "|---|---|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for _, r in m.iterrows():
        lines.append(
            f"| {r['mode']} | {r['effect_type']} | {int(r['n_gold_before'])} | "
            f"{int(r['observed_before'])} | {r['decoy_p_before']:.4g} | "
            f"{int(r['n_gold_after'])} | {int(r['observed_after'])} | "
            f"{r['decoy_p_after']:.4g} | {'YES' if r['crossed_p05'] else 'no'} |"
        )
    crossed = m[m["crossed_p05"]]
    lines += ["", "## Newly significant (decoy null, α=0.05)", ""]
    if len(crossed) == 0:
        lines.append(
            "No previously non-significant cell crossed p<0.05 after expansion "
            "(gained denominators unchanged if Source A/C additions remain excluded)."
        )
    else:
        for _, r in crossed.iterrows():
            lines.append(
                f"- **{r['mode']} {r['effect_type']}**: "
                f"decoy_p {r['decoy_p_before']:.4g} → {r['decoy_p_after']:.4g} "
                f"(obs {int(r['observed_before'])}/{int(r['n_gold_before'])} → "
                f"{int(r['observed_after'])}/{int(r['n_gold_after'])})"
            )
    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT}")
    print(f"Wrote {OUT_MD}")
    print(m[keys + ["n_gold_before", "n_gold_after", "decoy_p_before", "decoy_p_after", "crossed_p05"]].to_string(index=False))


if __name__ == "__main__":
    main()
