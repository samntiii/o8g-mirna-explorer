"""TargetScan de novo predictions for oxidized (and WT) seeds.

Bartel TargetScan Custom has no API. Offline ``targetscan_70.pl`` (vendored under
``third_party/targetscan/``) finds canonical sites by Watson–Crick seed matching.

Biology bridge for o8G
----------------------
o8G(syn)·A is WC-equivalent to replacing that seed G with **T/U** for site search:
TargetScan then looks for A opposite that position. We encode oxidized seeds that
way before calling the site finder.

Backends
--------
``python`` (default) — live ``TargetScanner`` / SeedState motifs (same site types
as TargetScanS: 8mer, 7mer-m8, 7mer-A1, 6mer). Fast; used in the Streamlit UI.

``perl`` — official ``targetscan_70.pl`` on a temp UTR table (human-only, no MSA
gaps). Slower; set ``O8G_TS_DENOVO_BACKEND=perl`` to force.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

import pandas as pd

from o8g_engine import SeedState, clean_seq

ROOT = Path(__file__).resolve().parent
TS_DIR = ROOT / "third_party" / "targetscan"
TS_PL = TS_DIR / "targetscan_70.pl"
HUMAN_TAXID = "9606"

SITE_RANK = {"8mer": 4, "7mer-m8": 3, "7mer-A1": 2, "6mer": 1}
# TargetScan 7 site_type strings in Perl output
_TS_TYPE_MAP = {
    "8mer": "8mer",
    "7mer-m8": "7mer-m8",
    "m8": "7mer-m8",
    "7mer-1a": "7mer-A1",
    "1a": "7mer-A1",
    "6mer": "6mer",
}


def wc_seed_for_oxidation(seed: str, oxidized_positions: tuple[int, ...] | list[int]) -> str:
    """Encode o8G positions as T so WC TargetScan finds o8G·A sites."""
    seed = clean_seq(seed)
    if len(seed) != 7:
        raise ValueError(f"seed must be 7 nt (pos 2–8), got {seed!r}")
    ox = set(int(p) for p in oxidized_positions)
    out = []
    for i, b in enumerate(seed):
        pos = i + 2
        if pos in ox:
            if b != "G":
                raise ValueError(f"cannot oxidize non-G at position {pos} in {seed}")
            out.append("T")
        else:
            out.append(b)
    return "".join(out)


def parse_ox_label(label: str) -> tuple[int, ...]:
    if not label or label == "none":
        return ()
    part = label.replace("o8G@", "")
    return tuple(int(x) for x in part.split(",") if x)


def backend_name() -> str:
    b = (os.environ.get("O8G_TS_DENOVO_BACKEND") or "python").strip().lower()
    return b if b in ("python", "perl") else "python"


def perl_available() -> bool:
    return TS_PL.is_file()


def targets_for_state(
    seed: str,
    label: str,
    *,
    scanner=None,
    min_rank: int = 3,
    family_id: str | None = None,
) -> pd.DataFrame:
    """De novo TargetScan-style targets for one seed oxidation state."""
    ox = parse_ox_label(label)
    state = SeedState(clean_seq(seed), ox)
    if backend_name() == "perl" and perl_available() and scanner is not None:
        try:
            return _targets_perl(state, scanner, min_rank=min_rank, family_id=family_id)
        except Exception:
            # fall through to python
            pass
    return _targets_python(state, scanner, min_rank=min_rank)


def _targets_python(state: SeedState, scanner, *, min_rank: int) -> pd.DataFrame:
    if scanner is None:
        raise RuntimeError(
            "TargetScan de novo (python) needs a TargetScanner — "
            "load utr3_human.parquet via get_scanner()."
        )
    df = scanner.scan_state(state)
    if df is None or df.empty:
        return pd.DataFrame(
            columns=["gene_id", "symbol", "site_rank", "site_type", "source"]
        )
    out = df[df["site_rank"] >= int(min_rank)].copy()
    out["source"] = "targetscan_denovo_python"
    return out.reset_index(drop=True)


def _targets_perl(
    state: SeedState,
    scanner,
    *,
    min_rank: int,
    family_id: str | None,
) -> pd.DataFrame:
    """Run vendored targetscan_70.pl on human UTRs (single-species, no gaps)."""
    fam = family_id or f"ox_{state.label}"
    wc_seed = wc_seed_for_oxidation(state.seed, state.oxidized_positions)
    with tempfile.TemporaryDirectory(prefix="o8g_ts_") as td:
        td_path = Path(td)
        mir_path = td_path / "mir.txt"
        utr_path = td_path / "utr.txt"
        out_path = td_path / "out.txt"
        mir_path.write_text(f"{fam}\t{wc_seed}\t{HUMAN_TAXID}\n")
        # TargetScan UTR format: GeneID, species_ID, sequence
        lines = []
        for gid, sym, u in zip(scanner.genes, scanner.symbols, scanner.utrs):
            # Prefer gene symbol for joinability with explorer tables
            name = str(sym) if sym is not None else str(gid)
            seq = str(u).upper().replace("U", "T")
            if not seq or set(seq) <= {"N", "-"}:
                continue
            lines.append(f"{name}\t{HUMAN_TAXID}\t{seq}\n")
        utr_path.write_text("".join(lines))
        cmd = ["perl", str(TS_PL), str(mir_path), str(utr_path), str(out_path)]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return _parse_ts_output(out_path, min_rank=min_rank)


def _parse_ts_output(path: Path, *, min_rank: int) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(
            columns=["gene_id", "symbol", "site_rank", "site_type", "source"]
        )
    df = pd.read_csv(path, sep="\t")
    # columns from targetscan_70.pl header
    gene_col = "a_Gene_ID" if "a_Gene_ID" in df.columns else df.columns[0]
    type_col = "Site_type" if "Site_type" in df.columns else None
    if type_col is None:
        return pd.DataFrame(
            columns=["gene_id", "symbol", "site_rank", "site_type", "source"]
        )
    rows = []
    best: dict[str, int] = {}
    best_type: dict[str, str] = {}
    for _, r in df.iterrows():
        sym = str(r[gene_col]).upper()
        raw = str(r[type_col]).strip().lower()
        # Perl may emit m8 / 1a / m8:1a group types — map site type only
        mapped = None
        for key, name in _TS_TYPE_MAP.items():
            if raw == key or raw.startswith(key):
                mapped = name
                break
        if mapped is None:
            if "8mer" in raw:
                mapped = "8mer"
            elif "m8" in raw:
                mapped = "7mer-m8"
            elif "1a" in raw or "a1" in raw:
                mapped = "7mer-A1"
            elif "6mer" in raw:
                mapped = "6mer"
        if mapped is None:
            continue
        rank = SITE_RANK[mapped]
        if rank > best.get(sym, 0):
            best[sym] = rank
            best_type[sym] = mapped
    for sym, rank in best.items():
        if rank < min_rank:
            continue
        rows.append(
            dict(
                gene_id=sym,
                symbol=sym,
                site_rank=rank,
                site_type=best_type[sym],
                source="targetscan_denovo_perl",
            )
        )
    return pd.DataFrame(rows)


@lru_cache(maxsize=256)
def _cached_python_key(seed: str, label: str, min_rank: int, scanner_id: int):
    return (seed, label, min_rank, scanner_id)


def partition_denovo(
    seed: str,
    ox_label: str,
    *,
    scanner,
    min_rank: int = 3,
) -> dict[str, set[str]]:
    """Lost/gained/shared using de novo TargetScan lists on both states."""
    u = targets_for_state(seed, "none", scanner=scanner, min_rank=min_rank)
    o = targets_for_state(seed, ox_label, scanner=scanner, min_rank=min_rank)
    su = set(u["symbol"].astype(str).str.upper()) if len(u) else set()
    so = set(o["symbol"].astype(str).str.upper()) if len(o) else set()
    return {
        "unmod": su,
        "oxid": so,
        "shared": su & so,
        "lost": su - so,
        "gained": so - su,
    }
