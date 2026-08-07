"""
o8g_genes.py
============
Offline gene ID resolution (symbol / Ensembl / Entrez) using locally cached
alias tables from scripts/fetch_gene_aliases.py.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
ALIAS_PATH = ROOT / "gene_aliases.parquet"
MAP_PATH = ROOT / "gene_id_map.parquet"
META_PATH = ROOT / "gene_aliases.meta.txt"

ID_TYPES = ("Auto", "Symbol", "Ensembl", "Entrez")

_ENSG_RE = re.compile(r"^ENSG\d+", re.I)
_ENTREZ_RE = re.compile(r"^\d+$")


@dataclass
class GeneHit:
    gene_idx: int
    ensembl: str
    entrez: int | None
    symbol: str
    matched_as: str  # which alias string matched
    id_type: str  # resolved type used


@dataclass
class ResolveResult:
    query: str
    id_type: str
    detected_type: str
    hits: list[GeneHit]
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.error is None and len(self.hits) >= 1

    @property
    def ambiguous(self) -> bool:
        return len(self.hits) > 1


def detect_id_type(query: str) -> str:
    q = query.strip()
    if not q:
        return "Symbol"
    if _ENSG_RE.match(q):
        return "Ensembl"
    if _ENTREZ_RE.match(q):
        return "Entrez"
    return "Symbol"


def alias_meta() -> dict[str, str]:
    if not META_PATH.exists():
        return {}
    out = {}
    for line in META_PATH.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k] = v
    return out


class GeneResolver:
    """Lookup gene_idx from symbol / Ensembl / Entrez using local parquets."""

    def __init__(
        self,
        alias_path: Path | str = ALIAS_PATH,
        map_path: Path | str = MAP_PATH,
    ):
        alias_path = Path(alias_path)
        map_path = Path(map_path)
        if not alias_path.exists() or not map_path.exists():
            raise FileNotFoundError(
                f"Missing {alias_path.name} / {map_path.name}. "
                "Run: python scripts/fetch_gene_aliases.py"
            )
        self.aliases = pd.read_parquet(alias_path)
        self.id_map = pd.read_parquet(map_path).set_index("gene_idx")
        self.meta = alias_meta()
        # index: alias_norm -> list of gene_idx (for symbol/general)
        self._by_alias: dict[str, list[int]] = {}
        for gene_idx, alias_norm in zip(
            self.aliases["gene_idx"].tolist(), self.aliases["alias_norm"].tolist()
        ):
            self._by_alias.setdefault(str(alias_norm), []).append(int(gene_idx))
        # dedicated maps
        self._by_ensembl = {
            str(e).upper(): int(i)
            for i, e in zip(self.id_map.index, self.id_map["ensembl"])
            if pd.notna(e)
        }
        self._by_entrez: dict[str, int] = {}
        for i, e in zip(self.id_map.index, self.id_map["entrez"]):
            if pd.notna(e):
                self._by_entrez[str(int(e))] = int(i)

    def _hit_from_idx(self, gene_idx: int, matched_as: str, id_type: str) -> GeneHit:
        row = self.id_map.loc[gene_idx]
        entrez = None if pd.isna(row["entrez"]) else int(row["entrez"])
        return GeneHit(
            gene_idx=int(gene_idx),
            ensembl=str(row["ensembl"]),
            entrez=entrez,
            symbol=str(row["symbol"]),
            matched_as=matched_as,
            id_type=id_type,
        )

    def resolve(self, query: str, id_type: str = "Auto") -> ResolveResult:
        q = (query or "").strip()
        if not q:
            return ResolveResult(query=q, id_type=id_type, detected_type="Symbol", hits=[], error="Empty query")
        if id_type not in ID_TYPES:
            return ResolveResult(query=q, id_type=id_type, detected_type="Symbol", hits=[], error=f"Unknown id_type {id_type}")

        detected = detect_id_type(q)
        use = detected if id_type == "Auto" else id_type

        idxs: list[int] = []
        if use == "Ensembl":
            key = q.split(".")[0].upper()
            if key in self._by_ensembl:
                idxs = [self._by_ensembl[key]]
        elif use == "Entrez":
            if q in self._by_entrez:
                idxs = [self._by_entrez[q]]
        else:  # Symbol — also allow alias / previous symbols
            idxs = list(dict.fromkeys(self._by_alias.get(q.upper(), [])))

        if not idxs:
            return ResolveResult(
                query=q,
                id_type=id_type,
                detected_type=detected,
                hits=[],
                error=f"No gene matched {use} query {q!r} in the local map "
                f"(gene may lack a longest 3′UTR in this database).",
            )

        hits = [self._hit_from_idx(i, matched_as=q, id_type=use) for i in idxs]
        return ResolveResult(query=q, id_type=id_type, detected_type=detected, hits=hits)
