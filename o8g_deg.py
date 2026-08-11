"""Parse user-uploaded DEG tables (CSV / TSV / Excel) into UP/DOWN gene sets."""
from __future__ import annotations

from dataclasses import dataclass, field
from io import BytesIO

import pandas as pd

MAX_ROWS = 50_000
MAX_BYTES = 20 * 1024 * 1024

SYMBOL_ALIASES = (
    "symbol",
    "gene",
    "gene_symbol",
    "gene_name",
    "genename",
    "hgnc_symbol",
    "external_gene_name",
)
LFC_ALIASES = (
    "log2foldchange",
    "log2fc",
    "logfc",
    "log2_fc",
    "lfc",
    "log2fold_change",
    "foldchange",
)
PADJ_ALIASES = (
    "padj",
    "p_adj",
    "padjust",
    "fdr",
    "qvalue",
    "q_value",
    "adj_pvalue",
    "adj.p.val",
    "p.adjust",
)


@dataclass
class DegSets:
    table: pd.DataFrame
    up: set[str]
    down: set[str]
    unmatched: int
    n_input: int
    symbol_col: str
    lfc_col: str
    padj_col: str
    notes: list[str] = field(default_factory=list)


def _norm_col(c: str) -> str:
    return str(c).strip().lower().replace(" ", "_").replace("-", "_")


def detect_columns(columns: list[str]) -> dict[str, str | None]:
    norms = {_norm_col(c): c for c in columns}
    out: dict[str, str | None] = {"symbol": None, "lfc": None, "padj": None}
    for a in SYMBOL_ALIASES:
        if a in norms:
            out["symbol"] = norms[a]
            break
    for a in LFC_ALIASES:
        if a in norms:
            out["lfc"] = norms[a]
            break
    for a in PADJ_ALIASES:
        if a in norms:
            out["padj"] = norms[a]
            break
    # index-as-symbol fallback: first column if still missing
    if out["symbol"] is None and columns:
        out["symbol"] = columns[0]
    return out


_DEG_SHEET_HINTS = (
    "deg",
    "degs",
    "differential",
    "deseq",
    "edger",
    "limma",
    "results",
    "result",
    "sig",
    "significant",
    "all_genes",
    "gene",
)


def is_excel_filename(filename: str) -> bool:
    return (filename or "").lower().endswith((".xlsx", ".xls"))


def list_excel_sheets(data: bytes, filename: str | None = None) -> list[str]:
    """Return sheet names for an Excel upload; empty list for non-Excel."""
    if filename is not None and not is_excel_filename(filename):
        return []
    if len(data) > MAX_BYTES:
        raise ValueError(f"File exceeds {MAX_BYTES // (1024 * 1024)} MB limit.")
    bio = BytesIO(data)
    try:
        xl = pd.ExcelFile(bio)
    except Exception as e:
        raise ValueError(f"Could not read Excel workbook: {e}") from e
    return [str(s) for s in xl.sheet_names]


def guess_deg_sheet(sheet_names: list[str]) -> str | None:
    """Pick a likely DEG results sheet; else first sheet."""
    if not sheet_names:
        return None
    norms = [(s, str(s).strip().lower().replace(" ", "_").replace("-", "_")) for s in sheet_names]
    for hint in _DEG_SHEET_HINTS:
        for original, n in norms:
            if hint in n:
                return original
    return sheet_names[0]


def read_upload(
    data: bytes,
    filename: str,
    *,
    sheet_name: str | int | None = None,
) -> pd.DataFrame:
    if len(data) > MAX_BYTES:
        raise ValueError(f"File exceeds {MAX_BYTES // (1024 * 1024)} MB limit.")
    name = (filename or "").lower()
    bio = BytesIO(data)
    if name.endswith((".xlsx", ".xls")):
        # Default: first sheet (or caller-selected sheet via dropdown)
        df = pd.read_excel(bio, sheet_name=0 if sheet_name is None else sheet_name)
        if isinstance(df, dict):
            # Defensive: pandas returns dict only if sheet_name=None with some versions
            first = next(iter(df.values()))
            df = first
    elif name.endswith(".tsv") or name.endswith(".txt"):
        df = pd.read_csv(bio, sep="\t")
    else:
        # try csv; fall back to tab
        try:
            df = pd.read_csv(bio)
        except Exception:
            bio.seek(0)
            df = pd.read_csv(bio, sep="\t")
    if df.shape[1] > 40 and df.shape[0] > 100:
        # Heuristic: raw count matrix (genes × many samples)
        raise ValueError(
            "Table looks like a wide expression/count matrix. "
            "Upload a differential-expression result with symbol, log2FC, and padj columns."
        )
    if len(df) > MAX_ROWS:
        df = df.head(MAX_ROWS)
    return df


def build_deg_sets(
    df: pd.DataFrame,
    *,
    symbol_col: str,
    lfc_col: str,
    padj_col: str,
    lfc_thr: float = 0.5,
    padj_thr: float = 0.05,
    universe: set[str] | None = None,
) -> DegSets:
    notes: list[str] = []
    work = df[[symbol_col, lfc_col, padj_col]].copy()
    work.columns = ["symbol", "lfc", "padj"]
    work["symbol"] = work["symbol"].astype(str).str.strip().str.upper()
    work["lfc"] = pd.to_numeric(work["lfc"], errors="coerce")
    work["padj"] = pd.to_numeric(work["padj"], errors="coerce")
    work = work.dropna(subset=["symbol", "lfc", "padj"])
    work = work[work["symbol"].str.len() > 0]
    work = work.drop_duplicates("symbol", keep="first")

    sig = work[work["padj"] < float(padj_thr)]
    up = set(sig.loc[sig["lfc"] >= float(lfc_thr), "symbol"])
    down = set(sig.loc[sig["lfc"] <= -float(lfc_thr), "symbol"])

    unmatched = 0
    if universe is not None:
        uni = {u.upper() for u in universe}
        n_before = len(up) + len(down)
        up &= uni
        down &= uni
        # genes significant but not in UTR index
        all_sig = set(sig["symbol"])
        unmatched = len(all_sig - uni)
        if unmatched:
            notes.append(
                f"{unmatched} significant genes not in the 3′UTR index were dropped from scoring."
            )
        if n_before and not (up or down):
            notes.append("All DEGs fell outside the UTR universe after filtering.")

    return DegSets(
        table=work,
        up=up,
        down=down,
        unmatched=unmatched,
        n_input=len(df),
        symbol_col=symbol_col,
        lfc_col=lfc_col,
        padj_col=padj_col,
        notes=notes,
    )


def parse_upload(
    data: bytes,
    filename: str,
    *,
    sheet_name: str | int | None = None,
    symbol_col: str | None = None,
    lfc_col: str | None = None,
    padj_col: str | None = None,
    lfc_thr: float = 0.5,
    padj_thr: float = 0.05,
    universe: set[str] | None = None,
) -> DegSets:
    df = read_upload(data, filename, sheet_name=sheet_name)
    detected = detect_columns(list(df.columns))
    sym = symbol_col or detected["symbol"]
    lfc = lfc_col or detected["lfc"]
    padj = padj_col or detected["padj"]
    if not sym or not lfc or not padj:
        raise ValueError(
            "Could not detect required columns (symbol, log2FC, padj). "
            f"Found columns: {list(df.columns)}"
        )
    return build_deg_sets(
        df,
        symbol_col=sym,
        lfc_col=lfc,
        padj_col=padj,
        lfc_thr=lfc_thr,
        padj_thr=padj_thr,
        universe=universe,
    )
