#!/usr/bin/env python3
"""Download NCBI (+ optional HGNC) gene ID maps and join onto o8g_targets.db genes.

Writes gene_aliases.parquet (long form: one row per alias string) and
gene_id_map.parquet (one row per gene_idx with canonical IDs).
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import io
import sqlite3
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw"

NCBI_GENE_INFO = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
)
NCBI_GENE2ENSEMBL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)


def download(url: str, dest: Path, force: bool = False) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        print(f"using cached {dest}", flush=True)
        return dest
    print(f"fetching {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "o8g-mirna-explorer/1.0"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        data = resp.read()
    dest.write_bytes(data)
    print(f"  wrote {dest} ({len(data)/1e6:.1f} MB)", flush=True)
    return dest


def _open_text(path: Path):
    if path.suffix == ".gz" or path.name.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, encoding="utf-8")


def load_gene2ensembl(path: Path) -> pd.DataFrame:
    """tax_id, GeneID, Ensembl_gene_identifier, ... — keep human Ensembl gene IDs."""
    rows = []
    with _open_text(path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            tax, gene_id, ens = parts[0], parts[1], parts[2]
            if tax != "9606":
                continue
            if not ens.startswith("ENSG"):
                continue
            # strip version if present
            ens = ens.split(".")[0]
            rows.append((int(gene_id), ens))
    df = pd.DataFrame(rows, columns=["entrez", "ensembl"])
    return df.drop_duplicates()


def load_gene_info(path: Path) -> pd.DataFrame:
    """NCBI Homo_sapiens.gene_info — Symbol + Synonyms."""
    # columns: tax_id GeneID Symbol LocusTag Synonyms ...
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip" if str(path).endswith(".gz") else None,
        dtype=str,
        low_memory=False,
    )
    df = df[df["#tax_id"] == "9606"] if "#tax_id" in df.columns else df[df["tax_id"] == "9606"]
    gid_col = "GeneID"
    sym_col = "Symbol"
    syn_col = "Synonyms"
    out = df[[gid_col, sym_col, syn_col]].copy()
    out.columns = ["entrez", "symbol", "synonyms"]
    out["entrez"] = out["entrez"].astype(int)
    return out


def load_hgnc(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    # ensembl_gene_id, entrez_id, symbol, alias_symbol, prev_symbol
    keep = [c for c in ["ensembl_gene_id", "entrez_id", "symbol", "alias_symbol", "prev_symbol"] if c in df.columns]
    return df[keep].copy()


def build_maps(db_path: Path, force: bool = False, skip_hgnc: bool = False):
    gene_info_path = download(NCBI_GENE_INFO, RAW / "Homo_sapiens.gene_info.gz", force=force)
    g2e_path = download(NCBI_GENE2ENSEMBL, RAW / "gene2ensembl.gz", force=force)
    hgnc_path = None
    if not skip_hgnc:
        try:
            hgnc_path = download(HGNC_URL, RAW / "hgnc_complete_set.txt", force=force)
        except Exception as exc:
            print(f"HGNC download failed ({exc}); continuing with NCBI only", flush=True)

    g2e = load_gene2ensembl(g2e_path)
    info = load_gene_info(gene_info_path)
    # prefer one entrez per ensembl (first)
    g2e = g2e.drop_duplicates("ensembl", keep="first")

    con = sqlite3.connect(str(db_path))
    genes = pd.read_sql("SELECT gene_idx, gene_id AS ensembl, symbol AS db_symbol FROM genes", con)
    con.close()
    genes["ensembl"] = genes["ensembl"].astype(str).str.split(".").str[0]

    m = genes.merge(g2e, on="ensembl", how="left")
    m = m.merge(info[["entrez", "symbol", "synonyms"]], on="entrez", how="left")
    # canonical symbol: NCBI symbol if present else DB symbol
    m["symbol"] = m["symbol"].fillna(m["db_symbol"])

    if hgnc_path is not None and hgnc_path.exists():
        h = load_hgnc(hgnc_path)
        h = h.rename(columns={"ensembl_gene_id": "ensembl", "entrez_id": "entrez_hgnc", "symbol": "hgnc_symbol"})
        h["ensembl"] = h["ensembl"].fillna("").astype(str).str.split(".").str[0]
        h = h[h["ensembl"] != ""]
        h = h.drop_duplicates("ensembl", keep="first")
        m = m.merge(h, on="ensembl", how="left")
        # fill missing entrez from HGNC
        miss = m["entrez"].isna() & m["entrez_hgnc"].notna() & (m["entrez_hgnc"] != "")
        m.loc[miss, "entrez"] = pd.to_numeric(m.loc[miss, "entrez_hgnc"], errors="coerce")
        # prefer HGNC symbol when available
        has_h = m["hgnc_symbol"].notna() & (m["hgnc_symbol"] != "")
        m.loc[has_h, "symbol"] = m.loc[has_h, "hgnc_symbol"]
    else:
        m["alias_symbol"] = None
        m["prev_symbol"] = None

    # long-form aliases
    rows = []
    for r in m.itertuples(index=False):
        gene_idx = int(r.gene_idx)
        ens = r.ensembl
        entrez = None if pd.isna(r.entrez) else int(r.entrez)
        sym = str(r.symbol) if pd.notna(r.symbol) else str(r.db_symbol)
        aliases = {sym, ens}
        if entrez is not None:
            aliases.add(str(entrez))
        # NCBI synonyms
        syn = getattr(r, "synonyms", None)
        if syn is not None and pd.notna(syn) and str(syn) not in ("-", ""):
            for a in str(syn).split("|"):
                a = a.strip()
                if a and a != "-":
                    aliases.add(a)
        # HGNC aliases / previous
        for field in ("alias_symbol", "prev_symbol"):
            val = getattr(r, field, None)
            if val is not None and pd.notna(val) and str(val) not in ("-", ""):
                for a in str(val).replace("|", ",").split(","):
                    a = a.strip()
                    if a and a != "-":
                        aliases.add(a)
        for a in aliases:
            rows.append(
                {
                    "gene_idx": gene_idx,
                    "ensembl": ens,
                    "entrez": entrez,
                    "symbol": sym,
                    "alias": a,
                    "alias_norm": a.upper(),
                }
            )

    aliases_df = pd.DataFrame(rows).drop_duplicates(["gene_idx", "alias_norm"])
    id_map = m[
        ["gene_idx", "ensembl", "entrez", "symbol", "db_symbol"]
    ].drop_duplicates("gene_idx")
    id_map["entrez"] = pd.to_numeric(id_map["entrez"], errors="coerce").astype("Int64")

    meta = {
        "built_at": dt.date.today().isoformat(),
        "n_genes": int(len(id_map)),
        "n_alias_rows": int(len(aliases_df)),
        "n_with_entrez": int(id_map["entrez"].notna().sum()),
        "sources": "NCBI gene_info + gene2ensembl"
        + (" + HGNC" if hgnc_path and hgnc_path.exists() else ""),
    }
    alias_out = ROOT / "gene_aliases.parquet"
    map_out = ROOT / "gene_id_map.parquet"
    aliases_df.to_parquet(alias_out, index=False)
    id_map.to_parquet(map_out, index=False)
    (ROOT / "gene_aliases.meta.txt").write_text(
        "\n".join(f"{k}={v}" for k, v in meta.items()) + "\n"
    )
    print(f"wrote {alias_out} ({len(aliases_df)} alias rows)")
    print(f"wrote {map_out} ({len(id_map)} genes, {meta['n_with_entrez']} with Entrez)")
    print(meta)
    return alias_out, map_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path, default=ROOT / "o8g_targets.db")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--skip-hgnc", action="store_true")
    args = ap.parse_args()
    build_maps(args.db, force=args.force, skip_hgnc=args.skip_hgnc)


if __name__ == "__main__":
    main()
