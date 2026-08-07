#!/usr/bin/env python3
"""Fetch human 3' UTRs from Ensembl BioMart; keep the longest UTR per gene."""
from __future__ import annotations

import argparse
import io
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BIOMART = "https://www.ensembl.org/biomart/martservice"

# Ensembl gene + 3' UTR sequence attributes
XML_QUERY = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="transcript_biotype" value="protein_coding"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="ensembl_transcript_id"/>
    <Attribute name="3utr"/>
  </Dataset>
</Query>
"""


def fetch_biomart(xml: str, retries: int = 6) -> str:
    payload = urllib.parse.urlencode({"query": xml}).encode()
    req = urllib.request.Request(
        BIOMART,
        data=payload,
        headers={"User-Agent": "o8g-mirna-explorer/1.0 (samnti.com)"},
    )
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            print(f"BioMart request attempt {attempt}/{retries}…", flush=True)
            with urllib.request.urlopen(req, timeout=1800) as resp:
                chunks = []
                while True:
                    block = resp.read(1024 * 1024)
                    if not block:
                        break
                    chunks.append(block)
                    if len(chunks) % 20 == 0:
                        n = sum(len(c) for c in chunks)
                        print(f"  downloaded {n/1e6:.1f} MB…", flush=True)
            text = b"".join(chunks).decode("utf-8", errors="replace")
            if text.startswith("Query ERROR") or "ERROR" in text[:200]:
                raise RuntimeError(text[:500])
            return text
        except Exception as exc:
            last_err = exc
            wait = min(60, 5 * attempt)
            print(f"  failed: {exc}\n  retry in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"BioMart fetch failed: {last_err}")


def longest_utr_per_gene(tsv_text: str) -> pd.DataFrame:
    df = pd.read_csv(io.StringIO(tsv_text), sep="\t")
    # BioMart header names vary slightly by release
    colmap = {}
    for c in df.columns:
        cl = c.lower().replace(" ", "_")
        if "ensembl_gene" in cl or cl == "geneid" or cl.endswith("gene_id"):
            colmap[c] = "gene_id"
        elif "external_gene" in cl or cl in {"hgnc_symbol", "gene_name", "genesymbol"}:
            colmap[c] = "symbol"
        elif "transcript" in cl:
            colmap[c] = "transcript_id"
        elif "3utr" in cl or "utr" in cl:
            colmap[c] = "utr3"
    df = df.rename(columns=colmap)
    need = {"gene_id", "symbol", "utr3"}
    missing = need - set(df.columns)
    if missing:
        raise RuntimeError(f"unexpected BioMart columns {list(df.columns)}; missing {missing}")

    df["utr3"] = df["utr3"].fillna("").astype(str).str.upper().str.replace("U", "T", regex=False)
    df = df[~df["utr3"].isin(["", "SEQUENCEUNAVAILABLE", "N/A", "NONE"])]
    df = df[df["utr3"].str.len() >= 20]
    df = df[~df["utr3"].str.contains(r"[^ACGTN]", regex=True)]
    df["utr_len"] = df["utr3"].str.len()
    df["symbol"] = df["symbol"].fillna("").replace("", pd.NA)
    df["symbol"] = df["symbol"].fillna(df["gene_id"])
    df = df.sort_values("utr_len", ascending=False).drop_duplicates("gene_id", keep="first")
    return df[["gene_id", "symbol", "utr3"]].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "utr3_human.parquet")
    ap.add_argument("--raw", type=Path, default=ROOT / "data_raw" / "ensembl_3utr.tsv")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.out.exists() and not args.force:
        df = pd.read_parquet(args.out)
        print(f"{args.out} exists ({len(df)} genes); use --force to refetch")
        return

    args.raw.parent.mkdir(parents=True, exist_ok=True)
    if args.raw.exists() and not args.force:
        text = args.raw.read_text()
    else:
        text = fetch_biomart(XML_QUERY)
        args.raw.write_text(text)
        print(f"wrote raw TSV {args.raw} ({len(text)/1e6:.1f} MB)", flush=True)

    df = longest_utr_per_gene(text)
    df.to_parquet(args.out, index=False)
    print(f"wrote {args.out}: {len(df)} genes with usable 3′UTRs")


if __name__ == "__main__":
    main()
