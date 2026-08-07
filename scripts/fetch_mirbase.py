#!/usr/bin/env python3
"""Download miRBase mature.fa and write hsa_mature.parquet."""
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from o8g_engine import clean_seq, extract_seed, g_positions

# miRBase occasionally renames well-known matures (same MIMAT / sequence).
# Keep the literature names searchable in the explorer.
NAME_ALIASES = {
    "hsa-miR-1-2-3p": ("hsa-miR-1-3p", "hsa-miR-1"),
}

MIRBASE_URLS = [
    "https://www.mirbase.org/download/mature.fa",
    "https://mirbase.org/ftp/CURRENT/mature.fa.gz",
    "https://www.mirbase.org/ftp/CURRENT/mature.fa.gz",
]


def _read_fasta_text(text: str) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    header, seq_parts = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header, seq_parts = line[1:].strip(), []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def parse_hsa_mature(fasta_text: str) -> pd.DataFrame:
    rows = []
    for header, seq in _read_fasta_text(fasta_text):
        # >hsa-miR-1-3p MIMAT0000416 Homo sapiens miR-1-3p
        name = header.split()[0]
        if not name.startswith("hsa-"):
            continue
        accession = ""
        for tok in header.split():
            if tok.startswith("MIMAT"):
                accession = tok
                break
        seq_dna = clean_seq(seq)
        if len(seq_dna) < 8:
            continue
        seed = extract_seed(seq_dna)
        rec = {
            "mirna": name,
            "accession": accession,
            "seq_dna": seq_dna,
            "seed": seed,
            "n_G": len(g_positions(seed)),
        }
        rows.append(rec)
        for alias in NAME_ALIASES.get(name, ()):
            alt = dict(rec)
            alt["mirna"] = alias
            rows.append(alt)
    return pd.DataFrame(rows).drop_duplicates("mirna").sort_values("mirna").reset_index(drop=True)


def download_mature_fa(dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    last_err = None
    for url in MIRBASE_URLS:
        try:
            print(f"fetching {url}", flush=True)
            with urllib.request.urlopen(url, timeout=120) as resp:
                data = resp.read()
            if url.endswith(".gz") or data[:2] == b"\x1f\x8b":
                import gzip

                text = gzip.decompress(data).decode("utf-8", errors="replace")
            else:
                text = data.decode("utf-8", errors="replace")
            dest.write_text(text)
            return text
        except Exception as exc:
            last_err = exc
            print(f"  failed: {exc}", flush=True)
    raise RuntimeError(f"could not download mature.fa: {last_err}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "hsa_mature.parquet")
    ap.add_argument("--fasta", type=Path, default=ROOT / "data_raw" / "mature.fa")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.out.exists() and not args.force:
        df = pd.read_parquet(args.out)
        print(f"{args.out} exists ({len(df)} miRNAs); use --force to refetch")
        return

    if args.fasta.exists() and not args.force:
        text = args.fasta.read_text()
    else:
        text = download_mature_fa(args.fasta)

    df = parse_hsa_mature(text)
    df.to_parquet(args.out, index=False)
    print(
        f"wrote {args.out}: {len(df)} human mature miRNAs, "
        f"{df['seed'].nunique()} unique seeds"
    )


if __name__ == "__main__":
    main()
