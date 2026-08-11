"""Fetch Enrichr TF / motif gene-set libraries into genesets/.

Usage:
  python scripts/fetch_tf_genesets.py
"""
from __future__ import annotations

import os
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "genesets")

# Motif + ChIP TF libraries used by the Transcription factors view.
LIBRARIES = [
    "TRANSFAC_and_JASPAR_PWMs",
    "JASPAR_PWM_Human_2025",
    "Genome_Browser_PWMs",  # UCSC/ENCODE-style PWM tracks (Ensembl Reg. Build–adjacent)
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
    "ChEA_2022",
    "ENCODE_TF_ChIP-seq_2015",
    "TRRUST_Transcription_Factors_2019",
    "TF_Perturbations_Followed_by_Expression",
]


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    for name in LIBRARIES:
        dest = os.path.join(OUT, f"{name}.gmt")
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            print(f"skip {name} ({os.path.getsize(dest)} bytes)")
            continue
        url = (
            "https://maayanlab.cloud/Enrichr/geneSetLibrary"
            f"?mode=text&libraryName={name}"
        )
        print(f"fetch {name} …")
        urllib.request.urlretrieve(url, dest)
        print(f"  -> {dest} ({os.path.getsize(dest)} bytes)")


if __name__ == "__main__":
    main()
