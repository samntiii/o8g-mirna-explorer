"""Shared Streamlit chrome for the o8G-miRNA Explorer."""
from __future__ import annotations

import os

import streamlit as st

from o8g_db import TargetDB

ROOT = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(ROOT, "o8g_targets.db")
RANK_LABEL = {1: "6mer", 2: "7mer-A1", 3: "7mer-m8", 4: "8mer"}

# Keep CSS scoped to our own classes only — never target Streamlit's [class*="css"]
# or blanket sidebar * rules (those break widgets / navigation).
CSS = """
.o8g-kicker {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 11px;
  letter-spacing: 0.18em;
  text-transform: uppercase;
  color: #b23a3a;
  margin-bottom: 0.35rem;
}
.o8g-hero { max-width: 760px; padding: 0.2rem 0 0.8rem; }
.o8g-hero h1 {
  font-size: 2.2rem;
  line-height: 1.15;
  margin: 0 0 0.6rem;
  color: #1c1915;
}
.o8g-lede {
  font-size: 1.1rem;
  line-height: 1.45;
  color: #3a342c;
  max-width: 42rem;
}
.o8g-rule {
  height: 1px;
  background: linear-gradient(90deg, #b23a3a 0 72px, #d7cec0 72px);
  border: 0;
  margin: 1rem 0 1.2rem;
}
.o8g-cardrow {
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: 0.85rem;
  margin: 0.4rem 0 1.2rem;
}
.o8g-card {
  background: #fffaf2;
  border: 1px solid #d7cec0;
  padding: 0.95rem 1rem;
}
.o8g-card strong {
  display: block;
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 12px;
  letter-spacing: 0.06em;
  text-transform: uppercase;
  color: #b23a3a;
  margin-bottom: 0.35rem;
}
.o8g-card p { margin: 0; color: #3a342c; font-size: 0.98rem; line-height: 1.4; }
.o8g-pair {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 1rem;
  margin: 0.6rem 0 1.2rem;
}
.o8g-pair .cell {
  border: 1px solid #d7cec0;
  background: #fffaf2;
  padding: 1rem 1.05rem;
}
.o8g-pair .cell .lab {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 11px;
  letter-spacing: 0.14em;
  text-transform: uppercase;
  margin-bottom: 0.35rem;
}
.o8g-pair .cell.g .lab { color: #2f4a5c; }
.o8g-pair .cell.ox .lab { color: #b23a3a; }
.o8g-pair .cell .big {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 1.3rem;
  font-weight: 600;
  color: #1c1915;
}
.o8g-pair .cell p { color: #3a342c; margin: 0.4rem 0 0; }
.o8g-footer {
  margin-top: 2rem;
  padding-top: 0.9rem;
  border-top: 1px solid #d7cec0;
  color: #6a6156;
  font-size: 0.88rem;
  line-height: 1.45;
}
.o8g-seed {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 26px;
  letter-spacing: 3px;
  line-height: 1.6;
  color: #1c1915;
}
.o8g-seed .pos {
  font-size: 11px;
  color: #8a8175;
  letter-spacing: 3px;
  margin-top: 2px;
}
.o8g-seed .g {
  background: #dbe4f0;
  color: #2f4a5c;
  padding: 3px 7px;
  font-weight: 700;
}
.o8g-seed .ox {
  background: #b23a3a;
  color: #fffaf2;
  padding: 3px 7px;
  font-weight: 700;
}
.o8g-seed .nt { padding: 3px 5px; color: #1c1915; }
.o8g-seed .end { color: #8a8175; font-size: 13px; letter-spacing: 0; }
@media (max-width: 900px) {
  .o8g-cardrow, .o8g-pair { grid-template-columns: 1fr; }
}
"""


def inject_css():
    st.markdown(f"<style>{CSS}</style>", unsafe_allow_html=True)


@st.cache_resource
def get_db() -> TargetDB:
    return TargetDB(DB_PATH)


def seed_html(seed: str, oxidized: set) -> str:
    cells = []
    for i, b in enumerate(seed):
        pos = i + 2
        if b == "G" and pos in oxidized:
            cells.append(
                f"<span class='ox' title='o8G at position {pos}'>G<sub>o8</sub></span>"
            )
        elif b == "G":
            cells.append(f"<span class='g' title='G at position {pos}'>{b}</span>")
        else:
            cells.append(f"<span class='nt'>{b}</span>")
    pos_row = "".join(f"&nbsp;{p}&nbsp;&nbsp;" for p in range(2, 9))
    return (
        "<div class='o8g-seed'>"
        "<span class='end'>5′-</span>"
        + "".join(cells)
        + "<span class='end'>-3′</span>"
        f"<div class='pos'>&nbsp;&nbsp;&nbsp;&nbsp;{pos_row}</div>"
        "</div>"
    )


def footer():
    st.markdown(
        """<div class="o8g-footer">
        Pairing rule: unmodified seed G pairs C; 8-oxoG (syn) Hoogsteen-pairs A.
        Targets are TargetScan-style 6mer / 7mer-A1 / 7mer-m8 / 8mer sites in the
        longest human 3′UTR per gene (Ensembl). Enrichment is hypergeometric ORA
        with BH correction. Large seed-match lists yield broad enrichment —
        treat pathway shifts as hypotheses and validate candidates experimentally.
        <br>oxomir.samnti.com · o8G-miRNA Retargeting Explorer
        </div>""",
        unsafe_allow_html=True,
    )
