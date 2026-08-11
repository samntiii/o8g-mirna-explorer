"""OBOE-inspired ranking of seed guanines likely to oxidize.

Live OBOE (Xia et al.; http://www.rnamd.org/o8GPredictor/) is a remote
RNABERT site predictor. Their public submit.php currently proxies to a local
FastAPI that is often down, so we do **not** hard-depend on it.

This module provides a transparent local prior grounded in the OBOE paper's
reported GC-rich / G-run sequence context preference, for ranking which seed
G positions (and therefore which ox states) to prioritize. Scores are
correlative priors — not causal o8G calls and not a substitute for OBOE.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from o8g_engine import clean_seq, extract_seed, g_positions, enumerate_states


@dataclass(frozen=True)
class GSiteScore:
    position: int  # mature 1-based
    base: str
    window: str
    score: float
    gc_frac: float
    gg_bonus: float


def _window(seq: str, pos1: int, flank: int = 3) -> str:
    """1-based position window on mature sequence (DNA alphabet)."""
    i = pos1 - 1
    a = max(0, i - flank)
    b = min(len(seq), i + flank + 1)
    return seq[a:b]


def score_g_sites(mature_seq: str, *, flank: int = 3) -> list[GSiteScore]:
    """Score every guanine on the mature (focus: seed Gs used for ox states)."""
    seq = clean_seq(mature_seq)
    out: list[GSiteScore] = []
    for i, b in enumerate(seq):
        if b != "G":
            continue
        pos = i + 1
        w = _window(seq, pos, flank=flank)
        gc = (w.count("G") + w.count("C")) / max(len(w), 1)
        # local G-run / GG preference highlighted by OBOE motif analyses
        left = seq[i - 1] if i > 0 else ""
        right = seq[i + 1] if i + 1 < len(seq) else ""
        gg = 0.20 * ((left == "G") + (right == "G"))
        # seed 3′ Gs (esp. 6–8) are repeatedly reported oxidized in wet-lab
        # oxomiR studies; blend lightly so AT-flanked pos7 is not zeroed out
        if 6 <= pos <= 8:
            lit = 0.22
        elif 2 <= pos <= 5:
            lit = 0.08
        else:
            lit = 0.0
        raw = 0.45 * gc + gg + lit
        score = max(0.0, min(1.0, raw))
        out.append(
            GSiteScore(
                position=pos,
                base=b,
                window=w,
                score=round(score, 4),
                gc_frac=round(gc, 4),
                gg_bonus=round(gg, 4),
            )
        )
    return out


def seed_g_table(mature_seq: str) -> pd.DataFrame:
    seed = extract_seed(mature_seq)
    seed_gs = set(g_positions(seed))
    rows = []
    for s in score_g_sites(mature_seq):
        rows.append(
            dict(
                mature_pos=s.position,
                in_seed=s.position in seed_gs,
                window=s.window,
                gc_frac=s.gc_frac,
                gg_bonus=s.gg_bonus,
                oboe_prior=s.score,
            )
        )
    return pd.DataFrame(rows)


def rank_oxidation_states(mature_seq: str) -> pd.DataFrame:
    """Rank non-'none' ox states by mean OBOE-style prior of oxidized seed Gs."""
    seed = extract_seed(mature_seq)
    priors = {s.position: s.score for s in score_g_sites(mature_seq)}
    rows = []
    for st in enumerate_states(seed):
        if st.label == "none":
            continue
        scores = [priors.get(p, 0.0) for p in st.oxidized_positions]
        mean_s = sum(scores) / len(scores) if scores else 0.0
        # prefer fewer oxidations when means tie (parsimony)
        rows.append(
            dict(
                ox_label=st.label,
                n_ox=len(st.oxidized_positions),
                positions=",".join(str(p) for p in st.oxidized_positions),
                mean_oboe_prior=round(mean_s, 4),
                min_oboe_prior=round(min(scores), 4) if scores else 0.0,
            )
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(
        ["mean_oboe_prior", "n_ox"], ascending=[False, True]
    ).reset_index(drop=True)


def try_remote_oboe(sequence: str, *, timeout: float = 8.0) -> dict | None:
    """Best-effort call to the public OBOE submit.php. Returns None on failure."""
    import json
    import urllib.error
    import urllib.request

    seq = clean_seq(sequence).replace("T", "U")
    boundary = "----O8GBoundary7"
    body = (
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="sequence"\r\n\r\n'
        f"{seq}\r\n"
        f"--{boundary}--\r\n"
    ).encode()
    req = urllib.request.Request(
        "http://www.rnamd.org/o8GPredictor/submit.php",
        data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read().decode("utf-8", "replace")
        # strip PHP warnings if present
        start = raw.find("{")
        if start < 0:
            return None
        return json.loads(raw[start:])
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, ValueError):
        return None
