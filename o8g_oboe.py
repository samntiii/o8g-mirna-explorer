"""OBOE-based ranking of seed guanines likely to oxidize.

Primary scorer: local fine-tuned RNABERT from Xia et al. OBOE
(Figshare DOI 10.6084/m9.figshare.29634239; see ``o8g_oboe_model`` and
``scripts/train_oboe_rnabert.py``). Falls back to a transparent GC-rich / G-run
motif prior if the checkpoint or torch stack is unavailable.

Scores are **site likelihood priors** for ranking ox states — not causal o8G
calls and not a substitute for wet-lab validation.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from o8g_engine import clean_seq, extract_seed, g_positions, enumerate_states

_MODEL_OK: bool | None = None
_MODEL_ERR: str | None = None


def model_status() -> tuple[bool, str]:
    """Return (available, message) for the local OBOE RNABERT checkpoint."""
    global _MODEL_OK, _MODEL_ERR
    if _MODEL_OK is not None:
        return bool(_MODEL_OK), _MODEL_ERR or ("OBOE RNABERT ready" if _MODEL_OK else "unavailable")
    try:
        from o8g_oboe_model import checkpoint_available, predict_o8g_prob

        if not checkpoint_available():
            _MODEL_OK = False
            _MODEL_ERR = "Checkpoint missing (run scripts/train_oboe_rnabert.py)"
            return False, _MODEL_ERR
        # smoke one short call so import/load failures surface once
        predict_o8g_prob("G" * 21)
        _MODEL_OK = True
        _MODEL_ERR = "Local OBOE RNABERT (fine-tuned on Xia et al. Figshare data)"
        return True, _MODEL_ERR
    except Exception as e:  # noqa: BLE001 — surface any optional-dep failure
        _MODEL_OK = False
        _MODEL_ERR = f"OBOE model unavailable ({e})"
        return False, _MODEL_ERR


@dataclass(frozen=True)
class GSiteScore:
    position: int  # mature 1-based
    base: str
    window: str
    score: float
    gc_frac: float
    gg_bonus: float
    source: str  # "oboe_rnabert" | "gc_prior"


def _window(seq: str, pos1: int, flank: int = 3) -> str:
    """1-based position window on mature sequence (DNA alphabet)."""
    i = pos1 - 1
    a = max(0, i - flank)
    b = min(len(seq), i + flank + 1)
    return seq[a:b]


def _gc_prior_score(seq: str, pos: int, flank: int = 3) -> tuple[float, float, float, str]:
    i = pos - 1
    w = _window(seq, pos, flank=flank)
    gc = (w.count("G") + w.count("C")) / max(len(w), 1)
    left = seq[i - 1] if i > 0 else ""
    right = seq[i + 1] if i + 1 < len(seq) else ""
    gg = 0.20 * ((left == "G") + (right == "G"))
    if 6 <= pos <= 8:
        lit = 0.22
    elif 2 <= pos <= 5:
        lit = 0.08
    else:
        lit = 0.0
    raw = 0.45 * gc + gg + lit
    score = max(0.0, min(1.0, raw))
    return score, gc, gg, w


def score_g_sites(mature_seq: str, *, flank: int = 3, prefer_model: bool = True) -> list[GSiteScore]:
    """Score every guanine on the mature (focus: seed Gs used for ox states).

    When the OBOE RNABERT checkpoint is available, combine model P(o8G|local
    window) with the GC/G-run motif prior so short-mature ranking stays
    discriminative (model alone is strong across sequences, flatter within one
    ~22 nt mature).
    """
    seq = clean_seq(mature_seq)
    g_pos = [i + 1 for i, b in enumerate(seq) if b == "G"]
    model_probs: dict[int, float] = {}
    if prefer_model and g_pos:
        ok, _ = model_status()
        if ok:
            try:
                from o8g_oboe_model import score_g_positions

                model_probs = score_g_positions(seq, g_pos)
            except Exception:
                model_probs = {}

    out: list[GSiteScore] = []
    for pos in g_pos:
        gc_score, gc, gg, w_short = _gc_prior_score(seq, pos, flank=flank)
        if pos in model_probs:
            try:
                from o8g_oboe_model import window_around

                w = window_around(seq, pos)
            except Exception:
                w = w_short
            # 70% OBOE RNABERT + 30% motif prior (within-mature contrast)
            score = 0.70 * float(model_probs[pos]) + 0.30 * float(gc_score)
            src = "oboe_rnabert"
        else:
            w = w_short
            score = gc_score
            src = "gc_prior"
        out.append(
            GSiteScore(
                position=pos,
                base="G",
                window=w,
                score=round(float(score), 4),
                gc_frac=round(gc, 4),
                gg_bonus=round(gg, 4),
                source=src,
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
                source=s.source,
            )
        )
    return pd.DataFrame(rows)


def rank_oxidation_states(mature_seq: str) -> pd.DataFrame:
    """Rank non-'none' ox states by mean OBOE prior of oxidized seed Gs."""
    seed = extract_seed(mature_seq)
    scored = score_g_sites(mature_seq)
    priors = {s.position: s.score for s in scored}
    src = scored[0].source if scored else "gc_prior"
    rows = []
    for st in enumerate_states(seed):
        if st.label == "none":
            continue
        scores = [priors.get(p, 0.0) for p in st.oxidized_positions]
        mean_s = sum(scores) / len(scores) if scores else 0.0
        rows.append(
            dict(
                ox_label=st.label,
                n_ox=len(st.oxidized_positions),
                positions=",".join(str(p) for p in st.oxidized_positions),
                mean_oboe_prior=round(mean_s, 4),
                min_oboe_prior=round(min(scores), 4) if scores else 0.0,
                source=src,
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
        start = raw.find("{")
        if start < 0:
            return None
        return json.loads(raw[start:])
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, ValueError):
        return None
