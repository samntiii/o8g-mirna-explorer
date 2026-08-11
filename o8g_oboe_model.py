"""Local OBOE RNABERT inference (Xia et al. Figshare 29634239).

Figshare ships training code + labeled windows, not weights. We fine-tune
``multimolecule/rnabert`` with ``scripts/train_oboe_rnabert.py`` and load the
checkpoint from ``models/oboe_rnabert/checkpoint``.

The published model is a **window-level** o8G classifier. For seed ranking on
short matures we score a local (~15 nt) window expanded asymmetrically into
real bases around each guanine (full-mature centering collapses all Gs to one
score).
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_CKPT = ROOT / "models" / "oboe_rnabert" / "checkpoint"
LOCAL_WINDOW = 15
MIN_WINDOW = 12  # shortest examples in OBOE train_0.9.csv


def checkpoint_available(path: Path | None = None) -> bool:
    p = Path(path) if path else DEFAULT_CKPT
    return (p / "model.safetensors").exists() or (p / "pytorch_model.bin").exists()


def _to_rna(seq: str) -> str:
    return "".join(c if c in "ACGU" else "N" for c in seq.upper().replace("T", "U"))


def window_around(seq_dna: str, pos1: int, width: int = LOCAL_WINDOW) -> str:
    """Symmetric clip around 1-based ``pos1`` (no asymmetric fill).

    Asymmetric expansion made every 5′ seed G share the same first 15 nt on a
    ~22 nt mature. Shorter edge windows are fine — OBOE training includes
    sequences down to 12 nt.
    """
    seq = _to_rna(seq_dna)
    n = len(seq)
    if n == 0:
        return "N" * MIN_WINDOW
    i = max(0, min(n - 1, pos1 - 1))
    half = width // 2
    start = max(0, i - half)
    end = min(n, i + half + 1)
    chunk = seq[start:end]
    # Only pad if the *entire mature* is shorter than MIN_WINDOW
    if n < MIN_WINDOW and len(chunk) < MIN_WINDOW:
        need = MIN_WINDOW - len(chunk)
        left = need // 2
        chunk = ("N" * left) + chunk + ("N" * (need - left))
    return chunk


@lru_cache(maxsize=1)
def _load_bundle(ckpt: str):
    import multimolecule  # noqa: F401
    import torch
    from transformers import AutoModelForSequenceClassification, AutoTokenizer

    path = Path(ckpt)
    tok = AutoTokenizer.from_pretrained(path, bos_token=None, eos_token=None)
    model = AutoModelForSequenceClassification.from_pretrained(path, num_labels=2)
    model.eval()
    device = "cpu"
    if torch.backends.mps.is_available():
        device = "mps"
    elif torch.cuda.is_available():
        device = "cuda"
    model.to(device)
    return tok, model, device


def predict_o8g_prob(sequence: str, *, ckpt: Path | None = None, max_length: int = 256) -> float:
    """P(o8G | sequence) for one window (OBOE inference.py semantics)."""
    return predict_o8g_probs([sequence], ckpt=ckpt, max_length=max_length)[0]


def predict_o8g_probs(
    sequences: list[str],
    *,
    ckpt: Path | None = None,
    max_length: int = 256,
    batch_size: int = 32,
) -> list[float]:
    if not sequences:
        return []
    path = str(Path(ckpt) if ckpt else DEFAULT_CKPT)
    tok, model, device = _load_bundle(path)
    import torch

    probs: list[float] = []
    with torch.no_grad():
        for i0 in range(0, len(sequences), batch_size):
            batch = [_to_rna(s) for s in sequences[i0 : i0 + batch_size]]
            inputs = tok(
                batch,
                return_tensors="pt",
                padding="max_length",
                truncation=True,
                max_length=max_length,
            )
            inputs = {k: v.to(device) for k, v in inputs.items()}
            logits = model(**inputs).logits
            p = torch.softmax(logits, dim=1)[:, 1].detach().cpu().numpy()
            probs.extend(float(x) for x in p)
    return probs


def score_g_positions(mature_dna: str, positions: list[int] | None = None) -> dict[int, float]:
    """Map 1-based mature G positions → OBOE P(o8G | local window)."""
    seq = mature_dna.upper().replace("U", "T")
    if positions is None:
        positions = [i + 1 for i, b in enumerate(seq) if b == "G"]
    if not positions:
        return {}
    windows = [window_around(seq, p, LOCAL_WINDOW) for p in positions]
    local = predict_o8g_probs(windows)
    return {p: round(float(lp), 4) for p, lp in zip(positions, local)}
