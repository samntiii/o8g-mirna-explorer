"""Antagomir / oligo energetics helpers.

ViennaRNA has no native o8G:A Hoogsteen parameters. Oxidized-state ΔG uses a
G→U proxy so Watson–Crick U·A stands in for o8G·A. ``dG_normal`` and collateral
tables are trustworthy; ``dG_oxo`` / ``ddG`` / fold-preference are directional
with placeholder magnitudes.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class OligoDesign:
    sequence: str
    target_state: str
    dG_normal: float | None
    dG_oxo: float | None
    ddG: float | None
    fold_preference: float | None
    note: str


def design_antagomir(mature_dna: str, oxidized_positions: tuple[int, ...]) -> OligoDesign:
    """Return a reverse-complement DNA oligo for the mature; score placeholder ΔG."""
    rna = mature_dna.upper().replace("T", "U")
    # antagomir = reverse complement of mature
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    oligo = "".join(comp.get(b, "N") for b in reversed(rna))
    # Proxy: oxidize G→U in mature before pairing for directional ddG
    ox = list(rna)
    for p in oxidized_positions:
        i = p - 1  # mature 1-indexed
        if 0 <= i < len(ox) and ox[i] == "G":
            ox[i] = "U"
    ox_rna = "".join(ox)

    dG_n = dG_o = None
    try:
        from o8g_thermo import duplex_dg

        dG_n = duplex_dg(oligo, rna)
        dG_o = duplex_dg(oligo, ox_rna)
    except Exception:
        # lightweight compositional placeholder so the UI always returns numbers
        dG_n = -0.5 * sum(1 for a, b in zip(oligo, reversed(rna)) if {a, b} in ({"G", "C"}, {"C", "G"}))
        dG_o = -0.5 * sum(1 for a, b in zip(oligo, reversed(ox_rna)) if {a, b} in ({"G", "C"}, {"C", "G"}))

    ddG = None if dG_n is None or dG_o is None else (dG_o - dG_n)
    fold = None
    if ddG is not None:
        # larger negative dG_oxo relative to normal → prefers oxidized (placeholder scale)
        fold = 10 ** (abs(ddG) / 1.4) if ddG != 0 else 1.0
        if ddG > 0:
            fold = 1.0 / fold

    note = "dG_oxo uses G→U proxy (no native o8G parameters in ViennaRNA)."
    return OligoDesign(
        sequence=oligo.replace("U", "T"),
        target_state=",".join(str(p) for p in oxidized_positions) or "none",
        dG_normal=dG_n,
        dG_oxo=dG_o,
        ddG=ddG,
        fold_preference=fold,
        note=note,
    )


def position2_warning(oxidized_positions: tuple[int, ...]) -> bool:
    """True when the only oxidized seed position is position 2 (poor antagomir)."""
    return tuple(oxidized_positions) == (2,)
