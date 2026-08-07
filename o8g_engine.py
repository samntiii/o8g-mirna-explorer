"""
o8g_engine.py
=============
Engine for predicting how 8-oxoguanine (o8G) oxidation of guanines in a miRNA
seed redirects its mRNA target repertoire.

Biology
-------
miRNAs recognize targets through Watson-Crick pairing of their *seed* region
(mature-miRNA positions 2-8) to complementary sites in target 3' UTRs.

An unmodified seed guanine (G) pairs a cytosine (C) in the target.
When that guanine is oxidized to 8-oxoguanine, it rotates to the *syn*
conformation and forms a Hoogsteen pair with adenine (o8G . A).  So an
oxidized seed position demands an **A** in the target where the unmodified
seed would have demanded a **C**.

For a seed carrying k guanines there are 2**k oxidation states (each G either
normal or o8G).  Each state defines its own target-site motif and therefore
its own predicted target-gene list.

Coordinate conventions
-----------------------
- Sequences are handled as DNA (U -> T) throughout.
- ``seed`` is the 7-nt string for mature-miRNA positions 2..8 (index 0 == pos 2).
- miRNA position 1 is the nucleotide 5' of the seed; an "A1" target site has an
  A opposite position 1.
- Target motifs are written 5'->3' on the mRNA (antiparallel to the miRNA).

TargetScan-style site types produced for every seed-state:
    8mer     : match to positions 2-8 + A opposite position 1   (revcomp(2-8)+A)
    7mer-m8  : match to positions 2-8                            (revcomp(2-8))
    7mer-A1  : match to positions 2-7 + A opposite position 1    (revcomp(2-7)+A)
    6mer     : match to positions 2-7                            (revcomp(2-7))
"""
from __future__ import annotations
from dataclasses import dataclass, field
from itertools import combinations
from typing import Iterable

# Watson-Crick complement for an *unmodified* base (miRNA base -> required target base)
_WC = {"A": "T", "C": "G", "G": "C", "T": "A"}

SITE_TYPES = ("8mer", "7mer-m8", "7mer-A1", "6mer")


def clean_seq(s: str) -> str:
    """Uppercase and RNA->DNA (U->T)."""
    return s.upper().replace("U", "T")


def extract_seed(mature_seq: str) -> str:
    """Seed = mature-miRNA positions 2-8 (7 nt).  Input may be RNA or DNA."""
    s = clean_seq(mature_seq)
    return s[1:8]


def g_positions(seed: str) -> list[int]:
    """miRNA positions (2-8, 1-based) of guanines in the seed, 5'->3'."""
    return [i + 2 for i, b in enumerate(seed) if b == "G"]


def _target_base(mirna_base: str, oxidized: bool) -> str:
    """Required target base opposite one seed nucleotide."""
    if oxidized:
        # only a guanine can be oxidized to o8G; o8G(syn) . A
        if mirna_base != "G":
            raise ValueError("only G can be oxidized")
        return "A"
    return _WC[mirna_base]


def state_label(oxidized_positions: Iterable[int]) -> str:
    """Human-readable oxidation-state label, e.g. 'none' or 'o8G@2,7'."""
    ox = sorted(oxidized_positions)
    return "none" if not ox else "o8G@" + ",".join(str(p) for p in ox)


@dataclass
class SeedState:
    """One oxidation state of a seed and the target motifs it recognizes."""
    seed: str                       # positions 2-8, DNA
    oxidized_positions: tuple[int, ...]   # miRNA positions (2-8) that are o8G
    label: str = field(init=False)
    motifs: dict[str, str] = field(init=False)   # site_type -> target motif (5'->3')

    def __post_init__(self):
        self.label = state_label(self.oxidized_positions)
        oxset = set(self.oxidized_positions)
        # per-position required target base, indexed by miRNA position 2..8
        req = [_target_base(b, (i + 2) in oxset) for i, b in enumerate(self.seed)]
        # target motif 5'->3' = reverse of complement list (antiparallel)
        m8_7 = "".join(reversed(req))            # positions 2-8  -> 7mer-m8
        core6 = "".join(reversed(req[:6]))       # positions 2-7  -> 6mer
        self.motifs = {
            "8mer": m8_7 + "A",
            "7mer-m8": m8_7,
            "7mer-A1": core6 + "A",
            "6mer": core6,
        }


def enumerate_states(seed: str, max_oxidized: int | None = None) -> list[SeedState]:
    """All 2**k oxidation states for a seed (k = #guanines).

    max_oxidized: optionally cap the number of simultaneously oxidized Gs
    (limits combinatorial blow-up for G-rich seeds; None = no cap).
    """
    seed = clean_seq(seed)
    gpos = g_positions(seed)
    k = len(gpos)
    hi = k if max_oxidized is None else min(max_oxidized, k)
    states: list[SeedState] = []
    for r in range(0, hi + 1):
        for combo in combinations(gpos, r):
            states.append(SeedState(seed=seed, oxidized_positions=tuple(combo)))
    return states


def enumerate_from_mature(mature_seq: str, max_oxidized: int | None = None) -> list[SeedState]:
    return enumerate_states(extract_seed(mature_seq), max_oxidized=max_oxidized)


# ---------------------------------------------------------------- self-test
if __name__ == "__main__":
    # miR-1-3p: TGGAATGTAAAG...  seed(2-8)=GGAATGT ; canonical 7mer-m8 = ACATTCC
    seed = extract_seed("UGGAAUGUAAAGAAGUAUGUAU")
    assert seed == "GGAATGT", seed
    assert g_positions(seed) == [2, 3, 7], g_positions(seed)
    unmod = SeedState(seed, ())
    assert unmod.motifs["7mer-m8"] == "ACATTCC", unmod.motifs
    assert unmod.motifs["8mer"] == "ACATTCCA"
    assert unmod.motifs["6mer"] == "CATTCC"
    # oxidize position 7 (the 7o8G-miR-1 species): target C->A at that position
    ox7 = SeedState(seed, (7,))
    assert ox7.motifs["7mer-m8"] == "AAATTCC", ox7.motifs
    # miR-124-3p: seed AAGGCAC, Gs at positions 4,5
    s124 = extract_seed("UAAGGCACGCGGUGAAUGCCAA")
    assert s124 == "AAGGCAC"
    assert g_positions(s124) == [4, 5]
    states124 = enumerate_states(s124)
    assert len(states124) == 4, len(states124)   # 2^2
    assert {st.label for st in states124} == {"none", "o8G@4", "o8G@5", "o8G@4,5"}
    print("all engine self-tests passed")
    print("miR-1 states:", len(enumerate_states(seed)))
    for st in enumerate_states(seed):
        print(f"  {st.label:12s} 7mer-m8={st.motifs['7mer-m8']}  8mer={st.motifs['8mer']}")
