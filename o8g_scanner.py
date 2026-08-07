"""
o8g_scanner.py
==============
Fast 3' UTR seed-match scanner for the o8G retargeting engine.

Design
------
TargetScan-style site classification for a seed needs only:
  * the 6mer core  = miRNA positions 2-7  (drawn 5'->3' on the mRNA)
  * the m8 flank   = base 5' of the core in the mRNA (opposite miRNA pos 8)
  * the A1 flank   = base 3' of the core in the mRNA (opposite miRNA pos 1)

So we index EVERY 6mer occurrence across all UTRs exactly once, recording
(gene index, 5'-flank base, is-A1).  Any seed-state query is then a
binary-search over the 6mer code plus a vectorised classification -- fast
(ms) and independent of the number of seed-states.

Site precedence (best kept per gene): 8mer > 7mer-m8 > 7mer-A1 > 6mer.

Base encoding: A=0 C=1 G=2 T=3 ; N / gene-separators = -1 (invalidates any
6mer window that contains them).

Context-style score (lightweight context++ analog)
--------------------------------------------------
Per-site features (best site per gene unless noted):
  - local AU fraction in ±30 nt flanks (Garcia et al. NSMB 2011; Agarwal eLife 2015)
  - 3′ supplementary pairing: miRNA positions 13–16 vs UTR opposite that register
    (Grimson et al. Mol Cell 2007; carried in TargetScan context++)
  - min distance to either UTR end (sites near ends often more accessible)
  - relative position along the UTR (0–1)

Linear combination into ``context_score`` (higher = more favorable). Coefficients
approximate published TargetScan context++ *signs* with simplified magnitudes
(not a byte-for-byte reimplementation of context++):

  CONTEXT_W_AU  = 1.0   # AU-rich flanks favor repression
  CONTEXT_W_3P  = 0.5   # supplementary 3′ pairing
  CONTEXT_W_END = 0.3   # proximity to UTR termini (scaled)
  CONTEXT_W_POS = -0.2  # mid-UTR slightly less favored vs terminal third

SITE_WEIGHT (multiplicity / confident_targets) remains separate from context_score.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

_B2I = {"A": 0, "C": 1, "G": 2, "T": 3}
_I2B = np.array(list("ACGT"))
SITE_RANK = {"8mer": 4, "7mer-m8": 3, "7mer-A1": 2, "6mer": 1}
RANK_SITE = {v: k for k, v in SITE_RANK.items()}

# --- context++ analog coefficients (see module docstring for citations) ---
CONTEXT_FLANK = 30  # nt on each side for local AU
CONTEXT_W_AU = 1.0
CONTEXT_W_3P = 0.5
CONTEXT_W_END = 0.3
CONTEXT_W_POS = -0.2
# Scale min-end distance: contribution saturates by ~CONTEXT_END_SCALE nt
CONTEXT_END_SCALE = 100.0


def _encode(seq: str) -> np.ndarray:
    a = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    out = np.full(a.shape, -1, dtype=np.int8)
    out[a == ord("A")] = 0
    out[a == ord("C")] = 1
    out[a == ord("G")] = 2
    out[a == ord("T")] = 3
    return out


def _code6(seq_i8: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rolling 6mer code (base-4) over an int8 sequence.

    Returns (code, valid) arrays of length len(seq)-5.  code is int32 in
    [0,4095] where valid; validity is False if the window has any base < 0.
    """
    L = seq_i8.shape[0]
    n = L - 5
    code = np.zeros(n, dtype=np.int32)
    valid = np.ones(n, dtype=bool)
    for j in range(6):
        b = seq_i8[j:j + n].astype(np.int32)
        valid &= (b >= 0)
        code += np.where(b >= 0, b, 0) * (4 ** (5 - j))
    return code, valid


def _au_fraction(seq: str, start: int, end: int) -> float:
    """AU fraction in [start, end) clipped to sequence bounds."""
    if end <= start or not seq:
        return 0.0
    start = max(0, start)
    end = min(len(seq), end)
    if end <= start:
        return 0.0
    window = seq[start:end]
    if not window:
        return 0.0
    return (window.count("A") + window.count("T")) / len(window)


def _supp_pair_score(mature_dna: str, utr: str, site_start_6mer: int) -> float:
    """Fraction of WC matches for miRNA pos 13–16 vs UTR 3′ of the seed site.

    Geometry (TargetScan / Grimson): seed match on mRNA; miRNA 3′ region pairs
    upstream on the mRNA relative to the seed. For a 6mer starting at ``i``
    (mRNA 5′→3′), the register opposite miRNA 13–16 is approximated as
    utr[i-6 : i-2] (immediately 5′ of the 6mer). Returns 0..1.
    """
    if not mature_dna or len(mature_dna) < 16:
        return 0.0
    mir = mature_dna[12:16]
    lo = site_start_6mer - 6
    hi = site_start_6mer - 2
    if lo < 0 or hi > len(utr) or hi <= lo:
        return 0.0
    utr_seg = utr[lo:hi][::-1]
    if len(utr_seg) != 4:
        return 0.0
    wc = {"A": "T", "T": "A", "G": "C", "C": "G"}
    hits = sum(1 for a, b in zip(mir, utr_seg) if wc.get(a) == b)
    return hits / 4.0


class TargetScanner:
    """Indexes all UTRs; classifies seed-state target sites per gene."""

    def __init__(self, genes: list[str], symbols: list[str]):
        self.genes = np.asarray(genes)
        self.symbols = np.asarray(symbols)
        self.utrs: list[str] = []
        self._sorted_code = None
        self._gene_idx = None
        self._flank5 = None
        self._a1 = None
        self._pos_in_utr = None

    @classmethod
    def from_parquet(cls, path: str) -> "TargetScanner":
        df = pd.read_parquet(path)
        sc = cls(df["gene_id"].tolist(), df["symbol"].tolist())
        sc.build(df["utr3"].tolist())
        return sc

    def build(self, utrs: list[str]):
        self.utrs = [u.upper().replace("U", "T") for u in utrs]
        SEP = np.array([-1], dtype=np.int8)
        parts, bounds, cum = [], [], 0
        for gi, u in enumerate(self.utrs):
            e = _encode(u)
            parts.append(e)
            parts.append(SEP)
            bounds.append((cum, cum + e.shape[0], gi))
            cum += e.shape[0] + 1
        big = np.concatenate(parts)
        code, valid = _code6(big)
        pos = np.nonzero(valid)[0]
        codes = code[pos]
        flank5 = np.full(pos.shape, -1, dtype=np.int8)
        m = pos >= 1
        flank5[m] = big[pos[m] - 1]
        a1pos = pos + 6
        a1 = np.zeros(pos.shape, dtype=bool)
        m2 = a1pos < big.shape[0]
        a1[m2] = big[a1pos[m2]] == 0
        starts = np.array([b[0] for b in bounds])
        gidx = np.searchsorted(starts, pos, side="right") - 1
        pos_in_utr = pos - starts[gidx]
        order = np.argsort(codes, kind="stable")
        self._sorted_code = codes[order]
        self._gene_idx = gidx[order].astype(np.int32)
        self._flank5 = flank5[order]
        self._a1 = a1[order]
        self._pos_in_utr = pos_in_utr[order].astype(np.int32)
        return self

    @staticmethod
    def _core_code(core6: str) -> int:
        c = 0
        for ch in core6:
            c = c * 4 + _B2I[ch]
        return c

    def scan_state(self, state) -> pd.DataFrame:
        """Best site type per gene for one SeedState."""
        core = state.motifs["6mer"]
        b8 = _B2I[state.motifs["7mer-m8"][0]]
        cc = self._core_code(core)
        lo = np.searchsorted(self._sorted_code, cc, "left")
        hi = np.searchsorted(self._sorted_code, cc, "right")
        if hi == lo:
            return pd.DataFrame(columns=["gene_id", "symbol", "site_type", "site_rank"])
        g = self._gene_idx[lo:hi]
        m8 = self._flank5[lo:hi] == b8
        a1 = self._a1[lo:hi]
        rank = np.where(m8 & a1, 4, np.where(m8 & ~a1, 3, np.where(~m8 & a1, 2, 1))).astype(
            np.int8
        )
        best = np.zeros(self.genes.shape[0], dtype=np.int8)
        np.maximum.at(best, g, rank)
        hit = np.nonzero(best > 0)[0]
        return pd.DataFrame(
            {
                "gene_id": self.genes[hit],
                "symbol": self.symbols[hit],
                "site_rank": best[hit],
                "site_type": [RANK_SITE[r] for r in best[hit]],
            }
        )

    SITE_WEIGHT = {4: 1.0, 3: 0.7, 2: 0.4, 1: 0.15}

    def scan_state_counts(self, state):
        """Per-gene site inventory + SITE_WEIGHT score."""
        core = state.motifs["6mer"]
        b8 = _B2I[state.motifs["7mer-m8"][0]]
        cc = self._core_code(core)
        lo = np.searchsorted(self._sorted_code, cc, "left")
        hi = np.searchsorted(self._sorted_code, cc, "right")
        ng = self.genes.shape[0]
        counts = np.zeros((ng, 4), dtype=np.int32)
        if hi > lo:
            g = self._gene_idx[lo:hi]
            m8 = self._flank5[lo:hi] == b8
            a1 = self._a1[lo:hi]
            rank = np.where(m8 & a1, 4, np.where(m8 & ~a1, 3, np.where(~m8 & a1, 2, 1)))
            for rk, col in ((4, 0), (3, 1), (2, 2), (1, 3)):
                np.add.at(counts[:, col], g[rank == rk], 1)
        n_sites = counts.sum(axis=1)
        hit = np.nonzero(n_sites > 0)[0]
        w = np.array([1.0, 0.7, 0.4, 0.15])
        score = counts[hit] @ w
        best_rank = np.where(
            counts[hit, 0] > 0,
            4,
            np.where(counts[hit, 1] > 0, 3, np.where(counts[hit, 2] > 0, 2, 1)),
        )
        return pd.DataFrame(
            {
                "gene_idx": hit.astype(np.int32),
                "symbol": self.symbols[hit],
                "n_8mer": counts[hit, 0],
                "n_7mer_m8": counts[hit, 1],
                "n_7mer_A1": counts[hit, 2],
                "n_6mer": counts[hit, 3],
                "n_sites": n_sites[hit],
                "best_rank": best_rank,
                "score": score,
            }
        )

    def scan_state_context(self, state, mature_dna: str | None = None, min_rank: int = 1):
        """Best site per gene with multiplicity score + lightweight context_score."""
        counts = self.scan_state_counts(state)
        if counts.empty:
            return counts
        counts = counts[counts["best_rank"] >= min_rank].copy()
        if counts.empty:
            return counts

        core = state.motifs["6mer"]
        b8 = _B2I[state.motifs["7mer-m8"][0]]
        cc = self._core_code(core)
        lo = np.searchsorted(self._sorted_code, cc, "left")
        hi = np.searchsorted(self._sorted_code, cc, "right")
        g = self._gene_idx[lo:hi]
        m8 = self._flank5[lo:hi] == b8
        a1 = self._a1[lo:hi]
        rank = np.where(m8 & a1, 4, np.where(m8 & ~a1, 3, np.where(~m8 & a1, 2, 1)))
        pos = self._pos_in_utr[lo:hi]

        best_pos: dict[int, int] = {}
        best_rk: dict[int, int] = {}
        for gi, rk, pu in zip(g.tolist(), rank.tolist(), pos.tolist()):
            if gi not in best_rk or rk > best_rk[gi]:
                best_rk[gi] = int(rk)
                best_pos[gi] = int(pu)

        mature = (mature_dna or "").upper().replace("U", "T")
        ctx_scores, au_vals, p3_vals, end_vals, pos_vals = [], [], [], [], []
        for gi in counts["gene_idx"].tolist():
            utr = self.utrs[int(gi)]
            L = len(utr)
            pu = int(best_pos.get(int(gi), 0))
            site_end = min(L, pu + 8)
            au = 0.5 * (
                _au_fraction(utr, pu - CONTEXT_FLANK, pu)
                + _au_fraction(utr, site_end, site_end + CONTEXT_FLANK)
            )
            p3 = _supp_pair_score(mature, utr, pu) if mature else 0.0
            min_end = float(min(pu, max(0, L - site_end)))
            end_term = min(1.0, min_end / CONTEXT_END_SCALE)
            rel = (pu / L) if L > 0 else 0.5
            mid_pen = 1.0 - abs(rel - 0.5) * 2.0
            ctx = (
                CONTEXT_W_AU * au
                + CONTEXT_W_3P * p3
                + CONTEXT_W_END * end_term
                + CONTEXT_W_POS * mid_pen
            )
            ctx_scores.append(ctx)
            au_vals.append(au)
            p3_vals.append(p3)
            end_vals.append(min_end)
            pos_vals.append(rel)

        counts["context_score"] = ctx_scores
        counts["au_flank"] = au_vals
        counts["supp3_score"] = p3_vals
        counts["min_end_dist"] = end_vals
        counts["site_position"] = pos_vals
        counts["site_start"] = [int(best_pos.get(int(gi), 0)) for gi in counts["gene_idx"].tolist()]
        counts["gene_id"] = self.genes[counts["gene_idx"].to_numpy()]
        counts["site_rank"] = counts["best_rank"]
        counts["site_type"] = [RANK_SITE[int(r)] for r in counts["best_rank"]]
        return counts

    def confident_targets(self, state, min_score: float = 1.0, min_sites: int = 1):
        """High-confidence rows: score>=min_score AND n_sites>=min_sites.

        Default min_score=1.0 ~ one 8mer or two 7mer-m8 (SITE_WEIGHT).
        """
        df = self.scan_state_counts(state)
        return df[(df["score"] >= min_score) & (df["n_sites"] >= min_sites)]

    def scan_state_arrays(self, state):
        """(gene_idx int32, best_rank int8) for genes with any site."""
        core = state.motifs["6mer"]
        b8 = _B2I[state.motifs["7mer-m8"][0]]
        cc = self._core_code(core)
        lo = np.searchsorted(self._sorted_code, cc, "left")
        hi = np.searchsorted(self._sorted_code, cc, "right")
        best = np.zeros(self.genes.shape[0], dtype=np.int8)
        if hi > lo:
            g = self._gene_idx[lo:hi]
            m8 = self._flank5[lo:hi] == b8
            a1 = self._a1[lo:hi]
            rank = np.where(
                m8 & a1, 4, np.where(m8 & ~a1, 3, np.where(~m8 & a1, 2, 1))
            ).astype(np.int8)
            np.maximum.at(best, g, rank)
        idx = np.nonzero(best > 0)[0].astype(np.int32)
        return idx, best[idx]

    def scan_state_arrays_full(self, state):
        """gene_idx, best_rank, n_8mer, n_7mer_m8 for genes with any site."""
        df = self.scan_state_counts(state)
        if df.empty:
            z = np.array([], dtype=np.int32)
            e = z.astype(np.int8)
            return z, e, e, e
        return (
            df["gene_idx"].to_numpy(dtype=np.int32),
            df["best_rank"].to_numpy(dtype=np.int8),
            np.clip(df["n_8mer"], 0, 127).to_numpy(dtype=np.int8),
            np.clip(df["n_7mer_m8"], 0, 127).to_numpy(dtype=np.int8),
        )

    def target_set(self, state, min_rank: int = 1) -> set:
        df = self.scan_state(state)
        return set(df.loc[df["site_rank"] >= min_rank, "symbol"])
