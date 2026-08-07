"""
o8g_db.py
=========
Read-only data layer over the precomputed o8G target database (o8g_targets.db).

Schema (v1 → v2 precision migration)
------------------------------------
states columns (v1): gene_blob, rank_blob  (rank>=3 strong sites)
optional (v2): n8_blob, n7m8_blob, cons_blob  (parallel int8 arrays)

Backward compatible: missing optional blobs → columns filled with defaults /
on-the-fly enrichment when a TargetScanner is supplied.
"""
from __future__ import annotations

import sqlite3
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

from o8g_precision import (
    PrecisionConfig,
    PrecisionMode,
    apply_precision_filter,
    partition_after_filter,
)

RANK_SITE = {4: "8mer", 3: "7mer-m8", 2: "7mer-A1", 1: "6mer"}

_ROOT = Path(__file__).resolve().parent
_DEFAULT_REVERSE = _ROOT / "o8g_reverse.db"
SCHEMA_VERSION_KEY = "schema_version"


def _primary_mirna(names: list[str]) -> str:
    if not names:
        return ""
    threep = [n for n in names if n.endswith("-3p")]
    pool = threep or names
    for pref in ("hsa-miR-1-3p", "hsa-miR-1", "hsa-miR-124-3p"):
        if pref in pool:
            return pref
    return sorted(pool)[0]


class TargetDB:
    def __init__(self, path: str = "o8g_targets.db", reverse_path: str | Path | None = None):
        self.path = path
        self._con = sqlite3.connect(path, check_same_thread=False)
        g = pd.read_sql("SELECT gene_idx, gene_id, symbol FROM genes", self._con)
        self.symbols = g.sort_values("gene_idx")["symbol"].to_numpy()
        self.gene_ids = g.sort_values("gene_idx")["gene_id"].to_numpy()
        self._gene_table = g.set_index("gene_idx")
        self._reverse_candidate = Path(reverse_path) if reverse_path else _DEFAULT_REVERSE
        self._rev: sqlite3.Connection | None = None
        self._state_cols = {
            r[1] for r in self._con.execute("PRAGMA table_info(states)").fetchall()
        }
        self.schema_version = self._read_schema_version()

    def _read_schema_version(self) -> int:
        tables = {
            r[0]
            for r in self._con.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
        }
        if "meta" not in tables:
            return 1
        row = self._con.execute(
            "SELECT value FROM meta WHERE key=?", [SCHEMA_VERSION_KEY]
        ).fetchone()
        return int(row[0]) if row else 1

    @property
    def reverse_path(self) -> Path | None:
        return self._reverse_candidate if self._reverse_candidate.exists() else None

    def _rev_con(self) -> sqlite3.Connection:
        if self.reverse_path is None:
            raise FileNotFoundError(
                "o8g_reverse.db not found. Run: python scripts/build_reverse_index.py"
            )
        if self._rev is None:
            self._rev = sqlite3.connect(str(self.reverse_path), check_same_thread=False)
        return self._rev

    def mirnas(self) -> pd.DataFrame:
        return pd.read_sql(
            "SELECT mirna, accession, seq_dna, seed, n_G FROM mirnas ORDER BY mirna",
            self._con,
        )

    def search_mirnas(self, query: str) -> pd.DataFrame:
        q = f"%{query}%"
        return pd.read_sql(
            "SELECT mirna, accession, seq_dna, seed, n_G FROM mirnas "
            "WHERE mirna LIKE ? ORDER BY mirna",
            self._con,
            params=[q],
        )

    def mirna_info(self, mirna: str) -> dict | None:
        r = self._con.execute(
            "SELECT mirna, accession, seq_dna, seed, n_G FROM mirnas WHERE mirna=?",
            [mirna],
        ).fetchone()
        if not r:
            return None
        return dict(zip(["mirna", "accession", "seq_dna", "seed", "n_G"], r))

    def states_for_seed(self, seed: str) -> pd.DataFrame:
        return pd.read_sql(
            "SELECT state_id, label, oxidized_positions, motif_6mer, motif_7mer_m8, motif_8mer, "
            "n_6mer, n_7mer_A1, n_7mer_m8, n_8mer, n_strong FROM states WHERE seed=? ORDER BY state_id",
            self._con,
            params=[seed],
        )

    def _decode_state_blobs(self, seed: str, label: str) -> pd.DataFrame:
        cols = ["gene_blob", "rank_blob"]
        optional = []
        for c in ("n8_blob", "n7m8_blob", "cons_blob"):
            if c in self._state_cols:
                optional.append(c)
                cols.append(c)
        sql = f"SELECT {', '.join(cols)} FROM states WHERE seed=? AND label=?"
        r = self._con.execute(sql, [seed, label]).fetchone()
        empty = pd.DataFrame(
            columns=[
                "gene_idx",
                "gene_id",
                "symbol",
                "site_type",
                "site_rank",
                "n_8mer",
                "n_7mer_m8",
                "n_sites",
                "score",
                "is_conserved",
            ]
        )
        if not r:
            return empty
        gidx = np.frombuffer(zlib.decompress(r[0]), dtype=np.int32)
        rnk = np.frombuffer(zlib.decompress(r[1]), dtype=np.int8)
        n = len(gidx)
        n8 = np.zeros(n, dtype=np.int8)
        n7 = np.zeros(n, dtype=np.int8)
        cons = np.zeros(n, dtype=np.int8)
        # Map optional blobs by name
        colmap = {c: r[i] for i, c in enumerate(cols)}
        if colmap.get("n8_blob"):
            n8 = np.frombuffer(zlib.decompress(colmap["n8_blob"]), dtype=np.int8)
        if colmap.get("n7m8_blob"):
            n7 = np.frombuffer(zlib.decompress(colmap["n7m8_blob"]), dtype=np.int8)
        if colmap.get("cons_blob"):
            cons = np.frombuffer(zlib.decompress(colmap["cons_blob"]), dtype=np.int8)
        # If multiplicity blobs missing, approximate from best rank:
        # 8mer → n8=1; 7mer-m8 → n7=1 (underestimates multi-site genes).
        if "n8_blob" not in colmap or colmap["n8_blob"] is None:
            n8 = (rnk == 4).astype(np.int8)
            n7 = (rnk == 3).astype(np.int8)
        score = n8.astype(np.float64) * 1.0 + n7.astype(np.float64) * 0.7
        n_sites = n8.astype(np.int32) + n7.astype(np.int32)
        return pd.DataFrame(
            {
                "gene_idx": gidx,
                "gene_id": self.gene_ids[gidx],
                "symbol": self.symbols[gidx],
                "site_rank": rnk,
                "site_type": [RANK_SITE[int(x)] for x in rnk],
                "n_8mer": n8.astype(int),
                "n_7mer_m8": n7.astype(int),
                "n_sites": n_sites,
                "score": score,
                "is_conserved": cons.astype(bool),
            }
        )

    def targets(self, seed: str, label: str) -> pd.DataFrame:
        """Strong-site target genes for one seed-state (backward-compatible API)."""
        df = self._decode_state_blobs(seed, label)
        return (
            df.drop(columns=["gene_idx"], errors="ignore")
            .sort_values(["site_rank", "symbol"], ascending=[False, True])
            .reset_index(drop=True)
        )

    def targets_enriched(
        self,
        seed: str,
        label: str,
        *,
        scanner=None,
        mature_dna: str | None = None,
        conserved_symbols: set[str] | None = None,
    ) -> pd.DataFrame:
        """Targets with score / multiplicity / optional live context + conservation."""
        from o8g_engine import SeedState

        df = self._decode_state_blobs(seed, label)
        if conserved_symbols is not None:
            df["is_conserved"] = df["symbol"].isin(conserved_symbols)
        if scanner is not None:
            # Live multiplicity + context from UTR index (authoritative when present)
            ox = []
            if label != "none":
                # parse o8G@2,7 → (2,7)
                part = label.replace("o8G@", "")
                ox = [int(x) for x in part.split(",") if x]
            state = SeedState(seed, tuple(ox))
            live = scanner.scan_state_context(state, mature_dna=mature_dna, min_rank=3)
            if not live.empty:
                keep = [
                    "symbol",
                    "n_8mer",
                    "n_7mer_m8",
                    "n_sites",
                    "score",
                    "context_score",
                    "site_rank",
                    "site_type",
                    "gene_id",
                ]
                live = live[keep]
                df = live.merge(
                    df[["symbol", "is_conserved"]],
                    on="symbol",
                    how="left",
                )
                df["is_conserved"] = df["is_conserved"].fillna(False)
        return df.sort_values(["site_rank", "symbol"], ascending=[False, True]).reset_index(
            drop=True
        )

    def targets_filtered(
        self,
        seed: str,
        label: str,
        cfg: PrecisionConfig | PrecisionMode | str,
        *,
        scanner=None,
        mature_dna: str | None = None,
        conserved_symbols: set[str] | None = None,
        mirna: str | None = None,
    ) -> pd.DataFrame:
        if not isinstance(cfg, PrecisionConfig):
            cfg = PrecisionConfig.from_mode(cfg)
        if conserved_symbols is None and cfg.mode == PrecisionMode.CONSENSUS and mirna:
            try:
                from conservation import get_conserved_index, build_seed_family_map

                idx = get_conserved_index()
                fam = build_seed_family_map(self.path).get(seed)
                conserved_symbols = idx.conserved_symbols_for_mirna(mirna, fam)
            except Exception:
                conserved_symbols = set()
        df = self.targets_enriched(
            seed,
            label,
            scanner=scanner,
            mature_dna=mature_dna,
            conserved_symbols=conserved_symbols,
        )
        return apply_precision_filter(
            df,
            cfg,
            conserved_symbols=conserved_symbols,
            is_unmodified_state=(label == "none"),
        )

    def retarget_partition(
        self,
        seed: str,
        ox_label: str,
        cfg: PrecisionConfig | PrecisionMode | str,
        **kwargs,
    ) -> dict[str, set[str]]:
        """Partition unmod vs oxidized after filtering both sides."""
        if not isinstance(cfg, PrecisionConfig):
            cfg = PrecisionConfig.from_mode(cfg)
        mirna = kwargs.get("mirna")
        conserved_symbols = kwargs.get("conserved_symbols")
        if conserved_symbols is None and cfg.mode == PrecisionMode.CONSENSUS and mirna:
            try:
                from conservation import get_conserved_index, build_seed_family_map

                idx = get_conserved_index()
                fam = build_seed_family_map(self.path).get(seed)
                conserved_symbols = idx.conserved_symbols_for_mirna(mirna, fam)
            except Exception:
                conserved_symbols = set()
        unmod = self.targets_enriched(
            seed,
            "none",
            scanner=kwargs.get("scanner"),
            mature_dna=kwargs.get("mature_dna"),
            conserved_symbols=conserved_symbols,
        )
        oxid = self.targets_enriched(
            seed,
            ox_label,
            scanner=kwargs.get("scanner"),
            mature_dna=kwargs.get("mature_dna"),
            conserved_symbols=conserved_symbols,
        )
        return partition_after_filter(
            unmod, oxid, cfg, conserved_symbols=conserved_symbols
        )

    def target_symbols(self, seed: str, label: str, min_rank: int = 3) -> list[str]:
        df = self.targets(seed, label)
        return df.loc[df["site_rank"] >= min_rank, "symbol"].tolist()

    def gene_info(self, gene_idx: int) -> dict | None:
        if gene_idx not in self._gene_table.index:
            return None
        row = self._gene_table.loc[gene_idx]
        return {
            "gene_idx": int(gene_idx),
            "gene_id": str(row["gene_id"]),
            "symbol": str(row["symbol"]),
        }

    def states_targeting_gene(self, gene_idx: int) -> pd.DataFrame:
        rev = self._rev_con()
        gt = pd.read_sql(
            "SELECT seed, label AS state_label, site_rank, oxidized_positions "
            "FROM gene_targets WHERE gene_idx=?",
            rev,
            params=[int(gene_idx)],
        )
        empty_cols = [
            "mirna",
            "all_mirnas",
            "seed",
            "state_label",
            "oxidized_positions",
            "site_rank",
            "site_type",
            "motif_7mer_m8",
            "motif_8mer",
            "vs_unmodified",
        ]
        if gt.empty:
            return pd.DataFrame(columns=empty_cols)

        seeds = gt["seed"].unique().tolist()
        ph = ",".join("?" * len(seeds))
        sm = pd.read_sql(
            f"SELECT seed, mirna FROM seed_mirnas WHERE seed IN ({ph})",
            rev,
            params=seeds,
        )
        mir_map: dict[str, list[str]] = {}
        for seed, mirna in sm.itertuples(index=False):
            mir_map.setdefault(seed, []).append(mirna)

        motifs = pd.read_sql(
            f"SELECT seed, label, motif_7mer_m8, motif_8mer FROM states WHERE seed IN ({ph})",
            self._con,
            params=seeds,
        ).rename(columns={"label": "state_label"})

        none_seeds = set(gt.loc[gt["state_label"] == "none", "seed"])

        def vs_unmod(seed: str, label: str) -> str:
            if label == "none":
                return "unmodified"
            if seed in none_seeds:
                return "also in unmodified"
            return "gained on oxidation"

        gt = gt.merge(motifs, on=["seed", "state_label"], how="left")
        gt["mirna"] = gt["seed"].map(lambda s: _primary_mirna(mir_map.get(s, [])))
        gt["all_mirnas"] = gt["seed"].map(lambda s: ";".join(sorted(mir_map.get(s, []))))
        gt["site_type"] = gt["site_rank"].map(RANK_SITE)
        gt["vs_unmodified"] = [
            vs_unmod(s, lab) for s, lab in zip(gt["seed"], gt["state_label"])
        ]
        gt["oxidized_positions"] = gt["oxidized_positions"].fillna("")
        return (
            gt[empty_cols]
            .sort_values(["mirna", "state_label", "site_rank"], ascending=[True, True, False])
            .reset_index(drop=True)
        )

    def close(self):
        self._con.close()
        if self._rev is not None:
            self._rev.close()
            self._rev = None
