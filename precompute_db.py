"""Precompute the o8G miRNA target-list database.

Outputs:
  o8g_targets.db   SQLite schema v2:
      genes, mirnas, states(+ optional n8/n7m8/cons blobs), meta
  o8g_states.parquet   per-state metadata + tier counts

See docs/MIGRATION_PRECISION.md for v1→v2 notes.
"""
import sqlite3, zlib, time, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, ".")
from o8g_engine import enumerate_states
from o8g_scanner import TargetScanner

STRONG = 3   # rank>=3 : 7mer-m8 + 8mer

# Write precision blobs (n8, n7m8, cons) into states. Requires TargetScan tables
# under paper/data for cons_blob; if missing, cons_blob is all zeros.
WRITE_PRECISION_BLOBS = True


def main():
    mir = pd.read_parquet("hsa_mature.parquet")
    utr = pd.read_parquet("utr3_human.parquet")
    sc = TargetScanner(utr["gene_id"].tolist(), utr["symbol"].tolist())
    sc.build(utr["utr3"].tolist())
    print(f"index built; {sc._sorted_code.shape[0]:,} 6mers", flush=True)

    conserved_by_seed = {}
    if WRITE_PRECISION_BLOBS:
        try:
            from conservation import get_conserved_index, build_seed_family_map, TARGETSCAN_RELEASE
            cons = get_conserved_index()
            cons.ensure_loaded()
            fam_map = build_seed_family_map("o8g_targets.db") if Path("o8g_targets.db").exists() else {
                "GGAATGT": "miR-1-3p/206",
                "AAGGCAC": "miR-124-3p.1",
                "GAGGTAG": "let-7-5p/98-5p",
                "GGAGTGT": "miR-122-5p",
            }
            seed_mirna = mir.drop_duplicates("seed").set_index("seed")["mirna"].to_dict()
            for seed, mname in seed_mirna.items():
                conserved_by_seed[seed] = cons.conserved_symbols_for_mirna(
                    mname, fam_map.get(seed)
                )
            ts_release = TARGETSCAN_RELEASE
        except Exception as e:
            print(f"WARN conservation unavailable ({e}); cons_blob=0", flush=True)
            ts_release = "unavailable"
    else:
        ts_release = "n/a"

    con = sqlite3.connect("o8g_targets.db")
    cur = con.cursor()
    cur.executescript("""
    DROP TABLE IF EXISTS genes; DROP TABLE IF EXISTS mirnas;
    DROP TABLE IF EXISTS states; DROP TABLE IF EXISTS meta;
    CREATE TABLE genes(gene_idx INTEGER PRIMARY KEY, gene_id TEXT, symbol TEXT);
    CREATE TABLE mirnas(mirna TEXT PRIMARY KEY, accession TEXT, seq_dna TEXT, seed TEXT, n_G INTEGER);
    CREATE TABLE states(state_id INTEGER PRIMARY KEY, seed TEXT, oxidized_positions TEXT, label TEXT,
        motif_6mer TEXT, motif_7mer_m8 TEXT, motif_8mer TEXT,
        n_6mer INT, n_7mer_A1 INT, n_7mer_m8 INT, n_8mer INT, n_strong INT,
        gene_blob BLOB, rank_blob BLOB,
        n8_blob BLOB, n7m8_blob BLOB, cons_blob BLOB);
    CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);
    CREATE INDEX idx_states_seed ON states(seed);
    """)
    cur.executemany("INSERT INTO genes VALUES (?,?,?)",
                    [(i, g, s) for i,(g,s) in enumerate(zip(sc.genes, sc.symbols))])
    cur.executemany("INSERT INTO mirnas VALUES (?,?,?,?,?)",
                    mir[["mirna","accession","seq_dna","seed","n_G"]].itertuples(index=False, name=None))

    seeds = mir.drop_duplicates("seed")["seed"].tolist()
    sid = 0; t0=time.time(); meta=[]
    for si, seed in enumerate(seeds):
        conserved = conserved_by_seed.get(seed, set())
        for st in enumerate_states(seed):
            if WRITE_PRECISION_BLOBS:
                idx, rank, n8a, n7a = sc.scan_state_arrays_full(st)
            else:
                idx, rank = sc.scan_state_arrays(st)
                n8a = (rank == 4).astype(np.int8)
                n7a = (rank == 3).astype(np.int8)
            n6=int((rank==1).sum()); n7a1=int((rank==2).sum()); n7m8=int((rank==3).sum()); n8=int((rank==4).sum())
            strong_mask = rank >= STRONG
            s_idx = idx[strong_mask]; s_rank = rank[strong_mask]
            s_n8 = n8a[strong_mask]; s_n7 = n7a[strong_mask]
            cons_arr = np.array(
                [1 if sc.symbols[int(i)] in conserved else 0 for i in s_idx],
                dtype=np.int8,
            )
            n_strong = int(strong_mask.sum())
            cur.execute(
                "INSERT INTO states VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (sid, seed, ",".join(map(str, st.oxidized_positions)), st.label,
                 st.motifs["6mer"], st.motifs["7mer-m8"], st.motifs["8mer"],
                 n6, n7a1, n7m8, n8, n_strong,
                 zlib.compress(s_idx.astype(np.int32).tobytes()),
                 zlib.compress(s_rank.astype(np.int8).tobytes()),
                 zlib.compress(s_n8.tobytes()),
                 zlib.compress(s_n7.tobytes()),
                 zlib.compress(cons_arr.tobytes())),
            )
            meta.append((sid, seed, st.label, st.motifs["6mer"], st.motifs["7mer-m8"],
                         n6, n7a1, n7m8, n8, n_strong))
            sid += 1
        if si % 200 == 0:
            con.commit()
            print(f"{si}/{len(seeds)} seeds, {sid} states, {time.time()-t0:.0f}s", flush=True)
    cur.executemany(
        "INSERT INTO meta VALUES (?,?)",
        [
            ("schema_version", "2"),
            ("targetscan_release", ts_release),
            ("strong_rank_min", str(STRONG)),
        ],
    )
    con.commit(); con.close()
    pd.DataFrame(meta, columns=["state_id","seed","label","motif_6mer","motif_7mer_m8",
                                "n_6mer","n_7mer_A1","n_7mer_m8","n_8mer","n_strong"]
                 ).to_parquet("o8g_states.parquet")
    print("done", sid, "states", flush=True)

if __name__ == "__main__":
    main()
