# Rebuild prompt — o8G-miRNA Retargeting Explorer

Paste everything in the fenced block below into a fresh Claude Science
session (or any capable coding agent with a Python environment and internet
access). It fully specifies the biology, the data sources, the algorithms,
and the deliverables so the project can be reconstructed from nothing.

---

```
Build a database and interactive web app that predicts how 8-oxoguanine
(o8G) oxidation of a human miRNA's seed rewires its mRNA target set.

BIOLOGY / RATIONALE
- A miRNA silences mRNAs by base-pairing its SEED (mature nucleotide
  positions 2-8, a 7-nt window) to complementary sites in target 3'UTRs.
- Normally guanine (G) Watson-Crick pairs cytosine (C). When a seed G is
  oxidized to 8-oxoguanine (o8G), it rotates to the syn conformation and
  Hoogsteen-pairs ADENINE (A) instead. So each oxidized seed G flips the
  target base it demands from C -> A, redirecting the miRNA to a DIFFERENT
  set of mRNAs.
- A seed containing k guanines therefore has 2^k possible oxidation states
  (each G independently normal or o8G), each with its own target list.
- Literature basis: 7-o8G-miR-1 drives cardiac hypertrophy (Nature 2020);
  seed oxidation of miR-124/miR-122 in tumors (Nature 2023).
- SCOPE: human only.

DATA SOURCES
1. Mature miRNA sequences: miRBase mature.fa, filter to human (hsa-). Parse
   to DNA (U->T). Seed = sequence[1:8] (0-based positions 1..7 = mature 2-8).
   ~2656 human mature miRNAs, ~2094 unique seeds.
2. 3'UTR sequences: Ensembl BioMart, human. Take the LONGEST 3'UTR per gene.
   ~19,159 genes with a usable UTR.

ENUMERATION ENGINE (module o8g_engine.py)
- extract_seed(seq), g_positions(seed) -> seed positions (in mature 2-8
  numbering) that are G.
- enumerate_states(seed): for each subset of G positions being oxidized,
  produce a state. state_label examples: "none", "o8G@7", "o8G@2,3,7".
- For each state build the required target-site motifs. The seed pairs the
  UTR antisense; at each seed position, a normal G requires C in the target,
  an o8G requires A. Build the four TargetScan site types from the resulting
  required 6-mer/7-mer/8-mer:
    8mer     = seed match positions 2-8 + A at target position 1  (rank 4)
    7mer-m8  = seed match positions 2-8                            (rank 3)
    7mer-A1  = seed match positions 2-7 + A at position 1          (rank 2)
    6mer     = seed match positions 2-7                            (rank 1)

TARGET SCANNER (module o8g_scanner.py)
- Build a fast index of all 6-mers across all 3'UTRs (~54.5M 6-mers).
- For a given state's motifs, find every gene whose UTR contains a matching
  site; record the best (highest-rank) site type per gene.
- "Strong-site" targets = site rank >= 3 (7mer-m8 and 8mer). Store these.

PRECOMPUTED DATABASE (script precompute_db.py -> o8g_targets.db, SQLite)
- Tables: genes (id, symbol), mirnas (id, name, seed, ...), and per
  seed x state a zlib-compressed blob of int32 gene indices + int8 site
  ranks (rank >= 3 only), plus per-tier counts.
- Enumerate every unique seed x oxidation state (~14,245) and store its
  strong-site target list. Also write o8g_states.parquet with per-state
  metadata + tier counts.
- Read layer o8g_db.py: TargetDB.targets(seed, state_label) -> per-gene
  gene_id/symbol/site_rank/site_type. RANK_SITE = {4:"8mer", 3:"7mer-m8",
  2:"7mer-A1", 1:"6mer"}.

PATHWAY ENRICHMENT (module o8g_enrich.py)
- Cache 6 Enrichr .gmt libraries locally: GO_Biological_Process_2023,
  GO_Molecular_Function_2023, GO_Cellular_Component_2023, KEGG_2021_Human,
  Reactome_2022, MSigDB_Hallmark_2020.
- enrich(query_genes, library, background=None, min_overlap=2): hypergeometric
  over-representation with Benjamini-Hochberg FDR (q-values).
- compare_states(...): per-state ORA compared across states for volcano/heatmap.
- HONESTY NOTE: seed-only prediction yields large, diffuse target lists
  (~3000-3900 strong-site genes per state), so few GO terms survive FDR on
  the full lists. Surface fold-enrichment and pathway shifts honestly;
  volcano compares confidence-filtered lists across states rather than
  forcing significance.

PLOTTING (module o8g_plots.py)
- Volcano: x = log2 odds-ratio (normal vs oxidized), y = -log10 q. Red =
  toward normal / lost targets; blue = toward o8G / gained.
- Heatmap: pathways x oxidation-states, -log10 q, magma, average linkage.
- Also: diverging bar (per-term two-state comparison) and dot plot
  (term x state, dot size = gene overlap).

WEB APP (app.py, Streamlit)
- Sidebar: search/select any human miRNA (all ~2656).
- Tab 1 "Targets for one state": pick an oxidation state, show target gene
  table + site types + download CSV + pathway enrichment.
- Tab 2 "Compare two states": Lost / Gained / Shared gene tables with CSV,
  plus differential pathway view (diverging bar default, volcano toggle).
- Tab 3 "All states": pathway comparison across every oxidation state of the
  seed (heatmap default, dot plot toggle).

DELIVERABLES
- The engine/scanner/db/enrich/plots modules, the precompute script, the
  SQLite database + parquet metadata, the gene-set libraries, the Streamlit
  app, a requirements.txt (streamlit, pandas, numpy, scipy, pyarrow, plotly,
  matplotlib), a Dockerfile, and a README + validation report.
- VALIDATION: recover known canonical targets. miR-1-3p (seed GGAATGT,
  Gs at 2,3,7) should recover canonical targets (HDAC4, TWF1, GJA1, KCNJ2,
  BDNF) in the unmodified state; o8G@7 should lose HDAC4 and redirect ~2800
  strong-site genes. miR-124-3p (seed AAGGCAC, Gs at 4,5) canonical targets
  should recover in the unmodified state.

Work in a dedicated conda env "o8g" (python=3.13, channels conda-forge +
bioconda). Package everything as a runnable bundle.
```

---

## Notes for whoever runs this
- The one-time Ensembl BioMart 3'UTR fetch is the slow part (it pulls
  ~644k UTR rows; took ~40 min in the original build). The target scan
  across all 14,245 states took ~13 min.
- If you only need the SAME database, don't rebuild — just copy
  `o8g_targets.db`, `o8g_states.parquet`, `utr3_human.parquet`, and the
  `genesets/` folder from the bundle. The rebuild is only for changing the
  inputs or the rules.
- Key validation numbers to check you reproduced it correctly: 2656 human
  mature miRNAs, 2094 unique seeds, 14,245 unique seed x state combinations,
  19,159 genes with 3'UTRs, miR-1 unmodified ~9043 total / ~3319 strong-site
  targets, HDAC4 present in miR-1 normal and absent in o8G@7.
