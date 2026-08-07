# o8G-miRNA Retargeting Database — Validation & Worked Examples

**What this is.** A precomputed database and interactive explorer that predicts how
8-oxoguanine (o8G) oxidation of guanines in a microRNA *seed* rewires its mRNA target
repertoire, with built-in pathway-enrichment comparison, volcano and heatmap views.

---

## 1. Biological premise

miRNAs silence mRNAs by Watson–Crick pairing of their **seed** (mature-miRNA
positions 2–8) to complementary sites in target 3′UTRs. An unmodified seed **guanine
pairs cytosine**. When that guanine is oxidized to 8-oxoguanine it rotates to the *syn*
conformation and forms a Hoogsteen pair with **adenine** (o8G·A). So every oxidized seed
position demands an **A** in the target where the unmodified seed demanded a **C** — a
one-base recoding of the targeting rule at that position.

A seed with *k* guanines therefore has **2^k oxidation states**, each defining its own
target-site motif and its own predicted target-gene list. This is the mechanism behind
7-o8G-miR-1 in cardiac hypertrophy (Nature 2020) and widespread seed oxidation of tumor
miRNAs such as miR-124 and miR-122 (Nature 2023).

## 2. Method

| Stage | Choice |
|---|---|
| miRNAs | 2,656 human mature miRNAs (miRBase), 2,094 unique seeds |
| Seed | positions 2–8 (7 nt), DNA alphabet (U→T) |
| Oxidation states | all 2^k combinations of seed-G oxidation (14,245 unique seed-states) |
| Targets | TargetScan-style 6mer / 7mer-A1 / 7mer-m8 / 8mer sites |
| 3′UTRs | longest UTR per gene, 19,159 human genes (Ensembl BioMart) |
| Pairing rule | normal G → C ; o8G → A (Hoogsteen) |
| Enrichment | hypergeometric over-representation, BH-corrected, across GO BP/MF/CC, KEGG, Reactome, MSigDB Hallmark |

For each seed-state the scanner records, per gene, the best site tier and a weighted
site score (8mer 1.0, 7mer-m8 0.7, 7mer-A1 0.4, 6mer 0.15). "Strong-site" target lists
(7mer-m8 + 8mer, rank ≥ 3) are stored in the database; "confident" lists apply a
site-score floor to focus enrichment.

## 3. Validation — canonical target recovery

The scanner was checked against experimentally established targets of two well-studied
miRNAs, using **only sequence** (no prior target knowledge). Recovery = a qualifying
seed-match site is present in the unmodified seed's target set.

**hsa-miR-1-3p** (seed `GGAATGT`, G at positions 2, 3, 7): **11/12** canonical targets
recovered — HDAC4, TWF1, GJA1, KCNJ2, BDNF, CDK6, PURB, SRSF9, PTBP1 as confident targets,
CNN3 and FBXO32 site-present; only MEF2A absent (its validated site is non-canonical /
outside the longest-UTR model).

**hsa-miR-124-3p** (seed `AAGGCAC`, G at positions 4, 5): **12/12** canonical targets
recovered — CDK6, SP1, PTBP1, ROCK1, VAMP3, RAB27A, SLC16A1, IQGAP1, CEBPA, STAT3, SNAI2
confident; EZH2 site-present.

Full table: `validation_canonical_targets.csv`.

## 4. Worked example — miR-1-3p

Seed `GGAATGT` has 3 guanines → **8 oxidation states**. The physiologically important
species is **7-o8G-miR-1** (oxidation at position 7), reported in cardiac hypertrophy.

- Unmodified seed: **9,043** target genes (3,319 strong-site).
- 7-o8G (o8G@7): **8,427** target genes (3,914 strong-site).
- Oxidation at position 7 flips the required target base at that position from C to A:
  7mer-m8 motif `ACATTCC` → `AAATTCC`.
- Canonical target **HDAC4 is lost** on oxidation; ~2,800 strong-site genes are gained.

**Retargeting at the pathway level.** Comparing per-state over-representation (confident
lists, score ≥ 2.0) between normal and 7-o8G finds 13 significantly shifted pathways:
7-o8G-miR-1 **gains** enrichment for MECP2-regulated neuronal receptors, NMDA-receptor /
postsynaptic signalling, and RB1-defect / mitotic-cell-cycle pathways, while **losing**
transmembrane receptor-Ser/Thr-kinase signalling and neuronal action-potential terms.

![miR-1 volcano]({{artifact:c9afb8aa-cd75-4219-881c-89c7c06ced52}})

![miR-1 heatmap]({{artifact:78d4e892-8ce6-418a-86b3-e3482e094502}})

The heatmap across all 8 states shows the shift is position-specific: o8G@7 drives the
NMDA/synaptic + mitotic signature, whereas o8G@3 uniquely lights up MECP2 regulation.

## 5. Worked example — miR-124-3p

Seed `AAGGCAC` has 2 guanines → **4 oxidation states**. miR-124 is a neuronal/tumor-
suppressor miRNA whose seed oxidation is documented in cancer.

- Unmodified: 1,276 confident targets; full oxidation (o8G@4,5): 1,559.
- Differential enrichment (normal vs o8G@4,5) recovers **4 significant retargeting
  pathways**, all *gained* on oxidation: Transcriptional Regulation by MECP2, Regulation
  of MECP2 Expression/Activity, and NMDA-receptor-mediated neuronal transmission.

![miR-124 volcano]({{artifact:aedb42a3-c70d-436e-804a-4f48ae1041cb}})

![miR-124 heatmap]({{artifact:2f9f1622-21fd-490d-93da-dfa05c06cc19}})

**Convergent theme.** Both miR-1 (o8G@7) and miR-124 (o8G@4,5) *gain* the MECP2 /
NMDA-receptor neuronal-signalling module on oxidation — a striking convergence given
their unrelated seeds, and a concrete, testable hypothesis the database surfaces.

## 6. Honest interpretation

Seed-match prediction yields **large target lists** (thousands of genes), so fold-
enrichments for individual pathways are real (5–50×) but the lists overlap heavily and
few pathways survive aggressive FDR when two big lists are compared directly (Fisher
differential is underpowered here). The database therefore frames the volcano/heatmap as
**per-state over-representation compared across states on confidence-filtered lists**,
which recovers genuine, position-specific pathway shifts. Treat the output as
**hypotheses about direction of retargeting**, not as a ranked list of validated targets;
candidate sites should be confirmed experimentally (reporter assays, o8G-specific CLIP).

## 7. Files

Engine: `o8g_engine.py` (state enumeration), `o8g_scanner.py` (UTR scanner),
`o8g_enrich.py` (enrichment), `o8g_plots.py` (figures), `o8g_db.py` (read-only data layer).
Data: `o8g_targets.db` (SQLite, 92 MB), `o8g_states.parquet`, `genesets.tar.gz`.
App: `app.py` (Streamlit) — run with `bash run_app.sh`.
