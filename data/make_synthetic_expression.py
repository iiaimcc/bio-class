#!/usr/bin/env python3
# make_synthetic_expression.py
# Generates a synthetic RNA-seq count matrix simulating a lung adenocarcinoma
# (LUAD) vs. normal study. Pure Python stdlib, no third-party dependencies.
#
# Design (teaching-oriented):
#   - 10,000 genes x 12 samples (6 LUAD / 6 normal, 2 batches)
#   - Negative-binomial counts with realistic library sizes and gene means
#   - ~3% differentially expressed genes; EGFR/MAPK signaling upregulated,
#     a few tumor suppressors downregulated (ground truth in *_truth.csv)
#   - A curated "druggable panel" of well-known drug targets is provided so
#     students can intersect DEGs with the panel to shortlist targets.
#
# Usage: python3 make_synthetic_expression.py [-o OUT_DIR]

import argparse
import csv
import math
import os
import random

OUT_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

# ---------------------------------------------------------------------------
# Curated lists of real human gene symbols (signaling, tumor, housekeeping...)
# ---------------------------------------------------------------------------

ONCO_SIGNAL = """EGFR ERBB2 ERBB3 ERBB4 KRAS NRAS HRAS BRAF RAF1 MAP2K1 MAP2K2 MAPK1 MAPK3
AKT1 AKT2 AKT3 MTOR RPTOR RICTOR TSC1 TSC2 RHEB PIK3CA PIK3CB PIK3CD PIK3R1 PTEN
ALK ROS1 RET MET NTRK1 NTRK2 NTRK3 FGFR1 FGFR2 FGFR3 FGFR4 VEGFA VEGFB KDR FLT1 FLT4
PDGFRA PDGFRB KIT FLT3 CSF1R JAK1 JAK2 JAK3 TYK2 STAT1 STAT3 STAT5A STAT5B
SRC YES1 FYN LYN LCK FGR HCK BLK ABL1 SYK ZAP70 BTK ITK TEC ROCK1 ROCK2
PRKCA PRKCB PRKCD PRKCE CDC42 RAC1 RHOA LIMK1 LIMK2 MAPK8 MAPK9 MAPK11 MAPK14""".split()

TUMOR_SUPPRESS = """TP53 MDM2 CDKN1A CDKN1B CDKN2A CDKN2B RB1 CCND1 CCNE1 CDK2 CDK4 CDK6
E2F1 MYC BCL2 BCL2L1 BAX BAK1 BAD BID MCL1 BBC3 PMAIP1 CASP3 CASP8 CASP9
TNFRSF1A TNFRSF1B NFKB1 NFKB2 RELA CHUK IKBKB NFKBIA TNFAIP3
BRCA1 BRCA2 RAD51 ATM ATR CHEK1 CHEK2 EPCAM""".split()

CYTOKINE_ISG = """TNF IL6 IL6R IL1B IL1A CXCL8 CXCL1 CXCL2 CCL2 CCL5 IL10 IL4 IL13 IFNG
IFNA1 IFNB1 OAS1 OAS2 MX1 ISG15 IRF1 IRF7 STAT2 TLR4 CD274 PDCD1 CTLA4 HAVCR2
TGFB1 TGFBR1 TGFBR2 SMAD2 SMAD3 SMAD4 WNT1 CTNNB1 APC GSK3B AXIN1 AXIN2
NOTCH1 NOTCH2 NOTCH3 HIF1A EPAS1 VHL ARNT TERT""".split()

MISC_TARGETS = """ACE ACE2 AGTR1 ADRB1 ADRB2 ADRA1A HTR1A HTR2A HTR3A DRD1 DRD2 DRD3
CHRM1 CHRM2 HIST1H3A HDAC1 HDAC2 HDAC3 SIRT1 TOP1 TOP2A DHFR TYMS TSHR
PDE5A PDE4B ALOX5 PTGS1 PTGS2 HMGCR PPARG PPARA NR3C1 NR3C2 ESR1 ESR2 AR
CYP3A4 CYP2D6 CYP2C19 CYP2C9 CYP1A2 UGT1A1 SLC6A4 SLC6A3 ABCB1 ABCG2 SLCO1B1""".split()

HOUSEKEEPING = """GAPDH ACTB TUBA1A TUBB B2M RPLP0 RPS18 HPRT1 TBP GUSB PPIA TFRC
PGK1 YWHAZ HMBS SDHA CYC1 ATP5FB NDUFA1 UQCRB EEF1A1 EEF2 RPL13A RPL32""".split()

EGFR_PATHWAY = """EGFR ERBB2 KRAS NRAS BRAF RAF1 MAP2K1 MAP2K2 MAPK1 MAPK3 AKT1 AKT2
PIK3CA PIK3R1 MTOR RPS6KB1 RPS6 TSC2 VEGFA EGF AREG HBEGF EREG TGFA NRG1 GRB2
SOS1 SHC1 JAK1 STAT3 MYC CCND1 MCL1 BCL2L1 BIRC5 SERPINE1 PLAU MMP9""".split()

DOWN_PATHWAY = """CDKN1A CDKN2A PTEN BCL2 BAX APAF1 BAD BIM""".split()
DOWN_PATHWAY = """CDKN1A CDKN2A PTEN BAX APAF1""".split()

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------

SEED = 20260901
N_GENES = 10000
SAMPLES = ["S{:02d}".format(i) for i in range(1, 13)]  # S01..S12
CONDITION = ["Normal"] * 6 + ["Tumor"] * 6              # S01-S06 normal, S07-S12 tumor
BATCH = ["B1", "B1", "B2", "B2", "B1", "B2"] * 2
BATCH_EFFECT = {"B1": 1.12, "B2": 0.88}                 # multiplicative library offset
LIBRARY_SIZE = [9.5e6] * 12                              # per-sample reads (scales counts)
N_DIFF = 300                                             # true DEG count (~3%)
UP_RATIO = 0.8

# log2FC ranges by class
LOG2FC = {
    "pathway": (1.2, 2.8),   # EGFR/MAPK core
    "other_up": (0.8, 1.8),
    "down": (-2.2, -1.0),
}


def log_normal_mean(rng):
    """Gene baseline expression: log-normal, mostly low-medium expressed."""
    return math.exp(rng.gauss(4.4, 1.35))  # median ~ 100 counts


def neg_binom(rng, lambda_mean, theta=0.4):
    """One negative-binomial count: NB(mean, theta), theta = 1/overdispersion.

    All draws must come from the seeded RNG; using the global random module
    here would make the matrix irreproducible even with a fixed SEED.
    """
    r = 1.0 / theta
    lam = rng.gammavariate(alpha=r, beta=lambda_mean / r)  # gamma mean = alpha*beta
    # Poisson sample (Knuth)
    if lam < 30:
        k, p = 0, 1.0
        exp_m = math.exp(-lam)
        while p > exp_m:
            k += 1
            p *= rng.random()
        return k - 1
    return int(round(lam + math.sqrt(lam) * rng.gauss(0.5, 0.8)))


def build_gene_table(rng):
    """Returns list of dicts: symbol, ensembl, chr, kind, baseline(per 1M reads)."""
    genes = []
    used = set()

    def add(sym, kind, base_scale=1.0):
        if sym in used or not sym:
            return
        used.add(sym)
        genes.append({
            "symbol": sym,
            "kind": kind,
            "scale": base_scale,
        })

    # Pathway members get a higher baseline so the designed scenario (a
    # strongly up-regulated EGFR/MAPK core) does not depend on a lucky draw.
    # EGFR is the anchor of the whole course chain (omics -> structure ->
    # design), so it also gets a higher baseline and a pinned effect below.
    for s in ONCO_SIGNAL + TUMOR_SUPPRESS + CYTOKINE_ISG + MISC_TARGETS:
        scale = 20.0 if s == "EGFR" else (4.0 if s in EGFR_PATHWAY else 1.0)
        add(s, "curated", base_scale=scale)
    for s in HOUSEKEEPING:
        add(s, "housekeeping", base_scale=3.0)

    # Fill the rest with HGNC-style simulated names (e.g. C1orf100, C10orf88)
    while len(genes) < N_GENES:
        chr_i = rng.randint(1, 22)
        num = rng.randint(4, 998)
        sym = "C{}orf{}".format(chr_i, num)
        kind = "lncRNA" if rng.random() < 0.12 else "protein_coding"
        add(sym, kind)

    # Assign baseline + fixed chromosome/ensembl ids
    for g in genes:
        g["chr"] = "{}x".format(rng.randint(1, 22)) if g["kind"] == "lncRNA" else str(rng.randint(1, 22))
        g["ensembl"] = "ENSG{:011d}".format(rng.randint(10**10, 9 * 10**10))
    return genes


def assign_differential_genes(genes, rng):
    """Pick N_DIFF genes as true DEGs; return dict symbol -> log2FC (tumor vs normal)."""
    # 1) guaranteed EGFR/MAPK pathway members (subset present in table)
    pathway = [g["symbol"] for g in genes if g["symbol"] in EGFR_PATHWAY and g["kind"] != "lncRNA"]
    down = [g["symbol"] for g in genes if g["symbol"] in DOWN_PATHWAY]
    rng.shuffle(pathway)
    rng.shuffle(down)

    n_pathway_up = min(30, len(pathway))
    n_down = 12
    deg = {s: rng.uniform(*LOG2FC["pathway"]) for s in pathway[:n_pathway_up]}
    deg.update({s: rng.uniform(*LOG2FC["down"]) for s in down[:n_down]})
    # Pin the anchor gene's effect: EGFR must stay the clear top hit no matter
    # which seed is used, otherwise the course narrative depends on luck.
    deg["EGFR"] = 3.5

    # 2) remaining DEGs picked at random
    eligible = [g["symbol"] for g in genes if g["symbol"] not in deg and g["kind"] != "lncRNA"]
    rng.shuffle(eligible)
    n_remaining = N_DIFF - len(deg)
    for s in eligible[:n_remaining]:
        upper = rng.random() < UP_RATIO
        deg[s] = rng.uniform(*LOG2FC["other_up" if upper else "down"])
    return deg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default=OUT_DEFAULT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rng = random.Random(SEED)
    genes = build_gene_table(rng)
    deg_map = assign_differential_genes(genes, rng)
    print("genes: {} | DEGs: {}".format(len(genes), len(deg_map)))

    # baseline per gene (counts per 1M reads); lncRNA-like genes kept low to
    # reproduce the "many silent genes" pattern of real RNA-seq
    for g in genes:
        base = log_normal_mean(rng) * g["scale"]
        if g["kind"] == "lncRNA":
            base *= 0.04
        g["base"] = base

    # build counts matrix
    counts = {}          # sample -> [gene x] counts
    for si, sid in enumerate(SAMPLES):
        cond = CONDITION[si]
        batch = BATCH[si]
        row = []
        for g in genes:
            mu = g["base"]
            if g["symbol"] in deg_map and cond == "Tumor":
                mu *= 2.0 ** deg_map[g["symbol"]]
            lam = mu * LIBRARY_SIZE[si] / 1e6 * BATCH_EFFECT[batch]
            row.append(neg_binom(rng, lam))
        counts[sid] = row

    # ---- write outputs ----
    def w(path, fn):
        with open(path, "w", newline="") as f:
            fn(f)

    # expression matrix: rows=genes, cols=samples
    w(os.path.join(args.out, "expression_matrix.csv"), lambda f: _w_matrix(f, genes, counts))
    # sample annotation
    w(os.path.join(args.out, "sample_meta.csv"), lambda f: _w_meta(f))
    # gene annotation
    w(os.path.join(args.out, "gene_annotation.csv"), lambda f: _w_annotation(f, genes, deg_map))
    # druggable panel (curated, for target shortlisting)
    w(os.path.join(args.out, "druggable_panel.csv"), lambda f: _w_panel(f))
    # ground truth (teacher only)
    w(os.path.join(args.out, "ground_truth_deg.csv"), lambda f: _w_truth(f, deg_map))

    print("wrote {} : {}".format(args.out,
          ", ".join(os.listdir(args.out))))


def _w_matrix(f, genes, counts):
    wr = csv.writer(f)
    wr.writerow([""] + SAMPLES)
    for gi, g in enumerate(genes):
        wr.writerow([g["symbol"]] + [counts[s][gi] for s in SAMPLES])


def _w_meta(f):
    wr = csv.writer(f)
    wr.writerow(["sample_id", "condition", "batch"])
    for i, s in enumerate(SAMPLES):
        wr.writerow([s, CONDITION[i], BATCH[i]])


def _w_annotation(f, genes, deg_map):
    wr = csv.writer(f)
    wr.writerow(["symbol", "ensembl_gene_id", "chr", "gene_type", "true_deg"])
    for g in genes:
        wr.writerow([g["symbol"], g["ensembl"], g["chr"], g["kind"],
                     "yes" if g["symbol"] in deg_map else "no"])


def _w_panel(f):
    panel = (set(ONCO_SIGNAL) | set(TUMOR_SUPPRESS) | set(CYTOKINE_ISG)
             | set(MISC_TARGETS) | set(HOUSEKEEPING)) - {"TFRC"}
    # annotate class for the panel
    class_map = {}
    for n, s in (("receptor_kinase", ONCO_SIGNAL), ("tumor_suppressor", TUMOR_SUPPRESS),
                 ("cytokine_immunity", CYTOKINE_ISG), ("pharma_gpcr_enzymes", MISC_TARGETS),
                 ("housekeeping", HOUSEKEEPING)):
        for sym in s:
            class_map.setdefault(sym, n)
    wr = csv.writer(f)
    wr.writerow(["symbol", "panel_class", "known_drug_class"])
    for sym in sorted(panel):
        wr.writerow([sym, class_map[sym], "approved_or_clinical"])


def _w_truth(f, deg_map):
    wr = csv.writer(f)
    wr.writerow(["symbol", "true_log2FC"])
    for sym, fc in sorted(deg_map.items(), key=lambda kv: kv[1], reverse=True):
        wr.writerow([sym, round(fc, 3)])


if __name__ == "__main__":
    main()