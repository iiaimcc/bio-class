#!/usr/bin/env Rscript
# =============================================================================
# 02_enrichment.R
# Course: "Omics-based disease mechanism dissection and target discovery" (Part 2)
#
# What this script does:
#   1. Reads the DEG table from 01_diff_analysis.R
#   2. Maps gene symbols to Entrez IDs
#   3. Runs GO (Biological Process) and Reactome pathway enrichment. Both are
#      given the tested genes as `universe`; without it they would fall back to
#      the full annotation background.
#      Both use local annotation packages (org.Hs.eg.db, reactome.db),
#      so the whole script works offline. KEGG is intentionally not used
#      here: it needs a live REST download that many classroom networks block.
#   4. Saves enrichment tables and dot plots
#
# How to run:
#   Rscript 02_enrichment.R  <data_dir>  <out_dir>
# =============================================================================

# ---- 0. Configuration --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "."
out_dir  <- if (length(args) >= 2) args[2] else "./results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# Annotation packages are local; keep a generous timeout anyway because the
# first ReactomePA call can be slow while R builds the pathway index.
options(timeout = 600)

# ---- 0.5 Required packages ----------------------------------------------------
for (p in c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "reactome.db",
            "ggplot2", "dplyr", "readr")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is missing. Install it first, see 环境与部署指南.md")
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)   # loads the org.Hs.eg.db annotation object
  library(ReactomePA)     # pathway enrichment
  library(reactome.db)    # local Reactome annotation, no download needed
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---- 1. Load DEG list (produced by 01_diff_analysis.R) ------------------------
deg_file <- file.path(out_dir, "deg_table_all.csv")
if (!file.exists(deg_file))
  stop("Missing ", deg_file, " - run 01_diff_analysis.R first (same <out_dir>)")
deg <- read_csv(deg_file, show_col_types = FALSE)
sig <- deg %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1,
         !is.na(symbol), symbol != "") %>%
  distinct(symbol, .keep_all = TRUE)
cat("Using", nrow(sig), "unique significant DEG symbols for enrichment\n")

# ---- 2. Symbol -> Entrez ID mapping ------------------------------------------
entrez <- bitr(sig$symbol, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat("Mapped", nrow(entrez), "/", nrow(sig), "genes to Entrez IDs\n")
if (nrow(entrez) < 0.5 * nrow(sig)) warning("Low mapping rate - check gene symbols")

# ---- 2.5 Background (universe) -------------------------------------------------
# The tested set is the genes that passed the low-count filter and actually
# entered the DESeq2 test. Leaving `universe` unset makes enrichGO fall back to
# every gene annotated in org.Hs.eg.db and enrichPathway to every gene in the
# Reactome annotation - i.e. tens of thousands of genes that this experiment
# could never have detected, which inflates significance. The tested set is the
# only defensible denominator, and mapping it to Entrez the same way as the DEGs
# keeps numerator and denominator on the same identifier space.
universe_sym <- deg %>%
  filter(!is.na(pvalue), !is.na(symbol), symbol != "") %>%
  pull(symbol) %>%
  unique()
universe <- bitr(universe_sym, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)$ENTREZID
universe <- unique(universe)
cat("Tested genes (universe):", length(universe_sym), "symbols ->",
    length(universe), "Entrez IDs\n")

# ---- 3. GO Biological Process enrichment --------------------------------------
ego <- enrichGO(
  gene          = entrez$ENTREZID,
  universe      = universe,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",            # Biological Process
  keyType       = "ENTREZID",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE             # show gene symbols instead of Entrez IDs
)
ego_df <- as.data.frame(ego)
write_csv(ego_df, file.path(out_dir, "go_bp_enrichment.csv"))
cat("GO enrichment terms:", nrow(ego_df), "\n")

if (nrow(ego_df) > 0) {
  p_go <- dotplot(ego, showCategory = 15) +
    ggtitle("GO Biological Process enrichment (TCGA-LUAD DEGs)")
  ggsave(file.path(out_dir, "go_bp_dotplot.png"), p_go,
         width = 9, height = 7, dpi = 300)
  cat("Saved GO dotplot\n")
}

# ---- 4. Reactome pathway enrichment (offline, via reactome.db) ----------------
ereact <- enrichPathway(
  gene          = entrez$ENTREZID,
  universe      = universe,
  organism      = "human",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
react_df <- as.data.frame(ereact)
write_csv(react_df, file.path(out_dir, "reactome_enrichment.csv"))
cat("Reactome terms:", nrow(react_df), "\n")

if (nrow(react_df) > 0) {
  p_react <- dotplot(ereact, showCategory = 15) +
    ggtitle("Reactome pathway enrichment (TCGA-LUAD DEGs)")
  ggsave(file.path(out_dir, "reactome_dotplot.png"), p_react,
         width = 9, height = 7, dpi = 300)
  cat("Saved Reactome dotplot\n")
}

# ---- 5. Console summary -------------------------------------------------------
cat("\n===== Top 15 GO terms =====\n")
if (nrow(ego_df) > 0) {
  print(ego_df %>% arrange(p.adjust) %>%
          select(ID, Description, Count, p.adjust) %>% head(15))
}
cat("\n===== Top 15 Reactome pathways =====\n")
if (nrow(react_df) > 0) {
  print(react_df %>% arrange(p.adjust) %>%
          select(ID, Description, Count, p.adjust) %>% head(15))
}
cat("\nAlternative online enrichment services (if local packages are missing):\n")
cat("Enrichr https://maayanlab.cloud/Enrichr , KEGG https://www.kegg.jp\n")

# Save the DEG symbol list for step 3 (target shortlist)
write_lines(sig$symbol, file.path(out_dir, "deg_symbols.txt"))
cat("\nSaved DEG symbol list (deg_symbols.txt). Next: run 03_target_shortlist.R\n")
