#!/usr/bin/env Rscript
# =============================================================================
# 02_enrichment.R
# Course: "Omics-based disease mechanism dissection and target discovery" (Part 2)
#
# What this script does:
#   1. Reads the DEG list from 01_diff_analysis.R
#   2. Maps gene symbols to Entrez IDs
#   3. Runs GO (Biological Process) and KEGG pathway enrichment
#      with clusterProfiler; KEGG needs internet, so it is wrapped to
#      fall back gracefully if the network is unavailable.
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

# ---- 0.5 Required packages ----------------------------------------------------
for (p in c("clusterProfiler", "org.Hs.eg.db", "ggplot2", "dplyr", "readr")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is missing. Install it first, see 环境与部署指南.md")
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)   # loads the org.Hs.eg.db annotation object
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---- 1. Load DEG list (produced by 01_diff_analysis.R) ------------------------
deg_file <- file.path(out_dir, "deg_table_all.csv")
if (!file.exists(deg_file))
  stop("Missing ", deg_file, " - run 01_diff_analysis.R first (same <out_dir>)")
deg <- read_csv(deg_file, show_col_types = FALSE)
sig <- deg %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
cat("Using", nrow(sig), "significant DEGs for enrichment\n")

# ---- 2. Symbol -> Entrez ID mapping ------------------------------------------
entrez <- bitr(sig$symbol, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat("Mapped", nrow(entrez), "/", nrow(sig), "genes to Entrez IDs\n")
if (nrow(entrez) < 0.5 * nrow(sig)) warning("Low mapping rate - check gene symbols")

# ---- 3. GO Biological Process enrichment --------------------------------------
ego <- enrichGO(
  gene          = entrez$ENTREZID,
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
    ggtitle("GO Biological Process enrichment (DEGs)")
  ggsave(file.path(out_dir, "go_bp_dotplot.png"), p_go,
         width = 9, height = 7, dpi = 300)
  cat("Saved GO dotplot\n")
}

# ---- 4. KEGG pathway enrichment (needs internet; tolerant to failure) ---------
kegg_ok <- TRUE
ekegg <- tryCatch({
  enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
             pvalueCutoff = 0.1, qvalueCutoff = 0.1)
}, error = function(e) {
  message("KEGG enrichment failed (network?): ", conditionMessage(e))
  kegg_ok <<- FALSE
  NULL
})

if (kegg_ok && !is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  kk_df <- as.data.frame(ekegg)
  write_csv(kk_df, file.path(out_dir, "kegg_enrichment.csv"))
  p_kegg <- dotplot(ekegg, showCategory = 15) +
    ggtitle("KEGG pathway enrichment (DEGs)")
  ggsave(file.path(out_dir, "kegg_dotplot.png"), p_kegg,
         width = 9, height = 7, dpi = 300)
  cat("KEGG enrichment terms:", nrow(kk_df), "\n")
} else {
  cat("KEGG step skipped (see message above). Alternative: use the public web\n")
  cat("service Enrichr (https://maayanlab.cloud/Enrichr) with the DEG list.\n")
}

# ---- 5. Console summary -------------------------------------------------------
cat("\n===== Top 15 GO terms =====\n")
if (nrow(ego_df) > 0) {
  print(ego_df %>% arrange(p.adjust) %>% select(ID, Description, Count, p.adjust) %>% head(15))
}

# Save the DEG symbol list for step 3 (target shortlist)
write_lines(sig$symbol, file.path(out_dir, "deg_symbols.txt"))
cat("\nSaved DEG symbol list (deg_symbols.txt). Next: run 03_target_shortlist.R\n")