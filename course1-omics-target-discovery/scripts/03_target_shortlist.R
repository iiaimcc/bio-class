#!/usr/bin/env Rscript
# =============================================================================
# 03_target_shortlist.R
# Course: "Omics-based disease mechanism dissection and target discovery" (Part 3)
#
# What this script does:
#   1. Reads the DEG table from step 1 and the curated druggable panel
#   2. Intersects DEGs with the panel -> "druggable candidate targets"
#   3. Ranks candidates by fold change and significance, writes a report
#   4. Prints a teaching-oriented summary table
#
# How to run:
#   Rscript 03_target_shortlist.R  <data_dir>  <out_dir>
# =============================================================================

# ---- 0. Configuration --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "."
out_dir  <- if (length(args) >= 2) args[2] else "./results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (p in c("dplyr", "readr")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' missing")
}
suppressPackageStartupMessages({ library(dplyr); library(readr) })

# ---- 1. Load inputs ------------------------------------------------------------
deg_file <- file.path(out_dir, "deg_table_all.csv")
if (!file.exists(deg_file))
  stop("Missing ", deg_file, " - run 01_diff_analysis.R first (same <out_dir>)")
panel_file <- file.path(data_dir, "druggable_panel.csv")
if (!file.exists(panel_file))
  stop("Missing ", panel_file, " - pass the course1 data dir as the first argument")
deg  <- read_csv(deg_file, show_col_types = FALSE)
panel <- read_csv(panel_file, show_col_types = FALSE)
cat("Deg table:", nrow(deg), "| druggable panel:", nrow(panel), "genes\n")

# ---- 2. Significant DEGs --------------------------------------------------------
sig <- deg %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
cat("Significant DEGs:", nrow(sig), "\n")

# ---- 3. Intersect with druggable panel: THE target shortlist --------------------
candidates <- inner_join(sig, panel, by = "symbol") %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(rank = row_number()) %>%
  select(rank, symbol, panel_class, log2FoldChange, padj, baseMean)

write_csv(candidates, file.path(out_dir, "target_shortlist.csv"))
cat("Saved target_shortlist.csv:", nrow(candidates), "candidate targets\n")

# ---- 4. Pathway members among DEGs (context for discussion) --------------------
pathway_genes <- c("EGFR", "ERBB2", "KRAS", "NRAS", "BRAF", "RAF1",
                   "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "AKT1", "AKT2",
                   "PIK3CA", "PIK3R1", "MTOR", "VEGFA", "MYC", "CCND1")
hits <- sig %>% filter(symbol %in% pathway_genes) %>%
  arrange(desc(log2FoldChange)) %>%
  select(symbol, log2FoldChange, padj)
cat("\n===== EGFR/MAPK pathway genes among DEGs =====\n")
print(hits)

# ---- 5. Console summary ---------------------------------------------------------
cat("\n===== TOP candidate targets (druggable, up-regulated) =====\n")
print(head(candidates, 10))
cat("\n===== Bottom-ranked candidate targets (smallest log2FC) =====\n")
print(tail(candidates, 5))

cat("\nInterpretation hint: an up-regulated druggable kinase (e.g. EGFR) with a\n")
cat("large fold change is a classic candidate: the tumour may depend on it.\n")
cat("Next step (course 2): validate EGFR at protein/structure level.\n")