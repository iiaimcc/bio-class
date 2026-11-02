#!/usr/bin/env Rscript
# =============================================================================
# 01_diff_analysis.R
# Course: "Omics-based disease mechanism dissection and target discovery" (Part 1)
#
# What this script does:
#   1. Loads the real TCGA-LUAD count matrix, sample annotation and gene annotation
#   2. Filters low-count genes, runs DESeq2 (Tumor vs adjacent Normal)
#   3. Writes the full result table and a volcano plot (ggplot2)
#
# Data: TCGA-LUAD STAR counts, 15 primary tumours + 15 solid tissue normals,
#       downloaded from the GDC API by data/fetch_tcga_luad.py.
#
# How to run (in RStudio, step through, or terminal):
#   Rscript 01_diff_analysis.R  <data_dir>  <out_dir>
#   e.g.  Rscript scripts/01_diff_analysis.R ../data ./results
#   With no arguments the script reads the current directory and writes ./results.
#
# Note: the script is written for teaching - read the comments line by line.
# =============================================================================

# ---- 0. Configuration (edit these two paths if needed) ----------------------
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "."
out_dir  <- if (length(args) >= 2) args[2] else "./results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 0.5 Required packages ---------------------------------------------------
pkgs <- c("DESeq2", "ggplot2", "dplyr", "tidyr", "readr", "tibble")
for (p in pkgs) {
  suppressPackageStartupMessages(if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is missing. Install it first, see 环境与部署指南.md"))
}
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ---- 1. Load data ------------------------------------------------------------
# expression_matrix.csv: rows = genes (Ensembl gene ID), columns = samples
counts <- read.csv(file.path(data_dir, "expression_matrix.csv"),
                   row.names = 1, check.names = FALSE)
# sample_meta.csv: one row per sample with condition (Normal/Tumor)
meta <- read.csv(file.path(data_dir, "sample_meta.csv"))
# gene_annotation.csv: Ensembl gene ID -> HGNC symbol, gene type
ann <- read.csv(file.path(data_dir, "gene_annotation.csv"))

stopifnot(all(colnames(counts) == meta$sample_id))
rownames(meta) <- meta$sample_id  # DESeq2 matches colData rows to count columns
cat("Loaded", nrow(counts), "genes x", ncol(counts), "samples\n")
print(table(meta$condition))

# ---- 2. Low-count filter -----------------------------------------------------
# A gene must reach >=10 counts in at least 6 samples (the size of one group).
# This removes genes that cannot support a stable dispersion estimate, and
# keeps the multiple-testing burden close to the expressed transcriptome.
keep <- rowSums(counts >= 10) >= 6
cat("Genes passing the low-count filter:", sum(keep), "/", nrow(counts), "\n")
counts <- counts[keep, , drop = FALSE]

# ---- 3. DESeq2 pipeline ------------------------------------------------------
# 3.1 Level "Normal" is the reference (baseline) group
meta$condition <- relevel(factor(meta$condition), ref = "Normal")

# 3.2 Design: condition only. The 30 samples come from different patients and
#     several tissue source sites, so a paired design is not available here and
#     TSS is not modelled (too many levels for 30 samples). This is a limitation
#     to discuss, not something the model fixes.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ condition
)

# 3.3 Estimate size factors, dispersions, fit the model (< 2 min typically).
dds <- DESeq(dds)

# 3.4 Extract the Tumor vs Normal contrast
res <- results(dds, contrast = c("condition", "Tumor", "Normal"),
               alpha = 0.05)

# 3.5 Attach gene symbols, tidy result table
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>%
  left_join(ann, by = "gene_id") %>%
  arrange(padj) %>%
  select(gene_id, symbol, gene_type, baseMean, log2FoldChange,
         lfcSE, stat, pvalue, padj)

# ---- 4. Write full result table ----------------------------------------------
deg_all <- file.path(out_dir, "deg_table_all.csv")
write_csv(res_df, deg_all)
cat("Saved full result table:", deg_all, "\n")

# ---- 5. Volcano plot (ggplot2) ----------------------------------------------
# A volcano plot shows significance (-log10 padj) vs. effect size (log2FC)
plt <- res_df %>%
  mutate(significance = case_when(
    padj < 0.05 & log2FoldChange > 1   ~ "up",
    padj < 0.05 & log2FoldChange < -1  ~ "down",
    TRUE                               ~ "ns"
  )) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c(down = "#2166ac", up = "#b2182b", ns = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "TCGA-LUAD: tumour vs adjacent normal",
       subtitle = "15 primary tumours vs 15 solid tissue normals (STAR counts)",
       x = "log2 fold change (Tumor/Normal)",
       y = "-log10(adjusted p-value)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "volcano_plot.png"), plt, width = 7, height = 5.5, dpi = 300)
cat("Saved volcano plot:", file.path(out_dir, "volcano_plot.png"), "\n")

# ---- 6. Quick summary printed to console -------------------------------------
sig <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
cat("\n===== Summary =====\n")
cat("Significant DEGs (padj<0.05, |log2FC|>1):", nrow(sig), "\n")
cat("  up-regulated:", sum(sig$log2FoldChange > 0),
    "| down-regulated:", sum(sig$log2FoldChange < 0), "\n")
cat("\nTop 10 up-regulated (by log2FC):\n")
print(head(sig %>% arrange(desc(log2FoldChange)) %>%
             select(symbol, log2FoldChange, padj), 10))
cat("\nTop 5 down-regulated:\n")
print(head(sig %>% arrange(log2FoldChange) %>%
             select(symbol, log2FoldChange, padj), 5))

# EGFR is the course target: print it even though it is not a DEG.
cat("\nEGFR (course target, for discussion):\n")
print(res_df %>% filter(symbol == "EGFR") %>%
        select(gene_id, baseMean, log2FoldChange, pvalue, padj))
cat("A large padj here is the lesson: expression change alone does not\n")
cat("nominate this target. See 03_target_shortlist.R and the mutation table.\n")
cat("\nDone. Next: run 02_enrichment.R with the DEG list.\n")
