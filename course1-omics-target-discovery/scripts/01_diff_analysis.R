#!/usr/bin/env Rscript
# =============================================================================
# 01_diff_analysis.R
# Course: "Omics-based disease mechanism dissection and target discovery" (Part 1)
#
# What this script does:
#   1. Loads the count matrix and sample annotation (CSV)
#   2. Runs DESeq2 differential expression analysis (Tumor vs Normal)
#   3. Writes the full result table and saves a volcano plot (ggplot2)
#
# How to run (in RStudio, step through, or terminal):
#   Rscript 01_diff_analysis.R  <data_dir>  <out_dir>
#   e.g.  Rscript scripts/01_diff_analysis.R . ./results   (from the course1 dir)
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
# expression_matrix.csv: rows = genes, columns = samples, first col = gene symbol
counts <- read.csv(file.path(data_dir, "expression_matrix.csv"),
                   row.names = 1, check.names = FALSE)
# sample_meta.csv: one row per sample with condition (Normal/Tumor) and batch
meta <- read.csv(file.path(data_dir, "sample_meta.csv"))

stopifnot(all(colnames(counts) == meta$sample_id))
rownames(meta) <- meta$sample_id  # DESeq2 matches colData rows to count columns
cat("Loaded", nrow(counts), "genes x", ncol(counts), "samples\n")
print(table(meta$condition, meta$batch))

# ---- 2. DESeq2 pipeline ------------------------------------------------------
# 2.1 Level "Normal" is the reference (baseline) group
meta$condition <- relevel(factor(meta$condition), ref = "Normal")

# 2.2 Build the DESeqDataSet. Design: condition is what we test;
#     batch is a known confounder that is modeled (not of interest).
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ batch + condition
)

# 2.3 Estimate size factors (library depth), dispersions, fit the model.
#     This is where most of the computation happens (< 2 min typically).
dds <- DESeq(dds)

# 2.4 Extract the Tumor vs Normal contrast
res <- results(dds, contrast = c("condition", "Tumor", "Normal"),
               alpha = 0.05)

# 2.5 Tidy result table
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "symbol") %>%
  arrange(padj)

# ---- 3. Write full result table ----------------------------------------------
deg_all <- file.path(out_dir, "deg_table_all.csv")
write_csv(res_df, deg_all)
cat("Saved full result table:", deg_all, "\n")

# ---- 4. Volcano plot (ggplot2) ----------------------------------------------
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
  labs(title = "Volcano plot: LUAD vs Normal (simulated dataset)",
       x = "log2 fold change (Tumor/Normal)",
       y = "-log10(adjusted p-value)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "volcano_plot.png"), plt, width = 7, height = 5.5, dpi = 300)
cat("Saved volcano plot:", file.path(out_dir, "volcano_plot.png"), "\n")

# ---- 5. Quick summary printed to console -------------------------------------
sig <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
cat("\n===== Summary =====\n")
cat("Significant DEGs (padj<0.05, |log2FC|>1):", nrow(sig), "\n")
cat("  up-regulated:", sum(sig$log2FoldChange > 0),
    "| down-regulated:", sum(sig$log2FoldChange < 0), "\n")
cat("\nTop 10 up-regulated (by log2FC):\n")
print(head(sig %>% arrange(desc(log2FoldChange)) %>% select(symbol, log2FoldChange, padj), 10))
cat("\nTop 5 down-regulated:\n")
print(head(sig %>% arrange(log2FoldChange) %>% select(symbol, log2FoldChange, padj), 5))
cat("\nDone. Next: run 02_enrichment.R with the DEG list.\n")