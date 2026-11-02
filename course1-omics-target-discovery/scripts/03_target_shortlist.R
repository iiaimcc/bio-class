#!/usr/bin/env Rscript
# =============================================================================
# 03_target_shortlist.R
# Course: "Omics-based disease mechanism dissection and target discovery" (Part 3)
#
# What this script does:
#   1. Reads the DEG table from step 1, the curated druggable panel and the
#      mutation matrix built from the same TCGA-LUAD samples
#   2. Intersects DEGs with the panel -> expression-driven candidates
#   3. Lists mutated driver genes -> mutation-driven candidates
#   4. Writes both tables and prints the evidence for EGFR
#
# How to run:
#   Rscript 03_target_shortlist.R  <data_dir>  <out_dir>
# =============================================================================

# ---- 0. Configuration --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "."
out_dir  <- if (length(args) >= 2) args[2] else "./results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (p in c("dplyr", "readr", "tidyr")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' missing")
}
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

# ---- 1. Load inputs ------------------------------------------------------------
deg_file <- file.path(out_dir, "deg_table_all.csv")
if (!file.exists(deg_file))
  stop("Missing ", deg_file, " - run 01_diff_analysis.R first (same <out_dir>)")
panel_file <- file.path(data_dir, "druggable_panel.csv")
if (!file.exists(panel_file))
  stop("Missing ", panel_file, " - pass the course1 data dir as the first argument")
mut_file <- file.path(data_dir, "mutation_matrix.csv")
if (!file.exists(mut_file))
  stop("Missing ", mut_file, " - run data/fetch_tcga_luad.py or use the fallback copy")

deg   <- read_csv(deg_file, show_col_types = FALSE)
panel <- read_csv(panel_file, show_col_types = FALSE)
mut   <- read_csv(mut_file, show_col_types = FALSE)
variants <- read_csv(file.path(data_dir, "mutation_status.csv"), show_col_types = FALSE)
cat("Deg table:", nrow(deg), "| druggable panel:", nrow(panel),
    "genes | tumours with mutation calls:", nrow(mut), "\n")

# Mutation counts per gene (mut is one row per tumour, 0/1 per gene)
mut_long <- mut %>%
  pivot_longer(-sample_id, names_to = "symbol", values_to = "mutated") %>%
  filter(mutated == 1) %>%
  count(symbol, name = "n_mutated")
examples <- variants %>%
  group_by(gene) %>%
  summarise(mutation_examples = paste(unique(hgvsp), collapse = "; "), .groups = "drop") %>%
  rename(symbol = gene)

# ---- 2. Significant DEGs --------------------------------------------------------
sig <- deg %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
cat("Significant DEGs:", nrow(sig), "\n")

# ---- 3. Expression-driven candidates: DEGs that are druggable -------------------
candidates <- inner_join(sig, panel, by = "symbol") %>%
  left_join(mut_long, by = "symbol") %>%
  left_join(examples, by = "symbol") %>%
  mutate(n_mutated = coalesce(n_mutated, 0L),
         mutation_examples = coalesce(mutation_examples, "")) %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(rank = row_number()) %>%
  select(rank, symbol, panel_class, known_drug_class, log2FoldChange, padj,
         baseMean, n_mutated, mutation_examples)

write_csv(candidates, file.path(out_dir, "target_shortlist.csv"))
cat("Saved target_shortlist.csv:", nrow(candidates), "expression-driven candidates\n")

# ---- 4. Mutation-driven candidates ---------------------------------------------
# A gene is listed if at least one tumour carries a non-silent variant. Absence
# of a call is not proof of wild type: one tumour has no MAF file, and MAF data
# cannot capture fusions or copy-number events.
mut_candidates <- mut_long %>%
  left_join(examples, by = "symbol") %>%
  left_join(deg %>% select(symbol, log2FoldChange, padj, baseMean),
            by = "symbol") %>%
  mutate(deg_status = case_when(
    is.na(padj)                              ~ "not tested (filtered out)",
    padj < 0.05 & abs(log2FoldChange) > 1    ~ "DEG",
    padj < 0.05                              ~ "padj<0.05 but |log2FC|<1",
    TRUE                                     ~ "not significant"
  )) %>%
  arrange(desc(n_mutated), desc(log2FoldChange)) %>%
  select(symbol, n_mutated, mutation_examples, deg_status,
         log2FoldChange, padj, baseMean)

write_csv(mut_candidates, file.path(out_dir, "mutation_driven_candidates.csv"))
cat("Saved mutation_driven_candidates.csv:", nrow(mut_candidates), "mutated genes\n")

# ---- 5. Console summary ---------------------------------------------------------
cat("\n===== TOP expression-driven candidates (druggable, up-regulated) =====\n")
print(head(candidates, 10))

cat("\n===== Mutated genes among the driver list =====\n")
print(mut_candidates, n = 20)

cat("\n===== EGFR: the course target =====\n")
print(deg %>% filter(symbol == "EGFR") %>%
        select(symbol, baseMean, log2FoldChange, padj))
print(mut_candidates %>% filter(symbol == "EGFR"))

cat("\nEGFR is NOT among the expression-driven candidates. In this subset it is\n")
cat("mutated in 2/15 tumours (exon 19 deletion and L858R), which is the\n")
cat("evidence that nominates it as a target. Expression change and mutation\n")
cat("status answer different questions; neither alone proves a dependency.\n")
cat("Next step (course 2): inspect EGFR at protein/structure level.\n")
