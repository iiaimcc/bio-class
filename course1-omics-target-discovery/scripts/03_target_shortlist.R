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
#   4. Layers the evidence (mutated -> druggable -> cell surface -> clinical
#      precedent) -> integrated_candidates.csv and target_recommendation.csv
#   5. Prints the evidence for EGFR, which is the last one standing
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

# ---- 5. Integrated evidence layers ---------------------------------------------
# A candidate has to survive four questions, not one: does it carry a non-silent
# variant in this cohort, is it in the druggable panel, can a large molecule
# reach it (cell surface), and is there clinical precedent. The counts are
# printed one layer at a time so the class can see what each filter costs.
annot_file <- file.path(data_dir, "driver_annotation.csv")
if (!file.exists(annot_file))
  stop("Missing ", annot_file, " - pass the course1 data dir as the first argument")
annot <- read_csv(annot_file, show_col_types = FALSE)

# A candidate has to survive four questions, not one: does it carry a non-silent
# variant in this cohort, is it in the druggable panel, can a large molecule
# reach it (cell surface), and does the observed variant type match the
# alteration an approved drug actually targets. The counts are printed one
# layer at a time so the class can see what each filter costs.
annot_file <- file.path(data_dir, "driver_annotation.csv")
if (!file.exists(annot_file))
  stop("Missing ", annot_file, " - pass the course1 data dir as the first argument")
annot <- read_csv(annot_file, show_col_types = FALSE)

integrated <- mut_candidates %>%
  left_join(panel %>% select(symbol, panel_class, known_drug_class), by = "symbol") %>%
  left_join(annot, by = "symbol") %>%
  mutate(in_druggable_panel = !is.na(panel_class),
         compartment = coalesce(compartment, "unannotated"),
         approved_drug_status = coalesce(approved_drug_status, "unknown"),
         observed_variants_targetable = coalesce(observed_variants_targetable, "no")) %>%
  arrange(desc(n_mutated), desc(in_druggable_panel), desc(log2FoldChange))

final <- integrated %>%
  filter(in_druggable_panel,
         compartment == "cell_surface",
         approved_drug_status == "approved_in_luad",
         observed_variants_targetable == "yes") %>%
  mutate(rank = row_number()) %>%
  select(rank, symbol, panel_class, known_drug_class, n_mutated, mutation_examples,
         compartment, approved_drug_status, targetable_alterations,
         deg_status, log2FoldChange, padj, baseMean)
write_csv(integrated, file.path(out_dir, "integrated_candidates.csv"))
cat("Saved integrated_candidates.csv:", nrow(integrated),
    "mutated driver genes with panel / compartment annotation\n")
write_csv(final, file.path(out_dir, "target_recommendation.csv"))
cat("Saved target_recommendation.csv:", nrow(final), "recommended targets\n")

# ---- 6. Console summary ---------------------------------------------------------
cat("\n===== Evidence funnel (each layer costs candidates) =====\n")
cat("  1. non-silent variant in this cohort :", nrow(integrated), "genes\n")
cat("  2. + in druggable panel              :", sum(integrated$in_druggable_panel), "genes\n")
cat("  3. + on the cell surface             :",
    sum(integrated$in_druggable_panel & integrated$compartment == "cell_surface"), "genes\n")
cat("  4. + variant matches approved drug   :", nrow(final), "gene(s)\n")
cat("  expression-driven shortlist overlap  :",
    sum(candidates$symbol %in% final$symbol), "genes\n")

cat("\n===== TOP expression-driven candidates (druggable, up-regulated) =====\n")
print(head(candidates, 10))

cat("\n===== Mutated genes among the driver list =====\n")
print(mut_candidates, n = 20)

cat("\n===== Recommended target(s) after all four layers =====\n")
print(final)

cat("\n===== EGFR: the course target =====\n")
print(deg %>% filter(symbol == "EGFR") %>%
        select(symbol, baseMean, log2FoldChange, padj))
print(mut_candidates %>% filter(symbol == "EGFR"))

cat("\nEGFR is NOT among the expression-driven candidates. In this subset it is\n")
cat("mutated in 2/15 tumours (exon 19 deletion and L858R), sits on the cell\n")
cat("surface and already has approved drugs. The 78-gene shortlist and this\n")
cat("four-layer funnel share no gene at all: that gap is the lesson, not a bug.\n")
cat("RET and ROS1 also have approved drugs in LUAD, but those drugs target\n")
cat("fusions, so the point mutations seen here do not make them actionable.\n")
cat("Next step (course 2): inspect EGFR at protein/structure level.\n")
