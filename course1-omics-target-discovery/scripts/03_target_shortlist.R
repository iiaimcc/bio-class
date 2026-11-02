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
#   4. Applies the biology-only funnel (variant present -> protein switched on
#      -> reachable by the large molecule we design in course 3) -> writes
#      integrated_candidates.csv and target_recommendation.csv
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

# ---- 5. Biology-only evidence layers -------------------------------------------
# The funnel deliberately uses ONLY properties of the tumour and of the molecule:
# does this gene carry a non-silent variant in this cohort, does that variant
# switch the protein ON (gain of function), and is the protein somewhere a large
# molecule can reach it. Drug-landscape facts (does an approved drug already
# exist, is the gene in the druggable panel) are NOT filters: using them would
# rebuild the answer we already know and call it a discovery. They are printed
# afterwards, as precedent, in section 6.
annot_file <- file.path(data_dir, "driver_annotation.csv")
if (!file.exists(annot_file))
  stop("Missing ", annot_file, " - pass the course1 data dir as the first argument")
annot <- read_csv(annot_file, show_col_types = FALSE)

integrated <- mut_candidates %>%
  left_join(panel %>% select(symbol, known_drug_class), by = "symbol") %>%
  left_join(annot, by = "symbol") %>%
  mutate(on_druggable_panel = !is.na(known_drug_class),
         compartment = coalesce(compartment, "unannotated"),
         variant_function = coalesce(variant_function, "unknown"),
         approved_drug_status = coalesce(approved_drug_status, "unknown"),
         gain_of_function = variant_function == "gain_of_function",
         reachable_by_large_molecule = compartment == "cell_surface") %>%
  arrange(desc(n_mutated), desc(gain_of_function), desc(log2FoldChange))

# Layer 3 encodes the modality we chose for course 3 (a designed mini-protein /
# antibody-like binder), not a statement about target quality: an intracellular
# protein is a perfectly good target if you are building a small molecule.
final <- integrated %>%
  filter(gain_of_function, reachable_by_large_molecule) %>%
  mutate(rank = row_number()) %>%
  select(rank, symbol, variant_function, compartment, n_mutated, mutation_examples,
         deg_status, log2FoldChange, padj, baseMean, known_drug_class,
         approved_drug_status, note)
write_csv(integrated, file.path(out_dir, "integrated_candidates.csv"))
cat("Saved integrated_candidates.csv:", nrow(integrated),
    "mutated genes with variant-function and compartment annotation\n")
write_csv(final, file.path(out_dir, "target_recommendation.csv"))
cat("Saved target_recommendation.csv:", nrow(final), "recommended targets\n")

# ---- 6. Console summary ---------------------------------------------------------
cat("\n===== Evidence funnel (biology only, each layer costs candidates) =====\n")
cat("  1. non-silent variant in this cohort  :", nrow(integrated), "genes\n")
cat("  2. + variant switches the protein ON  :", sum(integrated$gain_of_function), "genes\n")
cat("  3. + reachable by a large molecule    :", nrow(final), "gene(s)\n")
cat("  expression-driven shortlist overlap   :",
    sum(candidates$symbol %in% final$symbol), "genes\n")

cat("\n---- Not part of the funnel: what the same data looks like for the other modality ----\n")
small_mol <- integrated %>% filter(gain_of_function, compartment == "intracellular")
cat("Switched-on AND intracellular (= small-molecule route):",
    paste(small_mol$symbol, collapse = ", "), "\n")

cat("\n---- Not part of the funnel: drug-landscape facts, printed as precedent only ----\n")
print(integrated %>% filter(gain_of_function) %>%
        select(symbol, compartment, variant_function, approved_drug_status, known_drug_class))

cat("\n===== TOP expression-driven candidates (druggable, up-regulated) =====\n")
print(head(candidates, 10))

cat("\n===== Mutated genes among the driver list =====\n")
print(mut_candidates, n = 20)

cat("\n===== Recommended target(s) after all three biology layers =====\n")
print(final)

cat("\n===== EGFR: the course target =====\n")
print(deg %>% filter(symbol == "EGFR") %>%
        select(symbol, baseMean, log2FoldChange, padj))
print(mut_candidates %>% filter(symbol == "EGFR"))

cat("\nEGFR is NOT among the expression-driven candidates. Here it is mutated in\n")
cat("2/15 tumours, both textbook activating events (exon 19 deletion, L858R),\n")
cat("and it sits on the cell surface. The expression shortlist and this funnel\n")
cat("share no gene at all: that gap is the lesson, not a bug.\n")
cat("KRAS is switched on in 6/15 tumours too, but it works inside the cell, so\n")
cat("the same data points a small-molecule programme at KRAS instead. Which\n")
cat("gene wins depends on the modality you committed to, not on a ranking.\n")
cat("Next step (course 2): inspect the recommended target at protein level.\n")
