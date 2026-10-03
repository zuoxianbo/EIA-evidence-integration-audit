#!/usr/bin/env Rscript
# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Master script to regenerate all 12 EIA figures (Fig 1-6 + ED Fig 1-6) using
# R/ggplot2. Each figure script is self-contained and sources common.R.
# Output base names match the existing figures/ files.

base_dir <- file.path(ROOT, "code", "R")
setwd(base_dir)

cat("Regenerating all EIA figures (v39) with R/ggplot2...\n")
cat("Base directory:", base_dir, "\n")

fig_scripts <- c(
  "fig1_framework.R",
  "fig2_provenance.R",
  "fig3_representation_aggregation.R",
  "fig4_baseline_delta.R",
  "fig5_disease_alignment.R",
  "fig6_interpretation_matrix.R",
  "ED_Fig1_benchmark_construction.R",
  "ED_Fig2_provenance_graph.R",
  "ED_Fig3_representation_sensitivity.R",
  "ED_Fig4_missingness_structure.R",
  "ED_Fig5_aggregation_sensitivity.R",
  "ED_Fig6_evaluation_sensitivity.R"
)

for (s in fig_scripts) {
  cat("\n--- Running", s, "---\n")
  source(file.path(base_dir, s))
}

cat("\nAll 12 figures regenerated.\n")
cat("Outputs in:", file.path(dirname(dirname(base_dir)), "figures"), "\n")
