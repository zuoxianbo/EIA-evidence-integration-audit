#!/usr/bin/env Rscript
# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Supplementary Fig S1 (v25, new): Delta attribution map
# Per critical review section 9: unify the audit results into one attribution landscape.
# x = provenance independence (1 - label-input overlap score)
# y = representation sensitivity (available-case - sentinel AUROC shift, harmonic)
# colour = fair-baseline-adjusted Delta (harmonic - STRING centrality)
# shape = evaluation context (in-sample PDAC vs external)

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
source(file.path(base_dir, "common.R"))

library(ggplot2)

d <- read.csv(file.path(base_dir, "source_data", "supp1_attribution_map.csv"), stringsAsFactors = FALSE)

# small horizontal offsets within the same provenance level (readability only)
xoff <- c("E3-C" = -0.035, "E3" = 0.035, "E6" = 0.0,
          "E1" = -0.045, "E2" = -0.015, "E3-A" = 0.015, "E4" = 0.045, "E5" = 0.075)
d$x <- d$provenance_independence + xoff[as.character(d$endpoint)]
d$context <- factor(d$context, levels = c("in-sample", "external"))

# manual label nudges to avoid overlap
ynudge <- c("E3-C" = 0.006, "E3" = -0.008, "E6" = 0.007,
           "E1" = 0.006, "E2" = 0.007, "E3-A" = -0.008, "E4" = 0.006, "E5" = -0.008)
d$y_lab <- d$repr_shift + ynudge[as.character(d$endpoint)]
xnudge <- c("E3-C" = -0.01, "E3" = 0.012, "E6" = 0.015,
            "E1" = -0.01, "E2" = -0.012, "E3-A" = 0.012, "E4" = 0.012, "E5" = 0.012)
d$x_lab <- d$x + xnudge[as.character(d$endpoint)]

p <- ggplot(d, aes(x = x, y = repr_shift)) +
  # provenance zones
  annotate("rect", xmin = -0.06, xmax = 0.16, ymin = -0.06, ymax = 0.075, fill = "#c0392b", alpha = 0.06) +
  annotate("text", x = 0.05, y = 0.069, label = "provenance compromised\n(label reuses input)", size = 2.8, colour = "#c0392b", family = "sans") +
  annotate("text", x = 0.78, y = 0.069, label = "independent label provenance", size = 2.8, colour = "#404040", family = "sans") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0.33, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  # points
  geom_point(aes(colour = delta_auroc, shape = context), size = 4.2, stroke = 1.1) +
  scale_colour_gradient2(low = "#2e7bd6", mid = "grey85", high = "#c0392b", midpoint = 0,
                          limits = c(-0.13, 0.32),
                          name = expression(Delta * "AUROC")) +
  scale_shape_manual(values = c("in-sample" = 16, "external" = 17), name = NULL) +
  geom_text(aes(x = x_lab, y = y_lab, label = endpoint), size = 2.9, family = "sans", colour = "#202020") +
  # reading guide
  annotate("text", x = 0.20, y = 0.043, label = "positive Delta occurs only in the compromised-provenance zone",
           size = 2.9, family = "sans", colour = "#c0392b", hjust = 0) +
  scale_x_continuous(limits = c(-0.08, 1.13), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  scale_y_continuous(limits = c(-0.062, 0.078), breaks = seq(-0.05, 0.05, 0.025)) +
  labs(x = "provenance independence (1 - label-input overlap)",
       y = "representation sensitivity (AUROC shift under recoding)",
       title = "Supplementary Fig. 1  Delta attribution map across the provenance ladder") +
  theme_nature() +
  theme(plot.title = element_text(hjust = 0, size = 11.5, face = "bold"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(t = -3))

save_fig(p, "Supp1_attribution_map", w = 9.5, h = 6.0)
