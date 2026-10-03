#!/usr/bin/env Rscript
# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 7 (v26 numbering):: EIA Toolkit + REAL claim-level audit scorecard
# Per critical review Table 4: replace the schematic scorecard with the actual
# audit verdicts and quantitative findings of this study (fig6_scorecard.csv).

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
source(file.path(base_dir, "common.R"))

library(ggplot2)

# ---- Top: toolkit modules + audit call ----
modules <- data.frame(
  x = c(0.4, 2.9, 5.4, 7.9), xend = c(2.6, 5.1, 7.6, 10.1),
  y = 7.6, yend = 8.9,
  fill = c("#e9c46a", "#e07a5f", "#2a9d8f", "#4a90e2"),
  line1 = c("eia_toolkit.provenance", "eia_toolkit.representation", "eia_toolkit.aggregation", "eia_toolkit.generalization"),
  line2 = c("overlap / DAG / deletion", "sentinel / available-case", "order / monotonicity", "E4 / E5 / E6 held-out")
)

sc <- read.csv(file.path(base_dir, "source_data", "fig6_scorecard.csv"), stringsAsFactors = FALSE)
sc$axis <- factor(sc$axis, levels = rev(sc$axis))

verdict_col <- c("FAIL" = "#c0392b", "WARN" = "#e07a5f", "PASS (audit condition satisfied)" = "#2a9d8f")
sc$vcol <- verdict_col[sc$verdict]

p <- ggplot() +
  # module boxes
  geom_rect(data = modules, aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = modules$fill, colour = "black", linewidth = 0.3) +
  geom_text(data = modules, aes(x = (x + xend)/2, y = 8.5, label = line1),
            colour = "white", size = 3.3, fontface = "bold", family = "sans") +
  geom_text(data = modules, aes(x = (x + xend)/2, y = 8.0, label = line2),
            colour = "white", size = 2.7, family = "sans") +
  # audit call box
  annotate("rect", xmin = 1.0, xmax = 9.5, ymin = 6.3, ymax = 7.2, fill = "#eef3f9", colour = "#2e7bd6", linewidth = 0.5) +
  annotate("text", x = 5.25, y = 6.75, label = "eia_audit(data, endpoints)  ->  claim-level scorecard + report + checklist",
           size = 4.0, family = "sans", fontface = "bold", colour = "#1a3d5c") +
  annotate("segment", x = 5.25, xend = 5.25, y = 7.55, yend = 7.25,
           arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.6, colour = "black") +
  annotate("segment", x = 5.25, xend = 5.25, y = 6.25, yend = 5.85,
           arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.6, colour = "black") +
  # scorecard frame
  annotate("rect", xmin = 0.3, xmax = 10.2, ymin = 0.15, ymax = 5.75, colour = "black", linewidth = 0.5, fill = NA) +
  annotate("text", x = 5.25, y = 5.45, label = "AUDIT SCORECARD  -  PDAC / CRC case study (all numbers from this study; reproducible via example_pdac_audit.py)",
           size = 3.4, fontface = "bold", family = "sans", hjust = 0.5) +
  # verdict chips
  geom_rect(data = sc, aes(xmin = 0.55, xmax = 1.5, ymin = as.numeric(axis) - 0.32, ymax = as.numeric(axis) + 0.32),
            fill = sc$vcol, colour = "black", linewidth = 0.3) +
  geom_text(data = sc, aes(x = 1.025, y = as.numeric(axis),
                          label = c("FAIL", "WARN", "WARN", "PASS")),
            colour = "white", size = 3.0, fontface = "bold", family = "sans") +
  # axis names
  geom_text(data = sc, aes(x = 1.75, y = as.numeric(axis), label = axis),
            size = 3.5, fontface = "bold", family = "sans", hjust = 0) +
  # headlines
  geom_text(data = sc, aes(x = 4.1, y = as.numeric(axis) + 0.16, label = headline),
            size = 3.0, family = "sans", hjust = 0) +
  # details
  geom_text(data = sc, aes(x = 4.1, y = as.numeric(axis) - 0.22, label = detail),
            size = 2.5, family = "sans", hjust = 0, colour = "#404040") +
  coord_cartesian(xlim = c(0, 10.5), ylim = c(0, 9.3), clip = "off") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_nature() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  labs(title = "Extended Data Fig. 1 The EIA Toolkit produces a claim-level audit scorecard from the four audit modules")

save_fig(p, "ExtendedDataFig1_toolkit", w = 11, h = 6.0)
