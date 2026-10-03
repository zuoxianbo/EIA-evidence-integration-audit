#!/usr/bin/env Rscript
# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 1 (v25 revision): claim -> Delta -> alternative explanations -> audit decision
# Redesigned per critical review Table 4: central spine instead of checklist layout.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))

library(ggplot2)

claim_box <- data.frame(x = 1.9, xend = 8.1, y = 8.6, yend = 9.7)
delta_box <- data.frame(x = 1.9, xend = 8.1, y = 6.7, yend = 7.9)
alt_boxes <- data.frame(
  x = c(0.4, 2.85, 5.3, 7.75), xend = c(2.25, 4.7, 7.15, 9.6),
  y = 3.9, yend = 6.0,
  fill = c("#e9c46a", "#e07a5f", "#2a9d8f", "#4a90e2"),
  name = c("Provenance", "Representation", "Aggregation", "External evaluation"),
  q = c("label reuses an\ninput layer?", "missingness /\nscale coding?", "rule breaks\npairwise order?", "gain survives\nnew context?"),
  test = c("overlap / deletion /\ncircular control", "sentinel vs\navailable-case", "order-preservation\n/ monotonicity", "held-out endpoints\nE4 / E5 / E6")
)
decision_box <- data.frame(x = 1.2, xend = 8.8, y = 0.7, yend = 2.9)

p <- ggplot() +
  # claim box
  geom_rect(data = claim_box, aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "#264653", colour = "black", linewidth = 0.4) +
  annotate("text", x = 5.0, y = 9.35, label = "Integration claim   C = (X, A, Y, B, P, G)",
           colour = "white", size = 4.6, fontface = "bold", family = "sans") +
  annotate("text", x = 5.0, y = 8.9, label = "evidence layers X  |  integrator A  |  endpoint Y  |  baseline B  |  protocol P  |  context G",
           colour = "white", size = 3.2, family = "sans") +
  # delta estimand box
  geom_rect(data = delta_box, aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "#eef3f9", colour = "#2e7bd6", linewidth = 0.7) +
  annotate("text", x = 5.0, y = 7.55, label = "Estimand:   Delta = Performance(A(X)) - Performance(B)",
           size = 4.2, fontface = "bold", family = "sans", colour = "#1a3d5c") +
  annotate("text", x = 5.0, y = 7.08, label = "measured on the same independent evaluation;\nobserved Delta > 0 is not yet an identified information gain",
           size = 3.0, lineheight = 1.05, family = "sans", colour = "#404040") +
  # arrow claim -> delta
  annotate("segment", x = 5.0, xend = 5.0, y = 8.55, yend = 7.95,
           arrow = arrow(length = unit(0.3, "cm")), linewidth = 0.7, colour = "black") +
  # red provenance threat arrow (down the right margin, into the delta box corner)
  annotate("segment", x = 8.15, xend = 7.65, y = 9.1, yend = 8.0,
           arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.5, colour = "#c0392b", linetype = "dashed") +
  annotate("text", x = 8.85, y = 8.55, label = "input -> label\nprovenance", colour = "#c0392b",
           size = 3.0, fontface = "italic", family = "sans") +
  # alternative explanation boxes
  geom_rect(data = alt_boxes, aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = alt_boxes$fill, colour = "black", linewidth = 0.3) +
  geom_text(data = alt_boxes, aes(x = (x + xend)/2, y = 5.65, label = name),
            colour = "white", size = 4.3, fontface = "bold", family = "sans") +
  geom_text(data = alt_boxes, aes(x = (x + xend)/2, y = 4.95, label = q),
            colour = "white", size = 3.0, fontface = "italic", family = "sans", lineheight = 1.0) +
  geom_text(data = alt_boxes, aes(x = (x + xend)/2, y = 4.25, label = test),
            colour = "white", size = 2.8, family = "sans", lineheight = 1.0) +
  # arrows delta -> alternatives (drop from a horizontal connector under the delta box)
  annotate("segment", x = 1.3, xend = 8.65, y = 6.62, yend = 6.62, linewidth = 0.5, colour = "black") +
  geom_segment(data = data.frame(x = c(1.3, 3.75, 6.2, 8.65), xend = c(1.3, 3.75, 6.2, 8.65)),
               aes(x = x, xend = xend, y = 6.62, yend = 6.05),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.5, colour = "black") +
  annotate("text", x = 5.0, y = 6.35, label = "alternative explanations for Delta > 0",
           size = 3.3, fontface = "italic", family = "sans", colour = "#404040") +
  # decision box
  geom_rect(data = decision_box, aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "#f7f7f7", colour = "black", linewidth = 0.5) +
  annotate("text", x = 5.0, y = 2.62, label = "Audit decision",
           size = 4.0, fontface = "bold", family = "sans") +
  annotate("text", x = 5.0, y = 2.22, label = "audit-condition semantics (not model-performance grades)",
           size = 2.6, family = "sans", colour = "#404040") +
  annotate("text", x = 1.9, y = 1.55, label = "PASS", colour = "#2a9d8f", size = 4.2, fontface = "bold", family = "sans") +
  annotate("text", x = 1.9, y = 1.1, label = "audit condition\nsatisfied", size = 2.6, family = "sans") +
  annotate("text", x = 5.0, y = 1.55, label = "WARN", colour = "#e07a5f", size = 4.2, fontface = "bold", family = "sans") +
  annotate("text", x = 5.0, y = 1.1, label = "evidence insufficient\nor sensitive", size = 2.6, family = "sans") +
  annotate("text", x = 8.1, y = 1.55, label = "FAIL", colour = "#c0392b", size = 4.2, fontface = "bold", family = "sans") +
  annotate("text", x = 8.1, y = 1.1, label = "attribution condition\nviolated", size = 2.6, family = "sans") +
  # arrows alternatives -> decision
  geom_segment(data = data.frame(x = c(1.3, 3.75, 6.2, 8.65), xend = c(2.6, 4.6, 5.4, 7.4)),
               aes(x = x, xend = xend, y = 3.85, yend = 2.95),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.5, colour = "black") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.4, 10.0), clip = "off") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_nature() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  labs(title = "Fig. 1 Evidence Integration Audit: from an integration claim to an audit decision through the estimand Delta")

save_fig(p, "Fig1_EIA_framework", w = 10.5, h = 7.0)
