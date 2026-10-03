#!/usr/bin/env Rscript
# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 3: Provenance manufactures apparent gains (Principle 1)

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))

library(ggplot2)
library(patchwork)
library(grid)

# helper to turn CSV escaped \n into real newlines
unesc <- function(x) gsub("\\\\n", "\n", x)

# ---- Panel a: provenance gradient (dual y-axis line) ----
ga <- read.csv(file.path(dat_dir, "fig2a_provenance_gradient.csv"), stringsAsFactors = FALSE)
# map centrality AUROC to the left-axis scale so it can be drawn on the same plot
ga$centrality_left <- ga$centrality_auroc - 0.66

pa <- ggplot(ga, aes(x = theta)) +
  # left axis: delta AUROC
  geom_line(aes(y = delta_auroc, colour = "Delta AUROC"), linewidth = 1) +
  geom_point(aes(y = delta_auroc, colour = "Delta AUROC"), size = 3) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "black", linewidth = 0.5) +
  # right-axis: centrality AUROC (mapped to left scale)
  geom_line(aes(y = centrality_left, colour = "centrality AUROC"), linetype = "dashed", linewidth = 1) +
  geom_point(aes(y = centrality_left, colour = "centrality AUROC", shape = "centrality AUROC"), size = 3) +
  scale_y_continuous(
    name = expression(Delta * "AUROC"),
    limits = c(-0.1, 0.35),
    breaks = seq(-0.1, 0.3, 0.1),
    sec.axis = sec_axis(~ . + 0.66, name = "centrality AUROC",
                        breaks = seq(0.58, 0.70, 0.04))
  ) +
  scale_colour_manual(values = c("Delta AUROC" = "#2e7bd6", "centrality AUROC" = "#e07a5f"), name = NULL) +
  scale_shape_manual(values = c("centrality AUROC" = 15), guide = "none") +
  labs(x = expression(paste("provenance leak ", theta)),
       title = "a  Provenance gradient") +
  theme_nature() +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"),
        axis.title.y.left = element_text(colour = "#2e7bd6"),
        axis.text.y.left = element_text(colour = "#2e7bd6"),
        axis.title.y.right = element_text(colour = "#e07a5f"),
        axis.text.y.right = element_text(colour = "#e07a5f"),
        legend.position = "none")

# ---- Panel b: E3 layer-deletion ----
gb <- read.csv(file.path(dat_dir, "fig2b_e3_deletion.csv"), stringsAsFactors = FALSE)
gb$condition <- unesc(gb$condition)
gb$condition <- factor(gb$condition, levels = gb$condition)

delta_b <- round(gb$value[1] - gb$value[2], 3)
pb <- ggplot(gb, aes(x = condition, y = value, fill = condition)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("full (E3)" = "#2e7bd6",
                               "remove label-embedded druggability" = "#e07a5f"), guide = "none") +
  scale_x_discrete(labels = c("full (E3)" = "full\n(E3)",
                              "remove label-embedded druggability" = "remove\nlabel-embedded\ndruggability")) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "Harmonic AUROC", title = "b  E3 layer-deletion") +
  annotate("text", x = 1.95, y = 0.72, label = paste0("DeLong P < 1e-300 (exact) | Delta = ", delta_b), colour = "#c0392b", size = 2.7, lineheight = 1.0) +
  annotate("segment", x = 1.75, xend = 1.98, y = 0.64, yend = 0.60, arrow = arrow(length = unit(0.15, "cm")), colour = "#c0392b", linewidth = 0.5) +
  theme_nature() +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"),
        axis.text.x = element_text(size = 9.5, lineheight = 0.9))

# ---- Panel c: E3-C circular control ----
gc <- read.csv(file.path(dat_dir, "fig2c_e3c_control.csv"), stringsAsFactors = FALSE)
gc$condition <- unesc(gc$condition)
gc$condition <- factor(gc$condition, levels = gc$condition)

delta_c <- round(gc$value[1] - gc$value[2], 3)
pc <- ggplot(gc, aes(x = condition, y = value, fill = condition)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("Harmonic (E3-C)" = "#2e7bd6",
                               "Harmonic no-druggability" = "#e07a5f"), guide = "none") +
  scale_x_discrete(labels = c("Harmonic (E3-C)" = "Harmonic\n(E3-C)",
                              "Harmonic no-druggability" = "Harmonic\nno-druggability")) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "Harmonic AUROC", title = "c  Circular control (E3-C)") +
  annotate("text", x = 1.5, y = 0.78, label = paste0("Delta = +", delta_c), colour = "#c0392b", size = 3.2) +
  annotate("segment", x = 1.35, xend = 1.85, y = 0.70, yend = 0.60, arrow = arrow(length = unit(0.15, "cm")), colour = "#c0392b", linewidth = 0.5) +
  theme_nature() +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"),
        axis.text.x = element_text(size = 9.5, lineheight = 0.9))

# ---- Panel d: mechanism 2 (v26) - label reuse + constructed interaction label ----
gd <- read.csv(file.path(dat_dir, "fig2d_mechanism2.csv"), stringsAsFactors = FALSE)
gd$mech <- ifelse(grepl("E7", gd$endpoint), "E7  OT genetics\nlabel reuse", "E8  constructed\ninteraction label")

gd_long <- do.call(rbind, lapply(seq_len(nrow(gd)), function(i) {
  rows <- list(
    data.frame(mech = gd$mech[i], what = "best single layer",
               auroc = gd$best_single_auroc[i]),
    data.frame(mech = gd$mech[i], what = "harmonic mean",
               auroc = gd$harmonic_auroc[i]),
    data.frame(mech = gd$mech[i], what = "best fixed composite",
               auroc = gd$best_composite_auroc[i]))
  rfv <- gd$rf_oof_auroc[i]
  if (!is.na(rfv) && rfv != "" && nzchar(as.character(rfv))) {
    rows[[length(rows) + 1]] <- data.frame(mech = gd$mech[i],
                                           what = "random forest (OOF)",
                                           auroc = as.numeric(rfv))
  }
  do.call(rbind, rows)
}))
gd_long$what <- factor(gd_long$what, levels = c("best single layer", "harmonic mean",
                                                "best fixed composite", "random forest (OOF)"))
gd_long$lab <- sprintf("%.3f", gd_long$auroc)

pd <- ggplot(gd_long, aes(x = what, y = auroc, fill = what)) +
  geom_bar(stat = "identity", width = 0.62, colour = "black", linewidth = 0.25) +
  facet_wrap(~ mech, ncol = 2, scales = "free_x") +
  geom_text(aes(label = lab), vjust = -0.5, size = 2.6, family = "sans") +
  scale_fill_manual(values = c("best single layer" = "#b8c4d4",
                               "harmonic mean" = "#2e7bd6",
                               "best fixed composite" = "#5f7d95",
                               "random forest (OOF)" = "#c62828"), guide = "none") +
  scale_x_discrete(labels = c("best single layer" = "best single\nlayer",
                              "harmonic mean" = "harmonic\nmean",
                              "best fixed composite" = "best fixed\ncomposite",
                              "random forest (OOF)" = "random forest\n(OOF)")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, 0.2)) +
  labs(x = NULL, y = "AUROC", title = "d  Mechanism 2 + interaction control") +
  theme_nature() +
  theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"),
        axis.text.x = element_text(size = 8.5, lineheight = 0.85))

# ---- Combine ----
fig <- (pa + pb + pc) / pd + plot_layout(heights = c(1, 0.9))
fig <- fig + plot_annotation(
  title = "Fig. 2 Controlled provenance generates apparent integration gains (two mechanisms + positive control)",
  theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, family = "sans"))
)

save_fig(fig, "Fig2_controlled_provenance", w = 11, h = 7)
