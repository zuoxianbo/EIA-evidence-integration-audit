# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Extended Data Fig. 2 | Evidence provenance graph
# Conceptual: each evidence layer traces to primary source -> integrated DB -> derived feature.
# No scores/ranks. Editable PDF + 300 dpi PNG via ragg/Quartz.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))
library(ggplot2)
library(grid)
outd <- file.path(PKG, "figures")

# three columns: primary sources (x=1), integrated databases (x=2.5), derived features (x=4), layers (x=5.5)
src  <- data.frame(x = 1.0, y = c(4.6, 3.4, 2.2, 1.0), lab = c("STRING", "Open Targets", "COSMIC / gnomAD", "DepMap / HPA"))
db   <- data.frame(x = 2.5, y = c(3.8, 2.6, 1.4), lab = c("Protein network", "Target-disease", "Dependency / expression"))
feat <- data.frame(x = 4.0, y = c(3.8, 2.6, 1.4), lab = c("Centrality", "Assoc. score", "Rank / level"))
lay  <- data.frame(x = 5.5, y = c(3.8, 2.6, 1.4), lab = c("Layer A", "Layer B", "Layer C"))

p <- ggplot() +
  geom_rect(data = src,  aes(xmin = x - 0.55, xmax = x + 0.55, ymin = y - 0.42, ymax = y + 0.42), fill = "#e8f0fe", colour = "white", linewidth = 0.6) +
  geom_rect(data = db,   aes(xmin = x - 0.62, xmax = x + 0.62, ymin = y - 0.42, ymax = y + 0.42), fill = "#d9ead3", colour = "white", linewidth = 0.6) +
  geom_rect(data = feat, aes(xmin = x - 0.55, xmax = x + 0.55, ymin = y - 0.42, ymax = y + 0.42), fill = "#fff2cc", colour = "white", linewidth = 0.6) +
  geom_rect(data = lay,  aes(xmin = x - 0.55, xmax = x + 0.55, ymin = y - 0.42, ymax = y + 0.42), fill = "#fce5cd", colour = "white", linewidth = 0.6) +
  geom_text(data = src,  aes(x = x, y = y, label = lab), size = 2.8, family = FONT, colour = "#222", hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  geom_text(data = db,   aes(x = x, y = y, label = lab), size = 2.8, family = FONT, colour = "#222", hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  geom_text(data = feat, aes(x = x, y = y, label = lab), size = 2.8, family = FONT, colour = "#222", hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  geom_text(data = lay,  aes(x = x, y = y, label = lab), size = 2.8, family = FONT, colour = "#222", hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  # edges
  annotate("segment", x = 1.55, y = 4.6, xend = 1.88, yend = 3.8, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 1.55, y = 3.4, xend = 1.88, yend = 3.8, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 1.55, y = 3.4, xend = 1.88, yend = 2.6, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 1.55, y = 2.2, xend = 1.88, yend = 2.6, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 1.55, y = 2.2, xend = 1.88, yend = 1.4, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 1.55, y = 1.0, xend = 1.88, yend = 1.4, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 3.12, y = 3.8, xend = 3.45, yend = 3.8, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 3.12, y = 2.6, xend = 3.45, yend = 2.6, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 3.12, y = 1.4, xend = 3.45, yend = 1.4, colour = "#888", linewidth = 0.5) +
  annotate("segment", x = 4.55, y = 3.8, xend = 4.95, yend = 3.8, colour = "#888", linewidth = 0.5, arrow = arrow(length = unit(0.06, "in"))) +
  annotate("segment", x = 4.55, y = 2.6, xend = 4.95, yend = 2.6, colour = "#888", linewidth = 0.5, arrow = arrow(length = unit(0.06, "in"))) +
  annotate("segment", x = 4.55, y = 1.4, xend = 4.95, yend = 1.4, colour = "#888", linewidth = 0.5, arrow = arrow(length = unit(0.06, "in"))) +
  coord_fixed(ratio = 0.85) +
  labs(title = "Extended Data Fig. 2 | Evidence provenance graph",
       x = NULL, y = NULL) +
  theme_nature(base_size = 11) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

save_fig(p, "ED_Fig2_provenance_graph", w = 8.5, h = 5.2, out_dir = outd)
