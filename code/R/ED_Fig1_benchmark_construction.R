# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Extended Data Fig. 1 | Benchmark construction and candidate flow
# Conceptual schematic (no scores/ranks). Editable PDF + 300 dpi PNG via ragg/Quartz.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))
library(ggplot2)
library(grid)
outd <- file.path(PKG, "figures")

stages <- data.frame(
  x = c(0.5, 0.5, 1.5, 1.5, 2.5, 3.5),
  y = c(3.2, 1.2, 3.2, 1.2, 2.2, 2.2),
  lab = c("PDAC candidate\nspace", "CRC candidate\nspace", "PDAC eligible\nbenchmark",
          "CRC eligible\nbenchmark", "Merged\nbenchmark\nuniverse",
          "Final outcome\n(PDAC dependency /\nCRC alteration)"),
  fill = c("#dbe9f4", "#fde9d9", "#bcdfe8", "#fbd3b0", "#d9ead3", "#cfe2f3")
)
excl <- data.frame(x = 2.5, y = 3.95, lab = "Excluded: no outcome /\ninsufficient evidence (both cancers)")

p <- ggplot() +
  geom_rect(data = stages, aes(xmin = x - 0.44, xmax = x + 0.44,
                              ymin = y - 0.56, ymax = y + 0.56),
            fill = stages$fill, colour = "white", linewidth = 0.6) +
  geom_text(data = stages, aes(x = x, y = y, label = lab),
            size = 3.0, family = FONT, colour = "#222222", hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  geom_text(data = excl, aes(x = x, y = y, label = lab),
            size = 2.6, family = FONT, colour = "#990000", hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  annotate("segment", x = 0.94, y = 3.2, xend = 1.06, yend = 3.2,
           arrow = arrow(length = unit(0.08, "in")), colour = "#555555") +
  annotate("segment", x = 0.94, y = 1.2, xend = 1.06, yend = 1.2,
           arrow = arrow(length = unit(0.08, "in")), colour = "#555555") +
  annotate("segment", x = 1.94, y = 3.2, xend = 2.02, yend = 2.20,
           arrow = arrow(length = unit(0.08, "in")), colour = "#555555") +
  annotate("segment", x = 1.94, y = 1.2, xend = 2.02, yend = 2.20,
           arrow = arrow(length = unit(0.08, "in")), colour = "#555555") +
  annotate("segment", x = 2.94, y = 2.2, xend = 3.06, yend = 2.2,
           arrow = arrow(length = unit(0.08, "in")), colour = "#555555") +
  coord_fixed(ratio = 0.72) +
  labs(title = "Extended Data Fig. 1 | Benchmark construction and candidate flow") +
  theme_nature(base_size = 11) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

save_fig(p, "ED_Fig1_benchmark_construction", w = 8.5, h = 5.2, out_dir = outd)
