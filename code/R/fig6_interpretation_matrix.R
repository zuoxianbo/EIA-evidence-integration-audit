# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 6 (clean.docx numbering): EIA interpretation matrix
# Conceptual synthesis of provenance independence x external persistence.
# No data, no scores/ranks -- a qualitative interpretation matrix built directly
# from the Discussion's four-state hierarchy. Editable PDF + 300 dpi PNG via ragg/Quartz.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))
library(ggplot2)
library(grid)

outd <- file.path(PKG, "figures")

# 2x2 grid: x = provenance independence, y = external persistence
cells <- data.frame(
  x = c(1, 2, 1, 2),
  y = c(1, 1, 2, 2),
  prov = c("Not independent", "Independent", "Not independent", "Independent"),
  pers = c("Not persistent", "Not persistent", "Persistent", "Persistent"),
  fill = c("#f4cccc", "#fce5cd", "#fff2cc", "#d9ead3"),
  border = c("#c0392b", "#e69138", "#b8860b", "#2a9d8f")
)

# short state label (top of each cell, wrapped to 2 lines so it stays inside the box)
cells$state <- c(
  "1  Performance gain\nobserved",
  "2  Gain under independent\nprovenance",
  "3  Gain persistent externally\n(provenance unclear)",
  "4  Independent provenance\n+ external persistence"
)
cells$interp <- c(
  "Delta > 0 within the original ecology.\nWeakest claim; does not establish\ninformation attribution.",
  "Provenance control fails to remove the\ngain without independent evidence.\nStronger, but still ecology-bound.",
  "Gain survives an independent evaluation\nenvironment; source not yet separated\nfrom representation.",
  "Strongest: gain persists externally AND\nis not explained by dependent\nprovenance. Supports attribution."
)

p <- ggplot() +
  # cell fills
  geom_rect(data = cells, aes(xmin = x - 0.48, xmax = x + 0.48,
                              ymin = y - 0.48, ymax = y + 0.48, fill = fill),
            colour = "white", linewidth = 0.5) +
  # cell borders (claim strength)
  geom_rect(data = cells, aes(xmin = x - 0.48, xmax = x + 0.48,
                              ymin = y - 0.48, ymax = y + 0.48),
            fill = NA, colour = cells$border, linewidth = 1.4) +
  # state label near top of each cell (2 lines, kept inside the box)
  geom_text(data = cells, aes(x = x, y = y + 0.34, label = state),
            size = 3.1, fontface = "bold", family = FONT, colour = "#1a1a1a",
            hjust = 0.5, vjust = 1, lineheight = 0.95) +
  # interpretation inside cell
  geom_text(data = cells, aes(x = x, y = y - 0.05, label = interp),
            size = 2.7, family = FONT, colour = "#333333",
            hjust = 0.5, vjust = 0.5, lineheight = 0.95) +
  # axis labels (outside)
  scale_x_continuous(limits = c(0.4, 2.6), expand = c(0, 0),
                    breaks = c(1, 2), labels = c("Not independent", "Independent")) +
  scale_y_continuous(limits = c(0.4, 2.6), expand = c(0, 0),
                    breaks = c(1, 2), labels = c("Not persistent", "Persistent")) +
  scale_fill_identity() +
  labs(x = "Provenance independence", y = "External persistence",
       title = "Figure 6 | EIA interpretation matrix") +
  theme_nature(base_size = 11) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_text(size = 9.5, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA))

save_fig(p, "Fig6_interpretation_matrix", w = 7.2, h = 6.2, out_dir = outd)
