# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Extended Data Fig. 3 | Representation sensitivity
# Panel a: pairwise-inversion information index I across composites (fig3b).
# Panel b: coverage-stratified inversion (fig3b / fig3c). Real data, editable PDF + 300 dpi PNG.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))
library(ggplot2)
library(grid)
library(patchwork)
outd <- file.path(PKG, "figures")

d <- read.csv(file.path(dat_dir, "fig3b_pairwise_inversion.csv"), stringsAsFactors = FALSE)
d$composite <- factor(d$composite, levels = d$composite[order(d$I_pairwise_inversion)])

pa <- ggplot(d, aes(x = composite, y = I_pairwise_inversion)) +
  geom_col(width = 0.6, fill = "#4a78b5") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.18, colour = "#222") +
  geom_text(aes(label = sprintf("%.3f", I_pairwise_inversion)), vjust = -0.45, size = 3.0, family = FONT) +
  scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0)) +
  labs(title = "a | Pairwise-inversion information index I", x = NULL, y = "I (pairwise inversion)") +
  theme_nature(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 8.5),
        plot.title = element_text(size = 11, face = "bold", hjust = 0))

# coverage-stratified: long form
cov <- data.frame(
  composite = rep(d$composite, 3),
  stratum = rep(c("both high", "mixed", "both low"), each = nrow(d)),
  value = c(d$I_both_high_coverage, d$I_mixed_coverage, d$I_both_low_coverage)
)
pb <- ggplot(cov, aes(x = composite, y = value, fill = stratum)) +
  geom_col(position = "dodge", width = 0.72) +
  scale_fill_manual(values = c("#9fc5e8", "#f9cb9c", "#e06666")) +
  scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0)) +
  labs(title = "b | Coverage-stratified inversion", x = NULL, y = "I (stratified)", fill = "Coverage") +
  theme_nature(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 8.5),
        legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(t = 3, unit = "pt"),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 11, face = "bold", hjust = 0))

fig <- pa + pb + plot_layout(widths = c(1, 1.25))
save_fig(fig, "ED_Fig3_representation_sensitivity", w = 11.5, h = 5.0, out_dir = outd)
