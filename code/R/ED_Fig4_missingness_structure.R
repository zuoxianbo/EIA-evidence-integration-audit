# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Extended Data Fig. 4 | Missingness structure
# Panel a: missingness decomposition of AUROC (values+missingness vs values-only vs missingness-only) (fig3d).
# Panel b: sentinel-layer missingness coverage (fig3a). Real data, editable PDF + 300 dpi PNG.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))
library(ggplot2)
library(grid)
library(dplyr)
library(patchwork)
outd <- file.path(PKG, "figures")

dd <- read.csv(file.path(dat_dir, "fig3d_missingness_decomposition.csv"), stringsAsFactors = FALSE)
# aggregate mean auroc per configuration across endpoints
agg <- dd %>% group_by(scorer) %>% summarise(auroc = mean(auroc, na.rm = TRUE))
agg$scorer <- factor(agg$scorer, levels = agg$scorer[order(agg$auroc)])

pa <- ggplot(agg, aes(x = scorer, y = auroc, fill = scorer)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.4, size = 3.0, family = FONT) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(
    "per-layer coverage indicator" = "per-layer coverage\nindicator",
    "coverage count" = "coverage count",
    "values-only" = "values-only",
    "missingness-only" = "missingness-only",
    "values+missingness (frozen harmonic)" = "values+miss.\n(frozen harmonic)")) +
  scale_y_continuous(limits = c(0, 0.72), expand = c(0, 0)) +
  labs(title = "a | Missingness decomposition of AUROC", x = NULL, y = "Mean AUROC") +
  theme_nature(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 18, hjust = 1, vjust = 1, size = 8.5),
        plot.title = element_text(size = 11, face = "bold", hjust = 0))

mm <- read.csv(file.path(dat_dir, "fig3a_sentinel_missingness.csv"), stringsAsFactors = FALSE)
mm$layer_label <- gsub("\n", " ", mm$layer_label)
mm$layer_label <- factor(mm$layer_label, levels = mm$layer_label[order(mm$pct_genes_at_sentinel)])

pb <- ggplot(mm, aes(x = layer_label, y = pct_genes_at_sentinel)) +
  geom_col(width = 0.62, fill = "#b45f06") +
  geom_text(aes(label = sprintf("%.1f%%", pct_genes_at_sentinel)), vjust = -0.35, size = 2.9, family = FONT) +
  scale_y_continuous(limits = c(0, max(mm$pct_genes_at_sentinel) * 1.18), expand = c(0, 0)) +
  labs(title = "b | Genes at sentinel, by evidence layer", x = NULL, y = "% genes at sentinel") +
  theme_nature(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8.0),
        plot.title = element_text(size = 11, face = "bold", hjust = 0))

fig <- pa + pb + plot_layout(widths = c(1.1, 1.4))
save_fig(fig, "ED_Fig4_missingness_structure", w = 12.0, h = 5.0, out_dir = outd)
