# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Extended Data Fig. 6 | Evaluation sensitivity
# Distribution of Delta AUROC across evaluation conditions, by composite (fig6a).
# Real data, editable PDF + 300 dpi PNG. No overlap: dodged/violin with adequate width.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))
library(ggplot2)
library(grid)
library(dplyr)
library(patchwork)
outd <- file.path(PKG, "figures")

d <- read.csv(file.path(dat_dir, "fig6a_delta_distributions.csv"), stringsAsFactors = FALSE)
d$method <- factor(d$method, levels = sort(unique(d$method)))

# panel a: jittered points per method
pa <- ggplot(d, aes(x = method, y = delta_auroc, colour = method)) +
  geom_hline(yintercept = 0, colour = "#888", linewidth = 0.4) +
  geom_jitter(width = 0.18, size = 1.6, alpha = 0.8) +
  scale_colour_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(-0.5, 0.25), expand = c(0, 0)) +
  labs(title = "a | Delta AUROC by composite (per endpoint)", x = NULL, y = "Delta AUROC") +
  theme_nature(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 8.0),
        plot.title = element_text(size = 11, face = "bold", hjust = 0))

# panel b: mean +/- sd across endpoints per method
summ <- d %>% group_by(method) %>%
  summarise(m = mean(delta_auroc), s = sd(delta_auroc))
summ$method <- factor(summ$method, levels = levels(d$method))
pb <- ggplot(summ, aes(x = method, y = m, fill = method)) +
  geom_col(width = 0.62) +
  geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 0.18, colour = "#222") +
  geom_text(aes(label = sprintf("%.3f", m)), vjust = ifelse(summ$m >= 0, -0.4, 1.4), size = 3.0, family = FONT) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(-0.45, 0.2), expand = c(0, 0)) +
  labs(title = "b | Mean Delta AUROC (s.d.) across endpoints", x = NULL, y = "Mean Delta AUROC") +
  theme_nature(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 8.0),
        plot.title = element_text(size = 11, face = "bold", hjust = 0))

fig <- pa + pb + plot_layout(widths = c(1.3, 1))
save_fig(fig, "ED_Fig6_evaluation_sensitivity", w = 12.0, h = 5.2, out_dir = outd)
