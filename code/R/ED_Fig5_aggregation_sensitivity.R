# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Extended Data Fig. 5 | Aggregation sensitivity
# Panel a: rank correlation (Kendall tau, Spearman rho) across composites (fig3c).
# Panel b: pairwise-inversion information index I across composites (fig3c). Real data, editable PDF + 300 dpi PNG.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))
library(ggplot2)
library(grid)
library(dplyr)
library(tidyr)
library(patchwork)
outd <- file.path(PKG, "figures")

d <- read.csv(file.path(dat_dir, "fig3c_aggregation_robustness.csv"), stringsAsFactors = FALSE)
d$composite <- factor(d$composite, levels = d$composite[order(d$kendall_tau)])

long <- d %>% select(composite, kendall_tau, spearman_rho) %>% pivot_longer(-composite, names_to = "metric", values_to = "value")
long$metric <- factor(long$metric, levels = c("kendall_tau", "spearman_rho"),
                      labels = c("Kendall tau", "Spearman rho"))

pa <- ggplot(long, aes(x = composite, y = value, fill = metric)) +
  geom_col(position = "dodge", width = 0.72) +
  scale_fill_manual(values = c("#4a78b5", "#e69138")) +
  scale_y_continuous(limits = c(-0.2, 0.3), expand = c(0, 0)) +
  geom_hline(yintercept = 0, colour = "#888", linewidth = 0.4) +
  labs(title = "a | Composite rank correlation", x = NULL, y = "Correlation", fill = "Metric") +
  theme_nature(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 8.5),
        legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(t = 3, unit = "pt"),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 11, face = "bold", hjust = 0))

pb <- ggplot(d, aes(x = composite, y = I_pairwise_inversion)) +
  geom_point(size = 3.2, colour = "#6aa84f") +
  geom_segment(aes(x = as.numeric(composite), xend = as.numeric(composite),
                   y = inv_both_high_cov, yend = inv_both_low_cov),
               colour = "#6aa84f", linewidth = 1.2) +
  scale_y_continuous(limits = c(0.3, 0.7), expand = c(0, 0)) +
  labs(title = "b | Pairwise-inversion I (coverage range)", x = NULL, y = "I (pairwise inversion)") +
  theme_nature(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 8.5),
        plot.title = element_text(size = 11, face = "bold", hjust = 0))

fig <- pa + pb + plot_layout(widths = c(1.2, 1))
save_fig(fig, "ED_Fig5_aggregation_sensitivity", w = 11.5, h = 5.0, out_dir = outd)
