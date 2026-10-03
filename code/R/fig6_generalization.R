#!/usr/bin/env Rscript
# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 6 (v26 numbering):: independent Delta landscape on external evaluations
# Per critical review 6.7/6.8/6.10: report the FULL Delta distribution (not only the max),
# audit-condition semantics instead of PASS=success, no "prospective" wording.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
source(file.path(base_dir, "common.R"))

library(ggplot2)
library(patchwork)

# ---- Panel a: full Delta distributions on E4/E5/E6 (v33: native 7 composites) ----
fa <- read.csv(file.path(base_dir, "source_data", "fig6a_delta_distributions.csv"), stringsAsFactors = FALSE)
fa$endpoint <- factor(fa$endpoint, levels = c("E4 CRC", "E5 clinical", "E6 drug response"))
fa$method <- sub(" \\(multiplicative\\)", "", fa$method)
fa$method <- sub(" aggregation", "", fa$method)
fa$method <- sub(" mean", "", fa$method)
fa$is_max <- FALSE
for (ep in levels(fa$endpoint)) {
  idx <- fa$endpoint == ep
  fa$is_max[idx] <- fa$delta_auroc[idx] == max(fa$delta_auroc[idx])
}

pa <- ggplot(fa, aes(x = endpoint, y = delta_auroc)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.6) +
  geom_point(aes(colour = method, shape = is_max), size = 3.2, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c("Harmonic" = "#2e7bd6", "Arithmetic" = "#8e9aaf", "ECS" = "#c0392b",
                                "Additive" = "#e9c46a", "Geometric" = "#2a9d8f", "Rank" = "#6d597a",
                                "Weighted rank" = "#9c6644"),
                      name = NULL) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("composite", "best composite"), name = NULL) +
  scale_y_continuous(limits = c(-0.32, 0.03), breaks = seq(-0.3, 0, 0.1)) +
  labs(x = NULL, y = expression(Delta * "AUROC vs best single layer"),
       title = "a  External evaluations: all seven composites, full Delta distribution\n(audit condition satisfied: no composite shows positive Delta)") +
  theme_nature() +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = 4, unit = "pt"),
        legend.spacing.y = unit(3, "pt"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(7, "pt")) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 1))

# ---- Panel b: E6 discrimination vs ranking (centrality vs harmonic) ----
fb <- read.csv(file.path(base_dir, "source_data", "fig5b_e6_discrimination_ranking.csv"), stringsAsFactors = FALSE)
fb <- fb[fb$method %in% c("STRING centrality", "Harmonic mean"), ]
fb$method <- factor(fb$method, levels = c("STRING centrality", "Harmonic mean"))
fb$metric <- factor(fb$metric, levels = c("AUROC", "AUPRC"))

cols_b <- c("STRING centrality" = "#8e9aaf", "Harmonic mean" = "#2e7bd6")
pb <- ggplot(fb, aes(x = metric, y = value, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = cols_b, name = NULL) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "score",
       title = "b  E6 external drug-response evaluation:\ndiscrimination != ranking (AUPRC inverts vs AUROC)") +
  theme_nature() +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        legend.position = c(0.62, 0.85),
        legend.background = element_blank())

# ---- Panel c: trivial-transfer baselines vs supervised transfer (S1 self-audit) ----
fc <- read.csv(file.path(base_dir, "source_data", "fig6c_transfer_baselines.csv"), stringsAsFactors = FALSE)
fc <- data.frame(method = rep(fc$method, 2),
                 direction = rep(c("PDAC -> CRC", "CRC -> PDAC"), each = nrow(fc)),
                 auroc = c(fc$pdac_to_crc, fc$crc_to_pdac))
fc$method_orig <- fc$method
fc$is_rf <- fc$method_orig == "Random forest (5-fold ensemble, no refit)"
fc$method <- factor(fc$method_orig,
  levels = c("Random forest (5-fold ensemble, no refit)",
             "Label copy (binary essential flag)",
             "Dependency score (continuous, direct)"),
  labels = c("Random forest\n(no refit)",
             "Label copy\n(binary)",
             "Dependency\nscore"))
pc <- ggplot(fc, aes(x = method, y = auroc, fill = is_rf)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           width = 0.7, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#8e9aaf"), guide = "none") +
  facet_wrap(~ direction, nrow = 1) +
  geom_hline(data = data.frame(y = 0.929), aes(yintercept = y), linetype = "dashed",
             colour = "#404040", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", auroc)), position = position_dodge(width = 0.9),
            vjust = -0.4, size = 2.6, family = "sans") +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "AUROC",
       title = "c  Transfer claim below trivial baselines\n(RF 0.818 vs label-copy 0.929 / dependency 0.982)") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 9.0, face = "bold", lineheight = 1.05),
        axis.text.x = element_text(size = 7.0, lineheight = 0.92),
        strip.text = element_text(size = 8.5, face = "bold"))

# ---- Combine ----
fig <- pa + pb + pc + plot_layout(widths = c(1.35, 1.0, 1.2))
fig <- fig + plot_annotation(
  title = "Fig. 6 Independent external evaluations determine whether the gain generalizes",
  theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, family = "sans"))
)

save_fig(fig, "Fig6_generalization", w = 15.0, h = 5.6)
