#!/usr/bin/env Rscript
# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 4: Aggregation strategy changes the observed performance gain
# Rebuilt (v39) to six panels a-f per the locked manuscript legend, using only
# real data from scorer_manifest.csv (fixed-form vs supervised integrators,
# internal-benchmark performance, and Delta under alternative aggregation).
# (The old "fair-baseline / Delta landscape" panels are superseded by Fig. 2/5.)

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))

library(ggplot2)
library(patchwork)

manifest <- read.csv(file.path(dat_dir, "scorer_manifest.csv"), stringsAsFactors = FALSE)

internal_ep <- c("E1 PDAC pan-dependency", "E2 PDAC-enriched dependency",
                 "E3 conjunctive actionability", "E3-A leakage-controlled essentiality",
                 "E3-C out-of-evidence druggability")
ep_short <- c(
  "E1 PDAC pan-dependency"              = "E1",
  "E2 PDAC-enriched dependency"         = "E2",
  "E3 conjunctive actionability"        = "E3",
  "E3-A leakage-controlled essentiality" = "E3-A",
  "E3-C out-of-evidence druggability"   = "E3-C")
ff_short <- c(
  "Additive"            = "Add",
  "Arithmetic mean"     = "Arith",
  "ECS (multiplicative)" = "ECS",
  "Geometric mean"      = "Geo",
  "Harmonic mean"       = "Harm",
  "Rank aggregation"    = "Rank agg",
  "Weighted rank"       = "Wtd rank")
sup_short <- c(
  "Elastic net"        = "Elastic net",
  "Logistic regression" = "Logistic reg",
  "Random forest"      = "Random forest")

ff <- manifest[manifest$class == "fixed-form integration" & manifest$endpoint %in% internal_ep, ]
sup <- manifest[manifest$class == "supervised out-of-fold" & manifest$endpoint %in% internal_ep, ]
sl <- manifest[manifest$class == "single-layer evidence" & manifest$endpoint %in% internal_ep, ]

# ---- Panel a: fixed-form composite (mean AUROC over internal benchmark) ----
ff_mean <- aggregate(auroc ~ scorer, data = ff, mean)
ff_mean$scorer <- factor(ff_mean$scorer, levels = names(ff_short))
pa <- ggplot(ff_mean, aes(x = scorer, y = auroc, fill = scorer)) +
  geom_bar(stat = "identity", width = 0.7, colour = "black", linewidth = 0.2, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.4, size = 2.8, family = "sans") +
  scale_x_discrete(labels = ff_short) +
  scale_y_continuous(limits = c(0, 0.80), breaks = seq(0, 0.8, 0.2)) +
  labs(x = NULL, y = "Mean AUROC (internal benchmark)",
       title = "a  Fixed-form evidence composite (7 rules)") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8))

# ---- Panel b: supervised evidence integrator (mean over all evaluated endpoints) ----
sup_all <- manifest[manifest$class == "supervised out-of-fold", ]
sup_mean <- aggregate(auroc ~ scorer, data = sup_all, mean)
sup_mean$scorer <- factor(sup_mean$scorer, levels = names(sup_short))
pb <- ggplot(sup_mean, aes(x = scorer, y = auroc, fill = scorer)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.2, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.4, size = 2.8, family = "sans") +
  scale_x_discrete(labels = sup_short) +
  scale_y_continuous(limits = c(0, 0.80), breaks = seq(0, 0.8, 0.2)) +
  labs(x = NULL, y = "Mean AUROC (all evaluated endpoints)",
       title = "b  Supervised evidence integrator (out-of-fold)") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8.5, lineheight = 0.85))

# ---- Panel c: internal benchmark performance, fixed-form (heatmap) ----
ff$scorer <- factor(ff$scorer, levels = names(ff_short))
ff$ep <- factor(ff$endpoint, levels = internal_ep)
ff$lab <- ff_short[as.character(ff$scorer)]
pc <- ggplot(ff, aes(x = ep, y = lab, fill = auroc)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", auroc)), size = 2.6, family = "sans",
            colour = ifelse(ff$auroc > 0.62, "white", "black")) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2166ac", limits = c(0.40, 0.95),
                      oob = scales::squish, name = "AUROC") +
  scale_x_discrete(labels = ep_short) +
  labs(x = NULL, y = NULL, title = "c  Internal benchmark: fixed-form composites") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8.5), axis.text.y = element_text(size = 8.5),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 7.5))

# ---- Panel d: internal benchmark performance, supervised (bar; internal = E1) ----
# Supervised OOF integrators are defined on E1 only within the internal benchmark
# (their remaining evaluations are on external endpoints), so this panel is a bar.
sup$scorer <- factor(sup$scorer, levels = names(sup_short))
pd <- ggplot(sup, aes(x = scorer, y = auroc, fill = scorer)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.2, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.4, size = 2.8, family = "sans") +
  scale_x_discrete(labels = sup_short) +
  scale_y_continuous(limits = c(0, 0.80), breaks = seq(0, 0.8, 0.2)) +
  labs(x = NULL, y = "AUROC (E1 internal)",
       title = "d  Internal benchmark (E1): supervised integrators") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8.5, lineheight = 0.85))

# ---- best single-layer baseline per endpoint for Delta ----
best_sl <- aggregate(auroc ~ endpoint, data = sl, max)
colnames(best_sl)[2] <- "best_sl"

# ---- Panel e: Delta under alternative fixed-form aggregation ----
ff_d <- merge(ff, best_sl, by = "endpoint")
ff_d$delta <- ff_d$auroc - ff_d$best_sl
ff_d_agg <- aggregate(delta ~ scorer, data = ff_d, mean)
ff_d_agg$scorer <- factor(ff_d_agg$scorer, levels = names(ff_short))
pe <- ggplot(ff_d_agg, aes(x = scorer, y = delta, fill = delta > 0)) +
  geom_bar(stat = "identity", width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.3f", delta)), vjust = ifelse(ff_d_agg$delta > 0, -0.4, 1.4),
            size = 2.8, family = "sans") +
  scale_x_discrete(labels = ff_short) +
  scale_fill_manual(values = c("TRUE" = "#2a9d8f", "FALSE" = "#c0392b"), guide = "none") +
  scale_y_continuous(limits = c(-0.28, 0.10), breaks = seq(-0.25, 0.10, 0.05)) +
  labs(x = NULL, y = expression(Delta * " AUROC vs best single layer"),
       title = "e  Delta under alternative fixed-form aggregation") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8))

# ---- Panel f: Delta under alternative supervised aggregation ----
sup_d <- merge(sup, best_sl, by = "endpoint")
sup_d$delta <- sup_d$auroc - sup_d$best_sl
sup_d_agg <- aggregate(delta ~ scorer, data = sup_d, mean)
sup_d_agg$scorer <- factor(sup_d_agg$scorer, levels = names(sup_short))
pf <- ggplot(sup_d_agg, aes(x = scorer, y = delta, fill = delta > 0)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.3f", delta)), vjust = ifelse(sup_d_agg$delta > 0, -0.4, 1.4),
            size = 2.8, family = "sans") +
  scale_x_discrete(labels = sup_short) +
  scale_fill_manual(values = c("TRUE" = "#2a9d8f", "FALSE" = "#c0392b"), guide = "none") +
  scale_y_continuous(limits = c(-0.10, 0.10), breaks = seq(-0.10, 0.10, 0.05)) +
  labs(x = NULL, y = expression(Delta * " AUROC vs best single layer"),
       title = "f  Delta under alternative supervised aggregation") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8.5, lineheight = 0.85))

# ---- combine ----
fig <- (pa | pb) / (pc | pd) / (pe | pf)
fig <- fig + plot_annotation(
  title = "Fig. 4 Aggregation strategy changes the observed performance gain",
  theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, family = "sans"))
)

save_fig(fig, "Fig4_aggregation_strategy", w = 12.5, h = 11)
