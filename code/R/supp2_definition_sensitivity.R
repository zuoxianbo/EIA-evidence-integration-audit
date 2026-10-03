# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 2 panel d addition is produced by fig2_provenance.R itself in v26.
# This script makes the standalone mechanism-2 data check + Supp Fig 2 (definition sensitivity)
base_dir <- file.path(ROOT, "code", "R")
source(file.path(base_dir, "common.R"))
library(patchwork)

sd <- file.path(base_dir, "source_data")
outd <- file.path(PKG, "figures")

# ---------------- Supp Fig 2: definition sensitivity ----------------
d <- read.csv(file.path(sd, "supp2_definition_sensitivity.csv"), stringsAsFactors = FALSE)
d$definition <- factor(d$definition,
  levels = c("A_ratio", "B_zeffect", "C_lineage_adj", "D_mixed"),
  labels = c("A  ratio", "B  z-effect", "C  lineage-adj", "D  mixed"))
d$lab <- sprintf("Delta = %+.3f", d$delta)

p1 <- ggplot(d, aes(x = definition, y = harmonic_auroc, group = 1)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50") +
  geom_line(colour = "grey55", linewidth = 0.4) +
  geom_point(aes(fill = "Harmonic mean"), shape = 21, size = 2.4, stroke = 0.6) +
  geom_point(aes(y = rf_auroc, colour = "Random forest (OOF)"), shape = 16, size = 2.0) +
  geom_point(aes(y = best_single_auroc, shape = "Best single layer"), size = 2.2,
             colour = "#1a1a1a", stroke = 0.7) +
  scale_fill_manual(values = c("Harmonic mean" = "#2166ac"), guide = "none") +
  scale_colour_manual(values = c("Random forest (OOF)" = "#e07b39"),
                      name = NULL) +
  scale_shape_manual(values = c("Best single layer" = 17), name = NULL) +
  geom_text(aes(label = lab), vjust = -1.15, size = 2.5, family = "sans", colour = "#404040") +
  labs(x = "Selective-dependency endpoint definition (top 250 genes each)",
       y = "AUROC",
       title = "Supp Fig 2a  Estimand sensitivity to endpoint definition (PDAC selective dependency)") +
  theme_nature(base_size = 10) +
  theme(legend.position = c(0.72, 0.18),
        legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3))

# cross-definition concordance bar (from v26 JSON, hard numbers passed via CSV below)
# Jaccard values are written by extract step; kept in supp2b csv if present
f2 <- file.path(sd, "supp2b_definition_concordance.csv")
if (file.exists(f2)) {
  j <- read.csv(f2, stringsAsFactors = FALSE)
  p2 <- ggplot(j, aes(x = pair, y = jaccard)) +
    geom_col(fill = "#7e9db9", colour = "black", linewidth = 0.3, width = 0.6) +
    labs(x = NULL, y = "Jaccard (top-250 sets)",
         title = "Supp Fig 2b  Cross-definition concordance of the gold standard") +
    theme_nature(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  fig <- p1 / p2 + plot_layout(heights = c(1.4, 1))
} else {
  fig <- p1
}

# ---------------- Supp Fig 2c: strict-label sensitivity (S4) ----------------
sc <- read.csv(file.path(sd, "supp2c_strict_labels.csv"), stringsAsFactors = FALSE)
sc_long <- do.call(rbind, lapply(split(sc, sc$endpoint), function(d) {
  data.frame(endpoint = d$endpoint,
             metric = c("Best single layer", "Best composite", "Random forest (OOF)"),
             auroc = c(d$best_single_auroc, d$best_composite_auroc, d$rf_auroc))
}))
sc_long$endpoint <- factor(sc_long$endpoint, levels = c("E1-strict", "E4-strict"))
sc_long$metric <- factor(sc_long$metric,
  levels = c("Best single layer", "Best composite", "Random forest (OOF)"),
  labels = c("Best single\nlayer", "Best\ncomposite", "Random forest\n(OOF)"))
p3 <- ggplot(sc_long, aes(x = metric, y = auroc, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.65,
           colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("Best single\nlayer" = "#8e9aaf",
                               "Best\ncomposite" = "#2e7bd6",
                               "Random forest\n(OOF)" = "#c0392b"), guide = "none") +
  facet_wrap(~ endpoint, nrow = 1) +
  geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.4, size = 2.6, family = "sans") +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "AUROC",
       title = "Supp Fig 2c  Strict-label sensitivity (dependency_score > 0.5)\nRF still leads best single by +0.045 / +0.047") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 9.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8, lineheight = 0.85),
        strip.text = element_text(size = 8.5, face = "bold"))

# ---------------- Supp Fig 2d: annotation-poor subset (S5) ----------------
sd2 <- read.csv(file.path(sd, "supp2d_annotated_subset.csv"), stringsAsFactors = FALSE)
sd2$ci_lo <- suppressWarnings(as.numeric(sd2$delta_comp_ci_lo))
sd2$ci_hi <- suppressWarnings(as.numeric(sd2$delta_comp_ci_hi))
sd2$delta_comp <- suppressWarnings(as.numeric(sd2$delta_comp))
sd2$lab <- factor(sub(" PDAC-selective.*", "", sd2$endpoint),
                  levels = rev(c("E1", "SelA", "SelB", "SelD")))
p4 <- ggplot(sd2, aes(x = delta_comp, y = lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                 linewidth = 0.5, colour = "grey50", na.rm = TRUE) +
  geom_point(size = 2.2, colour = "#c0392b", stroke = 0.6) +
  scale_y_discrete(limits = rev(levels(sd2$lab))) +
  labs(x = expression(Delta * " (best composite - best single, AUROC)"),
       y = NULL,
       title = "Supp Fig 2d  Annotation-poor subset (coverage >= 4/9 layers)\nDelta flips negative: gain was a coverage artifact") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 9.5, face = "bold", lineheight = 1.1))

fig <- (p1 | p3) / (p2 | p4) + plot_layout(heights = c(1.4, 1.0))
save_fig(fig, "Supp2_definition_sensitivity", w = 13.5, h = 9.0, out_dir = outd)
