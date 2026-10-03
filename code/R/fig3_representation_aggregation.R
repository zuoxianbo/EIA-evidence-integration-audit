#!/usr/bin/env Rscript
# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 3: Evidence availability and representation affect the measured increment
# Rebuilt (v39) to six panels a-f per the locked manuscript legend, using only
# real source CSVs. Data sources are documented in FIGURE_FIX_REPORT.md.

rm(list = ls())
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))

library(ggplot2)
library(patchwork)

# ---- shared short labels -------------------------------------------------
layer_short <- c(
  "string_centrality"        = "STRING\ncentr.",
  "mutation_freq"            = "Mut.\nfreq.",
  "impc_animal_ko"           = "IMPC\nKO",
  "genetic_constraint"       = "Gen.\nconstraint",
  "cancer_driver"            = "Cancer\ndriver",
  "ot_genetics_pdac"         = "OT\n(PDAC)",
  "druggability"             = "Drug",
  "hpa_pdac_prognostic"      = "HPA\nprog",
  "hpa_rna_tissue_spec"      = "HPA\nRNA")
ep_levels <- c("E1 PDAC pan-dependency", "E4 CRC pan-dependency",
               "E5 clinical concordance", "SelA PDAC-selective (A_ratio)",
               "SelB PDAC-selective (B_zeffect==C)", "SelD PDAC-selective (D_mixed)")
ep_short <- c(
  "E1 PDAC pan-dependency"              = "E1 PDAC",
  "E4 CRC pan-dependency"               = "E4 CRC",
  "E5 clinical concordance"             = "E5 clin",
  "SelA PDAC-selective (A_ratio)"       = "SelA",
  "SelB PDAC-selective (B_zeffect==C)"  = "SelB",
  "SelD PDAC-selective (D_mixed)"       = "SelD")

# =====================================================================
# Panel a | candidate-by-layer evidence-availability matrix (REAL)
#   Built from evidence_layers.json per-gene per-layer present-flags:
#   a gene x 9 FEATURES-layer binary availability matrix (1 = annotated,
#   0 = missing). The two source CSVs (fig3a_candidate_layer_availability.csv,
#   fig3a_layer_coverage.csv) are computed by _build_fig3a_availability.py from
#   the real upstream JSON -- no proxy. Rendered as a binary raster (genes
#   sorted by total coverage descending => left dense, right sparse) because a
#   literal 20,751-gene x 9-tile heatmap is invisibly dense; every cell is
#   derived from the real gene x layer present flags.
# =====================================================================
# gd is still required by panels c/d/e below -- keep this read.
gd <- read.csv(file.path(dat_dir, "fig3d_missingness_decomposition.csv"),
               stringsAsFactors = FALSE)

avail <- read.csv(file.path(dat_dir, "fig3a_candidate_layer_availability.csv"),
                  stringsAsFactors = FALSE)
covsum <- read.csv(file.path(dat_dir, "fig3a_layer_coverage.csv"),
                   stringsAsFactors = FALSE)
feat <- names(layer_short)                       # 9 FEATURES layers, in order
M <- as.matrix(avail[, feat, drop = FALSE])
n_genes <- nrow(M)
ord <- order(rowSums(M), decreasing = TRUE)       # sort genes by total coverage
M <- M[ord, , drop = FALSE]
# long form for geom_raster: x = gene index, y = layer, fill = present(0/1)
xv <- rep(seq_len(n_genes), length(feat))
yv <- rep(seq_along(feat), each = n_genes)
fv <- as.vector(M)                                # column-major: layer1 all genes, ...
df_av <- data.frame(
  x = xv,
  y = factor(yv, levels = seq_along(feat), labels = feat),
  present = factor(fv, levels = c(0, 1), labels = c("0 missing", "1 annotated")))
cov_pct <- setNames(covsum$pct_present, covsum$layer)
ylab_a <- paste0(layer_short[feat], "\n(", sprintf("%.1f", cov_pct[feat]), "%)")

pa <- ggplot(df_av, aes(x = x, y = y, fill = present)) +
  geom_raster() +
  scale_fill_manual(values = c("0 missing" = "#f4f4f4", "1 annotated" = "#2166ac"),
                    name = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(expand = c(0, 0), labels = ylab_a) +
  labs(x = paste0("candidates (N = ", n_genes, ", sorted by total coverage, high to low)"),
       y = NULL,
       title = "a  Candidate-by-layer evidence availability\n(1 = annotated, 0 = missing)") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_blank(),            # 20,751 genes: no per-gene ticks
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7.5, lineheight = 0.85),
        legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7.5))

# =====================================================================
# Panel b | missingness patterns across layers (fig3a, real)
# =====================================================================
ga <- read.csv(file.path(dat_dir, "fig3a_sentinel_missingness.csv"),
               stringsAsFactors = FALSE)
ga$layer_id <- c("string_centrality","mutation_freq","impc_animal_ko","genetic_constraint",
                 "cancer_driver","ot_genetics_pdac","druggability","hpa_pdac_prognostic",
                 "hpa_rna_tissue_spec")
ga$layer_id <- factor(ga$layer_id, levels = names(layer_short))

pb <- ggplot(ga, aes(x = layer_id, y = pct_genes_at_sentinel)) +
  geom_bar(stat = "identity", width = 0.65, fill = "#e07a5f", colour = "black", linewidth = 0.25) +
  geom_hline(yintercept = 70, linetype = "dashed", colour = "#404040", linewidth = 0.4) +
  scale_x_discrete(labels = layer_short) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 110)) +
  labs(x = NULL, y = "% genes without sentinel",
       title = "b  Missingness patterns across layers\n(higher = less available)") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8))

# =====================================================================
# Panel c | alternative missingness conditions (fig3d main rows, real)
# =====================================================================
main <- gd[gd$scorer %in% c("values+missingness (frozen harmonic)", "values-only",
                            "missingness-only", "coverage count") &
            !grepl("\\[", gd$endpoint), ]
main$cond <- factor(main$scorer,
  levels = c("values-only", "missingness-only", "coverage count",
             "values+missingness (frozen harmonic)"),
  labels = c("values", "missing", "coverage", "frozen"))
main$ep <- factor(main$endpoint, levels = ep_levels)

pc <- ggplot(main, aes(x = ep, y = auroc, fill = cond)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("values" = "#bcb5b2",
                               "missing" = "#e76f51",
                               "coverage" = "#8e9aaf",
                               "frozen" = "#264653"),
                    name = "condition") +
  scale_x_discrete(labels = ep_short) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "AUROC",
       title = "c  Alternative missingness conditions\n(frozen harmonic family)") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8.5),
        legend.position = "bottom", legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 7.5),
        legend.margin = margin(t = 3, unit = "pt"),
        legend.box = "horizontal")

# =====================================================================
# Panel d | Delta under matched missingness: full minus values-only (real)
# =====================================================================
froz <- main[main$scorer == "values+missingness (frozen harmonic)", c("endpoint","auroc")]
vo   <- main[main$scorer == "values-only", c("endpoint","auroc")]
dd <- merge(froz, vo, by = "endpoint", suffixes = c("_f","_v"))
dd$delta <- dd$auroc_f - dd$auroc_v
dd$ep <- factor(dd$endpoint, levels = ep_levels)

pd <- ggplot(dd, aes(x = ep, y = delta, fill = delta > 0)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.3f", delta)), vjust = ifelse(dd$delta > 0, -0.4, 1.4),
            size = 2.8, family = "sans") +
  scale_fill_manual(values = c("TRUE" = "#2a9d8f", "FALSE" = "#c0392b"), guide = "none") +
  scale_x_discrete(labels = ep_short) +
  scale_y_continuous(limits = c(-0.08, 0.10), breaks = seq(-0.08, 0.10, 0.04)) +
  labs(x = NULL, y = expression(Delta * " AUROC (frozen - values-only)"),
       title = "d  Delta under matched condition 1:\nmissingness-aware vs values-only") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8.5))

# =====================================================================
# Panel e | Delta under matched missingness: values-only minus coverage (real)
# =====================================================================
cc <- main[main$scorer == "coverage count", c("endpoint","auroc")]
ee <- merge(vo, cc, by = "endpoint", suffixes = c("_v","_c"))
ee$delta <- ee$auroc_v - ee$auroc_c
ee$ep <- factor(ee$endpoint, levels = ep_levels)

pe <- ggplot(ee, aes(x = ep, y = delta, fill = delta > 0)) +
  geom_bar(stat = "identity", width = 0.6, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.3f", delta)), vjust = ifelse(ee$delta > 0, -0.4, 1.4),
            size = 2.8, family = "sans") +
  scale_fill_manual(values = c("TRUE" = "#2a9d8f", "FALSE" = "#c0392b"), guide = "none") +
  scale_x_discrete(labels = ep_short) +
  scale_y_continuous(limits = c(-0.06, 0.14), breaks = seq(-0.06, 0.14, 0.05)) +
  labs(x = NULL, y = expression(Delta * " AUROC (values-only - coverage count)"),
       title = "e  Delta under matched condition 2:\nactual values vs coverage only") +
  theme_nature(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        axis.text.x = element_text(size = 8.5))

# =====================================================================
# Panel f | representation audit: source provenance vs model-input
#          representation  (clearly-labelled conceptual schematic)
# =====================================================================
pf <- ggplot() +
  # two distinct audit axes (top row)
  annotate("rect", xmin = 0.3, xmax = 4.5, ymin = 5.8, ymax = 8.8,
           fill = "#dbe9f4", colour = "black", linewidth = 0.4) +
  annotate("text", x = 2.4, y = 8.1, label = "Source provenance",
           size = 3.1, fontface = "bold", family = "sans") +
  annotate("text", x = 2.4, y = 6.9, label = "independent source\nvs dependent / circular source",
           size = 2.7, family = "sans", lineheight = 1.0) +
  annotate("rect", xmin = 5.5, xmax = 9.7, ymin = 5.8, ymax = 8.8,
           fill = "#fce5cd", colour = "black", linewidth = 0.4) +
  annotate("text", x = 7.6, y = 8.1, label = "Model-input representation",
           size = 3.1, fontface = "bold", family = "sans") +
  annotate("text", x = 7.6, y = 6.9, label = "scale coding\nand missingness encoding",
           size = 2.7, family = "sans", lineheight = 1.0) +
  annotate("segment", x = 4.55, y = 7.3, xend = 5.45, yend = 7.3,
           arrow = arrow(length = unit(0.12, "cm")), colour = "#555") +
  # the audit point (bottom, full width)
  annotate("rect", xmin = 0.3, xmax = 9.7, ymin = 1.2, ymax = 4.6,
           fill = "#fff2cc", colour = "#b8860b", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 5.0, y = 3.9, label = "Two distinct audit questions",
           size = 3.0, fontface = "bold", family = "sans") +
  annotate("text", x = 5.0, y = 2.4,
           label = "Same provenance can give different representations,\nand different provenance can give the same representation;\nboth axes must be audited.",
           size = 2.6, family = "sans", lineheight = 1.05) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.6, 9.4), clip = "off") +
  labs(title = "f  Representation audit: provenance vs model-input",
       x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust = 0, size = 10.5, face = "bold", lineheight = 1.1),
        plot.margin = margin(2, 2, 2, 2, "pt"))

# ---- combine -------------------------------------------------------------
fig <- (pa + pb + pc) / (pd + pe + pf) + plot_layout(heights = c(1.15, 1))
fig <- fig + plot_annotation(
  title = "Fig. 3 Evidence availability and representation affect the measured increment",
  theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, family = "sans"))
)

save_fig(fig, "Fig3_evidence_representation", w = 14, h = 10)
