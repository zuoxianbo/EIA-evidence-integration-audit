# --- EIA reproducibility: paths are resolved relative to the repository root ---
loc <- function() {
  # this file lives in <repo>/code/R ; repo root is two levels up
  p <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
  if (is.na(p) || !nzchar(p)) p <- getwd()
  dirname(dirname(p))
}
ROOT <- loc()

# Fig 5 (v26 numbering): per-disease-type full-dimension alignment matrix (v33 pre-submission audit)
base_dir <- file.path(ROOT, "code", "R")
dat_dir <- file.path(ROOT, "data", "source_data")
source(file.path(base_dir, "common.R"))
library(patchwork)
library(ggrepel)

sd <- dat_dir
outd <- file.path(PKG, "figures")

# ---------------- panel a: verdict heatmap ----------------
m <- read.csv(file.path(sd, "fig5a_alignment_matrix.csv"), stringsAsFactors = FALSE)

verdict_provenance <- function(o) ifelse(o == 0, "PASS", ifelse(o >= 1, "FAIL", "WARN"))
verdict_representation <- function(s) ifelse(is.na(s), "N/A",
                                    ifelse(abs(s) <= 0.02, "PASS",
                                    ifelse(abs(s) <= 0.05, "WARN", "FAIL")))
verdict_aggregation <- function(d, circ) ifelse(d <= 0, "PASS",
                                                ifelse(circ == "TRUE" | circ == "True" | circ == TRUE, "FAIL", "GAIN?"))
verdict_baseline <- function(z) ifelse(z == "" | is.na(z), "N/A", ifelse(as.numeric(z) <= 2, "PASS",
                                                   ifelse(as.numeric(z) <= 3, "WARN", "FAIL")))
verdict_delta <- function(d, lo, hi, circ) {
  lo_n <- suppressWarnings(as.numeric(lo))
  hi_n <- suppressWarnings(as.numeric(hi))
  ifelse(is.na(lo_n), ifelse(d > 0 & (circ == "TRUE" | circ == "True" | circ == TRUE), "FAIL",
                             ifelse(d <= 0, "PASS", "WARN")),
  ifelse(hi_n <= 0, "PASS",
  ifelse(lo_n > 0 & (circ == "TRUE" | circ == "True" | circ == TRUE), "FAIL",
  ifelse(lo_n > 0, "GAIN?", "WARN"))))
}

short <- c(
  "E1 PDAC pan-dependency" = "E1 PDAC pan-dep",
  "E4 CRC pan-dependency" = "E4 CRC pan-dep\n(zero-shot)",
  "E5 clinical concordance" = "E5 clinical",
  "SelA PDAC-selective (A_ratio)" = "SelA A-ratio",
  "SelB z-effect (==C)" = "SelB z-effect",
  "SelC PDAC-selective (C_lineage_adj)" = "SelC lineage-adj",
  "SelD PDAC-selective (D_mixed)" = "SelD mixed",
  "E7 genetic-association reuse" = "E7 OT-reuse\n(mech 2a)",
  "E8 interaction label" = "E8 interaction\n(mech 2b, PC)")

vm <- data.frame(
  endpoint = factor(m$endpoint, levels = rev(m$endpoint)),
  provenance = verdict_provenance(m$overlap_score),
  representation = verdict_representation(m$harmonic_shift_repr),
  aggregation = verdict_aggregation(m$delta, as.character(m$circular)),
  baseline_family = verdict_baseline(as.character(m$degree_z_harmonic)),
  delta_identified = verdict_delta(m$delta, m$delta_ci_lo, m$delta_ci_hi, as.character(m$circular))
)

# long format
vml <- do.call(rbind, lapply(c("provenance", "representation", "aggregation",
                               "baseline_family", "delta_identified"), function(d) {
  data.frame(endpoint = vm$endpoint, dimension = d,
             verdict = vm[[d]])
}))
vml$dimension <- factor(vml$dimension,
  levels = c("provenance", "representation", "aggregation", "baseline_family", "delta_identified"),
  labels = c("Provenance\n(overlap)", "Representation\n(sentinel)", "Aggregation\n(Delta > 0)",
             "Baseline family\n(degree null)", "Delta estimand\n(CI)"))

vcol <- c("PASS" = "#2e7d32", "WARN" = "#f9a825", "FAIL" = "#c62828",
          "GAIN?" = "#1565c0", "N/A" = "#bdbdbd")

pa <- ggplot(vml, aes(x = dimension, y = endpoint, fill = verdict)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  scale_fill_manual(values = vcol, drop = FALSE,
    name = "Audit verdict") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "a  Per-disease-type full-dimension alignment: same audit, every endpoint") +
  theme_nature(base_size = 10) +
  theme(axis.text.x.top = element_text(angle = 0, size = 8.2, hjust = 0.5, vjust = 0),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.35, "cm"))

# ---------------- panel b: Delta volcano under Benjamini-Hochberg correction (v33 m=60) ----------------
f <- read.csv(file.path(sd, "fig5b_delta_forest.csv"), stringsAsFactors = FALSE)
f$ci_lo <- suppressWarnings(as.numeric(f$ci_lo))
f$ci_hi <- suppressWarnings(as.numeric(f$ci_hi))
f$delta <- suppressWarnings(as.numeric(f$delta))
f$bh_q <- suppressWarnings(as.numeric(as.character(f$bh_q)))
f$neglog10q <- -log10(pmax(f$bh_q, 1e-12))
# endpoint / scorer short labels for the three survivors
ep_short <- c(
  "E1 PDAC pan-dependency" = "E1",
  "E4 CRC pan-dependency" = "E4",
  "E5 clinical concordance" = "E5",
  "SelA PDAC-selective (A_ratio)" = "SelA",
  "SelB PDAC-selective (B_zeffect==C)" = "SelB",
  "SelD PDAC-selective (D_mixed)" = "SelD")
scorer_short <- c(
  "Harmonic mean" = "Harm", "ECS (multiplicative)" = "ECS", "Additive" = "Add",
  "Geometric mean" = "Geo", "Arithmetic mean" = "Ari", "Rank aggregation" = "Rank",
  "Weighted rank" = "WRank", "Random forest" = "RF",
  "Logistic regression" = "LR", "Elastic net" = "EN")
f$lab <- paste0(ep_short[as.character(f$endpoint)], " \u00b7 ",
                scorer_short[as.character(f$scorer)])
f$surv <- ifelse(f$bh_q < 0.05 & f$delta > 0, "survives BH (q<0.05)", "not significant after BH")
f$lbl <- ifelse(f$surv == "survives BH (q<0.05)", f$lab, "")
# manual vertical nudge for the two nearly-overlapping RF survivors
f$nudge_y <- 0
f$nudge_y[f$lbl == "E1 \u00b7 RF"] <- 0.22
f$nudge_y[f$lbl == "E4 \u00b7 RF"] <- -0.22
pb <- ggplot(f, aes(x = delta, y = neglog10q, colour = surv)) +
  geom_point(size = 1.9, stroke = 0.5, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_text_repel(aes(label = lbl), nudge_y = f$nudge_y, size = 3.2, colour = "#0d47a1",
                  min.segment.length = 0, box.padding = 0.22,
                  point.padding = 0.15, force = 2.5, force_pull = 0.3,
                  max.overlaps = Inf, seed = 42) +
  annotate("text", x = max(f$delta, na.rm = TRUE), y = -log10(0.05),
           vjust = -0.5, hjust = 1, size = 3.0, colour = "grey45",
           label = "BH 5% FDR (q = 0.05)") +
  scale_colour_manual(values = c("survives BH (q<0.05)" = "#1565c0",
                                 "not significant after BH" = "#c62828"),
                      name = "BH 5% FDR") +
  labs(x = expression(Delta * " (composite - best single layer, AUROC)"),
       y = expression(-log[10] * " (Benjamini-Hochberg q)"),
       title = "b  Full frozen Delta family: 10 integrators x 6 real endpoints (m = 60)\nunder Benjamini-Hochberg: only four survive 5% FDR") +
  theme_nature(base_size = 9) +
  theme(legend.position = "bottom", panel.grid = element_blank())

# ---------------- panel c: method x endpoint AUROC heatmap ----------------
h <- read.csv(file.path(sd, "fig5c_method_auroc_matrix.csv"), stringsAsFactors = FALSE,
              check.names = FALSE)
rownames(h) <- h$method
h <- h[, -1]
h <- as.matrix(h)
class(h) <- "numeric"

# order methods: single layers, composites, supervised
meth_order <- c("STRING centrality", "Mutation frequency", "IMPC animal KO",
                "Genetic constraint", "Cancer-driver annotation", "OT genetics (PDAC)",
                "Druggability", "HPA PDAC prognostic", "HPA RNA tissue spec",
                "ECS (multiplicative)", "Additive", "Geometric mean", "Harmonic mean",
                "Arithmetic mean", "Rank aggregation", "Weighted rank",
                "Logistic regression", "Elastic net", "Random forest")
meth_order <- meth_order[meth_order %in% rownames(h)]
h <- h[meth_order, , drop = FALSE]

hl <- do.call(rbind, lapply(seq_len(nrow(h)), function(i) {
  do.call(rbind, lapply(seq_len(ncol(h)), function(j) {
    data.frame(method = rownames(h)[i], endpoint = colnames(h)[j],
               auroc = h[i, j])
  }))
}))
hl$method <- factor(hl$method, levels = rev(meth_order))
hl$endpoint <- factor(hl$endpoint, levels = colnames(h))
hl$auroc_txt <- ifelse(is.na(hl$auroc), "", sprintf("%.2f", hl$auroc))

# x uses the endpoint FACTOR (correct CSV order) + a label lookup, so that no two
# endpoints with an unmapped short label can collapse into a single "NA" column.
pc <- ggplot(hl[!is.na(hl$auroc), ], aes(x = endpoint, y = method, fill = auroc)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = auroc_txt), size = 2.1, family = "sans",
            colour = ifelse(hl[!is.na(hl$auroc), "auroc"] > 0.82, "white", "black")) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2166ac", limits = c(0.5, 0.95),
                      oob = scales::squish,
                      name = "AUROC", na.value = "grey92") +
  scale_x_discrete(position = "top", labels = short) +
  labs(x = NULL, y = NULL,
       title = "c  Method x disease endpoint: unsupervised composites, single layers, supervised OOF") +
  theme_nature(base_size = 9) +
  theme(axis.text.x.top = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 7.5, lineheight = 0.85),
        panel.grid = element_blank())

fig <- pa + pb + pc + plot_layout(ncol = 1, heights = c(0.9, 1.4, 2.1))
save_fig(fig, "Fig5_external_evaluation", w = 11.5, h = 15, out_dir = outd)
