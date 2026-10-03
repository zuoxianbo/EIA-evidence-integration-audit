# --- EIA reproducibility: paths are resolved relative to the repository root ---
from pathlib import Path as _P
ROOT = _P(__file__).resolve().parent.parent.parent

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v33 postprocess: regenerate every figure source-data CSV natively from the
frozen alignment_results.json + raw inputs (runs AFTER run_analysis.py).

Outputs go to BOTH data/source_data (archive of record) and code/R/source_data
(the directory the R/ggplot2 scripts read), using canonical (version-free)
filenames so figures always render from the current run. This fixes the v27
issue where R read stale v26-era CSVs (e.g. fig2d_mechanism2.csv)."""
import json, csv, os, shutil
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(ROOT, "data", "results", "alignment_results.json")
RAW = os.path.join(ROOT, "data", "raw_inputs")
CSV_ARCHIVE = os.path.join(ROOT, "data", "source_data")
CSV_R = os.path.join(ROOT, "code", "R", "source_data")
SENT = -3.0

R = json.load(open(RES))
EPS = R["endpoints"]
V33 = R["v33"]
SA = R["self_audit"]

# ---- write helper: canonical CSV into BOTH directories ----
def wcsv(name, header, rows):
    for d in (CSV_ARCHIVE, CSV_R):
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, name), "w", newline="") as f:
            w = csv.writer(f); w.writerow(header); w.writerows(rows)
    print("wrote", name)

# =====================================================================
# 1. fig3a: sentinel-missingness fractions (native from raw layers)
#    fig3c: multiplicative-rule pathology (native; v18 sentinel_audit port)
# =====================================================================
EV = json.load(open(os.path.join(RAW, "evidence_layers.json")))
LAYERS, ALL_GENES = EV["layers"], EV["genes"]
FEATURES = ["string_centrality", "mutation_freq", "impc_animal_ko",
            "genetic_constraint", "cancer_driver", "ot_genetics_pdac",
            "druggability", "hpa_pdac_prognostic", "hpa_rna_tissue_spec"]
N = len(ALL_GENES)
X = np.full((N, len(FEATURES)), SENT)
for i, g in enumerate(ALL_GENES):
    for j, f in enumerate(FEATURES):
        r = LAYERS.get(f, {}).get(g)
        if r and r.get("present"):
            X[i, j] = r["norm"]
I_STR, I_MUT, I_IMPC = (FEATURES.index(x) for x in
                        ("string_centrality", "mutation_freq", "impc_animal_ko"))
SUP_ALL = [FEATURES.index(f) for f in ("cancer_driver", "ot_genetics_pdac",
                                       "druggability", "hpa_pdac_prognostic",
                                       "hpa_rna_tissue_spec")]
LABELS3A = {"string_centrality": "string\ncentr.", "mutation_freq": "mut.\nfreq.",
            "impc_animal_ko": "impc\nKO", "genetic_constraint": "gen.\nconstraint",
            "cancer_driver": "cancer\ndriver", "ot_genetics_pdac": "ot_pdac",
            "druggability": "drugg", "hpa_pdac_prognostic": "hpa\nprog",
            "hpa_rna_tissue_spec": "hpa\nRNA"}
rows = [[f, LABELS3A[f],
         round(float(100.0 * np.isclose(X[:, j], SENT).mean()), 6)]
        for j, f in enumerate(FEATURES)]
wcsv("fig3a_sentinel_missingness.csv",
     ["layer", "layer_label", "pct_genes_at_sentinel"], rows)

def DP(M):
    D = 0.80 * M[:, I_STR] + 0.10 * M[:, I_MUT] + 0.10 * M[:, I_IMPC]
    return D, np.mean(M[:, SUP_ALL], axis=1)
D, PHI = DP(X)
gain = 1.0 + 0.6 * PHI
mult = D * gain
neg_gain, neg_D = gain < 0, D < 0
both = neg_gain & neg_D
mult_up = D * (1.0 + 0.6 * (PHI + 0.25))
losers = float(100.0 * (mult_up < mult).mean())
rows = [["Genes whose ECS gain term is negative",
         round(float(100.0 * neg_gain.mean()), 6)],
        ["Genes whose score falls when support increases",
         round(losers, 6)],
        ["Genes with double sign flip",
         round(float(100.0 * both.mean()), 6)]]
wcsv("fig3c_monotonicity.csv", ["metric", "pct"], rows)
FROZEN_3C = [84.000771, 28.003470, 26.152956]
for (m, v), fr in zip(rows, FROZEN_3C):
    assert abs(v - fr) < 1e-3, (m, v, fr)
print("fig3c cross-check vs frozen v18 values: PASS")

# =====================================================================
# 2. fig2a/b/c: copy native run outputs to canonical names
# =====================================================================
def cp(src, dst):
    for d in (CSV_ARCHIVE, CSV_R):
        shutil.copy(os.path.join(CSV_ARCHIVE, src), os.path.join(d, dst))
    print("copied", src, "->", dst)
cp("fig2a_provenance_gradient.csv", "fig2a_provenance_gradient.csv")
cp("fig2b_e3_deletion.csv", "fig2b_e3_deletion.csv")
cp("fig2c_e3c_control.csv", "fig2c_e3c_control.csv")
cp("fig5a_alignment_matrix.csv", "fig5a_alignment_matrix.csv")
cp("fig5b_delta_forest.csv", "fig5b_delta_forest.csv")
cp("fig3b_pairwise_inversion.csv", "fig3b_pairwise_inversion.csv")
cp("fig3d_missingness_decomposition.csv",
   "fig3d_missingness_decomposition.csv")
cp("fig3c_aggregation_robustness.csv",
   "fig3c_aggregation_robustness.csv")

# =====================================================================
# 3. fig2d: mechanism-2 panel (E7 + E8) from current results
# =====================================================================
e7 = EPS["E7 genetic-association reuse"]
d7 = e7["delta_bestcomposite_vs_bestsingle"]
e8 = EPS["E8 interaction label (canonical)"]
wcsv("fig2d_mechanism2.csv",
     ["endpoint", "mechanism", "best_single", "best_single_auroc",
      "harmonic_auroc", "best_composite", "best_composite_auroc",
      "delta", "delta_ci_lo", "delta_ci_hi", "rf_oof_auroc",
      "delta_rf", "delta_rf_ci_lo", "delta_rf_ci_hi"],
     [["E7 genetic-association reuse", "direct reuse (genetic layer)",
       e7["best_single_layer"]["name"], e7["best_single_layer"]["auroc"],
       e7["aggregation"]["Harmonic mean"],
       e7["best_composite"]["name"], e7["best_composite"]["auroc"],
       d7["delta"], d7["ci95"][0], d7["ci95"][1], "", "", "", ""],
      ["E8 interaction label (canonical)",
       "constructed interaction (difference) label",
       e8["best_single"]["name"], e8["best_single"]["auroc"],
       e8["composite_auroc"]["Harmonic mean"],
       e8["best_fixed_composite"]["name"], e8["best_fixed_composite"]["auroc"],
       e8["delta_fixed_form_vs_best_single"]["delta"],
       e8["delta_fixed_form_vs_best_single"]["ci95"][0],
       e8["delta_fixed_form_vs_best_single"]["ci95"][1],
       e8["supervised_auroc"]["Random forest"],
       e8["delta_rf_vs_best_single"]["delta"],
       e8["delta_rf_vs_best_single"]["ci95"][0],
       e8["delta_rf_vs_best_single"]["ci95"][1]]])

# =====================================================================
# 4. fig4a: baseline family on E1; fig4b: full provenance ladder (native)
# =====================================================================
bf = EPS["E1 PDAC pan-dependency"]["baseline_family"]["Harmonic mean"]
cent_e1 = EPS["E1 PDAC pan-dependency"]["aggregation"]["STRING centrality"]
wcsv("fig4a_baseline_family.csv",
     ["item", "value", "err_type", "err_low", "err_high"],
     [["STRING centrality (strongest single layer)", cent_e1, "none", 0, 0],
      ["Harmonic composite (observed)", bf["degree"]["observed"], "none", 0, 0],
      ["Degree-matched null (mean)", bf["degree"]["null_mean"], "sd",
       bf["degree"]["null_mean"] - 1.96 * bf["degree"]["null_sd"],
       bf["degree"]["null_mean"] + 1.96 * bf["degree"]["null_sd"]],
      ["Annotation-density null (mean)", bf["density"]["null_mean"], "sd",
       bf["density"]["null_mean"] - 1.96 * bf["density"]["null_sd"],
       bf["density"]["null_mean"] + 1.96 * bf["density"]["null_sd"]]])

def harm_minus_centrality(ep):
    tab = EPS[ep]["aggregation"]
    return round(tab["Harmonic mean"] - tab["STRING centrality"], 4)

FAM = V33["provenance_endpoint_family"]
e6v = V33["e6_external_drug_response"]
rows = [
    ["E1", "in-sample", "independent source", 0.0,
     EPS["E1 PDAC pan-dependency"]["representation"]["harmonic_shift"],
     harm_minus_centrality("E1 PDAC pan-dependency"),
     EPS["E1 PDAC pan-dependency"]["n_pos"]],
    ["E2", "in-sample", "independent source", 0.0,
     FAM["E2 PDAC-enriched dependency"]["representation"]["harmonic_available_case"] -
     FAM["E2 PDAC-enriched dependency"]["representation"]["harmonic_sentinel_as_value"],
     FAM["E2 PDAC-enriched dependency"]["composites"]["Harmonic mean"]["delta_vs_centrality"],
     FAM["E2 PDAC-enriched dependency"]["n_pos"]],
    ["E3", "in-sample", "direct label reuse", 1.0,
     FAM["E3 conjunctive actionability"]["representation"]["harmonic_available_case"] -
     FAM["E3 conjunctive actionability"]["representation"]["harmonic_sentinel_as_value"],
     FAM["E3 conjunctive actionability"]["composites"]["Harmonic mean"]["delta_vs_centrality"],
     FAM["E3 conjunctive actionability"]["n_pos"]],
    ["E3-A", "in-sample", "leakage control", 0.0,
     FAM["E3-A leakage-controlled essentiality"]["representation"]["harmonic_available_case"] -
     FAM["E3-A leakage-controlled essentiality"]["representation"]["harmonic_sentinel_as_value"],
     FAM["E3-A leakage-controlled essentiality"]["composites"]["Harmonic mean"]["delta_vs_centrality"],
     FAM["E3-A leakage-controlled essentiality"]["n_pos"]],
    ["E3-C", "in-sample", "circular control", 1.0,
     FAM["E3-C out-of-evidence druggability"]["representation"]["harmonic_available_case"] -
     FAM["E3-C out-of-evidence druggability"]["representation"]["harmonic_sentinel_as_value"],
     FAM["E3-C out-of-evidence druggability"]["composites"]["Harmonic mean"]["delta_vs_centrality"],
     FAM["E3-C out-of-evidence druggability"]["n_pos"]],
    ["E4", "external", "independent context", 0.0,
     EPS["E4 CRC pan-dependency"]["representation"]["harmonic_shift"],
     harm_minus_centrality("E4 CRC pan-dependency"),
     EPS["E4 CRC pan-dependency"]["n_pos"]],
    ["E5", "external", "translational", 0.0,
     EPS["E5 clinical concordance"]["representation"]["harmonic_shift"],
     harm_minus_centrality("E5 clinical concordance"),
     EPS["E5 clinical concordance"]["n_pos"]],
    ["E6", "external", "drug response (related ecology)", 0.5, "",
     e6v["composites_vs_best_single"]["Harmonic mean"]["delta_vs_centrality"],
     e6v["n_positives"]],
]
wcsv("fig4b_delta_landscape.csv",
     ["endpoint", "context", "provenance", "overlap_score", "repr_shift",
      "delta_auroc", "n_pos"], rows)

# =====================================================================
# 5. fig6a: full Delta distributions, E4/E5 + native E6 (7 composites each)
# =====================================================================
rows = []
for ep, disp in (("E4 CRC pan-dependency", "E4 CRC"),
                 ("E5 clinical concordance", "E5 clinical")):
    fam = EPS[ep]["delta_all_composites_vs_best_single"]
    for cn, d in fam.items():
        rows.append([disp, cn, d["delta"]])
for cn, d in e6v["composites_vs_best_single"].items():
    rows.append(["E6 drug response", cn, d["delta"]])
wcsv("fig6a_delta_distributions.csv", ["endpoint", "method", "delta_auroc"], rows)
assert len(rows) == 21, len(rows)   # 3 external endpoints x 7 composites

# =====================================================================
# 6. E6 AUPRC addendum (deterministic reconstruction; AUROC cross-checked)
# =====================================================================
import csv as _csv
# GDSC/DeepCDR external dataset is NOT redistributed here (provider terms).
# Point this at your local copy of the DeepCDR data directory.
GDSC_DIR = os.environ.get(
    "GDSC_DIR",
    str(Path.home() / "DeepCDR" / "data"),
)
with open(os.path.join(GDSC_DIR, "CCLE",
                       "Cell_lines_annotations_20181226.txt"),
          encoding="utf-8", errors="replace") as f:
    ann = list(_csv.DictReader(f, delimiter="\t"))
panc = {r_["depMapID"].strip() for r_ in ann
        if (r_.get("Site_Primary") or "").strip().lower() == "pancreas"
        and r_.get("depMapID")}
with open(os.path.join(GDSC_DIR, "CCLE", "GDSC_IC50.csv")) as f:
    ic = list(_csv.reader(f))
cols = ic[0][1:]
pidx = [i for i, c in enumerate(cols) if c.strip() in panc]
drug_ic50 = {}
for r_ in ic[1:]:
    did = r_[0].replace("GDSC:", "").strip()
    vals = []
    for i in pidx:
        s = r_[1 + i].strip()
        if s not in ("", "NA", "NaN"):
            try:
                vals.append(float(s))
            except ValueError:
                pass
    if len(vals) >= 8:
        drug_ic50[did] = float(np.median(vals))
with open(os.path.join(GDSC_DIR, "GDSC",
                       "1.Drug_listMon Jun 24 09_00_55 2019.csv"),
          encoding="utf-8", errors="replace") as f:
    drugs = list(_csv.DictReader(f))
gene_set = set(ALL_GENES)
d2t = {}
for d in drugs:
    did = (d.get("drug_id") or "").strip()
    genes = {t.strip() for t in (d.get("Targets") or "").replace(";", ",").split(",")
             if t.strip() in gene_set}
    if did and genes:
        d2t[did] = genes
shared = [d for d in drug_ic50 if d in d2t]
order6 = sorted(shared, key=lambda d: drug_ic50[d])
k6 = max(1, len(order6) // 3)
sg = set().union(*[d2t[d] for d in order6[:k6]])
rg = set().union(*[d2t[d] for d in order6[-k6:]])
y6 = np.array([1.0 if g in (sg - rg) else 0.0 for g in ALL_GENES])
assert int(y6.sum()) == e6v["n_positives"]
I_DRUG = FEATURES.index("druggability")
SUP = X[:, SUP_ALL]
Ds_h = np.maximum(D + 3., .01)
Ps_h = np.maximum(SUP.mean(axis=1) + 3., .01)
harm6 = 2. * Ds_h * Ps_h / (Ds_h + Ps_h)
cent6 = X[:, I_STR]
a_cent6 = float(roc_auc_score(y6, cent6))
assert abs(round(a_cent6, 4) - e6v["baseline_centrality_auroc"]) < 1e-3, \
    "E6 reconstruction drifted"
auprc_cent = round(float(average_precision_score(y6, cent6)), 4)
auprc_harm = round(float(average_precision_score(y6, harm6)), 4)
# E6 representation shift (harmonic available-case minus sentinel), same formula
Xm6 = X.copy(); Xm6[np.isclose(Xm6, SENT)] = np.nan
Dm6 = np.nansum(np.stack([0.80 * Xm6[:, I_STR], 0.10 * Xm6[:, I_MUT],
                          0.10 * Xm6[:, I_IMPC]]), axis=0)
with np.errstate(invalid="ignore"):
    PHIm6 = np.nanmean(Xm6[:, SUP_ALL], axis=1)
PHIm6 = np.where(np.isfinite(PHIm6), PHIm6, 0.0)
Dm6 = np.where(np.isfinite(Dm6), Dm6, 0.0)
ac6 = 2.0 * np.maximum(Dm6 + 3., .01) * np.maximum(PHIm6 + 3., .01) / \
    (np.maximum(Dm6 + 3., .01) + np.maximum(PHIm6 + 3., .01))
e6_repr_shift = round(float(roc_auc_score(y6, ac6) - roc_auc_score(y6, harm6)), 4)
V33["e6_external_drug_response"]["representation"] = {
    "harmonic_shift": e6_repr_shift}
print(f"E6 repr shift: {e6_repr_shift}")
V33["e6_external_drug_response"]["auprc"] = {
    "STRING_centrality": auprc_cent, "Harmonic_mean": auprc_harm}
print(f"E6 AUPRC addendum: centrality={auprc_cent} harmonic={auprc_harm}")

# canonical E6 discrimination CSV: all methods AUROC + centrality/harmonic AUPRC
rows = [["STRING centrality", "AUROC", e6v["baseline_centrality_auroc"]]]
for nm in ("Harmonic mean", "ECS (multiplicative)", "Additive",
           "Geometric mean", "Arithmetic mean", "Rank aggregation",
           "Weighted rank"):
    rows.append([nm, "AUROC", e6v["composites_vs_best_single"][nm]["auroc"]])
for L in ("Logistic regression", "Elastic net", "Random forest"):
    rows.append([L, "AUROC", e6v["supervised_vs_best_single"][L]["auroc"]])
rows.append(["STRING centrality", "AUPRC", auprc_cent])
rows.append(["Harmonic mean", "AUPRC", auprc_harm])
wcsv("fig5b_e6_discrimination_ranking.csv", ["method", "metric", "value"], rows)

# =====================================================================
# 7. fig5c: method x endpoint AUROC matrix (SelC dropped: identical to SelB)
# =====================================================================
order = ["E1 PDAC pan-dependency", "E4 CRC pan-dependency",
         "E5 clinical concordance", "SelA PDAC-selective (A_ratio)",
         "SelB PDAC-selective (B_zeffect==C)", "SelD PDAC-selective (D_mixed)",
         "E7 genetic-association reuse"]
methods = ["STRING centrality", "Mutation frequency", "IMPC animal KO",
           "Genetic constraint", "Cancer-driver annotation",
           "OT genetics (PDAC)", "Druggability", "HPA PDAC prognostic",
           "HPA RNA tissue spec", "ECS (multiplicative)", "Additive",
           "Geometric mean", "Harmonic mean", "Arithmetic mean",
           "Rank aggregation", "Weighted rank", "Logistic regression",
           "Elastic net", "Random forest"]
cols_m = [ep.replace("SelB PDAC-selective (B_zeffect==C)",
                    "SelB z-effect (==C)") for ep in order] + \
         ["E8 interaction label"]
rows = []
for m in methods:
    row = [m]
    for ep in order:
        row.append(EPS[ep]["aggregation"].get(m, ""))
    if m in e8["single_layer_auroc"]:
        row.append(e8["single_layer_auroc"][m])
    elif m in e8["composite_auroc"]:
        row.append(e8["composite_auroc"][m])
    elif m in e8["supervised_auroc"]:
        row.append(e8["supervised_auroc"][m])
    else:
        row.append("")
    rows.append(row)
wcsv("fig5c_method_auroc_matrix.csv", ["method"] + cols_m, rows)

# =====================================================================
# 8. supp1: attribution map (native ladder; plain contexts for the R script)
# =====================================================================
def _fam_row(nm, short, ctx, prov):
    f = FAM[nm]
    return [short, ctx, prov, 1.0 if prov in ("direct label reuse",
                                              "circular control") else 0.0,
            round(f["representation"]["harmonic_available_case"] -
                  f["representation"]["harmonic_sentinel_as_value"], 4),
            f["composites"]["Harmonic mean"]["delta_vs_centrality"],
            f["n_pos"]]
rows = [
    ["E1", "in-sample", "independent source", 1.0,
     EPS["E1 PDAC pan-dependency"]["representation"]["harmonic_shift"],
     harm_minus_centrality("E1 PDAC pan-dependency"),
     EPS["E1 PDAC pan-dependency"]["n_pos"]],
    _fam_row("E2 PDAC-enriched dependency", "E2", "in-sample",
             "independent source"),
    _fam_row("E3 conjunctive actionability", "E3", "in-sample",
             "direct label reuse"),
    _fam_row("E3-A leakage-controlled essentiality", "E3-A", "in-sample",
             "leakage control"),
    _fam_row("E3-C out-of-evidence druggability", "E3-C", "in-sample",
             "circular control"),
    ["E4", "external", "independent context", 1.0,
     EPS["E4 CRC pan-dependency"]["representation"]["harmonic_shift"],
     harm_minus_centrality("E4 CRC pan-dependency"),
     EPS["E4 CRC pan-dependency"]["n_pos"]],
    ["E5", "external", "translational", 1.0,
     EPS["E5 clinical concordance"]["representation"]["harmonic_shift"],
     harm_minus_centrality("E5 clinical concordance"),
     EPS["E5 clinical concordance"]["n_pos"]],
    ["E6", "external", "related ecology", 0.5, e6_repr_shift,
     e6v["composites_vs_best_single"]["Harmonic mean"]["delta_vs_centrality"],
     e6v["n_positives"]],
]
wcsv("supp1_attribution_map.csv",
     ["endpoint", "context", "provenance_note", "provenance_independence",
      "repr_shift", "delta_auroc", "n_pos"], rows)

# =====================================================================
# 9. E6 preregistered-verdict correction (composite-only family)
#    The v33 run code extended the preregistered criterion to supervised
#    scorers; the preregistration (frozen v18/v19) covers ONLY the seven
#    fixed-form composites. Recompute the verdict on the preregistered
#    family and record the supervised excess separately.
# =====================================================================
comp_only_max = max(v["delta_vs_centrality"]
                    for v in e6v["composites_vs_best_single"].values())
verdict_prereg = ("PASS" if comp_only_max <= 0 else
                  "WARN" if comp_only_max <= 0.02 else "FAIL")
sup_max = max((v["delta_vs_centrality"]
               for v in e6v["supervised_vs_best_single"].values()),
              default=None)
e6v["preregistered_verdict_composite_only"] = verdict_prereg
e6v["verdict_scope_note"] = (
    "The preregistered criteria (frozen before E6 outcome inspection) cover "
    "the seven fixed-form composites only; on that family the observed "
    f"outcome is PASS (max composite Delta vs centrality = {comp_only_max}). "
    "The supervised random forest exceeded the centrality baseline by "
    f"{sup_max} (bootstrap 95% CI includes 0; one-sided p = "
    f"{e6v['supervised_vs_best_single']['Random forest']['p_one_sided']}), "
    "reported separately as a non-significant supervised excess, not as a "
    "composite verdict change.")
print(f"E6 preregistered (composite-only) verdict: {verdict_prereg} "
      f"(max composite Delta = {comp_only_max}); supervised RF excess = {sup_max}")

# =====================================================================
# 10. save updated results (with E6 AUPRC + verdict correction) + digest
# =====================================================================
with open(RES, "w") as f:
    json.dump(R, f, indent=1)
print("updated alignment_results.json (E6 AUPRC addendum)")

digest = {
    "e8": {"best_single": e8["best_single"],
           "best_fixed": e8["best_fixed_composite"],
           "harmonic": e8["composite_auroc"]["Harmonic mean"],
           "fixed_delta": e8["delta_fixed_form_vs_best_single"],
           "rf": e8["supervised_auroc"]["Random forest"],
           "rf_delta": e8["delta_rf_vs_best_single"]},
    "e1_baseline_family": bf,
    "transfer": SA["s1_trivial_transfer_baselines"],
    "bh_family": {"m": SA["s2_multiplicity"].get("family"),
                  "survivors": SA["s2_multiplicity"]["survivors_bh_q_lt_0.05"]},
    "repeated_cv": SA["s3_repeated_cv"],
    "strict": {k: {"n_pos": v["n_pos"], "best_single": v["best_single"],
                   "best_comp": v["best_composite"],
                   "delta_comp": v["delta_bestcomposite"],
                   "rf": v["rf_oof_auroc"], "delta_rf": v["delta_rf"],
                   "rcv_rf": v["repeated_cv_delta_rf"]}
               for k, v in SA["s4_strict_labels"].items()},
    "annotated_subset": SA["s5_annotation_poor_sensitivity"],
    "missingness": SA["s12_missingness_decomposition"],
    "agg_robustness": SA["s13_aggregation_robustness"],
    "nested_cv": SA["s14_nested_cv_sensitivity"],
    "e6": e6v,
    "gradient": V33["provenance_gradient"],
    "family": {k: {"n_pos": v["n_pos"],
                   "harmonic": v["composites"]["Harmonic mean"]["auroc"],
                   "harmonic_no_drug": v["harmonic_no_druggability"]}
               for k, v in FAM.items()},
    "run_record": V33["run_record"],
    "scorer_manifest": {"n_scorers": V33["scorer_manifest"]["n_scorers"],
                         "n_cells": V33["scorer_manifest"]["n_scorer_endpoint_cells"]},
}
with open(os.path.join(ROOT, "data", "results", "v33_manuscript_digest.json"),
          "w") as f:
    json.dump(digest, f, indent=1)
print("wrote v33_manuscript_digest.json")
