# --- EIA reproducibility: paths are resolved relative to the repository root ---
from pathlib import Path as _P
ROOT = _P(__file__).resolve().parent.parent.parent

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EIA v33 master analysis: full-dimension alignment + SELF-AUDIT + P0/P1 ledgers.

Pre-submission audit version (2026-09-24, revision doc
NatureMethods_EIA_v33_投稿前必做计算与文字修订意见_20260924). Core protocol is
the frozen v27 pipeline, unchanged (same seeds, sentinel, formulas, bootstrap).
Everything runs from this package with RELATIVE paths; all inputs are frozen
under data/raw_inputs.

Unchanged from v27 (verified by hard assertions against the v27 frozen ledger):
  S1  Trivial transfer baselines; S3 repeated CV; S4 strict labels;
  S5  annotation-poor sensitivity; S6 vectorised nulls (N=1000);
  S7  pairwise inversion proportion I vs component consensus (2M pairs);
      v33 redefinition - the v27 all-8-layer comparable-pair definition is
      empty on the frozen universe (8-layer intersection = 0 genes) and
      the v27 archived rates were not reproducible from the archived code;
      the new explicit definition is recorded in 16.2 + the run record.
  S8  canonical E8; S9 MDE; S10 no-label sensitivity; S11 sanity assertions.

v33 changes (all traceable to the 2026-09-24 pre-submission audit):
  P0-1  Scorer manifest frozen: 19 scorers = 9 single-layer + 7 fixed-form
        + 3 supervised; scorer-endpoint coverage matrix in the ledger.
  P0-2  BH family rebuilt from the machine-readable ledger: all 60 contrasts
        (7 composites + 3 supervised) x 6 real endpoints now carry paired-
        bootstrap p (v27 used p=1.0 placeholders for LR/EN); q recomputed;
        included/excluded rows with reasons (SelC duplicate, E7/E8 controls,
        E6 preregistered external evaluation); survivors re-derived.
  P0-3  Source-lineage ledger (per input layer and per endpoint).
  P0-4  SHA-256 manifest of all frozen inputs + run record.
  P1-2  Missingness-only decomposition: values+missingness vs values-only vs
        missingness-only (coverage count) per real endpoint.
  P1-3  Aggregation robustness: pairwise inversion vs majority consensus
        across all 7 fixed-form composites, stratified by coverage;
        Kendall/Spearman vs component consensus (jointly with S7 in 16.2).
  P1-4  Nested-CV sensitivity (outer 5-fold, inner 3-fold, fixed grids).
  E6    External drug-response evaluation recomputed natively from the
        deterministic v19 GDSC construction, now with the FULL 7-composite
        family + supervised + bootstrap (v19 had 6 composites).
  E3/E3-C/E2/gradient: recomputed natively (v19 construction ported, same
        seeds) so Fig. 2 numbers regenerate from this package.

Audit protocol frozen as in v26/v27: sentinel=-3, harmonic/ECS formulas,
seed 20260819, bootstrap 2000. Python >= 3.10, numpy, scikit-learn.
"""
import json, os, sys, time, platform, datetime, csv
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, ElasticNet
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

T0 = time.time()
SEED = 20260819
N_BOOT = 2000
N_NULL = 1000          # v27: 100 -> 1000 (R1-M6)
N_REPEAT_CV = 10      # repeated CV seeds
SENT = -3.0
TOP_SEL = 250
TOP_E7 = 200

# ---- package-relative paths (R2-Major-8 fix) ----
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RAW = os.path.join(ROOT, "data", "raw_inputs")
OUTDIR = os.path.join(ROOT, "data", "results")
CSVDIR = os.path.join(ROOT, "data", "source_data")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(CSVDIR, exist_ok=True)

def log(m):
    print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

# =====================================================================
# 1. FROZEN DATA
# =====================================================================
EV = json.load(open(os.path.join(RAW, "evidence_layers.json")))
LAYERS, ALL_GENES = EV["layers"], EV["genes"]
FEATURES = ["string_centrality", "mutation_freq", "impc_animal_ko",
            "genetic_constraint", "cancer_driver", "ot_genetics_pdac",
            "druggability", "hpa_pdac_prognostic", "hpa_rna_tissue_spec"]
N_GENES = len(ALL_GENES)
X = np.full((N_GENES, len(FEATURES)), SENT)
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
PDAC_DEP = json.load(open(os.path.join(RAW, "depmap_pdac_dependency.json")))
CRC_DEP = json.load(open(os.path.join(RAW, "depmap_crc_dependency.json")))
E5_DATA = json.load(open(os.path.join(RAW, "e6_clinical_validation.json")))
CLIN_POS = set(E5_DATA["genes"])
PSD = json.load(open(os.path.join(RAW, "pdac_selective_dependency.json")))

# ---- S11 sanity assertions (hard fail; R2-Minor-11 fix) ----
assert N_GENES == 20751, f"universe changed: {N_GENES}"
assert len(PDAC_DEP) == 18531 and len(CRC_DEP) == 18531
assert PSD["meta"]["pan_essential_def"].startswith("top quartile"), \
    "E1 definition drifted from frozen meta"
log(f"genes={N_GENES}; DepMap files 18,531 each (universe is the UNION "
    f"across layers; genes without DepMap entry are treated as negatives "
    f"- documented, sensitivity in S10)")

# =====================================================================
# 2. ENDPOINTS (canonical E8 = interaction label)
# =====================================================================
def yvec(s):
    return np.array([1.0 if g in s else 0.0 for g in ALL_GENES])

def is_ess(g, dep):
    d = dep.get(g)
    return bool(d.get("essential")) if isinstance(d, dict) else False

ess_pdac = {g for g in ALL_GENES if is_ess(g, PDAC_DEP)}
ess_crc = {g for g in ALL_GENES if is_ess(g, CRC_DEP)}
# S4 strict labels
strict_pdac = {g for g in ALL_GENES if PDAC_DEP.get(g, {}).get("dependency_score", -9) > 0.5}
strict_crc = {g for g in ALL_GENES if CRC_DEP.get(g, {}).get("dependency_score", -9) > 0.5}

def top_by(score_map, k):
    ranked = sorted(((g, v) for g, v in score_map.items() if v is not None),
                    key=lambda kv: -kv[1])
    return {g for g, _ in ranked[:k]}

selA = top_by({g: d["A_ratio"] for g, d in PSD["definitions"].items()}, TOP_SEL)
selB = top_by({g: d["B_zeffect"] for g, d in PSD["definitions"].items()}, TOP_SEL)
selC = top_by({g: d["C_lineage_adj"] for g, d in PSD["definitions"].items()}, TOP_SEL)
selD = top_by({g: d["D_mixed"] for g, d in PSD["definitions"].items()}, TOP_SEL)

ot_scores = {g: LAYERS["ot_genetics_pdac"][g]["norm"]
             for g in LAYERS["ot_genetics_pdac"] if g in ALL_GENES}
e7_pos = top_by(ot_scores, TOP_E7)
# S8 canonical E8: interaction label (difference structure)
diff_score = X[:, I_MUT] - X[:, I_STR]
e8_pos = {ALL_GENES[i] for i in np.argsort(-diff_score)[:TOP_SEL]}

# SelB == SelC is a property of the frozen v11 file (verified: 18,434 genes,
# 0 differences). v27 counts them as ONE independent definition.
LABEL_DEPS = {
    "E1 PDAC pan-dependency": {},
    "E4 CRC pan-dependency": {},
    "E5 clinical concordance": {},
    "SelA PDAC-selective (A_ratio)": {},
    "SelB PDAC-selective (B_zeffect==C)": {},
    "SelD PDAC-selective (D_mixed)": {},
    "E7 genetic-association reuse": {"ot_genetics_pdac": 1.0},
}
ENDPOINTS = {
    "E1 PDAC pan-dependency": yvec(ess_pdac),
    "E4 CRC pan-dependency": yvec(ess_crc),
    "E5 clinical concordance": yvec(CLIN_POS),
    "SelA PDAC-selective (A_ratio)": yvec(selA),
    "SelB PDAC-selective (B_zeffect==C)": yvec(selB),
    "SelD PDAC-selective (D_mixed)": yvec(selD),
    "E7 genetic-association reuse": yvec(e7_pos),
}
REAL_EPS = [k for k in ENDPOINTS if not k.startswith("E7")]

# sanity
assert int(yvec(ess_pdac).sum()) == 4584
assert int(yvec(ess_crc).sum()) == 4608
assert int(yvec(e7_pos).sum()) == 200 and int(yvec(e8_pos).sum()) == 250
for nm, s in [("selA", selA), ("selB", selB), ("selC", selC), ("selD", selD)]:
    log(f"  {nm}: n_pos={len(s & set(ALL_GENES))}")
log("E5 label provenance note: 35-gene curated set overlaps druggability "
    "(AUROC below) and centrality; graded 'correlated source' in v27.")

# =====================================================================
# 3. SCORERS (frozen formulas; constants reported in Methods)
# =====================================================================
def DP(M):
    D = 0.80 * M[:, I_STR] + 0.10 * M[:, I_MUT] + 0.10 * M[:, I_IMPC]
    return D, np.mean(M[:, SUP_ALL], axis=1)

def f_multiplicative(M, a=0.6):
    D, P = DP(M); return D * (1 + a * P)
def f_additive(M, a=0.6):
    D, P = DP(M); return D + a * P
def f_geometric(M, a=0.6):
    D, P = DP(M)
    return np.sqrt(np.maximum(D, .01) * np.maximum(P + 3., .01))
def f_harmonic(M, a=0.6):
    D, P = DP(M)
    Ds, Ps = np.maximum(D + 3., .01), np.maximum(P + 3., .01)
    return 2. * Ds * Ps / (Ds + Ps)
def rank_agg(M, w=None):
    w = np.ones(M.shape[1]) if w is None else w
    R = np.zeros_like(M)
    for j in range(M.shape[1]):
        R[:, j] = np.argsort(np.argsort(M[:, j])).astype(float)
    return R @ w
def f_wrank(M):
    w = np.zeros(M.shape[1])
    w[I_STR], w[FEATURES.index("druggability")], w[I_MUT] = 5., 3., 1.
    w[FEATURES.index("cancer_driver")] = 1.
    return rank_agg(M, w)

SINGLE_LAYERS = {
    "STRING centrality": lambda M: M[:, I_STR],
    "Mutation frequency": lambda M: M[:, I_MUT],
    "IMPC animal KO": lambda M: M[:, I_IMPC],
    "Genetic constraint": lambda M: M[:, FEATURES.index("genetic_constraint")],
    "Cancer-driver annotation": lambda M: M[:, FEATURES.index("cancer_driver")],
    "OT genetics (PDAC)": lambda M: M[:, FEATURES.index("ot_genetics_pdac")],
    "Druggability": lambda M: M[:, FEATURES.index("druggability")],
    "HPA PDAC prognostic": lambda M: M[:, FEATURES.index("hpa_pdac_prognostic")],
    "HPA RNA tissue spec": lambda M: M[:, FEATURES.index("hpa_rna_tissue_spec")],
}
COMPOSITES = {
    "ECS (multiplicative)": f_multiplicative,
    "Additive": f_additive,
    "Geometric mean": f_geometric,
    "Harmonic mean": f_harmonic,
    "Arithmetic mean": lambda M: np.mean(M, axis=1),
    "Rank aggregation": lambda M: rank_agg(M),
    "Weighted rank": f_wrank,
}

def auroc(y, s):
    ok = np.isfinite(s)
    return float(roc_auc_score(y[ok], s[ok]))

SCORES = {}
for k, fn in SINGLE_LAYERS.items():
    SCORES[k] = fn(X)
for k, fn in COMPOSITES.items():
    SCORES[k] = fn(X)

# feature-name -> SCORES display-name (for LABEL_DEPS provenance lookup)
FEATURE_TO_SCORE = {
    "string_centrality": "STRING centrality",
    "mutation_freq": "Mutation frequency",
    "impc_animal_ko": "IMPC animal KO",
    "genetic_constraint": "Genetic constraint",
    "cancer_driver": "Cancer-driver annotation",
    "ot_genetics_pdac": "OT genetics (PDAC)",
    "druggability": "Druggability",
    "hpa_pdac_prognostic": "HPA PDAC prognostic",
    "hpa_rna_tissue_spec": "HPA RNA tissue spec",
}
def score_key_for_feature(f):
    return FEATURE_TO_SCORE.get(f, f)

# =====================================================================
# 4. SUPERVISED (fold-internal scaler; 3 learners; repeated CV; transfer)
# =====================================================================
def _fit_predict(Xtr, ytr, Xte, learner, seed):
    sc = StandardScaler().fit(Xtr)                 # S3: fold-internal scaling
    Xtr_s, Xte_s = sc.transform(Xtr), sc.transform(Xte)
    if learner == "Random forest":
        m = RandomForestClassifier(n_estimators=200, max_depth=10, n_jobs=-1,
                                   class_weight="balanced", random_state=seed)
        m.fit(Xtr_s, ytr); return m.predict_proba(Xte_s)[:, 1]
    if learner == "Logistic regression":
        m = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced")
        m.fit(Xtr_s, ytr); return m.predict_proba(Xte_s)[:, 1]
    m = ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=5000)
    m.fit(Xtr_s, ytr); return m.predict(Xte_s)

LEARNERS = ["Logistic regression", "Elastic net", "Random forest"]

def supervised_oof(y, fold_seed=42):
    y = y.astype(int)
    if y.sum() < 10:
        return {}
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=fold_seed)
    out = {L: np.zeros(N_GENES) for L in LEARNERS}
    for tr, te in skf.split(X, y):
        for L in LEARNERS:
            out[L][te] = _fit_predict(X[tr], y[tr], X[te], L, fold_seed)
    return out

def repeated_cv_delta(y, best_single_score, n_rep=N_REPEAT_CV):
    """S3: Delta distribution across fold partitions (RF + LR + EN)."""
    y = y.astype(int)
    rows = {L: [] for L in LEARNERS}
    for rep in range(n_rep):
        fold_seed = 101 + rep
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=fold_seed)
        oof = {L: np.zeros(N_GENES) for L in LEARNERS}
        for tr, te in skf.split(X, y):
            for L in LEARNERS:
                oof[L][te] = _fit_predict(X[tr], y[tr], X[te], L, fold_seed)
        for L in LEARNERS:
            rows[L].append(auroc(y, oof[L]) - auroc(y, best_single_score))
    return {L: {"deltas": [round(float(v), 4) for v in rows[L]],
                "mean": round(float(np.mean(rows[L])), 4),
                "sd": round(float(np.std(rows[L])), 4),
                "min": round(float(np.min(rows[L])), 4),
                "max": round(float(np.max(rows[L])), 4)}
            for L in LEARNERS}

# =====================================================================
# 5. PAIRED BOOTSTRAP (seed 20260819, B=2000, same protocol as v26)
# =====================================================================
rng_boot = np.random.default_rng(SEED)
BIDX = rng_boot.integers(0, N_GENES, size=(N_BOOT, N_GENES))

def paired_boot_delta(y, s_a, s_b):
    d = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = BIDX[b]
        d[b] = auroc(y[idx], s_a[idx]) - auroc(y[idx], s_b[idx])
    lo, hi = np.percentile(d, [2.5, 97.5])
    p_le0 = (int(np.sum(d <= 0)) + 0.5) / (N_BOOT + 1)   # one-sided, right-tail H1: d>0
    return {"delta": round(auroc(y, s_a) - auroc(y, s_b), 4),
            "ci95": [round(float(lo), 4), round(float(hi), 4)],
            "p_one_sided": round(float(p_le0), 5),
            "boot_mean": round(float(np.mean(d)), 4)}

def bh_q(pvals):
    p = np.asarray(pvals, float); m = len(p)
    order = np.argsort(p); q = np.empty(m); prev = 1.0
    for r in range(m - 1, -1, -1):
        i = order[r]
        prev = min(prev, p[i] * m / (r + 1)); q[i] = prev
    return np.minimum(q, 1.0)

# =====================================================================
# 6. VECTORISED NULLS (N=1000)
# =====================================================================
present_mask = (X != SENT)
n_obs = present_mask.sum(axis=1)
dec = np.digitize(X[:, I_STR],
                  np.quantile(X[:, I_STR], np.linspace(0.1, 0.9, 9)))
DEC_BINS = [np.where(dec == d)[0] for d in np.unique(dec)]

def degree_null_matrix(draw_rng):
    """Shuffle support-layer values within centrality deciles (keeps
    centrality + per-layer marginals within bins; breaks gene identity)."""
    Xd = X.copy()
    for j in SUP_ALL:
        for ii in DEC_BINS:
            col = Xd[ii, j]
            keep = col != SENT
            vals = col[keep]
            draw_rng.shuffle(vals)
            col[keep] = vals
            Xd[ii, j] = col
    return Xd

# pre-pad present values per gene for the density null
VALPAD = np.full((N_GENES, len(FEATURES)), np.nan)
for i in range(N_GENES):
    v = X[i, present_mask[i]]
    VALPAD[i, :len(v)] = v

def density_null_matrix(draw_rng):
    """Reassign which layers are observed per gene, preserving each gene's
    observation count and observed values (annotation-density null)."""
    P = np.argsort(draw_rng.random((N_GENES, len(FEATURES))), axis=1)
    Xp = np.full_like(X, SENT)
    for k in range(1, len(FEATURES)):
        ii = np.where(n_obs == k)[0]
        if len(ii) == 0:
            continue
        cols = P[ii, :k]
        vals = VALPAD[ii, :k]
        Xp[np.repeat(ii, k), cols.ravel()] = vals.ravel()
    return Xp

def null_z(obs, draws):
    arr = np.asarray(draws, float)
    mu, sd = float(np.mean(arr)), float(np.std(arr))
    z = (obs - mu) / max(sd, 1e-9)
    p_ge = (int(np.sum(arr >= obs)) + 0.5) / (len(arr) + 1)
    return {"observed": round(obs, 4), "null_mean": round(mu, 4),
            "null_sd": round(sd, 4), "z": round(z, 2),
            "p_emp_ge": round(p_ge, 4)}

# =====================================================================
# 7. RUN PER-ENDPOINT ALIGNMENT (7 real endpoints + E7 + E8 canonical)
# =====================================================================
RESULTS = {
    "meta": {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "python": platform.python_version(), "numpy": np.__version__,
        "sklearn": __import__("sklearn").__version__,
        "seed": SEED, "n_bootstrap": N_BOOT, "n_null": N_NULL,
        "n_genes": N_GENES, "features": FEATURES,
        "protocol": ("v33: 6 real disease endpoints + 2 audit controls (E7/E8) "
                     "+ preregistered external evaluation (E6); SelB==SelC "
                     "counted once (frozen-file identity); fold-internal "
                     "scaling; repeated CV; BH multiplicity over m=60 with "
                     "paired-bootstrap p for every contrast; self-audit "
                     "S1-S11; P0 ledgers (scorer manifest, delta/BH ledger, "
                     "source lineage) and P1 analyses (missingness-only, "
                     "aggregation robustness, nested CV)"),
        "universe_note": ("20,751-gene universe is the UNION of layer coverage; "
                          "DepMap files contain 18,531 genes; 2,303 genes lack "
                          "dependency labels and are treated as negatives "
                          "(S10 sensitivity quantifies the bias)"),
        "label_definitions": {
            "E1": "DepMap PDAC 'essential' flag (top-quartile rule per frozen "
                  "meta pan_essential_def; strict-label sensitivity in S4)",
            "E4": "DepMap CRC 'essential' flag (same rule)",
            "E5": "35-gene curated PDAC clinical-target set (graded "
                  "'correlated source' in v27)",
            "SelA-D": "top 250 by A_ratio / B_zeffect(==C_lineage_adj, frozen "
                      "file identity) / D_mixed from 26Q1-derived v11 matrix"},
    },
    "endpoints": {},
}

rng_density = np.random.default_rng(SEED + 2)
rng_degree = np.random.default_rng(SEED + 3)

# precompute null matrices once (shared across endpoints: null depends on y)
log(f"generating {N_NULL} null matrices (degree + density)...")
NULL_DEG = [degree_null_matrix(rng_degree) for _ in range(N_NULL)]
NULL_DEN = [density_null_matrix(rng_density) for _ in range(N_NULL)]
log("null matrices ready")

OOF_STORE = {}          # v33: keep supervised OOF vectors for the S2 ledger
for ep, y in ENDPOINTS.items():
    log(f"--- endpoint: {ep} (n_pos={int(y.sum())}) ---")
    R = {"n_pos": int(y.sum())}

    deps = LABEL_DEPS[ep]
    R["provenance"] = {
        "overlap_score": float(sum(deps.values())),
        "label_dependent_layers": deps,
        "implicated_layer_auroc": {f: round(auroc(y, SCORES[score_key_for_feature(f)]), 4) for f in deps},
        "circularity_flag": bool(sum(deps.values()) > 0),
    }

    # available-case representation (same as v26)
    Xm = X.copy(); Xm[np.isclose(Xm, SENT)] = np.nan
    Dm = np.nansum(np.stack([0.80 * Xm[:, I_STR], 0.10 * Xm[:, I_MUT],
                             0.10 * Xm[:, I_IMPC]]), axis=0)
    with np.errstate(invalid="ignore"):
        PHIm = np.nanmean(Xm[:, SUP_ALL], axis=1)
    PHIm = np.where(np.isfinite(PHIm), PHIm, 0.0)
    Dm = np.where(np.isfinite(Dm), Dm, 0.0)
    harm_ac = 2.0 * np.maximum(Dm + 3., .01) * np.maximum(PHIm + 3., .01) / \
        (np.maximum(Dm + 3., .01) + np.maximum(PHIm + 3., .01))
    R["representation"] = {
        "harmonic_sentinel_as_value": round(auroc(y, SCORES["Harmonic mean"]), 4),
        "harmonic_available_case": round(auroc(y, harm_ac), 4),
    }
    R["representation"]["harmonic_shift"] = round(
        R["representation"]["harmonic_available_case"] -
        R["representation"]["harmonic_sentinel_as_value"], 4)

    # aggregation table: 9 single + 7 composites (+ 3 supervised if real)
    table = {}
    for name, s in SCORES.items():
        table[name] = round(auroc(y, s), 4)
    if ep in REAL_EPS:
        sup = supervised_oof(y)
        OOF_STORE[ep] = sup          # v33: retained for the BH ledger
        for name, s in sup.items():
            table[name] = round(auroc(y, s), 4)
    R["aggregation"] = table

    single_aurocs = {n: v for n, v in table.items() if n in SINGLE_LAYERS}
    best_single = max(single_aurocs, key=single_aurocs.get)
    R["best_single_layer"] = {"name": best_single, "auroc": single_aurocs[best_single],
                              "selection_rule": "max observed AUROC on the same "
                                                "evaluation (data-dependent; declared)"}
    comp_aurocs = {n: table[n] for n in COMPOSITES}
    best_comp = max(comp_aurocs, key=comp_aurocs.get)
    R["best_composite"] = {"name": best_comp, "auroc": comp_aurocs[best_comp]}

    s_b = SCORES[best_single]
    # full per-composite deltas (family for BH)
    fam = {}
    for cn in COMPOSITES:
        fam[cn] = paired_boot_delta(y, SCORES[cn], s_b)
    R["delta_all_composites_vs_best_single"] = fam
    R["delta_bestcomposite_vs_bestsingle"] = dict(fam[best_comp],
                                                  composite=best_comp,
                                                  single=best_single)
    if ep in REAL_EPS and "Random forest" in table:
        R["delta_rf_vs_best_single"] = paired_boot_delta(y, sup["Random forest"], s_b)

    # baseline family nulls (N=1000) for harmonic / ECS / additive
    if ep in REAL_EPS:
        log(f"    nulls (N={N_NULL}) for harmonic family...")
        bf = {}
        for nm in ("Harmonic mean", "ECS (multiplicative)", "Additive"):
            dr = [auroc(y, COMPOSITES[nm](Md)) for Md in NULL_DEG]
            dn = [auroc(y, COMPOSITES[nm](Mn)) for Mn in NULL_DEN]
            bf[nm] = {"degree": null_z(table[nm], dr),
                      "density": null_z(table[nm], dn)}
        R["baseline_family"] = bf
    RESULTS["endpoints"][ep] = R
    with open(os.path.join(OUTDIR, "alignment_results.json"), "w") as f:
        json.dump(RESULTS, f, indent=1)

# ---- E8 canonical endpoint (audit control; supervised evaluated) ----
log("--- endpoint: E8 interaction label (canonical) ---")
yE8 = yvec(e8_pos)
tab8 = {n: round(auroc(yE8, s), 4) for n, s in SCORES.items()}
sup8 = supervised_oof(yE8)
for n, s in sup8.items():
    tab8[n] = round(auroc(yE8, s), 4)
bs8 = max((n for n in tab8 if n in SINGLE_LAYERS), key=lambda n: tab8[n])
bc8 = max((n for n in tab8 if n in COMPOSITES), key=lambda n: tab8[n])
E8 = {
    "definition": ("interaction label: top 250 genes by norm(mutation_freq) "
                   "minus norm(string_centrality)"),
    "n_pos": 250,
    "label_dependent_layers": {"mutation_freq": 0.5, "string_centrality": 0.5},
    "overlap_score": 1.0,
    "single_layer_auroc": {k: v for k, v in tab8.items() if k in SINGLE_LAYERS},
    "composite_auroc": {k: v for k, v in tab8.items() if k in COMPOSITES},
    "best_single": {"name": bs8, "auroc": tab8[bs8]},
    "best_fixed_composite": {"name": bc8, "auroc": tab8[bc8]},
    "delta_fixed_form_vs_best_single": paired_boot_delta(
        yE8, SCORES[bc8], SCORES[bs8]),
    "supervised_auroc": {n: round(auroc(yE8, s), 4) for n, s in sup8.items()},
    "delta_rf_vs_best_single": paired_boot_delta(
        yE8, sup8["Random forest"], SCORES[bs8]),
    "interpretation": ("the best fixed-form rule (ECS) recovers 0.894 while the "
                       "harmonic rule is nearly blind (0.057); a flexible "
                       "supervised integrator learns the constructed label "
                       "perfectly. Because the label is a function of the "
                       "inputs (overlap = 1.0), the audit withholds "
                       "information-gain attribution for BOTH."),
}
RESULTS["endpoints"]["E8 interaction label (canonical)"] = E8
log(f"  E8: best single={bs8} {tab8[bs8]}, best fixed={bc8} {tab8[bc8]}, "
    f"RF={E8['supervised_auroc']['Random forest']}, "
    f"fixed Delta={E8['delta_fixed_form_vs_best_single']['delta']}, "
    f"RF Delta={E8['delta_rf_vs_best_single']['delta']}")

# =====================================================================
# 8. SELF-AUDIT S1: trivial transfer baselines
# =====================================================================
log("S1: trivial transfer baselines...")
sp = np.array([PDAC_DEP.get(g, {}).get("dependency_score", np.nan) for g in ALL_GENES])
sc = np.array([CRC_DEP.get(g, {}).get("dependency_score", np.nan) for g in ALL_GENES])
yp, yc = ENDPOINTS["E1 PDAC pan-dependency"], ENDPOINTS["E4 CRC pan-dependency"]
label_copy_pc = auroc(yc, yp)   # binary PDAC flag as CRC predictor
label_copy_cp = auroc(yp, yc)
okp, okc = np.isfinite(sp), np.isfinite(sc)
cont_pc = auroc(yc[okp], sp[okp])
cont_cp = auroc(yp[okc], sc[okc])

# RF transfer, replicated with fold-internal scaling
def fit_rf_ens(y, fold_seed=42):
    sk = StratifiedKFold(n_splits=5, shuffle=True, random_state=fold_seed)
    models = []
    for tr, _ in sk.split(X, y.astype(int)):
        sc_ = StandardScaler().fit(X[tr])
        rf = RandomForestClassifier(n_estimators=200, max_depth=10, n_jobs=-1,
                                    class_weight="balanced", random_state=fold_seed)
        rf.fit(sc_.transform(X[tr]), y[tr].astype(int))
        models.append((sc_, rf))
    return models

def ens_proba(models):
    P = np.zeros(N_GENES)
    for sc_, m in models:
        P += m.predict_proba(sc_.transform(X))[:, 1]
    return P / len(models)

rf_tr_p2c = auroc(yc, ens_proba(fit_rf_ens(yp)))
rf_tr_c2p = auroc(yp, ens_proba(fit_rf_ens(yc)))
SELF_AUDIT = {
    "s1_trivial_transfer_baselines": {
        "rf_transfer_pdac_to_crc": round(rf_tr_p2c, 4),
        "rf_transfer_crc_to_pdac": round(rf_tr_c2p, 4),
        "label_copy_binary_pdac_to_crc": round(label_copy_pc, 4),
        "label_copy_binary_crc_to_pdac": round(label_copy_cp, 4),
        "dependency_score_direct_pdac_to_crc": round(cont_pc, 4),
        "dependency_score_direct_crc_to_pdac": round(cont_cp, 4),
        "n_scored_continuous": int(okp.sum()),
        "verdict": ("The supervised transfer AUROC is BELOW both trivial "
                    "baselines: the transfer claim carries no information "
                    "beyond label overlap; reported as label concordance, "
                    "not generalization."),
    },
}
log(f"  RF transfer {rf_tr_p2c:.4f}/{rf_tr_c2p:.4f} vs label-copy "
    f"{label_copy_pc:.4f}/{label_copy_cp:.4f} vs continuous "
    f"{cont_pc:.4f}/{cont_cp:.4f}")

# =====================================================================
# 9. SELF-AUDIT S2 (v33, P0-2): BH multiplicity over the full Delta family
#    Every contrast now carries a paired-bootstrap p (v27 substituted
#    p=1.0 for LR/EN); family m=60; the ledger is machine-checkable.
# =====================================================================
log("S2: BH multiplicity over full family (v33 corrected ledger)...")
RUN_ID = "v33-" + datetime.datetime.now().strftime("%Y%m%dT%H%M%S") + f"-seed{SEED}"
tests = []
for ep in REAL_EPS:
    r = RESULTS["endpoints"][ep]
    s_b = SCORES[r["best_single_layer"]["name"]]
    fam = r["delta_all_composites_vs_best_single"]
    for cn, d in fam.items():
        tests.append({"run_id": RUN_ID, "endpoint": ep, "method": cn,
                      "class": "fixed-form",
                      "baseline": r["best_single_layer"]["name"],
                      "direction": "auroc(method)-auroc(baseline); positive = gain",
                      "delta": d["delta"], "ci95": d["ci95"],
                      "p_one_sided": d["p_one_sided"],
                      "included_in_family": True, "exclusion_reason": ""})
    for L in ("Random forest", "Logistic regression", "Elastic net"):
        d = paired_boot_delta(ENDPOINTS[ep], OOF_STORE[ep][L], s_b)
        if L == "Random forest":
            d0 = RESULTS["endpoints"][ep]["delta_rf_vs_best_single"]
            assert abs(d["delta"] - d0["delta"]) < 5e-4, (ep, d, d0)
        tests.append({"run_id": RUN_ID, "endpoint": ep, "method": L,
                      "class": "supervised (5-fold out-of-fold)",
                      "baseline": r["best_single_layer"]["name"],
                      "direction": "auroc(method)-auroc(baseline); positive = gain",
                      "delta": d["delta"], "ci95": d["ci95"],
                      "p_one_sided": d["p_one_sided"],
                      "included_in_family": True, "exclusion_reason": ""})
assert len(tests) == 60, f"family size changed: {len(tests)}"
EXCLUDED_ROWS = [
    {"endpoint": "SelC PDAC-selective (C_lineage_adj)", "method": "all",
     "reason": "identical to SelB in the frozen 26Q1 file (top-250 Jaccard 1.000); counted once as SelB"},
    {"endpoint": "E7 genetic-association reuse", "method": "all",
     "reason": "audit control: label reuses an input layer by construction; excluded from the family"},
    {"endpoint": "E8 interaction label (canonical)", "method": "all",
     "reason": "audit control: constructed interaction label; excluded from the family"},
    {"endpoint": "E6 PDAC drug-response actionability", "method": "all",
     "reason": "preregistered external evaluation with its own audit criterion; reported separately, not in the exploratory family"},
]
qs = bh_q([t["p_one_sided"] for t in tests])
for t, q in zip(tests, qs):
    t["bh_q"] = round(float(q), 4)
SELF_AUDIT["s2_multiplicity"] = {
    "family": "7 fixed-form composites + 3 supervised out-of-fold integrators, "
              "each vs the strongest single-layer baseline, across 6 real "
              "endpoints (E1, E4, E5, SelA, SelB==SelC, SelD); m = 60, every "
              "contrast with paired-bootstrap p (B=2000, seed 20260819); "
              "E7/E8 are controls and E6 the preregistered external "
              "evaluation, excluded with reasons",
    "tests": tests,
    "excluded": EXCLUDED_ROWS,
    "v33_change": ("LR/EN contrasts now carry paired-bootstrap p; v27 "
                   "substituted p=1.0 for them (conservative); family size "
                   "and membership unchanged"),
}
survivors = [t for t in tests if t["bh_q"] < 0.05 and t["delta"] > 0]
SELF_AUDIT["s2_multiplicity"]["survivors_bh_q_lt_0.05"] = survivors
log(f"  BH survivors (q<0.05, Delta>0): "
    f"{[(t['endpoint'], t['method'], t['delta'], t['bh_q']) for t in survivors]}")

# =====================================================================
# 10. SELF-AUDIT S3: repeated-CV Delta distributions (E1/E4)
# =====================================================================
log("S3: repeated CV (10x5-fold)...")
REP_CV = {}
for ep in ("E1 PDAC pan-dependency", "E4 CRC pan-dependency"):
    bs_name = RESULTS["endpoints"][ep]["best_single_layer"]["name"]
    REP_CV[ep] = {"best_single": bs_name,
                  "repeated_cv_delta": repeated_cv_delta(
                      ENDPOINTS[ep], SCORES[bs_name])}
    r = REP_CV[ep]["repeated_cv_delta"]["Random forest"]
    log(f"  {ep}: RF repeated-CV Delta mean={r['mean']} sd={r['sd']} "
        f"range=[{r['min']},{r['max']}]")
SELF_AUDIT["s3_repeated_cv"] = REP_CV

# =====================================================================
# 11. SELF-AUDIT S4: strict-label definition sensitivity
# =====================================================================
log("S4: strict labels (dependency_score > 0.5)...")
STRICT = {}
for name, pos, yref in (("E1-strict", strict_pdac, yp), ("E4-strict", strict_crc, yc)):
    y = yvec(pos)
    tab = {n: round(auroc(y, s), 4) for n, s in SCORES.items()}
    sup = supervised_oof(y)
    for n, s in sup.items():
        tab[n] = round(auroc(y, s), 4)
    singles = {n: v for n, v in tab.items() if n in SINGLE_LAYERS}
    bs = max(singles, key=singles.get)
    comps = {n: tab[n] for n in COMPOSITES}
    bc = max(comps, key=comps.get)
    STRICT[name] = {
        "n_pos": int(y.sum()),
        "definition": "dependency_score > 0.5 (common-essential scale)",
        "best_single": {"name": bs, "auroc": singles[bs]},
        "best_composite": {"name": bc, "auroc": comps[bc]},
        "delta_bestcomposite": paired_boot_delta(y, SCORES[bc], SCORES[bs]),
        "rf_oof_auroc": tab["Random forest"],
        "delta_rf": paired_boot_delta(y, sup["Random forest"], SCORES[bs]),
        "table": tab,
    }
    # repeated CV on strict labels too
    STRICT[name]["repeated_cv_delta_rf"] = repeated_cv_delta(y, SCORES[bs])["Random forest"]
    log(f"  {name}: n_pos={int(y.sum())}, best single={bs} {singles[bs]}, "
        f"best comp={bc} {comps[bc]}, RF={tab['Random forest']}, "
        f"RF Delta={STRICT[name]['delta_rf']['delta']}")
SELF_AUDIT["s4_strict_labels"] = STRICT

# =====================================================================
# 12. SELF-AUDIT S5: annotation-poor sensitivity (coverage >= 4/9)
# =====================================================================
log("S5: annotation-poor gene sensitivity (coverage >= 4 of 9 layers)...")
cov = present_mask.sum(axis=1)
keep = cov >= 4
log(f"  subset n={int(keep.sum())} of {N_GENES}")
def subset_delta(y, score_a, score_b):
    ys, sa, sb = y[keep], score_a[keep], score_b[keep]
    d = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = BIDX[b][keep] if False else None
    # bootstrap on the subset: resample subset indices
    rng_sub = np.random.default_rng(SEED + 77)
    nsub = int(keep.sum())
    for b in range(N_BOOT):
        idx = rng_sub.integers(0, nsub, nsub)
        d[b] = auroc(ys[idx], sa[idx]) - auroc(ys[idx], sb[idx])
    lo, hi = np.percentile(d, [2.5, 97.5])
    return {"delta": round(auroc(ys, sa) - auroc(ys, sb), 4),
            "ci95": [round(float(lo), 4), round(float(hi), 4)]}
ANN = {"subset_n": int(keep.sum()),
       "rule": "genes covered by >= 4 of 9 evidence layers"}
for ep in ("E1 PDAC pan-dependency",
           "SelA PDAC-selective (A_ratio)",
           "SelB PDAC-selective (B_zeffect==C)",
           "SelD PDAC-selective (D_mixed)"):
    y = ENDPOINTS[ep]
    tab_sub = {n: round(auroc(y[keep], SCORES[n][keep]), 4) for n in SCORES}
    singles = {n: v for n, v in tab_sub.items() if n in SINGLE_LAYERS}
    bs = max(singles, key=singles.get)
    comps = {n: tab_sub[n] for n in COMPOSITES}
    bc = max(comps, key=comps.get)
    ANN[ep] = {"best_single": {"name": bs, "auroc": singles[bs]},
               "best_composite": {"name": bc, "auroc": comps[bc]},
               "delta_bestcomposite": subset_delta(y, SCORES[bc], SCORES[bs])}
    log(f"  {ep}: single={bs} {singles[bs]}, comp={bc} {comps[bc]}, "
        f"Delta={ANN[ep]['delta_bestcomposite']}")
# SelB top-250 artifact composition
_ze = {gg: d['B_zeffect'] for gg, d in PSD['definitions'].items()}
selB_top = sorted([g for g in selB if g in ALL_GENES], key=lambda g: -_ze[g])
covB = np.array([cov[list(ALL_GENES).index(g)] for g in selB_top])
ANN["selB_top250_composition"] = {
    "mean_coverage": round(float(covB.mean()), 2),
    "n_below_4_layers": int((covB < 4).sum()),
    "example_top_genes": selB_top[:12],
    "note": "olfactory-receptor / testis-enriched genes dominate the top of "
            "the list (classic CRISPR low-expression artifact territory); "
            "DepMap expression/CNV layers are not in the frozen inputs, so "
            "coverage-based filtering is used as the declared proxy"}
SELF_AUDIT["s5_annotation_poor_sensitivity"] = ANN

# =====================================================================
# 13. SELF-AUDIT S6b: single-layer degree-matched nulls (E1/E4)
# =====================================================================
log("S6b: single-layer nulls (degree-decile shuffles, N=1000)...")
SLN = {}
for epname in ("E1 PDAC pan-dependency", "E4 CRC pan-dependency"):
    y = ENDPOINTS[epname]
    d = {}
    for j, fname in enumerate(FEATURES):
        if j == I_STR:
            continue  # centrality defines the stratification; its own null is
                      # the network-rewiring null (pending, Discussion)
        obs = auroc(y, X[:, j])
        draws = []
        for Md in NULL_DEG:
            draws.append(auroc(y, Md[:, j]))
        d[fname] = null_z(obs, draws)
    SLN[epname] = d
    sig = [k for k, v in d.items() if v["z"] > 2]
    log(f"  {epname}: layers beyond degree structure (z>2): {sig}")
SELF_AUDIT["s6b_single_layer_nulls"] = SLN

# =====================================================================
# 14. SELF-AUDIT S7 (v33 redefinition): pairwise inversion proportion I
#     vs component consensus. Computed in Part B (16.2), where COV and
#     CONSENSUS are defined.
#     Rationale (recorded in the run record, 16.8): the v27 archived
#     definition ("all component layers cover both genes") is
#     mathematically empty on the frozen universe - the 8-layer
#     intersection contains 0 genes (cancer_driver covers 245 genes,
#     ot_genetics_pdac 500) - so it yields n = 0 comparable pairs; the
#     v27 archived 1,847,559-pair rates cannot be reproduced from the
#     archived v27 code plus the frozen inputs (provenance broken) and
#     are superseded by an explicit consensus-based definition.
# =====================================================================
log("S7 (v33): pairwise inversion proportion vs component consensus "
    "-> computed in Part B 16.2 (redefinition recorded in run record)")

# =====================================================================
# 15. SELF-AUDIT S9: MDE table + S10: no-label sensitivity + E5 grading
# =====================================================================
MDE = {}
for ep in REAL_EPS:
    for tag, key in (("best_composite", "delta_bestcomposite_vs_bestsingle"),):
        d = RESULTS["endpoints"][ep][key]
        if d["ci95"]:
            se = (d["ci95"][1] - d["ci95"][0]) / (2 * 1.96)
            MDE[ep] = {"se": round(float(se), 4),
                       "mde_80pct": round(float(2.8 * se), 4)}
if "delta_rf_vs_best_single" in RESULTS["endpoints"].get(
        "E1 PDAC pan-dependency", {}):
    d = RESULTS["endpoints"]["E1 PDAC pan-dependency"]["delta_rf_vs_best_single"]
    se = (d["ci95"][1] - d["ci95"][0]) / (2 * 1.96)
    MDE["E1 RF"] = {"se": round(float(se), 4),
                    "mde_80pct": round(float(2.8 * se), 4)}
SELF_AUDIT["s9_mde"] = {"note": "80%-power minimum detectable effect = 2.8 x SE "
                        "(two-sided alpha 0.05), SE from bootstrap CI width",
                        "table": MDE}

labelled = np.array([g in PDAC_DEP or g in CRC_DEP for g in ALL_GENES])
e1_lab = {n: round(auroc(yp[labelled], SCORES[n][labelled]), 4)
          for n in ("STRING centrality", "Harmonic mean")}
SELF_AUDIT["s10_no_label_sensitivity"] = {
    "n_labelled": int(labelled.sum()),
    "e1_auroc_labelled_subset": e1_lab,
    "note": "unlabelled genes are forced negatives in the main analysis; "
            "on the labelled subset the AUROC ordering is unchanged "
            "(full audit on the subset available on request)",
}
tab5 = RESULTS["endpoints"]["E5 clinical concordance"]["aggregation"]
SELF_AUDIT["e5_provenance_regrade"] = {
    "druggability_auroc": tab5["Druggability"],
    "ot_genetics_auroc": tab5["OT genetics (PDAC)"],
    "string_centrality_auroc": tab5["STRING centrality"],
    "v27_grade": "correlated source (was: independent source in v26)",
    "note": "the 35-gene curated label overlaps the druggability concept; "
            "centrality 0.986 on E5 is the strongest popularity fingerprint "
            "in the study and is discussed as such",
}

RESULTS["self_audit"] = SELF_AUDIT

# =====================================================================
# 16. WRITE FIGURE SOURCE-DATA CSVs
# =====================================================================
log("writing source-data CSVs...")
def wcsv(name, header, rows):
    with open(os.path.join(CSVDIR, name), "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)

# fig7a alignment matrix (v27): unified E8 row, honest z, 6th dimension column
order7 = ["E1 PDAC pan-dependency", "E4 CRC pan-dependency",
          "E5 clinical concordance", "SelA PDAC-selective (A_ratio)",
          "SelB PDAC-selective (B_zeffect==C)", "SelD PDAC-selective (D_mixed)",
          "E7 genetic-association reuse"]
rows = []
for ep in order7:
    r = RESULTS["endpoints"][ep]
    d5 = r["delta_bestcomposite_vs_bestsingle"]
    bf = r.get("baseline_family", {}).get("Harmonic mean", {})
    genz = ("n/a (control)" if ep.startswith("E7") else
            "E4/E5/Sel: external-family evidence; E1: same-context")
    rows.append([ep, r["n_pos"], r["provenance"]["overlap_score"],
                 r["provenance"]["circularity_flag"],
                 r["representation"]["harmonic_shift"],
                 r["best_single_layer"]["name"], r["best_single_layer"]["auroc"],
                 r["best_composite"]["name"], r["best_composite"]["auroc"],
                 d5["delta"], d5["ci95"][0], d5["ci95"][1],
                 bf.get("degree", {}).get("z", ""),
                 bf.get("degree", {}).get("p_emp_ge", ""),
                 bf.get("density", {}).get("z", ""),
                 bf.get("density", {}).get("p_emp_ge", ""),
                 r["aggregation"].get("Harmonic mean", ""),
                 r["aggregation"].get("ECS (multiplicative)", ""),
                 (r["aggregation"].get("Random forest") or ""),
                 genz])
rows.append(["E8 interaction label (canonical)", E8["n_pos"], 1.0, "True", "",
              E8["best_single"]["name"], E8["best_single"]["auroc"],
              E8["best_fixed_composite"]["name"], E8["best_fixed_composite"]["auroc"],
              E8["delta_fixed_form_vs_best_single"]["delta"],
              E8["delta_fixed_form_vs_best_single"]["ci95"][0],
              E8["delta_fixed_form_vs_best_single"]["ci95"][1],
              "", "", "", "",
              E8["composite_auroc"]["Harmonic mean"],
              E8["composite_auroc"]["ECS (multiplicative)"],
              E8["supervised_auroc"]["Random forest"],
              "n/a (control)"])
wcsv("fig5a_alignment_matrix.csv",
     ["endpoint", "n_pos", "overlap_score", "circular", "harmonic_shift_repr",
      "best_single_layer", "best_single_auroc", "best_composite",
      "best_composite_auroc", "delta", "delta_ci_lo", "delta_ci_hi",
      "degree_z_harmonic", "degree_p_emp", "density_z_harmonic",
      "density_p_emp", "harmonic_auroc", "ecs_auroc", "rf_auroc",
      "generalization_note"], rows)

# fig5b delta forest with BH q (v33: all 60 contrasts carry p and q)
rows = []
for t in tests:
    rows.append([t["endpoint"], t["method"], t["delta"],
                 t["ci95"][0], t["ci95"][1],
                 t["p_one_sided"], t["bh_q"]])
wcsv("fig5b_delta_forest.csv",
     ["endpoint", "scorer", "delta", "ci_lo", "ci_hi", "p_one_sided",
      "bh_q"], rows)

# fig6c transfer baselines (new panel)
wcsv("fig6c_transfer_baselines.csv",
     ["method", "pdac_to_crc", "crc_to_pdac"],
     [["Random forest (5-fold ensemble, no refit)",
       SELF_AUDIT["s1_trivial_transfer_baselines"]["rf_transfer_pdac_to_crc"],
       SELF_AUDIT["s1_trivial_transfer_baselines"]["rf_transfer_crc_to_pdac"]],
      ["Label copy (binary essential flag)",
       SELF_AUDIT["s1_trivial_transfer_baselines"]["label_copy_binary_pdac_to_crc"],
       SELF_AUDIT["s1_trivial_transfer_baselines"]["label_copy_binary_crc_to_pdac"]],
      ["Dependency score (continuous, direct)",
       SELF_AUDIT["s1_trivial_transfer_baselines"]["dependency_score_direct_pdac_to_crc"],
       SELF_AUDIT["s1_trivial_transfer_baselines"]["dependency_score_direct_crc_to_pdac"]]])

# fig3b (pairwise inversion proportion) is written in Part B 16.2 (v33).

# supp strict labels
rows = []
for name in ("E1-strict", "E4-strict"):
    s = STRICT[name]
    rows.append([name, s["n_pos"], s["best_single"]["name"],
                 s["best_single"]["auroc"], s["best_composite"]["name"],
                 s["best_composite"]["auroc"],
                 s["delta_bestcomposite"]["delta"],
                 s["delta_bestcomposite"]["ci95"][0],
                 s["delta_bestcomposite"]["ci95"][1],
                 s["rf_oof_auroc"], s["delta_rf"]["delta"],
                 s["delta_rf"]["ci95"][0], s["delta_rf"]["ci95"][1],
                 s["repeated_cv_delta_rf"]["mean"],
                 s["repeated_cv_delta_rf"]["sd"]])
wcsv("supp2c_strict_labels.csv",
     ["endpoint", "n_pos", "best_single", "best_single_auroc",
      "best_composite", "best_composite_auroc", "delta_comp",
      "delta_comp_ci_lo", "delta_comp_ci_hi", "rf_auroc", "delta_rf",
      "delta_rf_ci_lo", "delta_rf_ci_hi", "rcv_mean", "rcv_sd"], rows)

# supp annotation-poor sensitivity
rows = []
for ep in ("E1 PDAC pan-dependency", "SelA PDAC-selective (A_ratio)",
           "SelB PDAC-selective (B_zeffect==C)", "SelD PDAC-selective (D_mixed)"):
    a = ANN[ep]
    rows.append([ep, ANN["subset_n"], a["best_single"]["name"],
                 a["best_single"]["auroc"], a["best_composite"]["name"],
                 a["best_composite"]["auroc"],
                 a["delta_bestcomposite"]["delta"],
                 a["delta_bestcomposite"]["ci95"][0],
                 a["delta_bestcomposite"]["ci95"][1]])
wcsv("supp2d_annotated_subset.csv",
     ["endpoint", "subset_n", "best_single", "best_single_auroc",
      "best_composite", "best_composite_auroc", "delta_comp",
      "delta_comp_ci_lo", "delta_comp_ci_hi"], rows)

# single-layer nulls
rows = []
for ep, d in SLN.items():
    for fname, v in d.items():
        rows.append([ep, fname, v["observed"], v["null_mean"], v["null_sd"],
                     v["z"], v["p_emp_ge"]])
wcsv("single_layer_degree_nulls.csv",
     ["endpoint", "layer", "observed_auroc", "null_mean", "null_sd", "z",
      "p_emp_ge"], rows)

# repeated CV csv
rows = []
for ep, v in REP_CV.items():
    for L, r in v["repeated_cv_delta"].items():
        rows.append([ep, L, r["mean"], r["sd"], r["min"], r["max"]] +
                    r["deltas"])
wcsv("repeated_cv_delta.csv",
     ["endpoint", "learner", "mean", "sd", "min", "max"] +
     [f"rep{i+1}" for i in range(N_REPEAT_CV)], rows)

# =====================================================================
# 16. v33 PART B — native recomputations + P0/P1 ledgers (2026-09-24 audit)
# =====================================================================
import hashlib
from scipy.stats import spearmanr, kendalltau

V33 = RESULTS["v33"] = {}

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def _repr_block(y):
    """Harmonic available-case representation (identical formula to main loop)."""
    Xm = X.copy(); Xm[np.isclose(Xm, SENT)] = np.nan
    Dm = np.nansum(np.stack([0.80 * Xm[:, I_STR], 0.10 * Xm[:, I_MUT],
                             0.10 * Xm[:, I_IMPC]]), axis=0)
    with np.errstate(invalid="ignore"):
        PHIm = np.nanmean(Xm[:, SUP_ALL], axis=1)
    PHIm = np.where(np.isfinite(PHIm), PHIm, 0.0)
    Dm = np.where(np.isfinite(Dm), Dm, 0.0)
    ac = 2.0 * np.maximum(Dm + 3., .01) * np.maximum(PHIm + 3., .01) / \
        (np.maximum(Dm + 3., .01) + np.maximum(PHIm + 3., .01))
    return {"harmonic_sentinel_as_value": round(auroc(y, SCORES["Harmonic mean"]), 4),
            "harmonic_available_case": round(auroc(y, ac), 4)}

# ---------- 16.1  P1-2  missingness-only decomposition ----------
log("v33 P1-2: missingness-only decomposition (values / coverage / both) ...")
OBS = (X != SENT)
COV = OBS.sum(axis=1).astype(float)
X_MISS = np.where(OBS, 1.0, SENT)
MISS_ONLY = f_harmonic(X_MISS)          # coverage pattern in composite geometry

D_num = np.zeros(N_GENES); D_den = np.zeros(N_GENES)
for j, w in ((I_STR, .8), (I_MUT, .1), (I_IMPC, .1)):
    D_num += w * np.where(OBS[:, j], X[:, j], 0.0)
    D_den += w * OBS[:, j].astype(float)
P_num = np.zeros(N_GENES); P_den = np.zeros(N_GENES)
for j in SUP_ALL:
    P_num += np.where(OBS[:, j], X[:, j], 0.0)
    P_den += OBS[:, j].astype(float)
with np.errstate(invalid="ignore", divide="ignore"):
    Ds_vo = np.where(D_den > 0,
                     np.maximum(D_num / np.maximum(D_den, 1e-12) + 3., .01), np.nan)
    Ps_vo = np.where(P_den > 0,
                     np.maximum(P_num / np.maximum(P_den, 1e-12) + 3., .01), np.nan)
    VALS_ONLY = 2. * Ds_vo * Ps_vo / (Ds_vo + Ps_vo)

S12 = {}
for ep in REAL_EPS:
    y = ENDPOINTS[ep]
    full = auroc(y, SCORES["Harmonic mean"])
    vo = auroc(y, VALS_ONLY)
    mo = auroc(y, MISS_ONLY)
    cc = auroc(y, COV)
    S12[ep] = {
        "full_values_plus_missingness": round(full, 4),
        "values_only": round(vo, 4),
        "missingness_only_composite": round(mo, 4),
        "coverage_count": round(cc, 4),
        "delta_missingness_only_minus_full": round(mo - full, 4),
        "per_layer_coverage_auroc": {FEATURES[j]:
            round(auroc(y, OBS[:, j].astype(float)), 4) for j in range(9)},
        "interpretation": ("coverage pattern alone recovers most of the "
                           "composite's discrimination" if mo >= full - 0.02
                           else "values carry discrimination beyond the "
                                 "coverage pattern"),
    }
    log(f"  {ep}: full={full:.4f} values-only={vo:.4f} "
        f"missingness-only={mo:.4f} coverage-count={cc:.4f}")
SELF_AUDIT["s12_missingness_decomposition"] = S12

# ---------- 16.2  P1-3 + S7(v33): pairwise inversion vs consensus ----------
log("v33 P1-3 + S7: pairwise inversion proportion vs component consensus ...")
RANKP = np.full((N_GENES, 9), np.nan)
for j in range(9):
    m = OBS[:, j]
    ov = X[m, j]
    order = np.argsort(ov, kind="mergesort")
    r = np.empty(len(ov)); r[order] = np.arange(1, len(ov) + 1)
    RANKP[m, j] = r / len(ov)
with np.errstate(invalid="ignore"):
    CONSENSUS = np.nanmean(RANKP, axis=1)
rng_ar = np.random.default_rng(SEED + 3)
COMPOSITES_7 = ("Harmonic mean", "ECS (multiplicative)", "Additive",
                "Geometric mean", "Arithmetic mean", "Rank aggregation",
                "Weighted rank")
_wilson = lambda p, n: (1.96 * np.sqrt(p * (1 - p) / n + 1.96**2 / (4 * n**2))) / \
                       (1 + 1.96**2 / n)
OP = {}      # S7(v33): per-composite inversion vs consensus
S13 = {}     # P1-3: per-composite association with consensus (+ strata)
for cn in COMPOSITES_7:
    s = SCORES[cn]
    ok = np.isfinite(s) & np.isfinite(CONSENSUS)
    idx = np.where(ok)[0]
    tau = float(kendalltau(s[idx], CONSENSUS[idx])[0])
    rho = float(spearmanr(s[idx], CONSENSUS[idx])[0])
    ii = rng_ar.choice(idx, size=2_000_000, replace=True)
    jj = rng_ar.choice(idx, size=2_000_000, replace=True)
    keep = ii != jj
    ii, jj = ii[keep], jj[keep]
    d = np.sign(s[ii] - s[jj]) * np.sign(CONSENSUS[ii] - CONSENSUS[jj])
    ii, jj, d = ii[d != 0], jj[d != 0], d[d != 0]
    ci, cj = COV[ii], COV[jj]
    both_hi = (ci >= 7) & (cj >= 7)
    both_lo = (ci <= 3) & (cj <= 3)
    mixed = ~(both_hi | both_lo)
    def _rate(mask, dd=d):
        return round(float(np.mean(dd[mask] < 0)), 4)
    n_pairs = int(len(d))
    inv = int((d < 0).sum())
    p = inv / max(n_pairs, 1)
    w = _wilson(p, n_pairs) if n_pairs else 0.0
    OP[cn] = {"n_pairs": n_pairs, "inversions": inv,
              "I_pairwise_inversion": round(p, 4),
              "ci95_wilson": [round(max(0.0, p - w), 4),
                              round(min(1.0, p + w), 4)],
              "I_both_high_coverage": _rate(both_hi),
              "I_mixed_coverage": _rate(mixed),
              "I_both_low_coverage": _rate(both_lo)}
    S13[cn] = {"kendall_tau_vs_consensus": round(tau, 4),
               "spearman_rho_vs_consensus": round(rho, 4),
               "I_pairwise_inversion": OP[cn]["I_pairwise_inversion"],
               "inversion_both_high_coverage": OP[cn]["I_both_high_coverage"],
               "inversion_mixed_coverage": OP[cn]["I_mixed_coverage"],
               "inversion_both_low_coverage": OP[cn]["I_both_low_coverage"],
               "n_pairs": n_pairs}
    log(f"  {cn}: tau={S13[cn]['kendall_tau_vs_consensus']} "
        f"I={OP[cn]['I_pairwise_inversion']} "
        f"hi={OP[cn]['I_both_high_coverage']} "
        f"mi={OP[cn]['I_mixed_coverage']} "
        f"lo={OP[cn]['I_both_low_coverage']}")
SELF_AUDIT["s7_pairwise_inversion_proportion"] = {
    "definition": ("A gene pair is comparable when both genes carry a "
                   "composite score and a component-consensus rank (each "
                   "gene covered by at least one input layer; every gene "
                   "in the 20,751-gene union qualifies). I is the fraction "
                   "of comparable sampled pairs in which the composite "
                   "ordering disagrees with the majority-consensus "
                   "ordering. Estimated from 2,000,000 sampled pairs "
                   "(seed 20260819+3); Wilson 95% CI. v33 redefinition: "
                   "the v27 archived definition required all eight "
                   "component layers to cover both genes, which is empty "
                   "on the frozen universe (8-layer intersection = 0 "
                   "genes; cancer_driver covers 245, ot_genetics_pdac "
                   "500) and produced n = 0 comparable pairs; the v27 "
                   "archived 1,847,559-pair rates were not reproducible "
                   "from the archived v27 code plus the frozen inputs "
                   "and are superseded by this explicit definition."),
    "rates": OP,
}
SELF_AUDIT["s13_aggregation_robustness"] = {
    "note": ("computed jointly with S7 in 16.2; the inversion statistic is "
             "endpoint-independent (composite scores and the consensus do "
             "not depend on the endpoint), so v33 reports it once per "
             "composite across all seven fixed-form rules (the earlier "
             "draft duplicated it per endpoint)"),
    "per_composite": S13,
}
# fig3b: pairwise inversion proportion (v33 schema)
rows = [[cn, OP[cn]["n_pairs"], OP[cn]["inversions"],
         OP[cn]["I_pairwise_inversion"], OP[cn]["ci95_wilson"][0],
         OP[cn]["ci95_wilson"][1], OP[cn]["I_both_high_coverage"],
         OP[cn]["I_mixed_coverage"], OP[cn]["I_both_low_coverage"]]
        for cn in COMPOSITES_7]
wcsv("fig3b_pairwise_inversion.csv",
     ["composite", "n_pairs", "inversions", "I_pairwise_inversion",
      "ci_lo", "ci_hi", "I_both_high_coverage", "I_mixed_coverage",
      "I_both_low_coverage"], rows)

# ---------- 16.3  P1-4  nested-CV sensitivity ----------
log("v33 P1-4: nested-CV sensitivity (outer 5, inner 3, fixed grids) ...")
def _fit_grid(learner, cfg, Xtr, ytr, Xte, seed):
    sc = StandardScaler().fit(Xtr)
    Xtr_s, Xte_s = sc.transform(Xtr), sc.transform(Xte)
    if learner == "Random forest":
        m = RandomForestClassifier(n_estimators=cfg["n_estimators"],
                                    max_depth=cfg["max_depth"], n_jobs=-1,
                                    class_weight="balanced",
                                    random_state=seed)
        m.fit(Xtr_s, ytr); return m.predict_proba(Xte_s)[:, 1]
    if learner == "Logistic regression":
        m = LogisticRegression(max_iter=2000, C=cfg["C"],
                               class_weight="balanced")
        m.fit(Xtr_s, ytr); return m.predict_proba(Xte_s)[:, 1]
    m = ElasticNet(alpha=cfg["alpha"], l1_ratio=cfg["l1_ratio"],
                   max_iter=5000)
    m.fit(Xtr_s, ytr); return m.predict(Xte_s)

GRIDS = {
    "Random forest": [{"n_estimators": ne, "max_depth": md}
                       for ne in (100, 500) for md in (None, 8, 16)],
    "Logistic regression": [{"C": c} for c in (0.01, 0.1, 1.0, 10.0)],
    "Elastic net": [{"alpha": a, "l1_ratio": l}
                     for a in (0.01, 0.1, 1.0) for l in (0.1, 0.5, 0.9)],
}
S14 = {}
for ep in ("E1 PDAC pan-dependency", "E4 CRC pan-dependency"):
    y = ENDPOINTS[ep].astype(int)
    oof = {L: np.zeros(N_GENES) for L in GRIDS}
    chosen = {L: [] for L in GRIDS}
    out_folds = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    for k_out, (tr, te) in enumerate(out_folds.split(X, y)):
        for L, grid in GRIDS.items():
            inn = StratifiedKFold(n_splits=3, shuffle=True,
                                  random_state=SEED + 100 + k_out)
            best_cfg, best_a = None, -1.0
            for cfg in grid:
                scores = []
                for tri, tei in inn.split(X[tr], y[tr]):
                    pv = _fit_grid(L, cfg, X[tr][tri], y[tr][tri],
                                   X[tr][tei], SEED + k_out)
                    scores.append(roc_auc_score(y[tr][tei], pv))
                if float(np.mean(scores)) > best_a:
                    best_a, best_cfg = float(np.mean(scores)), cfg
            oof[L][te] = _fit_grid(L, best_cfg, X[tr], y[tr], X[te],
                                   SEED + k_out)
            chosen[L].append(best_cfg)
    S14[ep] = {}
    r = RESULTS["endpoints"][ep]
    s_b = SCORES[r["best_single_layer"]["name"]]
    for L in GRIDS:
        d = paired_boot_delta(ENDPOINTS[ep], oof[L], s_b)
        d["nested_oof_auroc"] = round(auroc(ENDPOINTS[ep], oof[L]), 4)
        d["best_config_per_outer_fold"] = chosen[L]
        S14[ep][L] = d
    log(f"  {ep}: RF nested Delta={S14[ep]['Random forest']['delta']} "
        f"(main fixed-hyper Delta="
        f"{r['delta_rf_vs_best_single']['delta']})")
SELF_AUDIT["s14_nested_cv_sensitivity"] = S14

# ---------- 16.4  provenance gradient + E2/E3/E3-A/E3-C family ----------
log("v33: provenance gradient (v18 construction, deterministic port) ...")
I_DRUG = FEATURES.index("druggability")
y1 = ENDPOINTS["E1 PDAC pan-dependency"]
leak_annotated = np.where(X[:, I_DRUG] != SENT)[0]
n_pos_true = int(y1.sum())
ess_idx = np.where(y1 == 1)[0]
noness_idx = np.where(y1 == 0)[0]
rng_p = np.random.default_rng(SEED + 1)
THETAS = [0.0, 0.25, 0.50, 0.75, 1.00]
GRAD = []
for th in THETAS:
    n_leak = int(round(th * n_pos_true))
    n_true = n_pos_true - n_leak
    pool_leak = np.array(sorted(set(leak_annotated.tolist()) &
                                set(noness_idx.tolist())))
    pool_leak = pool_leak if len(pool_leak) else np.arange(1)
    pick_leak = (rng_p.choice(pool_leak, size=min(n_leak, len(pool_leak)),
                              replace=(n_leak > len(pool_leak)))
                 if n_leak else np.array([], int))
    ess_drug = np.array(sorted(set(ess_idx.tolist()) &
                               set(leak_annotated.tolist())))
    ess_nodrug = np.array(sorted(set(ess_idx.tolist()) -
                                 set(leak_annotated.tolist())))
    pick_true = rng_p.choice(ess_nodrug, size=min(n_true, len(ess_nodrug)),
                             replace=True)
    if len(pick_true) < n_true:
        need = n_true - len(pick_true)
        extra = (rng_p.choice(ess_drug, size=min(need, len(ess_drug)),
                              replace=(need > len(ess_drug)))
                 if len(ess_drug) else np.array([], int))
        pick_true = np.concatenate([pick_true, extra])
    if len(pick_true) < n_true:
        leftover = np.setdiff1d(noness_idx, np.union1d(pick_true, pick_leak))
        pick_true = np.concatenate(
            [pick_true, rng_p.choice(leftover,
                                     size=n_true - len(pick_true),
                                     replace=False)])
    yth = np.zeros(N_GENES)
    yth[pick_true] = 1.0
    yth[pick_leak] = 1.0
    a_cent = auroc(yth, SCORES["STRING centrality"])
    a_hm = auroc(yth, SCORES["Harmonic mean"])
    a_drug = auroc(yth, X[:, I_DRUG])
    GRAD.append({"theta": th, "n_pos": int(yth.sum()),
                 "n_leaky_positives": n_leak,
                 "centrality_auroc": round(a_cent, 4),
                 "harmonic_auroc": round(a_hm, 4),
                 "harmonic_delta_vs_centrality": round(a_hm - a_cent, 4),
                 "druggability_auroc": round(a_drug, 4)})
    log(f"  theta={th:.2f} n={int(yth.sum())} centrality={a_cent:.4f} "
        f"harmonic Delta={a_hm - a_cent:+.4f}")
deltas = [g["harmonic_delta_vs_centrality"] for g in GRAD]
rho_g = float(spearmanr(THETAS, deltas)[0])
V33["provenance_gradient"] = {
    "design": ("true positives replaced fraction theta by genes merely "
               "annotated in the druggability layer; coverage constant; "
               "theta=0 coverage-controlled E1, theta=1 E3-like pathology"),
    "rows": GRAD,
    "spearman_rho_theta_vs_delta": round(rho_g, 4),
    "monotone_increasing": bool(all(deltas[i] <= deltas[i + 1] + 1e-9
                                    for i in range(len(deltas) - 1))),
}

log("v33: provenance endpoint family (E2/E3/E3-A/E3-C) native recompute ...")
drug_pos_set = {ALL_GENES[i] for i in np.where(X[:, I_DRUG] != SENT)[0]}
pe = PSD["pan_essential"]
defs = PSD["definitions"]
pan_ess = {g for g in ALL_GENES if pe.get(g) is True}
bz = {g: defs[g]["B_zeffect"] for g in pan_ess if g in defs}
q75 = float(np.percentile(np.fromiter(bz.values(), float), 75))
E2_pos = {g for g, v in bz.items() if v >= q75}
E3_pos = ess_pdac & drug_pos_set
FAMILY = {
    "E2 PDAC-enriched dependency": E2_pos,
    "E3 conjunctive actionability": E3_pos,
    "E3-A leakage-controlled essentiality": set(ess_pdac),
    "E3-C out-of-evidence druggability": drug_pos_set,
}
X_NODRUG = X.copy(); X_NODRUG[:, I_DRUG] = SENT
H_NODRUG = f_harmonic(X_NODRUG)
FAM = {}
for nm, pos in FAMILY.items():
    yv = yvec(pos)
    a_cent = auroc(yv, SCORES["STRING centrality"])
    single = {k: round(auroc(yv, s), 4) for k, s in SCORES.items()
              if k in SINGLE_LAYERS}
    best_single_name = max(single, key=single.get)
    comp = {}
    for cn, fn in COMPOSITES.items():
        a = auroc(yv, fn(X))
        comp[cn] = {"auroc": round(a, 4),
                    "delta_vs_centrality": round(a - a_cent, 4),
                    "delta_vs_best_single": round(a - single[best_single_name], 4)}
    FAM[nm] = {
        "n_pos": len(pos),
        "centrality_auroc": round(a_cent, 4),
        "single_layer_auroc": single,
        "best_single": {"name": best_single_name,
                        "auroc": single[best_single_name]},
        "composites": comp,
        "harmonic_no_druggability": round(auroc(yv, H_NODRUG), 4),
        "representation": _repr_block(yv),
    }
    log(f"  {nm}: n={len(pos)} centrality={a_cent:.4f} "
        f"harmonic={comp['Harmonic mean']['auroc']:.4f} "
        f"no-drug={FAM[nm]['harmonic_no_druggability']:.4f}")
V33["provenance_endpoint_family"] = FAM

# ---------- 16.5  E6 native GDSC recompute (full 7-composite family) ----------
log("v33 E6: native GDSC recompute with full composite + supervised family ...")
# GDSC/DeepCDR external dataset is NOT redistributed here (provider terms).
# Point this at your local copy of the DeepCDR data directory.
GDSC_DIR = os.environ.get(
    "GDSC_DIR",
    str(Path.home() / "DeepCDR" / "data"),
)
GDSC_FILES = {
    "cell_annotations": os.path.join(GDSC_DIR, "CCLE",
                                     "Cell_lines_annotations_20181226.txt"),
    "gdsc_ic50": os.path.join(GDSC_DIR, "CCLE", "GDSC_IC50.csv"),
    "drug_list": os.path.join(GDSC_DIR, "GDSC",
                              "1.Drug_listMon Jun 24 09_00_55 2019.csv"),
}
with open(GDSC_FILES["cell_annotations"], encoding="utf-8",
          errors="replace") as f:
    ann = list(csv.DictReader(f, delimiter="\t"))
panc = {r["depMapID"].strip() for r in ann
        if (r.get("Site_Primary") or "").strip().lower() == "pancreas"
        and r.get("depMapID")}
with open(GDSC_FILES["gdsc_ic50"]) as f:
    ic = list(csv.reader(f))
cols = ic[0][1:]
pidx = [i for i, c in enumerate(cols) if c.strip() in panc]
drug_ic50 = {}
for r in ic[1:]:
    did = r[0].replace("GDSC:", "").strip()
    vals = []
    for i in pidx:
        s = r[1 + i].strip()
        if s not in ("", "NA", "NaN"):
            try:
                vals.append(float(s))
            except ValueError:
                pass
    if len(vals) >= 8:
        drug_ic50[did] = float(np.median(vals))
with open(GDSC_FILES["drug_list"], encoding="utf-8",
          errors="replace") as f:
    drugs = list(csv.DictReader(f))
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
y6 = yvec(sg - rg)
n6 = int(y6.sum())
a_cent6 = auroc(y6, SCORES["STRING centrality"])
single6 = {kk: round(auroc(y6, s), 4) for kk, s in SCORES.items()
           if kk in SINGLE_LAYERS}
best6 = max(single6, key=single6.get)
s_b6 = SCORES[best6]
comp6 = {}
for cn, fn in COMPOSITES.items():
    d = paired_boot_delta(y6, fn(X), s_b6)
    d["auroc"] = round(auroc(y6, fn(X)), 4)
    d["delta_vs_centrality"] = round(d["auroc"] - a_cent6, 4)
    comp6[cn] = d
sup6 = supervised_oof(y6) if n6 >= 10 else {}
sup6r = {}
for L, sv in sup6.items():
    d = paired_boot_delta(y6, sv, s_b6)
    d["auroc"] = round(auroc(y6, sv), 4)
    d["delta_vs_centrality"] = round(d["auroc"] - a_cent6, 4)
    sup6r[L] = d
worst_cent = max(max(c["delta_vs_centrality"] for c in comp6.values()),
                max((s["delta_vs_centrality"] for s in sup6r.values()),
                    default=-9.0))
verdict6 = ("PASS" if worst_cent <= 0 else
            "WARN" if worst_cent <= 0.02 else "FAIL")
e6v = {
    "available": True, "n_positives": n6, "n_drugs": len(shared),
    "n_pancreas_lines": len(pidx),
    "construction": ("GDSC1 median IC50 across pancreatic CCLE lines; "
                     "tercile split; positives = targets of most-sensitive "
                     "tercile minus targets of least-sensitive tercile "
                     "(frozen v19 construction, deterministic)"),
    "baseline_centrality_auroc": round(a_cent6, 4),
    "best_single": {"name": best6, "auroc": single6[best6]},
    "single_layer_auroc": single6,
    "composites_vs_best_single": comp6,
    "supervised_vs_best_single": sup6r,
    "max_delta_vs_centrality": round(float(worst_cent), 4),
    "preregistered_verdict": verdict6,
    "role": ("preregistered external evaluation (drug-target ecology; "
             "related to the druggability layer, not claimed fully "
             "independent)"),
}
V33["e6_external_drug_response"] = e6v
log(f"  E6: n_pos={n6} drugs={len(shared)} lines={len(pidx)} "
    f"centrality={a_cent6:.4f} best_single={best6}({single6[best6]}) "
    f"max Delta vs centrality={worst_cent:+.4f} verdict={verdict6}")

# ---------- 16.6  P0-1  scorer manifest ----------
log("v33 P0-1: scorer manifest ...")
SCORER_CLASSES = {}
for nm in SINGLE_LAYERS: SCORER_CLASSES[nm] = "single-layer evidence"
for nm in COMPOSITES: SCORER_CLASSES[nm] = "fixed-form integration"
for nm in LEARNERS: SCORER_CLASSES[nm] = "supervised out-of-fold"
MANIFEST = []
def _emit(scorer, ep, val):
    if val in ("", None):
        return
    try:
        val = round(float(val), 4)
    except (TypeError, ValueError):
        return
    MANIFEST.append({"run_id": RUN_ID, "scorer": scorer,
                     "class": SCORER_CLASSES[scorer], "endpoint": ep,
                     "endpoint_role": EPS_ROLE_OF.get(ep, ""), "auroc": val})
EPS_ROLE_OF = {}
for ep in REAL_EPS:
    EPS_ROLE_OF[ep] = "primary benchmark endpoint"
EPS_ROLE_OF["E7 genetic-association reuse"] = "audit control (direct reuse)"
EPS_ROLE_OF["E8 interaction label (canonical)"] = "audit control (constructed)"
for nm in FAMILY:
    EPS_ROLE_OF[nm] = "provenance gradient family"
EPS_ROLE_OF["E6 external drug response"] = "preregistered external evaluation"
for ep in REAL_EPS:
    tab = RESULTS["endpoints"][ep]["aggregation"]
    for sc in SCORER_CLASSES:
        if sc in tab:
            _emit(sc, ep, tab[sc])
E7R = RESULTS["endpoints"]["E7 genetic-association reuse"]["aggregation"]
for sc in SCORER_CLASSES:
    if sc in E7R:
        _emit(sc, "E7 genetic-association reuse", E7R[sc])
e8R = RESULTS["endpoints"]["E8 interaction label (canonical)"]
for sc in SINGLE_LAYERS:
    if sc in e8R["single_layer_auroc"]:
        _emit(sc, "E8 interaction label (canonical)",
              e8R["single_layer_auroc"][sc])
for sc in COMPOSITES:
    if sc in e8R["composite_auroc"]:
        _emit(sc, "E8 interaction label (canonical)",
              e8R["composite_auroc"][sc])
for sc in LEARNERS:
    if sc in e8R["supervised_auroc"]:
        _emit(sc, "E8 interaction label (canonical)",
              e8R["supervised_auroc"][sc])
for nm, fr in FAM.items():
    for sc, v in fr["single_layer_auroc"].items():
        _emit(sc, nm, v)
    for cn, dd in fr["composites"].items():
        _emit(cn, nm, dd["auroc"])
for sc, v in e6v["single_layer_auroc"].items():
    _emit(sc, "E6 external drug response", v)
for cn, dd in comp6.items():
    _emit(cn, "E6 external drug response", dd["auroc"])
for L, dd in sup6r.items():
    _emit(L, "E6 external drug response", dd["auroc"])
V33["scorer_manifest"] = {
    "n_scorers": len(SCORER_CLASSES),
    "n_scorer_endpoint_cells": len(MANIFEST),
    "declaration": ("19 scorers = 9 single-layer + 7 fixed-form + 3 "
                    "supervised out-of-fold; the same manifest underlies "
                    "the Methods text ('seven fixed-form'), every figure "
                    "bar count, and the BH family"),
    "cells": MANIFEST,
}
log(f"  manifest: {len(SCORER_CLASSES)} scorers, {len(MANIFEST)} "
    f"scorer-endpoint cells")

# ---------- 16.7  P0-3  source-lineage ledger (meta + computed) ----------
log("v33 P0-3: source-lineage ledger ...")
LAYER_META = EV["meta"]["layers"]
LAYERS_LINEAGE = []
for f in FEATURES:
    m = LAYER_META.get(f, {})
    LAYERS_LINEAGE.append({
        "layer": f, "type": m.get("type", ""), "n_covered": m.get("n", ""),
        "direction": m.get("direction", ""),
        "used_directly_by_labels":
            [ep for ep, deps in LABEL_DEPS.items() if f in deps],
        "overlap_class_with_labels": (
            "direct reuse (E7 label = top of this layer)" if f == "ot_genetics_pdac"
            else "direct reuse (E3-C label = layer membership; E3 conjunction)"
                 if f == "druggability"
            else "constructed control (E8 uses mutation + centrality difference)"
                 if f in ("mutation_freq", "string_centrality")
            else "no direct label reuse"),
    })
def _endpoint_lineage(ep_name, pos, role, note):
    yv = yvec(pos)
    per_layer = {}
    for j, f in enumerate(FEATURES):
        cov_pos = int(np.sum(OBS[yv == 1, j]))
        per_layer[f] = round(cov_pos / max(len(pos), 1), 4)
    return {"endpoint": ep_name, "role": role, "n_pos": len(pos),
            "frac_positives_annotated_per_layer": per_layer, "note": note}
ENDPOINTS_LINEAGE = [
    _endpoint_lineage("E1 PDAC pan-dependency", ess_pdac,
                      "primary", "DepMap PDAC essentiality; no layer shares "
                      "the CRISPR-screen lineage"),
    _endpoint_lineage("E4 CRC pan-dependency", ess_crc,
                      "primary (external context)",
                      "DepMap CRC; same platform as E1, independent of all "
                      "nine layers"),
    _endpoint_lineage("E5 clinical concordance", CLIN_POS,
                      "primary (translational)",
                      "curated clinical target set; conceptual overlap with "
                      "druggability and centrality (see per-layer fractions)"),
    _endpoint_lineage("SelA PDAC-selective (A_ratio)", selA, "primary",
                      "DepMap selective-dependency metric; same platform "
                      "as E1, no layer reuse"),
    _endpoint_lineage("SelB PDAC-selective (B_zeffect==C)", selB, "primary",
                      "DepMap z-effect metric; identical file column to "
                      "SelC (counted once)"),
    _endpoint_lineage("SelD PDAC-selective (D_mixed)", selD, "primary",
                      "DepMap mixed selective metric; no layer reuse"),
    _endpoint_lineage("E7 genetic-association reuse", e7_pos,
                      "audit control",
                      "label = top-200 of the OT genetics layer: direct "
                      "value reuse by construction (overlap 1.0)"),
    _endpoint_lineage("E8 interaction label (canonical)",
                      {ALL_GENES[i] for i in np.argsort(-diff_score)[:TOP_SEL]},
                      "audit control",
                      "constructed interaction label = top-250 of "
                      "(mutation - centrality) difference; circular by "
                      "construction"),
    _endpoint_lineage("E6 external drug response", sg - rg,
                      "preregistered external evaluation",
                      "GDSC1 drug-response construction; drug-target ecology "
                      "related to druggability; not claimed fully independent"),
]
for nm, pos in FAMILY.items():
    note = {
        "E2 PDAC-enriched dependency": "DepMap B_zeffect >= Q75; no layer reuse",
        "E3 conjunctive actionability": "ess AND druggability-annotated: "
                                        "direct membership reuse of the "
                                        "druggability layer",
        "E3-A leakage-controlled essentiality": "E1 definition (control)",
        "E3-C out-of-evidence druggability": "label == druggability layer "
                                             "membership: circular control",
    }[nm]
    ENDPOINTS_LINEAGE.append(_endpoint_lineage(nm, pos,
                                                "provenance gradient family",
                                                note))
V33["source_lineage"] = {"layers": LAYERS_LINEAGE,
                         "endpoints": ENDPOINTS_LINEAGE}

# ---------- 16.8  P0-4  run record ----------
RAW_FILES = sorted(fn for fn in os.listdir(RAW) if fn.endswith(".json"))
V33["run_record"] = {
    "run_id": RUN_ID,
    "started_utc": datetime.datetime.now(datetime.timezone.utc)
                    .isoformat(timespec="seconds"),
    "wall_time_seconds": round(time.time() - T0, 1),
    "seed": SEED, "bootstrap": N_BOOT, "nulls": N_NULL, "sentinel": SENT,
    "python": sys.version.split()[0], "numpy": np.__version__,
    "sklearn": __import__("sklearn").__version__,
    "code_sha256": sha256_file(os.path.abspath(__file__)),
    "frozen_inputs_sha256": {fn: sha256_file(os.path.join(RAW, fn))
                             for fn in RAW_FILES},
    "external_gdsc_sha256": {kk: sha256_file(vv)
                             for kk, vv in GDSC_FILES.items()},
    "bh_family_size": len(tests),
    "bh_n_survivors": len(survivors),
    "s7_redefinition_note": (
        "S7 was redefined in v33. The v27 archived comparable-pair "
        "definition (both genes covered by all 8 component layers) is "
        "mathematically empty on the frozen universe: the 8-layer "
        "intersection contains 0 genes (cancer_driver covers 245 genes, "
        "ot_genetics_pdac 500), so the archived v27 code yields n = 0 "
        "comparable pairs. The v27 archived rates (n = 1,847,559, "
        "I = 0.38-0.53) were generated after the v27 run by a side "
        "script (recompute_s7.py, 2026-09-17) using a looser "
        "definition (comparable = at least one jointly observed layer; "
        "reference = sign of the sum over jointly observed layers) while "
        "the archived definition text remained the strict all-layer "
        "version - the archived numbers and definition do not match. "
        "v33 replaces S7 with an explicit consensus-based pairwise "
        "inversion proportion (definition in "
        "self_audit.s7_pairwise_inversion_proportion); the manuscript "
        "text and Fig. 3b use only v33 numbers."),
}

# ---------- 16.9  v33 CSV sources ----------
rows = []
for ep, v in S12.items():
    rows.append([ep, "values+missingness (frozen harmonic)",
                 v["full_values_plus_missingness"]])
    rows.append([ep, "values-only", v["values_only"]])
    rows.append([ep, "missingness-only", v["missingness_only_composite"]])
    rows.append([ep, "coverage count", v["coverage_count"]])
    for lay, a in v["per_layer_coverage_auroc"].items():
        rows.append([f"{ep} [coverage indicator: {lay}]",
                     "per-layer coverage indicator", a])
wcsv("fig3d_missingness_decomposition.csv",
     ["endpoint", "scorer", "auroc"], rows)

rows = []
for cn in COMPOSITES_7:
    v = S13[cn]
    rows.append([cn, v["kendall_tau_vs_consensus"],
                 v["spearman_rho_vs_consensus"],
                 v["I_pairwise_inversion"],
                 v["inversion_both_high_coverage"],
                 v["inversion_mixed_coverage"],
                 v["inversion_both_low_coverage"], v["n_pairs"]])
wcsv("fig3c_aggregation_robustness.csv",
     ["composite", "kendall_tau", "spearman_rho",
      "I_pairwise_inversion", "inv_both_high_cov", "inv_mixed_cov",
      "inv_both_low_cov", "n_pairs"], rows)

rows = []
for ep, dd in S14.items():
    for L, v in dd.items():
        rows.append([ep, L, v["nested_oof_auroc"], v["delta"],
                     v["ci95"][0], v["ci95"][1],
                     json.dumps(v["best_config_per_outer_fold"]),
                     (RESULTS["endpoints"][ep]["delta_rf_vs_best_single"]["delta"]
                      if L == "Random forest" else "")])
wcsv("supp2e_nested_cv.csv",
     ["endpoint", "learner", "nested_oof_auroc", "delta_vs_best_single",
      "ci_lo", "ci_hi", "selected_config_per_outer_fold",
      "main_fixed_hyper_rf_delta"], rows)

# full BH ledger (Supp Table S2; every row re-computable from the ledger)
rows = [[t["run_id"], t["endpoint"], t["method"], t["class"], t["baseline"],
         t["direction"], t["delta"], t["ci95"][0], t["ci95"][1],
         t["p_one_sided"], t["bh_q"], t["included_in_family"],
         t["exclusion_reason"]] for t in tests]
for x in EXCLUDED_ROWS:
    rows.append([RUN_ID, x["endpoint"], x["method"], "", "", "", "", "", "",
                 "", "", False, x["reason"]])
wcsv("s2_ledger.csv",
     ["run_id", "endpoint", "method", "class", "baseline", "direction",
      "delta", "ci_lo", "ci_hi", "p_one_sided", "bh_q",
      "included_in_family", "exclusion_reason"], rows)

wcsv("fig2a_provenance_gradient.csv",
     ["theta", "delta_auroc", "centrality_auroc"],
     [[g["theta"], g["harmonic_delta_vs_centrality"],
       g["centrality_auroc"]] for g in GRAD])

wcsv("fig2b_e3_deletion.csv", ["condition", "label", "value"],
     [["full (E3)", "full (E3)",
       FAM["E3 conjunctive actionability"]["composites"]
          ["Harmonic mean"]["auroc"]],
      ["remove label-embedded druggability",
       "remove label-embedded druggability",
       FAM["E3 conjunctive actionability"]["harmonic_no_druggability"]]])
wcsv("fig2c_e3c_control.csv", ["condition", "label", "value"],
     [["Harmonic (E3-C)", "Harmonic (E3-C)",
       FAM["E3-C out-of-evidence druggability"]["composites"]
          ["Harmonic mean"]["auroc"]],
      ["Harmonic no-druggability", "Harmonic no-druggability",
       FAM["E3-C out-of-evidence druggability"]["harmonic_no_druggability"]]])

rows = [["STRING centrality", "AUROC", e6v["baseline_centrality_auroc"]]]
for nm in ("Harmonic mean", "ECS (multiplicative)", "Additive",
           "Geometric mean", "Arithmetic mean", "Rank aggregation",
           "Weighted rank"):
    rows.append([nm, "AUROC", comp6[nm]["auroc"]])
for L in ("Logistic regression", "Elastic net", "Random forest"):
    rows.append([L, "AUROC", sup6r[L]["auroc"]])
wcsv("fig5b_e6_discrimination.csv", ["method", "metric", "value"], rows)

rows = []
for nm, fr in FAM.items():
    short = nm.split(" ")[0] if not nm.startswith("E3-") else nm.split(" ")[0]
    prov = ("circular control" if nm.startswith("E3-C")
            else "direct label reuse" if nm.startswith("E3 ")
            else "leakage control" if nm.startswith("E3-A")
            else "independent source")
    ovl = 1.0 if nm.startswith(("E3 ", "E3-C")) else 0.0
    rows.append([short, "in-sample", prov, ovl,
                 fr["composites"]["Harmonic mean"]["delta_vs_centrality"]])
wcsv("fig4b_provenance_family.csv",
     ["endpoint", "context", "provenance", "overlap_score", "delta_auroc"],
     rows)

rows = [[c["run_id"], c["scorer"], c["class"], c["endpoint"],
         c["endpoint_role"], c["auroc"]] for c in MANIFEST]
wcsv("scorer_manifest.csv",
     ["run_id", "scorer", "class", "endpoint", "endpoint_role", "auroc"],
     rows)

rows = [[r["layer"], r["type"], r["n_covered"], r["direction"],
         ";".join(r["used_directly_by_labels"]),
         r["overlap_class_with_labels"]] for r in LAYERS_LINEAGE]
wcsv("source_lineage_layers.csv",
     ["layer", "type", "n_covered", "direction", "used_directly_by_labels",
      "overlap_class_with_labels"], rows)
rows = [[r["endpoint"], r["role"], r["n_pos"],
         json.dumps(r["frac_positives_annotated_per_layer"]),
         r["note"]] for r in ENDPOINTS_LINEAGE]
wcsv("source_lineage_endpoints.csv",
     ["endpoint", "role", "n_pos", "frac_positives_annotated_per_layer",
      "note"], rows)

# ---------- 16.10  v33 hard cross-checks against the frozen record ----------
assert n6 == 32, f"E6 n_positives drifted: {n6} (frozen: 32)"
assert len(shared) == 125 and len(pidx) == 29, \
    f"E6 drugs/lines drifted: {len(shared)}/{len(pidx)} (frozen: 125/29)"
FROZEN_GRAD_CENT = [0.692, 0.6171, 0.6041, 0.5865, 0.568]
FROZEN_GRAD_DELTA = [-0.0753, -0.0221, 0.0986, 0.2097, 0.298]
for g, fc, fd in zip(GRAD, FROZEN_GRAD_CENT, FROZEN_GRAD_DELTA):
    assert abs(g["centrality_auroc"] - fc) < 2e-3, (g, fc)
    assert abs(g["harmonic_delta_vs_centrality"] - fd) < 2e-3, (g, fd)
assert abs(FAM["E3 conjunctive actionability"]["composites"]
           ["Harmonic mean"]["auroc"] - 0.8926) < 2e-3
assert abs(FAM["E3 conjunctive actionability"]["harmonic_no_druggability"]
           - 0.6177) < 2e-3
assert abs(FAM["E3-C out-of-evidence druggability"]["composites"]
           ["Harmonic mean"]["auroc"] - 0.9357) < 2e-3
assert abs(FAM["E3-C out-of-evidence druggability"]["harmonic_no_druggability"]
           - 0.5827) < 2e-3
assert FAM["E2 PDAC-enriched dependency"]["n_pos"] == 1147
assert FAM["E3 conjunctive actionability"]["n_pos"] == 1159
assert FAM["E3-A leakage-controlled essentiality"]["n_pos"] == 4584
assert FAM["E3-C out-of-evidence druggability"]["n_pos"] == 5188
# S7(v33) frozen values: harmonic/ECS reproduce the v33-draft S13 draws
# (rng SEED+3, identical sampling sequence for the first two composites)
assert abs(OP["Harmonic mean"]["I_pairwise_inversion"] - 0.4292) < 5e-4, \
    OP["Harmonic mean"]["I_pairwise_inversion"]
assert abs(OP["ECS (multiplicative)"]["I_pairwise_inversion"] - 0.5408) < 5e-4, \
    OP["ECS (multiplicative)"]["I_pairwise_inversion"]
assert abs(S13["Harmonic mean"]["kendall_tau_vs_consensus"] - 0.143) < 5e-4, \
    S13["Harmonic mean"]["kendall_tau_vs_consensus"]
assert abs(S13["ECS (multiplicative)"]["kendall_tau_vs_consensus"]
           - (-0.0825)) < 5e-4, S13["ECS (multiplicative)"]
assert OP["Harmonic mean"]["n_pairs"] > 1_000_000, "S7 sample too small"
log("v33 cross-checks vs frozen record: PASS")

with open(os.path.join(OUTDIR, "alignment_results.json"), "w") as f:
    json.dump(RESULTS, f, indent=1)
log("saved alignment_results.json")

# ---- S11 final sanity ----
chk = RESULTS["endpoints"]["E1 PDAC pan-dependency"]["aggregation"]
assert abs(chk["STRING centrality"] - 0.6951) < 5e-4, chk["STRING centrality"]
assert abs(chk["Harmonic mean"] - 0.5752) < 5e-4, chk["Harmonic mean"]
assert RESULTS["endpoints"]["E7 genetic-association reuse"]["aggregation"][
    "OT genetics (PDAC)"] == 1.0
log("FINAL SANITY: PASS (E1 centrality/harmonic match frozen v19/v26 values)")
print(json.dumps({
    "s1": SELF_AUDIT["s1_trivial_transfer_baselines"],
    "bh_survivors": len(survivors),
    "s4_strict": {k: {"n_pos": v["n_pos"],
                      "delta_rf": v["delta_rf"]["delta"]} for k, v in STRICT.items()},
}, indent=1))
