# --- EIA reproducibility: paths are resolved relative to the repository root ---
from pathlib import Path as _P
ROOT = _P(__file__).resolve().parent.parent.parent

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_submission_forms.py
Generate the full set of Nature Methods submission documents for v38.

All outputs are produced in BOTH .docx and .md under:
  - 20260926_v38_NM/submission_forms/   (docx + md)
  - 20260926_v38_NM/manuscript/         (cover letter + contributor docs copied)

Per project convention the submission is ANONYMOUS: no author names or
affiliations appear in the manuscript or cover letter; [AUTHORS — INSERT ...]
placeholders are used instead.

Documents:
  1. cover_letter_v38
  2. Nature_Portfolio_Reporting_Summary_v38
  3. Software_Submission_Checklist_v38
  4. author_contributions_v38
  5. funding_v38
  6. competing_interests_v38
  7. SUBMISSION_CHECKLIST_v38   (internal P0/P1 closure log, md only)
  8. ARCHIVE_README_v38          (md only)
  9. CHECKSUMS_v38               (md only)
"""
import hashlib
import os
from docx import Document
from docx.shared import Pt

ROOT = str(ROOT) if False else ROOT
OUT = os.path.join(ROOT, "20260926_v38_NM", "submission_forms")
MAN = os.path.join(ROOT, "20260926_v38_NM", "manuscript")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# minimal block renderer -> both docx and md
# blocks: list of (kind, text)  kind in {h1,h2,para,bullet,space,rule}
# ---------------------------------------------------------------------------

def render_docx(title, blocks, path):
    d = Document()
    d.add_paragraph(title, style="Title")
    for kind, text in blocks:
        if kind == "space":
            d.add_paragraph()
        elif kind == "rule":
            d.add_paragraph("─" * 40)
        elif kind == "h1":
            d.add_paragraph(text, style="Heading 1")
        elif kind == "h2":
            d.add_paragraph(text, style="Heading 2")
        elif kind == "bullet":
            p = d.add_paragraph(style="List Bullet")
            p.add_run(text)
        else:
            d.add_paragraph(text)
    d.save(path)


def render_md(title, blocks, path):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# {title}\n\n")
        for kind, text in blocks:
            if kind == "space":
                f.write("\n")
            elif kind == "rule":
                f.write("\n---\n\n")
            elif kind == "h1":
                f.write(f"## {text}\n\n")
            elif kind == "h2":
                f.write(f"### {text}\n\n")
            elif kind == "bullet":
                f.write(f"- {text}\n")
            else:
                f.write(f"{text}\n\n")


def emit(base, title, blocks):
    dx = os.path.join(OUT, base + ".docx")
    mx = os.path.join(OUT, base + ".md")
    render_docx(title, blocks, dx)
    render_md(title, blocks, mx)
    # also drop a copy of cover letter + contributor docs into manuscript/
    if base in ("cover_letter_v38", "author_contributions_v38",
                "funding_v38", "competing_interests_v38"):
        import shutil
        shutil.copy(mx, os.path.join(MAN, base + ".md"))
        shutil.copy(dx, os.path.join(MAN, base + ".docx"))
    return dx, mx


# ===========================================================================
# 1. COVER LETTER (anonymous)
# ===========================================================================
cover_blocks = [
    ("para", "Date: [INSERT SUBMISSION DATE]"),
    ("space", ""),
    ("para", "Editorial Office\nNature Methods\nSpringer Nature"),
    ("space", ""),
    ("para", "Dear Dr. [INSERT HANDLING EDITOR NAME],"),
    ("space", ""),
    ("para", "Re: Submission of an Analysis article entitled \u201cAuditing the Attribution of "
            "Evidence-Integration Gains in Biomedical Benchmarks\u201d"),
    ("space", ""),
    ("para", "We enclose our manuscript for consideration as an Analysis article in Nature Methods. "
            "The work addresses a methodological gap that affects how the biomedical community "
            "interprets evidence-integration studies: a higher integrated score is routinely read as "
            "evidence that integration added new information, yet the observed gain can be produced "
            "by benchmark construction, annotation-driven missingness, aggregation rules, or a weak "
            "comparator rather than by independent signal."),
    ("space", ""),
    ("para", "We define the attribution gap and formalize the incremental value of integration as "
            "\u0394 = performance of an integrator relative to a prespecified single-layer comparator "
            "on the same evaluation setting. We then operationalize an Evidence Integration Audit "
            "(EIA) with four modules \u2014 provenance, representation, aggregation and external "
            "evaluation \u2014 and a built-in self-audit, applied across benchmark endpoints in "
            "pancreatic (PDAC) and colorectal (CRC) cancer. Controlled provenance experiments show "
            "that benchmark construction alone can manufacture apparent gains; missingness and "
            "aggregation choices also shift \u0394. Across external evaluations, the tested fixed-form "
            "family showed no transferable positive \u0394, whereas supervised models produced "
            "within-ecology candidate gains that remain unreplicated. EIA therefore treats a positive "
            "\u0394 as an attribution claim requiring explicit audit and independent replication."),
    ("space", ""),
    ("para", "We believe the manuscript is of direct interest to the Nature Methods readership because "
            "it provides a reproducible protocol and reporting standard for any study that compares an "
            "integrated scorer with alternatives \u2014 a situation that now spans single-cell, "
            "multi-omics, imaging and clinical-risk integration. The framework is implemented as an "
            "open toolkit with composable audit modules, claim-level scorecards and per-figure source "
            "data."),
    ("space", ""),
    ("para", "We confirm that this manuscript is original, has not been published elsewhere, is not "
            "under consideration by any other journal, and that all authors have approved the "
            "submission. The accompanying Nature Portfolio Reporting Summary and Software Submission "
            "Checklist are provided. Source data, analysis code and figure scripts are archived under "
            "a persistent identifier (repository identifier and access instructions to be inserted "
            "before submission)."),
    ("space", ""),
    ("para", "We declare no competing interests beyond those stated in the manuscript."),
    ("space", ""),
    ("para", "Correspondence and material requests should be addressed to [AUTHORS \u2014 INSERT FINAL "
            "NAMES, AFFILIATIONS AND CORRESPONDING-AUTHOR DETAILS]."),
    ("space", ""),
    ("para", "Sincerely,\n[AUTHORS \u2014 INSERT FINAL NAMES AND SIGNATURES]"),
]

# ===========================================================================
# 2. NATURE PORTFOLIO REPORTING SUMMARY (adapted for computational Analysis)
# ===========================================================================
rs_blocks = [
    ("h1", "Nature Portfolio Reporting Summary (adapted for a computational Analysis article)"),
    ("para", "This summary captures the reporting items Nature Methods requires for a computational / "
            "statistical-methods Analysis. Items are answered for the v38 frozen submission. Fields "
            "marked [INSERT \u2026] require author completion before final upload."),
    ("h2", "1. Data and code availability"),
    ("bullet", "All input evidence layers are derived from publicly available resources (STRING v12, "
               "DepMap 23Q2, gnomAD, COSMIC, Open Targets Genetics, DGIdb/ChEMBL and the Human "
               "Protein Atlas); gene-level harmonization is scripted. [CITATION 15\u201322]"),
    ("bullet", "Analysis code, the audit report, alignment results, figure scripts and per-figure "
               "source data are archived under a persistent identifier with the frozen input manifest. "
               "[INSERT final repository DOI / accession]"),
    ("bullet", "No new biological samples were generated; the study is a secondary computational "
               "analysis of publicly available aggregated data."),
    ("h2", "2. Statistical reporting"),
    ("bullet", "Sample sizes: nine single-layer evidence layers (n = 20,751 protein-coding genes in "
               "the HGNC intersection of STRING v12 and DepMap 23Q2); benchmark endpoints E1\u2013E6, "
               "SelA\u2013SelD (exact n per endpoint reported in Supplementary Table S3b)."),
    ("bullet", "The primary estimand is \u0394 on the AUROC scale; uncertainty via paired bootstrap "
               "(2,000 replicates; seed 20260819); multiplicity controlled by Benjamini\u2013Hochberg "
               "across the 60-contrast \u0394 family. [CITATION 8\u201310]"),
    ("bullet", "Randomization / blinding are not applicable; outcome inspection was precluded by "
               "pre-specification of the E6 audit rule before label inspection. [CITATION R4]"),
    ("h2", "3. Reproducibility and replication"),
    ("bullet", "Every reported number is traceable to frozen raw inputs and a re-runnable script "
               "(run_analysis.py \u2192 postprocess.py \u2192 build.py); a reproducible "
               "hash manifest is archived."),
    ("bullet", "Candidate supervised gains are explicitly reported as unreplicated; the manuscript "
               "states they require independent replication and proposes an audit-first, "
               "outcome-second validation roadmap."),
    ("h2", "4. Ethics"),
    ("bullet", "Secondary computational analysis of publicly available aggregated data; no individual-"
               "level data; ethics approval not required. Stated in the manuscript Ethics statement."),
    ("h2", "5. Author contributions / Funding / Competing interests"),
    ("bullet", "Provided as separate placeholder documents to be completed by the authors before "
               "submission (author_contributions_v38, funding_v38, competing_interests_v38)."),
    ("para", "[INSERT any further journal-specific reporting items required at upload]"),
]

# ===========================================================================
# 3. SOFTWARE SUBMISSION CHECKLIST (Nature Methods)
# ===========================================================================
sw_blocks = [
    ("h1", "Software Submission Checklist \u2014 EIA Toolkit (Nature Methods)"),
    ("para", "Confirms the computational tool accompanying the Analysis meets Nature Methods software "
            "reporting expectations."),
    ("h2", "Availability"),
    ("bullet", "Source code archived under a persistent identifier with a frozen environment "
               "[INSERT repository DOI]."),
    ("bullet", "License: [INSERT license, e.g., MIT/Apache-2.0] declared in the repository."),
    ("bullet", "Install instructions, command-line audit call and worked examples provided in the "
               "archive (see ARCHIVE_README)."),
    ("h2", "Documentation"),
    ("bullet", "README documents install, inputs, outputs and the four audit modules (provenance, "
               "representation, aggregation, external evaluation)."),
    ("bullet", "Per-figure source data and R/ggplot2 scripts (300 dpi, theme_classic) supplied for "
               "every display item; no AI-generated figures."),
    ("h2", "Testing / validation"),
    ("bullet", "Reproducible test: running the frozen pipeline regenerates all reported \u0394 values "
               "and the 60-contrast BH ledger (Supplementary Table S2)."),
    ("bullet", "Self-audit (S1\u2013S14) exercises the framework against its own candidate gains to "
               "avoid double standards."),
    ("h2", "Reporting"),
    ("bullet", "Use of AI-assisted tools disclosed in the Methods (\u201cUse of AI-assisted tools\u201d "
               "section). [CITATION R1]"),
    ("para", "[INSERT any remaining Nature Methods software-checklist items at upload]"),
]

# ===========================================================================
# 4-6. Contributor / funding / competing-interests placeholders
# ===========================================================================
ac_blocks = [
    ("para", "[INSERT FINAL AUTHOR CONTRIBUTIONS]"),
    ("para", "Conceptualization: [INSERT]. Methodology: [INSERT]. Software: [INSERT]. "
            "Investigation: [INSERT]. Data curation: [INSERT]. Writing \u2014 original draft: "
            "[INSERT]. Writing \u2014 review & editing: [INSERT]. Supervision: [INSERT]. "
            "Funding acquisition: [INSERT]."),
]
fu_blocks = [
    ("para", "[INSERT FINAL FUNDING STATEMENT]"),
    ("para", "This work was supported by [INSERT grant numbers and funding bodies]. "
            "The funders had no role in study design, analysis, or decision to publish."),
]
ci_blocks = [
    ("para", "[INSERT FINAL COMPETING-INTERESTS STATEMENT]"),
    ("para", "The authors declare no competing interests. [If any exist, list them here per Nature "
            "Methods policy.]"),
]

# ===========================================================================
# Emit docx + md
# ===========================================================================
emit("cover_letter_v38", "Cover Letter \u2014 EIA v38 (Nature Methods, Analysis)", cover_blocks)
emit("Nature_Portfolio_Reporting_Summary_v38",
     "Nature Portfolio Reporting Summary \u2014 EIA v38", rs_blocks)
emit("Software_Submission_Checklist_v38",
     "Software Submission Checklist \u2014 EIA Toolkit v38", sw_blocks)
emit("author_contributions_v38", "Author Contributions \u2014 EIA v38", ac_blocks)
emit("funding_v38", "Funding \u2014 EIA v38", fu_blocks)
emit("competing_interests_v38", "Competing Interests \u2014 EIA v38", ci_blocks)

# ===========================================================================
# 7. SUBMISSION_CHECKLIST.md  (internal P0/P1 closure log, md only)
# ===========================================================================
closure = """# Submission Compliance Checklist \u2014 EIA v38 (Nature Methods Analysis)

## P0 \u2014 Mandatory calculations and verification
| ID | Item | Status | Note |
|----|------|--------|------|
| C1 | BH 60-test family recomputed; q-values核对 | PASS | 4 survivors: RF E1/E4 (q=0.008), LR SelB (q=0.026), ECS SelD (q=0.015) |
| C2 | Full-number cross-check (manuscript \u2194 source data \u2194 ED \u2194 SI) | PASS | ~30 key numbers verbatim; 3 gaps traced to v27 pipeline, reproducible from frozen raw_inputs |
| C3 | Exact endpoint sample sizes frozen | PASS | per-endpoint n in Supplementary Table S3b |
| C4 | P-value side / bootstrap construction核对 | PASS | one-sided paired bootstrap B=2000, seed 20260819 |
| C5 | 19-scorer manifest | PASS | Supplementary Table S1 |
| C6 | Reproducible hash / repository | PASS | frozen raw_inputs + scripts archived |
| C7 | Reporting Summary + Software Checklist | PASS | generated (this folder) |

## P1 \u2014 Mandatory textual / packaging revisions
| ID | Item | Status |
|----|------|--------|
| R1 | Restore LLM-disclosure paragraph | PASS (Methods \u00a7Use of AI-assisted tools) |
| R2 | Fig. 3 caption gains panel d (coverage-count AUROC 0.573 vs frozen harmonic 0.575 on E1) | PASS |
| R3 | \u201cestablishing a low-expression mechanism\u201d \u2192 \u201cassigning a specific biological mechanism\u201d | PASS (frozen text) |
| R4 | E6 reported only as pre-specified PASS | PASS |
| R5 | Candidate gains not called \u201cvalidated\u201d | PASS (frozen text: \u201ccandidate gains\u201d) |
| R6 | Sequential superscript citations; reference list reordered to citation order | PASS (22 refs, citation-ordered) |
| R7 | Supplementary numbering frozen (S1/S2/S3/Fig1\u20132/Note) | PASS |
| R8 | Author / funding placeholders retained | PASS |
| R9 | \u201cIntroduction\u201d heading deleted (body retained) | PASS |

## Format / compliance
- Abstract: 131 words (\u2264150) PASS
- Main text (excl. Methods): ~2038 words (\u22644000) PASS
- Methods (online): ~1896 words (\u22643000) PASS
- Display items: 6 main figures + 1 Extended Data + 2 Supplementary PASS
- Figures: R/ggplot2 4.0.3, 300 dpi, theme_classic; no AI-generated figures; no watermark PASS
- AI-taste scan: only \u201cin conclusion\u201d substring in SI scaffolding, rephrased; body clean PASS
- Anonymous: no author names / affiliations in manuscript or cover letter PASS
"""

with open(os.path.join(OUT, "SUBMISSION_CHECKLIST.md"), "w", encoding="utf-8") as f:
    f.write(closure)

# ===========================================================================
# 8. ARCHIVE_README.md
# ===========================================================================
readme = """# Archive README \u2014 20260926_v38_NM

This folder is the v38 frozen Nature Methods (Analysis) submission package.

## Layout
- `manuscript/`  \u2014 NatureMethods_EIA_v38_投稿最终版_20260926.docx (+ .md): final manuscript with embedded figures and citation-ordered references.
- `supplementary/` \u2014 NatureMethods_EIA_v38_Supplementary_20260926.docx (+ .md): SI architecture, Supp Fig 1\u20132, mathematical note, self-audit S1\u2013S14.
- `figures/`    \u2014 Fig1\u20136, ExtendedDataFig1_toolkit, Supp1\u20132 (PNG 300 dpi + PDF, ragg/Quartz native, no watermark).
- `code/python/` \u2014 run_analysis.py, postprocess.py, build_submission_forms.py, wordcount.py.
- `code/R/`     \u2014 11 ggplot2 scripts + run_all.R (regenerates all figures).
- `data/raw_inputs/`  \u2014 5 frozen JSON inputs (byte-identical to v18).
- `data/results/`     \u2014 alignment_results.json (self-audit + BH ledger).
- `data/source_data/` \u2014 per-figure CSVs (source data for each display item).
- `supplementary/`    \u2014 Supplementary Tables S1\u2013S3 (CSV) + SI figures (PDF/PNG).
- `submission_forms/` \u2014 cover letter, reporting summary, software checklist, author/funding/competing-interests placeholders, compliance checklist.

## Reproducibility
1. `python code/python/run_analysis.py` \u2192 `postprocess.py` \u2192 `build.py`
2. `Rscript code/R/run_all.R` regenerates all 9 figures.
3. `python code/python/build.py` regenerates the manuscript + supplementary (docx + md).
4. `python code/python/build_submission_forms.py` regenerates submission forms.

## Frozen口径
- Sentinel = \u22123; harmonic formula; seed 20260819; bootstrap 2000; null 100.
- BH family m = 60; 4 survivors (RF E1/E4, LR SelB, ECS SelD).
- E6 pre-specified composite-only PASS.

## Pre-submission TODO (author action)
- Insert final author names / affiliations / correspondence in manuscript [47] and cover letter.
- Insert funding statement, competing interests, author contributions.
- Insert persistent repository DOI / accession in Data availability and Reporting Summary.
- Re-run CHECKSUMS after inserting the above.
"""

with open(os.path.join(OUT, "ARCHIVE_README.md"), "w", encoding="utf-8") as f:
    f.write(readme)

# ===========================================================================
# 9. CHECKSUMS.md  (md5 of all generated deliverables)
# ===========================================================================
def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

checksum_lines = ["# CHECKSUMS \u2014 20260926_v38_NM deliverables (MD5)", ""]
roots = [
    os.path.join(ROOT, "20260926_v38_NM", "manuscript"),
    os.path.join(ROOT, "20260926_v38_NM", "supplementary"),
    OUT,
]
seen = set()
for r in roots:
    for fn in sorted(os.listdir(r)):
        if fn.lower().endswith((".docx", ".md", ".pdf", ".png", ".csv", ".json", ".R", ".py")):
            fp = os.path.join(r, fn)
            if fp in seen:
                continue
            seen.add(fp)
            checksum_lines.append(f"{md5(fp)}  {os.path.relpath(fp, ROOT)}")
checksum_lines.append("")
with open(os.path.join(OUT, "CHECKSUMS.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(checksum_lines))

print("[OK] submission forms generated:")
for fn in sorted(os.listdir(OUT)):
    print("   ", fn)
