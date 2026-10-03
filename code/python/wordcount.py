# --- EIA reproducibility: paths are resolved relative to the repository root ---
from pathlib import Path as _P
ROOT = _P(__file__).resolve().parent.parent.parent

# -*- coding: utf-8 -*-
"""Structure-aware word-count audit for the v33 manuscript (Nature Methods
Analysis quotas: Abstract <= 200, main text (Introduction..end of Results/
Discussion, excluding figure legends, Methods, references) <= 4,000,
Online Methods budget tracked separately <= 3,000).

Word convention (frozen in this project): a word matches [A-Za-z][A-Za-z'-]*;
numbers, statistics and symbols are not counted as words.

Usage: python wordcount.py [path/to/manuscript.docx]
"""
import re
import sys
from docx import Document

WORD = re.compile(r"[A-Za-z][A-Za-z'\-]*")
CAPTION = re.compile(r"^(Fig\. \d|Supplementary Fig|Extended Data)")


def words(text):
    return len(WORD.findall(text))


def main(path):
    doc = Document(path)
    paras = [(p.style.name if p.style is not None else "", p.text.strip())
             for p in doc.paragraphs]

    sections = {}
    current = None
    for style, text in paras:
        if style.startswith("Heading"):
            current = text
            sections[current] = []
        elif current is not None:
            sections[current].append(text)

    # locate section boundaries by heading text
    def find(*needles):
        for name in sections:
            if any(n.lower() in name.lower() for n in needles):
                return name
        return None

    h_abs = find("Abstract")
    h_intro = find("Introduction")
    h_methods = find("Methods")
    h_refs = find("References")

    abstract = " ".join(sections.get(h_abs, []))

    def main_text_words():
        names = list(sections)
        i0 = names.index(h_intro)
        i1 = names.index(h_methods) if h_methods in names else len(names)
        total, caps = 0, 0
        for name in names[i0:i1]:
            for t in sections[name]:
                w = words(t)
                total += w
                if CAPTION.match(t):
                    caps += w
        return total, caps, total - caps

    def methods_words():
        if h_methods not in sections:
            return 0, 0, 0
        names = list(sections)
        i0 = names.index(h_methods)
        i1 = names.index(h_refs) if h_refs in names else len(names)
        total, caps = 0, 0
        for name in names[i0:i1]:
            for t in sections[name]:
                w = words(t)
                total += w
                if CAPTION.match(t):
                    caps += w
        return total, caps, total - caps

    a = words(abstract)
    mt, mtc, mtn = main_text_words()
    me, mec, men = methods_words()

    print(f"file: {path}")
    print(f"Abstract                    : {a:6d}  (limit 200)  "
          f"{'PASS' if a <= 200 else 'FAIL'}")
    print(f"Main text (incl. legends)   : {mt:6d}")
    print(f"  figure legend words       : {mtc:6d}")
    print(f"Main text (excl. legends)   : {mtn:6d}  (limit 4,000)  "
          f"{'PASS' if mtn <= 4000 else 'FAIL'}")
    print(f"Online Methods              : {men:6d}  (online budget 3,000)  "
          f"{'PASS' if men <= 3000 else 'OVER'}")
    if a > 200 or mtn > 4000:
        sys.exit(1)


if __name__ == "__main__":
    p = (sys.argv[1] if len(sys.argv) > 1 else
         str(ROOT / "manuscript") + "/"
         "NatureMethods_EIA_v33_投稿前审计修订_20260924.docx")
    main(p)
