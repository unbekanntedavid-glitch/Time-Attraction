# SPDX-License-Identifier: CC-BY-NC-4.0
# Copyright (c) 2026 Szabolcsi David (ORCID 0009-0005-6300-074X). Commercial use only with the author's written permission.
"""Check that the EN and HU article agree in numbers, units, signs and caveats.

usage (from the article folder): python3 src/check_en_hu.py

  1  the same multiset of value macros \\V{id} (numbers and their signs come only from these)
  2  the same counts of every unit macro
  3  the same caveat markers (%CAV:name, one per stated caveat)
  4  the same section, appendix, table, figure and equation labels
  5  the same multiset of numeric tokens in the PDF text (HU decimal comma read as point; footer lines ignored),
     with signs kept
Exit code 0 only if every check passes.
"""
import collections
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import texscan  # noqa: E402
from build_article import PDF  # noqa: E402

UNITS = ["um", "mmu", "cmu", "kmu", "mU", "kgu", "mss", "uGal", "kgm", "Gu", "pct", "degu"]
FOOTER = re.compile(r"ORCID|creativecommons|commercial use|kereskedelmi|page \d+ of|\d+/\d+\. oldal")
results = []


def check(name, ok, detail=""):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))


def diff(a, b):
    ca, cb = collections.Counter(a), collections.Counter(b)
    only_a = ca - cb
    only_b = cb - ca
    return dict(only_a), dict(only_b)


def pdf_numbers(pdf, lang):
    txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True).stdout
    toks = []
    for line in txt.splitlines():
        if FOOTER.search(line):
            continue
        line = line.replace("\u2212", "-")
        if lang == "hu":
            line = re.sub(r"(\d),(\d)", r"\1.\2", line)
        for m in re.finditer(r"-?\d+(?:\.\d+)*", line):
            toks.append(m.group(0).rstrip("."))
    return toks


def main():
    src = {}
    for lang in ("en", "hu"):
        src[lang] = texscan.expand_tables(open(f"tex/article_{lang}.tex", encoding="utf-8").read(), "generated")
    a, b = diff(texscan.macros(src["en"]), texscan.macros(src["hu"]))
    check(f"1 same value macros ({len(texscan.macros(src['en']))} uses)", not a and not b, f"EN only {a}, HU only {b}")
    bad = []
    for u in UNITS:
        n = [len(re.findall(r"\\" + u + r"(?![A-Za-z])", re.sub(r"(?<!\\)%.*", "", src[l]))) for l in ("en", "hu")]
        if n[0] != n[1]:
            bad.append(f"\\{u}: EN {n[0]}, HU {n[1]}")
    check("2 same unit-macro counts", not bad, "; ".join(bad))
    ce, ch = texscan.caveats("tex/article_en.tex"), texscan.caveats("tex/article_hu.tex")
    check(f"3 same caveat markers ({len(ce)})", ce == ch, f"{set(ce) ^ set(ch)}")
    le = re.findall(r"\\label\{([^}]*)\}", src["en"]); lh = re.findall(r"\\label\{([^}]*)\}", src["hu"])
    check(f"4 same labels ({len(le)})", le == lh, f"{set(le) ^ set(lh)}")
    ne, nh = pdf_numbers(PDF["en"], "en"), pdf_numbers(PDF["hu"], "hu")
    a, b = diff(ne, nh)
    check(f"5 same numeric tokens in the PDF text ({len(ne)} tokens)", not a and not b, f"EN only {a}; HU only {b}")
    print(f"{sum(results)}/{len(results)} checks passed")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
