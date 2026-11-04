# SPDX-License-Identifier: CC-BY-NC-4.0
# Copyright (c) 2026 Szabolcsi David (ORCID 0009-0005-6300-074X). Commercial use only with the author's written permission.
"""Helpers shared by the checks: read the article sources and find hand-typed numbers."""
import os
import re

TABLE_RE = re.compile(r"\\tabinput\s+\.\./generated/tab_(\w+)\.tex")

# Digits that are not quantitative claims: identifiers, citation years, source figure/table/equation numbers, and
# structural integers of formulas. Everything else that contains a digit must be a value macro \V{id}.
WHITELIST = [
    r"\b(P[123]|S[12]|K(?:89|90)|B(?:89|90)|V(?:20|26|28)|H1|r4)\b",
    r"Step~27", r"GLO-90", r"8662/8663", r"8663", r"CODATA~?\s?2018", r"python3",
    r"Fig\.~\d+", r"Table~11\.4", r"eq\.~\d", r"11\.4\.~?táblázat\w*", r"\d+\.~?(?:ábr|táblázat|egyenlet)\w*",
    r"(?:HUST|Eöt-Wash|Lee|AHN|Aeppli|Bothwell|Schlamminger|Takamoto|Grotti|Touboul|Kasameyer|Kammeraad|"
    r"Adelberger et al\.\\ |Lee et al\.\\ )[~ ]?(?:19|20)\d\d(?:/20\d\d)?",
    r"1/\(4\\pi\)", r"G=c=1", r"SHA-256", r"CC BY-NC 4\.0", r"method-1", r"\(0,", r"(?:=|\\ne|>|<)\s*0(?![\d.,])",
    r"(?<![\d.,])1\s*[+-](?=\s*\\?[a-zA-Z(])", r"(?<![\d.,])10\^",
    r"\d*\\pi", r"\^\{?-?\d+(?:/\d+)?\}?", r"_\{?\d+\}?", r"\{2\\alpha_Y\}", r"\\frac\{1\}", r"1-2z/a", r"\b2z\b",
    r"\(1\s*[+-]", r"\[1\s*\+", r"1/\(4\\pi\)", r"1\\sigma", r"(?:\\alpha_Y|R|k_c|\|\\alpha\|)\s*=\s*1(?![\d.,{])",
    r"\b10\^",
]
STRIP = [
    r"(?<!\\)%.*", r"\\V\{[^}]*\}", r"\\label\{[^}]*\}", r"\\(?:eq)?ref\{[^}]*\}", r"\\doi\{[^}]*\}", r"\\url\{[^}]*\}",
    r"\\file\{[^}]*\}", r"\\includegraphics(?:\[[^]]*\])?\{[^}]*\}", r"\\tabinput\s+\S+", r"\\begin\{tabular\}\{[^}]*\}",
    r"\\setlength\{[^}]*\}\{[^}]*\}", r"\\multicolumn\{\d\}\{[^}]*\}", r"\\texttt\{[^}]*\}", r"width=[0-9.]+\\(?:text|line)width",
]


def body(path):
    t = open(path, encoding="utf-8").read()
    b = t.split("\\begin{document}", 1)[1]
    for marker in ("\\section*{References}", "\\section*{Hivatkozások}"):
        if marker in b:
            b = b.split(marker, 1)[0]
    return b


def expand_tables(text, gen_dir):
    def rep(m):
        return open(os.path.join(gen_dir, f"tab_{m.group(1)}.tex"), encoding="utf-8").read()
    return TABLE_RE.sub(rep, text)


def hand_typed(text):
    """Return (line number, token, line) for digits that are neither value macros nor whitelisted."""
    out = []
    for i, line in enumerate(text.splitlines(), 1):
        s = line
        for p in STRIP:
            s = re.sub(p, " ", s)
        for p in WHITELIST:
            s = re.sub(p, " ", s)
        for m in re.finditer(r"\d[\d.,]*", s):
            out.append((i, m.group(0), line.strip()[:140]))
    return out


def macros(text):
    return re.findall(r"\\V\{([^}]+)\}", re.sub(r"(?<!\\)%.*", "", text))


def caveats(path):
    return sorted(re.findall(r"^%CAV:(\w+)", open(path, encoding="utf-8").read(), re.M))
