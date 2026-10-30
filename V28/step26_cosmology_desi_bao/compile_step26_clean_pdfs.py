#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Self-contained LaTeX/PDF builder for V27 (Step 26) HU and EN academic reports.
Uses relative directory paths (os.path.dirname(os.path.abspath(__file__))) and
inlines all table/texttt/equation formatting helpers so it has zero external
local module dependencies.

Exit status (review item F11):
  default (PDF mode): pandoc and pdflatex are checked BEFORE any work; if either is
    missing, or a PDF is not newly written by this run, or pdflatex fails, the script
    exits with a nonzero code (2 = missing tool, 3 = PDF not newly written, 4 = pdflatex error).
  --tex-only: only the .tex files are written (needs pandoc); status "TEX-ONLY",
    which is reported separately and is never a PDF build success.
"""
import sys
import time
import argparse

import os
import re
import shutil
import subprocess
import hashlib

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def strip_latex_for_len(s: str) -> str:
    s = re.sub(r"\\allowbreak\{\}", "", s)
    s = re.sub(r"\\(?:textbf|textit|emph|texttt|text|mathrm|mathbf|mathcal|frac|sqrt|underline)\b", "", s)
    s = re.sub(r"\\[a-zA-Z]+", " ", s)
    s = re.sub(r"[{}^_$~]", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def make_texttt_breakable_in_chunk(inner: str) -> str:
    s = inner.replace(r"\allowbreak{}", "")
    s = s.replace(r"\_", r"\_\allowbreak{}")
    s = re.sub(r"([/.,=:;+])", r"\1\\allowbreak{}", s)
    s = re.sub(r"(?<!-)-(?!-)", r"-\\allowbreak{}", s)
    s = s.replace(r"\ ", r"\ \allowbreak{}")

    def break_long_alnum(m):
        if m.start() > 0 and s[m.start() - 1] == "\\":
            return m.group(0)
        word = m.group(0)
        return r"\allowbreak{}".join(word[i:i+7] for i in range(0, len(word), 7))

    s = re.sub(r"[0-9a-zA-Z]{13,}", break_long_alnum, s)
    return s


def process_all_texttt(tex: str) -> str:
    out = []
    i = 0
    n = len(tex)
    target = r"\texttt{"
    while i < n:
        if tex.startswith(r"\href{", i):
            out.append(r"\href{")
            i += len(r"\href{")
            depth = 1
            while i < n and depth > 0:
                if tex[i] == "{" and tex[i - 1] != "\\":
                    depth += 1
                elif tex[i] == "}" and tex[i - 1] != "\\":
                    depth -= 1
                out.append(tex[i])
                i += 1
            continue
        if tex.startswith(r"\texorpdfstring{", i):
            out.append(r"\texorpdfstring{")
            i += len(r"\texorpdfstring{")
            depth = 1
            first_start = i
            while i < n and depth > 0:
                if tex[i] == "{" and tex[i - 1] != "\\":
                    depth += 1
                elif tex[i] == "}" and tex[i - 1] != "\\":
                    depth -= 1
                i += 1
            first_arg = tex[first_start:i - 1]
            out.append(process_all_texttt(first_arg))
            out.append("}")
            if i < n and tex[i] == "{":
                out.append("{")
                i += 1
                depth = 1
                while i < n and depth > 0:
                    if tex[i] == "{" and tex[i - 1] != "\\":
                        depth += 1
                    elif tex[i] == "}" and tex[i - 1] != "\\":
                        depth -= 1
                    out.append(tex[i])
                    i += 1
            continue
        if tex.startswith(target, i):
            out.append(target)
            i += len(target)
            depth = 1
            start = i
            while i < n and depth > 0:
                if tex[i] == "{" and tex[i - 1] != "\\":
                    depth += 1
                elif tex[i] == "}" and tex[i - 1] != "\\":
                    depth -= 1
                i += 1
            inner = tex[start:i - 1]
            out.append(make_texttt_breakable_in_chunk(inner))
            out.append("}")
        else:
            out.append(tex[i])
            i += 1
    return "".join(out)


def split_top_level_ampersands(row_str: str) -> list:
    cells = []
    cur = []
    depth = 0
    i = 0
    n = len(row_str)
    while i < n:
        ch = row_str[i]
        if ch == "\\" and i + 1 < n:
            cur.append(row_str[i:i + 2])
            i += 2
            continue
        if ch == "{":
            depth += 1
        elif ch == "}" and depth > 0:
            depth -= 1
        elif ch == "&" and depth == 0:
            cells.append("".join(cur))
            cur = []
            i += 1
            continue
        cur.append(ch)
        i += 1
    cells.append("".join(cur))
    return cells


def split_top_level_row_breaks(body_str: str) -> list:
    rows = []
    cur = []
    depth = 0
    i = 0
    n = len(body_str)
    while i < n:
        if body_str[i] == "\\" and i + 1 < n:
            nxt = body_str[i + 1]
            if nxt == "\\" and depth == 0:
                rows.append("".join(cur))
                cur = []
                i += 2
                continue
            else:
                cur.append(body_str[i:i + 2])
                i += 2
                continue
        if body_str[i] == "{":
            depth += 1
        elif body_str[i] == "}" and depth > 0:
            depth -= 1
        cur.append(body_str[i])
        i += 1
    if "".join(cur).strip():
        rows.append("".join(cur))
    return rows


def relax_inline_math_in_cell(cell: str, num_cols: int = 4) -> str:
    def fix_math(m):
        content = m.group(1)
        content = re.sub(r"\\left\s*([(\[\{|])", r"\1", content)
        content = re.sub(r"\\right\s*([)\]\}\|])", r"\1", content)
        content = content.replace(r"\left\lvert", r"\lvert").replace(r"\right\rvert", r"\rvert")
        content = content.replace(r"\left\langle", r"\langle").replace(r"\right\rangle", r"\rangle")
        content = content.replace(r"\left.", "").replace(r"\right.", "")
        content = re.sub(r"\\(?:big|Big|bigg|Bigg)[lr]?\s*([(\[\{\)\|\]\}])", r"\1", content)
        if num_cols >= 8 and r"\frac{" in content and len(content) > 35:
            brace_pat = r"(?:[^{}]|\{(?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*\})+"
            content = re.sub(
                rf"\\frac\{{({brace_pat})\}}\{{({brace_pat})\}}",
                lambda m2: f"({m2.group(1)})/\\allowbreak{{}}({m2.group(2)})",
                content,
            )
        content = content.replace(r"\;=\;", "=").replace(r"\;\equiv\;", r"\equiv ").replace(r"\;\approx\;", r"\approx ")
        content = content.replace(r"\;\implies\;", r"\implies ").replace(r"\;\Longrightarrow\;", r"\Longrightarrow ")
        content = content.replace(r"\text{–}", r"\text{–}\allowbreak{}")
        content = content.replace(r"\text{--}", r"\text{--}\allowbreak{}")
        content = content.replace("/", r"/\allowbreak{}")
        if num_cols >= 6:
            content = re.sub(r"\\!\s*\(", "(", content)
            content = content.replace("(", r"\allowbreak{}(\allowbreak{}")
            content = content.replace(r"\lvert", r"\allowbreak{}\lvert\allowbreak{}")
        content = re.sub(r",\s*", r",\\allowbreak ", content)
        return r"\(" + content + r"\)"

    return re.sub(r"\\\((.*?)\\\)", fix_math, cell, flags=re.S)


def process_longtable(match: re.Match) -> str:
    lt = match.group(0)
    col_specs = re.findall(r">\{[^}]*\}p\{\(\\linewidth\s*-\s*\d+\\tabcolsep\)\s*\*\s*\\real\{([0-9.]+)\}\}", lt)
    num_cols = len(col_specs)
    if num_cols == 0 or "\\endlastfoot" not in lt:
        return lt

    lt = re.sub(
        r"\\begin\{minipage\}\[b\]\{\\linewidth\}\s*\\(?:raggedleft|raggedright|centering)",
        r"\\begin{minipage}[b]{\\linewidth}\\raggedright\\hspace{0pt}",
        lt,
    )
    pre_foot, post_foot = lt.split("\\endlastfoot", 1)
    body, _ = post_foot.rsplit("\\end{longtable}", 1)

    raw_rows = split_top_level_row_breaks(body)
    parsed_rows = []
    for r in raw_rows:
        r_clean = r.strip()
        if not r_clean or r_clean == r"\bottomrule\noalign{}":
            continue
        if r"\bottomrule" in r_clean:
            r_clean = r_clean.replace(r"\bottomrule\noalign{}", "").replace(r"\bottomrule", "").strip()
        if r_clean:
            cells = split_top_level_ampersands(r_clean)
            cells = [relax_inline_math_in_cell(c.strip(), num_cols) for c in cells]
            parsed_rows.append(cells)

    if not parsed_rows:
        return lt

    col_max_len = [0] * num_cols
    col_sum_len = [0] * num_cols
    header_part = pre_foot.split("\\toprule")[-1] if "\\toprule" in pre_foot else pre_foot
    header_cells = split_top_level_ampersands(header_part)
    for c_idx in range(min(num_cols, len(header_cells))):
        st = strip_latex_for_len(header_cells[c_idx])
        col_max_len[c_idx] = max(col_max_len[c_idx], len(st))

    for cells in parsed_rows:
        for c_idx in range(min(num_cols, len(cells))):
            st = strip_latex_for_len(cells[c_idx])
            L = len(st)
            col_max_len[c_idx] = max(col_max_len[c_idx], L)
            col_sum_len[c_idx] += L

    n_rows = max(1, len(parsed_rows))
    weights = []
    for c_idx in range(num_cols):
        avg_l = col_sum_len[c_idx] / n_rows
        max_l = col_max_len[c_idx]
        score = (avg_l + 0.35 * max_l + 4.0) ** 0.52
        if max_l <= 7:
            min_w, max_w = 0.055, 0.075
        elif max_l <= 14:
            min_w, max_w = 0.075, 0.11
        elif max_l <= 24:
            min_w, max_w = 0.10, 0.16
        else:
            min_w, max_w = 0.11, 0.56
        weights.append((score, min_w, max_w))

    total_Target = 0.98
    raw_scores = [w[0] for w in weights]
    sum_scores = sum(raw_scores) or 1.0
    reals = [total_Target * s / sum_scores for s in raw_scores]

    for _ in range(5):
        fixed_sum = 0.0
        free_indices = []
        for idx, val in enumerate(reals):
            mn, mx = weights[idx][1], weights[idx][2]
            if val < mn:
                reals[idx] = mn
                fixed_sum += mn
            elif val > mx:
                reals[idx] = mx
                fixed_sum += mx
            else:
                free_indices.append(idx)
        rem = total_Target - fixed_sum
        if free_indices and rem > 0:
            free_sum = sum(reals[i] for i in free_indices) or 1.0
            for i in free_indices:
                reals[i] = reals[i] * (rem / free_sum)
        else:
            s_all = sum(reals) or 1.0
            reals = [v * (total_Target / s_all) for v in reals]
            break

    spec_idx = [0]
    def repl_spec(m):
        idx = spec_idx[0]
        r_val = reals[idx] if idx < len(reals) else (0.98 / num_cols)
        spec_idx[0] += 1
        tab_mult = 2 * num_cols - 2
        return f">{{\\raggedright\\arraybackslash\\hspace{{0pt}}}}p{{(\\linewidth - {tab_mult}\\tabcolsep) * \\real{{{r_val:.4f}}}}}"

    pre_foot_new = re.sub(
        r">\{[^}]*\}p\{\(\\linewidth\s*-\s*\d+\\tabcolsep\)\s*\*\s*\\real\{[0-9.]+\}\}",
        repl_spec,
        pre_foot,
    )
    pre_foot_new = relax_inline_math_in_cell(pre_foot_new, num_cols)

    new_body_rows = []
    for cells in parsed_rows:
        if len(cells) < num_cols:
            cells = cells + [""] * (num_cols - len(cells))
        new_body_rows.append(" & ".join(cells) + r" \\")

    rebuilt_lt = (
        pre_foot_new
        + "\\endlastfoot\n"
        + "\n".join(new_body_rows)
        + "\n\\bottomrule\\noalign{}\n\\end{longtable}"
    )

    if num_cols >= 12:
        wrapper_open = "{\\tiny\\setlength{\\tabcolsep}{1.5pt}\\renewcommand{\\arraystretch}{1.10}\n"
    elif num_cols >= 6:
        wrapper_open = "{\\scriptsize\\setlength{\\tabcolsep}{2.0pt}\\renewcommand{\\arraystretch}{1.12}\n"
    else:
        wrapper_open = "{\\footnotesize\\setlength{\\tabcolsep}{3.0pt}\\renewcommand{\\arraystretch}{1.13}\n"
    return wrapper_open + rebuilt_lt + "\n}\n"


def wrap_display_equations(tex: str) -> str:
    def repl_eq(m):
        inner = m.group(1).strip()
        return f"\\fitdisplay{{{inner}}}"
    return re.sub(r"\\\[(.*?)\\\]", repl_eq, tex, flags=re.S)


def find_pdflatex():
    cand = shutil.which("pdflatex")
    if os.environ.get("TAL_PDFLATEX_FALLBACK", "1") == "0":
        return cand  # PATH only (used to test the missing-compiler exit code)
    if cand:
        return cand
    for p in [
        "/home/box/.TinyTeX/bin/x86_64-linux/pdflatex",
        "/root/.TinyTeX/bin/x86_64-linux/pdflatex",
        os.path.expanduser("~/.TinyTeX/bin/x86_64-linux/pdflatex"),
    ]:
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return None


def build_pdf(lang: str, tex_only: bool = False):
    pandoc_bin = shutil.which("pandoc")
    pdflatex_bin = find_pdflatex()
    if not pandoc_bin:
        print(f"[{lang}] ERROR: pandoc not found; no TeX/PDF written.")
        sys.exit(2)
    if not tex_only and not pdflatex_bin:
        print(f"[{lang}] ERROR: pdflatex not found; PDF build requested but impossible (use --tex-only for TeX only).")
        sys.exit(2)

    if lang == "HU":
        md_file = "step26_desi_bao_planck_sne_rsd_bbn_nogo_hu.md"
        base_name = "step26_desi_bao_planck_sne_rsd_bbn_nogo_hu"
        babel_line = r"\usepackage[magyar]{babel}"
        pdf_title = "Idővonzástörvény (TAL) — V27 26. lépés: Kozmológiai H(z), DESI DR1 BAO, Planck 2018 CMB, Szintetikus SNe Ia Alakkontroll, RSD f*sigma8(z), Algebrai BBN és Inverz Távolságlétra No-Go Audit"
        pdf_subj = "26. lépés (V27): Kozmológiai audit és No-Go függelék L = 10 km, alpha_Y = 1 mellett"
        head_l = r"Id\H{o}vonz\'ast\"orv\'eny --- Skal\'ar--tenzor Yukawa-illeszt\'es, V27 (26. l\'ep\'es)"
        head_r = r"V26 DOI: \href{https://doi.org/10.5281/zenodo.23727585}{\texttt{10.5281/zenodo.23727585}} (\texttt{2026-10-28})"
        foot_l = r"\textcopyright\ 2026 Szabolcsi D\'avid (ORCID: 0009-0005-6300-074X) $\cdot$ Licenc: CC BY-NC 4.0 International"
        page_lbl = r"Oldal \thepage"
        box_top = r"Akad\'emiai B\'{\i}r\'alati Kiad\'as $\cdot$ Kozmol\'ogiai \'es Hubble-fesz\"ults\'eg Audit $\cdot$ 26. Modul (V27)"
        box_title = r"Id\H{o}vonz\'ast\"orv\'eny --- Skal\'ar--tenzor Yukawa-illeszt\'es (V27)"
        box_sub = r"26. l\'ep\'es: Kozmol\'ogiai $H(z)$ ($0 < z < 2{,}5$), DESI DR1 BAO, Planck 2018 CMB, Szintetikus SNe Ia $\mu(z)$,\\Vari\'aci\'os RSD $f\sigma_8(z)$, Algebrai BBN \'es Inverz T\'avols\'agl\'etra No-Go Audit ($L = 10\;\text{km}$, $\alpha_Y = 1$)"
        author_line = r"\textbf{Szabolcsi D\'avid} \quad (\texttt{ORCID: 0009-0005-6300-074X}) --- F\"uggetlen kutat\'o, Unterterzen, Sv\'ajc"
        meta_rows = r"""
      \textbf{Zenodo V26 DOI:} & \href{https://doi.org/10.5281/zenodo.23727585}{\texttt{10.5281/zenodo.23727585}} (\textbf{V26}, 1--25. l\'ep\'es; GitHub: \href{https://github.com/unbekanntedavid-glitch/Time-Attraction/releases/tag/V26}{\texttt{Time-Attraction/releases/tag/V26}}) \\
      \textbf{Kor\'abbi DOI-k:} & \textbf{V25}: \texttt{10.5281/zenodo.23699685} (szigor\'uan V25); \textbf{V24}: \texttt{10.5281/zenodo.23680715} (szigor\'uan V24) \\
      \textbf{Modul (\texttt{step26}):} & DESI DR1 BAO ($12\times 12$ kovariancia + profilkontroll), modellspecifikus $r_s, r_d$, szintetikus SNe $\mu(z)$, vari\'aci\'os RSD $f\sigma_8(z)$, algebrai BBN \'es No-Go F\"uggel\'ek \\
      \textbf{D\'atum \'es Licenc:} & 2026-10-28 $\cdot$ CC BY-NC 4.0 International
"""
    else:
        md_file = "step26_desi_bao_planck_sne_rsd_bbn_nogo_en.md"
        base_name = "step26_desi_bao_planck_sne_rsd_bbn_nogo_en"
        babel_line = r"\usepackage[english]{babel}"
        pdf_title = "Time Attraction Law (TAL) — V27 Step 26: Cosmological H(z), DESI DR1 BAO, Planck 2018 CMB, Synthetic SNe Ia Shape Control, Variational RSD f*sigma8(z), Algebraic BBN, and Inverse Distance Ladder No-Go Audit"
        pdf_subj = "Step 26 (V27): Cosmological Audit and No-Go Appendix at fixed L = 10 km, alpha_Y = 1"
        head_l = r"Time Attraction Law --- Scalar--Tensor Yukawa Fit, V27 (Step 26)"
        head_r = r"V26 DOI: \href{https://doi.org/10.5281/zenodo.23727585}{\texttt{10.5281/zenodo.23727585}} (\texttt{2026-10-28})"
        foot_l = r"\textcopyright\ 2026 D\'avid Szabolcsi (ORCID: 0009-0005-6300-074X) $\cdot$ License: CC BY-NC 4.0 International"
        page_lbl = r"Page \thepage"
        box_top = r"Academic Review Edition $\cdot$ Cosmological \& Hubble-Tension Audit $\cdot$ Module 26 (V27)"
        box_title = r"Time Attraction Law --- Scalar--Tensor Yukawa Fit (V27)"
        box_sub = r"Step 26: Cosmological $H(z)$ ($0 < z < 2.5$), DESI DR1 BAO, Planck 2018 CMB, Synthetic SNe Ia $\mu(z)$,\\Variational RSD $f\sigma_8(z)$, Algebraic BBN, and Inverse Distance Ladder No-Go Audit ($L = 10\;\text{km}$, $\alpha_Y = 1$)"
        author_line = r"\textbf{D\'avid Szabolcsi} \quad (\texttt{ORCID: 0009-0005-6300-074X}) --- Independent Researcher, Unterterzen, Switzerland"
        meta_rows = r"""
      \textbf{Zenodo V26 DOI:} & \href{https://doi.org/10.5281/zenodo.23727585}{\texttt{10.5281/zenodo.23727585}} (\textbf{V26}, Steps 1--25; GitHub: \href{https://github.com/unbekanntedavid-glitch/Time-Attraction/releases/tag/V26}{\texttt{Time-Attraction/releases/tag/V26}}) \\
      \textbf{Earlier DOIs:} & \textbf{V25}: \texttt{10.5281/zenodo.23699685} (strictly V25); \textbf{V24}: \texttt{10.5281/zenodo.23680715} (strictly V24) \\
      \textbf{Module (\texttt{step26}):} & DESI DR1 BAO ($12\times 12$ covariance + profile control), model-specific $r_s, r_d$, synthetic SNe $\mu(z)$, variational RSD $f\sigma_8(z)$, algebraic BBN, and No-Go Appendix \\
      \textbf{Date \& License:} & 2026-10-28 $\cdot$ CC BY-NC 4.0 International
"""

    raw_tex_body = subprocess.check_output(
        [pandoc_bin, os.path.join(OUT_DIR, md_file), "-f", "markdown", "-t", "latex", "--no-highlight"],
        cwd=OUT_DIR
    ).decode("utf-8")

    raw_tex_body = re.sub(r"\\pandocbounded\{(.*?)\}", r"\1", raw_tex_body, flags=re.S)
    raw_tex_body = re.sub(r"\\includegraphics\[[^\]]*\]", r"\\includegraphics[width=0.95\\linewidth]", raw_tex_body)

    preamble = rf"""\documentclass[10pt,a4paper]{{article}}
\usepackage[utf8]{{inputenc}}
\usepackage[T1]{{fontenc}}
{babel_line}
\usepackage{{lmodern}}
\usepackage{{textcomp}}
\usepackage{{microtype}}
\usepackage[a4paper,top=16mm,bottom=18mm,left=14mm,right=14mm,headheight=14pt,headsep=5mm,footskip=9mm]{{geometry}}
\usepackage{{amsmath,amssymb,amsfonts,mathtools}}
\usepackage{{booktabs,longtable,array,calc,graphicx,float}}
\usepackage[table]{{xcolor}}
\usepackage{{fancyhdr,upquote,xurl,fvextra}}
\usepackage[hidelinks,pdfusetitle=false]{{hyperref}}
\usepackage{{bookmark}}

\definecolor{{navy}}{{HTML}}{{12284C}}
\definecolor{{slate}}{{HTML}}{{334E68}}
\definecolor{{boxbg}}{{HTML}}{{F4F7FA}}
\definecolor{{slatebg}}{{HTML}}{{F8FAFC}}

\hypersetup{{
  colorlinks=true,
  linkcolor=navy,
  urlcolor=navy,
  citecolor=navy,
  pdftitle={{{pdf_title}}},
  pdfauthor={{Dávid Szabolcsi (ORCID: 0009-0005-6300-074X)}},
  pdfsubject={{{pdf_subj}}},
  pdfcreator={{LaTeX with hyperref}}
}}
\pdfsuppressptexinfo=-1 \pdftrailerid{{}}

\pagestyle{{fancy}}
\fancyhf{{}}
\fancyhead[L]{{\footnotesize\scshape\color{{slate}}{head_l}}}
\fancyhead[R]{{\footnotesize\color{{slate}}{head_r}}}
\fancyfoot[L]{{\scriptsize\color{{slate}}{foot_l}}}
\fancyfoot[R]{{\footnotesize\bfseries\color{{navy}}{page_lbl}}}
\renewcommand{{\headrulewidth}}{{0.4pt}}
\renewcommand{{\footrulewidth}}{{0.4pt}}
\fancypagestyle{{plain}}{{}}

\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{5pt plus 2pt minus 1pt}}
\setlength{{\emergencystretch}}{{4em}}
\setlength{{\LTpre}}{{4pt}}
\setlength{{\LTpost}}{{4pt}}
\providecommand{{\tightlist}}{{\setlength{{\itemsep}}{{3pt}}\setlength{{\parskip}}{{1pt}}}}
\makeatletter\def\fps@figure{{htbp}}\makeatother
\renewcommand{{\arraystretch}}{{1.15}}

\setcounter{{secnumdepth}}{{0}}
\allowdisplaybreaks[4]
\relpenalty=50
\binoppenalty=50
\tolerance=2000
\hbadness=10000
\vbadness=10000
\hfuzz=1pt
\vfuzz=1pt

\newsavebox{{\eqbox}}
\newcommand{{\fitdisplay}}[1]{{%
  \[%
    \sbox{{\eqbox}}{{$\displaystyle #1$}}%
    \ifdim\wd\eqbox>\linewidth
      \resizebox{{\linewidth}}{{!}}{{\usebox{{\eqbox}}}}%
    \else
      \usebox{{\eqbox}}%
    \fi
  \]%
}}

\begin{{document}}

\begin{{center}}
\fcolorbox{{navy}}{{boxbg}}{{%
  \begin{{minipage}}{{0.96\linewidth}}
    \vspace{{3.5mm}}
    \centering
    {{\footnotesize\scshape\color{{slate}}{box_top}\par}}
    \vspace{{1.5mm}}
    {{\LARGE\bfseries\color{{navy}}{box_title}\par}}
    \vspace{{1.5mm}}
    {{\normalsize\bfseries\color{{slate}}{box_sub}\par}}
    \vspace{{2.5mm}}
    {{\small {author_line}\par}}
    \vspace{{2mm}}
    \begin{{tabular}}{{@{{}}r >{{\raggedright\arraybackslash}}p{{0.68\linewidth}}@{{}}}}
{meta_rows.strip()}
    \end{{tabular}}
    \vspace{{2.5mm}}
  \end{{minipage}}%
}}
\end{{center}}
\vspace{{2mm}}
"""

    raw_tex_body = process_all_texttt(raw_tex_body)
    raw_tex_body = re.sub(r"\\begin\{longtable\}.*?\\end\{longtable\}", process_longtable, raw_tex_body, flags=re.S)
    raw_tex_body = wrap_display_equations(raw_tex_body)
    full_tex = preamble + "\n" + raw_tex_body + "\n\\end{document}\n"

    tex_path = os.path.join(OUT_DIR, f"{base_name}.tex")
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write(full_tex)

    if tex_only:
        print(f"[{lang}] STATUS TEX-ONLY: wrote {base_name}.tex; no PDF compiled (not a PDF build success).")
        return "TEX-ONLY"

    pdf_path = os.path.join(OUT_DIR, f"{base_name}.pdf")
    t_start_ns = time.time_ns()
    for _ in range(2):
        r = subprocess.run([pdflatex_bin, "-interaction=nonstopmode", f"{base_name}.tex"], cwd=OUT_DIR, stdout=subprocess.DEVNULL)
        if r.returncode != 0:
            print(f"[{lang}] ERROR: pdflatex exit code {r.returncode}")
            sys.exit(4)
    if not os.path.exists(pdf_path) or os.stat(pdf_path).st_mtime_ns < t_start_ns:
        print(f"[{lang}] ERROR: {base_name}.pdf was not newly written by this run.")
        sys.exit(3)

    log_path = os.path.join(OUT_DIR, f"{base_name}.log")
    log_txt = open(log_path, encoding="utf-8", errors="ignore").read()
    ov_hbox = [m.group(0) for m in re.finditer(r"Overfull \\hbox \(([0-9.]+)pt too wide\)", log_txt) if float(m.group(1)) > 10.0]
    ov_vbox = re.findall(r"Overfull \\vbox.*", log_txt)
    for ext in [".aux", ".out", ".toc", ".log"]:
        p_ext = os.path.join(OUT_DIR, f"{base_name}{ext}")
        if os.path.exists(p_ext):
            os.remove(p_ext)

    sha = hashlib.sha256(open(pdf_path, "rb").read()).hexdigest()
    sz = os.path.getsize(pdf_path)
    n_pages = -1
    if shutil.which("pdfinfo"):
        pages = subprocess.check_output(["pdfinfo", pdf_path]).decode("utf-8")
        p_m = re.search(r"Pages:\s+(\d+)", pages)
        n_pages = int(p_m.group(1)) if p_m else -1
    print(f"[{lang}] STATUS PDF-OK: {base_name}.pdf: {sz} B, {n_pages} pp, SHA256={sha}, Overfull hbox>10pt={len(ov_hbox)}, Overfull vbox={len(ov_vbox)}")
    return "PDF-OK"


def main():
    ap = argparse.ArgumentParser(description="Build the Step 26 HU/EN TeX and PDF files.")
    ap.add_argument("--tex-only", action="store_true", help="write .tex only (separate TEX-ONLY status, not a PDF build)")
    a = ap.parse_args()
    st = [build_pdf("HU", a.tex_only), build_pdf("EN", a.tex_only)]
    print(f"compile_step26_clean_pdfs: status {'/'.join(st)}")


if __name__ == "__main__":
    main()
