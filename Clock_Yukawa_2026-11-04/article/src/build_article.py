# SPDX-License-Identifier: CC-BY-NC-4.0
# Copyright (c) 2026 Szabolcsi David (ORCID 0009-0005-6300-074X). Commercial use only with the author's written permission.
"""Build the article: figures copied byte-for-byte from the module archives, values and tables generated, pdfLaTeX x2.

usage: python3 src/build_article.py <modules_dir>        (run from the article folder)
"""
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from claimlib import Store  # noqa: E402
import make_values  # noqa: E402

FIGURES = [  # (module, path inside the module, file name in figures/)
    ("clock", "docs/generated/fig_delta_vs_lambda.png", "fig_delta_vs_lambda.png"),
    ("bren", "outputs/fig_residual_curve.png", "fig_residual_curve.png"),
    ("subcm", "outputs/figures/subcm_window_en.png", "subcm_window_en.png"),
    ("subcm", "outputs/figures/subcm_window_hu.png", "subcm_window_hu.png"),
]
PDF = {"en": "TAL_clock_accel_yukawa_article_EN.pdf", "hu": "TAL_clock_accel_yukawa_article_HU.pdf"}


def main(modules_dir):
    art = os.getcwd()
    st = Store(modules_dir, art)
    os.makedirs("figures", exist_ok=True)
    for mod, rel, name in FIGURES:
        with open(os.path.join("figures", name), "wb") as f:
            f.write(st.raw(mod, rel))
    make_values.main(modules_dir, art)
    for lang in ("en", "hu"):
        for _ in range(2):
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"article_{lang}.tex"],
                               cwd="tex", capture_output=True, text=True)
            if r.returncode != 0:
                sys.stdout.write(r.stdout[-4000:])
                raise SystemExit(f"pdflatex failed ({lang})")
        log = open(os.path.join("tex", f"article_{lang}.log"), encoding="latin-1").read()
        for bad in ("undefined references", "Rerun to get", "Overfull \\hbox"):
            if bad in log:
                print(f"WARNING {lang}: '{bad}' in log")
        shutil.copyfile(os.path.join("tex", f"article_{lang}.pdf"), PDF[lang])
        print("built", PDF[lang])
    # claim table again, now with the section usage of the compiled sources
    make_values.main(modules_dir, art)


if __name__ == "__main__":
    main(sys.argv[1])
