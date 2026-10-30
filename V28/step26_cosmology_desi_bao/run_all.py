#!/usr/bin/env python3
"""
Step 26 (V27) reproduction runner.

Order:
  0. SHA-256 check of the two DESI DR1 input files (Cobaya bao_data mean and covariance).
  1. Block 1 via hubble_desi_bao_tal_audit.py (run_step1_1, run_step1_2_and_1_3, run_step1_4_to_1_6,
     then step 2 and the re-audit controls)            -- skipped with --skip-block1
  2. run_step2_nogo_and_report.py                       -- M1-M4 branches, inverse ladder CSV and figures
  3. run_step26_revision_controls.py                    -- re-audit controls R1-R6, 5.1-5.3, strict JSON merge
  4. generate_step26_markdown.py                        -- HU/EN Markdown and README from JSON/CSV
     test_step26_text_regression.py                     -- HU/EN text regression test (review T01), nonzero exit on failure
  5. compile_step26_clean_pdfs.py                       -- only with --pdf (needs pandoc + pdflatex)
                                                           or --tex-only (TeX files only, separate status)

Exit status (review item F11): with --pdf a dependency pre-check runs first; a missing
pandoc/pdflatex, a failed compile or a PDF that is not newly written gives a nonzero exit.
The numerical run status and the PDF build status are printed separately.

Usage: python3 run_all.py [--skip-block1] [--pdf | --tex-only]
"""
import os
import sys
import json
import hashlib
import argparse
import subprocess

D = os.path.dirname(os.path.abspath(__file__))

DATA_SHA256 = {
    "desi_2024_gaussian_bao_ALL_GCcomb_mean.txt": "dd2873a0b88459a491af3c0c0307ba059f62df9211d5b976760f310565a1be68",
    "desi_2024_gaussian_bao_ALL_GCcomb_cov.txt": "bbafa9074b51cf1a45e0d10e4f37db8c0e80a5d1d1788857abb7fc49fb21abcc",
}


def check_data():
    ok = True
    for name, ref in DATA_SHA256.items():
        p = os.path.join(D, name)
        if not os.path.exists(p):
            print(f"[data] MISSING {name}")
            ok = False
            continue
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        status = "OK" if h == ref else "MISMATCH"
        ok &= (h == ref)
        print(f"[data] {name}: {status} ({h})")
    return ok


def run(script, *args):
    print(f"\n=== python3 {script} {' '.join(args)} ===", flush=True)
    r = subprocess.run([sys.executable, script, *args], cwd=D)
    if r.returncode != 0:
        sys.exit(f"{script} failed with exit code {r.returncode}")


def check_json():
    def bad(name):
        raise ValueError(f"non-strict JSON constant {name}")
    with open(os.path.join(D, "tal_hubble_desi_bao_summary.json"), encoding="utf-8") as fh:
        json.loads(fh.read(), parse_constant=bad)
    print("[json] tal_hubble_desi_bao_summary.json: strict JSON OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-block1", action="store_true", help="do not run Block 1 (use the shipped CSV/JSON)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--pdf", action="store_true", help="also compile the PDFs (pandoc + pdflatex); nonzero exit if impossible")
    g.add_argument("--tex-only", action="store_true", help="write the .tex files only (status TEX-ONLY, not a PDF build)")
    a = ap.parse_args()
    if a.pdf:
        sys.path.insert(0, D)
        import shutil
        from compile_step26_clean_pdfs import find_pdflatex
        missing = [n for n, p in (("pandoc", shutil.which("pandoc")), ("pdflatex", find_pdflatex())) if not p]
        if missing:
            print(f"[pdf] ERROR: --pdf requested but missing: {', '.join(missing)}")
            sys.exit(2)
    if not check_data():
        sys.exit("DESI input hash check failed")
    if not a.skip_block1:
        run("hubble_desi_bao_tal_audit.py")
    run("run_step2_nogo_and_report.py")
    run("run_step26_revision_controls.py")
    check_json()
    run("generate_step26_markdown.py")
    run("test_step26_text_regression.py")
    print("\nrun_all: numerical run and Markdown generation OK")
    if a.pdf:
        run("compile_step26_clean_pdfs.py")
        print("run_all: PDF build status PDF-OK")
    elif a.tex_only:
        run("compile_step26_clean_pdfs.py", "--tex-only")
        print("run_all: PDF build status TEX-ONLY (no PDF compiled)")
    else:
        print("run_all: PDF build status NOT-REQUESTED")
    print("run_all: done")


if __name__ == "__main__":
    main()
