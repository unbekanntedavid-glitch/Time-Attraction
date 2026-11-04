# Test log (article and publication set, 2026-11-04)

All times are Europe/Zurich (CET, UTC+1).

## Environment
Linux 6.12.94+ x86_64, Python 3.13.5, NumPy 2.5.4, SciPy 1.18.1, matplotlib 3.11.2, Pillow 13.0.0,
pdfTeX 3.141592653-2.6-1.40.29 (TeX Live 2026), poppler-utils (pdftotext, pdfinfo, pdftoppm), 8 CPUs.

## 1. Fresh module runs (inputs of the article)
Each r4 module archive was unpacked into its own empty folder and run once. Logs: `logs/fresh_run_*.txt`
(see `logs/README.md` for the archive SHA-256 values and commands).

| module | command | result |
|---|---|---|
| clock (`tal_branch_clock_accel_test_r4.zip`) | `python3 run_check.py --out runs/fresh_article` | `RESULT: PASS`, exit 0 |
| BREN (`bren_continuation_zenodo_r4.zip`) | `python3 run_check.py` | `all steps passed`, exit 0 |
| sub-cm (`subcm_yukawa_window_zenodo_r4.zip`) | `python3 run_check.py` | `all steps passed (frozen-digitization mode)`, exit 0 |

## 2. Article checks
- `python3 src/check_claims.py`: 14/14 checks passed, exit 0 (295 claims re-read from the archives; generated
  files reproduced byte-for-byte; no hand-typed numbers in EN and HU; 7 cross-module equalities; bound files and
  figures byte-identical to the module files; only project DOIs; licence footer on every page; fresh-run logs).
- `python3 src/check_en_hu.py`: 5/5 checks passed, exit 0 (374 value-macro uses, 28 caveat markers, 31 labels,
  942 numeric tokens in the PDF text, identical in EN and HU).

## 3. Negative tests (the checks must fail on a broken copy)
Eight deliberately broken copies; every one made the relevant check fail with exit 1:

| broken copy | failing check |
|---|---|
| hand-typed number in the EN body | check_claims 4 |
| edited value in `generated/values_en.tex` | check_claims 2 |
| edited value in the claim table | check_claims 1 and 2 |
| HU uses another value macro | check_en_hu 1 |
| HU caveat marker removed | check_en_hu 3 |
| one figure changed | check_claims 7 |
| foreign DOI added | check_claims 8 |
| HU value edited and HU PDF recompiled | check_en_hu 5 |

## 4. Set verification from a fresh unzip
- `sha256sum -c SHA256SUMS.txt`: all sums OK; no file outside the list (except `SHA256SUMS.txt` itself).
- check_claims 14/14 and check_en_hu 5/5 in the unpacked set; both PDFs have 10 pages.
- Rebuild in a copy (`python3 src/build_article.py ..`): checks pass again and every generated file, `CLAIMS_SOURCES.csv`
  and `CLAIMS_SOURCES.md` is identical to the shipped copy.
- The three module archives of the set were run again, each in its own empty folder: clock exit 0 / `RESULT: PASS`;
  BREN exit 0 / `all steps passed`; sub-cm exit 0 / `all steps passed (frozen-digitization mode)`.

## 5. Visual check
Every page of both PDFs (10 EN + 10 HU) was rendered at 150 dpi (`pdftoppm -r 150`) and inspected: no overfull
lines, tables and figures readable, footer with licence and page number present on each page.
