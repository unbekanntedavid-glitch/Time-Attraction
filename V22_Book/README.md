# Time Attraction Law — V22 Book Edition (Research Monograph)

- **Zenodo DOI (V22 Book Edition):** [10.5281/zenodo.23604423](https://doi.org/10.5281/zenodo.23604423)
- **Title on Zenodo:** *Time Attraction Law - V22 Book Edition: Monograph (15 pages, 12 chapters + Appendix revision3) - C = 4πG/c² - line-averaged time attraction*
- **Author:** David Szabolcsi (ORCID: [0009-0005-6300-074X](https://orcid.org/0009-0005-6300-074X))
- **Publication Date:** 2026-10-24
- **License:** Creative Commons Attribution-NonCommercial 4.0 International (`CC BY-NC 4.0`). All commercial, industrial, and technological rights are reserved by David Szabolcsi.
- **Previous Versions:**
  - Zenodo V21 (= GitHub V19/V21): [10.5281/zenodo.23602501](https://doi.org/10.5281/zenodo.23602501)
  - Zenodo V20 (= GitHub V18): [10.5281/zenodo.23540741](https://doi.org/10.5281/zenodo.23540741)
  - Zenodo V19 (7th model): [10.5281/zenodo.23497128](https://doi.org/10.5281/zenodo.23497128)
  - Zenodo V17 (SPARC & H.E.S.S. audit): [10.5281/zenodo.19666220](https://doi.org/10.5281/zenodo.19666220)
  - Zenodo V14: [10.5281/zenodo.23217189](https://doi.org/10.5281/zenodo.23217189)

## Overview & Scope

This directory contains the 15-page, 12-chapter (+ 2-page Round-3 technical reproduction audit appendix) research monograph on the **Time Attraction Law** across Versions 1 through 21:
- **Primary English Monograph (`tal_book_EN_2026-10-24.pdf` & `.tex`):** Byte-for-byte identical to `tal_book_EN_2026-10-24 (1).pdf` on Zenodo record `23604423` (MD5 `8af3aac2ee04538dd40b6907f8bcc8c0`, SHA-256 `7e36bcc669840d9165c2d857b9851885408f370a56e2c310eb41ad175f38873a`).
- **Hungarian Edition (`tal_book_HU_2026-10-24.pdf` & `.tex`):** Byte-for-byte identical to `tal_book_HU_2026-10-24.pdf` on Zenodo record `23604423` (MD5 `ae47aaa9b2a25f6e252ff01cc9c3a0e8`, SHA-256 `cf18bef5718ba995471be8c7e57505524710f5091a9037d5ee1163ce48c371fc`).
- **German Edition (`tal_book_DE_2026-10-24.pdf` & `.tex`):** 15-page German translation matching the English and Hungarian monographs page for page (MD5 `42fc33f6ecec559c7d264976715fe052`, SHA-256 `9474f16f08f742cfcc3e9c28db5a92e061198f757271abef3e82e775ad5520b8`).

The calculations summarized in this monograph test selected equilibrium and observational couplings; they do not establish a closed covariant theory or an experimental detection of a new clock effect. Algebraic results are separated from model assumptions, and numerical checks from experimental evidence: the archived calculations document what has been tested, while open tasks (a)–(d) (specifically covariant action closure, task `50aa2db0`) state what is still required to link those results within a single dynamics.

## Contents of `V22_Book/`

| File | Description | SHA-256 |
|---|---|---|
| `tal_book_EN_2026-10-24.pdf` | Primary 15-page English Research Monograph (PDF) | `7e36bcc669840d9165c2d857b9851885408f370a56e2c310eb41ad175f38873a` |
| `tal_book_EN_2026-10-24.tex` | LaTeX source of the 15-page English Monograph | `0d273486120cf22247183634427f22c9e6c7a4682dfafc946290fc15d1c149e4` |
| `tal_book_HU_2026-10-24.pdf` | 15-page Hungarian Research Monograph (PDF) | `cf18bef5718ba995471be8c7e57505524710f5091a9037d5ee1163ce48c371fc` |
| `tal_book_HU_2026-10-24.tex` | LaTeX source of the 15-page Hungarian Monograph | `ee37fa95fd8011242c608d2f07e921182f5f6c53fbccd93cc9da307799aaa1db` |
| `tal_book_DE_2026-10-24.pdf` | 15-page German Research Monograph (PDF) | `9474f16f08f742cfcc3e9c28db5a92e061198f757271abef3e82e775ad5520b8` |
| `tal_book_DE_2026-10-24.tex` | LaTeX source of the 15-page German Monograph | `c97825ee312ce2a15d2d4d402d8f0a9b9aa346ad39730a01b15feb67c253aaad` |
| `generate_figures.py` | Reproducible Python script generating all 12 figures (`EN`, `HU`, `DE`) | See `SHA256SUMS.txt` |
| `figures/*.png` | 12 high-resolution figures (SPARC, TOV, radial stability, 8-coil & granite lab) | See `SHA256SUMS.txt` |
| `revision3_review_EN.md` / `HU.md` | Verbatim Round-3 AI-assisted technical reproduction audit texts | See `SHA256SUMS.txt` |

## Verification & Reproduction

```bash
cd V22_Book
sha256sum -c SHA256SUMS.txt
python3 generate_figures.py
pdflatex -interaction=nonstopmode tal_book_EN_2026-10-24.tex
pdflatex -interaction=nonstopmode tal_book_EN_2026-10-24.tex
```
