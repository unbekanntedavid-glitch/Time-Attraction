# TAL publication set, 2026-11-04, r2

Article *Clock and acceleration predictions of a fixed scalar–tensor branch: calibration and Yukawa constraints*
with the three runnable modules it is computed from. Preliminary research. Method of the checks: numerical
cross-checks with Academia Co-scientist (AI-assisted). The checks are technical (scripted reproduction, numeric
comparison, negative controls); a passing check is not an empirical confirmation.

## Contents

| item | what it is |
|---|---|
| `article/` | the article (EN and HU PDF), LaTeX sources, claim–source table, check scripts, fresh-run logs, change and test logs (see `article/README.md`) |
| `tal_branch_clock_accel_test_r4.zip` | clock and acceleration module (main analysis) |
| `bren_continuation_zenodo_r4.zip` | BREN tower continuation module (bound module, Zenodo variant) |
| `subcm_yukawa_window_zenodo_r4.zip` | sub-centimetre Yukawa window module (bound module, Zenodo variant) |
| `LICENSE.txt` | licence, including the note on third-party material |
| `SHA256SUMS.txt` | SHA-256 of every file of the set |

## How the modules relate

The clock module consumes outputs of the two bound modules: its `inputs/bounds/` folder holds byte-identical copies of
the BREN results file (under the historical name `bren_continuation_results_r2.json`) and of three sub-cm outputs,
recorded with their SHA-256 in `inputs/PROVENANCE_bounds.json`. The three modules are therefore one chain, not
independent evidence. Each archive must be unpacked into its own folder; do not unpack them on top of each other.

## Main commands

```
sha256sum -c SHA256SUMS.txt

mkdir clock bren subcm
unzip -q tal_branch_clock_accel_test_r4.zip -d clock
unzip -q bren_continuation_zenodo_r4.zip    -d bren
unzip -q subcm_yukawa_window_zenodo_r4.zip  -d subcm
(cd clock/tal_branch_clock_accel_test && python3 run_check.py)   # RESULT: PASS
(cd bren/bren_continuation            && python3 run_check.py)   # all steps passed
(cd subcm/subcm_yukawa_window         && python3 run_check.py)   # all steps passed (frozen-digitization mode)

cd article
python3 src/check_claims.py      # every number of the article against the archives
python3 src/check_en_hu.py       # EN/HU agreement
```

The module checks need Python 3 with NumPy, SciPy and matplotlib (versions in each module's `requirements.txt`;
the versions used for the logged runs are in `article/logs/`). The article checks need only Python 3 and poppler
(`pdftotext`, `pdfinfo`). The sub-cm run uses the frozen digitization data; regenerating the BREN DEM from raw
Copernicus tiles is a separate task that the module does not perform.

## Identifiers

| module | DOI |
|---|---|
| TAL work, concept DOI | 10.5281/zenodo.22983681 |
| clock and acceleration | V20 10.5281/zenodo.23540741 (source geometry), V26 10.5281/zenodo.23727585 (branch and action), V28 10.5281/zenodo.23814053 (calibration convention, Step 27) |
| BREN continuation | V28 10.5281/zenodo.23814053 (refinement of Step 27) |
| sub-cm Yukawa window | V28 10.5281/zenodo.23814053 (candidate of Step 27) |

EN
© 2026 Dávid Szabolcsi (ORCID 0009-0005-6300-074X)
License: Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
https://creativecommons.org/licenses/by-nc/4.0/
All commercial, industrial and technological implementation rights are expressly retained by the author.
Commercial use only with the author's written permission. Third-party material remains under its own terms (see LICENSE.txt).
