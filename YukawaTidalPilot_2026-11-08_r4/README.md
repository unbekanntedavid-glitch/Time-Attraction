# Static tensor tidal response in a local massive-scalar neutron-star model: a fixed-laboratory-G SLy/APR4 pilot (revision r4)

Revision r4. Publication date: 2026-11-08. Author: Szabolcsi Dávid (ORCID [0009-0005-6300-074X](https://orcid.org/0009-0005-6300-074X)).

Zenodo DOI: [10.5281/zenodo.24057132](https://doi.org/10.5281/zenodo.24057132)  
Concept DOI: [10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681) · Zenodo version index 38.

This record is the static tidal pilot. It is not the binding-energy revision. DOI [10.5281/zenodo.24053942](https://doi.org/10.5281/zenodo.24053942) remains YBE-2026-11-08-a. DOI [10.5281/zenodo.24054760](https://doi.org/10.5281/zenodo.24054760) remains the separate binding-energy revision b.

The deposited archive does not name this DOI. `handoff.json` binds the zip and `MANIFEST.json` hashes and was written as delivery metadata. Those bytes stay unchanged. This note records the DOI.

Numerical cross-checks with Academia Co-scientist (AI-assisted); not peer review.

## Scope

The pilot is an action-derived static quadrupole for a coupled metric, scalar, and fluid, with exponential coupling to the bare energy density. The GR Love-number formula is not placed on a modified stellar background. The applied tide is tensor-only: there is no independent applied growing scalar mode, and the induced exterior scalar response is kept.

SLy and APR4 are compared at the same distant Kepler mass M_K = 1.4 solar masses. For a = 1 the laboratory-G calibration assumes unscreened dust and uses G_*/G_lab = 0.9262883177508389. The modified ADM energy mass in the archived table is 1.511408460164327 solar masses, distinct from that Kepler mass. Ranges are L = 1 km and 10 km.

Λ_metric is the dimensionless metric response of this calculation. It is not automatically the tidal parameter inferred by a GR gravitational-wave likelihood. The shifts belong to the calibrated model as a whole, not to an isolated extra scalar force. Two equations of state do not separate modified gravity from EOS uncertainty. The result is not an observation, a stability conclusion, a priority claim, or an empirical validation of the model.

The deposited record describes r4 as output guards, comparison alignment, and a notation correction. It states that the solver, EOS, and scientific table are unchanged relative to r3, that 70 regression tests passed (6 accepting cases and 64 rejecting cases, 35 of them under `python -O`), that six principal states were reproduced without stored backgrounds, and that the 18-script `code/run_all.py` suite was not rerun for r4. The optional a = 10^-8 APR4, L = 10 km background is a rejected node-limit case and has no Λ. Agreement with the historical bridge source SHA-256 `8e1f370db232d71d8813054f011c84db47ab6d1d2acfbcb03f5fc0c3e9533f32` is only the attributed Academia Co-scientist report in the deposited notes. That source file is not in this archive, and this import did not rerun it.

## Files

| file | role |
|---|---|
| `Yukawa_tidal_pilot_20261108_r4.zip` | Deposited archive, 2440905 bytes. SHA-256 `8001ccbfbb304db1b8506001b2b3b835b57543496f5b1699f42511a0cce720d4`. No PDF is inside it. |
| `handoff.json` | Deposited handoff. Zenodo stores the same bytes as `handoff (1).json`. |
| `LICENSE.txt` | License notice copied from `Yukawa_tidal_pilot_20261108_r4/LICENSE.txt` inside the archive. |
| `SHA256SUMS` | SHA-256 manifest of these imported files. |

`MANIFEST.json` inside the archive lists 160 payload files and does not hash itself. Its SHA-256 is `63020aa16ef67e7efa96a4b1a97fe4a1e2fbd389acd6ccd9c974fb57506e3caa`. The handoff names `delivery/archive_verification.json` and `internal/r4_audit_response_HU.md`. Those two files are not in the zip and are not separate Zenodo files, so they are not added here.

Verification: `cd YukawaTidalPilot_2026-11-08_r4 && sha256sum -c SHA256SUMS`.

## Safe extraction and reproduction

Reject members whose names are absolute or contain `..` before writing them. The archive has no such members. `ZipFile.testzip()` reports no bad CRC. After a safe extract, check every `MANIFEST.json` payload hash. From the package root, the deposited instructions are:

```sh
OPENBLAS_NUM_THREADS=1 python3 code/run_all.py
OPENBLAS_NUM_THREADS=1 python3 code/reproduce_pilot.py
python3 code/build_outputs.py
python3 code/test_publication_guards.py
```

`run_all.py` is the 18-script physics suite. `reproduce_pilot.py` rebuilds the six baseline states. `build_outputs.py` rewrites figures and status text from the tables. `test_publication_guards.py` is the postprocessing regression and rewrites `checks/r4_publication_regressions.json` in that tree. This import ran only the guard suite, on a scratch copy, and did not run `run_all.py`, `reproduce_pilot.py`, or `build_outputs.py` against the deposited tree. The package records Python 3.13 and `requirements.txt`.

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International · https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission.

## Magyar

Statikus tenzori árapályválasz lokális tömeges skalármezős neutroncsillagmodellben: rögzített laboratóriumi G melletti SLy/APR4 pilot (r4). Zenodo DOI: [10.5281/zenodo.24057132](https://doi.org/10.5281/zenodo.24057132). Verzióindex: 38.

Ez nem a kötésienergia-revízió. A hatásból származó csatolt metrika–skalár–folyadék sztatikus kvadrupólusát számolja, csupasz energiasűrűség exponenciális csatolásával. A GR Love-képlet nem kerül módosított háttérre. Közös távoli Kepler-tömeg: M_K = 1,4 naptömeg. a = 1 mellett G_*/G_lab = 0,9262883177508389, árnyékolatlan porfeltevés mellett. Az ADM-energiatömeg ettől külön érték. Λ_metric nem automatikusan a GR hullámforma-következtetés Λ-ja. Két EOS nem választja szét a módosított gravitációt az EOS-bizonytalanságtól. Nincs észlelési, stabilitási, elsőbbségi vagy modellvalidációs állítás.

Numerikus keresztellenőrzés az Academia Co-scientist segítségével (MI-alapú); nem szakmai lektorálás.

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International · https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission.
