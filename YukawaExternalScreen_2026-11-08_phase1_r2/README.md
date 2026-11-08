# Weak-field forces and external constraints in a local massive-scalar model: finite-source geometry and laboratory-G calibration (phase 1, revision r2)

Phase 1, revision r2. Publication date: 2026-11-08. Author: Szabolcsi Dávid (ORCID [0009-0005-6300-074X](https://orcid.org/0009-0005-6300-074X)).

Zenodo DOI: [10.5281/zenodo.24060636](https://doi.org/10.5281/zenodo.24060636)  
Concept DOI: [10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681) · Zenodo version index 39.

This record is the weak-field external-screen package. It is not the binding-energy study and not the static tidal pilot. The deposited archive does not name this DOI. `handoff.json` binds the zip and `MANIFEST.json` hashes. Those bytes stay unchanged. This note records the DOI.

Numerical cross-checks with Academia Co-scientist (AI-assisted); not peer review.

## Scope

The calculation is the ε-source, universal unscreened dust branch. In the stated weak-field convention, q = a/(4π) and λ = L. Laboratory calibration is geometry-dependent: G_*/G_fit = 1/(1 + q κ_cal). The ideal short-distance ratio G_*/G_lab = 1/(1 + q) is not applied to every experiment. This branch is not the separate clock-based H1 hypothesis.

For a = 1 and L = 1 km, q ≈ 0.07957747. Under the standard universal Yukawa reading, that value is incompatible with the publisher-deposited Cruz (1991) abstract bound |α| < 0.001 at λ = 1000 m. The archived ratio to that bound is 79.57747154594767, about 79.58. That ratio is not a statistical significance and not a confidence-level exclusion. The abstract states no confidence level. The record’s necessary condition a < 0.01256637 is not a demonstrated allowed region.

For L = 10 km, a verified primary experimental response and confidence level remain unresolved. That benchmark is not classified as experimentally allowed. Neither the 1 km bound nor the synthetic finite-source rows replace a 10 km experimental constraint.

No experimental detection, empirical model validation, or exclusion of every TAL variant is claimed. GR–EOS mimicry and radial stellar stability are not calculated here.

Revision r2 changes the floating-point grid match and automatic creation of the delivery directory. The deposited preservation record says the scientific CSV tables, the two figures, `code/weak_field.py`, `code/reproduce.py`, and the derivation are byte-identical to the earlier phase-1 archive.

## Files

| file | role |
|---|---|
| `Yukawa_external_screen_20261108_phase1_r2.zip` | Deposited archive, 257885 bytes. SHA-256 `d6dadb52c60390eb26ffca5445cd158d12a406f14913df82b47f3c13338090e4`. No PDF is inside it. |
| `handoff.json` | Deposited handoff. Zenodo stores the same bytes as `handoff (3).json`. |
| `LICENSE.txt` | License notice copied from `LICENSE.txt` inside the archive. |
| `SHA256SUMS` | SHA-256 manifest of these imported files. |

`MANIFEST.json` lists 35 payload files and does not hash itself. Its SHA-256 is `2b165fe4d75dd34a89704f7e330a60a733787e3980c4c74b5ba36815e3ce436e`.

Verification: `cd YukawaExternalScreen_2026-11-08_phase1_r2 && sha256sum -c SHA256SUMS`.

## Safe extraction and checks

Reject members whose names are absolute or contain `..`. This archive has none. `ZipFile.testzip()` reports no bad CRC. After a safe extract, `python3 code/verify_manifest.py` checks the 35 payload hashes. From the package root, the deposited commands are:

```sh
python3 code/verify_manifest.py
python3 code/reproduce.py
python3 code/data_gate.py
python3 code/regress_guard.py
python3 code/test_grid_roundoff.py
python3 code/test_archive_workflow.py path/to/Yukawa_external_screen_20261108_phase1_r2.zip delivery/workflow.json
```

`reproduce.py` rewrites `data/`, `figures/`, and `checks/validation.json` in that tree. `regress_guard.py` and `test_grid_roundoff.py` rewrite their own check JSON files and create `delivery/`. Run them on a copy if the deposited bytes must stay unchanged. The archived environment is Python 3.13.5 with numpy 2.5.4, scipy 1.18.1, and matplotlib 3.11.2.

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International · https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission.

## Magyar

Gyenge téri erők és külső korlátok lokális tömeges skalármezős modellben: véges forrásgeometria és laboratóriumi G kalibráció (1. fázis, r2). Zenodo DOI: [10.5281/zenodo.24060636](https://doi.org/10.5281/zenodo.24060636). Verzióindex: 39.

Ez az ε-forrású, univerzális, árnyékolatlan dust-ág. q = a/(4π), λ = L. A G_*/G_fit = 1/(1 + q κ_cal) geometriafüggő. Az a = 1, L = 1 km pont a Cruz (1991) absztrakt |α| < 0,001 korlátjával a megadott univerzális értelmezésben összeférhetetlen. A mintegy 79,58-szoros arány nem szignifikancia és nem konfidenciaszintű kizárás. Az L = 10 km leképezés rendezetlen, és nem minősül kísérletileg megengedettnek. Nincs kísérleti észlelés, GR-utánzó EOS-keresés vagy radiális stabilitásszámítás. Az órás H1-ág külön hipotézis.

Numerikus keresztellenőrzés az Academia Co-scientist segítségével (MI-alapú); nem szakmai lektorálás.

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International · https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission.
