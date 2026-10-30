# Time Attraction Law — Scalar-Tensor Yukawa Fit, V28

**Author:** Dávid Szabolcsi (ORCID: 0009-0005-6300-074X)
**Date:** 2026-10-30 (revision r3)
**Status:** preliminary research package. The Step 27 conclusion is conditional; it is not a verified empirical exclusion. The functional checks (checksums, reruns, the 117 simulator comparisons) are not a physics audit of every simulator panel.

## Extraction (two ZIP files)

Extract both ZIP files (`TAL_V28_r3_part1_step26.zip`, `TAL_V28_r3_part2_simulator_step27.zip`) into the same directory; they merge into one `TAL_V28/` folder. The main `SHA256SUMS.txt` is valid only when both parts are extracted together; each part can be checked alone with `SHA256SUMS_part1.txt` / `SHA256SUMS_part2.txt`. From inside `TAL_V28/`: `sha256sum -c SHA256SUMS.txt`. Details: `README_EXTRACT_en.md`.

**DOIs**

- V28: assigned by Zenodo on publication (not yet available).
- V28 is a new version of V27: `10.5281/zenodo.23761574` (strictly V27, Step 26).
- V26: `10.5281/zenodo.23727585` (strictly V26, steps 1–25).
- V25: `10.5281/zenodo.23699685`; V24: `10.5281/zenodo.23680715`.

## Model

Scalar-tensor Yukawa branch with L = λ = 10 km, α_Y = 1, β = 1/(4π), adiabatic attractor C = 0. The effective local coupling is α_eff = q = α_Y²β = 0.0795775, positive, universal and unscreened. G dictionary: G_lab ≈ G_∞(1 + q), G_eff/G_lab (r ≫ L) = 0.9262883.

V28 contains three linked studies of this model.

## 1. Step 26 — cosmology audit against DESI DR1 BAO (`step26_cosmology_desi_bao/`)

The V27 Step 26 audit, carried over byte-identical from the published V27 release (GitHub `V27/`, 48 files, its `SHA256SUMS.txt` verifies 46/46). The Academia versions of the HU/EN reports are in `docs/`.

Question: is the rigid TAL background with H0 = 72.8 km/s/Mpc compatible with the DESI DR1 BAO distances (12×12 covariance)?

- Rigid background, H0 = 72.8: χ²_BAO = 194.342.
- M4 shape with its own re-integrated horizon (r_s = 144.443 Mpc, r_d = 147.092 Mpc): χ²_BAO = 20.785.
- Two-parameter profile (Ω_m, P = H0·r_d): χ²_BAO = 12.7405; at H0 = 72.8 the data require r_d = 140.02 ± 1.74 Mpc.
- M3 control with its own horizons: r_s = 144.366 Mpc, r_d = 147.012 Mpc, χ²_BAO = 2201.3.
- Joint BAO + Planck r_d prior + assumed local H0 prior (72.8 ± 1.0): χ² = 19.713, Δχ² = 6.973 against the BAO-only optimum.

Conclusion: the required r_d ≈ 140 Mpc lies well below the ≈147 Mpc horizons that the control models re-integrate for themselves. The rigid background is therefore not compatible with DESI BAO unless the early-universe sound horizon is modified. Reproduce with `python3 run_all.py` (see the Step 26 README).

## 2. Interactive simulator (`simulator/`)

`index.html` is a single, self-contained HU/EN page (open it in a browser, no server needed; 278,402 B, SHA-256 `46379f6323f7e9bdc342b95f0033d52c6c77fb5ffb1d466c676aabd69fae4226`). It differs from the V27 GitHub Pages build only by the point-source note below and the data-attribution line in the footer. Six panels: Yukawa force law and G, SPARC rotation curves, GW170817 and lensing, Cosmology V27 · DESI BAO, further probes, and info. Every number is computed in the browser from the package formulas. `src/` holds the sources and the build script (`build.py`), `shot_*.png` are screenshots of the panels (`shot_yukawa.png` HU and `shot_yukawa_en.png` EN show the r2 state; `shot_footer_hu.png` / `shot_footer_en.png` show the r3 footer data-attribution line).

In the Yukawa panel the Sun / Earth / neutron-star selection evaluates the scalar field with the point-source formula φ(r) = −α_Y β GM/(c² r) e^(−r/L) (`phiPoint`). The panel states this visibly: point-source, weak-field illustration; not the integrated field of the body's finite mass distribution.

Regenerating the embedded data (`src/build_data.py`) is a separate step that needs external V26 inputs (V26, doi:10.5281/zenodo.23727585); input paths are portable (`--v27`, `--v26`, `--sparc`). See `simulator/README_DATA_en.md`.

Reproducibility check: `node verify_node.js` recomputes 117 quantities with the page's own JavaScript and compares them with the published Step 26 JSON values (output: `verify_node_out.json`). 113 rows agree to a relative deviation below 10⁻⁶. 4 rows lie above 10⁻⁶ (largest 8.0×10⁻⁵, NGC 3198 C_req); the same 4 rows as in the reviewed V27 simulator state (rounded reference values, e.g. NGC 3198 χ² = 32639.3, and optimizer tolerances).

## 3. Step 27 — test of the fixed, unscreened local Yukawa branch (`step27_local_yukawa_branch/`)

Question: is the fixed coupling α = 0.0795775 at λ = 10 km compatible with geophysical (ocean and tower) gravity measurements?

R is the **ratio to bound** = 0.0795775 / α_lim(λ = 10 km). It is only meaningful under the same potential convention; R > 1 means the model lies outside the allowed region of the bound in question. R is **not** a significance in σ. Every α_lim below is our own re-derivation in idealized layer/half-space models, with the quoted errors read as Gaussian 1σ (×1.96 for 95%), and the G calibration held consistently fixed (G_∞ = G_lab/(1 + α), limits found by exact inversion). All R values are therefore **CONDITIONAL**.

| Record | α_lim (10 km) | R |
|---|---|---|
| Ocean (Zumberge 1991 via March 2019), exact inversion — decisive | 0.00225606 | 35.27 |
| Ocean, March ±1× interval | 0.000910599 | 87.39 |
| Ocean, March eq. (81) as printed (scenario) | 0.00350120 | 22.73 |
| Tower: Thomas 1989 (ρ = 2000 / 2670 kg/m³) | 0.00650259 / 0.00486292 | 12.24 / 16.36 |
| Tower: Romaides 1994 | 0.00228048 / 0.00170725 | 34.90 / 46.61 |
| Tower: Romaides 1997 | 0.00066375 / 0.00049711 | 119.89 / 160.08 |
| Tower: BREN (OSTI 6190822) | 0.00755665 / 0.00564970 | 10.53 / 14.09 |

- Applied R range: 10.53–160.08. The range collects scenarios (error readings, layer formulas, assumed densities); it is not a combined confidence interval.
- Back-substitution: 11/11 tests pass (each α_lim returns the edge of the allowed residual).
- Jekeli 1990: no 10 km value is applied (λ range unknown from the abstract). Cruz 1991 (Erie): only a conditional estimate (≈0.00088), not used in the range.

Verdict: our own re-derivation indicates strong tension with the fixed, unscreened Yukawa branch; confirming an exclusion requires verification of the measurement response and the uncertainty treatment. "Strong tension" holds only in this sense: in the idealized layer/half-space models, under the assumed error reading, a large numerical discrepancy results. It is not an empirical exclusion with a confidence level. The scope is only the fixed, unscreened local branch, not all TAL variants, and Step 26 is not reopened.

HU: „A saját újraszármaztatás erős feszültséget jelez a rögzített, árnyékolatlan Yukawa-ággal; a kizárás megerősítéséhez a mérési válasz és a bizonytalanságkezelés ellenőrzése szükséges.”

Data status:

- The Zumberge 1991 ocean inputs are secondary-source data taken from March et al. 2019 (PRD 100, 042002), not directly verified primary data.
- Missing primary sources: Zumberge et al. 1991 (PRL 67, 3051; not obtained), Cruz et al. 1991 (JGR 96, 20073; full text not obtained), Jekeli, Eckhardt, Romaides 1990 (PRL 64, 1204; full text not obtained). Further details: `step27_yukawa_fixed_branch_test/step27_results.json` (`missing_primary_sources`) and `sources/README_sources.md`.
- Third-party documents are not redistributed: `sources/README_sources.md` cites each with DOI, landing and download URL (checked) and the pages/equations used. Distribution rights of third-party documents not verified; documents are referenced by DOI/link instead of redistributed. The Adelberger 2003 Fig. 4 values in `reference_inputs/` are our own digitization, with attribution; the figure crop is not included.

Untested candidates (new branches, not verified escape routes): screening for Earth-scale sources; α ≲ 10⁻³ at 10 km; a different λ (e.g. λ < 1 cm) with unchanged α, to be re-tested.

Reproduce: `python3 step27_yukawa_fixed_branch_test/run_step27.py` (numpy, matplotlib) regenerates `step27_results.json` and the four figures. The bundle's own checksum list is `step27_local_yukawa_branch/SHA256SUMS.txt`. The review history (S01–S04) is in `step27_review1_changes_en.md`.

## Data sources and licenses / Adatforrások és licencek

The author's own content is CC BY-NC 4.0; embedded third-party data keep their own terms. Sources, verbatim citation and acknowledgement texts and the unclear points: `DATA_LICENSES.md`.

*HU: A csomag saját tartalmára a CC BY-NC 4.0 vonatkozik; a beágyazott külső adatok megtartják a saját feltételeiket. Részletek: `DATA_LICENSES.md`.*

| Forrás / Source | Licenc / License or terms | Hivatkozás / Citation | Újraterjesztés / Redistribution |
|---|---|---|---|
| DESI DR1 BAO mean + covariance (`step26_cosmology_desi_bao/desi_2024_*.txt`, `simulator`) | CC BY 4.0 — https://data.desi.lbl.gov/doc/acknowledgments/ | DR1 paper arXiv:2503.14745; DESI 2024 III/VI (arXiv:2404.03000, 2404.03002); DESI acknowledgement text (verbatim in `DATA_LICENSES.md` and in the simulator footer) | allowed with attribution and change notice / hivatkozással és a változtatások jelzésével megengedett |
| SPARC rotation curves (`simulator`) | no explicit license stated / nincs kifejezett licenc — https://astroweb.case.edu/SPARC/ | Lelli, McGaugh & Schombert 2016, AJ 152, 157, doi:10.3847/0004-6256/152/6/157 | not stated; unclear / nincs megadva; tisztázatlan |
| RSD fσ8 (12 values) | no data license stated / nincs adatlicenc | 8 source papers (list in `DATA_LICENSES.md`) | not stated / nincs megadva |
| Planck 2018 parameters | no explicit data license; acknowledgement requested — https://www.cosmos.esa.int/web/planck/planck-data-use | Planck 2018 VI, doi:10.1051/0004-6361/201833910 + ESA Planck acknowledgement | not stated / nincs megadva |
| BBN: PDG 2024 Y_p; Cooke+2018 D/H | PDG 2024 (PRD 110, 030001): CC BY 4.0; Cooke+2018: no data license stated | doi:10.1103/PhysRevD.110.030001; doi:10.3847/1538-4357/aaab53 | PDG: allowed with attribution; Cooke: not stated |
| `PANTHEON_GRID` | own synthetic grid, not third-party data / saját szintetikus rács | — | — |

Unclear: the official SPARC page states no license (it only asks for a citation); a 2025 Zenodo copy (doi:10.5281/zenodo.16284118) is labeled CC BY 4.0, but it is not confirmed to come from the authors.

## Changes in r3 (2026-10-30)

1. New `DATA_LICENSES.md` (HU/EN): license or stated terms of the embedded third-party data (DESI DR1 BAO, SPARC, RSD fσ8, Planck 2018, BBN), with checked URLs, verbatim citation and acknowledgement texts and redistribution status; where no explicit license exists, it says so.
2. New section in this README, "Data sources and licenses / Adatforrások és licencek"; pointer in `simulator/README_DATA_hu.md` / `_en.md`.
3. Simulator: visible HU/EN data-attribution line in the footer and a collapsible block with the verbatim DESI and Planck acknowledgements. `index.html` rebuilt (278,402 B, SHA-256 `46379f6323f7e9bdc342b95f0033d52c6c77fb5ffb1d466c676aabd69fae4226`); `verify_node.js` gives the same 117 rows (113 below 10⁻⁶, the same 4 above). No computation changed.
4. Step 26 and Step 27 files unchanged.

## Changes in r2 (2026-10-30)

1. Third-party PDFs (March 2019, OSTI 6190822, OSTI 7173840, Fischbach & Talmadge 1996, Adelberger et al. 2003) and the Adelberger figure crop removed; replaced by DOI/URL references in `step27_local_yukawa_branch/sources/README_sources.md`. Step 27 results unchanged (rerun gives an identical `step27_results.json`). The Step 26 part contains only the author's own PDFs.
2. `simulator/src/build_data.py`: portable input paths (`--v27`, `--v26`, `--sparc`, `--out`, or environment variables) instead of fixed local paths; data sources documented in `simulator/README_DATA_hu.md` / `_en.md`. Regeneration reproduces the shipped `data.js` byte for byte.
3. Simulator Yukawa panel: visible HU/EN note that the Sun / Earth / neutron-star values are a point-source, weak-field illustration. `index.html` rebuilt; `verify_node.js` gives the same 117 rows (113 below 10⁻⁶, the same 4 above).
4. Extraction guide (`README_EXTRACT_hu.md` / `_en.md`, in both ZIPs) and per-part manifests `SHA256SUMS_part1.txt` / `SHA256SUMS_part2.txt`.

## Integrity

`SHA256SUMS.txt` at the top level lists the SHA-256 of every file in `TAL_V28/` (both parts). Check: `sha256sum -c SHA256SUMS.txt`.
