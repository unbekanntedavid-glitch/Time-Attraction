# Simulator data generation (`src/build_data.py`)

The page `index.html` runs as shipped: all data are already embedded (`src/data.js`, inlined by `src/build.py`). Regenerating `src/data.js` from the underlying inputs is a **separate step that requires external V26 inputs**, which are not part of the V28 package.

## Inputs

| Input | Location | Inside V28? |
|---|---|---|
| `tal_cosmology_common.py`, `run_step1_2_and_1_3.py`, `run_step1_4_to_1_6.py`, `tal_hubble_desi_bao_summary.json` (and the DESI mean/cov files read by `tal_cosmology_common.py`) | V27 Step 26 tree; in V28: `../step26_cosmology_desi_bao/` (default). Also: V27, doi:[10.5281/zenodo.23761574](https://doi.org/10.5281/zenodo.23761574); GitHub `V27/`: https://github.com/unbekanntedavid-glitch/Time-Attraction/tree/main/V27 | yes |
| `gw170817_and_lensing_summary.json`, `lensing_deflection_and_ppn_slip.csv`, `gw170817_dispersion_and_cT.csv` (V26 Step 23 outputs) | root of the V26 tree: V26, doi:[10.5281/zenodo.23727585](https://doi.org/10.5281/zenodo.23727585) (folder `scalar-tensor-yukawa-fit/` in the V26 ZIP), or GitHub `V26/`: https://github.com/unbekanntedavid-glitch/Time-Attraction/tree/main/V26 | **no (external)** |
| SPARC rotation curves `*_rotmod.dat` (175 files) | `<V26>/data/Rotmod_LTG/` (same V26 sources). Original data: Lelli, McGaugh & Schombert 2016, AJ 152, 157, doi:[10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157); https://astroweb.case.edu/SPARC/ | **no (external)** |

## Command

From inside the extracted `TAL_V28/` folder:

```
python3 simulator/src/build_data.py --v26 /path/to/V26 --out /tmp/data.js
cmp /tmp/data.js simulator/src/data.js     # identical bytes expected
python3 simulator/src/build.py              # inlines src/* into simulator/index.html
node simulator/verify_node.js               # 117 recomputed quantities vs the published JSON
```

Options (or environment variables): `--v27` / `TAL_V27_DIR` (default `../../step26_cosmology_desi_bao` relative to `src/`), `--v26` / `TAL_V26_DIR` (required), `--sparc` / `TAL_SPARC_DIR` (default `<V26>/data/Rotmod_LTG`), `--out` (default `src/data.js`). Missing inputs stop the script with a list of the missing files. Inputs are only read; nothing is written into the input trees.

Test (2026-10-30): from an extracted `TAL_V28/` tree with `--v26` pointing to the GitHub `V26/` folder, the regenerated `data.js` (150,822 B, 175 galaxies) is byte-identical to the shipped `src/data.js`.

Note: `src/data.js` (and so `index.html`) embeds the SPARC rotation-curve tables and the DESI DR1 BAO mean/covariance values as numeric inputs; their original sources are cited above and in the Step 26 reports.

Data sources and licenses (DESI CC BY 4.0, SPARC: no explicit license stated, RSD, Planck, BBN), with verbatim citation and acknowledgement texts: `../DATA_LICENSES.md`. The simulator footer states the same briefly in HU/EN.
