# The TAL factor f_L = exp(−γ_L ρ̄) as a seventh model in the hash-locked v17 SPARC runner

**Author of the law:** David Szabolcsi (ORCID 0009-0005-6300-074X)
**Analysis date:** 2026-10-20
**Status:** supplementary analysis. Nothing in the v17 deposit, the repository or the manuscript was changed.
**Licence:** CC BY-NC 4.0, as for the v17 deposit

---

## 1. Purpose

The v17 SPARC audit (Zenodo [10.5281/zenodo.23479940](https://doi.org/10.5281/zenodo.23479940)) compares six purely baryonic baselines on a hash-locked blind sample of 123 galaxies (2596 radii). It uses one error metric for all of them. The v17 audit also says that applying the TAL factor itself was outside its scope (v17 PDF, §2 and §8).

This note adds the published TAL factor

  f_L(ρ̄_L) = exp(−γ_L ρ̄_L),  γ_L = G_N L² / c²

as a **seventh model**. It uses the same runner, the same input hashes, the same sample and the same σ and χ²_red definitions. It then compares that model with the six baselines using paired per-galaxy tests.

## 2. Data and integrity

| Item | Value |
|---|---|
| Rotation curves | SPARC `Rotmod_LTG.zip` (Lelli, McGaugh & Schombert 2016), SHA256 `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588` (verified) |
| Blind sample list | `blind_123.txt`, SHA256 `3637d2f5bda6b2b97e645c5fba5cf39d310a1826bee388286c2c74b9171aa283` (verified) |
| Sample | 123 galaxies, 2596 radius points (the runner's built-in check passed) |
| Runner | `sparc_full_runner.py` from `v17_package/github/`, SHA256 `eb40036842162dc9a42216bd7a4d35c7f2ce8addcf7e3ed8ffd3626427b22313`. An unmodified copy was used. |
| Error model | σ = max(0.05·V_obs, 20 km/s); χ² = sum over radii of [(V_pred − V_obs)/σ]²; χ²_red = χ²/(N − 1) |

**Baseline reproduction.** The six baselines were re-run from the copied runner. All 12 CSVs (`points_*.csv`, `summary_*.csv`) are **byte-for-byte identical** to `v17_package/github/run_output/` (checked with `cmp`). The medians are vbar 3.1433, vbar_095 3.1995, vbar_105 2.9555, no_bulge 3.3131, no_gas 4.0474, disk_only 4.2623.

**Report referred to as "SPARC_zarojelentes.pdf".** No file of that name exists anywhere on the analysis machine (searched with `find /`). Its content was taken from the "Akadémia zárójelentés (kivonat)" section of the v17 PDF (§8). That section states that the six baselines are reproducible and that the TAL factor was not part of the audit.

## 3. Formula and sources

### 3.1 Published elements used

| Element | Published statement | Source |
|---|---|---|
| TAL factor | f_L(ρ̄_L) = exp(−γ_L ρ̄_L), dτ_L = f_L dt, f_L ∈ (0, 1] | `Time_Attraction.pdf` (repository root; formalism of v1–v7, header DOI [10.5281/zenodo.23066993](https://doi.org/10.5281/zenodo.23066993)), Def. 3, Eq. (7); also `Time_Attraction_LabxPitch.pdf`; v17 PDF §1 (axiomatic form f(ρ) = exp(−γρ)) |
| Coupling constant | γ_L = G_N L²/c² [m³ kg⁻¹] | `Time_Attraction.pdf` Def. 2, Eq. (2); v17 PDF §2 |
| Averaged density | ρ̄_L(x) = (3/4πL³) ∫ over \|x′−x\| ≤ L of ρ(x′) d³x′ (top-hat); ρ̄_L → ρ as L → 0 | `Time_Attraction.pdf` Def. 1, Eq. (1) |
| Favoured scale | L₀ = 1.2×10⁻¹¹ m, γ₀ = 1.07×10⁻⁴⁹ m³ kg⁻¹ | `Time_Attraction.pdf` Eqs. (4)–(5) |
| Clock bound | γ < 10⁻²² m³ kg⁻¹, L < 367 m | `Time_Attraction.pdf` Eq. (3) |
| Matter coupling | G_μν = 8πG_N f_L T_μν, T_eff = f_L T; TOV: dm/dr = 4πr²ρ f_L | `Time_Attraction.pdf` Eqs. (11)–(12) and §8 |
| Galactic density proxy | "ρ_eff is the baryonic (disk + bulge) density proxy" | v17 PDF §2 |
| Disk density | ρ_disk(R,0) = Υ_disk·SB_disk(R)·M_⊙ pc⁻²/(2h_z), with h_z = 0.2 kpc, Υ_disk = 0.5, Υ_bulge = 0.7 | `Time-Attraction_v16.0.pdf`, constants line and Eq. (5) |
| Bulge density | "ρ_bulge proxy = 3 Υ_bulge V_bul²/(4πG r³)" (spherical enclosed-mass proxy) | `Time-Attraction_v16.0.pdf`, line following Eq. (6) |
| Other published scales | L₀ = 2.36 kpc (median R_disk of the 123-galaxy sample); "γ_L = 10²³ conflicts with the laboratory" | v17 PDF §3 and §7 |
| Author's own expectation | the rotation-curve correction at γ_lab is ~10⁻³⁴, which is negligible; L ~ 10²⁶ m would be needed for e^{γρ} ~ 5 | `Time_Attraction.pdf` §8, Eq. (14) |

### 3.2 Seventh model (primary, parameter-free)

  **V_pred²(R) = f_L(ρ̄(R)) · V_bar²(R),  f_L = exp(−γ₀ ρ̄),  γ₀ = G_N L₀²/c², L₀ = 1.2×10⁻¹¹ m**

  ρ̄(R) = Υ_disk·SB_disk(R)·M_⊙ pc⁻²/(2h_z) + 3Υ_bulge V_bul(R)²/(4πG R²)

Here V_bar = (V_gas² + V_disk² + V_bul²)^{1/2}, exactly as in the runner's `vbar` baseline. With these constants, γ₀ = 1.0694×10⁻⁴⁹ m³ kg⁻¹ (`arith_checks.py`).

### 3.3 Implementation choices (not specified in the published material)

1. **How f_L enters V(R)** (*implementation choice, not specified in the published material*). The published material gives no rotation-curve prediction rule for f_L. The v17 PDF mentions only the phenomenological rule V² = V_bar²·C, as distinct from f_L. The primary model takes the Newtonian limit of the published field equation T_eff = f_L T, which also matches the published TOV form dm/dr = 4πr²ρ f_L. That limit is an effective source density f_L ρ, so V² = f_L V_bar² for locally constant f_L. As a sensitivity check, the opposite coupling V² = V_bar²/f_L ("inverse", a time-rate rescaling reading) was also run. Because f_L ≤ 1, the primary coupling can only lower V_pred, and only the inverse coupling can raise it.
2. **ρ̄ from SPARC data** (*implementation choice*). The definition is a top-hat average over radius L. For the published L (≤ 367 m), the average reduces to the local density, by Eq. (1) in the limit L ≪ disk scale height. The local density is evaluated in the disk midplane at each radius, from the published v16 disk formula plus a bulge term. The gas term is omitted, following the v17 wording "disk + bulge". At 141 of the 2596 points, SB_disk = 0 and V_bul = 0, so ρ̄ = 0 there and f_L = 1.
3. **Bulge term dimension** (*implementation choice*). As printed in v16, the bulge proxy 3Υ_bulge V_bul²/(4πG r³) does not have units of density (it gives kg m⁻⁴). The dimensionally consistent mean enclosed density 3Υ_bulge V_bul²/(4πG r²) was used instead.
4. **Υ values** (*implementation choice*). Υ_disk = 0.5 and Υ_bulge = 0.7 are used **only inside ρ̄**, as published in v16. V_bar itself is kept identical to the runner's `vbar` baseline, which uses the tabulated SPARC components unscaled. This makes the seventh model a strict modification of `vbar`.
5. **Alternative density proxy** (sensitivity only). ρ_enc = 3V_bar²/(4πG R²), the mean enclosed baryonic density.
6. **Order-one scales.** L_unity = c/√(G ρ̄) is the L at which γ_L ρ̄ = 1. Its sample median is computed from the data (positive-ρ̄ points only). That makes it data-derived and **not blind**. The values 64.2 and 85.3 Mpc come from the author's earlier audit, with a different density proxy. They are re-run here with the present proxy.
7. **Fitted L** (**fitted, not blind**). One global L was chosen to minimise the median χ²_red over the same 123 galaxies. The grid ran from 0.1 to 1000 Mpc in 0.025-dex steps (`fit_scan_L.csv`). The optimum is in-sample and cannot count as a blind test.

## 4. Results

### 4.1 All seven models

| # | Model | V_pred | Median χ²_red | Mean χ²_red | Galaxies |
|---|---|---|---|---|---|
| 1 | vbar | V_bar | 3.1433 | 4.2161 | 123 |
| 2 | vbar_095 | 0.95 V_bar | 3.1995 | 4.3725 | 123 |
| 3 | vbar_105 | 1.05 V_bar | 2.9555 | 4.3224 | 123 |
| 4 | no_bulge | (V_gas² + V_disk²)^{1/2} | 3.3131 | 6.3928 | 123 |
| 5 | no_gas | (V_disk² + V_bul²)^{1/2} | 4.0474 | 5.2533 | 123 |
| 6 | disk_only | \|V_disk\| | 4.2623 | 7.5707 | 123 |
| **7** | **tal_fL (primary)** | **f_L^{1/2} V_bar, L₀ = 1.2×10⁻¹¹ m** | **3.1433** | **4.2161** | **123** |

For the primary model, max γ₀ρ̄ = 6.3×10⁻⁶⁷ and **max |f_L − 1| = 0** over all 2596 points. In IEEE double precision, exp(−6×10⁻⁶⁷) is exactly 1. All 2596 V_pred values are identical to `points_vbar.csv`. All 123 per-galaxy χ²_red values in `summary_tal_fL.csv` are identical to `summary_vbar.csv`; only the model-name column differs.

### 4.2 Paired tests: seventh model against each baseline

The tests compare per-galaxy χ²_red values (n = 123), and a lower value counts as a win. The Wilcoxon signed-rank test is two-sided and drops zero differences. The sign test is an exact binomial test with ties excluded.

| Seventh model vs | Wins (tal_fL) | Wins (baseline) | Ties | Wilcoxon p | Sign-test p | Median Δχ²_red (7th − baseline) |
|---|---|---|---|---|---|---|
| vbar | 0 | 0 | 123 | undefined (all differences 0) | undefined | 0 |
| vbar_095 | 97 | 26 | 0 | 2.57×10⁻⁴ | 7.97×10⁻¹¹ | −0.2455 |
| vbar_105 | 32 | 91 | 0 | 1.46×10⁻² | 9.74×10⁻⁸ | +0.1905 |
| no_bulge | 22 | 4 | 97 | 2.68×10⁻⁴ | 5.34×10⁻⁴ | 0 |
| no_gas | 110 | 13 | 0 | 9.22×10⁻¹⁹ | 2.61×10⁻²⁰ | −0.7127 |
| disk_only | 111 | 12 | 0 | 3.45×10⁻¹⁷ | 3.03×10⁻²¹ | −1.1301 |

The seventh model is identical to `vbar`. Every row therefore reproduces the vbar-vs-baseline comparison. In particular, vbar_105 beats it in 91 of 123 galaxies (Wilcoxon p = 0.015).

### 4.3 Sensitivity variants (compared with vbar)

| Variant | γ_L [m³ kg⁻¹] | Coupling | max \|f_L − 1\| | Median χ²_red | Mean χ²_red | Wins vs vbar (variant/vbar/ties) | Wilcoxon p | Sign p |
|---|---|---|---|---|---|---|---|---|
| L = 367 m (clock bound) | 1.0002×10⁻²² | forward | 0 | 3.1433 | 4.2161 | 0/0/123 | undefined | undefined |
| γ = 10⁻²² (clock bound) | 1.0000×10⁻²² | forward | 0 | 3.1433 | 4.2161 | 0/0/123 | undefined | undefined |
| L = 2.36 kpc (v17 L₀) | 3.938×10¹² | forward | 2.33×10⁻⁵ | 3.1433 | 4.2161 | 52/71/0 | 0.087 | 0.104 |
| L = 2.36 kpc (v17 L₀) | 3.938×10¹² | inverse | 2.33×10⁻⁵ | 3.1433 | 4.2161 | 71/52/0 | 0.087 | 0.104 |
| L_unity = 27.64 Mpc (data-derived, not blind) | 5.40×10²⁰ | forward | ≈1 | 7.2578 | 29.29 | 2/121/0 | 7.3×10⁻²² | 1.4×10⁻³³ |
| L_unity = 27.64 Mpc (data-derived, not blind) | 5.40×10²⁰ | inverse | ≈1 | 1935.2 | ∞ (overflow) | 24/99/0 | n/a | 5.2×10⁻¹² |
| ρ_enc, L_unity = 42.24 Mpc (not blind) | 1.26×10²¹ | forward | ≈1 | 8.4003 | 27.61 | 1/122/0 | 6.6×10⁻²² | 2.3×10⁻³⁵ |
| ρ_enc, L_unity = 42.24 Mpc (not blind) | 1.26×10²¹ | inverse | ≈1 | 7230.7 | ∞ (overflow) | 22/101/0 | n/a | 2.7×10⁻¹³ |
| L = 64.2 Mpc | 2.914×10²¹ | forward | ≈1 | 12.620 | 39.62 | 0/123/0 | 6.3×10⁻²² | 1.9×10⁻³⁷ |
| L = 64.2 Mpc | 2.914×10²¹ | inverse | ≈1 | 1.9×10¹⁹ | ∞ (overflow) | 7/116/0 | n/a | 1.4×10⁻²⁶ |
| L = 85.3 Mpc | 5.145×10²¹ | forward | ≈1 | 14.409 | 42.53 | 0/123/0 | 6.3×10⁻²² | 1.9×10⁻³⁷ |
| L = 85.3 Mpc | 5.145×10²¹ | inverse | ≈1 | 1.2×10³⁴ | ∞ (overflow) | 4/119/0 | n/a | 1.8×10⁻³⁰ |
| ρ_enc, L = 64.2 Mpc | 2.914×10²¹ | inverse | ≈1 | 1.14×10⁸ | ∞ (overflow) | 10/113/0 | n/a | 3.1×10⁻²³ |
| ρ_enc, L = 85.3 Mpc | 5.145×10²¹ | inverse | ≈1 | 5.37×10¹³ | ∞ (overflow) | 6/117/0 | n/a | 8.4×10⁻²⁸ |
| γ = 10²³ (v17 §7 value) | 1×10²³ | forward | ≈1 | 21.505 | 52.13 | 0/123/0 | 6.3×10⁻²² | 1.9×10⁻³⁷ |
| γ = 10²³ (v17 §7 value) | 1×10²³ | inverse | ≈1 | ∞ (overflow) | ∞ | 0/123/0 | n/a | 1.9×10⁻³⁷ |
| **FITTED L = 2.371 Mpc (not blind)** | 3.976×10¹⁸ | forward | ≈1 (max γρ̄ = 23.5) | **2.7083** | 4.1761 | 42/81/0 | 0.249 | 5.6×10⁻⁴ |
| **FITTED L = 0.473 Mpc (not blind)** | 1.583×10¹⁷ | inverse | 0.608 (max γρ̄ = 0.935) | **3.1425** | 4.3265 | 70/53/0 | 0.053 | 0.149 |

Notes on the table:
- "≈1" means that max |f_L − 1| rounds to 1, i.e. f_L → 0 at the densest points.
- "∞ (overflow)": exp(γρ̄/2) exceeds the double-precision range at some points, so χ²_red is infinite for at least one galaxy. The Wilcoxon test is then not computed (n/a). The sign test is still valid.
- The data-derived L_unity medians are 27.64 Mpc (5–95 %: 3.0–489 Mpc) for the primary proxy and 42.24 Mpc (3.8–204 Mpc) for ρ_enc.
- For the fitted rows, γ_L = G_N L²/c² with the optimal L from `run_meta.json` (2.3714 Mpc forward, 0.47315 Mpc inverse).

**Fitted forward variant against vbar_105:** 29 wins vs 94; Wilcoxon p = 3.0×10⁻³; sign p = 3.4×10⁻⁹.

## 5. Interpretation

1. **The published TAL factor is numerically identical to the pure-baryon baseline.** With the published favoured scale (L₀ = 1.2×10⁻¹¹ m), or with anything up to the published clock bound (L < 367 m, γ < 10⁻²²), γ_L ρ̄ is at most ~6×10⁻⁴⁰. At that size f_L equals 1 to machine precision. The seventh model returns the same V_pred at every radius, the same χ²_red in every galaxy and the same median (3.1433) as `vbar`. No paired test can separate them. It ranks exactly where `vbar` ranks: better than vbar_095, no_bulge, no_gas and disk_only, and worse than vbar_105. This agrees with the author's own published statement that the galactic correction is ~10⁻³⁴, i.e. negligible (`Time_Attraction.pdf` §8).
2. **Larger published scales also change nothing measurable.** At the v17 L₀ = 2.36 kpc, |f_L − 1| ≤ 2.3×10⁻⁵, and the median χ²_red changes only in the eighth significant digit.
3. **Order-one scales make the fit worse, not better.** At L ~ 28–85 Mpc, or at γ = 10²³ (the value v17 §7 itself says conflicts with the laboratory), f_L departs strongly from 1. The forward coupling then suppresses V_pred in the dense inner regions and loses to vbar in 121–123 of 123 galaxies. The inverse coupling inflates V_pred exponentially and diverges. Both directions are rejected at p < 10⁻¹¹.
4. **Fitted L is not evidence for the law.** With the forward coupling, one global L ≈ 2.37 Mpc tuned on the same 123 galaxies lowers the median χ²_red to 2.708. The per-galaxy comparison goes the other way: the variant is worse than vbar in 81 of 123 galaxies (sign p = 5.6×10⁻⁴), and worse than the simple 5 % rescaling vbar_105 in 94 of 123. The lower median comes from damping a minority of dense inner points. It is not a systematic improvement, and it was obtained in-sample. The best inverse-coupling L gives a median (3.1425) that matches vbar to three decimals, with no significant paired difference.
5. **The order-one scales assume a local density.** At L of tens of Mpc, the published top-hat definition, Eq. (1), averages over a sphere of that radius. Such a sphere is dominated by near-cosmic-mean density, not by the galaxy's local density. With the critical density quoted in `Time_Attraction.pdf` §8 (9×10⁻²⁷ kg m⁻³), γ_L ρ̄ is only 5×10⁻⁶ to 5×10⁻⁵ for L = 27.6–85.3 Mpc (`arith_checks.py`). Under a strict reading of the definition, f_L would again be ≈ 1 for galaxies. The Mpc variants above therefore overstate the possible effect.

**Conclusion.** Implemented as published, the TAL factor leaves the SPARC rotation-curve fit exactly at the pure-baryon level (median χ²_red = 3.1433 on the locked 123-galaxy sample). It neither improves nor worsens any of the six baryonic baselines. Scales large enough to give f_L a visible effect conflict with the published laboratory bound and worsen the per-galaxy fits.

## 6. Reproduction

Python 3.13 was used. The comparison step needs `scipy` (version 1.18.1 here) and `numpy`.

```bash
mkdir -p tal_7th_model_work/data && cd tal_7th_model_work
cp <repo>/v17_package/github/sparc_full_runner.py .
cp <repo>/v17_package/github/data/blind_123.txt data/
cp <path>/Rotmod_LTG.zip data/                 # SHA256 0a80cc90…f588
python3 sparc_full_runner.py                   # six baselines -> output/
for f in <repo>/v17_package/github/run_output/*.csv; do cmp "$f" output/$(basename "$f"); done
python3 tal_7th_model.py                        # 7th model + variants + L scan -> output_7th/
python3 compare_7th.py                          # paired tests -> output_7th/comparison.json
python3 arith_checks.py                         # gamma_L arithmetic quoted in §5
```

Work folder: `/workspace/artifacts/tal_7th_model_work/`. Logs: `run_log.txt`, `compare_log.txt`, `arith_log.txt`. The full hash list is in `SHA256SUMS.txt`.

| File | SHA256 |
|---|---|
| `tal_7th_model.py` | `73962aed3380230f91c898db2d1876cb7c4984c4fe01523bdcf6ef1727f572ef` |
| `compare_7th.py` | `af1c9501e775a724b60ad228753d164bacd71d4fdb3fc71eff9df095490882d7` |
| `arith_checks.py` | `c7de8b01f42a03f3659b4419ead7feb3172adc2a538040773b914f16c0bee7aa` |
| `sparc_full_runner.py` (unmodified copy) | `eb40036842162dc9a42216bd7a4d35c7f2ce8addcf7e3ed8ffd3626427b22313` |
| `output_7th/points_tal_fL.csv` | `9a4ed57b0fd0b078a462b2a55180b41d1d35e01fdd1328f47e36318d9d733800` |
| `output_7th/summary_tal_fL.csv` | `121cc97209f4e2f4a1aef97a43d9f1bb8803af6254c1ee23be62f3702f7ae724` |
| `output_7th/points_tal_fL_FITTED_forward.csv` | `2b4022fef37a705aa9c6a963ba80528b733667102d43ad0c8c352a0d243ec79d` |
| `output_7th/summary_tal_fL_FITTED_forward.csv` | `2a1927fb7359b1c560bcdc864c1c6254ca5e633b765eea9c90a66ea59375216e` |
| `output_7th/points_tal_fL_FITTED_inverse.csv` | `a99ad9207ec364771d2dfb7c633a848cb43ff8d3ccf37b2536fb89029b5909d7` |
| `output_7th/summary_tal_fL_FITTED_inverse.csv` | `07b4d765919b9ab7cfb505e73e0a8ee7a93c0e037a177c85c67626492a983db8` |
| `output_7th/comparison.json` | `795ff27185a7490a403e3979a1b9a2e8eb4df432e05b4882f3aee61107098ce2` |
| `output_7th/run_meta.json` | `2981af691b76de2c88f9cf6205a30de8789a813f6af4ccea12fe552c001f3768` |

## References

1. Lelli, F., McGaugh, S. S., Schombert, J. M. (2016). SPARC: Mass models for 175 disk galaxies with Spitzer photometry and accurate rotation curves. *AJ* 152, 157. [https://doi.org/10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157)
2. Szabolcsi, D. (2026). Time Attraction Law v17 – SPARC Audit. Zenodo. [https://doi.org/10.5281/zenodo.23479940](https://doi.org/10.5281/zenodo.23479940) (concept DOI [10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681))
3. Szabolcsi, D. The Law of Time Attraction: Density-Dependent Local Time Scale (v1–v7 formalism). Zenodo. [https://doi.org/10.5281/zenodo.23066993](https://doi.org/10.5281/zenodo.23066993)
4. Repository: [https://github.com/unbekanntedavid-glitch/Time-Attraction](https://github.com/unbekanntedavid-glitch/Time-Attraction) (`v17_package/`, `Time_Attraction.pdf`, `Time-Attraction_v16.0.pdf`, `Time_Attraction_LabxPitch.pdf`)
