# Step 27 (part of V28): re-derivation method, ocean and BREN

Series: Time Attraction Law — Scalar-Tensor Yukawa Fit, V28 (no V28 DOI yet; earlier DOIs: V24 10.5281/zenodo.23680715, V25 10.5281/zenodo.23699685, V26 10.5281/zenodo.23727585, V27 10.5281/zenodo.23761574).

This document shows how a published gravity residual is turned into an α limit at λ = 10 km for the two re-derivations a reviewer asked about: the ocean experiment (Zumberge et al. 1991, inputs via March et al. 2019) and the BREN tower (LLNL report OSTI 6190822). Every α_lim here is **our own, unverified re-derivation**. Every R built on it is **conditional**.

Convention: V(r) = −G_∞ m₁m₂/r · (1 + α e^{−r/λ}), so a laboratory measurement (r ≪ λ) gives G_lab = G_N = G_∞(1+α), i.e. **G_∞ = G_N/(1+α)**. Model point: α = 0.0795775, λ = 10 km, α > 0, universal, unscreened.

Code map (`run_step27.py`):

| Step | Function | What it does |
|---|---|---|
| 1 | `yukawa_slab_field_infinite_layers` | Yukawa field of horizontal layers (volume integral) |
| 2 | `ocean_model_residual_slope` | Model residual with the experiment's definition, and its fitted slope |
| 3 | `ocean_alpha_lim_from_residual` | α where the model slope reaches the edge of the measured interval (exact inversion edge/(S − edge)) |
| 4 | `march_eq81_F_difference` | Cross-check of step 1 against March eq. (81) and a brute-force integral |
| 5 | `tower_model_residual_uGal` | Half-space tower residual (observed − continued prediction) |
| 6 | `tower_alpha_lim_from_residual` | α where the model residual reaches μ − kσ (exact root E/(K − E)) |
| 7 | `cruz_conditional_scaling_factor` | Conditional 1 km → 10 km scaling for Erie (not a bound) |
| 8 | `back_substitution_tests` | Checks that every limit point returns the allowed residual edge (pass/fail in JSON) |

G handling in every α search (review r1, S01): the measured G_lab is held fixed, so each trial α uses G_∞(α) = G_lab/(1+α); α → 0 sensitivities use G_lab.

All intermediate numbers are written to `step27_results.json` → `review_trace`.

## A. Ocean: Zumberge et al. 1991, via March et al. 2019

**Data status: secondary-source data (not directly verified primary data).** The Zumberge paper itself (PRL 67, 3051) was **not obtained**. Every input below is as quoted by March et al., Phys. Rev. D 100, 042002 (2019), arXiv:1904.12789v2. Page numbers refer to that PDF.

### A1. Residual definition (what is subtracted, which G)
- March p. 11, eq. (72): corrected gravity differences g_c(z_N) − g_c(z_N,1). Raw gravimeter readings are corrected for instrument drift and for (i) temporal effects (Eötvös effect, Earth and ocean tides, vertical acceleration of the submersible) and (ii) gravity disturbances δg (isostatically compensated local seafloor topography, regional density inhomogeneities beneath the ocean). March refers to Zumberge for the size and uncertainty of these corrections.
- March p. 11, eq. (73): **Δg_obs(z_N) = g_c(z_N) − g_c(z_N,1) − β_N (z_N − z_N,1)**. This is the excess over the Newtonian value computed with **G = G_N (laboratory G)**. β_N is the Newtonian vertical gradient inside seawater (March p. 9: β_N = 2.226×10⁻⁶ s⁻² in oceanographic use; p. 10, eq. (67): c₀ = β_N + 4πG_Nρ_w).
- Profiles cover z_N,1 = 500 m to z_N,2 = 4800 m (March p. 11).

### A2. Quoted value and uncertainty
- March p. 11, eq. (74): the fit to the average slope of Δg_obs over the four profiles is **dΔg_obs/dz_N = −0.060 ± 0.172 mGal/km**.
- **The statistical meaning of ±0.172 is not stated in March, and we could not check Zumberge.** March eqs. (75)–(76) use it as a hard interval, −0.232 ≤ slope ≤ +0.112 mGal/km (±1 × 0.172).
- **Our assumption:** ±0.172 is a 1σ Gaussian uncertainty.

### A3. Geophysical background
- March p. 7, eq. (46): water 0–z_w with z_w = 5000 m; oceanic crust z_w–z_c with z_c − z_w = 7000 m and ρ_c = 2860 kg/m³; mantle below z_c with ρ_m = 3400 kg/m³.
- March p. 12, eq. (79) and text: seawater density 1023.6 → 1050.5 kg/m³; linear-profile mean ρ̄_w = 1037.05 kg/m³.
- Layers are infinite horizontal slabs. March's lateral terms ΔY and ΔG are neglected (p. 12, before eq. 81). The 36 m sediment layer is negligible (p. 7).
- Terrain, isostatic and regional effects are handled only through Zumberge's corrections inside g_c (A1). We add no further terrain model.
- Free air: the free-air part comes from GM_⊕, which is fixed by orbits. It is identical in data and reference and cancels.
- Regional correction: it comes from the source's g_c/δg preprocessing, which is not reconstructed here. The linear fit estimates the remaining slope and does not separate a remaining linear regional error from the Yukawa signal.

### A4. Formula: residual → α
Model residual, the same definition as eq. (73) (March eqs. 63–64 give the same structure):

  Δg_m(z) = −4π (G_∞ − G_N) ρ̄_w (z − z₁)  +  2π G_∞ α [F(z) − F(z₁)],   with G_∞ = G_N/(1+α)

  F(z) = ∫₀^∞ ρ(z′) sign(z′ − z) e^{−|z′−z|/λ} dz′   (volume integral over infinite slabs; March eqs. 47 and 66, F = 𝓕)

- The first term is the G renormalisation. The data reference uses G_N, but the true Newtonian Bouguer-type gradient uses G_∞. For small α it equals +4πG_N ρ̄_w α (z − z₁).
- The model slope s(α) is the least-squares slope of Δg_m over z ∈ [500, 4800] m.
- With G_lab fixed, both terms scale exactly as α/(1+α): **s(α) = S·α/(1+α)**.
- α → 0 sensitivity: S = **123.108 mGal/km per unit α**. Endpoint chords over 0–5000 m, computed with G_lab: **86.97917** (G renormalisation, Bouguer-type) and **34.78485** (Yukawa volume integral). These are chord values, not direct components of the fitted S. (The previous version gave 32.2 for the Yukawa chord because G_∞ was frozen at q.)
- Model point: s(q) = **+9.07 mGal/km** (= S·q/(1+q)).
- Edge: s_obs + k·σ with k = 1.95996 gives 0.27711 mGal/km.
- **Exact inversion (primary): α_lim = edge/(S − edge) = 0.00225606**, R = 35.27 (34.6 with the rounded 0.0023).
- Linear approximation, reported only as such: α_lim ≈ edge/S = 0.00225098, R = 35.35 (0.23% difference; it does not change the rounded ≈0.0023, ≈35). **R is conditional.**

### A5. Statistical meaning of α_lim
With σ = 0.172 read as 1σ Gaussian and k = 1.96, the edge is the upper end of the central 95% interval of the slope. For α > 0 (which predicts a positive slope), that makes α_lim a **one-sided 97.5% upper limit**. A one-sided 95% limit (k = 1.645) would give 0.00181. These numbers mean what they say only if the uncertainty is Gaussian, 1σ, and dominated by statistics. None of that is verified.

### A6. Our assumptions that the source does not make
1. Pure, universal, unscreened Yukawa. March's NMC "extra force" (pressure-to-depth term) is dropped.
2. ±0.172 read as 1σ Gaussian (March uses it as a ±1× interval).
3. Mantle half-space extends to infinity. March p. 7 instead assumes the Yukawa force is suppressed below the upper mantle, "to be verified a posteriori".
4. The slope is fitted to our model over 500–4800 m on a uniform grid. Zumberge's fit used the measured profiles.
5. Depth is used directly as z. March p. 10 notes a pressure-to-depth conversion effect of order cm at α ~ 10⁻³. At α ≈ 0.08 this was not re-evaluated.
6. λ = 10 km lies outside the 1–5000 m range of Zumberge's published bound. Our extension relies only on the layered model above.

### A7. Cross-check (from the code)
- Our layer integral F(z_w) − F(0) = 4.147×10⁶ kg/m². An independent brute-force numerical integral gives 4.148×10⁶ (agreement better than 10⁻⁴).
- **March eq. (81) as printed** (p. 12) gives −9.006×10⁵ kg/m². The difference, 5.048×10⁶ kg/m², is exactly λ(1 − e^{−z_w/λ})(2ρ_c − ρ_w − ρ_m). Our direct layer integral agrees with March's input profile, eqs. (46)–(47); the printed reduced expressions (48)/(50)/(81) do not give the same algebraic layer result. The discrepancy **cannot be resolved by a sign or G-convention change**; the authors' intent or any correction was **not verified**.
- **Separate scenario** (not a statistical error band): if eq. (81) as printed were used, the endpoint-chord sensitivity is 86.98 − 7.55 = 79.4 mGal/km per α, giving α_lim = 0.00350 (exact inversion; linear 0.00349) and R = 22.7.
- Using March's own ±1× interval, the edge is 0.112, giving α_lim = 0.00091 (exact) and R = 87.4.
- Plausibility check at λ ≲ 100 m: the model gives α_lim ≈ 0.0032 at the ×1.96 edge, of the same order as Zumberge's published 0.002 (quoted by March p. 6).
- Taken together, the ocean α_lim at 10 km spans **≈0.0009–0.0035, R ≈ 23–87**, depending on choices we cannot verify. These are **scenarios**, not one combined confidence interval; R is never a σ significance.

### A8. Open questions (not reconstructed)
- The original uncertainty of the slope, the weighting of the measured profiles and the model dependence of the pressure-derived depth are not reconstructed.
- The analytic layer check only verifies the idealized response, not the experiment's actual response.

## B. BREN tower: LLNL report OSTI 6190822

**Data status: primary document, preliminary results** (Thomas, Kasameyer et al., "A test of Newton's Law of Gravity using the BREN Tower, Nevada", LLNL conference report, https://www.osti.gov/biblio/6190822; 16 pp.). The journal version, Thomas et al., PRL 63, 1902 (1989), was available as abstract only.

### B1. Residual definition
- Report pp. 6–7 (eqs. 4–7): the surface gravity anomaly Δg is upward-continued with Poisson's integral, in the spherical approximation with all surface data treated as if at one elevation. The prediction on the tower is compared with the observed Δg(r) differences. Centrifugal and tidal effects are removed from all data (p. 7).
- Residual = observed − predicted change of Δg with height. **No laboratory G enters the leading term**: the prediction is a harmonic continuation of measured gravity.

### B2. Quoted value and uncertainty
- Report p. 11 (section 4): "Our predictions are within 93 µGal of observations at all heights on the tower, and the uncertainties at each height are as large as 95 µGal. The error bars in Figure 2 indicate the single standard deviation uncertainty from sampling errors."
- The abstract (p. 3) says "consistent with the Newtonian hypothesis to within 93 ± 95 µGal at the top of the tower".
- **The sign of the 93 µGal is not stated.** We use the conservative **μ = −93 µGal** (the same sign as the α > 0 signal) with **σ = 95 µGal (1σ, sampling only)** at the top, h = 454 m.

### B3. Geophysical background
- Report p. 11: 91 stations within 2.5 km plus 60,000 USGS/NOAA points out to 300 km. Sector means and a "regional" planar trend are removed per ring before the uncertainty is estimated (p. 13).
- Report p. 12: the spherical approximation is "potentially our most significant source of systematic error". A terrain model (density 1.8 g/cm³) shows that it **underestimates Δg(r₁) − Δg(r₀) by about 100 µGal**. That systematic is **not** included in the 95 µGal.
- Our model ignores terrain entirely (flat half-space).

### B4. Formula: residual → α
Half-space of uniform density ρ (assumed 2000 or 2670 kg/m³; not from the report). The Yukawa part of the surface field, 2πG_∞αρλ, is continued as a constant, while the true field at height h is 2πG_∞αρλe^{−h/λ}:

  Δg(h) = −2π G_∞(α) α ρ λ (1 − e^{−h/λ}) = −K α/(1+α),   G_∞(α) = G_lab/(1+α),   K = 2π G_lab ρ λ (1 − e^{−h/λ})

G_lab is held fixed for every trial α (review r1, S01). At h = 454 m and λ = 10 km, K = 37,226 µGal (ρ = 2000) or 49,697 µGal (ρ = 2670); this is the α → 0 sensitivity. (At α = 1 the residual is −K/2, i.e. computed with G_N/2.)

With k = 1.96 the edge is μ − kσ = −279.2 µGal, E = 279.2 µGal, and the exact root is **α_lim = E/(K − E) = 0.00756 (ρ = 2000) or 0.00565 (ρ = 2670), so R ≈ 10.5–14.1 (conditional)**. The previous version froze G_∞ at q and gave 0.0081/0.0061. If the sign were +93 µGal, the result would be 0.0025 (ρ = 2000) or 0.0019 (ρ = 2670), R ≈ 32–42.

### B5. Statistical meaning
One-sided 97.5% upper limit, but only if σ is Gaussian and 1σ and the systematics are negligible. The report's own ~100 µGal spherical-approximation systematic contradicts the last condition.

### B6. Our assumptions that the source does not make
1. Pure universal Yukawa; flat half-space; uniform assumed density; terrain ignored.
2. The sign of the 93 µGal residual (conservative choice).
3. Only the top-of-tower point is used. The full height profile of Figure 2 is not digitized.
4. The ~100 µGal systematic is not added. It needs a separate nuisance treatment and cannot be added automatically as an independent Gaussian σ to the 95 µGal. Illustration only: a bounded ±100 µGal nuisance at its worst-case edge would move the ρ = 2000 edge from 0.00756 to 0.01029 (ρ = 2670: 0.00565 → 0.00769); this is not a new data fit and not an empirical bound.
5. The surface Yukawa field is treated as uniform over the survey region, even though the data extend to 300 km, which is ≫ λ.

## C. Known gaps
- Zumberge 1991 not obtained, so the meaning of ±0.172, the corrections and their uncertainties are unchecked.
- The disagreement between our layer integral and March eq. (81) as printed is unresolved.
- No re-derivation includes the experiment's own response to a Yukawa field in the real geometry (finite survey area, terrain, how the reference model is built). This is the "measurement response" that must be verified before any exclusion can be claimed.
- Towers: the actual surface sampling, terrain, depth-dependent density and the joint Yukawa response of the Poisson continuation (report eq. 7) are not computed. Tower densities are assumed. The towers' published per-height residuals are not used.
- Every limit point passes the back-substitution test (it returns the allowed residual edge to a relative deviation ≤ 1×10⁻¹⁰); this checks the inversion only, not the physical model.
- Jekeli et al. 1990: the abstract value α < 0.001 has an unknown λ range. It is a text-only reference and is not applied (not plotted, not in the 10 km R summary).
- The Erie (Cruz 1991) 1 km bound cannot be transferred to 10 km without that experiment's response function, so 0.00088 is only a conditional estimate (`cruz_conditional_scaling_factor`).
- Missing primary sources (listed separately; Zumberge values used here are secondary-source data via March et al. 2019): Zumberge et al. 1991; full texts of Cruz et al. 1991, Jekeli et al. 1990, Thomas et al. 1989, Romaides et al. 1994/1997; not reached: Stacey et al. 1987, Moore et al. 1988, Ander et al. 1989, Bartlett & Tew 1989, Achilli et al. 1997, Rapp 1987, Fischbach & Talmadge 1999 (book).
