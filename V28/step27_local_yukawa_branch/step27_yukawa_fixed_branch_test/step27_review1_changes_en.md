# Step 27 (part of V28): changes in review round 1

Basis: external review of `step27_V28_review_bundle.zip` (2026-10-29), points S01–S04 and open questions §5. No V28 DOI exists yet. Scope is unchanged: only the fixed, unscreened, universal local Yukawa branch (α = q = 0.0795775, λ = 10 km); this is not a re-opening of Step 26.

## S01: G calibration in the α-limit search (`run_step27.py`)
- New `g_inf_of(α) = G_lab/(1+α)`; every trial α uses it. The global `G_INF` (frozen at q) is now used only at α = q.
- Towers: Δg(α) = −Kα/(1+α), K = 2πG_lab ρλ(1 − e^{−h/λ})·10⁸ µGal (`tower_K_uGal`); exact root α_lim = E/(K − E), E = −(μ − 1.96σ).
- `tower_model_residual_uGal(1.0, …)` now uses G_N/2, consistent with the documented G(α).
- Ocean: exact inversion α_lim = edge/(S − edge) is primary; the linear edge/S is reported as an approximation. S is the α → 0 sensitivity with G_lab.
- March endpoint chords use G_lab: Yukawa 32.22081 → **34.78485** mGal/km; Bouguer **86.97917** (unchanged).
- New `back_substitution_tests()`: each limit point must return the allowed residual edge; result in JSON → `back_substitution_tests` (tolerance 1×10⁻⁹ relative). **All 11 tests pass**; maximum relative deviation 8.8×10⁻¹¹ (ocean, finite-difference S); towers ≤ 3.5×10⁻¹⁶.
- Model residuals at α = q are unchanged (difference 0.0 for every tower and density; ocean slope +9.07449725 mGal/km unchanged).

| Bound | Old α_lim | New α_lim (exact) | Old R | New R |
|---|--:|--:|--:|--:|
| ocean (primary) | 0.00225098 (linear) | **0.00225606** | 35.352 | **35.273** |
| ocean, linear approximation | — | 0.00225098 | — | 35.352 |
| ocean, March ±1× interval | 0.00090977 | 0.00091060 | 87.47 | 87.39 |
| ocean, March eq. (81) scenario | 0.0034647 | 0.0035012 (linear 0.0034890) | 22.97 | 22.73 |
| thomas1989 / 2000 | 0.00697470 | 0.00650259 | 11.409 | 12.238 |
| thomas1989 / 2670 | 0.00522449 | 0.00486292 | 15.232 | 16.364 |
| romaides1994_493m / 2000 | 0.00245635 | 0.00228048 | 32.397 | 34.895 |
| romaides1994_493m / 2670 | 0.00183996 | 0.00170725 | 43.250 | 46.612 |
| romaides1997_568m / 2000 | 0.00071609 | 0.00066375 | 111.127 | 119.891 |
| romaides1997_568m / 2670 | 0.00053640 | 0.00049711 | 148.355 | 160.081 |
| bren / 2000 | 0.00809680 | 0.00755665 | 9.828 | 10.531 |
| bren / 2670 | 0.00606502 | 0.00564970 | 13.121 | 14.085 |
| bren, sign + variant / 2000; 2670 | 0.00270; 0.00202 | 0.00251; 0.00188 | 29.4; 39.3 | 31.7; 42.4 |

The eq. (81) scenario changed because the March endpoint sensitivity now uses G_lab (S81 = 79.98 → 79.43 mGal/km) and the exact inversion. Overall applied R range: 9.83–148.4 → **10.53–160.1** (scenarios, not a confidence interval).

## S02: Jekeli 1990
- JSON: `alpha_lim_10km_published` = null; new fields `abstract_alpha_value` = 0.001, `lambda_status` = "unknown", `applied` = false; `ratio_to_bound_R` = null (was 79.58).
- Removed the 10 km plot point and the q/0.001 term from the 10 km R summary (`verdict.min_R/max_R`). Kept as a text-only, non-applied reference.
- Erie (Cruz 1991): 1 km value and its conditional 10 km estimate (0.00088) remain separate and outside the decision.

## S03: regional trend sentence (method A3, HU and EN)
- Old: "Regional trend: removed only through the linear slope fit."
- New: "Regional correction: it comes from the source's g_c/δg preprocessing, which is not reconstructed here. The linear fit estimates the remaining slope and does not separate a remaining linear regional error from the Yukawa signal." No similar claim found elsewhere.

## S04: figure labels
- Ocean band: "observed slope band (95%)" → "quoted error ×1.96; assumed Gaussian 1σ" (HU: „idézett hiba ×1.96; feltételezett Gauss 1σ”).
- Own ocean bound curve: marked CONDITIONAL. Model curves: "idealized layer/half-space model" / "idealized half-space". Tower panel title uses the same error wording. Adelberger 95% label unchanged.

## §5 open questions (summaries section 2c, method A8/B6/C, JSON)
- Ocean: original uncertainty, profile weighting, pressure-derived depth not reconstructed (JSON `open_questions`).
- Towers: real surface sampling, terrain, depth density, Poisson continuation response not computed.
- BREN ~100 µGal systematic: separate nuisance treatment needed; illustrative bounded ±100 µGal nuisance moves the ρ = 2000 edge 0.00756 → 0.01029 (ρ = 2670: 0.00565 → 0.00769), illustrative only.
- λ < 1 cm: reworded in text and JSON (`verdict.untested_candidates_to_retest`, replacing `what_would_relieve_the_tension`) as an untested candidate to re-test, not a verified escape range.
- R ranges are scenarios, not a combined confidence interval, never σ.
- March eq. (81): separate scenario; the discrepancy cannot be resolved by a sign or G-convention change; author intent not verified.

## Verdict wording
- Kept: "Our own re-derivation indicates strong tension with the fixed, unscreened Yukawa branch; confirming an exclusion requires verification of the measurement response and the uncertainty treatment."
- Added: strong tension only in the sense that, in the idealized layer/half-space models under the assumed error reading, a large numerical discrepancy results; not an empirical exclusion with a confidence level. R = ratio to bound, stays CONDITIONAL. Zumberge data = secondary-source data via March 2019. Missing primary sources listed separately (Zumberge 1991, Cruz 1991 full text, Jekeli 1990 full text).
