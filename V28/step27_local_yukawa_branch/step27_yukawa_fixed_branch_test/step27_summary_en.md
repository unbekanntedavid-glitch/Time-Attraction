# Step 27 (part of V28): local gravity test of the fixed Yukawa branch

Series: Time Attraction Law — Scalar-Tensor Yukawa Fit, V28. V28 covers the V27 Step 26 cosmology package, the interactive simulator and this Step 27. No V28 DOI exists yet. Earlier DOIs: V24 10.5281/zenodo.23680715, V25 10.5281/zenodo.23699685, V26 10.5281/zenodo.23727585, V27 10.5281/zenodo.23761574.

Scope: only the **fixed, unscreened, universal local Yukawa branch** of the TAL scalar-tensor Yukawa fit (L = 10 km, α_Y = 1, β_Yuk = 1/(4π), C = 0 attractor). The result does not apply automatically to other TAL variants.

## 1. What the package predicts

Convention: V(r) = −G_∞ m₁m₂/r · (1 + α e^{−r/λ}). A laboratory measurement (r ≪ λ) gives G_lab = G_∞(1+α).

- Action in `hubble_desi_bao_test/tal_cosmology_common.py`, lines 50–64: the kinetic term is c⁴/(4πGβ_Yuk), and matter couples universally to A²(φ)g with A = e^{α_Y φ}. For a point source, φ = −α_Y β_Yuk GM e^{−r/L}/(c²r). The potential on a test mass is c²α_Y φ = −α_Y²β_Yuk GM e^{−r/L}/r.
  - So **α = q = α_Y²β_Yuk = 1/(4π) = 0.0795775** and **λ = L = 10 km**.
  - The sign is positive, so the extra force is attractive.
- The same result appears in three other places:
  - `step26_md_templates.py`, line 103: G_lab ≈ G_*(1+q) "for r ≪ L and unscreened small bodies".
  - `step5_coupling_dictionary_ppn_screening_en.md`, lines 54, 73 and 101: α_std = α_Y/(4π), and G_eff(r) = G[1+α_std(1+r/L)e^{−r/L}].
  - `run_step26_revision_controls.py`, line 244: G_eff/G_lab = 1/(1+q) = 0.9262883.
- Screening: **none**.
  - `step26_md_templates.py`, line 106, and `run_step26_revision_controls.py`, line 253, both say "no screening branch is introduced".
  - The chameleon/symmetron mechanism in Step 11 (`step11_..._en.md`, line 18) is a separate nonlinear extension hypothesis, not part of the fixed branch.
  - In Step 5, line 100, the word "screened" only means the e^{−r/L} decay.
- Composition dependence: the coupling is universal, because all matter couples to A²g. There is no equivalence-principle signal, so the test reduces to the inverse-square law.

Note on CODATA G: every laboratory G measurement has r ≪ 10 km, so they all measure the same G_∞(1+α). The spread between them therefore carries no signal. The signal appears only on scales of order 1–10 km.


## 2. Bounds and model observables (λ = 10 km)

R = 0.0795775 / α_lim(10 km) is the **ratio to bound**.
- It is meaningful only in the convention above and under the stated assumptions.
- R > 1 means the model point lies outside the allowed region *of that record*.
- R is **not** a significance and is not quoted in σ.
- **Every R in this table is conditional**: no denominator is a verified, directly citable primary bound at 10 km.

Data labels:
- **secondary-source data (not directly verified primary data)**: numbers taken from a later paper that quotes the original (Zumberge inputs via March et al. 2019).
- **primary abstract text**: numbers from the original paper's abstract; full text not obtained.
- **primary report**: the original document itself was obtained (BREN, OSTI 6190822).

Column guide: **Model (a)** is the fixed model point, α = 0.0796 at λ = 10 km. **Model (b)** is the experiment's own observable computed for the model in that setup.

| Source | Data label | Observable | Assumptions / CL | α_lim at 10 km | Model (a) | Model (b) | R (conditional) |
|---|---|---|---|---|---|---|---|
| Zumberge et al., PRL 67, 3051 (1991), doi:10.1103/PhysRevLett.67.3051, inputs via March et al., PRD 100, 042002 (2019), arXiv:1904.12789 | secondary-source data (not directly verified primary data) | Slope of the gravity residual in the water column, 500–4800 m: −0.060 ± 0.172 mGal/km (March eq. 74) | α > 0, universal, no screening; idealized layer model; quoted error ×1.96 under an assumed Gaussian 1σ reading (one-sided 97.5% only under that assumption) | Published 0.002 only for 1–5000 m; **no 10 km value**. Re-derived (exact inversion): **0.00226** (≈0.0023); linear approximation 0.00225 | 0.0796 | Slope **+9.1 mGal/km** | **≈35** (35.27; linear approx. 35.35; 34.6 with rounded 0.0023) |
| Jekeli, Eckhardt, Romaides, PRL 64, 1204 (1990), doi:10.1103/PhysRevLett.64.1204 | primary abstract text | ~600 m tower vs upward continuation | Sign, λ range and CL not stated | **Not applied.** Abstract value α < 0.001, λ range unknown; no 10 km value (text-only reference) | 0.0796 | not computed | not applied |
| Thomas et al., PRL 63, 1902 (1989), doi:10.1103/PhysRevLett.63.1902 | primary abstract text | Residual at 465 m: (−60 ± 95) µGal | Idealized half-space, ρ assumed 2000–2670 kg/m³, terrain ignored; ± read as 1σ | Re-derived: 0.0049–0.0065 | 0.0796 | −2.8 to −3.7 mGal | ≈12–16 |
| Romaides et al., PRD 50, 3608 (1994) | primary abstract text | Residual at 493 m: (−33 ± 30) µGal | Same | Re-derived: 0.0017–0.0023 | 0.0796 | −3.0 to −4.0 mGal | ≈35–47 |
| Romaides et al., PRD 55, 4532 (1997), doi:10.1103/PhysRevD.55.4532 | primary abstract text | Residual at 568 m: (32 ± 32) µGal | Same | Re-derived: 0.00050–0.00066 | 0.0796 | −3.4 to −4.6 mGal | ≈120–160 |
| Thomas, Kasameyer et al., LLNL report, OSTI 6190822 (preliminary BREN analysis) | primary report | "within 93 µGal of observations at all heights", uncertainties up to 95 µGal (1σ, sampling); sign not stated | Conservative μ = −93 µGal at 454 m; ρ assumed 2000–2670 | Re-derived: 0.0056–0.0076 (if the sign were +: 0.0019–0.0025); the ~100 µGal systematic is not included (see 2c) | 0.0796 | −2.7 to −3.7 mGal | ≈11–14 |
| Cruz, Harrison, Speake, Niebauer, McHugh, Keyser, Faller, Mäkinen, Beruff, J. Geophys. Res. 96(B12), 20073 (1991), doi:10.1029/91JB01756 | primary abstract text | 300 m Erie tower; rms difference 10 µGal | Bound on \|α\|; CL not stated | **\|α\| < 0.001 at λ = 1000 m**. At 10 km: ≈0.00088 is only a **conditional estimate** (see below) | 0.0796 | not computed | 80 at 1 km; not used at 10 km |
| Context: Adelberger, Heckel, Nelson, Annu. Rev. Nucl. Part. Sci. 53, 77 (2003), Fig. 4 | digitized compilation | — | 95% | ≈4×10⁻⁴, digitized; reading uncertainty ≈ ±0.1 dex | 0.0796 | not computed | context only |
| Context: Fischbach & Talmadge, arXiv:hep-ph/9606249 (1996), Fig. 1 | digitized compilation | — | 2σ; axis is −α, i.e. **α < 0 only** | ≈4×10⁻⁴ (3–5×10⁻⁴), digitized; ≈ ±0.15 dex | 0.0796 | not computed | not used (opposite sign) |

**Erie (Cruz et al. 1991), 10 km value.** The 0.00088 figure is a conditional estimate, not a bound. A 1 km limit cannot be transferred to 10 km from the Yukawa exponential alone. That needs the full response function of the tower measurement difference, including the source geometry and the G normalization, and we do not have it (the publisher full text returned HTTP 403). Our scaling assumed: (i) the 1 km limit is set by the top-of-tower residual at h = 300 m; (ii) a flat half-space of uniform density; (iii) the Yukawa part of the surface field is continued as a constant; (iv) the same G normalization at both λ. Under those assumptions the factor is S(1 km)/S(10 km) with S(λ) = λ(1 − e^{−h/λ}), i.e. 259.2/295.5 = 0.877.

**G normalization of the bounds (review r1, S01).** The measured G_lab is held fixed for every trial α, so G_∞(α) = G_lab/(1+α). Towers: Δg(α) = −Kα/(1+α), K = 2πG_lab ρλ(1 − e^{−h/λ}), and the exact root is α_lim = E/(K − E) with E = −(μ − 1.96σ). Ocean: s(α) = Sα/(1+α), α_lim = edge/(S − edge). The previous version froze G_∞ at the model point q while searching α; its tower limits were 7–8% too large (e.g. BREN ρ = 2000: 0.00810 → 0.00756). Model residuals at α = q are unchanged. A back-substitution test in the code confirms that every limit point returns the allowed residual edge (relative deviation ≤ 1×10⁻¹⁰).

The re-derivation method, with formulas, source equation numbers and known gaps, is in `step27_rederivation_method_en.md`.

## 2a. Decisive bound and its status

- **Decisive record: the ocean re-derivation**, α_lim(10 km) ≈ 0.0023, **R ≈ 35 (conditional)**. R_status: *conditional: denominator α_lim ≈ 0.0023 is an unverified re-derivation*. Inputs are secondary-source data (not directly verified primary data).
- Sensitivity of that number to choices we could not check:
  - reading ±0.172 as the bounds of an interval, as March et al. do in eqs. (75)–(76): α_lim ≈ 0.00091, R ≈ 87;
  - using March eq. (81) as printed for the layer integral: α_lim ≈ 0.0035, R ≈ 23 (exact 0.00350, R = 22.7). This is a **separate scenario**, not a statistical error band. Our direct layer integral and a brute-force integral disagree with the printed eq. (81); the discrepancy cannot be resolved by a sign or G-convention change, and the authors' intent or any correction was not verified.
- The ocean figures 0.0009–0.0035 (R ≈ 23–87) are **scenarios** (error reading, layer formula), not one combined confidence interval. R is never a σ significance.
- The ocean α_lim is the exact inversion edge/(S − edge) = 0.00225606 (R = 35.27); the linear approximation edge/S = 0.00225098 (R = 35.35) is reported only as an approximation.
- No primary source found states a numeric α limit at exactly λ = 10 km for α > 0.

## 2b. Missing primary sources

Listed separately from the data used. The Zumberge values used here are **secondary-source data via March et al. 2019**, not directly verified primary data. Key missing primary sources: **Zumberge et al. 1991** (not obtained), **Cruz et al. 1991 full text**, **Jekeli et al. 1990 full text**.

- Zumberge et al., PRL 67, 3051 (1991): **not obtained**; ocean inputs come from March et al. 2019.
- Cruz et al., JGR 96, 20073 (1991): abstract only (publisher PDF HTTP 403).
- Jekeli et al., PRL 64, 1204 (1990): abstract only.
- Thomas et al., PRL 63, 1902 (1989): abstract only; the preliminary LLNL report was obtained.
- Romaides et al. 1994 and 1997: abstracts only.
- Not reached: Stacey et al. RMP 59, 157 (1987); Moore et al. PRD 38, 1023 (1988); Ander et al. PRL 62, 985 (1989); Bartlett & Tew PRL 63, 1531 (1989); Achilli et al. 1997; Rapp 1987; Fischbach & Talmadge 1999 (book).

## 2c. Open questions (not reconstructed in Step 27)

- **Ocean:** the original uncertainty, the weighting of the measured profiles and the model dependence of the pressure-derived depth are not reconstructed. The analytic layer check only verifies the idealized response.
- **Towers:** the actual surface sampling, terrain, depth-dependent density and the joint Yukawa response of the Poisson continuation are not computed.
- **BREN:** the report's ~100 µGal systematic needs a separate nuisance treatment; it cannot be added automatically as an independent Gaussian σ to the 95 µGal. Illustration only: a bounded ±100 µGal nuisance taken at its worst-case edge would move the ρ = 2000 half-space edge from 0.00756 to 0.01029 (ρ = 2670: 0.00565 → 0.00769). This is not a new data fit and not an empirical bound.
- **R ranges** across error readings, layer formulas and assumed densities are scenarios, not one combined confidence interval; R is never a σ significance.
- **March eq. (81):** kept as a separate scenario; the discrepancy is not resolvable by a sign or G-convention change; the authors' intent was not verified.
- **λ < 1 cm:** only an untested candidate to be re-tested (see section 3); Step 27 did not examine whether that range is allowed.

## 3. Verdict

**Our own re-derivation indicates strong tension with the fixed, unscreened Yukawa branch; confirming an exclusion requires verification of the measurement response and the uncertainty treatment.**

Meaning of "strong tension": **in the idealized layer/half-space models, under the assumed error reading, a large numerical discrepancy results.** This is not an empirical exclusion with a confidence level: the actual measurement response, the background corrections and the systematic uncertainties are not reconstructed.

Scope: only the fixed, unscreened, universal local Yukawa branch (α = 0.0796, λ = 10 km); not all TAL variants, and not a re-opening of Step 26. R is the **ratio to bound** and stays **CONDITIONAL**; all applied conditional R values are above 1 (≈11–160, scenarios, not a confidence interval). R is a ratio, not a σ value.

Untested candidates (each would define a new branch that must be re-tested; none is verified):
- a screening mechanism effective for Earth-scale and laboratory sources (the package has none);
- α ≲ 10⁻³ at λ = 10 km, i.e. α_Y²β_Yuk reduced by a factor of about 80 or more;
- a different λ (for example λ < 1 cm) with unchanged α: an **untested candidate** to be re-tested, **not** a verified escape range; Step 27 did not examine whether that range is allowed.

Changes in review round 1: `step27_review1_changes_en.md`.

Reproduce: `python run_step27.py` → `step27_results.json` (strict JSON, including `back_substitution_tests`), plus the figures `step27_alpha_lambda_{en,hu}.png` and `step27_observables_{en,hu}.png`.
