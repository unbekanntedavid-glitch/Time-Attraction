# Covariant Action Candidates and Nonlocal Static Consistency Tests for the Time Attraction Law (Task `50aa2db0`)
# Kovariáns hatásjelöltek és nemlokális statikus konzisztenciavizsgálatok az Idővonzási Törvényhez (`50aa2db0` feladat)

**Author:** David Szabolcsi (ORCID: [0009-0005-6300-074X](https://orcid.org/0009-0005-6300-074X))  
**Date:** 2026-10-24 (audit revisions 2026-10-25, incl. follow-up audit)  
**License:** CC BY-NC 4.0  
**Related Archival DOIs:** V14 (`10.5281/zenodo.23217189`), V17 (`10.5281/zenodo.19666220`), V19 (`10.5281/zenodo.23497128`), V20 (`10.5281/zenodo.23540741`), V21 (`10.5281/zenodo.23602501`), V22 Book Edition (`10.5281/zenodo.23604423`)

---

## English Overview (audit revision)
1. **Action candidates & normalization** (`closure_tov_comparison.csv`): Formulation I (nonlocal 3D top-hat, V21 algebraic closure), a flat-kernel fixed-point equilibrium with bilateral reciprocal force (`bilateral_eom`, unweighted; `bilateral_lapse_weighted`, $V_{\text{eff},N} = f_L - (\gamma_L/N)\int N K f_L\varepsilon\,dV$ with $P_{\text{matter}} = f_L P$; SLy/APR4/H4, $L = 1, 10$ km; not a full covariant EL solution, since the metric variation of kernel and volume element is absent; 2-cell check (0.74881594, 0.72084293) vs. 0.73160178 in the JSON), Formulation II (Yukawa) with explicit `alpha_Y_unit` and `alpha_std_poisson` ($\alpha_Y = 4\pi\alpha_{\text{std}}$; $\alpha_{\text{std}} = 1/(4\pi)$ reproduces the homogeneous top-hat limit), `yukawa_flat_phenom_Tphi` (phenomenological TOV model with flat-space Yukawa field and added scalar stress) and `yukawa_action_el` (action-consistent: $f_L$-weighted source, $P_{\text{matter}}=f_LP$, $h'=-\nu'+\phi'$, curved operator with $B=1-2m/r$, exterior match $\phi'(R)=-(1/R+1/L)\phi(R)$, tail mass `Delta_M_phi_ext_*`). SymPy: `sym_yukawa_el_source_factor_residual`, `sym_yukawa_el_hydrostatic_residual` = 0.
2. **Symbolic checks** (`conservation_Residuals.json`): conservation identities, `sym_disformal_g00_residual` ($\tilde g_{\mu\nu} = g_{\mu\nu} + (1-f_L^2)u_\mu u_\nu/c^2 \Rightarrow \tilde g_{00} = -f_L^2$), `sym_yukawa_center_limit_residual` ($\phi_L(0) = \alpha_Y (G/c^2)\int r' e^{-r'/L}\rho\,dr'$, no extra $L$) — all exactly 0.
3. **Three-tier residuals**: symbolic 0; algebraic X-sector quadrature $4.44\times10^{-16}$ (mass, $/|\rho_X\nu'|$) and $4.39\times10^{-14}$ (gradient, $/|P'_{r,X}|$); full spline divergence $/|P'|$: $2.70\times10^{-5}$ (mass), $2.36\times10^{-5}$ (gradient). Figure 3 shows both (y-range $10^{-18}$–$10^{-3}$).
4. **Static nonlocal tangent model and quasi-dynamic BVP test on the V21 background** (`coupled_stability_results.csv`; linearized around V21 `dense_profile`, not the new branches): selected BVP-branch zero at $+6.97\times10^{-5}$ dex (mass) and $+7.39\times10^{-5}$ dex (gradient) from $dM/d\rho_c=0$. This verifies the static linearization; $W\to4W$ gives $\omega^2$ ratio 0.2500 with identical $F(0)$, so the zero is inertia-invariant. Symmetrized Galerkin kinetic form has positive diagonal. Column `omega0_sq_sym_tangent_s2`: first-order tangent estimate. The curved-space window variation remains open.
5. **Clock coupling** (`clock_coupling_H0_H1.csv`): V20 granite chamber IN ($\mathbf r = 0$) and OUT ($D_{\text{out}} = 2.0$ m): `GR_newtonian_IN` $=-9.269\times10^{-24}$, `GR_newtonian_IN_minus_OUT` $=-6.094\times10^{-24}$, `H0_TAL_excess_IN` $=1.331\times10^{-47}$, `H1_shift_IN` $=-1.436\times10^{-24}$, `H1_shift_IN_minus_OUT` $=-1.436\times10^{-24}$ ($L=1$ m). Sirius B: H0 and H1 on fixed-$M$ and on a common effective-mass proxy (fixed $n=3$ Lane–Emden profile with the H0 effective-mass shift, not a newly solved equilibrium). Granite $\sigma_y=10^{-18}$ is an illustrative benchmark.

## Magyar nyelvű összefoglaló (audit utáni változat)
1. **Hatásjelöltek és normálás**: I. megfogalmazás (V21 algebrai lezárás), bilaterális reciprok erővel kiegészített, síkbeli kernelű fixpontos egyensúly (`bilateral_eom`, illetve lapse-súlyozott `bilateral_lapse_weighted`), II. megfogalmazás explicit `alpha_Y_unit` / `alpha_std_poisson` oszlopokkal, síkbeli Yukawa-mezővel és hozzáadott skalárfeszültséggel számolt fenomenológiai TOV-modell (`yukawa_flat_phenom_Tphi`), valamint a hatással konzisztens `yukawa_action_el` ág.
2. **Szimbolikus ellenőrzések**: a diszform metrika helyes előjele ($\tilde g_{00} = -f_L^2$) és a Yukawa középponti határérték (nincs többlet $L$) — maradék 0.
3. **Háromszintű maradékok**: szimbolikus 0; algebrai X-szektor $4{,}44\times10^{-16}$ / $4{,}39\times10^{-14}$; spline-divergencia $2{,}70\times10^{-5}$ / $2{,}36\times10^{-5}$.
4. **Statikus érintőoperátor**: a $\sim7\times10^{-5}$ dex egyezés a statikus linearizációt igazolja; $W\to4W$ mellett a zérushely változatlan.
5. **Óratáblázat**: gránitkamra IN–OUT differenciális jel (σ_y = 10⁻¹⁸ szemléltető érték), Szíriusz B közös effektívtömeg-proxyn is.

## Reproducibility
Run `python3 covariant_closure_50aa2db0.py` inside the unpacked folder (numpy, scipy, pandas, sympy, matplotlib). All helper modules (`tal_profiles.py`, `stellar_tal_v21.py`, `radial_modes.py`) and the V21 reference table (`v21_radial_summary_ref.csv`) are included; paths are relative to the script.

## Directory Contents / Fájlok jegyzéke
- `covariant_closure_50aa2db0.py` — calculation and verification script
- `tal_profiles.py`, `stellar_tal_v21.py`, `radial_modes.py` — V21 helper modules
- `v21_radial_summary_ref.csv` — V21 radial turning-point reference table
- `tal_covariant_closure_50aa2db0_EN_2026-10-24.pdf` & `.tex` — 4-page English report
- `tal_covariant_closure_50aa2db0_HU_2026-10-24.pdf` & `.tex` — 4-page Hungarian report
- `closure_tov_comparison.csv`, `coupled_stability_results.csv`, `conservation_Residuals.json`, `clock_coupling_H0_H1.csv` — data tables
- `fig1_action_profiles.png`, `fig2_coupled_stability.png`, `fig3_conservation_check.png` — figures
- `SHA256SUMS.txt` — SHA-256 checksum manifest (17 files)
