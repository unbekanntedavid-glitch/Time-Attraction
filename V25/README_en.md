# Scalar-Tensor and Yukawa/Disformal Action Mapping (V23 $\leftrightarrow$ Jordan–Brans–Dicke / Dilaton Framework) — Steps 1–22 Derivation and Audit Reproduction Package

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Publication status and version delineation (V24 vs. V25):** The earlier Zenodo publication (**V24**, DOI: `10.5281/zenodo.23680715`) refers exclusively to the earlier V24 release. The sign fixes in `curved_yukawa_adm_and_stability.py` (the exterior metric coupling term `a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0` at line 593 of `linearized_exterior_response`, and setting the exact left-hand side sign `scalar_m_sign = -1.0`, i.e., `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n`, as the default in the $L_u u = R_h\,\delta h$ operator of `solve_yukawa_coupled_omega2_gl2`), the three new diagnostics in `curved_yukawa_audit_summary.json` (`full_jacobian_multi_eos_and_points`, `baryonic_mass_and_first_law_audit`, `eigenmode_xi_nodes_and_static_tangent_audit`), and the closure of the theoretical verification belong exclusively to the new **V25** version — the earlier V24 DOI must never be cited as if these fixes and closures were already contained in it.

**V25 diagnostics summary (exact data from `curved_yukawa_audit_summary.json` and `curved_yukawa_coupled_stability.csv`):**
- Derivative verification (SLy, $\log_{10}\rho_{0c}=15.0$, $L=10$ km): `exterior_dphi_prime_dm_R` $=-1.683651\times10^{-4}\;\text{km}^{-2}$ (rel. diff $3.98\times10^{-10}$); `interior_r6km_a_m_sc` ($\partial\phi_L''/\partial m$) $=-1.728334\times10^{-4}\;\text{km}^{-3}$ (rel. diff $6.08\times10^{-9}$).
- Full 4×4 Jacobian: 4 configurations (SLy and APR4), 108 checks, max relative error $2.74\times10^{-5}$.
- Baryonic mass / First Law: `delta_dex_ADM_vs_baryon_turning` between $M_{\text{ADM}}$ and the adaptive Gauss–Kronrod baryonic turning point ($M_B^{\text{quad}}$) is $-4.58\times10^{-6}$ to $+8.43\times10^{-6}$ dex; non-conformal excess $8.23\times10^{-5}\,M_\odot$/dex.
- Main stability zero crossing and eigenmode audit (with exact LHS sign `scalar_m_sign = -1.0`): the main dynamical $\omega_0^2=0$ zero crossing in `curved_yukawa_coupled_stability.csv` for SLy $L=10\;\text{km},\;\alpha_Y=1$ is **`15.316410`** ($+5.80\times10^{-6}$ dex above the ADM peak `15.316404` and $+1.33\times10^{-6}$ dex above the $M_B^{\text{quad}}$ peak `15.316409`; across all 5 sequences $\le 6.9\times 10^{-6}$ dex offset, nodeless $n=0$ fundamental mode, and $(0.6\text{--}5.0)\times 10^{-4}$ relative $L^2$ agreement with the static equilibrium tangent; the diagnostic `scalar_m_sign = +1.0` root is `15.319129`, $+2.725\times 10^{-3}$ dex — details: `step22`, Sections 4 and 7).

**Date:** 2026-10-27  
**Working Directory:** `scalar-tensor-yukawa-fit`  
**Self-Contained Audit Archive:** `scalar-tensor-yukawa-fit_audit_package.zip`  
**License:** CC BY-NC 4.0 — Commercial use permitted only with the prior written consent of the author.  
**Official Archival DOIs:**  
- **V14:** `10.5281/zenodo.23217189`  
- **V17 (SPARC 155 / 123 blind sample):** `10.5281/zenodo.19666220`  
- **V19 (3D spherical top-hat SPARC, 175 galaxies):** `10.5281/zenodo.23497128`  
- **V20 (Two-clock granite chamber and 8-coil magnetic test):** `10.5281/zenodo.23540741`  
- **V21 (Compact stars, TOV, and Sirius B audit):** `10.5281/zenodo.23602501`  
- **V22 Book Edition (Trilingual monograph, HU/EN/DE):** `10.5281/zenodo.23604423`  
- **V23 (Covariant action closure and coupled radial stability, `50aa2db0`):** `10.5281/zenodo.23639084`  

---

## 1. Research Question, Theoretical Framework & Branch Taxonomy

This audit package contains the exact mathematical correspondence and numerical verification mapping the covariant action candidates developed in the V19–V23 manuscripts and monograph edition of the **Time Attraction Law (TAL)** (`10.5281/zenodo.23604423`, `10.5281/zenodo.23639084`) — **Formulation I** (nonlocal, positive-definite bilateral top-hat kernel action with auxiliary fields) and **Formulation II** (local, screened Helmholtz–Yukawa scalar-field action principle) — onto general **scalar–tensor theories of gravity** (Jordan–Brans–Dicke family, low-energy string-theory **dilaton** coupling, and Bekenstein-type **disformal** clock metrics).

### 1.1. Core Equations & Canonical Normalization
1. **Exponential Time Factor and Scalar/Dilaton Field (`step1_dilaton_form.md` / `step1_dilaton_form_en.md`):**
   $$f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(\mathbf{x})\bigr) \;\equiv\; e^{-\phi_L(\mathbf{x})}, \qquad \phi_L(\mathbf{x}) \;=\; \gamma_L\,\bar{\rho}_L(\mathbf{x}) \;\ge\; 0,$$
   where the canonical coupling coefficient originates from the gravitational volumetric capacity ($C = \frac{4\pi G}{c^2} = 9.33306\times 10^{-27}\;\text{m/kg}$):
   $$\gamma_L \;=\; \frac{G L^2}{c^2} \;=\; \frac{C}{4\pi}\,L^2, \qquad m_\phi \;=\; \frac{\hbar}{L c}.$$
2. **Conformal Jordan Metric vs. Medium-Anchored Disformal Clock Metric (`step2`, `step5`, `step8`, `step14`):**
   $$\tilde{g}_{\mu\nu}^{\text{conf}} \;=\; A^2(\phi_L)\,g_{\mu\nu} \;=\; f_L^2\,g_{\mu\nu} \qquad \text{vs.} \qquad \tilde{g}_{\mu\nu}^{\text{disf}} \;=\; g_{\mu\nu} + \frac{1 - f_L^2}{c^2}\,u_\mu u_\nu \;\Longrightarrow\; \tilde{g}_{00}^{\text{disf}} \;=\; f_L^2\,g_{00}\;(\text{in the comoving frame}).$$
   In a comoving ($u^\mu = (c/\sqrt{-g_{00}},\mathbf{0})$, $u_\mu = (-\sqrt{-g_{00}}\,c,\mathbf{0}) \implies u_0 u_0/c^2 = -g_{00}$) baryonic medium and in atomic clocks, both metrics yield the same proper-time factor ($d\tilde{\tau} = f_L\,d\tau_g$). For ultrarelativistic light rays ($g_{\mu\nu}k^\mu k^\nu = 0$) and PPN tests, however, the two branches separate sharply:
   - In the universal **conformal Jordan branch (`Branch II-Conf`)**, $\sqrt{-\tilde{g}^{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$, the Maxwell action is conformally invariant, and the standard Brans–Dicke / PPN dictionary applies:
     $$\alpha_{\text{std}} \;=\; 2\beta^2 \;=\; 2\tilde{\alpha}_L \;=\; \frac{\alpha_Y}{4\pi}, \qquad \omega_{\text{BD}} \;=\; \frac{1 - 3\alpha_{\text{std}}}{2\alpha_{\text{std}}} \;=\; \frac{4\pi}{\alpha_Y} - \frac{3}{2}, \qquad \gamma_{\text{PPN}}^{\text{conf}}(r) \;=\; \frac{1 - \alpha_{\text{std}}(1+r/L)e^{-r/L}}{1 + \alpha_{\text{std}}(1+r/L)e^{-r/L}}.$$
   - In the medium-anchored **disformal branch (`Branch II-Disf`)**, $\sqrt{-\tilde{g}^{\text{disf}}} = f_L\sqrt{-g}$, magnetic energy ($\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$) directly sources $\phi_L$ (`step15`), and light is deflected by the achromatic refractive index $n_t(s) = f_L(s)^{-1} = e^{+\phi_L(s)}$ (`step9`, `step14`), so the conformal Brans–Dicke $\gamma_{\text{PPN}}^{\text{conf}}$ formula does not apply to this branch.
3. **V23 Local Helmholtz–Yukawa Action Principle, Exact Normalization, and Curved Exterior Boundary-Value Problem (`step3`, `step5`, `step7`, `step22`):**
   The exact SI-unit coefficient of the scalar kinetic and mass terms (with $\alpha_Y \equiv 8\pi\tilde{\alpha}_L = 4\pi\alpha_{\text{std}}$) is:
   $$\boxed{\;S_{\text{II}}[g_{\mu\nu},\phi_L,\Psi] \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}\left(R \;-\; \frac{1}{\tilde{\alpha}_L}\Bigl(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{\phi_L^2}{L^2}\Bigr)\right) \;+\; e^{-\phi_L}\,\mathcal{L}_{\text{matter}}\right] \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}R \;-\; \frac{c^4}{2 G\,\alpha_Y}\Bigl(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{\phi_L^2}{L^2}\Bigr) \;+\; e^{-\phi_L}\,\mathcal{L}_{\text{matter}}\right],\;}$$
   which in geometrized ($G = c = 1$) units yields the V23 code form $\mathcal{L}_\phi = -\frac{1}{2\alpha_Y}\bigl((\nabla\phi_L)^2 + \phi_L^2/L^2\bigr)$ and the exact curved-spacetime Euler–Lagrange equation (`yukawa_action_el`, $B(r) = 1 - 2m(r)/r$):
   $$B\,\phi_L'' \;+\; \left[\frac{2B}{r} + B\nu' + \frac{B'}{2}\right]\phi_L' \;-\; \frac{\phi_L}{L^2} \;=\; -\frac{\alpha_Y G}{c^2}\,e^{-\phi_L}\,\varepsilon(r), \qquad \rho_\phi(r) \;=\; \frac{c^2}{2 G\,\alpha_Y}\left[B\bigl(\phi_L'\bigr)^2 + \frac{\phi_L^2}{L^2}\right].$$
   In the curved exterior vacuum outside the stellar surface ($r \ge R_*$, where $\varepsilon = P = 0$, $B\nu' + B'/2 = \frac{2m}{r^2} - \frac{4\pi r}{\alpha_Y L^2}\phi_L^2$), the coupled exterior Einstein–scalar system (`step22`, `curved_yukawa_adm_and_stability.py`) is:
   $$\left(1 - \frac{2m(r)}{r}\right)\phi_L'' \;+\; \left[\frac{2}{r}\left(1 - \frac{m(r)}{r}\right) - \frac{4\pi r}{\alpha_Y L^2}\phi_L^2\right]\phi_L' \;-\; \frac{\phi_L}{L^2} \;=\; 0, \qquad \frac{dm(r)}{dr} \;=\; 4\pi r^2\rho_\phi(r),$$
   which uniquely determines the genuine curved Dirichlet-to-Neumann slope $\phi_{L,\text{curved}}'(R_*) = -k_{\text{curved}}\phi_L(R_*)$ (for SLy $\log_{10}\rho_{0c}=15.0$, $L=10\;\text{km}$, $\alpha_Y=1$, $k_{\text{curved}} = 0.233542\;\text{km}^{-1}$, which is $+24.74\%$ steeper than the flat-space Robin approximation $k_{\text{flat}} = 1/R_* + 1/L = 0.187216\;\text{km}^{-1}$) and the asymptotic **ADM mass**:
   $$M_{\text{ADM}} \;\equiv\; \lim_{r\to\infty} m(r) \;=\; M_{\text{surf}} \;+\; \Delta M_{\phi,\text{ext}}^{\text{curved}}, \qquad M_{\text{surf}} \;\equiv\; m(R_*).$$

### 1.2. Explicit Branch Taxonomy Table

The assumptions, kernel moments, and domain of validity of the six theoretical branches examined in this package are separated in the table below (`step2`, `step5`, `step10`, `step11`, `step21`, `step22`):

| Branch | Action Principle / Metric Definition | Kernel Moments (Expansion of $\bar{\rho}_L(\mathbf{x})$) | Interior Hydrostatic Balance | Exterior ($r > R_*$) Vacuum and Mass | Domain of Validity & Steps |
|---|---|---|---|---|---|
| **Branch I-A** (V21 nonlocal 3D top-hat: `tophat_mass`, `tophat_gradient`) | Phenomenological 3D spherical top-hat average + algebraic Bianchi X-sector | $\left(1 + \frac{L^2}{10}\nabla^2 + \mathcal{O}(L^4)\right)\rho$ | `mass`: $h' = -\nu'$; `gradient`: $h' = -\nu' = -f_L A_{\text{GR}}$ | Compact support $[R_*, R_*+L]$; $m(r) = M_{\text{surf}}$ constant for $r \ge R_*$ | V21 769 models, V23 Eqs. (8)–(9) static tangent BVP (`step7`, `step16`–`step18`, `step21`) |
| **Branch I-B** (V23 Formulation I: `bilateral_lapse_weighted`) | Auxiliary-field ($\phi_L,\lambda_L$) bilateral kernel action: $S_m = -\int\sqrt{-g}\,e^{-\phi_L}\varepsilon$ | $\left(1 + \frac{L^2}{10}\nabla^2 + \mathcal{O}(L^4)\right)\rho$ | $h' = -\nu' - (\ln f_{\text{recip}})'$ (reciprocal lapse-weighted force) | Compact support $[R_*, R_*+L]$ | V23 reciprocal bilateral potential $V_{\text{eff},N}(r)$ and Synge boundary variation (`step10`, `step13`) |
| **Branch II-Fluid** (V23 Formulation II: `yukawa_action_el` / `yukawa_action_el_curved_adm`) | Local Helmholtz–Yukawa action with fluid-density coupling: $\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$ | $(1 - L^2\nabla^2)^{-1}\rho = \left(1 + L^2\nabla^2 + \mathcal{O}(L^4)\right)\rho$ | $h'(r) = -\nu'(r) + \phi_L'(r)$ (`sym_yukawa_el_hydrostatic = 0`) | Curved Einstein–scalar tail $r \in [R_*,\infty)$; $M_{\text{ADM}} = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$ | Curved exterior ADM mass, mass–radius curves, and coupled $(\delta h,\delta\phi_L,\delta m)$ stability (`step3`, `step4`, `step7`, `step12`, `step22`) |
| **Branch II-Conf** (Universal conformal Jordan branch) | $\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$, $\sqrt{-\tilde{g}^{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$, $S_m[\tilde{g}^{\text{conf}},\Psi]$ | $(1 - L^2\nabla^2)^{-1}\rho_{\text{trace}}$ | $h' = -\nu' + \phi_L'$ (trace-free EM does not source $\phi_L$) | Yukawa scalar tail; conformally invariant null geodesics | Standard Brans–Dicke / PPN dictionary ($\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$, $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2}$, $\gamma_{\text{PPN}}^{\text{conf}}(r)$; `step1`, `step2`, `step5`, `step14`) |
| **Branch II-Disf** (Medium-anchored purely temporal disformal branch) | $\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$, $\sqrt{-\tilde{g}^{\text{disf}}} = f_L\sqrt{-g}$ | Top-hat or Yukawa kernel on comoving energy density | Comoving clocks $d\tilde{\tau} = f_L d\tau_g$; $u_\mu$ variation yields extra constraint | Magnetic energy ($\rho_B = B^2/(2\mu_0 c^2)$) sources $\phi_L$; refractive index $n_t = e^{+\phi_L}$ | V20 8-coil magnetic test, V20 granite chamber H1, Sirius B H1, VHE gamma dispersion (`step2`, `step8`, `step9`, `step14`, `step15`, `step19`, `step20`) |
| **Branch II-Cham** (Nonlinear chameleon/symmetron extension) | Hypothetical density-dependent $\beta_{\text{eff}}(\rho)$ and self-interaction potential $V_0(\phi_L)$ | Density-dependent $L_{\text{eff}}(\rho) = [\partial^2 V_{\text{eff}}/\partial\phi_L^2]^{-1/2}$ | Screened (thin shell) at high density, amplified at low density | Environment-dependent Yukawa range | Separate theoretical extension hypothesis beyond the constant-$\alpha_Y$ V23 action (`step6`, `step11`) |

*(Note on the top-hat $\leftrightarrow$ Yukawa correspondence: the second-order Taylor moment of the 3D spherical top-hat kernel is $\frac{1}{10}L^2\nabla^2$, whereas that of the Yukawa Green's function is $L^2\nabla^2$; hence choosing $\alpha_Y = 1$ matches the zeroth-order homogeneous [$\nabla^2\rho \to 0$] limit, while matching the gradient-sensitive second moment requires $L_{\text{Yukawa}} = L_{\text{top-hat}}/\sqrt{10}$ or rescaling $\alpha_Y$ — see `step10`.)*

---

## 2. Step-by-Step Summary of Hungarian Files (`step1`–`step22`, `README.md`, `MANIFEST.txt`) and English Files (`step1_*_en.md`–`step22_*_en.md`, `README_en.md`, `MANIFEST_en.txt`)

### 2.1. Hungarian Derivation Chapters (`step1_*.md`–`step22_*.md`)

| Step | File Name | Size (B) | SHA-256 (prefix) | One-Line Gist (HU / EN) |
|---:|---|---:|---|---|
| **1** | `step1_dilaton_form.md` | 11 20511 205 ecd894436341cbdcecd894436341cbdc
| **2** | `step2_conformal_disformal.md` | 18 40218 402 7a3761a4adcca6fd7a3761a4adcca6fd
| **3** | `step3_yukawa_helmholtz_eom.md` | 12 64712 647 77fa5ad7646bad0b77fa5ad7646bad0b
| **4** | `step4_covariant_conservation_hydrostatic.md` | 12 35612 356 52b1cab653da99d652b1cab653da99d6
| **5** | `step5_coupling_dictionary_ppn_screening.md` | 14 94914 949 d9c9d69ff93e6625d9c9d69ff93e6625
| **6** | `step6_sparc_galactic_ngc3198_consistency.md` | 14 61514 615 8fbdd2157d0a4aa98fbdd2157d0a4aa9
| **7** | `step7_neutron_star_tov_stability.md` | 21 21521 215 d656ac9249a495aed656ac9249a495ae
| **8** | `step8_microscopic_clock_coupling_siriusb_lab_v20.md` | 16 30616 306 f4d790f47b5ab8d4f4d790f47b5ab8d4
| **9** | `step9_vhe_photon_propagation_hess_dispersion.md` | 12 16412 164 af5b8d0ab66ff02caf5b8d0ab66ff02c
| **10** | `step10_bilateral_kernel_reciprocity_eft_synthesis.md` | 13 42013 420 3cfe4aaa90e66f3f3cfe4aaa90e66f3f
| **11** | `step11_chameleon_screening_multiscale_eft_synthesis.md` | 13 29713 297 69c908b620c7562669c908b620c75626
| **12** | `step12_second_variation_coupled_pulsation_hessian.md` | 19 54519 545 66e0fdb1964263f166e0fdb1964263f1
| **13** | `step13_curved_window_metric_variation_synge_world_function.md` | 11 59511 595 261d16eadbf89cf4261d16eadbf89cf4
| **14** | `step14_gravitational_lensing_null_geodesics_cosmological_background.md` | 16 32416 324 d377321aba8ddfffd377321aba8ddfff
| **15** | `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md` | 13 34013 340 75332aacb1edc9f375332aacb1edc9f3
| **16** | `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md` | 14 17714 177 54fd2d65267c463354fd2d65267c4633
| **17** | `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md` | 14 38314 383 fce1edf9f49c8d62fce1edf9f49c8d62
| **18** | `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance.md` | 21 17821 178 6ba1c66a313683626ba1c66a31368362
| **19** | `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy.md` | 18 88218 882 bf6eea32f1a730c0bf6eea32f1a730c0
| **20** | `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature.md` | 17 36717 367 bedb9d8511f3c287bedb9d8511f3c287
| **21** | `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse.md` | 18 43118 431 29eb46950ed70c5729eb46950ed70c57
| **22** | `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md` | 35 77131 713 0c09a01ca8973d539013d9fd621f9235

### 2.2. English Translation Files (`step1_*_en.md`–`step22_*_en.md`, `README_en.md`, `MANIFEST_en.txt`)

| Step | English File Name | Size (B) | SHA-256 (prefix) | One-Line Description |
|---:|---|---:|---|---|
| **1 (EN)** | `step1_dilaton_form_en.md` | 10 89010 890 83a2e8e88808373983a2e8e888083739
| **2 (EN)** | `step2_conformal_disformal_en.md` | 18 01418 014 abd9db8d104e3bcfabd9db8d104e3bcf
| **3 (EN)** | `step3_yukawa_helmholtz_eom_en.md` | 12 25812 258 8ca1db6f110e26178ca1db6f110e2617
| **4 (EN)** | `step4_covariant_conservation_hydrostatic_en.md` | 11 98611 986 1deaba07182e7cfa1deaba07182e7cfa
| **5 (EN)** | `step5_coupling_dictionary_ppn_screening_en.md` | 14 70114 701 e57d506cdd3ba02be57d506cdd3ba02b
| **6 (EN)** | `step6_sparc_galactic_ngc3198_consistency_en.md` | 13 87013 870 6f75c1f0ccba91de6f75c1f0ccba91de
| **7 (EN)** | `step7_neutron_star_tov_stability_en.md` | 20 54620 546 e219a8cfaf4c81a2e219a8cfaf4c81a2
| **8 (EN)** | `step8_microscopic_clock_coupling_siriusb_lab_v20_en.md` | 15 62015 620 b6caf38ba27692eab6caf38ba27692ea
| **9 (EN)** | `step9_vhe_photon_propagation_hess_dispersion_en.md` | 11 45411 454 dabd6f01ef17584bdabd6f01ef17584b
| **10 (EN)** | `step10_bilateral_kernel_reciprocity_eft_synthesis_en.md` | 12 89012 890 748620c6254121e1748620c6254121e1
| **11 (EN)** | `step11_chameleon_screening_multiscale_eft_synthesis_en.md` | 12 50812 508 375040d0e3dc34c7375040d0e3dc34c7
| **12 (EN)** | `step12_second_variation_coupled_pulsation_hessian_en.md` | 19 01019 010 798f21c32075f3a7798f21c32075f3a7
| **13 (EN)** | `step13_curved_window_metric_variation_synge_world_function_en.md` | 11 23311 233 dddf81914b76e505dddf81914b76e505
| **14 (EN)** | `step14_gravitational_lensing_null_geodesics_cosmological_background_en.md` | 15 51515 515 43051f7ff852c02543051f7ff852c025
| **15 (EN)** | `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20_en.md` | 12 94712 947 8b537673b566d4cb8b537673b566d4cb
| **16 (EN)** | `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress_en.md` | 13 53013 530 f3886101a664595bf3886101a664595b
| **17 (EN)** | `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness_en.md` | 13 76513 765 672dc1adfa28d46f672dc1adfa28d46f
| **18 (EN)** | `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance_en.md` | 20 83220 832 ad7267dc00729932ad7267dc00729932
| **19 (EN)** | `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy_en.md` | 18 06718 067 9fab3349db68588b9fab3349db68588b
| **20 (EN)** | `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature_en.md` | 16 82016 820 800ff0bb69664079800ff0bb69664079
| **21 (EN)** | `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse_en.md` | 17 54217 542 106b9ee27d97f4a4106b9ee27d97f4a4
| **22 (EN)** | `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability_en.md` | 36 89631 863dc84080449d2fcbac89fccc6bcfe0ffe
| **README (EN)** | `README_en.md` | — | — | Full English edition of `README.md` listing all Hungarian and English chapters, numerical tables, and reproduction instructions. |
| **MANIFEST (EN)** | `MANIFEST_en.txt` | — | — | Full English manifest listing all tracked files (including the 22 English files) with exact SHA-256 checksums and byte sizes. |

---

## 3. Summary of Key Numerical Results

| Physical Domain | Quantity | Exact Numerical Value | Source Step & Script |
|---|---|---|---|
| **Scalar–Tensor / Yukawa Dictionary (`Branch II-Conf`)** | Canonical Yukawa strength ($\alpha_Y = 1$) and Brans–Dicke $\omega_{\text{BD}}$ | $\alpha_{\text{std}} = 2\beta^2 = 2\tilde{\alpha}_L = \frac{\alpha_Y}{4\pi} = \frac{1}{4\pi} \approx 0.0795775$; $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2} \approx 11.06637 \approx 11.07$ ($\alpha_Y=1$) or $36.1991$ ($\alpha_Y=1/3$) | `step5`, `covariant_closure_50aa2db0.py` |
| **SPARC 175 Galaxies ($\rho < 0$ positivity test)** | Number of negative-density points (across 3391 radial points) | **0 points** (positivity check of photometric + 3D top-hat + 3D Yukawa kernels) vs. **239 points** (140 non-holdout galaxies, $\Upsilon_d=0.5$) / **301–303 points** (175 galaxies, spherical derivative formula) | `step6`, `SPARC_mode.py`, `sparc_tophat.py` |
| **NGC 3198 Negative Control ($R = 30.08\;\text{kpc}$)** | Observed $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$ vs. isolated constant-coupling $C_{\text{pred}}^{\text{isol}}$ | $V_{\text{obs}} = 146 \pm 2\;\text{km/s}$; $\Upsilon_d=1.0$: $V_{\text{bar}} = 88.53\;\text{km/s}$, **$C_{\text{req}} = 2.720$**, whereas for the isolated constant-coupling model **$C_{\text{pred}}^{\text{isol}} \le 1 + \alpha_{\text{std}} = 1.0796$** ($V_{\text{pred}}^{\text{isol}} = 88.53\text{--}91.99\;\text{km/s}$; $\Upsilon_d=0.5$: $V_{\text{bar}}=70.32\;\text{km/s}$, $C_{\text{req}}=4.311$, $V_{\text{pred}}^{\text{isol}} = 70.32\text{--}73.06\;\text{km/s}$) | `step6`, `SPARC_mode.py`, `data/NGC3198_table.csv` |
| **Cosmological Horizon Scales ($\Omega_b$ vs. $\Omega_m$)** | $L_{\text{hor}} = \frac{c}{\sqrt{G\rho_{\text{bg}}}} = \sqrt{\frac{8\pi}{3\Omega}}\,\frac{c}{H_0}$ ($H_0 = 67.4\;\text{km/s/Mpc}$) | **Baryonic background ($\Omega_b = 0.0493$, $\rho_{b,0} = 4.2067\times 10^{-28}\;\text{kg/m}^3$):** $\sqrt{8\pi/(3\Omega_b)} = 13.031 \implies L_{\text{hor}}^{(b)} = 1.7892\times 10^{27}\;\text{m} \approx \mathbf{57.98\;\text{Gpc}}$;<br>**Matter background ($\Omega_m = 0.315$, $\rho_{m,0} = 2.686\times 10^{-27}\;\text{kg/m}^3$):** $\sqrt{8\pi/(3\Omega_m)} = 5.160 \implies L_{\text{hor}}^{(m)} = 7.081\times 10^{26}\;\text{m} \approx \mathbf{22.94\;\text{Gpc}}$ | `step14`, `sparc_tophat.py` |
| **Bianchi / Noether Conservation** | $\max_r |r\cdot\nabla_\mu(T^{\mu}_{\;r,\text{matter}} + T^{\mu}_{\;r,X})|/\varepsilon_c$ ($s=r^2$ PCHIP) | Symbolic: `0.0`; continuum: $1.82\times 10^{-11}$ (`tophat_mass`), $1.59\times 10^{-11}$ (`tophat_gradient`), $3.65\times 10^{-10}$ (`yukawa_action_el`) | `step4`, `step16`, `step21`, `covariant_closure_50aa2db0.py` |
| **Curved Exterior Yukawa BVP and ADM Mass (`Branch II-Fluid`)** | $R_{\max}$ convergence ($[R_*, R_*+N_{\text{mult}}L]$, $N_{\text{mult}}\in\{10,12,14,16\}$), $k_{\text{curved}}$ vs. $k_{\text{flat}}$, and $M_{\text{ADM},\max}$ | $R_{\max}$ spread: $\Delta k_{\text{curved}}/k_{\text{curved}} = 6.17\times 10^{-12}$, $\Delta M_{\text{ADM}}/M_{\text{ADM}} = 5.54\times 10^{-13}$; SLy $15.0$, $L=10\;\text{km}$, $\alpha_Y=1$: $k_{\text{curved}} = 0.233542\;\text{km}^{-1}$ vs. $k_{\text{flat}} = 0.187216\;\text{km}^{-1}$ (**$+24.74\%$**); $M_{\text{ADM},\max}^{\text{SLy}} = \mathbf{2.03321\,M_\odot}$ ($M_{\text{surf},\max} = 2.02911\,M_\odot$, $M_{B,\max}^{\text{quad}} = \mathbf{2.42339\,M_\odot}$, $\Delta M_{\phi,\text{ext}}^{\text{curved}} = 0.00411\,M_\odot$), $R_{1.4}^{\text{ADM}} = \mathbf{11.5359\;\text{km}}$ ($M_{\text{ADM},\max}^{\text{APR4}} = 2.18068\,M_\odot$, $M_{\text{ADM},\max}^{\text{H4}} = 1.97809\,M_\odot$) | `step7`, `step22`, `curved_yukawa_adm_and_stability.py` |
| **Coupled $(\delta h,\delta\phi_L,\delta m)$ Stability on the Curved `yukawa_action_el` Branch** | $\omega_0^2 = 0$ zero crossing vs. $\frac{dM_{\text{ADM}}}{d\rho_c} = 0$ and $\frac{dM_B}{d\rho_c} = 0$ (`scalar_m_sign = -1.0`) | SLy $L=10\;\text{km}$, $\alpha_Y=1$: $\log_{10}\rho_{0c}(M_{\text{ADM},\max}) = 15.316404$, $\log_{10}\rho_{0c}(M_{B,\max}^{\text{quad}}) = 15.316409$, coupled $\omega_0^2=0$: **$15.316410$** ($\mathbf{+5.80\times 10^{-6}\;\text{dex}}$, **$2097\times$** reduction over frozen-$\phi_L$ shift $+1.22\times 10^{-2}\;\text{dex}$); $\alpha_Y=1/3$: **$15.307170$** (**$+4.20\times 10^{-6}\;\text{dex}$**, **$995\times$** reduction); $L=1\;\text{km}$: **$15.303080$** (**$-6.86\times 10^{-6}\;\text{dex}$**, **$165\times$** reduction); APR4 ($L=10\;\text{km}$): **$15.293025$** (**$+3.57\times 10^{-6}\;\text{dex}$**, **$2392\times$** reduction); H4 ($L=10\;\text{km}$): **$15.227783$** (**$-1.34\times 10^{-6}\;\text{dex}$**, **$17593\times$** reduction) | `step12`, `step22`, `curved_yukawa_adm_and_stability.py` |
| **V21 Nonlocal 3D Top-Hat Tangent BVP (`Branch I-A`)** | SLy $L = 10\;\text{km}$ `mass` / `gradient` $\omega_0^2 = 0$ vs. $\partial M/\partial\rho_c = 0$ and $W\to 4W$ scaling | In `mass` mode: $15.29931$ vs. $15.29924$ ($+6.97\times 10^{-5}\;\text{dex}$, $58.3\times$ improvement); in `gradient` mode: $15.30113$ vs. $15.30106$ ($+7.39\times 10^{-5}\;\text{dex}$, $59.3\times$ improvement); inertia scaling: $\omega^2_{4W}/\omega^2_W = 0.2500$ | `step7`, `step17`, `step18`, `covariant_closure_50aa2db0.py` |
| **Sirius B H0 vs. H1 Redshift** | $v_{\text{obs}} = 80.65\pm 0.77\;\text{km/s}$ ($M=1.018\,M_\odot$, $R=5846.6\;\text{km}$) | H0: $|\Delta v|\le 1.27\;\text{m/s}$ (unconstrained); conditional H1 $2\sigma$ ($\Delta\chi^2=4$) excludes: $L\in[4.08\times 10^6, 3.37\times 10^7]\;\text{m}$ ($n=1.5$, $\max\Delta\chi^2=52.3$ / $72.0$) | `step8`, `step19`, `covariant_closure_50aa2db0.py` |
| **Granite Chamber ($8558\;\text{kg}$) IN–OUT Test** | Conditional H1 differential clock signal $\Delta y_{\text{H1}}(\text{IN}-\text{OUT})$ | At $L=1.00\;\text{m}$: Top-hat $-1.436\times 10^{-24}$ ($\equiv 0$ for $L\ge 2.944\;\text{m}$); Yukawa $-3.781\times 10^{-25}$ ($-4.818\times 10^{-25}$ at $L=10\;\text{m}$) | `step8`, `step19`, `covariant_closure_50aa2db0.py` |
| **8-Coil Solenoid ($N=2000$, $14\;\text{mA}$)** | Magnetic peak signal ($L_p=0.630957\;\text{m}$), $r=L/t$ quadrature, and $L=0.95\;\text{m}$ | $\Delta W(L_p)=1.96114452\;\text{J/A}^2$, $k_g = -6.13117767\times 10^{-45}\;\text{A}^{-2}$, $\Delta y(L_p) = -1.2017108\times 10^{-48} \sim -10^{-48}$; Path A vs. B ($L=3\;\text{m}$): $7.44\times 10^{-16}$; $L=0.95\;\text{m}$ identity: $8.88\times 10^{-16}\;\text{J/A}^2$ | `step15`, `step20`, `make_revision3_checks.py` |

---

### 3.1. V25 Fixes and Diagnostics: 4-Component Exterior BVP, Exact Exterior Linearization, Dynamical Scalar Wave, and GL2 N-Convergence

*(Source: `curved_yukawa_adm_and_stability.py`, `curved_yukawa_audit_summary.json`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_adm_sequences.csv`.)*

- **Exterior Lapse (Issue 5):** the exterior $\nu$ is co-integrated as the fourth ODE component $(\phi_L,\psi_L,m,\nu)$ with `solve_bvp` `tol=1e-12` and outer boundary condition $\nu(R_{\max})=\tfrac12\ln(1-2m(R_{\max})/R_{\max})$. Spreads for SLy, $\log_{10}\rho_{0c}=15.0$, $L=10\;\text{km}$, $N_{\text{mult}}\in\{10,12,14,16\}$: $\delta k/k=6.17\times 10^{-12}$, $\delta M_{\text{ADM}}/M_{\text{ADM}}=5.54\times 10^{-13}$, $\delta\nu_R=1.63\times 10^{-14}$.

| $N_{\text{mult}}$ | $R_{\max}$ (km) | $k_{\text{curved}}$ (km$^{-1}$) | $M_{\text{ADM}}$ ($M_\odot$) | $\nu_R$ |
|---|---|---|---|---|
| `10` | `111.4657` | `0.233542329280` | `1.509998138706` | `-0.246306687228` |
| `12` | `131.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `14` | `151.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `16` | `171.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |

- **Exterior Linearization and Finite-Difference Derivative Tests (SLy, $\log_{10}\rho_{0c}=15.0$, $L=10\;\text{km}$, $\alpha_Y=1$):** `linearized_exterior_response` (line 593: `a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0`) and the left-hand side ($L_u u = R_h\,\delta h$) interior scalar operator coefficient `a_m_sc` (`scalar_m_sign = -1.0`: `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n`) verified against central finite differences:
  - Exterior response $\partial\phi_L'(R_*)/\partial m(R_*)$: linearized BVP **`-1.68365059e-04`** $\text{km}^{-2}$ vs. central finite difference of the full 4-component exterior BVP **`-1.68365059e-04`** $\text{km}^{-2}$ (relative difference: **$3.98\times 10^{-10}$**).
  - Interior $r = 6.0\;\text{km}$ scalar–mass derivative ($\partial\phi_L''/\partial m$): analytical **`-1.72833429e-04`** $\text{km}^{-3}$ vs. central finite difference of the interior scalar RHS **`-1.72833428e-04`** $\text{km}^{-3}$ (relative difference: **$6.08\times 10^{-9}$**).
- **Dynamical Scalar Block:** the scalar block includes the $-\omega^2 e^{-2\nu}\delta\phi_L$ term and the exterior response is supplied by the $\omega$-dependent BVP (`dynamical_scalar_wave=True`); solver: GL2 piecewise-segment Volterra discretization, $N=480$, with exact LHS sign `scalar_m_sign = -1.0`. The quasi-static and dynamical roots agree to all 6 reported decimal places.
- **Turning Points and Zero Crossings:** `minimize_scalar` + 7-point quadratic fit over $\pm0.002$ dex; zero crossings from a 3-point quadratic fit over $\pm0.0005$ dex around the root.

| Sequence | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dynamical $\omega_0^2=0$ | $\Delta_{\text{ADM}}$ Dynamical (dex) | $\Delta_{M_B}$ Dynamical (dex) | Quasi-Static $\omega_0^2=0$ | Frozen $\delta\phi_L=0$ | $\Delta_{\text{ADM}}$ Frozen (dex) |
|---|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316409` | **`15.316410`** | **`+5.80e-06`** | `+1.33e-06` | `15.316410` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.307169` | **`15.307170`** | **`+4.20e-06`** | `+6.41e-07` | `15.307170` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303082` | **`15.303080`** | **`-6.86e-06`** | `-2.61e-06` | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.293025` | **`15.293025`** | **`+3.57e-06`** | `-1.40e-07` | `15.293025` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227787` | **`15.227783`** | **`-1.34e-06`** | `-4.05e-06` | `15.227783` | `15.251384` | `+2.36e-02` |

**GL2 Convergence on the $\omega_0^2=0$ Root (across $N$ and $R_{\max}$) and at $\log_{10}\rho_{0c}=15.311404$** (SLy, $L=10\;\text{km}$, $\alpha_Y=1$): on the root, $p_{\text{eff,root}}=4.248$, Richardson root $\log_{10}\rho_{0c,\text{Rich}}(\omega_0^2=0)=`15.31640981`$ ($N=480$ discretization error $6.72\times 10^{-8}$ dex; spread across $R_{\max}\in[R_*+10L,R_*+20L]$ is $2.42\times 10^{-13}$ dex); at $\log_{10}\rho_{0c}=15.311404$, $p_{\text{eff}}=5.223$, Richardson value $\omega_0^2=6369417.611\;\text{s}^{-2}$ ($N=480$ relative error $1.29\times 10^{-5}$).

| $N$ | $\log_{10}\rho_{0c}(\omega_0^2=0)$ | $|\Delta\log_{10}\rho_{0c}|$ vs. Richardson Root (dex) | $\Delta_{\text{ADM}}$ (dex) | $\Delta_{M_B}$ (dex) | $\omega_0^2(15.311404)$ ($\text{s}^{-2}$) | Rel. Diff. ($\omega_0^2$) |
|---|---|---|---|---|---|---|
| `80` | `15.31645808` | `4.83e-05` | `+5.40e-05` | `+4.95e-05` | `6432300.503336` | `9.87e-03` |
| `120` | `15.31642223` | `1.24e-05` | `+1.82e-05` | `+1.37e-05` | `6385729.603235` | `2.56e-03` |
| `160` | `15.31641426` | `4.46e-06` | `+1.02e-05` | `+5.72e-06` | `6375276.486537` | `9.20e-04` |
| `240` | `15.31641079` | `9.80e-07` | `+6.71e-06` | `+2.24e-06` | `6370703.105233` | `2.02e-04` |
| `360` | `15.31641001` | `2.06e-07` | `+5.94e-06` | `+1.47e-06` | `6369680.451299` | `4.13e-05` |
| `480` | `15.31640987` | `6.72e-08` | `+5.80e-06` | `+1.33e-06` | `6369499.539632` | `1.29e-05` |
| `600` | `15.31640983` | `2.56e-08` | `+5.76e-06` | `+1.29e-06` | `6369446.532399` | `4.54e-06` |
| `800` | `15.31640982` | `7.55e-09` | `+5.74e-06` | `+1.27e-06` | `6369424.048161` | `1.01e-06` |
| `1200` | `15.31640981` | `1.35e-09` | `+5.73e-06` | `+1.26e-06` | `6369418.385568` | `1.22e-07` |

**Interpretive Note (Boundary of Theoretical Closure and Conditional Stability Closing Statement):** With the exact left-hand side sign (`scalar_m_sign = -1.0`), the root of the GL2 operator converges with respect to $N$ and $R_{\max}$ to the reported precision ($<10^{-7}\;\text{dex}$, `15.31640981`), and lies within $\le 6.9\times 10^{-6}\;\text{dex}$ of the ADM mass maximum across all 5 sequences (for SLy $L=10\;\text{km}$, $\alpha_Y=1$: $15.316410$ vs. $15.316404$, $+5.80\times 10^{-6}\;\text{dex}$; the diagnostic comparison with `scalar_m_sign = +1.0` is documented in Section 7.4 of `step22`). For SLy $L=10\;\text{km}$, $\alpha_Y=1$, the relative $L^2$ differences between the eigenfunction and the static tangent are $\delta h = 5.04\times 10^{-4}$, $\delta m = 1.77\times 10^{-4}$, $\delta\phi_L = 2.27\times 10^{-4}$, $\delta\phi_L' = 3.09\times 10^{-4}$, while the relative difference of the free-surface displacements is $3.54\%$ (`rel_diff_xi_R_exact_minus_vs_static_tangent = 0.035406`). While the minus-sign root and grid convergence represent a major step forward, because of the $3.54\%$ free-surface difference and the ODE vs. quadrature baryonic-mass turning-point spread they do not by themselves constitute an independent proof of the full radial fundamental spectrum and free-surface variation; therefore, the stability closing statement is formulated conditionally until the independent verification of the coupled variational derivation and boundary-condition matching is concluded.

---

## 4. Local Archival Manuscripts Used (`manuscripts/`)

All referenced V19–V23 LaTeX (`.tex`) source files are included in the `manuscripts/` subdirectory:

1. **`manuscripts/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`** and **`manuscripts/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex`**  
   - *Version & DOI:* V23 (`50aa2db0`), Zenodo DOI `10.5281/zenodo.23639084`.  
   - *Role in the Derivation:* Primary source for the two covariant action candidates (Formulation I bilateral kernel action and Formulation II local Helmholtz–Yukawa scalar–tensor action), Bianchi residuals, coupled radial stability BVP (`eq:factored_diff`, `eq:coupled_bvp`), and H0/H1 clock coupling (Sirius B and $8558\;\text{kg}$ granite chamber).
2. **`manuscripts/hu_p1.tex`**, **`manuscripts/hu_p2.tex`**, **`manuscripts/hu_p3.tex`**, **`manuscripts/en_part1.tex`**, **`manuscripts/en_part2.tex`**  
   - *Version & DOI:* V22 Book Edition chapter sources, Zenodo DOI `10.5281/zenodo.23604423`.  
   - *Role in the Derivation:* Source of verbatim citations for the Cauchy derivation (Chapter 2), 3D top-hat and SPARC 175-galaxy analysis (Chapters 4–5), V17 H.E.S.S./MAGIC gamma-ray timing audit (Chapter 6), 7-segment piecewise-polytrope TOV and $s=r^2$ regularization (Chapter 7), coupled radial pulsation (Chapter 8), Sirius B redshift (Chapter 9), $8558\;\text{kg}$ granite chamber and 8-coil solenoid experiment (Chapter 10), and the numerical Appendix.
3. **`manuscripts/tal_book_HU_2026-10-24.tex`**, **`manuscripts/tal_book_EN_2026-10-24.tex`**, **`manuscripts/tal_book_DE_2026-10-24.tex`**  
   - *Version & DOI:* V22 Book Edition unified Hungarian, English, and German LaTeX monograph files, Zenodo DOI `10.5281/zenodo.23604423`.
4. **`manuscripts/sparc_tophat_record_2026-10-22_v5.tex`**  
   - *Version & DOI:* V19 SPARC 3D spherical top-hat manuscript, Zenodo DOI `10.5281/zenodo.23497128`.
5. **`manuscripts/tal_v20_fop.tex`**  
   - *Version & DOI:* V20 laboratory two-clock and 8-coil magnetic manuscript (*Foundations of Physics* format), Zenodo DOI `10.5281/zenodo.23540741`.
6. **`manuscripts/tal_v21_fop.tex`**  
   - *Version & DOI:* V21 compact stellar TOV and Sirius B manuscript, Zenodo DOI `10.5281/zenodo.23602501`.

---

## 5. Python Programs & SPARC Data Files

### 5.1. Python Programs (in the Root Directory and `code/` Subdirectory)
1. **`curved_yukawa_adm_and_stability.py`** — Module solving the coupled interior + curved exterior ($r \in [R_*, R_*+12L]$) Einstein–scalar boundary-value problem of `yukawa_action_el` (**Branch II-Fluid**), verifying $R_{\max}$-independence, curved vs. flat Robin slope, recomputed mass–radius sequences in $M_{\text{ADM}}$ and $M_B$ (`curved_yukawa_adm_sequences.csv`), and coupled $(\delta h, \delta\phi_L, \delta m)$ dynamical radial pulsation eigenvalues (`curved_yukawa_coupled_stability.csv`, `curved_yukawa_audit_summary.json`).
2. **`covariant_closure_50aa2db0.py`** — Main computational and symbolic verification module of V23 (`10.5281/zenodo.23639084`: SymPy Bianchi identities, Formulation I and II TOV comparison — including `yukawa_action_el_curved_adm` rows —, coupled `coupled_bvp` / `sym_galerkin` radial stability, $W\to 4W$ scaling invariance, Sirius B and $8558\;\text{kg}$ granite chamber H0/H1 clock coupling).
3. **`radial_modes.py`** — Finite-element / banded-matrix eigenvalue solver for the Chandrasekhar Sturm–Liouville radial pulsation operator (`std`, `modA`, `modB`, and `fluid` vs. `metric` inertia weight).
4. **`tal_profiles.py`** — Profile module implementing $s = r^2$ even-parity central PCHIP interpolation (`Audit Fix O1`), the 3D top-hat shell integral, and the `gradient`-mode exterior vacuum-halo lapse correction $\Delta\nu_{\text{ext}}$ on $r\in[R_*, R_*+L]$.
5. **`stellar_tal_v21.py`** (and `code/stellar_tal.py`, `code/diagnostics.py`) — V21/V23 module containing the Read et al. (2009) 7-segment piecewise-polytrope equations of state (SLy, APR4, H4), the log-enthalpy ($h$) TOV integrator (769 models, 33 causal endpoints), and the white-dwarf (Sirius B) integrator.
6. **`make_revision3_checks.py`** (along with `code/energy_integrator.py`, `code/fields.py`, `code/tophat_finite_source.py`, `code/v20_fixed_case_runner.py`, `code/build_pdf.py`) — Scripts running the V20 (`10.5281/zenodo.23540741`) 8-coil solenoid elliptic Biot–Savart check, $L_p = 0.630957\;\text{m}$ peak recalculation, $r=L/t$ compactified exterior integral (Paths A and B), and $L=0.95\;\text{m}$ interior/exterior switching identity.
7. **`sparkto_b_raw.py`** and **`sparc_tophat.py`** — V19 program computing the 3D spherical top-hat mean density ($\bar{\rho}_L(\mathbf{x})$) and $f_L = \exp(-\gamma_L\bar{\rho}_L)$ factor of the 3D axisymmetric baryonic disk+gas+bulge model for all 175 SPARC galaxies.
8. **`SPARC_mode.py`**, **`ta_pipeline.py`**, **`ta_foundation.py`**, **`sparc_full_runner.py`**, **`tal_factor_7th_model.py`**, **`ta_gammaL.py`**, **`ta_run.py`**, **`ta_run_l0.py`**, **`ta_run_l0_v2.py`** — Programs running the SPARC 175 / 155 / 123 galaxy photometric vs. derivative-based (`foundation`) density comparison (0 vs. 239 / 301–303 negative-density points), the NGC 3198 ($R=30.08\;\text{kpc}$) control point, and V17/V19 7th-model rotation-curve fits.

### 5.2. SPARC and V19–V23 Data Files (`data/`, `csv/`, `run_output/`)
- **`data/Rotmod_LTG.zip`** (SHA-256: `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`) and **`data/Rotmod_LTG/*_rotmod.dat`** — Original Lelli, McGaugh & Schombert (2016) rotation-curve data files for all 175 SPARC galaxies (including `data/Rotmod_LTG/NGC3198_rotmod.dat`).
- **`data/SPARC_Lelli2016c.mrt`** (SHA-256: `5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243`) — SPARC Table 1 (distances, inclinations, $L_{3.6}$, $R_{\text{disk}}$, $M_{\text{HI}}$, $R_{\text{HI}}$, quality flags).
- **`data/blind_123.txt`** (and `v17_package/github/data/blind_123.txt`, SHA-256: `3637d2f5bda6b2b97e645c5fba5cf39d310a1826bee388286c2c74b9171aa283`) — Hash-locked 123-galaxy blind sample list of V17/V19.
- **`data/NGC3198_table.csv`** and **`data/NGC3198_table.md`** — Detailed table of all 43 radial points of NGC 3198 (highlighting the $R=30.08\;\text{kpc}$ control point).
- **`data/sparc_tophat_per_galaxy_max.csv`** and **`data/sparc_tophat_summary.json`** — V19 3D spherical top-hat peak values and Monte Carlo cross-checks across all 175 SPARC galaxies.
- **`data/v19_7th_model/`** and **`run_output/`** — Reference outputs of the V17/V19 baryonic baselines (`points_vbar.csv`, `summary_vbar.csv`) and the 7th model (`points_tal_fL_FITTED_forward.csv`, `fit_scan_L.csv`, `comparison.json`, `run_meta.json`).
- **`csv/`** and **`data/v23/`** — Reference files of the V20 laboratory granite chamber and 8-coil solenoid, together with the V23 (`50aa2db0`) and Step 22 curved exterior ADM mass and coupled stability output files (`curved_yukawa_adm_sequences.csv`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_audit_summary.json`, `conservation_Residuals.json`, `closure_tov_comparison.csv`, `coupled_stability_results.csv`, `clock_coupling_H0_H1.csv`).
- **`data/v21/`** and **`observational_inputs.csv`** — V21/V23 769 stellar models (`mr_sequences.csv`, `stellar_summary.csv`, `causal_endpoint_convergence.csv`, `profiles_L10km.csv`, `radius_shifts.csv`, `sirius_redshift.csv`, `numerical_diagnostics.json`, `wd_summary.json`) and NICER/LIGO/pulsar observational inputs.

---

## 6. Exact Reproduction Instructions

### 6.1. Software Dependencies
- **Operating System:** Linux / macOS / POSIX-compatible environment  
- **Python Version:** `Python >= 3.10` (tested under Python 3.11–3.13; on the local workstation, `/workspace/scratch/venv_openitems/bin/python3` has all scientific packages pre-installed)  
- **Required Python Packages:** `numpy`, `scipy`, `pandas`, `sympy`, `matplotlib`  
  Installation (if needed):
  ```bash
  pip install numpy scipy pandas sympy matplotlib
  ```

### 6.2. Running Each Program & Expected Outputs
Every command below can be executed directly from the unpacked `scalar-tensor-yukawa-fit` working directory without network access:

1. **Curved Exterior `yukawa_action_el` Einstein–Scalar BVP, ADM Mass ($M_{\text{ADM}}$), Mass–Radius Curves, and Coupled $(\delta h, \delta\phi_L, \delta m)$ Radial Stability (`curved_yukawa_adm_and_stability.py`):**
   ```bash
   python3 curved_yukawa_adm_and_stability.py
   ```
   - *Expected Output (~2–3 minutes):*  
     - Generates `curved_yukawa_adm_sequences.csv`, `curved_yukawa_coupled_stability.csv`, and `curved_yukawa_audit_summary.json` (and updates `data/v23/`).  
     - Verifies $R_{\max}$-independence ($\Delta k_{\text{curved}}/k_{\text{curved}} = 6.17\times 10^{-12}$, $\Delta M_{\text{ADM}}/M_{\text{ADM}} = 5.54\times 10^{-13}$), curved vs. flat Robin slope ($+24.74\%$ for SLy $\log_{10}\rho_{0c}=15.0$, $L=10\;\text{km}$, $\alpha_Y=1$), ADM mass peaks ($M_{\text{ADM},\max}^{\text{SLy},10\text{km}} = 2.03321\,M_\odot$, $R_{1.4}^{\text{ADM}} = 11.5359\;\text{km}$), finite-difference derivative tests ($\partial\phi_L'(R_*)/\partial m(R_*)=-1.68365\times 10^{-4}\;\text{km}^{-2}$, interior $r=6\;\text{km}$ `a_m_sc` $=-1.72833\times 10^{-4}\;\text{km}^{-3}$), and coupled stability zero crossings (`scalar_m_sign = -1.0`: $\log_{10}\rho_{0c}(\omega_0^2=0)=15.316410$, $\Delta_{\text{ADM}} = +5.80\times 10^{-6}\;\text{dex}$, and $-6.86\times 10^{-6}\;\text{dex}$ at $L=1\;\text{km}$).

2. **SPARC 175 Galaxies ($\rho < 0$ Audit), NGC 3198 ($R = 30.08\;\text{kpc}$) Control Point, and V19 3D Top-Hat Peaks (`SPARC_mode.py`):**
   ```bash
   python3 SPARC_mode.py --mode summary
   ```
   - *Expected Output (JSON on stdout, ~0.2 s):*  
     - `"photometric_and_3D_tophat_yukawa_rho_le_0_points": 0` vs. `"derivative_formula_rho_lt_0_points_visible140_Ydisk_0p5": 239` (and `301` across all 175 galaxies).  
     - `"NGC3198_control_point"`: `"R_kpc": 30.08`, `"V_obs_km_s": 146.0`, `"V_bar_Ydisk_1p0_km_s": 88.53`, `"C_req_Ydisk_1p0": 2.72`, `"V_bar_Ydisk_0p5_km_s": 70.32`, `"C_req_Ydisk_0p5": 4.311`.

3. **SPARC 3D Spherical Top-Hat Mean Density and $f_L$ Factor Calculation (`sparkto_b_raw.py` / `sparc_tophat.py`):**
   ```bash
   python3 sparkto_b_raw.py --galaxies NGC3198,UGC02487,UGC03546 --skip_checks --out ./output_tophat_check
   ```

4. **V20 Eight-Coil Magnetic-Field Experiment, Elliptic Biot–Savart, $r=L/t$ Exterior Quadrature, and $L=0.95\;\text{m}$ Identity (`make_revision3_checks.py`):**
   ```bash
   python3 make_revision3_checks.py . ./csv/revision3_checks_reproduced.json
   ```

5. **V23 Covariant Action Closure, Bianchi Residuals, Coupled TOV Stability BVP, Sirius B, and $8558\;\text{kg}$ Granite Chamber (`covariant_closure_50aa2db0.py`):**
   ```bash
   python3 covariant_closure_50aa2db0.py
   ```

6. **SHA-256 Integrity Verification of Steps 1–22 and All Package Files:**
   ```bash
   sha256sum -c SHA256SUMS.txt
   ```
