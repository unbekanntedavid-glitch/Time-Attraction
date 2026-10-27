# Scalar-Tensor and Yukawa/Disformal Action Mapping (V23 $\leftrightarrow$ Jordan–Brans–Dicke / Dilaton Framework) — Steps 1–22 Derivation and Audit Reproduction Package

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**Working Directory:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit`  
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

## 2. Step-by-Step Summary of Hungarian Files (`step1`–`step22`, `README.md`, `MANIFEST.txt`) and English Files (`step1_*_en.md`–`step20_*_en.md`, `README_en.md`, `MANIFEST_en.txt`)

### 2.1. Hungarian Derivation Chapters (`step1_*.md`–`step22_*.md`)

| Step | File Name | Size (B) | SHA-256 (prefix) | One-Line Gist (HU / EN) |
|---:|---|---:|---|---|
| **1** | `step1_dilaton_form.md` | 11 282 | `6d12c0fd49bb3a50` | Derives the exact equivalence between the Cauchy-derived exponential factor $f_L = \exp(-\gamma_L\bar{\rho}_L) = e^{-\phi_L}$ and the string-theory dilaton matter coupling $A(\phi)=e^{\beta\phi}$. |
| **2** | `step2_conformal_disformal.md` | 18 425 | `5a9fe423508e3329` | Derives the Jacobian, inverse metric, comoving-clock equivalence, and explicit branch taxonomy of the conformal and disformal metric frames. |
| **3** | `step3_yukawa_helmholtz_eom.md` | 12 670 | `9ca162a9375eb71c` | Variationally derives the curved-spacetime massive Klein–Gordon / Helmholtz–Yukawa field equation, Green function, and flat-space Robin approximation at $r=R_*$. |
| **4** | `step4_covariant_conservation_hydrostatic.md` | 12 379 | `8f5446c737f3c39f` | Proves exact covariant stress-energy conservation $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$ and the modified hydrostatic Euler–TOV balance $h' = -\nu' + \phi_L'$. |
| **5** | `step5_coupling_dictionary_ppn_screening.md` | 14 972 | `33ce437c1eec32ba` | Establishes the exact $-\frac{c^4}{2G\alpha_Y}$ normalized coupling dictionary ($\alpha_{\text{std}} = 1/(4\pi)$, $\omega_{\text{BD}} \approx 11.07$) and separates conformal vs. disformal PPN behavior. |
| **6** | `step6_sparc_galactic_ngc3198_consistency.md` | 14 692 | `1d3c17156dbf5ed4` | Verifies 0 negative-density points on 175 SPARC galaxies and contrasts NGC 3198's required $C_{\text{req}}=2.720$ against the isolated constant-coupling prediction $C_{\text{pred}}^{\text{isol}} \le 1.0796$ as a negative control. |
| **7** | `step7_neutron_star_tov_stability.md` | 20 067 | `ce336cc89fa20b2b` | Derives the curved-spacetime TOV equations, compares flat-Robin surface mass vs. curved-exterior ADM mass, and separates Branch I-A vs. Branch II-Fluid stability. |
| **8** | `step8_microscopic_clock_coupling_siriusb_lab_v20.md` | 16 329 | `8e3645451604bbc8` | Derives H0 vs. conditional H1 clock coupling from the Dirac–Maxwell action and quantifies Sirius B redshift exclusions, the $8558\;\text{kg}$ granite chamber ($\sim 10^{-24}$), and the 8-coil solenoid ($\sim 10^{-48}$). |
| **9** | `step9_vhe_photon_propagation_hess_dispersion.md` | 12 187 | `9e735e9bfd15e9bc` | Derives the achromatic refractive index $n_t(s)=f_L(s)^{-1}$ in the disformal metric and the multi-source H.E.S.S./MAGIC VHE gamma-ray timing audit. |
| **10** | `step10_bilateral_kernel_reciprocity_eft_synthesis.md` | 13 443 | `4e336d901146cfea` | Derives Formulation I's auxiliary-field bilateral kernel action, reciprocal potential $V_{\text{eff},N}(r)$, and zeroth- vs. second-moment correspondence to Formulation II. |
| **11** | `step11_chameleon_screening_multiscale_eft_synthesis.md` | 13 320 | `3d7262216b373998` | Formulates environment-dependent chameleon/symmetron screening (`Branch II-Cham`) as an explicit nonlinear extension hypothesis beyond the constant-$\alpha_Y$ V23 action. |
| **12** | `step12_second_variation_coupled_pulsation_hessian.md` | 18 743 | `1c2f9401c56e140b` | Explicitly derives $\mathcal{Q}_{\xi\xi}, \mathcal{Q}_{\phi\phi}, \mathcal{C}_{\xi\phi}, \mathcal{C}_{\phi\xi}$ after metric constraint elimination, moving-surface matching, and numerical stability on the curved `yukawa_action_el` sequence. |
| **13** | `step13_curved_window_metric_variation_synge_world_function.md` | 11 618 | `ec8dcd74bd038ec2` | Derives the exact metric variation $\delta[\sqrt{h}K_L]/\delta g^{\mu\nu}$ of the curved-spacetime top-hat window via Synge's world function and the $\delta_D(L-s_{\text{geo}})$ surface stress. |
| **14** | `step14_gravitational_lensing_null_geodesics_cosmological_background.md` | 16 347 | `05725c4639528c5f` | Derives conformal vs. disformal light deflection, FLRW cancellation in $\nabla\ln f_L$, and exact horizon scales for both $\Omega_b=0.0493$ ($57.98\;\text{Gpc}$) and $\Omega_m=0.315$ ($22.94\;\text{Gpc}$). |
| **15** | `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md` | 13 386 | `ff2ed9a3ea6f3080` | Proves that while $T^\mu_{\;\mu,\text{EM}}=0$ suppresses conformal trace coupling, disformal coupling $u_\mu u_\nu T^{\mu\nu}_{\text{EM}}$ sources $\phi_L$ via $\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$ with exact polarity-even cross-term cancellation. |
| **16** | `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md` | 14 200 | `0af30fe0290528d5` | Derives the V21/V23 algebraic X-sector anisotropic stress reconstruction from the Bianchi identity ($\max|p_{tX}/P_c|=0.140$) and compares it with the dynamical Yukawa stress tensor. |
| **17** | `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md` | 14 406 | `b94fc6996982b636` | Analytically derives the V23 Eq. (8) factorized Chandrasekhar–hydrostatic difference identity and proves the non-positivity of the `gradient`-mode `metric` inertia weight for $L \ge 367\;\text{m}$. |
| **18** | `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance.md` | 20 862 | `0ebe52095efe2a94` | Derives the V23 Eq. (9) nonlocal static tangent BVP on the V21 background, the 3D spherical-shell kernel $K_{ij}$, and the exact inertia-scaling invariance $\omega^2_{4W}/\omega^2_W = 0.2500$. |
| **19** | `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy.md` | 18 905 | `4aeec2bc59136cb8` | Derives the 3D IN–OUT differential quadrature of the $8558\;\text{kg}$ granite chamber across three geometric $L$-regimes and the Sirius B–Sirius A `common_background` vs. `fixedM` redshift proxy. |
| **20** | `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature.md` | 17 390 | `f373f8084f4257e2` | Derives the 8-coil elliptic Biot–Savart field, peak signal $\Delta y(L_p) = -1.2017\times 10^{-48}$, $r=L/t$ compactified exterior quadrature ($7.44\times 10^{-16}$ agreement), and the $L=0.95\;\text{m}$ switching identity. |
| **21** | `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse.md` | 18 011 | `94ff920e1b578faf` | Derives the 7-segment piecewise-polytrope log-enthalpy TOV transformation, $s=r^2$ even-parity PCHIP regularization, exterior vacuum-halo lapse $\Delta\nu_{\text{ext}}$, and contrasts V21 $h'=-\nu'$ against Yukawa $h'=-\nu'+\phi_L'$. |
| **22** | `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md` | 16 654 | `dce145dee013eb92` | Solves the curved exterior Einstein–scalar vacuum BVP of `yukawa_action_el`, determines $M_{\text{ADM}}$, recomputes all mass–radius curves, and solves the coupled $(\delta h,\delta\phi_L,\delta m)$ radial stability problem. |

### 2.2. English Translation Files (`step1_*_en.md`–`step20_*_en.md`, `README_en.md`, `MANIFEST_en.txt`)

| Step | English File Name | Size (B) | SHA-256 (prefix) | One-Line Description |
|---:|---|---:|---|---|
| **1 (EN)** | `step1_dilaton_form_en.md` | 10 967 | `2c9a146a645f68da` | Derives the exact equivalence between the Cauchy-derived exponential factor $f_L = \exp(-\gamma_L\bar{\rho}_L) = e^{-\phi_L}$ and the string-theory dilaton matter coupling $A(\phi)=e^{\beta\phi}$. |
| **2 (EN)** | `step2_conformal_disformal_en.md` | 18 037 | `11e46f3c0e249a66` | Derives the Jacobian, inverse metric, comoving-clock equivalence, and explicit branch taxonomy of the conformal and disformal metric frames. |
| **3 (EN)** | `step3_yukawa_helmholtz_eom_en.md` | 12 281 | `39d4c08f70704587` | Variationally derives the curved-spacetime massive Klein–Gordon / Helmholtz–Yukawa field equation, Green function, and flat-space Robin approximation at $r=R_*$. |
| **4 (EN)** | `step4_covariant_conservation_hydrostatic_en.md` | 12 009 | `3f3ac1ebe1b6a984` | Proves exact covariant stress-energy conservation $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$ and the modified hydrostatic Euler–TOV balance $h' = -\nu' + \phi_L'$. |
| **5 (EN)** | `step5_coupling_dictionary_ppn_screening_en.md` | 14 724 | `1c281eb0af68813f` | Establishes the exact $-\frac{c^4}{2G\alpha_Y}$ normalized coupling dictionary ($\alpha_{\text{std}} = 1/(4\pi)$, $\omega_{\text{BD}} \approx 11.07$) and separates conformal vs. disformal PPN behavior. |
| **6 (EN)** | `step6_sparc_galactic_ngc3198_consistency_en.md` | 13 947 | `ea83e05ff331e3ef` | Verifies 0 negative-density points on 175 SPARC galaxies and contrasts NGC 3198's required $C_{\text{req}}=2.720$ against the isolated constant-coupling prediction $C_{\text{pred}}^{\text{isol}} \le 1.0796$ as a negative control. |
| **7 (EN)** | `step7_neutron_star_tov_stability_en.md` | 19 506 | `d0be25723a09fdcb` | Derives the curved-spacetime TOV equations, compares flat-Robin surface mass vs. curved-exterior ADM mass, and separates Branch I-A vs. Branch II-Fluid stability. |
| **8 (EN)** | `step8_microscopic_clock_coupling_siriusb_lab_v20_en.md` | 15 643 | `3332edd07261c63f` | Derives H0 vs. conditional H1 clock coupling from the Dirac–Maxwell action and quantifies Sirius B redshift exclusions, the $8558\;\text{kg}$ granite chamber ($\sim 10^{-24}$), and the 8-coil solenoid ($\sim 10^{-48}$). |
| **9 (EN)** | `step9_vhe_photon_propagation_hess_dispersion_en.md` | 11 477 | `003941bcc65c727f` | Derives the achromatic refractive index $n_t(s)=f_L(s)^{-1}$ in the disformal metric and the multi-source H.E.S.S./MAGIC VHE gamma-ray timing audit. |
| **10 (EN)** | `step10_bilateral_kernel_reciprocity_eft_synthesis_en.md` | 12 913 | `399f6fce0c03dc32` | Derives Formulation I's auxiliary-field bilateral kernel action, reciprocal potential $V_{\text{eff},N}(r)$, and zeroth- vs. second-moment correspondence to Formulation II. |
| **11 (EN)** | `step11_chameleon_screening_multiscale_eft_synthesis_en.md` | 12 531 | `4ad3ed4131afd2be` | Formulates environment-dependent chameleon/symmetron screening (`Branch II-Cham`) as an explicit nonlinear extension hypothesis beyond the constant-$\alpha_Y$ V23 action. |
| **12 (EN)** | `step12_second_variation_coupled_pulsation_hessian_en.md` | 18 198 | `919d54ad8e3455a7` | Explicitly derives $\mathcal{Q}_{\xi\xi}, \mathcal{Q}_{\phi\phi}, \mathcal{C}_{\xi\phi}, \mathcal{C}_{\phi\xi}$ after metric constraint elimination, moving-surface matching, and numerical stability on the curved `yukawa_action_el` sequence. |
| **13 (EN)** | `step13_curved_window_metric_variation_synge_world_function_en.md` | 11 256 | `e2ae56ad1b197b06` | Derives the exact metric variation $\delta[\sqrt{h}K_L]/\delta g^{\mu\nu}$ of the curved-spacetime top-hat window via Synge's world function and the $\delta_D(L-s_{\text{geo}})$ surface stress. |
| **14 (EN)** | `step14_gravitational_lensing_null_geodesics_cosmological_background_en.md` | 15 538 | `632373a149b21224` | Derives conformal vs. disformal light deflection, FLRW cancellation in $\nabla\ln f_L$, and exact horizon scales for both $\Omega_b=0.0493$ ($57.98\;\text{Gpc}$) and $\Omega_m=0.315$ ($22.94\;\text{Gpc}$). |
| **15 (EN)** | `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20_en.md` | 12 993 | `d038c9b3de65558a` | Proves that while $T^\mu_{\;\mu,\text{EM}}=0$ suppresses conformal trace coupling, disformal coupling $u_\mu u_\nu T^{\mu\nu}_{\text{EM}}$ sources $\phi_L$ via $\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$ with exact polarity-even cross-term cancellation. |
| **16 (EN)** | `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress_en.md` | 13 553 | `24622e0f5d74bbe4` | Derives the V21/V23 algebraic X-sector anisotropic stress reconstruction from the Bianchi identity ($\max|p_{tX}/P_c|=0.140$) and compares it with the dynamical Yukawa stress tensor. |
| **17 (EN)** | `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness_en.md` | 13 788 | `9a04ce99cd38bbe2` | Analytically derives the V23 Eq. (8) factorized Chandrasekhar–hydrostatic difference identity and proves the non-positivity of the `gradient`-mode `metric` inertia weight for $L \ge 367\;\text{m}$. |
| **18 (EN)** | `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance_en.md` | 20 532 | `d04a05e7c0507687` | Derives the V23 Eq. (9) nonlocal static tangent BVP on the V21 background, the 3D spherical-shell kernel $K_{ij}$, and the exact inertia-scaling invariance $\omega^2_{4W}/\omega^2_W = 0.2500$. |
| **19 (EN)** | `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy_en.md` | 18 090 | `0b4953fec2eb9fce` | Derives the 3D IN–OUT differential quadrature of the $8558\;\text{kg}$ granite chamber across three geometric $L$-regimes and the Sirius B–Sirius A `common_background` vs. `fixedM` redshift proxy. |
| **20 (EN)** | `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature_en.md` | 16 843 | `ec747ce9b2db70b9` | Derives the 8-coil elliptic Biot–Savart field, peak signal $\Delta y(L_p) = -1.2017\times 10^{-48}$, $r=L/t$ compactified exterior quadrature ($7.44\times 10^{-16}$ agreement), and the $L=0.95\;\text{m}$ switching identity. |
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
| **Curved Exterior Yukawa BVP and ADM Mass (`Branch II-Fluid`)** | $R_{\max}$ convergence ($[R_*, R_*+N_{\text{mult}}L]$, $N_{\text{mult}}\in\{10,12,14,16\}$), $k_{\text{curved}}$ vs. $k_{\text{flat}}$, and $M_{\text{ADM},\max}$ | $R_{\max}$ spread: $\Delta k_{\text{curved}}/k_{\text{curved}} = 6.17\times 10^{-12}$, $\Delta M_{\text{ADM}}/M_{\text{ADM}} = 5.54\times 10^{-13}$; SLy $15.0$, $L=10\;\text{km}$, $\alpha_Y=1$: $k_{\text{curved}} = 0.233542\;\text{km}^{-1}$ vs. $k_{\text{flat}} = 0.187216\;\text{km}^{-1}$ (**$+24.74\%$**); $M_{\text{ADM},\max}^{\text{SLy}} = \mathbf{2.03321\,M_\odot}$ ($M_{\text{surf},\max} = 2.02911\,M_\odot$, $\Delta M_{\phi,\text{ext}}^{\text{curved}} = 0.00411\,M_\odot$), $R_{1.4}^{\text{ADM}} = \mathbf{11.5359\;\text{km}}$ ($M_{\text{ADM},\max}^{\text{APR4}} = 2.18068\,M_\odot$, $M_{\text{ADM},\max}^{\text{H4}} = 1.97809\,M_\odot$) | `step7`, `step22`, `curved_yukawa_adm_and_stability.py` |
| **Coupled $(\delta h,\delta\phi_L,\delta m)$ Stability on the Curved `yukawa_action_el` Branch** | $\omega_0^2 = 0$ zero crossing vs. $\frac{dM_{\text{ADM}}}{d\rho_c} = 0$ and $\frac{dM_B}{d\rho_c} = 0$ | SLy $L=10\;\text{km}$, $\alpha_Y=1$: $\log_{10}\rho_{0c}(M_{\text{ADM},\max}) = 15.316404$, $\log_{10}\rho_{0c}(M_{B,\max}) = 15.316388$, coupled $\omega_0^2=0$: **$15.316425$** ($\mathbf{+2.04\times 10^{-5}\;\text{dex}}$, **$595\times$** improvement over frozen-$\phi_L$ error $+1.22\times 10^{-2}\;\text{dex}$); $\alpha_Y=1/3$: **$+5.84\times 10^{-6}\;\text{dex}$** (**$715\times$** improvement); $L=1\;\text{km}$: **$-6.86\times 10^{-6}\;\text{dex}$** (**$165\times$** improvement); APR4 ($L=10\;\text{km}$): **$+2.07\times 10^{-5}\;\text{dex}$** (**$413\times$** improvement); H4 ($L=10\;\text{km}$): **$+5.71\times 10^{-6}\;\text{dex}$** (**$4137\times$** improvement) | `step12`, `step22`, `curved_yukawa_adm_and_stability.py` |
| **V21 Nonlocal 3D Top-Hat Tangent BVP (`Branch I-A`)** | SLy $L = 10\;\text{km}$ `mass` / `gradient` $\omega_0^2 = 0$ vs. $\partial M/\partial\rho_c = 0$ and $W\to 4W$ scaling | In `mass` mode: $15.29931$ vs. $15.29924$ ($+6.97\times 10^{-5}\;\text{dex}$, $58.3\times$ improvement); in `gradient` mode: $15.30113$ vs. $15.30106$ ($+7.39\times 10^{-5}\;\text{dex}$, $59.3\times$ improvement); inertia scaling: $\omega^2_{4W}/\omega^2_W = 0.2500$ | `step7`, `step17`, `step18`, `covariant_closure_50aa2db0.py` |
| **Sirius B H0 vs. H1 Redshift** | $v_{\text{obs}} = 80.65\pm 0.77\;\text{km/s}$ ($M=1.018\,M_\odot$, $R=5846.6\;\text{km}$) | H0: $|\Delta v|\le 1.27\;\text{m/s}$ (unconstrained); conditional H1 $2\sigma$ ($\Delta\chi^2=4$) excludes: $L\in[4.08\times 10^6, 3.37\times 10^7]\;\text{m}$ ($n=1.5$, $\max\Delta\chi^2=52.3$ / $72.0$) | `step8`, `step19`, `covariant_closure_50aa2db0.py` |
| **Granite Chamber ($8558\;\text{kg}$) IN–OUT Test** | Conditional H1 differential clock signal $\Delta y_{\text{H1}}(\text{IN}-\text{OUT})$ | At $L=1.00\;\text{m}$: Top-hat $-1.436\times 10^{-24}$ ($\equiv 0$ for $L\ge 2.944\;\text{m}$); Yukawa $-3.781\times 10^{-25}$ ($-4.818\times 10^{-25}$ at $L=10\;\text{m}$) | `step8`, `step19`, `covariant_closure_50aa2db0.py` |
| **8-Coil Solenoid ($N=2000$, $14\;\text{mA}$)** | Magnetic peak signal ($L_p=0.630957\;\text{m}$), $r=L/t$ quadrature, and $L=0.95\;\text{m}$ | $\Delta W(L_p)=1.96114452\;\text{J/A}^2$, $k_g = -6.13117767\times 10^{-45}\;\text{A}^{-2}$, $\Delta y(L_p) = -1.2017108\times 10^{-48} \sim -10^{-48}$; Path A vs. B ($L=3\;\text{m}$): $7.44\times 10^{-16}$; $L=0.95\;\text{m}$ identity: $8.88\times 10^{-16}\;\text{J/A}^2$ | `step15`, `step20`, `make_revision3_checks.py` |

---

### 3.1. V27 Re-Review Fixes: 4-Component Exterior BVP, Exact Exterior Linearization, Dynamical Scalar Wave, and GL2 N-Convergence

*(Source: `curved_yukawa_adm_and_stability.py`, `curved_yukawa_audit_summary.json`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_adm_sequences.csv`.)*

- **Exterior Lapse (Issue 5):** the exterior $\nu$ is co-integrated as the fourth ODE component $(\phi_L,\psi_L,m,\nu)$ with `solve_bvp` `tol=1e-12` and outer boundary condition $\nu(R_{\max})=\tfrac12\ln(1-2m(R_{\max})/R_{\max})$. Spreads for SLy, $\log_{10}\rho_{0c}=15.0$, $L=10\;\text{km}$, $N_{\text{mult}}\in\{10,12,14,16\}$: $\delta k/k=6.17\times 10^{-12}$, $\delta M_{\text{ADM}}/M_{\text{ADM}}=5.54\times 10^{-13}$, $\delta\nu_R=1.63\times 10^{-14}$.

| $N_{\text{mult}}$ | $R_{\max}$ (km) | $k_{\text{curved}}$ (km$^{-1}$) | $M_{\text{ADM}}$ ($M_\odot$) | $\nu_R$ |
|---|---|---|---|---|
| `10` | `111.4657` | `0.233542329280` | `1.509998138706` | `-0.246306687228` |
| `12` | `131.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `14` | `151.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `16` | `171.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |

- **Exterior Linearization (Issue 2):** `linearized_exterior_response` obtains $\delta\phi_L'(R_*)$, $\delta m(R_*)$, and $\delta M_{\text{ADM}}=\delta m(\infty)$ from the exact first-order variational BVP in $(\delta\phi_L(R_*),\delta m(R_*),\omega^2)$; the eigenvalue condition is $\delta M_{\text{ADM}}=0$.
- **Dynamical Scalar Block (Issues 1 and 3):** the scalar block includes the $-\omega^2 e^{-2\nu}\delta\phi_L$ term and the exterior response is supplied by the $\omega$-dependent BVP (`dynamical_scalar_wave=True`); solver: GL2 piecewise-segment Volterra discretization, $N=480$. The quasi-static and dynamical roots agree to all 6 reported decimal places.
- **Turning Points:** `minimize_scalar` + 7-point quadratic fit over $\pm0.002$ dex; zero crossings from a 3-point quadratic fit.

| Sequence | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dynamical $\omega_0^2=0$ | $\Delta$ Dynamical (dex) | Quasi-Static $\omega_0^2=0$ | Frozen $\delta\phi_L=0$ | $\Delta$ Frozen (dex) |
|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316388` | **`15.316425`** | **`+2.04e-05`** | `15.316425` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.306886` | **`15.307172`** | **`+5.84e-06`** | `15.307172` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303330` | **`15.303080`** | **`-6.86e-06`** | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.292989` | **`15.293042`** | **`+2.07e-05`** | `15.293042` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227981` | **`15.227790`** | **`+5.71e-06`** | `15.227790` | `15.251385` | `+2.36e-02` |

**GL2 $N$-Convergence** (SLy, $L=10\;\text{km}$, $\alpha_Y=1$, $\log_{10}\rho_{0c}=15.311404$): $p_{\text{eff}}=5.23$ (from the three finest levels), Richardson value $\omega_0^2=6387684.936\;\text{s}^{-2}$; at $N=480$ the relative error is $\approx1.28\times 10^{-5}$. The sequence converges monotonically without sign reversal (the $N$-dependent sign flip of the earlier FD solver does not occur).

| $N$ | $\omega_0^2$ ($\text{s}^{-2}$) | Relative Difference from Richardson Value |
|---|---|---|
| `80` | `6450570.762819` | `9.84e-03` |
| `120` | `6403998.296759` | `2.55e-03` |
| `160` | `6393544.377742` | `9.17e-04` |
| `240` | `6388970.545464` | `2.01e-04` |
| `360` | `6387947.756886` | `4.11e-05` |
| `480` | `6387766.842088` | `1.28e-05` |
| `600` | `6387713.828704` | `4.52e-06` |
| `800` | `6387691.354875` | `1.00e-06` |
| `1200` | `6387685.705914` | `1.21e-07` |

**Interpretive Scope:** the agreement between $\omega_0^2=0$ and the $M_{\text{ADM}}$ peak is an internal consistency check of the same system of equations; the $\sim10^{-5}$ dex offset lies within the resolution of the 3-point root fit and 7-point turning-point fit, and is not an independent physical test.

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
Every command below can be executed directly from the unpacked `/workspace/ta-v27-work/scalar-tensor-yukawa-fit` working directory without network access:

1. **Curved Exterior `yukawa_action_el` Einstein–Scalar BVP, ADM Mass ($M_{\text{ADM}}$), Mass–Radius Curves, and Coupled $(\delta h, \delta\phi_L, \delta m)$ Radial Stability (`curved_yukawa_adm_and_stability.py`):**
   ```bash
   python3 curved_yukawa_adm_and_stability.py
   ```
   - *Expected Output (~2–3 minutes):*  
     - Generates `curved_yukawa_adm_sequences.csv`, `curved_yukawa_coupled_stability.csv`, and `curved_yukawa_audit_summary.json` (and updates `data/v23/`).  
     - Verifies $R_{\max}$-independence ($\Delta k_{\text{curved}}/k_{\text{curved}} = 6.17\times 10^{-12}$, $\Delta M_{\text{ADM}}/M_{\text{ADM}} = 5.54\times 10^{-13}$), curved vs. flat Robin slope ($+24.74\%$ for SLy $\log_{10}\rho_{0c}=15.0$, $L=10\;\text{km}$, $\alpha_Y=1$), ADM mass peaks ($M_{\text{ADM},\max}^{\text{SLy},10\text{km}} = 2.03321\,M_\odot$, $R_{1.4}^{\text{ADM}} = 11.5359\;\text{km}$), and coupled stability zero crossings ($\Delta\log_{10}\rho_{0c} = 2.04\times 10^{-5}\;\text{dex}$, $595\times$ improvement).

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
