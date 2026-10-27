# Step 5: Dictionary of Coupling Constants, Weak-Field Yukawa Fifth Force, and Conformal vs. Disformal PPN Light Deflection

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step5_coupling_dictionary_ppn_screening_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`), `step2_conformal_disformal_en.md` (`step2_conformal_disformal.md`), `step3_yukawa_helmholtz_eom_en.md` (`step3_yukawa_helmholtz_eom.md`), `step4_covariant_conservation_hydrostatic_en.md` (`step4_covariant_conservation_hydrostatic.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` Hungarian manuscript): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 1.2, Eqs. 6–7, lines 56–65; Section 3, lines 161–167; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (Chapter 3, lines 245–255, Eqs. 4–5) and `hu_p2.tex` (Chapter 4, lines 23–38; DOI: `10.5281/zenodo.23604423`)

---

## 1. Verbatim Citations from the V22 and V23 Manuscripts

In Step 2 (`step2_conformal_disformal.md`), we established the taxonomy of the V19–V23 model branches and showed that one must distinguish between the universal conformal Jordan metric (**Branch II-Conf**), the comoving timelike disformal metric (**Branch II-Disf**), and the effective fluid-density coupling (**Branch II-Fluid**).

1. **The Coupling Constant $C$, $\gamma_L$, and the Compactness Scale in the V22 Monograph (`hu_p1.tex`, lines 247–254, Eqs. 4–5):**
   $$\begin{equation}
     C \;\equiv\; \frac{4\pi G}{c^2} \;=\; 9.33306\times 10^{-27}\;\text{m}\,\text{kg}^{-1}, \qquad
     \gamma_L \;=\; \frac{G\,L^2}{c^2} \;=\; \frac{C\,L^2}{4\pi} \;=\; 7.425915\times 10^{-28}\,L^2\;\text{m}^3\,\text{kg}^{-1}. \tag{V22 Eqs. 4–5}
   \end{equation}$$
   > *"When a point mass $M$ is averaged over a sphere of radius $L$, its mean density is $\bar{\rho}_L = \frac{3M}{4\pi L^3}$, so the exponent is $\gamma_L \bar{\rho}_L = \frac{3GM}{4\pi c^2 L} = \frac{3\,C\,M}{(4\pi)^2 L}$, which up to the geometric factor $3/(4\pi)$ recovers the dimensionless compactness scale $GM/(c^2 L)$."*

2. **Solar System Screening in Chapter 4 of V22 (`hu_p2.tex`, lines 32 and 38):**
   > *"For subgalactic $L \gtrsim 1\,\text{pc}$, $\gamma_L\bar{\rho}_L^\odot = \frac{3GM_\odot}{4\pi c^2 L} \lesssim 10^{-14}$; planetary orbits remain untouched! [...] when $L$ exceeds the size of planetary orbits, the solar mass $M_\odot$ is diluted over the volume $\frac{4}{3}\pi L^3$, so gradients of $f_L$ along planetary orbits are suppressed below the Cassini and ephemeris bounds."*

3. **The V23 Yukawa Normalization Constants (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 65):**
   > *"The unit-normalized kernel corresponds to $\alpha_Y = 1$, i.e., $\alpha_{\text{std}} = 1/(4\pi)\approx 0.0796$, which recovers the homogeneous top-hat limit ($\phi_L = \gamma_L\rho_{\text{src}}$); $\alpha_Y = 1/3$ ($\alpha_{\text{std}} = 1/(12\pi)$) is the volume-matched variant."*

---

## 2. Weak-Field Effective Potential, Yukawa Fifth Force, and the Effective Gravitational Constant ($G_{\text{eff}}(r)$)

In Step 4 (`step4_covariant_conservation_hydrostatic.md`), we derived from the covariant conservation law ($\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$) the modified hydrostatic / Euler equation:
$$h'(r) \;=\; -\nu'(r) \;+\; \phi_L'(r), \qquad \mathbf{a}_{\text{eff}}(\mathbf{r}) \;=\; -\nabla\Phi_{\text{eff}}(\mathbf{r}) \;=\; -\nabla\Bigl(\Phi_N(\mathbf{r}) \;-\; c^2\,\phi_L(\mathbf{r})\Bigr).$$

### 2.1. Exterior Yukawa Potential Around a Compact Source ($r > R_*$)
In the exterior of a spherically symmetric body of mass $M$ and radius $R_*$ ($r > R_*$, in the weak-field approximation), the flat-space Green's function solution written in Step 3 (`step3_yukawa_helmholtz_eom.md`, V23 Eq. 7) takes the form:
$$\phi_L^{\text{Yuk}}(r) \;=\; \alpha_Y\,\frac{G}{c^2}\,\frac{L\,e^{-r/L}}{r}\int_0^{R_*} r'\,\sinh\!\left(\frac{r'}{L}\right)\rho_{\text{src}}(r')\,dr' \;=\; \frac{\alpha_{\text{std}}\,G M_{\text{eff}}(L)}{c^2\,r}\,e^{-r/L},$$
where $\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$, and in the point-source limit ($R_* \ll L$, where $\sinh(r'/L) \simeq r'/L$):
$$M_{\text{eff}}(L) \;=\; 4\pi\int_0^{R_*} r'^2\left[\frac{L}{r'}\sinh\!\left(\frac{r'}{L}\right)\right]\rho_{\text{src}}(r')\,dr' \;\xrightarrow{R_*\ll L}\; \int_0^{R_*} 4\pi r'^2\rho_{\text{src}}(r')\,dr' \;=\; M.$$
Substituting this alongside the Newtonian potential $\Phi_N(r) = -\frac{GM}{r}$ into the effective potential gives:
$$\Phi_{\text{eff}}(r) \;=\; \Phi_N(r) \;-\; c^2\,\phi_L^{\text{Yuk}}(r) \;=\; -\frac{G M}{r}\Bigl(1 \;+\; \alpha_{\text{std}}\,e^{-r/L}\Bigr).$$
Differentiating with respect to $r$ yields the distance-dependent effective gravitational constant:
$$\boxed{\;a_{\text{eff}}(r) \;=\; -\frac{G_{\text{eff}}(r)\,M}{r^2}, \qquad G_{\text{eff}}(r) \;=\; G\left[1 \;+\; \alpha_{\text{std}}\left(1 + \frac{r}{L}\right)e^{-r/L}\right].\;}$$

### 2.2. Canonical Scalar–Tensor Normalization and the Relation $\alpha_{\text{std}} = 2\beta^2$
Let us rewrite the V23 Eq. (8) scalar action:
$$S_\phi \;=\; -\int d^4x\,\sqrt{-g}\,\frac{c^4}{16\pi G}\frac{1}{\tilde{\alpha}_L}\left(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{\phi_L^2}{L^2}\right) \;=\; -\int d^4x\,\sqrt{-g}\,\frac{c^4}{2 G\,\alpha_Y}\left(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{\phi_L^2}{L^2}\right)$$
(where $\alpha_Y = 8\pi\tilde{\alpha}_L = 4\pi\alpha_{\text{std}}$, so $\tilde{\alpha}_L = \frac{\alpha_{\text{std}}}{2} = \frac{\alpha_Y}{8\pi}$, and in units $G=c=1$ the scalar Lagrangian is precisely $-\frac{1}{2\alpha_Y}[(\nabla\phi_L)^2 + \phi_L^2/L^2]$) in terms of the canonical Einstein-frame scalar field $\varphi$ of Damour–Esposito-Farèse (1992) / Damour–Polyakov (1994), whose kinetic term is $-\frac{c^4}{4\pi G}\,\frac{1}{2}(\nabla\varphi)^2 = -\frac{c^4}{8\pi G}(\nabla\varphi)^2$:
$$\frac{c^4}{2 G\,\alpha_Y}(\nabla\phi_L)^2 \;=\; \frac{c^4}{8\pi G\,\alpha_{\text{std}}}(\nabla\phi_L)^2 \;\equiv\; \frac{c^4}{8\pi G}\,2(\nabla\varphi)^2 \quad\implies\quad \phi_L \;=\; \sqrt{2\alpha_{\text{std}}}\,\varphi \;=\; 2\beta\,\varphi.$$
Since the matter-coupling factor is $f_L = e^{-\phi_L} = e^{-2\beta\varphi}$, the exact relation between the dimensionless scalar–matter coupling constant ($\beta$) and the V23 parameters $\alpha_{\text{std}}, \alpha_Y, \tilde{\alpha}_L$ is:
$$\boxed{\;\alpha_{\text{std}} \;=\; 2\beta^2 \;=\; 2\tilde{\alpha}_L \;=\; \frac{\alpha_Y}{4\pi}, \qquad \beta \;=\; \sqrt{\frac{\alpha_{\text{std}}}{2}} \;=\; \sqrt{\frac{\alpha_Y}{8\pi}}.\;}$$
- For the **unit-normalized top-hat counterpart** ($\alpha_Y = 1$):
  $$\alpha_{\text{std}} \;=\; \frac{1}{4\pi} \;\approx\; 0.079577, \qquad \beta \;=\; \frac{1}{\sqrt{8\pi}} \;\approx\; 0.199471, \qquad 2\omega_{\text{BD}} + 3 \;=\; \frac{1}{\beta^2} \;=\; 8\pi \implies \omega_{\text{BD}} \;=\; 4\pi - \frac{3}{2} \;\approx\; 11.066.$$
- For the **volume-matched variant** ($\alpha_Y = 1/3$):
  $$\alpha_{\text{std}} \;=\; \frac{1}{12\pi} \;\approx\; 0.026526, \qquad \beta \;=\; \frac{1}{\sqrt{24\pi}} \;\approx\; 0.115165, \qquad \omega_{\text{BD}} \;=\; 12\pi - \frac{3}{2} \;\approx\; 36.199.$$

---

## 3. Strict Separation of the Three Coupling Branches (Branch II-Conf, Branch II-Fluid, Branch II-Disf) in the PPN and Light-Deflection Sector

In accordance with peer-review requirements, we record that the conformal coupling $\tilde{g}_{\mu\nu}^{\text{conf}}$, the disformal coupling $\tilde{g}_{\mu\nu}^{\text{disf}}$, and the effective fluid coupling $S_m = -\int d^4x\sqrt{-g}\,e^{-\phi_L}\varepsilon(n)$ represent **three distinct physical hypotheses** regarding the electromagnetic/null-geodesic sector, and consequently their PPN behaviors also differ:

1. **Branch II-Conf — Universal Conformal Jordan Coupling ($\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$, $\sqrt{-\tilde{g}_{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$):**
   - In this branch, every matter and electromagnetic field couples universally to the metric $\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$.
   - Because the Maxwell action is conformally invariant in $D=4$ dimensions ($\sqrt{-\tilde{g}_{\text{conf}}}\,\tilde{g}_{\text{conf}}^{\mu\alpha}\tilde{g}_{\text{conf}}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} = \sqrt{-g}\,g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}$), photon null geodesics ($\tilde{g}_{\mu\nu}^{\text{conf}}dx^\mu dx^\nu = 0 \iff g_{\mu\nu}dx^\mu dx^\nu = 0$) sense only the Einstein-frame metric $g_{\mu\nu}$, and the electromagnetic field does not source the scalar $\phi_L$ ($T^\mu_{\;\mu}(\text{EM}) = 0$).
   - The standard scalar–tensor / Brans–Dicke PPN parameter (Will, *Living Rev. Relativ.* **17**, 4, 2014) holds **exclusively for this Branch II-Conf**:
     $$\boxed{\;\gamma_{\text{PPN}}^{\text{conf}}(r) \;=\; \frac{1 - \alpha_{\text{std}}\,e^{-r/L}}{1 + \alpha_{\text{std}}\,e^{-r/L}},\qquad \gamma_{\text{PPN}}^{\text{conf}}(r \ll L) \;=\; \frac{2\omega_{\text{BD}} + 1}{2\omega_{\text{BD}} + 3} \;=\; \frac{4\pi - 1}{4\pi + 1} \;\approx\; 0.8526\;(\alpha_Y=1),\qquad \gamma_{\text{PPN}}^{\text{conf}}(r \gg L) \;\to\; 1.\;}$$

2. **Branch II-Fluid + $H_0$ — Effective Baryonic Fluid Coupling with Minimal Photon/Atomic Metric ($S_m = -\int d^4x\sqrt{-g}\,e^{-\phi_L}\varepsilon(n)$, $\Psi_{\text{atomic/EM}}$ Coupled to $g_{\mu\nu}$):**
   - Macroscopic baryonic fluid experiences the fifth force ($h' = -\nu' + \phi_L'$), so the effective gravitational constant between two macroscopic bodies is $G_{\text{eff}}(r) = G\bigl[1 + \alpha_{\text{std}}(1 + r/L)e^{-r/L}\bigr]$.
   - However, in the regime $r \gg L$ (for example when $L$ is submillimeter or when on astronomical scales $r_{\text{AU}} \gg L$), the Yukawa term decays exponentially ($e^{-r/L} \to 0$), the exterior spacetime is asymptotically Schwarzschild geometry (with mass $M_{\text{ADM}}$), and $\gamma_{\text{PPN}}(r \gg L) = 1$, $\beta_{\text{PPN}}(r \gg L) = 1$.

3. **Branch II-Disf + $H_1$ — Comoving Timelike Disformal Coupling ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$, $\sqrt{-\tilde{g}_{\text{disf}}} = f_L\sqrt{-g}$):**
   - In this branch, the volume element $\sqrt{-\tilde{g}_{\text{disf}}} = e^{-\phi_L}\sqrt{-g}$ coincides exactly with the $e^{-\phi_L}$ weight factor of Branch II-Fluid, but the four-velocity $u^\mu$ enters the metric explicitly (which either requires a comoving khronon/aether fixing or yields an additional constraint equation under free $u^\mu$ variation).
   - Because the disformal transformation is non-conformal, it breaks conformal invariance in the Maxwell sector: on the one hand, a static magnetic field sources $\phi_L$ with effective mass density $\rho_B = \frac{|\mathbf{B}|^2}{2\mu_0 c^2}$ (Steps 15 and 20), and on the other hand, light rays propagate on the disformal light cone $\tilde{g}_{\mu\nu}^{\text{disf}}k^\mu k^\nu = 0$ with an effective refractive index $n_{\text{refr}}(r) = f_L(r)^{-1} = e^{+\phi_L(r)}$ (Step 14).
   - Consequently, light deflection in Branch II-Disf is **not** described by the conformal Brans–Dicke formula $\gamma_{\text{PPN}}^{\text{conf}} = (1 - \alpha_{\text{std}}e^{-r/L})/(1 + \alpha_{\text{std}}e^{-r/L})$, but rather by the refractive deflection angle $\hat{\alpha}_{\text{disf}}(b) = \frac{4GM}{c^2 b} - \int_{-\infty}^\infty \frac{\partial_\perp\phi_L}{c^2}dz$ derived in Step 14.

---

## 4. Summary Coupling Dictionary (TAL $\leftrightarrow$ Scalar–Tensor / Dilaton)

| Physical Quantity | Time Attraction Law (V19–V23) | Scalar–Tensor / Dilaton Counterpart | Exact Conversion and Branch-Specific Note |
| :--- | :--- | :--- | :--- |
| **Dimensionless scalar field** | $\phi_L(\mathbf{x})$ | $2\beta\varphi(\mathbf{x})$ (or $2\phi_{\text{dil}}(\mathbf{x})$) | $f_L(\mathbf{x}) = e^{-\phi_L} \leftrightarrow e^{-2\beta\varphi} \leftrightarrow e^{-2\phi_{\text{dil}}}$ |
| **Fundamental gravitational–temporal constant** | $C \equiv \frac{4\pi G}{c^2} = 9.33306\times 10^{-27}\;\text{m/kg}$ | $8\pi G / c^4$ ($\times c^2/2$) | $\frac{G}{c^2} = \frac{C}{4\pi} = 7.425915\times 10^{-28}\;\text{m/kg}$ |
| **Scale-dependent density coupling** | $\gamma_L = \frac{G L^2}{c^2} = \frac{C L^2}{4\pi}$ ($\text{m}^3/\text{kg}$) | $\frac{G}{c^2\mu_\phi^2}$ | $[\gamma_L] = \text{m}^3\,\text{kg}^{-1} = [\text{density}]^{-1}$ |
| **Smoothing radius / Range** | $L$ (3D top-hat / Yukawa scale) | Compton wavelength $\lambda_\phi = \mu_\phi^{-1}$ | $m_\phi = \frac{\hbar}{L c}$ (scalar/dilaton mass) |
| **Action normalization constants** | $\alpha_Y = 4\pi\alpha_{\text{std}} = 8\pi\tilde{\alpha}_L$ | $8\pi\beta^2 = \frac{4\pi}{2\omega_{\text{BD}} + 3}$ | Scalar prefactor in action: $-\frac{c^4}{16\pi G}\frac{1}{\tilde{\alpha}_L} = -\frac{c^4}{2 G\alpha_Y}$ ($-\frac{1}{2\alpha_Y}$ in $G=c=1$) |
| **Conformal branch (Branch II-Conf)** | $\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$ | Jordan–Brans–Dicke metric | $\sqrt{-\tilde{g}_{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$; $\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$, $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2}$ |
| **Disformal branch (Branch II-Disf)** | $\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$ | Bekenstein disformal metric | $\tilde{g}_{00} = f_L^2 g_{00},\; \tilde{g}_{ij} = g_{ij},\; \sqrt{-\tilde{g}_{\text{disf}}} = f_L\sqrt{-g}$; magnetic source $\rho_B$ and refractive index $n=e^{\phi_L}$ |

---

## 5. Checks

1. **Yukawa Limits ($r \ll L$ and $r \gg L$):**
   - $r \gg L$ (strongly screened regime, $e^{-r/L} \to 0$): $G_{\text{eff}}(r) \to G$ and $\gamma_{\text{PPN}}^{\text{conf}}(r) \to 1$, so the model recovers general relativity at all distances $r \gg L$ $\checkmark$.
   - $r \ll L$ (unscreened short distance, $e^{-r/L} \to 1$): $G_{\text{eff}}(r) \to G(1 + \alpha_{\text{std}}) = G(1 + 2\beta^2)$, which matches the standard scalar–tensor gravity two-body interaction formula (Damour & Esposito-Farèse 1992, Eq. 2.12) $\checkmark$.
2. **Solar System Dilution Bound of the 3D Top-Hat (V22 `hu_p2.tex`, line 32):**
   For the solar mass ($M_\odot = 1.98847\times 10^{30}\;\text{kg}$) at $L = 1\;\text{pc} = 3.0857\times 10^{16}\;\text{m}$:
   $$\gamma_L\bar{\rho}_L^\odot \;=\; \frac{3 G M_\odot}{4\pi c^2 L} \;=\; \frac{3 \times 1.4766\times 10^3\;\text{m}}{4\pi \times 3.0857\times 10^{16}\;\text{m}} \;\approx\; 1.14\times 10^{-14} \;\ll\; 10^{-5} \quad (\text{Cassini PPN bound}) \quad \checkmark$$
