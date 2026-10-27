# Step 2: Conformal, Disformal, and Effective Fluid-Coupled Jordan Metrics ($\tilde{g}_{\mu\nu}$) and the Taxonomy of V19–V23 Model Branches

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step2_conformal_disformal_en.md`  
**Previous Step:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` Hungarian manuscript): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Sections 1, 1.1, 1.2, and 3, Eqs. 1–12; DOI: `10.5281/zenodo.23639084`)  
- V23 (`50aa2db0` English manuscript): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex`  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` and `hu_p2.tex` (DOI: `10.5281/zenodo.23604423`)  
- V23 SymPy verification script: `covariant_closure_50aa2db0.py` (lines 105–138)

---

## 1. Verbatim Citations from the V23 Manuscript (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`)

### 1.1. Formulation of the Covariant Action (V23 EN/HU Section 1, Eq. (1), lines 33–35)
Section 1 of the V23 manuscript (*"Covariant Action Principle: Two Explicit Formulations"*) establishes the following covariant action using the $(-+++)$ metric signature convention:
$$S\bigl[g_{\mu\nu},\,\phi_L,\,\Psi_{\text{matter}}\bigr] \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}\,R(g) \;+\; \mathcal{L}_{\text{matter}}\!\bigl(\Psi_{\text{matter}},\,\tilde{g}_{\mu\nu}(g,\phi_L)\bigr) \;+\; \mathcal{L}_{\phi}\!\bigl(g_{\mu\nu},\,\phi_L\bigr)\right]. \tag{V23 Eq. 1}$$

### 1.2. Verbatim Definition of the Microscopic $H_0$ vs. $H_1$ Coupling and the Two Effective Jordan Metrics (V23 Section 3, Eqs. (11)–(12), lines 95–101)
Section 3 of the V23 manuscript (*"Microscopic Clock and Spectral-Line Coupling from the Action ($H_0$ vs. $H_1$)"*) defines the two coupling branches and the two Jordan metrics verbatim as follows:
> *"In the covariant action (1), the microscopic coupling of atomic clocks and electromagnetic spectral lines is uniquely determined by which metric contracts the indices of the Dirac + Maxwell Lagrangian $\mathcal{L}_{\text{atom}}(\psi_e, A_\mu)$:*
> - **Minimal physical-metric coupling ($H_0$ — default covariant closure):** *$\mathcal{L}_{\text{atom}}(\psi_e, A_\mu; g_{\mu\nu})$ couples directly to the physical spacetime metric $g_{\mu\nu}$, while the macroscopic baryonic energy–momentum source is weighted by $f_L$ through Formulation I or II. Then the proper time of a stationary atom is $d\tau_{H_0} = \sqrt{-g_{00}}\,dt = e^{\nu(r)}dt$, the local transition frequencies $\omega_{\text{atom}} = m_e c^2 \alpha_{\text{EM}}^2/(2\hbar)$ are spacetime constants, and atomic clocks measure precisely the metric lapse function arising from the $f_L$-modified Einstein equations. The macroscopic clock anomaly appears exclusively through the gravitational potential change $\Delta\Phi_N$ originating from the modified mass:*
>   $$\left.\frac{\Delta\nu}{\nu}\right|_{H_0} \;=\; \frac{\Delta\Phi_{\text{metric}}}{c^2}, \qquad \Delta z_{H_0}(\text{Sirius B}) \;=\; z_{\text{GR}}\!\left(\frac{\Delta M_B}{M_B}\right)_{\text{TAL}}. \tag{V23 Eq. 11}$$
> - **Disformal / conformal proper-time coupling ($H_1$ — direct clock-coupling hypothesis):** *$\mathcal{L}_{\text{atom}}(\psi_e, A_\mu; \tilde{g}_{\mu\nu})$ couples to the effective metric*
>   $$\tilde{g}_{\mu\nu}^{\text{disf}} \;=\; g_{\mu\nu} \;+\; \frac{1 - f_L(\phi_L)^2}{c^2}\,u_\mu u_\nu \qquad\text{or}\qquad \tilde{g}_{\mu\nu}^{\text{conf}} \;=\; f_L(\phi_L)^2\,g_{\mu\nu} \tag{V23 Eq. 12}$$
>   *(with fixed $m_e$), so that the proper-time element is $d\tilde{\tau}_{H_1} = \sqrt{-\tilde{g}_{00}}\,dt = f_L(\mathbf{x})\,e^{\nu(\mathbf{x})}dt$. This induces a direct non-metric frequency shift:*
>   $$\left.\frac{\Delta\nu}{\nu}\right|_{H_1}(\mathbf{x}) \;=\; \frac{\Delta\Phi_{\text{metric}}(\mathbf{x})}{c^2} \;+\; \bigl(f_L(\mathbf{x}) - 1\bigr) \;\simeq\; \frac{\Delta\Phi_{\text{metric}}(\mathbf{x})}{c^2} \;-\; \phi_L(\mathbf{x}).\text{"}$$

---

## 2. Derivation of the Conformal ($\tilde{g}_{\mu\nu}^{\text{conf}}$) Transformation and the Brans–Dicke Relation

In Step 1 (`step1_dilaton_form.md`), we proved that $f_L(\phi_L) = e^{-\phi_L}$. Substituting this into the second equation of V23 Eq. (12):
$$\boxed{\;\tilde{g}_{\mu\nu}^{\text{conf}}(x) \;=\; f_L(\phi_L(x))^2\,g_{\mu\nu}(x) \;=\; e^{-2\phi_L(x)}\,g_{\mu\nu}(x),\qquad \Omega(\phi_L) \;\equiv\; f_L(\phi_L) \;=\; e^{-\phi_L}.\;}$$

### 2.1. Inverse Metric, Determinant, and Proper-Time Element
1. **Inverse conformal metric:** Since $\tilde{g}_{\mu\alpha}^{\text{conf}}\tilde{g}_{\text{conf}}^{\alpha\nu} = \delta_\mu^\nu$,
   $$\tilde{g}_{\text{conf}}^{\mu\nu} \;=\; f_L(\phi_L)^{-2}\,g^{\mu\nu} \;=\; e^{+2\phi_L}\,g^{\mu\nu}.$$
2. **Four-volume element (determinant):** In $D=4$ spacetime dimensions, $\det(e^{-2\phi_L}g_{\mu\nu}) = (e^{-2\phi_L})^4 \det(g_{\mu\nu})$, therefore:
   $$\sqrt{-\tilde{g}_{\text{conf}}} \;=\; f_L(\phi_L)^4\,\sqrt{-g} \;=\; e^{-4\phi_L}\,\sqrt{-g}.$$
3. **Proper-time element of a static observer ($dx^i = 0$):** In a static spherically symmetric physical metric ($g_{00} = -e^{2\nu(r)}$):
   $$d\tilde{\tau}_{\text{conf}}^2 \;=\; -\frac{1}{c^2}\tilde{g}_{00}^{\text{conf}}\,dt^2 \;=\; -f_L^2\frac{g_{00}}{c^2}\,dt^2 \;=\; f_L(r)^2\,e^{2\nu(r)}\,dt^2 \quad\implies\quad \boxed{\;d\tilde{\tau}_{\text{conf}} \;=\; f_L(r)\,e^{\nu(r)}\,dt \;=\; e^{\nu(r) - \phi_L(r)}\,dt.\;}$$

### 2.2. Conformal Transformation of the Ricci Scalar (Jordan Frame $\leftrightarrow$ Einstein Frame)
From standard $D=4$ conformal geometry (e.g., Wald 1984, Appendix D; Fujii & Maeda 2003), if $g_{\mu\nu} = \Omega^{-2}\tilde{g}_{\mu\nu}^{\text{conf}} = e^{+2\phi_L}\tilde{g}_{\mu\nu}^{\text{conf}}$ (that is, $\tilde{g}_{\mu\nu}^{\text{conf}} = \Omega^2 g_{\mu\nu}$ with $\Omega = e^{-\phi_L}$), the relation between the two Ricci scalars is:
$$R(g) \;=\; \Omega^2\left[\tilde{R}(\tilde{g}_{\text{conf}}) \;+\; 6\,\tilde{\Box}\ln\Omega \;-\; 6\,\tilde{g}_{\text{conf}}^{\mu\nu}\,\partial_\mu\ln\Omega\,\partial_\nu\ln\Omega\right].$$
Accordingly, rewriting the Einstein–Hilbert term in the Jordan frame ($\sqrt{-g} = \Omega^{-4}\sqrt{-\tilde{g}_{\text{conf}}}$) yields:
$$\int d^4x\,\sqrt{-g}\,\frac{c^4}{16\pi G}\,R(g) \;=\; \int d^4x\,\sqrt{-\tilde{g}_{\text{conf}}}\,\frac{c^4}{16\pi G}\left[\Omega^{-2}\tilde{R}(\tilde{g}_{\text{conf}}) \;-\; 6\,\Omega^{-4}\tilde{g}_{\text{conf}}^{\mu\nu}\partial_\mu\Omega\,\partial_\nu\Omega\right],$$
where the boundary divergence term $\nabla_\mu(\sqrt{-g}\,6\Omega\nabla^\mu\Omega)$ has been dropped via integration by parts. Introducing the **Brans–Dicke scalar field**:
$$\boxed{\;\Phi_{\text{BD}}(x) \;\equiv\; \Omega(x)^{-2} \;=\; f_L(x)^{-2} \;=\; e^{+2\phi_L(x)},\;}$$
for which $\frac{\partial_\mu\Phi_{\text{BD}}}{\Phi_{\text{BD}}} = +2\partial_\mu\phi_L$, the local scalar kinetic term of V23 Formulation II ($\mathcal{L}_\phi = -\frac{c^4}{8\pi G\tilde{\alpha}_L}\bigl[\frac{1}{2}g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{1}{2L^2}\phi_L^2\bigr]$, V23 Eq. 8) in the Jordan frame becomes:
$$-\frac{1}{2}g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L \;\longrightarrow\; -\frac{\omega_{\text{BD}}}{\Phi_{\text{BD}}}\,\tilde{g}_{\text{conf}}^{\mu\nu}\partial_\mu\Phi_{\text{BD}}\,\partial_\nu\Phi_{\text{BD}}.$$
Matching the coefficients yields the standard scalar–tensor identity:
$$2\omega_{\text{BD}} + 3 \;=\; \frac{1}{\tilde{\alpha}_L} \;=\; \frac{1}{2\beta^2} \quad\Longleftrightarrow\quad \boxed{\;\omega_{\text{BD}} \;=\; \frac{1}{2\tilde{\alpha}_L} - \frac{3}{2} \;=\; \frac{1}{4\beta^2} - \frac{3}{2}.\;}$$
*(The exact numerical values of $\tilde{\alpha}_L$, $\alpha_Y = 8\pi\tilde{\alpha}_L$, $\alpha_{\text{std}} = 2\tilde{\alpha}_L = \frac{\alpha_Y}{4\pi}$, and $\omega_{\text{BD}}$, together with the Cassini PPN bound, are analyzed in detail in Step 5.)*

---

## 3. The Disformal ($\tilde{g}_{\mu\nu}^{\text{disf}}$) Coupling and Algebraic Projector Decomposition

The first equation of V23 Eq. (12) defines a purely timelike **Bekenstein-type disformal metric** (Bekenstein 1993, *Phys. Rev. D* 48, 3641) using the four-velocity $u^\mu$ of the medium ($g_{\mu\nu}u^\mu u^\nu = -c^2$, $u_\mu = g_{\mu\nu}u^\nu$):
$$\boxed{\;\tilde{g}_{\mu\nu}^{\text{disf}} \;=\; g_{\mu\nu} \;+\; \frac{1 - f_L(\phi_L)^2}{c^2}\,u_\mu u_\nu \;=\; g_{\mu\nu} \;+\; \frac{1 - e^{-2\phi_L}}{c^2}\,u_\mu u_\nu.\;}$$

### 3.1. Exact Projector Decomposition (SymPy-Verified Identity, `covariant_closure_50aa2db0.py` lines 108–125)
Let us introduce the spatial projector of the spacetime metric orthogonal to $u^\mu$:
$$h_{\mu\nu} \;\equiv\; g_{\mu\nu} + \frac{1}{c^2}u_\mu u_\nu \qquad\left(h_{\mu\nu}u^\nu = 0,\quad h^\mu_{\;\alpha}h^\alpha_{\;\nu} = h^\mu_{\;\nu},\quad h^\mu_{\;\mu} = 3\right).$$
Then $g_{\mu\nu} = h_{\mu\nu} - \frac{1}{c^2}u_\mu u_\nu$, so the separated projector form of the disformal metric is:
$$\boxed{\;\tilde{g}_{\mu\nu}^{\text{disf}} \;=\; h_{\mu\nu} \;-\; f_L(\phi_L)^2\,\frac{u_\mu u_\nu}{c^2} \;=\; \left(g_{\mu\nu} + \frac{u_\mu u_\nu}{c^2}\right) \;-\; e^{-2\phi_L}\,\frac{u_\mu u_\nu}{c^2}.\;}$$
From this, the physical interpretation follows immediately:
1. **The disformal coupling leaves the spatial geometry untouched:** for every spatial displacement $\delta x^\mu$ orthogonal to $u^\mu$ ($u_\mu \delta x^\mu = 0$), $\tilde{g}_{\mu\nu}^{\text{disf}}\delta x^\mu\delta x^\nu = g_{\mu\nu}\delta x^\mu\delta x^\nu$.
2. **It scales the timelike (lapse) component by precisely the factor $f_L(\phi_L) = e^{-\phi_L}$:** the proper time of a clock comoving with the medium ($\delta x^\mu = \frac{u^\mu}{c}d\tau$) is:
   $$d\tilde{\tau}_{\text{disf}}^2 \;=\; -\frac{1}{c^2}\tilde{g}_{\mu\nu}^{\text{disf}}\,dx^\mu dx^\nu \;=\; f_L(\phi_L)^2\left(-\frac{1}{c^2}g_{\mu\nu}dx^\mu dx^\nu\right) \;=\; f_L(\phi_L)^2\,d\tau^2.$$
3. **Inverse disformal metric and determinant:** Since $h^{\mu\alpha}u_\alpha = 0$ and $\bigl(-\frac{u_\mu u_\alpha}{c^2}\bigr)\bigl(-\frac{u^\alpha u^\nu}{c^2}\bigr) = -\frac{u_\mu u^\nu}{c^2}$, the inverse metric and volume element are exactly:
   $$\boxed{\;\tilde{g}_{\text{disf}}^{\mu\nu} \;=\; h^{\mu\nu} \;-\; f_L(\phi_L)^{-2}\,\frac{u^\mu u^\nu}{c^2} \;=\; g^{\mu\nu} \;-\; \frac{f_L(\phi_L)^{-2} - 1}{c^2}\,u^\mu u^\nu,\qquad \sqrt{-\tilde{g}_{\text{disf}}} \;=\; f_L(\phi_L)\,\sqrt{-g} \;=\; e^{-\phi_L}\sqrt{-g}.\;}$$

### 3.2. General Static Spherically Symmetric Background and `SymPy` Verification
In a static spherically symmetric physical background metric:
$$ds^2 \;=\; g_{\mu\nu}dx^\mu dx^\nu \;=\; -e^{2\nu(r)}c^2 dt^2 + e^{2\lambda(r)}dr^2 + r^2 d\Omega^2, \qquad g_{00} = -e^{2\nu(r)},$$
the four-velocity of a stationary (comoving) observer is $u^\mu = \bigl(c/\sqrt{-g_{00}},\,0,\,0,\,0\bigr) = \bigl(c\,e^{-\nu(r)},\,0,\,0,\,0\bigr)$, and its covariant form is $u_\mu = g_{\mu\nu}u^\nu = \bigl(-\sqrt{-g_{00}}\,c,\,0,\,0,\,0\bigr) = \bigl(-c\,e^{\nu(r)},\,0,\,0,\,0\bigr)$. Thus:
$$\frac{u_0 u_0}{c^2} \;=\; -g_{00} \;=\; +e^{2\nu(r)},$$
and substituting this yields the $(0,0)$-component of the disformal metric:
$$\tilde{g}_{00}^{\text{disf}} \;=\; g_{00} \;+\; \bigl(1 - f_L(r)^2\bigr)\frac{u_0 u_0}{c^2} \;=\; g_{00} \;+\; \bigl(1 - f_L(r)^2\bigr)\bigl(-g_{00}\bigr) \;=\; f_L(r)^2\,g_{00} \;=\; -f_L(r)^2\,e^{2\nu(r)}.$$
In `covariant_closure_50aa2db0.py` (lines 108–125), within the local orthonormal tetrad ($g_{\hat{a}\hat{b}} = \eta_{\hat{a}\hat{b}} = \text{diag}(-1,1,1,1)$, $u^{\hat{a}} = (c,0,0,0)$, $u_{\hat{a}} = (-c,0,0,0)$), the same identity is verified: $\tilde{g}_{\hat{0}\hat{0}}^{\text{disf}} = -1 + (1 - f_L^2) = -f_L^2$, so the general coordinate expression is $\tilde{g}_{00}^{\text{disf}} = e_{\;0}^{\hat{0}} e_{\;0}^{\hat{0}} \tilde{g}_{\hat{0}\hat{0}}^{\text{disf}} = e^{2\nu}(-f_L^2) = f_L^2 g_{00}$.

---

## 4. Explicit Taxonomy of V19–V23 Model Branches (Branch Taxonomy Table)

For peer-review transparency, the table below unambiguously separates the model branches appearing across the V19–V23 manuscripts and the present derivation chain, including their action principle/definition, volume element, hydrostatic equation, clock/light coupling, and domain of validity:

| Branch Code | Name and Manuscript Source | Scalar / Kernel Definition and 2nd Moment | Matter Coupling / Volume Element | Internal Hydrostatic Equilibrium | Electromagnetic / Clock / PPN Behavior |
|---|---|---|---|---|---|
| **Branch I-A** (`tophat_mass`, `tophat_gradient`) | V21 phenomenological 3D spherical top-hat TOV (`stellar_tal_v21.py`, V23 Eqs. 4–5) | Nonlocal 3D spherical window: $\phi_L = \gamma_L(K_L*\varepsilon)$, $K_L = \frac{3}{4\pi L^3}\Theta(L-r)$, $1 + \frac{L^2}{10}\nabla^2 + \dots$ | Phenomenological $f_L$ multiplier in $m'$ (`mass`) or $\nu'$ (`gradient`), with algebraic $T^\mu_{\;\nu}(X)$ closure | $h'(r) = -\nu'(r)$ (`mass`: $\nu'=A_{\text{geo}}$; `gradient`: $\nu'=f_L A_{\text{geo}}$) | $H_0$: minimal physical metric; $H_1$: direct $f_L$ clock shift (excluded by Sirius B for $L\ge 1\;\text{km}$) |
| **Branch I-B** (`bilateral_eom`, `bilateral_lapse_weighted`) | V23 Formulation I: bilateral nonlocal kernel action (V23 Eqs. 2–3, 6–7) | Symmetric two-point kernel $K_L(x,y) = K_L(y,x)$, $1 + \frac{L^2}{10}\nabla^2 + \dots$ | $S_{\text{I},m} = -\int d^4x\sqrt{-g}\,e^{-\phi_L}\varepsilon$ (in fixed-metric Euclidean kernel approximation) | $V_{\text{eff}}(r) = f_L(r) - \gamma_L\int K_L f_L\varepsilon$ with reciprocity force | Macroscopic hydrostatic test branch; metric variation inside the kernel is not closed |
| **Branch II-Fluid** (`yukawa_action_el`) | V23 Formulation II: local Yukawa-screened Einstein–scalar–fluid action (V23 Eqs. 8–10) | Local Helmholtz–Yukawa PDE: $(-\Box_g + L^{-2})\phi_L = \frac{\alpha_Y G}{c^2}e^{-\phi_L}\varepsilon$, $1 + L^2\nabla^2 + \dots$ | $S_m = -\int d^4x\sqrt{-g}\,e^{-\phi_L}\varepsilon(n)$, $P_{\text{matter}} = e^{-\phi_L}P$, $\rho_\phi = \frac{B(\phi_L')^2 + \phi_L^2/L^2}{2\alpha_Y}$ | $h'(r) = -\nu'(r) + \phi_L'(r)$ (fifth force + curved exterior scalar tail $M_{\text{ADM}}$) | Under $H_0$, atomic clocks couple to physical $g_{\mu\nu}$; turning points of $M_{\text{ADM}}$ and $M_B$ coincide |
| **Branch II-Conf** ($\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$) | Universal conformal Jordan–Brans–Dicke / dilaton metric (V23 Eq. 12b, `step1`, `step5`) | Same local Yukawa PDE; $\Phi_{\text{BD}} = e^{+2\phi_L}$ | $\sqrt{-\tilde{g}_{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$; conformally invariant Maxwell action ($\tilde{F}_{\mu\nu}\tilde{F}^{\mu\nu}\sqrt{-\tilde{g}} = F_{\mu\nu}F^{\mu\nu}\sqrt{-g}$) | In the non-relativistic limit ($\varepsilon \simeq \rho_0 c^2 \gg 3P$), coincides with the $\beta = \sqrt{\alpha_{\text{std}}/2}$ fifth force | Standard conformal Brans–Dicke / PPN dictionary applies to this branch: $\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$, $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2} \approx 11.07$ ($\alpha_Y=1$), light deflection is conformally invariant |
| **Branch II-Disf** ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$) | Comoving timelike disformal clock and EM metric (V23 Eq. 12a, V20/V22 solenoid, `step14`–`step15`) | Same local Yukawa PDE, augmented with magnetic source $\rho_B = \frac{|\mathbf{B}|^2}{2\mu_0 c^2}$ | $\sqrt{-\tilde{g}_{\text{disf}}} = f_L\sqrt{-g} = e^{-\phi_L}\sqrt{-g}$ (coincides with the $e^{-\phi_L}$ weight of Branch II-Fluid!) | If $u^\mu$ is fixed by a comoving khronon/aether field, energy weight is $e^{-\phi_L}$; free $u^\mu$ variation yields an extra constraint | Refractive index $n_{\text{refr}} = e^{+\phi_L}$ (refractive light deflection), magnetic solenoid source ($\rho_B$); **not** described by the conformal Brans–Dicke PPN formula |
| **Branch II-Cham** ($\beta_{\text{eff}}(\rho)$ extension) | Density-dependent chameleon / symmetron nonlinear EFT hypothesis (`step11`) | Density-dependent effective mass $m_{\text{eff}}^2(\rho)$ and coupling $\beta_{\text{eff}}(\rho)$ | Nonlinear $V(\phi) + A(\phi)\rho$ extension beyond the quadratic constant-$\alpha_Y$ V23 action | Environment-dependent fifth force (strong screening in laboratory/NS, enhancement in dilute galactic medium) | **Standalone extension hypothesis**, not a consequence of the constant-$\alpha_Y$ V23 `yukawa_action_el` equations |

---

## 5. Checks

1. **Limit Check ($\phi_L \to 0$, i.e., $f_L \to 1$):**
   - When $\phi_L = 0$ ($f_L = 1$), $\tilde{g}_{\mu\nu}^{\text{conf}} = g_{\mu\nu}$ and $\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + 0 = g_{\mu\nu}$, so both Jordan metrics reduce to the physical Einstein metric $g_{\mu\nu}$, and $d\tilde{\tau}_{H_1} = d\tau_{H_0}$ $\checkmark$.
2. **Proper-Time Agreement for Stationary Clocks (along $u^\mu$):**
   - For a clock comoving with the baryonic medium ($dx^\mu \parallel u^\mu$), the conformal and disformal metrics yield the **identical** proper-time element and frequency shift:
     $$d\tilde{\tau}_{\text{conf}} \;=\; d\tilde{\tau}_{\text{disf}} \;=\; f_L(x)\,d\tau \;=\; e^{-\phi_L(x)}\,\sqrt{-g_{00}}\,dt,$$
     which is also verified symbolically by `covariant_closure_50aa2db0.py` (`sym_res["proper_time_ratio_H1_over_H0"] = "f_L"`, `sym_res["weak_field_H1_clock_shift"] = "Phi_N/c^2 - phi_L"`) $\checkmark$.
3. **Difference in the Spatial Sector, Volume Element, and Light Propagation:**
   - Whereas $\sqrt{-\tilde{g}_{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$ and the conformal metric leaves null geodesics (light rays) invariant ($\tilde{g}_{\mu\nu}^{\text{conf}}k^\mu k^\nu = 0 \iff g_{\mu\nu}k^\mu k^\nu = 0$), the disformal metric satisfies $\sqrt{-\tilde{g}_{\text{disf}}} = e^{-\phi_L}\sqrt{-g}$ and induces an effective refractive index ($n = f_L^{-1} = e^{+\phi_L}$) for light while leaving spatial distances unscaled $\checkmark$.
