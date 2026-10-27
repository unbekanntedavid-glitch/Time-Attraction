# Step 7: Compact Objects, Curved-Spacetime Scalar-Tensor / Yukawa TOV Equations (Flat-Robin Approximation vs. Curved Exterior ADM Mass), and Static Nonlocal Tangent Stability Consistency (V21 / V23)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step7_neutron_star_tov_stability_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step6_sparc_galactic_ngc3198_consistency_en.md` (`step6_sparc_galactic_ngc3198_consistency.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p2.tex` (Chapter 7, lines 102–124, Eq. 10) and `hu_p3.tex` (Chapter 8, lines 1–28, Eqs. 11–12; DOI: `10.5281/zenodo.23604423`)  
- V21 Zenodo record: DOI `10.5281/zenodo.23602501` (Neutron Star Conditional TOV Sensitivity and Radial Pulsation Study)  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Sections 1.1–1.2 and 2, lines 37–159, Eqs. 4–9, Tables 1–2; DOI: `10.5281/zenodo.23639084`)  
- Curved exterior `yukawa_action_el` ADM mass and coupled stability: `curved_yukawa_adm_and_stability.py` and `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`

---

## 1. Verbatim Citations from the V22 Monograph and the V23 (`50aa2db0`) Manuscript

In this step, we derive the equations of **compact objects and neutron-star TOV stability (V21 / V23)**, strictly separating:
1. the nonlocal static tangent operator V23 Eq. (9) formulated on the **Branch I-A** (V21 nonlocal 3D top-hat `tophat_mass` and `tophat_gradient`) background, and
2. the **Branch II-Fluid** (V23 local `yukawa_action_el`) coupled Einstein–scalar–fluid system, for which we separately compare the **flat-space Robin approximation and surface mass** ($M_{\text{surf}} = m(R_*)$) used in Table 1 of the V23 manuscript with the **exact Dirichlet–Neumann slope and asymptotic ADM mass determined from the coupled curved exterior vacuum solution ($r \in [R_*, \infty)$)** ($M_{\text{ADM}} = \lim_{r\to\infty} m(r) = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$, detailed derivation in Step 22, `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`).

1. **The Two Conditional TOV Modes in V21, the $1.214\times 10^{-6}$ Kernel Accuracy, the 769 Stellar Models (Including the 33 Causality Endpoints), and the Methodological Stipulation on $e^\nu F$ (`hu_p2.tex`, Chapter 7, lines 107–124):**
   > *"In our audited code, the measured maximum relative discrepancy between the analytical shell kernel and independent 3D quadrature across all profiles is **$1.214\times 10^{-6}$**. I then integrated two conditional Tolman–Oppenheimer–Volkoff (TOV) closures:"*
   $$\begin{equation}
     \textbf{(A) \texttt{mass} mode:}\quad \frac{dm_{\text{eff}}}{dr} = 4\pi r^2 \rho(r)\,f_L(r), \qquad
     \textbf{(B) \texttt{gradient} mode:}\quad \frac{d\Phi_{\text{eff}}}{dr} = f_L(r)\,\frac{G\bigl(m(r) + 4\pi r^3 P(r)/c^2\bigr)}{c^2 r\bigl(r - 2Gm(r)/c^2\bigr)}, \tag{V22 Eq. 10}
   \end{equation}$$
   > *"I integrated 33 equilibrium sequences across three nuclear equations of state — **SLy** (nucleonic), **APR4** (variational nucleonic; reaching the causal sound-speed limit $c_s = c$ at rest-mass density **$\rho_0 \approx 1.4458\times 10^{15}\,\text{g/cm}^3$**), and **H4** (stiff hyperonic) — over 11 logarithmically spaced smoothing scales $L \in [1\,\text{m}, 100\,\text{km}]$. The V21 archival baseline table contains **182 equilibrium models**, while the dense diagnostic stability grid consists of **769 stellar models in total** (including all **33 exact causality endpoints** where $c_s/c = 1$)."*
   > *"When evaluating the TOV equations, I explicitly state: **the appearance of the modified hydrostatic factor $e^\nu F$ (or $\Phi_{\text{eff}}'$) in the conditional pressure-gradient equation is an algebraic equilibrium assumption, and by itself does not derive how atomic clocks, nuclear transitions, or surface photon redshifts couple to the factor $f_L$.**"*

2. **Algebraic $X$-Sector Stress Reconstruction and the Limitations of Frozen-$f_L$ Chandrasekhar Operators (`hu_p3.tex`, Chapter 8, lines 8–22 and `tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, lines 48–52):**
   > *"Across the 769 audited neutron-star models, the relative tangential stress remains bounded, with a maximum of **$\max |p_{tX}/P_c| = 0.140$** ($14.0\%$ of the central pressure). I emphasize that this $X$-sector reconstruction is an algebraic bookkeeping closure of the background Bianchi identity, not yet a dynamical field equation derived from an action principle."*
   $$\begin{align}
     \textbf{Mass mode } (m' = 4\pi r^2 f_L\varepsilon,\; \nu' = A_{\text{geo}}): &\quad \varepsilon_X = -(1-f_L)\varepsilon,\quad P_{r,X} = 0,\quad P_{t,X} = -\frac{r}{2}(1-f_L)\varepsilon\nu', \tag{V23 Eq. 4}\\
     \textbf{Gradient mode } (m' = 4\pi r^2\varepsilon,\; \nu' = f_L A_{\text{geo}}): &\quad \varepsilon_X = 0,\quad P_{r,X} = -(1-f_L)\!\left(P + \frac{m}{4\pi r^3}\right),\; P_{t,X} = P_{r,X} + \frac{r}{2}\bigl(P'_{r,X} + P_{r,X}\nu'\bigr). \tag{V23 Eq. 5}
   \end{align}$$

3. **Curved-Spacetime Yukawa Action-EL (`yukawa_action_el`) Equations of V23 (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 65):**
   > *"The `yukawa_action_el` rows follow from the action $\mathcal{L}=-[(\nabla\phi)^2+\phi^2/L^2]/(2\alpha_Y)-f_L\varepsilon$ ($G=c=1$; `SymPy`: `sym_yukawa_el_source_factor` and `_hydrostatic` residual $=0$): $B\phi''+(2B/r+B\nu'+B'/2)\phi'-\phi/L^2=-\alpha_Y f_L\varepsilon$, $B=1-2m/r$, $P_{\text{matter}}=f_LP$, $h'=-\nu'+\phi'$, $\rho_\phi=(B\phi'^2+\phi^2/L^2)/(2\alpha_Y)$, exterior matching $\phi'(R)=-(1/R+1/L)\phi(R)$; the exterior tail $\Delta M_{\phi,\text{ext}}$ is listed separately (SLy, $L=10\text{ km}$: $3.4\times 10^{-3}\,M_\odot$ at $\log_{10}\rho_{0c}=15$, $9.0\times 10^{-3}\,M_\odot$ at $M_{\max}$)."*

4. **The Frozen-$f_L$ Sturm–Liouville Residual and the Nonlocal Static Tangent Operator (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, Section 2, lines 122–135, Eqs. 8–9):**
   $$\begin{equation}
     \mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} = \frac{4\pi r^{5/2}(f_L - 1)\Gamma P \varepsilon e^{2\nu}}{(r - 2m)^{5/2}}\Big[(4\pi r^3 P - 3r + 7m)\xi - r(r-2m)\xi'\Big]. \tag{V23 Eq. 8}
   \end{equation}$$
   $$\begin{equation}
     \delta f_L(r_i) = -\gamma_L f_L(r_i)\sum_j K_{ij}\,\frac{\varepsilon(r_j)+P(r_j)}{c_s^2(r_j)}\,\delta h(r_j), \qquad
     \begin{cases}
       \dfrac{d(\delta h)}{dr} = -\delta\nu'\!\left[\delta h, \delta m, \delta f_L[\delta h]\right] + \omega_0^2\,e^{\lambda-2\nu}\xi, \\[4pt]
       \dfrac{d(\delta m)}{dr} = 4\pi r^2\!\left[f_{\text{mass}}\frac{\varepsilon+P}{c_s^2}\delta h + \varepsilon\,\delta f_{\text{mass}}[\delta h]\right].
     \end{cases} \tag{V23 Eq. 9}
   \end{equation}$$

---

## 2. Derivation of the Curved-Spacetime Scalar-Tensor / Yukawa TOV System (`yukawa_action_el`)

### 2.1. The Curved-Spacetime d'Alembertian Operator and the Yukawa Field Equation
In a static, spherically symmetric spacetime (in geometrized units $G = c = 1$):
$$ds^2 \;=\; -e^{2\nu(r)}dt^2 + B(r)^{-1}dr^2 + r^2(d\theta^2 + \sin^2\theta\,d\varphi^2), \qquad B(r) \;\equiv\; e^{-2\lambda(r)} \;=\; 1 - \frac{2m(r)}{r},$$
the metric determinant is $\sqrt{-g} = e^{\nu(r)} B(r)^{-1/2} r^2 \sin\theta$. For a scalar field $\phi_L(r)$ depending only on $r$, the covariant Laplace–Beltrami / d'Alembert operator is:
$$\Box_g \phi_L \;=\; \frac{1}{\sqrt{-g}}\partial_r\!\left(\sqrt{-g}\,g^{rr}\,\partial_r\phi_L\right) \;=\; \frac{B^{1/2}}{e^\nu r^2}\frac{d}{dr}\!\left(e^\nu B^{1/2} r^2 \phi_L'\right) \;=\; B\phi_L'' + \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L'.$$
For the Lagrangian $\mathcal{L}_{\phi+m} = -\frac{1}{2\alpha_Y}\left[g^{\mu\nu}\nabla_\mu\phi_L\nabla_\nu\phi_L + \frac{\phi_L^2}{L^2}\right] - f_L(\phi_L)\varepsilon(n)$ ($f_L = e^{-\phi_L}$, $\alpha_Y = 4\pi\alpha_{\text{std}} = 8\pi\tilde{\alpha}_L = 2\beta^2\cdot 4\pi$), the Euler–Lagrange equation with respect to $\phi_L$ ($\frac{1}{\sqrt{-g}}\partial_\mu(\sqrt{-g}\frac{\partial\mathcal{L}}{\partial(\partial_\mu\phi_L)}) - \frac{\partial\mathcal{L}}{\partial\phi_L} = 0$, taking into account $\frac{\partial(-e^{-\phi_L}\varepsilon)}{\partial\phi_L} = +f_L\varepsilon$) is:
$$\boxed{\;B\phi_L'' \;+\; \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L' \;-\; \frac{\phi_L}{L^2} \;=\; -\alpha_Y\,f_L(\phi_L)\,\varepsilon(r).\;}$$
This is precisely the `yukawa_action_el` equation appearing in Section 1.2 of V23 (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 65), where the source term is multiplied by the dilaton-like self-consistent weight factor $f_L = e^{-\phi_L}$ (`sym_yukawa_el_source_factor = 0`).

### 2.2. The Coupled Einstein–Yukawa–Euler (TOV) System, the Curved Exterior Boundary Condition, and the ADM Mass ($M_{\text{ADM}}$ vs. $M_{\text{surf}}$)
From the tensors $T^\mu_{\;\nu}(\text{matter}) = \text{diag}(-f_L\varepsilon, f_L P, f_L P, f_L P)$ and $T^\mu_{\;\nu}(\phi)$ derived in Step 4 (`step4_covariant_conservation_hydrostatic.md`), the $(tt)$ and $(rr)$ Einstein equations together with the conservation law $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$ yield the coupled, action-closed TOV system:
$$\begin{aligned}
  \rho_\phi(r) &\;=\; -T^t_{\;t}(\phi) \;=\; \frac{1}{2\alpha_Y}\!\left(B(r)\bigl[\phi_L'(r)\bigr]^2 + \frac{\phi_L(r)^2}{L^2}\right), \qquad
  p_{r,\phi}(r) \;=\; T^r_{\;r}(\phi) \;=\; \frac{1}{2\alpha_Y}\!\left(B(r)\bigl[\phi_L'(r)\bigr]^2 - \frac{\phi_L(r)^2}{L^2}\right), \\[4pt]
  \frac{dm(r)}{dr} &\;=\; 4\pi r^2\bigl[f_L(r)\varepsilon(r) + \rho_\phi(r)\bigr], \\[4pt]
  \frac{d\nu(r)}{dr} &\;=\; \frac{m(r) + 4\pi r^3\bigl[f_L(r)P(r) + p_{r,\phi}(r)\bigr]}{r\bigl(r - 2m(r)\bigr)}, \\[4pt]
  \frac{dh(r)}{dr} &\;=\; \frac{P'(r)}{\varepsilon(r) + P(r)} \;=\; -\frac{d\nu(r)}{dr} \;+\; \frac{d\phi_L(r)}{dr}.
\end{aligned}$$

**Critical Methodological Clarification (Point 1 and Point 2 — detailed derivation in Step 22):**
1. **Difference Between the Flat-Space Robin Approximation and the True Curved Exterior Dirichlet–Neumann Slope:**
   - The V23 manuscript (the original `yukawa_action_el` row of `covariant_closure_50aa2db0.py`) used as boundary condition at the stellar surface $r = R_*$ the derivative of the flat-space Helmholtz solution ($\phi_L^{\text{flat}}(r) = \phi_L(R_*)\frac{R_*}{r}e^{-(r-R_*)/L}$):
     $$\phi_{L,\text{flat}}'(R_*) \;=\; -k_{\text{flat}}\,\phi_L(R_*), \qquad k_{\text{flat}} \;=\; \frac{1}{R_*} + \frac{1}{L}.$$
   - In reality, however, in the exterior vacuum outside the stellar surface ($r \ge R_*$, where $\varepsilon = P = 0$), spacetime is strongly curved ($B(R_*) = 1 - 2M_{\text{surf}}/R_* \approx 0.50\text{--}0.62$). In the exterior vacuum, the exact identity is $B\nu' + B'/2 = \frac{2m(r)}{r^2} - \frac{4\pi r}{\alpha_Y L^2}\phi_L(r)^2$, so the exterior coupled Einstein–scalar system is:
     $$\left(1 - \frac{2m(r)}{r}\right)\phi_L''(r) + \left[\frac{2}{r}\left(1 - \frac{m(r)}{r}\right) - \frac{4\pi r}{\alpha_Y L^2}\phi_L(r)^2\right]\phi_L'(r) - \frac{\phi_L(r)}{L^2} \;=\; 0, \qquad m'(r) \;=\; 4\pi r^2\rho_\phi(r).$$
     Because $B(R_*) < 1$, the effective local decay rate of the mass term $\phi_L/L^2$ at the surface is $\sim B(R_*)^{-1/2}/L > 1/L$. As proven in Step 22 (`curved_yukawa_adm_and_stability.py`), for SLy at $\log_{10}\rho_{0c} = 15.0$, $L = 10\;\text{km}$, $\alpha_Y = 1$, the true curved Dirichlet–Neumann slope is $k_{\text{curved}} = -\phi_{L,\text{curved}}'(R_*)/\phi_L(R_*) = 0.233542\;\text{km}^{-1}$, which is **$+24.74\%$ larger** than the flat Robin approximation $k_{\text{flat}} = 0.187216\;\text{km}^{-1}$!
2. **Surface Mass ($M_{\text{surf}}$) vs. Asymptotic ADM Mass ($M_{\text{ADM}}$):**
   - In Table 1 of the V23 manuscript, the $M_{\max}$ and $R_{1.4}$ columns of the `yukawa_action_el` row were optimized for the surface mass $m(R_*) = M_{\text{surf}}$, while the exterior scalar mass $\Delta M_{\phi,\text{ext}}$ was listed in a separate column.
   - For a distant observer (and for the Bardeen–Carter–Hawking / Harrison–Wheeler turning-point theorem), the physical gravitational mass is the asymptotic **ADM mass**:
     $$\boxed{\;M_{\text{ADM}} \;\equiv\; \lim_{r \to \infty} m(r) \;=\; M_{\text{surf}} \;+\; \Delta M_{\phi,\text{ext}}^{\text{curved}}, \qquad \Delta M_{\phi,\text{ext}}^{\text{curved}} \;=\; \int_{R_*}^\infty 4\pi r^2 \rho_\phi(r)\,dr.\;}$$
   - On the `yukawa_action_el_curved_adm` sequences recomputed with the curved exterior solution (SLy, $L = 10\;\text{km}$, $\alpha_Y = 1$):
     - $M_{\text{ADM},\max} = \mathbf{2.03321\,M_\odot}$ ($\log_{10}\rho_{0c}(M_{\text{ADM},\max}) = \mathbf{15.316404}$, $R(M_{\text{ADM},\max}) = 9.8853\;\text{km}$, $\Delta M_{\phi,\text{ext}}^{\text{curved}} = 0.00411\,M_\odot$, $M_{\text{surf}} = 2.02911\,M_\odot$),
     - the maximum of the conserved baryonic rest mass ($M_B = \int_0^{R_*} 4\pi r^2 B(r)^{-1/2}\rho_0(r)\,dr$) occurs at $\log_{10}\rho_{0c}(M_{B,\max}) = \mathbf{15.316388}$ ($|\Delta\log_{10}\rho_{0c}| = 1.56\times 10^{-5}\;\text{dex}$ agreement with the ADM-mass turning point!),
     - $R(M_{\text{ADM}} = 1.4\,M_\odot) = \mathbf{11.5359\;\text{km}}$ (whereas $R(M_{\text{surf}} = 1.4\,M_\odot) = 11.5348\;\text{km}$),
     - and the zero crossing of the coupled $(\delta h, \delta\phi_L)$ dynamical radial pulsation eigenvalue solved on the `yukawa_action_el` background is $\log_{10}\rho_{0c}(\omega_0^2 = 0) = \mathbf{15.316425}$ (a discrepancy of $\mathbf{2.04\times 10^{-5}\;\text{dex}}$ from the turning point $dM_{\text{ADM}}/d\rho_c = 0$; see Steps 12 and 22).

---

## 3. Stability Analysis of the V21 Nonlocal 3D Top-Hat Background (**Branch I-A**): Frozen-$f_L$ Residual vs. $\delta f_L[\delta h]$ Tangent Operator

*(Note on separating model branches — Point 4: Subsections 3.1–3.2 below apply specifically to the V21 nonlocal 3D spherical top-hat `tophat_mass` and `tophat_gradient` [**Branch I-A**] backgrounds, whereas the coupled stability of `yukawa_action_el` [**Branch II-Fluid**] is presented in Steps 12 and 22.)*

### 3.1. Algebraic Residual Term of the Frozen-$f_L$ Chandrasekhar Operator (V23 Eq. 8)
In general relativity ($f_L \equiv 1$), taking the derivative $\partial/\partial\rho_c$ along the static equilibrium sequence generates a zero-eigenfrequency ($\omega_0^2 = 0$) radial trial perturbation at the extremal point $dM/d\rho_c = 0$. However, the derivation of the classical second-order Chandrasekhar Sturm–Liouville operator ($\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi]$) substitutes during integration by parts the GR background identity
$$\lambda' + \nu' \;=\; 4\pi r e^{2\lambda} (\varepsilon + P) \qquad (\text{only when } f_L \equiv 1\text{!}).$$
When the background is modified according to the V21 `mass` mode, the difference between the Chandrasekhar operator and the linearized hydrostatic operator evaluated at $\delta f_L = 0$ ($\mathcal{H}_{\text{lin}}[\xi]\big|_{\delta f_L=0}$) does not vanish, but produces precisely the residual proportional to $(f_L - 1)$ given in V23 Eq. (8):
$$\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} \;=\; \frac{4\pi r^{5/2}(f_L - 1)\Gamma P \varepsilon e^{2\nu}}{(r - 2m)^{5/2}}\Big[(4\pi r^3 P - 3r + 7m)\xi - r(r-2m)\xi'\Big].$$

### 3.2. Incorporating the Nonlocal $\delta f_L[\delta h]$ Tangent Perturbation on the V21 Background (V23 Eq. 9 and Table 2)
When the central density is varied by $\rho_{0c} \to \rho_{0c} + \delta\rho_{0c}$ on the V21 3D top-hat background, the Eulerian perturbation of the energy density is $\delta\varepsilon(r) = \frac{\varepsilon(r)+P(r)}{c_s^2(r)}\,\delta h(r)$, which also perturbs the factor $f_L$ through the nonlocal kernel:
$$\delta f_L(r_i) \;=\; -\gamma_L\,f_L(r_i)\,\sum_j K_{ij}\,\frac{\varepsilon(r_j)+P(r_j)}{c_s^2(r_j)}\,\delta h(r_j).$$
Substituting this operator $\delta f_L[\delta h]$ back into the linearized $(\delta h, \delta m)$ boundary-value problem around the V21 background (V23 Eq. 9), at $\omega_0^2 = 0$ the system coincides exactly with the tangent operator of the V21 static equilibrium equations with respect to $\log_{10}\rho_{0c}$. As a result, on the SLy ($L = 10\text{ km}$) V21 sequence, the discrepancy between the zero crossing $\omega_0^2 = 0$ and the static turning point $dM/d\rho_c = 0$ (V23 Table 2) is:
- **GR baseline ($L = 0\text{ km}$):** $\Delta\log_{10}\rho_{0c} = +7.12\times 10^{-5}\;\text{dex}$ (discretization baseline),
- **V21 SLy `mass` mode ($L = 10\text{ km}$, turning point $\log_{10}\rho_{0c} = 15.29924$):** improves from the frozen-$f_L$ BVP discrepancy of $-4.06\times 10^{-3}\;\text{dex}$ (and `matter` SL $-1.34\times 10^{-3}\;\text{dex}$) to **$\mathbf{+6.97\times 10^{-5}\;\text{dex}}$** (a **$58.3\times$** and **$19.2\times$** improvement, respectively),
- **V21 SLy `gradient` mode ($L = 10\text{ km}$, turning point $\log_{10}\rho_{0c} = 15.30106$):** improves from the frozen-$f_L$ BVP discrepancy of $-4.39\times 10^{-3}\;\text{dex}$ (and `matter` SL $-3.57\times 10^{-3}\;\text{dex}$) to **$\mathbf{+7.39\times 10^{-5}\;\text{dex}}$** (a **$59.3\times$** and **$48.3\times$** improvement, respectively).

---

## 4. Checks

1. **$f_L \to 1$ ($L \to 0$ or $\alpha_Y \to 0$) GR Limit:**
   - $\phi_L \to 0 \implies \rho_\phi \to 0,\; p_{r,\phi} \to 0,\; f_L \to 1$, the `yukawa_action_el` TOV system recovers the standard GR TOV equations ($M_{\max}^{\text{SLy}} = 2.04841\,M_\odot$, $M_{\max}^{\text{APR4}} = 2.18800\,M_\odot$, $M_{\max}^{\text{H4}} = 2.01385\,M_\odot$), and the V23 Eq. (8) residual term vanishes identically because $(f_L - 1) = 0$ $\checkmark$.
2. **Agreement Between Curved Exterior `yukawa_action_el` ADM Mass and Coupled Stability Turning Point (see `step22`):**
   - SLy $L = 10\;\text{km}$, $\alpha_Y = 1$: $\log_{10}\rho_{0c}(dM_{\text{ADM}}/d\rho_c = 0) = 15.316404$, $\log_{10}\rho_{0c}(dM_B/d\rho_c = 0) = 15.316388$, coupled $(\delta h, \delta\phi_L)$ zero crossing $\log_{10}\rho_{0c}(\omega_0^2 = 0) = 15.316425$ ($|\Delta\log_{10}\rho_{0c}| = 3.07\times 10^{-5}\;\text{dex}$, a $439\times$ improvement compared to the frozen-$\phi_L$ shift of $+1.22\times 10^{-2}\;\text{dex}$) $\checkmark$.
3. **Inertia-Scaling Invariance at the V21 Zero Crossing ($\omega_0^2 = 0$):**
   - Scaling $W \to 4W$ at $\log_{10}\rho_{0c} = 15.30427$ yields the exact ratio $\omega^2_{4W}/\omega^2_W = 0.2500$ ($\omega^2_W = 4.897\times 10^6\;\text{s}^{-2}$), whereas at the zero crossing $\omega_0^2 = 0$, $F(0)$ is invariant, confirming that the $+6.97\times 10^{-5}\;\text{dex}$ and $+7.39\times 10^{-5}\;\text{dex}$ agreement in V21 is a consequence of the exact consistency of the static nonlocal tangent operator $\checkmark$.

> **V27 update (curved `yukawa_action_el` branch):** 4-component exterior BVP ($\delta\nu_R=1.63\times 10^{-14}$), exact exterior linearization, dynamical scalar wave, GL2 $N$-convergence ($p_{\text{eff}}=5.23$). SLy $L=10$ km, $\alpha_Y=1$: $M_{\text{ADM}}$ turning point $15.316404$, dynamical $\omega_0^2=0$ root $15.316425$ ($2.04\times 10^{-5}$ dex). Detailed tables: `step22`, Section 4.
