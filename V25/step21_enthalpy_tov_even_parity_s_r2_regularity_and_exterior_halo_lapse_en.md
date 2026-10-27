# Step 21: Exact Derivation of the 7-Segment Piecewise-Polytrope Log-Enthalpy ($h$) TOV Transformation, the $s = r^2$ Even-Parity Central Regularization (`Audit Fix O1`), and the `gradient`-Mode Exterior Vacuum-Halo Lapse Integral ($\Delta\nu_{\text{ext}}$) over $r \in [R_*, R_* + L]$

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-27  
**File Path:** `scalar-tensor-yukawa-fit/step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse_en.md`  
**Hungarian Original:** `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse.md`  
**Previous Steps:** `step1_dilaton_form_en.md` – `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature_en.md`  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p2.tex` (Chapter 7, lines 102–120; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 1.3 and Table 1, lines 109–156; DOI: `10.5281/zenodo.23639084`)  
- V21/V23 profile and TOV integrator: `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_profiles.py` (lines 14–108) and `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/stellar_tal_v21.py` (lines 33–136; DOI: `10.5281/zenodo.23602501`)

---

## 1. Verbatim Citations from the V22 and V23 Manuscripts and the `tal_profiles.py` Module

In Steps 7, 16, 17, and 18 we used the V21/V23 neutron-star equilibrium profiles (`dense_profile`, SLy, APR4, and H4 equations of state, continuum Bianchi residuals $1.82\times 10^{-11}$ and $1.59\times 10^{-11}$), without yet deriving the three analytical transformations in `stellar_tal_v21.py` and `tal_profiles.py` that guarantee this accuracy and central regularity:

1. **The 3 Equations of State (SLy, APR4, H4), the Causal $c_s = c$ Bound, and the Two Conditional TOV Modes (`hu_p2.tex`, Chapter 7, lines 108–115):**
   $$\begin{equation}
     \textbf{(A) \texttt{mass} mode:}\quad \frac{dm_{\text{eff}}}{dr} = 4\pi r^2 \rho(r)\,f_L(r), \qquad
     \textbf{(B) \texttt{gradient} mode:}\quad \frac{d\Phi_{\text{eff}}}{dr} = f_L(r)\,\frac{G\bigl(m(r) + 4\pi r^3 P(r)/c^2\bigr)}{c^2 r\bigl(r - 2Gm(r)/c^2\bigr)}, \tag{V22 Ch. 7}
   \end{equation}$$
   > *"I integrated 33 equilibrium sequences across three nuclear equations of state — **SLy** (nucleonic), **APR4** (variational nucleonic; reaches the causal sound-speed threshold $c_s = c$ at rest-mass density **$\rho_0 \approx 1.4458\times 10^{15}\,\text{g/cm}^3$**), and **H4** (stiff hyperonic) — across 11 logarithmically spaced smoothing scales $L \in [1\,\text{m}, 100\,\text{km}]$. The V21 archival baseline table contains **182 equilibrium models**, whereas the dense diagnostic stability grid comprises **769 stellar models** in total (including all **33 exact causal endpoints** where $c_s/c = 1$)."*

2. **The $s = r^2$ Even-Parity Central Mapping and the `gradient`-Mode Exterior Lapse Integral (`tal_profiles.py`, lines 53–56 and 91–93):**
   - Lines 53–56:
     > *"`# Profile mapping, regular at the centre (audit fix O1, run 2): h and mu = m/r^3 are even, analytic functions of r, so they are interpolated in s = r^2 with their exact centre values h(0)=hc, mu(0)=(4pi/3) eps_c f_m. Then h'(r)=2r dh/ds -> 0, m=r^3 mu, A=r(mu+4pi P)/(1-2 mu r^2) are regular. (Run 1 interpolated h(r), m(r) in r with a single node at r=0 -> wrong slopes near the centre.)`"*
   - Lines 91–93:
     > *"`def exterior_lapse(cur, L, M, R, nq=160): \"\"\"nu(R) = -int_R^inf f(r) M/(r(r-2M)) dr with f=exp(-L^2 rhobar_L) from the converged profile; f=1 beyond R+L. Returns nu(R) and diagnostic dict.\"\"\"`"*

3. **Continuum Bianchi Conservation Residuals in V23 Table 1 (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, lines 144–147):**
   - `tophat_mass`: $1.82\times 10^{-11}$, `tophat_gradient`: $1.59\times 10^{-11}$ (exact symbolic residual `0.0`).

---

## 2. Read et al. (2009) 7-Segment Piecewise-Polytrope Log-Enthalpy ($h$) Closure and the $h$-Transformation of the TOV System

### 2.1. Thermodynamic Consistency, $a_i$ Continuity Recursion, and the $c_s^2 / c^2 = 1$ Causal Threshold (`stellar_tal_v21.py`, lines 33–69)
The Read et al. (2009) parameterization consists of 4 crust segments ($\Gamma_0,\dots,\Gamma_3$) and 3 core segments ($\Gamma_4,\Gamma_5,\Gamma_6 = g_1,g_2,g_3$, with dividing densities $\rho_{0,4} = \rho_{\text{match}}$, $\rho_{0,5} = 10^{14.7}\;\text{g/cm}^3$, $\rho_{0,6} = 10^{15}\;\text{g/cm}^3$). In the $i$-th segment ($\rho_0 \in [\rho_{0,i-1}, \rho_{0,i}]$), the pressure (in $c^2$-scaled CGS units) is $P(\rho_0) = K_i \rho_0^{\Gamma_i}$. Pressure continuity at the segment boundaries enforces:
$$K_{i+1} \;=\; K_i\,\rho_{0,i}^{\Gamma_i - \Gamma_{i+1}}, \qquad \rho_{\text{match}} \;=\; \left(\frac{K_3}{K_4}\right)^{\frac{1}{\Gamma_4 - \Gamma_3}}.$$
Integrating the first law of thermodynamics at constant entropy per baryon ($d(\varepsilon / \rho_0) = -P\,d(1/\rho_0) = \frac{P}{\rho_0^2}d\rho_0$) yields the total mass-energy density (`stellar_tal_v21.py`, lines 47–49 and 67):
$$\varepsilon(\rho_0) \;=\; (1 + a_i)\,\rho_0 \;+\; \frac{P(\rho_0)}{\Gamma_i - 1}, \qquad a_{i+1} \;=\; a_i \;+\; \frac{K_i}{\Gamma_i - 1}\rho_{0,i}^{\Gamma_i - 1} \;-\; \frac{K_{i+1}}{\Gamma_{i+1} - 1}\rho_{0,i}^{\Gamma_{i+1} - 1} \quad (a_0 = 0).$$
The specific relativistic log-enthalpy ($dh \equiv \frac{dP}{\varepsilon + P}$) can then be integrated and inverted analytically (`stellar_tal_v21.py`, lines 50 and 60–68):
$$\boxed{\;h(\rho_0) \;=\; \ln\!\left(\frac{\varepsilon + P}{\rho_0}\right) \;=\; \ln\!\left(1 + a_i + \frac{\Gamma_i}{\Gamma_i - 1}K_i\,\rho_0^{\Gamma_i - 1}\right) \;\Longleftrightarrow\; \rho_0(h) \;=\; \left[\frac{\bigl(e^h - 1 - a_i\bigr)(\Gamma_i - 1)}{\Gamma_i K_i}\right]^{\frac{1}{\Gamma_i - 1}}.\;}$$
The squared adiabatic sound speed in each segment is:
$$\frac{c_s^2(h)}{c^2} \;=\; \left(\frac{\partial P}{\partial\varepsilon}\right)_s \;=\; \frac{\Gamma_i\,P(h)}{\varepsilon(h) + P(h)} \;=\; (\Gamma_i - 1)\,\frac{e^h - 1 - a_i}{e^h}.$$
Because $\Gamma_i > 2$ in the high-density core (for APR4, $\Gamma_5 = 3.445$, $\Gamma_6 = 3.348$), $c_s^2(h)/c^2$ increases strictly monotonically with $h$, and the causal threshold ($c_s^2/c^2 = 1$) is located by `brentq` to $10^{-15}$ accuracy at $h_{\text{causal}}$ (for APR4, $\rho_0(h_{\text{causal}}) \approx 1.4458\times 10^{15}\;\text{g/cm}^3$, `hu_p2.tex`, line 115).

### 2.2. Rewriting the TOV Equations with $h$ as Independent Variable and Central Taylor Startup (`tal_profiles.py`, lines 38–50)
Because the hydrostatic balance equation is $P'(r) = (\varepsilon + P)h'(r) = -(\varepsilon + P)f_{\text{grad}}(r)A_{\text{GR}}(r)$, the $(\varepsilon + P)$ factor cancels identically:
$$\frac{dh}{dr} \;=\; -f_{\text{grad}}(r)\,A_{\text{GR}}(r) \;=\; -f_{\text{grad}}(r)\,\frac{m(r) + 4\pi r^3 P(h)}{r\bigl(r - 2m(r)\bigr)},$$
where $f_{\text{grad}}(r) = f_L(r)$ in `gradient` mode and $1$ in `mass` (as well as GR) mode, while $f_{\text{mass}}(r) = f_L(r)$ in `mass` mode and $1$ in `gradient` mode. Since $h(r)$ decreases strictly monotonically inside the star from central $h_c$ down to surface $h(R_*) = 0$, the ODE system for $(r(h), m(h))$ reads (`tal_profiles.py`, lines 43–48):
$$\boxed{\;\frac{dr}{dh} \;=\; -\frac{r(r - 2m)}{\bigl(m + 4\pi r^3 P(h)\bigr)\,f_{\text{grad}}(r)}, \qquad \frac{dm}{dh} \;=\; 4\pi r^2 \varepsilon(h)\,f_{\text{mass}}(r)\,\frac{dr}{dh}.\;}$$
**Why is integrating in $h$ advantageous?** When $r$ is the independent variable, the stellar surface $P(R_*) = 0$ must be found by event detection where $\varepsilon(P) \to 0$ varies steeply. With $h$ as the independent variable, the surface is a fixed integration endpoint ($h_{\text{surf}} = 10^{-10}$) at which $R_* = r(h_{\text{surf}})$ and $M = m(h_{\text{surf}})$ are read out directly without root-finding overshoot! To avoid the $r = 0$ coordinate singularity, the integration starts at $r_0 = 10^{-4}\;\text{km}$ using the second-order central Taylor expansion (`tal_profiles.py`, lines 41–42):
$$m(r_0) \;=\; \frac{4\pi}{3}\varepsilon_c f_m\,r_0^3 + \mathcal{O}(r_0^5), \qquad h(r_0) \;=\; h_c - \delta h, \qquad \delta h \;=\; \frac{2\pi}{3}(\varepsilon_c f_m + 3 P_c)\,f_g\,r_0^2 + \mathcal{O}(r_0^4).$$

---

## 3. Even-Parity $s = r^2$ Central Regularization (`Audit Fix O1`, `tal_profiles.py`, lines 53–79)

### 3.1. Why Direct Spline Interpolation in $r$ Causes Central Singularities, and How $s = r^2$ Resolves It
Any spherically symmetric scalar field that is smooth at the origin ($C^\infty(\mathbb{R}^3)$) — such as the log-enthalpy $h(\mathbf{x})$, the mean interior mass density $\mu(\mathbf{x}) \equiv \frac{m(r)}{r^3}$, or the 3D top-hat exponent $x(\mathbf{x}) \equiv \gamma_L\bar{\rho}_L(r)$ — can contain only powers of $|\mathbf{x}|^2 = x^2 + y^2 + z^2 = r^2 \equiv s$ in its Taylor expansion, i.e., it is an **even function of $r$**:
$$h(r) \;=\; h_s(r^2), \qquad \mu(r) \;\equiv\; \frac{m(r)}{r^3} \;=\; \mu_s(r^2), \qquad x(r) \;=\; x_s(r^2),$$
with $h_s(0) = h_c$ and $\mu_s(0) = \frac{4\pi}{3}\varepsilon_c f_m$. If $h(r)$ and $m(r)$ are interpolated directly in $r$ via PCHIP splines (`Run 1`), the cubic Hermite polynomial at the single node $r = 0$ enforces neither the exact value of $h''(0)$ nor $m'(0) = m''(0) = 0$, injecting numerical noise and $\mathcal{O}(1/r)$ errors into $\frac{m(r)}{r^3}$ and $\frac{h'(r)}{r}$ near the center!

By contrast, `tal_profiles.py` (lines 64–79, `Audit Fix O1`) constructs `PchipInterpolator` objects on the grid $s = r^2 \in [0, R_*^2]$:
$$h_s(s), \qquad \mu_s(s) \;=\; \frac{m(\sqrt{s})}{s^{3/2}} \quad (\mu_s(0) = \mu_c), \qquad x_s(s) \;=\; x(\sqrt{s}).$$
By the chain rule ($\frac{d}{dr} = 2r\frac{d}{ds}$), every central derivative and curvature factor is **analytically regular** as $r \to 0^+$ (`tal_profiles.py`, lines 66–79):
$$\boxed{\;\begin{aligned}
  h'(r) &\;=\; 2r\,\frac{dh_s}{ds}\!(r^2) \;=\; \mathcal{O}(r), \qquad &m(r) &\;=\; r^3\,\mu_s(r^2) \;=\; \mathcal{O}(r^3), \\[4pt]
  \frac{d\mu}{dr}(r) &\;=\; 2r\,\frac{d\mu_s}{ds}\!(r^2) \;=\; \mathcal{O}(r), \qquad &\frac{d\ln f_L}{dr}(r) &\;=\; -2r\,\frac{dx_s}{ds}\!(r^2) \;=\; \mathcal{O}(r), \\[4pt]
  A_{\text{GR}}(r) &\;=\; \frac{m(r) + 4\pi r^3 P(r)}{r\bigl(r - 2m(r)\bigr)} \;=\; \frac{r\bigl(\mu_s(r^2) + 4\pi P(r)\bigr)}{1 - 2r^2\,\mu_s(r^2)} \;=\; \mathcal{O}(r).
\end{aligned}\;}$$
Consequently, in the ratios $\frac{h'(r)}{r} = 2\frac{dh_s}{ds}(r^2)$ and $\frac{A_{\text{GR}}(r)}{r} = \frac{\mu_s(r^2) + 4\pi P(r)}{1 - 2r^2\mu_s(r^2)}$ entering the Chandrasekhar operator and the Bianchi identity, the factor $r$ **cancels algebraically** at $r = 0$! This regularization reduces the continuum Bianchi conservation residual in V23 Table 1 to $1.82\times 10^{-11}$ (`tophat_mass`) and $1.59\times 10^{-11}$ (`tophat_gradient`).

---

## 4. `gradient`-Mode Exterior Vacuum-Halo Lapse Integral ($\Delta\nu_{\text{ext}}$) over $r \in [R_*, R_* + L]$ and the Interior Bernoulli Integral (`tal_profiles.py`, lines 80–106)

### 4.1. Why the $f_L(r) < 1$ Domain Extends Beyond the Stellar Surface $R_*$ out to $R_* + L$
Outside the star ($r > R_*$), the local matter density vanishes identically ($\varepsilon(r) = 0$, $P(r) = 0$), so $m(r) \equiv M$ and $A_{\text{GR}}(r) = \frac{M}{r(r - 2M)}$ for all $r \ge R_*$.
- **In `mass` mode (and in GR):** the exterior lapse gradient $\nu'(r) = A_{\text{GR}}(r) = \frac{M}{r(r-2M)}$ carries no $f_L(r)$ factor. Thus the exterior spacetime is exact Schwarzschild geometry, and the surface lapse is (`tal_profiles.py`, line 81):
  $$\nu_{\text{mass}}(R_*) \;=\; \nu_{\text{GR}}(R_*) \;=\; -\int_{R_*}^\infty \frac{M}{r(r - 2M)}\,dr \;=\; \frac{1}{2}\ln\!\left(1 - \frac{2M}{R_*}\right).$$
- **In `gradient` mode:** the lapse gradient is $\nu'(r) = f_L(r)\,A_{\text{GR}}(r)$ everywhere. Although $\varepsilon(r) = 0$ outside the star, a 3D top-hat ball $B(\mathbf{x}, L)$ of radius $L$ centered at an exterior point $r \in (R_*, R_* + L)$ **still intersects the outer layer $[r - L, R_*]$ of the star** (via the spherical-cap intersection formula derived in Step 18)! Consequently, a **nonlocal vacuum-halo shell** of thickness $L$ ($r \in [R_*, R_* + L]$) forms around the star, within which:
  $$x(r) \;=\; L^2\bar{\rho}_L(r) \;=\; \frac{3}{4r\,L}\int_{r-L}^{R_*} s\,\varepsilon(s)\,\left[L^2 - (r - s)^2\right]\,ds \;>\; 0 \quad\Longrightarrow\quad f_L(r) \;=\; e^{-x(r)} \;<\; 1 \qquad (R_* \le r < R_* + L),$$
  and only for $r \ge R_* + L$ does $x(r) \equiv 0 \implies f_L(r) \equiv 1$!

### 4.2. Derivation of the Finite-Range Integral $\Delta\nu_{\text{ext}}$ and the Interior $\nu(r)$ Profile (`tal_profiles.py`, lines 86–105)
Integrating the `gradient`-mode equation $\nu'(r) = f_L(r)\,A_{\text{GR}}(r)$ from asymptotically flat infinity ($\nu(\infty) = 0$) inward to the stellar surface $r = R_*$, and using $f_L(r) = 1 - (1 - e^{-x(r)})$ where $1 - e^{-x(r)} \equiv 0$ for all $r \ge R_* + L$, splits the improper integral into the analytical Schwarzschild value plus a finite correction over the compact interval $[R_*, R_* + L]$ (`tal_profiles.py`, lines 91–105):
$$\boxed{\;\begin{aligned}
  \nu_{\text{grad}}(R_*) &\;=\; -\int_{R_*}^\infty e^{-x(r)}\,\frac{M}{r(r - 2M)}\,dr \;=\; -\int_{R_*}^\infty \frac{M}{r(r - 2M)}\,dr \;+\; \int_{R_*}^{R_* + L}\bigl(1 - e^{-x(r)}\bigr)\,\frac{M}{r(r - 2M)}\,dr \\[6pt]
  &\;=\; \frac{1}{2}\ln\!\left(1 - \frac{2M}{R_*}\right) \;+\; \Delta\nu_{\text{ext}}, \qquad \Delta\nu_{\text{ext}} \;\equiv\; \int_{R_*}^{R_* + L}\bigl(-\operatorname{expm1}(-x(r))\bigr)\,\frac{M}{r(r - 2M)}\,dr \;>\; 0,
\end{aligned}\;}$$
which `exterior_lapse()` evaluates via 200-point Gauss–Legendre quadrature on $[R_*, R_* + L]$.
In the stellar interior ($0 \le r \le R_*$), both `mass` and `gradient` modes satisfy the relativistic Bernoulli identity $\nu'(r) = -h'(r)$, yielding the interior lapse profile without any additional numerical quadrature (`tal_profiles.py`, lines 86–87):
$$\boxed{\;\nu(r) \;=\; \nu(R_*) \;-\; \bigl[h(r) - h(R_*)\bigr], \qquad \forall\,r \in [0, R_*].\;}$$

---

## 5. Checks

1. **Physical Verification of the Sign and Limits of $\Delta\nu_{\text{ext}} > 0$:**
   - Because $x(r) > 0$ on $r \in [R_*, R_* + L)$ and $r > 2M$, the integrand $(1 - e^{-x(r)})\frac{M}{r(r-2M)} > 0$, so $\Delta\nu_{\text{ext}} > 0$, meaning `gradient` mode softens the surface gravitational potential depth ($|\nu_{\text{grad}}(R_*)| < |\nu_{\text{GR}}(R_*)|$). As $L \to 0^+$, the integration interval length $L \to 0$ and $\Delta\nu_{\text{ext}} \to 0$, continuously recovering the Schwarzschild limit $\checkmark$.
2. **Verification of the Derivative of the Interior Bernoulli Integral ($\nu(r) = \nu(R_*) - [h(r) - h(R_*)]$)**:
   - In both modes, $\frac{d\nu}{dr} = -\frac{dh}{dr} = -2r\frac{dh_s}{ds}(r^2) = f_{\text{grad}}(r)A_{\text{GR}}(r)$, which gives $\nu'(0) = 0$ at $r = 0$ and $\nu(R_*)$ at $r = R_*$ $\checkmark$.

---

## 6. Explicit Comparison with the V23 Local Yukawa (`yukawa_action_el`, Branch II-Fluid) Branch (Point 5)

For audit clarity, we strictly separate the properties of the V21 phenomenological 3D top-hat branch (`tophat_mass`, `tophat_gradient` — **Branch I-A**) derived above from the V23 local Yukawa-screened Einstein–scalar–fluid branch (`yukawa_action_el` — **Branch II-Fluid**, see `step4`, `step7`, and `step22`):

| Physical Property | V21 `tophat_mass` / `tophat_gradient` (**Branch I-A**, Step 21) | V23 `yukawa_action_el` (**Branch II-Fluid**, Steps 4, 7, and 22) |
|---|---|---|
| **Interior Hydrostatic / Bernoulli Equation ($0 \le r \le R_*$)** | $h'(r) = -\nu'(r) \implies \nu(r) = \nu(R_*) - [h(r) - h(R_*)]$ ($f_L$ enters directly into $\nu' = f_L A_{\text{GR}}$ or $m' = 4\pi r^2 f_L\varepsilon$) | $h'(r) = -\nu'(r) + \phi_L'(r) \implies \nu(r) - \phi_L(r) = \nu(R_*) - \phi_L(R_*) - [h(r) - h(R_*)]$ (distinct scalar fifth force $\phi_L'$) |
| **Structure of Exterior Domain ($r > R_*$)** | Compact-support 3D top-hat halo on $r \in [R_*, R_* + L]$; for $r \ge R_* + L$, $x(r) \equiv 0$, $f_L(r) \equiv 1$, exact Schwarzschild vacuum | Exponentially decaying gravitating scalar tail over $r \in [R_*, \infty)$: $\rho_\phi = \frac{1}{2\alpha_Y}[B(\phi_L')^2 + \phi_L^2/L^2] > 0$ |
| **Surface Mass ($M_{\text{surf}}$) vs. ADM Mass ($M_{\text{ADM}}$)** | $m'(r > R_*) = 0 \implies M_{\text{ADM}} \equiv M_{\text{surf}} = m(R_*)$ | $m'(r > R_*) = 4\pi r^2\rho_\phi(r) > 0 \implies M_{\text{ADM}} = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}} > M_{\text{surf}}$ (see `step22`) |

> **V25 Update (Curved `yukawa_action_el` Branch):** 4-component exterior BVP ($\delta\nu_R=1.63\times 10^{-14}$), exact exterior linearization (`a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0`) and interior left-hand side scalar coupling (`scalar_m_sign = -1.0`: `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n`), dynamical scalar wave, and GL2 root convergence ($p_{\text{eff,root}}=4.248$). For SLy $L=10\;\text{km}$, $\alpha_Y=1$: $M_{\text{ADM}}$ turning point $15.316404$, adaptive Gauss–Kronrod $M_B$ turning point $15.316409$, main GL2 operator root (`scalar_m_sign=-1.0`) $15.316410$ (within $+5.80\times 10^{-6}\;\text{dex}$ of the $M_{\text{ADM}}$ peak; with diagnostic `scalar_m_sign=+1.0`, $15.319129$, $+2.725\times 10^{-3}\;\text{dex}$). Detailed tables and diagnostics: `step22`, Sections 4 and 7.
