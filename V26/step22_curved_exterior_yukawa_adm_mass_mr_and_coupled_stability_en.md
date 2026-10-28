# Step 22 (V25): Curved Exterior Vacuum Boundary-Value Problem of `yukawa_action_el`, Exact ADM Mass ($M_{\text{ADM}}$), Recomputed Mass–Radius ($M_{\text{ADM}}$–$R_*$) Sequences, and Coupled $(\delta h, \delta\phi_L, \delta m)$ Radial Pulsation Stability

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-27  
**File Path:** `scalar-tensor-yukawa-fit/step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` – `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse_en.md`  
**License:** CC BY-NC 4.0  
**Local Sources and Runnable Code:**  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` (Sections 1.2 and 2, lines 65–153; DOI: `10.5281/zenodo.23639084`)  
- Runnable solver code: `scalar-tensor-yukawa-fit/curved_yukawa_adm_and_stability.py` (and integrated into `covariant_closure_50aa2db0.py`)  
- Generated output files: `curved_yukawa_adm_sequences.csv`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_audit_summary.json`, `closure_tov_comparison.csv`, `conservation_Residuals.json`

---

## 1. Exact Exterior Coupled Einstein–Scalar Vacuum Equations ($r \ge R_*$) and Curved Robin Matching (Point 3)

In the earlier V23 (`50aa2db0`) `yukawa_action_el` computation (`covariant_closure_50aa2db0.py`, `solve_yukawa_to_surface`), the interior boundary-value problem was closed at the stellar surface $r = R_*$ using the flat-space Yukawa decay relation $\phi_{L,\text{flat}}'(R_*) = -\left(\frac{1}{R_*} + \frac{1}{L}\right)\phi_L(R_*)$, and the mass–radius sequences were optimized for the surface mass $M_{\text{surf}} \equiv m(R_*)$.

In reality, in the exterior vacuum outside the star ($r \ge R_*$, where $\varepsilon(r) = P(r) = 0$), the gravitational field is strongly curved ($B(R_*) \equiv e^{-2\lambda(R_*)} = 1 - \frac{2M_{\text{surf}}}{R_*} \approx 0.50\text{--}0.62$), and the energy density and radial pressure of the scalar field $\phi_L(r)$ continue to source the metric components $(m(r), \nu(r))$:
$$\boxed{\;\rho_\phi(r) \;=\; \frac{1}{2\alpha_Y}\left[B(r)\bigl(\phi_L'(r)\bigr)^2 + \frac{\phi_L(r)^2}{L^2}\right], \qquad p_{r,\phi}(r) \;=\; \frac{1}{2\alpha_Y}\left[B(r)\bigl(\phi_L'(r)\bigr)^2 - \frac{\phi_L(r)^2}{L^2}\right], \qquad B(r) \;=\; 1 - \frac{2m(r)}{r}.\;}$$
Substituting the exact algebraic identity of the vacuum Einstein equations ($\varepsilon = P = 0$),
$$B(r)\nu'(r) + \frac{B'(r)}{2} \;=\; \frac{2m(r)}{r^2} + 4\pi r\bigl(p_{r,\phi}(r) - \rho_\phi(r)\bigr) \;=\; \frac{2m(r)}{r^2} - \frac{4\pi r}{\alpha_Y L^2}\,\phi_L(r)^2,$$
into the radial Klein–Gordon equation $B\phi_L'' + \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L' - \frac{\phi_L}{L^2} = 0$, we obtain the **exact coupled exterior Einstein–scalar ODE system on $r \in [R_*, \infty)$**:
$$\boxed{\;\begin{aligned}
  \left(1 - \frac{2m(r)}{r}\right)\phi_L''(r) \;+\; \left[\frac{2}{r}\left(1 - \frac{m(r)}{r}\right) - \frac{4\pi r}{\alpha_Y L^2}\,\phi_L(r)^2\right]\phi_L'(r) \;-\; \frac{\phi_L(r)}{L^2} &\;=\; 0, \\[4pt]
  \frac{dm}{dr} &\;=\; 4\pi r^2\,\rho_\phi(r) \;=\; \frac{2\pi r^2}{\alpha_Y}\left[\left(1 - \frac{2m(r)}{r}\right)\bigl(\phi_L'(r)\bigr)^2 + \frac{\phi_L(r)^2}{L^2}\right], \\[4pt]
  \frac{d\nu}{dr} &\;=\; \frac{m(r) + 4\pi r^3\,p_{r,\phi}(r)}{r\bigl(r - 2m(r)\bigr)}, \qquad \lim_{r\to\infty}\nu(r) \;=\; 0.
\end{aligned}\;}$$

### Numerical Verification (`curved_yukawa_adm_and_stability.py`): Curved vs. Flat Exterior Robin Slope and $R_{\max}$-Convergence
The script `curved_yukawa_adm_and_stability.py` solves the interior + exterior two-domain boundary-value problem ($s = r^2$ regular core up to $r = R_*$, coupled continuously in $C^1$ to the 4-component exterior curved vacuum $y = [\phi_L, \phi_L', m, \nu]$ on $r \in [R_*, R_{\max}]$ with asymptotic WKB boundary condition $\phi_L'(R_{\max}) + \left(\frac{1}{L} + \frac{1 + m(R_{\max})/L}{R_{\max}}\right)\phi_L(R_{\max}) = 0$ and $\nu(R_{\max}) = \frac{1}{2}\ln(1 - 2m(R_{\max})/R_{\max})$). For the SLy benchmark configuration ($\log_{10}\rho_{0c} = 15.0$, $L = 10\;\text{km}$, $\alpha_Y = 1$, $R_* = 11.46572423\;\text{km}$, $M_{\text{surf}} = 1.50766122\,M_\odot$, $\phi_L(R_*) = 0.0054086823$, $B(R_*) = 0.6119$; `curved_yukawa_adm_sequences.csv` and `curved_yukawa_audit_summary.json`):
- **$R_{\max}$-independence across $N_{\text{mult}} \in \{10, 12, 14, 16\}$ ($R_{\max} = R_* + N_{\text{mult}}L$):**
  - Relative spread of the curved Dirichlet–Neumann slope: $\delta k_{\text{curved}}/k_{\text{curved}} = \mathbf{6.17\times 10^{-12}}$,
  - Relative spread of the ADM mass: $\delta M_{\text{ADM}}/M_{\text{ADM}} = \mathbf{5.54\times 10^{-13}}$,
  - Spread of the surface lapse $\nu(R_*)$: $\delta\nu_R = \mathbf{1.63\times 10^{-14}}$.
- **Curved vs. flat surface slope ($\log_{10}\rho_{0c} = 15.0$, $R_* = 11.46572423\;\text{km}$):**
  - Exact curved exterior Dirichlet–Neumann slope: $k_{\text{curved}} = -\phi_{L,\text{curved}}'(R_*)/\phi_L(R_*) = \mathbf{0.23354233\;\text{km}^{-1}}$ ($\phi_{L,\text{curved}}'(R_*) = -1.26315627\times 10^{-3}\;\text{km}^{-1}$),
  - Flat Robin approximation: $k_{\text{flat}} = \frac{1}{R_*} + \frac{1}{L} = \mathbf{0.18721647\;\text{km}^{-1}}$ ($\phi_{L,\text{flat}}'(R_*) = -1.01259441\times 10^{-3}\;\text{km}^{-1}$),
  - Relative difference in the surface derivative: **$\mathbf{+24.74\%}$** (spacetime curvature $B(R_*)^{-1/2} \approx 1.278$ steepens the exterior Yukawa decay rate).
  - Consequently, the exact curved exterior scalar halo mass is $\Delta M_{\phi,\text{ext}}^{\text{curved}} = \mathbf{0.00233692\,M_\odot}$ at $\log_{10}\rho_{0c} = 15.0$ (giving $M_{\text{ADM}}(15.0) = \mathbf{1.50999814\,M_\odot}$, whereas the flat Robin estimate gave $\Delta M_{\phi,\text{ext}}^{\text{flat}} = 0.003416\,M_\odot$), and $\Delta M_{\phi,\text{ext}}^{\text{curved}} = \mathbf{0.00410774\,M_\odot}$ at the ADM mass maximum (whereas the flat estimate gave $0.00898\,M_\odot$).

| $N_{\text{mult}}$ | $R_{\max}$ (km) | $k_{\text{curved}}$ ($\text{km}^{-1}$) | $\phi_{L,\text{curved}}'(R_*)$ ($\text{km}^{-1}$) | $\Delta M_{\phi,\text{ext}}^{\text{curved}}$ ($M_\odot$) | $M_{\text{ADM}}$ ($M_\odot$) | $\nu_R \equiv \nu_{\text{exact}}(R_*)$ |
|---|---|---|---|---|---|---|
| `10` | `111.4657` | `0.233542329280` | `-0.001263156266` | `0.002336923550` | `1.509998138706` | `-0.246306687228` |
| `12` | `131.4657` | `0.233542329278` | `-0.001263156266` | `0.002336923550` | `1.509998138707` | `-0.246306687228` |
| `14` | `151.4657` | `0.233542329278` | `-0.001263156266` | `0.002336923550` | `1.509998138707` | `-0.246306687228` |
| `16` | `171.4657` | `0.233542329278` | `-0.001263156266` | `0.002336923550` | `1.509998138707` | `-0.246306687228` |

---

## 2. Exact Exterior Lapse ($\nu(r)$) and Asymptotic Redshift Synchronization at $r = R_*$ (Point 5)

In the `yukawa_action_el` model (**Branch II-Fluid**), the exterior lapse $\nu(r)$ is **not** equal to the vacuum Schwarzschild formula $\frac{1}{2}\ln\bigl(1 - 2M_{\text{surf}}/r\bigr)$, because in the exterior domain $r \in [R_*, \infty)$ both $m(r)$ and $p_{r,\phi}(r)$ vary with radius. Integrating the 4-component exterior BVP $y = [\phi_L, \phi_L', m, \nu]$ with outer boundary condition $\nu(R_{\max}) = \frac{1}{2}\ln\bigl(1 - 2m(R_{\max})/R_{\max}\bigr)$ from $R_{\max} = R_* + N_{\text{mult}}L$ down to the stellar surface $r = R_*$ (defined by $h(R_*) = P(R_*) = 0$), the exact surface lapse is:
$$\boxed{\;\nu_{\text{exact}}(R_*) \;\equiv\; \nu_R \;=\; \frac{1}{2}\ln\!\left(1 - \frac{2m(R_{\max})}{R_{\max}}\right) \;-\; \int_{R_*}^{R_{\max}} \frac{m(r) + 4\pi r^3\,p_{r,\phi}(r)}{r\bigl(r - 2m(r)\bigr)}\,dr.\;}$$
For the SLy ($\log_{10}\rho_{0c} = 15.0$, $L = 10\;\text{km}$, $\alpha_Y = 1$) star at the exact surface boundary point $R_* = 11.46572423\;\text{km}$ ($M_{\text{surf}} = 1.50766122\,M_\odot$, $M_{\text{ADM}} = 1.50999814\,M_\odot$; `curved_yukawa_adm_sequences.csv` column `nu_R_15p0` and `curved_yukawa_audit_summary.json` block `rmax_independence_checks_SLy_15p0_L10km`), the three definitions at $r = R_*$ yield:
- $\nu_{\text{Schw}}(M_{\text{surf}}) = \frac{1}{2}\ln\bigl(1 - 2M_{\text{surf}}/R_*\bigr) = -0.2457823130$,
- $\nu_{\text{Schw}}(M_{\text{ADM}}) = \frac{1}{2}\ln\bigl(1 - 2M_{\text{ADM}}/R_*\bigr) = -0.2462745915$,
- $\nu_{\text{exact}}(R_*) \equiv \nu_R(15.0) = \mathbf{-0.2463066872}$ (the curved exterior scalar tail shifts the surface lapse at $r = R_*$ by $-5.2437\times 10^{-4}$ relative to $\nu_{\text{Schw}}(M_{\text{surf}})$ and by $-3.2096\times 10^{-5}$ relative to $\nu_{\text{Schw}}(M_{\text{ADM}})$, and the 4-component exterior BVP matches backward quadrature to $\delta\nu_R = 1.63\times 10^{-14}$).

Inside the star ($0 \le r \le R_*$), the modified Euler equation $h'(r) = -\nu'(r) + \phi_L'(r)$ integrates directly to the exact modified Bernoulli first integral:
$$\boxed{\;h(r) + \nu(r) - \phi_L(r) \;=\; \nu_{\text{exact}}(R_*) - \phi_L(R_*) \;=\; \text{const.}\;}$$

---

## 3. Recomputing the $M_{\text{ADM}}\text{–}R_*$ Sequences with the Asymptotic ADM Mass ($M_{\text{ADM}}$) and Baryonic Rest Mass ($M_B$)

The table below (`curved_yukawa_adm_sequences.csv` and `closure_tov_comparison.csv`) reports the **ADM-mass ($M_{\text{ADM}}$) optimized** quantities across all 8 `yukawa_action_el` sequences, together with the surface mass ($M_{\text{surf}}$) and the maximum and turning point of the conserved baryonic rest mass ($M_B^{\text{quad}} = \int_0^{R_*} 4\pi r^2 B(r)^{-1/2}\rho_0(r)\,dr$, evaluated via PCHIP + adaptive Gauss–Kronrod quadrature in $s=r^2$ with `epsrel=1e-12`):

| EOS | $L$ (km) | $\alpha_Y$ | $M_{\text{ADM},\max}$ ($M_\odot$) | $M_{\text{surf},\max}$ ($M_\odot$) | $M_{B,\max}^{\text{quad}}$ ($M_\odot$) | $\Delta M_{\phi,\text{ext}}^{\text{curved}}(M_{\max})$ ($M_\odot$) | $\log_{10}\rho_{0c}(M_{\text{ADM},\max})$ | $\log_{10}\rho_{0c}(M_{B,\max}^{\text{quad}})$ | $R(M_{\text{ADM},\max})$ (km) | $R(M_{\text{ADM}}=1.4\,M_\odot)$ (km) | $R(M_{\text{surf}}=1.4\,M_\odot)$ (km) | $M_{\text{ADM}}(15.0)$ ($M_\odot$) | $f_L(0; 15.0)$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **SLy** | `1.0` | `1.000` ($1/4\pi$) | `2.04684` | `2.04684` | `2.42634` | `0.000002` | `15.303086` | `15.303082` | `9.9689` | `11.6933` | `11.6933` | `1.56171` | `0.999198` |
| **SLy** | `1.0` | `0.333` ($1/12\pi$) | `2.04789` | `2.04789` | `2.42744` | `0.000001` | `15.302656` | `15.302660` | `9.9733` | `11.7018` | `11.7018` | `1.56386` | `0.999733` |
| **SLy** | `10.0` | `1.000` ($1/4\pi$) | **`2.03321`** | `2.02911` | **`2.42339`** | `0.004108` | **`15.316404`** | **`15.316409`** | `9.8853` | **`11.5359`** | `11.5348` | `1.51000` | `0.982378` |
| **SLy** | `10.0` | `0.333` ($1/12\pi$) | `2.04321` | `2.04183` | `2.42630` | `0.001381` | `15.307166` | `15.307169` | `9.9446` | `11.6475` | `11.6471` | `1.54638` | `0.993967` |
| **APR4** | `1.0` | `1.000` ($1/4\pi$) | `2.18693` | `2.18693` | `2.65730` | `0.000003` | `15.282808` | `15.282805` | `9.8947` | `11.3085` | `11.3085` | `1.56316` | `0.999201` |
| **APR4** | `10.0` | `1.000` ($1/4\pi$) | `2.18068` | `2.17624` | `2.66733` | `0.004449` | `15.293022` | `15.293025` | `9.8462` | `11.1830` | `11.1828` | `1.51271` | `0.982098` |
| **H4** | `1.0` | `1.000` ($1/4\pi$) | `2.01139` | `2.01138` | `2.29535` | `0.000001` | `15.209751` | `15.209739` | `11.6917` | `13.9388` | `13.9388` | `1.88120` | `0.999166` |
| **H4** | `10.0` | `1.000` ($1/4\pi$) | `1.97809` | `1.97502` | `2.26306` | `0.003077` | `15.227785` | `15.227787` | `11.5289` | `13.7179` | `13.7169` | `1.82871` | `0.981076` |

**Baryonic Mass Definition and Turning-Point Precision (`baryonic_mass_and_first_law_audit`):** In `curved_yukawa_adm_sequences.csv` and `closure_tov_comparison.csv`, the columns `M_baryon_max_Msun` and `log10_rho0c_at_M_baryon_max` are defined via adaptive Gauss–Kronrod quadrature on $s=r^2$ PCHIP-interpolated interior profiles ($M_B^{\text{quad}}$, `epsrel=1e-12`) with a 7-point quadratic fit over $\pm 0.002\;\text{dex}$. For SLy $L=10\;\text{km}$, $\alpha_Y=1$, $M_{B,\max}^{\text{quad}} = \mathbf{2.42338503\,M_\odot}$ at $\log_{10}\rho_{0c}(M_{B,\max}^{\text{quad}}) = 15.31640854$ (converging across $N_{\text{grid}} = 301\to 2401$ within $1.2\times 10^{-6}\;\text{dex}$), refined trapezoidal quadrature at $N_{\text{grid}}=2401$ gives $M_B^{\text{trapz}} = 2.42338667\,M_\odot$ with a 7-point peak at $15.31641814$ ($+9.6\times 10^{-6}\;\text{dex}$), whereas single-pass integration of $M_B^{\text{ODE}}$ as a 5th ODE state variable gives $2.42203137\,M_\odot$ and a 3-point peak at $15.31652192$ (shifted by $-0.001354\,M_\odot$ and $+1.13\times 10^{-4}\;\text{dex}$ due to step-size sensitivity of the steep $\sqrt{s}\,\rho_0(h)$ integrand near the crust edge $h\to 0$).

*(Comparison with the earlier flat-Robin `yukawa_action_el` values for SLy $L = 10\;\text{km}$, $\alpha_Y = 1$: under the flat-Robin approximation, $M_{\text{surf},\max}^{\text{flat}} = 2.03583\,M_\odot$ and $M_{\text{surf},\max}^{\text{flat}} + \Delta M_{\phi,\text{ext}}^{\text{flat}} = 2.04481\,M_\odot$ at $\log_{10}\rho_{0c} = 15.31307$; with the exact curved exterior BVP solution, the surface mass peak is $M_{\text{surf},\max}^{\text{curved}} = 2.02911\,M_\odot$ and the exact ADM mass peak is $M_{\text{ADM},\max}^{\text{curved}} = \mathbf{2.03321\,M_\odot}$ at $\log_{10}\rho_{0c} = \mathbf{15.316404}$.)*

---

## 4. Coupled $(\delta h, \delta\phi_L, \delta m)$ Radial Pulsation Eigenvalue Problem on the Curved `yukawa_action_el` Background (Point 4)

To ensure that the radial stability analysis covers not only the V21 nonlocal 3D top-hat background (`step18`), we solved the coupled $(\delta h, \delta\phi_L, \delta m)$ linear perturbation system derived in Step 12 directly on the curved-exterior `yukawa_action_el` backgrounds (`curved_yukawa_adm_and_stability.py`, `solve_yukawa_coupled_omega2_gl2`):
- **Interior coupled equations ($0 \le r \le R_*$):** $(\delta h, \delta m, \delta\phi_L, \delta\phi_L')$ integrated with regular center expansion ($\delta h(r) = \delta h_0 + \mathcal{O}(r^2)$, $\delta\phi_L(r) = \delta\phi_c + \mathcal{O}(r^2)$, $\delta m(r) = \mathcal{O}(r^3)$), where the fluid inertia term $\omega^2 e^{2(\lambda-\nu)}\xi(r)$ and the dynamical scalar wave term $+\omega^2 e^{-2\nu}\delta\phi_L$ are closed via the exact adiabatic relation $(\varepsilon+P)r\bigl(\delta\nu' - \delta\phi_L' + \delta h'\bigr) = \bigl[\text{static background derivative}\bigr]\,\xi(r)$.
- **Moving surface and curved exterior tail ($R_* \le r \le R_{\max}$):** at $r = R_*^-$, because $\varepsilon(R_*^-) = P(R_*^-) = 0$ on the polytropic crust, the jump in every perturbation variable vanishes ($\llbracket\delta\phi_L\rrbracket = \llbracket\delta\phi_L'\rrbracket = \llbracket\delta m\rrbracket = 0$). In the exterior domain $r \in [R_*, R_{\max}]$, the linearized exterior $(\delta\phi_L, \delta m)$ system is integrated out to $R_{\max} = R_* + 14L$, and the eigenvalue $\omega_0^2$ is determined by the vanishing of the asymptotic ADM mass perturbation $\delta M_{\text{ADM}}(\omega^2) = \delta m(R_{\max}; \omega^2) = 0$.

### Numerical Results (`curved_yukawa_coupled_stability.csv` and `curved_yukawa_audit_summary.json`): Dynamical vs. Quasi-Static vs. Frozen-$\phi_L$ ($\delta\phi_L = 0$) Zero Crossings

**Methodological Summary of `curved_yukawa_adm_and_stability.py`:**
- **4-component exterior BVP (`solve_curved_yukawa_star`, lines 234–364):** $y = [\phi_L, \phi_L', m, \nu]$ integrated simultaneously on $r \in [R_*, R_{\max}]$ with boundary conditions $\phi_L(R_*)=\phi_{\text{surf}}$, $m(R_*)=M_{\text{surf}}$, $\phi_L'(R_{\max}) + \bigl(1/R_{\max} + 1/(L\sqrt{B(R_{\max})})\bigr)\phi_L(R_{\max}) = 0$, and $\nu(R_{\max}) = \frac{1}{2}\ln\bigl(1 - 2m(R_{\max})/R_{\max}\bigr)$. Consistency between the backward quadrature and the 4th BVP component at $r=R_*$ is $\delta\nu_R = 1.63\times 10^{-14}$.
- **Exterior linearization and finite-difference derivative tests (SLy, $\log_{10}\rho_{0c}=15.0$, $L=10\;\text{km}$, $\alpha_Y=1$):** `linearized_exterior_response` (line 593: `a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0`) and the left-hand side ($L_u u = R_h\,\delta h$) interior scalar operator coefficient `a_m_sc` (`scalar_m_sign = -1.0`: `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n`, obtained by moving the right-hand side derivative $\partial(\phi_L'')/\partial m = +(2\phi_L''/r + 2\phi_L'/r^2)/B$ to the left-hand side) verified against central finite differences:
  - Exterior response $\partial\phi_L'(R_*)/\partial m(R_*)$: linearized BVP **`-1.68365059e-04`** $\text{km}^{-2}$ vs. 4-component exterior BVP central finite difference **`-1.68365059e-04`** $\text{km}^{-2}$ (relative difference: **$3.98\times 10^{-10}$**; $K_\phi = 0.23358323\;\text{km}^{-1}$, rel. diff. $1.42\times 10^{-10}$; $\partial M_{\text{ADM}}/\partial m(R_*) = 0.99947576$, rel. diff. $4.26\times 10^{-10}$).
  - Interior $r = 6.0\;\text{km}$ scalar-mass derivative ($\partial\phi_L''/\partial m$): analytic **`-1.72833429e-04`** $\text{km}^{-3}$ vs. interior scalar RHS central finite difference **`-1.72833428e-04`** $\text{km}^{-3}$ (relative difference: **$6.08\times 10^{-9}$**).
- **Dynamical scalar block:** the scalar block includes the $-\omega^2 e^{-2\nu}\delta\phi_L$ term and the exterior response is obtained from an $\omega$-dependent BVP (`dynamical_scalar_wave=True`); solver: GL2 piecewise segment Volterra discretization, $N=480$, with the exact left-hand side sign `scalar_m_sign = -1.0`. The quasi-static and dynamical roots agree to all 6 reported decimal places (since at $\omega^2=0$ the two operators are identical).
- **Turning points and zero crossings:** `minimize_scalar` + 7-point quadratic fit over $\pm 0.002\;\text{dex}$; zero crossings from a 3-point quadratic fit over $\pm 0.0005\;\text{dex}$ around the root.

| Sequence | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dynamical $\omega_0^2=0$ | $\Delta_{\text{ADM}}$ Dynamical (dex) | $\Delta_{M_B}$ Dynamical (dex) | Quasi-Static $\omega_0^2=0$ | Frozen $\delta\phi_L=0$ | $\Delta_{\text{ADM}}$ Frozen (dex) |
|---|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316409` | **`15.316410`** | **`+5.80e-06`** | `+1.33e-06` | `15.316410` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.307169` | **`15.307170`** | **`+4.20e-06`** | `+6.41e-07` | `15.307170` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303082` | **`15.303080`** | **`-6.86e-06`** | `-2.61e-06` | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.293025` | **`15.293025`** | **`+3.57e-06`** | `-1.40e-07` | `15.293025` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227787` | **`15.227783`** | **`-1.34e-06`** | `-4.05e-06` | `15.227783` | `15.251384` | `+2.36e-02` |

**Convergence Study on the $\omega_0^2=0$ Zero Crossing Itself (with respect to $N$ and $R_{\max}$)** (SLy, $L=10\;\text{km}$, $\alpha_Y=1$; together with the convergence of $\omega_0^2$ at $\log_{10}\rho_{0c}=15.311404$):
- **Zero-crossing $N$-convergence ($N_{\text{mult}}=14$):** $p_{\text{eff,root}}=4.248$ (from the three finest levels), Richardson-extrapolated root $\log_{10}\rho_{0c,\text{Rich}}(\omega_0^2=0) = `15.31640981`$; at $N=480$, the discretization error relative to the Richardson root is $6.72\times 10^{-8}\;\text{dex}$ (and between $N=240\to 800$ across all 5 sequences $\le 2.02\times 10^{-5}\;\text{dex}$).
- **Zero-crossing $R_{\max}$-convergence ($N=480$, $R_{\max}=R_*+N_{\text{mult}}L$, $N_{\text{mult}}\in\{10,12,14,16,20\}$, $R_{\max}\in[109.89,\;209.89]\;\text{km}$):** the spread of the root is $2.42\times 10^{-13}\;\text{dex}$ (across all 5 sequences $\le 2.72\times 10^{-13}\;\text{dex}$).
- **Fixed-point ($\log_{10}\rho_{0c}=15.311404$) $\omega_0^2$ convergence:** $p_{\text{eff}}=5.223$, Richardson value $\omega_0^2=6369417.611\;\text{s}^{-2}$ (at $N=480$, relative error $1.29\times 10^{-5}$).

| $N$ | $\log_{10}\rho_{0c}(\omega_0^2=0)$ | $|\Delta\log_{10}\rho_{0c}|$ from Richardson root (dex) | $\Delta_{\text{ADM}}$ (dex) | $\Delta_{M_B}$ (dex) | $\omega_0^2(15.311404)$ ($\text{s}^{-2}$) | Rel. error ($\omega_0^2$) |
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

**Interpretive Note:** With the exact left-hand side sign (`scalar_m_sign = -1.0`), the root of the GL2 operator converges with respect to $N$ and $R_{\max}$ to better than $<10^{-7}\;\text{dex}$ (`15.31640981`), the $\omega_0^2=0$ zero crossing lies within $\le 6.9\times 10^{-6}\;\text{dex}$ of the ADM mass maximum across all 5 sequences (for SLy $L=10\;\text{km}$, $\alpha_Y=1$: $15.316410$ vs. $15.316404$, $\Delta_{\text{ADM}} = +5.80\times 10^{-6}\;\text{dex}$, $\Delta_{M_B^{\text{quad}}} = +1.33\times 10^{-6}\;\text{dex}$), and the interior relative $L^2$ differences between the eigenfunction and the static equilibrium tangent for SLy $L=10\;\text{km}$, $\alpha_Y=1$ are $\delta h = 5.04\times 10^{-4}$, $\delta m = 1.77\times 10^{-4}$, $\delta\phi_L = 2.27\times 10^{-4}$, $\delta\phi_L' = 3.09\times 10^{-4}$ (while the relative difference of the free-surface displacements is $\Delta\xi(R_*)/\xi = 3.54\%$; see Section 7.4). Because the $3.54\%$ free-surface displacement difference and the ODE vs. quadrature baryonic-mass turning-point spread show that the minus-sign root and grid convergence do not by themselves constitute an independent proof of the full radial fundamental spectrum and free-surface variation, **the stability closing statement is formulated conditionally until the independent verification of the coupled second-variation derivation and interior–exterior boundary-condition matching is concluded.**

---

## 5. Checks

1. **Comparison of the ADM Mass ($M_{\text{ADM}}$), Baryonic Mass ($M_B$), and Coupled $\omega_0^2 = 0$ Zero Crossing:**
   - For SLy, $L = 10\;\text{km}$, $\alpha_Y = 1$:
     $$\log_{10}\rho_{0c}\!\left(\frac{dM_{\text{ADM}}}{d\rho_c}=0\right) = 15.316404, \quad \log_{10}\rho_{0c}\!\left(\frac{dM_B}{d\rho_c}=0\right) = 15.316409, \quad \log_{10}\rho_{0c}\bigl(\omega_0^2=0\bigr) = 15.316410,$$
     meaning the coupled dynamical stability zero crossing coincides with the ADM mass peak within $+5.80\times 10^{-6}\;\text{dex}$ (and with the baryonic mass peak within $+1.33\times 10^{-6}\;\text{dex}$), while across the remaining 4 sequences the offset is $+4.20\times 10^{-6}\;\text{dex}$ ($\alpha_Y=1/3$), $-6.86\times 10^{-6}\;\text{dex}$ ($L=1\;\text{km}$), $+3.57\times 10^{-6}\;\text{dex}$ (APR4), and $-1.34\times 10^{-6}\;\text{dex}$ (H4) $\checkmark$.
2. **Physical Effect of Omitting the Scalar-Field Perturbation ($\delta\phi_L$):**
   - If the scalar field is frozen on the `yukawa_action_el` background ($\delta\phi_L \equiv 0$), the zero crossing for SLy $L=10\;\text{km},\;\alpha_Y=1$ shifts to $\log_{10}\rho_{0c} = 15.328565$ ($+1.22\times 10^{-2}\;\text{dex}$ above the ADM peak), because without the compressional perturbation of the attractive scalar self-energy ($\delta\phi_L > 0$) the star appears artificially stiffer. Including the coupled $(\delta h, \delta\phi_L, \delta m)$ system reduces this offset by a factor of **$2097\times$** (to $+5.80\times 10^{-6}\;\text{dex}$, and by factors of $165\times\text{--}17593\times$ across all 5 sequences) $\checkmark$.

---

## 6. Reproducibility and Generated Files

All results of Step 22 can be reproduced with a single command:
```bash
python3 curved_yukawa_adm_and_stability.py
```
Running the script generates and updates the following files (in the root directory and in the `csv/` and `data/v23/` subdirectories):
- `curved_yukawa_adm_sequences.csv`: summary table of the 8 mass–radius sequences recomputed for the ADM mass ($M_{\text{ADM}}$), surface mass ($M_{\text{surf}}$), and three baryonic rest-mass definitions ($M_B^{\text{quad}}$, $M_B^{\text{ODE}}$, $M_B^{\text{trapz}}$).
- `curved_yukawa_coupled_stability.csv`: coupled $(\delta h, \delta\phi_L, \delta m)$ dynamical, quasi-static, and frozen-$\phi_L$ ($\delta\phi_L = 0$) radial eigenvalue scans and zero crossings across the 5 `yukawa_action_el` sequences (generated with the `scalar_m_sign = -1.0` main operator).
- `curved_yukawa_audit_summary.json`: machine-readable audit log of $R_{\max}$-convergence, curved vs. flat Robin slope, finite-difference derivative tests, the 108-point $4\times 4$ Jacobian verification, the baryonic mass / First-Law audit, and the eigenmode vs. static tangent comparison.
- `closure_tov_comparison.csv` and `conservation_Residuals.json`: V23 reference files updated with the `yukawa_action_el_curved_adm` rows and the `curved_exterior_yukawa_adm_and_coupled_stability` block.

---

## 7. V25 Diagnostic Analysis: Sign Fixes, the Three New Diagnostics, and the Conditional Stability Closing Statement (`curved_yukawa_audit_summary.json`)

**Publication Status and Version Delineation (V24 vs. V25):** The earlier Zenodo publication (**V24**, DOI: `10.5281/zenodo.23680715`) refers exclusively to the earlier V24 release. The sign fixes in `curved_yukawa_adm_and_stability.py` (`a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0` in `linearized_exterior_response` and making the exact left-hand side sign `scalar_m_sign = -1.0`, i.e., `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n`, the default in the $L_u u = R_h\,\delta h$ operator of `solve_yukawa_coupled_omega2_gl2`), the three new diagnostics in `curved_yukawa_audit_summary.json` (`full_jacobian_multi_eos_and_points`, `baryonic_mass_and_first_law_audit`, `eigenmode_xi_nodes_and_static_tangent_audit`), and the closure of the theoretical verification belong exclusively to the new **V25** version — the earlier V24 DOI must never be cited as if these fixes and closures had already been included in it.

### 7.1. Finite-Difference Derivative Verification (`derivative_verification_tests_SLy_15p0_L10km`)
Finite-difference verification of the exterior `linearized_exterior_response` and interior scalar-mass coupling (SLy, $\log_{10}\rho_{0c}=15.0$, $L=10\;\text{km}$, $\alpha_Y=1$):

| Quantity | Linearized BVP / Analytic | Central Finite Difference (FD) | Relative Difference |
|---|---|---|---|
| `exterior_dphi_prime_dm_R` ($\partial\phi_L'(R_*)/\partial m(R_*)$, $\text{km}^{-2}$) | `-1.68365059e-04` | `-1.68365059e-04` | `3.98e-10` |
| `exterior_K_phi_R` ($-\partial\phi_L'(R_*)/\partial\phi_L(R_*)$, $\text{km}^{-1}$) | `2.33583231e-01` | `2.33583231e-01` | `1.42e-10` |
| `exterior_M_m_R` ($\partial M_{\text{ADM}}/\partial m(R_*)$) | `9.99475763e-01` | `9.99475764e-01` | `4.26e-10` |
| `exterior_M_phi_R` ($\partial M_{\text{ADM}}/\partial\phi_L(R_*)$, $\text{km}$) | `1.27572889e+00` | `1.27572919e+00` | `2.31e-07` |
| `interior_r6km_a_m_sc` ($\partial\phi_L''/\partial m$ at $r=6\;\text{km}$, $\text{km}^{-3}$) | `-1.72833429e-04` | `-1.72833428e-04` | `6.08e-09` |

### 7.2. Multi-EOS, Multi-Density, and Multi-Radius Full $4\times 4$ Jacobian Verification (`full_jacobian_multi_eos_and_points`)
The function `verify_full_jacobian_multi_eos_and_points()` verifies the full interior $4\times 4$ Jacobian matrix ($\delta h'$, $\delta\phi_L''$, $\delta m'$, $\delta\nu'$, and the coefficients of the adiabatic displacement $\xi$) against central finite differences across 4 configurations (SLy at $\log_{10}\rho_{0c}=15.0$ and the turning point $15.316404$; APR4 at $\log_{10}\rho_{0c}=15.0$ and the turning point $15.293022$), at 3 interior radii per configuration ($r \in \{3.0,\;6.0,\;8.5\}\;\text{km}$) plus surface matching. Across all **108 checks**, the maximum relative error is `max_relative_error_across_all_108_checks = 2.74e-05` (`all_checks_pass_tol_1e_4 = true`):

| EOS | $\log_{10}\rho_{0c}$ | $K_m$ rel. err | $K_\phi$ rel. err | $M_m$ rel. err | $M_\phi$ rel. err | Interior radii $r$ (km) | Max interior rel. err |
|---|---|---|---|---|---|---|---|
| SLy | `15.000000` | `3.98e-10` | `1.42e-10` | `4.26e-10` | `2.31e-07` | `3.0, 6.0, 8.5` | `2.74e-05` |
| SLy | `15.316404` | `7.12e-10` | `4.46e-10` | `4.64e-10` | `7.40e-08` | `3.0, 6.0, 8.5` | `1.88e-05` |
| APR4 | `15.000000` | `1.42e-09` | `6.91e-10` | `5.42e-10` | `8.70e-08` | `3.0, 6.0, 8.5` | `1.94e-05` |
| APR4 | `15.293022` | `5.79e-10` | `4.84e-10` | `6.21e-10` | `1.12e-07` | `3.0, 6.0, 8.5` | `1.59e-05` |

### 7.3. Baryonic Rest Mass ($M_B$) and First-Law Audit (`baryonic_mass_and_first_law_audit`)
The function `verify_baryonic_mass_and_first_law()` computes the baryonic rest mass ($M_B = \int_0^{R_*} 4\pi r^2 B(r)^{-1/2}\rho_0(r)\,dr$) via three independent methods: (i) post-facto trapezoidal rule (`M_B_trapz`), (ii) as the 5th state variable of the $(s=r^2)$ TOV ODE system (`dM_B/ds = 2*pi*sqrt(s)*rho0/sqrt(B)`, `rtol=1e-12`, `atol=1e-14`, `M_B_ode`), and (iii) adaptive Gauss–Kronrod quadrature on $s=r^2$ PCHIP-interpolated profiles (`quad`, `epsrel=1e-12`, `M_B_quad`). Grid convergence (SLy, $L=10\;\text{km}$, $\alpha_Y=1$):

| `ngrid` | 7-point peak $M_B^{\text{trapz}}$ | 7-point peak $M_B^{\text{quad}}$ | 7-point peak $M_{\text{ADM}}$ | 3-point peak $M_B^{\text{ODE}}$ |
|---|---|---|---|---|
| `301` | `15.316389` | `15.316409` | `15.316404` | `15.316522` |
| `601` | `15.316471` | `15.316410` | `15.316414` | `15.316522` |
| `1201` | `15.316433` | `15.316410` | `15.316406` | `15.316523` |
| `2401` | `15.316418` | `15.316410` | `15.316415` | `15.316522` |

Across all 8 sequences, the turning points of $M_B^{\text{quad}}$ (computed with adaptive Gauss–Kronrod quadrature in $s=r^2$) and $M_{\text{ADM}}$ agree to within $\le 8.43\times 10^{-6}\;\text{dex}$ (`delta_dex_ADM_vs_baryon_turning` between $-4.58\times 10^{-6}$ and $+8.43\times 10^{-6}\;\text{dex}$). The 9-point First-Law test around the turning point (`dM_ADM/dx - e^{\nu(R_*) - \phi_L(R_*)} dM_B/dx`) quantifies the exact residual arising from the non-conformal ($\mathcal{L}_m = -e^{-\phi_L}\varepsilon$, `USE_TRACE_SOURCE=False`) coupling:
$$\frac{dM_{\text{ADM}}}{dx} - e^{\nu(R_*) - \phi_L(R_*)}\frac{dM_B}{dx} \;=\; \int_0^{R_*} \mathcal{O}_{\text{nc}}(r)\,\frac{\partial\phi_L}{\partial x}\,dr \;\approx\; +8.23\times 10^{-5}\,M_\odot/\text{dex}.$$

### 7.4. Eigenmode Reconstruction ($\xi(r)$, Node Count) and Comparison with the Static Equilibrium Tangent (`eigenmode_xi_nodes_and_static_tangent_audit`)
The function `verify_eigenmode_xi_nodes_and_static_tangent()` reconstructs the coupled fundamental radial displacement profile $\xi(r)$, verifies the interior node count (`n_nodes = 0`) and the free-surface Lagrangian pressure condition ($\Delta P(R_*) = -(\varepsilon+P)r\,\delta h(R_*) = 0$), and compares the $\omega_0^2=0$ eigenmode against the static equilibrium tangent vector $\partial_{\log_{10}\rho_{0c}}(h, m, \phi_L, \phi_L')$ for both the main (`scalar_m_sign = -1.0`, exact LHS) and plus-sign (`scalar_m_sign = +1.0`) operator conventions:

| Sequence | Main Zero Crossing (`scalar_m_sign = -1.0`, exact LHS) | $\Delta_{\text{ADM}}$ (dex) | Diagnostic Zero Crossing (`scalar_m_sign = +1.0`) | $\Delta_{\text{ADM}}$ (dex) | $\xi(r)$ nodes (`-`) | $\xi(R_*)$ free surface vs. static tangent (km) | Surface $\xi(R_*)$ rel. diff. | $\xi/r$ central reg. |
|---|---|---|---|---|---|---|---|---|
| `SLy_L10km_alphaY_1.0000` | **`15.316410`** | **`+5.813e-06`** | `15.319129` | `+2.725e-03` | `0` | `-2.8994` vs. `-2.8003` | `3.541e-02` (`3.54%`) | `1.08e-03` |
| `SLy_L10km_alphaY_0.3333` | **`15.307170`** | **`+4.089e-06`** | `15.308114` | `+9.479e-04` | `0` | `-2.9237` vs. `-2.8986` | `8.666e-03` (`0.87%`) | `1.03e-03` |
| `SLy_L1km_alphaY_1.0000` | **`15.303080`** | **`-6.412e-06`** | `15.303088` | `+2.212e-06` | `0` | `-2.9406` vs. `-2.9411` | `1.767e-04` (`0.018%`) | `1.01e-03` |
| `APR4_L10km_alphaY_1.0000` | **`15.293025`** | **`+3.065e-06`** | `15.294175` | `+1.153e-03` | `0` | `-2.1155` vs. `-2.0681` | `2.289e-02` (`2.29%`) | `1.29e-03` |
| `H4_L10km_alphaY_1.0000` | **`15.227783`** | **`-1.658e-06`** | `15.234109` | `+6.324e-03` | `0` | `-5.5394` vs. `-5.5623` | `4.116e-03` (`0.41%`) | `5.30e-04` |

Relative $L^2$-norm difference between the static equilibrium tangent ($\partial_{\log_{10}\rho_{0c}}$) and the $\omega_0^2=0$ eigenfunction across the entire stellar volume:

| Sequence | Operator Variant | $\|\delta h - \partial h\|_2$ | $\|\delta m - \partial m\|_2$ | $\|\delta\phi_L - \partial\phi_L\|_2$ | $\|\delta\phi_L' - \partial\phi_L'\|_2$ |
|---|---|---|---|---|---|
| `SLy_L10km_alphaY_1.0000` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, main result) | `5.043e-04` | `1.767e-04` | `2.273e-04` | `3.093e-04` |
| `SLy_L10km_alphaY_1.0000` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnostic) | `5.026e-03` | `4.099e-03` | `9.347e-01` | `5.974e-01` |
| `SLy_L10km_alphaY_0.3333` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, main result) | `1.319e-04` | `1.588e-04` | `5.441e-04` | `3.635e-04` |
| `SLy_L10km_alphaY_0.3333` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnostic) | `1.525e-03` | `1.341e-03` | `9.397e-01` | `5.967e-01` |
| `SLy_L1km_alphaY_1.0000` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, main result) | `5.177e-05` | `1.638e-04` | `1.709e-04` | `2.478e-03` |
| `SLy_L1km_alphaY_1.0000` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnostic) | `4.052e-05` | `1.468e-04` | `2.243e-02` | `2.143e-02` |
| `APR4_L10km_alphaY_1.0000` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, main result) | `3.048e-04` | `2.108e-04` | `1.239e-03` | `6.328e-04` |
| `APR4_L10km_alphaY_1.0000` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnostic) | `4.077e-02` | `6.161e-03` | `1.112e+00` | `7.202e-01` |
| `H4_L10km_alphaY_1.0000` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, main result) | `7.340e-05` | `1.885e-04` | `6.043e-05` | `1.400e-04` |
| `H4_L10km_alphaY_1.0000` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnostic) | `1.059e-01` | `1.964e-03` | `5.713e-01` | `3.420e-01` |

> **Summary and Boundary of Theoretical Closure (`mathematical_origin_of_0p002725_dex_offset`):** In the nonlinear interior scalar equation ($B\phi_L'' + C_0\phi_L' = \phi_L/L^2 - \alpha_Y e^{-\phi_L}\varepsilon$), the partial derivative of the right-hand side (RHS) with respect to $m$ at fixed $(h, \phi_L, \phi_L')$ is $\partial(\phi_L'')/\partial m = +(2\phi_L''/r + 2\phi_L'/r^2)/B$ (at $r=6\;\text{km}$, $-1.728334\times 10^{-4}\;\text{km}^{-3}$). When rearranged onto the left-hand side ($L_u u = R_h\,\delta h$) of the Volterra system in `solve_yukawa_coupled_omega2_gl2`, the exact left-hand side sign is `scalar_m_sign = -1.0` (`a_m_sc = -\partial(\phi_L'')/\partial m`), with which all V25 main results (`curved_yukawa_coupled_stability.csv`, `curved_yukawa_audit_summary.json`) are generated: across all 5 sequences it yields a nodeless ($n=0$) $\xi(r)$ fundamental mode and an $\omega_0^2=0$ zero crossing (`15.316410`) lying within $\le 6.9\times 10^{-6}\;\text{dex}$ of the ADM mass peak (whereas `scalar_m_sign = +1.0` shifts the root by $+0.002725\;\text{dex}$ to `15.319129`). For SLy $L=10\;\text{km}$, $\alpha_Y=1$, the relative $L^2$ differences between the eigenfunction and the static equilibrium tangent are $\delta h = 5.04\times 10^{-4}$, $\delta m = 1.77\times 10^{-4}$, $\delta\phi_L = 2.27\times 10^{-4}$, $\delta\phi_L' = 3.09\times 10^{-4}$, while the relative difference of the free-surface displacements is $3.54\%$ (`rel_diff_xi_R_exact_minus_vs_static_tangent = 0.035406`). While the minus-sign root and grid convergence represent a major step forward, they do **not by themselves constitute an independent proof of the full radial fundamental spectrum and free-surface variation**, and therefore **the stability closing statement is formulated conditionally until the independent verification of the coupled variational derivation and boundary-condition matching is concluded.**
