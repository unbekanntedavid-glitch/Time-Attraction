# Step 23 (Module 1): GW170817 Tensor Wave Speed ($c_T = 1$), Massive Scalar Dispersion, and Gravitational Lensing ($\Phi + \Psi$, $\gamma_{\text{PPN}}(r)$, Curved-Spacetime Light Deflection $\hat{\alpha}(b)$, Shapiro Time Delay, and Weak-Lensing Convergence $\kappa$) at $L = 10\;\text{km}$, $\alpha_Y = 1$

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-27  
**Version:** V26  
**File Path:** `scalar-tensor-yukawa-fit/step23_gw170817_and_lensing_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` – `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability_en.md` (V24 Zenodo DOI: `10.5281/zenodo.23680715`; V25 Zenodo DOI: `10.5281/zenodo.23699685`)  
**License:** CC BY-NC 4.0  
**Local Sources and Runnable Code:**  
- Runnable solver code: `scalar-tensor-yukawa-fit/gw170817_and_lensing_audit.py` (synced to `code/gw170817_and_lensing_audit.py`)  
- Generated output files: `gw170817_dispersion_and_cT.csv`, `lensing_deflection_and_ppn_slip.csv`, `gw170817_and_lensing_summary.json` (synced to `csv/` and `data/v23/`)

---

## 1. Fundamental Action Principle and Coupling Dictionary ($L = 10\;\text{km}$, $\alpha_Y = 1$)

In the covariant scalar–tensor formulation of the Time Attraction Law (TAL, `step1`–`step22`, V25 DOI: `10.5281/zenodo.23699685`), the Einstein-frame gravitational and massive Yukawa scalar action reads:

$$S \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}\,R \;-\; \frac{c^4}{16\pi G\,\tilde{\alpha}_L}\left(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + 2\,V(\phi_L)\right)\right] \;+\; S_m\!\left[\tilde{g}_{\mu\nu},\,\Psi_m\right],$$

where the massive Yukawa screening potential is:

$$V(\phi_L) \;=\; \frac{1}{2}\,m_\phi^2\,\phi_L^2, \qquad m_\phi \;=\; \frac{1}{L} \;=\; 0.1\;\text{km}^{-1} \;=\; 10^{-4}\;\text{m}^{-1},$$

and the coupling dictionary (`step1`–`step5`, `step14`, `step22`) is:

$$\alpha_Y \;=\; 4\pi\,\alpha_{\text{std}} \;=\; 8\pi\,\tilde{\alpha}_L \;=\; 8\pi\,\beta^2, \qquad \gamma_L \;=\; \frac{G L^2}{c^2} \;=\; 7.42616027\times 10^{-20}\;\text{m}^3\,\text{kg}^{-1}.$$

**Single action, single field variable (V26 normalization fix).** All statements of Steps 23–25 are derived from the one action $S[\phi_L]$ above with $\tilde{\alpha}_L = \alpha_Y/(8\pi)$ and the matter coupling $S_m[A^2(\phi_L)g_{\mu\nu},\Psi_m]$, $A(\phi_L) = \exp(-\alpha_Y\phi_L)$, $d\ln A/d\phi_L = -\alpha_Y$ ($\phi_L > 0$ near matter, $g^{(J)}_{\mu\nu} = \exp(-2\alpha_Y\phi_L)g_{\mu\nu}$). If the kinetic prefactor inside $\frac{c^4}{16\pi G}[\dots]$ is changed from $1/\tilde{\alpha}_L$ to $2$ (Damour–Esposito-Farèse form), the field **must** be rescaled; the same symbol $\phi_L$ with unchanged $A = \exp(-\alpha_Y\phi_L)$ cannot be kept:

$$\frac{1}{\tilde{\alpha}_L}(\partial\phi_L)^2 \;=\; 2(\partial\phi_{\text{canon}})^2 \quad\Longleftrightarrow\quad \phi_{\text{canon}} \;=\; \frac{\phi_L}{\sqrt{2\tilde{\alpha}_L}} \;=\; \sqrt{\frac{4\pi}{\alpha_Y}}\,\phi_L, \qquad \phi_L \;=\; \sqrt{\frac{\alpha_Y}{4\pi}}\,\phi_{\text{canon}},$$

$$S_{\text{can}} \;=\; \int d^4x\,\sqrt{-g}\,\frac{c^4}{16\pi G}\Bigl[R - 2\,g^{\mu\nu}\partial_\mu\phi_{\text{canon}}\partial_\nu\phi_{\text{canon}} - 2\,m_\phi^2\phi_{\text{canon}}^2\Bigr] \;+\; S_m\!\left[A^2 g_{\mu\nu},\,\Psi_m\right], \qquad A \;=\; \exp\!\Bigl(-\alpha_Y\sqrt{2\tilde{\alpha}_L}\,\phi_{\text{canon}}\Bigr),$$

$$\frac{d\ln A}{d\phi_{\text{canon}}} \;=\; \sqrt{2\tilde{\alpha}_L}\,\frac{d\ln A}{d\phi_L} \;=\; -\alpha_Y\sqrt{\frac{\alpha_Y}{4\pi}} \;=\; -\frac{1}{\sqrt{4\pi}} \;=\; -0.2820947918 \quad (\alpha_Y = 1).$$

Equivalently, with the $4\pi$-normalized variable $\hat{\phi}_{\text{canon}} \equiv \phi_{\text{canon}}/\sqrt{4\pi} = \phi_L/\sqrt{\alpha_Y}$ (kinetic term $-\frac{c^4}{2G}(\partial\hat{\phi}_{\text{canon}})^2$), $d\ln A/d\hat{\phi}_{\text{canon}} = -\alpha_Y^{3/2}$. The fifth-force strength is therefore

$$\boxed{\;\beta_{\text{Yuk}} \;=\; \Bigl(\frac{d\ln A}{d\phi_{\text{canon}}}\Bigr)^{2} \;=\; \frac{(d\ln A/d\hat{\phi}_{\text{canon}})^2}{4\pi} \;=\; 2\tilde{\alpha}_L\Bigl(\frac{d\ln A}{d\phi_L}\Bigr)^{2} \;=\; 2\tilde{\alpha}_L\,\alpha_Y^2 \;=\; \frac{\alpha_Y^3}{4\pi} \;=\; \frac{1}{4\pi} \;\approx\; 0.0795774715 \quad (\alpha_Y = 1).\;}$$

Varying $S[\phi_L]$ in the non-relativistic dust / weak-field limit ($\tilde{P} \ll \tilde{\varepsilon} \simeq \rho_m c^2$, $|\phi_L| \ll 1$) gives the linearized static source equation and the Jordan-frame Newtonian potential:

$$\bigl(\nabla^2 - L^{-2}\bigr)\phi_L \;=\; -8\pi\tilde{\alpha}_L\,\alpha_Y\,\frac{G\rho_m}{c^2} \;=\; -\alpha_Y^2\,\frac{G\rho_m}{c^2}, \qquad \tilde{\Psi}(r) \;=\; -\frac{GM}{r}\bigl(1 + \beta_{\text{Yuk}}e^{-r/L}\bigr), \qquad \frac{G_{\text{eff}}(k)}{G} \;=\; 1 + \beta_{\text{Yuk}}\frac{k^2L^2}{1 + k^2L^2}.$$

At $\alpha_Y = 1$ the value $\beta_{\text{Yuk}} = 1/(4\pi)$ is unchanged from V25; for $\alpha_Y \neq 1$ the consistent conformal scaling is $\beta_{\text{Yuk}} = 2\tilde{\alpha}_L\alpha_Y^2 = \alpha_Y^3/(4\pi)$ (`BETA_YUK = 2*ALPHA_TILDE_L*ALPHA_Y**2`), with effective Brans–Dicke parameter $\omega_{\text{BD}} = \frac{1}{2\beta_{\text{Yuk}}} - \frac{3}{2} = \frac{2\pi}{\alpha_Y^3} - \frac{3}{2} = 4.7831853072$ at $\alpha_Y = 1$ (distinguished from the legacy `step5` $4\pi$-convention formula $\frac{1 - 3\alpha_{\text{std}}}{2\alpha_{\text{std}}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2} = 11.0663706144$). Here $\alpha_{\text{DEF}} \equiv d\ln A/d\phi_{\text{canon}} = -1/\sqrt{4\pi} = -0.2820947918$ denotes the first-order Damour–Esposito-Farèse coupling derivative (with second derivative $\beta_{\text{DEF}} \equiv d^2\ln A/d\phi_{\text{canon}}^2 = 0$), while $\beta_{\text{Yuk}} \equiv \alpha_{\text{DEF}}^2 = 1/(4\pi)$ is the dimensionless fifth-force ratio.

### 1.1. Exact Comparison of `Branch II-Conf` (Step 23) vs. `Branch II-Fluid` (Step 22): Interior Equations, Frame Mapping, and Exterior Vacuum Identity

While `Branch II-Fluid` (`step3`, `step4`, `step7`, `step12`, `step22`) and `Branch II-Conf` (`step1`, `step2`, `step5`, `step14`, `step23`–`step25`) share the exact same Einstein-frame gravitational and massive scalar action $S_{g+\phi}[g_{\mu\nu},\phi_L]$, they use two distinct interior matter couplings for a barotropic fluid with equation of state $\tilde{\varepsilon}(h), \tilde{P}(h), \tilde{\rho}_0(h)$ ($dh = d\tilde{P}/(\tilde{\varepsilon} + \tilde{P})$):

1. **`Branch II-Fluid` (`yukawa_action_el` / Step 22 direct fluid-density coupling):**
   - *Matter Lagrangian density:* $\mathcal{L}_m^{\text{(II-Fluid)}} = -e^{-\phi_L}\,\varepsilon(n)$, with $P(n) = n\,d\varepsilon/dn - \varepsilon(n)$ and Einstein-frame fluid stress-energy tensor $T^{\mu\nu}_{\text{matter,(II-Fluid)}} = e^{-\phi_L}[(\varepsilon + P)u^\mu u^\nu + P g^{\mu\nu}]$, so $\rho_m^{(E)} = e^{-\phi_L}\varepsilon(h)$ and $P_m^{(E)} = e^{-\phi_L}P(h)$.
   - *Einstein ($m', \nu'$), scalar ($\phi_L''$), and Euler–TOV ($h'$) equations ($B(r) = 1 - 2m(r)/r$, $G=c=1$):*
     $$\frac{dm}{dr} = 4\pi r^2\bigl(e^{-\phi_L}\varepsilon + \rho_\phi\bigr), \qquad \frac{d\nu}{dr} = \frac{m + 4\pi r^3\bigl(e^{-\phi_L}P + p_{r,\phi}\bigr)}{r(r - 2m)}, \qquad \frac{dh}{dr} = -\nu'(r) + \phi_L'(r),$$
     $$B\,\phi_L'' + \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L' - \frac{\phi_L}{L^2} \;=\; -\frac{\alpha_Y G}{c^2}\,e^{-\phi_L}\,\varepsilon(r), \qquad \rho_\phi = \frac{B(\phi_L')^2 + \phi_L^2/L^2}{2\alpha_Y}, \quad p_{r,\phi} = \frac{B(\phi_L')^2 - \phi_L^2/L^2}{2\alpha_Y}.$$

2. **`Branch II-Conf` (Step 23 universal conformal Jordan coupling):**
   - *Matter action and frame mapping:* $S_m^{\text{(II-Conf)}}[\tilde{g}_{\mu\nu},\Psi_m]$ with $\tilde{g}_{\mu\nu} = A^2(\phi_L)g_{\mu\nu} = e^{-2\alpha_Y\phi_L}g_{\mu\nu}$ and $\sqrt{-\tilde{g}} = A^4(\phi_L)\sqrt{-g} = e^{-4\alpha_Y\phi_L}\sqrt{-g}$. A fluid with proper Jordan-frame energy density $\tilde{\varepsilon}(h)$, pressure $\tilde{P}(h)$, and 4-velocity $\tilde{u}^\mu = A^{-1}u^\mu$ has Jordan stress-energy $\tilde{T}^{\mu\nu} = (\tilde{\varepsilon} + \tilde{P})\tilde{u}^\mu\tilde{u}^\nu + \tilde{P}\tilde{g}^{\mu\nu}$ (trace $\tilde{T} = \tilde{g}_{\mu\nu}\tilde{T}^{\mu\nu} = -\tilde{\varepsilon} + 3\tilde{P}$) and Einstein-frame stress-energy tensor $T^{\mu\nu}_{\text{matter,(II-Conf)}} = \frac{2}{\sqrt{-g}}\frac{\delta S_m}{\delta g_{\mu\nu}} = A^6(\phi_L)\tilde{T}^{\mu\nu}$, i.e., $T^\mu{}_{\nu,\text{matter,(II-Conf)}} = e^{-4\alpha_Y\phi_L}\,\text{diag}(-\tilde{\varepsilon}, \tilde{P}, \tilde{P}, \tilde{P})$, so $\rho_m^{(E)} = e^{-4\alpha_Y\phi_L}\tilde{\varepsilon}(h)$, $P_m^{(E)} = e^{-4\alpha_Y\phi_L}\tilde{P}(h)$, and $T_{\text{matter}}^{(E)} = -e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$.
   - *Einstein ($m', \nu'$), scalar ($\phi_L''$), and Euler–TOV ($h'$) equations:*
     $$\frac{dm}{dr} = 4\pi r^2\bigl(e^{-4\alpha_Y\phi_L}\tilde{\varepsilon} + \rho_\phi\bigr), \qquad \frac{d\nu}{dr} = \frac{m + 4\pi r^3\bigl(e^{-4\alpha_Y\phi_L}\tilde{P} + p_{r,\phi}\bigr)}{r(r - 2m)}, \qquad \frac{dh}{dr} = -\nu'(r) + \alpha_Y\phi_L'(r),$$
     $$B\,\phi_L'' + \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L' - \frac{\phi_L}{L^2} \;=\; -\frac{\alpha_Y^2 G}{c^2}\,e^{-4\alpha_Y\phi_L}\bigl(\tilde{\varepsilon}(r) - 3\tilde{P}(r)\bigr).$$

**Structural consequences for Steps 22–25:**
- **Exterior vacuum ($r > R_*$, where $\tilde{\varepsilon} = \tilde{P} = 0$):** Both branches satisfy the **exact same 4-component curved exterior Einstein–scalar vacuum ODE system** (`step22`), differing only in the surface boundary values $(R_*, M_{\text{surf}}, \phi_L(R_*))$ supplied by the interior integration.
- **Weak-field non-relativistic dust limit ($\tilde{P} \ll \tilde{\varepsilon} \simeq \rho_m c^2$, $|\phi_L| \ll 1$):** At the baseline $\alpha_Y = 1$ ($\alpha_Y^2 = \alpha_Y = 1$), the two linearized sources coincide; for $\alpha_Y \neq 1$, the `Branch II-Conf` source scales as $\alpha_Y^2$ whereas the `Branch II-Fluid` source scales as $\alpha_Y$.
- **Relativistic neutron-star interior ($0 \le r \le R_*$):** Because relativistic pressure reduces the trace ($3\tilde{P}_c/\tilde{\varepsilon}_c = 0.566065$, so $1 - 3\tilde{P}_c/\tilde{\varepsilon}_c = 0.433935$ at the center of the SLy $\log_{10}\rho_{0c}=15.0$ star) and the conformal volume weight is $e^{-4\alpha_Y\phi_L}$ rather than $e^{-\phi_L}$, the two interior BVPs are **distinct**. Accordingly, in Section 3.3 below, the ray-tracing table on the V25 Step 22 background (`Table 2`, $M_{\text{ADM}} = 1.50999814\;M_\odot$) is explicitly labeled as a **conditional cross-branch diagnostic evaluated on the `Branch II-Fluid` background**, and is accompanied by the **self-consistent `Branch II-Conf` trace-sourced interior + exterior BVP solution** (`Table 2B`, $M_{\text{ADM}}^{\text{(II-Conf)}} = 1.56203491\;M_\odot$).

At the fixed baseline parameters $L = 10\;\text{km}$, $\alpha_Y = 1$, the numerical values of the coupling dictionary are:

| Quantity | Formula | Numerical Value ($L = 10\;\text{km}$, $\alpha_Y = 1$) | Unit |
|---|---|---|---|
| Yukawa screening length ($L$) | $L = 1/m_\phi$ | `10.0` (`1.0e4`) | $\text{km}$ ($\text{m}$) |
| Yukawa inverse length ($m_\phi$) | $m_\phi = 1/L$ | `0.1` (`1.0e-4`) | $\text{km}^{-1}$ ($\text{m}^{-1}$) |
| Scalar Compton mass ($m_\phi \hbar c$) | $\hbar c / L$ | `1.9732698e-11` | $\text{eV}$ |
| Scalar cutoff angular frequency ($\omega_{\text{cut}}$) | $c / L$ | `29979.2458` | $\text{rad}\,\text{s}^{-1}$ |
| Scalar cutoff linear frequency ($f_{\text{cut}}$) | $c / (2\pi L)$ | `4771.345159` | $\text{Hz}$ |
| Yukawa source strength ($\alpha_Y$) | $4\pi\alpha_{\text{std}} = 8\pi\tilde{\alpha}_L$ | `1.0` | dimensionless |
| Standard fifth-force coupling ($\alpha_{\text{std}}$) | $\alpha_Y / (4\pi)$ | `0.0795774715` | dimensionless |
| Action normalization constant ($\tilde{\alpha}_L$) | $\alpha_Y / (8\pi)$ | `0.0397887358` | dimensionless |
| Canonical rescaling ($\phi_{\text{canon}}/\phi_L$) | $1/\sqrt{2\tilde{\alpha}_L} = \sqrt{4\pi/\alpha_Y}$ | `3.5449077018` | dimensionless |
| First-order DEF derivative ($\lvert\alpha_{\text{DEF}}\rvert = \lvert d\ln A/d\phi_{\text{canon}}\rvert$; $\beta_{\text{DEF}} = d^2\ln A/d\phi_{\text{canon}}^2 = 0$) | $\alpha_Y\sqrt{2\tilde{\alpha}_L} = 1/\sqrt{4\pi}$ | `0.2820947918` | dimensionless |
| Fifth-force strength ($\beta_{\text{Yuk}}$) | $(d\ln A/d\phi_{\text{canon}})^2 = (d\ln A/d\hat{\phi}_{\text{canon}})^2/(4\pi) = 2\tilde{\alpha}_L\alpha_Y^2$ | `0.0795774715` | dimensionless |
| Effective Brans–Dicke parameter ($\omega_{\text{BD}}$) | $1/(2\beta_{\text{Yuk}}) - 3/2 = 2\pi/\alpha_Y^3 - 3/2$ | `4.7831853072` | dimensionless |

---

## 2. Part A: GW170817 / GRB 170817A and Tensor Gravitational-Wave Speed ($c_T = 1$)

### 2.1. Covariant Perturbation Proof in Both Frames (Einstein & Jordan / Horndeski Mapping)

1. **Einstein-Frame Tensor Modes ($g_{\mu\nu} = \bar{g}_{\mu\nu} + h_{\mu\nu}$):**
   In the Einstein frame, the gravitational kinetic term is purely the Einstein–Hilbert term $\frac{c^4}{16\pi G}R$; the scalar field $\phi_L$ has no derivative couplings to the Riemann, Ricci, or Einstein tensors. In transverse-traceless (TT) gauge ($h^\mu{}_\mu = 0$, $\bar{\nabla}^\mu h_{\mu\nu}^{\text{TT}} = 0$), the quadratic tensor action reads:

   $$S_T^{(2)} \;=\; \frac{c^4}{64\pi G}\int d^4x\,\sqrt{-\bar{g}}\left[-\bar{g}^{\mu\nu}\partial_\mu h_{ij}^{\text{TT}}\,\partial_\nu h^{ij\,\text{TT}}\right] \quad\Longrightarrow\quad \Box_{\bar{g}}\,h_{ij}^{\text{TT}} \;=\; 0.$$

   Hence the propagation speed of tensor gravitational waves in the Einstein frame is identically:

   $$c_{T,\text{Einstein}}^2 \;=\; c^2 \quad\Longleftrightarrow\quad \frac{c_{T,\text{Einstein}} - c}{c} \;=\; 0.$$

2. **Jordan-Frame Horndeski Functions ($\tilde{g}_{\mu\nu} = A^2(\phi_L) g_{\mu\nu}$):**
   Transforming to the physical Jordan-frame metric $\tilde{g}_{\mu\nu} = A^2(\phi_L) g_{\mu\nu}$ (in the locked V26 convention $A(\phi_L) = \exp(-\alpha_Y\phi_L)$, $\phi_L > 0$), the scalar–tensor action maps into the Horndeski class with curvature coupling functions:

   $$\tilde{G}_4(\phi_L) \;=\; \frac{c^4}{16\pi G}\,A^{-2}(\phi_L), \qquad \tilde{G}_{4,X} \;\equiv\; \frac{\partial\tilde{G}_4}{\partial X} \;=\; 0, \qquad \tilde{G}_5(\phi_L, X) \;=\; 0,$$

   where $X \equiv -\frac{1}{2}\tilde{g}^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L$. Because $\tilde{G}_4$ is a purely algebraic function of $\phi_L$ and independent of the kinetic invariant $X$ ($\tilde{G}_{4,X} = 0$, $\tilde{G}_5 = 0$), the general Horndeski tensor speed formula yields:

   $$\boxed{\;c_{T,\text{Jordan}}^2 \;=\; \frac{\tilde{G}_4}{\tilde{G}_4 - 2 X \tilde{G}_{4,X}}\,c^2 \;=\; c^2 \quad\Longrightarrow\quad \frac{c_T - c}{c} \;=\; 0 \quad \text{identically for all }(L,\,\alpha_Y)\text{!}\;}$$

   **Scope (V26):** $c_T = c$ is proven on the purely conformal, minimal Einstein–Hilbert + canonical scalar + universal Maxwell coupling branch ($G_{4,X} = 0$, $G_5 = 0$); it does not automatically prove $c_T = c$ on a disformal/refractive branch ($D(\phi_L, X) \ne 0$).

3. **Conformal Invariance of Null Cones:**
   The conformal transformation $\tilde{g}_{\mu\nu} = A^2(\phi_L) g_{\mu\nu}$ preserves null cones identically ($\tilde{g}_{\mu\nu} k^\mu k^\nu = A^2(\phi_L) g_{\mu\nu} k^\mu k^\nu = 0 \Longleftrightarrow g_{\mu\nu} k^\mu k^\nu = 0$). Because four-dimensional vacuum Maxwell electrodynamics ($\sqrt{-\tilde{g}}\,\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} = \sqrt{-g}\,g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}$) is conformally invariant, electromagnetic waves (photons) and tensor gravitational waves $h_{+,\times}^{\text{TT}}$ propagate along the exact same null cone.
   - **Comparison with the GW170817 / GRB 170817A Observational Bound** ($\Delta t = 1.74 \pm 0.05\;\text{s}$ over $D = 40\;\text{Mpc} = 1.2343\times 10^{20}\;\text{km}$, light travel time $t_{\text{light}} = 4.1171\times 10^{15}\;\text{s} \approx 1.3046\times 10^8\;\text{yr}$, Abbott et al. 2017):
     $$-3\times 10^{-15} \;\le\; \frac{c_T - c}{c} \;\le\; +7\times 10^{-15}, \qquad \left(\frac{c_T - c}{c}\right)_{\text{TAL}} \;=\; 0.$$
     Both symbolic (`SymPy`) and numerical checks (`gw170817_and_lensing_summary.json`) confirm a residual of `0` in both frames.

---

### 2.2. Massive Scalar Breathing Mode Dispersion and the $4.77\;\text{kHz}$ Cutoff

Linearizing the vacuum scalar wave equation ($\Box\phi_L - m_\phi^2\phi_L = 0$, $m_\phi = 1/L$) for plane waves ($\delta\phi_L \propto e^{i(\mathbf{k}\cdot\mathbf{x} - \omega t)}$) yields the Klein–Gordon dispersion relation:

$$\omega^2 \;=\; c^2 k^2 + \omega_{\text{cut}}^2, \qquad \omega_{\text{cut}} \;=\; \frac{c}{L} \;=\; 29\,979.2458\;\text{rad}\,\text{s}^{-1}, \qquad f_{\text{cut}} \;=\; \frac{c}{2\pi L} \;=\; 4771.345159\;\text{Hz}.$$

1. **Sub-Cutoff (Evanescent / Yukawa-Screened) Regime ($f < f_{\text{cut}} = 4771.345\;\text{Hz}$):**
   At all frequencies below the $4.771\;\text{kHz}$ Compton cutoff — including binary-pulsar orbital frequencies ($f \sim 10^{-4}\;\text{Hz}$), the LISA band ($10^{-2}\;\text{Hz}$), the entire GW170817 LIGO inspiral band ($f_{\text{GW}} \in [24, 600]\;\text{Hz}$), and even the post-merger neutron-star fundamental quadrupole $f$-mode ($f \sim 1500\text{--}3000\;\text{Hz}$) — the wavenumber is purely imaginary:

   $$k \;=\; i\,\kappa_{\text{ev}}(f), \qquad \kappa_{\text{ev}}(f) \;=\; \frac{1}{L}\sqrt{1 - \left(\frac{f}{f_{\text{cut}}}\right)^2}, \qquad \lambda_{\text{ev}}(f) \;=\; \frac{1}{\kappa_{\text{ev}}(f)} \;=\; \frac{L}{\sqrt{1 - (f/f_{\text{cut}})^2}} \;\approx\; 10\;\text{km}.$$

   The scalar breathing mode does not propagate as a wave; it is an evanescent near-zone field decaying exponentially over $\lambda_{\text{ev}}(f) \approx 10\text{--}12.86\;\text{km}$. Over the GW170817 propagation distance $D = 40\;\text{Mpc} = 1.2343\times 10^{20}\;\text{km}$, the scalar amplitude suppression is:

   $$\exp\!\left(-\frac{D}{\lambda_{\text{ev}}(f)}\right) \;\sim\; 10^{-5.36\times 10^{19}} \;=\; 0.$$

2. **Super-Cutoff (Propagating Dispersive) Regime ($f > f_{\text{cut}} = 4771.345\;\text{Hz}$):**
   If a high-frequency post-merger burst excited modes with $f > 4771.345\;\text{Hz}$, their group and phase velocities would be:

   $$v_g(f) \;=\; \frac{\partial\omega}{\partial k} \;=\; c\sqrt{1 - \left(\frac{f_{\text{cut}}}{f}\right)^2} \;<\; c, \qquad v_p(f) \;=\; \frac{\omega}{k} \;=\; \frac{c^2}{v_g(f)} \;>\; c.$$

   Over $D = 40\;\text{Mpc}$, the arrival time delay of a super-cutoff scalar mode relative to the $c_T = c$ tensor wave is:
   - At $f = 5000\;\text{Hz}$ ($v_g/c = 0.298949$): $\Delta t_{\text{scalar}} = 9.655\times 10^{15}\;\text{s} = 3.059\times 10^8\;\text{yr}$ ($305.9$ million years),
   - At $f = 10\,000\;\text{Hz}$ ($v_g/c = 0.878830$): $\Delta t_{\text{scalar}} = 5.676\times 10^{14}\;\text{s} = 1.799\times 10^7\;\text{yr}$ ($18.0$ million years).
   Consequently, terrestrial detectors at $D = 40\;\text{Mpc}$ observe purely the $c_T = c$ tensor polarizations $h_+, h_\times$ with zero scalar contamination.

#### Table 1: GW170817 Tensor Wave Speed and Massive Scalar Mode Dispersion ($L = 10\;\text{km}$, $\alpha_Y = 1$; `gw170817_dispersion_and_cT.csv`)

| Mode / Band (`mode_label`) | $f$ ($\text{Hz}$) | $f / f_{\text{cut}}$ | Scalar Regime | $c_T / c$ (Einstein & Jordan) | $(c_T - c)/c$ | $k_{\text{real}}$ ($\text{km}^{-1}$) | $\kappa_{\text{ev}}$ ($\text{km}^{-1}$) | $\lambda_{\text{ev}}$ ($\text{km}$) | $\log_{10}\!\left(A_{\text{scalar}}(40\;\text{Mpc})\right)$ | $v_g / c$ | $\Delta t_{\text{scalar}}(40\;\text{Mpc})$ ($\text{yr}$) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `PSR_B1913p16_orbital` | `3.583e-5` | `7.509e-9` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.100000` | `10.000000` | `-5.360e+19` | `0.0` | `0.0` |
| `PSR_J0737m3039_orbital` | `1.132e-4` | `2.372e-8` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.100000` | `10.000000` | `-5.360e+19` | `0.0` | `0.0` |
| `LISA_mHz_band` | `1.000e-2` | `2.096e-6` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.100000` | `10.000000` | `-5.360e+19` | `0.0` | `0.0` |
| `GW170817_LIGO_entry_24Hz` | `24.0` | `0.005030` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.099999` | `10.000127` | `-5.360e+19` | `0.0` | `0.0` |
| `GW170817_LIGO_100Hz` | `100.0` | `0.020958` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.099978` | `10.002197` | `-5.359e+19` | `0.0` | `0.0` |
| `GW170817_LIGO_300Hz` | `300.0` | `0.062875` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.099802` | `10.019825` | `-5.350e+19` | `0.0` | `0.0` |
| `GW170817_LIGO_600Hz` | `600.0` | `0.125751` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.099206` | `10.080016` | `-5.318e+19` | `0.0` | `0.0` |
| `NS_post_merger_1500Hz` | `1500.0` | `0.314377` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.094930` | `10.534097` | `-5.089e+19` | `0.0` | `0.0` |
| `NS_f_mode_2500Hz` | `2500.0` | `0.523961` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.085174` | `11.740643` | `-4.566e+19` | `0.0` | `0.0` |
| `NS_high_3000Hz` | `3000.0` | `0.628754` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.077760` | `12.860005` | `-4.168e+19` | `0.0` | `0.0` |
| `Sub_cutoff_4500Hz` | `4500.0` | `0.943130` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.033242` | `30.082119` | `-1.782e+19` | `0.0` | `0.0` |
| `Yukawa_cutoff_f_cut` | `4771.345159` | `1.000000` | `cutoff_threshold` | `1.0` | `0.0` | `0.0` | `0.000000` | `inf` | `0.0` | `0.0` | `inf` |
| `Super_cutoff_4800Hz` | `4800.0` | `1.006006` | `propagating_dispersive` | `1.0` | `0.0` | `0.010976` | `0.0` | `0.0` | `0.0` | `0.109105` | `1.065291e+09` |
| `Super_cutoff_5000Hz` | `5000.0` | `1.047923` | `propagating_dispersive` | `1.0` | `0.0` | `0.031328` | `0.0` | `0.0` | `0.0` | `0.298949` | `3.059412e+08` |
| `Super_cutoff_7500Hz` | `7500.0` | `1.571884` | `propagating_dispersive` | `1.0` | `0.0` | `0.121277` | `0.0` | `0.0` | `0.0` | `0.771541` | `3.863088e+07` |
| `Super_cutoff_10000Hz` | `10000.0` | `2.095845` | `propagating_dispersive` | `1.0` | `0.0` | `0.184189` | `0.0` | `0.0` | `0.0` | `0.878830` | `1.798767e+07` |

---

## 3. Part B: Gravitational Lensing ($\Phi + \Psi$, $\gamma_{\text{PPN}}(r)$, Curved-Spacetime $\hat{\alpha}(b)$, and Convergence $\kappa$) at $L = 10\;\text{km}$, $\alpha_Y = 1$

### 3.1. Weak-Field Potentials in the Einstein and Jordan Frames: Conformal Cancellation in the Weyl Potential

In Newtonian gauge, the weak-field Einstein-frame metric reads:

$$ds^2 \;=\; -\left(1 + \frac{2\Psi_E(r)}{c^2}\right)c^2 dt^2 \;+\; \left(1 - \frac{2\Phi_E(r)}{c^2}\right)d\mathbf{x}^2, \qquad \Psi_E(r) \;=\; \Phi_E(r) \;=\; -\frac{G M_{\text{ADM}}}{r}.$$

In the physical Jordan-frame metric ($\tilde{g}_{\mu\nu} = A^2(\phi_L) g_{\mu\nu} = \exp(-2\alpha_Y\phi_L) g_{\mu\nu}$, $A(\phi_L) = \exp(-\alpha_Y\phi_L)$, $\phi_L > 0$), the temporal ($\tilde{\Psi}$) and spatial ($\tilde{\Phi}$) metric potentials are shifted at linear order with opposite signs:

$$\tilde{\Psi}(r) \;=\; \Psi_E(r) - \alpha_Y c^2 \phi_L(r), \qquad \tilde{\Phi}(r) \;=\; \Phi_E(r) + \alpha_Y c^2 \phi_L(r).$$

Consequently, **in the Weyl lensing potential governing null geodesics (gravitational lensing and Shapiro time delay at fixed $M_{\text{ADM}}$), the conformal factor cancels identically at linear order**:

$$\boxed{\;\Phi_{\text{lens}}(r) \;\equiv\; \frac{1}{2}\bigl(\tilde{\Psi}(r) + \tilde{\Phi}(r)\bigr) \;=\; \frac{1}{2}\bigl(\Psi_E(r) + \Phi_E(r)\bigr) \;=\; -\frac{G M_{\text{ADM}}}{r}\quad (\text{outside the Yukawa energy cloud}).\;}$$

---

### 3.2. Scale-Dependent PPN Slip ($\gamma_{\text{PPN}}(r)$) and the Lensing-to-Dynamical Mass Ratio

Whereas relativistic photons respond to the Weyl combination $\Phi_{\text{lens}}(r) = \frac{1}{2}(\tilde{\Psi} + \tilde{\Phi}) = -G M_{\text{ADM}}/r$, non-relativistic test bodies (planets, stars, gas clouds) respond to the dynamical Newtonian potential:

$$\tilde{\Psi}(r) \;=\; -\frac{G M_{\text{ADM}}}{r}\bigl(1 + \beta_{\text{Yuk}}\,e^{-r/L}\bigr), \qquad \tilde{\Phi}(r) \;=\; -\frac{G M_{\text{ADM}}}{r}\bigl(1 - \beta_{\text{Yuk}}\,e^{-r/L}\bigr).$$

This yields the scale-dependent Eddington–Robertson PPN slip parameter $\gamma_{\text{PPN}}(r) \equiv \tilde{\Phi}(r)/\tilde{\Psi}(r)$ and the lensing-to-dynamical potential ratio:

$$\boxed{\;\gamma_{\text{PPN}}(r) \;=\; \frac{1 - \beta_{\text{Yuk}}\,e^{-r/L}}{1 + \beta_{\text{Yuk}}\,e^{-r/L}}, \qquad \frac{\Phi_{\text{lens}}(r)}{\tilde{\Psi}(r)} \;=\; \frac{1 + \gamma_{\text{PPN}}(r)}{2} \;=\; \frac{1}{1 + \beta_{\text{Yuk}}\,e^{-r/L}}.\;}$$

The V26 package uses a single normalization (`lensing_deflection_and_ppn_slip.csv`): $\beta_{\text{Yuk}}(\alpha_Y) = 2\tilde{\alpha}_L\alpha_Y^2 = \alpha_Y^3/(4\pi) = 1/(4\pi) \approx 0.0795774715$ at $\alpha_Y = 1$ ($\tilde{\alpha}_L = \frac{1}{8\pi}$). No parallel $\beta_{\text{Yuk}} = 2.0$ column is reported.

At the Cassini 2003 solar conjunction impact parameter ($b = 1.6\,R_\odot = 1.11312\times 10^6\;\text{km}$, $b/L = 1.11312\times 10^5$), the exponential Yukawa factor is $e^{-b/L} = 10^{-48\,342.19} = 0$, giving:

$$\gamma_{\text{PPN}}(1.6\,R_\odot) \;=\; 1, \qquad \bigl|1 - \gamma_{\text{PPN}}(1.6\,R_\odot)\bigr| \;<\; 10^{-48\,342} \;\ll\; 2.3\times 10^{-5}\;\text{(Bertotti et al. 2003)}.$$

---

### 3.3. Exact Curved-Spacetime Deflection Angle ($\hat{\alpha}(b)$) and Shapiro Time Delay: Conditional Diagnostic on the V25 Step 22 (`Branch II-Fluid`) Background vs. Self-Consistent `Branch II-Conf` Trace-Sourced Solution

As established in Section 1.1, the V25 Step 22 neutron-star background (`curved_yukawa_adm_and_stability.py`) integrates the interior with the `Branch II-Fluid` coupling $\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$ (source $-\alpha_Y e^{-\phi_L}\varepsilon$), whereas `Branch II-Conf` sources $\phi_L$ inside a relativistic star via the trace $-\alpha_Y^2 e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$. Evaluating null-geodesic light deflection on the V25 Step 22 metric is therefore a **conditional cross-branch diagnostic on the `Branch II-Fluid` background**, which we present in **Table 2** alongside the **self-consistent `Branch II-Conf` trace-sourced interior + exterior BVP solution** in **Table 2B** (`gw170817_and_lensing_summary.json`, key `branch_ii_fluid_vs_branch_ii_conf_interior_bvp_comparison`).

For the V25 (`step22`, DOI: `10.5281/zenodo.23699685`, `Branch II-Fluid` background) reference SLy neutron star ($\log_{10}(\rho_{0c}/[\text{g cm}^{-3}]) = 15.0$, $L = 10\;\text{km}$, $\alpha_Y = 1$), the exact 4-component exterior BVP parameters are:
- Stellar radius: $R_* = 11.46572423\;\text{km}$,
- Coordinate-radius surface mass: $M_{\text{surf}} = m(R_*) = 1.50766122\;M_\odot$ ($2.22625027\;\text{km}$),
- Asymptotically flat ADM mass: $M_{\text{ADM}} = \lim_{r\to\infty} m(r) = 1.50999814\;M_\odot$ ($2.22970113\;\text{km}$),
- Exterior scalar cloud mass contribution: $\Delta M_{\phi,\text{ext}}^{\text{curved}} = +0.00233692\;M_\odot$ ($+0.155003\%$ relative to $M_{\text{surf}}$),
- Surface scalar field and curved derivative: $\phi_L(R_*) = 0.0054086823$, $\phi_{L,\text{curved}}'(R_*) = -0.0012631563\;\text{km}^{-1}$ ($k_{\text{curved}} = 0.23354233\;\text{km}^{-1}$),
- Surface metric potential: $\nu_{\text{exact}}(R_*) = -0.2463066872$ (compared with the vacuum Schwarzschild values $\nu_{\text{Schw}}(M_{\text{surf}}) = -0.2457823129$ and $\nu_{\text{Schw}}(M_{\text{ADM}}) = -0.2462745914$).

In the spherically symmetric exterior spacetime ($ds^2 = -e^{2\nu(r)}c^2 dt^2 + B(r)^{-1}dr^2 + r^2 d\Omega^2$, $B(r) = 1 - 2m(r)/r$), the exact curved-spacetime deflection angle for a light ray with impact parameter $b = r_0\,e^{-\nu_{\text{eff}}(r_0)}$ (evaluated via cancellation-free 128-point Gauss–Legendre quadrature with $u = r_0/r = \sin\theta$, $\theta \in [0, \pi/2]$) is:

$$\hat{\alpha}(b) \;=\; 2\int_0^{\pi/2}\left[\frac{\cos\theta}{\sqrt{B(r_0/\sin\theta)}\sqrt{\exp\!\bigl(2[\nu_{\text{eff}}(r_0) - \nu_{\text{eff}}(r_0/\sin\theta)]\bigr) - \sin^2\theta}} \;-\; 1\right]d\theta,$$

where on the conformal Jordan branch $\nu_{\text{eff}}(r) = \nu(r)$ (since $d\tilde{s}^2 = e^{-2\alpha_Y\phi_L}ds^2 = 0 \Longleftrightarrow ds^2 = 0$), while on the disformal refractive branch (`step2`, `step14`, $n_{\text{disf}}(r) = e^{\phi_L(r)}$) $\nu_{\text{eff}}(r) = \nu(r) - \phi_L(r)$.

**Three Physical Consequences from the Stellar Surface ($b = 14.749\;\text{km}$) to Galaxy Clusters ($b = 1\;\text{Mpc}$):**
1. **Birkhoff-Violating Effect of the Exterior Yukawa Energy Cloud near the Neutron Star ($b \in [14.75, 100]\;\text{km}$):**
   - For the surface-grazing ray ($b = 14.749032\;\text{km}$, closest approach $r_0 = 11.558475\;\text{km}$), the enclosed mass is $m(r_0) = 1.50772597\;M_\odot < M_{\text{ADM}}$.
   - The conformal Yukawa deflection angle is $\hat{\alpha}_{\text{Yukawa,conf}} = 1.29085207\;\text{rad}$ ($73.96037^\circ$), which is **$+0.36521\%$ larger** than the vacuum Schwarzschild deflection with the surface mass ($\hat{\alpha}_{\text{Schw}}(M_{\text{surf}}) = 1.28615489\;\text{rad}$), but **$-0.02867\%$ smaller** than the vacuum Schwarzschild deflection with the asymptotic ADM mass ($\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}}) = 1.29122231\;\text{rad}$) because the light ray dives inside the exterior scalar energy halo!
   - On the disformal refractive branch ($n_{\text{disf}} = e^{\phi_L}$), the grazing deflection angle is $\hat{\alpha}_{\text{Yukawa,disf}} = 1.34618158\;\text{rad}$ ($+4.28628\%$ above the conformal branch), decaying exponentially to $+5.253\times 10^{-6}$ at $b = 100\;\text{km}$ and $0$ at $b \ge 1000\;\text{km}$.
2. **Asymptotic ADM Mass Calibration ($b \ge 200\;\text{km}$):**
   - For $b \ge 200\;\text{km} = 20\,L$, the scalar field vanishes ($\phi_L(200\;\text{km}) = 1.069\times 10^{-12}$), so $\hat{\alpha}_{\text{Yukawa,conf}}(b)$ coincides with $\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}}, b)$ to machine precision ($\le 10^{-12}$ relative difference, and $0.0$ for $b \ge R_\odot$), while maintaining a permanent **$+0.155003\% = \Delta M_{\phi,\text{ext}}^{\text{curved}}/M_{\text{surf}}$** excess over the naive surface-mass Schwarzschild curve $\hat{\alpha}_{\text{Schw}}(M_{\text{surf}}, b)$ at all asymptotic impact parameters.

#### Table 2: Scale-Dependent PPN Slip $\gamma_{\text{PPN}}(b)$, Curved-Spacetime Light Deflection $\hat{\alpha}(b)$, and Shapiro Time Delay from $b = 14.75\;\text{km}$ to $1\;\text{Mpc}$ (`lensing_deflection_and_ppn_slip.csv`)

| Scale (`scale_label`) | $b$ ($\text{km}$) | $r_0$ ($\text{km}$) | $b / L$ | $\log_{10}(e^{-b/L})$ | $m(r_0)$ ($M_\odot$) | $\gamma_{\text{PPN}}^{\text{TAL}}(b)$ | $\lvert 1 - \gamma_{\text{PPN}}^{\text{TAL}}\rvert$ | $\hat{\alpha}_{\text{Yuk,conf}}$ ($\text{rad}$) | $\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})$ ($\text{rad}$) | $\hat{\alpha}_{\text{Schw}}(M_{\text{surf}})$ ($\text{rad}$) | $\frac{\hat{\alpha}_{\text{Yuk,conf}} - \hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})}{\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})}$ | $\frac{\hat{\alpha}_{\text{Yuk,conf}} - \hat{\alpha}_{\text{Schw}}(M_{\text{surf}})}{\hat{\alpha}_{\text{Schw}}(M_{\text{surf}})}$ | $\Delta t_{\text{Shapiro}}^{\text{Yuk,conf}}$ ($\mu\text{s}$) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `NS_surface_grazing` | `14.749032` | `11.558475` | `1.4749` | `-0.6405` | `1.507726` | `0.964236` | `3.576e-2` | `1.29085207` | `1.29122231` | `1.28615489` | `-2.867e-4` | `+3.652e-3` | `420.139048` |
| `NS_near_halo_15km` | `15.000000` | `11.843795` | `1.5000` | `-0.6514` | `1.507913` | `0.965107` | `3.489e-2` | `1.23768742` | `1.23800029` | `1.23332361` | `-2.527e-4` | `+3.538e-3` | `417.501587` |
| `NS_Yukawa_cloud_20km` | `20.000000` | `17.215854` | `2.0000` | `-0.8686` | `1.509517` | `0.978690` | `2.131e-2` | `0.70230838` | `0.70233272` | `0.70052511` | `-3.465e-5` | `+2.546e-3` | `387.157454` |
| `NS_Yukawa_cloud_30km` | `30.000000` | `27.455794` | `3.0000` | `-1.3029` | `1.509957` | `0.992107` | `7.893e-3` | `0.38749878` | `0.38749944` | `0.38670428` | `-1.700e-6` | `+2.055e-3` | `361.964526` |
| `NS_Yukawa_cloud_50km` | `50.000000` | `47.600307` | `5.0000` | `-2.1715` | `1.509998` | `0.998928` | `1.072e-3` | `0.20641763` | `0.20641763` | `0.20604614` | `-1.039e-8` | `+1.803e-3` | `339.066547` |
| `NS_Yukawa_tail_100km` | `100.000000` | `97.690938` | `10.000` | `-4.3429` | `1.509998` | `0.999993` | `7.226e-6` | `0.09556441` | `0.09556441` | `0.09540578` | `-2.099e-12` | `+1.663e-3` | `313.665881` |
| `NS_Yukawa_tail_200km` | `200.000000` | `197.731861` | `20.000` | `-8.6859` | `1.509998` | `1.000000` | `3.280e-10` | `0.04612015` | `0.04612015` | `0.04604631` | `+8.027e-13` | `+1.603e-3` | `290.894079` |
| `NS_far_1000km` | `1000.00000` | `997.762797` | `100.00` | `-43.429` | `1.509998` | `1.000000` | `5.921e-45` | `0.00897785` | `0.00897785` | `0.00896386` | `+2.771e-12` | `+1.560e-3` | `241.379351` |
| `Sirius_B_radius_5844km` | `5844.00000` | `5841.76902` | `584.40` | `-253.80` | `1.509998` | `1.000000` | `2.513e-255` | `0.00152786` | `0.00152786` | `0.00152550` | `+1.062e-12` | `+1.552e-3` | `188.460119` |
| `NS_far_field_b_6p957e5_km` | `6.95700e+5` | `6.95698e+5` | `6.957e+4` | `-3.021e+4` | `1.509998` | `1.000000` | `0.0` | `1.282002e-5` | `1.282002e-5` | `1.280018e-5` | `0.0` | `+1.550e-3` | `102.503515` |
| `NS_far_field_b_1p113e6_km` | `1.11312e+6` | `1.11312e+6` | `1.113e+5` | `-4.834e+4` | `1.509998` | `1.000000` | `0.0` | `8.012485e-6` | `8.012485e-6` | `8.000084e-6` | `0.0` | `+1.550e-3` | `102.503351` |
| `NS_far_field_b_1AU` | `1.49598e+8` | `1.49598e+8` | `1.496e+7` | `-6.497e+6` | `1.509998` | `1.000000` | `0.0` | `5.961853e-8` | `5.961853e-8` | `5.952626e-8` | `0.0` | `+1.550e-3` | `102.503080` |
| `Outer_Solar_10_AU` | `1.49598e+9` | `1.49598e+9` | `1.496e+8` | `-6.497e+7` | `1.509998` | `1.000000` | `0.0` | `5.961852e-9` | `5.961852e-9` | `5.952626e-9` | `0.0` | `+1.550e-3` | `102.503078` |
| `Subparsec_1_pc` | `3.08568e+13` | `3.08568e+13` | `3.086e+12` | `-1.340e+12` | `1.509998` | `1.000000` | `0.0` | `2.890388e-13` | `2.890388e-13` | `2.885914e-13` | `0.0` | `+1.550e-3` | `102.503078` |
| `Galactic_1_kpc` | `3.08568e+16` | `3.08568e+16` | `3.086e+15` | `-1.340e+15` | `1.509998` | `1.000000` | `0.0` | `2.890388e-16` | `2.890388e-16` | `2.885914e-16` | `0.0` | `+1.550e-3` | `102.503078` |
| `NGC3198_30p08_kpc` | `9.28172e+17` | `9.28172e+17` | `9.282e+16` | `-4.031e+16` | `1.509998` | `1.000000` | `0.0` | `9.609001e-18` | `9.609001e-18` | `9.594130e-18` | `0.0` | `+1.550e-3` | `102.503078` |
| `Bullet_Cluster_250_kpc` | `7.71419e+18` | `7.71419e+18` | `7.714e+17` | `-3.350e+17` | `1.509998` | `1.000000` | `0.0` | `1.156155e-18` | `1.156155e-18` | `1.154366e-18` | `0.0` | `+1.550e-3` | `102.503078` |
| `Cluster_virial_1_Mpc` | `3.08568e+19` | `3.08568e+19` | `3.086e+18` | `-1.340e+18` | `1.509998` | `1.000000` | `0.0` | `2.890388e-19` | `2.890388e-19` | `2.885914e-19` | `0.0` | `+1.550e-3` | `102.503078` |

**Note (V26):** every row of Table 2 is evaluated on the $1.50999814\;M_\odot$ Step 22 `Branch II-Fluid` neutron-star background (`source_body = NS_SLy_1p51Msun`); the `NS_far_field_b_*` rows are far-field benchmark impact parameters only, not solar data. Table 2B below compares this `Branch II-Fluid` diagnostic background directly against the self-consistent `Branch II-Conf` trace-sourced ($\tilde{\varepsilon} - 3\tilde{P}$) BVP solution for the same SLy $\log_{10}\rho_{0c} = 15.0$, $L = 10\;\text{km}$, $\alpha_Y = 1$ star, and Table 3 gives the dedicated $M = 1.000000\;M_\odot$ Solar–Cassini benchmark.

#### Table 2B: Direct Comparison of the Step 22 `Branch II-Fluid` Background vs. the Self-Consistent Step 23 `Branch II-Conf` Trace-Sourced Interior + Curved-Exterior BVP (SLy, $\log_{10}\rho_{0c} = 15.0$, $L = 10\;\text{km}$, $\alpha_Y = 1$)

| Quantity / Observable | `Branch II-Fluid` (Step 22 Background Diagnostic) | `Branch II-Conf` (Self-Consistent Trace-Sourced BVP) | Physical Origin of Difference |
|---|---|---|---|
| Interior matter action | $\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$ | $S_m[e^{-2\alpha_Y\phi_L}g_{\mu\nu},\Psi_m]$ ($\sqrt{-\tilde{g}} = e^{-4\alpha_Y\phi_L}\sqrt{-g}$) | Direct energy-density weight vs. universal Jordan metric |
| Interior scalar source RHS | $-(\alpha_Y G/c^2)\,e^{-\phi_L}\varepsilon(r)$ | $-(\alpha_Y^2 G/c^2)\,e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$ | Relativistic pressure trace reduction ($1 - 3\tilde{P}_c/\tilde{\varepsilon}_c = 0.433935$) |
| Interior Euler–TOV balance | $h'(r) = -\nu'(r) + \phi_L'(r)$ | $h'(r) = -\nu'(r) + \alpha_Y\phi_L'(r)$ | Identical form at $\alpha_Y = 1$, evaluated on respective $\phi_L(r)$ |
| Exterior ($r > R_*$) vacuum ODEs | 4-component curved Einstein–scalar BVP | Identical 4-component curved Einstein–scalar BVP | Vacuum equations coincide identically ($\tilde{\varepsilon} = \tilde{P} = 0$) |
| Stellar radius $R_*$ ($\text{km}$) | `11.46572423` | `11.66269842` | Smaller $e^{-4\phi_L}$ weight and trace source shift equilibrium radius |
| Surface mass $M_{\text{surf}} = m(R_*)$ ($M_\odot$) | `1.50766122` | `1.56063550` | Interior mass integral with $e^{-4\alpha_Y\phi_L}\tilde{\varepsilon}$ |
| Asymptotic ADM mass $M_{\text{ADM}}$ ($M_\odot$) | `1.50999814` | `1.56203491` | $M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$ |
| Exterior scalar cloud mass $\Delta M_{\phi,\text{ext}}^{\text{curved}}$ ($M_\odot$) | `+0.00233692` (`+0.155003%`) | `+0.00139941` (`+0.089669%`) | Lower surface scalar field due to $-3\tilde{P}$ in trace source |
| Central / surface scalar field $\phi_L(0)$ / $\phi_L(R_*)$ | `0.01777896` / `0.00540868` | `0.01157688` / `0.00414318` | Trace factor $1 - 3\tilde{P}/\tilde{\varepsilon} < 1$ inside core |
| Curved surface slope $\phi_{L,\text{curved}}'(R_*)$ ($\text{km}^{-1}$) / $k_{\text{curved}}$ ($\text{km}^{-1}$) | `-0.00126316` / `0.23354233` | `-0.00096508` / `0.23293146` | Curved Dirichlet-to-Neumann exterior matching |
| Grazing impact parameter $b_{\text{graz}}$ ($\text{km}$) / $r_0$ ($\text{km}$) | `14.749032` / `11.558475` | `15.064914` / `11.736304` | Evaluated at $1.0001\,R_* e^{-(\nu(R_*) - \phi_L(R_*))}$ |
| Grazing deflection $\hat{\alpha}_{\text{Yuk,conf}}(b_{\text{graz}})$ ($\text{rad}$) | `1.29085207` | `1.33402365` | Conformal null geodesic at each branch's own grazing ray |
| Fixed-$b$ deflection $\hat{\alpha}_{\text{Yuk,conf}}(b = 15.064914\;\text{km})$ ($\text{rad}$) | `1.22474419` ($r_0 = 11.917144\;\text{km}$) | `1.33402365` ($r_0 = 11.736304\;\text{km}$) | Evaluated at shared exterior $b = b_{\text{graz}}^{\text{(II-Conf)}}$ ($+8.92\%$ from larger $M_{\text{ADM}}$) |
| Grazing relative shift vs. $\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})$ / $\hat{\alpha}_{\text{Schw}}(M_{\text{surf}})$ | `-2.867e-4` / `+3.652e-3` | `-1.671e-4` / `+2.184e-3` | Ray dives inside exterior scalar cloud ($m(r_0) < M_{\text{ADM}}$) |
| Far-field relative shift vs. $\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})$ ($b \ge 200\;\text{km}$) | `<= 8.03e-13` (`0.0` at $b \ge R_\odot$) | `<= 1.79e-13` (`0.0` at $b \ge R_\odot$) | Exponential screening $\exp(-b/L) \to 0$ on both branches |
| Interior–exterior surface slope mismatch $\lvert\phi'_{L,\text{int}}(R_*) - \phi'_{L,\text{ext}}(R_*)\rvert$ ($\text{km}^{-1}$) | `1.88e-09` (4-step Robin) | Robin: `1.465e-09` (4-step) / `3.937e-15` (every-step); FD2 ($N=4000$): `4.723e-10` | Enforced Robin BC vs. independent one-sided FD2 ($\mathcal{O}(h^2)$) |
| Converged interior scalar ODE residual $\max_{r\in(0,R_*)}\lvert\text{res}_{\text{ODE}}\rvert$ ($\text{km}^{-2}$) | `1.42e-10` (iter) / `1.95e-14` (banded) | `1.078e-10` (iter, `tol=2e-8`) / `6.958e-14` (`tol=1e-11`) / `1.480e-14` (banded) | Discrete residual of 2nd-order banded Helmholtz–Yukawa system |

#### 2B.1. Targeted Numerical Verification of the `Branch II-Conf` SLy $\log_{10}(\rho_{0c}) = 15.0$ BVP: Interior Residual, Surface $\phi'_L(R_*)$ Matching, and $N = 2000$ vs. $N = 4000$ Mesh Refinement

The script `gw170817_and_lensing_audit.py` (function `compute_branch_ii_conf_bvp_convergence_and_mesh_refinement`, stored under `gw170817_and_lensing_summary.json` $\to$ `bvp_convergence_and_mesh_refinement`) directly logs three numerical verification diagnostics for the single examined $L = 10\;\text{km}$, $\alpha_Y = 1$, SLy $\log_{10}(\rho_{0c}) = 15.0$ `Branch II-Conf` trace-sourced reference star BVP:

1. **Converged Interior Discrete Equation Residual ($\max_{r\in(0,R_*)}\lvert\text{res}_{\text{ODE}}\rvert$, $\text{km}^{-2}$) and Exterior Normalized Collocation Residual:**
   - On the converged $N = 601$ grid, the discrete residual of the interior scalar ODE $B\phi_L'' + (2B/r + B\nu' + B'/2)\phi_L' - \phi_L/L^2 + \alpha_Y^2 e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P}) = 0$ is **$\max\lvert\text{res}_{\text{ODE,banded}}\rvert = 1.480\times 10^{-14}\;\text{km}^{-2}$** for the direct tridiagonal banded solve (`_solve_curved_scalar_with_k`, which measures the accuracy of solving its own discrete linear system on that grid), **$\max\lvert\text{res}_{\text{ODE,iter}}\rvert = 1.078\times 10^{-10}\;\text{km}^{-2}$** for the relaxed fixed-point iterate (`tol_iter = 2e-8`, 19 iterations, $\max\lvert\phi_L - \phi_{L,\text{banded}}\rvert = 3.065\times 10^{-9}$), and **$\max\lvert\text{res}_{\text{ODE,iter}}\rvert = 6.958\times 10^{-14}\;\text{km}^{-2}$** on the separate `tol_iter = 1e-11` run with every-step 4-component exterior Robin updates ($N = 601$, 36 iterations). The exterior 4-component BVP (`scipy.integrate.solve_bvp`, `tol = 1e-12`) has maximum normalized RMS collocation residual **$\max(\text{rms\_residuals}) = 9.998\times 10^{-13}$** on $[R_*, R_{\max}]$ (a dimensionless normalized collocation residual on the exterior mesh, distinct from the interior $\text{km}^{-2}$ discrete ODE residual).
2. **Distinction Between the Enforced Robin Derivative ($\phi'_{L,\text{int,Robin}}(R_*)$) and the Independent One-Sided Finite-Difference Derivative ($\phi'_{L,\text{int,FD2}}(R_*)$) at $r = R_*$:**
   - **Enforced Robin boundary derivative ($\phi'_{L,\text{int,Robin}}(R_*) = -k_{\text{int}}\phi_L(R_*)$):** Under the 4-step Robin update schedule (`it in (1, 3, 6, 10)`, `tol_iter = 2e-8`, $N = 601$), the interior Robin coefficient is $k_{\text{int}} = 0.232931108388\;\text{km}^{-1}$ vs. the exterior 4-component curved Dirichlet-to-Neumann coefficient $k_{\text{curved}} = 0.232931461902\;\text{km}^{-1}$ ($\lvert\Delta k\rvert = 3.535\times 10^{-7}\;\text{km}^{-1}$, relative difference $1.518\times 10^{-6}$); the enforced Robin derivative $\phi'_{L,\text{int,Robin}}(R_*) = -0.000965074633\;\text{km}^{-1}$ and exterior $\phi'_{L,\text{ext,curved}}(R_*) = -0.000965076098\;\text{km}^{-1}$ differ by **$\lvert\Delta\phi'_{L,\text{Robin}}(R_*)\rvert = 1.465\times 10^{-9}\;\text{km}^{-1}$**. When the curved Robin coefficient is updated in **every outer iteration** (`it >= 1`) from the 4-component exterior BVP (`tol_iter = 1e-11`, $N = 601$, 36 iterations), $k_{\text{int}} = 0.232931462480\;\text{km}^{-1}$ and $k_{\text{curved}} = 0.232931462481\;\text{km}^{-1}$ (**$\lvert\Delta k\rvert = 9.502\times 10^{-13}\;\text{km}^{-1}$**), so the **enforced Robin derivative** matches to **$\lvert\phi'_{L,\text{int,Robin}}(R_*) - \phi'_{L,\text{ext,curved}}(R_*)\rvert = 3.937\times 10^{-15}\;\text{km}^{-1}$** (and $\lvert\phi'_{L,\text{int,ghost}}(R_*) - \phi'_{L,\text{ext,curved}}(R_*)\rvert = 1.430\times 10^{-14}\;\text{km}^{-1}$ for the ghost-point centered interior derivative), while shifting $R_*$ and $M_{\text{ADM}}$ by only $\lvert\Delta R_*\rvert = 2.06\times 10^{-7}\;\text{km}$ and $\lvert\Delta M_{\text{ADM}}\rvert = 6.07\times 10^{-8}\;M_\odot$.
   - **Independent 3-point one-sided second-order finite-difference derivative ($\phi'_{L,\text{int,FD2}}(R_*) = (3\phi_{N-1} - 4\phi_{N-2} + \phi_{N-3})/(2h)$):** Separated from the enforced Robin condition, the one-sided FD2 derivative differs from $\phi'_{L,\text{ext,curved}}(R_*)$ by $2.1213\times 10^{-8}\;\text{km}^{-1}$ at $N = 601$, $1.9048\times 10^{-9}\;\text{km}^{-1}$ at $N = 2000$, and **$4.7228\times 10^{-10}\;\text{km}^{-1}$ at $N = 4000$** ($\approx 4.89\times 10^{-7}$ relative to the slope $\lvert\phi'_L(R_*)\rvert = 9.6507\times 10^{-4}\;\text{km}^{-1}$), and converges in second order ($\mathcal{O}(h^2)$) as the mesh is refined ($N = 2000 \to 4000$ empirical improvement ratio: `4.0332`, consistent with the expected second-order mesh scaling $(h_{2000}/h_{4000})^2 = 4.0020$).
3. **Mesh Refinement Test ($N = 601$ vs. $N = 2000$ vs. $N = 4000$, `run_variant`: `DOP853`, `rtol_ivp = 1e-11`, `atol = 1e-13`, `max_step = 0.02 km`, `tol_iter = 1e-10`, Every-Step Robin Update):**

| Grid $N$ | $h = \Delta r$ ($\text{km}$) | $R_*$ ($\text{km}$) | $M_{\text{surf}}$ ($M_\odot$) | $M_{\text{ADM}}$ ($M_\odot$) | $\phi_L(0)$ | $\phi_L(R_*)$ | $\phi'_{L,\text{int,Robin}}(R_*)$ ($\text{km}^{-1}$) | $\phi'_{L,\text{ext,curved}}(R_*)$ ($\text{km}^{-1}$) | $\lvert\phi'_{L,\text{int,FD2}} - \phi'_{L,\text{ext}}\rvert$ ($\text{km}^{-1}$) | $\max\lvert\text{res}_{\text{ODE,iter}}\rvert$ ($\text{km}^{-2}$) | $\hat{\alpha}_{\text{Yuk,conf}}(b_{\text{graz}})$ ($\text{rad}$) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `601` | `0.019438` | `11.662698212` | `1.560635444` | `1.562034848` | `0.011576871` | `0.004143168` | `-0.000965074210` | `-0.000965074204` | `2.1213e-08` | `7.568e-13` | `1.334023631` |
| `2000` | `0.005834` | `11.662697490` | `1.560635270` | `1.562034651` | `0.011576835` | `0.004143134` | `-0.000965066327` | `-0.000965066321` | `1.9048e-09` | `7.844e-13` | `1.334023603` |
| `4000` | `0.002916` | `11.662697448` | `1.560635260` | `1.562034640` | `0.011576833` | `0.004143132` | `-0.000965065881` | `-0.000965065874` | `4.7228e-10` | `7.430e-13` | `1.334023601` |
| **$\lvert\Delta\rvert_{2000\to 4000}$** | — | **`4.213e-08`** (`3.61e-9`) | **`9.670e-09`** (`6.20e-9`) | **`1.097e-08`** (`7.03e-9`) | **`2.042e-09`** (`1.76e-7`) | **`1.920e-09`** (`4.63e-7`) | **`4.463e-10`** | **`4.463e-10`** | **Scaling: `4.0332`** (expected $\mathcal{O}(h^2)$: `4.0020`) | — | **`1.055e-09`** (`7.91e-10`) |

#### Table 3: Dedicated $1.000000\;M_\odot$ Solar–Cassini Benchmark ($GM_\odot/c^2 = 1.4766250385\;\text{km}$, $R_\odot = 695700.0\;\text{km}$, $r_E = 1\;\text{AU}$, $r_S = 8.43\;\text{AU}$; `source_body = Sun_1p00Msun`)

| Row (`scale_label`) | $b$ ($\text{km}$) | $b/R_\odot$ | $\hat{\alpha}_{\text{exact}}$ ($''$) | $\hat{\alpha}_{\text{1PN}}$ ($''$) | $\Delta t_{\text{Shapiro}}(10\,r_0)$ ($\mu\text{s}$) | $\Delta t_{\text{Cassini,rt}}$ ($\mu\text{s}$) | $\Delta t_{\text{Cassini,rt}}$ with geometric $-1$ ($\mu\text{s}$) | $\log_{10}\lvert\gamma_{\text{PPN}}-1\rvert$ |
|---|---|---|---|---|---|---|---|---|
| `Solar_1Msun_limb_1_Rsun` | `6.957000e+05` | `1.0000` | `1.751201` | `1.751190` | `67.8831` | `280.9434` | `261.2415` | `-30214.67` |
| `Solar_1Msun_Cassini_1p6_Rsun` | `1.113120e+06` | `1.6000` | `1.094498` | `1.094494` | `67.8830` | `262.4234` | `242.7215` | `-48342.99` |
| `Solar_1Msun_10_Rsun` | `6.957000e+06` | `10.0000` | `0.175119` | `0.175119` | `67.8829` | — | — | `-302139.47` |
| `Solar_1Msun_Earth_1_AU` | `1.495979e+08` | `215.0322` | `0.008144` | `0.008144` | `67.8829` | — | — | `-6496953.77` |

---

### 3.4. Weak-Lensing Convergence ($\kappa(\mathbf{\theta})$) on Galactic and Cluster Scales: Bridge to Module 2 (Bullet Cluster)

The weak-lensing convergence ($\kappa(\mathbf{\theta}) = \frac{1}{2}\nabla_{\mathbf{\theta}}^2\psi(\mathbf{\theta})$) is governed by the line-of-sight integral of the Weyl potential ($\Phi_{\text{lens}} = \frac{1}{2}(\tilde{\Psi} + \tilde{\Phi})$):

$$\kappa(\mathbf{b}_\perp) \;=\; \frac{1}{2\,c^2}\,\frac{D_L D_{LS}}{D_S}\int_{-\infty}^{+\infty}\nabla_\perp^2\bigl(\tilde{\Psi}(\mathbf{b}_\perp, z) + \tilde{\Phi}(\mathbf{b}_\perp, z)\bigr)\,dz \;=\; \frac{\Sigma_{\text{lens}}(\mathbf{b}_\perp)}{\Sigma_{\text{crit}}}, \qquad \Sigma_{\text{crit}} \;=\; \frac{c^2}{4\pi G}\,\frac{D_S}{D_L D_{LS}}.$$

Let us evaluate what the fixed parameter pair $L = 10\;\text{km}$, $\alpha_Y = 1$ predicts on galactic ($\sim 30\;\text{kpc}$, SPARC NGC 3198) and galaxy-cluster ($\sim 250\;\text{kpc}$, Bullet Cluster 1E 0657-558) scales (`gw170817_and_lensing_summary.json`):

1. **Behavior of the Local Fixed-$L = 10\;\text{km}$ Yukawa Component on kpc–Mpc Scales:**
   - **Exterior field of compact sources (stars, remnants):** At $b = 1\;\text{kpc}$ ($b/L = 3.086\times 10^{15}$) and $b = 250\;\text{kpc}$ ($b/L = 7.714\times 10^{17}$), the exponential screening factor $e^{-b/L} = 10^{-3.35\times 10^{17}} = 0$ exponentially suppresses to zero any local $L = 10\;\text{km}$ Yukawa tail.
   - **Equilibrium local field inside diffuse baryonic media (galactic disk or X-ray intracluster plasma, $\rho_{\text{ICM}} \sim 10^{-22}\;\text{kg/m}^3$):** In a quasi-homogeneous background of density $\rho$, the local Helmholtz–Yukawa equilibrium field is $\phi_L^{\text{local}}(L = 10\;\text{km}) = \alpha_Y \gamma_L \rho$:
     - In the outer disk of NGC 3198 ($\rho_{\text{gal}} \sim 10^{-24}\;\text{kg/m}^3$, $R = 30.08\;\text{kpc}$):
       $$\phi_L^{\text{local}}(10\;\text{km}) \;=\; 7.426\times 10^{-44}, \qquad \frac{|\Phi_N|}{c^2} \;\approx\; \frac{V_{\text{bar}}^2}{c^2} \;=\; 5.499\times 10^{-8} \quad\Longrightarrow\quad \frac{c^2\phi_L^{\text{local}}(10\;\text{km})}{|\Phi_N|} \;=\; 1.350\times 10^{-36}.$$
     - In the hot intracluster medium of the Bullet Cluster ($\rho_{\text{ICM}} \sim 10^{-22}\;\text{kg/m}^3$, $b = 250\;\text{kpc}$):
       $$\phi_L^{\text{local}}(10\;\text{km}) \;=\; 7.426\times 10^{-42} \quad\Longrightarrow\quad \kappa_{\text{lens}}^{\text{local }L=10\text{ km}}(\mathbf{b}_\perp) \;=\; \kappa_{\text{baryon}}(\mathbf{b}_\perp) \quad \left(\frac{\kappa_{\text{lens}}^{\text{local }L=10\text{ km}}}{\kappa_{\text{baryon}}} = 1.000000000000\right).$$

2. **Physical Conclusion and Boundary Condition for Module 2 (Bullet Cluster 1E 0657-558):**
   - The local $L = 10\;\text{km}$, $\alpha_Y = 1$ Yukawa interaction is strictly a **strong-field / compact-object ($r \lesssim 100\;\text{km}$) fifth force** that:
     - satisfies the GW170817 $c_T = c$ bound identically ($(c_T - c)/c = 0$),
     - suppresses scalar dipolar and monopolar radiation across binary-pulsar and LIGO inspiral frequencies via the $4771.35\;\text{Hz}$ mass cutoff,
     - satisfies the Solar-System Cassini bound at $b = 1.6\,R_\odot$ with $\lvert 1 - \gamma_{\text{PPN}}\rvert = 0 < 2.3\times 10^{-5}$,
     - and simultaneously produces **zero non-GR contribution** on kiloparsec and megaparsec scales ($b/L \sim 10^{16}\text{--}10^{18}$), neither to galactic rotation curves (where NGC 3198 at $R = 30.08\;\text{kpc}$ requires a dynamical boost $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2 = 4.311$ at $\Upsilon_* = 0.5$; `step6`) nor to the Bullet Cluster weak-lensing convergence map (where the galaxy subcluster peak exhibits $\kappa_{\text{obs}}/\kappa_{\text{baryon,gal}} \approx 0.232 / 0.0258 \approx 9.0$ at $\Sigma_{\text{crit}} \approx 3.10\times 10^9\;M_\odot\,\text{kpc}^{-2}$).
   - For a macroscopic clock field ($\phi_{\text{macro}}$, such as the chameleon/nonlocal $L_{\text{eff}}(\rho) \sim 10\text{--}300\;\text{kpc}$ branch analyzed in `step6`, `step11`, `step14`) to enhance weak-lensing convergence ($\kappa > \kappa_{\text{baryon}}$), **two conditions must hold simultaneously**:
     1. The metric coupling cannot be purely conformal ($\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu}$), because the conformal factor cancels in the Weyl potential $(\tilde{\Psi} + \tilde{\Phi})/2$; it requires a **disformal / refractive coupling** ($\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu} + B(\phi)\partial_\mu\phi\partial_\nu\phi$ or effective optical metric potential $\Phi_{\text{lens}}^{\text{eff}} = \Phi_N - \frac{1}{2}c^2\phi_{\text{macro}}$, `step14`); and
     2. The effective range of the scalar field must be environment-dependent ($L_{\text{eff}}(\rho) \propto \rho^{-1/2}$) or nonlocal, since the fixed $L = 10\;\text{km}$ Yukawa kernel vanishes identically for $b \gg 10\;\text{km}$.
