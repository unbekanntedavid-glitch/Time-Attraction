# Step 12: Second Variation ($\delta^2 S$) of the Scalar-Tensor / Yukawa Action, Explicit Coefficients of the Coupled $(\xi, \delta\phi_L)$ Radial Pulsation Operator, and Derivation of the Positive-Definite Physical Kinetic Hessian

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step12_second_variation_coupled_pulsation_hessian_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step11_chameleon_screening_multiscale_eft_synthesis_en.md` (`step11_chameleon_screening_multiscale_eft_synthesis.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p3.tex` (Chapter 8, lines 13–23, Eq. 12, and Chapter 11 task **(b)**, line 102; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 1.2, line 65 and Section 2, lines 120–153, Eqs. 8–9, Table 2; DOI: `10.5281/zenodo.23639084`)  
- V21 Zenodo record: DOI `10.5281/zenodo.23602501`  
- Numerical implementation and verification on the `yukawa_action_el` sequence: `curved_yukawa_adm_and_stability.py`, `curved_yukawa_coupled_stability.csv`, and `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`

---

## 1. Verbatim Citations from the V22 Monograph and the V23 (`50aa2db0`) Manuscript

In Step 7 (`step7_neutron_star_tov_stability.md`), we derived the curved-spacetime `yukawa_action_el` TOV background equations as well as the **static ($\omega_0^2 = 0$) nonlocal tangent model** ($\delta f_L[\delta h]$) around the V21 background. Both open task **(b)** of Chapter 11 of the V22 monograph (`hu_p3.tex`, line 102) and Section 2 of the V23 manuscript (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, lines 134–136) explicitly distinguish the static tangent model (and its Gram matrix) from the **genuine coupled dynamical pulsation system and physical kinetic Hessian derived from the second variation $\delta^2 S$ of the action**:

1. **Theoretical Task (b) of Chapter 11 of the V22 Monograph (`hu_p3.tex`, line 102):**
   > *"**(b) Genuine Modified Radial Stability Analysis Derived from the Coupling:** Once task (a) provides the dynamical field equations of the perturbations $\delta\phi$ and $\delta T^{\mu\nu}_X$, the Chandrasekhar radial pulsation operator must be re-derived variationally from the second variation $\delta^2 S$ of the action, eliminating the non-variational metric operator pathology flagged above $L \ge 367\,\text{m}$ and lifting the frozen-$f$ ($\delta f_L = 0$) approximation."*

2. **Non-Variational Character of the V21 `gradient`/`metric` Operator Above $L \ge 367\text{ m}$ (`hu_p3.tex`, Chapter 8, line 21):**
   > *"**`gradient`/`metric` operator (flagged as non-variational above $L \ge 367\,\text{m}$):** Substitutes the metric lapse derivative $\nu' = 2\Phi_{\text{eff}}'$ into the standard GR curvature terms of $\mathcal{Q}(r)$ without accounting for the perturbation of the $X$-sector tangential stress $\delta p_{tX}$. My audit openly flags this operator as **non-positive-definite and non-variational above $L \ge 367\,\text{m}$ across 12 sequences (609 models)**, where uncompensated residual terms produce unphysical eigenvalues $\omega_0^2 < 0$ even along the rising branch $\partial M/\partial\rho_c > 0$."*

3. **Inertia Invariance of the V23 Static Tangent Model and the Difference Between the Gram Matrix and the Physical Kinetic Hessian (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, Section 2, lines 134–136):**
   > *"Because at $\omega_0^2=0$ the problem reduces to $\mathcal{H}_{\text{lin}}[\delta h]=0$, this agreement verifies static nonlocal linearization rather than dynamical inertia: scaling $W\to 4W$ at $\log_{10}\rho_{0c}=15.30427$ yields a ratio $\omega^2_{4W}/\omega^2_W = 0.2500$ ($\omega^2_W=4.897\times 10^{6}\,\text{s}^{-2}$) with the identical shooting function $F(0)=0.023252$, so the zero crossing is inertia-invariant. The raw BVP $W$ matrix is not symmetric; the symmetrized Galerkin form $\mathbf{M}_{\text{sym}} = \mathbf{D}^{T}\mathrm{diag}(w\,\Delta r)\mathbf{D}$ is symmetric with $\min\mathrm{diag}\,\mathbf{M}_{\text{sym}} = 5.68\times 10^{-9}>0$; this is a Gram-matrix property, not a proof of positive-definiteness of the physical kinetic Hessian. The `omega0_sq_sym_tangent_s2` column gives the first-order tangent estimate ($\omega^2\simeq -F(0)/F'(0)$)."*

---

## 2. Mathematical Derivation: The Second Variation $\delta^2 S$, Elimination of Metric Constraints, and Explicit Coefficients

### 2.1. Time-Dependent Spherically Symmetric Perturbations and the Metric Constraints $(\delta\lambda, \delta m, \delta\nu')$
Consider the covariant action of V23 Formulation II (`yukawa_action_el`, **Branch II-Fluid**, $G = c = 1$) in a time-dependent, spherically symmetric spacetime ($ds^2 = -e^{2\nu(r,t)}dt^2 + e^{2\lambda(r,t)}dr^2 + r^2 d\Omega^2$, where on the static background $B_0(r) = e^{-2\lambda_0(r)} = 1 - 2m_0(r)/r$):
$$S[g_{\mu\nu}, \phi_L, \Psi_{\text{fluid}}] \;=\; \int dt\int_0^\infty dr\,4\pi r^2 e^{\nu+\lambda}\left[\frac{R}{16\pi} - \frac{1}{2\alpha_Y}\left(-e^{-2\nu}(\dot{\phi}_L)^2 + e^{-2\lambda}(\phi_L')^2 + \frac{\phi_L^2}{L^2}\right) - e^{-\phi_L}\varepsilon(n)\right].$$
Let us introduce small radial Lagrangian fluid displacements $\xi(r,t)$ (for which the enthalpy perturbation is $\delta h(r,t) = \frac{\delta P(r,t)}{\varepsilon_0(r) + P_0(r)}$) and the Eulerian scalar-field perturbation $\delta\phi_L(r,t)$ around the static equilibrium background ($\nu_0(r), \lambda_0(r), \phi_{L,0}(r), \varepsilon_0(r), P_0(r)$):
$$r_{\text{fluid}}(r,t) \;=\; r + \xi(r,t), \qquad \phi_L(r,t) \;=\; \phi_{L,0}(r) + \delta\phi_L(r,t).$$

1. **Metric Constraints from the $(tr)$ and $(tt)$ Einstein Equations ($\delta\lambda$ and $\delta m$):**
   Integrating the spherically symmetric $(tr)$ Einstein equation ($\delta G_t^{\;r} = 8\pi\,\delta T_t^{\;r}$) with respect to time in the interior ($0 < r < R_*$) gives:
   $$\frac{\partial(\delta\lambda)}{\partial t} \;=\; -4\pi r e^{2\lambda_0}\left[f_{L,0}(\varepsilon_0 + P_0)\,\dot{\xi} - \frac{1}{\alpha_Y}\phi_{L,0}'\,\delta\dot{\phi}_L\right],$$
   whereas from the $(tt)$ Einstein equation ($\delta m'(r,t) = 4\pi r^2\,\delta(-T_t^{\;t})$, where $\delta B = -2\delta m / r = -2 B_0\,\delta\lambda$) at fixed baryon number ($\delta N_B = 0$), the radial derivative of the mass perturbation in variables $(\delta h, \delta\phi_L, \delta\phi_L', \delta m)$ is:
   $$\boxed{\;\frac{\partial(\delta m)}{\partial r} \;=\; 4\pi r^2\left[f_{L,0}\frac{\varepsilon_0+P_0}{c_s^2}\,\delta h - f_{L,0}\varepsilon_0\,\delta\phi_L + \frac{1}{\alpha_Y}\left(B_0\phi_{L,0}'\,\delta\phi_L' + \frac{\phi_{L,0}}{L^2}\delta\phi_L - \frac{(\phi_{L,0}')^2}{r}\delta m\right)\right].\;}$$
   Equivalently, expressed in terms of the Lagrangian displacement $\xi(r)$ (by substituting the baryon-number conservation constraint $\delta n = -\frac{1}{r^2 e^{\lambda_0}}(r^2 e^{\lambda_0} n_0\xi)' - n_0\delta\lambda$ and integrating by parts), the mass perturbation is:
   $$\delta m(r) \;=\; -4\pi r^2 f_{L,0}(\varepsilon_0 + P_0)\,\xi(r) \;+\; \int_0^r 4\pi \tilde{r}^2\left[\frac{B_0\phi_{L,0}'}{\alpha_Y}\delta\phi_L' + \left(\frac{\phi_{L,0}}{\alpha_Y L^2} - f_{L,0}\varepsilon_0\right)\delta\phi_L + f_{L,0}(\varepsilon_0+P_0)\phi_{L,0}'\,\xi - \frac{(\phi_{L,0}')^2}{\alpha_Y\tilde{r}}\delta m\right]d\tilde{r}.$$
2. **Lapse-Gradient Constraint from the $(rr)$ Einstein Equation ($\delta\nu'$):**
   Because $\nu_0'(r) = \frac{m_0 + 4\pi r^3 P_{\text{tot},0}}{r(r-2m_0)}$, where $P_{\text{tot},0} = f_{L,0}P_0 + p_{r,\phi,0}$ and $p_{r,\phi,0} = \frac{1}{2\alpha_Y}[B_0(\phi_{L,0}')^2 - \phi_{L,0}^2/L^2]$, the total pressure perturbation is:
   $$\delta P_{\text{tot}} \;=\; f_{L,0}(\varepsilon_0 + P_0)\,\delta h - f_{L,0}P_0\,\delta\phi_L + \frac{1}{\alpha_Y}\left(B_0\phi_{L,0}'\,\delta\phi_L' - \frac{\phi_{L,0}}{L^2}\delta\phi_L - \frac{(\phi_{L,0}')^2}{r}\delta m\right),$$
   from which the perturbation of the lapse gradient follows:
   $$\boxed{\;\delta\nu'(r) \;=\; \frac{1 + 2 r \nu_0'(r)}{r\bigl(r - 2m_0(r)\bigr)}\,\delta m(r) \;+\; \frac{4\pi r^3}{r\bigl(r - 2m_0(r)\bigr)}\,\delta P_{\text{tot}}(r).\;}$$

### 2.2. Second-Order Action $\delta^2 S$ and the Positive-Definite Physical Kinetic Hessian ($\mathbf{T}_{\text{kin}}$)
After substituting the metric constraints back into the quadratic second-order action ($\delta^2 S = \int dt\,(T_{\text{kin}}[\dot{\xi}, \delta\dot{\phi}_L] - U_{\text{pot}}[\xi, \delta\phi_L])$), the kinetic part under harmonic time dependence ($\xi(r,t) = \xi(r)e^{i\omega t}$, $\delta\phi_L(r,t) = \delta\phi_L(r)e^{i\omega t}$) is:
$$\boxed{\;T_{\text{kin}}\bigl[\xi, \delta\phi_L\bigr] \;=\; \frac{\omega^2}{2}\int_0^{R_*} dr\,\mathcal{W}_\xi(r)\,\xi(r)^2 \;+\; \frac{\omega^2}{2}\int_0^\infty dr\,\mathcal{W}_\phi(r)\,\bigl[\delta\phi_L(r)\bigr]^2,\;}$$
where the two physical inertia (kinetic) weight functions are:
$$\boxed{\;\mathcal{W}_\xi(r) \;=\; 4\pi r^2\,e^{3\lambda_0(r) - \nu_0(r)}\,f_{L,0}(r)\,\bigl(\varepsilon_0(r) + P_0(r)\bigr) \;>\; 0 \quad (0 < r < R_*), \qquad \mathcal{W}_\phi(r) \;=\; \frac{4\pi r^2}{\alpha_Y}\,e^{\lambda_0(r) - \nu_0(r)} \;>\; 0 \quad (r > 0).\;}$$

### 2.3. Explicit Coefficients of the Coupled $(\xi, \delta\phi_L)$ and $(\delta h, \delta\phi_L, \delta m)$ Eigenvalue System (Point 3)
The coupled radial pulsation system can be written in two strictly equivalent forms — the second-order Lagrangian $(\xi, \delta\phi_L)$ form and the first-/second-order canonical $(\delta h, \delta\phi_L, \delta m)$ form:

#### (A) Explicit Second-Order Coupled $(\xi, \delta\phi_L)$ Sturm–Liouville / Klein–Gordon Form:
$$\begin{cases}
  -\dfrac{d}{dr}\!\left(\mathcal{P}_\xi(r)\dfrac{d\xi}{dr}\right) + \mathcal{Q}_{\xi\xi}(r)\,\xi(r) + \mathcal{C}_{\xi\phi}\!\left[\delta\phi_L\right](r) \;=\; \omega^2\,\mathcal{W}_\xi(r)\,\xi(r), \\[8pt]
  -\dfrac{d}{dr}\!\left(\mathcal{P}_\phi(r)\dfrac{d(\delta\phi_L)}{dr}\right) + \mathcal{Q}_{\phi\phi}(r)\,\delta\phi_L(r) + \mathcal{C}_{\phi\xi}\!\left[\xi\right](r) \;=\; \omega^2\,\mathcal{W}_\phi(r)\,\delta\phi_L(r),
\end{cases}$$
where the principal-part coefficients are:
$$\mathcal{P}_\xi(r) \;=\; 4\pi r^2\,e^{\nu_0(r) + \lambda_0(r)}\,f_{L,0}(r)\,\Gamma_1(r)\,P_0(r) \;>\; 0, \qquad \mathcal{P}_\phi(r) \;=\; \frac{4\pi r^2}{\alpha_Y}\,e^{\nu_0(r) - \lambda_0(r)} \;>\; 0,$$
the fluid–fluid potential coefficient ($\mathcal{Q}_{\xi\xi}$, with the modified hydrostatic gradient $P_0' = -(\varepsilon_0+P_0)(\nu_0' - \phi_{L,0}')$ and the combined background pressure/energy $\rho_{\text{tot},0} = f_{L,0}\varepsilon_0 + \rho_{\phi,0}$, $P_{\text{tot},0} = f_{L,0}P_0 + p_{r,\phi,0}$) is:
$$\boxed{\;\begin{aligned}
  \mathcal{Q}_{\xi\xi}(r) &\;=\; 4\pi r^2 e^{\nu_0+\lambda_0}f_{L,0}\left\{\frac{4}{r}P_0' + \frac{(P_0')^2}{\varepsilon_0 + P_0} - 8\pi e^{2\lambda_0}(\varepsilon_0 + P_0)\bigl(f_{L,0}P_0 + p_{r,\phi,0}\bigr)\right. \\
  &\qquad\left. -\; \frac{d}{dr}\!\left(\frac{2\Gamma_1 P_0}{r}\right) + (\varepsilon_0 + P_0)\left[\phi_{L,0}'' + \bigl(\nu_0' + \lambda_0' - \phi_{L,0}'\bigr)\phi_{L,0}'\right]\right\},
\end{aligned}\;}$$
the scalar–scalar potential coefficient ($\mathcal{Q}_{\phi\phi}$, including the metric back-reaction from $\delta\nu'$ and $\delta B$) is:
$$\boxed{\;\mathcal{Q}_{\phi\phi}(r) \;=\; \frac{4\pi r^2}{\alpha_Y}\,e^{\nu_0+\lambda_0}\left[\frac{1}{L^2} + \alpha_Y\,f_{L,0}(r)\,\varepsilon_0(r) + \frac{8\pi r\,e^{2\lambda_0}}{\alpha_Y L^2}\phi_{L,0}(r)\phi_{L,0}'(r)\right] \;+\; \mathcal{Q}_{\phi\phi}^{\text{metric}}[\delta m],\;} $$
and the mutually adjoint fluid–scalar cross-coupling operators ($\langle \eta, \mathcal{C}_{\xi\phi}[\delta\phi_L]\rangle = \langle \mathcal{C}_{\phi\xi}[\eta], \delta\phi_L\rangle$) are:
$$\boxed{\;\begin{aligned}
  \mathcal{C}_{\xi\phi}[\delta\phi_L](r) &\;=\; 4\pi r^2 e^{\nu_0+\lambda_0}f_{L,0}(\varepsilon_0+P_0)\left[-\delta\phi_L'(r) + \delta\nu'\big|_{\delta\phi_L, \delta m[\delta\phi_L]}(r)\right], \\[4pt]
  \mathcal{C}_{\phi\xi}[\xi](r) &\;=\; -\frac{d}{dr}\!\Bigl(4\pi r^2 e^{\nu_0+\lambda_0}f_{L,0}(\varepsilon_0+P_0)\,\xi(r)\Bigr) \;+\; 4\pi r^2 e^{\nu_0+\lambda_0}f_{L,0}\frac{\varepsilon_0+P_0}{c_s^2}\,\delta h[\xi](r) \;+\; \mathcal{C}_{\phi\xi}^{\text{metric}}[\delta m[\xi]](r).
\end{aligned}\;}$$

#### (B) Exact Canonical $(\delta h, \delta\phi_L, \delta m)$ System and Moving-Surface Matching Conditions (Implemented in `curved_yukawa_adm_and_stability.py`):
For numerical integration and exact handling of polytrope transitions ($r_{\text{trans}}$), the canonical variables $(\delta h, \delta\phi_L, \delta m)$ are optimal (since $\delta h$ and $\delta m$ remain continuous even across jump points of $c_s^2$!):
$$\boxed{\;\begin{aligned}
  \frac{d(\delta h)}{dr} &\;=\; -\delta\nu'(r) \;+\; \delta\phi_L'(r) \;+\; \omega^2\,e^{2(\lambda_0-\nu_0)}\,\xi(r), \\[4pt]
  \frac{d(\delta m)}{dr} &\;=\; 4\pi r^2\left[f_{L,0}\frac{\varepsilon_0+P_0}{c_s^2}\,\delta h - f_{L,0}\varepsilon_0\,\delta\phi_L + \frac{1}{\alpha_Y}\left(B_0\phi_{L,0}'\,\delta\phi_L' + \frac{\phi_{L,0}}{L^2}\delta\phi_L - \frac{(\phi_{L,0}')^2}{r}\delta m\right)\right], \\[4pt]
  B_0\,\delta\phi_L'' &\;+\; \left(\frac{2B_0}{r} + B_0\nu_0' + \frac{B_0'}{2}\right)\delta\phi_L' - \left(\frac{1}{L^2} + \alpha_Y f_{L,0}\varepsilon_0 - \omega^2 e^{-2\nu_0}\right)\delta\phi_L \\
  &\;=\; -\alpha_Y f_{L,0}\frac{\varepsilon_0+P_0}{c_s^2}\,\delta h + \frac{2\,\delta m}{r}\phi_{L,0}'' - \left[\delta(B\nu') + \frac{\delta B'}{2} - \frac{2\,\delta m}{r^2}\right]\phi_{L,0}'.
\end{aligned}\;}$$
- **Moving-Surface Boundary and Matching Conditions ($r = R_*$):**
  1. At the free surface, the Lagrangian pressure perturbation vanishes: $\Delta P(R_*) = \delta P(R_*) + P_0'(R_*)\xi(R_*) = 0 \implies \delta h(R_*) = -h_0'(R_*)\xi(R_*)$, which together with integrated mass conservation ($\delta m(R_*^-) = -4\pi R_*^2 f_{L,0}(R_*)(\varepsilon_0+P_0)\big|_{R_*^-}\xi(R_*) = 0$, since for a polytropic crust $\varepsilon_0(R_*^-) = P_0(R_*^-) = 0$) guarantees $\delta m(R_*) \in C^0$.
  2. Across the moving surface $r = R_* + \xi(R_*)$, the scalar field and its normal derivative are continuous:
     $$\llbracket \delta\phi_L + \xi\phi_{L,0}' \rrbracket_{R_*} \;=\; 0, \qquad \llbracket \delta\phi_L' + \xi\phi_{L,0}'' \rrbracket_{R_*} \;=\; 0.$$
     Since $\phi_{L,0} \in C^1([0,\infty))$, $\llbracket\phi_{L,0}'\rrbracket_{R_*} = 0 \implies \llbracket\delta\phi_L\rrbracket_{R_*} = 0$; furthermore, $\llbracket\phi_{L,0}''\rrbracket_{R_*} = \frac{\alpha_Y}{B_0(R_*)}f_{L,0}(R_*)\varepsilon_0(R_*^-) = 0$ at the surface of the polytropic crust ($\varepsilon_0(R_*^-) = 0$), so $\llbracket\delta\phi_L'\rrbracket_{R_*} = 0$ holds as well!
  3. In the curved exterior vacuum ($r \in [R_*, R_{\max}]$), $\varepsilon_0 = P_0 = 0$, so the $(\delta\phi_L, \delta m)$ system is integrated in coupled form out to $R_{\max} = R_* + 12 L$, where the asymptotic Yukawa decay condition closes the boundary-value problem, and the asymptotic mass perturbation $\delta M_{\text{ADM}}(\omega^2) = \delta m(R_{\max}; \omega^2) = 0$ determines the eigenvalue!

---

## 3. Numerical Verification on the Curved Exterior `yukawa_action_el` Sequence (Point 4)

The script `curved_yukawa_adm_and_stability.py` solves the coupled $(\delta h, \delta\phi_L, \delta m)$ eigenvalue problem ($r \in [0, R_* + 14L]$, with the exact left-hand side sign `scalar_m_sign = -1.0`) directly on the `yukawa_action_el` equilibrium sequences equipped with the curved exterior solution. The results (`curved_yukawa_coupled_stability.csv` and `curved_yukawa_audit_summary.json`) are:

| Sequence | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dynamical $\omega_0^2=0$ | $\Delta_{\text{ADM}}$ Dynamical (dex) | $\Delta_{M_B}$ Dynamical (dex) | Quasi-Static $\omega_0^2=0$ | Frozen $\delta\phi_L=0$ | $\Delta_{\text{ADM}}$ Frozen (dex) |
|---|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316409` | **`15.316410`** | **`+5.80e-06`** | `+1.33e-06` | `15.316410` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.307169` | **`15.307170`** | **`+4.20e-06`** | `+6.41e-07` | `15.307170` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303082` | **`15.303080`** | **`-6.86e-06`** | `-2.61e-06` | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.293025` | **`15.293025`** | **`+3.57e-06`** | `-1.40e-07` | `15.293025` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227787` | **`15.227783`** | **`-1.34e-06`** | `-4.05e-06` | `15.227783` | `15.251384` | `+2.36e-02` |

*(GL2 Volterra, $N=480$, `scalar_m_sign = -1.0`, dynamical scalar wave, finite-difference verified exterior $\partial\phi_L'(R_*)/\partial m(R_*)=-1.68365\times 10^{-4}\;\text{km}^{-2}$ and interior $r=6\;\text{km}$ $\partial\phi_L''/\partial m=-1.72833\times 10^{-4}\;\text{km}^{-3}$ derivatives; $N$- and $R_{\max}$-convergence on the root: $p_{\text{eff,root}}=4.248$, $\log_{10}\rho_{0c,\text{Rich}}(\omega_0^2=0)=15.31640981$, see `step22`, Sections 4 and 7.)*

---

## 4. Checks

1. **$f_{L,0} \to 1$, $\alpha_Y \to 0$ (GR Limit):**
   - For $\phi_{L,0} \to 0$ and $\delta\phi_L \to 0$, $\mathcal{W}_\xi(r) \to 4\pi r^2 e^{3\lambda_0-\nu_0}(\varepsilon_0 + P_0)$, $\mathcal{P}_\xi(r) \to 4\pi r^2 e^{\nu_0+\lambda_0}\Gamma_1 P_0$, and $\mathcal{Q}_{\xi\xi}(r)$ recovers the classical Chandrasekhar (1964) general-relativistic radial pulsation Sturm–Liouville operator $\checkmark$.
2. **Comparison Between the Coupled $(\delta h, \delta\phi_L, \delta m)$ Eigenvalue Zero Crossing and the ADM/Baryonic Mass Peaks on the `yukawa_action_el` Branch:**
   - For SLy $L = 10\;\text{km}$, $\alpha_Y = 1$, the frozen-$\phi_L$ ($\delta\phi_L = 0$) approximation shifts the zero crossing by $+1.22\times 10^{-2}\;\text{dex}$ from the $dM_{\text{ADM}}/d\rho_c = 0$ turning point ($15.328565$ vs. $15.316404$), whereas the coupled $(\delta h, \delta\phi_L, \delta m)$ operator (`scalar_m_sign = -1.0`) yields $\log_{10}\rho_{0c}(\omega_0^2=0)=15.316410$ ($+5.80\times 10^{-6}\;\text{dex}$ above the ADM mass peak and $+1.33\times 10^{-6}\;\text{dex}$ above the baryonic mass peak, a $2097\times$ reduction relative to the frozen-$\phi_L$ shift; across all 5 sequences the agreement is within $\le 6.9\times 10^{-6}\;\text{dex}$, a $165\times\text{--}17593\times$ reduction) $\checkmark$.
3. **Positivity of the Kinetic Weights Above $L \ge 367\text{ m}$:**
   - Because $f_{L,0}(r) = e^{-\phi_{L,0}(r)} > 0$ and $\alpha_Y > 0$, $\mathcal{W}_\xi(r) > 0$ and $\mathcal{W}_\phi(r) > 0$ hold across every scale $L \in (0, \infty)$ without sign reversal $\checkmark$.
