# Step 4: Covariant Stress-Energy Conservation, Scalar-Matter Momentum Exchange, and Modified Hydrostatic Equilibrium ($h' = -\nu' + \phi_L'$)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step4_covariant_conservation_hydrostatic_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`), `step2_conformal_disformal_en.md` (`step2_conformal_disformal.md`), `step3_yukawa_helmholtz_eom_en.md` (`step3_yukawa_helmholtz_eom.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` Hungarian manuscript): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Subsections 1.1 and 1.2, lines 47–65, Eqs. 4–5 and Eqs. 6–7; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (Chapter 2, lines 205–207; DOI: `10.5281/zenodo.23604423`)  
- V23 SymPy verification script: `covariant_closure_50aa2db0.py` (lines 72–98 and 143–160) and `conservation_Residuals.json`

---

## 1. Verbatim Citations from the V22 and V23 Manuscripts

In Step 3 (`step3_yukawa_helmholtz_eom.md`), we derived the Euler–Lagrange equation of the scalar/dilaton field $\phi_L$. The next necessary theoretical step is the variation with respect to the metric $g^{\mu\nu}$ and the derivation of **covariant energy–momentum conservation following from the Bianchi identity ($\nabla_\mu G^{\mu\nu} \equiv 0$)** ($\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_X) = 0$), together with the resulting **modified hydrostatic equilibrium (fifth force / scalar-gradient force)**.

1. **Formulation of the Physical Problem in the V22 Monograph (`hu_p1.tex`, line 206):**
   > *"At the same time, inserting $f_L(\mathbf{x})$ as an external multiplier into the Tolman–Oppenheimer–Volkoff equations or the clock-frequency formulas leaves behind a residual $\nabla_\mu T^{\mu\nu}_{\text{matter}} \neq 0$ unless balanced by a dynamical field sector. Therefore, I formulated the dynamical 'why' as task `50aa2db0`: elevating the nonlocal kernel into a covariant scalar–tensor or nonlocal action principle $S = \int d^4x\sqrt{-g}\,\bigl[\frac{c^4}{16\pi G}R + \mathcal{L}_{\text{matter}}(\psi, \tilde{g}_{\mu\nu}) + \mathcal{L}_\phi(\phi,\nabla\phi)\bigr]$ whose Euler–Lagrange equation recovers the profile $\phi(\mathbf{x}) \to \gamma_L\bar{\rho}_L(\mathbf{x})$ in the quasi-static limit while guaranteeing the conservation law $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_X) = 0$."*

2. **Algebraic $X$-Sector Reconstruction in Formulation I of V23 (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, lines 47–52, Eqs. 4–5):**
   > *"In a static spherically symmetric interior, if the auxiliary sector $T^\mu_{\;\nu}(X) = \text{diag}(-\varepsilon_X, P_{r,X}, P_{t,X}, P_{t,X})$ absorbs the window stress on the right-hand side of $G^\mu_{\;\nu} = \frac{8\pi G}{c^4}(T^\mu_{\;\nu,\text{matter}} + T^\mu_{\;\nu,X})$ while the fluid elements satisfy $\nabla_\mu T^{\mu\nu}_{\text{matter}} = 0$ ($P' = -(\varepsilon+P)\nu'$ in geometrized units, $G=c=1$, $\nu = \Phi_{\text{eff}}/c^2$, $g_{tt}=-e^{2\nu}=-e^{2\Phi_{\text{eff}}/c^2}$), then the Bianchi identity $\nabla_\mu G^\mu_{\;r} \equiv 0$ uniquely fixes $T^\mu_{\;\nu}(X)$ in both V21 modes:"*
   $$\begin{aligned}
     \textbf{Mass mode } (m' = 4\pi r^2 f_L\varepsilon,\; \nu' = A_{\text{geo}}): &\quad \varepsilon_X = -(1-f_L)\varepsilon,\quad P_{r,X} = 0,\quad P_{t,X} = -\frac{r}{2}(1-f_L)\varepsilon\nu', \tag{V23 Eq. 4}\\
     \textbf{Gradient mode } (m' = 4\pi r^2\varepsilon,\; \nu' = f_L A_{\text{geo}}): &\quad \varepsilon_X = 0,\quad P_{r,X} = -(1-f_L)\!\left(P + \frac{m}{4\pi r^3}\right),\; P_{t,X} = P_{r,X} + \frac{r}{2}\bigl(P'_{r,X} + P_{r,X}\nu'\bigr), \tag{V23 Eq. 5}
   \end{aligned}$$
   > *"where $A_{\text{geo}}(r) = (m + 4\pi r^3 P)/[r(r-2m)]$. Symbolic verification performed in `SymPy` (`conservation_Residuals.json`) confirms that $\nabla_\mu T^\mu_{\;r}(\text{matter}) \equiv 0$ and $\nabla_\mu T^\mu_{\;r}(X) \equiv 0$ hold identically."*

3. **Local Scalar–Tensor (Yukawa) Stress–Energy Tensor and Action-Consistent Hydrostatic Equilibrium in Formulation II (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 65):**
   > *"and the variation $\delta S_\phi/\delta g^{\mu\nu}$ provides the local covariant stress–energy tensor $T_{\mu\nu}^{(\phi)} = \frac{c^4}{8\pi G\tilde{\alpha}_L}\bigl[\nabla_\mu\phi_L\nabla_\nu\phi_L - \frac{1}{2}g_{\mu\nu}\bigl((\nabla\phi_L)^2 + L^{-2}\phi_L^2\bigr)\bigr]$. [...] The `yukawa_action_el` rows follow from the action $\mathcal{L}=-[(\nabla\phi)^2+\phi^2/L^2]/(2\alpha_Y)-f_L\varepsilon$ ($G=c=1$; `SymPy`: `sym_yukawa_el_source_factor` and `_hydrostatic` residual $=0$): $B\phi''+(2B/r+B\nu'+B'/2)\phi'-\phi/L^2=-\alpha_Y f_L\varepsilon$, $B=1-2m/r$, $P_{\text{matter}}=f_LP$, $h'=-\nu'+\phi'$, $\rho_\phi=(B\phi'^2+\phi^2/L^2)/(2\alpha_Y)$, exterior matching $\phi'(R)=-(1/R+1/L)\phi(R)$."*

---

## 2. Derivation of Scalar–Matter Momentum Exchange and Modified Hydrostatic Equilibrium

### 2.1. Derivation of the Effective Matter Pressure ($P_{\text{matter}} = f_L P$)
In the V23 `yukawa_action_el` branch (and in general Einstein-frame scalar–tensor theory where the matter Lagrangian is $\mathcal{L}_{\text{fluid}}(n,\phi_L) = -f_L(\phi_L)\,\varepsilon(n) = -e^{-\phi_L}\varepsilon(n)$, with $n$ denoting the baryon number density), the thermodynamic pressure is defined via the Legendre transformation under baryon-number conservation:
$$P_{\text{matter}} \;=\; n\,\frac{\partial\bigl(f_L(\phi_L)\,\varepsilon(n)\bigr)}{\partial n} \;-\; f_L(\phi_L)\,\varepsilon(n) \;=\; f_L(\phi_L)\Bigl(n\,\frac{d\varepsilon}{dn} - \varepsilon(n)\Bigr) \;=\; f_L(\phi_L)\,P(n).$$
Thus the stress–energy tensor of the $f_L$-weighted ideal fluid is:
$$T^{\mu\nu}_{\text{matter}} \;=\; f_L(\phi_L)\,\Bigl[\bigl(\varepsilon + P\bigr)\frac{u^\mu u^\nu}{c^2} + P\,g^{\mu\nu}\Bigr] \quad\implies\quad \varepsilon_{\text{matter}} = f_L\varepsilon, \qquad P_{\text{matter}} = f_L P.$$

### 2.2. Divergence of the Scalar-Field Stress–Energy Tensor ($\nabla^\mu T_{\mu\nu}^{(\phi)}$)
In geometrized units ($G = c = 1$, where $\alpha_Y = 8\pi\tilde{\alpha}_L$, so that $\frac{c^4}{8\pi G\tilde{\alpha}_L} = \frac{1}{\alpha_Y}$), the stress–energy tensor of the scalar field is:
$$T_{\mu\nu}^{(\phi)} \;=\; \frac{1}{\alpha_Y}\left[\nabla_\mu\phi_L\nabla_\nu\phi_L - \frac{1}{2}g_{\mu\nu}\left(g^{\alpha\beta}\nabla_\alpha\phi_L\nabla_\beta\phi_L + \frac{1}{L^2}\phi_L^2\right)\right].$$
Let us compute its covariant divergence ($\nabla^\mu T_{\mu\nu}^{(\phi)}$):
$$\nabla^\mu T_{\mu\nu}^{(\phi)} \;=\; \frac{1}{\alpha_Y}\left[\left(\Box_g\phi_L\right)\nabla_\nu\phi_L + \left(\nabla^\mu\phi_L\right)\nabla_\mu\nabla_\nu\phi_L - \left(\nabla^\alpha\phi_L\right)\nabla_\nu\nabla_\alpha\phi_L - \frac{1}{L^2}\phi_L\nabla_\nu\phi_L\right].$$
Since second covariant derivatives commute on a scalar field ($\nabla_\mu\nabla_\nu\phi_L = \nabla_\nu\nabla_\mu\phi_L$), the second and third terms cancel:
$$\nabla^\mu T_{\mu\nu}^{(\phi)} \;=\; \frac{1}{\alpha_Y}\left(\Box_g\phi_L - \frac{1}{L^2}\phi_L\right)\nabla_\nu\phi_L.$$
Now substitute the exact `yukawa_action_el` Euler–Lagrange field equation derived in Step 3 (`step3_yukawa_helmholtz_eom.md`):
$$\left(-\Box_g + \frac{1}{L^2}\right)\phi_L \;=\; \alpha_Y\,f_L(\phi_L)\,\varepsilon \quad\Longleftrightarrow\quad \Box_g\phi_L - \frac{1}{L^2}\phi_L \;=\; -\alpha_Y\,f_L\,\varepsilon.$$
With this, the divergence of the scalar sector in closed form is:
$$\boxed{\;\nabla^\mu T_{\mu\nu}^{(\phi)} \;=\; -f_L(\phi_L)\,\varepsilon\,\nabla_\nu\phi_L \;=\; \varepsilon\,\nabla_\nu f_L(\phi_L).\;}$$

### 2.3. Joint Conservation ($\nabla^\mu(T_{\mu\nu}^{\text{matter}} + T_{\mu\nu}^{(\phi)}) = 0$) and the Modified Euler/TOV Equation
Because of the contracted Bianchi identity ($\nabla^\mu G_{\mu\nu} \equiv 0$) on the left-hand side of the Einstein equations ($G_{\mu\nu} = 8\pi(T_{\mu\nu}^{\text{matter}} + T_{\mu\nu}^{(\phi)})$), the joint divergence of the right-hand side vanishes identically:
$$\nabla^\mu T_{\mu\nu}^{\text{matter}} \;+\; \nabla^\mu T_{\mu\nu}^{(\phi)} \;=\; 0 \quad\implies\quad \nabla^\mu T_{\mu\nu}^{\text{matter}} \;=\; +f_L(\phi_L)\,\varepsilon\,\nabla_\nu\phi_L.$$
In a static spherically symmetric spacetime ($ds^2 = -e^{2\nu(r)}dt^2 + B(r)^{-1}dr^2 + r^2 d\Omega^2$, $B(r) = 1 - 2m(r)/r$), the radial ($\nu = r$) covariant divergence of the matter tensor is:
$$\nabla_\mu T^\mu_{\;r}(\text{matter}) \;=\; \frac{d}{dr}\bigl(f_L P\bigr) \;+\; f_L(\varepsilon + P)\,\nu' \;=\; f_L\,P' \;-\; f_L\,P\,\phi_L' \;+\; f_L(\varepsilon + P)\,\nu',$$
where we used $f_L' = \frac{d}{dr}(e^{-\phi_L}) = -f_L\,\phi_L'$. Adding the divergence of the scalar sector $\nabla_\mu T^\mu_{\;r}(\phi) = -f_L\,\varepsilon\,\phi_L'$ gives:
$$\nabla_\mu T^\mu_{\;r}(\text{matter}) + \nabla_\mu T^\mu_{\;r}(\phi) \;=\; f_L\Bigl[P' + (\varepsilon + P)\nu' - (\varepsilon + P)\phi_L'\Bigr] \;=\; 0.$$
Since $f_L = e^{-\phi_L} > 0$, we can divide by $f_L(\varepsilon + P)$, and introducing the specific enthalpy ($dh = \frac{dP}{\varepsilon + P}$) yields:

$$\boxed{\;h'(r) \;\equiv\; \frac{P'(r)}{\varepsilon(r) + P(r)} \;=\; -\nu'(r) \;+\; \phi_L'(r) \;=\; -\nu'(r) \;-\; \frac{d\ln f_L(r)}{dr}.\;}$$

---

## 3. Physical Interpretation in the Scalar–Tensor / Dilaton Framework

The relation $h' = -\nu' + \phi_L'$ derived above is the direct dynamical bridge between scalar–tensor gravity (Damour & Esposito-Farèse 1992; Fujii & Maeda 2003) and the Time Attraction Law:
1. **Decomposition of the Effective Potential:** The net equilibrium force acting on a baryonic fluid element consists of two parts:
   - the standard metric gravitational attraction $-\nu'(r) = -\frac{1}{c^2}\frac{d\Phi_{\text{eff}}}{dr}$, and
   - the **scalar / time-coherence gradient force** $+\phi_L'(r) = -\frac{d\ln f_L}{dr}$ (referred to in the scalar–tensor literature as the *fifth force* / dilaton-gradient force).
2. **Weak-Field / Galactic Limit:**
   In a non-relativistic medium ($P \ll \varepsilon \simeq \rho c^2$), the effective acceleration of a test body or fluid element is:
   $$\mathbf{a}_{\text{eff}} \;=\; -\nabla\Phi_N \;+\; c^2\,\nabla\phi_L \;=\; -\nabla\left(\Phi_N - c^2\phi_L\right) \;=\; -\nabla\left(\Phi_N + c^2\ln f_L\right).$$
   Because the density around a mass concentration increases inward ($\frac{d\bar{\rho}_L}{dr} < 0$), $\phi_L(r)$ increases inward ($\phi_L'(r) < 0$, or $\frac{d\ln f_L}{dr} > 0$), so the $+\phi_L'(r)$ term provides a density-gradient-driven contribution in both hydrostatic equilibrium and orbital dynamics.

---

## 4. Checks

1. **SymPy Symbolic Identities in the V23 Package (`covariant_closure_50aa2db0.py`, lines 72–98 and 150–160; `conservation_Residuals.json`):**
   - `sym_yukawa_el_matter_pressure_residual = 0`: verifies the identity $P_{\text{matter}} = n\frac{\partial(f_L\varepsilon)}{\partial n} - f_L\varepsilon = f_L P$ $\checkmark$.
   - `sym_yukawa_el_hydrostatic_residual = 0`: verifies that substituting the Yukawa EL equation into $\nabla_\mu T^\mu_{\;r}(\text{matter}) + \nabla_\mu T^\mu_{\;r}(\phi) = 0$ yields $h' = -\nu' + \phi_L'$ exactly $\checkmark$.
   - `div_matter_mass_symbolic = 0`, `div_X_mass_symbolic = 0`, `div_matter_grad_symbolic = 0`, `div_X_grad_symbolic = 0`: verifies conservation of the algebraic $X$-sector in V23 Formulation I Eqs. (4)–(5) $\checkmark$.
2. **General Relativistic (GR) Limit ($\phi_L \to 0$, $f_L \to 1$):**
   When $\phi_L \equiv 0$ ($\phi_L' = 0$, $f_L = 1$), the scalar stress–energy tensor vanishes ($T_{\mu\nu}^{(\phi)} = 0$), the effective pressure reduces to $P_{\text{matter}} \to P$, and hydrostatic equilibrium returns to the standard GR Tolman–Oppenheimer–Volkoff equation ($h' = -\nu'$, i.e., $P' = -(\varepsilon + P)\nu'$) $\checkmark$.
