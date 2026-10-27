# Step 16: Derivation of the V21/V23 Algebraic X-Sector Anisotropic Stress Reconstruction (Mass and Gradient Modes, V23 Eqs. 4–5) from the Bianchi Identity and Comparison with the Dynamical Yukawa Scalar Stress Tensor ($T^\mu_{\;\nu}(\phi)$)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20_en.md` (`step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 1.1, lines 47–53, Eqs. 4–5 and Fig. 1(a–b); DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (Chapter 8, point 1, lines 400–408; DOI: `10.5281/zenodo.23604423`) and `hu_p3.tex` (Chapter 8, point 1, lines 6–12)

---

## 1. Verbatim Citations from the V22 and V23 Manuscripts

In Steps 4 and 7 (`step4_covariant_conservation_hydrostatic.md`, `step7_neutron_star_tov_stability.md`), we derived the local covariant stress–energy tensor ($T_{\mu\nu}^{(\phi)}$) of Formulation II (`yukawa_action_el`), but we had not yet derived step by step **Equations (4) and (5)** of the V23 manuscript (`eq:mass_X` and `eq:grad_X`), the non-zero central limit ($p_{tX}(0) = p_{rX}(0) = -(1-f_c)(\varepsilon_c/3 + P_c)$) in Chapter 8 of the V22 monograph, and the three-tier covariant divergence residual hierarchy (`conservation_Residuals.json`):

1. **Equations (4)–(5) of V23 and the Three-Tier Bianchi Residual (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, Section 1.1, lines 47–53):**
   > *"In a static spherically symmetric interior, if an auxiliary sector $T^\mu_{\;\nu}(X) = \text{diag}(-\varepsilon_X, P_{r,X}, P_{t,X}, P_{t,X})$ absorbs the window stress on the right-hand side of $G^\mu_{\;\nu} = \frac{8\pi G}{c^4}(T^\mu_{\;\nu,\text{matter}} + T^\mu_{\;\nu,X})$ while fluid elements satisfy $\nabla_\mu T^{\mu\nu}_{\text{matter}} = 0$ ($P' = -(\varepsilon+P)\nu'$ in geometrized units, $G=c=1$, $\nu = \Phi_{\text{eff}}/c^2$, $g_{tt}=-e^{2\nu}=-e^{2\Phi_{\text{eff}}/c^2}$), then the Bianchi identity $\nabla_\mu G^\mu_{\;r} \equiv 0$ uniquely fixes $T^\mu_{\;\nu}(X)$ in both V21 modes:"*
   $$\begin{align}
     \textbf{Mass mode } (m' = 4\pi r^2 f_L\varepsilon,\; \nu' = A_{\text{geo}}): &\quad \varepsilon_X = -(1-f_L)\varepsilon,\quad P_{r,X} = 0,\quad P_{t,X} = -\frac{r}{2}(1-f_L)\varepsilon\nu', \tag{V23 Eq. 4}\\
     \textbf{Gradient mode } (m' = 4\pi r^2\varepsilon,\; \nu' = f_L A_{\text{geo}}): &\quad \varepsilon_X = 0,\quad P_{r,X} = -(1-f_L)\!\left(P + \frac{m}{4\pi r^3}\right),\; P_{t,X} = P_{r,X} + \frac{r}{2}\bigl(P'_{r,X} + P_{r,X}\nu'\bigr), \tag{V23 Eq. 5}
   \end{align}$$
   > *"where $A_{\text{geo}}(r) = (m + 4\pi r^3 P)/[r(r-2m)]$. Symbolic verification in `SymPy` (`conservation_Residuals.json`) confirms that $\nabla_\mu T^\mu_{\;r}(\text{matter}) \equiv 0$ and $\nabla_\mu T^\mu_{\;r}(X) \equiv 0$ hold identically. In Table 1 and Fig. 1(b), evaluating relations (4)–(5) on numerical TOV grids ($L=10\,\text{km}$, $\log_{10}\rho_{0c}=15.0$) yields three distinct tiers (`conservation_Residuals.json`): (i) symbolic identities, residual exactly $0$; (ii) algebraic X-sector quadrature residual $\max|\nabla_\mu T^\mu_{\;r}(X)|/|\rho_X\nu'| = 4.44\times 10^{-16}$ (mass) and $/|P'_{r,X}| = 4.39\times 10^{-14}$ (gradient); (iii) divergence of the numerically (spline-) differentiated background $\max|\nabla_\mu T^\mu_{\;r}|/|P'| = 2.70\times 10^{-5}$ (mass) and $2.36\times 10^{-5}$ (gradient), governed by the finite-difference derivative of $h(r)$, not by the action."*

2. **Central Limit of the Gradient Mode in the V22 Monograph (`tal_book_HU_2026-10-24.tex`, Chapter 8, point 1, lines 404–408):**
   > *"The background X-sector differs between the two insertion modes:*
   > $$\textbf{Gradient mode:}\quad \rho_X = 0,\quad p_{rX} = -(1-f)\!\left(P + \frac{m}{4\pi r^3}\right),\quad p_{tX} = p_{rX} + \frac{r}{2}\bigl(p'_{rX} + p_{rX}\nu'\bigr),$$
   > *whose central limit $p_{tX}(0) = p_{rX}(0) = -(1-f_c)(\varepsilon_c/3 + P_c)$ is generally non-zero. In the **mass-source mode**, by contrast, $\rho_X = -(1-f)\varepsilon$, $p_{rX} = 0$, and $p_{tX} = -\frac{r}{2}(1-f)\varepsilon\nu'$ ($p_{tX}(0)=0$). The value $0.140$ in the 36-configuration stress-diagnostics table is a background bookkeeping diagnostic of that configuration, which does not replace a covariant action principle (`50aa2db0` task)."*

---

## 2. Exact Derivation of V23 Equations (4) and (5) from the Einstein Equations and the Bianchi Identity

In a static, spherically symmetric spacetime ($ds^2 = -e^{2\nu(r)}dt^2 + e^{2\lambda(r)}dr^2 + r^2 d\Omega^2$, in geometrized units $G = c = 1$, where $e^{-2\lambda(r)} \equiv 1 - \frac{2m(r)}{r}$), the non-vanishing mixed components of the Einstein tensor are:
$$G^t_{\;t} \;=\; -\frac{2m'(r)}{r^2}, \qquad G^r_{\;r} \;=\; -\frac{2m(r)}{r^3} \;+\; \frac{2}{r}\left(1 - \frac{2m(r)}{r}\right)\nu'(r).$$
Let us write the Einstein equations for the sum of the baryonic perfect fluid ($T^\mu_{\;\nu,\text{matter}} = \text{diag}(-\varepsilon, P, P, P)$) and the anisotropic $X$-sector ($T^\mu_{\;\nu,X} = \text{diag}(-\varepsilon_X, P_{r,X}, P_{t,X}, P_{t,X})$):
$$G^t_{\;t} \;=\; -8\pi\bigl(\varepsilon + \varepsilon_X\bigr), \qquad G^r_{\;r} \;=\; 8\pi\bigl(P + P_{r,X}\bigr), \qquad G^\theta_{\;\theta} \;=\; G^\phi_{\;\phi} \;=\; 8\pi\bigl(P + P_{t,X}\bigr).$$
If the fluid satisfies $P'(r) = -(\varepsilon + P)\nu'(r)$, then $\nabla_\mu T^\mu_{\;r}(\text{matter}) = P' + (\varepsilon + P)\nu' \equiv 0$. Because $\nabla_\mu G^\mu_{\;r} \equiv 0$ by the contracted Bianchi identity, the $X$-sector by itself must also satisfy the radial covariant conservation equation:
$$\boxed{\;\nabla_\mu T^\mu_{\;r}(X) \;=\; P'_{r,X}(r) \;+\; \bigl[\varepsilon_X(r) + P_{r,X}(r)\bigr]\nu'(r) \;+\; \frac{2}{r}\bigl[P_{r,X}(r) - P_{t,X}(r)\bigr] \;=\; 0.\;}$$

### 2.1. Derivation in the V21/V23 `mass` (Mass-Source) Mode (V23 Eq. 4)
In the `mass` mode, the mass function and lapse gradient are defined by:
$$m'(r) \;=\; 4\pi r^2 f_L(r)\,\varepsilon(r), \qquad \nu'(r) \;=\; A_{\text{geo}}(r) \;\equiv\; \frac{m(r) + 4\pi r^3 P(r)}{r\bigl(r - 2m(r)\bigr)}.$$
1. **From the $(tt)$ Einstein Equation ($\varepsilon_X$):**
   $$-\frac{2m'}{r^2} = -8\pi f_L\varepsilon = -8\pi(\varepsilon + \varepsilon_X) \quad\implies\quad \varepsilon_X(r) \;=\; -\bigl(1 - f_L(r)\bigr)\,\varepsilon(r).$$
2. **From the $(rr)$ Einstein Equation ($P_{r,X}$):**
   Since $\nu' = A_{\text{geo}}$, the expression for $G^r_{\;r}$ yields $8\pi P(r)$ directly, hence:
   $$8\pi P(r) = 8\pi\bigl(P(r) + P_{r,X}(r)\bigr) \quad\implies\quad P_{r,X}(r) \;\equiv\; 0.$$
3. **From the Radial Bianchi Identity ($P_{t,X}$):**
   Substituting $P_{r,X} = 0$ and $P'_{r,X} = 0$ into $\nabla_\mu T^\mu_{\;r}(X) = 0$ gives:
   $$\varepsilon_X(r)\,\nu'(r) \;-\; \frac{2}{r}P_{t,X}(r) \;=\; 0 \quad\implies\quad \boxed{\;P_{t,X}(r) \;=\; \frac{r}{2}\,\varepsilon_X(r)\,\nu'(r) \;=\; -\frac{r}{2}\bigl(1 - f_L(r)\bigr)\,\varepsilon(r)\,\nu'(r),\;}$$
   which is precisely **V23 Equation (4)** (`eq:mass_X`), and since $\nu'(r) = \mathcal{O}(r)$ at the origin, $P_{t,X}(r) = \mathcal{O}(r^2) \xrightarrow{r\to 0} 0$.

### 2.2. Derivation in the V21/V23 `gradient` Mode (V23 Eq. 5) and Proof of Central Isotropy
In the `gradient` mode, the mass function and lapse gradient are defined by:
$$m'(r) \;=\; 4\pi r^2\,\varepsilon(r), \qquad \nu'(r) \;=\; f_L(r)\,A_{\text{geo}}(r) \;=\; f_L(r)\,\frac{m(r) + 4\pi r^3 P(r)}{r\bigl(r - 2m(r)\bigr)}.$$
1. **From the $(tt)$ Einstein Equation ($\varepsilon_X$):**
   $$-\frac{2m'}{r^2} = -8\pi\varepsilon = -8\pi(\varepsilon + \varepsilon_X) \quad\implies\quad \varepsilon_X(r) \;\equiv\; 0.$$
2. **From the $(rr)$ Einstein Equation ($P_{r,X}$):**
   Substituting $\nu'(r) = f_L(r)\frac{m + 4\pi r^3 P}{r(r-2m)}$ into $G^r_{\;r} = 8\pi(P + P_{r,X})$ gives:
   $$P(r) + P_{r,X}(r) \;=\; \frac{G^r_{\;r}}{8\pi} \;=\; -\frac{m(r)}{4\pi r^3} \;+\; f_L(r)\left(\frac{m(r)}{4\pi r^3} + P(r)\right) \quad\implies\quad \boxed{\;P_{r,X}(r) \;=\; -\bigl(1 - f_L(r)\bigr)\left(P(r) + \frac{m(r)}{4\pi r^3}\right).\;}$$
3. **From the Radial Bianchi Identity ($P_{t,X}$):**
   Since $\varepsilon_X = 0$, solving $\nabla_\mu T^\mu_{\;r}(X) = P'_{r,X} + P_{r,X}\nu' + \frac{2}{r}(P_{r,X} - P_{t,X}) = 0$ for $P_{t,X}(r)$ yields:
   $$\boxed{\;P_{t,X}(r) \;=\; P_{r,X}(r) \;+\; \frac{r}{2}\Bigl[P'_{r,X}(r) \;+\; P_{r,X}(r)\,\nu'(r)\Bigr],\;}$$
   which is precisely **V23 Equation (5)** (`eq:grad_X`).
4. **Central Limit and Isotropy ($r \to 0$):**
   Near the stellar center ($r \to 0$), $m(r) = \frac{4\pi}{3}\varepsilon_c r^3 + \mathcal{O}(r^5)$, so $\lim_{r\to 0}\frac{m(r)}{4\pi r^3} = \frac{\varepsilon_c}{3}$. Because the derivative of every regular scalar profile vanishes at the origin ($P'(0) = 0$, $f_L'(0) = 0$, $\nu'(0) = 0$), we have $P'_{r,X}(r) = \mathcal{O}(r)$, so the anisotropy term $\frac{r}{2}[P'_{r,X} + P_{r,X}\nu'] = \mathcal{O}(r^2)$ vanishes at the origin, and:
   $$\boxed{\;\lim_{r\to 0} P_{t,X}(r) \;=\; \lim_{r\to 0} P_{r,X}(r) \;=\; -\bigl(1 - f_c\bigr)\left(\frac{\varepsilon_c}{3} + P_c\right) \;\neq\; 0, \qquad \lim_{r\to 0}\bigl[P_{t,X}(r) - P_{r,X}(r)\bigr] \;=\; 0.\;}$$
   Hence no geometric conical singularity arises at the origin: although the central pressure is non-zero, the stress of the $X$-sector at the origin is **isotropic** ($P_{r,X}(0) = P_{t,X}(0)$)!

---

## 3. Comparison with the Dynamical Scalar Stress Tensor of Formulation II (`yukawa_action_el`)

In Formulation II of V23 (`yukawa_action_el`, `step4` and `step7`), the $X$-sector is determined not by algebraic reconstruction, but directly from the covariant action ($S = \int d^4x\sqrt{-g}\bigl[\frac{R}{16\pi} - \frac{(\nabla\phi_L)^2 + \phi_L^2/L^2}{2\alpha_Y} - f_L(\phi_L)\varepsilon(n)\bigr]$). Relative to the standard GR source $(\varepsilon, P)$, the resulting effective $X$-sector consists of two dynamical parts:
$$T^\mu_{\;\nu,X,\text{dyn}} \;=\; \bigl(T^\mu_{\;\nu,\text{matter}} - T^\mu_{\;\nu,\text{GR}}\bigr) \;+\; T^\mu_{\;\nu}(\phi),$$
that is, component by component ($B(r) = 1 - 2m(r)/r$):
$$\begin{aligned}
  \varepsilon_{X,\text{dyn}}(r) &\;=\; -\bigl(1 - f_L(r)\bigr)\,\varepsilon(r) \;+\; \frac{1}{2\alpha_Y}\left(B(r)\bigl[\phi_L'(r)\bigr]^2 + \frac{\phi_L(r)^2}{L^2}\right), \\
  P_{r,X,\text{dyn}}(r) &\;=\; -\bigl(1 - f_L(r)\bigr)\,P(r) \;+\; \frac{1}{2\alpha_Y}\left(B(r)\bigl[\phi_L'(r)\bigr]^2 - \frac{\phi_L(r)^2}{L^2}\right), \\
  P_{t,X,\text{dyn}}(r) &\;=\; -\bigl(1 - f_L(r)\bigr)\,P(r) \;-\; \frac{1}{2\alpha_Y}\left(B(r)\bigl[\phi_L'(r)\bigr]^2 + \frac{\phi_L(r)^2}{L^2}\right).
\end{aligned}$$
Three deep structural relationships between Formulations I and II follow immediately:
1. **Origin of the $-(1-f_L)\varepsilon$ Term in the Mass Mode:** The first term of $\varepsilon_{X,\text{dyn}}(r)$ is the exact dilution term $\varepsilon_X = -(1-f_L)\varepsilon$ of the V21/V23 `mass` mode, supplemented by the positive self-energy density $\rho_\phi = \frac{1}{2\alpha_Y}(B[\phi_L']^2 + \phi_L^2/L^2)$ of the scalar field.
2. **Purely Gradient Origin of the Anisotropy ($\Pi_\phi \equiv P_{r,X,\text{dyn}} - P_{t,X,\text{dyn}}$):**
   $$\boxed{\;\Pi_\phi(r) \;\equiv\; P_{r,X,\text{dyn}}(r) - P_{t,X,\text{dyn}}(r) \;=\; p_{r,\phi}(r) - p_{t,\phi}(r) \;=\; \frac{B(r)}{\alpha_Y}\bigl[\phi_L'(r)\bigr]^2 \;\ge\; 0.\;}$$
   Because $\phi_L'(r) = \mathcal{O}(r)$ near the origin, the dynamical anisotropy also obeys $\Pi_\phi(r) = \mathcal{O}(r^2) \xrightarrow{r\to 0} 0$, so at the origin the pressure is strictly isotropic: $P_{r,X,\text{dyn}}(0) = P_{t,X,\text{dyn}}(0) = -(1-f_c)P_c - \frac{\phi_L(0)^2}{2\alpha_Y L^2} < 0$.
3. **Why the Euler–Lagrange Equation for $\phi_L$ Replaces Separate Bianchi Reconstruction:**
   Whereas in the algebraic V21 modes $P_{t,X}(r)$ had to be reconstructed a posteriori from $\nabla_\mu T^\mu_{\;r}(X) = 0$, in the `yukawa_action_el` system the Noether identity $\nabla_\mu\bigl(T^\mu_{\;r,\text{matter}} + T^\mu_{\;r}(\phi)\bigr) \equiv 0$ proved in Step 4 (`step4_covariant_conservation_hydrostatic.md`) is satisfied automatically by virtue of the field equation for $\phi_L$!

---

## 4. Checks

1. **Symbolic Bianchi Identity Check for Both V23 Algebraic Modes:**
   - **Mass mode:** $P'_{r,X} + (\varepsilon_X + P_{r,X})\nu' + \frac{2}{r}(P_{r,X} - P_{t,X}) = 0 + \varepsilon_X\nu' + \frac{2}{r}\bigl(0 - \frac{r}{2}\varepsilon_X\nu'\bigr) = \varepsilon_X\nu' - \varepsilon_X\nu' \equiv 0$ $\checkmark$.
   - **Gradient mode:** $P'_{r,X} + (0 + P_{r,X})\nu' + \frac{2}{r}\Bigl(P_{r,X} - \bigl[P_{r,X} + \frac{r}{2}(P'_{r,X} + P_{r,X}\nu')\bigr]\Bigr) = P'_{r,X} + P_{r,X}\nu' - (P'_{r,X} + P_{r,X}\nu') \equiv 0$ $\checkmark$.
2. **$f_L \to 1$ ($\phi_L \to 0$) GR Limit:**
   - In the `mass` mode ($\varepsilon_X, P_{r,X}, P_{t,X} \to 0$), in the `gradient` mode ($\varepsilon_X, P_{r,X}, P_{t,X} \to 0$), and in the dynamical `yukawa_action_el` case ($\varepsilon_{X,\text{dyn}}, P_{r,X,\text{dyn}}, P_{t,X,\text{dyn}} \to 0$), the $X$-sector vanishes identically, recovering the isotropic GR TOV equations $\checkmark$.
