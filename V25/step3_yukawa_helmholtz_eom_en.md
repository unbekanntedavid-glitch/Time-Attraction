# Step 3: Mapping the V23 Local Screened Helmholtz–Yukawa Action and Field Equation to the Massive Scalar-Tensor / Dilaton Klein–Gordon Equation

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step3_yukawa_helmholtz_eom_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`), `step2_conformal_disformal_en.md` (`step2_conformal_disformal.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` Hungarian manuscript): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Subsection 1.2, lines 56–65, Eqs. 6–7; DOI: `10.5281/zenodo.23639084`)  
- V23 (`50aa2db0` English manuscript): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` (Sec. 1.2, lines 56–65, Eqs. 6–7)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (lines 200–206 and 245–254; DOI: `10.5281/zenodo.23604423`)  
- V23 SymPy verification script: `covariant_closure_50aa2db0.py` (lines 126–149) and `conservation_Residuals.json`

---

## 1. Verbatim Citations from the V23 (`50aa2db0`) and V22 Manuscripts

In Step 1 (`step1_dilaton_form.md`), we recorded as an open question how the nonlocal exponent $\phi_L(\mathbf{x}) = \gamma_L\bar{\rho}_L(\mathbf{x})$ becomes a local dynamical scalar field endowed with its own kinetic term. **Subsection 1.2** of the V23 manuscript (*Formulation II: Local Screened Scalar-Tensor (Helmholtz–Yukawa) Action*, lines 56–65, **Eqs. 6–7**) defines this verbatim as follows:

- **Lagrangian and Field Equation of the Scalar Field (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, lines 57–60, Eq. 6):**
  > *"In Formulation II, I replace the nonlocal integral kernel with a local canonical screened scalar action:"*
  $$\begin{equation} \mathcal{L}_\phi = -\frac{c^4}{8\pi G \tilde{\alpha}_L}\!\left[\frac{1}{2}g^{\mu\nu}\nabla_\mu\phi_L\nabla_\nu\phi_L + \frac{1}{2L^2}\phi_L^2\right], \quad \left(-\Box_g + \frac{1}{L^2}\right)\phi_L = \frac{4\pi G}{c^2}\alpha_{\text{std}}\,\rho_{\text{src}} \iff (-L^2\nabla^2+1)\phi_L = \alpha_Y\gamma_L\rho_{\text{src}}. \tag{V23 Eq. 6} \end{equation}$$

- **Static Spherically Symmetric Green's Function Solution and Scalar Stress–Energy Tensor (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, lines 61–65, Eq. 7):**
  > *"In a static spherical system, the exact positive-definite Green's function solution of the radial Helmholtz operator $-\frac{1}{r^2}\frac{d}{dr}(r^2\phi_L') + L^{-2}\phi_L = \frac{4\pi G}{c^2}\alpha_{\text{std}}\rho_{\text{src}}$ ($\alpha_Y = 4\pi\alpha_{\text{std}}$) is"*
  $$\begin{equation} \phi_L^{\text{Yuk}}(r) = \alpha_Y\,\frac{G}{c^2}\int_0^{R_*} dr'\,\frac{r' L}{r}\,\sinh\!\frac{r_<}{L}\,e^{-r_>/L}\rho_{\text{src}}(r'), \qquad \phi_L^{\text{Yuk}}(0) = \alpha_Y\,\frac{G}{c^2}\int_0^{R_*} r' e^{-r'/L}\rho_{\text{src}}(r')\,dr', \tag{V23 Eq. 7} \end{equation}$$
  > *"and the variation $\delta S_\phi/\delta g^{\mu\nu}$ provides the local covariant stress–energy tensor $T_{\mu\nu}^{(\phi)} = \frac{c^4}{8\pi G\tilde{\alpha}_L}\bigl[\nabla_\mu\phi_L\nabla_\nu\phi_L - \frac{1}{2}g_{\mu\nu}\bigl((\nabla\phi_L)^2 + L^{-2}\phi_L^2\bigr)\bigr]$. (`SymPy`: `sym_yukawa_center_limit_residual` $=0$, no extra $L$-factor). The unit-normalized kernel corresponds to $\alpha_Y = 1$, i.e., $\alpha_{\text{std}} = 1/(4\pi)\approx 0.0796$, which recovers the homogeneous top-hat limit ($\phi_L = \gamma_L\rho_{\text{src}}$); $\alpha_Y = 1/3$ ($\alpha_{\text{std}} = 1/(12\pi)$) is the volume-matched variant."*

- **Action-Consistent (`yukawa_action_el`) Euler–Lagrange Branch (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 65):**
  > *"The `yukawa_action_el` rows follow from the action $\mathcal{L}=-[(\nabla\phi)^2+\phi^2/L^2]/(2\alpha_Y)-f_L\varepsilon$ ($G=c=1$; `SymPy`: `sym_yukawa_el_source_factor` and `_hydrostatic` residual $=0$): $B\phi''+(2B/r+B\nu'+B'/2)\phi'-\phi/L^2=-\alpha_Y f_L\varepsilon$, $B=1-2m/r$, $P_{\text{matter}}=f_LP$, $h'=-\nu'+\phi'$, $\rho_\phi=(B\phi'^2+\phi^2/L^2)/(2\alpha_Y)$, exterior matching $\phi'(R)=-(1/R+1/L)\phi(R)$."*

---

## 2. Derivation of the Euler–Lagrange Field Equation and Mapping to a Massive Dilaton / Scalar-Tensor Field

### 2.1. Variation with Respect to $\phi_L$ (with Matter-Sector Coupling $f_L(\phi_L) = e^{-\phi_L}$)
Consider the local action of V23 Formulation II, where the matter energy density couples with the factor $f_L(\phi_L) = e^{-\phi_L}$ derived in Steps 1 and 2 ($\mathcal{L}_m = -f_L(\phi_L)\,\rho c^2 = -e^{-\phi_L}\varepsilon$):
$$S_{\text{II}}[g_{\mu\nu},\phi_L,\Psi] \;=\; \int d^4x\sqrt{-g}\left[\frac{c^4}{16\pi G}R \;-\; \frac{c^4}{8\pi G\tilde{\alpha}_L}\left(\frac{1}{2}g^{\mu\nu}\nabla_\mu\phi_L\nabla_\nu\phi_L + \frac{1}{2L^2}\phi_L^2\right) \;-\; e^{-\phi_L}\,\rho c^2\right].$$
Varying the action with respect to $\phi_L$ ($\delta S_{\text{II}} / \delta\phi_L = 0$):
$$\frac{1}{\sqrt{-g}}\frac{\delta S_{\text{II}}}{\delta\phi_L} \;=\; \frac{c^4}{8\pi G\tilde{\alpha}_L}\left(\Box_g\phi_L - \frac{1}{L^2}\phi_L\right) \;+\; e^{-\phi_L}\,\rho c^2 \;=\; 0,$$
where $\Box_g\phi_L = \frac{1}{\sqrt{-g}}\partial_\mu\bigl(\sqrt{-g}\,g^{\mu\nu}\partial_\nu\phi_L\bigr)$. Rearranging and multiplying by $-\frac{8\pi G\tilde{\alpha}_L}{c^4}$ yields:
$$\left(-\Box_g + \frac{1}{L^2}\right)\phi_L \;=\; \frac{8\pi G}{c^2}\,\tilde{\alpha}_L\,\rho_{\text{src}}, \qquad \rho_{\text{src}} \;=\; \begin{cases} \rho & \text{(weak-field / linear source, } f_L \simeq 1\text{)}, \\[3pt] f_L(\phi_L)\,\rho = e^{-\phi_L}\rho & \text{(exact } \texttt{yukawa\_action\_el} \text{ source)}. \end{cases}$$

### 2.2. Identification of the Constants ($\tilde{\alpha}_L$, $\alpha_{\text{std}}$, $\alpha_Y$, $\gamma_L$, and $C$)
Comparing term by term with V23 Eq. (6):
1. **Relation between the Poisson and Yukawa normalizations:**
   $$\frac{8\pi G}{c^2}\,\tilde{\alpha}_L \;=\; \frac{4\pi G}{c^2}\,\alpha_{\text{std}} \;=\; C\,\alpha_{\text{std}} \quad\implies\quad \tilde{\alpha}_L \;=\; \frac{\alpha_{\text{std}}}{2},$$
   where $C = \frac{4\pi G}{c^2}$ is the fundamental coupling constant $C$ of the V22 monograph (`hu_p1.tex`, Eq. 4).
2. **Factoring out $\gamma_L = \frac{G L^2}{c^2}$ ($-L^2\Box_g + 1$ form):**
   Multiplying both sides of the field equation by $L^2$ (in the quasi-static limit $\Box_g \to \nabla^2$):
   $$(-L^2\nabla^2 + 1)\phi_L \;=\; 4\pi\alpha_{\text{std}}\,\left(\frac{G L^2}{c^2}\right)\rho_{\text{src}} \;=\; \alpha_Y\,\gamma_L\,\rho_{\text{src}}, \qquad \boxed{\;\alpha_Y \;\equiv\; 4\pi\,\alpha_{\text{std}} \;=\; 8\pi\,\tilde{\alpha}_L.\;}$$
   Or in an equivalent notation (as written in Section 1.1 of `README.md`):
   $$(\nabla^2 - L^{-2})\phi_L \;=\; -\frac{4\pi G}{c^2}\alpha_{\text{std}}\,\rho_{\text{src}} \qquad \left(\text{or in the }\phi_L \to -\phi_L\text{ sign convention, } +\frac{\alpha_{\text{std}}G}{c^2}\rho_{\text{src}}\text{ with the } 4\pi\text{ factor absorbed}\right).$$

### 2.3. Mapping to the Massive Scalar-Tensor / Dilaton Klein–Gordon Equation
In string-theoretic dilaton gravity or general scalar–tensor gravity (e.g., Damour & Polyakov 1994; Fujii & Maeda 2003; Khoury & Weltman 2004), the field equation of a canonical scalar field of mass $m_\phi$ with a quadratic self-interaction potential ($V(\varphi) = \frac{1}{2}\mu_\phi^2\varphi^2$, where $\mu_\phi = \frac{m_\phi c}{\hbar}$) is:
$$\left(-\Box_g + \mu_\phi^2\right)\phi_L \;=\; \frac{4\pi G}{c^2}\,\alpha_{\text{std}}\,\rho_{\text{src}}.$$
From this, the exact dictionary correspondence between the TAL smoothing scale ($L$) and the effective Compton wavelength and mass of the scalar/dilaton field is:
$$\boxed{\;\mu_\phi \;\equiv\; \frac{1}{L}, \qquad m_\phi \;=\; \frac{\hbar}{L\,c}\;}$$
- **Physical interpretation:** What the nonlocal formulation of V19/V21 describes as a 3D spherical top-hat windowing of radius $L$ corresponds in V23 Formulation II to the exchange of a **Yukawa-screened scalar/dilaton field with range $\lambda_{\text{Compton}} = L$ (i.e., mass $m_\phi = \hbar/(Lc)$)**.

---

## 3. How Does the Yukawa Solution Recover the 3D Top-Hat Law $\phi_L = \gamma_L\bar{\rho}_L$?

1. **The 3D Yukawa Green's Function:**
   The fundamental solution of $(-L^2\nabla^2 + 1)G_L(\mathbf{x}-\mathbf{y}) = \delta^{(3)}(\mathbf{x}-\mathbf{y})$ in three-dimensional Euclidean space is:
   $$G_L(\mathbf{x}-\mathbf{y}) \;=\; \frac{\exp\!\left(-\frac{|\mathbf{x}-\mathbf{y}|}{L}\right)}{4\pi L^2\,|\mathbf{x}-\mathbf{y}|}, \qquad \int_{\mathbb{R}^3} G_L(\mathbf{r})\,d^3r \;=\; \frac{4\pi}{4\pi L^2}\int_0^\infty r\,e^{-r/L}\,dr \;=\; \frac{L^2}{L^2} \;=\; 1.$$
   Thus, the solution of V23 Eq. (6) for an arbitrary distribution $\rho_{\text{src}}(\mathbf{y})$ is a **positive-definite Yukawa smoothing**:
   $$\phi_L^{\text{Yuk}}(\mathbf{x}) \;=\; \alpha_Y\,\gamma_L\int_{\mathbb{R}^3} G_L(\mathbf{x}-\mathbf{y})\,\rho_{\text{src}}(\mathbf{y})\,d^3y \;=\; \alpha_Y\,\frac{G}{c^2}\int_{\mathbb{R}^3}\frac{e^{-|\mathbf{x}-\mathbf{y}|/L}}{4\pi\,|\mathbf{x}-\mathbf{y}|}\,\rho_{\text{src}}(\mathbf{y})\,d^3y.$$
2. **Homogeneous / Slowly Varying Density Limit ($\rho_{\text{src}} = \text{const.}$):**
   Since $\int_{\mathbb{R}^3} G_L(\mathbf{r})\,d^3r = 1$, in a homogeneous medium:
   $$\phi_L^{\text{Yuk}} \;=\; \alpha_Y\,\gamma_L\,\rho_{\text{src}}.$$
   Choosing the **unit-normalized kernel** ($\alpha_Y = 1 \iff \alpha_{\text{std}} = \frac{1}{4\pi} \approx 0.079577$):
   $$\phi_L^{\text{Yuk}} \;=\; \gamma_L\,\rho_{\text{src}} \quad\implies\quad f_L \;=\; e^{-\phi_L^{\text{Yuk}}} \;=\; \exp(-\gamma_L\rho_{\text{src}}),$$
   which recovers precisely the TAL exponential law examined in Step 1 (`step1_dilaton_form.md`)!
3. **Point-Mass / Compact-Source Limit ($R_* \ll r$):**
   Outside a compact source of mass $M = \int \rho_{\text{src}}\,d^3y$ ($r > R_*$), the solution V23 Eq. (7) simplifies to a classical fifth-force Yukawa potential:
   $$\phi_L^{\text{Yuk}}(r) \;\propto\; \alpha_Y\,\frac{G M}{c^2 r}\,e^{-r/L}.$$

---

## 4. Checks

### 4.1. Dimensional Analysis of V23 Eqs. (6) and (7)
- **Left-hand side of V23 Eq. (6):** $\phi_L$ is dimensionless ($[\phi_L] = 1$), the dimension of the d'Alembert/Laplace operator is $[\Box_g] = [\nabla^2] = \text{m}^{-2}$, and $[L^{-2}] = \text{m}^{-2}$, therefore:
  $$\left[\left(-\Box_g + \frac{1}{L^2}\right)\phi_L\right] \;=\; \text{m}^{-2}.$$
- **Right-hand side of V23 Eq. (6):** $[G/c^2] = \text{m}\,\text{kg}^{-1}$, $[\alpha_{\text{std}}] = 1$, $[\rho_{\text{src}}] = \text{kg}\,\text{m}^{-3}$, therefore:
  $$\left[\frac{4\pi G}{c^2}\alpha_{\text{std}}\,\rho_{\text{src}}\right] \;=\; \left(\text{m}\,\text{kg}^{-1}\right)\cdot \left(\text{kg}\,\text{m}^{-3}\right) \;=\; \text{m}^{-2} \quad \checkmark$$
- **Central integral of V23 Eq. (7) ($\phi_L^{\text{Yuk}}(0)$):**
  $$\left[\alpha_Y\,\frac{G}{c^2}\int_0^{R_*} r'\,e^{-r'/L}\rho_{\text{src}}(r')\,dr'\right] \;=\; 1 \cdot \left(\text{m}\,\text{kg}^{-1}\right)\cdot \text{m}\cdot 1\cdot \left(\text{kg}\,\text{m}^{-3}\right)\cdot \text{m} \;=\; \text{m}^3\,\text{kg}^{-1}\cdot \text{kg}\,\text{m}^{-3} \;=\; 1 \quad (\text{dimensionless}) \quad \checkmark$$

### 4.2. Agreement of SymPy Verifications with the V23 Package (`conservation_Residuals.json`)
In V23, `covariant_closure_50aa2db0.py` (lines 126–148) runs three symbolic checks on this sector, each yielding an exact residual of $0$:
1. `sym_yukawa_center_limit_residual = 0`: $\lim_{r\to 0} \frac{L}{r}\sinh(r/L) = 1$, so $\phi_L^{\text{Yuk}}(0)$ contains no extra $L$-multiplier $\checkmark$.
2. `sym_yukawa_alphaY_eq_4pi_alphastd_residual = 0`: $\alpha_{\text{std}} = \frac{1}{4\pi} \implies \alpha_Y = 4\pi\alpha_{\text{std}} = 1 \implies \phi_L = \gamma_L\rho$ in a homogeneous medium $\checkmark$.
3. `sym_yukawa_el_source_factor_residual = 0` and `sym_yukawa_el_matter_pressure_residual = 0`: variation of the matter term $\mathcal{L} \supset -e^{-\phi_L}\varepsilon(n)$ yields the source term $-\alpha_Y f_L\varepsilon$, the effective pressure $P_{\text{matter}} = f_L P$, and the modified hydrostatic equilibrium $h' = -\nu' + \phi_L'$ $\checkmark$.
