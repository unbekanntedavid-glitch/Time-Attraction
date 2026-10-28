# Step 1: Mathematical Equivalence of the Form of the TAL Factor $f_L = \exp(-\gamma_L \bar{\rho}_L)$ and the String-Theory Dilaton Coupling

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step1_dilaton_form_en.md` (English translation of `step1_dilaton_form.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0`): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` & `tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` & `tal_book_HU_2026-10-24.tex` (DOI: `10.5281/zenodo.23604423`)  
- V19 (model 7): 3D spherical top-hat law (DOI: `10.5281/zenodo.23497128`)  
- Local working note: `docs/C_derivation_status.md` (line 15: `Time_Attraction.pdf (v1–v7)`)

---

## 1. Verbatim definitions of the Time Attraction Law (TAL) from the local manuscripts

In the local V19–V23 manuscripts, the positive-definite exponential time-flow / damping factor $f_L(\mathbf{x})$, the 3D spherical top-hat mean density $\bar{\rho}_L(\mathbf{x})$, and the scale-dependent coupling coefficient $\gamma_L$ (together with the coupling constant $C$) appear verbatim in the following form:

### 1.1. The exponential law and the definition of the constants (V22 `hu_p1.tex`, lines 141–144 and 245–253)
- **Exponential time-flow factor:**
  $$f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(\mathbf{x})\bigr), \qquad 0 < f_L(\mathbf{x}) \le 1$$
- **3D spherical top-hat mean density (with smoothing radius $L$ and volume $V_L = \frac{4}{3}\pi L^3$):**
  $$\bar{\rho}_L(\mathbf{x}) \;=\; \frac{3}{4\pi L^3}\int_{|\mathbf{y}-\mathbf{x}|\le L}\rho(\mathbf{y})\,d^3y \;\ge\; 0$$
- **Fundamental gravitational–temporal constants (in SI units, CODATA $G = 6{,}67430\times 10^{-11}\;\text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}$, $c = 2{,}99792458\times 10^{8}\;\text{m}\,\text{s}^{-1}$):**
  $$C \;\equiv\; \frac{4\pi G}{c^2} \;=\; 9{,}33306\times 10^{-27}\;\text{m}\,\text{kg}^{-1}, \qquad \frac{G}{c^2} \;=\; \frac{C}{4\pi} \;=\; 7{,}425915\times 10^{-28}\;\text{m}\,\text{kg}^{-1}$$
  $$\gamma_L \;=\; \frac{G\,L^2}{c^2} \;=\; \frac{C\,L^2}{4\pi} \;=\; 7{,}425915\times 10^{-28}\,L^2\;\text{m}^3\,\text{kg}^{-1}$$

### 1.2. The covariant auxiliary scalar and the matter-Lagrangian coupling in V23 (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex`, lines 32–46 and 54; `C_derivation_status.md`, line 15)
- The early `Time_Attraction.pdf` (v1–v7, cited in line 15 of `docs/C_derivation_status.md`) and the V23 `50aa2db0` manuscript (Section 1, Eq. 1–3 and line 54) weight the matter action by the exponential factor $f_L$:
  $$S \;=\; \int d^4x\sqrt{-g}\left[\frac{c^4}{16\pi G}R + \mathcal{L}_X(g_{\mu\nu}, \phi_L, \lambda_L) + \mathcal{L}_{\text{matter}}\right], \qquad \mathcal{L}_{\text{fluid}}(n,\phi_L) \;=\; -f_L(\phi_L)\,\varepsilon(n)$$
  or, in the formulation of versions 1–7:
  $$S \;=\; \int d^4x\sqrt{-g}\left[\frac{c^4}{16\pi G}R + f_L\,\mathcal{L}_m\right]$$
- In Formulation I of V23 (Eq. 2–3), the dimensionless scalar $\phi_L(x)$ is tied by the Lagrange multiplier $\lambda_L(x)$ to the 3D spherical top-hat average taken on the spatial hypersurface $\Sigma_u$ ($K_L(x,y) = \delta(\tau_u(x)-\tau_u(y))\,\Theta(L - d_\Sigma(x,y))/V_L(x)$):
  $$\frac{\delta S}{\delta \lambda_L(x)} = 0 \quad\implies\quad \phi_L(x) \;=\; \gamma_L\,\bar{\rho}_L(x), \qquad f_L(x) \;=\; e^{-\phi_L(x)} \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(x)\bigr).$$

---

## 2. The low-energy effective action of string theory and the conventions of the dilaton coupling

In the low-energy effective field theory of string theory (bosonic and superstring), reduced to 4 dimensions, the gravitational sector necessarily contains a scalar field, the **dilaton ($\phi$)** (see e.g. Polchinski: *String Theory*, Vol. I, Sec. 3.7; Gasperini & Veneziano, *Phys. Rept.* **373**, 1, 2003; Damour & Polyakov, *Nucl. Phys. B* **423**, 532, 1994).

1. **String frame:**
   At tree level (weighted by the square of the string coupling constant $g_s = e^{\phi}$, i.e. by $g_s^{-2} = e^{-2\phi}$), the effective action of the gravitational and closed/open string sector is:
   $$S_{\text{string}} \;=\; \int d^4x\sqrt{-g_{(s)}}\,e^{-2\phi}\left[\frac{c^4}{16\pi G_0}\Big(R_{(s)} + 4\,g_{(s)}^{\mu\nu}\partial_\mu\phi\,\partial_\nu\phi\Big) - \dots\right] + \int d^4x\sqrt{-g_{(s)}}\,e^{-\eta\phi}\,\mathcal{L}_m\!\left(g_{(s)\mu\nu},\Psi\right),$$
   where:
   - in the closed-string / tree-level sector $\eta = 2$, i.e. the matter Lagrangian is multiplied directly by the exponential factor **$e^{-2\phi}\,\mathcal{L}_m$**;
   - in the open-string / D-brane (Born–Infeld) sector $\eta = 1$, i.e. the coupling factor is **$e^{-\phi}\,\mathcal{L}_m$**.

2. **Einstein frame — exponential matter coupling alongside a canonical $R$ term:**
   In order for the factor in front of the curvature scalar $R$ (as in TAL V23 Eq. (1): $\frac{c^4}{16\pi G}R$) to take a canonical, $\phi$-independent form, we perform the conformal rescaling $g_{(s)\mu\nu} = e^{2\phi} g_{\mu\nu}$ (Damour & Polyakov 1994; Fujii & Maeda: *The Scalar-Tensor Theory of Gravitation*, Cambridge Univ. Press, 2003). The curvature term then becomes purely Einstein–Hilbert in form, while an exponential dilaton multiplier appears in the matter sector:
   $$S_{\text{Einstein}} \;=\; \int d^4x\sqrt{-g}\left[\frac{c^4}{16\pi G}\Big(R - 2\,g^{\mu\nu}\partial_\mu\phi\,\partial_\nu\phi\Big) \;+\; e^{-2\beta_{\text{dil}}\phi}\,\mathcal{L}_m(\Psi, g_{\mu\nu})\right],$$
   (where, depending on the convention and the matter sector, the usual notation of the exponential multiplier is $e^{-2\phi}$, $e^{-\phi}$, or in general scalar–tensor notation $e^{-2\beta_{\text{dil}}\phi}$).

---

## 3. The formal mathematical correspondence (Step 1 Mapping)

Comparing the TAL matter-coupling term ($f_L(\mathbf{x})\,\mathcal{L}_m$, where $f_L(\mathbf{x}) = e^{-\phi_L(\mathbf{x})} = \exp(-\gamma_L\bar{\rho}_L(\mathbf{x}))$) with the exponential matter coupling of the string-theory dilaton ($e^{-2\phi}\,\mathcal{L}_m$ in the $2\phi$ convention, or $e^{-\phi}\,\mathcal{L}_m$ in the $1\phi$ convention), the **mathematical functional form of the two expressions is exactly the same**:

$$\boxed{\;f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(\mathbf{x})\bigr) \;=\; e^{-\phi_L(\mathbf{x})} \quad\Longleftrightarrow\quad e^{-2\phi(\mathbf{x})} \;\;\text{(or } e^{-\phi(\mathbf{x})}\text{)}\;}$$

i.e. the direct correspondence of the dimensionless exponents:

$$\boxed{\;\phi_L(\mathbf{x}) \;\equiv\; \gamma_L\,\bar{\rho}_L(\mathbf{x}) \;=\; \frac{G L^2}{c^2}\,\bar{\rho}_L(\mathbf{x}) \quad\Longleftrightarrow\quad 2\phi(\mathbf{x}) \;\;\text{(in the } e^{-2\phi} \text{ convention), or } \phi(\mathbf{x}) \;\;\text{(in the } e^{-\phi} \text{ convention).}\;}$$

### Open question recorded for later steps (not resolved in this step)
- **The similarity:** In both theories, a negative exponential function of a dimensionless scalar quantity ($e^{-\phi_L}$ and $e^{-2\phi}$, respectively) multiplies the matter Lagrangian.
- **The difference (open question for Step 2 and beyond):** In the basic definition of the Time Attraction Law (and in Formulation I of V23), the exponent $\phi_L(\mathbf{x}) = \gamma_L \bar{\rho}_L(\mathbf{x})$ is determined algebraically / via a nonlocal integral directly by the mass-energy density $\bar{\rho}_L$ averaged over the scale $L$, whereas in string theory the dilaton $\phi$ is an independent dynamical field with its own kinetic term ($\propto (\partial\phi)^2$). The connection between this dynamical kinetic term and the local Helmholtz–Yukawa operator of Formulation II of V23, $(-L^2\Box_g + 1)\phi_L = \alpha_Y \gamma_L \rho_{\text{src}}$, is examined in later steps; we do not touch it here.

---

## 4. Checks

### 4.1. Dimensional analysis (inverse-density dimension of $\gamma_L$)
We check whether the exponent $\gamma_L \bar{\rho}_L$ is dimensionless in SI base units ($\text{m}, \text{kg}, \text{s}$):
- Dimension of Newton's gravitational constant:
  $$[G] \;=\; \text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}$$
- Dimension of the square of the speed of light:
  $$[c^2] \;=\; \text{m}^2\,\text{s}^{-2}$$
- Dimension of the mass–length ratio ($G/c^2$, or $C = 4\pi G/c^2$):
  $$\left[\frac{G}{c^2}\right] \;=\; \frac{\text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}}{\text{m}^2\,\text{s}^{-2}} \;=\; \text{m}\,\text{kg}^{-1}$$
- Dimension of the square of the smoothing scale:
  $$[L^2] \;=\; \text{m}^2$$
- Consequently, the dimension of the coupling coefficient $\gamma_L = \frac{G L^2}{c^2}$:
  $$[\gamma_L] \;=\; \left[\frac{G}{c^2}\right]\cdot [L^2] \;=\; \left(\text{m}\,\text{kg}^{-1}\right)\cdot \text{m}^2 \;=\; \text{m}^3\,\text{kg}^{-1} \;=\; \frac{1}{\text{kg}\,\text{m}^{-3}} \;=\; \frac{1}{[\text{density}]} \quad \checkmark$$
- Since the dimension of the 3D volumetric mass density is $[\bar{\rho}_L] = \text{kg}\,\text{m}^{-3}$, the product is:
  $$[\gamma_L\,\bar{\rho}_L] \;=\; \left(\text{m}^3\,\text{kg}^{-1}\right)\cdot \left(\text{kg}\,\text{m}^{-3}\right) \;=\; 1 \quad (\text{dimensionless}) \quad \checkmark$$

*(Note for the written-out check: the ratio $[G/c^2]$ by itself is $\text{m}\,\text{kg}^{-1}$ [the reciprocal of a linear mass density], and multiplication by the factor $L^2$ [$\text{m}^2$] makes the dimension of $\gamma_L = G L^2/c^2$ equal to $\text{m}^3\,\text{kg}^{-1}$, which is exactly the reciprocal of the 3D volumetric density [$\text{kg}\,\text{m}^{-3}$], consistent with Eq. (3) and Eq. (5) of V22 `hu_p1.tex`.)*

### 4.2. Checks of limits and signs
1. **Vacuum limit ($\bar{\rho}_L \to 0$):**
   When the averaged mass-energy density in the neighbourhood of radius $L$ vanishes ($\bar{\rho}_L \to 0$), the exponent $\phi_L = \gamma_L \bar{\rho}_L \to 0$, so:
   $$\lim_{\bar{\rho}_L \to 0} f_L(\bar{\rho}_L) \;=\; \exp(0) \;=\; 1 \quad \checkmark$$
   This exactly recovers the unmodified standard matter coupling ($\mathcal{L}_{\text{matter}} \to 1\cdot \mathcal{L}_m$), which in the dilaton picture corresponds to the asymptotic vacuum $\phi \to 0$ ($e^{-2\phi} \to 1$).
2. **Positive-density domain ($\bar{\rho}_L > 0$):**
   Since $G > 0$, $c > 0$, $L > 0 \implies \gamma_L > 0$, and the 3D spherical top-hat average is strictly positive ($\bar{\rho}_L > 0$) under the non-negative energy condition ($\rho(\mathbf{y}) \ge 0$), we have $\gamma_L \bar{\rho}_L > 0$, therefore:
   $$0 \;<\; f_L(\bar{\rho}_L) \;=\; \exp(-\gamma_L \bar{\rho}_L) \;<\; 1 \qquad (\forall\,\bar{\rho}_L > 0) \quad \checkmark$$
   and the function is strictly monotonically decreasing:
   $$\frac{df_L}{d\bar{\rho}_L} \;=\; -\gamma_L\,\exp(-\gamma_L \bar{\rho}_L) \;<\; 0 \quad \checkmark$$
