# Step 10: The Nonlocal Bilateral Kernel Action (Formulation I), Lapse-Weighted Reciprocity Potential ($V_{\text{eff},N}$), and the Moment/EFT Bridge Between Formulations I and II (Yukawa Scalar-Tensor)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step10_bilateral_kernel_reciprocity_eft_synthesis_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step9_vhe_photon_propagation_hess_dispersion_en.md` (`step9_vhe_photon_propagation_hess_dispersion.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Sections 1.1 and 1.2, lines 37–66, Eqs. 2–3, Table 1; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p3.tex` (Chapters 11 and 12, lines 94–120; DOI: `10.5281/zenodo.23604423`)

---

## 1. Verbatim Citations from the V23 (`50aa2db0`) Manuscript and the V22 Monograph

In Steps 1–9, we developed Formulation II (the local screened Helmholtz–Yukawa scalar–tensor action) together with its observational and laboratory consequences. In this step, we derive **Formulation I (the nonlocal bilateral kernel action with auxiliary fields $\phi_L, \lambda_L$, unweighted `bilateral_eom`, and lapse-weighted `bilateral_lapse_weighted` reciprocal equilibrium)** as well as the **exact moment expansion connecting Formulation I (3D spherical top-hat, Branch I-A/I-B) and Formulation II (Yukawa scalar–tensor, Branch II-Fluid, $\alpha_Y = 1$ vs. $\alpha_Y = 1/3$)**, strictly recording that the second-order gradient moments of the two kernels ($\frac{1}{10}L^2\nabla^2$ vs. $L^2\nabla^2$) differ from each other.

1. **Formulation I of V23: Nonlocal Positive-Definite Kernel Action and Lagrange Multiplier (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, Section 1.1, lines 38–46, Eqs. 2–3):**
   > *"In Formulation I, a Lagrange multiplier field $\lambda_L(x)$ couples the auxiliary scalar $\phi_L(x)$ to the spatial top-hat average taken on the hypersurface $\Sigma_u$ orthogonal to the fluid four-velocity $u^\mu$ ($h_{\mu\nu} = g_{\mu\nu} + u_\mu u_\nu/c^2$):"*
   $$\begin{equation}
     S_{\text{I},X}[\phi_L, \lambda_L, g_{\mu\nu}] = -\frac{c^4}{8\pi G \gamma_L}\int d^4x\sqrt{-g}\,\lambda_L(x)\!\left[\phi_L(x) - \gamma_L \int_{\Sigma_u} d^4y\sqrt{-g(y)}\,K_L(x,y)\,\rho_{\text{src}}(y)\right], \tag{V23 Eq. 2}
   \end{equation}$$
   > *"where $K_L(x,y) = \delta(\tau_u(x)-\tau_u(y))\,\Theta(L - d_\Sigma(x,y))/V_L(x)$ is the positive-definite 3D spherical kernel of radius $L$ and proper volume $V_L = \frac{4}{3}\pi L^3$. Variation with respect to $\lambda_L(x)$, $\phi_L(x)$, and $\rho_{\text{src}}(x)$ yields the following Euler–Lagrange closure:"*
   $$\begin{equation}
     \frac{\delta S}{\delta \lambda_L(x)} = 0 \implies \phi_L(x) = \gamma_L\bar{\rho}_L(x),\quad f_L(x) = e^{-\phi_L(x)}; \qquad
     \frac{\delta S}{\delta \phi_L(x)} = 0 \implies \lambda_L(x) = \frac{8\pi G\gamma_L}{c^4}\frac{\partial \mathcal{L}_{\text{eff}}}{\partial \phi_L}(x). \tag{V23 Eq. 3}
   \end{equation}$$

2. **Bilateral Kernel Reciprocity ($x \leftrightarrow y$), Lapse-Weighted $V_{\text{eff},N}$, and the Two-Cell Variational Test (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 54):**
   > *"If $\lambda_L$ is absorbed into the fluid in the form $\mathcal{L}_{\text{fluid}}(n,\phi_L) = -f_L(\phi_L)\,\varepsilon(n)$ and $\varepsilon$ is varied both explicitly and inside the kernel, then $V_{\text{eff}}(x) = f_L(x) - \gamma_L\!\int K_L(y,x) f_L(y)\varepsilon(y)\,d^3y$ and the modified equilibrium is $dh/dr = -A(r) - dV_{\text{recip}}/dr$, $V_{\text{recip}} = \ln V_{\text{eff}}$. I solve this as a self-consistent fixed point (`bilateral_eom` rows in Table 1); it differs from the algebraically closed V21 mass mode. Varying the lapse-weighted integral $\int dV\,N f_L\varepsilon$ with $N=e^\nu$ (at fixed kernel and metric) gives $V_{\text{eff},N}(x)=f_L-\frac{\gamma_L}{N(x)}\!\int dV_y\,N(y)K(y,x)f_L\varepsilon$, which I solve in the lapse equation with $P_{\text{matter}}=f_LP$ (`bilateral_lapse_weighted` rows; two-cell test $N=(1;0.8)$, $\varepsilon=(1;2)$, $\gamma=0.1$: $(0.74882;0.72084)$, unweighted $0.73160$). Covariant closure still requires metric variation of the kernel and volume element; neither row is the EL solution of action (2)."*

3. **The $\alpha_Y = 1$ (Unit-Normalized) and $\alpha_Y = 1/3$ (Volume-Matched) Yukawa Normalizations (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 65):**
   > *"The unit-normalized kernel corresponds to $\alpha_Y = 1$, i.e., $\alpha_{\text{std}} = 1/(4\pi)\approx 0.0796$, which recovers the homogeneous top-hat limit ($\phi_L = \gamma_L\rho_{\text{src}}$); $\alpha_Y = 1/3$ ($\alpha_{\text{std}} = 1/(12\pi)$) is the volume-matched variant."*

---

## 2. Mathematical Derivation

### 2.1. Derivation of Bilateral Kernel Reciprocity ($x \leftrightarrow y$) and Lapse-Weighted $V_{\text{eff},N}(x)$
Consider the static energy functional in which the matter Lagrangian $\mathcal{L}_{\text{fluid}}(n, \phi_L) = -f_L(\phi_L)\,\varepsilon(n)$ ($f_L(x) = \exp[-\gamma_L\int dV_y\,K_L(x,y)\varepsilon(y)]$) appears weighted by the lapse factor $N(x) = e^{\nu(x)}$ (with fixed background metric and fixed kernel $K_L(x,y)$):
$$E_{\text{matter}}[\varepsilon] \;=\; \int dV_x\,N(x)\,f_L(x)\,\varepsilon(x), \qquad \phi_L(x) \;=\; \gamma_L\int dV_y\,K_L(x,y)\,\varepsilon(y).$$
Let us vary $E_{\text{matter}}[\varepsilon]$ with respect to the energy density $\varepsilon(x)$! Because $\varepsilon$ appears in two places — once as an explicit factor at point $x$ and once in the exponent of $f_L(y) = \exp[-\gamma_L\int dV_x\,K_L(y,x)\varepsilon(x)]$ at every point $y$ whose smoothing ball contains $x$ — the functional derivative by the chain rule and by exchanging the dummy integration variables $x \leftrightarrow y$ (*bilateral kernel reciprocity*) is:
$$\delta E_{\text{matter}} \;=\; \int dV_x\,N(x)\,f_L(x)\,\delta\varepsilon(x) \;+\; \int dV_y\,N(y)\,\varepsilon(y)\,\frac{\partial f_L(y)}{\partial\phi_L(y)}\,\gamma_L\!\int dV_x\,K_L(y,x)\,\delta\varepsilon(x).$$
Since $\frac{\partial f_L(y)}{\partial\phi_L(y)} = -f_L(y)$, exchanging the integration order of $dV_y$ and $dV_x$ in the second integral and factoring out $N(x)\,\delta\varepsilon(x)$ yields:
$$\delta E_{\text{matter}} \;=\; \int dV_x\,N(x)\,\underbrace{\left[f_L(x) \;-\; \frac{\gamma_L}{N(x)}\int dV_y\,N(y)\,K_L(y,x)\,f_L(y)\,\varepsilon(y)\right]}_{\equiv\; V_{\text{eff},N}(x)}\delta\varepsilon(x).$$
This derives the two bilateral potentials appearing in Section 1.1 (line 54) of V23:
1. **Lapse-weighted bilateral potential (`bilateral_lapse_weighted`):**
   $$\boxed{\;V_{\text{eff},N}(x) \;=\; f_L(x) \;-\; \frac{\gamma_L}{N(x)}\int dV_y\,N(y)\,K_L(y,x)\,f_L(y)\,\varepsilon(y).\;}$$
2. **Unweighted bilateral potential (`bilateral_eom`, $N(y)/N(x) \to 1$):**
   $$\boxed{\;V_{\text{eff}}(x) \;=\; f_L(x) \;-\; \gamma_L\int dV_y\,K_L(y,x)\,f_L(y)\,\varepsilon(y).\;}$$
In the hydrostatic Euler equation ($\delta\varepsilon/\delta n = (\varepsilon + P)/n = \mu e^h$), the condition of constant chemical potential ($N(r)\,V_{\text{eff}}(r)\,e^{h(r)} = \text{const.}$) after logarithmic differentiation directly gives the modified equilibrium equation:
$$\frac{dh}{dr} \;=\; -\frac{d\nu}{dr} \;-\; \frac{d V_{\text{recip}}}{dr}, \qquad V_{\text{recip}}(r) \;\equiv\; \ln V_{\text{eff}}(r).$$
Since $\gamma_L\int K_L f_L\varepsilon\,dV_y > 0$, the bilateral reciprocal term **reduces** $V_{\text{eff}}$ in the dense stellar core ($V_{\text{eff}}(0) < f_L(0)$), producing an inward additional compressive force that explains why the `bilateral_eom` and `bilateral_lapse_weighted` branches yield a lower maximum mass than GR ($M_{\max}^{\text{SLy}} = 1.93214\,M_\odot$ and $1.90217\,M_\odot$ at $L = 10\text{ km}$, compared to $2.09822\,M_\odot$ in the V21 `mass` mode; V23 Table 1).

### 2.2. Exact Moment-Expansion Comparison: Where Do the 3D Top-Hat (Branch I) and Yukawa (Branch II) Formulations Agree and Where Do They Differ?
Expanding the source density $\rho_{\text{src}}(\mathbf{x} + \mathbf{u})$ in a Taylor series around point $\mathbf{x}$ inside the 3D spherical top-hat average ($|\mathbf{u}| \le L$) gives:
$$\bar{\rho}_L(\mathbf{x}) \;=\; \frac{3}{4\pi L^3}\int_{|\mathbf{u}|\le L}\!\left[\rho_{\text{src}}(\mathbf{x}) + u_i\partial_i\rho_{\text{src}}(\mathbf{x}) + \frac{1}{2}u_i u_j\partial_i\partial_j\rho_{\text{src}}(\mathbf{x}) + \dots\right]d^3u.$$
By spherical symmetry, $\int u_i\,d^3u = 0$ and $\int u_i u_j\,d^3u = \frac{1}{3}\delta_{ij}\int_0^L 4\pi u^4\,du = \frac{4\pi}{15}L^5\,\delta_{ij}$, so the gradient expansion of the 3D top-hat operator is:
$$\boxed{\;\phi_L^{\text{top-hat}}(\mathbf{x}) \;=\; \gamma_L\,\bar{\rho}_L(\mathbf{x}) \;=\; \gamma_L\left(1 \;+\; \frac{L^2}{10}\nabla^2 \;+\; \mathcal{O}(L^4\nabla^4)\right)\rho_{\text{src}}(\mathbf{x}).\;}$$
On the other hand, the formal inverse of the Helmholtz–Yukawa equation of Formulation II ($(-L^2\nabla^2 + 1)\phi_L^{\text{Yuk}} = \alpha_Y\gamma_L\rho_{\text{src}}$) for small gradients is:
$$\boxed{\;\phi_L^{\text{Yuk}}(\mathbf{x}) \;=\; \alpha_Y\,\gamma_L\,\bigl(1 - L^2\nabla^2\bigr)^{-1}\rho_{\text{src}}(\mathbf{x}) \;=\; \alpha_Y\,\gamma_L\left(1 \;+\; L^2\nabla^2 \;+\; \mathcal{O}(L^4\nabla^4)\right)\rho_{\text{src}}(\mathbf{x}).\;}$$

**Critical Theoretical Distinction (Point 5):**
1. **Zeroth-Moment (Homogeneous $\nabla^2\rho_{\text{src}} \to 0$) Limit:**
   Under the normalization $\alpha_Y = 1$ ($\alpha_{\text{std}} = \frac{1}{4\pi} \approx 0.07958$), the 3D top-hat and the Yukawa Green's function agree **exclusively in the zeroth moment** ($\int K_L(\mathbf{u})\,d^3u = \int G_L(\mathbf{u})\,d^3u = 1$, i.e., in a homogeneous medium: $\phi_L^{\text{Yuk}} = \gamma_L\rho_{\text{src}} = \phi_L^{\text{top-hat}}$)!
2. **Difference in the Second-Order Gradient Moment ($\frac{1}{10}L^2\nabla^2$ vs. $L^2\nabla^2$) and Compact Support:**
   In an inhomogeneous medium ($\nabla^2\rho_{\text{src}} \neq 0$), the two models are **not** the same effective field theory:
   - at fixed $L$ and $\alpha_Y = 1$, the second-order $\nabla^2\rho_{\text{src}}$ coefficient differs by a factor of 10 ($\frac{1}{10}L^2$ for the 3D top-hat kernel vs. $1\cdot L^2$ for the Yukawa Green's function; matching the second moment would require $\ell_Y = L/\sqrt{10}$);
   - the 3D top-hat kernel has compact support ($r \le L$, so the exterior stellar halo terminates sharply at $r = R_* + L$), whereas the Yukawa Green's function generates an exponentially decaying ($\propto e^{-r/L}/r$) curved exterior scalar tail over an infinite domain and an ADM mass contribution $\Delta M_{\phi,\text{ext}}^{\text{curved}}$;
   - the $\alpha_Y = \frac{1}{3}$ ($\alpha_{\text{std}} = \frac{1}{12\pi} \approx 0.02653$) variant matches the ratio of the spherical volume $\frac{4}{3}\pi L^3$ to the Yukawa normalization volume $4\pi L^3$ ($\frac{4\pi/3}{4\pi} = \frac{1}{3}$) at the stellar center ($r \to 0$).

---

## 3. Checks

1. **Exact Reproduction of the V23 Two-Cell Variational Test (`bilateral_lapse_weighted` vs. `bilateral_eom`):**
   - Consider two cells of equal volume ($\Delta V_1 = \Delta V_2 = 1/2$, symmetric overlapping kernel $K_{ij} = 1$, so $\bar{\varepsilon} = \frac{1}{2}(\varepsilon_1 + \varepsilon_2)$), where according to line 54 of V23, $N = (1.0;\; 0.8)$, $\varepsilon = (1.0;\; 2.0)$, and $\gamma_L = 0.1$:
     - Mean density: $\bar{\varepsilon} = \frac{1 + 2}{2} = 1.5 \implies \phi_L = \gamma_L\bar{\varepsilon} = 0.15 \implies f_L(1) = f_L(2) = e^{-0.15} = 0.86070798$.
     - **Unweighted (`bilateral_eom`):**
       $$V_{\text{eff}}(1) = V_{\text{eff}}(2) = f_L\bigl(1 - \gamma_L\bar{\varepsilon}\bigr) = e^{-0.15}\times (1 - 0.15) = 0.85\,e^{-0.15} = \mathbf{0.73160178} \approx \mathbf{0.73160} \quad \checkmark$$
     - **Lapse-weighted (`bilateral_lapse_weighted`):**
       The lapse-weighted source average is $\overline{N\varepsilon} = \frac{1}{2}(N_1\varepsilon_1 + N_2\varepsilon_2) = \frac{1}{2}(1.0\times 1.0 + 0.8\times 2.0) = \frac{2.6}{2} = 1.3$:
       $$V_{\text{eff},N}(1) = f_L\left(1 - \frac{\gamma_L\,\overline{N\varepsilon}}{N_1}\right) = e^{-0.15}\left(1 - \frac{0.13}{1.0}\right) = 0.87\,e^{-0.15} = \mathbf{0.74881594} \approx \mathbf{0.74882} \quad \checkmark$$
       $$V_{\text{eff},N}(2) = f_L\left(1 - \frac{\gamma_L\,\overline{N\varepsilon}}{N_2}\right) = e^{-0.15}\left(1 - \frac{0.13}{0.8}\right) = 0.8375\,e^{-0.15} = \mathbf{0.72084293} \approx \mathbf{0.72084} \quad \checkmark$$
   - All three values agree digit for digit with the text of the V23 manuscript (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 54).
