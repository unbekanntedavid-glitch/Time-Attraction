# Step 13: Metric Variation of the Curved-Spacetime Bilateral Top-Hat Window and Volume Element ($\delta[\sqrt{h}\,K_L(x,y)]/\delta g^{\mu\nu}$), Synge's World Function, and the $\delta_D(L - d_\Sigma)$ Boundary Surface Stress

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step13_curved_window_metric_variation_synge_world_function_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step12_second_variation_coupled_pulsation_hessian_en.md` (`step12_second_variation_coupled_pulsation_hessian.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Abstract, line 28 and Section 1.1, lines 38–55, Eqs. 2–3; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p2.tex` (Chapters 4 and 7) and `hu_p3.tex` (Chapters 8 and 11; DOI: `10.5281/zenodo.23604423`)

---

## 1. Verbatim Citations from the V23 (`50aa2db0`) Manuscript

In Step 10 (`step10_bilateral_kernel_reciprocity_eft_synthesis.md`), we derived the unweighted (`bilateral_eom`) and lapse-weighted (`bilateral_lapse_weighted`) reciprocal equilibrium of V23 Formulation I **at fixed kernel and metric**. However, the abstract (line 28) and Section 1.1 (lines 38–54) of the V23 manuscript explicitly state that one open mathematical step remains for the covariant Euler–Lagrange closure of Formulation I (Eq. 2): **the variation of the curved-space window function $K_L(x,y)$ and the volume element with respect to the metric $g_{\mu\nu}$**:

1. **Covariant Action of V23 Formulation I and the $\Sigma_u$ Hypersurface Kernel (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, Section 1.1, lines 38–46, Eqs. 2–3):**
   > *"In Formulation I, a Lagrange multiplier field $\lambda_L(x)$ couples the auxiliary scalar $\phi_L(x)$ to the spatial top-hat average taken on the hypersurface $\Sigma_u$ orthogonal to the fluid four-velocity $u^\mu$ ($h_{\mu\nu} = g_{\mu\nu} + u_\mu u_\nu/c^2$):"*
   $$\begin{equation}
     S_{\text{I},X}[\phi_L, \lambda_L, g_{\mu\nu}] = -\frac{c^4}{8\pi G \gamma_L}\int d^4x\sqrt{-g}\,\lambda_L(x)\!\left[\phi_L(x) - \gamma_L \int_{\Sigma_u} d^4y\sqrt{-g(y)}\,K_L(x,y)\,\rho_{\text{src}}(y)\right], \tag{V23 Eq. 2}
   \end{equation}$$
   > *"where $K_L(x,y) = \delta(\tau_u(x)-\tau_u(y))\,\Theta(L - d_\Sigma(x,y))/V_L(x)$ is the positive-definite 3D spherical kernel of radius $L$ and proper volume $V_L = \frac{4}{3}\pi L^3$."*

2. **Curved-Space Window and Volume-Element Variation Flagged as Open in the V23 Manuscript (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, lines 28 and 54):**
   > - **Abstract (line 28):** *"[T]his is not a dynamical proof of **curved-space window variation**, which remains an open question."*
   > - **Section 1.1 (line 54):** *"Varying the lapse-weighted integral $\int dV\,N f_L\varepsilon$ with $N=e^\nu$ (at fixed kernel and metric) gives $V_{\text{eff},N}(x)=f_L-\frac{\gamma_L}{N(x)}\!\int dV_y\,N(y)K(y,x)f_L\varepsilon$, which I solve in the lapse equation with $P_{\text{matter}}=f_LP$ (`bilateral_lapse_weighted` rows; two-cell test $N=(1;0.8)$, $\varepsilon=(1;2)$, $\gamma=0.1$: $(0.74882;0.72084)$, unweighted $0.73160$). **Covariant closure still requires metric variation of the kernel and volume element; neither row is the EL solution of action (2).**"*

---

## 2. Mathematical Derivation: Metric Variation of the Curved-Space Volume Element and the Geodesic Window $d_\Sigma(x,y)$

### 2.1. Variation of the Hypersurface Integration Measure on $\Sigma_u$ and the Proper Ball Volume $V_L(x)$
In a static spacetime ($ds^2 = -N(x)^2 c^2 dt^2 + h_{ij}(x)dx^i dx^j$, where $N(x) = e^{\nu(x)}$ is the lapse function and $h_{ij}(x)$ is the induced 3-metric on the spatial hypersurface $\Sigma_u$, with $\sqrt{-g(x)} = N(x)\sqrt{h(x)}$), the kernel integral restricted to the simultaneous hypersurface $\tau_u(x) = \tau_u(y)$ is:
$$\bar{\rho}_L(x) \;=\; \int_{\Sigma_u} d^3y\,\sqrt{h(y)}\,\mathcal{K}_L(x,y)\,\varepsilon(y), \qquad \mathcal{K}_L(x,y) \;\equiv\; \frac{\Theta\!\bigl(L - d_\Sigma(x,y)\bigr)}{V_L(x)},$$
where $d_\Sigma(x,y)$ is the minimal spatial geodesic arc length connecting points $x$ and $y$ on $\Sigma_u$, and the proper volume of the geodesic ball $B_\Sigma(x,L) = \{y \in \Sigma_u : d_\Sigma(x,y) \le L\}$ is:
$$V_L(x) \;=\; \int_{\Sigma_u} d^3y\,\sqrt{h(y)}\,\Theta\!\bigl(L - d_\Sigma(x,y)\bigr).$$
By V23 Eq. (3), the Lagrange multiplier evaluated on the solution ($\phi_L = \gamma_L\bar{\rho}_L$, $\mathcal{L}_{\text{eff}} = -f_L(\phi_L)\varepsilon$) is:
$$\lambda_L(x) \;=\; \frac{8\pi G\gamma_L}{c^4}\,\frac{\partial\bigl(-e^{-\phi_L(x)}\varepsilon(x)\bigr)}{\partial\phi_L(x)} \;=\; \frac{8\pi G\gamma_L}{c^4}\,f_L(x)\,\varepsilon(x).$$
When the action $S_{\text{I}} = \int dt\int d^3x\,N(x)\sqrt{h(x)}\,\bigl[-f_L(x)\varepsilon(x)\bigr]$ is varied with respect to the spatial metric $h_{ij}(z)$ (at fixed $N(x)$ and $\varepsilon(x)$), three distinct metric contributions appear in the variation $\delta\phi_L(x) = \gamma_L\,\delta\bar{\rho}_L(x)$:
$$\delta\bar{\rho}_L(x) \;=\; \int_{\Sigma_u} d^3y\,\delta\Bigl[\sqrt{h(y)}\,\mathcal{K}_L(x,y)\Bigr]\varepsilon(y) \;=\; \delta\bar{\rho}_L^{(\text{vol})}(x) \;+\; \delta\bar{\rho}_L^{(\text{norm})}(x) \;+\; \delta\bar{\rho}_L^{(\text{bdry})}(x).$$

1. **Direct Variation of the Volume Element $\sqrt{h(y)}$ ($\delta\bar{\rho}_L^{(\text{vol})}$):**
   Since $\delta\sqrt{h(y)} = \frac{1}{2}\sqrt{h(y)}\,h^{ij}(y)\,\delta h_{ij}(y)$:
   $$\delta\bar{\rho}_L^{(\text{vol})}(x) \;=\; \frac{1}{2}\int_{\Sigma_u} d^3y\,\sqrt{h(y)}\,\mathcal{K}_L(x,y)\,\varepsilon(y)\,h^{ij}(y)\,\delta h_{ij}(y).$$
2. **Variation of the Normalization Volume $V_L(x)$ ($\delta\bar{\rho}_L^{(\text{norm})}$):**
   Since $\delta\bigl(V_L(x)^{-1}\bigr) = -V_L(x)^{-2}\,\delta V_L(x)$, the volume-element variation from the normalization contributes a subtraction weighted by the local mean $\bar{\rho}_L(x)$:
   $$\delta\bar{\rho}_L^{(\text{norm})}\big|_{\text{vol}}(x) \;=\; -\bar{\rho}_L(x)\,\frac{1}{2}\int_{\Sigma_u} d^3y\,\sqrt{h(y)}\,\mathcal{K}_L(x,y)\,h^{ij}(y)\,\delta h_{ij}(y).$$
   The sum of the two volume-element terms is remarkable in its own right: **in a homogeneous medium ($\varepsilon(y) = \bar{\rho}_L(x)$), the two volume-element variations ($\delta\bar{\rho}_L^{(\text{vol})} + \delta\bar{\rho}_L^{(\text{norm})}\big|_{\text{vol}}$) cancel each other identically!**

### 2.2. Variation of the Geodesic Distance $\delta d_\Sigma(x,y)$ (Synge's World Function) and the Boundary Surface Stress $\delta_D(L - d_\Sigma)$
Let $\gamma_{xy}: [0, d_\Sigma(x,y)] \to \Sigma_u$, $l \mapsto z^k(l)$ denote the minimal spatial geodesic curve connecting points $x$ and $y$, parameterized by proper arc length $l$, with unit tangent vector $\hat{t}^k(l) = \frac{dz^k}{dl}$ ($h_{ij}\hat{t}^i\hat{t}^j = 1$). Because of the geodesic equation, the first variation of Synge's spatial world function ($\sigma_\Sigma(x,y) = \frac{1}{2}d_\Sigma(x,y)^2$) at fixed endpoints ($x, y$) is independent of the path displacement $\delta z^k(l)$ of the geodesic, and reduces purely to the line integral of the metric perturbation $\delta h_{ij}$ along the geodesic:
$$\boxed{\;\delta d_\Sigma(x,y) \;=\; \frac{1}{2}\int_0^{d_\Sigma(x,y)} dl\,\hat{t}^i(l)\,\hat{t}^j(l)\,\delta h_{ij}\bigl(\gamma_{xy}(l)\bigr) \;=\; \frac{1}{2}\int_{\Sigma_u} d^3z\,\delta^{(3)}\!\bigl(z - \gamma_{xy}(l)\bigr)\,\hat{t}^i(l)\hat{t}^j(l)\,\delta h_{ij}(z)\,dl.\;}$$
Because the derivative of the Heaviside step function is the Dirac delta ($\frac{\partial}{\partial d_\Sigma}\Theta(L - d_\Sigma) = -\delta_D(L - d_\Sigma)$), the displacement of the window boundary generates a **spherical surface ($d_\Sigma(x,y) = L$) boundary stress**:
$$\boxed{\;\delta\bar{\rho}_L^{(\text{bdry})}(x) \;=\; -\frac{1}{2 V_L(x)}\oint_{d_\Sigma(x,y)=L} dA_y\,\bigl[\varepsilon(y) - \bar{\rho}_L(x)\bigr]\int_0^L dl\,\hat{t}^i(l)\hat{t}^j(l)\,\delta h_{ij}\bigl(\gamma_{xy}(l)\bigr),\;}$$
where $dA_y$ is the proper surface element of the geodesic sphere $d_\Sigma(x,y) = L$, and the factor $\varepsilon(y) - \bar{\rho}_L(x)$ arises from combining the boundary variation of $\Theta(L - d_\Sigma)$ in the numerator with that of $V_L(x)$ in the denominator.

### 2.3. Evaluation Inside a Spherically Symmetric Neutron Star ($h_{rr} = e^{2\lambda(r)} = B(r)^{-1}$) and Connection to Formulation II (Yukawa)
In a spherically symmetric metric ($ds_\Sigma^2 = e^{2\lambda(r)}dr^2 + r^2 d\Omega^2$, where $e^{-2\lambda(r)} = B(r) = 1 - 2m(r)/r$), the radial metric variation $\delta h_{rr}(r) = 2 e^{2\lambda(r)}\delta\lambda(r)$ in the $(rr)$ Einstein equation ($G^r_{\;r} = 8\pi T^r_{\;r,\text{tot}}$) adds the following nonlocal radial window pressure ($P_{r,\text{window}}(r)$) to the `bilateral_lapse_weighted` equations:
$$P_{r,\text{window}}(r) \;=\; -\frac{1}{4\pi r^2 N(r) e^{\lambda(r)}}\,\frac{\delta S_{\text{I},X}}{\delta\lambda(r)}.$$
This expression illuminates two fundamental statements of the V23 manuscript:
1. **Why does `bilateral_lapse_weighted` computed with a flat-space kernel differ from the EL solution of action (2)?**
   Because the flat-space kernel ($d_{\text{flat}}(x,y) = |\mathbf{x} - \mathbf{y}|$) omits the curved volume-element factor $e^{\lambda(r)} = (1 - 2m/r)^{-1/2}$ and the geodesic boundary line integral $\delta d_\Sigma(x,y)$.
2. **Why does Formulation II (Helmholtz–Yukawa, `yukawa_action_el`) eliminate the distributional surface boundary stress $\delta_D(L - d_\Sigma)$?**
   When the sharp top-hat window $\Theta(L - d_\Sigma)$ is replaced by the smooth 3D Helmholtz–Yukawa resolvent ($(-\Box_\Sigma + L^{-2})^{-1}$), the nonlocal boundary stress arising from geodesic line integrals transforms into the **local, continuously differentiable scalar pressure** derived in Step 7 (`step7_neutron_star_tov_stability.md`):
   $$p_{r,\phi}(r) \;=\; \frac{1}{2\alpha_Y}\left(e^{-2\lambda(r)}\bigl[\phi_L'(r)\bigr]^2 - \frac{\phi_L(r)^2}{L^2}\right),$$
   in which the variation of the factor $g^{rr} = e^{-2\lambda(r)} = B(r)$ ($\delta g^{rr} = -2e^{-2\lambda}\delta\lambda$) automatically and algebraically incorporates all metric back-reactions of the curved-space distance measure!

---

## 3. Checks

1. **Homogeneous-Density Limit ($\varepsilon(y) \equiv \varepsilon_0 = \text{const.}$):**
   - If $\varepsilon(y) = \varepsilon_0$ is constant inside the ball $d_\Sigma(x,y) \le L$, then $\bar{\rho}_L(x) = \varepsilon_0$, so both the volume-element variation ($\varepsilon(y) - \bar{\rho}_L(x) = 0$) and the $\delta_D(L - d_\Sigma)$ surface boundary stress ($\varepsilon(y)\big|_{d_\Sigma=L} - \bar{\rho}_L(x) = 0$) **vanish identically**: $\delta\bar{\rho}_L(x) \equiv 0$ $\checkmark$.
2. **Flat-Spacetime ($\lambda(r) \to 0$, $N(r) \to 1$) Limit:**
   - $h_{ij} \to \delta_{ij} \implies d_\Sigma(x,y) \to |\mathbf{x} - \mathbf{y}|$, $V_L(x) \to \frac{4}{3}\pi L^3$, and the variational equation recovers the `bilateral_eom` / `bilateral_lapse_weighted` expressions derived in Step 10 $\checkmark$.
