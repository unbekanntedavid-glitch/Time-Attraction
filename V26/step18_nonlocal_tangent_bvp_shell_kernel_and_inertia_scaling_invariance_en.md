# Step 18: Exact Derivation of the V23 Eq. (9) Nonlocal Static Tangent Operator (`eq:coupled_bvp`), the 3D Spherical-Shell Intersection Kernel $K_{ij}$, and the $W \to \kappa W$ Inertia-Scaling Invariance

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness_en.md` (`step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 2, lines 127–153, Eq. 9 and Table 2; DOI: `10.5281/zenodo.23639084`)  
- V23 verification script and results: `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/covariant_closure_50aa2db0.py` (lines 179–191 and 583–695) and `coupled_stability_results.csv`  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (Chapter 7, point 1, line 374 and Chapter 8, point 2, lines 409–419; DOI: `10.5281/zenodo.23604423`)

---

## 1. Verbatim Citations from the V22 and V23 Manuscripts

In Step 17 (`step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md`), we derived V23 Equation (8) (`eq:factored_diff`) and the surface sign reversal $\mathcal{W}_{\text{metric}}(r) \le 0$ of the `gradient`-mode `metric` operator. Immediately following this in Section 2 of the V23 manuscript comes **V23 Equation (9) (`eq:coupled_bvp`)**, the nonlocal kernel perturbation $K_{ij}$ aligned to the piecewise-polytrope transition radii $r_{\text{trans}}$, the $W \to 4W$ inertia-scaling invariance, the first-order tangent estimate $\omega^2 \simeq -F(0)/F'(0)$, and the symmetrized Galerkin matrix $\mathbf{M}_{\text{sym}} = \mathbf{D}^T\mathrm{diag}(w\,\Delta r)\mathbf{D}$:

1. **V23 Eq. (9) Coupled Static Nonlocal Tangent Model and Quasi-Dynamical BVP Test (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, Section 2, lines 127–136, Eq. 9):**
   > *"I linearize the coupled system (around the V21 `dense_profile` background, not around the new branches) in terms of the enthalpy perturbation $\delta h(r) = \delta P/(\varepsilon+P)$ and mass perturbation $\delta m(r)$ on a grid aligned with the exact piecewise-polytrope transition radii $r_{\text{trans}}$, incorporating the nonlocal kernel perturbation operator:"*
   $$\begin{equation}
     \delta f_L(r_i) = -\gamma_L f_L(r_i)\sum_j K_{ij}\,\frac{\varepsilon(r_j)+P(r_j)}{c_s^2(r_j)}\,\delta h(r_j), \qquad
     \begin{cases}
       \dfrac{d(\delta h)}{dr} = -\delta\nu'\!\left[\delta h, \delta m, \delta f_L[\delta h]\right] + \omega_0^2\,e^{\lambda-2\nu}\xi, \\[4pt]
       \dfrac{d(\delta m)}{dr} = 4\pi r^2\!\left[f_{\text{mass}}\frac{\varepsilon+P}{c_s^2}\delta h + \varepsilon\,\delta f_{\text{mass}}[\delta h]\right].
     \end{cases} \tag{V23 Eq. 9}
   \end{equation}$$
   > *"At $\omega_0^2 = 0$, system (9) is the exact tangent operator of the static equilibrium equations with respect to $\log_{10}\rho_{0c}$. As shown in Table 2 and Fig. 3, solving the coupled eigenvalue problem (9) along the SLy $L=10\,\text{km}$ sequences yields a zero crossing within $\sim 7\times 10^{-5}\,\text{dex}$ of $dM/d\rho_c=0$ (the pure GR discretization baseline is $7.12\times 10^{-5}\,\text{dex}$). Because at $\omega_0^2=0$ the problem reduces to $\mathcal{H}_{\text{lin}}[\delta h]=0$, this agreement verifies static nonlocal linearization rather than dynamical inertia: scaling $W\to 4W$ at $\log_{10}\rho_{0c}=15.30427$ yields a ratio $\omega^2_{4W}/\omega^2_W = 0.2500$ ($\omega^2_W=4.897\times 10^{6}\,\text{s}^{-2}$) with the identical shooting function $F(0)=0.023252$, so the zero crossing is inertia-invariant. The raw BVP $W$ matrix is not symmetric; the symmetrized Galerkin form $\mathbf{M}_{\text{sym}} = \mathbf{D}^{T}\mathrm{diag}(w\,\Delta r)\mathbf{D}$ is symmetric with $\min\mathrm{diag}\,\mathbf{M}_{\text{sym}} = 5.68\times 10^{-9}>0$; this is a Gram-matrix property, not a proof of positive-definiteness of the physical kinetic Hessian. The `omega0_sq_sym_tangent_s2` column gives the first-order tangent estimate ($\omega^2\simeq -F(0)/F'(0)$)."*

2. **3D Spherical-Shell Intersection Integral (`tal_book_HU_2026-10-24.tex`, Chapter 7, point 1, line 374):**
   > *"Inside a spherically symmetric star of radius $R_*$, the exact 3D top-hat average $\bar{\rho}_L(r)$ at radial coordinate $r$ can be evaluated analytically over the intersection of the ball $B(\mathbf{x},L)$ with spherical shells $r' \in [|r-L|, r+L]$, yielding $f_L(r) = \exp(-\gamma_L\bar{\rho}_L(r))$."*

---

## 2. Exact Derivation of the 3D Spherical-Shell Intersection Quadrature Kernel $K_{ij}$ and the $r_{\text{trans}}$ Transition Grid

### 2.1. Solid-Angle Fraction of the 3D Spherical Top-Hat Shell Intersection (`build_tophat_matrix`)
Let $\mathbf{x}_i = (0,0,r_i)$ be the query point ($r_i > 0$), and consider a concentric spherical shell of radius $r_j > 0$ and thickness $\Delta r_j$ around the stellar center. By the law of cosines, the squared Euclidean distance of a point $\mathbf{y} = (r_j\sin\theta\cos\varphi,\,r_j\sin\theta\sin\varphi,\,r_j\cos\theta)$ on the shell from $\mathbf{x}_i$ is:
$$|\mathbf{y} - \mathbf{x}_i|^2 \;=\; r_i^2 \;+\; r_j^2 \;-\; 2\,r_i\,r_j\cos\theta.$$
The condition $|\mathbf{y} - \mathbf{x}_i| \le L$ is equivalent to:
$$\cos\theta \;\ge\; \mu_{\min}(r_i, r_j; L) \;\equiv\; \frac{r_i^2 + r_j^2 - L^2}{2\,r_i\,r_j}.$$
The surface fraction of the spherical shell of radius $r_j$ falling inside the ball $B(\mathbf{x}_i, L)$ ($\mathrm{frac}(r_i, r_j; L) = \frac{1}{4\pi}\int_{|\mathbf{y}-\mathbf{x}_i|\le L}\sin\theta\,d\theta\,d\varphi = \frac{1}{2}\int_{\max(-1, \mu_{\min})}^{1} d(\cos\theta)$) splits into three geometric regimes:
1. **Fully Enclosed Shell ($r_i + r_j \le L$):** Here $r_i^2 + r_j^2 - L^2 \le -2 r_i r_j \implies \mu_{\min} \le -1$, so the entire spherical shell lies inside $B(\mathbf{x}_i, L)$: $\mathrm{frac}(r_i, r_j; L) = 1$.
2. **Partially Intersected Shell ($|r_i - r_j| < L < r_i + r_j$):** Here $-1 < \mu_{\min} < +1$, so the solid-angle fraction of the spherical cap is:
   $$\mathrm{frac}(r_i, r_j; L) \;=\; \frac{1 - \mu_{\min}}{2} \;=\; \frac{2 r_i r_j - (r_i^2 + r_j^2 - L^2)}{4\,r_i\,r_j} \;=\; \frac{L^2 - (r_i - r_j)^2}{4\,r_i\,r_j}.$$
3. **Exterior Shell ($|r_i - r_j| \ge L$):** Here $\mu_{\min} \ge +1$, so the intersection is empty: $\mathrm{frac}(r_i, r_j; L) = 0$.

Because the volume of the $j$-th spherical shell is $4\pi r_j^2\,\Delta r_j$ and the volume of the ball $B(\mathbf{x}_i, L)$ is $V_L = \frac{4}{3}\pi L^3$, the entry of the discrete 3D top-hat kernel matrix is (`covariant_closure_50aa2db0.py`, lines 188–190):
$$\boxed{\;K_{ij} \;=\; \frac{4\pi r_j^2\,\mathrm{frac}(r_i, r_j; L)\,\Delta r_j}{\frac{4}{3}\pi L^3} \;=\; \frac{3\,r_j^2\,\Delta r_j}{L^3} \begin{cases}
  1, & r_i + r_j \le L, \\[4pt]
  \dfrac{L^2 - (r_i - r_j)^2}{4\,r_i\,r_j}, & |r_i - r_j| < L < r_i + r_j, \\[6pt]
  0, & |r_i - r_j| \ge L.
\end{cases}\;}$$
Consequently, the first variation of the nonlocal average $\bar{\rho}_L(r_i) = \sum_j K_{ij}\,\varepsilon(r_j)$ and of the factor $f_L(r_i) = \exp(-\gamma_L\bar{\rho}_L(r_i))$ with respect to the enthalpy perturbation $\delta\varepsilon(r_j) = \frac{d\varepsilon}{dh}\big|_{r_j}\delta h(r_j) = \frac{\varepsilon(r_j)+P(r_j)}{c_s^2(r_j)}\delta h(r_j)$ gives the first part of V23 Equation (9) identically:
$$\delta f_L(r_i) \;=\; -\gamma_L\,f_L(r_i)\sum_{j=1}^{N_{\text{tot}}} K_{ij}\,\frac{\varepsilon(r_j) + P(r_j)}{c_s^2(r_j)}\,\delta h(r_j).$$

### 2.2. Why Must Cell Boundaries Be Aligned with the Transition Radii $r_{\text{trans}}$?
In piecewise-polytrope equations of state (SLy, APR4, H4), the pressure $P(h)$ and energy density $\varepsilon(h)$ are continuous across the polytrope boundaries $h_{\text{tr}} \in \text{EOS}.h_b$, but the adiabatic index $\Gamma_k$ and hence the squared sound speed $c_s^2(h) = \frac{dP}{d\varepsilon}$ undergo a **discontinuous jump** at the transition radii $h(r_{\text{trans}}) = h_{\text{tr}}$! If a quadrature cell straddled a transition radius $r_{\text{trans}}$, the jump in the integrand $\frac{\varepsilon+P}{c_s^2}\delta h$ would introduce a first-order $\mathcal{O}(\Delta r)$ quadrature error, spoiling the turning-point agreement at $dM/d\rho_c = 0$. Therefore, `solve_coupled_radial_omega2` (`covariant_closure_50aa2db0.py`, lines 604–616) uses Brent's method (`brentq(lambda x: p["h"](x) - h_tr, 0.0, R)`) to locate all transition radii $0 < r_{\text{trans}} < R_*$ to machine precision and embeds them as fixed cell boundaries in the grid.

---

## 3. Matrix Form of the V23 Eq. (9) Coupled BVP, the $W \to \kappa W$ Scaling Invariance, and Galerkin Symmetrization

### 3.1. Operator Matrices ($\mathbf{S}$, $\mathbf{D}_M$, $\mathbf{D}_A$, $\mathbf{W}$) and the Fredholm Shooting Function $F(\omega^2)$
In geometrized units ($G = c = 1$, where $\gamma_L = L^2$), the two first-order equations of system V23 (9) for the vector of cell-center enthalpy perturbations $\delta\mathbf{h} = (\delta h(r_1),\dots,\delta h(r_{N_{\text{tot}}}))^T$ (`covariant_closure_50aa2db0.py`, lines 629–660) take the form:
1. **Linearized Mass-Source Matrix ($\mathbf{S}$) and Integrated Mass Perturbation ($\delta\mathbf{m} = \mathbf{D}_M\,\delta\mathbf{h}$):**
   $$\frac{d(\delta m)}{dr} \;=\; (\mathbf{S}\,\delta\mathbf{h})_i, \qquad \mathbf{S} \;=\; \begin{cases}
     \mathrm{diag}\!\bigl(4\pi r^2 f_L\bigr)\,\bigl[\mathbf{I} - L^2\,\mathrm{diag}(\varepsilon)\,\mathbf{K}\bigr]\,\mathrm{diag}\!\left(\dfrac{\varepsilon+P}{c_s^2}\right), & \text{(\texttt{mass} mode)}, \\[8pt]
     \mathrm{diag}\!\left(4\pi r^2 \dfrac{\varepsilon+P}{c_s^2}\right), & \text{(\texttt{gradient} mode)}.
   \end{cases}$$
   With the trapezoidal lower-triangular integration matrix $\bigl(\mathbf{M}_{\text{int}}\bigr)_{ij} = \Delta r_j\,\Theta(i - j) + \frac{1}{2}\Delta r_i\,\delta_{ij}$, the mass-perturbation operator is $\mathbf{D}_M \equiv \mathbf{M}_{\text{int}}\mathbf{S}$, i.e., $\delta\mathbf{m} = \mathbf{D}_M\,\delta\mathbf{h}$ (automatically satisfying the central boundary condition $\delta m(0) = 0$).
2. **Linearized Hydrostatic Gradient Matrix ($\mathbf{D}_A$):**
   Since $A(r) = \frac{m + 4\pi r^3 P}{r(r-2m)}$, the partial derivatives are $c_{dm}(r) \equiv \frac{\partial A}{\partial m} = \frac{1 + 2 r A}{r(r-2m)}$ and $c_{dh}(r) \equiv (\varepsilon+P)\frac{\partial A}{\partial P} = \frac{4\pi r^3(\varepsilon+P)}{r(r-2m)}$, so the base perturbation is $\mathbf{D}_{A,\text{base}} = \mathrm{diag}(c_{dm})\,\mathbf{D}_M + \mathrm{diag}(c_{dh})$, and:
   $$\delta\boldsymbol{\nu}' \;=\; \mathbf{D}_A\,\delta\mathbf{h}, \qquad \mathbf{D}_A \;=\; \begin{cases}
     \mathbf{D}_{A,\text{base}}, & \text{(\texttt{mass} mode)}, \\[6pt]
     \mathrm{diag}(f_L)\,\mathbf{D}_{A,\text{base}} \;-\; L^2\,\mathrm{diag}(f_L A)\,\mathbf{K}\,\mathrm{diag}\!\left(\dfrac{\varepsilon+P}{c_s^2}\right), & \text{(\texttt{gradient} mode)}.
   \end{cases}$$
3. **Adiabatic Mass-Conservation Relation and the Inertia Matrix ($\mathbf{W}$):**
   For an adiabatic radial displacement, Lagrangian mass conservation $\Delta m = 0 \implies \delta m(r) = -4\pi r^2 f_{\text{mass}}(\varepsilon + P)\,r\xi(r)$ ($f_{\text{mass}} = f_L$ in `mass` mode and $f_{\text{mass}} = 1$ in `gradient` mode), allowing the displacement to be expressed in terms of the mass perturbation:
   $$-\omega^2 e^{2(\lambda-\nu)} r\xi(r) \;=\; +\omega^2\,w_{\text{inert}}(r)\,\delta m(r), \qquad w_{\text{inert}}(r) \;\equiv\; \frac{e^{2(\lambda-\nu)}}{4\pi r^2 f_{\text{mass}}(r)\bigl(\varepsilon(r) + P(r)\bigr)} \;>\; 0.$$
   Therefore, the inertia matrix is $\mathbf{W}(\kappa) \equiv \kappa\,\mathrm{diag}(w_{\text{inert}})\,\mathbf{D}_M$ (where $\kappa \equiv W_{\text{scale}} > 0$ is the scaling factor), and integrating the V23 Eq. (9) enthalpy equation with central normalization $\delta h(0) = 1$ yields the discrete Fredholm system (`covariant_closure_50aa2db0.py`, lines 664–665):
   $$\boxed{\;\Bigl[\mathbf{I} \;+\; \mathbf{M}_{\text{int}}\bigl(\mathbf{D}_A \;+\; \omega^2\,\mathbf{W}(\kappa)\bigr)\Bigr]\,\delta\mathbf{h}(\omega^2; \kappa) \;=\; \mathbf{1}, \qquad F(\omega^2; \kappa) \;\equiv\; \delta m(R_*; \omega^2, \kappa) \;=\; \boldsymbol{\Delta r}^T \mathbf{S}\,\delta\mathbf{h}(\omega^2; \kappa) \;=\; 0.\;}$$
   The surface Birkhoff boundary condition $F(\omega^2; \kappa) = \delta m(R_*; \omega^2, \kappa) = 0$ is equivalent to the free-surface boundary condition $\delta m(R_*) = -4\pi R_*^3 f_{\text{mass}}(\varepsilon + P)\xi(R_*)$ and $\Delta P(R_*) = 0$.

### 3.2. Proof of the $W \to \kappa W$ Inertia-Scaling Invariance and the First-Order Tangent Theorem
Three exact algebraic consequences follow from the Fredholm equation above:
1. **Exact Scaling Identity ($\forall\,\kappa > 0$):**
   Because $\mathbf{W}(\kappa) = \kappa\,\mathbf{W}(1)$, the combination $\mathbf{D}_A + \omega^2\mathbf{W}(\kappa) = \mathbf{D}_A + (\kappa\omega^2)\mathbf{W}(1)$ depends on the eigenvalue $\omega^2$ and inertia scale $\kappa$ solely through the product $\tilde{\omega}^2 \equiv \kappa\,\omega^2$:
   $$\delta\mathbf{h}(\omega^2; \kappa) \;\equiv\; \delta\mathbf{h}(\kappa\omega^2; 1) \quad\Longrightarrow\quad F(\omega^2; \kappa) \;\equiv\; F(\kappa\omega^2; 1).$$
   Consequently, if $\omega^2_W$ is the eigenvalue at $\kappa = 1$ ($F(\omega^2_W; 1) = 0$), then under $\kappa$-fold inertia scaling the exact eigenvalue is:
   $$\boxed{\;\omega^2_{\kappa W} \;=\; \frac{\omega^2_W}{\kappa}, \qquad \frac{\omega^2_{4W}}{\omega^2_W} \;=\; \frac{1}{4} \;=\; 0.2500\quad\text{identically!}\;}$$
   Specifically, at $\log_{10}\rho_{0c} = 15.30427$ (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 135), $F(0; \kappa) \equiv F(0; 1) = 0.023252$ (independent of $\kappa$), $\omega^2_W = 4.897\times 10^6\;\text{s}^{-2}$, and $\omega^2_{4W}/\omega^2_W = 0.2500$.
2. **Inertia Invariance of the Static Turning Point ($\omega_0^2 = 0$):**
   At $\omega^2 = 0$, the term $\omega^2\mathbf{W}(\kappa)$ drops out, leaving $(\mathbf{I} + \mathbf{M}_{\text{int}}\mathbf{D}_A)\delta\mathbf{h}(0) = \mathbf{1}$. Because $\mathbf{D}_A$ and $\mathbf{S}$ are the exact Fréchet derivatives of the nonlocal static TOV equations with respect to the central enthalpy perturbation $\delta h(0)$, the solution $\delta\mathbf{h}(0)$ is proportional to the static sequence derivative $\frac{\partial h(r; \rho_{0c})}{\partial\log_{10}\rho_{0c}}$, and the shooting function is:
   $$F(0) \;=\; \delta m(R_*; 0) \;=\; \left(\frac{\partial h_c}{\partial\log_{10}\rho_{0c}}\right)^{-1}\frac{dM}{d\log_{10}\rho_{0c}}.$$
   Therefore $F(0) = 0 \iff \frac{dM}{d\rho_c} = 0$, explaining why the turning-point offset drops to $\Delta\log_{10}\rho_{0c} = +6.97\times 10^{-5}\;\text{dex}$ (`mass`, a $19.2\times$/$58.3\times$ improvement) and $+7.39\times 10^{-5}\;\text{dex}$ (`gradient`, a $48.3\times$/$59.3\times$ improvement), matching the pure GR discretization baseline ($+7.12\times 10^{-5}\;\text{dex}$, V23 Table 2).
3. **First-Order Tangent Estimate (`omega0_sq_sym_tangent_s2`) and the Galerkin Gram Matrix ($\mathbf{M}_{\text{sym}}$):**
   First-order Taylor expansion of the root $F(\omega^2) = 0$ around $\omega^2 = 0$ gives $F(0) + F'(0)\,\omega^2 + \mathcal{O}(\omega^4) = 0$, i.e. (`covariant_closure_50aa2db0.py`, lines 680–690):
   $$\omega^2_{\text{tangent}} \;\simeq\; -\frac{F(0)}{F'(0)}.$$
   Although the raw BVP matrix $\mathbf{W} = \mathrm{diag}(w_{\text{inert}})\mathbf{D}_M$ is not symmetric (since $\mathbf{D}_M = \mathbf{M}_{\text{int}}\mathbf{S}$ is multiplied by a lower-triangular integrator), the symmetrized Galerkin matrix of the quadratic kinetic form expressed in terms of the mass perturbations $\delta\mathbf{m} = \mathbf{D}_M\delta\mathbf{h}$ (`covariant_closure_50aa2db0.py`, lines 685–686) is:
   $$\boxed{\;\mathbf{M}_{\text{sym}} \;\equiv\; \mathbf{D}_M^T\,\mathrm{diag}\!\bigl(w_{\text{inert}}\,\Delta r\bigr)\,\mathbf{D}_M \;=\; \mathbf{M}_{\text{sym}}^T.\;}$$
   Because $w_{\text{inert}}(r_i)\,\Delta r_i > 0$ in every cell $i = 1,\dots,N_{\text{tot}}$ (in both `mass` and `gradient` modes, since the denominator contains the matter sum $\varepsilon + P > 0$!), for any vector $\mathbf{v} \neq \mathbf{0}$ we have $\mathbf{v}^T\mathbf{M}_{\text{sym}}\mathbf{v} = \sum_{i=1}^{N_{\text{tot}}} w_{\text{inert}}(r_i)\Delta r_i\,(\mathbf{D}_M\mathbf{v})_i^2 \ge 0$, with minimum diagonal entry $\min\mathrm{diag}\,\mathbf{M}_{\text{sym}} = 5.68\times 10^{-9} > 0$. As emphasized in the V23 manuscript (lines 135 and 152), this is by construction a **Gram-matrix property** ($\mathbf{B}^T\mathbf{B}$ form with $\mathbf{B} = \mathrm{diag}(\sqrt{w_{\text{inert}}\Delta r})\mathbf{D}_M$), whereas the positive-definiteness of the genuine dynamical kinetic Hessian was derived in Step 12 (`step12_second_variation_coupled_pulsation_hessian.md`) from the second variation $\delta^2 S$ of the action.

---

## 4. Checks

1. **Normalization of the 3D Spherical-Shell Intersection Kernel $K_{ij}$ ($\sum_j K_{ij} = 1$):**
   - When the stellar radius satisfies $R_* > r_i + L$ and the density is homogeneous ($\varepsilon(r_j) \equiv 1$), in the continuum limit:
     $$\sum_j K_{ij} \;\longrightarrow\; \frac{3}{L^3}\left[\int_0^{L - r_i} r'^2\,dr' \;+\; \int_{|L - r_i|}^{L + r_i} r'^2\,\frac{L^2 - (r_i - r')^2}{4\,r_i\,r'}\,dr'\right] \;=\; 1 \qquad \forall\,r_i \ge 0 \;\checkmark.$$
2. **Agreement Between $\omega^2_{\text{tangent}} \simeq -F(0)/F'(0)$ and $\omega^2_{\text{BVP}}$ Near the Turning Point (`coupled_stability_results.csv`):**
   - In row 23 of `coupled_stability_results.csv` (SLy GR, $\log_{10}\rho_{0c} = 15.302445$), $\omega^2_{\text{BVP}} = 1.04011\times 10^5\;\text{s}^{-2}$ and $\omega^2_{\text{tangent}} = 1.03992\times 10^5\;\text{s}^{-2}$ ($0.018\%$ difference); in row 56 (SLy `mass`, $L = 10\;\text{km}$, $\log_{10}\rho_{0c} = 15.308269$), $\omega^2_{\text{BVP}} = 9.84022\times 10^4\;\text{s}^{-2}$ and $\omega^2_{\text{tangent}} = 9.83845\times 10^4\;\text{s}^{-2}$ ($0.018\%$ difference), confirming the accuracy of first-order tangent linearization in the neighborhood of the $\omega_0^2 = 0$ turning point $\checkmark$.

---

## 5. Explicit Separation Between the V21 Nonlocal 3D Top-Hat Tangent Model (`step18`) and the V23 Local Curved `yukawa_action_el` Coupled Stability (`step12`, `step22`)

To distinguish the model branches unambiguously, we record the exact scope of the two stability calculations:

| Property | **Branch I-A**: V21 Nonlocal 3D Top-Hat (`step18`, `covariant_closure_50aa2db0.py`) | **Branch II-Fluid**: V23 Curved Exterior `yukawa_action_el` (`step12`, `step22`, `curved_yukawa_adm_and_stability.py`) |
|---|---|---|
| **Equilibrium Background** | V21 `tophat_mass` and `tophat_gradient` (`tp.dense_profile`) | Coupled interior + curved exterior ($r \in [0, R_*+12L]$) `yukawa_action_el` |
| **Perturbation Variables** | $(\delta h, \delta m)$ with nonlocal shell kernel $\delta f_L = -\gamma_L f_L K_{ij}\delta\varepsilon_j$ | Coupled $(\delta h, \delta\phi_L, \delta m)$ in the interior and curved exterior vacuum tail |
| **Mass Turning Point (SLy, $L=10\;\text{km}$)** | $\log_{10}\rho_{0c}(M_{\max}) = 15.29924$ (`mass`), $15.30106$ (`gradient`) | $\log_{10}\rho_{0c}(M_{\text{ADM},\max}) = 15.316404$ ($\alpha_Y=1$), $15.307166$ ($\alpha_Y=1/3$) |
| **$\omega_0^2 = 0$ Zero Crossing** | $15.29931$ (`mass`), $15.30113$ (`gradient`) | $15.316410$ ($\alpha_Y=1$), $15.307170$ ($\alpha_Y=1/3$) |
| **Turning-Point Offset $\Delta\log_{10}\rho_{0c}$** | $+6.97\times 10^{-5}\;\text{dex}$ (`mass`), $+7.39\times 10^{-5}\;\text{dex}$ (`gradient`) | $+5.80\times 10^{-6}\;\text{dex}$ ($\alpha_Y=1$, $2097\times$ reduction), $+4.20\times 10^{-6}\;\text{dex}$ ($\alpha_Y=1/3$, $995\times$ reduction; $L=1\;\text{km}$: $-6.86\times 10^{-6}\;\text{dex}$, $165\times$) |

> **V25 Update (curved `yukawa_action_el` branch):** 4-component exterior BVP ($\delta\nu_R=1.63\times 10^{-14}$), finite-difference verified exterior ($\partial\phi_L'(R_*)/\partial m(R_*)=-1.68365\times 10^{-4}\;\text{km}^{-2}$) and interior ($r=6\;\text{km}$ $\partial\phi_L''/\partial m=-1.72833\times 10^{-4}\;\text{km}^{-3}$, `scalar_m_sign = -1.0` main operator) derivatives, dynamical scalar wave, GL2 root convergence ($p_{\text{eff,root}}=4.248$). SLy $L=10$ km, $\alpha_Y=1$: $M_{\text{ADM}}$ turning point $15.316404$, main dynamical $\omega_0^2=0$ root $15.316410$ ($+5.80\times 10^{-6}$ dex). Detailed tables: `step22`, Sections 4 and 7.
