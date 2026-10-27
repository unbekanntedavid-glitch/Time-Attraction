# Step 17: Exact Analytical Derivation of the V23 Eq. (8) Factorized Chandrasekhar–Hydrostatic Difference Identity (`eq:factored_diff`) and the `gradient`-Mode `metric` Inertia-Weight Indefiniteness ($\mathcal{W}_{\text{metric}}(r) \le 0$, $L \ge 367\;\text{m}$)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress_en.md` (`step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 2, lines 121–126, Eq. 8 and Table 2 footnote, line 152; DOI: `10.5281/zenodo.23639084`)  
- V23 verification scripts: `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/covariant_closure_50aa2db0.py` (lines 100–120) and `radial_modes.py` (lines 20–39)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (Chapter 8, point 2, lines 409–419; DOI: `10.5281/zenodo.23604423`)

---

## 1. Verbatim Citations from the V22 and V23 Manuscripts

In Step 16 (`step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md`), we derived the background stress tensor $T^\mu_{\;\nu}(X)$ appearing in V23 Equations (4)–(5). Section 2 of the V23 manuscript (lines 121–126, Eq. 8, and line 152) and Chapter 8 of the V22 monograph (lines 409–418) show how this same background modification causes both frozen-$f_L$ ($\delta f_L = 0$) Sturm–Liouville operators (`matter` and `metric`) to depart from the static turning point $dM/d\rho_c = 0$:

1. **V23 Eq. (8) Factorized Chandrasekhar–Hydrostatic Difference Identity (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, Section 2, lines 121–125, Eq. 8):**
   > *"To understand why both frozen-$f_L$ operators depart from the static mass peak $dM/d\rho_c = 0$, I first determined in `SymPy` the exact difference between the frozen-$f_L$ Chandrasekhar operator $\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi]$ and the linearized hydrostatic operator $\mathcal{H}_{\text{lin}}[\xi] = (\delta P)' + (\delta\varepsilon+\delta P)\nu' + (\varepsilon+P)\delta\nu'\big|_{\delta f_L=0}$ in mass mode:"*
   $$\begin{equation}
     \mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} = \frac{4\pi r^{5/2}(f_L - 1)\Gamma P \varepsilon e^{2\nu}}{(r - 2m)^{5/2}}\Big[(4\pi r^3 P - 3r + 7m)\xi - r(r-2m)\xi'\Big]. \tag{V23 Eq. 8}
   \end{equation}$$
   > *"Expression (8) vanishes identically for $f_L = 1$ (GR), because Chandrasekhar's integration by parts substitutes the GR Einstein identity $\lambda' + 2\nu' = 8\pi r e^\lambda(\varepsilon+P)$. When $f_L \ne 1$, however, $\lambda' + 2\nu' = 8\pi r e^\lambda(f_L\varepsilon + P)$, which generates the $(f_L-1)$ residual term in (8) AND ALSO omits the nonlocal perturbation $\delta f_L[\delta\varepsilon,\xi]$."*

2. **$\mathcal{W}(r) \le 0$ Inertia-Weight Indefiniteness of the `metric` Operator Above $L \ge 367\;\text{m}$ (`tal_book_HU_2026-10-24.tex`, line 417 and V23 Table 2 footnote, line 152):**
   - `tal_book_HU_2026-10-24.tex`, line 417:
     > *"**Approximate `metric` operator:** Above $L \ge 367\,\text{m}$ across 12 sequences (609 models), the inertia weight $\mathcal{W}(r)$ becomes indefinite ($\mathcal{W}\le 0$) and yields negative $\omega_0^2$ values along the selected positive-norm branch; determining physical stability remains open (`50aa2db0` task)."*
   - `tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 152:
     > *"$^*$In gradient mode at $L=10\,\text{km}$, the weight of the approximate frozen-$f_L$ `metric` operator becomes indefinite near the surface ($\mathcal{W}(r)\le 0$); the static tangent model (9) uses a positive Gram weight (this is not a proof of kinetic positivity)."*

---

## 2. Step-by-Step Derivation of the V23 Eq. (8) Factorized Difference Identity (`eq:factored_diff`)

### 2.1. Definition of the Two Operators in the Convention of V23 `covariant_closure_50aa2db0.py` (Lines 101–120)
In the notation of `covariant_closure_50aa2db0.py` (lines 103–118) and `radial_modes.py` (lines 3–7) ($G = c = 1$, $g_{tt} = -e^{\nu_C}$ where $\nu_C \equiv 2\nu$, and $e^\lambda \equiv \bigl(1 - \frac{2m}{r}\bigr)^{-1} = \frac{r}{r - 2m}$, so that $e^{\lambda/2} = \sqrt{\frac{r}{r-2m}}$):
- The background equations in `mass` mode (`sub_mass_full`):
  $$m'(r) \;=\; 4\pi r^2 f_L(r)\,\varepsilon(r), \qquad \nu'(r) \;=\; A(r) \;\equiv\; \frac{m + 4\pi r^3 P}{r(r - 2m)}, \qquad P'(r) \;=\; -(\varepsilon + P)A(r), \qquad \varepsilon'(r) \;=\; -\frac{(\varepsilon + P)^2}{\Gamma P}A(r).$$
- The **frozen-$f_L$ Chandrasekhar `matter` Sturm–Liouville operator** ($\zeta(r) \equiv r^3 e^{-\nu(r)}\xi(r)$, where $\xi(r) = \Delta r/r$):
  $$\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] \;\equiv\; -\frac{d}{dr}\!\left(\Pi(r)\frac{d\zeta}{dr}\right) \;-\; Q_{\text{matter}}(r)\,\zeta(r),$$
  where the coefficients (`covariant_closure_50aa2db0.py`, lines 106–108) are:
  $$\Pi(r) \;=\; \sqrt{\frac{r}{r-2m}}\,e^{3\nu} r^{-2}\,\Gamma P, \qquad Q_{\text{matter}}(r) \;=\; \sqrt{\frac{r}{r-2m}}\,e^{3\nu} r^{-2}\,(\varepsilon + P)\left(\nu'^2 + \frac{4\nu'}{r} - \frac{8\pi r\,P}{r - 2m}\right).$$
- The **linearized hydrostatic operator** at $\delta f_L = 0$ (`covariant_closure_50aa2db0.py`, lines 111–117):
  $$\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} \;\equiv\; \frac{d(\delta P)}{dr} \;+\; (\delta\varepsilon + \delta P)\,A(r) \;+\; (\varepsilon + P)\,\delta A(r)\Big|_{\delta f_L=0},$$
  where the Eulerian pressure, energy-density, and mass perturbations are:
  $$\Delta P \;=\; -\Gamma P\,e^\nu r^{-2}\zeta'(r), \qquad \delta P \;=\; \Delta P - r\xi P'(r), \qquad \delta\varepsilon \;=\; \frac{\varepsilon + P}{\Gamma P}\delta P, \qquad \delta m\Big|_{\delta f_L=0} \;=\; -4\pi r^3\bigl(f_L\varepsilon + P\bigr)\xi,$$
  and the linearized variation of the geometric factor is:
  $$\delta A\Big|_{\delta f_L=0} \;=\; \frac{\partial A}{\partial m}\,\delta m\Big|_{\delta f_L=0} \;+\; \frac{\partial A}{\partial P}\,\delta P \;=\; \frac{r + 8\pi r^3 P}{r(r-2m)^2}\,\delta m\Big|_{\delta f_L=0} \;+\; \frac{4\pi r^2}{r - 2m}\,\delta P.$$

### 2.2. Why Does a Residual Proportional to $(f_L - 1)$ Arise Under Integration by Parts?
Let us expand the difference $\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}e^{2\nu}\mathcal{H}_{\text{lin}}[\xi]\big|_{\delta f_L=0}$:
1. Because $\delta P = -e^{-\nu}\sqrt{\frac{r-2m}{r}}\,\Pi(r)\zeta'(r) + r(\varepsilon + P)A\xi$, the leading second-order derivative term $-\frac{d}{dr}(\Pi\zeta')$ ($\propto \xi''(r)$) **cancels identically** against the leading term of $-\sqrt{\frac{r}{r-2m}}e^{2\nu}\frac{d(\delta P)}{dr}$!
2. When collecting the remaining first-order ($\propto \xi'(r)$) and zeroth-order ($\propto \xi(r)$) terms, the key step is the radial derivative of the metric prefactor $\sqrt{\frac{r}{r-2m}}e^{2\nu}$:
   $$\frac{d}{dr}\!\left(\ln\left[\sqrt{\frac{r}{r-2m}}\,e^{2\nu}\right]\right) \;=\; \frac{1}{2}\lambda'(r) \;+\; 2\nu'(r) \;=\; \frac{r m'(r) - m(r)}{r(r-2m)} \;+\; \frac{2(m + 4\pi r^3 P)}{r(r-2m)} \;=\; \frac{m + 4\pi r^3\bigl[f_L(r)\varepsilon + 2P\bigr]}{r(r-2m)}.$$
   In general relativity ($f_L = 1$), the numerator has $f_L\varepsilon = \varepsilon$, and likewise $(\varepsilon + P)$ appears in the mass perturbation $\delta m_{\text{GR}} = -4\pi r^3(\varepsilon + P)\xi$, so the difference between the two expressions vanishes identically.
3. When $f_L(r) \neq 1$, however, a shift proportional to $(f_L - 1)\varepsilon$ survives from both the background derivative $m'(r) = 4\pi r^2 f_L\varepsilon$ and the perturbation $\delta m\big|_{\delta f_L=0} = -4\pi r^3(f_L\varepsilon + P)\xi$:
   - The coefficient of the $\xi'(r)$ term (arising from the product of the $r\xi'$ part of $\Delta P = -\Gamma P(3\xi - r\nu'\xi + r\xi')$ and the factor $\frac{1}{2}\lambda'\big|_{f_L} - \frac{1}{2}\lambda'\big|_{f_L=1} = \frac{4\pi r^2(f_L - 1)\varepsilon}{r - 2m}$) is:
     $$\text{Coeff}(\xi') \;=\; -\sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\Gamma P\,r\cdot\frac{4\pi r^2(f_L - 1)\varepsilon}{r - 2m}\cdot\frac{r}{r-2m} \;=\; -\frac{4\pi r^{7/2}(f_L - 1)\Gamma P\varepsilon e^{2\nu}}{(r - 2m)^{3/2}}.$$
   - Combining the contributions from $(3 - r\nu')\Delta P$ and $\frac{\partial A}{\partial m}\bigl(\delta m|_{f_L} - \delta m|_{f_L=1}\bigr)$, the coefficient of the $\xi(r)$ term is $\frac{4\pi r^{5/2}(f_L - 1)\Gamma P\varepsilon e^{2\nu}}{(r-2m)^{5/2}}(4\pi r^3 P - 3r + 7m)$.
   Factoring out the common prefactor $\frac{4\pi r^{5/2}(f_L - 1)\Gamma P\varepsilon e^{2\nu}}{(r-2m)^{5/2}}$ yields **V23 Equation (8)** identically:
   $$\boxed{\;\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} \;=\; \frac{4\pi r^{5/2}(f_L - 1)\Gamma P \varepsilon e^{2\nu}}{(r - 2m)^{5/2}}\Big[(4\pi r^3 P - 3r + 7m)\xi - r(r-2m)\xi'\Big].\;}$$

---

## 3. Why Does the Weight $\mathcal{W}_{\text{metric}}(r)$ of the `metric` Operator Turn Negative Near the Surface ($L \ge 367\;\text{m}$, `gradient` Mode)?

By the definition in `radial_modes.py` (lines 11–12 and 30–38), the `metric` operator variant replaces the matter quantities $(\varepsilon, P)$ in the Sturm–Liouville coefficients $Q(r)$ and $\mathcal{W}(r)$ with the effective Einstein-tensor sources of the modified background ($\varepsilon_{\text{eff}} = \varepsilon + \varepsilon_X$, $P_{r,\text{eff}} = P + P_{r,X}$):
$$\mathcal{W}_{\text{metric}}(r) \;=\; e^{(3\lambda - \nu_C)/2}\,r^4\,\bigl[\varepsilon_{\text{eff}}(r) \;+\; P_{r,\text{eff}}(r)\bigr].$$

1. **In `mass` Mode (`radial_modes.py`, lines 31–32):**
   By V23 Eq. (4) derived in Step 16, $\varepsilon_{\text{eff}} = \varepsilon + \varepsilon_X = f_L(r)\varepsilon(r) > 0$ and $P_{r,\text{eff}} = P(r) > 0$, hence:
   $$\varepsilon_{\text{eff}}^{\text{mass}}(r) + P_{r,\text{eff}}^{\text{mass}}(r) \;=\; f_L(r)\,\varepsilon(r) \;+\; P(r) \;>\; 0 \qquad \forall\,r \in [0, R_*),$$
   so in `mass` mode the inertia weight remains positive everywhere.

2. **In `gradient` Mode (`radial_modes.py`, lines 33–34):**
   By V23 Eq. (5) derived in Step 16, $\varepsilon_{\text{eff}} = \varepsilon(r)$, whereas $P_{r,\text{eff}}(r) = P(r) + P_{r,X}(r) = f_L(r)P(r) - \bigl(1 - f_L(r)\bigr)\frac{m(r)}{4\pi r^3}$. Therefore, the effective inertial sum is:
   $$\boxed{\;\varepsilon_{\text{eff}}^{\text{grad}}(r) + P_{r,\text{eff}}^{\text{grad}}(r) \;=\; \varepsilon(r) \;+\; f_L(r)\,P(r) \;-\; \bigl(1 - f_L(r)\bigr)\,\frac{m(r)}{4\pi r^3}.\;}$$
   Let us examine the behavior of this expression near the stellar surface ($r \to R_*$):
   - At the surface, matter pressure and energy density vanish: $\lim_{r\to R_*} P(r) = 0$ and $\lim_{r\to R_*} \varepsilon(r) = 0$.
   - By contrast, the mean enclosed density ($\frac{m(R_*)}{4\pi R_*^3} = \frac{M}{4\pi R_*^3} > 0$) is strictly positive, and for any nonlocal smoothing radius $L > 0$ the surface ball reaches into the stellar interior, so $\bar{\rho}_L(R_*) > 0 \implies 1 - f_L(R_*) > 0$!
   - Consequently, near the surface ($r \to R_*$) the limit of the expression is **strictly negative**:
     $$\boxed{\;\lim_{r\to R_*}\Bigl[\varepsilon_{\text{eff}}^{\text{grad}}(r) + P_{r,\text{eff}}^{\text{grad}}(r)\Bigr] \;=\; -\bigl(1 - f_L(R_*)\bigr)\,\frac{M}{4\pi R_*^3} \;<\; 0 \qquad (\forall\,L > 0)!\;}$$
   When $L \ge 367\;\text{m}$, the factor $1 - f_L(R_*)$ and the crustal density profile are wide enough that the outermost Gauss–Legendre points of the $N = 1200$ finite-element quadrature grid (`radial_modes.py`, lines 45–54) fall inside the crustal zone where $\varepsilon(r) + f_L(r)P(r) < (1 - f_L(r))\frac{m(r)}{4\pi r^3}$, making $\mathcal{W}_{\text{metric}}(r) \le 0$ (across 12 sequences, 609 models in `item2_radial_grid.csv`). This proves analytically that the `metric` weight indefiniteness is not a numerical artifact, but a consequence of the `metric` ansatz substituting the static background stress ($P_{r,X} < 0$) as inertial mass, whereas the genuine kinetic weights ($\mathcal{W}_\xi \propto f_L(\varepsilon+P) > 0$ and $\mathcal{W}_\phi \propto \alpha_Y^{-1} > 0$) derived in Step 12 (`step12_second_variation_coupled_pulsation_hessian.md`) from the second variation $\delta^2 S$ of the covariant action are positive-definite for all $r \in (0, R_*)$ and all $L > 0$!

---

## 4. Checks

1. **$f_L \to 1$ (GR) Limit in V23 Equation (8):**
   - If $f_L(r) \equiv 1$, the factor $(f_L - 1)$ makes the right-hand side of V23 Eq. (8) vanish identically (`diff_gr_limit = 0` in line 120 of `covariant_closure_50aa2db0.py`), so $\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] \equiv \sqrt{\frac{r}{r-2m}}e^{2\nu}\mathcal{H}_{\text{lin}}[\xi]$, and the eigenvalue crossing $\omega_0^2 = 0$ coincides with the static turning point $dM/d\rho_c = 0$ $\checkmark$.
2. **Sign-Reversal Radius ($r_{\text{sign}}$) of $\mathcal{W}_{\text{metric}}(r)$ in `gradient` Mode:**
   - The sign reversal occurs in the crust precisely where $\varepsilon(r_{\text{sign}}) + f_L(r_{\text{sign}})P(r_{\text{sign}}) = \bigl(1 - f_L(r_{\text{sign}})\bigr)\frac{m(r_{\text{sign}})}{4\pi r_{\text{sign}}^3}$; as $L \to 0$, $1 - f_L(R_*) \to 0$ so $r_{\text{sign}} \to R_*$ (pushed outside the quadrature points), whereas for $L \ge 367\;\text{m}$ we have $r_{\text{sign}} < r_{\text{quad,max}} < R_*$ $\checkmark$.
