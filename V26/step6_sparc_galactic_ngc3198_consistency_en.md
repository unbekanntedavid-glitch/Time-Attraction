# Step 6: Galactic Rotation Curves (155 / 175 SPARC Galaxies), Positive-Definite 3D Top-Hat / Yukawa Kernel Positivity, and the NGC 3198 ($R = 30.08\text{ kpc}$) Negative Control Check

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step6_sparc_galactic_ngc3198_consistency_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step5_coupling_dictionary_ppn_screening_en.md` (`step5_coupling_dictionary_ppn_screening.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (Chapter 3, lines 256–265, Eq. 6) and `hu_p2.tex` (Chapters 4–5, lines 1–61, Eqs. 7–8; DOI: `10.5281/zenodo.23604423`)  
- V17 and V19 Zenodo records: DOI `10.5281/zenodo.19666220` (V17 SPARC 155 / 123 blind sample) and DOI `10.5281/zenodo.23497128` (V19 3D spherical top-hat, 175 SPARC galaxies)  
- Local SPARC and NGC 3198 working documents: `docs/README_draft.md` (lines 6–24), `negative_density.md` (lines 1–12 and 55–63), `foundation_mode.md` (lines 13–15 and 50–63)

---

## 1. Verbatim Citations from the V22 Monograph and Local SPARC Documents

In Steps 1–5, we constructed the V23 covariant scalar–tensor / Yukawa framework and verified its weak-field / Solar System (PPN) limit. In this step, we examine **the mathematical positivity of the 3D kernel on the galactic scale (exclusion of negative-density points across all 175 SPARC galaxies)** and the **NGC 3198 ($R = 30.08\text{ kpc}$) negative control point**, strictly separating the boost required by observation ($C_{\text{req}}$) from the boost actually predicted by the isolated model with fixed coupling $\alpha_{\text{std}} = 1/(4\pi)$ ($C_{\text{pred}}^{\text{isol}}$).

1. **Definition of the 3D Spherical Top-Hat and Exclusion of Negative Density ($\rho < 0$) (`hu_p2.tex`, lines 12–21, Eqs. 7–8):**
   > *"Instead of pointwise derivatives or linear approximations, I define the source of temporal damping at any point $\mathbf{x} \in \mathbb{R}^3$ as the unweighted three-dimensional top-hat average of the physical mass–energy density $\rho(\mathbf{y}) = T_{00}(\mathbf{y})/c^2 \ge 0$ inside a ball $B(\mathbf{x},L)$ of radius $L > 0$:"*
   $$\begin{equation}
     \bar{\rho}_L(\mathbf{x}) \;\equiv\; \frac{3}{4\pi L^3}\int_{|\mathbf{y}-\mathbf{x}|\le L} \rho(\mathbf{y})\,d^3y \;=\; \frac{M_{\text{enc}}(\mathbf{x};L)}{\frac{4}{3}\pi L^3}, \qquad
     f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(\mathbf{x})\bigr). \tag{V22 Eq. 7}
   \end{equation}$$
   > *"Consequently, the exponential damping factor satisfies the exact analytical bounds: $0 < f_L(\mathbf{x}) \le 1 \quad \forall\,\mathbf{x}\in\mathbb{R}^3,\;\forall\,L > 0$. Neither negative density, nor sign reversal, nor a metric horizon inversion arising from $f_L$ can occur."*

2. **The SPARC 175 / 155 / 123 Galaxy Sample and the 3D Axisymmetric Disk Model (`hu_p2.tex`, lines 47–55):**
   > *"The main SPARC catalog contains **175 late-type and dwarf galaxies** with $3.6\,\mu\text{m}$ surface photometry and resolved $\text{H\,\textsc{i}}/\text{H}\alpha$ rotation curves. During my V17 processing, after filtering out 20 systems exhibiting severe kinematic asymmetry or inclination quality flags ($i < 30^\circ$ or $Q=3$), the **processed rotation-curve sample of 155 galaxies** was obtained. Strict data-hygiene cuts (at least 5 radial measurement points and uninhibited $V_{\text{gas}}, V_{\text{disk}}, V_{\text{bulge}}$ baryonic decomposition) defined the **123-galaxy blind test sample**, encompassing **2596 radial velocity measurements**."*
   > *"In Version 19, I replaced the thin-disk midplane approximation $\rho_{\text{mid}}(r) = \Sigma_{\text{bary}}(r)/(2h_z)$ with the exact 3D axisymmetric mass distribution $\rho_{\text{bary}}(R,z) = \frac{\Sigma_{\text{bary}}(R)}{2h_z}\exp(-|z|/h_z)$ (with scale height $h_z = 0.196\,R_d^{0.633}$) [...]. Across all 175 SPARC galaxies and over the range $L \in [10^{16}, 10^{23}]\,\text{m}$, the maximum isolated galactic damping is $\max |1 - f_L| = 2.284\times 10^{-7}$ on the logarithmic grid (and **$2.300\times 10^{-7}$** under continuous optimization at $L \approx 23.09\,\text{kpc}$ in the massive spiral galaxy **UGC02487**)."*

3. **Exact Clarification of the NGC 3198 ($R = 30.08\text{ kpc}$) Velocity Values (`docs/README_draft.md`, lines 6–19 and `foundation_mode.md`, lines 13–14):**
   > *"At $R = 30.08\text{ kpc}$ (`NGC3198_rotmod.dat`, $D = 13.8\text{ Mpc}$):*
   > - **`Vobs` = `146 ± 2 km/s`**: *Observed rotation speed. '145' in earlier text is consistent with this within $1\sigma$. Neighbouring points: 148 at 28.16 kpc, 147 at 32.14 kpc.*
   > - **`Vbar` ($\Upsilon_{\text{disk}} = 0.5, \Upsilon_{\text{bul}} = 0.7$, signed gas) = `70.32 km/s`**: *Baryonic speed with the documented v16 mass-to-light ratio.*
   > - **`Vbar` ($\Upsilon_{\text{disk}} = 1.0$) = `88.53 km/s`**: *Baryonic speed computed from the raw SPARC columns $\sqrt{V_{\text{gas}}^2 + V_{\text{disk}}^2}$. **'88.5' is a baryonic value with $\Upsilon_{\text{disk}} = 1$, not an observed speed.***
   > *The required boost is $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$, which is $4.311$ for $\Upsilon = 0.5$ and $2.720$ for $\Upsilon = 1.0$."*

---

## 2. Derivation of the Galactic Circular Velocity and Positivity Check of the 3D Kernel

### 2.1. Circular Velocity from the Modified Euler Equation ($h' = -\nu' + \phi_L'$)
In Steps 4 and 5 (`step4_covariant_conservation_hydrostatic.md`, `step5_coupling_dictionary_ppn_screening.md`), we derived from the covariant conservation law ($\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$) that the radial effective acceleration acting on a baryonic element orbiting in the midplane ($z = 0$) of a galactic disk is:
$$a_{r,\text{eff}}(r) \;=\; -\frac{d\Phi_N(r)}{dr} \;+\; c^2\,\frac{d\phi_L(r)}{dr} \;=\; -\frac{V_{\text{bar}}^2(r)}{r} \;-\; c^2\,\frac{d\ln f_L(r)}{dr}.$$
The equilibrium circular velocity ($V_{\text{rot}}^2(r) = -r\,a_{r,\text{eff}}(r)$) can therefore be written in two ways:
1. **Additive scalar-gradient (V19–V23 canonical scalar–tensor / Yukawa) form:**
   $$\boxed{\;V_{\text{rot}}^2(r) \;=\; V_{\text{bar}}^2(r) \;-\; c^2\,r\,\frac{d\phi_L(r)}{dr} \;=\; V_{\text{bar}}^2(r) \;+\; c^2\,r\,\gamma_L\,\left|\frac{d\bar{\rho}_L(r)}{dr}\right|\;}$$
   (since $\frac{d\bar{\rho}_L}{dr} < 0 \implies \frac{d\phi_L}{dr} < 0$, the scalar gradient provides an inward additional attraction).
2. **Multiplicative $C(\rho)$ force-multiplier (V14–V16 phenomenological) form (`hu_p1.tex`, Eq. 6):**
   $$V_{\text{rot}}^2(r) \;=\; C(\rho(r))\,V_{\text{bar}}^2(r), \qquad C(\rho) \;=\; 1 + C_0\exp\!\left(-\frac{\rho}{\rho_c}\right).$$
   The exact algebraic relation between the two formulations at any radius $r$ is:
   $$C_{\text{eff}}(r) \;\equiv\; \frac{V_{\text{rot}}^2(r)}{V_{\text{bar}}^2(r)} \;=\; 1 \;-\; \frac{c^2\,r}{V_{\text{bar}}^2(r)}\,\frac{d\phi_L(r)}{dr}.$$

### 2.2. What Does the Result "0 Negative-Density Points Across 175 SPARC Galaxies" Prove, and What Does It Not Prove?
For peer-review clarity, we strictly state the exact meaning of the V19 (`sparc_tophat.py`) test performed across the 175 SPARC galaxies:
- **Failure of the Earlier Spherical Derivative Formula (`negative_density.md`, line 11; `foundation_mode.md`, lines 16–18 and 55–63):**
  When density is estimated from the numerical derivative of the rotation curve ($\rho_{\text{deriv}}(r) = \frac{1}{4\pi G r^2}\frac{d(r V_{\text{bar}}^2)}{dr}$), non-monotonic gas and disk profiles cause $\rho_{\text{deriv}}(r) < 0$ at **239** ($\Upsilon_{\text{disk}}=0.5$) and **303** ($\Upsilon_{\text{disk}}=1.0$) interior points across the 140 non-holdout SPARC galaxies (of which 22 points fall in the $29.5\text{--}31.5\;\text{kpc}$ band, closest to $30.37\;\text{kpc}$ in **NGC 4013** at $R = 30.31\;\text{kpc}$ with $\rho = -1.56\times 10^{-24}\;\text{kg/m}^3$).
- **Mathematical Positivity of the V19 3D Top-Hat and V23 Yukawa Kernels:**
  Both the 3D spherical top-hat kernel ($K_L(\mathbf{x},\mathbf{y}) = \frac{\Theta(L - |\mathbf{x}-\mathbf{y}|)}{\frac{4}{3}\pi L^3} \ge 0$) and the 3D Yukawa Green's function derived in Step 3 ($G_L(\mathbf{x}-\mathbf{y}) = \frac{e^{-|\mathbf{x}-\mathbf{y}|/L}}{4\pi L^2|\mathbf{x}-\mathbf{y}|} > 0$) **integrate the photometric surface-brightness density $\rho_{\text{bary}}(R,z) \ge 0$ with a strictly positive kernel**. Thus, at every radial point of every galaxy:
  $$\bar{\rho}_L(r) > 0, \qquad \phi_L^{\text{Yuk}}(r) > 0, \qquad 0 < f_L(r) \le 1 \quad \checkmark$$
- **Critical Methodological Distinction:** This result is a **proof of the mathematical positivity of the 3D integral kernel** on the positive photometric mass model of all 175 SPARC galaxies (`n_negative_density_points = 0`), **not** a rotation-curve fitting success of the constant-$\alpha_Y = 1$ ($\alpha_{\text{std}} = 1/(4\pi)$) Yukawa action! As demonstrated numerically in Section 3 below, the isolated model with fixed coupling $\alpha_{\text{std}} = 1/(4\pi)$ by itself does not reproduce flat galactic rotation curves.

---

## 3. Exact Numerical Negative Control: Isolated Fixed-Coupling ($\alpha_{\text{std}} = 1/(4\pi)$) Model vs. NGC 3198 ($R = 30.08\text{ kpc}$)

In accordance with peer-review requirements, the table below lists in separate columns for the control point `NGC3198_rotmod.dat` ($R = 30.08\;\text{kpc}$) the observed velocity ($V_{\text{obs}}$), the baryonic velocity ($V_{\text{bar}}$), the boost **required** by observation ($C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$), and the boost ($C_{\text{pred}}^{\text{isol}}$) and rotation velocity ($V_{\text{pred}}^{\text{isol}} = \sqrt{C_{\text{pred}}^{\text{isol}}}\,V_{\text{bar}}$) **actually predicted** by the isolated model with fixed $\gamma_L = G L^2/c^2$ ($\alpha_Y = 1$, $\alpha_{\text{std}} = 1/(4\pi) \approx 0.07958$):

| $\Upsilon_{\text{disk}},\,\Upsilon_{\text{bul}}$ Choice | Observed $V_{\text{obs}}$ ($\text{km/s}$) | Baryonic $V_{\text{bar}}$ ($\text{km/s}$) | Required Boost $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$ | Predicted Boost $C_{\text{pred}}^{\text{isol}}$ (Isolated 3D Top-Hat, $\gamma_L = GL^2/c^2$) | Predicted Boost $C_{\text{pred}}^{\text{isol,max}}$ (Point-Source Yukawa Upper Bound, $1+\alpha_{\text{std}}$) | Predicted Velocity $V_{\text{pred}}^{\text{isol}}$ (Isolated $\alpha_{\text{std}} = 1/(4\pi)$, $\text{km/s}$) | Observed Minus Predicted Residual $\Delta V = V_{\text{obs}} - V_{\text{pred}}^{\text{isol}}$ |
|---|---|---|---|---|---|---|---|
| **$\Upsilon_{\text{disk}} = 1.0,\;\Upsilon_{\text{bul}} = 0$** (raw SPARC columns) | $146.0 \pm 2.0$ | $88.53$ | **$2.720$** | $1 + \mathcal{O}(10^{-7}) \approx 1.0000002$ | $\le 1 + \frac{1}{4\pi} = \mathbf{1.07958}$ | **$88.53\text{--}91.99$** | **$-54.0\text{ to }-57.5\;\text{km/s}$** (does not fit the flat curve!) |
| **$\Upsilon_{\text{disk}} = 0.5,\;\Upsilon_{\text{bul}} = 0.7$** (V16/V17 mass-to-light ratio) | $146.0 \pm 2.0$ | $70.32$ | **$4.311$** | $1 + \mathcal{O}(10^{-7}) \approx 1.0000002$ | $\le 1 + \frac{1}{4\pi} = \mathbf{1.07958}$ | **$70.32\text{--}73.06$** | **$-72.9\text{ to }-75.7\;\text{km/s}$** (does not fit the flat curve!) |

### Physical Conclusion (Negative Control):
1. In an isolated galaxy, the depth of the gravitational potential is $\frac{G M_{\text{gal}}}{c^2 R} \sim \frac{V_{\text{bar}}^2}{c^2} \sim 10^{-7}$, and the 3D top-hat damping computed with the canonical coupling $\gamma_L = \frac{G L^2}{c^2}$ across all 175 SPARC galaxies satisfies $\max |1 - f_L| \le 2.30\times 10^{-7}$ (`sparc_tophat_summary.json`).
2. Even in the optimal transition zone $r \sim L$, the Yukawa fifth force with fixed $\alpha_Y = 1$ ($\alpha_{\text{std}} = \frac{1}{4\pi} \approx 0.07958$) increases $V_{\text{bar}}^2$ by at most a factor of $1 + \alpha_{\text{std}} \approx 1.0796$ (yielding $91.99\;\text{km/s}$ for $\Upsilon_{\text{disk}}=1.0$ and $73.06\;\text{km/s}$ for $\Upsilon_{\text{disk}}=0.5$, instead of the observed $146\;\text{km/s}$).
3. Consequently, the NGC 3198 ($R = 30.08\;\text{kpc}$) control point proves as an **exact negative control** that the isolated V23 Yukawa model with constant coupling $\alpha_{\text{std}} = 1/(4\pi)$ by itself **does not** explain flat galactic rotation curves; doing so requires either the V16/V17 phenomenological midplane boost $C(\rho) = 1 + C_0 e^{-\rho/\rho_c}$ ($C_0 \approx 3.34\text{--}5.09$) or a chameleon/symmetron extension with density-dependent coupling ($\beta_{\text{eff}}(\rho)$, $\alpha_{\text{eff}}^{\text{gal}} = 2\beta_{\text{eff}}^2 \approx 1.72\text{--}3.31$) going beyond the constant-$\alpha_Y$ V23 action (see Step 11, **Branch II-Cham**).

---

## 4. Checks

1. **NGC 3198 $R = 30.08\text{ kpc}$ Numerical Identities ($C_{\text{req}}$ vs. $C_{\text{pred}}^{\text{isol}}$):**
   - $V_{\text{bar}}(\Upsilon_{\text{disk}}=1.0) = 88.53\;\text{km/s} \implies C_{\text{req}} = (146 / 88.53)^2 = 2.7199 \approx 2.720$, whereas $C_{\text{pred}}^{\text{isol}} \le 1 + \frac{1}{4\pi} = 1.07958 \implies V_{\text{pred}}^{\text{isol}} \le 91.99\;\text{km/s} < 146\;\text{km/s}$ $\checkmark$.
   - $V_{\text{bar}}(\Upsilon_{\text{disk}}=0.5) = 70.32\;\text{km/s} \implies C_{\text{req}} = (146 / 70.32)^2 = 4.3107 \approx 4.311$, whereas $C_{\text{pred}}^{\text{isol}} \le 1.07958 \implies V_{\text{pred}}^{\text{isol}} \le 73.06\;\text{km/s} < 146\;\text{km/s}$ $\checkmark$.
2. **Positivity Check Across All (175) SPARC Galaxies:**
   - Photometric 3D disk model ($\rho_{\text{bary}}(R,z) \ge 0$) + 3D spherical top-hat or 3D Yukawa Green's function: number of points with $\rho \le 0$ is **0** across all 175 galaxies (`negative_density.md`, lines 5–7 and `hu_p2.tex`, lines 17–21) $\checkmark$.
