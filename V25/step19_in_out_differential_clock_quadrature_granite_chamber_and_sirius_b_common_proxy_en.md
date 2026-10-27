# Step 19: Exact Derivation of the V20 $8558\;\text{kg}$ Granite-Chamber 3D IN–OUT Differential Quadrature (Three Geometric $L$-Regimes, Top-Hat vs. Yukawa) and the Sirius B–Sirius A Common Effective-Mass Proxy (`common_background` vs. `fixedM`)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance_en.md` (`step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 3 and Table 3 + footnote, lines 161–189; DOI: `10.5281/zenodo.23639084`)  
- V23 verification script and dataset: `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/covariant_closure_50aa2db0.py` (Part 4, lines 700–885) and `clock_coupling_H0_H1.csv` (rows 1–13)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (Chapter 9, lines 426–455 and Chapter 10, point 1, lines 463–468; DOI: `10.5281/zenodo.23604423`)

---

## 1. Verbatim Citations from the V22 and V23 Manuscripts

In Step 8 (`step8_microscopic_clock_coupling_siriusb_lab_v20.md`), we derived the Dirac–Maxwell atomic-clock action ($S_{\text{atomic}}[\tilde{g}_{\mu\nu}, \psi_e, A_\mu]$) under H0 ($\tilde{g}_{\mu\nu} = g_{\mu\nu}$) and conditional H1 ($\tilde{g}_{00} = f_L^2 g_{00}$). Table 3 of Section 3 of the V23 manuscript (lines 170–189) and Part 4 of `covariant_closure_50aa2db0.py` (lines 700–885), however, contain two additional exact mathematical structures that we now derive step by step:

1. **Footnote to V23 Table 3 on the Granite-Chamber IN–OUT Differential Signal and the Sirius B Common Effective-Mass Proxy (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 188):**
   > *"V20 rows: H0 column = TAL excess over Newtonian GR at the IN ($\mathbf{r}=\mathbf{0}$) position; H1 columns = IN shift (IN$-$OUT differential signal, $D_{\text{out}}=2.0\,\text{m}$). Newtonian GR: $\Phi_N(\text{IN})/c^2 = -9.269\times 10^{-24}$, IN$-$OUT $= -6.094\times 10^{-24}$; $\sigma_y=10^{-18}$ is an illustrative reference value. On the Sirius B common effective-mass proxy (fixed $n=3$ Lane–Emden profile with the H0 effective-mass shift, not a newly solved equilibrium) at $L=10^3\,\text{km}$: H0 $80.69839$, H1 top-hat $80.69878\,\text{km/s}$."*

2. **Geometry and Mass Budget of the $8558\;\text{kg}$ Granite Chamber (`tal_book_HU_2026-10-24.tex`, Chapter 10, point 1, line 466):**
   > *"**Exact Geometry and Mass Budget:** The open cylindrical barrel of dimensions $R_{\text{in}} = 0.25\,\text{m}$, $R_{\text{out}} = 0.85\,\text{m}$, $H = 1.20\,\text{m}$ has volume $V_{\text{cyl}} = \pi(0.85^2 - 0.25^2)\times 1.20 = 2.48814\,\text{m}^3$ and mass **$M_{\text{cyl}} = 6842.39\,\text{kg}$**; the nominal $V_{\text{tot}} = 3.1102\,\text{m}^3$ (**$M_{\text{tot}} = 8558\,\text{kg}$**) includes an additional $0.6220\,\text{m}^3$ ($1710.6\,\text{kg}$) of endcap/support granite elements (or an outer radius $R_{\text{out}} = 0.942\,\text{m}$)."*

---

## 2. 3D IN–OUT Quadrature of the V20 $8558\;\text{kg}$ Granite Chamber and Derivation of the Three Geometric $L$-Regimes

### 2.1. Domain of the Cylindrical Barrel + Endcaps ($\Omega_{\text{granite}}$) and the Off-Axis Distance Function
In `covariant_closure_50aa2db0.py` (lines 802–822), the granite chamber of density $\rho_{\text{granite}} = 2750\;\text{kg/m}^3$ and mass $M_{\text{tot}} = 8558\;\text{kg}$ ($V_{\text{tot}} = M_{\text{tot}}/\rho_{\text{granite}} = 8558/2750 = 3.112000\;\text{m}^3$) is constructed in cylindrical coordinates ($\varrho \in [0, R_{\text{out}}]$, $z \in [-z_{\text{out}}, +z_{\text{out}}]$, $\varphi \in [0, 2\pi)$), exploiting $z \leftrightarrow -z$ mirror symmetry, as the union of:
- **Open Cylindrical Barrel:** $R_{\text{in}} = 0.25\;\text{m} \le \varrho \le R_{\text{out}} = 0.85\;\text{m}$, $|z| \le z_{\text{in}} = 0.60\;\text{m}$ ($H_{\text{in}} = 2 z_{\text{in}} = 1.20\;\text{m}$), whose volume and mass are:
  $$V_{\text{cyl}} \;=\; \pi\bigl(R_{\text{out}}^2 - R_{\text{in}}^2\bigr)\,(2 z_{\text{in}}) \;=\; 2.488141\;\text{m}^3, \qquad M_{\text{cyl}} \;=\; \rho_{\text{granite}}\,V_{\text{cyl}} \;=\; 6842.39\;\text{kg}.$$
- **Two Solid Endcaps:** $0 \le \varrho \le R_{\text{out}} = 0.85\;\text{m}$, $z_{\text{in}} \le |z| \le z_{\text{out}}$, where the endcap thickness $\Delta z_{\text{cap}} = z_{\text{out}} - z_{\text{in}}$ (`covariant_closure_50aa2db0.py`, line 808) is:
  $$\Delta z_{\text{cap}} \;=\; \frac{V_{\text{tot}} - V_{\text{cyl}}}{2\pi R_{\text{out}}^2} \;=\; \frac{3.112000 - 2.488141}{2\pi\times 0.85^2} \;=\; 0.137426\;\text{m} \quad\Longrightarrow\quad z_{\text{out}} \;=\; 0.60 + 0.137426 \;=\; 0.737426\;\text{m}.$$
  Thus, the characteristic domain $\Omega_{\text{granite}}$ occupied by granite on the half-meridian plane ($\varrho \in [0, R_{\text{out}}]$, $z \in [0, z_{\text{out}}]$) is:
  $$\Omega_{\text{granite}} \;=\; \bigl\{(\varrho, z) \in [0, R_{\text{out}}]\times [0, z_{\text{out}}] \;:\; (\varrho \ge R_{\text{in}}) \;\lor\; (z \ge z_{\text{in}})\bigr\}, \qquad dV \;=\; 2\cdot\bigl(2\pi\varrho\,d\varrho\,dz\bigr),$$
  where the factor of $2$ accounts for the sum over the $\pm z$ half-spaces (giving $M_{\text{granite,num}} = 8558.23\;\text{kg}$ on an $N_\varrho \times N_z = 500\times 500$ grid).

The Euclidean distance of a source point $(\varrho, \varphi, \pm z) \in \Omega_{\text{granite}}$ from the two clock positions (`covariant_closure_50aa2db0.py`, lines 820 and 834) is:
1. **IN Position ($\mathbf{r}_{\text{in}} = (0,0,0)$, Geometric Center of the Chamber):** independent of azimuth:
   $$d_{\text{in}}(\varrho, z) \;=\; \sqrt{\varrho^2 + z^2}.$$
2. **OUT Position ($\mathbf{r}_{\text{out}} = (D_{\text{out}}, 0, 0)$, $D_{\text{out}} = 2.0\;\text{m}$ in the Midplane):** by the law of cosines:
   $$d_{\text{out}}(\varrho, z, \varphi) \;=\; \sqrt{\varrho^2 + D_{\text{out}}^2 - 2\varrho D_{\text{out}}\cos\varphi + z^2},$$
   integrated over $N_\varphi = 256$ uniform azimuthal points ($\varphi_k = \frac{2\pi(k + 1/2)}{256}$, $dm_p = dm_g / N_\varphi$). From this, the Newtonian gravitational potential and standard GR redshift difference (`clock_coupling_H0_H1.csv`, rows 8–13) are:
   $$\frac{\Phi_N(\text{IN})}{c^2} \;=\; -\frac{G}{c^2}\int_{\Omega_{\text{granite}}}\frac{dm_g}{d_{\text{in}}} \;=\; -9.268639\times 10^{-24}, \qquad \frac{\Phi_N(\text{OUT})}{c^2} \;=\; -\frac{G}{c^2}\int_{\Omega_{\text{granite}}}\frac{dm_p}{d_{\text{out}}} \;=\; -3.174686\times 10^{-24},$$
   $$\Delta y_{\text{GR}}(\text{IN}-\text{OUT}) \;=\; \frac{\Phi_N(\text{IN}) - \Phi_N(\text{OUT})}{c^2} \;=\; -6.093953\times 10^{-24}.$$

### 2.2. Exact Proof of the Four Geometric Threshold Distances and Top-Hat vs. Yukawa IN–OUT Behavior
The extremal Euclidean distances of the chamber domain $\Omega_{\text{granite}}$ from the IN and OUT points are:
- **At the IN Position ($\mathbf{r}_{\text{in}} = \mathbf{0}$):**
  $$d_{\min}^{\text{IN}} \;=\; \min_{(\varrho,z)\in\Omega_{\text{granite}}}\sqrt{\varrho^2 + z^2} \;=\; \min(R_{\text{in}}, z_{\text{in}}) \;=\; \min(0.25,\,0.60) \;=\; 0.25\;\text{m},$$
  $$d_{\max}^{\text{IN}} \;=\; \max_{(\varrho,z)\in\Omega_{\text{granite}}}\sqrt{\varrho^2 + z^2} \;=\; \sqrt{R_{\text{out}}^2 + z_{\text{out}}^2} \;=\; \sqrt{0.85^2 + 0.737426^2} \;=\; 1.12530\;\text{m}.$$
- **At the OUT Position ($D_{\text{out}} = 2.00\;\text{m}$):**
  $$d_{\min}^{\text{OUT}} \;=\; \min_{\Omega_{\text{granite}}} d_{\text{out}} \;=\; D_{\text{out}} - R_{\text{out}} \;=\; 2.00 - 0.85 \;=\; 1.15\;\text{m} \quad (\varphi = 0,\; \varrho = R_{\text{out}},\; z = 0),$$
  $$d_{\max}^{\text{OUT}} \;=\; \max_{\Omega_{\text{granite}}} d_{\text{out}} \;=\; \sqrt{(D_{\text{out}} + R_{\text{out}})^2 + z_{\text{out}}^2} \;=\; \sqrt{2.85^2 + 0.737426^2} \;=\; 2.94385\;\text{m} \quad (\varphi = \pi,\; \varrho = R_{\text{out}},\; z = z_{\text{out}}).$$

From these four threshold distances ($0.25\;\text{m} < 1.1253\;\text{m} < 1.15\;\text{m} < 2.9439\;\text{m}$), the sharp physical contrast between the 3D spherical top-hat ($\phi_L^{\text{top-hat}}(\mathbf{r}) = \frac{3 G}{4\pi c^2 L}\int_{d \le L} dm$) and the local Yukawa kernel ($\phi_L^{\text{Yuk}}(\mathbf{r}) = \frac{G}{4\pi c^2}\int \frac{e^{-d/L}}{d}\,dm$ at $\alpha_Y = 1$) follows immediately (`clock_coupling_H0_H1.csv`, rows 8–13):

1. **Cavity Regime ($L < d_{\min}^{\text{IN}} = 0.25\;\text{m}$, e.g. $L = 0.20\;\text{m}$):**
   Both balls $B(\mathbf{0}, L)$ and $B(\mathbf{r}_{\text{out}}, L)$ remain outside the granite walls (inside the inner cavity and outside the chamber, respectively), so $\phi_L^{\text{top-hat}}(\text{IN}) = \phi_L^{\text{top-hat}}(\text{OUT}) = 0$, giving $\Delta y_{\text{H1}}^{\text{top-hat}} = 0$. By contrast, the Yukawa Green's function $e^{-d/L}/d$ senses the walls across the cavity: $\Delta y_{\text{H1}}^{\text{Yuk}}(\text{IN}) = -3.690\times 10^{-26}$ and $\Delta y_{\text{H1}}^{\text{Yuk}}(\text{IN}-\text{OUT}) = -3.683\times 10^{-26}$.
2. **Purely Interior Top-Hat Regime ($d_{\min}^{\text{IN}} < L < d_{\min}^{\text{OUT}} = 1.15\;\text{m}$, e.g. $L = 0.50,\;0.85,\;1.00\;\text{m}$ in V23 Table 3):**
   Because $L < 1.15\;\text{m}$, the OUT ball $B(\mathbf{r}_{\text{out}}, L)$ **does not intersect the granite chamber** ($d_{\text{out}} \ge 1.15\;\text{m} > L$), so:
   $$\boxed{\;\phi_L^{\text{top-hat}}(\text{OUT}) \;\equiv\; 0 \quad\Longrightarrow\quad \Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}-\text{OUT}) \;\equiv\; \Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}) \qquad \forall\,L < 1.15\;\text{m}!\;}$$
   This is why V23 Table 3 (lines 183–185) has the exact same number before and inside the parentheses in the Top-Hat column: $-3.318\times 10^{-25}\;(-3.318\times 10^{-25})$ ($L=0.50\;\text{m}$), $-1.303\times 10^{-24}\;(-1.303\times 10^{-24})$ ($L=0.85\;\text{m}$), and $-1.436\times 10^{-24}\;(-1.436\times 10^{-24})$ ($L=1.00\;\text{m}$). By contrast, the exponential tail of the Yukawa kernel is non-zero even at $D_{\text{out}} = 2.0\;\text{m}$ ($\phi_L^{\text{Yuk}}(\text{OUT}) = +3.735\times 10^{-26}$ at $L = 1.00\;\text{m}$), so the Yukawa differential signal is smaller in magnitude than the absolute IN signal: $-3.781\times 10^{-25}\;(-3.407\times 10^{-25})$ at $L = 1.00\;\text{m}$!
3. **Partial Exterior Overlap ($d_{\min}^{\text{OUT}} < L < d_{\max}^{\text{OUT}}$, e.g. $L = 1.50\;\text{m}$):**
   Here the IN ball already encloses the entire chamber ($L = 1.50\;\text{m} > d_{\max}^{\text{IN}} = 1.1253\;\text{m}$), while the OUT ball partially intersects the chamber: $\Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}) = -1.0115\times 10^{-24}$, $\Delta y_{\text{H1}}^{\text{top-hat}}(\text{OUT}) = -1.0064\times 10^{-25}$, yielding $\Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}-\text{OUT}) = -9.1086\times 10^{-25}$.
4. **Common Enclosure Regime ($L \ge d_{\max}^{\text{OUT}} = 2.94385\;\text{m}$, e.g. $L = 10.0\;\text{m}$):**
   When $L \ge 2.944\;\text{m}$, both balls $B(\mathbf{0}, L)$ and $B(\mathbf{r}_{\text{out}}, L)$ enclose the entire $M_{\text{tot}} = 8558\;\text{kg}$ granite chamber! Consequently, the 3D spherical top-hat average at both positions is $\bar{\rho}_L^{\text{top-hat}}(\text{IN}) = \bar{\rho}_L^{\text{top-hat}}(\text{OUT}) = \frac{3 M_{\text{tot}}}{4\pi L^3}$, so the **Top-Hat IN$-$OUT differential signal vanishes identically**:
   $$\boxed{\;\Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}-\text{OUT})\Big|_{L \ge 2.944\;\text{m}} \;=\; 0\quad(\text{numerically } -5.28\times 10^{-40}\text{ rounding noise at } L=10\;\text{m}),\;}$$
   whereas the Yukawa kernel, owing to its $1/d$ weighting, still produces a distinct differential signal at $L = 10.0\;\text{m}$: $\Delta y_{\text{H1}}^{\text{Yuk}}(\text{IN}-\text{OUT}) = -4.818\times 10^{-25}$!
5. **Second-Order H0 TAL Excess in the Cavity (`dy_H0_cavity`):**
   Under H0 there is no direct $f_L$ clock prefactor, so the excess over Newtonian GR arises from the $1 - f_L \simeq \phi_L$ modification of the effective gravitating mass of the chamber (`covariant_closure_50aa2db0.py`, line 852):
   $$\Delta y_{\text{H0}}^{\text{excess}}(\text{IN}) \;=\; \left|\frac{\Phi_N(\text{IN})}{c^2}\right|\,\phi_L^{\text{top-hat}}(\text{IN}) \;=\; 9.268639\times 10^{-24}\times 1.436124\times 10^{-24} \;=\; +1.331091\times 10^{-47} \quad (L = 1.00\;\text{m}).$$

---

## 3. Derivation of the Sirius B–Sirius A Differential Redshift and the `common_background` vs. `fixedM` Proxy

In `covariant_closure_50aa2db0.py` (lines 710–799), the white dwarf Sirius B ($M_B = 1.018 \pm 0.011\,M_\odot$, $R_B = 0.00803 \pm 0.00011\,R_\odot = 5586.471\;\text{km}$) and the main-sequence reference star Sirius A ($M_A = 2.042\,M_\odot$, $R_A = 1.711\,R_\odot = 1.19034\times 10^6\;\text{km}$) are modeled as $n = 3$ Lane–Emden polytropes:
1. **GR Baseline and Combined Uncertainty (`covariant_closure_50aa2db0.py`, lines 719–723):**
   The surface compactness is $\mathcal{C}_B = \frac{GM_B}{R_B c^2} = 2.690878\times 10^{-4}$, giving the exact Schwarzschild redshift and velocity equivalent:
   $$z_{\text{GR},B} \;=\; \bigl(1 - 2\mathcal{C}_B\bigr)^{-1/2} - 1 \;=\; 2.6918805\times 10^{-4}, \qquad v_{\text{GR}} \;=\; c\,z_{\text{GR},B} \;=\; 80.700550\;\text{km/s}.$$
   The parameter uncertainty is $\sigma_{\text{GR}} = v_{\text{GR}}\sqrt{(\sigma_M/M_B)^2 + (\sigma_R/R_B)^2} = 1.4080\;\text{km/s}$, which combined in quadrature with the observational uncertainty $\sigma_{\text{obs}} = 0.77\;\text{km/s}$ yields $\sigma_{\text{comb}} = \sqrt{\sigma_{\text{obs}}^2 + \sigma_{\text{GR}}^2} = 1.604806\;\text{km/s}$ (`clock_coupling_H0_H1.csv`, row 2).
2. **Sirius B–Sirius A Differential H1 Clock Shift (`covariant_closure_50aa2db0.py`, lines 767–773):**
   Because the astrometric gravitational redshift is measured from the spectral-line difference between Sirius B and Sirius A, and under H1 the surface redshift is $1 + z_{\text{H1}} = (1 + z_{\text{metric}})\,f_L(R_*)^{-1} = (1 + z_{\text{metric}})\,e^{+x_*}$ (where $x_* \equiv \gamma_L\bar{\rho}_L(R_*)$), the direct H1 excess velocity between the two stars is:
   $$\boxed{\;\Delta v_{\text{H1,extra}}(L) \;=\; \Delta v_B(L) \;-\; \Delta v_A(L) \;=\; c\,(1 + z_{\text{GR},B})\bigl(e^{x_B(L)} - 1\bigr) \;-\; c\,\bigl(e^{x_A(L)} - 1\bigr).\;}$$
   - At $L = 1000\;\text{km} = 10^6\;\text{m}$: $\Delta v_{\text{H1,extra}}^{\text{top-hat}} = +0.000392\;\text{km/s}$ and $\Delta v_{\text{H1,extra}}^{\text{Yuk}} = +0.044215\;\text{km/s}$.
   - At $L = 5000\;\text{km} = 5\times 10^6\;\text{m}$: $\Delta v_{\text{H1,extra}}^{\text{top-hat}} = +4.688557\;\text{km/s}$ ($v_{\text{H1}}^{\text{top-hat}} = 85.38911\;\text{km/s}$, $2.95\sigma$) and $\Delta v_{\text{H1,extra}}^{\text{Yuk}} = +2.150821\;\text{km/s}$ ($v_{\text{H1}}^{\text{Yuk}} = 82.85137\;\text{km/s}$, $1.37\sigma$).
3. **Difference Between Fixed-Mass (`fixedM`) and the Common Effective-Mass Proxy (`common_background`) (`covariant_closure_50aa2db0.py`, lines 758–766 and 794–798; V23 line 188):**
   - If the gravitating mass $M_B = 1.018\,M_\odot$ of the star is held fixed from orbital astrometry (`fixedM`), the exterior Schwarzschild metric is unchanged, so $v_{\text{H0,fixedM}} \equiv v_{\text{GR}} = 80.70055\;\text{km/s}$, and $v_{\text{H1,fixedM}}(L) = v_{\text{GR}} + \Delta v_{\text{H1,extra}}(L)$ ($80.70094\;\text{km/s}$ at $L = 10^3\;\text{km}$).
   - If instead one evaluates the `mass`-mode effective mass reduction at a fixed $n = 3$ Lane–Emden baryonic profile $\rho_B(r)$ (`common_background` proxy):
     $$\frac{\Delta M_{\text{eff}}(L)}{M_B} \;=\; \frac{\int_0^{R_B} 4\pi r^2\rho_B(r)\,\bigl[f_L(r) - 1\bigr]\,dr}{\int_0^{R_B} 4\pi r^2\rho_B(r)\,dr} \;<\; 0,$$
     then the H0 prediction is $v_{\text{H0,mass}}(L) = v_{\text{GR}}\bigl(1 + \frac{\Delta M_{\text{eff}}(L)}{M_B}\bigr)$ ($80.69839\;\text{km/s}$ at $L = 10^3\;\text{km}$ and $80.69487\;\text{km/s}$ at $L = 5\times 10^3\;\text{km}$), while the H1 value built on that same common effective-mass proxy is:
     $$\boxed{\;v_{\text{H1,common}}(L) \;=\; v_{\text{H0,mass}}(L) \;+\; \Delta v_{\text{H1,extra}}(L) \;=\; 80.698385 + 0.000392 \;=\; 80.69878\;\text{km/s}\quad(\text{at } L = 10^3\;\text{km, Top-Hat}).\;}$$

---

## 4. Checks

1. **Ratio of Top-Hat and Yukawa $\phi_L(\text{IN})$ at Large $L$ ($L = 10.0\;\text{m} \gg d_{\max}^{\text{IN}} = 1.125\;\text{m}$):**
   - When $L \gg d_{\max}^{\text{IN}}$, the entire mass $M_{\text{tot}}$ lies inside the ball of volume $\frac{4}{3}\pi L^3$ in the 3D spherical top-hat average, so $\phi_L^{\text{top-hat}}(\text{IN}) = \frac{G L^2}{c^2}\frac{3 M_{\text{tot}}}{4\pi L^3} = \frac{3 G M_{\text{tot}}}{4\pi c^2 L} = 1.517259\times 10^{-25}$ (at $L = 10\;\text{m}$), whereas in the Yukawa integral $e^{-d_{\text{in}}/L} \approx 1 - d_{\text{in}}/L$, yielding:
     $$\phi_L^{\text{Yuk}}(\text{IN})\Big|_{L \gg d_{\max}^{\text{IN}}} \;=\; \frac{G}{4\pi c^2}\int_{\Omega_{\text{granite}}}\frac{e^{-d_{\text{in}}/L}}{d_{\text{in}}}\,dm_g \;\approx\; \frac{1}{4\pi}\left|\frac{\Phi_N(\text{IN})}{c^2}\right| - \frac{G M_{\text{tot}}}{4\pi c^2 L} \;=\; 7.3757\times 10^{-25} - 0.5058\times 10^{-25} \;\approx\; 6.87\times 10^{-25},$$
     in close agreement with the exact quadrature value ($6.888137\times 10^{-25}$) in row 13 of `clock_coupling_H0_H1.csv` $\checkmark$.
2. **Verification of the Sirius B `common_background` vs. `fixedM` Difference at $L = 10^3\;\text{km}$:**
   - $v_{\text{H1,fixedM}}^{\text{top-hat}} - v_{\text{H1,common}}^{\text{top-hat}} = 80.70094147 - 80.69877704 = +0.00216443\;\text{km/s} \equiv v_{\text{GR}} - v_{\text{H0,mass}}$ ($80.70054975 - 80.69838532 = +0.00216443\;\text{km/s}$) to machine precision $\checkmark$.
