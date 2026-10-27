# Step 14: Conformal vs. Disformal Gravitational Lensing (Null-Geodesic Light Deflection Angle $\hat{\alpha}_{\text{def}}(b)$), NICER Ray-Tracing, and Cosmological Background Embedding ($L_{\text{hor}}^{(m)} \approx 22.94\;\text{Gpc}$ vs. $L_{\text{hor}}^{(b)} \approx 57.98\;\text{Gpc}$)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step14_gravitational_lensing_null_geodesics_cosmological_background_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step13_curved_window_metric_variation_synge_world_function_en.md` (`step13_curved_window_metric_variation_synge_world_function.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (Chapter 2, line 199), `hu_p2.tex` (Chapter 4, line 31; Chapter 5, lines 52–55), and `hu_p3.tex` (Chapter 9, lines 59–60; Chapter 11, line 104; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 3, lines 161–186; DOI: `10.5281/zenodo.23639084`)

---

## 1. Verbatim Citations from the V22–V23 Manuscripts

In Step 5 (`step5_coupling_dictionary_ppn_screening.md`), we derived the weak-field PPN parameter $\gamma_{\text{PPN}}$, and in Step 9 (`step9_vhe_photon_propagation_hess_dispersion.md`) we analyzed straight line-of-sight VHE gamma-ray time-of-flight dispersion. In this step, we derive three closely connected results: **(A)** what difference the purely temporal disformal clock metric ($\tilde{g}_{00} = f_L^2 g_{00}$, $\tilde{g}_{ij} = g_{ij}$, **Branch II-Disf**) induces in the **light deflection angle (gravitational lensing, $\hat{\alpha}_{\text{def}}(b)$)** compared to the conformal case $\tilde{g}_{\mu\nu}^{\text{conf}} = f_L^2 g_{\mu\nu}$ (**Branch II-Conf**); **(B)** how **cosmological background embedding** ($\gamma_L\rho_{\text{bg}} = 1$) enters the exponential factor, computing separately the horizon scales associated with the baryonic density ($\Omega_b = 0.0493 \implies L_{\text{hor}}^{(b)} \approx 57.98\;\text{Gpc}$) and the total cosmic matter density ($\Omega_m \approx 0.315 \implies L_{\text{hor}}^{(m)} \approx 22.94\;\text{Gpc}$); and **(C)** what this implies for **NICER pulsar ray-tracing likelihoods** (PSR J0030+0451, PSR J0740+6620, PSR J0437-4715).

1. **Temporal Refractive Index and Cosmological Background Embedding in the V22 Monograph (`hu_p1.tex` line 199 and `hu_p2.tex` line 54):**
   - `hu_p1.tex`, line 199:
     > *"I interpret the flow of time as possessing an effective local carrier density (analogous to an optical medium with refractive index $n_t(\mathbf{x}) = f_L(\mathbf{x})^{-1} = \exp(+\gamma_L \bar{\rho}_L(\mathbf{x}))$ or a compressible time fluid)."*
   - `hu_p2.tex`, Chapter 5, line 54:
     > *"**Cosmological Background Embedding ($L \approx 22.94\,\text{Gpc}$):** For $\bar{\rho}_L$ on megaparsecs scales to reach $\gamma_L\bar{\rho}_L \sim 1$, the ball $B(\mathbf{x},L)$ must include the ambient intergalactic and cosmological background density $\rho_{\text{bg}} \approx \Omega_b\rho_{\text{crit}}$. At uniform cosmic baryon density, the condition $\gamma_L\bar{\rho}_{\text{bg}} = 1$ occurs at the formal horizon scale $L \approx 22.94\,\text{Gpc}$ ($7.08\times 10^{26}\,\text{m}$)."*

2. **NICER Pulsar Likelihoods in the V22 Monograph (`hu_p3.tex`, lines 59–60 and 104):**
   - `hu_p3.tex`, Chapter 9, line 60:
     > *"The joint log-likelihood computed for the three NICER pulsars (PSR J0030+0451, PSR J0740+6620, PSR J0437-4715) under causal endpoint ($M_{\text{max}}(\text{causal})$) mass normalization shows that at $L = 1\,\text{km}$ the shift relative to GR is at most **$4.8\times 10^{-4}\,\text{dex}$** ($\Delta\ln\mathcal{L} \approx 1.1\times 10^{-3}$). Consequently, **current NICER mass–radius data cannot statistically distinguish the scale $L$ from GR under either hypothesis**."*

---

## 2. Conformal vs. Purely Temporal Disformal Gravitational Lensing ($\hat{\alpha}_{\text{def}}(b)$)

Consider a static, weak gravitational field in which the line element of the Einstein-frame gravitational metric ($g_{\mu\nu}$) in isotropic coordinates is:
$$ds^2 \;=\; g_{\mu\nu}dx^\mu dx^\nu \;=\; -\left(1 + \frac{2\Phi_N(\mathbf{x})}{c^2}\right)c^2 dt^2 \;+\; \left(1 - \frac{2\Psi_N(\mathbf{x})}{c^2}\right)d\mathbf{x}^2, \qquad \Phi_N(\mathbf{x}) = \Psi_N(\mathbf{x}) = -\frac{G M(<\!r)}{r}.$$

### 2.1. Light Deflection under Purely Conformal Jordan Coupling (Branch II-Conf: $\tilde{g}_{\mu\nu}^{\text{conf}} = f_L^2 g_{\mu\nu}$)
If the Maxwell field couples to the conformally rescaled metric $\tilde{g}_{\mu\nu}^{\text{conf}} = f_L(\mathbf{x})^2 g_{\mu\nu} = e^{-2\phi_L(\mathbf{x})}g_{\mu\nu}$, then in $D=4$ spacetime dimensions the conformal invariance $\sqrt{-\tilde{g}^{\text{conf}}}\tilde{g}_{\text{conf}}^{\mu\alpha}\tilde{g}_{\text{conf}}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} = \sqrt{-g}g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}$ ensures that null geodesics of light rays ($\tilde{g}_{\mu\nu}^{\text{conf}}dx^\mu dx^\nu = 0 \iff g_{\mu\nu}dx^\mu dx^\nu = 0$) are **independent of the scalar field $\phi_L(\mathbf{x})$**. The effective optical refractive index entering Fermat's principle ($\delta\int n_{\text{conf}}(\mathbf{x})\,dl = 0$) is:
$$n_{\text{conf}}(\mathbf{x}) \;=\; \sqrt{\frac{g_{rr}}{-g_{00}}} \;\simeq\; 1 \;-\; \frac{\Phi_N(\mathbf{x}) + \Psi_N(\mathbf{x})}{c^2} \;=\; 1 \;-\; \frac{2\Phi_N(\mathbf{x})}{c^2}.$$
Consequently, the deflection angle of a light ray with impact parameter $b$ ($\hat{\alpha}_{\text{def}}(b) = -\int_{-\infty}^{+\infty}\nabla_\perp\ln n\,dz$) under purely conformal coupling **senses only the baryonic Weyl potential ($\Phi_{\text{lens}}^{\text{conf}} = \frac{1}{2}(\Phi_N + \Psi_N) = \Phi_N$)**:
$$\boxed{\;\hat{\alpha}_{\text{def}}^{\text{conf}}(b) \;=\; \frac{2}{c^2}\int_{-\infty}^{+\infty}\nabla_\perp\Phi_N\bigl(\sqrt{b^2+z^2}\bigr)\,dz \;=\; \frac{4 G M_{\text{bary}}(<\!b)}{c^2 b}.\;}$$
That is, in purely conformal coupling the scalar field ($\phi_L$) enhances the circular orbital velocity of non-relativistic stars and gas ($V_{\text{rot}}^2(r) = V_{\text{bar}}^2(r) - c^2 r\phi_L'(r)$), **but does not contribute to gravitational lensing**!

### 2.2. Light Deflection in the V23 Purely Temporal Disformal Clock Metric (Branch II-Disf: $\tilde{g}_{00} = f_L^2 g_{00}$, $\tilde{g}_{ij} = g_{ij}$)
By contrast, if electromagnetic waves propagate on the V23 disformal physical metric ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$), then the null-cone equation ($d\tilde{s}_{\text{disf}}^2 = 0$) is:
$$0 \;=\; \tilde{g}_{00}^{\text{disf}}\,c^2 dt^2 \;+\; \tilde{g}_{ij}^{\text{disf}}\,dx^i dx^j \;=\; -f_L(\mathbf{x})^2\left(1 + \frac{2\Phi_N(\mathbf{x})}{c^2}\right)c^2 dt^2 \;+\; \left(1 - \frac{2\Psi_N(\mathbf{x})}{c^2}\right)d\mathbf{x}^2.$$
From this, the effective optical refractive index associated with Fermat's principle (`hu_p1.tex`, line 199: $n_t(\mathbf{x}) = f_L(\mathbf{x})^{-1} = e^{+\phi_L(\mathbf{x})}$) is:
$$n_{\text{disf}}(\mathbf{x}) \;=\; \frac{c}{v_{\text{phase}}(\mathbf{x})} \;=\; f_L(\mathbf{x})^{-1}\sqrt{\frac{1 - 2\Psi_N/c^2}{1 + 2\Phi_N/c^2}} \;\simeq\; 1 \;-\; \frac{\Phi_N(\mathbf{x}) + \Psi_N(\mathbf{x})}{c^2} \;+\; \phi_L(\mathbf{x}).$$
Substituting the transverse gradient $\nabla_\perp\ln n_{\text{disf}} \simeq -\frac{1}{c^2}\nabla_\perp(\Phi_N + \Psi_N - c^2\phi_L)$ into the Fermat integral, the disformal lensing potential and light deflection angle are:
$$\boxed{\;\Phi_{\text{lens}}^{\text{disf}}(\mathbf{x}) \;=\; \Phi_N(\mathbf{x}) \;-\; \frac{1}{2}c^2\,\phi_L(\mathbf{x}), \qquad \hat{\alpha}_{\text{def}}^{\text{disf}}(b) \;=\; \frac{4 G M_{\text{bary}}(<\!b)}{c^2 b} \;-\; \int_{-\infty}^{+\infty}\frac{b}{r}\,\frac{d\phi_L(r)}{dr}\,dz.\;}$$
Since $\frac{d\phi_L}{dr} < 0$ for an attractive time potential (density decreases outward), the second term produces a **positive (inward-directed) excess light deflection**: the disformal temporal refractive index ($n_t = e^{+\phi_L}$) around galaxies and galaxy clusters amplifies not only rotation curves, but gravitational lensing as well!

---

## 3. Cosmological Background Embedding: Exact Arithmetic for $\Omega_b = 0.0493$ ($57.98\;\text{Gpc}$) and $\Omega_m = 0.315$ ($22.94\;\text{Gpc}$)

Let us decompose the cosmic mass–energy density into a spatially homogeneous FLRW background density ($\rho_{\text{bg}}(t) = \Omega_{\text{bg}}\,\rho_{\text{crit},0}$) and the local overdensity ($\delta\rho(\mathbf{x})$):
$$\rho(\mathbf{x}) \;=\; \rho_{\text{bg}} \;+\; \delta\rho(\mathbf{x}) \quad\implies\quad \bar{\rho}_L(\mathbf{x}) \;=\; \frac{3}{4\pi L^3}\int_{|\mathbf{y}-\mathbf{x}|\le L}\bigl[\rho_{\text{bg}} + \delta\rho(\mathbf{y})\bigr]d^3y \;=\; \rho_{\text{bg}} \;+\; \delta\bar{\rho}_L(\mathbf{x}).$$

### 3.1. General Formula and Exact Separation of the Two Cosmological Density Parameters ($\Omega_b$ vs. $\Omega_m$)
Using the Planck 2018 (TT,TE,EE+lowE+lensing) parameters:
$$H_0 \;=\; 67.4\;\text{km\,s}^{-1}\text{Mpc}^{-1} \;=\; 2.184285\times 10^{-18}\;\text{s}^{-1}, \qquad R_H \;\equiv\; \frac{c}{H_0} \;=\; 1.37250\times 10^{26}\;\text{m} \;=\; 4.4480\;\text{Gpc},$$
$$\rho_{\text{crit},0} \;=\; \frac{3 H_0^2}{8\pi G} \;=\; 8.53287\times 10^{-27}\;\text{kg/m}^3.$$
The exact general solution of the unit-background-exponent condition ($\gamma_L\rho_{\text{bg}} = 1$, where $\gamma_L = \frac{G L_{\text{hor}}^2}{c^2}$) for any density fraction $\Omega_{\text{bg}} \equiv \rho_{\text{bg}}/\rho_{\text{crit},0}$ is:
$$\gamma_L\,\rho_{\text{bg}} \;=\; \frac{G L_{\text{hor}}^2}{c^2}\,\Omega_{\text{bg}}\,\frac{3 H_0^2}{8\pi G} \;=\; \frac{3\,\Omega_{\text{bg}}}{8\pi}\left(\frac{L_{\text{hor}}}{c/H_0}\right)^2 \;=\; 1 \quad\Longleftrightarrow\quad \boxed{\;L_{\text{hor}}(\Omega_{\text{bg}}) \;=\; \sqrt{\frac{8\pi}{3\,\Omega_{\text{bg}}}}\,\frac{c}{H_0} \;=\; \frac{c}{\sqrt{G\,\rho_{\text{bg}}}}.\;}$$
Observe that Newton's constant $G$ cancels identically from the product $\gamma_L\rho_{\text{bg}}$. To clarify the numerical values cited in line 54 of V22 `hu_p2.tex`, it is essential to separate the purely baryonic ($\Omega_b$) and combined cosmic matter ($\Omega_m$) cases:

| Background Density Model | Density Parameter $\Omega_{\text{bg}}$ | Background Density $\rho_{\text{bg}} = \Omega_{\text{bg}}\rho_{\text{crit},0}$ ($\text{kg/m}^3$) | Dimensionless Factor $\sqrt{\frac{8\pi}{3\Omega_{\text{bg}}}}$ | Horizon Scale $L_{\text{hor}}$ ($\text{m}$) | Horizon Scale $L_{\text{hor}}$ ($\text{Gpc}$) | Note / Source |
|---|---|---|---|---|---|---|
| **Purely Baryonic Cosmic Background ($\Omega_b$)** | $\Omega_b = 0.0493$ | $4.2067\times 10^{-28}$ | **$13.0306$** | **$1.7885\times 10^{27}$** | **$57.96\text{--}57.98\;\text{Gpc}$** | Exact $\gamma_L\rho_{b,0} = 1$ solution for $\Omega_b = 0.0493$ |
| **Cosmic Matter Density Background ($\Omega_m$)** | $\Omega_m = 0.3147\text{--}0.3150$ | $2.6855\text{--}2.6879\times 10^{-27}$ | **$5.1583\text{--}5.1605$** | **$7.080\text{--}7.083\times 10^{26}$** | **$22.94\text{--}22.95\;\text{Gpc}$** | Value $7.08\times 10^{26}\;\text{m} \approx 22.94\;\text{Gpc}$ in V22 `hu_p2.tex` line 54 (corresponds to $\Omega_m \approx 0.315$!) |

*(Note on line 54 of V22 `hu_p2.tex`: the manuscript prose verbally mentions $\Omega_b\rho_{\text{crit}}$, but the reported numerical value $7.08\times 10^{26}\;\text{m} = 22.94\;\text{Gpc}$ arises arithmetically from the cosmic matter density $\Omega_m \approx 0.315$, $\rho_{\text{bg}} = \frac{1}{\gamma_L(L_{\text{hor}})} = 2.686\times 10^{-27}\;\text{kg/m}^3$ and the factor $\sqrt{8\pi/(3\Omega_m)} = 5.160$, whereas for $\Omega_b = 0.0493$ one obtains $\sqrt{8\pi/(3\Omega_b)} = 13.031$ and $L_{\text{hor}}^{(b)} = 1.789\times 10^{27}\;\text{m} \approx 57.98\;\text{Gpc}$.)*

### 3.2. Factorization of the Homogeneous Background Factor from Local Dynamical and Lensing Gradients
By Cauchy's multiplicative property (`hu_p1.tex`, Eq. 2: $f(\rho_1+\rho_2) = f(\rho_1)f(\rho_2)$), the exponential time factor factorizes as a product:
$$f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L[\rho_{\text{bg}} + \delta\bar{\rho}_L(\mathbf{x})]\bigr) \;=\; \underbrace{e^{-\gamma_L\rho_{\text{bg}}}}_{f_{L,\text{bg}}\;(\text{global constant})}\,\cdot\,\underbrace{e^{-\gamma_L\,\delta\bar{\rho}_L(\mathbf{x})}}_{f_{L,\text{loc}}(\mathbf{x})}.$$
Since $f_{L,\text{bg}}$ is spatially constant ($\nabla f_{L,\text{bg}} = 0$), rescaling the global time coordinate $d\tilde{t} = f_{L,\text{bg}}\,dt$ absorbs the homogeneous background factor into the cosmological background clock, while the logarithmic gradient appearing in local dynamical acceleration ($\mathbf{a}_{\text{eff}} = -\nabla\Phi_N + c^2\nabla\phi_L$) and disformal light deflection ($\hat{\alpha}_{\text{def}}^{\text{disf}}$) is **independent of the homogeneous background $\rho_{\text{bg}}$**:
$$\nabla\ln f_L(\mathbf{x}) \;=\; -\nabla\phi_L(\mathbf{x}) \;=\; -\gamma_L\,\nabla\delta\bar{\rho}_L(\mathbf{x}).$$

---

## 4. Implication for NICER Pulsar Ray-Tracing (PSR J0030+0451, PSR J0740+6620, PSR J0437-4715)

In ray-tracing of X-ray photons emitted from the surface ($r \ge R_*$) of a neutron star of radius $R_*$ (Beloborodov / Pechenick–Ftaclas–Cohen formula) in the exterior vacuum ($\rho(r > R_*) = 0$):
- **For a 3D top-hat kernel ($L \le 1\;\text{km} \ll R_*$):** Outside the stellar surface ($r > R_* + L$), the top-hat sphere does not intersect the star, so $\bar{\rho}_L(r > R_* + L) = 0 \implies f_L = 1$; exterior light deflection follows exact Schwarzschild geodesics, and at $L = 1\;\text{km}$ the relative shift of the mass–radius curve is only $4.8\times 10^{-4}\;\text{dex}$ ($\Delta\ln\mathcal{L} \approx 1.1\times 10^{-3}$, V22 `hu_p3.tex` line 60).
- **For a Yukawa kernel (`yukawa_action_el`):** The exterior scalar profile in curved exterior spacetime decays exponentially as $r \to \infty$ ($\phi_L(r) \propto r^{-(1+M_{\text{ADM}}/L)}e^{-r/L}$), and the asymptotic gravitational mass of the star is absorbed into the ADM mass $M_{\text{ADM}} = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$ derived in Steps 7 and 22 (`step7_neutron_star_tov_stability.md`, `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`).

---

## 5. Checks

1. **Vacuum and GR Limit ($\phi_L \to 0 \implies f_L \to 1$):**
   - As $\phi_L(\mathbf{x}) \to 0$, $\frac{d\phi_L}{dr} \to 0$, so both the conformal ($\hat{\alpha}_{\text{def}}^{\text{conf}}$) and disformal ($\hat{\alpha}_{\text{def}}^{\text{disf}}$) light deflection angles recover the Einstein value $\hat{\alpha}_{\text{GR}}(b) = \frac{4GM}{c^2 b}$ $\checkmark$.
2. **Numerical Verification of $L_{\text{hor}}^{(m)} \approx 22.94\;\text{Gpc}$ and $L_{\text{hor}}^{(b)} \approx 57.98\;\text{Gpc}$:**
   - $\Omega_m = 0.3147 \implies \rho_{m,0} = 2.6855\times 10^{-27}\;\text{kg/m}^3 \implies \sqrt{8\pi/(3\Omega_m)} = 5.1605 \implies L_{\text{hor}}^{(m)} = 7.083\times 10^{26}\;\text{m} = 22.95\;\text{Gpc}$ (and for $\rho_{\text{bg}} = 2.686\times 10^{-27}\;\text{kg/m}^3$, $7.08\times 10^{26}\;\text{m} = 22.94\;\text{Gpc}$) $\checkmark$.
   - $\Omega_b = 0.0493 \implies \rho_{b,0} = 4.2067\times 10^{-28}\;\text{kg/m}^3 \implies \sqrt{8\pi/(3\Omega_b)} = 13.0306 \implies L_{\text{hor}}^{(b)} = 1.7885\times 10^{27}\;\text{m} = 57.96\text{--}57.98\;\text{Gpc}$ $\checkmark$.
