# Step 24 — Bullet Cluster (`1E 0657-558`) Gravitational Lensing Convergence $\kappa(\mathbf{\theta})$ and Gas vs. Galaxy Peak Offset Analysis at $L = 10\;\text{km}$, $\alpha_Y = 1$ (`step24_en`, V26)

## I. Physical and Methodological Framework: Synthetic 1D Illustrative Collision-Axis Test of the Bullet Cluster (`1E 0657-558`)

In `step23_gw170817_and_lensing_en.md`, we proved that within the canonical V24 (`10.5281/zenodo.23680715`) and V25 (`10.5281/zenodo.23699685`) short-range scalar–tensor Yukawa baseline ($L = 10\;\text{km} = 10^4\;\text{m}$, $\alpha_Y = 1$, $\beta_{\text{Yuk}} = 2\tilde{\alpha}_L\alpha_Y^2 = \alpha_Y^3/(4\pi) = 1/(4\pi) \approx 0.0795775$), the tensor gravitational-wave propagation speed on the purely conformal minimal branch is identically $c_T = c$ ($\Delta c_T/c = 0$), and under the conformal Jordan metric ($\tilde{g}_{\mu\nu} = A^2(\phi_L)g_{\mu\nu} = e^{-2\alpha_Y\phi_L}g_{\mu\nu}$) the linear Weyl lensing potential satisfies $\Phi_{\text{lens}} = \frac{1}{2}(\tilde{\Psi} + \tilde{\Phi}) = \Phi_N$. In Step 24 (`step24`), we apply this framework to a **synthetic 1D illustrative collision-axis model** of the **Bullet Cluster (`1E 0657-558`)** (Markevitch et al. 2004, *ApJ* 606, 819; Clowe et al. 2006, *ApJ* 648, L109; Bradač et al. 2006, *ApJ* 652, 937).

### 1.1. Cosmological and Lensing Parameters (`1E 0657-558`)
- **Lens redshift:** $z_L = 0.296$, angular diameter distance $D_L = 910.0\;\text{Mpc}$ ($1' \leftrightarrow 264.7\;\text{kpc}$).
- **Characteristic background source redshift:** $z_S \approx 0.8\text{–}1.0$ (fiducial value $z_S = 0.85$).
- **Critical surface mass density:**
  $$\Sigma_{\text{crit}} = \frac{c^2 D_S}{4\pi G D_L D_{LS}} = 3.10\times 10^9\;M_\odot\,\text{kpc}^{-2} = 0.6474\;\text{g}\,\text{cm}^{-2}.$$
- **Supersonic merger shock velocity:** $v_{\text{coll}} \approx 4500\;\text{km/s} \approx 0.01501\,c$, driving ram-pressure stripping that spatially separates the collisional X-ray intracluster medium (ICM) from the collisionless galaxies.

### 1.2. Four-Component Baryonic Model Along the Collision Axis ($x \in [-800, +800]\;\text{kpc}$)
Motivated by the characteristic mass and scale ratios from Chandra X-ray surface brightness maps and HST/VLT/Magellan optical imaging, the baryonic mass distribution along the collision axis ($y = 0$, $x \in [-800, +800]\;\text{kpc}$) is represented by four analytic components with aperture-normalized ($R_{\text{ap}} = 650.0\;\text{kpc}$) King / modified-Hubble ($\beta = 1$) surface density profiles $\Sigma_i(R) = \Sigma_{0,i}(1 + R^2/R_{c,i}^2)^{-1}$, where $\Sigma_{0,i} = \frac{M_i}{\pi R_{c,i}^2 \ln(1 + R_{\text{ap}}^2/R_{c,i}^2)}$:

| Component ID | Physical Description | Centroid $x_{c,i}$ ($\text{kpc}$) | Aperture Mass $M_i$ ($M_\odot$) | Core Radius $R_{c,i}$ ($\text{kpc}$) | Central $\Sigma_{0,i}$ ($M_\odot/\text{kpc}^2$) | Central $\kappa_{0,i} = \Sigma_{0,i}/\Sigma_{\text{crit}}$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `gal_main` | Main Cluster Collisionless Galaxies (Stars + Compact Remnants) | `-350.0` | `5.50e12` | `120.0` | `3.9329e7` | `0.012687` |
| `gas_main` | Main Cluster Collisional X-ray ICM Gas (Chandra, ram-pressure stripped) | `-170.0` | `2.60e13` | `150.0` | `1.3638e8` | `0.043995` |
| `gas_sub` | Subcluster ("Bullet") Shocked X-ray ICM Gas | `+140.0` | `6.50e12` | `75.0` | `8.4628e7` | `0.027300` |
| `gal_sub` | Subcluster ("Bullet") Collisionless Galaxies ($160\;\text{kpc}$ ahead of shock) | `+300.0` | `1.80e12` | `65.0` | `2.9296e7` | `0.009450` |

- **Total baryonic mass budget in the aperture:**
  $$M_{\text{gas,tot}} = 3.25\times 10^{13}\;M_\odot, \qquad M_{\text{gal,tot}} = 7.30\times 10^{12}\;M_\odot, \qquad M_{\text{bary,tot}} = 3.98\times 10^{13}\;M_\odot, \qquad \frac{M_{\text{gas,tot}}}{M_{\text{gal,tot}}} = 4.4521.$$
- **Observationally motivated peak positions and synthetic 1D reference curve (`kappa_synth_target_1D`):**
  In weak + strong lensing reconstructions (Clowe et al. 2006; Bradač et al. 2006), the convergence peaks lie near the collisionless galaxy centroids ($x_{\kappa,\text{main}}^{\text{obs}} \approx -350.0\;\text{kpc}$, $\kappa_{\text{peak,main}}^{\text{obs}} \approx 0.320$; $x_{\kappa,\text{sub}}^{\text{obs}} \approx +300.0\;\text{kpc}$, $\kappa_{\text{peak,sub}}^{\text{obs}} \approx 0.220$) rather than the dominant X-ray gas centroids ($x = -170\;\text{kpc}$ and $+140\;\text{kpc}$).
  **Methodological Scope and Caveat (V26):** The `kappa_synth_target_1D` column in `bullet_cluster_convergence_profile.csv` is **not an observational data series**, but a **self-constructed 1D synthetic reference curve**:
  $$\kappa_{\text{synth}}(x) = \kappa_{\text{GR,bary}}(x) + \frac{w_{\text{synth,main}}\,\Sigma_{\text{gal,main}}(x) + w_{\text{synth,sub}}\,\Sigma_{\text{gal,sub}}(x)}{\Sigma_{\text{crit}}}, \qquad (w_{\text{synth,main}} = 22.9765,\; w_{\text{synth,sub}} = 21.3095),$$
  built from the model's own four King profiles with two galaxy weights tuned to reproduce the illustrative peak heights at $x = -350\;\text{kpc}$ ($\kappa = 0.320$) and $x = +300\;\text{kpc}$ ($\kappa = 0.220$). Accordingly, the metric reported in `synth_mse_norm_vs_target`,
  $$M_{\text{synth}} \equiv \frac{1}{N}\sum_{i=1}^{N}\left[\frac{\kappa(x_i) - \kappa_{\text{synth}}(x_i)}{\sigma_{\text{ref}}}\right]^2 \qquad (\sigma_{\text{ref}} = 0.020,\; N = 321),$$
  is a **normalized mean-squared deviation from this synthetic 1D target curve**, **NOT** an observational $\chi^2_\nu$ fit to Clowe et al. (2006) or Bradač et al. (2006) data. A genuine observational $\chi^2$ fit requires a published 2D convergence map, full error covariance matrix, pre-registered spatial selection, and proper subtraction of fitted degrees of freedom ($N_{\text{dof}} = N_{\text{bins}} - p$).

---

## II. Test 1 — Pure Local Single-Scale Yukawa Model at Fixed $L = 10\;\text{km}$, $\alpha_Y = 1$ (Strict Cluster-Scale Falsification / No-Go Result)

We first evaluate the exact 3D scalar field $\phi_L^{\text{local}}(\mathbf{r})$ and projected 2D lensing convergence $\kappa_{\text{local}}(x)$ produced by the V24/V25 local Yukawa model at fixed $L = 10\;\text{km} = 3.24077929\times 10^{-16}\;\text{kpc}$ and $\alpha_Y = 1$.

### 2.1. Algebraic Contact Limit on Cluster Scales ($R_c \sim 65\text{–}150\;\text{kpc} \sim (2.0\text{–}4.6)\times 10^{17}\,L$)
In the V23–V25 sign convention ($\phi_L > 0$ near matter sources), the static 3D Helmholtz–Yukawa field equation is:
$$\left(\nabla^2 - \frac{1}{L^2}\right)\phi_L(\mathbf{r}) = -\alpha_{\text{std}}\,\frac{4\pi G}{c^2}\,\rho_b(\mathbf{r}), \qquad \alpha_{\text{std}} = \frac{\alpha_Y}{4\pi} \approx 0.0795775.$$
Because the cluster core radii satisfy $R_{c,i} \in [65, 150]\;\text{kpc}$, the dimensionless scale ratio is:
$$\frac{R_{c,i}}{L} \in \left[2.006\times 10^{17},\; 4.628\times 10^{17}\right] \implies \left(\frac{L}{R_{c,i}}\right)^2 \in \left[4.668\times 10^{-36},\; 2.486\times 10^{-35}\right].$$
Spatial derivatives ($|\nabla^2\phi_L| \sim R_c^{-2}|\phi_L|$) are suppressed by more than 35 orders of magnitude relative to the mass term $-L^{-2}\phi_L$, reducing the 3D scalar field to the algebraic contact limit:
$$\left|\phi_L^{\text{local}}(\mathbf{r})\right| \simeq \alpha_{\text{std}}\,L^2\,\frac{4\pi G}{c^2}\,\rho_b(\mathbf{r}) \in \left[10^{-44},\; 3.086\times 10^{-42}\right].$$
For individual compact objects and stars inside the galaxies, each star's $L = 10\;\text{km}$ Yukawa halo decays as $e^{-r/(10\;\text{km})}$, falling below $e^{-10} \approx 4.5\times 10^{-5}$ within $100\;\text{km}$ ($< 1.5\times 10^{-4}\,R_\odot$, or $< 3\times 10^{-18}\;\text{kpc}$) and vanishing as $\exp(-10^{17}) = 0$ on kiloparsec scales.

### 2.2. Conformal vs. Disformal Lensing Convergence at $L = 10\;\text{km}$
1. **Conformal Jordan Branch (`step14` Branch II-Conf, `step23`):**
   Under $\tilde{g}_{\mu\nu} = A^2(\phi_L)g_{\mu\nu} = e^{-2\alpha_Y\phi_L}g_{\mu\nu}$, the scalar field cancels identically out of the linear Weyl potential ($\Phi_{\text{lens}}^{\text{conf}} = \Phi_N$), while the exterior scalar energy cloud of compact objects ($\Delta M_{\phi,\text{ext}}^{\text{curved}}/M_{\text{surf}} = +1.55\times 10^{-3}$ for neutron stars, and $\lesssim 10^{-16}$ for main-sequence stars) is already absorbed into the stars' asymptotic ADM masses $M_{\text{gal}}$. Consequently:
   $$\Delta\kappa_{\text{local }\,L=10\text{ km}}^{\text{conf}}(x) = 0 \quad (\text{identically}).$$
2. **Disformal Refractive Clock Branch (`step14` Branch II-Disf):**
   If on a separate disformal optical branch $\Phi_{\text{lens}}^{\text{disf}} = \Phi_N - \frac{1}{2}c^2\phi_L$, integrating along the line of sight $z \in (-\infty, +\infty)$ and applying the 2D transverse Laplacian yields the 2D Fourier filter $\frac{k_\perp^2}{k_\perp^2 + L^{-2}} \simeq L^2 k_\perp^2 \leftrightarrow -L^2\nabla_\perp^2$, giving:
   $$\Delta\kappa_{\text{local }\,L=10\text{ km}}^{\text{disf}}(x) = -\frac{\alpha_{\text{std}}}{2}\,L^2\,\nabla_\perp^2 \kappa_{\text{GR,bary}}(x), \qquad \max_{x}\left|\frac{\Delta\kappa_{\text{local}}^{\text{disf}}}{\kappa_{\text{GR,bary}}}\right| = 2.171\times 10^{-36} \quad (\text{and } < 10^{-40}\text{ in the outer halo}).$$

### 2.3. Transparent Scientific Verdict for Test 1
On both metric coupling branches:
$$\kappa_{\text{local }\,L=10\text{ km}}(x) = \kappa_{\text{GR,bary}}(x) = \frac{\Sigma_{\text{gas}}(x) + \Sigma_{\text{gal}}(x)}{\Sigma_{\text{crit}}}.$$
- The convergence peaks remain locked to the **X-ray gas centroids**: in the 4-component global superposition they sit at $x_{\text{peak,main}}^{\text{local}} = -175.44\;\text{kpc}$ ($\kappa_{\text{peak,main}}^{\text{local}} = 0.04503$) and $x_{\text{peak,sub}}^{\text{local}} = +137.01\;\text{kpc}$ ($\kappa_{\text{peak,sub}}^{\text{local}} = 0.03697$), or in the shock-separated / steeper-outer-gas subcluster-pair decomposition at $x_{\text{peak,main}}^{\text{local}} = -170.7\;\text{kpc}$ and $x_{\text{peak,sub}}^{\text{local}} = +142.3\;\text{kpc}$ ($\kappa_{\text{peak,main}}^{\text{local}} = 0.0450$).
- The main peak amplitude is **$7.11\times$ smaller** than the $\kappa = 0.320$ reference peak, and the peaks are **displaced by $157.7\text{–}179.3\;\text{kpc}$** from the galaxy centroids ($-350\;\text{kpc}$ and $+300\;\text{kpc}$) toward the collisional X-ray gas ($M_{\text{synth}} = 32.089$).
- **Conclusion:** **At fixed $L = 10\;\text{km}$ and $\alpha_Y = 1$, the short-range Yukawa fifth force is purely a compact-object/neutron-star effect and produces zero macroscopic modification to Bullet Cluster lensing ($\Delta\kappa/\kappa_{\text{bary}} \le 2.17\times 10^{-36}$; $< 10^{-40}$ in the outer halo).**

---

## III. Test 2 — Phenomenological Macroscopic Extensions on the Synthetic 1D Test: What It Takes to Move the $\kappa(x)$ Peak from the X-ray Gas to the Galaxies

We now examine quantitatively within the synthetic 1D King model (`step11`, `step14` Branch II-Disf: $\Phi_{\text{lens}} = \Phi_N - \frac{1}{2}c^2\phi_{\text{macro}}$, with illustrative cluster Helmholtz scale $L_{\text{macro}} = 300.0\;\text{kpc}$) the two macroscopic sourcing cases.

### 3.1. Case 2A (`Macro_unscreened_bulk_baryon`) — Static Macroscopic Kernel Depending Only on Smoothed Bulk Baryon Density ($\bar{\rho}_b = \rho_{\text{gas}} + \rho_{\text{gal}}$)
Suppose a macroscopic scalar field is sourced universally by total baryon density with a uniform coupling $\eta_{\text{bulk}}$:
$$\left(\nabla^2 - L_{\text{macro}}^{-2}\right)\phi_{\text{macro}}(\mathbf{r}) = -\frac{4\pi G}{c^2}\,\eta_{\text{bulk}}\left[\rho_{\text{gas}}(\mathbf{r}) + \rho_{\text{gal}}(\mathbf{r})\right], \qquad \widetilde{\Sigma}_{\phi,\text{bulk}}(\mathbf{k}_\perp) = \eta_{\text{bulk}}\,\frac{k_\perp^2}{k_\perp^2 + L_{\text{macro}}^{-2}}\,\widetilde{\Sigma}_b(\mathbf{k}_\perp).$$
Calibrating $\eta_{\text{bulk}} = 8.6523$ boosts the main convergence peak by $7.11\times$ to match $\kappa_{\text{peak,main}} = 0.3200$. **However**, because $M_{\text{gas,tot}}/M_{\text{gal,tot}} = 4.4521$, the convergence peaks **remain locked to the X-ray gas**:
$$x_{\text{peak,main}}^{\text{bulk}} = -176.02\;\text{kpc} \quad (\text{offset from galaxies: } +173.98\;\text{kpc}), \qquad x_{\text{peak,sub}}^{\text{bulk}} = +138.35\;\text{kpc} \quad (\text{offset from galaxies: } -161.65\;\text{kpc}),$$
($M_{\text{synth}} = 20.723$). This demonstrates that **no static modified-gravity or scalar–tensor kernel sourced solely by smoothed bulk baryon density ($\rho_{\text{gas}} + \rho_{\text{gal}}$) can shift the convergence peaks from the dominant X-ray gas to the collisionless galaxies.**

### 3.2. Case 2B (`Macro_collisionless_compact_clock`) — Phenomenological Collisionless / Compact-Well Sourcing Hypothesis
For a macroscopic field to center its lensing convergence peaks on the collisionless galaxies without cold dark matter (CDM) particles, the macroscopic sector would need one of two physical properties:
1. **Hypothesis 2B-1 — Compactness / Deep-Potential-Well Clock Sourcing:** In the diffuse X-ray ICM plasma ($n_e \sim 10^{-3}\;\text{cm}^{-3}$), small-scale metric-gradient variance vanishes ($\langle (\nabla\nu)^2 \rangle_{\text{micro}} \approx 0$), whereas in the galaxy component baryons are bound inside stellar and compact-remnant potential wells. If the effective source of a macroscopic clock field is weighted by small-scale gradient energy or local compactness, it couples predominantly to the collisionless ($\sigma_{\text{coll}} \approx 0$) stellar/compact-remnant population.
2. **Hypothesis 2B-2 — Non-Equilibrium / Retarded Scalar Memory:** If the macroscopic field has a dynamical relaxation timescale $\tau_{\text{relax}} \gtrsim 150\;\text{Myr}$ longer than the core-passage time $t_{\text{cross}} = \Delta x / v_{\text{coll}} \approx 180\;\text{kpc} / (4500\;\text{km/s}) \approx 39.1\;\text{Myr}$, the scalar halo can travel ballistically with the collisionless galaxies during the supersonic merger.

### 3.3. Model-Dependent Threshold Calculation for $\eta_{\text{gal}} / \eta_{\text{gas}}$ and Methodological Caveats
With phenomenological macroscopic refractive coupling weights $\eta_{\text{gal}}$ and $\eta_{\text{gas}}$ in the chosen 1D King geometry ($R_{c,i} \in [65, 150]\;\text{kpc}$) at $L_{\text{macro}} = 300.0\;\text{kpc}$:
1. **Peak-Dominance Threshold ($\kappa(-350\;\text{kpc}) \ge \kappa(-170\;\text{kpc})$ and $\kappa(+300\;\text{kpc}) \ge \kappa(+140\;\text{kpc})$):**
   - Direct surface-density ratio (without Helmholtz filtering):
     $$\frac{\Sigma_{\text{gas,main}}(-170)}{\Sigma_{\text{gal,main}}(-350)} = 3.4592, \qquad \frac{\Sigma_{\text{gas,sub}}(+140)}{\Sigma_{\text{gal,sub}}(+300)} = 2.8895.$$
   - Including $L_{\text{macro}} = 300\;\text{kpc}$ Helmholtz–Yukawa filtering (at reference $\eta_{\text{gas}} = 1$):
     $$\left(\frac{\eta_{\text{gal}}}{\eta_{\text{gas}}}\right)_{\text{dom,main}} = 5.1204, \qquad \left(\frac{\eta_{\text{gal}}}{\eta_{\text{gas}}}\right)_{\text{dom,sub}} = 5.0972.$$
2. **$< 5.0\;\text{kpc}$ Centroid-Offset Threshold ($|x_{\text{peak,main}} - (-350\;\text{kpc})| \le 5.0\;\text{kpc}$):**
   - To prevent the X-ray gas tail in this 1D King profile from pulling the main convergence peak by more than $5.0\;\text{kpc}$ away from the galaxy centroid at $-350.0\;\text{kpc}$, the required ratio is:
     $$\left(\frac{\eta_{\text{gal}}}{\eta_{\text{gas}}}\right)_{|\Delta x| \le 5\;\text{kpc}} \ge 25.6650.$$
3. **Calibrated Case 2B Illustrative Model (`Macro_collisionless_compact_clock`, $\eta_{\text{gas}} = 0$, $\eta_{\text{gal,main}} = 31.7448$, $\eta_{\text{gal,sub}} = 24.1556$):**
   - Main cluster peak: $x_{\text{peak,main}} = -347.86\;\text{kpc}$ (centroid offset from galaxies: **$+2.14\;\text{kpc} < 5\;\text{kpc}$**), $\kappa_{\text{peak,main}} = 0.32011$ ($\kappa(-350\;\text{kpc}) = 0.32000$).
   - Subcluster peak: $x_{\text{peak,sub}} = +299.45\;\text{kpc}$ (centroid offset from galaxies: **$-0.55\;\text{kpc} < 5\;\text{kpc}$**), $\kappa_{\text{peak,sub}} = 0.22002$ ($\kappa(+300\;\text{kpc}) = 0.22000$).
   - Deviation from the synthetic 1D target curve (`kappa_synth_target_1D`): $\text{RMS}(\Delta\kappa) = 0.02559$, normalized synthetic mean-squared deviation $M_{\text{synth}} = 1.637$ (compared with $M_{\text{synth}} = 32.089$ for the local $L = 10\;\text{km}$ model and $M_{\text{synth}} = 20.723$ for the Case 2A bulk-baryon model).
4. **Strict Theoretical and Methodological Delineation (V26):**
   - The peak-dominance ratio $\eta_{\text{gal}}/\eta_{\text{gas}} \ge 5.1204$ and the $< 5\;\text{kpc}$ centroid-offset ratio $\eta_{\text{gal}}/\eta_{\text{gas}} \ge 25.6650$ are **model-dependent thresholds specific to the chosen 1D King geometry ($R_{c,i} = 65\text{–}150\;\text{kpc}$) and Helmholtz scale $L_{\text{macro}} = 300\;\text{kpc}$**, not a general model-independent physical theorem.
   - Setting $\eta_{\text{gas}} = 0$ and fitting two separate galactic weights ($\eta_{\text{gal,main}} = 31.7448 \neq \eta_{\text{gal,sub}} = 24.1556$) are **phenomenological, illustrative assumptions that are not derived from the covariant scalar–tensor action**.

---

## IV. Synthetic 1D Illustrative Peak-Offset and Convergence Profile Tables

### 4.1. Summary of Peak Positions, Amplitudes, and Deviations from the Synthetic 1D Target Curve Across All Four Models (`bullet_cluster_peak_offsets_summary.csv` — Synthetic 1D Illustrative Test)

| Model (`model_name`) | $x_{\text{peak,main}}$ ($\text{kpc}$) | $\kappa_{\text{peak,main}}$ | $x_{\text{peak,sub}}$ ($\text{kpc}$) | $\kappa_{\text{peak,sub}}$ | $\Delta x_{\text{main,gal}}$ ($\text{kpc}$) | $\Delta x_{\text{main,gas}}$ ($\text{kpc}$) | $\Delta x_{\text{sub,gal}}$ ($\text{kpc}$) | $\Delta x_{\text{sub,gas}}$ ($\text{kpc}$) | RMS $\Delta\kappa_{\text{synth}}$ | $M_{\text{synth}}$ (`synth_mse_norm_vs_target`) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GR_baryon_only` | `-175.44` | `0.04503` | `+137.01` | `0.03697` | `+174.56` | `  -5.44` | `-162.99` | `  -2.99` | `0.11329` | `32.089` |
| `TAL_local_Yukawa_L10km_alphaY1` | `-175.44` | `0.04503` | `+137.01` | `0.03697` | `+174.56` | `  -5.44` | `-162.99` | `  -2.99` | `0.11329` | `32.089` |
| `Macro_unscreened_bulk_baryon` | `-176.02` | `0.32002` | `+138.35` | `0.26531` | `+173.98` | `  -6.02` | `-161.65` | `  -1.65` | `0.09104` | `20.723` |
| `Macro_collisionless_compact_clock` | `-347.86` | `0.32011` | `+299.45` | `0.22002` | `  +2.14` | `-177.86` | `  -0.55` | `+159.45` | `0.02559` | `1.637` |

*(Note: $M_{\text{synth}} = \frac{1}{N}\sum_{i=1}^N [(\kappa(x_i) - \kappa_{\text{synth}}(x_i))/\sigma_{\text{ref}}]^2$ is the normalized mean-squared deviation from the self-constructed 1D `kappa_synth_target_1D` reference curve with $\sigma_{\text{ref}} = 0.020$, not an observational $\chi^2_\nu$.)*

### 4.2. Projected Surface Density and $\kappa(x)$ Profile at Representative Collision-Axis Positions (`bullet_cluster_convergence_profile.csv` — Synthetic 1D Illustrative Test)

| $x$ ($\text{kpc}$) | $\Sigma_{\text{gas}}$ ($M_\odot/\text{kpc}^2$) | $\Sigma_{\text{gal}}$ ($M_\odot/\text{kpc}^2$) | $\kappa_{\text{GR,bary}}$ | $\kappa_{\text{local }\,L=10\text{ km}}$ | $\Delta\kappa_{\text{local,disf}}$ | $\kappa_{\text{Macro,2A (bulk)}}$ | $\kappa_{\text{Macro,2B (clock)}}$ | $\kappa_{\text{synth,1D}}$ (`kappa_synth_target_1D`) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `-600.0` | `1.4233e+07` | `6.8239e+06` | `0.006793` | `0.006793` | `-5.740e-40` | `0.010277` | `0.028049` | `0.062017` |
| `-450.0` | `2.8832e+07` | `2.1245e+07` | `0.016154` | `0.016154` | `+2.251e-42` | `0.076874` | `0.165259` | `0.188492` |
| `-350.0` | `5.2453e+07` | `3.5918e+07` | `0.028507` | `0.028507` | `+1.240e-38` | `0.178665` | `0.320000` | `0.320000` |
| `-260.0` | `9.3503e+07` | `2.3192e+07` | `0.037644` | `0.037644` | `+8.951e-39` | `0.254340` | `0.203216` | `0.225535` |
| `-175.0` | `1.2766e+08` | `1.1934e+07` | `0.045031` | `0.045031` | `+2.869e-38` | `0.320000` | `0.104232` | `0.141178` |
| `-170.0` | `1.2794e+08` | `1.1514e+07` | `0.044984` | `0.044984` | `+2.882e-38` | `0.319487` | `0.100405` | `0.137693` |
| ` -80.0` | `9.9460e+07` | `6.7119e+06` | `0.034249` | `0.034249` | `+6.544e-39` | `0.219456` | `0.048876` | `0.087474` |
| `  +0.0` | `7.2882e+07` | `5.0651e+06` | `0.025144` | `0.025144` | `-3.132e-39` | `0.138610` | `0.027722` | `0.064214` |
| ` +70.0` | `7.9995e+07` | `4.8621e+06` | `0.027373` | `0.027373` | `+3.408e-40` | `0.166192` | `0.029235` | `0.063429` |
| `+140.0` | `1.0829e+08` | `6.1784e+06` | `0.036924` | `0.036924` | `+8.015e-38` | `0.265185` | `0.049209` | `0.080507` |
| `+220.0` | `5.5598e+07` | `1.3197e+07` | `0.022192` | `0.022192` | `-2.748e-39` | `0.134507` | `0.088322` | `0.110826` |
| `+300.0` | `2.6688e+07` | `3.0558e+07` | `0.018466` | `0.018466` | `+3.559e-38` | `0.110815` | `0.220000` | `0.220000` |
| `+450.0` | `1.1510e+07` | `5.4294e+06` | `0.005464` | `0.005464` | `-9.488e-40` | `0.006652` | `0.021648` | `0.042188` |
| `+600.0` | `6.7045e+06` | `1.8771e+06` | `0.002768` | `0.002768` | `-1.602e-40` | `-0.004653` | `0.000158` | `0.015925` |

---

## V. Reproducibility and Generated Artifacts

- **Numerical audit script:** `bullet_cluster_lensing_audit.py` (synced to `code/`).
- **1D collision-axis synthetic convergence profile (321 grid points, $x \in [-800, +800]\;\text{kpc}$):** `bullet_cluster_convergence_profile.csv` (synced to root, `code/`, `csv/`, and `data/v23/`).
- **Peak-offset and normalized synthetic deviation ($M_{\text{synth}}$) summary table:** `bullet_cluster_peak_offsets_summary.csv` (synced to root, `code/`, `csv/`, and `data/v23/`).
- **Structured JSON summary:** `bullet_cluster_audit_summary.json` (synced to root, `code/`, `csv/`, and `data/v23/`).
