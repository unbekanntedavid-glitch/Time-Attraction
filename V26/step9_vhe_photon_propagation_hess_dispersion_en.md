# Step 9: High-Energy Photon Propagation in the Disformal/Conformal Scalar-Tensor Metric, the Effective Refractive Index $n_t(s) = f_L(s)^{-1}$, and the Multi-Source V17 H.E.S.S./MAGIC Gamma-Ray Timing Audit

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step9_vhe_photon_propagation_hess_dispersion_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step8_microscopic_clock_coupling_siriusb_lab_v20_en.md` (`step8_microscopic_clock_coupling_siriusb_lab_v20.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p2.tex` (Chapter 6, lines 63–99, Eq. 9 and the 5-source H.E.S.S./MAGIC control table; DOI: `10.5281/zenodo.23604423`)  
- V17 Zenodo record: DOI `10.5281/zenodo.19666220` (SPARC 155 galaxies and multi-source H.E.S.S. VHE gamma-ray timing audit)  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Sections 1 and 3, Eq. 1 and lines 161–166; DOI: `10.5281/zenodo.23639084`)

---

## 1. Verbatim Citations from Chapter 6 of the V22 Monograph (`hu_p2.tex`)

In Step 8 (`step8_microscopic_clock_coupling_siriusb_lab_v20.md`), we derived the stationary proper-time and redshift coupling of a comoving atomic clock ($\Psi_{\text{atomic}} = (\psi_e, A_\mu)$) on the effective metric $\tilde{g}_{\mu\nu}$. The next logical step is to derive the **propagation of the Maxwell photon field ($F_{\mu\nu}$) across astronomical and cosmological baseline distances ($s \in [0, D]$)** from the scalar–tensor / disformal action principle, comparing it with the multi-source V17 H.E.S.S./MAGIC gamma-ray timing audit documented in Chapter 6 of the V22 monograph:

1. **Local Refractive Index of the Temporal Carrier ($n_t(s) = f_L(s)^{-1}$), Achromatic Delay, and Dispersion Parameterization (`hu_p2.tex`, Chapter 6, lines 68–74, Eq. 9):**
   > *"When a VHE gamma-ray photon traverses cosmic voids, filaments, host-galaxy halos, and galactic interstellar structures along a line of sight $s \in [0, D]$, the local refractive index of the temporal carrier is $n_t(s) = f_L(s)^{-1} = \exp(+\gamma_L \bar{\rho}_L(s))$. In the purely metric (geometric) limit, where all photon energies couple identically to the effective lapse factor $f_L(s)$, the arrival time experiences an achromatic Shapiro-like excess delay: $\Delta t_{\text{geom}} = c^{-1}\int_0^D \bigl[\exp(\gamma_L\bar{\rho}_L(s)) - 1\bigr]\,ds$. If, however, high-energy ($E$) photon wave packets undergo finite-wavelength dispersion by scattering across spatial gradients $\nabla\bar{\rho}_L$ on scale $L$, the parameterized form of the energy-dependent time lag between energy channels $E_1$ and $E_2$ is:"*
   $$\begin{equation}
     \Delta t(E_1, E_2) \;=\; \tau_n\,\frac{E_2^n - E_1^n}{E_{\text{ref}}^n}\int_0^D \gamma_L\,\bar{\rho}_L(s)\,(1+z(s))^n\,ds, \qquad n \in \{1, 2\}, \tag{V22 Eq. 9}
   \end{equation}$$
   > *"where $\tau_n$ is the dispersion coupling coefficient and $D(z)$ is the redshift-weighted comoving path length."*

2. **The Five-Source V17 H.E.S.S. / MAGIC Control Table (`hu_p2.tex`, Chapter 6, lines 78–91):**

   | VHE Source | Redshift / Dist. | Energy Band ($\text{TeV}$) | Observed Lag / Bound | V17 Audit Result |
   | :--- | :--- | :--- | :--- | :--- |
   | **PKS 2155-304 (H.E.S.S.)** | $z = 0.116$ | $0.2\text{--}4.0\,\text{TeV}$ | $\Delta t \approx +28 \pm 12\,\text{s}$ | Marginal $\sim 2\sigma$ lag; dominated by intrinsic SSC emission profile. |
   | **Mrk 501 (1997 / H.E.S.S.)** | $z = 0.034$ | $0.15\text{--}10\,\text{TeV}$ | $\Delta t \approx +240 \pm 60\,\text{s}$ | Larger lag at $3.4\times$ lower $z$ violates monotonic distance scaling! |
   | **GRB 190114C (MAGIC)** | $z = 0.4245$ | $0.2\text{--}1.0\,\text{TeV}$ | $|\Delta t / \Delta E| \le 19\,\text{s/TeV}$ | Stringent cosmological bound ruling out large linear vacuum dispersion. |
   | **Crab Pulsar (H.E.S.S./MAGIC)** | $d \approx 2.0\,\text{kpc}$ | $0.1\text{--}1.5\,\text{TeV}$ | $|\Delta t| < 0.4\,\text{ms}$ | Tight Galactic-disk bound across dense ISM ($\bar{n}_H \sim 1\,\text{cm}^{-3}$). |
   | **Vela Pulsar (H.E.S.S. II)** | $d \approx 0.287\,\text{kpc}$ | $0.02\text{--}0.1\,\text{TeV}$ | $|\Delta t| < 0.2\,\text{ms}$ | Sub-millisecond phase alignment across the dense Vela SNR shell. |

3. **Two Main Conclusions of the V17 Multi-Source Audit (`hu_p2.tex`, Chapter 6, lines 95–98):**
   > - **Rejection of Naive Distance-Proportional Dispersion:** *"Because Mrk 501 ($z=0.034$) exhibits a larger apparent flare lag than PKS 2155-304 ($z=0.116$) while GRB 190114C ($z=0.4245$) and the Crab/Vela pulsars yield tight null bounds, explaining the $28\,\text{s}$ lag of PKS 2155-304 via a uniform intergalactic dispersion constant is refuted. Intrinsic blazar synchrotron self-Compton (SSC) acceleration timescales must be modeled first."*
   > - **Unique Environmental Overdensity Prediction vs. $\Lambda$CDM:** *"Unlike Lorentz-invariance-violating (LIV) models that predict a purely redshift-dependent delay $\Delta t \propto \int_0^z \frac{(1+z')\,dz'}{H(z')}$, the Time Attraction Law predicts that any genuine non-geometric timing residual must correlate with the baryonic overdensity integral ($\int \bar{\rho}_L(s)\,ds$) along lines of sight passing through rich galaxy clusters (such as Coma or Perseus) compared to cosmic voids at the same redshift."*

---

## 2. Mathematical Derivation: Light Cone, Effective Refractive Index, and EFT Dispersion in the Scalar-Tensor Framework

### 2.1. Conformal Invariance vs. Disformal Light Cone: Derivation of the Refractive Index $n_t(s) = f_L(s)^{-1}$
Let us examine the free Maxwell action ($S_{\text{EM}} = -\frac{1}{4\mu_0}\int d^4x\sqrt{-\tilde{g}}\,\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}$) under the two metric transformations derived in Step 2 (`step2_conformal_disformal.md`):

1. **Pure Conformal Transformation ($\tilde{g}_{\mu\nu}^{\text{conf}} = f_L^2 g_{\mu\nu}$):**
   In $D = 4$ spacetime dimensions, $\sqrt{-\tilde{g}^{\text{conf}}} = f_L^4\sqrt{-g}$ and $\tilde{g}_{\text{conf}}^{\mu\alpha}\tilde{g}_{\text{conf}}^{\nu\beta} = f_L^{-4}g^{\mu\alpha}g^{\nu\beta}$, so the factor $f_L^4 f_L^{-4} = 1$ cancels:
   $$S_{\text{EM}}\!\left[f_L^2 g_{\mu\nu}, A_\mu\right] \;=\; S_{\text{EM}}\!\left[g_{\mu\nu}, A_\mu\right].$$
   Therefore, in a purely conformal scalar–tensor coupling, the light cone ($\tilde{g}_{\mu\nu}^{\text{conf}}k^\mu k^\nu = 0 \iff g_{\mu\nu}k^\mu k^\nu = 0$) is identical to that of the background metric.

2. **Medium-Anchored Disformal Clock Metric of V23 ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$):**
   If electromagnetic waves couple to the disformal metric $\tilde{g}_{\mu\nu}^{\text{disf}}$ defined by the four-velocity $u^\mu$ of the baryonic medium, then in the local rest frame ($g_{\mu\nu} \simeq \eta_{\mu\nu}$, $u^\mu = (c, \mathbf{0})$) the line element is:
   $$d\tilde{s}_{\text{disf}}^2 \;=\; -f_L(s)^2 c^2 dt^2 + ds^2, \qquad ds^2 = dx^2 + dy^2 + dz^2.$$
   In the geometric-optics (eikonal / null-geodesic) limit, $d\tilde{s}_{\text{disf}}^2 = 0$, from which the coordinate speed of light of the photon and the effective temporal refractive index of the medium follow:
   $$v_{\text{eff}}(s) \;=\; \frac{ds}{dt} \;=\; c\,f_L(s) \;=\; c\,\exp\!\bigl(-\phi_L(s)\bigr) \quad\Longrightarrow\quad \boxed{\;n_t(s) \;\equiv\; \frac{c}{v_{\text{eff}}(s)} \;=\; f_L(s)^{-1} \;=\; \exp\!\bigl(+\phi_L(s)\bigr) \;=\; \exp\!\bigl(+\gamma_L\,\bar{\rho}_L(s)\bigr).\;}$$
   Integrating over the path $s \in [0, D]$ directly yields the achromatic Shapiro-like travel-time excess given in line 69 of V22 `hu_p2.tex`:
   $$\Delta t_{\text{geom}} \;=\; \int_0^D \left(\frac{1}{v_{\text{eff}}(s)} - \frac{1}{c}\right)ds \;=\; \frac{1}{c}\int_0^D \Bigl[\exp\!\bigl(\gamma_L\bar{\rho}_L(s)\bigr) - 1\Bigr]\,ds \;\simeq\; \frac{1}{c}\int_0^D \phi_L(s)\,ds \qquad (\phi_L \ll 1).$$
   Because $\Delta t_{\text{geom}}$ does not depend on the photon energy $E = \hbar\omega$, in the pure eikonal limit the relative delay between photons of energies $E_1$ and $E_2$ emitted simultaneously vanishes: $\Delta t_{\text{geom}}(E_1, E_2) = 0$.

### 2.2. Finite-Wavelength Scalar–Photon EFT Operators and Derivation of the Dispersion Formula V22 Eq. (9)
When the wavelength of high-energy (VHE, $E \sim 0.1\text{--}10\;\text{TeV}$) gamma-ray photons interacts with gradients of the scalar field $\phi_L = \gamma_L\bar{\rho}_L$ (or $\phi_L^{\text{Yuk}}$), the leading dimension-$4+n$ ($n \in \{1,2\}$) scalar–Maxwell operator in the effective field theory (EFT) expansion with the medium four-velocity $u^\mu$ is:
$$\Delta\mathcal{L}_{\text{disp}}^{(n)} \;=\; -\frac{c\,\tau_n}{4\mu_0}\,\frac{\phi_L(x)}{E_{\text{ref}}^n}\,F_{\mu\nu}\,\bigl(i\hbar\,u^\alpha\nabla_\alpha\bigr)^n F^{\mu\nu}.$$
In the plane-wave approximation ($A_\mu \propto e^{i k_\alpha x^\alpha}$, where the local comoving photon energy is $E_{\text{loc}}(s) = -\hbar u^\alpha k_\alpha = E_{\text{obs}}\,(1 + z(s))$), the modified dispersion relation is:
$$\omega^2 \;=\; c^2 k^2\left[1 - \frac{2c\,\tau_n}{n+1}\,\phi_L(s)\,\left(\frac{E_{\text{obs}}(1+z(s))}{E_{\text{ref}}}\right)^n\right].$$
Integrating the reciprocal of the group velocity ($v_g = \partial\omega/\partial k$) between two distinct observed energy channels ($E_1 < E_2$) yields V22 Eq. (9):
$$\boxed{\;\Delta t(E_1, E_2) \;=\; t_{\text{arr}}(E_2) - t_{\text{arr}}(E_1) \;=\; \tau_n\,\frac{E_2^n - E_1^n}{E_{\text{ref}}^n}\int_0^D \phi_L(s)\,\bigl(1 + z(s)\bigr)^n\,ds \;=\; \tau_n\,\frac{E_2^n - E_1^n}{E_{\text{ref}}^n}\int_0^D \gamma_L\,\bar{\rho}_L(s)\,\bigl(1 + z(s)\bigr)^n\,ds.\;}$$

---

## 3. Checks and Physical Implications

1. **Why Does the 5-Source V17 Audit Exclude Universal Vacuum Dispersion ($\phi_L = \text{const.}$) and Support Local $\phi_L(\mathbf{x})$ Coupling?**
   - If dispersion were caused by density-independent vacuum LIV ($\phi_L(s) \equiv \text{const.}$), then $\Delta t / \Delta E$ would increase strictly monotonically with redshift $z$. Instead, **Mrk 501** ($z = 0.034$, $\Delta t \approx +240 \pm 60\;\text{s}$) shows an $\sim 8.6\times$ larger lag at $3.4\times$ lower redshift than **PKS 2155-304** ($z = 0.116$, $\Delta t \approx +28 \pm 12\;\text{s}$), while **GRB 190114C** at $z = 0.4245$ gives $|\Delta t/\Delta E| \le 19\;\text{s/TeV}$ and the Galactic **Crab** ($|\Delta t| < 0.4\;\text{ms}$) / **Vela** ($|\Delta t| < 0.2\;\text{ms}$) pulsars impose tight upper bounds $\checkmark$.
2. **Discriminating Prediction of Scalar–Tensor / Yukawa TAL Against LIV Models:**
   - Because the scalar field $\phi_L(s)$ is proportional to the baryonic mass density along the line of sight ($\phi_L = \gamma_L\bar{\rho}_L$ or $\phi_L^{\text{Yuk}}$), for sources at identical redshift $z$ the line-of-sight integral $\mathcal{I}_L \equiv \int_0^D \phi_L(s)(1+z(s))^n\,ds$ is proportional to the baryonic overdensity of intervening galaxy clusters (e.g., Coma, Perseus) relative to cosmic voids $\checkmark$.
3. **Vacuum Limit ($\bar{\rho}_L \to 0 \implies \phi_L \to 0$):**
   - $\phi_L(s) \to 0 \implies n_t(s) = e^0 = 1$, $\Delta t_{\text{geom}} = 0$, and $\Delta t(E_1, E_2) = 0$: light propagation recovers the special/general relativistic vacuum light cone $\checkmark$.
