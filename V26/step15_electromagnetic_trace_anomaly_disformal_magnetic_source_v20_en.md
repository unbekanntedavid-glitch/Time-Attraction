# Step 15: Electromagnetic Trace-Freeness ($T^\mu_{\;\mu,\text{EM}} = 0$) vs. Disformal Magnetic Source Coupling ($\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$), V20 Polarity-Even Cross-Term Cancellation, and the Methodological Boundary ($\Phi_t \equiv 0$)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step14_gravitational_lensing_null_geodesics_cosmological_background_en.md` (`step14_gravitational_lensing_null_geodesics_cosmological_background.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (Chapter 2, points 3–4, lines 208–221) and `hu_p3.tex` (Chapter 10, point 2, lines 75–85; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Sections 1 and 3, lines 33–35, 58–66, and 161–166; DOI: `10.5281/zenodo.23639084`)  
- `README.md` (`scalar-tensor-yukawa-fit/README.md`, Section 3, point 5, lines 60–61)

---

## 1. Verbatim Citations from the V22 (`hu_p1.tex`, `hu_p3.tex`) and V23 Manuscripts

In Step 8 (`step8_microscopic_clock_coupling_siriusb_lab_v20.md`), we derived how a given scalar profile $\phi_L$ shifts the frequencies of atomic clocks and quartz oscillators (H0 vs. H1), and evaluated the numerical value for the V20 8-coil solenoid ($\Delta y_{\text{peak}} = -1.2017108\times 10^{-48}$ at $L = 0.630957\;\text{m}$). Within the scalar–tensor framework, however, there is a fundamental theoretical question highlighted in `hu_p1.tex` (lines 208–221), `hu_p3.tex` (lines 75–85), and `README.md` (line 61): **how can pure magnetic energy density ($\rho_B = \frac{|\mathbf{B}|^2}{2\mu_0 c^2}$) act as a source for the scalar field $\phi_L$ in a covariant action principle when the trace of the 4-dimensional Maxwell stress–energy tensor vanishes identically ($T^\mu_{\;\mu,\text{EM}} = g^{\mu\nu}T_{\mu\nu}^{\text{EM}} = 0$)?**

1. **Methodological Boundary Between the Audited Physical Core ($\rho = T_{00}/c^2$, $\Phi_t \equiv 0$) and the V13–V16 Perspective Layer (`hu_p1.tex`, Chapter 2, points 3–4, lines 209–221):**
   $$\begin{equation}
     \rho_{\text{eff}} \;=\; \rho \;+\; \Xi\,\varepsilon_0\,\Phi_t \;+\; \beta\,\Phi_B \;+\; \gamma\,\kappa, \qquad
     T_{\mu\nu}^{\text{total}} \;=\; T_{\mu\nu}^{\text{matter}} \;+\; T_{\mu\nu}^{\text{flux}} \;+\; T_{\mu\nu}^{\text{thought}}, \tag{V22 Eq. 2}
   \end{equation}$$
   > *"When presenting this book to the academic community, I draw a clear boundary between the two layers:*  
   > - ***Audited Physical Core (V19–V21):** Every numerical result in the SPARC galactic rotation analysis, the H.E.S.S. timing bounds, the neutron-star TOV and pulsation sequences, and the V20 granite-chamber and 8-coil quartz-oscillator laboratory protocol is built exclusively on the standard baryonic and electromagnetic mass–energy density $\rho = T_{00}/c^2$, with $\Phi_t \equiv 0$.*  
   > - ***Perspective Bio-Coherence Layer (V13–V16):** The conscious-coherence terms ($\Phi_t$, $T_{\mu\nu}^{\text{thought}}$, EEG/HRV bio-coherence) are retained as an open, not-yet-audited theoretical perspective representing part of the intellectual history of the Time Attraction Law."*

2. **V20 Eight-Coil Magnetic Energy-Density Source and Background Cross-Term Cancellation (`hu_p3.tex`, Chapter 10, point 2, lines 76 and 84):**
   > - *"To test whether pure magnetic energy density $u_B(\mathbf{x}) = \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0}$ slows the flow of time via $\rho_B(\mathbf{x}) = u_B(\mathbf{x})/c^2$, I designed a symmetric eight-coil arrangement..."*
   > - *"**Cancellation of the Background Cross-Term via Polarity Reversal:** Since $|\mathbf{B}_{\text{tot}}|^2 = |\mathbf{B}_{\text{bg}}|^2 + 2\mathbf{B}_{\text{bg}}\cdot\mathbf{B}_{\text{coil}} + |\mathbf{B}_{\text{coil}}|^2$, reversing the current direction ($+I \leftrightarrow -I$) reverses $\mathbf{B}_{\text{coil}}$, so both the linear Zeeman sensitivity of the quartz oscillator and the cross-term $2\mathbf{B}_{\text{bg}}\cdot\mathbf{B}_{\text{coil}}$ cancel out of the polarity-even average $\Delta y_{\text{even}}$."*

---

## 2. Why Does a Magnetic Field Not Source a Purely Conformal Brans–Dicke / Dilaton Scalar, and Why Does It Source the V23 Disformal Coupling?

### 2.1. Electromagnetic Blindness of Purely Conformal Coupling ($T^\mu_{\;\mu,\text{EM}} = 0$ in $D=4$)
Consider the standard Maxwell action of the free electromagnetic field ($F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$) on a general Jordan metric $\tilde{g}_{\mu\nu}$:
$$S_{\text{EM}}\bigl[\tilde{g}_{\mu\nu}, A_\mu\bigr] \;=\; -\frac{1}{4\mu_0}\int d^4x\,\sqrt{-\tilde{g}}\,\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}.$$
If the metric depends on the scalar field $\phi_L$ in a purely **conformal** manner ($\tilde{g}_{\mu\nu}^{\text{conf}} = \Omega^2(\phi_L)\,g_{\mu\nu} = f_L^2\,g_{\mu\nu}$, so that $\tilde{g}_{\text{conf}}^{\mu\alpha} = f_L^{-2}g^{\mu\alpha}$ and in $D=4$ dimensions $\sqrt{-\tilde{g}^{\text{conf}}} = f_L^4\sqrt{-g}$), then the factors of $f_L$ cancel identically:
$$\sqrt{-\tilde{g}^{\text{conf}}}\,\tilde{g}_{\text{conf}}^{\mu\alpha}\tilde{g}_{\text{conf}}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} \;=\; \bigl(f_L^4\sqrt{-g}\bigr)\bigl(f_L^{-2}g^{\mu\alpha}\bigr)\bigl(f_L^{-2}g^{\nu\beta}\bigr)F_{\mu\nu}F_{\alpha\beta} \;=\; \sqrt{-g}\,g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}.$$
Accordingly, the variation of the Maxwell action with respect to $\phi_L$ vanishes identically in purely conformal coupling (which is equivalent to the vanishing of the trace of the Maxwell stress–energy tensor in $D=4$, $T^\mu_{\;\mu,\text{EM}} = g^{\mu\nu}T_{\mu\nu}^{\text{EM}} = -u_E - u_B + 3\cdot\frac{1}{3}(u_E + u_B) = 0$):
$$\boxed{\;\frac{\delta S_{\text{EM}}\bigl[\tilde{g}_{\mu\nu}^{\text{conf}}, A_\mu\bigr]}{\delta\phi_L(x)} \;=\; \frac{\partial\tilde{g}_{\mu\nu}^{\text{conf}}}{\partial\phi_L}\,\frac{\delta S_{\text{EM}}}{\delta\tilde{g}_{\mu\nu}^{\text{conf}}} \;=\; \sqrt{-g}\,f_L^2\,\frac{\partial\ln f_L}{\partial\phi_L}\,g^{\mu\nu}T_{\mu\nu}^{\text{EM}} \;=\; 0.\;}$$
**Physical consequence:** A purely conformal Jordan–Brans–Dicke scalar field is not sourced at tree level by a static magnetic field $\mathbf{B}$ at all!

### 2.2. Exact Magnetic Source Term of the V23 Disformal Clock Metric ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$)
Now let us substitute the V23 **temporal disformal metric** (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, line 165):
$$\tilde{g}_{\mu\nu}^{\text{disf}} \;=\; g_{\mu\nu} \;+\; \frac{1 - f_L^2(\phi_L)}{c^2}\,u_\mu u_\nu, \qquad \tilde{g}_{\text{disf}}^{\mu\nu} \;=\; g^{\mu\nu} \;-\; \frac{1 - f_L^{-2}(\phi_L)}{c^2}\,u^\mu u^\nu, \qquad \sqrt{-\tilde{g}^{\text{disf}}} \;=\; f_L(\phi_L)\,\sqrt{-g}.$$
Decomposing the field-strength tensor $F_{\mu\nu}$ into electric ($E_\mu = F_{\mu\nu}u^\nu/c$) and magnetic ($B_\mu = \frac{1}{2}\epsilon_{\mu\nu\alpha\beta}u^\nu F^{\alpha\beta}/c$) field four-vectors comoving with the laboratory four-velocity $u^\mu$:
$$\frac{1}{4}\tilde{g}_{\text{disf}}^{\mu\alpha}\tilde{g}_{\text{disf}}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} \;=\; \frac{1}{2}\left[\tilde{g}_{\text{disf}}^{ij}\tilde{g}_{\text{disf}}^{kl}F_{ik}F_{jl} \;+\; 2\,\tilde{g}_{\text{disf}}^{00}\tilde{g}_{\text{disf}}^{ij}F_{0i}F_{0j}\right] \;=\; \frac{|\mathbf{B}|^2}{2\mu_0}\cdot\mu_0 \;-\; f_L^{-2}(\phi_L)\,\frac{\varepsilon_0|\mathbf{E}|^2}{2}\cdot\mu_0.$$
Substituting this together with the volume element $\sqrt{-\tilde{g}^{\text{disf}}} = f_L(\phi_L)\sqrt{-g}$ into the Maxwell action gives:
$$\boxed{\;S_{\text{EM}}\bigl[\tilde{g}_{\mu\nu}^{\text{disf}}, A_\mu\bigr] \;=\; \int d^4x\,\sqrt{-g}\left[-\,f_L(\phi_L)\,\frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0} \;+\; f_L(\phi_L)^{-1}\,\frac{\varepsilon_0|\mathbf{E}(\mathbf{x})|^2}{2}\right].\;}$$
In a static magnetic field ($\mathbf{E} = \mathbf{0}$, as in the V20 8-coil solenoid), the Lagrangian density $\mathcal{L}_{\text{EM}}^{\text{disf}}(\mathbf{E}=\mathbf{0}) = -f_L(\phi_L)\,u_B(\mathbf{x})$ takes the exact same $-f_L(\phi_L)\varepsilon$ form as a static baryonic fluid (`step4_covariant_conservation_hydrostatic.md`: $\mathcal{L}_{\text{matter}} = -f_L(\phi_L)\varepsilon_{\text{baryon}}$)! Varying the V23 Formulation II (Helmholtz–Yukawa) or Formulation I (bilateral top-hat) action with respect to $\phi_L$ ($f_L = e^{-\phi_L}$, $\partial f_L/\partial\phi_L = -f_L$) yields:
$$\frac{\partial\bigl(\mathcal{L}_{\text{matter}} + \mathcal{L}_{\text{EM}}^{\text{disf}}\bigr)}{\partial\phi_L} \;=\; f_L(\phi_L)\,\Bigl(\varepsilon_{\text{baryon}}(\mathbf{x}) + u_B(\mathbf{x})\Bigr) \;=\; f_L(\phi_L)\,c^2\Bigl(\rho_{\text{baryon}}(\mathbf{x}) + \rho_B(\mathbf{x})\Bigr), \qquad \rho_B(\mathbf{x}) \;\equiv\; \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0 c^2}.$$
We have thus derived from a covariant action principle that the V23 disformal clock metric ($\sqrt{-\tilde{g}} = f_L\sqrt{-g}$, $\tilde{g}_{ij} = g_{ij}$) automatically selects the magnetic mass–energy density $T_{00}^{\text{EM}}/c^2 = \rho_B = \frac{|\mathbf{B}|^2}{2\mu_0 c^2}$ as the source of the scalar field $\phi_L$!

---

## 3. Covariant Derivation of V20 Polarity-Even ($\Delta y_{\text{even}}$) Background Cross-Term Cancellation

In the laboratory, the net magnetic field is the superposition of the static background field $\mathbf{B}_{\text{bg}}(\mathbf{x})$ of the Earth (and laboratory environment) and the 8-coil field $\mathbf{B}_{\text{coil}}(\mathbf{x}; I) = I\,\mathbf{b}_{\text{coil}}(\mathbf{x})$ proportional to the current $I$ (`hu_p3.tex`, line 84):
$$\rho_B(\mathbf{x}; I) \;=\; \frac{|\mathbf{B}_{\text{bg}}(\mathbf{x}) + I\,\mathbf{b}_{\text{coil}}(\mathbf{x})|^2}{2\mu_0 c^2} \;=\; \underbrace{\frac{|\mathbf{B}_{\text{bg}}(\mathbf{x})|^2}{2\mu_0 c^2}}_{\rho_{B,\text{bg}}(\mathbf{x})\;(I^0)} \;+\; \underbrace{I\,\frac{\mathbf{B}_{\text{bg}}(\mathbf{x})\cdot\mathbf{b}_{\text{coil}}(\mathbf{x})}{\mu_0 c^2}}_{\rho_{B,\text{cross}}(\mathbf{x}; I)\;(I^1,\;\text{odd})} \;+\; \underbrace{I^2\,\frac{|\mathbf{b}_{\text{coil}}(\mathbf{x})|^2}{2\mu_0 c^2}}_{\rho_{B,\text{coil}}(\mathbf{x}; I)\;(I^2,\;\text{even})}.$$
Because both the 3D top-hat ($\mathcal{K}_L^{\text{TH}}$) and 3D Yukawa ($\mathcal{G}_L^{\text{Yuk}}$) operators are linear in the source density, the ON ($I = \pm I_0$) minus OFF ($I = 0$) difference of the scalar potential at the oscillator position $\mathbf{x}_1 = \mathbf{0}$ is:
$$\Delta\phi_L(\pm I_0) \;\equiv\; \phi_L(\mathbf{0}; \pm I_0) - \phi_L(\mathbf{0}; 0) \;=\; \pm I_0\,\gamma_L\bar{\rho}_{B,\text{cross}}(\mathbf{0}; 1) \;+\; I_0^2\,\gamma_L\bar{\rho}_{B,\text{coil}}(\mathbf{0}; 1).$$
Even if the oscillator possesses a parasitic linear Zeeman sensitivity ($\kappa_Z\,\mathbf{B}_{\text{coil}}(\mathbf{0}; \pm I_0) = \pm\kappa_Z I_0\mathbf{b}_{\text{coil}}(\mathbf{0})$), the **polarity-even average** of the H1 fractional frequency shift ($\Delta y(\pm I_0) \simeq -\Delta\phi_L(\pm I_0) \pm \kappa_Z I_0 b_{\text{coil}}$) (`hu_p3.tex`, lines 79 and 84) is:
$$\boxed{\;\Delta y_{\text{even}}(L) \;\equiv\; \frac{\Delta y(+I_0) + \Delta y(-I_0)}{2} \;=\; -I_0^2\,\gamma_L\bar{\rho}_{B,\text{coil}}(\mathbf{0}; 1) \;=\; -\frac{G L^2}{c^2}\int d^3y\,K_L(\mathbf{0},\mathbf{y})\,\frac{|\mathbf{B}_{\text{coil}}(\mathbf{y}; I_0)|^2}{2\mu_0 c^2},\;}$$
from which both the linear Zeeman term $\pm\kappa_Z I_0$ and the background cross-term $\pm 2\mathbf{B}_{\text{bg}}\cdot\mathbf{B}_{\text{coil}}$ **cancel identically**!

---

## 4. Checks

1. **Conformal vs. Disformal Maxwell-Trace Check:**
   - Conformal metric ($\tilde{g}_{\mu\nu} = f_L^2 g_{\mu\nu}$): $\sqrt{-\tilde{g}}\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta} = f_L^4 f_L^{-2}f_L^{-2}\sqrt{-g}g^{\mu\alpha}g^{\nu\beta} = \sqrt{-g}g^{\mu\alpha}g^{\nu\beta}$, independent of $\phi_L \implies \delta S_{\text{EM}}^{\text{conf}}/\delta\phi_L \equiv 0$ $\checkmark$.
   - Disformal metric ($\tilde{g}_{00} = f_L^2 g_{00}$, $\tilde{g}_{ij} = g_{ij}$): $\sqrt{-\tilde{g}}\tilde{g}^{ik}\tilde{g}^{jl} = f_L\cdot 1\cdot 1\,\sqrt{-g}g^{ik}g^{jl} \implies -\partial\mathcal{L}_{\text{EM}}^{\text{disf}}/\partial\phi_L\big|_{\mathbf{E}=\mathbf{0}} = -f_L u_B$ $\checkmark$.
2. **Polarity-Reversal ($\pm I_0$) Parity Check:**
   - $\Delta\phi_L(+I_0) + \Delta\phi_L(-I_0) = (+I_0 - I_0)\gamma_L\bar{\rho}_{B,\text{cross}} + 2 I_0^2\gamma_L\bar{\rho}_{B,\text{coil}} = 2\,\gamma_L\bar{\rho}_{B,\text{coil}}(I_0)$, which at $L = 0.630957\;\text{m}$ and $I_0 = 14.0\;\text{mA}$ reproduces the audited V20 value $\Delta y_{\text{peak}} = -1.2017108\times 10^{-48}$ $\checkmark$.
