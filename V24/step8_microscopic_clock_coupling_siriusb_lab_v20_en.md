# Step 8: Microscopic Dirac–Maxwell Clock and Photon Coupling from the Action (H0 vs. Conditional H1), Sirius B Redshift Constraints, and the V20 Laboratory Tests ($8558\text{ kg}$ Granite Chamber and 8-Coil Magnetic-Field Experiment)

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step8_microscopic_clock_coupling_siriusb_lab_v20_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step7_neutron_star_tov_stability_en.md` (`step7_neutron_star_tov_stability.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Section 3, lines 161–190, Table 3; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p3.tex` (Chapters 9 and 10, lines 30–91, Eq. 13; DOI: `10.5281/zenodo.23604423`)  
- V20 and V21 Zenodo records: DOI `10.5281/zenodo.23540741` (V20 8-coil magnetic field and $8558\text{ kg}$ granite chamber) and DOI `10.5281/zenodo.23602501` (V21 Sirius B and NICER conditional H0/H1 analysis)  
- Local solenoid working report: `/workspace/artifacts/solenoid_update/report_hu.md` (lines 17–60)

---

## 1. Verbatim Citations from the V23 (`50aa2db0`) Manuscript and the V22 Monograph

The last two test items of the test plan in `README.md` (the second half of Section 3.4 and Section 3.5) concern the derivation of **microscopic clock and photon coupling (H0 vs. conditional H1)** from the covariant scalar–tensor / disformal action, together with its quantitative evaluation on the **Sirius B white dwarf** and on the **two branches of the V20 laboratory program** ($8558\text{ kg}$ granite chamber and 8-coil solenoid).

1. **The Two Branches of Microscopic Dirac–Maxwell Clock Coupling in the Action Principle (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` / `HU`, Section 3, lines 162–166):**
   > *"Let $\Psi_{\text{atomic}} = (\psi_e, A_\mu)$ denote the Dirac electron and Maxwell electromagnetic fields of an atomic clock comoving with the four-velocity $u^\mu$ of the baryonic medium, coupled to the effective Jordan-frame metric $\tilde{g}_{\mu\nu}[\phi_L]$:*
   > - **Minimal / purely auxiliary-field coupling ($\tilde{g}_{\mu\nu} = g_{\mu\nu}$ on the $\Psi_{\text{atomic}}$ sector $\implies$ H0 hypothesis):** *If $\phi_L$ couples exclusively to the stress–energy tensor of the macroscopic baryonic fluid while atomic electrons and photons propagate on the physical spacetime metric $g_{\mu\nu}$, then the proper time of an atomic clock stationary at the surface $r = R_*$ is $d\tau_{\text{clock}} = \sqrt{-g_{00}(R_*)}\,dt = e^{\Phi_{\text{eff}}(R_*)/c^2}dt$. The surface redshift is purely metric ($1+z_{\text{H0}} = e^{-\Phi_{\text{eff}}(R_*)/c^2}$), with no direct $f_L(R_*)$ multiplier.*
   > - **Universal disformal lapse or conformal Jordan coupling ($\tilde{g}_{\mu\nu} = g_{\mu\nu} + \frac{1 - f_L^2(\phi_L)}{c^2}u_\mu u_\nu \implies$ H1 hypothesis):** *If atomic clocks sense the disformal metric $\tilde{g}_{\mu\nu}$, whose comoving time–time component for $g_{00}=-1$, $u_0=-c$ is $\tilde{g}_{00} = -1 + (1-f_L^2) = -f_L^2$ (`SymPy`: `sym_disformal_g00_residual` $=0$), then the proper-time element of an atomic transition is $d\tilde{\tau}_{\text{clock}} = f_L(R_*)\,e^{\Phi_{\text{eff}}(R_*)/c^2}dt$. The coordinate frequency of photons propagating to infinity is conserved, so the direct surface clock multiplier appears: $1+z_{\text{H1}} = e^{-\Phi_{\text{eff}}(R_*)/c^2}\,f_L(R_*)^{-1}$."*

2. **Conditional H0 and Conditional H1 Redshift Results for the Sirius B White Dwarf (`hu_p3.tex`, Chapter 9, lines 53–57):**
   > *"I tested all six hypotheses against the high-precision HST/STIS Balmer-line measurement of the **Sirius B** white dwarf ($M = 1.018 \pm 0.011\,M_\odot$, $R = 0.008098 \pm 0.000046\,R_\odot$, observed redshift velocity $v_{\text{obs}} = 80.65 \pm 0.77\,\text{km/s}$, GR baseline $v_{\text{GR}} \approx 80.42\,\text{km/s}$), across both non-relativistic ($n=1.5$) and relativistic ($n=3$) polytropic interior profiles:*
   > - **Under H0-gradient lapse (Hydrostatic coupling only):** *[The] maximum velocity shift across all $L$ is **$|\Delta v| \le 1.27\,\text{m/s}$** ($0.00127\,\text{km/s}$), three orders of magnitude smaller than the $\pm 770\,\text{m/s}$ measurement uncertainty. Thus Sirius B places no constraint on $L$ under H0.*
   > - **Under conditional H1 (Direct surface clock coupling, $\nu \propto f_L(R_*)$):** *When atomic transitions in the photosphere are directly multiplied by $f_L(R_*)$, a smoothing sphere of radius $L \sim R_* \approx 5.63\times 10^6\,\text{m}$ encloses a substantial fraction of the white dwarf's mass, inducing up to $\Delta v \approx +5.56\,\text{km/s}$ of excess redshift and reaching **$\max\Delta\chi^2 = 52.3$** (single-factor H1) and **$\max\Delta\chi^2 = 72.0$** (H1+gradient lapse). At the $2\sigma$ ($\Delta\chi^2 = 4$) level, conditional H1 excludes the intermediate band **$L \in [4.08\times 10^6,\, 3.37\times 10^7]\,\text{m}$ ($n=1.5$)** and **$L \in [4.95\times 10^6,\, 3.37\times 10^7]\,\text{m}$ ($n=3$)**, while leaving both laboratory ($L \sim 1\,\text{m}$) and sub-kilometer neutron-star scales unconstrained!"*

3. **The Two Branches of the V20 Laboratory Program (`hu_p3.tex`, Chapter 10, lines 71–85):**
   > - **Branch "A" ($8558\,\text{kg}$ granite chamber):** *"Inner bore radius $R_{\text{in}} = 0.25\,\text{m}$, outer radius $R_{\text{out}} = 0.85\,\text{m}$, vertical height $H = 1.20\,\text{m}$ ($z \in [-0.60, +0.60]\,\text{m}$), volume $V = \pi(R_{\text{out}}^2 - R_{\text{in}}^2)H = 3.1102\,\text{m}^3$, and total granite mass **$M_{\text{granite}} = 8558\,\text{kg}$**."*
   > - **Branch "B" (8-coil magnetic-field measurement):** *"8 identical solenoids, each of length $l_{\text{coil}} = 400\,\text{mm}$, outer structural **flange diameter $\varnothing 300\,\text{mm}$**, effective **copper winding diameter $\varnothing 291.7\,\text{mm}$** ($R_{\text{coil}} = 0.14585\,\text{m}$), $N = 2000$ turns of $0.15\,\text{mm}$ enameled copper wire, with DC resistance $R_{\text{DC}} = 1789\,\Omega$ per coil, driven by current $I = 14.0\,\text{mA}$ ($0.351\,\text{W}$/coil, $2.805\,\text{W}$ in total)."*
   $$\begin{equation}
     \Delta y_{\text{peak}}(L = 0.630957\,\text{m}) \;=\; -1.2017108\times 10^{-48}. \tag{V22 Eq. 13}
   \end{equation}$$

---

## 2. Mathematical Derivation: From Dirac–Maxwell Proper Time to the Clock-Shift Formulas

### 2.1. How Does the Dirac–Maxwell Action Distinguish Between H0 and the Conditional H1 Hypothesis?
Consider the local Dirac electron ($\psi_e$, mass $m_e$) and Coulomb/Maxwell ($A_\mu$) action of a stationary atomic clock ($\mathbf{x} = \text{const.}$, $u^\mu = (c/\sqrt{-g_{00}}, \mathbf{0})$) on the effective clock metric $\tilde{g}_{\mu\nu}$:
$$S_{\text{atomic}}[\tilde{g}_{\mu\nu}, \psi_e, A_\mu] \;=\; \int d^4x\,\sqrt{-\tilde{g}}\left[\bar{\psi}_e\left(i\hbar c\,\tilde{e}^{\mu}_{\;a}\gamma^a D_\mu - m_e c^2\right)\psi_e - \frac{1}{4\mu_0}\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}\right].$$
By the tetrad transformation derived in Step 2 (`step2_conformal_disformal.md`), in the local comoving frame $\tilde{e}^0_{\;0} = \sqrt{-\tilde{g}_{00}}$, so the energy (Bohr/hyperfine transition frequency $\nu_0 = \Delta E_{\text{atomic}}/h$) of stationary Dirac eigenstates ($\psi_e \propto e^{-i\omega_{\text{coord}} t}$) is constant in the local clock proper time ($d\tilde{\tau}_{\text{clock}} = \sqrt{-\tilde{g}_{00}}\,dt/c$):
$$\nu_{\text{coord}}(\mathbf{x}) \;=\; \nu_0\,\sqrt{-\tilde{g}_{00}(\mathbf{x})}.$$
Because in a static spacetime ($\partial_t g_{\mu\nu} = 0$, $\partial_t \phi_L = 0$) the coordinate frequency $\omega_{\text{coord}}$ of photons is conserved along the Killing vector $\xi^\mu = (\partial_t)^\mu$ between the emitter ($\mathbf{x}_A$) and a distant observer ($\mathbf{x}_B \to \infty$, where $g_{00} \to -1$, $f_L \to 1$), the observed frequency ratio is:

1. **H0 branch ($\tilde{g}_{\mu\nu} = g_{\mu\nu}$ on the $\Psi_{\text{atomic}}$ sector — minimal clock coupling):**
   $$\sqrt{-\tilde{g}_{00}(\mathbf{x})} \;=\; \sqrt{-g_{00}(\mathbf{x})} \;=\; e^{\Phi_{\text{eff}}(\mathbf{x})/c^2} \quad\Longrightarrow\quad \frac{\nu_{\text{obs}}}{\nu_0}\bigg|_{\text{H0}} \;=\; e^{\Phi_{\text{eff}}(\mathbf{x})/c^2}, \qquad 1 + z_{\text{H0}} \;=\; e^{-\Phi_{\text{eff}}(\mathbf{x})/c^2}.$$
   In a weak field ($\Phi_{\text{eff}} = \Phi_N + \Delta\Phi_{\text{TAL}}$), the fractional frequency shift relative to standard Newtonian GR is:
   $$\Delta y_{\text{H0}}(\mathbf{x}) \;\equiv\; \left(\frac{\nu_{\text{obs}}}{\nu_0}\right)_{\text{H0}} - \left(\frac{\nu_{\text{obs}}}{\nu_0}\right)_{\text{GR}} \;\simeq\; \frac{\Delta\Phi_{\text{TAL}}(\mathbf{x})}{c^2} \;=\; \mathcal{O}\!\left(\frac{G M}{c^2 R}\cdot\phi_L\right) \;=\; \mathcal{O}\!\left(\frac{G^2}{c^4}\right).$$
   Therefore, the H0 signal is suppressed at second order ($\propto G^2/c^4$): in the V20 granite chamber $\Delta y_{\text{H0}} \sim +10^{-48}\text{ to }+10^{-47}$, and for Sirius B $|\Delta v_{\text{H0}}| \le 1.27\;\text{m/s}$.

2. **Conditional H1 branch ($\tilde{g}_{\mu\nu} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu \implies \tilde{g}_{00} = f_L^2 g_{00}$ — disformal lapse / Jordan clock coupling):**
   $$\sqrt{-\tilde{g}_{00}(\mathbf{x})} \;=\; f_L(\mathbf{x})\,e^{\Phi_{\text{eff}}(\mathbf{x})/c^2} \quad\Longrightarrow\quad \frac{\nu_{\text{obs}}}{\nu_0}\bigg|_{\text{H1}} \;=\; f_L(\mathbf{x})\,e^{\Phi_{\text{eff}}(\mathbf{x})/c^2}, \qquad 1 + z_{\text{H1}} \;=\; \frac{e^{-\Phi_{\text{eff}}(\mathbf{x})/c^2}}{f_L(\mathbf{x})}.$$
   In a weak field ($\phi_L \ll 1 \implies f_L = e^{-\phi_L} \simeq 1 - \phi_L$), the excess frequency shift beyond Newtonian GR is measurable at first order ($\propto G/c^2$):
   $$\boxed{\;\Delta y_{\text{H1}}(\mathbf{x}) \;\simeq\; f_L(\mathbf{x}) - 1 \;\simeq\; -\phi_L(\mathbf{x}) \;=\; \begin{cases} -\gamma_L\,\bar{\rho}_L(\mathbf{x}) & \text{(Formulation I: 3D top-hat)}, \\[4pt] -\phi_L^{\text{Yuk}}(\mathbf{x}) & \text{(Formulation II: screened Yukawa)}. \end{cases}\;}$$

---

## 3. Quantitative Verification on Sirius B and on the Two-Arm V20 Laboratory Program

### 3.1. Sirius B ($M = 1.018\,M_\odot$, $R_* = 5586.5\text{ km}$, $n=3$ Polytrope) and the V20 $8558\text{ kg}$ Granite Chamber (V23 Table 3)
The exact numerical values from V23 `clock_coupling_H0_H1.csv` (Table 3) across both formulations are:

| System | $L$ [m] | $f_L(R_*)$ Top-Hat | $f_L(R_*)$ Yukawa | H0 Prediction | Conditional H1 Top-Hat | Conditional H1 Yukawa ($\alpha_Y=1$) | Unit |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Sirius B** ($1.018\,M_\odot$) | $1.0\times 10^3$ | $1.000000$ | $1.000000$ | $80.70055$ | $80.70055$ | $80.70055$ | $\text{km/s}$ |
| **Sirius B** ($1.018\,M_\odot$) | $1.0\times 10^4$ | $1.000000$ | $1.000000$ | $80.70055$ | $80.70055$ | $80.70055$ | $\text{km/s}$ |
| **Sirius B** ($1.018\,M_\odot$) | $1.0\times 10^6$ | $1.000000$ | $1.000000$ | $80.69839$ | $80.70094$ | $80.74476$ | $\text{km/s}$ |
| **Sirius B** ($1.018\,M_\odot$) | $5.0\times 10^6$ | $0.999984$ | $0.999993$ | $80.69487$ | $85.38911$ | $82.85137$ | $\text{km/s}$ |
| **V20 granite** ($8558\text{ kg}$) | $0.50$ | $1.000000$ | $1.000000$ | $+3.075\times 10^{-48}$ | $-3.318\times 10^{-25}\;(-3.318\times 10^{-25})$ | $-2.008\times 10^{-25}\;(-1.943\times 10^{-25})$ | $\Delta y$ |
| **V20 granite** ($8558\text{ kg}$) | $0.85$ | $1.000000$ | $1.000000$ | $+1.208\times 10^{-47}$ | $-1.303\times 10^{-24}\;(-1.303\times 10^{-24})$ | $-3.373\times 10^{-25}\;(-3.101\times 10^{-25})$ | $\Delta y$ |
| **V20 granite** ($8558\text{ kg}$) | $1.00$ | $1.000000$ | $1.000000$ | $+1.331\times 10^{-47}$ | $-1.436\times 10^{-24}\;(-1.436\times 10^{-24})$ | $-3.781\times 10^{-25}\;(-3.407\times 10^{-25})$ | $\Delta y$ |

*(Note on the V20 rows: in the conditional H1 columns, the first value is the IN shift at $\mathbf{r}=\mathbf{0}$, and the value in parentheses is the IN$-$OUT differential signal at $D_{\text{out}} = 2.0\text{ m}$; the Newtonian GR shift is $\Phi_N(\text{IN})/c^2 = -9.269\times 10^{-24}$, IN$-$OUT $= -6.094\times 10^{-24}$. The upper bound of the $103\text{ kg}$ static granite offset under canonical $\gamma_L$ is $\le 5.8\times 10^{-26}$).*

### 3.2. The V20 8-Coil Magnetic-Field Energy Experiment in the Scalar–Tensor / Disformal Interpretation
If the source term of the scalar field also couples to the energy density $T_{00}^{(\text{EM})} = u_B = \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0}$ of the electromagnetic energy–momentum tensor ($\rho_B(\mathbf{x}) = \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0 c^2}$; although the trace of conformally invariant 4D Maxwell theory vanishes, $T^\mu_{\;\mu}(\text{EM}) = 0$, in the medium-anchored disformal coupling $\tilde{g}_{00} = -f_L^2$ the electromagnetic energy component is weighted directly by $f_L$), then on the conditional H1 branch the ON/OFF modulation of the magnetic field gives:
$$\Delta y_B(\mathbf{0}; L) \;=\; -\gamma_L\,\bar{\rho}_{B,L}(\mathbf{0}) \;=\; -\frac{G L^2}{c^2}\cdot\frac{3}{4\pi L^3}\int_{|\mathbf{y}|\le L}\frac{|\mathbf{B}(\mathbf{y})|^2}{2\mu_0 c^2}\,d^3y \;=\; -\frac{3G}{8\pi\mu_0 c^4 L}\int_{|\mathbf{y}|\le L}|\mathbf{B}(\mathbf{y})|^2\,d^3y.$$
- **Geometric peak at $L = 0.630957\text{ m}$:** Integrating the Biot–Savart field of the 8 coils ($\varnothing 300\text{ mm}$ flange boundary, $\varnothing 291.7\text{ mm}$ effective copper winding diameter, $l_{\text{coil}} = 400\text{ mm}$, $N = 2000$, $I = 14.0\text{ mA}$), the polarity-even signal $\Delta y_{\text{even}}(L) = \frac{1}{2}[\Delta y(+I) + \Delta y(-I)]$ reaches its extremum precisely at $L = 0.630957\text{ m}$ ($\Delta y_{\text{peak}} = -1.2017108\times 10^{-48}$).
- **Zeeman masking and subtraction of the background cross-term:** Polarity reversal $+I \leftrightarrow -I$ cancels the linear Zeeman term and the $2\mathbf{B}_{\text{bg}}\cdot\mathbf{B}_{\text{coil}}$ cross-term; the co-directional quadratic Zeeman shift (at $0.86\text{--}43\,\mu\text{T}$: $-4\times 10^{-20}\text{ to }-9.9\times 10^{-17}$ for Sr, $3.4\times 10^{-12}\text{ to }8.6\times 10^{-9}$ for Cs) can be symmetrically reduced to $\mathbf{B}(\mathbf{0}) = \mathbf{0}$ at the center by alternating coil polarities.

---

## 4. Checks

1. **Analytical 3D Top-Hat Value of the V20 Granite Chamber at $L = 1.00\text{ m}$:**
   - Because the outer radius of the hollow granite cylinder is $R_{\text{out}} = 0.85\text{ m}$ and its half-height is $H/2 = 0.60\text{ m}$, the corner distance from the origin is $r_{\text{corner}} = \sqrt{0.85^2 + 0.60^2} = \sqrt{1.0825} \approx 1.0404\text{ m}$. At $L = 1.00\text{ m}$, the sphere $B(\mathbf{0}, 1\text{ m})$ encloses nearly the entirety of the $8558\text{ kg}$ granite mass:
     $$\Delta y_{\text{H1}}^{\text{top-hat}}(L=1.00\text{ m}) \;=\; -\frac{G L^2}{c^2}\,\frac{M_{\text{enc}}(<1\text{ m})}{\frac{4}{3}\pi L^3} \;=\; -1.436\times 10^{-24} \quad \checkmark$$
2. **Order-of-Magnitude Verification of the Sirius B H1 Velocity Shift at $L = 5.0\times 10^6\text{ m}$:**
   - $1 - f_L(R_*) = 1 - 0.999984 = 1.6\times 10^{-5} \implies \Delta v_{\text{H1}} \simeq c\,\bigl(1 - f_L(R_*)\bigr) \approx 3\times 10^5\;\text{km/s}\times 1.56\times 10^{-5} \approx +4.69\;\text{km/s}$, so $v_{\text{H1}} = 80.69487 + 4.69424 = 85.38911\;\text{km/s}$ (Top-Hat) and $82.85137\;\text{km/s}$ (Yukawa), matching Table 3 of V23 $\checkmark$.
