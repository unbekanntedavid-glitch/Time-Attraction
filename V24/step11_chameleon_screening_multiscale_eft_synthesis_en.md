# Step 11: Environment-Dependent Scalar-Tensor Screening Extension Hypothesis (Branch II-Cham: $V_{\text{eff}}(\phi_L; \rho)$, Chameleon/Symmetron $\beta_{\text{eff}}(\rho)$), Correlation Length $L_{\text{eff}}(\rho)$, and Multi-Scale EFT Synthesis

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step11_chameleon_screening_multiscale_eft_synthesis_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step10_bilateral_kernel_reciprocity_eft_synthesis_en.md` (`step10_bilateral_kernel_reciprocity_eft_synthesis.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (Chapters 1 and 3, lines 171–180 and 256–265, Eq. 6), `hu_p2.tex` (Chapters 4–5, lines 24–61), and `hu_p3.tex` (Chapters 11–12, lines 99–120; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Sections 1–3; DOI: `10.5281/zenodo.23639084`)

---

## 1. Verbatim Citations and Methodological Delimitation of the Branch II-Cham Extension

In Steps 1–10, under a fixed correlation/Yukawa length $L$ and fixed coupling $\alpha_Y \in \{1, 1/3\}$ ($\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$), we derived the equations of V23 Formulation I (3D spherical top-hat, **Branch I-A/I-B**) and Formulation II (local quadratic Helmholtz–Yukawa scalar–tensor, **Branch II-Fluid / II-Conf / II-Disf**). As shown in Step 6 (`step6_sparc_galactic_ngc3198_consistency.md`) via the NGC 3198 ($R = 30.08\;\text{kpc}$) negative control check, the isolated V23 model with constant coupling $\alpha_{\text{std}} = 1/(4\pi) \approx 0.07958$ yields at most a squared-velocity boost of $C_{\text{pred}}^{\text{isol}} \le 1 + \alpha_{\text{std}} = 1.0796$, compared to the observed $C_{\text{req}} = 2.720\text{--}4.311$.

**Methodological Delimitation (Point 5):** The density-dependent $\beta_{\text{eff}}(\rho)$ and $V(\phi_L)$ chameleon/symmetron mechanism examined in this step (**Branch II-Cham** in the taxonomy of Step 2) is a **standalone nonlinear extension hypothesis** that goes beyond the quadratic, constant-$\alpha_Y$ Yukawa action of V23 (`50aa2db0`), providing a candidate 4D effective field theory (EFT) framework for open task **(d)** of Chapter 11 of the V22 monograph (`hu_p3.tex`, line 104: *"together with a physical scale-selection mechanism for $L$"*).

1. **The V5–V16 Chameleon and Horndeski Precedents, and the Force Multiplier $C(\rho)$ (`hu_p1.tex`, lines 172–174 and 256–265, Eq. 6):**
   > *"V5 – V7: Weak-Field and Chameleon Explorations — I investigated scalar–tensor and chameleon-like screening analogies to connect local density gradients with low-density spatial effective acceleration."*
   > *"In Versions 14–16 (DOI `10.5281/zenodo.23217189`), during my first phenomenological studies of galactic rotation curves, I parameterized the effective gravitational boost using the following density-dependent force multiplier:"*
   $$\begin{equation}
     C(\rho) \;=\; 1 + C_0\,\exp\!\left(-\frac{\rho}{\rho_c}\right), \qquad V_{\text{rot}}^2(r) \;=\; C(\rho(r))\,V_{\text{bar}}^2(r), \tag{V22 Eq. 6}
   \end{equation}$$
   > *"where $\rho_c$ is the critical transition density, and in the early works I set $C_0 = 5$ (so that for $\rho \ll \rho_c$, $C(\rho) \to 6$, whereas at high density $\rho \gg \rho_c \implies C(\rho) \to 1$, recovering the Newtonian limit in the Solar System)."*

2. **Open Theoretical Task (d) of Chapter 11 of the V22 Monograph (`hu_p3.tex`, line 104):**
   > *"**(d) Self-Consistent Multi-Channel Observational Evaluation:** Once the action principle and matter coupling in (a)–(c) are fixed, galactic rotation curves (including cosmic-background embedding $\rho_{\text{bg}}$), VHE gamma-ray timing across cluster overdensities, Sirius B, and NICER posteriors must be jointly fitted **together with a physical scale-selection mechanism for $L$**."*

---

## 2. Mathematical Derivation: The Environment-Dependent Effective Potential $V_{\text{eff}}(\phi_L; \rho)$ and Scale Selection $L_{\text{eff}}(\rho)$

### 2.1. General Self-Interacting Scalar-Tensor Action and the Effective Potential
Let us generalize the quadratic mass term ($\frac{1}{2L_0^2}\phi_L^2$) of V23 Formulation II (Eq. 8) to a general self-interaction potential $V(\phi_L)$ in the 4-dimensional scalar–tensor action:
$$S \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}R \;-\; \frac{\kappa_\phi}{2}g^{\mu\nu}\nabla_\mu\phi_L\nabla_\nu\phi_L \;-\; \kappa_\phi\,V(\phi_L) \;+\; f_L(\phi_L)\,\mathcal{L}_m\right], \qquad \kappa_\phi \;=\; \frac{c^4}{4\pi G\alpha_{\text{std}}} \;=\; \frac{c^4}{G\alpha_Y}, \quad f_L(\phi_L) = e^{-\phi_L}.$$
*(Note on normalization: $\kappa_\phi/2 = \frac{c^4}{2 G\alpha_Y} = \frac{c^4}{16\pi G}\frac{1}{\tilde{\alpha}_L}$, agreeing with the V23 Eq. 8 normalization used in Steps 3 and 7!)* In a non-relativistic baryonic medium ($\mathcal{L}_m = -\rho c^2$), the Euler–Lagrange equation for the field $\phi_L$ is:
$$\Box_g\phi_L \;=\; \frac{\partial V_{\text{eff}}(\phi_L; \rho)}{\partial\phi_L}, \qquad \boxed{\;V_{\text{eff}}(\phi_L; \rho) \;\equiv\; V(\phi_L) \;+\; \frac{\rho c^2}{\kappa_\phi}\,e^{-\phi_L} \;=\; V(\phi_L) \;+\; \frac{\alpha_Y G}{c^2}\,\rho\,e^{-\phi_L}\;}$$
(where $\alpha_Y = 4\pi\alpha_{\text{std}} = 8\pi\tilde{\alpha}_L$).

### 2.2. Equilibrium Minimum $\phi_{\min}(\rho)$ and the Density-Dependent Yukawa Length $L_{\text{eff}}(\rho)$
In a locally quasi-homogeneous region of density $\rho$, the scalar field takes the minimum $\frac{\partial V_{\text{eff}}}{\partial\phi_L}\big|_{\phi_{\min}} = 0$, and small perturbations $\delta\phi_L = \phi_L - \phi_{\min}(\rho)$ propagate with an **environment-dependent Yukawa length** ($L_{\text{eff}}(\rho)$) determined by the local curvature:
$$\frac{\partial V_{\text{eff}}}{\partial\phi_L}\bigg|_{\phi_{\min}} \;=\; V'\bigl(\phi_{\min}(\rho)\bigr) \;-\; \frac{\alpha_Y G}{c^2}\,\rho\,e^{-\phi_{\min}(\rho)} \;=\; 0,$$
$$\boxed{\;L_{\text{eff}}^{-2}(\rho) \;\equiv\; \mu_{\text{eff}}^2(\rho) \;\equiv\; \frac{\partial^2 V_{\text{eff}}(\phi_L; \rho)}{\partial\phi_L^2}\bigg|_{\phi_L = \phi_{\min}(\rho)} \;=\; V''\bigl(\phi_{\min}(\rho)\bigr) \;+\; \frac{\alpha_Y G}{c^2}\,\rho\,e^{-\phi_{\min}(\rho)}.\;}$$
Let us examine the two fundamental physical cases of this formula:

1. **Quadratic (Constant Bare Mass) Limit of V23 (**Branch II-Fluid**, $V(\phi_L) = \frac{\phi_L^2}{2L_0^2}$):**
   Here $V'(\phi_L) = \phi_L / L_0^2$ and $V''(\phi_L) = L_0^{-2}$. The minimum condition can be written in terms of the Lambert $W$ function ($\phi_{\min}e^{\phi_{\min}} = \alpha_Y\gamma_{L_0}\rho$, so $\phi_{\min}(\rho) \simeq \alpha_Y\gamma_{L_0}\rho$ for small $\phi_L$), and the local effective Yukawa length **decreases with increasing density** even from the bare $e^{-\phi_L}\rho c^2$ coupling alone:
   $$L_{\text{eff}}^{-2}(\rho) \;=\; L_0^{-2}\Bigl(1 + \phi_{\min}(\rho)\Bigr) \;\simeq\; L_0^{-2}\Bigl(1 + \alpha_Y\,\gamma_{L_0}\,\rho\Bigr) \quad\Longrightarrow\quad L_{\text{eff}}(\rho) \;=\; \frac{L_0}{\sqrt{1 + \phi_{\min}(\rho)}}.$$

2. **Nonlinear Chameleon / Symmetron Extension Hypothesis (**Branch II-Cham**, `hu_p1.tex` V5–V7):**
   If, going beyond the V23 quadratic model, the bare potential is $V(\phi_L) = \Lambda_0^4\,\exp(\sigma\phi_L)$ or $V(\phi_L) = M_0^4\phi_L^{-n}$ ($n > 0$) and the coupling function has a general nonlinear form $f_{\text{eff}}(\phi_L) = \exp[-A(\phi_L)]$ whose effective slope at the local minimum is $\beta_{\text{eff}}(\rho) = \beta_0\exp\!\left(-\frac{\rho}{2\rho_c}\right)$, then the effective gravitational multiplier ($G_{\text{eff}}/G = 1 + 2\beta_{\text{eff}}^2(\rho)$) reproduces the V14–V16 phenomenological force multiplier $C(\rho)$ (`hu_p1.tex`, Eq. 6):
   $$\boxed{\;\frac{G_{\text{eff}}(\rho)}{G} \;=\; 1 + 2\beta_{\text{eff}}^2(\rho) \;=\; 1 + 2\beta_0^2\exp\!\left(-\frac{\rho}{\rho_c}\right) \;\equiv\; C(\rho), \qquad C_0 \;\equiv\; 2\beta_0^2 \;=\; \alpha_{\text{std},0}.\;}$$
   - **Galactic outer regime ($\rho \ll \rho_c \sim 10^{-23}\text{--}10^{-21}\;\text{kg/m}^3$):** in this nonlinear extension hypothesis, screening is lifted and $2\beta_{\text{eff}}^2(\rho) \to C_0$ (at the NGC 3198 $R = 30.08\text{ kpc}$ point, $C_{\text{req}} - 1 = 1.720$ [$\Upsilon=1.0$] and $3.311$ [$\Upsilon=0.5$]).
   - **High-density Solar System and terrestrial regime ($\rho \gg \rho_c$):** the exponential factor suppresses the long-range scalar force multiplier ($C(\rho) - 1 = C_0 e^{-\rho/\rho_c} \to 0$), leaving only the short-scale ($L$) canonical $\gamma_L = G L^2/c^2$ ($\alpha_{\text{std}} = 1/(4\pi)$) V23 mode.

---

## 3. Multi-Scale EFT Overview Across the V17–V23 Investigated Regimes (with Branch Labels)

The table below specifies for every physical system examined which model branch (**Branch I-A**, **Branch II-Fluid**, **Branch II-Disf**, **Branch II-Cham**) governs the corresponding result:

| Physical Scale and System | Model Branch (`step2` Taxonomy) | Characteristic Density $\rho$ [$\text{kg/m}^3$] | Investigated $L$ Scale | $\phi_L = \gamma_L\bar{\rho}_L$ (or $\phi_L^{\text{Yuk}}$) | Primary Observable and V17–V23 Audit Outcome |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **V20 8-coil magnetic field** (`step8`, `step20`) | **Branch II-Disf / I-A** ($H_1$) | $\rho_B = \frac{B^2}{2\mu_0 c^2} \sim 10^{-22}$ | $L_{\text{peak}} = 0.630957\text{ m}$ | $1.2017108\times 10^{-48}$ | Polarity-even $\Delta y_{\text{peak}} = -1.2017108\times 10^{-48}$; Zeeman + cross-term canceled. |
| **V20 $8558\text{ kg}$ granite chamber** (`step8`, `step19`) | **Branch I-A / II-Fluid** ($H_0$ vs. $H_1$) | $2.75\times 10^3$ | $L = 0.50\text{--}1.00\text{ m}$ | $1.436\times 10^{-24}$ ($1\text{ m}$) | $H_0$: $+1.331\times 10^{-47}$; conditional $H_1$: $-1.436\times 10^{-24}$ (top-hat), $-3.781\times 10^{-25}$ (Yukawa). |
| **Sirius B white dwarf** (`step8`, `step19`) | **Branch I-A / II-Fluid** ($H_0$ vs. $H_1$) | $\bar{\rho} \sim 3\times 10^9$ | $L = 1\text{--}5000\text{ km}$ | $1.6\times 10^{-5}$ ($5000\text{ km}$) | $H_0$: $|\Delta v|\le 1.27\text{ m/s}$ (compatible); $H_1$: excluded for $L \ge 1\text{ km}$ ($L \lesssim 10^6\text{ m}$ for local kernel). |
| **Neutron stars** (`step7`, `step18`, `step22`) | **Branch I-A** vs. **Branch II-Fluid** | $10^{17}\text{--}1.45\times 10^{18}$ | $L = 1\text{--}10\text{ km}$ | $0.016\text{--}0.050$ ($10\text{ km}$) | Curved exterior `yukawa_action_el` $M_{\text{ADM}}$ and coupled $(\delta h,\delta\phi_L)$ stability (`step22`); V21 nonlocal tangent BVP (`step18`). |
| **Solar System (PPN / Cassini)** (`step5`) | **Branch II-Conf** vs. **II-Disf / I-A** | $\bar{\rho}_L^\odot \propto L^{-3}$ | $L = 1\text{ AU}\text{--}1\text{ pc}$ | $1.14\times 10^{-14}$ ($1\text{ pc}$) | 3D top-hat dilution: $|1-f_L^\odot| \approx 1.14\times 10^{-14} \ll 10^{-5}$; for $r \gg L$, $\gamma_{\text{PPN}} \to 1$. |
| **SPARC (175/155) and NGC 3198** (`step6`) | **Branch I-A / II-Fluid** (negative control) vs. **II-Cham** | $10^{-24}\text{--}10^{-20}$ | $L \sim 1\text{--}23.09\text{ kpc}$ | $\le 2.300\times 10^{-7}$ (isolated) | 0 points with $\rho \le 0$ (kernel positivity); NGC 3198 ($30.08\text{ kpc}$): $C_{\text{req}} = 2.720 / 4.311$, whereas constant $\alpha_{\text{std}}$: $C_{\text{pred}}^{\text{isol}} \le 1.0796$. |
| **V17 H.E.S.S./MAGIC gamma timing** (`step9`) | **Branch II-Disf** ($H_1$ EM) | $\rho_{\text{IGM/ISM}} \sim 10^{-27}\text{--}10^{-21}$ | $z = 0.034\text{--}0.4245$ | $\int_0^D \phi_L(s)\,ds$ | Naive $z$-proportional LIV excluded; prediction correlates with cluster overdensity ($\int\bar{\rho}_L\,ds$). |

---

## 4. Checks

1. **Limiting Values of the Correspondence $C(\rho) = 1 + 2\beta_{\text{eff}}^2(\rho)$ (Branch II-Cham):**
   - $\rho \to 0$ (deep galactic outer regime): $2\beta_{\text{eff}}^2(0) = 2\beta_0^2 = C_0 \implies C(0) = 1 + C_0$ $\checkmark$.
   - $\rho \gg \rho_c$ (stellar interior, Solar System, laboratory): $2\beta_{\text{eff}}^2(\rho) = C_0 e^{-\rho/\rho_c} \to 0 \implies C(\rho) \to 1$ $\checkmark$.
2. **Verification of the Quadratic Yukawa Limit $L_{\text{eff}}(\rho) = L_0 / \sqrt{1 + \phi_{\min}(\rho)}$ (Branch II-Fluid):**
   - At low density ($\phi_{\min}(\rho) \to 0$), $L_{\text{eff}}(\rho) \to L_0$, recovering the constant Yukawa length $L_0$ of V23 Formulation II $\checkmark$.
