# Step 20: Exact Derivation of the V20 Eight-Coil Elliptic Biot–Savart Loop Integral, the $r = L/t$ Compactified Exterior Quadrature ($3.7\times 10^{-16}$ / $7.4\times 10^{-16}$ Agreement), and the $L = 0.95\;\text{m}$ Switching Identity

**Author:** Dávid Szabolcsi (ORCID: `0009-0005-6300-074X`)  
**Date:** 2026-10-26  
**File Path:** `scalar-tensor-yukawa-fit/step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature_en.md`  
**Previous Steps:** `step1_dilaton_form_en.md` (`step1_dilaton_form.md`) – `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy_en.md` (`step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy.md`)  
**License:** CC BY-NC 4.0  
**Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (Chapter 10, point 2, lines 470–485 and Appendix, lines 557–618; DOI: `10.5281/zenodo.23604423`)  
- V20 (`10.5281/zenodo.23540741`) verification script and dataset: `/workspace/ta-v18/V20/make_revision3_checks.py` (lines 1–77) and `/workspace/ta-v18/V20/revision3_checks.json` (lines 123–197)  
- V23 (`50aa2db0` task): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (line 28; DOI: `10.5281/zenodo.23639084`)

---

## 1. Verbatim Citations from the V22 and V20 Manuscripts

In Step 15 (`step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md`), we derived why the trace-free Maxwell tensor generates a source $\rho_B(\mathbf{x}) = \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0 c^2}$ in the V23 disformal metric, and in Step 19 we closed the 3D IN–OUT quadrature of V20 Branch "A" ($8558\;\text{kg}$ granite chamber). Point 2 of Chapter 10 (lines 470–485) and the Appendix (lines 557–618) of the V22 monograph, together with V20 `make_revision3_checks.py` (lines 25–75), contain the exact elliptic-integral magnetic-field evaluation of V20 Branch "B" (the eight-coil magnetic-energy experiment), the two exterior quadrature paths (Path A and Path B), and the $L = 0.95\;\text{m}$ interior/exterior switching identity:

1. **Parameters of the 8-Coil Arrangement and the $L = 0.630957\;\text{m}$ Grid Maximum (`tal_book_HU_2026-10-24.tex`, Chapter 10, point 2, lines 473–478):**
   > *"**Audited 8-Coil Hardware Parameters:** 8 identical solenoids, each of length $l_{\text{coil}} = 400\,\text{mm}$, outer structural **$\varnothing 300\,\text{mm}$ flange diameter**, effective **$\varnothing 291.7\,\text{mm}$ copper winding diameter** ($R_{\text{coil}} = 0.14585\,\text{m}$), $N = 2000$ turns of $0.15\,\text{mm}$ enameled copper wire, DC resistance $R_{\text{DC}} = 1789\,\Omega$ per coil, driven by a current $I = 14.0\,\text{mA}$ ($0.351\,\text{W}$/coil, $2.805\,\text{W}$ in all)."*
   > *"**Maximum on the Sampled 25-Point $L$-Grid ($L = 0.630957\,\text{m}$):** Evaluating the magnetic energy density of the 8 coils on the 25-point logarithmic grid $L \in [0.1, 3.0]\,\text{m}$, the conditional polarity-even shift $\Delta y(L) = \frac{1}{2}[\Delta y(+I) + \Delta y(-I)]$ attains its discrete maximum on the grid at **$L = 0.630957\,\text{m}$** (not a continuous extremum fit):*
   $$\begin{equation}
     \Delta y_{\text{grid max}}(L = 0.630957\,\text{m}) \;=\; -1.2017108\times 10^{-48}. \tag{V22 Eq. 14}
   \end{equation}$$
   > *Here I separate: (i) the $3.7\times 10^{-16}$ agreement of the two deterministic integration routines at $L=3.0\,\text{m}$, (ii) the $-0.0722\%$ Sobol quadrature cross-check, (iii) the $1.23\times 10^{-3}$ matching residual at $L=0.95\,\text{m}$, and (iv) the actual measurement sensitivity."*

2. **The $r = L/t$ Compactified Exterior Integral (Paths A and B) and the $L = 0.95\;\text{m}$ Interior/Exterior Switch (`tal_book_HU_2026-10-24.tex`, Appendix, lines 581 and 615):**
   - Line 581:
     > *"The new Path A uses panel integration plus an analytic tail to infinity, while Path B uses the direct substitution $r=L/t$. In the $2.6\text{--}3.01\,\text{m}$ range, the targeted rerun shows the two paths to agree to a **relative precision of about $10^{-15}$**."*
   - Line 615:
     > *"**B. The new, smaller matching residual is reported and genuine.** The report already notes that at $L=0.95\,\text{m}$ the interior/exterior switch is not strictly continuous. The archived values are **$2.17738947$ and $2.17597607\,\text{J/A}^2$**, with difference **$0.00141340\,\text{J/A}^2$**, i.e. **$0.0649\%$** of the interior value. This agrees algebraically with the difference of the total-energy estimates; it is not a persistence of the old large $3\,\text{m}$ error."*

---

## 2. Exact Elliptic Biot–Savart Field of Circular Current Loops and the 8-Coil Superposition

### 2.1. Off-Axis Magnetic Field of a Single Circular Current Loop of Radius $a = D_{\text{turn}}/2$ (`make_revision3_checks.py`, Lines 30–35)
Consider a circular current loop of radius $a = D_{\text{turn}}/2 = 0.2917 / 2 = 0.14585\;\text{m}$ lying in the $xy$-plane and carrying current $I$. In cylindrical coordinates ($\varrho = \sqrt{x^2+y^2}$, $z$), the azimuthal vector potential $A_\varphi(\varrho, z)$ from the Biot–Savart law is:
$$A_\varphi(\varrho, z) \;=\; \frac{\mu_0 I\,a}{4\pi}\int_0^{2\pi}\frac{\cos\phi'\,d\phi'}{\sqrt{a^2 + \varrho^2 + z^2 - 2 a\varrho\cos\phi'}}.$$
With the substitution $\phi' = \pi - 2\theta$ ($\cos\phi' = 2\sin^2\theta - 1$) and the elliptic parameter (`make_revision3_checks.py`, lines 31–32):
$$m(\varrho, z; a) \;\equiv\; \frac{4 a\varrho}{(a + \varrho)^2 + z^2} \;\in\; [0, 1), \qquad d(\varrho, z; a) \;\equiv\; \sqrt{(a + \varrho)^2 + z^2}, \qquad q(\varrho, z; a) \;\equiv\; (a - \varrho)^2 + z^2,$$
the components of $\mathbf{B} = \nabla\times(A_\varphi\hat{\boldsymbol{\varphi}}) = -\frac{\partial A_\varphi}{\partial z}\hat{\boldsymbol{\varrho}} + \frac{1}{\varrho}\frac{\partial(\varrho A_\varphi)}{\partial\varrho}\hat{\mathbf{z}}$ can be expressed in closed form in terms of the elliptic integrals of the first and second kind, $K(m) = \int_0^{\pi/2}(1 - m\sin^2\theta)^{-1/2}d\theta$ and $E(m) = \int_0^{\pi/2}(1 - m\sin^2\theta)^{1/2}d\theta$ (`make_revision3_checks.py`, lines 33–34):
$$\boxed{\;\begin{aligned}
  B_z(\varrho, z; a, I) &\;=\; \frac{\mu_0 I}{2\pi\,d(\varrho,z;a)}\left[K(m) \;+\; \frac{a^2 - \varrho^2 - z^2}{q(\varrho,z;a)}\,E(m)\right], \\[6pt]
  B_\varrho(\varrho, z; a, I) &\;=\; \frac{\mu_0 I\,z}{2\pi\,\varrho\,d(\varrho,z;a)}\left[-K(m) \;+\; \frac{a^2 + \varrho^2 + z^2}{q(\varrho,z;a)}\,E(m)\right] \qquad (\varrho > 0,\;\; B_\varrho(0,z) = 0).
\end{aligned}\;}$$

### 2.2. Superposition of the $N_c = 8 \times N = 2000$ Loops and the Derby–Olbert Current-Sheet Limit (`make_revision3_checks.py`, Lines 25–51)
In the V20 layout, the centers of the $N_c = 8$ solenoids lie on a circle of radius $R_{\text{ax}} = 0.45\;\text{m}$ at vertical coordinate $z_c = 0.200 - 0.300 = -0.100\;\text{m}$:
$$\mathbf{c}_k \;=\; \left(R_{\text{ax}}\cos\frac{2\pi k}{8},\; R_{\text{ax}}\sin\frac{2\pi k}{8},\; z_c\right), \qquad k = 0, 1, \dots, 7.$$
Each coil has length $\ell = 0.400\;\text{m}$ and contains $N = 2000$ turns at midpoint positions $z_n = \bigl(\frac{n + 1/2}{N} - \frac{1}{2}\bigr)\ell$ ($n = 0,\dots,N-1$). At an arbitrary point $\mathbf{P} = (x,y,z)$, the net magnetic field is:
$$\mathbf{B}_{\text{loops}}(\mathbf{P}; I_0) \;=\; \sum_{k=0}^{7}\sum_{n=0}^{N-1}\left[B_\varrho(\varrho_k, z - z_c - z_n; a, I_0)\,\frac{(x - c_{k,x})\hat{\mathbf{x}} + (y - c_{k,y})\hat{\mathbf{y}}}{\varrho_k} \;+\; B_z(\varrho_k, z - z_c - z_n; a, I_0)\,\hat{\mathbf{z}}\right],$$
where $\varrho_k = \sqrt{(x - c_{k,x})^2 + (y - c_{k,y})^2}$. As $N \to \infty$ at fixed surface current density $N I_0 / \ell$, the midpoint sum $\frac{1}{N}\sum_{n=0}^{N-1} f(z_n) = \frac{1}{\ell}\int_{-\ell/2}^{+\ell/2} f(z')\,dz' + \mathcal{O}(N^{-2})$ converges to second order to the continuous Derby–Olbert cylindrical current-sheet integral (`fields.coils_B`). Accordingly, at $I_0 = 14.0\;\text{mA}$, the two independent evaluations (`revision3_checks.json`, lines 124–136) give:
- **At the Main Clock ($\mathbf{P}_{\text{main}} = (0,0,0)$):** $|\mathbf{B}_{\text{loops}}| = 11.63290030\;\mu\text{T}$, $|\mathbf{B}_{\text{sheet}}| = 11.63289980\;\mu\text{T}$ (relative difference $-4.28\times 10^{-8}$).
- **At the Reference Clock ($\mathbf{P}_{\text{ref}} = (1.60, 0, 0)\;\text{m}$):** $|\mathbf{B}_{\text{loops}}| = 0.430539231\;\mu\text{T}$, $|\mathbf{B}_{\text{sheet}}| = 0.430539228\;\mu\text{T}$ (relative difference $-7.56\times 10^{-9}$).

---

## 3. Differential Magnetic Energy ($\Delta W(L)$), the $r = L/t$ Compactified Exterior Quadrature, and the $L = 0.95\;\text{m}$ Identity

### 3.1. Expression of the Polarity-Even Fractional Frequency Shift in Terms of $\Delta W(L)$ (`make_revision3_checks.py`, Lines 53–58)
By definition, the magnetic energy normalized to unit current ($I = 1\;\text{A}$) inside a ball $B(\mathbf{P}, L)$ of radius $L$ centered at $\mathbf{P}$ is:
$$W_{\text{in}}(\mathbf{P}; L) \;\equiv\; \frac{1}{2\mu_0 I_0^2}\int_{|\mathbf{x}-\mathbf{P}|\le L}|\mathbf{B}(\mathbf{x}; I_0)|^2\,d^3x \qquad [\text{J/A}^2],$$
and the differential enclosed energy between the main clock ($\mathbf{P}_{\text{main}} = \mathbf{0}$) and the reference clock ($\mathbf{P}_{\text{ref}}$) is $\Delta W(L) \equiv W_{\text{in}}(\mathbf{P}_{\text{main}}; L) - W_{\text{in}}(\mathbf{P}_{\text{ref}}; L)$. Because the 3D top-hat averaged magnetic mass-density difference is $\Delta\bar{\rho}_{B,L} = \frac{I_0^2\,\Delta W(L)}{c^2\,V_L}$ (where $V_L = \frac{4}{3}\pi L^3$), the conditional H1 clock-shift sensitivity coefficient ($k_g(L)$, $[\text{A}^{-2}]$) and the signal at $I_0 = 14\;\text{mA}$ (`make_revision3_checks.py`, lines 55–57; `revision3_checks.json`, lines 139–145) are:
$$\boxed{\;k_g(L) \;=\; -\frac{G L^2}{c^2}\frac{\Delta W(L)}{c^2\,\frac{4}{3}\pi L^3} \;=\; -\frac{3 G\,\Delta W(L)}{4\pi c^4 L}, \qquad \Delta y(L) \;=\; k_g(L)\,I_0^2.\;}$$
At the grid peak on the 25-point logarithmic grid ($L_p = 10^{-0.2}\;\text{m} = 0.6309573444801932\;\text{m}$):
$$\Delta W(L_p) \;=\; 1.9611445207162488\;\text{J/A}^2 \quad\Longrightarrow\quad k_g(L_p) \;=\; -6.131177672606294\times 10^{-45}\;\text{A}^{-2},$$
$$\Delta y(L_p) \;=\; k_g(L_p)\times (0.014\;\text{A})^2 \;=\; -1.2017108238308338\times 10^{-48}.$$

### 3.2. Why Is the Exterior Complementary Integral Advantageous at Large $L$, and How Does the $r = L/t$ Compactification Work?
When $L > L_{\text{jump}} = 0.95\;\text{m}$ is large enough that both the main-clock ball $B(\mathbf{P}_{\text{main}}, L)$ and the reference-clock ball $B(\mathbf{P}_{\text{ref}}, L)$ enclose the intense near-coil zone, directly subtracting the two interior integrals ($W_{\text{in,main}}(L) - W_{\text{in,ref}}(L)$) would subtract two nearly equal numbers. In Euclidean space, however, the magnetic energy integrated over all of $\mathbb{R}^3$ is **translation-invariant**:
$$W_{\text{tot}} \;\equiv\; \frac{1}{2\mu_0 I_0^2}\int_{\mathbb{R}^3}|\mathbf{B}(\mathbf{x}; I_0)|^2\,d^3x \;=\; W_{\text{in}}(\mathbf{P}_{\text{main}}; L) + W_{\text{ext}}(\mathbf{P}_{\text{main}}; L) \;=\; W_{\text{in}}(\mathbf{P}_{\text{ref}}; L) + W_{\text{ext}}(\mathbf{P}_{\text{ref}}; L),$$
where $W_{\text{ext}}(\mathbf{P}; L) \equiv \frac{1}{2\mu_0 I_0^2}\int_{|\mathbf{x}-\mathbf{P}| \ge L}|\mathbf{B}(\mathbf{x}; I_0)|^2\,d^3x$ is the exterior (complementary) energy outside the ball. Because the exact total energy cancels ($W_{\text{tot,main}}^{\text{exact}} - W_{\text{tot,ref}}^{\text{exact}} \equiv 0$), the interior energy difference equals the negative difference of the smooth, singularity-free exterior integrals:
$$\boxed{\;\Delta W_{\text{out}}(L) \;\equiv\; -\Bigl[W_{\text{ext}}(\mathbf{P}_{\text{main}}; L) \;-\; W_{\text{ext}}(\mathbf{P}_{\text{ref}}; L)\Bigr] \;=\; \frac{1}{2\mu_0 I_0^2}\int_{r \ge L} r^2\,dr\int_{S^2} d\Omega\,\Bigl[|\mathbf{B}(\mathbf{P}_{\text{ref}} + r\hat{\mathbf{n}})|^2 - |\mathbf{B}(\mathbf{P}_{\text{main}} + r\hat{\mathbf{n}})|^2\Bigr].\;}$$

V20 `energy_integrator.py` evaluates this exterior integral over the half-line $[L, \infty)$ along two independent paths (`make_revision3_checks.py`, lines 62–68; `tal_book_HU_2026-10-24.tex`, line 581):
- **Path "B" (`dW_B`): Exact $r = L / t$ Compactification Without Finite Cutoff.**
  Introducing the dimensionless variable $t \equiv L / r \in (0, 1]$, the differential is $dr = -\frac{L}{t^2}\,dt$ and the volume measure is $r^2\,dr = -\frac{L^3}{t^4}\,dt$. Thus the infinite domain $[L, \infty)$ maps onto the compact interval $(0, 1]$:
  $$\boxed{\;\Delta W_{\text{B}}(L) \;=\; \frac{L^3}{2\mu_0 I_0^2}\int_0^1 \frac{dt}{t^4}\int_{S^2} d\Omega\,\Bigl[\bigl|\mathbf{B}(\mathbf{P}_{\text{ref}} + \tfrac{L}{t}\hat{\mathbf{n}})\bigr|^2 \;-\; \bigl|\mathbf{B}(\mathbf{P}_{\text{main}} + \tfrac{L}{t}\hat{\mathbf{n}})\bigr|^2\Bigr].\;}$$
  Let us inspect the behavior of the integrand as $t \to 0^+$ ($r = L/t \to \infty$): because the net charge of the 8 closed current loops is zero, the far magnetic field is dipolar at leading order: $|\mathbf{B}(r\hat{\mathbf{n}})| = \mathcal{O}(r^{-3}) = \mathcal{O}(t^3/L^3)$, so the energy density is $|\mathbf{B}|^2 = \mathcal{O}(t^6/L^6)$, and the difference between the two shifted balls is $|\mathbf{B}(\mathbf{P}_{\text{ref}} + \frac{L}{t}\hat{\mathbf{n}})|^2 - |\mathbf{B}(\mathbf{P}_{\text{main}} + \frac{L}{t}\hat{\mathbf{n}})|^2 = \mathcal{O}(t^7/L^7)$. Multiplied by the Jacobian factor $\frac{L^3}{t^4}$, the integrand vanishes at $t \to 0^+$ **as a regular polynomial ($\mathcal{O}(t^2)$–$\mathcal{O}(t^3)$)**! Consequently, $n_t = 32$ (or $64$) point Gauss–Legendre quadrature on $(0,1]$ achieves machine precision while simultaneously proving the far-field scaling:
  $$\Delta W(L) \;\propto\; L^{-3}\;\text{(or } L^{-5}\text{ under quadrupole/shift symmetry)} \quad\Longrightarrow\quad |\Delta y(L)| \;\propto\; \frac{\Delta W(L)}{L} \;\propto\; L^{-4}\text{--}L^{-6}.$$
- **Path "A" (`dW_A`): Panel Integral + Analytic Asymptotic Multipole Tail.**
  At $L = 3.0\;\text{m}$, the relative difference between Path "A" (`dW_A = 0.0023313149428584538 J/A^2`) and Path "B" (`dW_B = 0.0023313149428584555 J/A^2`) is $-7.44\times 10^{-16}$ (`revision3_checks.json`, lines 164–165; the $3.7\times 10^{-16}$ agreement of the two deterministic routines in `tal_book_HU_2026-10-24.tex`, line 478).

### 3.3. Exact Algebraic Identity of the $L_{\text{jump}} = 0.95\;\text{m}$ Interior/Exterior Switching Jump (`make_revision3_checks.py`, Lines 70–74)
When the interior and exterior quadratures evaluated at finite angular resolution $n_{\text{hi}} = 24000$ are added together, the numerical total energies centered on the two clock positions (`revision3_checks.json`, lines 186–195) are:
$$W_{\text{tot,main}}^{\text{num}} \;=\; W_{\text{in,main}}^{\text{num}}(L) + W_{\text{ext,main}}^{\text{num}}(L) \;=\; 2.237861258956951\;\text{J/A}^2,$$
$$W_{\text{tot,ref}}^{\text{num}} \;=\; W_{\text{in,ref}}^{\text{num}}(L) + W_{\text{ext,ref}}^{\text{num}}(L) \;=\; 2.2364478629240048\;\text{J/A}^2.$$
Subtracting the two equations yields an exact algebraic identity valid for every $L$ (including the switching point $L_{\text{jump}} = 0.95\;\text{m}$):
$$\boxed{\;\Delta W_{\text{in}}^{\text{num}}(0.95\;\text{m}) \;-\; \Delta W_{\text{out}}^{\text{num}}(0.95\;\text{m}) \;\equiv\; W_{\text{tot,main}}^{\text{num}} \;-\; W_{\text{tot,ref}}^{\text{num}} \;=\; 0.001413396032946\;\text{J/A}^2\;}$$
($2.17738947 - 2.17597607 = 0.00141340\;\text{J/A}^2$, i.e., $0.06491\%$ of the interior value; the algebraic difference check residual is `identity_check_abs` $= 8.88\times 10^{-16}\;\text{J/A}^2$).

---

## 4. Checks

1. **Consistency of the DC Coil Resistance with the Measured Winding Diameter $D_{\text{turn}} = 291.7\;\text{mm}$ (`tal_book_HU_2026-10-24.tex`, Line 577):**
   - For $N = 2000$ turns of copper wire of diameter $d_{\text{wire}} = 0.15\;\text{mm}$ ($A_{\text{wire}} = \frac{\pi}{4}d_{\text{wire}}^2 = 1.767146\times 10^{-8}\;\text{m}^2$, $\rho_{\text{Cu}}(20^\circ\text{C}) = 1.7241\times 10^{-8}\;\Omega\,\text{m}$) and mean winding diameter $D_{\text{turn}} = 0.2917\;\text{m}$: $L_{\text{wire}} = \pi D_{\text{turn}} N = 1832.805\;\text{m}$, giving $R_{\text{DC,calc}} = 1788.160\;\Omega$, which differs by only $-0.0470\%$ from the nominal $1789\;\Omega$ $\checkmark$.
2. **Verification of the Peak Signal via Formula (1) at $L_p = 0.63095734448\;\text{m}$:**
   - $k_g(L_p) = -\frac{3\times 6.67430\times 10^{-11}\times 1.961144520716}{4\pi\times (299792458)^4\times 0.63095734448} = -6.1311776726\times 10^{-45}\;\text{A}^{-2}$, and $\Delta y(L_p) = k_g(L_p)\times 0.014^2 = -1.2017108238\times 10^{-48}$, reproducing the main V20/V22 result $\checkmark$.
