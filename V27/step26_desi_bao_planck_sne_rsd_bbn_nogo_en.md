# Step 26 (V27): Cosmological $H(z)$, DESI DR1 BAO, Planck 2018 Acoustic Scale, Synthetic SNe Ia Shape Control, Variational RSD $f\sigma_8$, Algebraic BBN Benchmark and Inverse Distance Ladder No-Go Audit

**Author:** Dávid Szabolcsi (`ORCID: 0009-0005-6300-074X`)  
**Date:** 2026-10-28 (text revised 2026-10-29 after the Academia Co-scientist review, the subsequent re-audit and the second review F01–F12 and the third re-audit N01, N02, T01, see Sections 11, 12 and 13)  
**Preceding version DOIs:**

- **V24:** Zenodo DOI `10.5281/zenodo.23680715` (strictly V24)
- **V25:** Zenodo DOI `10.5281/zenodo.23699685` (strictly V25, `step1`–`step22`)
- **V26:** Zenodo DOI `10.5281/zenodo.23727585` (strictly V26, `step1`–`step25`)

**Source of the numbers:** computed results are generated from the JSON/CSV sources (`tal_hubble_desi_bao_summary.json` and the Step 26 CSV tables, produced with the shared background module `tal_cosmology_common.py`) by `generate_step26_markdown.py`. Fixed input constants (data values, priors, literature benchmarks, external inputs) are stated in the code and cited in the text.

---

## 0. Scope and summary

This module tests one specific branch: the V26 `Step 25` M4 background at fixed $L=10\;\text{km}$, $\alpha_Y=1$, $\beta_{\rm Yuk}=1/(4\pi)=0.079577$, on the adiabatic particular solution of the scalar equation with negligible free scalar oscillation energy. The narrower, defensible statement is: **in this branch the bare coupling $\beta_{\rm Yuk}$ inserted as a multiplier of $H_0$ is not a derived mechanism for the shift $67.4\to 72.8\;\text{km/s/Mpc}$.** This is not an exclusion of every possible cosmological solution of TAL.

The revision implements the review items as follows:

- one action, one scalar equation (SI units, cosmic seconds), one attractor formula and one $G_{\rm eff}$ normalization, identical in text and code (Section 1);
- explicit branch selection (adiabatic attractor, free oscillation amplitude $C=0$) instead of a claim of uniqueness (Section 1.3);
- DESI DR1 BAO with the full-precision Cobaya 12-vector and $12\times 12$ covariance, plus profiled controls (Section 2);
- model-specific re-integrated sound horizons $r_s(z_*)$ and $r_d(z_{\rm drag})$ for every CMB model and for the M3 control; the earlier $\theta_*$ and $r_d$ sigma exclusion labels are removed (Sections 3 and 8);
- the supernova part is renamed a **synthetic distance-modulus shape control** (Section 4);
- RSD: diagonal benchmark plus a first-order variational $\delta f\sigma_8$ (Section 5);
- BBN: algebraic order-of-magnitude benchmark only (Section 6);
- local environments: constant-density bulk benchmarks with `expm1` (Section 7);
- the re-audit items R1–R6 and 5.1–5.3 are mapped to the changes in the revision log (Section 11); the second review items F01–F12 are mapped in Revision log 2 (Section 12), the third re-audit items N01, N02, T01 in Revision log 3 (Section 13), and the open limitations P01–P09 are listed in Section 9.1.

---

## 1. Unified action, scalar equation and $G_{\rm eff}$ normalization

### 1.1 Action and sign convention

With signature $(-,+,+,+)$, Einstein-frame metric $g_{\mu\nu}$ and Jordan-frame matter metric $\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}$, $A(\phi)=e^{+\alpha_Y\phi}$:

$$
S=\int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}R-\frac{c^4}{4\pi G\beta_{\rm Yuk}}\left(\frac{1}{2}g^{\mu\nu}\partial_\mu\phi\,\partial_\nu\phi+\frac{\phi^2}{2L^2}\right)-\rho_{\rm vac}^{\rm eff}c^2\right]+S_m\!\left[\Psi_m,\,A^2(\phi)g_{\mu\nu}\right].
$$

With $T^{\mu\nu}=(2/\sqrt{-g})\,\delta S_m/\delta g_{\mu\nu}$ and $T=-\rho_m c^2+3p_m$, variation with respect to $\phi$ gives (units m$^{-2}$):

$$
\Box\phi-\frac{\phi}{L^2}=-\frac{4\pi G}{c^4}\,\alpha_Y\beta_{\rm Yuk}\,T=+\frac{4\pi G}{c^2}\,\alpha_Y\beta_{\rm Yuk}\left(\rho_m-\frac{3p_m}{c^2}\right).
$$

### 1.2 Homogeneous FLRW equation (cosmic seconds, SI)

With $\Box\phi=-c^{-2}(\ddot\phi+3H\dot\phi)$ for a homogeneous field, both sides now have units s$^{-2}$:

$$
\ddot{\bar\phi}+3H\dot{\bar\phi}+\frac{c^2}{L^2}\bar\phi=-4\pi G\,\alpha_Y\beta_{\rm Yuk}\left(\rho_m-\frac{3p_m}{c^2}\right).
$$

The scalar Compton rate is $c/L=29979.2458\;\text{s}^{-1}$ and $H_0L/c=7.286\times 10^{-23}$.

### 1.3 Attractor and branch selection

Because $4\pi\beta_{\rm Yuk}=1$, the adiabatic particular solution ($\ddot{\bar\phi},\,3H\dot{\bar\phi}\ll (c/L)^2\bar\phi$) is

$$
\bar\phi_{\rm attr}(z)=-\alpha_Y\gamma_L\left(\rho_m(z)-\frac{3p_m(z)}{c^2}\right),\qquad \gamma_L=\frac{4\pi\beta_{\rm Yuk}GL^2}{c^2}=\frac{GL^2}{c^2}=7.42616\times 10^{-20}\;\text{m}^3/\text{kg},
$$

so that $\bar\phi_{\rm attr}(0)=-1.9949\times 10^{-46}$ and $\bar\phi_{\rm attr}(z_*)\approx-2.589\times 10^{-37}$. The general homogeneous solution also contains a free oscillation

$$
\delta\phi\simeq C\,a^{-3/2}\cos\!\left[(c/L)\,t+\vartheta\right],\qquad \langle\rho_{\delta\phi}\rangle\propto a^{-3}.
$$

$H/(c/L)\ll 1$ justifies the adiabatic approximation but does not remove this mode. **This analysis therefore selects the branch $C=0$ ($\Omega_{\delta\phi,0}=0$): the adiabatic particular solution with negligible free oscillation energy.** No initial-value integration of the scalar equation is performed, and the conclusions below hold only for this branch.

### 1.4 Scalar energy density

From the action, $\rho_\phi=\dfrac{c^2}{8\pi G\beta_{\rm Yuk}}\left[(\dot\phi/c)^2+(\phi/L)^2\right]$ (mass density), hence on the attractor

$$
\Omega_{\phi,\rm attr}(0)=\frac{(c/L)^2\,\bar\phi_{\rm attr}(0)^2}{3\beta_{\rm Yuk}H_0^2}=3.140\times 10^{-47}.
$$

### 1.5 $G_{\rm eff}$ normalization

For a static point mass, $\phi(r)=-\alpha_Y\beta_{\rm Yuk}\,\dfrac{GM}{c^2r}\,e^{-r/L}<0$, and the matter coupling $A=e^{\alpha_Y\phi}$ yields the fifth-force strength $\alpha_Y^2\beta_{\rm Yuk}$. For linear perturbations with comoving wavenumber $k$:

$$
\frac{G_{\rm eff}(k,a)}{G_N}=1+\alpha_Y^2\beta_{\rm Yuk}\,\frac{(kL/a)^2}{1+(kL/a)^2}.
$$

This canonical normalization is used throughout. The $2\alpha_Y^2\beta_{\rm Yuk}$ row is a separate, double-strength sensitivity control; it is not the growth law of the published V26 Step25 (Section 5). Local proper time is $d\tau_J=A(\phi)\,d\tau_E=e^{-\alpha_Y|\phi|}d\tau_E$.

**G dictionary (review F02).** Four symbols are kept apart:

| Symbol | Meaning |
|:-----|:-----------------------------|
| $G_*$ | the $G$ in the action of Section 1.1 (Einstein frame) |
| $G_T$ | Jordan-frame tensor coupling: $g_J=A^2g_E$ gives the Ricci term $A^{-2}R/(16\pi G_*)$, so $1/G_T=A^{-2}/G_*$, $G_T=G_*A^2$ ($\phi<0$ decreases $G_T$) |
| $G_{\rm lab}$ | operationally measured $G$: for $r\ll L$ and unscreened small bodies $\Phi=-G_*M(1+q\,e^{-r/L})/r$, so $G_{\rm lab}\approx G_*(1+q)$, $q=\alpha_Y^2\beta_{\rm Yuk}=0.0795775$ |
| $G_N$ | the symbol in the perturbation equation of Section 5: $G_N:=G_*$ |

On scales $kL/a\ll1$ the growth equation therefore uses $G_{\rm eff}/G_{\rm lab}=1/(1+q)=0.9262883$ if the coupling is calibrated by a laboratory value; the $10^{-41}$-level statement of Section 5 is a formal control at the same bare $G_*$. No $f\sigma_8$ exclusion claim is made and no screening branch is introduced.

### 1.6 Background normalization and observer frame

All scripts use one background module (`tal_cosmology_common.py`). The smooth massive-neutrino interpolation alone gives $E^2(0)=1.00000013546$ (offset $1.3546\times 10^{-7}$); the dark-energy closure absorbs this offset so that $E^2(0)=1$ exactly and $H(0)=H_0$. This offset ($1.3546\times 10^{-7}$ in $\Omega_\nu$) is the $z=0$ normalization drift of the interpolation, not its approximation error. Against the Fermi–Dirac energy integral of the massive species (matched to the same relativistic and non-relativistic limits), the $p=1.83$ interpolation has a relative error of $+1.149\times 10^{-3}$ in the massive-neutrino density and $+1.795\times 10^{-5}$ in the total $H$ at $z=1000$; the largest absolute values on the grid $0\le z\le10^4$ are $3.435\times 10^{-3}$ (massive component, $z\approx33.5$) and $1.795\times 10^{-5}$ (total $H$, $z\approx1000$). This is an approximation error of the $\Lambda$CDM-like background, common to all compared models, and separate from the theoretical TAL correction ($\sim 10^{-46}$).

| $z$ | $\rho_\nu^{\rm FD}$ [$E^2$ units] | $\rho_\nu^{p=1.83}$ [$E^2$ units] | rel. error (massive) | rel. error ($H$) |
|-----:|-------:|-------:|------:|------:|
| 0 | $1.41991\times 10^{-3}$ | $1.41998\times 10^{-3}$ | $+4.449\times 10^{-5}$ | $+3.159\times 10^{-8}$ |
| 1 | $1.13610\times 10^{-2}$ | $1.13626\times 10^{-2}$ | $+1.355\times 10^{-4}$ | $+2.402\times 10^{-7}$ |
| 10 | $1.90137\times 10^{0}$ | $1.90428\times 10^{0}$ | $+1.529\times 10^{-3}$ | $+3.454\times 10^{-6}$ |
| 100 | $2.02529\times 10^{3}$ | $2.02475\times 10^{3}$ | $-2.670\times 10^{-4}$ | $-8.115\times 10^{-7}$ |
| 1000 | $1.27145\times 10^{7}$ | $1.27291\times 10^{7}$ | $+1.149\times 10^{-3}$ | $+1.795\times 10^{-5}$ |
| 3000 | $1.01910\times 10^{9}$ | $1.01945\times 10^{9}$ | $+3.388\times 10^{-4}$ | $+1.083\times 10^{-5}$ |

 The computations use $\Lambda$CDM-like $H(z)$ functions with a separate clock factor; the Jordan-frame mapping ($dt_J=A\,dt_E$, $a_J=A\,a_E$, $H_J=(H_E+\alpha_Y\dot\phi_E)/A$) is a leading-order identity at this $|\phi|$, not a numerically solved coupled Friedmann–scalar–matter system.

### 1.7 Version reconciliation of the kinetic normalization (re-audit R1)

Paths below are relative to the V26 tree `scalar-tensor-yukawa-fit/` (V26, DOI `10.5281/zenodo.23727585`).

- The first Step 26 draft wrote the kinetic term as $\frac{c^4}{16\pi G}\frac{1}{2\beta_{\rm Yuk}}(\partial\phi)^2=\frac{c^4}{32\pi G\beta_{\rm Yuk}}(\partial\phi)^2$, while the action of Section 1.1 has $\frac{c^4}{8\pi G\beta_{\rm Yuk}}(\partial\phi)^2$: $K_{\rm old}/K_{\rm new}=0.25$.
- V26 uses $\frac{c^4}{16\pi G\tilde\alpha_L}(\partial\phi_L)^2$ with $\tilde\alpha_L=\alpha_Y/(8\pi)$, i.e. $c^4/(2G\alpha_Y)$ (`step23_gw170817_and_lensing_en.md:19` and `:29`, `README_en.md:50`, `gw170817_and_lensing_audit.py:104–113`). At $\alpha_Y=1$, $\beta_{\rm Yuk}=1/(4\pi)$ this equals the new coefficient ($K_{\rm V26}/K_{\rm new}=1.000000$), and the V26 fifth-force strength $2\tilde\alpha_L\alpha_Y^2=0.079577$ equals $\alpha_Y^2\beta_{\rm Yuk}$. The V26 scalar equation $(\Box-L^{-2})\phi=-(\alpha_Y^2G/c^4)\varepsilon_m$ (`step25_cmb_and_cosmology_en.md:37`) and $\phi_{\rm attr}=\alpha_Y^2\gamma_L\rho$ (`cmb_and_cosmology_audit.py:327`, `:620`) agree with Section 1.3 up to the sign convention $\phi_{\rm V26}=-\phi_{\rm V27}$. **Conclusion:** V26 used the new normalization; the $c^4/(32\pi G\beta)$ form was a transcription inconsistency of the first draft, and the computed numbers already followed the new normalization. For $\alpha_Y\neq1$ the two versions differ, since V26 ties $\beta_{\rm Yuk}=\alpha_Y^3/(4\pi)$ whereas here $\beta_{\rm Yuk}$ is an independent parameter.
- The field redefinition $\phi_{\rm old}=2\phi_{\rm new}$ maps the old kinetic term onto the new one, but $A=e^{\alpha_{\rm old}\phi_{\rm old}}=e^{2\alpha_{\rm old}\phi_{\rm new}}$, so it requires $\alpha_{\rm new}=2\alpha_{\rm old}$. With the same $\alpha_Y=1$ the old form would give the fifth-force strength $4\alpha_Y^2\beta_{\rm Yuk}=0.318310$ instead of $\alpha_Y^2\beta_{\rm Yuk}=0.079577$.
- Carried over unchanged: all background, BAO, CMB, SNe-shape and RSD numbers (the background scalar-energy correction is $\Omega_\phi\sim10^{-47}$; the direct Yukawa perturbation of RSD is a separate negligible term quantified in Section 5: $\delta G_{\rm eff}/G_*\approx3.8\times 10^{-41}$, $\delta f\sigma_8\approx1.67\times 10^{-40}$), $\bar\phi_{\rm attr}(0)$, $\Omega_{\phi,\rm attr}(0)$ and the bulk local benchmarks at $\alpha_Y=1$. The V26 background additionally used a $(\phi/4\pi)^2$ energy convention (`cmb_and_cosmology_audit.py:330`, `:623`), giving $\Omega_\phi=1.582\times 10^{-50}$ instead of the action-consistent $3.140\times 10^{-47}$; both are negligible.
- The factor-2 growth convention $2\alpha_Y^2\beta_{\rm Yuk}$ is relabelled a **sensitivity control**, not an equivalent convention.
- The V21 top-hat interior external benchmark $|\phi|=0.0412$ (Section 7) is an interior value of the V21 top-hat fixed-point profile $x=L^2(G/c^2)\bar\rho$ (`stellar_tal_v21.py:140–160`, `data/v21/profiles_L10km.csv`: SLy, mode `gradient`, $r=2.0046409680935633$ km, $x=0.041239199439254325$; APR4, mode `mass`, $r=2.75417748649486$ km, $x=0.041233030079658145$), not a surface value and **not automatically a solution of the new Helmholtz–Yukawa action**; it remains an external input.

---

## 2. DESI DR1 BAO: $12\times 12$ covariance and profiled controls

The 12 primary observables are the Cobaya DESI DR1 vector and covariance (`desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`, `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`). The separate LRG3 and ELG1 single-tracer rows are diagnostics only and are not added to the likelihood; the $z=0.845$ ELG1 row ($D_V/r_d=18.33\pm0.62$) is a **chosen illustration without a verified primary source** and not a DESI measurement (it is not in DESI III Table 18 or DESI VI Table 1). Per-bin $\chi^2$ for primary bins use the $12\times 12$ covariance; diagnostic rows (italic) use their rounded table errors. Predictions in the table are $D_M/r_d$; $D_H/r_d$ for anisotropic bins and $D_V/r_d$ for isotropic bins.

| Tracer | $z_{\rm eff}$ | DESI DR1 ($X/r_d$) | M4 prediction | $\chi^2_{\rm bin}$ M4 | $\chi^2_{\rm bin}$ 2P-prof. | $\chi^2_{\rm bin}$ rigid 72.8 |
|:-----------|-----:|:------------------------------------------|:--------------------|-----:|-----:|-----:|
| **BGS** | 0.295 | $D_V$: $7.93\pm0.15$ | $8.053$ | 0.724 | 0.090 | 9.69 |
| **LRG1** | 0.510 | $D_M$: $13.62\pm0.25$; $D_H$: $20.98\pm0.61$ ($r=-0.445$) | $13.495$; $22.733$ | 8.963 | 6.010 | 24.43 |
| **LRG2** | 0.706 | $D_M$: $16.85\pm0.32$; $D_H$: $20.08\pm0.60$ ($r=-0.420$) | $17.694$; $20.164$ | 8.967 | 3.097 | 12.92 |
| **LRG3+ELG1** | 0.930 | $D_M$: $21.71\pm0.28$; $D_H$: $17.88\pm0.35$ ($r=-0.389$) | $21.917$; $17.608$ | 0.826 | 2.031 | 75.07 |
| *LRG3 (sep.) (diagn.)* | 0.922 | $D_M$: $21.81\pm0.31$; $D_H$: $17.83\pm0.38$ ($r=-0.393$) | $21.775$; $17.692$ | *0.209* | *2.931* | *69.61* |
| *ELG1 (z=0.955) (diagn.)* | 0.955 | $D_V$: $20.01\pm0.41$ | $20.229$ | *0.285* | *0.004* | *9.77* |
| *ELG1 (chosen illustration, not a DESI measurement)* | 0.845 | $D_V$: $18.33\pm0.62$ | $18.666$ | *0.294* | *0.021* | *2.86* |
| **ELG2** | 1.317 | $D_M$: $27.79\pm0.69$; $D_H$: $13.82\pm0.42$ ($r=-0.444$) | $28.017$; $14.094$ | 0.883 | 0.380 | 18.57 |
| **QSO** | 1.491 | $D_V$: $26.07\pm0.67$ | $26.026$ | 0.005 | 0.141 | 8.74 |
| **Ly$\alpha$ QSO** | 2.330 | $D_M$: $39.71\pm0.94$; $D_H$: $8.52\pm0.17$ ($r=-0.477$) | $39.173$; $8.615$ | 0.418 | 0.990 | 44.93 |
| **Sum, 12 primary observables (12×12 cov.)** | | | | **20.785** | **12.741** | **194.34** |

- Fixed M4 point ($H_0=67.40$, $r_d=147.0922\;\text{Mpc}$): $\chi^2=20.7852$ for 12 observables, $p=0.0536$ (rounded Table 1 input: $20.6996$).
- Rigid $H_0=72.8$ with unchanged shape and $r_d$: $\chi^2=194.342$, $p=4.79\times 10^{-35}$ (rounded Table 1 input: $195.178$).

**Profiled controls** (local Fisher errors; $p$-values are conditional frequentist values for the stated Gaussian likelihood, not probabilities of TAL being true):

| Control | $\Omega_m$ | $P=H_0r_d$ [km/s] | $\chi^2$ (12×12) | dof | $p$ | $r_d$ needed at $H_0=72.8$ [Mpc] |
|:------------------------------------|-----:|------------:|------:|---:|------:|------------:|
| Fixed M4 point (no free BAO parameter) | 0.3148 | 9914.02 | 20.7852 | 12 | 0.0536 | 136.18 |
| M4 shape, rigid $H_0=72.8$, unchanged $r_d$ | 0.3148 | 10708.32 | 194.3420 | 12 | $4.79\times 10^{-35}$ | — (147.09 fixed) |
| M4 shape, $P=H_0r_d$ profiled | 0.3148 | 10029.65 $\pm$ 47.43 | 14.7015 | 11 | 0.1966 | 137.77 $\pm$ 0.65 |
| Flat $\Lambda$CDM shape, $\Omega_m$ and $P$ profiled | 0.2939 $\pm$ 0.0145 | 10193.67 $\pm$ 126.47 | 12.7405 | 10 | 0.2385 | 140.02 $\pm$ 1.74 |

- At fixed $E(z)$, BAO constrains only $P=H_0r_d$ (rank-1 information in the $(H_0,r_d)$ plane).
- Two-parameter control: $\Omega_m=0.2939\pm0.0145$, $P=10193.67\pm126.47\;\text{km/s}$, $\chi^2_{\min}=12.7405$ (10 dof, $p=0.2385$).
- $r_d$ required at $H_0=72.8$ exactly: $140.02\pm1.74\;\text{Mpc}$; propagating $H_0=72.8\pm1.0$: $\pm2.59\;\text{Mpc}$. Here $72.8\pm1.0$ is an **assumed local $H_0$ benchmark prior** (a chosen target value, review F05); it is not the Riess et al. 2022 value $73.04\pm1.04$ and no specific SH0ES release is claimed for it.
- With the model-dependent Planck prior $r_d=147.09\pm0.26\;\text{Mpc}$: $H_0=69.30\pm0.87\;\text{km/s/Mpc}$.
- Joint illustrative test (BAO + Planck $r_d$ prior + local $H_0=72.8\pm1.0$, profiling $\Omega_m,H_0,r_d$): best fit $\Omega_m=0.2719$, $H_0=70.806$, $r_d=147.025\;\text{Mpc}$, $\chi^2_{\rm joint}=19.713$, split as $\chi^2_{\rm BAO}=15.67431585$, $\chi^2_{r_d\,\rm prior}=0.06234684$, $\chi^2_{H_0\,\rm prior}=3.97660606$; $\Delta\chi^2=6.973$ relative to the BAO-only minimum ($\sqrt{\Delta\chi^2}=2.64$). This is not a Planck+TAL exclusion significance.

![DESI DR1 BAO, Step 25 M4 and profiled controls (`step1_1_desi_dr1_bao_plot_en.png`)](step1_1_desi_dr1_bao_plot_en.png){ width=96% }

---

## 3. Planck 2018 acoustic scale with model-specific sound horizons

For every model, $r_s(z_*=1089.80)$ and $r_d(z_{\rm drag}=1059.94)$ are re-integrated with that model's own $H(z)$ and $c_s(z)$ (fixed $z_*$, $z_{\rm drag}$ as a control approximation). The 3D Gaussian prior on $(\omega_b,\omega_{cb},100\theta_*)$ is an **illustrative assumed prior** (Planck 2018 marginal means and errors with an approximate correlation matrix whose chain provenance is not documented here); it is not the Planck TT,TE,EE+lowE+lensing likelihood. Reference $100\theta_*=1.04109\pm0.00030$.

| Model | $H_0$ | $\Omega_m$ | $\Omega_K$ | $\omega_{cb}$ | $r_s(z_*)$ [Mpc] | $r_d(z_{\rm drag})$ [Mpc] | $100\theta_*$ | $\chi^2_{\rm 3D,illustr.}$ | $\chi^2_{\rm BAO}$ ($r_d$ re-int.) |
|:------------------------------|----:|------:|------:|------:|-------:|-------:|-------:|---------:|-------:|
| Step 25 M4 / TAL | 67.40 | 0.3148 | 0.0000 | 0.14237 | 144.443 | 147.092 | 1.04132 | 0.757 | 20.79 |
| Planck 2018 reference | 67.36 | 0.3152 | 0.0000 | 0.14237 | 144.443 | 147.092 | 1.04120 | 0.173 | 21.19 |
| All-$z$ rigid $H\to sH$ (kinematic control) | 72.80 | 0.3148 | 0.0000 | 0.14237 | 133.729 | 136.182 | 1.04132 | 0.757 | 20.79 |
| *$D_M/s$ with fixed early $r_s,r_d$ (inconsistent geometric illustration, not a model)* | 72.80 | 0.3148 | 0.0000 | 0.14237 | 144.443 | 147.092 | 1.12475 | 101040.4 | 194.34 |
| Flat $H_0=72.8$, same $\Omega_m$ | 72.80 | 0.3148 | 0.0000 | 0.16620 | 138.879 | 141.375 | 1.08075 | 24800.4 | 43.73 |
| Flat $H_0=72.8$, $\Omega_m=0.244828$ | 72.80 | 0.2448 | 0.0000 | 0.12911 | 147.903 | 150.651 | 1.04174 | 196.000 | 70.30 |
| Flat $H_0=72.8$, $\theta_*$-matched $\Omega_m\approx0.2442$ | 72.80 | 0.2442 | 0.0000 | 0.12877 | 147.996 | 150.747 | 1.04132 | 209.592 | 70.84 |
| Open $H_0=72.8$, fixed $\omega_{cb}$ | 72.80 | 0.2698 | 0.0098 | 0.14237 | 144.442 | 147.092 | 1.04132 | 0.727 | 61.21 |

- Step 25 M4: $r_s=144.4428\;\text{Mpc}$, $r_d=147.0922\;\text{Mpc}$, $100\theta_*=1.041320$ (pull $0.77$), $\chi^2_{\rm 3D}=0.757$.
- All-$z$ rigid rescaling $H\to sH$ with unchanged $c_s$: $r_s\to r_s/s=133.7286\;\text{Mpc}$, $r_d=136.1816\;\text{Mpc}$ and $D_M\to D_M/s$, so $100\theta_*=1.041320$ is unchanged. Scaling $D_M$ alone would give $1.124750$; that inconsistent combination is no longer used as a test.
- Flat $H_0=72.8$, $\Omega_m=0.244828$: re-integration gives $r_s=147.9028$, $r_d=150.6513\;\text{Mpc}$ and $100\theta_*=1.041745$ (with the fixed early $r_s$ it would be $1.017374$). This does **not** mean the model passes Planck: $\omega_{cb}=0.12911$ is offset from the prior and the BAO $\chi^2$ with the re-integrated $r_d$ is $70.30$ (with the old fixed $r_d$: $29.30$). Early and late parts must be refitted together.
- **Two $z_*$-shift diagnostics (review F12).** With the Hu–Sugiyama fit shift $z_*\to z_*+\Delta z_*$ ($z_*+\Delta z_*=1088.6345$ for this model), the CSV keeps two separately named values: `theta_star_100_rs_only_shifted_fixed_DM` $=100\,r_s(z_*+\Delta z_*)/D_M(1089.80)=1.04248582$ (only the $r_s$ endpoint shifted) and `theta_star_100_both_endpoints_shifted` $=100\,r_s(z_*+\Delta z_*)/D_M(z_*+\Delta z_*)=1.04250287$. Neither is a recombination solution.
- **Late transition (re-audit R2).** The earlier "late-time-only" row combined $D_M(z_*)/s$ with unchanged early $r_s,r_d$ ($100\theta_*=1.124750$); it has no transition redshift and is **reclassified as an inconsistent geometric illustration** (italic row in the table), removed from the consistent-model statements. A defined control is the piecewise $H(z)=s\,H_{\rm M4}(z)$ for $z<z_t$ and $H_{\rm M4}(z)$ for $z\ge z_t$, $s=72.8/67.4=1.080119$, with $D_M(z_*)$ integrated numerically and $r_s=144.442796\;\text{Mpc}$ unchanged (M4: $D_M(z_*)=13871.1192\;\text{Mpc}$). $H$ jumps at $z_t$, so this is a kinematic control, not a dynamical model; the pull column is a formal 1D distance from the Planck value, not a likelihood significance. For $z_t\ge2.5$ all DESI bins lie in the rescaled region, so the BAO $\chi^2$ equals the rigid value $194.342$.

| $z_t$ | $D_M(z_*)$ [Mpc] | $100\theta_*$ | pull vs $100\theta_*^{\rm Pl}$ (formal) |
|-----:|-----------:|---------:|--------:|
| 2.5 | 13428.2963 | 1.075660 | 115.2 |
| 10 | 13156.6704 | 1.097867 | 189.3 |
| 100 | 12920.1717 | 1.117963 | 256.2 |
| 1000 | 12843.5875 | 1.124630 | 278.5 |
| $D_M/s$ (no $z_t$; inconsistent illustration) | 12842.2175 | 1.124750 | — |

---

## 4. Synthetic distance-modulus shape control (not a Pantheon+ / Union3 likelihood)

This section uses **no** observed SN catalogue: 18 redshift nodes ($z=0.015$–$2.26$) with representative errors $\sigma_\mu$ and a synthetic reference curve from a flat model with $H_0=72.8$, $\Omega_m=0.334$. There are no SN identifiers, observed $\mu$ values or statistical+systematic covariance. The $\chi^2$ values below are model-versus-model shape metrics, not observational goodness of fit.

- For identical $E(z)$ geometry the offset would be the constant identity $5\log_{10}(72.8/67.4)=+0.167357\;\text{mag}$ (positive sign: M4 minus 72.8 reference).
- The tabulated lowest-node difference $0.16782511\;\text{mag}$ at $z=0.015$ is a model difference (the reference has $\Omega_m=0.334$), not this identity.
- **Symmetric zero-point treatment (re-audit R5):** both models receive the same profiled constant offset (absolute-magnitude marginalization), weights $1/\sigma_\mu^2$:

| Model | best offset [mag] | $\chi^2$ without offset | $\chi^2$ with offset (17 dof) |
|:--------------------|------:|------:|------:|
| M4 ($H_0=67.4$) | 0.17595029 | 522.00704 | 0.916614 |
| Flat $H_0=72.8$, $\Omega_m=0.2448$ | 0.04159799 | 51.90332 | 22.777527 |

- The large drop for M4 only reflects that its shape is close to the reference after removing the constant; node positions and representative errors are synthetic, so changes of nodes or errors do not imply a physical improvement.
- The figure shows both models with and without the profiled constant zero-point offset (low-$\Omega_m$ with offset: $\chi^2=22.777527$); $N=18$ nodes, $\nu=17$ degrees of freedom after one profiled offset. The offset is a profiled nuisance constant, not an observational calibration (review F09).
- Uncalibrated SNe alone do not measure $H_0$ ($H_0$–$M_B$ degeneracy). A real catalogue statement requires reading the Pantheon+ or Union3 likelihood and covariance, which this package does not do.

![Synthetic SNe Ia distance-modulus shape control (`step1_3_pantheon_sne_plot_en.png`)](step1_3_pantheon_sne_plot_en.png){ width=96% }

| $z$ node | $\sigma_\mu$ [mag] | $\mu_{\rm ref}$ (72.8, 0.334) | $\mu_{\rm M4}$ (67.4) | $\Delta\mu_{\rm M4-ref}$ | $\Delta\mu_{\rm M4,shift-ref}$ | $\Delta\mu_{0.2448-{\rm ref}}$ | $\Delta\mu_{0.2448,{\rm shift-ref}}$ |
|------:|------:|----------:|----------:|---------:|---------:|---------:|---------:|
| 0.015 | 0.042 | 33.9780 | 34.1459 | +0.1678 | -0.0081 | +0.0022 | -0.0394 |
| 0.023 | 0.038 | 34.9189 | 35.0869 | +0.1681 | -0.0079 | +0.0033 | -0.0383 |
| 0.035 | 0.035 | 35.8493 | 36.0177 | +0.1684 | -0.0075 | +0.0051 | -0.0365 |
| 0.050 | 0.032 | 36.6467 | 36.8156 | +0.1689 | -0.0070 | +0.0073 | -0.0343 |
| 0.075 | 0.030 | 37.5641 | 37.7338 | +0.1697 | -0.0063 | +0.0109 | -0.0307 |
| 0.100 | 0.028 | 38.2245 | 38.3950 | +0.1704 | -0.0055 | +0.0145 | -0.0271 |
| 0.150 | 0.026 | 39.1726 | 39.3446 | +0.1720 | -0.0040 | +0.0217 | -0.0199 |
| 0.200 | 0.025 | 39.8603 | 40.0337 | +0.1734 | -0.0025 | +0.0287 | -0.0129 |
| 0.250 | 0.025 | 40.4036 | 40.5784 | +0.1749 | -0.0011 | +0.0357 | -0.0059 |
| 0.300 | 0.026 | 40.8543 | 41.0306 | +0.1763 | +0.0003 | +0.0425 | +0.0009 |
| 0.400 | 0.028 | 41.5785 | 41.7575 | +0.1789 | +0.0030 | +0.0556 | +0.0140 |
| 0.500 | 0.030 | 42.1508 | 42.3322 | +0.1814 | +0.0055 | +0.0679 | +0.0263 |
| 0.650 | 0.033 | 42.8339 | 43.0187 | +0.1848 | +0.0088 | +0.0849 | +0.0433 |
| 0.800 | 0.036 | 43.3808 | 43.5685 | +0.1877 | +0.0118 | +0.0999 | +0.0583 |
| 1.000 | 0.042 | 43.9724 | 44.1634 | +0.1910 | +0.0151 | +0.1173 | +0.0757 |
| 1.300 | 0.055 | 44.6699 | 44.8649 | +0.1950 | +0.0191 | +0.1384 | +0.0968 |
| 1.700 | 0.075 | 45.3817 | 45.5806 | +0.1989 | +0.0229 | +0.1593 | +0.1177 |
| 2.260 | 0.110 | 46.1321 | 46.3347 | +0.2026 | +0.0267 | +0.1796 | +0.1380 |

---

## 5. RSD $f\sigma_8(z)$: diagonal benchmark and variational TAL correction

The growth equation $\delta_m''+\left(2+\frac{d\ln H}{d\ln a}\right)\delta_m'-\frac{3}{2}\Omega_{cb}(a)\frac{G_{\rm eff}(k,a)}{G_N}\delta_m=0$ is integrated with DOP853 (`rtol=1e-12`, `atol=1e-14`, dense output); $\sigma_{8,0}=0.8111$ is a fixed Planck input, not a TAL-derived amplitude.

- At $k=0.1\,h\,\text{Mpc}^{-1}$: $kL=2.18429\times 10^{-20}$ ($(kL)^2=4.7711\times 10^{-40}$) at $z=0$; $\Delta G_{\rm eff}/G_N=3.7967\times 10^{-41}$ (canonical $\alpha_Y^2\beta_{\rm Yuk}$) or $7.5934\times 10^{-41}$ (factor-2 convention).
- In float64, $1+O(10^{-41})$ equals 1, so the direct difference $f\sigma_8^{\rm TAL}-f\sigma_8^{\rm GR}$ is $0$ by floating-point identity, not by a 40-digit measurement. The first-order variational equation gives $\max\lvert\delta f\sigma_8\rvert=1.666\times 10^{-40}$ (canonical) and $3.333\times 10^{-40}$ (factor-2 sensitivity control) on the 12 points; an independent rescaled check in the re-audit confirmed $1.666\times10^{-40}$ (maximum $1.6663681521\times10^{-40}$, relative deviation $\le5.1\times10^{-6}$). The float64 branch uses $G_{\rm eff}/G_N=1+\Delta G_{\rm eff,can}/G_N$ (re-audit 5.1). The figure legend shows the canonical maximum $1.6664\times 10^{-40}$; the factor-2 value is only the separate sensitivity control (review F07). Here $G_N:=G_*$ (G dictionary, Section 1.5).
- Diagonal $\chi^2$ (no bin covariances, overlaps or AP treatment): M4 with clustering source $\Omega_{cb}$: $10.1233$; with all matter: $10.2176$; flat $H_0=72.8$, $\Omega_m=0.2448$: $9.4350$ (12 points). This is a benchmark, not a covariance likelihood.

![RSD $f\sigma_8(z)$ diagonal benchmark (`step1_4_rsd_fsigma8_plot_en.png`)](step1_4_rsd_fsigma8_plot_en.png){ width=96% }

| Survey | $z$ | $f\sigma_8$ obs. | $f\sigma_8$ M4 | pull M4 | $f\sigma_8$ (72.8, 0.2448) | $\lvert\delta f\sigma_8\rvert$ variational | $\Delta G_{\rm eff}/G_N$ |
|:------------------------------|-----:|-----------:|-------:|------:|-------:|----------:|----------:|
| 6dFGS (Beutler+ 2012) | 0.067 | $0.423\pm0.055$ | 0.4424 | +0.35 | 0.3887 | $7.57\times 10^{-41}$ | $4.32\times 10^{-41}$ |
| SDSS MGS (Howlett+ 2015) | 0.150 | $0.530\pm0.160$ | 0.4574 | -0.45 | 0.4083 | $8.30\times 10^{-41}$ | $5.02\times 10^{-41}$ |
| BOSS DR12 (Alam+ 2017 z1) | 0.380 | $0.497\pm0.045$ | 0.4750 | -0.49 | 0.4411 | $1.01\times 10^{-40}$ | $7.23\times 10^{-41}$ |
| WiggleZ (Blake+ 2012 z1) | 0.440 | $0.413\pm0.080$ | 0.4750 | +0.77 | 0.4451 | $1.06\times 10^{-40}$ | $7.87\times 10^{-41}$ |
| BOSS DR12 (Alam+ 2017 z2) | 0.510 | $0.459\pm0.038$ | 0.4732 | +0.37 | 0.4478 | $1.11\times 10^{-40}$ | $8.66\times 10^{-41}$ |
| WiggleZ (Blake+ 2012 z2) | 0.600 | $0.390\pm0.063$ | 0.4686 | +1.25 | 0.4487 | $1.17\times 10^{-40}$ | $9.72\times 10^{-41}$ |
| BOSS DR12 (Alam+ 2017 z3) | 0.610 | $0.436\pm0.034$ | 0.4680 | +0.94 | 0.4486 | $1.17\times 10^{-40}$ | $9.84\times 10^{-41}$ |
| eBOSS DR16 LRG (Bautista+21) | 0.700 | $0.473\pm0.041$ | 0.4612 | -0.29 | 0.4467 | $1.23\times 10^{-40}$ | $1.10\times 10^{-40}$ |
| WiggleZ (Blake+ 2012 z3) | 0.730 | $0.437\pm0.072$ | 0.4586 | +0.30 | 0.4457 | $1.25\times 10^{-40}$ | $1.14\times 10^{-40}$ |
| eBOSS DR16 ELG (Tamone+ 20) | 0.850 | $0.315\pm0.095$ | 0.4471 | +1.39 | 0.4395 | $1.32\times 10^{-40}$ | $1.30\times 10^{-40}$ |
| FastSound (Okumura+ 2016) | 1.360 | $0.482\pm0.116$ | 0.3896 | -0.80 | 0.3960 | $1.60\times 10^{-40}$ | $2.11\times 10^{-40}$ |
| eBOSS DR16 QSO (Neveux+ 20) | 1.480 | $0.462\pm0.045$ | 0.3762 | -1.91 | 0.3844 | $1.67\times 10^{-40}$ | $2.34\times 10^{-40}$ |

---

## 6. Algebraic BBN order-of-magnitude benchmark

No nuclear reaction network (e.g. PArthENoPE) is run. The abundances are linearized standard values at $\eta_{10}=6.1245$ with an algebraic TAL shift; $H(T)$ is a simple extrapolation without $e^\pm$ thermodynamics or variable $g_*$, and the scalar source contains only dust matter, not the thermal-plasma trace. The pulls divide only by the quoted observational errors.

**Proxies and sources (re-audit R6).** The clock proxy $|\delta G/G|\approx2|\phi|$ and the coefficients $0.08$ ($\delta Y_p\approx0.08\,\delta H/H$) and $1.6$ ($\delta\ln({\rm D/H})\approx1.6\,\delta H/H$) are **assumed order-of-magnitude proxies**, not derived here. Frame and sign (review F01): with $g_J=A^2g_E$ the Jordan-frame Ricci term carries $A^{-2}/G_*$ ($A^{-2}$ multiplies the Planck mass squared), so the tensor coupling is $1/G_T=A^{-2}/G_*$, i.e. $G_T=G_*A^2=G_*e^{2\alpha_Y\phi}$, and $\phi<0$ **decreases** $G_T$. The tensor coupling $G_T$, the operationally measured $G_{\rm lab}$ (Section 1.5) and the bare $2|\phi|$ proxy are three different quantities; only the magnitude $2\alpha_Y|\phi|$ is used, and the sign and reference state are not propagated. For comparison, the Kneller–Steigman fits (New J. Phys. 6, 117, 2004) give $\delta Y_p\approx0.16\,(S-1)$ and $\delta\ln({\rm D/H})\approx1.6\,(S-1)$ with $S=H'/H$: the coefficient 1.6 matches, 0.08 is a factor 2 smaller (immaterial at $|\delta H/H|\lesssim10^{-25}$). Input values: $Y_p^{\rm std}=0.24709$ is the PRIMAT central value (Pitrou, Coc, Uzan & Vangioni 2018, Phys. Rep. 754, 1); $10^5\,{\rm D/H}^{\rm std}=2.5147$ is an **arithmetic proxy** (review P02): a linearized benchmark value without verified single-paper provenance, not a sourced BBN prediction (PRIMAT 2018 gives 2.459); $Y_p^{\rm obs}=0.245\pm0.003$ is the PDG 2024 BBN review recommendation (Fields, Molaro & Sarkar), consistent with Aver et al. 2021 (JCAP 03, 027: $0.2453\pm0.0034$); $10^5\,{\rm D/H}^{\rm obs}=2.527\pm0.030$ is Cooke, Pettini & Steidel 2018 (ApJ 855, 102) — the CSV field name `DH_e5_obs_PDG` is a legacy label (PDG 2024 recommends $2.547\pm0.025$).

| Epoch | $z$ | $H/(c/L)$ | $\rho_m$ [kg/m$^3$] | $\lvert\bar\phi_{\rm attr}\rvert$ | $\delta H/H$ (scalar energy) | $\delta G/G$ (clock) | $\lvert\delta Y_p\rvert$ |
|:------------------------------|------:|------:|------:|------:|------:|------:|------:|
| end of BBN, $T\approx 0.02$ MeV | $1.0\times 10^{8}$ | $6.992\times 10^{-9}$ | $2.686\times 10^{-3}$ | $1.995\times 10^{-22}$ | $1.705\times 10^{-27}$ | $3.990\times 10^{-22}$ | $1.596\times 10^{-23}$ |
| D bottleneck, $T\approx 0.09$ MeV | $4.0\times 10^{8}$ | $1.119\times 10^{-7}$ | $1.719\times 10^{-1}$ | $1.277\times 10^{-20}$ | $2.728\times 10^{-26}$ | $2.553\times 10^{-20}$ | $1.021\times 10^{-21}$ |
| $z=10^9$, $T\approx 0.23$ MeV | $1.0\times 10^{9}$ | $6.992\times 10^{-7}$ | $2.686\times 10^{0}$ | $1.995\times 10^{-19}$ | $1.705\times 10^{-25}$ | $3.990\times 10^{-19}$ | $1.596\times 10^{-20}$ |

- $Y_p=0.24709$ versus $0.245\pm0.003$ (pull $0.70$); $10^5\,{\rm D/H}=2.5147$ versus $2.527\pm0.030$ (pull $-0.41$).
- Maximum values over the three epochs: $H/(c/L)\le6.992\times 10^{-7}$, $|\bar\phi_{\rm attr}|\le1.995\times 10^{-19}$, $\delta H/H\le1.705\times 10^{-25}$, $|\delta Y_p|\le1.596\times 10^{-20}$. The qualitative expectation of a tiny adiabatic correction survives; this is not a general BBN bound.

---

## 7. Local Yukawa constant-density bulk benchmarks (`expm1`)

The rows are bulk products $|\phi|_{\rm bulk}=\alpha_Y\gamma_L\rho$ for constant density. The half-space edge value $|\phi|_{\rm bulk}/2$ for $R\gg L$ is given only where a boundary geometry is defined (mean solar density, mean Earth density, white dwarf interior; status `half-space boundary`). The homogeneous FLRW, cosmic void and Milky Way disk rows have no boundary geometry here: their edge cell is n.a. (CSV/JSON value null, `boundary_status` = `n.a. (no boundary geometry)`, edge $\Delta H_0$ null); their bulk value is neither copied into the edge column nor halved (review N01). The values are not solved surface fields of the Sun or the Earth (that needs the real source and density profile). Clock factors use `expm1` to avoid round-off zeros; the apparent shift is $\Delta H_0^{\rm app}=67.40\,{\rm expm1}(\alpha_Y|\phi|)$. The row "V21 top-hat interior external benchmark" ($|\phi|=0.0412$) is an **external input, not reproduced here**: an interior value of the V21 top-hat fixed-point profile (SLy `gradient`, $r=2.0046$ km, $x=0.041239$; APR4 `mass`, $r=2.7542$ km, $x=0.041233$; Section 1.7), not a surface value (edge column n.a.) and not a solution of the present action.

| Benchmark | $\rho$ [kg/m$^3$] | $\lvert\phi\rvert_{\rm bulk}=\alpha_Y\gamma_L\rho$ | $\lvert\phi\rvert$ ($R\gg L$ edge) | ${\rm expm1}(-\alpha_Y\lvert\phi\rvert)$ | $\Delta H_0^{\rm app}$ [km/s/Mpc] | fraction of 5.4 gap |
|:------------------------------|------:|------:|------:|------:|------:|------:|
| Homogeneous FLRW, $z=0$ | $2.686\times 10^{-27}$ | $1.995\times 10^{-46}$ | n.a. | $-1.995\times 10^{-46}$ | $1.345\times 10^{-44}$ | $2.490\times 10^{-45}$ |
| Cosmic void ($0.2\rho_{m,0}$) | $5.373\times 10^{-28}$ | $3.990\times 10^{-47}$ | n.a. | $-3.990\times 10^{-47}$ | $2.689\times 10^{-45}$ | $4.980\times 10^{-46}$ |
| Milky Way disk (average) | $1.000\times 10^{-20}$ | $7.426\times 10^{-40}$ | n.a. | $-7.426\times 10^{-40}$ | $5.005\times 10^{-38}$ | $9.269\times 10^{-39}$ |
| Mean solar density (bulk) | $1408.0000$ | $1.046\times 10^{-16}$ | $5.228\times 10^{-17}$ | $-1.046\times 10^{-16}$ | $7.047\times 10^{-15}$ | $1.305\times 10^{-15}$ |
| Mean Earth density (bulk) | $5514.0000$ | $4.095\times 10^{-16}$ | $2.047\times 10^{-16}$ | $-4.095\times 10^{-16}$ | $2.760\times 10^{-14}$ | $5.111\times 10^{-15}$ |
| White dwarf interior (bulk) | $1.000\times 10^{9}$ | $7.426\times 10^{-11}$ | $3.713\times 10^{-11}$ | $-7.426\times 10^{-11}$ | $5.005\times 10^{-9}$ | $9.269\times 10^{-10}$ |
| V21 top-hat interior external benchmark | n.a. | $0.0412$ | n.a. | $-0.0404$ | $2.8349$ | $0.5250$ |

- Mean solar density bulk: $\Delta H_0^{\rm app}=7.047\times 10^{-15}\;\text{km/s/Mpc}$ (not a mathematical zero).
- Only the external V21 interior benchmark reaches $\Delta H_0^{\rm app}=2.835\;\text{km/s/Mpc}$ (fraction $0.525$ of the 5.4 km/s/Mpc gap); this is a conditional illustration; identifying a local clock factor with the cosmological $H_0$ needs a separate derivation with a defined observer convention.

---

## 8. No-go appendix

### 8.1 M2 / M3 branches

- `M2` (flat, uncompensated negative vacuum density): $E^2(0)=-0.37020$, so the declared normalization to a positive $H_0$ today is unphysical. The root of $E^2(z)=0$ is $z\approx0.295737$ in the full background (matter-only approximation: $0.295878$). A negative vacuum energy is not excluded in general; this specific normalized point is.
- `M3` (open, $\Omega_K=1.3702$), with its **own** sound horizons re-integrated with the M3 $H(z)$ (review F03): $r_s(z_*)=144.36622618\;\text{Mpc}$, $r_d(z_{\rm drag})=147.01159289\;\text{Mpc}$, $D_M(z_*)=34916.4\;\text{Mpc}$, $100\theta_*=0.41346227$ (formal geometric pull $-2092.09$ on the 1D $\theta_*$ prior, not a likelihood significance) and BAO $\chi^2=2201.29927$: a very poor acoustic-geometry benchmark. (The superseded values with the frozen M4 horizons were $0.41368156$ and $2209.46156$.)

### 8.2 The $67.36\,(1+1/4\pi)$ coincidence

$67.36\,(1+1/4\pi)=72.7203\;\text{km/s/Mpc}$ and $67.36\,e^{1/4\pi}=72.9394\;\text{km/s/Mpc}$ are correct arithmetic identities, but no dynamical $H_0$ correction follows from them: the clock factor is set by the solved potential $|\phi|$ (Sections 1 and 7), not by the bare coupling $\beta_{\rm Yuk}$.

### 8.3 Inverse distance ladder with profiled statistics

| Scenario | $H_0$ | $\Omega_m$ | $\Omega_K$ | $r_d$ [Mpc] | $P=H_0r_d$ [km/s] | $\chi^2_{\rm BAO}$ only (12×12) | $\Delta r_d$ vs Planck [Mpc] |
|:------------------------------------|-----:|------:|------:|----------:|----------:|-------:|-------:|
| 1. Fixed Step 25 M4 baseline | 67.40 | 0.3148 | 0.0000 | 147.09 | 9914.02 | 20.785 | +0.00 |
| 2. M4 shape, $P$ profiled | 72.80 | 0.3148 | 0.0000 | 137.77 $\pm$ 0.65 | 10029.65 $\pm$ 47.43 | 14.701 | -9.32 |
| 3. $\Omega_m$ and $P$ profiled | 72.80 | 0.2939 | 0.0000 | 140.02 $\pm$ 1.74 | 10193.67 $\pm$ 126.47 | 12.741 | -7.07 |
| 4. Joint BAO + Planck $r_d$ + assumed local $H_0$ prior | 70.81 | 0.2719 | 0.0000 | 147.03 (prior $\sigma$ 0.26) | 10410.24 | 15.674 | -0.06 |
| 5. Rigid all-$z$ $H_0=72.8$, unchanged $r_d$ | 72.80 | 0.3148 | 0.0000 | 147.09 | 10708.32 | 194.342 | +0.00 |
| 6. Fixed-$P$ illustration | 72.80 | 0.3148 | 0.0000 | 136.18 | 9914.02 | 20.785 | -10.91 |
| 7. Flat 72.8, $\Omega_m=0.2448$, re-int. $r_d$ | 72.80 | 0.2448 | 0.0000 | 150.65 | 10967.41 | 70.304 | +3.56 |
| 8. Open 72.8, fixed $\omega_{cb}$ | 72.80 | 0.2698 | 0.0098 | 147.09 | 10708.27 | 61.206 | +0.00 |

- The $\chi^2_{\rm BAO}$ column contains only the DESI $12\times12$ term. Row 4 adds the priors: $\chi^2_{\rm joint}=\chi^2_{\rm BAO}+\chi^2_{r_d\,\rm prior}+\chi^2_{H_0\,\rm prior}=15.67431585+0.06234684+3.97660606=19.71326876$; its $r_d$ entry quotes the **input prior error** $\pm0.26$ Mpc, not a fitted uncertainty. Row 5 is the rigid all-$z$ rescaling with unchanged $r_d$, not a late-transition model (Section 3).
- The fixed-$P$ illustration (row 6: $P=9914.02\;\text{km/s}$, $r_d=136.18\;\text{Mpc}$) is a conditional arithmetic example; $P=9914.02$ is not the BAO best-fit product, and its earlier sigma label is withdrawn.
- The profiled BAO requirement at $H_0=72.8$ is $r_d=140.02\pm1.74\;\text{Mpc}$ (two-parameter control), and the joint BAO + Planck $r_d$ + assumed local $H_0$ benchmark prior control gives $\Delta\chi^2=6.973$.
- The rigid rescaling at unchanged shape and $r_d$ is strongly disfavoured by BAO ($\chi^2=194.342$, $p=4.79\times 10^{-35}$).

![Inverse distance ladder in the $(H_0,r_d)$ plane (`step2_3_inverse_ladder_rd_H0_exclusion_en.png`)](step2_3_inverse_ladder_rd_H0_exclusion_en.png){ width=96% }

---

## 9. What remains and what is not established

**Remains:**

- A rigid 8% increase of $H(z)$ at all redshifts with fixed $E(z)$ shape and unchanged $r_d$ fits BAO very poorly.
- The bare $\beta_{\rm Yuk}$ is not equal to a solved scalar potential or a clock calibration shift.
- In the short-range adiabatic branch with negligible oscillation energy the late-time cosmological scalar correction is extremely small.
- The specific M2 normalization is unphysical and the M3 curvature is a poor acoustic-geometry benchmark.
- The profiled constant SNe zero-point offset and the $H_0$–$M_B$ degeneracy are correct geometric illustrations.

**Not established in this package:**

- a joint Planck + Pantheon+/Union3 + RSD + BBN data validation;
- general physical exclusion strengths of the kind formerly suggested by sigma labels;
- a cosmological closure of the scalar–tensor theory for arbitrary allowed initial conditions;
- a new TAL-specific empirical advantage over a properly fitted $\Lambda$CDM.

**Closing statement.** In the examined M4 branch with $L=10$ km, fixed coupling and negligible free scalar oscillation energy, the local Yukawa factor provides no derived mechanism for a global Hubble constant of 72.8 km/s/Mpc. A simple global rescaling of a fixed early background and a fixed late-time shape is strongly disfavoured by BAO. The result is not a comprehensive joint Planck–SNe–RSD–BBN likelihood analysis and not a no-go theorem for all cosmological branches of the theory.

### 9.1 Nyitott korlátok / Open limitations

| Item | Open limitation | Status in this package |
|:-----|:-----------------------------|:-----------------------------|
| P01 | Verified primary source of the $z=0.845$ ELG1 point $18.33\pm0.62$ | chosen illustration without a verified primary source; not a DESI measurement; not in the 12-point likelihood |
| P02 | Single-paper source of $10^5\,({\rm D/H})_{\rm std}=2.5147$ | arithmetic proxy; an empirical standard-model test needs a verified input and theory error |
| P03 | Coupled scalar–Friedmann–matter initial-value problem and free homogeneous amplitudes | branch $C=0$ is selected explicitly; no statement about all initial conditions or all cosmological TAL branches |
| P04 | Planck TT/TE/EE + lensing likelihood and a recombination calculation | 3D prior assumed, $z_*$ and $z_{\rm drag}$ fixed control parameters; geometric benchmark, not a CAMB/CLASS validation |
| P05 | Real Pantheon+/Union3 data vector and covariance | synthetic 18-node metric, not an observational model preference |
| P06 | RSD covariance, overlaps, AP mapping, derived amplitude, $k$ spectrum | $k=0.1\,h$/Mpc variation with fixed $\sigma_{8,0}$; $1.666\times10^{-40}$ refers to this ODE, normalization and initial conditions |
| P07 | BBN thermal trace, variable $g_*$, $e^\pm$ thermodynamics, reaction network | three sample epochs with assumed sensitivities; not an abundance likelihood |
| P08 | Real solar/terrestrial density profiles, surface Yukawa solution, link of local clocks to the distance ladder | bulk and half-space limits; the $\Delta H_0^{\rm app}$ table is an arithmetic illustration |
| P09 | The V26 book as a whole, all historical conventions and all stellar BVPs | this review concerns Step 26; the V26 normalization identity is a targeted check only |

---

## 10. Files and reproducibility

- `desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`: 376 B, SHA-256 `dd2873a0b88459a491af3c0c0307ba059f62df9211d5b976760f310565a1be68`
- `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`: 2,170 B, SHA-256 `bbafa9074b51cf1a45e0d10e4f37db8c0e80a5d1d1788857abb7fc49fb21abcc`
- `tal_cosmology_common.py`: 24,847 B, SHA-256 `54c5636b35631af314f501060e39adf723a46cada78db15526b0591a87c2c109`
- `tal_hubble_desi_bao_summary.json`: 52,790 B, SHA-256 `94ba04e4c09fc8a7d6a94eed04eecf03c9f794444a3a9f646a58f80aa1e3ca53`

- Generators: `tal_cosmology_common.py` (shared background), `run_step1_1.py`, `run_step1_2_and_1_3.py`, `run_step1_4_to_1_6.py`, `run_step2_nogo_and_report.py`, `run_step26_revision_controls.py` (re-audit controls), `hubble_desi_bao_tal_audit.py`; text: `generate_step26_markdown.py` + `step26_md_templates.py`, text regression test `test_step26_text_regression.py`; PDF: `compile_step26_clean_pdfs.py`; runner: `run_all.py` (order and data hash check, `--skip-block1`), dependencies: `requirements.txt`.
- The package is self-contained: both DESI input files, all CSV and PNG outputs, the strict JSON (`allow_nan=False`, NaN written as `null`) and `SHA256SUMS.txt` are included.
- External sources: DESI VI (arXiv:2404.03002), DESI III (arXiv:2404.03000), Planck 2018 VI (arXiv:1807.06209), Cobaya `bao_data` repository (DESI 2024 Gaussian BAO mean and covariance).

---

## 11. Revision log

| Item | Change | Where |
|:-----|:-------|:------|
| R1 | Kinetic normalization reconciled with V26 (file:line), $\phi_{\rm old}=2\phi_{\rm new}$ needs $\alpha_{\rm new}=2\alpha_{\rm old}$; carried-over results listed; compact-star $\lvert\phi\rvert$ marked as external V21 interior value; $2\alpha^2\beta$ relabelled sensitivity control | 1.5, 1.7, 7 |
| R2 | Piecewise $H=sH_{\rm M4}$ ($z<z_t$) control for $z_t\in\{2.5,10,100,1000\}$; $D_M/s$ row reclassified as inconsistent geometric illustration | 3, `step26_r2_piecewise_transition.csv` |
| R3 | Self-contained package: DESI inputs with hashes, all CSV/PNG, `requirements.txt`, `run_all.py`, strict JSON; old `build_step26_pdfs.py` moved out of the package | 10, README |
| R4 | Joint $\chi^2$ split into BAO, $r_d$-prior and $H_0$-prior terms; table column is BAO-only; $r_d$ $\pm0.26$ labelled input prior error | 2, 8.3 |
| R5 | Symmetric zero-point offset for M4 and the low-$\Omega_m$ model; lowest-node value separated from the $5\log_{10}$ identity | 4, `step26_r5_sne_symmetric_offsets.csv` |
| R6 | BBN clock proxy and coefficients labelled assumed proxies; frame and sign discussed; sources given | 6 |
| 5.1 | Float64 branch uses $1+\Delta G_{\rm eff,can}/G_N$; absolute-value headers; independent rescaled check quoted | 5 |
| 5.2 | $z=0$ drift separated from the Fermi–Dirac interpolation error, which is computed | 1.6 |
| 5.3 | $\theta_*$-matched model relabelled to $\Omega_m=0.24419$, $\omega_{cb}=0.12877$; inverse-ladder figure uses $\sqrt{\Delta\chi^2}$ wording; provenance statement reworded | 3, 8.3, header |


---

## 12. Revíziós napló 2 / Revision log 2 (second review, F01–F12)

| Item | Change | Targeted check |
|:-----|:-----------------------------|:-----------------------------|
| F01 | $G_J=GA^{-2}$ replaced by $1/G_T=A^{-2}/G_*$, $G_T=G_*A^2$; $\phi<0$ decreases $G_T$; $G_T$, $G_{\rm lab}$ and the bare $2\lvert\phi\rvert$ proxy separated (Section 6, code comments, R6 metadata) | text and JSON `R6_bbn_proxy_provenance.clock_proxy` contain the new relation |
| F02 | G dictionary after Section 1.5 ($G_*$, $G_T$, $G_{\rm lab}\approx G_*(1+q)$, $G_N:=G_*$) | $1/(1+q)=0.9262883$ computed in `run_step26_revision_controls.py` |
| F03 | M3 uses its own $r_s$, $r_d$ | $r_s=144.36622618$, $r_d=147.01159289$ Mpc, $100\theta_*=0.41346227$, $\chi^2_{\rm BAO}=2201.29927$ (assertion against the review values) |
| F04 | step2 status strings and CLI text rewritten: "BAO fit", "deviation on illustrative prior", "pull vs fixed target"; no Planck exclusion or curvature-likelihood wording | automatic string scan in `run_step26_revision_controls.py` |
| F05 | $72.8\pm1.0$ labelled assumed local $H_0$ benchmark prior; JSON key `H0_local_benchmark_prior_km_s_Mpc`; figure labels changed | labels no longer attribute the value to SH0ES (CSV/JSON/PNG text) |
| F06 | Section 1.5 sentence on the $2\alpha_Y^2\beta_{\rm Yuk}$ row | text check |
| F07 | RSD legend shows the canonical maximum | $1.6664\times 10^{-40}$ in the legend and JSON |
| F08 | lowest-node sentence removed from Section 4, identity kept | text check |
| F09 | SNe figure: low-$\Omega_m$ offset curve added, "profiled constant zero-point offset", $N=18$, $\nu=17$ | low-$\Omega_m$ with offset $\chi^2=22.777527$ |
| F10 | row renamed "V21 top-hat interior external benchmark", edge column n.a., provenance added | SLy $r=2.0046409680935633$ km, $x=0.041239199439254325$; APR4 $r=2.75417748649486$ km, $x=0.041233030079658145$ |
| F11 | PDF builder and `run_all.py --pdf` exit nonzero without pandoc/pdflatex or without a newly written PDF; `--tex-only` has its own status | test with a `PATH` lacking pdflatex |
| F12 | `theta_star_100_rs_only_shifted_fixed_DM` (renamed) and `theta_star_100_both_endpoints_shifted` (new) | low-$\Omega_m$: $1.04248582$ and $1.04250287$ |

## 13. Revíziós napló 3 / Revision log 3 (third re-audit, N01, N02, T01)

| Item | Change | Targeted check |
|:-----|:-----------------------------|:-----------------------------|
| N01 | Section 7 table: the "$R\gg L$ edge" cell of the homogeneous FLRW, cosmic void and Milky Way disk rows is n.a. (no boundary geometry); `step1_6_local_yukawa_environments.csv` and the JSON carry null in `abs_phi_surface_R_gg_L` and `Delta_H0_apparent_surf_km_s_Mpc` and the new per-row field `boundary_status`; bulk values unchanged (not halved) | assertion in `run_step26_revision_controls.py` (JSON `review3_N01_N02_T01`); text regression test |
| N02 | Section 1.7 "Carried over unchanged": the parenthesis now separates the background scalar-energy correction $\Omega_\phi\sim10^{-47}$ from the direct Yukawa perturbation of RSD ($\delta G_{\rm eff}/G_*\approx3.8\times 10^{-41}$, $\delta f\sigma_8\approx1.67\times 10^{-40}$, Section 5); numbers taken from the JSON | text regression test (HU and EN sentence) |
| T01 | New `test_step26_text_regression.py`, run by `run_all.py` after the Markdown generation, nonzero exit on failure: checks the HU/EN F06 sentence, the F08 identity wording, the absence of the superseded wrong sentences (V26 factor-2 growth-law attribution; the identity value attributed to the lowest table node), the N01 table cells, the N02 sentence and the canonical legend string of the RSD figure source (F07); the comment in `run_step26_revision_controls.py` updated | mutation test on a scratch copy: reinserting a wrong sentence gives a nonzero exit |
