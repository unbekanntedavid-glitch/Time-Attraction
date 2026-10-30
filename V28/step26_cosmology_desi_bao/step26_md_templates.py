# -*- coding: utf-8 -*-
"""Text templates (EN/HU/README) for generate_step26_markdown.py.
Placeholders <<key>> are filled from tal_hubble_desi_bao_summary.json and the Step 26 CSVs (computed results)."""

TEMPLATE_EN = r"""# Step 26 (V27): Cosmological $H(z)$, DESI DR1 BAO, Planck 2018 Acoustic Scale, Synthetic SNe Ia Shape Control, Variational RSD $f\sigma_8$, Algebraic BBN Benchmark and Inverse Distance Ladder No-Go Audit

**Author:** Dávid Szabolcsi (`ORCID: 0009-0005-6300-074X`)  
**Date:** 2026-10-28 (text revised 2026-10-29 after the Academia Co-scientist review, the subsequent re-audit and the second review F01–F12 and the third re-audit N01, N02, T01, see Sections 11, 12 and 13)  
**Preceding version DOIs:**

- **V24:** Zenodo DOI `10.5281/zenodo.23680715` (strictly V24)
- **V25:** Zenodo DOI `10.5281/zenodo.23699685` (strictly V25, `step1`–`step22`)
- **V26:** Zenodo DOI `10.5281/zenodo.23727585` (strictly V26, `step1`–`step25`)

**Source of the numbers:** computed results are generated from the JSON/CSV sources (`tal_hubble_desi_bao_summary.json` and the Step 26 CSV tables, produced with the shared background module `tal_cosmology_common.py`) by `generate_step26_markdown.py`. Fixed input constants (data values, priors, literature benchmarks, external inputs) are stated in the code and cited in the text.

---

## 0. Scope and summary

This module tests one specific branch: the V26 `Step 25` M4 background at fixed $L=<<L_km>>\;\text{km}$, $\alpha_Y=1$, $\beta_{\rm Yuk}=1/(4\pi)=<<beta>>$, on the adiabatic particular solution of the scalar equation with negligible free scalar oscillation energy. The narrower, defensible statement is: **in this branch the bare coupling $\beta_{\rm Yuk}$ inserted as a multiplier of $H_0$ is not a derived mechanism for the shift $67.4\to 72.8\;\text{km/s/Mpc}$.** This is not an exclusion of every possible cosmological solution of TAL.

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

The scalar Compton rate is $c/L=<<mphic>>\;\text{s}^{-1}$ and $H_0L/c=<<H0overm>>$.

### 1.3 Attractor and branch selection

Because $4\pi\beta_{\rm Yuk}=1$, the adiabatic particular solution ($\ddot{\bar\phi},\,3H\dot{\bar\phi}\ll (c/L)^2\bar\phi$) is

$$
\bar\phi_{\rm attr}(z)=-\alpha_Y\gamma_L\left(\rho_m(z)-\frac{3p_m(z)}{c^2}\right),\qquad \gamma_L=\frac{4\pi\beta_{\rm Yuk}GL^2}{c^2}=\frac{GL^2}{c^2}=<<gammaL>>\;\text{m}^3/\text{kg},
$$

so that $\bar\phi_{\rm attr}(0)=<<phi0>>$ and $\bar\phi_{\rm attr}(z_*)\approx<<phistar>>$. The general homogeneous solution also contains a free oscillation

$$
\delta\phi\simeq C\,a^{-3/2}\cos\!\left[(c/L)\,t+\vartheta\right],\qquad \langle\rho_{\delta\phi}\rangle\propto a^{-3}.
$$

$H/(c/L)\ll 1$ justifies the adiabatic approximation but does not remove this mode. **This analysis therefore selects the branch $C=0$ ($\Omega_{\delta\phi,0}=0$): the adiabatic particular solution with negligible free oscillation energy.** No initial-value integration of the scalar equation is performed, and the conclusions below hold only for this branch.

### 1.4 Scalar energy density

From the action, $\rho_\phi=\dfrac{c^2}{8\pi G\beta_{\rm Yuk}}\left[(\dot\phi/c)^2+(\phi/L)^2\right]$ (mass density), hence on the attractor

$$
\Omega_{\phi,\rm attr}(0)=\frac{(c/L)^2\,\bar\phi_{\rm attr}(0)^2}{3\beta_{\rm Yuk}H_0^2}=<<Omphi0>>.
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
| $G_{\rm lab}$ | operationally measured $G$: for $r\ll L$ and unscreened small bodies $\Phi=-G_*M(1+q\,e^{-r/L})/r$, so $G_{\rm lab}\approx G_*(1+q)$, $q=\alpha_Y^2\beta_{\rm Yuk}=<<qval>>$ |
| $G_N$ | the symbol in the perturbation equation of Section 5: $G_N:=G_*$ |

On scales $kL/a\ll1$ the growth equation therefore uses $G_{\rm eff}/G_{\rm lab}=1/(1+q)=<<Gratio>>$ if the coupling is calibrated by a laboratory value; the $10^{-41}$-level statement of Section 5 is a formal control at the same bare $G_*$. No $f\sigma_8$ exclusion claim is made and no screening branch is introduced.

### 1.6 Background normalization and observer frame

All scripts use one background module (`tal_cosmology_common.py`). The smooth massive-neutrino interpolation alone gives $E^2(0)=<<E2unnorm>>$ (offset $<<E2offset>>$); the dark-energy closure absorbs this offset so that $E^2(0)=1$ exactly and $H(0)=H_0$. This offset ($<<NUdrift>>$ in $\Omega_\nu$) is the $z=0$ normalization drift of the interpolation, not its approximation error. Against the Fermi–Dirac energy integral of the massive species (matched to the same relativistic and non-relativistic limits), the $p=1.83$ interpolation has a relative error of $<<NUrel1000>>$ in the massive-neutrino density and $<<NUH1000>>$ in the total $H$ at $z=1000$; the largest absolute values on the grid $0\le z\le10^4$ are $<<NUmax>>$ (massive component, $z\approx<<NUmaxZ>>$) and $<<NUHmax>>$ (total $H$, $z\approx<<NUHmaxZ>>$). This is an approximation error of the $\Lambda$CDM-like background, common to all compared models, and separate from the theoretical TAL correction ($\sim 10^{-46}$).

<<TABLE_NU>>

 The computations use $\Lambda$CDM-like $H(z)$ functions with a separate clock factor; the Jordan-frame mapping ($dt_J=A\,dt_E$, $a_J=A\,a_E$, $H_J=(H_E+\alpha_Y\dot\phi_E)/A$) is a leading-order identity at this $|\phi|$, not a numerically solved coupled Friedmann–scalar–matter system.

### 1.7 Version reconciliation of the kinetic normalization (re-audit R1)

Paths below are relative to the V26 tree `scalar-tensor-yukawa-fit/` (V26, DOI `10.5281/zenodo.23727585`).

- The first Step 26 draft wrote the kinetic term as $\frac{c^4}{16\pi G}\frac{1}{2\beta_{\rm Yuk}}(\partial\phi)^2=\frac{c^4}{32\pi G\beta_{\rm Yuk}}(\partial\phi)^2$, while the action of Section 1.1 has $\frac{c^4}{8\pi G\beta_{\rm Yuk}}(\partial\phi)^2$: $K_{\rm old}/K_{\rm new}=<<R1Kratio>>$.
- V26 uses $\frac{c^4}{16\pi G\tilde\alpha_L}(\partial\phi_L)^2$ with $\tilde\alpha_L=\alpha_Y/(8\pi)$, i.e. $c^4/(2G\alpha_Y)$ (`step23_gw170817_and_lensing_en.md:19` and `:29`, `README_en.md:50`, `gw170817_and_lensing_audit.py:104–113`). At $\alpha_Y=1$, $\beta_{\rm Yuk}=1/(4\pi)$ this equals the new coefficient ($K_{\rm V26}/K_{\rm new}=<<R1Kv26>>$), and the V26 fifth-force strength $2\tilde\alpha_L\alpha_Y^2=<<R1betaV26>>$ equals $\alpha_Y^2\beta_{\rm Yuk}$. The V26 scalar equation $(\Box-L^{-2})\phi=-(\alpha_Y^2G/c^4)\varepsilon_m$ (`step25_cmb_and_cosmology_en.md:37`) and $\phi_{\rm attr}=\alpha_Y^2\gamma_L\rho$ (`cmb_and_cosmology_audit.py:327`, `:620`) agree with Section 1.3 up to the sign convention $\phi_{\rm V26}=-\phi_{\rm V27}$. **Conclusion:** V26 used the new normalization; the $c^4/(32\pi G\beta)$ form was a transcription inconsistency of the first draft, and the computed numbers already followed the new normalization. For $\alpha_Y\neq1$ the two versions differ, since V26 ties $\beta_{\rm Yuk}=\alpha_Y^3/(4\pi)$ whereas here $\beta_{\rm Yuk}$ is an independent parameter.
- The field redefinition $\phi_{\rm old}=2\phi_{\rm new}$ maps the old kinetic term onto the new one, but $A=e^{\alpha_{\rm old}\phi_{\rm old}}=e^{2\alpha_{\rm old}\phi_{\rm new}}$, so it requires $\alpha_{\rm new}=2\alpha_{\rm old}$. With the same $\alpha_Y=1$ the old form would give the fifth-force strength $4\alpha_Y^2\beta_{\rm Yuk}=<<R1aold>>$ instead of $\alpha_Y^2\beta_{\rm Yuk}=<<R1anew>>$.
- Carried over unchanged: all background, BAO, CMB, SNe-shape and RSD numbers (the background scalar-energy correction is $\Omega_\phi\sim10^{-47}$; the direct Yukawa perturbation of RSD is a separate negligible term quantified in Section 5: $\delta G_{\rm eff}/G_*\approx<<N02dG>>$, $\delta f\sigma_8\approx<<N02dfs>>$), $\bar\phi_{\rm attr}(0)$, $\Omega_{\phi,\rm attr}(0)$ and the bulk local benchmarks at $\alpha_Y=1$. The V26 background additionally used a $(\phi/4\pi)^2$ energy convention (`cmb_and_cosmology_audit.py:330`, `:623`), giving $\Omega_\phi=<<R1Omv26>>$ instead of the action-consistent $<<Omphi0>>$; both are negligible.
- The factor-2 growth convention $2\alpha_Y^2\beta_{\rm Yuk}$ is relabelled a **sensitivity control**, not an equivalent convention.
- The V21 top-hat interior external benchmark $|\phi|=<<ns_phi>>$ (Section 7) is an interior value of the V21 top-hat fixed-point profile $x=L^2(G/c^2)\bar\rho$ (`stellar_tal_v21.py:140–160`, `data/v21/profiles_L10km.csv`: SLy, mode `gradient`, $r=2.0046409680935633$ km, $x=0.041239199439254325$; APR4, mode `mass`, $r=2.75417748649486$ km, $x=0.041233030079658145$), not a surface value and **not automatically a solution of the new Helmholtz–Yukawa action**; it remains an external input.

---

## 2. DESI DR1 BAO: $12\times 12$ covariance and profiled controls

The 12 primary observables are the Cobaya DESI DR1 vector and covariance (`desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`, `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`). The separate LRG3 and ELG1 single-tracer rows are diagnostics only and are not added to the likelihood; the $z=0.845$ ELG1 row ($D_V/r_d=18.33\pm0.62$) is a **chosen illustration without a verified primary source** and not a DESI measurement (it is not in DESI III Table 18 or DESI VI Table 1). Per-bin $\chi^2$ for primary bins use the $12\times 12$ covariance; diagnostic rows (italic) use their rounded table errors. Predictions in the table are $D_M/r_d$; $D_H/r_d$ for anisotropic bins and $D_V/r_d$ for isotropic bins.

<<TABLE_BAO>>

- Fixed M4 point ($H_0=67.40$, $r_d=<<rd>>\;\text{Mpc}$): $\chi^2=<<chi2M4>>$ for 12 observables, $p=<<pM4>>$ (rounded Table 1 input: $<<chi2M4t>>$).
- Rigid $H_0=72.8$ with unchanged shape and $r_d$: $\chi^2=<<chi2R>>$, $p=<<pR>>$ (rounded Table 1 input: $<<chi2Rt>>$).

**Profiled controls** (local Fisher errors; $p$-values are conditional frequentist values for the stated Gaussian likelihood, not probabilities of TAL being true):

<<TABLE_PROF>>

- At fixed $E(z)$, BAO constrains only $P=H_0r_d$ (rank-1 information in the $(H_0,r_d)$ plane).
- Two-parameter control: $\Omega_m=<<Om2>>\pm<<Om2s>>$, $P=<<P2>>\pm<<P2s>>\;\text{km/s}$, $\chi^2_{\min}=<<chi2P2>>$ (10 dof, $p=<<pP2>>$).
- $r_d$ required at $H_0=72.8$ exactly: $<<rdP2>>\pm<<rdP2s>>\;\text{Mpc}$; propagating $H_0=72.8\pm1.0$: $\pm<<rdP2sH>>\;\text{Mpc}$. Here $72.8\pm1.0$ is an **assumed local $H_0$ benchmark prior** (a chosen target value, review F05); it is not the Riess et al. 2022 value $73.04\pm1.04$ and no specific SH0ES release is claimed for it.
- With the model-dependent Planck prior $r_d=147.09\pm0.26\;\text{Mpc}$: $H_0=<<H0Pl>>\pm<<H0Pls>>\;\text{km/s/Mpc}$.
- Joint illustrative test (BAO + Planck $r_d$ prior + local $H_0=72.8\pm1.0$, profiling $\Omega_m,H_0,r_d$): best fit $\Omega_m=<<OmJ>>$, $H_0=<<H0J>>$, $r_d=<<rdJ>>\;\text{Mpc}$, $\chi^2_{\rm joint}=<<chi2J>>$, split as $\chi^2_{\rm BAO}=<<R4bao>>$, $\chi^2_{r_d\,\rm prior}=<<R4rd>>$, $\chi^2_{H_0\,\rm prior}=<<R4h0>>$; $\Delta\chi^2=<<dchi2J>>$ relative to the BAO-only minimum ($\sqrt{\Delta\chi^2}=<<sigJ>>$). This is not a Planck+TAL exclusion significance.

![DESI DR1 BAO, Step 25 M4 and profiled controls (`step1_1_desi_dr1_bao_plot_en.png`)](step1_1_desi_dr1_bao_plot_en.png){ width=96% }

---

## 3. Planck 2018 acoustic scale with model-specific sound horizons

For every model, $r_s(z_*=1089.80)$ and $r_d(z_{\rm drag}=1059.94)$ are re-integrated with that model's own $H(z)$ and $c_s(z)$ (fixed $z_*$, $z_{\rm drag}$ as a control approximation). The 3D Gaussian prior on $(\omega_b,\omega_{cb},100\theta_*)$ is an **illustrative assumed prior** (Planck 2018 marginal means and errors with an approximate correlation matrix whose chain provenance is not documented here); it is not the Planck TT,TE,EE+lowE+lensing likelihood. Reference $100\theta_*=<<th100ref>>\pm<<th100sig>>$.

<<TABLE_CMB>>

- Step 25 M4: $r_s=<<rs>>\;\text{Mpc}$, $r_d=<<rd>>\;\text{Mpc}$, $100\theta_*=<<th100M4>>$ (pull $<<pullthM4>>$), $\chi^2_{\rm 3D}=<<chi23DM4>>$.
- All-$z$ rigid rescaling $H\to sH$ with unchanged $c_s$: $r_s\to r_s/s=<<rsRig>>\;\text{Mpc}$, $r_d=<<rdRig>>\;\text{Mpc}$ and $D_M\to D_M/s$, so $100\theta_*=<<thRig>>$ is unchanged. Scaling $D_M$ alone would give $<<thRigFix>>$; that inconsistent combination is no longer used as a test.
- Flat $H_0=72.8$, $\Omega_m=0.244828$: re-integration gives $r_s=<<rsLow>>$, $r_d=<<rdLow>>\;\text{Mpc}$ and $100\theta_*=<<thLow>>$ (with the fixed early $r_s$ it would be $<<thLowFix>>$). This does **not** mean the model passes Planck: $\omega_{cb}=<<wcbLow>>$ is offset from the prior and the BAO $\chi^2$ with the re-integrated $r_d$ is $<<chi2Low>>$ (with the old fixed $r_d$: $<<chi2LowFix>>$). Early and late parts must be refitted together.
- **Two $z_*$-shift diagnostics (review F12).** With the Hu–Sugiyama fit shift $z_*\to z_*+\Delta z_*$ ($z_*+\Delta z_*=<<zsLowShift>>$ for this model), the CSV keeps two separately named values: `theta_star_100_rs_only_shifted_fixed_DM` $=100\,r_s(z_*+\Delta z_*)/D_M(1089.80)=<<thLowRsOnly>>$ (only the $r_s$ endpoint shifted) and `theta_star_100_both_endpoints_shifted` $=100\,r_s(z_*+\Delta z_*)/D_M(z_*+\Delta z_*)=<<thLowBoth>>$. Neither is a recombination solution.
- **Late transition (re-audit R2).** The earlier "late-time-only" row combined $D_M(z_*)/s$ with unchanged early $r_s,r_d$ ($100\theta_*=<<R2illus>>$); it has no transition redshift and is **reclassified as an inconsistent geometric illustration** (italic row in the table), removed from the consistent-model statements. A defined control is the piecewise $H(z)=s\,H_{\rm M4}(z)$ for $z<z_t$ and $H_{\rm M4}(z)$ for $z\ge z_t$, $s=72.8/67.4=<<R2s>>$, with $D_M(z_*)$ integrated numerically and $r_s=<<R2rs>>\;\text{Mpc}$ unchanged (M4: $D_M(z_*)=<<R2DMM4>>\;\text{Mpc}$). $H$ jumps at $z_t$, so this is a kinematic control, not a dynamical model; the pull column is a formal 1D distance from the Planck value, not a likelihood significance. For $z_t\ge2.5$ all DESI bins lie in the rescaled region, so the BAO $\chi^2$ equals the rigid value $<<chi2R>>$.

<<TABLE_R2>>

---

## 4. Synthetic distance-modulus shape control (not a Pantheon+ / Union3 likelihood)

This section uses **no** observed SN catalogue: <<nn>> redshift nodes ($z=<<zmin>>$–$<<zmax>>$) with representative errors $\sigma_\mu$ and a synthetic reference curve from a flat model with $H_0=72.8$, $\Omega_m=0.334$. There are no SN identifiers, observed $\mu$ values or statistical+systematic covariance. The $\chi^2$ values below are model-versus-model shape metrics, not observational goodness of fit.

- For identical $E(z)$ geometry the offset would be the constant identity $5\log_{10}(72.8/67.4)=+<<dmu>>\;\text{mag}$ (positive sign: M4 minus 72.8 reference).
- The tabulated lowest-node difference $<<R5node>>\;\text{mag}$ at $z=<<zmin>>$ is a model difference (the reference has $\Omega_m=0.334$), not this identity.
- **Symmetric zero-point treatment (re-audit R5):** both models receive the same profiled constant offset (absolute-magnitude marginalization), weights $1/\sigma_\mu^2$:

<<TABLE_R5>>

- The large drop for M4 only reflects that its shape is close to the reference after removing the constant; node positions and representative errors are synthetic, so changes of nodes or errors do not imply a physical improvement.
- The figure shows both models with and without the profiled constant zero-point offset (low-$\Omega_m$ with offset: $\chi^2=<<chi2LowSh>>$); $N=18$ nodes, $\nu=17$ degrees of freedom after one profiled offset. The offset is a profiled nuisance constant, not an observational calibration (review F09).
- Uncalibrated SNe alone do not measure $H_0$ ($H_0$–$M_B$ degeneracy). A real catalogue statement requires reading the Pantheon+ or Union3 likelihood and covariance, which this package does not do.

![Synthetic SNe Ia distance-modulus shape control (`step1_3_pantheon_sne_plot_en.png`)](step1_3_pantheon_sne_plot_en.png){ width=96% }

<<TABLE_SNE>>

---

## 5. RSD $f\sigma_8(z)$: diagonal benchmark and variational TAL correction

The growth equation $\delta_m''+\left(2+\frac{d\ln H}{d\ln a}\right)\delta_m'-\frac{3}{2}\Omega_{cb}(a)\frac{G_{\rm eff}(k,a)}{G_N}\delta_m=0$ is integrated with DOP853 (`rtol=1e-12`, `atol=1e-14`, dense output); $\sigma_{8,0}=0.8111$ is a fixed Planck input, not a TAL-derived amplitude.

- At $k=<<kcom>>\,h\,\text{Mpc}^{-1}$: $kL=<<kL>>$ ($(kL)^2=<<kL2>>$) at $z=0$; $\Delta G_{\rm eff}/G_N=<<dGcan>>$ (canonical $\alpha_Y^2\beta_{\rm Yuk}$) or $<<dG2x>>$ (factor-2 convention).
- In float64, $1+O(10^{-41})$ equals 1, so the direct difference $f\sigma_8^{\rm TAL}-f\sigma_8^{\rm GR}$ is $0$ by floating-point identity, not by a 40-digit measurement. The first-order variational equation gives $\max\lvert\delta f\sigma_8\rvert=<<dfsCan>>$ (canonical) and $<<dfs2x>>$ (factor-2 sensitivity control) on the 12 points; an independent rescaled check in the re-audit confirmed $1.666\times10^{-40}$ (maximum $1.6663681521\times10^{-40}$, relative deviation $\le5.1\times10^{-6}$). The float64 branch uses $G_{\rm eff}/G_N=1+\Delta G_{\rm eff,can}/G_N$ (re-audit 5.1). The figure legend shows the canonical maximum $<<dfsCanF>>$; the factor-2 value is only the separate sensitivity control (review F07). Here $G_N:=G_*$ (G dictionary, Section 1.5).
- Diagonal $\chi^2$ (no bin covariances, overlaps or AP treatment): M4 with clustering source $\Omega_{cb}$: $<<chi2RSD>>$; with all matter: $<<chi2RSDm>>$; flat $H_0=72.8$, $\Omega_m=0.2448$: $<<chi2RSDlow>>$ (12 points). This is a benchmark, not a covariance likelihood.

![RSD $f\sigma_8(z)$ diagonal benchmark (`step1_4_rsd_fsigma8_plot_en.png`)](step1_4_rsd_fsigma8_plot_en.png){ width=96% }

<<TABLE_RSD>>

---

## 6. Algebraic BBN order-of-magnitude benchmark

No nuclear reaction network (e.g. PArthENoPE) is run. The abundances are linearized standard values at $\eta_{10}=<<eta10>>$ with an algebraic TAL shift; $H(T)$ is a simple extrapolation without $e^\pm$ thermodynamics or variable $g_*$, and the scalar source contains only dust matter, not the thermal-plasma trace. The pulls divide only by the quoted observational errors.

**Proxies and sources (re-audit R6).** The clock proxy $|\delta G/G|\approx2|\phi|$ and the coefficients $0.08$ ($\delta Y_p\approx0.08\,\delta H/H$) and $1.6$ ($\delta\ln({\rm D/H})\approx1.6\,\delta H/H$) are **assumed order-of-magnitude proxies**, not derived here. Frame and sign (review F01): with $g_J=A^2g_E$ the Jordan-frame Ricci term carries $A^{-2}/G_*$ ($A^{-2}$ multiplies the Planck mass squared), so the tensor coupling is $1/G_T=A^{-2}/G_*$, i.e. $G_T=G_*A^2=G_*e^{2\alpha_Y\phi}$, and $\phi<0$ **decreases** $G_T$. The tensor coupling $G_T$, the operationally measured $G_{\rm lab}$ (Section 1.5) and the bare $2|\phi|$ proxy are three different quantities; only the magnitude $2\alpha_Y|\phi|$ is used, and the sign and reference state are not propagated. For comparison, the Kneller–Steigman fits (New J. Phys. 6, 117, 2004) give $\delta Y_p\approx0.16\,(S-1)$ and $\delta\ln({\rm D/H})\approx1.6\,(S-1)$ with $S=H'/H$: the coefficient 1.6 matches, 0.08 is a factor 2 smaller (immaterial at $|\delta H/H|\lesssim10^{-25}$). Input values: $Y_p^{\rm std}=0.24709$ is the PRIMAT central value (Pitrou, Coc, Uzan & Vangioni 2018, Phys. Rep. 754, 1); $10^5\,{\rm D/H}^{\rm std}=2.5147$ is an **arithmetic proxy** (review P02): a linearized benchmark value without verified single-paper provenance, not a sourced BBN prediction (PRIMAT 2018 gives 2.459); $Y_p^{\rm obs}=0.245\pm0.003$ is the PDG 2024 BBN review recommendation (Fields, Molaro & Sarkar), consistent with Aver et al. 2021 (JCAP 03, 027: $0.2453\pm0.0034$); $10^5\,{\rm D/H}^{\rm obs}=2.527\pm0.030$ is Cooke, Pettini & Steidel 2018 (ApJ 855, 102) — the CSV field name `DH_e5_obs_PDG` is a legacy label (PDG 2024 recommends $2.547\pm0.025$).

<<TABLE_BBN>>

- $Y_p=<<Yp>>$ versus $<<YpO>>\pm<<YpE>>$ (pull $<<pYp>>$); $10^5\,{\rm D/H}=<<DH>>$ versus $<<DHO>>\pm<<DHE>>$ (pull $<<pDH>>$).
- Maximum values over the three epochs: $H/(c/L)\le<<HbbnMax>>$, $|\bar\phi_{\rm attr}|\le<<phibbnMax>>$, $\delta H/H\le<<dHbbnMax>>$, $|\delta Y_p|\le<<dYpMax>>$. The qualitative expectation of a tiny adiabatic correction survives; this is not a general BBN bound.

---

## 7. Local Yukawa constant-density bulk benchmarks (`expm1`)

The rows are bulk products $|\phi|_{\rm bulk}=\alpha_Y\gamma_L\rho$ for constant density. The half-space edge value $|\phi|_{\rm bulk}/2$ for $R\gg L$ is given only where a boundary geometry is defined (mean solar density, mean Earth density, white dwarf interior; status `half-space boundary`). The homogeneous FLRW, cosmic void and Milky Way disk rows have no boundary geometry here: their edge cell is n.a. (CSV/JSON value null, `boundary_status` = `n.a. (no boundary geometry)`, edge $\Delta H_0$ null); their bulk value is neither copied into the edge column nor halved (review N01). The values are not solved surface fields of the Sun or the Earth (that needs the real source and density profile). Clock factors use `expm1` to avoid round-off zeros; the apparent shift is $\Delta H_0^{\rm app}=67.40\,{\rm expm1}(\alpha_Y|\phi|)$. The row "V21 top-hat interior external benchmark" ($|\phi|=<<ns_phi>>$) is an **external input, not reproduced here**: an interior value of the V21 top-hat fixed-point profile (SLy `gradient`, $r=2.0046$ km, $x=0.041239$; APR4 `mass`, $r=2.7542$ km, $x=0.041233$; Section 1.7), not a surface value (edge column n.a.) and not a solution of the present action.

<<TABLE_LOC>>

- Mean solar density bulk: $\Delta H_0^{\rm app}=<<sun_dH0>>\;\text{km/s/Mpc}$ (not a mathematical zero).
- Only the external V21 interior benchmark reaches $\Delta H_0^{\rm app}=<<ns_dH0>>\;\text{km/s/Mpc}$ (fraction $<<ns_frac>>$ of the 5.4 km/s/Mpc gap); this is a conditional illustration; identifying a local clock factor with the cosmological $H_0$ needs a separate derivation with a defined observer convention.

---

## 8. No-go appendix

### 8.1 M2 / M3 branches

- `M2` (flat, uncompensated negative vacuum density): $E^2(0)=<<E2M2>>$, so the declared normalization to a positive $H_0$ today is unphysical. The root of $E^2(z)=0$ is $z\approx<<zr_full>>$ in the full background (matter-only approximation: $<<zr_mat>>$). A negative vacuum energy is not excluded in general; this specific normalized point is.
- `M3` (open, $\Omega_K=<<OmKM3>>$), with its **own** sound horizons re-integrated with the M3 $H(z)$ (review F03): $r_s(z_*)=<<rsM3>>\;\text{Mpc}$, $r_d(z_{\rm drag})=<<rdM3>>\;\text{Mpc}$, $D_M(z_*)=<<DMM3>>\;\text{Mpc}$, $100\theta_*=<<thM3x>>$ (formal geometric pull $<<pullM3x>>$ on the 1D $\theta_*$ prior, not a likelihood significance) and BAO $\chi^2=<<chi2M3x>>$: a very poor acoustic-geometry benchmark. (The superseded values with the frozen M4 horizons were $<<thM3old>>$ and $<<chi2M3old>>$.)

### 8.2 The $67.36\,(1+1/4\pi)$ coincidence

$67.36\,(1+1/4\pi)=<<H0beta>>\;\text{km/s/Mpc}$ and $67.36\,e^{1/4\pi}=<<H0exp>>\;\text{km/s/Mpc}$ are correct arithmetic identities, but no dynamical $H_0$ correction follows from them: the clock factor is set by the solved potential $|\phi|$ (Sections 1 and 7), not by the bare coupling $\beta_{\rm Yuk}$.

### 8.3 Inverse distance ladder with profiled statistics

<<TABLE_NOGO>>

- The $\chi^2_{\rm BAO}$ column contains only the DESI $12\times12$ term. Row 4 adds the priors: $\chi^2_{\rm joint}=\chi^2_{\rm BAO}+\chi^2_{r_d\,\rm prior}+\chi^2_{H_0\,\rm prior}=<<R4bao>>+<<R4rd>>+<<R4h0>>=<<R4joint>>$; its $r_d$ entry quotes the **input prior error** $\pm0.26$ Mpc, not a fitted uncertainty. Row 5 is the rigid all-$z$ rescaling with unchanged $r_d$, not a late-transition model (Section 3).
- The fixed-$P$ illustration (row 6: $P=<<PM4>>\;\text{km/s}$, $r_d=<<rdfixP>>\;\text{Mpc}$) is a conditional arithmetic example; $P=<<PM4>>$ is not the BAO best-fit product, and its earlier sigma label is withdrawn.
- The profiled BAO requirement at $H_0=72.8$ is $r_d=<<rdP2>>\pm<<rdP2s>>\;\text{Mpc}$ (two-parameter control), and the joint BAO + Planck $r_d$ + assumed local $H_0$ benchmark prior control gives $\Delta\chi^2=<<dchi2J>>$.
- The rigid rescaling at unchanged shape and $r_d$ is strongly disfavoured by BAO ($\chi^2=<<chi2R>>$, $p=<<pR>>$).

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

<<DATA_FILES>>

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
| 5.3 | $\theta_*$-matched model relabelled to $\Omega_m=<<MTom>>$, $\omega_{cb}=<<MTwcb>>$; inverse-ladder figure uses $\sqrt{\Delta\chi^2}$ wording; provenance statement reworded | 3, 8.3, header |


---

## 12. Revíziós napló 2 / Revision log 2 (second review, F01–F12)

| Item | Change | Targeted check |
|:-----|:-----------------------------|:-----------------------------|
| F01 | $G_J=GA^{-2}$ replaced by $1/G_T=A^{-2}/G_*$, $G_T=G_*A^2$; $\phi<0$ decreases $G_T$; $G_T$, $G_{\rm lab}$ and the bare $2\lvert\phi\rvert$ proxy separated (Section 6, code comments, R6 metadata) | text and JSON `R6_bbn_proxy_provenance.clock_proxy` contain the new relation |
| F02 | G dictionary after Section 1.5 ($G_*$, $G_T$, $G_{\rm lab}\approx G_*(1+q)$, $G_N:=G_*$) | $1/(1+q)=<<Gratio>>$ computed in `run_step26_revision_controls.py` |
| F03 | M3 uses its own $r_s$, $r_d$ | $r_s=<<rsM3>>$, $r_d=<<rdM3>>$ Mpc, $100\theta_*=<<thM3x>>$, $\chi^2_{\rm BAO}=<<chi2M3x>>$ (assertion against the review values) |
| F04 | step2 status strings and CLI text rewritten: "BAO fit", "deviation on illustrative prior", "pull vs fixed target"; no Planck exclusion or curvature-likelihood wording | automatic string scan in `run_step26_revision_controls.py` |
| F05 | $72.8\pm1.0$ labelled assumed local $H_0$ benchmark prior; JSON key `H0_local_benchmark_prior_km_s_Mpc`; figure labels changed | labels no longer attribute the value to SH0ES (CSV/JSON/PNG text) |
| F06 | Section 1.5 sentence on the $2\alpha_Y^2\beta_{\rm Yuk}$ row | text check |
| F07 | RSD legend shows the canonical maximum | $<<dfsCanF>>$ in the legend and JSON |
| F08 | lowest-node sentence removed from Section 4, identity kept | text check |
| F09 | SNe figure: low-$\Omega_m$ offset curve added, "profiled constant zero-point offset", $N=18$, $\nu=17$ | low-$\Omega_m$ with offset $\chi^2=<<chi2LowSh>>$ |
| F10 | row renamed "V21 top-hat interior external benchmark", edge column n.a., provenance added | SLy $r=2.0046409680935633$ km, $x=0.041239199439254325$; APR4 $r=2.75417748649486$ km, $x=0.041233030079658145$ |
| F11 | PDF builder and `run_all.py --pdf` exit nonzero without pandoc/pdflatex or without a newly written PDF; `--tex-only` has its own status | test with a `PATH` lacking pdflatex |
| F12 | `theta_star_100_rs_only_shifted_fixed_DM` (renamed) and `theta_star_100_both_endpoints_shifted` (new) | low-$\Omega_m$: $<<thLowRsOnly>>$ and $<<thLowBoth>>$ |

## 13. Revíziós napló 3 / Revision log 3 (third re-audit, N01, N02, T01)

| Item | Change | Targeted check |
|:-----|:-----------------------------|:-----------------------------|
| N01 | Section 7 table: the "$R\gg L$ edge" cell of the homogeneous FLRW, cosmic void and Milky Way disk rows is n.a. (no boundary geometry); `step1_6_local_yukawa_environments.csv` and the JSON carry null in `abs_phi_surface_R_gg_L` and `Delta_H0_apparent_surf_km_s_Mpc` and the new per-row field `boundary_status`; bulk values unchanged (not halved) | assertion in `run_step26_revision_controls.py` (JSON `review3_N01_N02_T01`); text regression test |
| N02 | Section 1.7 "Carried over unchanged": the parenthesis now separates the background scalar-energy correction $\Omega_\phi\sim10^{-47}$ from the direct Yukawa perturbation of RSD ($\delta G_{\rm eff}/G_*\approx<<N02dG>>$, $\delta f\sigma_8\approx<<N02dfs>>$, Section 5); numbers taken from the JSON | text regression test (HU and EN sentence) |
| T01 | New `test_step26_text_regression.py`, run by `run_all.py` after the Markdown generation, nonzero exit on failure: checks the HU/EN F06 sentence, the F08 identity wording, the absence of the superseded wrong sentences (V26 factor-2 growth-law attribution; the identity value attributed to the lowest table node), the N01 table cells, the N02 sentence and the canonical legend string of the RSD figure source (F07); the comment in `run_step26_revision_controls.py` updated | mutation test on a scratch copy: reinserting a wrong sentence gives a nonzero exit |
"""

TEMPLATE_HU = r"""# 26. lépés (V27): Kozmológiai $H(z)$, DESI DR1 BAO, Planck 2018 akusztikus skála, szintetikus SNe Ia alakkontroll, variációs RSD $f\sigma_8$, algebrai BBN-benchmark és inverz távolságlétra no-go audit

**Szerző:** Szabolcsi Dávid (`ORCID: 0009-0005-6300-074X`)  
**Dátum:** 2026-10-28 (a szöveg az Academia Co-scientist bírálata, az azt követő újrabírálat és az F01–F12 második bírálat, valamint a harmadik újraaudit N01, N02, T01 pontjai után átdolgozva: 2026-10-29, lásd 11., 12. és 13. fejezet)  
**Előzmény DOI-k:**

- **V24:** Zenodo DOI `10.5281/zenodo.23680715` (szigorúan V24)
- **V25:** Zenodo DOI `10.5281/zenodo.23699685` (szigorúan V25, `step1`–`step22`)
- **V26:** Zenodo DOI `10.5281/zenodo.23727585` (szigorúan V26, `step1`–`step25`)

**A számok forrása:** a számítási eredményeket a JSON/CSV-forrásokból generáljuk (`tal_hubble_desi_bao_summary.json` és a 26. lépés CSV-táblázatai, a közös `tal_cosmology_common.py` háttérmodullal készítve) a `generate_step26_markdown.py` segítségével. A rögzített bemeneti állandók (adatértékek, priorok, irodalmi benchmarkok, külső bemenetek) a kódban szerepelnek, és a szöveg hivatkozik rájuk.

---

## 0. Hatókör és összefoglaló

A modul egyetlen konkrét ágat vizsgál: a V26 `Step 25` M4 hátteret rögzített $L=<<L_km>>\;\text{km}$, $\alpha_Y=1$, $\beta_{\rm Yuk}=1/(4\pi)=<<beta>>$ mellett, a skaláregyenlet adiabatikus partikuláris megoldásán, elhanyagolható szabad skaláros oszcillációs energiával. A szűkebb, védhető állítás: **ebben az ágban a $\beta_{\rm Yuk}$ csatolás puszta behelyettesítése $H_0$ szorzójaként nem levezetett mechanizmus a $67.4\to 72.8\;\text{km/s/Mpc}$ eltolásra.** Ez nem a TAL minden lehetséges kozmológiai megoldásának kizárása.

Az átdolgozás a bírálati pontokat így valósítja meg:

- egyetlen hatás, egyetlen skaláregyenlet (SI egységek, kozmikus másodperc), egyetlen attraktorformula és egyetlen $G_{\rm eff}$-normalizálás, a szövegben és a kódban azonosan (1. fejezet);
- explicit ágkiválasztás (adiabatikus attraktor, $C=0$ szabad oszcillációs amplitúdó) a kizárólagosság állítása helyett (1.3);
- DESI DR1 BAO a kerekítetlen Cobaya 12-elemű vektorral és $12\times 12$ kovarianciával, profilozott kontrollokkal (2. fejezet);
- minden CMB-modellhez és az M3-kontrollhoz modellspecifikus, újraintegrált $r_s(z_*)$ és $r_d(z_{\rm drag})$; a korábbi $\theta_*$- és $r_d$-alapú szigma-kizárási címkék törölve (3. és 8. fejezet);
- a szupernóva-rész új neve: **szintetikus távolságmodulus-alakkontroll** (4. fejezet);
- RSD: diagonális benchmark és elsőrendű variációs $\delta f\sigma_8$ (5. fejezet);
- BBN: csak algebrai nagyságrendi benchmark (6. fejezet);
- lokális környezetek: állandó sűrűségű bulk-benchmarkok `expm1`-gyel (7. fejezet);
- az R1–R6 és 5.1–5.3 újrabírálati pontok és a változtatások megfeleltetése a revíziós naplóban található (11. fejezet); a második bírálat F01–F12 pontjait a Revíziós napló 2 (12. fejezet), a harmadik újraaudit N01, N02, T01 pontjait a Revíziós napló 3 (13. fejezet), a P01–P09 nyitott korlátokat a 9.1 fejezet tartalmazza.

---

## 1. Egységes hatás, skaláregyenlet és $G_{\rm eff}$-normalizálás

### 1.1 Hatás és előjelkonvenció

$(-,+,+,+)$ szignatúrával, $g_{\mu\nu}$ Einstein-keretbeli metrikával és $\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}$ Jordan-keretbeli anyagmetrikával, $A(\phi)=e^{+\alpha_Y\phi}$:

$$
S=\int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}R-\frac{c^4}{4\pi G\beta_{\rm Yuk}}\left(\frac{1}{2}g^{\mu\nu}\partial_\mu\phi\,\partial_\nu\phi+\frac{\phi^2}{2L^2}\right)-\rho_{\rm vac}^{\rm eff}c^2\right]+S_m\!\left[\Psi_m,\,A^2(\phi)g_{\mu\nu}\right].
$$

$T^{\mu\nu}=(2/\sqrt{-g})\,\delta S_m/\delta g_{\mu\nu}$ és $T=-\rho_m c^2+3p_m$ mellett a $\phi$ szerinti variáció (egység: m$^{-2}$):

$$
\Box\phi-\frac{\phi}{L^2}=-\frac{4\pi G}{c^4}\,\alpha_Y\beta_{\rm Yuk}\,T=+\frac{4\pi G}{c^2}\,\alpha_Y\beta_{\rm Yuk}\left(\rho_m-\frac{3p_m}{c^2}\right).
$$

### 1.2 Homogén FLRW-egyenlet (kozmikus másodperc, SI)

Homogén térre $\Box\phi=-c^{-2}(\ddot\phi+3H\dot\phi)$, így mindkét oldal egysége s$^{-2}$:

$$
\ddot{\bar\phi}+3H\dot{\bar\phi}+\frac{c^2}{L^2}\bar\phi=-4\pi G\,\alpha_Y\beta_{\rm Yuk}\left(\rho_m-\frac{3p_m}{c^2}\right).
$$

A skalár Compton-rátája $c/L=<<mphic>>\;\text{s}^{-1}$, és $H_0L/c=<<H0overm>>$.

### 1.3 Attraktor és ágkiválasztás

Mivel $4\pi\beta_{\rm Yuk}=1$, az adiabatikus partikuláris megoldás ($\ddot{\bar\phi},\,3H\dot{\bar\phi}\ll (c/L)^2\bar\phi$):

$$
\bar\phi_{\rm attr}(z)=-\alpha_Y\gamma_L\left(\rho_m(z)-\frac{3p_m(z)}{c^2}\right),\qquad \gamma_L=\frac{4\pi\beta_{\rm Yuk}GL^2}{c^2}=\frac{GL^2}{c^2}=<<gammaL>>\;\text{m}^3/\text{kg},
$$

így $\bar\phi_{\rm attr}(0)=<<phi0>>$ és $\bar\phi_{\rm attr}(z_*)\approx<<phistar>>$. Az általános homogén megoldás szabad oszcillációt is tartalmaz:

$$
\delta\phi\simeq C\,a^{-3/2}\cos\!\left[(c/L)\,t+\vartheta\right],\qquad \langle\rho_{\delta\phi}\rangle\propto a^{-3}.
$$

A $H/(c/L)\ll 1$ feltétel az adiabatikus közelítést támasztja alá, de ezt a módust nem tünteti el. **Az elemzés ezért a $C=0$ ($\Omega_{\delta\phi,0}=0$) ágat választja: az adiabatikus partikuláris megoldást elhanyagolható szabad oszcillációs energiával.** A skaláregyenlet kezdetiérték-integrálása nem történik meg, és az alábbi következtetések csak erre az ágra érvényesek.

### 1.4 A skalár energiasűrűsége

A hatásból $\rho_\phi=\dfrac{c^2}{8\pi G\beta_{\rm Yuk}}\left[(\dot\phi/c)^2+(\phi/L)^2\right]$ (tömegsűrűség), így az attraktoron:

$$
\Omega_{\phi,\rm attr}(0)=\frac{(c/L)^2\,\bar\phi_{\rm attr}(0)^2}{3\beta_{\rm Yuk}H_0^2}=<<Omphi0>>.
$$

### 1.5 A $G_{\rm eff}$ normalizálása

Statikus pontszerű tömegre $\phi(r)=-\alpha_Y\beta_{\rm Yuk}\,\dfrac{GM}{c^2r}\,e^{-r/L}<0$, és az $A=e^{\alpha_Y\phi}$ anyagcsatolás $\alpha_Y^2\beta_{\rm Yuk}$ ötödikerő-erősséget ad. Lineáris perturbációkra, $k$ együttmozgó hullámszámmal:

$$
\frac{G_{\rm eff}(k,a)}{G_N}=1+\alpha_Y^2\beta_{\rm Yuk}\,\frac{(kL/a)^2}{1+(kL/a)^2}.
$$

Végig ezt a kanonikus normalizálást használjuk. A $2\alpha_Y^2\beta_{\rm Yuk}$ sor külön, kétszeres erősségű érzékenységi kontroll; nem a publikált V26 Step25 növekedési törvénye (5. fejezet). A lokális sajátidő $d\tau_J=A(\phi)\,d\tau_E=e^{-\alpha_Y|\phi|}d\tau_E$.

**G-szótár (F02 bírálati pont).** Négy jelölést választunk szét:

| Jelölés | Jelentés |
|:-----|:-----------------------------|
| $G_*$ | az 1.1 fejezet hatásában szereplő $G$ (Einstein-keret) |
| $G_T$ | Jordan-kereti tenzorcsatolás: $g_J=A^2g_E$ mellett a Ricci-tag $A^{-2}R/(16\pi G_*)$, így $1/G_T=A^{-2}/G_*$, $G_T=G_*A^2$ ($\phi<0$ esetén $G_T$ csökken) |
| $G_{\rm lab}$ | operatívan mért $G$: $r\ll L$ és árnyékolatlan kis testek esetén $\Phi=-G_*M(1+q\,e^{-r/L})/r$, így $G_{\rm lab}\approx G_*(1+q)$, $q=\alpha_Y^2\beta_{\rm Yuk}=<<qval>>$ |
| $G_N$ | az 5. fejezet perturbációs egyenletének jelölése: $G_N:=G_*$ |

$kL/a\ll1$ skálákon ezért, ha a csatolást laboratóriumi érték kalibrálja, a növekedési egyenletben $G_{\rm eff}/G_{\rm lab}=1/(1+q)=<<Gratio>>$ szerepel; az 5. fejezet $10^{-41}$-es állítása formális kontroll ugyanazon csupasz $G_*$ mellett. $f\sigma_8$-kizárási állítást nem teszünk, árnyékolási ágat nem vezetünk be.

### 1.6 Háttérnormalizálás és megfigyelői keret

Minden program egyetlen háttérmodult használ (`tal_cosmology_common.py`). A sima tömeges-neutrínó interpoláció önmagában $E^2(0)=<<E2unnorm>>$ értéket ad (eltérés: $<<E2offset>>$); a sötétenergia-lezárás ezt elnyeli, így $E^2(0)=1$ pontosan és $H(0)=H_0$. Ez az eltérés ($\Omega_\nu$-ben $<<NUdrift>>$) az interpoláció $z=0$-beli normalizálási sodródása, nem a közelítés hibája. A tömeges neutrínófaj Fermi–Dirac energiaintegráljához képest (azonos relativisztikus és nemrelativisztikus határesetekkel) a $p=1.83$ interpoláció relatív hibája $z=1000$-nél $<<NUrel1000>>$ a tömeges neutrínók sűrűségében és $<<NUH1000>>$ az összesített $H$-ban; a $0\le z\le10^4$ rácson a legnagyobb abszolút értékek $<<NUmax>>$ (tömeges komponens, $z\approx<<NUmaxZ>>$) és $<<NUHmax>>$ (összesített $H$, $z\approx<<NUHmaxZ>>$). Ez a $\Lambda$CDM-szerű háttér közelítési hibája, minden összevetett modellnél közös, és elkülönül az elméleti TAL-korrekciótól ($\sim 10^{-46}$).

<<TABLE_NU>>

 A számítások $\Lambda$CDM-szerű $H(z)$ függvényeket és külön órafaktort használnak; a Jordan-keretbe való leképezés ($dt_J=A\,dt_E$, $a_J=A\,a_E$, $H_J=(H_E+\alpha_Y\dot\phi_E)/A$) ennél a $|\phi|$-nél vezetőrendű azonosság, nem numerikusan megoldott, csatolt Friedmann–skalár–anyag rendszer.

### 1.7 A kinetikus normalizálás verzióegyeztetése (R1 újrabírálati pont)

Az alábbi útvonalak a V26 `scalar-tensor-yukawa-fit/` könyvtárfához viszonyítottak (V26, DOI `10.5281/zenodo.23727585`).

- A 26. lépés első vázlata a kinetikus tagot $\frac{c^4}{16\pi G}\frac{1}{2\beta_{\rm Yuk}}(\partial\phi)^2=\frac{c^4}{32\pi G\beta_{\rm Yuk}}(\partial\phi)^2$ alakban írta, az 1.1 hatásában viszont $\frac{c^4}{8\pi G\beta_{\rm Yuk}}(\partial\phi)^2$ áll: $K_{\rm régi}/K_{\rm új}=<<R1Kratio>>$.
- A V26 a $\frac{c^4}{16\pi G\tilde\alpha_L}(\partial\phi_L)^2$ tagot használja $\tilde\alpha_L=\alpha_Y/(8\pi)$ mellett, azaz $c^4/(2G\alpha_Y)$-t (`step23_gw170817_and_lensing_en.md:19` és `:29`, `README_en.md:50`, `gw170817_and_lensing_audit.py:104–113`). $\alpha_Y=1$, $\beta_{\rm Yuk}=1/(4\pi)$ mellett ez megegyezik az új együtthatóval ($K_{\rm V26}/K_{\rm új}=<<R1Kv26>>$), és a V26 ötödikerő-erőssége, $2\tilde\alpha_L\alpha_Y^2=<<R1betaV26>>$, egyenlő $\alpha_Y^2\beta_{\rm Yuk}$-kal. A V26 skaláregyenlete, $(\Box-L^{-2})\phi=-(\alpha_Y^2G/c^4)\varepsilon_m$ (`step25_cmb_and_cosmology_en.md:37`), és $\phi_{\rm attr}=\alpha_Y^2\gamma_L\rho$ (`cmb_and_cosmology_audit.py:327`, `:620`) a $\phi_{\rm V26}=-\phi_{\rm V27}$ előjelkonvenciótól eltekintve egyezik az 1.3 fejezettel. **Következtetés:** a V26 az új normalizálást használta; a $c^4/(32\pi G\beta)$ alak az első vázlat átírási inkonzisztenciája volt, a számított értékek már az új normalizálást követték. $\alpha_Y\neq1$ esetén a két verzió eltér, mert a V26 $\beta_{\rm Yuk}=\alpha_Y^3/(4\pi)$ kötést használ, itt pedig $\beta_{\rm Yuk}$ független paraméter.
- A $\phi_{\rm régi}=2\phi_{\rm új}$ térátdefiniálás a régi kinetikus tagot az újba viszi, de $A=e^{\alpha_{\rm régi}\phi_{\rm régi}}=e^{2\alpha_{\rm régi}\phi_{\rm új}}$, ezért $\alpha_{\rm új}=2\alpha_{\rm régi}$ szükséges. Azonos $\alpha_Y=1$ mellett a régi alak $4\alpha_Y^2\beta_{\rm Yuk}=<<R1aold>>$ ötödikerő-erősséget adna $\alpha_Y^2\beta_{\rm Yuk}=<<R1anew>>$ helyett.
- Változatlanul átvihető: minden háttér-, BAO-, CMB-, SNe-alak- és RSD-szám (a háttér skalárenergia-korrekciója $\Omega_\phi\sim10^{-47}$; az RSD közvetlen Yukawa-perturbációja külön, az 5. fejezetben számszerűsített elhanyagolható tag: $\delta G_{\rm eff}/G_*\approx<<N02dG>>$, $\delta f\sigma_8\approx<<N02dfs>>$), $\bar\phi_{\rm attr}(0)$, $\Omega_{\phi,\rm attr}(0)$ és a bulk lokális benchmarkok $\alpha_Y=1$ mellett. A V26 háttere emellett egy $(\phi/4\pi)^2$ energiakonvenciót is használt (`cmb_and_cosmology_audit.py:330`, `:623`), amely $\Omega_\phi=<<R1Omv26>>$ értéket ad a hatással konzisztens $<<Omphi0>>$ helyett; mindkettő elhanyagolható.
- A kétszeres növekedési konvenció ($2\alpha_Y^2\beta_{\rm Yuk}$) új címkéje **érzékenységi kontroll**, nem egyenértékű konvenció.
- A V21 top-hat belső külső benchmark, $|\phi|=<<ns_phi>>$ (7. fejezet), a V21 top-hat fixpont-profil ($x=L^2(G/c^2)\bar\rho$) belső értéke (`stellar_tal_v21.py:140–160`, `data/v21/profiles_L10km.csv`: SLy, `gradient` mód, $r=2.0046409680935633$ km, $x=0.041239199439254325$; APR4, `mass` mód, $r=2.75417748649486$ km, $x=0.041233030079658145$), nem felszíni érték, és **nem automatikusan megoldása az új Helmholtz–Yukawa hatásnak**; külső bemenet marad.

---

## 2. DESI DR1 BAO: $12\times 12$ kovariancia és profilozott kontrollok

A 12 elsődleges megfigyelés a Cobaya DESI DR1 vektora és kovarianciája (`desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`, `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`). A külön LRG3 és ELG1 egytracer-sorok csak diagnosztikák, a likelihoodba nem kerülnek; a $z=0.845$ ELG1 sor ($D_V/r_d=18.33\pm0.62$) **igazolt elsődleges forrás nélküli, választott illusztráció**, nem DESI-mérés (nem szerepel a DESI III Table 18-ban és a DESI VI Table 1-ben). Az elsődleges binek $\chi^2$-e a $12\times 12$ kovarianciával, a (dőlt) diagnosztikai soroké a kerekített táblázati hibákkal számolt. A táblázat jóslatai anizotrop binekre $D_M/r_d$; $D_H/r_d$, izotrop binekre $D_V/r_d$.

<<TABLE_BAO>>

- Rögzített M4-pont ($H_0=67.40$, $r_d=<<rd>>\;\text{Mpc}$): $\chi^2=<<chi2M4>>$ 12 megfigyelésre, $p=<<pM4>>$ (kerekített Table 1 bemenettel: $<<chi2M4t>>$).
- Merev $H_0=72.8$ változatlan alakkal és $r_d$-vel: $\chi^2=<<chi2R>>$, $p=<<pR>>$ (kerekített Table 1 bemenettel: $<<chi2Rt>>$).

**Profilozott kontrollok** (helyi Fisher-hibák; a $p$-értékek a megadott Gauss-likelihood feltételes gyakorisági értékei, nem a TAL igazságának valószínűségei):

<<TABLE_PROF>>

- Rögzített $E(z)$ mellett a BAO csak a $P=H_0r_d$ szorzatot korlátozza (rang-1 információ a $(H_0,r_d)$ síkon).
- Kétparaméteres kontroll: $\Omega_m=<<Om2>>\pm<<Om2s>>$, $P=<<P2>>\pm<<P2s>>\;\text{km/s}$, $\chi^2_{\min}=<<chi2P2>>$ (10 dof, $p=<<pP2>>$).
- Szükséges $r_d$ pontosan $H_0=72.8$ mellett: $<<rdP2>>\pm<<rdP2s>>\;\text{Mpc}$; a $H_0=72.8\pm1.0$ bizonytalanságot is propagálva: $\pm<<rdP2sH>>\;\text{Mpc}$. A $72.8\pm1.0$ **feltételezett lokális $H_0$ benchmark-prior** (választott célérték, F05 bírálati pont); nem a Riess és mtsai. 2022-es $73.04\pm1.04$ értéke, és konkrét SH0ES-kiadást nem tulajdonítunk neki.
- A modellfüggő Planck $r_d=147.09\pm0.26\;\text{Mpc}$ priorral: $H_0=<<H0Pl>>\pm<<H0Pls>>\;\text{km/s/Mpc}$.
- Szemléltető közös teszt (BAO + Planck $r_d$ prior + lokális $H_0=72.8\pm1.0$, $\Omega_m,H_0,r_d$ profilozva): legjobb illesztés $\Omega_m=<<OmJ>>$, $H_0=<<H0J>>$, $r_d=<<rdJ>>\;\text{Mpc}$, $\chi^2_{\rm közös}=<<chi2J>>$, felbontva: $\chi^2_{\rm BAO}=<<R4bao>>$, $\chi^2_{r_d\,\rm prior}=<<R4rd>>$, $\chi^2_{H_0\,\rm prior}=<<R4h0>>$; $\Delta\chi^2=<<dchi2J>>$ a csak-BAO minimumhoz képest ($\sqrt{\Delta\chi^2}=<<sigJ>>$). Ez nem Planck+TAL kizárási szignifikancia.

![DESI DR1 BAO, Step 25 M4 és profilozott kontrollok (`step1_1_desi_dr1_bao_plot.png`)](step1_1_desi_dr1_bao_plot.png){ width=96% }

---

## 3. Planck 2018 akusztikus skála modellspecifikus hanghorizontokkal

Minden modellnél az $r_s(z_*=1089.80)$ és $r_d(z_{\rm drag}=1059.94)$ a modell saját $H(z)$ és $c_s(z)$ függvényével újraintegrált (rögzített $z_*$, $z_{\rm drag}$ kontrollközelítésként). Az $(\omega_b,\omega_{cb},100\theta_*)$ 3D Gauss-prior **szemléltető, feltételezett prior** (Planck 2018 marginális középértékek és hibák, közelítő korrelációs mátrixszal, amelynek láncproveniencája itt nincs dokumentálva); nem a Planck TT,TE,EE+lowE+lensing likelihood. Referencia: $100\theta_*=<<th100ref>>\pm<<th100sig>>$.

<<TABLE_CMB>>

- Step 25 M4: $r_s=<<rs>>\;\text{Mpc}$, $r_d=<<rd>>\;\text{Mpc}$, $100\theta_*=<<th100M4>>$ (pull $<<pullthM4>>$), $\chi^2_{\rm 3D}=<<chi23DM4>>$.
- Minden $z$-n merev $H\to sH$ átskálázás változatlan $c_s$ mellett: $r_s\to r_s/s=<<rsRig>>\;\text{Mpc}$, $r_d=<<rdRig>>\;\text{Mpc}$ és $D_M\to D_M/s$, így $100\theta_*=<<thRig>>$ változatlan. Csak $D_M$ skálázása $<<thRigFix>>$ értéket adna; ezt az inkonzisztens kombinációt már nem használjuk tesztként.
- Sík $H_0=72.8$, $\Omega_m=0.244828$: az újraintegrálás $r_s=<<rsLow>>$, $r_d=<<rdLow>>\;\text{Mpc}$ és $100\theta_*=<<thLow>>$ értéket ad (rögzített korai $r_s$-sel $<<thLowFix>>$ lenne). Ez **nem** jelenti, hogy a modell átmegy a Planck-teszten: $\omega_{cb}=<<wcbLow>>$ eltér a priortól, és a BAO $\chi^2$ az újraintegrált $r_d$-vel $<<chi2Low>>$ (a régi rögzített $r_d$-vel: $<<chi2LowFix>>$). A korai és késői részt együtt kell újrailleszteni.
- **Két $z_*$-eltolási diagnosztika (F12 bírálati pont).** A Hu–Sugiyama-illesztés $z_*\to z_*+\Delta z_*$ eltolásával (ennél a modellnél $z_*+\Delta z_*=<<zsLowShift>>$) a CSV két külön nevű értéket tárol: `theta_star_100_rs_only_shifted_fixed_DM` $=100\,r_s(z_*+\Delta z_*)/D_M(1089.80)=<<thLowRsOnly>>$ (csak az $r_s$ végpontja eltolva) és `theta_star_100_both_endpoints_shifted` $=100\,r_s(z_*+\Delta z_*)/D_M(z_*+\Delta z_*)=<<thLowBoth>>$. Egyik sem rekombinációs megoldás.
- **Késői átmenet (R2 újrabírálati pont).** A korábbi „csak késői” sor a $D_M(z_*)/s$ értéket változatlan korai $r_s,r_d$-vel kombinálta ($100\theta_*=<<R2illus>>$); nincs átmeneti vöröseltolódása, ezért **inkonzisztens geometriai szemléltetésként sorolódik át** (dőlt sor a táblázatban), és kikerül a konzisztens modellekre vonatkozó állításokból. Definiált kontroll a szakaszos $H(z)=s\,H_{\rm M4}(z)$ ($z<z_t$), illetve $H_{\rm M4}(z)$ ($z\ge z_t$), $s=72.8/67.4=<<R2s>>$, numerikusan integrált $D_M(z_*)$-gal és változatlan $r_s=<<R2rs>>\;\text{Mpc}$-vel (M4: $D_M(z_*)=<<R2DMM4>>\;\text{Mpc}$). $H$ ugrik $z_t$-nél, így ez kinematikai kontroll, nem dinamikai modell; a pull oszlop formális 1D távolság a Planck-értéktől, nem likelihood-szignifikancia. $z_t\ge2.5$ esetén minden DESI-bin az átskálázott tartományba esik, így a BAO $\chi^2$ a merev értékkel egyezik: $<<chi2R>>$.

<<TABLE_R2>>

---

## 4. Szintetikus távolságmodulus-alakkontroll (nem Pantheon+ / Union3 likelihood)

Ez a rész **nem** használ megfigyelt SN-katalógust: <<nn>> vöröseltolódás-csomópont ($z=<<zmin>>$–$<<zmax>>$) reprezentatív $\sigma_\mu$ hibákkal, és egy $H_0=72.8$, $\Omega_m=0.334$ sík modellből számolt szintetikus referenciagörbe. Nincs SN-azonosító, megfigyelt $\mu$ és statisztikai+szisztematikus kovariancia. Az alábbi $\chi^2$ értékek modell–modell alakmetrikák, nem megfigyelési illeszkedési mértékek.

- Azonos $E(z)$ geometria esetén az eltolás a konstans $5\log_{10}(72.8/67.4)=+<<dmu>>\;\text{mag}$ azonosság lenne (pozitív előjel: M4 mínusz a 72.8-as referencia).
- A táblázat $z=<<zmin>>$-beli $<<R5node>>\;\text{mag}$ különbsége modellkülönbség (a referencia $\Omega_m=0.334$), nem ez az azonosság.
- **Szimmetrikus nullpontkezelés (R5 újrabírálati pont):** mindkét modell ugyanazt a profilozott konstans eltolást kapja (abszolút fényesség marginalizálása), $1/\sigma_\mu^2$ súlyokkal:

<<TABLE_R5>>

- Az M4 nagy csökkenése csak azt mutatja, hogy a konstans levonása után alakja közel van a referenciához; a csomópontok és a reprezentatív hibák szintetikusak, ezért a csomópontok vagy hibák változtatása nem jelent fizikai javulást.
- Az ábra mindkét modellt profilozott konstans nullpont-eltolással és anélkül is mutatja (az alacsony $\Omega_m$-ű modell eltolással: $\chi^2=<<chi2LowSh>>$); $N=18$ csomópont, egy profilozott eltolás után $\nu=17$ szabadsági fok. Az eltolás profilozott zavaró konstans, nem megfigyelési kalibráció (F09 bírálati pont).
- Kalibrálatlan SNe-adatokból önmagában nincs $H_0$-mérés ($H_0$–$M_B$ degeneráció). Valódi katalógusállításhoz a Pantheon+ vagy Union3 likelihoodot és kovarianciát kellene beolvasni; ez a csomag ezt nem teszi.

![Szintetikus SNe Ia távolságmodulus-alakkontroll (`step1_3_pantheon_sne_plot.png`)](step1_3_pantheon_sne_plot.png){ width=96% }

<<TABLE_SNE>>

---

## 5. RSD $f\sigma_8(z)$: diagonális benchmark és variációs TAL-korrekció

A $\delta_m''+\left(2+\frac{d\ln H}{d\ln a}\right)\delta_m'-\frac{3}{2}\Omega_{cb}(a)\frac{G_{\rm eff}(k,a)}{G_N}\delta_m=0$ növekedési egyenletet DOP853 integrálja (`rtol=1e-12`, `atol=1e-14`, dense output); $\sigma_{8,0}=0.8111$ rögzített Planck-bemenet, nem TAL-ból levezetett amplitúdó.

- $k=<<kcom>>\,h\,\text{Mpc}^{-1}$ mellett: $kL=<<kL>>$ ($(kL)^2=<<kL2>>$) $z=0$-nál; $\Delta G_{\rm eff}/G_N=<<dGcan>>$ (kanonikus $\alpha_Y^2\beta_{\rm Yuk}$) vagy $<<dG2x>>$ (kétszeres konvenció).
- float64-ben $1+O(10^{-41})$ pontosan 1, ezért az $f\sigma_8^{\rm TAL}-f\sigma_8^{\rm GR}$ közvetlen különbség lebegőpontos azonosság miatt $0$, nem 40 jegyű mérés eredménye. Az elsőrendű variációs egyenlet a 12 ponton $\max\lvert\delta f\sigma_8\rvert=<<dfsCan>>$ (kanonikus) és $<<dfs2x>>$ (kétszeres érzékenységi kontroll) értéket ad; az újrabírálat független, átskálázott ellenőrzése $1.666\times10^{-40}$-t erősített meg (maximum $1.6663681521\times10^{-40}$, relatív eltérés $\le5.1\times10^{-6}$). A float64-ág $G_{\rm eff}/G_N=1+\Delta G_{\rm eff,kan}/G_N$ alakot használ (5.1 újrabírálati pont). Az ábra jelmagyarázata a kanonikus $<<dfsCanF>>$ maximumot mutatja; a kétszeres érték csak a külön érzékenységi kontroll (F07 bírálati pont). Itt $G_N:=G_*$ (G-szótár, 1.5 fejezet).
- Diagonális $\chi^2$ (bin-kovarianciák, átfedések és AP-kezelés nélkül): M4 $\Omega_{cb}$ csomósodási forrással: $<<chi2RSD>>$; az egész anyaggal: $<<chi2RSDm>>$; sík $H_0=72.8$, $\Omega_m=0.2448$: $<<chi2RSDlow>>$ (12 pont). Ez benchmark, nem kovarianciás likelihood.

![RSD $f\sigma_8(z)$ diagonális benchmark (`step1_4_rsd_fsigma8_plot.png`)](step1_4_rsd_fsigma8_plot.png){ width=96% }

<<TABLE_RSD>>

---

## 6. Algebrai BBN nagyságrendi benchmark

Nukleáris reakcióhálózat (pl. PArthENoPE) nem fut. Az elemgyakoriságok linearizált standard értékek $\eta_{10}=<<eta10>>$ mellett, algebrai TAL-eltolással; a $H(T)$ egyszerű extrapoláció $e^\pm$-termodinamika és változó $g_*$ nélkül, a skalárforrás pedig csak a poranyagot tartalmazza, a termikus plazma trace-ét nem. A pullok csak a megadott megfigyelési hibával osztanak.

**Proxyk és források (R6 újrabírálati pont).** A $|\delta G/G|\approx2|\phi|$ óraproxy és a $0.08$ ($\delta Y_p\approx0.08\,\delta H/H$), illetve $1.6$ ($\delta\ln({\rm D/H})\approx1.6\,\delta H/H$) együtthatók **feltételezett nagyságrendi proxyk**, itt nem levezetettek. Keret és előjel (F01 bírálati pont): $g_J=A^2g_E$ mellett a Jordan-kereti Ricci-tag $A^{-2}/G_*$ szorzót hordoz ($A^{-2}$ a Planck-tömeg négyzetét szorozza), így a tenzorcsatolás $1/G_T=A^{-2}/G_*$, azaz $G_T=G_*A^2=G_*e^{2\alpha_Y\phi}$, és $\phi<0$ esetén $G_T$ **csökken**. A $G_T$ tenzorcsatolás, az operatívan mért $G_{\rm lab}$ (1.5 fejezet) és a csupasz $2|\phi|$ proxy három különböző mennyiség; csak a $2\alpha_Y|\phi|$ nagyságot használjuk, az előjelet és a referenciaállapotot nem propagáljuk. Összevetésként a Kneller–Steigman-illesztések (New J. Phys. 6, 117, 2004) $\delta Y_p\approx0.16\,(S-1)$ és $\delta\ln({\rm D/H})\approx1.6\,(S-1)$ értéket adnak, $S=H'/H$: az 1.6-os együttható egyezik, a 0.08 kétszer kisebb (ez $|\delta H/H|\lesssim10^{-25}$ mellett lényegtelen). Bemeneti értékek: $Y_p^{\rm std}=0.24709$ a PRIMAT központi értéke (Pitrou, Coc, Uzan & Vangioni 2018, Phys. Rep. 754, 1); $10^5\,{\rm D/H}^{\rm std}=2.5147$ **számtani proxy** (P02 bírálati pont): linearizált benchmarkérték igazolt egyedi forrás nélkül, nem forrásolt BBN-jóslat (a PRIMAT 2018 2.459-et ad); $Y_p^{\rm obs}=0.245\pm0.003$ a PDG 2024 BBN-áttekintés ajánlott értéke (Fields, Molaro & Sarkar), összhangban Aver és mtsai. 2021 eredményével (JCAP 03, 027: $0.2453\pm0.0034$); $10^5\,{\rm D/H}^{\rm obs}=2.527\pm0.030$ Cooke, Pettini & Steidel 2018 (ApJ 855, 102) értéke — a `DH_e5_obs_PDG` CSV-mezőnév örökölt címke (a PDG 2024 ajánlása $2.547\pm0.025$).

<<TABLE_BBN>>

- $Y_p=<<Yp>>$, összevetve: $<<YpO>>\pm<<YpE>>$ (pull $<<pYp>>$); $10^5\,{\rm D/H}=<<DH>>$, összevetve: $<<DHO>>\pm<<DHE>>$ (pull $<<pDH>>$).
- Maximumok a három korszakon: $H/(c/L)\le<<HbbnMax>>$, $|\bar\phi_{\rm attr}|\le<<phibbnMax>>$, $\delta H/H\le<<dHbbnMax>>$, $|\delta Y_p|\le<<dYpMax>>$. A kis adiabatikus korrekció kvalitatív várakozása fennmarad; ez nem általános BBN-korlát.

---

## 7. Lokális Yukawa állandó sűrűségű bulk-benchmarkok (`expm1`)

A sorok állandó sűrűségre vett $|\phi|_{\rm bulk}=\alpha_Y\gamma_L\rho$ szorzatok. Az $R\gg L$ félteres peremérték ($|\phi|_{\rm bulk}/2$) csak ott szerepel, ahol peremgeometria definiált (Nap átlagsűrűség, Föld átlagsűrűség, fehér törpe belseje; státusz: `half-space boundary`). A homogén FLRW-, a kozmikus üreg- és a Tejútrendszer-korong sornak itt nincs peremgeometriája: a peremcellájuk n.a. (CSV/JSON-ban null, `boundary_status` = `n.a. (no boundary geometry)`, a perem-$\Delta H_0$ null); a bulk értéküket nem másoljuk a peremoszlopba, és nem felezzük (N01 bírálati pont). Az értékek nem a Nap vagy a Föld megoldott felszíni terei (ahhoz a valódi forrás és sűrűségprofil kell). Az órafaktorok `expm1`-gyel készülnek, a kerekítési nullák elkerülésére; a látszólagos eltolás $\Delta H_0^{\rm app}=67.40\,{\rm expm1}(\alpha_Y|\phi|)$. A „V21 top-hat belső külső benchmark” sor ($|\phi|=<<ns_phi>>$) **külső bemenet, itt nem reprodukált**: a V21 top-hat fixpont-profil belső értéke (SLy `gradient`, $r=2.0046$ km, $x=0.041239$; APR4 `mass`, $r=2.7542$ km, $x=0.041233$; 1.7 fejezet), nem felszíni érték (a peremoszlop n.a.), és nem a jelen hatás megoldása.

<<TABLE_LOC>>

- Nap átlagsűrűségű bulk: $\Delta H_0^{\rm app}=<<sun_dH0>>\;\text{km/s/Mpc}$ (nem matematikai nulla).
- Csak a külső V21 belső benchmark éri el a $\Delta H_0^{\rm app}=<<ns_dH0>>\;\text{km/s/Mpc}$ értéket (az 5.4 km/s/Mpc rés $<<ns_frac>>$ hányada); ez feltételes szemléltetés; egy lokális órafaktor és a kozmológiai $H_0$ azonosítása külön levezetést igényel, definiált megfigyelői konvencióval.

---

## 8. No-go függelék

### 8.1 Az M2 / M3 ágak

- `M2` (sík, kompenzálatlan negatív vákuumsűrűség): $E^2(0)=<<E2M2>>$, így a deklarált, ma pozitív $H_0$-ra vett normalizálás nem fizikai. Az $E^2(z)=0$ gyöke a sugárzást és neutrínót is tartalmazó háttérben $z\approx<<zr_full>>$ (csak anyagos közelítés: $<<zr_mat>>$). A negatív vákuumenergia általában nincs kizárva; ez a konkrét normalizált pont igen.
- `M3` (nyílt, $\Omega_K=<<OmKM3>>$), az M3 $H(z)$-vel újraintegrált **saját** hanghorizontokkal (F03 bírálati pont): $r_s(z_*)=<<rsM3>>\;\text{Mpc}$, $r_d(z_{\rm drag})=<<rdM3>>\;\text{Mpc}$, $D_M(z_*)=<<DMM3>>\;\text{Mpc}$, $100\theta_*=<<thM3x>>$ (formális geometriai pull az 1D $\theta_*$ prioron: $<<pullM3x>>$, nem likelihood-szignifikancia) és BAO $\chi^2=<<chi2M3x>>$: nagyon rossz akusztikus geometriai benchmark. (A felülírt, rögzített M4-horizontos értékek: $<<thM3old>>$ és $<<chi2M3old>>$.)

### 8.2 A $67.36\,(1+1/4\pi)$ numerikus egybeesés

$67.36\,(1+1/4\pi)=<<H0beta>>\;\text{km/s/Mpc}$ és $67.36\,e^{1/4\pi}=<<H0exp>>\;\text{km/s/Mpc}$ helyes számtani azonosságok, de dinamikai $H_0$-korrekció nem következik belőlük: az órafaktort a megoldott $|\phi|$ potenciál adja (1. és 7. fejezet), nem a puszta $\beta_{\rm Yuk}$ csatolás.

### 8.3 Inverz távolságlétra profilozott statisztikával

<<TABLE_NOGO>>

- A $\chi^2_{\rm BAO}$ oszlop csak a DESI $12\times12$ tagot tartalmazza. A 4. sor a priorokat is hozzáadja: $\chi^2_{\rm közös}=\chi^2_{\rm BAO}+\chi^2_{r_d\,\rm prior}+\chi^2_{H_0\,\rm prior}=<<R4bao>>+<<R4rd>>+<<R4h0>>=<<R4joint>>$; $r_d$-bejegyzése a **bemeneti prior hibáját** ($\pm0.26$ Mpc) közli, nem illesztett bizonytalanságot. Az 5. sor a minden $z$-n merev átskálázás változatlan $r_d$-vel, nem késői átmeneti modell (3. fejezet).
- A rögzített-$P$ szemléltetés (6. sor: $P=<<PM4>>\;\text{km/s}$, $r_d=<<rdfixP>>\;\text{Mpc}$) feltételes számtani példa; a $P=<<PM4>>$ nem a BAO legjobb illesztett szorzata, korábbi szigma-címkéjét visszavontuk.
- A profilozott BAO-igény $H_0=72.8$ mellett $r_d=<<rdP2>>\pm<<rdP2s>>\;\text{Mpc}$ (kétparaméteres kontroll), a közös BAO + Planck $r_d$ + feltételezett lokális $H_0$ benchmark-prior kontroll pedig $\Delta\chi^2=<<dchi2J>>$ értéket ad.
- A változatlan alakú és $r_d$-jű merev átskálázást a BAO erősen elutasítja ($\chi^2=<<chi2R>>$, $p=<<pR>>$).

![Inverz távolságlétra a $(H_0,r_d)$ síkon (`step2_3_inverse_ladder_rd_H0_exclusion.png`)](step2_3_inverse_ladder_rd_H0_exclusion.png){ width=96% }

---

## 9. Mi marad meg, és mi nem igazolt

**Megmarad:**

- A rögzített $E(z)$ alak és változatlan $r_d$ mellett minden vöröseltolódáson 8%-kal megemelt $H(z)$ nagyon rosszul illik a BAO-hoz.
- A puszta $\beta_{\rm Yuk}$ nem egyenlő a forrásból megoldott skalárpotenciállal vagy egy órakalibrációs eltolódással.
- A rövid hatótávú, adiabatikus, elhanyagolható oszcillációs energiájú ágban a késői kozmológiai skalárkorrekció rendkívül kicsi.
- A konkrét M2 normalizálás nem fizikai, az M3 görbülete pedig rossz akusztikus geometriai benchmark.
- A profilozott konstans SNe-nullpont-eltolás és a $H_0$–$M_B$ degeneráció geometriai szemléltetése helyes.

**Nem igazolt ebben a csomagban:**

- közös Planck + Pantheon+/Union3 + RSD + BBN adatvalidáció;
- a korábbi szigma-címkék által sugallt általános fizikai kizárási erő;
- a skalár–tenzor elmélet kozmológiai lezárása tetszőleges megengedett kezdeti feltételekkel;
- új, TAL-specifikus empirikus előny a megfelelően illesztett $\Lambda$CDM-hez képest.

**Zárómondat.** A vizsgált L=10 km-es, rögzített csatolású, elhanyagolható szabad skaláros oszcillációs energiájú M4 ágban a lokális Yukawa-faktor nem szolgáltat levezetett mechanizmust a 72.8 km/s/Mpc globális Hubble-állandóhoz. A rögzített korai háttér és rögzített késői alak egyszerű globális átskálázását a BAO erősen elutasítja. Az eredmény nem átfogó, közös Planck–SNe–RSD–BBN likelihood-elemzés, és nem az elmélet összes kozmológiai ágára vonatkozó no-go tétel.

### 9.1 Nyitott korlátok / Open limitations

| Pont | Nyitott korlát | Állapot ebben a csomagban |
|:-----|:-----------------------------|:-----------------------------|
| P01 | A $z=0.845$ ELG1 pont ($18.33\pm0.62$) igazolt elsődleges forrása | igazolt elsődleges forrás nélküli, választott illusztráció; nem DESI-mérés; nincs a 12 pontos likelihoodban |
| P02 | A $10^5\,({\rm D/H})_{\rm std}=2.5147$ egyedi forrása | számtani proxy; empirikus standardmodell-teszthez igazolt bemenet és elméleti hiba kell |
| P03 | Csatolt skalár–Friedmann–anyag kezdetiérték-probléma és szabad homogén amplitúdók | a $C=0$ ág expliciten választott; nincs állítás minden kezdeti feltételről vagy minden kozmológiai TAL-ágról |
| P04 | Planck TT/TE/EE + lensing likelihood és rekombinációs számítás | a 3D prior feltételezett, $z_*$ és $z_{\rm drag}$ rögzített kontrollparaméter; geometriai benchmark, nem CAMB/CLASS-validáció |
| P05 | Valódi Pantheon+/Union3 adatvektor és kovariancia | szintetikus 18 csomópontos metrika, nem megfigyelési modellpreferencia |
| P06 | RSD-kovariancia, átfedések, AP-leképezés, levezetett amplitúdó, $k$-spektrum | $k=0.1\,h$/Mpc variáció rögzített $\sigma_{8,0}$-val; az $1.666\times10^{-40}$ erre az ODE-re, normalizálásra és kezdeti feltételekre vonatkozik |
| P07 | BBN termikus trace, változó $g_*$, $e^\pm$-termodinamika, reakcióhálózat | három mintaepocha feltételezett érzékenységekkel; nem elemgyakoriság-likelihood |
| P08 | Valódi nap- és földi sűrűségprofil, felszíni Yukawa-megoldás, a lokális órák távolságlétrába illesztése | bulk és félteres limit; a $\Delta H_0^{\rm app}$ táblázat számtani szemléltetés |
| P09 | A V26 könyv egésze, minden történeti konvenció és minden csillag-BVP | ez a bírálat a 26. lépésre vonatkozik; a V26 normalizálási azonossága csak célzott ellenőrzés |

---

## 10. Fájlok és reprodukálhatóság

<<DATA_FILES>>

- Generátorok: `tal_cosmology_common.py` (közös háttér), `run_step1_1.py`, `run_step1_2_and_1_3.py`, `run_step1_4_to_1_6.py`, `run_step2_nogo_and_report.py`, `run_step26_revision_controls.py` (újrabírálati kontrollok), `hubble_desi_bao_tal_audit.py`; szöveg: `generate_step26_markdown.py` + `step26_md_templates.py`, szöveges regressziós teszt `test_step26_text_regression.py`; PDF: `compile_step26_clean_pdfs.py`; futtató: `run_all.py` (sorrend és adat-hash ellenőrzés, `--skip-block1`), függőségek: `requirements.txt`.
- A csomag önállóan reprodukálható: mindkét DESI bemeneti fájl, az összes CSV- és PNG-kimenet, a szigorú JSON (`allow_nan=False`, a NaN `null`-ként írva) és a `SHA256SUMS.txt` benne van.
- Külső források: DESI VI (arXiv:2404.03002), DESI III (arXiv:2404.03000), Planck 2018 VI (arXiv:1807.06209), Cobaya `bao_data` tároló (DESI 2024 Gauss BAO középérték és kovariancia).

---

## 11. Revíziós napló

| Pont | Változtatás | Hol |
|:-----|:-------|:------|
| R1 | A kinetikus normalizálás egyeztetése a V26-tal (fájl:sor); $\phi_{\rm régi}=2\phi_{\rm új}$ esetén $\alpha_{\rm új}=2\alpha_{\rm régi}$ szükséges; az átvihető eredmények felsorolva; a kompaktcsillag-$\lvert\phi\rvert$ külső V21 belső értékként jelölve; a $2\alpha^2\beta$ új címkéje érzékenységi kontroll | 1.5, 1.7, 7 |
| R2 | Szakaszos $H=sH_{\rm M4}$ ($z<z_t$) kontroll $z_t\in\{2.5,10,100,1000\}$ értékekre; a $D_M/s$ sor inkonzisztens geometriai szemléltetésként átsorolva | 3, `step26_r2_piecewise_transition.csv` |
| R3 | Önállóan reprodukálható csomag: DESI-bemenetek hash-sel, minden CSV/PNG, `requirements.txt`, `run_all.py`, szigorú JSON; a régi `build_step26_pdfs.py` kikerült a csomagból | 10, README |
| R4 | A közös $\chi^2$ felbontása BAO-, $r_d$-prior- és $H_0$-prior-tagra; a táblázat oszlopa csak BAO; az $r_d$ $\pm0.26$ bemeneti prior hibaként jelölve | 2, 8.3 |
| R5 | Szimmetrikus nullpont-eltolás az M4 és az alacsony $\Omega_m$-ű modell számára; a legalsó csomópont értéke elválasztva az $5\log_{10}$ azonosságtól | 4, `step26_r5_sne_symmetric_offsets.csv` |
| R6 | A BBN óraproxy és együtthatók feltételezett proxyként jelölve; keret és előjel tárgyalva; források megadva | 6 |
| 5.1 | A float64-ág $1+\Delta G_{\rm eff,kan}/G_N$ alakot használ; abszolútérték-fejlécek; a független átskálázott ellenőrzés idézve | 5 |
| 5.2 | A $z=0$ sodródás elválasztva a Fermi–Dirac interpolációs hibától, amelyet kiszámoltunk | 1.6 |
| 5.3 | A $\theta_*$-hez illesztett modell átcímkézve: $\Omega_m=<<MTom>>$, $\omega_{cb}=<<MTwcb>>$; az inverz-létra ábra $\sqrt{\Delta\chi^2}$ megfogalmazást használ; a provenienciaállítás átfogalmazva | 3, 8.3, fejléc |


---

## 12. Revíziós napló 2 / Revision log 2 (második bírálat, F01–F12)

| Pont | Változtatás | Célzott ellenőrzés |
|:-----|:-----------------------------|:-----------------------------|
| F01 | $G_J=GA^{-2}$ helyett $1/G_T=A^{-2}/G_*$, $G_T=G_*A^2$; $\phi<0$ esetén $G_T$ csökken; $G_T$, $G_{\rm lab}$ és a csupasz $2\lvert\phi\rvert$ proxy szétválasztva (6. fejezet, kódkommentek, R6-metaadat) | a szöveg és a JSON `R6_bbn_proxy_provenance.clock_proxy` mezője az új összefüggést tartalmazza |
| F02 | G-szótár az 1.5 fejezet után ($G_*$, $G_T$, $G_{\rm lab}\approx G_*(1+q)$, $G_N:=G_*$) | $1/(1+q)=<<Gratio>>$, számolva a `run_step26_revision_controls.py`-ban |
| F03 | Az M3 saját $r_s$, $r_d$ értékeket használ | $r_s=<<rsM3>>$, $r_d=<<rdM3>>$ Mpc, $100\theta_*=<<thM3x>>$, $\chi^2_{\rm BAO}=<<chi2M3x>>$ (ellenőrzés a bírálati értékekkel szemben) |
| F04 | A step2 státuszszövegek és a CLI-szöveg átírva: „BAO fit”, „deviation on illustrative prior”, „pull vs fixed target”; nincs Planck-kizárási vagy görbületi-likelihood megfogalmazás | automatikus szövegellenőrzés a `run_step26_revision_controls.py`-ban |
| F05 | A $72.8\pm1.0$ feltételezett lokális $H_0$ benchmark-priorként címkézve; JSON-kulcs `H0_local_benchmark_prior_km_s_Mpc`; ábracímkék módosítva | a címkék már nem tulajdonítják az értéket a SH0ES-nek (CSV/JSON/PNG szöveg) |
| F06 | Az 1.5 fejezet mondata a $2\alpha_Y^2\beta_{\rm Yuk}$ sorról | szövegellenőrzés |
| F07 | Az RSD-jelmagyarázat a kanonikus maximumot mutatja | $<<dfsCanF>>$ a jelmagyarázatban és a JSON-ban |
| F08 | A legalsó csomópontra vonatkozó mondat törölve a 4. fejezetből, az azonosság megmarad | szövegellenőrzés |
| F09 | SNe-ábra: az alacsony $\Omega_m$-ű eltolt görbe hozzáadva, „profilozott konstans nullpont-eltolás”, $N=18$, $\nu=17$ | alacsony $\Omega_m$ eltolással $\chi^2=<<chi2LowSh>>$ |
| F10 | A sor új neve „V21 top-hat belső külső benchmark”, a peremoszlop n.a., proveniencia hozzáadva | SLy $r=2.0046409680935633$ km, $x=0.041239199439254325$; APR4 $r=2.75417748649486$ km, $x=0.041233030079658145$ |
| F11 | A PDF-építő és a `run_all.py --pdf` nem nulla kóddal lép ki pandoc/pdflatex vagy újonnan írt PDF hiányában; a `--tex-only` mód külön státuszt ad | teszt pdflatex nélküli `PATH`-szal |
| F12 | `theta_star_100_rs_only_shifted_fixed_DM` (átnevezve) és `theta_star_100_both_endpoints_shifted` (új) | alacsony $\Omega_m$: $<<thLowRsOnly>>$ és $<<thLowBoth>>$ |

## 13. Revíziós napló 3 / Revision log 3 (harmadik újraaudit, N01, N02, T01)

| Pont | Változtatás | Célzott ellenőrzés |
|:-----|:-----------------------------|:-----------------------------|
| N01 | 7. fejezet táblázata: a homogén FLRW-, a kozmikus üreg- és a Tejútrendszer-korong sor „$R\gg L$ perem” cellája n.a. (nincs peremgeometria); a `step1_6_local_yukawa_environments.csv` és a JSON az `abs_phi_surface_R_gg_L` és a `Delta_H0_apparent_surf_km_s_Mpc` mezőben null értéket, valamint új, soronkénti `boundary_status` mezőt tartalmaz; a bulk értékek változatlanok (nincsenek felezve) | ellenőrzés a `run_step26_revision_controls.py`-ban (JSON `review3_N01_N02_T01`); szöveges regressziós teszt |
| N02 | 1.7 fejezet, „Változatlanul átvihető”: a zárójel szétválasztja a háttér skalárenergia-korrekcióját ($\Omega_\phi\sim10^{-47}$) és az RSD közvetlen Yukawa-perturbációját ($\delta G_{\rm eff}/G_*\approx<<N02dG>>$, $\delta f\sigma_8\approx<<N02dfs>>$, 5. fejezet); a számok a JSON-ból származnak | szöveges regressziós teszt (HU és EN mondat) |
| T01 | Új `test_step26_text_regression.py`, amelyet a `run_all.py` a Markdown-generálás után futtat, hiba esetén nem nulla kilépési kóddal: ellenőrzi a HU/EN F06-mondatot, az F08 azonosság megfogalmazását, a felülírt hibás mondatok hiányát (V26 kétszeres növekedési törvényének tulajdonítás; az azonossági érték a táblázat legalsó csomópontjához rendelve), az N01 táblázatcelláit, az N02 mondatot és az RSD-ábra forrásának kanonikus jelmagyarázat-szövegét (F07); a `run_step26_revision_controls.py` kommentje frissítve | mutációs teszt scratch-másolaton: egy hibás mondat visszaírása nem nulla kilépési kódot ad |
"""

README_TPL = r"""# Time Attraction Law — Scalar-Tensor Yukawa Fit (V27): Step 26 Cosmological Audit / 26. lépés kozmológiai audit

**Author / Szerző:** Dávid Szabolcsi (`ORCID: 0009-0005-6300-074X`)  
**Date / Dátum:** 2026-10-28 (revised after the Academia Co-scientist review, the re-audit and the second review F01–F12 and the third re-audit N01, N02, T01 / a bírálat, az újrabírálat, az F01–F12 második bírálat és a harmadik újraaudit N01, N02, T01 után átdolgozva: 2026-10-29)

**Preceding DOIs / Előzmény DOI-k:**

- V24: `10.5281/zenodo.23680715` (strictly V24 / szigorúan V24)
- V25: `10.5281/zenodo.23699685` (strictly V25, `step1`–`step22`)
- V26: `10.5281/zenodo.23727585` (strictly V26, `step1`–`step25`)

## Reports / Jelentések

- EN: `step26_desi_bao_planck_sne_rsd_bbn_nogo_en.md` / `.pdf` (copy: `v27_cosmology_desi_planck_sne_rsd_audit_en.md`)
- HU: `step26_desi_bao_planck_sne_rsd_bbn_nogo_hu.md` / `.pdf` (copy: `v27_cosmology_desi_planck_sne_rsd_audit_hu.md`)

Computed results are generated from the JSON/CSV sources (`tal_hubble_desi_bao_summary.json`, Step 26 CSVs) by `generate_step26_markdown.py`; fixed input constants are stated in the code. / A számítási eredményeket a JSON/CSV-forrásokból generáljuk (`tal_hubble_desi_bao_summary.json`, a 26. lépés CSV-i) a `generate_step26_markdown.py` segítségével; a rögzített bemeneti állandók a kódban szerepelnek.

## Key results / Fő eredmények (branch: adiabatic attractor, $C=0$)

- Unified normalization: $\gamma_L=GL^2/c^2=<<gammaL>>$ m^3/kg, $\bar\phi_{\rm attr}(0)=<<phi0>>$, $\Omega_{\phi,\rm attr}(0)=<<Omphi0>>$, $E^2(0)=1$ (DE closure).
- DESI DR1 BAO (12×12): M4 $\chi^2=<<chi2M4>>$ ($p=<<pM4>>$); rigid 72.8 $\chi^2=<<chi2R>>$; 2-param profile $\Omega_m=<<Om2>>$, $P=<<P2>>$ km/s, $\chi^2=<<chi2P2>>$; $r_d$ needed at 72.8: $<<rdP2>>\pm<<rdP2s>>$ Mpc; $H_0$ with Planck $r_d$: $<<H0Pl>>\pm<<H0Pls>>$; joint $\chi^2=<<R4bao>>+<<R4rd>>+<<R4h0>>=<<R4joint>>$ (BAO + $r_d$ prior + $H_0$ prior), $\Delta\chi^2=<<dchi2J>>$.
- Planck acoustic scale (re-integrated $r_s$, $r_d$): M4 $100\theta_*=<<th100M4>>$.
- Late transition (piecewise $H=sH_{\rm M4}$ for $z<z_t$): see Section 3 table; the $D_M/s$ value $100\theta_*=<<R2illus>>$ is an inconsistent illustration.
- Synthetic SNe shape control (not a catalogue likelihood), symmetric offsets: M4 $<<R5offM4>>$ mag, $\chi^2$ $<<R5c0M4>>\to<<R5c1M4>>$; low-$\Omega_m$ $<<R5offLow>>$ mag, $\chi^2$ $<<R5c0Low>>\to<<R5c1Low>>$.
- RSD diagonal benchmark $\chi^2=<<chi2RSD>>$; variational $\max|\delta f\sigma_8|=<<dfsCan>>$.
- Algebraic BBN benchmark ($10^5\,{\rm D/H}^{\rm std}=2.5147$ is an arithmetic proxy); bulk local Yukawa benchmarks with `expm1`; the V21 top-hat interior external benchmark $|\phi|=0.0412$ is an external input.
- M3 control with its own horizons: $r_s=<<rsM3>>$ Mpc, $r_d=<<rdM3>>$ Mpc, $100\theta_*=<<thM3x>>$, $\chi^2_{\rm BAO}=<<chi2M3x>>$.
- $72.8\pm1.0$ km/s/Mpc is an assumed local $H_0$ benchmark prior (not a sourced SH0ES value). / A $72.8\pm1.0$ km/s/Mpc feltételezett lokális $H_0$ benchmark-prior.
- Second review F01–F12: Revision log 2 (Section 12 of the reports); open limitations P01–P09: Section 9.1; third re-audit N01, N02, T01: Revision log 3 (Section 13); text regression test: `test_step26_text_regression.py` (run by `run_all.py`). / Második bírálat: Revíziós napló 2 (12. fejezet); nyitott korlátok: 9.1 fejezet; harmadik újraaudit: Revíziós napló 3 (13. fejezet); szöveges regressziós teszt: `test_step26_text_regression.py`.

## Build / Fordítás

Dependencies / Függőségek: `pip install -r requirements.txt` (numpy, scipy, pandas, matplotlib); PDF: pandoc + pdflatex.

Recommended / Ajánlott: `python3 run_all.py` (full run including Block 1, a few seconds / az 1. blokkal együtt, néhány másodperc) or `python3 run_all.py --skip-block1` — checks the SHA-256 of the DESI input files, then runs in order / ellenőrzi a DESI-bemenetek SHA-256 értékét, majd sorban futtatja:

1. `run_step2_nogo_and_report.py` (M1–M4 branches, inverse ladder CSV and figures)
2. `run_step26_revision_controls.py` (re-audit controls R1–R6, 5.1–5.3; merges into the strict JSON)
3. `generate_step26_markdown.py` (Markdown from JSON/CSV)
   then `test_step26_text_regression.py` (HU/EN text regression test, review T01; nonzero exit on failure / szöveges regressziós teszt, hiba esetén nem nulla kilépési kód)
4. `compile_step26_clean_pdfs.py` (only with `--pdf`; pandoc + pdflatex; with `--pdf` a missing pandoc/pdflatex or a PDF that is not newly written gives a nonzero exit code; `--tex-only` writes only the `.tex` files with the separate status TEX-ONLY / `--pdf` esetén hiányzó pandoc/pdflatex vagy nem újonnan írt PDF nem nulla kilépési kódot ad; a `--tex-only` csak `.tex` fájlt ír, külön TEX-ONLY státusszal)

Without `--skip-block1`, `run_all.py` first runs `hubble_desi_bao_tal_audit.py` (Block 1: `run_step1_1.py`, `run_step1_2_and_1_3.py`, `run_step1_4_to_1_6.py`, then step 2 and the re-audit controls). / `--skip-block1` nélkül előbb a `hubble_desi_bao_tal_audit.py` fut (1. blokk, majd a 2. lépés és az újrabírálati kontrollok).

Input data / Bemeneti adatok:

<<DATA_FILES>>

## Academia Co-scientist review (round 4, 2026-10-29)

The `academia_review/` folder contains the independent reviewer's closing report on the package with ZIP SHA-256 `e08c8dfcf5454c318762a9d154ad219ad2fc0b4907b465c6dafc8f3638c58150`, reproduced verbatim:

- `step26_academia_coscientist_review4_lezarasi_jelentes_hu.md`: closing report (HU). N01, N02 and T01 are closed, the F01–F12 closures remain valid, and no new required change is listed.
- `step26_academia_coscientist_review4_ellenorzesi_jegyzokonyv.json`: the reviewer's verification record, with input hashes, the run log, the numerical comparison, independent controls and the 16 mutation tests.

These two files are the reviewer's text and are not edited by the author.

License / Licenc: CC BY-NC 4.0 International
"""
