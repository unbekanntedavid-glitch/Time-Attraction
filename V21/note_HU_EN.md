# TAL with a self-consistent stellar top-hat: conditional compact-star sensitivity

23 October 2026. Inputs: `tal_manuscript_v17.md`, Eq. (1), §3.2 and Appendix A.5; the prior SPARC true-top-hat audit and its saved project result. The input manuscript and previous audit were not modified.

## Rövid összefoglaló (HU)

A valódi, minden sugárpont köré helyezett top-hat gömbátlaggal a TAL-kitevő neutroncsillagban elérheti a százalékos időfaktor-változás tartományát; fehér törpében és a Napban nem. A 33 feltételes egyensúlyi sorozatban az 1,4 M☉-nál kapott legnagyobb sugáreltolódás **0,164 km**, miközben a három GR-EOS sugárterjedelme **2,632 km**. 2,0 M☉-nál ugyanez **0,623 km**, illetve **1,557 km**. A sugárváltozás tehát ebben az összevetésben EOS-degenerált. A Sirius B additív óracsatolása feltételesen kizárhat egy köztes, csillagméretű L-sávot, de nem ad monoton felső L-határt; az 1 m–100 km értékeket a vizsgált profilokkal nem korlátozza. Nincs detektálási állítás. A csatolás, a nagy sűrűségű EOS kauzalitása és a módosított sorozatok stabilitása külön feltevés; az 50aa2db0 konzisztenciafeladat nyitva marad.

## Summary (EN)

A true local top-hat produces percent-level clock factors in neutron stars, but not in white dwarfs or the Sun. Across 33 conditional equilibrium sequences, the largest radius shift at 1.4 M☉ is **0.164 km**, compared with a **2.632 km** three-EOS GR spread; at 2.0 M☉ the corresponding values are **0.623 km** and **1.557 km**. The radius effect is therefore EOS-degenerate in this comparison. Sirius B can conditionally exclude an intermediate, stellar-scale L band under an additive clock coupling, not impose a monotonic upper limit on L. The requested 1 m–100 km values are unconstrained by the tested WD profiles. These are sensitivity calculations, not detections, and they do not settle the pending coupling-consistency task.

## 1. Finite-source bound: where a percent-level factor is possible

Write \(x_L=\gamma_L\bar\rho_L\), \(f_L=e^{-x_L}\). For the manuscript's positive-density, Euclidean averaging ball,

\[
x_L(\mathbf x)=\frac{3G}{4\pi c^2 L}M_{\mathrm{ball}(\mathbf x,L)}.
\]

For a spherically symmetric, radially decreasing stellar density, a ball of a given radius contains most mass when centred at the stellar centre. Consequently,

\[
\sup_{\mathbf x,L}x_L(\mathbf x)\leq\frac{3}{4\pi}\max_L\frac{G M_0(<L)}{Lc^2},\qquad
0\leq1-f_L\leq1-e^{-x_{\max}}\leq x_{\max}.
\]

Here \(M_0(<L)=4\pi\int_0^L\rho(s)s^2ds\) is the **unweighted coordinate-density integral**. In the mass-source-modified model it is not the modified gravitational mass \(m(L)\). The identity involving \(\Phi_{\rm tophat}=-GM_{\rm ball}/L\) uses this monopole-at-the-ball-edge proxy; it is not an identity for the actual potential at the ball centre, nor a general-relativistic potential scalar.

Using the compactness scales specified in the request gives:

| Object | Compactness scale \(GM/(Rc^2)\) | \(x\sim(3/4\pi)C\) | \(1-e^{-x}\) |
|---|---:|---:|---:|
| Neutron star | 0.15–0.25 | 0.0358–0.0597 | **3.52–5.79%** |
| Typical white dwarf | \(10^{-4}\) | \(2.39\times10^{-5}\) | 0.00239% |
| Sun | \(2\times10^{-6}\) | \(4.77\times10^{-7}\) | 0.0000477% |

These are scale estimates, not rigorous ceilings obtained by replacing the interior maximum with surface compactness. Central concentration can increase \(\max GM_0(<L)/(Lc^2)\) relative to the surface value. The exact bound retains that maximum. **Among ordinary stars, white dwarfs and neutron stars, only the strongly compact/neutron-star regime can reach a percent-level TAL factor. White dwarfs cannot.** This statement concerns finite stellar sources, not an arbitrarily extended homogeneous cosmological background.

For a finite source, \(x_L\propto L^2\) at sufficiently small \(L\) in smooth interior matter, but \(x_L\propto1/L\) once the whole source is enclosed. An ever larger \(L\) does not create an arbitrarily large stellar correction. This is the same geometric point established in the SPARC audit.

The manuscript's \(x=0.0296\) gives \(f=0.970834\), consistent with its 0.9708 rounding. At \(L=10\) km the constants used here give \(\gamma=7.4262\times10^{-20}\,\mathrm{m^3\,kg^{-1}}\); the manuscript's \(7.41\times10^{-20}\) is an approximate value. Its fixed-density example is not a stellar-profile calculation. The quoted 367 m laboratory value is not imposed here as an established experimental exclusion: it follows only if the assumed finite-source density contrast and achieved clock sensitivity are justified.

## 2. Exactly what was solved

Let \(\rho_0\) be EOS rest-mass density, \(\epsilon\) total energy density and \(\rho=\epsilon/c^2\). **The latter is used both in TOV's mass equation and in the TAL averaging integral.** This convention recovers GR at \(f=1\). The manuscript does not unambiguously distinguish baryon density from mass-equivalent energy density; using \(\rho_0\) instead would be another model, not a negligible notational substitution.

Define

\[
A(r)=\frac{G[m(r)+4\pi r^3P(r)/c^2]}{c^2r^2[1-2Gm(r)/(rc^2)]}.
\]

Three equilibrium prescriptions were evaluated:

| Label in files | Mass equation | Pressure/lapse-gradient prescription |
|---|---|---|
| `GR` | \(m'=4\pi r^2\epsilon/c^2\) | \(P'=-(\epsilon+P)A\), \(\nu'=A\) |
| `mass` — published mass-source insertion | \(m'=4\pi r^2(\epsilon/c^2)f_L\) | Retain \(P'=-(\epsilon+P)A\), \(\nu'=A\) |
| `gradient` — labelled alternative | \(m'=4\pi r^2\epsilon/c^2\) | \(P'=-(\epsilon+P)f_LA\), \(\nu'=f_LA\) |

The metric lapse notation is \(g_{tt}=-e^{2\nu}c^2\). The alternative multiplies the **gradient**, not the lapse itself. Multiplying \(e^\nu\) by \(f_L\) would introduce \(d\ln f_L/dr\) and is not the same prescription. No model simultaneously applies both modifications.

**These are conditional ansätze, not solutions of a completed TAL field theory.** Backlog consistency task **50aa2db0** may change which insertion, if any, is admissible. This note does not resolve that task, derive a covariant action, or supply missing effective stresses. In particular, an extremum of these modified equilibrium sequences is not a demonstrated radial-stability boundary. `Mmax` below means the maximum mass on the computed sequence; the GR stability interpretation must not automatically be transferred to the modified models.

### True top-hat, centred locally

The density is zero outside the stellar surface. At each stellar radius \(r\), the averaging sphere is centred at that point, not at the stellar centre. For \(r>0\), full shells \(0<s<L-r\), when present, contribute \(4\pi\int\rho s^2ds\). Partial shells contribute

\[
M_{\rm cap}(r,L)=\frac{\pi}{r}\int_{|r-L|}^{\min(r+L,R)}s\rho(s)\,[L^2-(s-r)^2]ds,
\]

with empty intervals omitted. At \(r=0\), \(M=4\pi\int_0^{\min(L,R)}\rho s^2ds\). This is divided by the **whole** ball volume \(4\pi L^3/3\), including any vacuum portion. Neither central density nor \(m(r)/(4\pi r^3/3)\) is substituted for the local ball average.

Integration uses logarithmic relativistic enthalpy \(h=\int dP/(\epsilon+P)\), a central expansion at \(r=10^{-4}\) km, and an outer endpoint \(h=10^{-10}\). Each solution is iterated: solve the equilibrium ODEs, rebuild the full density profile, recompute the top-hat, and repeat until the relative mass/radius change and maximum exponent change are both below \(2\times10^{-7}\). Smooth interpolation is used for the exponent feedback. Production runs use ODE relative tolerance \(10^{-10}\), 801 radial profile points, and 80-point Gauss shell quadrature. The convolution is the manuscript's Euclidean coordinate-volume prescription in areal coordinates, **not** a new proper-volume/geodesic-ball definition. At neutron-star compactness that choice is physically consequential and remains conditional.

### EOSs and scan

The three standard EOSs are represented by the **Read et al. (2009) piecewise-polytropic fits**, not by newly downloaded microscopic EOS tables. These were explicitly allowed in the request. Pressure and energy density are matched continuously, including the four-piece SLy crust.

| EOS | \(\log_{10}P_1\), dyne cm\(^{-2}\) | \(\Gamma_1\) | \(\Gamma_2\) | \(\Gamma_3\) |
|---|---:|---:|---:|---:|
| SLy | **34.384** | 3.005 | 2.988 | 2.851 |
| APR4 (AP4 in Read's table) | 34.269 | 2.830 | 3.445 | 3.348 |
| H4, stiff-radius comparison | 34.669 | 2.909 | 2.246 | 2.144 |

Core boundaries are \(10^{14.7}\) and \(10^{15}\) g cm\(^{-3}\). The low-density coefficients printed in Read's Table II are interpreted as \(K/c^2\) in cgs mass-density units, consistent with their numerical magnitude and the GR benchmarks; pressure is restored in physical units. The later three crust coefficients are fixed by continuity, avoiding small discontinuities from rounded table entries.

For every EOS, GR and both modifications were run for \(L=1,367,1000,10000,100000\) m: **33 sequences in total**. The initial 53-point central-rest-density grid spans \(\log_{10}\rho_{0c}[\mathrm{g\,cm^{-3}}]=14.55\)–15.65. A local bounded search refines each maximum. Radii at fixed masses are found by root solving on the lower-central-density branch, rather than reading off a coarse plotting grid. A blank fixed-mass entry means the sequence cannot reach that target on the selected branch; it is not extrapolated.

Some standard fits become superluminal at high density. `max_cs2_over_c2` explicitly flags this in the summary table. No causal completion was imposed, so especially the high-density APR4 maximum is a **benchmark of the published fit**, not an admitted causal-EOS prediction.

## 3. Numerical results and EOS degeneracy

### GR baseline

| EOS | Maximum sequence mass, M☉ | R(1.4 M☉), km | R(2.0 M☉), km |
|---|---:|---:|---:|
| SLy | 2.048411 | 11.7061 | 10.6266 |
| APR4 | 2.188001 | 11.3197 | 10.8667 |
| H4 | 2.013851 | 13.9522 | 12.1834 |

### Full EOS × L comparison

Masses are in M☉ and radii in km. The two radius columns are at **fixed gravitational mass 1.4 M☉**. All digits are numerical outputs, not observational precision; differences at 1 m are below the tested numerical error floor.

| EOS | L, m | Mmax: mass-source | Mmax: gradient | R1.4: mass-source | R1.4: gradient |
|---|---:|---:|---:|---:|---:|
| SLy | 1 | 2.048411 | 2.048411 | 11.7061 | 11.7061 |
| SLy | 367 | 2.048591 | 2.048567 | 11.7065 | 11.7063 |
| SLy | 1000 | 2.049739 | 2.049558 | 11.7091 | 11.7074 |
| SLy | 10000 | 2.098219 | 2.108482 | 11.8703 | 11.7920 |
| SLy | 100000 | 2.053290 | 2.058655 | 11.7356 | 11.7297 |
| APR4 | 1 | 2.188001 | 2.188001 | 11.3197 | 11.3197 |
| APR4 | 367 | 2.188182 | 2.188160 | 11.3201 | 11.3199 |
| APR4 | 1000 | 2.189339 | 2.189168 | 11.3227 | 11.3209 |
| APR4 | 10000 | 2.241521 | 2.248766 | 11.4757 | 11.3937 |
| APR4 | 100000 | 2.193438 | 2.198771 | 11.3471 | 11.3393 |
| H4 | 1 | 2.013851 | 2.013851 | 13.9522 | 13.9522 |
| H4 | 367 | 2.013994 | 2.013975 | 13.9525 | 13.9523 |
| H4 | 1000 | 2.014904 | 2.014767 | 13.9545 | 13.9532 |
| H4 | 10000 | 2.059898 | 2.069632 | 14.0952 | 14.0267 |
| H4 | 100000 | 2.018823 | 2.025597 | 13.9875 | 13.9821 |

All sequences pass the **nominal 2 M☉ mass-capability check**, before imposing EOS causality or modified radial stability. This is not evidence for TAL: the GR baselines already pass it. The mass-source insertion reduces the source at a fixed density profile, but the star readjusts; its maximum equilibrium mass need not decrease.

Only APR4's GR sequence reaches the **central** J0740 mass 2.073 M☉. SLy's 10 km modified sequences also reach it; H4's maximum is 2.0599 M☉ for the mass-source ansatz and 2.0696 M☉ for the gradient ansatz, just below that central value. The measured mass has a 0.069 M☉ posterior uncertainty, so these statements are not individual EOS exclusions. The APR4 fit has maximum sound-speed squared about **1.38–1.40 c²** near its mass peak. Several mass-source SLy peaks are slightly superluminal, up to about **1.011 c²**. H4 remains below c² in these peak models. These flags prevent interpreting the uncompleted fit maxima as universal physical mass limits.

### Does the radius shift exceed EOS uncertainty?

| Fixed M, M☉ | GR EOS spread, km | Largest mass-source ΔR, km | Largest gradient ΔR, km |
|---|---:|---:|---:|
| 1.4 | 2.6325 | 0.1642 | 0.0859 |
| 2.0 | 1.5569 | 0.6234 | 0.5475 |

The largest shifts occur at L=10 km in the requested set. At 1.4 M☉, L=367 m gives only **0.306–0.412 m** of mass-source radius increase and **0.139–0.178 m** for the gradient ansatz. L=1 km gives about 2.26–3.04 m and 1.03–1.32 m, respectively. L=100 km is weaker again: approximately 27.4–35.3 m and 19.6–29.9 m.

At 2.073 M☉ there is **no meaningful three-EOS GR spread**: only one GR EOS reaches that mass. `radius_shifts.csv` therefore leaves that spread blank instead of reporting a misleading zero. The APR4 10 km mass-source radius shift there is about 0.325 km, but this cannot be promoted into a multi-EOS comparison at that mass.

The profile diagnostics demonstrate the difference from a central-density proxy. At 1.4 M☉ and L=10 km, the true centre exponents are **0.0321–0.0469**, while extending the central density uniformly across the averaging ball would give **0.0458–0.0738**. The surface exponents are only **0.00616–0.01181**. The factors vary through the star and are recomputed on the modified solution.

The chosen three-EOS radius range is a comparison envelope, not a probability interval or the full set of possible EOSs. A shift smaller than this envelope is not a proof of exact functional degeneracy at every mass; it means radius alone does not identify TAL in this comparison. No EOS parameters or \(L\) were fit to NICER.

### Numerical checks

Every EOS/coupling/L case was re-solved at its peak central density and near 1.4 M☉ with **1601 radial points, 160 quadrature points and ODE rtol=10⁻¹¹**: 66 refinement comparisons. The largest differences were **6.43×10⁻⁷ M☉** in mass and **3.92×10⁻⁶ km** in radius; the estimated fixed-mass radius difference near 1.4 M☉ was at most **4.10×10⁻⁶ km**. The refined exponent fixed-point residual was at most **4.41×10⁻⁸**. These are tested changes, not rigorous global error bounds. The tiny 1 m radius shifts are numerically unresolved and their displayed signs have no physical meaning.

A separate, polynomial-exact integral on the piecewise-linear profile agrees with the 80-point production convolution to within **9.0×10⁻⁵ relative** in the 15 checked geometry cases, including all five L values at centre, half radius and surface. GR maxima and R1.4 agree with the Read table's quoted *fit* predictions to within **0.051%**; the small common offset is retained, not renormalized away.

Tests also include exact uniform-sphere limits, the surface overlap \(\bar\rho_L/\rho=5/16\) for \(L=R\), full enclosure, EOS enthalpy inversion, and a separate piecewise-polynomial-exact shell integration of selected profiles. Diagnostics retain errors rather than rounding them away. Numerical precision is not EOS accuracy; the Read fits themselves differ from the microscopic tabulations.

## 4. NICER reference values: verified, but not a TAL likelihood

Primary text was read for all values below. Uncertainties are marginal 16th–84th-percentile credible intervals, approximately 68%.

| Pulsar and selected analysis | Mass, \(M_\odot\) | Equatorial circumferential radius, km |
|---|---:|---:|
| J0740+6620, Salmi et al. 2024 | \(2.073\pm0.069\) | \(12.49^{+1.28}_{-0.88}\) |
| J0437−4715, Choudhury et al. 2024 | \(1.418\pm0.037\) | \(11.36^{+0.95}_{-0.63}\) |
| J0030+0451, **Kini et al. 2026**, preferred PDT-U | \(1.43^{+0.20}_{-0.17}\) | \(12.68^{+1.31}_{-1.04}\) |
| J0030+0451, same 2026 data, alternative ST+PDT | \(1.29^{+0.13}_{-0.12}\) | \(11.63^{+0.95}_{-0.89}\) |

The J0030 update uses six years of NICER observations and is not the older 2024 \(1.70\,M_\odot,14.44\) km preferred-model result. The verified 2026 primary source is cited as arXiv:2602.23743; its retrieved arXiv record labels it submitted to ApJ.

The requested 2024 J0437 posterior mass error is **0.037**, not the radio-prior error **0.044**. An important current sensitivity is Miller et al.'s 2025 reanalysis including modulated nonthermal emission: its radius interval is approximately **11.8–15.1 km**, with a bimodal radius posterior. Its Table 5 gives median 13.892 km and bounds 11.823 and 15.142 km. This is not silently combined with, or substituted for, the specifically requested Choudhury result. It further weakens any claim that a sub-km theory shift is identifiable from the 2024 result alone.

**Figure prescription:** each plotted reference is an *uncorrelated, symmetrized Gaussian approximation*. The radius/mass standard deviations are the arithmetic mean of the published lower and upper errors. Ellipses satisfy \((\Delta R/\sigma_R)^2+(\Delta M/\sigma_M)^2=-2\ln(1-p)\), for **joint** \(p=0.68,0.95\). They are not 1D error bars relabelled as joint contours, and they are not the papers' posterior samples. Actual mass–radius correlation, skewness and hot-spot multimodality are absent. The figure labels this explicitly.

All NICER posteriors were inferred with GR-based photon propagation, atmosphere models and rotational geometry. Their comparison with nonrotating conditional TAL sequences is an **illustrative overlay**, not a self-consistent alternative-gravity likelihood. Rotation changes equatorial radius and maximum mass; neither rotating stellar models nor TAL pulse-profile inference were done. The horizontal \(2\,M_\odot\) line is a conventional necessary mass-capability check, not a sharp observational cutoff or proof of stability.

## 5. Sirius B: additive-clock redshift and a finite excluded band

Joyce et al. (2018), abstract, §3.1.7 and Table 3, give the preferred E1 result

\[
v_{\rm gr}=80.65\pm0.77\ \mathrm{km\,s^{-1}},\quad
M_{\rm redshift}=1.017\pm0.025\,M_\odot,
\quad R=(0.00803\pm0.00011)R_\odot.
\]

Bond et al. (2017), primary abstract, give the independent dynamical mass \(M_B=1.018\pm0.011\,M_\odot\). The mass agreement must not be tested against a second mass inferred from the same redshift as if it were independent. Instead use the dynamical mass and independently estimated radius to predict GR redshift.

The **additional clock-only assumption** here is

\[
1+z=(1+z_{\rm GR})\exp(x_{\rm em}-x_{\rm obs}),\qquad
\delta v_B=c(1+z_{\rm GR,B})[e^{x_B}-1].
\]

This retains a GR exterior and multiplies its local clock factor by TAL. It is **not uniquely implied** by either equilibrium insertion in §2. In particular, the pressure-gradient prescription does not itself imply this multiplicative lapse law. The WD test therefore remains separately conditional on the clock coupling. The differential Sirius A/B observing strategy cancels the observer factor; the modelled additional Sirius A contribution is subtracted. No atmospheric line-transfer model or refitted orbital dynamics is supplied.

Using R☉=6.957×10⁸ m gives **R_B=5.586471×10⁶ m** and compactness **2.69079×10⁻⁴**. The exact exterior Schwarzschild expression gives **v_GR=80.70055 km/s**. With independent mass/radius uncertainties, σ_GR=1.4080 km/s and the measured-minus-GR residual is **−0.05055±1.60481 km/s** after adding the redshift measurement error. The approximate two-sided-95% positive extra-shift allowance is therefore **3.0949 km/s**. Using the paper's larger quoted redshift-mass uncertainty instead gives a more conservative allowance **4.1645 km/s**. Both are approximate Gaussian consistency allowances, not a covariance-aware posterior on L; mass and radius share astrometric inputs.

A particularly simple result does **not** need a WD density model. For a surface emitter and **L≥2R_B=1.1172942×10⁷ m**, the local ball contains the entire WD:

\[
x_B=\frac{3GM_B}{4\pi Lc^2},\qquad
\delta v_B=c(1+z_{\rm GR,B})\left[\exp\!\left(\frac{3GM_B}{4\pi Lc^2}\right)-1\right].
\]

At L=2R_B this adds approximately **9.63 km/s**, larger than either allowance. Holding the GR exterior/dynamical interpretation fixed and neglecting the negligible Sirius A correction over this interval, the all-enclosed conditional exclusion is

- **1.117×10⁷ ≲ L ≲ 3.477×10⁷ m** with direct independent-error propagation;
- **1.117×10⁷ ≲ L ≲ 2.584×10⁷ m** with the more conservative quoted-mass allowance.

These are a profile-independent *subset* of the possible excluded band under this specific additive-clock model, not global bounds on L. Below 2R_B the result depends on the density profile; above the upper endpoint the signal has fallen below the allowance.

For the illustrative WD profiles, the extra surface redshift at the requested scales is:

| L, m | δv_B, km/s, n=1.5 | δv_B, km/s, n=3 |
|---|---:|---:|
| 1 | 2.4×10⁻²³ | 3.1×10⁻³⁴ |
| 367 | 2.3×10⁻¹⁴ | 2.0×10⁻²¹ |
| 1000 | 7.6×10⁻¹³ | 3.0×10⁻¹⁹ |
| 10000 | 2.4×10⁻⁹ | 3.0×10⁻¹⁴ |
| 100000 | 7.6×10⁻⁶ | 3.1×10⁻⁹ |

The very small-L entries are idealized zero-density-surface numbers, not atmosphere-resolved predictions. At 100 km even the larger illustrative shift is only **0.0076 m/s**, compared with the measured **770 m/s** redshift error. The sampled peak extra shifts are about **11.56 km/s** at L≈8.25×10⁶ m for n=1.5 and **13.57 km/s** at L≈7.26×10⁶ m for n=3. These redshift velocities are small dimensionless clock corrections, of order 4×10⁻⁵, not percent-level WD factors.

For illustrative \(L<2R_B\) sensitivity only, normalized Lane–Emden WD profiles with \(n=1.5\) and \(n=3\) were used at the observed mass/radius, with an \(n=3\) Sirius A profile. These are **bracketing illustrative polytropes, not a tabulated finite-temperature Sirius B EOS or guaranteed bounds on its profile**. Very small-\(L\) photospheric estimates are especially atmosphere-dependent because ideal polytropes terminate at zero density. The all-enclosed band above does not have that dependence. Sirius A's additional correction is negligible over that band; its physical parameters were retained from Joyce's discussion for the differential-reference calculation.

The sampled WD curves do **not** constrain any of the requested \(L\le100\) km values at the reported precision in these profiles. Once the whole WD is enclosed, the extra shift decreases as \(1/L\), so large \(L\) becomes allowed again. **There is no monotonic upper limit on \(L\) from this finite-source argument.** A different redshift/dynamical coupling, correlated astrometric errors, atmospheric systematics or a replacement rather than additive law could change the inferred band. This is a conditional consistency check, not a TAL detection or a final universal exclusion.

## 6. Reproduction and interpretation limits

Extract `tal_compact_stars_code_results.zip` into `user/workspace/co-scientist-outputs/tal-compact-stars/` to restore the code and supporting numerical files, then run from the workspace root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 user/workspace/co-scientist-outputs/tal-compact-stars/stellar_tal.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 user/workspace/co-scientist-outputs/tal-compact-stars/diagnostics.py
```

Requires Python, NumPy, SciPy, pandas and matplotlib; exact versions are in `numerical_diagnostics.json`. No network, external EOS library or fitted TAL parameter is needed. The main script regenerates the curves, tables and figures. The diagnostic script adds detailed convergence tests, EOS causal flags and profile checks. The ZIP contains code and all numerical products; loose summary/plot files are supplied for quick inspection.

The following choices are not hidden uncertainty bars: which density enters TAL; Euclidean versus proper-volume averaging; the mass-source versus gradient coupling; causal EOS completion; modified radial stability; and the mapping from the model to GR-derived NICER or spectroscopic inferences. **Task 50aa2db0 remains unresolved. No detection claim or unconditional bound on \(L\) is made.**

## 7. V21 addendum: CSV-driven NICER overlay, causal onset and turning-point heuristic

Prepared for Zenodo V21 (GitHub V19), 23 October 2026. Sections 1–6 above are unchanged. Everything below was computed by the scripts in this package. The numbers come from `stellar_summary.csv` and `causality_stability_log.txt`.

### 7.1 Rövid összefoglaló (HU)

**NICER-bemenetek.** A `plots()` függvény a NICER-referenciaértékeket most az `observational_inputs.csv` fájlból olvassa. A kódban nem maradt beégetett NICER-szám. A CSV-ben szerepelnek a szimmetrizált hibák (az alsó és felső hiba számtani közepe), a 68/95%-os ellipszisek valószínűségei és átlátszóságai, a címkék és a színek. Benne van az RA/Dec (J2000) is, az ATNF Pulsar Catalogue v2.6.5 alapján. Az auditjavítások előtti CSV-alapú rajzoló pixelre és bájtra azonos ábrát adott az eredeti `mr_nicer.png`-vel azonos adatokból. A jelenlegi ábra az új stílusok, címkék és végpontjelölők miatt szándékosan eltér (7.4).

**Hangsebesség.** A program a \(c_s^2=dP/d\epsilon\) értéket minden elfogadott ODE-lépésben, minden kimeneti csomópontban és minden EOS-szakaszhatár két oldalán kiszámolja. Az első szuperluminális sűrűség (nyugalmi tömegsűrűség \(\rho_0\), illetve \(\epsilon/c^2\)):

- **APR4:** \(\rho_0=1{,}4458\times10^{15}\) g cm⁻³ (\(\log_{10}=15{,}1601\)), \(\epsilon/c^2=1{,}8013\times10^{15}\) g cm⁻³.
- **SLy:** \(\rho_0=2{,}0071\times10^{15}\) g cm⁻³ (\(\log_{10}=15{,}3026\)), \(\epsilon/c^2=2{,}8626\times10^{15}\) g cm⁻³.
- **H4:** \(\rho_0=1{,}0461\times10^{16}\) g cm⁻³ (\(\log_{10}=16{,}0196\)). Ez minden számolt csillag központi sűrűsége fölött van.

**APR4-nél \(\rho_0\approx1{,}446\times10^{15}\) g cm⁻³ (\(\epsilon/c^2\approx1{,}801\times10^{15}\) g cm⁻³) felett a GR alap is szuperluminális, a módosított ágak ezt öröklik. A kauzális tartomány végpontjait külön közöljük, a teljes egyensúlyi sorozatokat pedig kifejezetten megjelölt, potenciálisan akauzális benchmarkként megtartjuk.** Az APR4 mind a 11 maximumcsillagában a kauzalitássértés a központi sűrűség 74,5–76,3%-ánál kezdődik, és egy 5,244–5,411 km sugarú magot érint (a \(c_s^2=c^2\) határ a sűrű ODE-megoldásból, 7.4). A kauzális maximális tömeg GR-ben 2,0781 M☉ (a teljes sorozatmaximum 2,1880 M☉). A módosított ágakban 2,0781–2,1414 M☉.

Az SLy GR-maximuma épp kauzális (\(c_s^2/c^2=0{,}99976\)). A tömegforrás-ansatz 1 km-es, 10 km-es és 100 km-es maximuma enyhén szuperluminális (legfeljebb 1,0112). Ezekben a sértés a központi sűrűség 98,7–99,9%-ánál kezdődik, és 0,230–1,067 km sugarú magra terjed ki (a \(c_s^2=c^2\) határ; a korábbi „legfeljebb 1,06 km” mintavételi alsó becslésen alapult, lásd 7.4). A kauzális maximális tömeg legfeljebb 0,00015 M☉-mel kisebb a sorozatmaximumnál. A gradiens-ansatz SLy-maximumai és mind a 11 H4-maximum kauzális (H4: \(c_s^2/c^2\le0{,}522\)).

Mind a 33 sorozatban az 1,4 és a 2,0 M☉-es összehasonlító pont a kauzális tartományban van, így az R(1,4) és R(2,0) értékek (és a fenti 0,164/2,632 km, illetve 0,623/1,557 km) kauzális csillagokra vonatkoznak.

**Stabilitás (heurisztika).** A fordulópont-kritériumot (\(dM/d\rho_c>0\)) mind a 33 sorozatra alkalmaztuk. Az 53 pontos rácson minden sorozatban pontosan egy előjelváltás van, a tömegmaximumnál. A stabil ág \(\log_{10}\rho_{0c}=14{,}55\)-től (a rács alja) a maximumig tart, és az 1,4 és a 2,0 M☉-es pont mindenhol ezen az ágon van. A finomító pontok közötti további előjelváltások mind a csúcs ±2×10⁻⁴ dex-es környezetében vannak, és ≤2,1×10⁻⁷ M☉ nagyságúak. Ez a numerikus finomítási szint (6,4×10⁻⁷ M☉) alatt van. **A módosított egyenletekre ez csak heurisztika, nem bizonyított stabilitási eredmény.** A teljes radiális módusanalízis nyitva marad, és az 50aa2db0 feladat is nyitva van.

### 7.2 Summary (EN)

**NICER inputs.** `plots()` now starts with `pd.read_csv(OUT/'observational_inputs.csv')`, resolved relative to the script directory, and draws every reference ellipse, marker, label, colour and the colour caption from that file. No NICER number is hard-coded any more. The CSV keeps the published asymmetric errors and adds the following:

- symmetrized errors (arithmetic mean of lower and upper);
- the joint-probability levels (0.68, 0.95) and fill alphas of the two ellipses;
- labels and colours;
- J2000 RA/Dec with position epoch, taken from the ATNF Pulsar Catalogue v2.6.5.

Rows with `used_in_plot=False` (the J0030 ST+PDT alternative and the Miller et al. J0437 sensitivity) are kept for reference and are not drawn. Before the audit fixes, with identical sequence data, the CSV-driven `--plots-only` output was byte-identical to the original `mr_nicer.png` (SHA-256 `1f2e4e84…a51e`, 0 of 4,924,800 pixels differ; the PDF differed only in its `/CreationDate` stamp). The present figure differs intentionally (new styles, `plot_label` legend, endpoint markers; §7.4).

**Sound speed.** `integrate()` evaluates \(c_s^2/c^2=dP/d\epsilon\) (\(=\Gamma P/(\epsilon+P)\) on each polytropic piece) at every accepted ODE step, every one of the 2401 output nodes and both sides of every EOS piece boundary. Each EOS has a single acausal interval, beginning at the following density:

| EOS | onset \(\rho_0\), g cm⁻³ | \(\log_{10}\rho_0\) | onset \(\epsilon/c^2\), g cm⁻³ | onset as % of \(\rho_{0c}\) at Mmax | max causal mass, M☉ |
|---|---:|---:|---:|---:|---:|
| APR4 | 1.4458×10¹⁵ | 15.1601 | 1.8013×10¹⁵ | 74.49–76.28 (all 11 sequences) | GR 2.0781; modified 2.0781–2.1414 |
| SLy | 2.0071×10¹⁵ | 15.3026 | 2.8626×10¹⁵ | 98.69–99.94 (mass-source L=1, 10, 100 km only) | at most 1.5×10⁻⁴ below Mmax |
| H4 | 1.0461×10¹⁶ | 16.0196 | 4.3012×10¹⁶ | never reached (all Mmax stars causal) | = Mmax |

**For APR4, above \(\rho_0\approx1.446\times10^{15}\) g cm⁻³ (\(\epsilon/c^2\approx1.801\times10^{15}\) g cm⁻³) the GR baseline is already superluminal. The modified branches inherit this. We report causal-range endpoints separately while retaining the full equilibrium sequences as explicitly labelled, potentially acausal benchmarks.** All 11 APR4 maximum-mass stars contain an acausal core whose \(c_s^2=c^2\) boundary lies at r = 5.244–5.411 km, with \(c_s^2/c^2\) up to 1.370–1.402. This refines the rounded 1.38–1.40 quoted in §3.

The SLy GR maximum is just causal (0.99976). The mass-source maxima at L=1, 10 and 100 km are slightly acausal (up to 1.0112), with acausal-core boundaries at r = 0.230–1.067 km (1.067 km at L=10 km). The earlier "≤1.06 km" was based on the outermost sampled acausal point, a lower estimate; see §7.4. The gradient-ansatz SLy maxima and all H4 maxima are causal (H4 ≤0.522).

In all 33 sequences the 1.4 and 2.0 M☉ comparison stars have central densities below the onset, so every R(1.4) and R(2.0) in this package, and the 0.164/2.632 km and 0.623/1.557 km comparison, lies in the causal range. New summary columns:

- `causal_violation_flag`, `cs2_max_over_c2_in_Mmax_star`;
- `cs2_onset_log10_rho0_cgs`, `cs2_onset_log10_eps_over_c2_cgs`, `cs2_onset_percent_of_rho0c_at_Mmax`, `r_acausal_core_at_Mmax_km` (core boundary), `r_acausal_core_sampled_lower_at_Mmax_km` (sampled lower estimate);
- `M_max_causal_Msun`, `R_at_M_max_causal_km`, `log10_rho0c_at_M_max_causal`;
- `log10_rho0c_at_*`, `causal_at_1p4`, `causal_at_2p0`.

`mr_sequences.csv` gains the per-star columns `cs2_max_over_c2_all_steps`, `acausal_point`, `r_acausal_core_km` (core boundary) and `r_acausal_core_sampled_lower_km`. No causal completion of the EOS was imposed. "Max causal mass" is the mass of the star whose centre sits exactly at the onset density on that sequence.

**Turning-point heuristic.** \(dM/d\rho_c>0\) was evaluated on the 53-point scan grid of all 33 sequences:

- Every sequence has exactly one sign change, at its mass maximum.
- The stable branch runs from the grid floor \(\log_{10}\rho_{0c}=14.55\) (M = 0.23 M☉ APR4, 0.29 M☉ SLy, 0.70 M☉ H4; lower densities were not computed) to the refined maximum.
- The 1.4 and 2.0 M☉ points are on that branch in all 33 sequences.

When the optimizer and root-finder points are included, additional sign flips appear. All of them lie within ±2×10⁻⁴ dex of the peak, with reversal steps ≤2.1×10⁻⁷ M☉, below the 6.4×10⁻⁷ M☉ refinement floor of §3. They are numerical noise on a flat maximum, not extra turning points. The peak location itself is therefore resolved only to about 10⁻⁴ dex. **For the modified (mass-source and gradient) equations this criterion is a GR-motivated heuristic, not a proven radial-stability result.** The full radial-mode (Sturm–Liouville) analysis with the TAL-modified perturbation equations remains open.

### 7.3 Reproduction in this V21 run and open items

- The unchanged original scripts were rerun (per-EOS `--eos` runs, then `--combine` and `diagnostics.py`) with Python 3.13.5, NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.0 and matplotlib 3.11.2. The original run used Python 3.12.3 and pandas 3.0.6.
  - **Not byte-identical, but numerically reproduced.** Mmax relative change ≤1.7×10⁻⁷. Fixed-mass radii ≤1.2×10⁻⁷ relative (≤1.4×10⁻⁶ km). Radius shifts ≤1.8×10⁻⁶ km.
  - The peak central density moved by ≤1.4×10⁻⁴ dex on the flat maximum, so the radius at Mmax moved by ≤9.4×10⁻⁴ km.
  - `mr_sequences.csv` has 3113 rather than 3103 rows, because the bounded optimizer took a few more evaluations.
  - All headline numbers in §§3–5 are unchanged at their quoted precision.
  - `sirius_redshift.png` and `geometry_tests.json` were byte-identical. `wd_summary.json` agrees to ≤5×10⁻¹⁶ relative and `sirius_redshift.csv` to ≤6.3×10⁻⁹ relative.
- The V21 code changes do not alter any previously reported value: all pre-existing columns are identical to the unchanged-code rerun in the same environment.
- Still open:
  - task **50aa2db0** (coupling consistency);
  - a full radial-mode stability analysis of the modified equations;
  - causal completion of the APR4 (and marginally SLy) high-density fits;
  - densities below \(\log_{10}\rho_{0c}=14.55\), which are not scanned.
- No detection claim or unconditional L bound is made.

**Additional reference (V21).** Manchester, R. N., Hobbs, G. B., Teoh, A., & Hobbs, M. (2005), AJ 129, 1993, *The ATNF Pulsar Catalogue*; catalogue v2.6.5, https://www.atnf.csiro.au/research/pulsar/psrcat/ (RAJ/DECJ queried 23 October 2026). SIMBAD (https://simbad.cds.unistra.fr) was queried as a cross-check; its listed positions agree to within the proper-motion/epoch differences.

### 7.4 Audit fixes (Academia audit of this V21 package, 23 October 2026)

**Rövid összefoglaló (HU).** Az audit négy hibáját (H1–H4) és a lefedettségi hiányt javítottuk; a csomag továbbra is „prepared for Zenodo V21 (GitHub V19)”.
- **H1:** az ábrán a TAL-görbe folytonos, a GR-görbe szaggatott a kauzális végpontig; a sorozatmaximumig tartó folytatás pontozott, halványított, potenciálisan akauzális benchmark. Üres kör: TAL kauzális végpont, üres négyzet: GR kauzális végpont. Érintett: APR4 11/11, SLy 3/11 sorozat. A 7.2 mondatát átírtuk (fent). A csomagban tárolt adatokon a meghatározott számlálási szabály 182 egyedi akauzális táblázatsort ad (173 APR4, 9 SLy): az `acausal_point=True` sorok közül azokat számoljuk, amelyek `log10_rho0c_cgs` értéke nem nagyobb a megfelelő EOS–csatolás–L sorozat `stellar_summary.csv` fájlban tárolt maximumának központi sűrűségénél. Az összehasonlítás a CSV-mezők pontos decimális értékén történik, hozzáadott tolerancia vagy újabb kerekítés nélkül; a különböző paneleken ismétlődő előfordulásokat nem számoljuk újra.
- **H2:** a magsugár most a \(c_s^2/c^2=1\) határ területi sugara: a küszöbentalpiát közrefogott gyökkereséssel (brentq) határozzuk meg, a sugarat a végső ODE-megoldás sűrű interpolációjából vesszük. A régi mintavételi érték külön oszlopban, alsó becslésként maradt. APR4: 5,244–5,411 km; SLy: 0,230–1,067 km. A „≤” jelölést töröltük.
- **H3:** a `plot_label` a NICER-jelmagyarázatban és a színmagyarázatban szerepel; módosítása megváltoztatja az ábrát. Új ellenőrzés: a szimmetrizált hiba = (alsó+felső)/2, eltérésnél hibával leáll; mind a 10 forrás/mennyiség sor PASS.
- **H4:** a README-ben öt végrehajtható parancssor szerepel.
- **Lefedettség:** a `diagnostics.py` mind a 33 kauzális végpontot és mind a 14 akauzális maximum maghatárát finomított beállítással újraszámolja (`causal_endpoint_convergence.csv`).

**Summary (EN).**

- **H1 (figure vs. causal-range wording).** In `mr_nicer.png/pdf` each TAL curve is solid and each GR curve dashed up to its causal endpoint (the star whose centre sits at the onset density). The continuation to the sequence maximum is dotted and faded and is labelled as a potentially acausal benchmark of the published fit; it contains the stored acausal points selected for the drawn branches. The archive point count is **182 unique table rows: 173 APR4 and 9 SLy**, not the number of repeated point occurrences across figure panels. The counting rule selects `acausal_point=True` rows in `mr_sequences.csv` whose `log10_rho0c_cgs` is less than or equal to the peak central density stored in `stellar_summary.csv` for the same EOS, coupling and L. This count compares the CSV fields as exact decimal values, with no additional rounding or inclusion tolerance. It is an archive bookkeeping count, not a tolerance-expanded selection of near-peak points. Open circles mark TAL causal endpoints and open squares GR causal endpoints; markers appear only on sequences that cross the onset (APR4: 11 of 11 sequences; SLy: 3 of 11; H4: none). The caption states this. The §7.2 sentence now reads: "We report causal-range endpoints separately while retaining the full equilibrium sequences as explicitly labelled, potentially acausal benchmarks." The HU text was changed consistently.
- **H2 (core radius).** `EOS.h_causal()` finds the threshold enthalpy \(h_*\) with \(c_s^2/c^2=1\) by a bracketed root (brentq, xtol=rtol=10⁻¹⁵; residual \(|c_s^2/c^2-1|\le2.2\times10^{-16}\)). `core_boundary()` evaluates the areal radius at \(h_*\) from the dense output of the converged ODE solution (for \(h_*\) inside the analytic central step, the series \(r_0\sqrt{(h_c-h_*)/\Delta h}\) is used). This is applied to every acausal star in `mr_sequences.csv` (654 points) and every acausal maximum in `stellar_summary.csv` and the log. The old sampled values are kept only in the new columns `r_acausal_core_sampled_lower_km` / `r_acausal_core_sampled_lower_at_Mmax_km`, labelled as lower estimates; the boundary exceeds them by up to 10.68 m at Mmax (17.6 m over all sequence points). Results at Mmax: APR4 5.244–5.411 km (11 stars), SLy mass-source 0.237, 1.067 and 0.230 km for L=1, 10, 100 km. Check against the audit's dense-solution values:

| Configuration | Audit (dense solution) | This package (boundary) | Sampled lower estimate |
|---|---:|---:|---:|
| APR4, GR | 5.246704 km | 5.246704 km | 5.246477 km |
| APR4, mass-source, L=10 km | 5.411172 km | 5.411172 km | 5.409369 km |
| SLy, mass-source, L=10 km | 1.067392 km | 1.067392 km | 1.056709 km |

- **H3 (labels and errors).** `plot_label` now feeds the "NICER (GR-derived)" legend (lower-left panel) and the colour caption ("Purple: J0740 (2024); brown: J0437 (2024); gray: J0030 (2026, PDT-U)."). A scratch test that changed one active `plot_label` changed the PNG hash. `check_observational_inputs()` verifies that `M_sym_err_Msun` and `R_sym_err_km` equal the mean of the lower and upper errors (tolerance 10⁻⁹) and raises `ValueError` naming the row otherwise; a scratch test with a deliberately wrong J0437 R error stopped the plot with that error. In this package all 10 source/quantity rows pass (largest difference 1.1×10⁻¹⁶, J0030 ST+PDT R). Results are printed per source and stored in `numerical_diagnostics.json`.
- **H4 (README).** The split-run line was replaced by five executable commands: `python3 stellar_tal.py --eos SLy`, `python3 stellar_tal.py --eos APR4`, `python3 stellar_tal.py --eos H4`, `python3 stellar_tal.py --combine`, `python3 diagnostics.py`.
- **Endpoint and core-boundary refinement tests.** `diagnostics.py` recomputes all 33 causal endpoints (14 of them below the sequence maximum) and the core boundary of all 14 acausal maxima with standard vs finer settings (rtol 10⁻¹¹, 1601 grid nodes, 160 quadrature nodes); results in `causal_endpoint_convergence.csv` and `numerical_diagnostics.json`. Maximum absolute differences: endpoint mass 8.5×10⁻⁷ M☉, endpoint radius 3.8×10⁻⁶ km, core boundary 7.0×10⁻⁷ km.
- **Unchanged values.** The full run was repeated (`--eos` ×3, `--combine`, `diagnostics.py`). Apart from the core-radius columns, every pre-existing column of `stellar_summary.csv`, `mr_sequences.csv`, `convergence.csv`, `radius_shifts.csv`, `profiles_L10km.csv` and `sirius_redshift.csv` is numerically identical to the pre-fix V21 outputs, and `wd_summary.json` and `geometry_tests.json` are byte-identical. The causal flags and the radius-shift conclusions are unaffected.

**Still open (not resolved by these fixes; listed in the audit's §5):**

1. **Coupling consistency (task 50aa2db0):** open. The two ansätze are not derived from a complete TAL field equation in this package.
2. **True radial stability:** the \(dM/d\rho_c>0\) test is a heuristic. The perturbation operator, boundary conditions and mode frequencies of the modified equations are not computed; `stable_branch_*` fields must be read with that restriction.
3. **Causal EOS completion:** not done. `M_max_causal_Msun` is the endpoint of the published fit's causal range on the followed branch, not the maximum mass of a causally completed EOS.
4. **NICER inference:** the overlay compares with approximate ellipses of GR-derived posteriors; it is not a TAL pulse-profile likelihood or a Bayes factor.
5. **Low densities:** the equilibrium branch below \(\log_{10}\rho_{0c}=14.55\) is not scanned; the "start" of the stable branch in the log is the grid limit, not a physical stability boundary.
6. **Sirius B:** the additive clock coupling and the retained GR exterior remain separate assumptions; the V21 causality extension does not establish them.

No detection claim or unconditional L bound is made.

**Changelog (audit fixes, V21).** H1 figure styles/endpoint markers and §7.2 wording; H2 core boundary from bracketed threshold + dense ODE output, new sampled-lower-estimate columns, "≤" removed; H3 `plot_label` in legend/caption and fail-loud symmetrized-error check; H4 executable README commands; new causal-endpoint/core-boundary refinement tests (`causal_endpoint_convergence.csv`); tables, figure, log, README, note and SHA256SUMS regenerated.

**Documentation-only clarification, 23 October 2026.** The speculative 182/183 explanation was replaced by the explicit archive-counting rule above, verified using exact decimal CSV values. Only `note_HU_EN.md`, `README.txt` and `SHA256SUMS.txt` were updated for this clarification. The Python code, numerical tables, figures and scientific conclusions are unchanged; no new full equilibrium-sequence run is claimed for this documentation update.

## References and verification trail

External references below were retrieved independently; their specific numbers/parameters were checked in primary text or page images, rather than inferred from search summaries. CiteRight also identified the principal sources; the numerical checks above rely on the primary documents.

1. **Read, J. S., Lackey, B. D., Owen, B. J., & Friedman, J. L. (2009).** *Constraints on a phenomenologically parametrized neutron-star equation of state*. Physical Review D **79**, 124032. [DOI](https://doi.org/10.1103/PhysRevD.79.124032); [primary PDF](https://arxiv.org/pdf/0812.2163). §III–IV, Appendix A, Table II p.19 and Table III p.20 were checked; the table page images confirm SLy \(\log P_1=34.384\).
2. **Salmi, T., et al. (2024).** *The Radius of the High-mass Pulsar PSR J0740+6620 with 3.6 yr of NICER Data*. ApJ **974**, 294. [DOI](https://doi.org/10.3847/1538-4357/ad5f1f); [primary text](https://arxiv.org/html/2406.14466v1). Abstract and §III–IV confirm all quoted headline errors and the distinction between radio prior and posterior.
3. **Choudhury, D., et al. (2024).** *A NICER View of the Nearest and Brightest Millisecond Pulsar: PSR J0437−4715*. ApJL **971**, L20. [DOI](https://doi.org/10.3847/2041-8213/ad5a6f); [primary text](https://arxiv.org/html/2407.06789). Abstract, §III.1, §VI.1 and §VII checked.
4. **Kini, Y., et al. (2026).** *A NICER View of PSR J0030+0451: Updated Constraints from Six Years of NICER Observations*. [arXiv:2602.23743](https://arxiv.org/abs/2602.23743); [primary text](https://arxiv.org/html/2602.23743). Abstract and §III.1 give both hot-spot-model results. The first author is Kini, not Vinciguerra.
5. **Joyce, S. R. G., et al. (2018).** *The gravitational redshift of Sirius B*. MNRAS **481**, 2361–2370. [DOI](https://doi.org/10.1093/mnras/sty2404); [primary PDF](https://arxiv.org/pdf/1809.01240). Abstract, §3.1.7 and Table 3 p.8 checked visually. The paper contains an incidental 78.31 versus 78.21 km/s raw-shift inconsistency; this note uses its consistent, final corrected 80.65±0.77 result.
6. **Bond, H. E., et al. (2017).** *The Sirius System and Its Astrophysical Puzzles: Hubble Space Telescope and Ground-based Astrometry*. [DOI](https://doi.org/10.3847/1538-4357/aa6af8); [primary abstract](https://arxiv.org/abs/1703.10625). It directly gives \(M_B=1.018\pm0.011\,M_\odot\). Joyce's quoted Sirius A mass differs slightly from Bond's abstract; no Sirius B dynamical mass was taken from that inconsistent companion value.
7. **Miller, M. C., et al. (2025).** *The Radius of PSR J0437−4715 from NICER Data*. [arXiv:2512.08790](https://arxiv.org/abs/2512.08790); [primary text](https://arxiv.org/html/2512.08790). Abstract and Table 5 supply the wider, nonthermal-emission sensitivity discussed above.
8. **User's workspace sources:** `user/workspace/uploads/tal_manuscript_v17.md`; `user/workspace/co-scientist-outputs/sparc-tophat-review/independent_review_HU.md`; project memory `tal-sparc-true-tophat-result`. These specify the TAL law and prior finite-source context, not external observational confirmation.
