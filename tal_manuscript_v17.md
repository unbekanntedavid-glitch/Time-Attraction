# The Time Attraction Law (TAL), version 17: formal framework, a locked blind audit of SPARC baryonic baselines, and a falsifiable laboratory test

**David Szabolcsi**

Independent researcher, Unterterzen, Switzerland · ORCID [0009-0005-6300-074X](https://orcid.org/0009-0005-6300-074X)

Manuscript version 17 · 19 October 2026

Version DOI: [https://doi.org/10.5281/zenodo.23479940](https://doi.org/10.5281/zenodo.23479940) · Concept DOI: [https://doi.org/10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681)

Code and data: [https://github.com/unbekanntedavid-glitch/Time-Attraction](https://github.com/unbekanntedavid-glitch/Time-Attraction) · License: CC BY-NC 4.0

---

## Abstract

The Time Attraction Law (TAL) is an alternative gravity proposal developed by a private researcher. Its premise is that time is the primary quantity and that matter density acts as its carrier. Four axioms (normalisation, monotonicity, multiplicativity and continuity) lead through Cauchy's functional equation to a local time factor f(ρ) = exp(−γ_L ρ). The scale is set by γ_L = G_N L²/c², which has dimensions of m³ kg⁻¹. Versions 14–16 added a phenomenological rotation-curve multiplier, C(ρ) = 1 + C₀ exp(−ρ/ρ_c).

This manuscript does four things:

1. It states every formula explicitly, with its dimensions and status (derived, postulated, fitted, or undefined in the source material).
2. It documents version 17, a locked, hash-verified blind audit of 123 SPARC galaxies (2596 radii). The audit uses six purely baryonic baselines and a single error model, σ = max(0.05·V_obs, 20 km s⁻¹). The median reduced χ² ranges from 2.9555 (V_bar × 1.05) to 4.2623 (disk only); for unscaled V_bar it is 3.1433. All twelve output files reproduce byte-for-byte.
3. On the same points and error model, it reports a comparison computed for this manuscript: the McGaugh–Lelli–Schombert radial acceleration relation (RAR) gives a median χ²_red of 1.0941.
4. It specifies a double-blind, pre-registered clock-comparison test of the published law Δf/f ≈ γ_L Δρ̄_L, with the falsification bound γ_L < k·σ_y/Δρ̄_L.

Version 17 contains no TAL density term. In earlier audits by the author, density-dependent multipliers gave no advantage over the RAR on SPARC. The constants C₀ and ρ_c used in versions 14–16 were postulated and never derived. Several author-reported results, including a GPS discrepancy factor, a dwarf-galaxy test and a laboratory apparatus, are labelled as not independently verified. TAL is therefore a testable hypothesis, not an established result.

**Keywords:** alternative gravity; density-dependent time; SPARC; rotation curves; blind analysis; reproducibility; optical clocks; falsifiability

---

## 1. Introduction

General relativity (GR) ties the rate of a clock to the gravitational potential. In the weak field, dτ/dt ≈ 1 + Φ/c². This dependence has been confirmed from the Pound–Rebka experiment to the daily operation of the Global Positioning System, whose satellite clocks run fast by about 38 μs per day relative to ground clocks. That net figure combines a gravitational blueshift of about +45 μs/day and a velocity time dilation of about −7 μs/day, and is standard published physics [5]. TAL asks a complementary question: can the local rate of time depend *directly* on a coarse-grained matter density, rather than only on the potential sourced by distant masses? In TAL's own wording, "the time does not just curve; the time is attracted."[^motto]

The project has gone through seventeen public versions in about three weeks, archived on Zenodo [3] and GitHub [4]. Over that period the formalism grew from a one-parameter clock law into a phenomenological rotation-curve multiplier and a cosmological reinterpretation, and then contracted again. Version 17 deliberately narrows the claim. It does **not** assert that dark matter is absent. Instead it establishes a reproducible *baryonic baseline* on a locked blind SPARC sample, against which any future TAL term must be judged. In the words of the v17 README: "Dark matter implications follow only if independent replication confirms the baseline."

The manuscript is organised as follows:

- §2 states the framework and gives a complete dimensional table.
- §3 goes through the derivations one step at a time and labels each step as derived or postulated.
- §4 documents the SPARC blind audit and adds comparisons computed for this manuscript.
- §5 specifies a laboratory test with explicit falsification criteria.
- §6 lists the limitations.
- §7 gives the version history.
- §8 is the reproduction guide.
- Appendix A records, version by version, every base formula, multiplier, axiom and substitution found in the sources.

A word on sources. Every number in this manuscript comes from one of four places: the v17 package, the earlier versions' PDFs and records, the author's own audit files, or computations on the public SPARC data documented in Appendix B. Where a quantity could not be checked, it is marked *author-reported, not independently verified*. Where the sources define no value or dimension, it is marked *not yet defined*. Only material published in the Zenodo records v1–v17 and the GitHub repository is used for the theory.

## 2. Framework and notation

### 2.1 Core law (axiomatic layer, present since v1)

Let ρ̄_L(x) be the density averaged over a sphere of radius L:

```
ρ̄_L(x) = (3 / 4πL³) ∫_{|x′−x|≤L} ρ(x′) d³x′                                  (1)
```

The time factor f_L ≡ dτ/dt is assumed to satisfy four axioms:

- (A1) f_L = f(ρ̄_L);
- (A2) f(0) = 1;
- (A3) f is monotonically decreasing;
- (A4) f(ρ₁+ρ₂) = f(ρ₁) f(ρ₂), with f continuous.

The v17 PDF restates these axioms and concludes:

```
f(ρ) = exp(−γ ρ),          γ_L = G_N L² / c²                                    (2)
```

In the v17 PDF the exponential form is called "the axiomatic form, independent of the scale of γ_L and of the phenomenological rule V² = V_bar²·C". In other words, the functional form is derived (§3.1), but the value of γ_L depends on a length L that the theory does not fix.

### 2.2 Extended density and multipliers (published in v1–v16)

The first manuscript already extended the density to include mass-equivalent energy densities:

```
ρ_tot = ρ_b + ρ_DM + B²/(2µ₀c²) + c²(R − 4Λ)/(8πG_N) · (1/c²) ,     f = exp(−γ_L ρ_tot)        (3)
```

Version 8.1 rewrote this as a sum of mass-equivalent terms with dimensionless "calibration factors" M_ρ ≈ 1.2 and M_B ≈ 1.1, and with ρ_GR = (c²/8πG)(R − 4Λ). Versions 14–16 replaced the time-factor approach for galaxies with a force multiplier:

```
V² = C(ρ_eff) · V_bar² ,     C(ρ) = 1 + C₀ exp(−ρ/ρ_c)                                             (4)
```

with C₀ = 5, ρ_c = 10⁻²³ kg m⁻³ and L₀ = 1 kpc. In v14.1, ρ_eff = ρ_local + ρ_integrated, where ρ_integrated = (1/L₀)∫₀^d ρ(l) dl. In v16 it becomes a pair-symmetric functional, given in Appendix A.13. Version 17 keeps only the axiomatic law (2) and evaluates baryonic baselines.

### 2.3 Dimensional table

**Table 1.** Symbols, dimensions and status, for all quantities that appear in the published versions.

| Symbol | Meaning | Dimension (SI) | Origin / status |
|---|---|---|---|
| ρ, ρ̄_L | mass density; L-averaged density | kg m⁻³ | defined (v1) |
| L, L₀ | averaging length; "fundamental" length (v1–v7) | m | free; values used: 1 m, 1.3×10⁷ m, 1.2×10⁻¹¹ m, 1 kpc, 2.36 kpc |
| γ_L = G_N L²/c² | time-attraction coefficient (alias G_p, v2) | m³ kg⁻¹ | dimensional construction; scale free |
| f_L = dτ/dt | time factor | 1 | derived from A1–A4 |
| γ_L0 = 1.27×10⁻¹³ | Earth-surface calibration (v1) | m³ kg⁻¹ | one-point calibration; excluded in the laboratory (v1 corrected) |
| γ_lab bound 10⁻²² | laboratory upper bound (v2) | m³ kg⁻¹ | derived from Δf/f < 10⁻¹⁸, Δρ ≈ 10⁴ kg m⁻³ |
| ρ_B = B²/(2µ₀c²) | magnetic mass-equivalent density | kg m⁻³ | standard (v1, v8.1) |
| ρ_GR = (c²/8πG)(R − 4Λ) | curvature mass-equivalent density | kg m⁻³ | v8.1 form correct; the v1 form carries an extra 1/c² |
| M_ρ, M_B; M_ρ = (L/L₀)^(D_f−3) | "calibration factors" (v8.1) | 1 | postulated values "~1.2", "~1.1" |
| γ(ρ) = γ₀/(1+(ρ/ρ_c)ⁿ) | chameleon coefficient (v3) | m³ kg⁻¹ | postulated; ρ_c "~10⁻¹⁸ kg m⁻³", n "1–2" |
| k | gradient coefficient in the v3 Lagrangian | m⁵ kg⁻¹ s⁻² | postulated |
| β (chameleon) | matter coupling in e^{βφ}L_m (v3 HTML) | 1 | postulated |
| κ | coefficient in r(a) = 1 − κ γ_L ρ_L (v8.2) | 1 | postulated (sign only discussed) |
| r, R | α_M/α_B; ΔΣ/Δµ (Horndeski subclass, v8.2–v13) | 1 | cited [6]; algebra derived |
| C₀, ρ_c | multiplier amplitude, critical density (v14–v16) | 1; kg m⁻³ | postulated (v14.0 "Postulate 3"); "fitted" in v14.1 without data or code |
| h_z, Υ_disk, Υ_bul | scale height 0.2 kpc; 0.5; 0.7 (v14.1, v16) | kpc; M_☉/L_☉ | assumed |
| λ₀ = 1.3×10²⁶ m | redshift length (v14.1, v16) | m | assumed ≈ c/H₀ |
| σ = max(0.05 V_obs, 20 km s⁻¹) | v17 velocity uncertainty | km s⁻¹ | chosen by the author (v17) |

## 3. Step-by-step derivations and their status

### 3.1 Exponential form (derived)

Let g(ρ) = −ln f(ρ). Axiom A4 gives g(ρ₁+ρ₂) = g(ρ₁) + g(ρ₂), and A2 gives g(0) = 0. The only continuous solution of Cauchy's additive equation is linear, g(ρ) = γρ. Monotonicity (A3) requires γ > 0, so f(ρ) = exp(−γρ). This step is mathematically sound. It fixes only the *form* of f. The coefficient γ is a free constant.

### 3.2 Choice of γ_L = G_N L²/c² (dimensional construction; L postulated)

The product γρ must be dimensionless. The combination G_N L²/c² has dimensions m³ kg⁻¹, so it qualifies. Early versions add a motivation:[^vr] for a uniform sphere of radius L, |Φ|/c² = (4π/3)·G_N ρ̄_L L²/c² = (4π/3)γ_L ρ̄_L. Hence γ_L ρ̄_L is, up to the factor 4π/3, the familiar weak-field potential ratio. This identity is exact, but it also shows that the weak-field law reproduces GR's potential dependence only when L is the size of the source. It does not supply an independent density effect.

### 3.3 Scale of the effect in galaxies (derived)

Take L ≈ r and a dynamically inferred density ρ ~ v²/(Gr²). Then γ_L ρ = G r² ρ/c² ~ (v/c)² ~ 10⁻⁷ for v ~ 100–300 km s⁻¹. The author's audit on SPARC confirms this: the largest value is γ_Lρ ≈ 9.5×10⁻⁸ for Υ_disk = 0.5 and 1.9×10⁻⁷ for Υ_disk = 1.0. An effect of order unity (γ_L ρ = 1) would need L_unity = c/√(Gρ), whose median is 85.3 Mpc for Υ_disk = 0.5 and 64.2 Mpc for Υ_disk = 1.0. That is about 10⁴ times the outermost measured radius. The v1–v7 documents had already concluded that the exponential law "alone does not replace dark matter".

### 3.4 Rotation-curve multiplier C(ρ) (postulated, never derived)

Version 14.0 introduced, as "Postulate 3":

```
C(ρ) = 1 + C₀ exp(−ρ/ρ_c),   V² = C·V_bar²,   C₀ = 5,  ρ_c = 10⁻²³ kg m⁻³        (5)
```

Version 14.1 described these constants as "fitted" (C₀ = 5.0 ± 0.4, ρ_c = (1.0 ± 0.2)×10⁻²³ kg m⁻³, median χ²_red = 1.18). The promised data file ("v14.1-data.zip") does not exist. The accompanying derivation, from an ordinary differential equation dX/dl = −X ρ/(ρ_c L₀) with X = 1/C, yields an exponential of an optical-depth integral, not 1/[1 + 5 exp(−ρ/ρ_c)]. Version 16 itself retracted that step ("ODE solution is exponential of optical-depth integral, not algebraic") and called C a "force multiplier ansatz". The critical density has appeared with three different values: 10⁻¹⁸ kg m⁻³ (v3), 10³ kg m⁻³ (the early chameleon document) and 10⁻²³ kg m⁻³ (v14–v16). The author's refit on SPARC training galaxies gave C₀ = 2.776 and ρ_c = 3.562×10⁻²¹ kg m⁻³, not the published values (§4.7). **Status: postulate with unverified constants.**

### 3.5 Action and conservation (open)

The v1–v7 formalism writes a genuine action, S = ∫d⁴x √−g [R/16πG + f_L L_m]. It gives the field equation only as G_µν = 8πG_N f_L T_µν + O(δf_L/δg), and it does not show that the effective stress–energy tensor is conserved, which a non-constant f_L generally violates. The v3 Lagrangian density, −ρ_L c² exp(−γρ_L) − ½k(∇ρ_L)², is postulated and has no stated field equation. The rotation-curve multiplier of v14–v16 has no action at all. Version 16 states this itself: "Force multiplier ansatz remains unless derived from field dynamics." It does, however, construct a conservative pair potential U_r(r) = −A∫_r^∞ C(s)/s² ds. Any future TAL term needs an action from which ∇_µT^µν = 0 follows.

## 4. SPARC blind audit (v17)

### 4.1 Data and sample

The input is the public SPARC archive of rotation curves and mass models [1], `Rotmod_LTG.zip`, with SHA-256 `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`. The blind list `blind_123.txt` has SHA-256 `3637d2f5bda6b2b97e645c5fba5cf39d310a1826bee388286c2c74b9171aa283`. It contains 123 galaxies with 2596 radii, and the runner stops with an error if the point count differs. The v17 PDF describes the split as SHA256("sparc-holdout-20261010:" + name), with 30 holdout galaxies "earlier explored, presented separately". The threshold of the hash split is not stated.

As a consistency check made for this manuscript, the SPARC catalogue [1] contains 153 galaxies with quality flag Q ∈ {1, 2} and inclination ≥ 30°. The 123 blind galaxies are a subset of these 153, and the 30 that remain include NGC 3198. The blind sample holds 78 Q = 1 and 45 Q = 2 galaxies, and its median disk scale length is R_disk = 2.36 kpc. Both values match those printed in the v17 PDF ("L₀ = 2.36 kpc", and 0.83 kpc for the 23 galaxies with Hubble type T ≥ 10).

### 4.2 Models

The runner `sparc_full_runner.py` reads the columns R, V_obs, V_gas, V_disk and V_bul and computes

```
V_bar = √(V_gas² + V_disk² + V_bul²)                                            (6)
```

V_disk and V_bul are used as tabulated in SPARC, i.e. for a mass-to-light ratio of 1 M_☉/L_☉ at 3.6 μm [1]. V_gas is squared, so the sign that SPARC uses to mark central gas depressions is lost. This affects 298 of the 2596 points in 33 galaxies. Six predictions are evaluated (Table 2). All six are **baryonic baselines**: none contains a TAL density term, a fitted parameter or a dark component. The v17 editorial letter states that "the phenomenological application of the TAL factor and dark-matter replacement fall outside the scope of this audit."

### 4.3 Error model and statistic

```
σ_i = max(0.05 · V_obs,i , 20 km s⁻¹)
χ² = Σ_i [(V_pred,i − V_obs,i)/σ_i]² ,     χ²_red = χ² / (N − 1)                (7)
```

Here N is the number of radii per galaxy, and the summary statistic is the median of χ²_red over galaxies. The largest V_obs in the sample is 383 km s⁻¹, below the 400 km s⁻¹ at which the 5 % term would exceed the floor. The 20 km s⁻¹ floor therefore applies at **all** 2596 points, as the v17 PDF states. In practice σ is a constant, and χ²_red is a scaled mean-square residual. This choice is robust to underestimated published errors, but it reduces the power to discriminate between models, especially for low-velocity galaxies.[^errV]

### 4.4 Results

**Table 2.** v17 baryonic baselines on the locked 123-galaxy blind set. The medians are those printed in the v17 PDF and reproduced exactly. The quartiles, means and fractions were computed for this manuscript from the v17 `summary_*.csv` files.

| Model | V_pred | Median χ²_red | Q1–Q3 | Mean | Fraction χ²_red < 1 |
|---|---|---|---|---|---|
| vbar | V_bar | **3.1433** | 1.461–4.934 | 4.216 | 17.9 % |
| vbar_095 | 0.95 V_bar | 3.1995 | 1.671–4.820 | 4.373 | 17.1 % |
| vbar_105 | 1.05 V_bar | 2.9555 | 1.296–5.512 | 4.322 | 17.9 % |
| no_bulge | √(V_gas² + V_disk²) | 3.3131 | 1.478–7.051 | 6.393 | 17.9 % |
| no_gas | √(V_disk² + V_bul²) | 4.0474 | 2.243–6.638 | 5.253 | 10.6 % |
| disk_only | \|V_disk\| | 4.2623 | 2.295–8.012 | 7.571 | 10.6 % |

Over all points, V_bar < V_obs at 66.8 % of the radii, and the median normalised residual (V_bar − V_obs)/σ is −0.81: the baryons alone under-predict the observed velocities, as expected. Section 6 of the v17 PDF states that varying Υ, V_disk, the gas profile, h_z, the bulge or L₀ "did not give a better fit than pure V_bar". On the six released baselines this does not hold exactly. The model vbar_105 has a lower median than vbar, and it has a lower χ²_red in 91 of 123 galaxies (median difference 0.190; two-sided Wilcoxon signed-rank p = 0.015). The editorial statement should be read as "no *TAL* variant improved on the baselines". The improvement from a global 5 % increase in V_bar is itself an expected sign of missing mass or acceleration.

### 4.5 Comparison with the radial acceleration relation (computed for this manuscript)

To give the baselines an external reference, the RAR of McGaugh, Lelli & Schombert [2] was evaluated on the **same 2596 points, with the same V_bar and the same error model**:

```
g_bar = V_bar²/R,    g_obs = g_bar / (1 − exp(−√(g_bar/a₀))),    a₀ = 1.2×10⁻¹⁰ m s⁻²,    V_RAR = √(g_obs R)      (8)
```

**Table 3.** Paired comparison on the 123 blind galaxies (computed for this manuscript; scripts in Appendix B). Here d = χ²_red(vbar) − χ²_red(RAR).

| Quantity | vbar | RAR |
|---|---|---|
| Median χ²_red | 3.1433 | **1.0941** |
| Q1–Q3 | 1.461–4.934 | 0.399–6.023 |
| Fraction χ²_red < 1 | 17.9 % | 48.0 % |
| Galaxies in which the model is better | 53 | 70 |
| One-sided Wilcoxon (H1: RAR better), p | — | 0.176 |
| One-sided sign test, p | — | 0.074 |

The RAR has a markedly lower median. Because the 20 km s⁻¹ floor dampens per-galaxy differences, the paired tests do not reach p < 0.05 on the full sample. The RAR is a one-parameter empirical relation with a₀ fixed from the literature, and it was not refitted here.

### 4.6 Dwarf subset (computed for this manuscript)

The author reports that a dwarf-galaxy test found TAL to be 240 times worse than the RAR. Neither that test's data, its code nor its model definition was supplied, so this figure is **author-reported, not independently verified**. The comparison that can be computed from the released data uses the baselines, not a TAL term (Table 4).

**Table 4.** Baryons-only (vbar) versus the RAR on dwarf subsets of the blind sample, using the same error model. Dwarf criteria are taken from the SPARC catalogue [1].

| Subset | N | Median χ²_red vbar | Median χ²_red RAR | RAR better / worse | Two-sided Wilcoxon p |
|---|---|---|---|---|---|
| T ≥ 10 (Im/BCD; median R_disk 0.83 kpc) | 23 | 0.520 | 0.588 | 12 / 11 | 0.33 |
| T ≥ 9, or 0 < V_flat < 60 km s⁻¹, or DDO | 42 | 1.537 | 0.428 | 28 / 14 | 4.7×10⁻⁴ |
| Complement (non-dwarfs) | 81 | 3.757 | 3.479 | 42 / 39 | 0.57 |

For the 23 latest-type galaxies, the 20 km s⁻¹ floor exceeds the typical velocities so much that both models give χ²_red < 1 and cannot be told apart. This is a property of the error model, not evidence for either model. With the broader definition, the RAR is significantly better than baryons alone.

### 4.7 Earlier density-term audits (author's audit; different sample and error model)

Before v17, the author audited density-dependent multipliers with a clean Python pipeline. That audit used a temporary blind set of N = 114 galaxies (Q ∈ {1, 2}, inclination > 30°, at least 5 radii), the published errV as the only error term, and Υ_disk = 0.5. The results cannot be compared numerically with Tables 2–4, but they bear directly on the TAL term:

- **Photometric density proxy (v16 Eq. 5–6).** The median χ²_red was 78.44 for C = 1 and 76.86 for the v16 C(ρ). A refit on training galaxies (C₀ = 2.776, ρ_c = 3.562×10⁻²¹ kg m⁻³) gave 21.70. The RAR gave 11.19.
- **NGC 3198 (holdout).** At R = 30.08 kpc, V_obs = 146 ± 2 km s⁻¹, and V_bar = 70.32 km s⁻¹ (Υ = 0.5) or 88.53 km s⁻¹ (Υ = 1.0). The galaxy's χ²_red (Υ = 0.5) was 776.82 for C = 1, 271.99 for the v16 local variant, 776.82 for the integrated variant and 13.05 for the RAR.
- **Dynamical density** ρ = d(rV_bar²)/dr/(4πGr²). The "integrated" v16 structure gave C ≈ 1 everywhere, identical to C = 1. The local variant improved on C = 1 but was worse than the RAR in 105 of 114 galaxies (medians 87.57 versus 12.26).
- **Negative densities.** The derivative-based density was negative at 239 interior points in 53 galaxies (Υ = 0.5). These points were flagged and excluded, never clipped. For C = 1 + 5 exp(−ρ/ρ_c), a negative ρ gives C > 6 rather than C = 1.

### 4.8 Statement on predictive advantage

On the evidence available, **no TAL density term has shown a predictive advantage on SPARC.** Version 17 does not evaluate one. Every earlier evaluation the author could reproduce performed worse than the RAR. The contribution of v17 is a locked, reproducible baseline and test harness, not support for TAL.

## 5. Laboratory proposal: clock comparison with a movable mass

### 5.1 Published prediction

All laboratory proposals in the published versions test the weak-field form of the core law. They are the two-clock lead protocol of v2, the Yb⁺ clock and tungsten-cube protocol of v1–v7, and the Sr clock and 10 t lead pitch of v16. That weak-field form is

```
Δf/f ≈ γ_L Δρ̄_L ,     γ_L = G_N L²/c²                                         (9)
```

Here Δρ̄_L is the change of the L-averaged density at the probe clock when a dense mass is moved between a near and a far position. For L = 1 m, γ_L = 7.4×10⁻²⁸ m³ kg⁻¹, and v2 estimates a predicted shift of about 3.7×10⁻²⁴ for a 10 t lead tower. That is far below present clock instability, which is why v1–v7 and v2 frame the experiment as a *bound*.

### 5.2 Falsification bound

Let σ_y be the fractional frequency instability (Allan deviation) at the integration time used, and k a pre-registered coverage factor (k = 2 for about 95 %). A null result implies

```
γ_L < k · σ_y / Δρ̄_L ,      equivalently   L < √( k σ_y c² / (G_N Δρ̄_L) )        (10)
```

With σ_y = 10⁻¹⁸ [7], Δρ̄_L = 10⁴ kg m⁻³ and k = 1, Eq. (10) gives γ_L < 10⁻²² m³ kg⁻¹ and L < 367 m. These are exactly the bounds published in v2 and v1–v7. The ordinary gravitational effect of the moved mass is a confound and must be subtracted. For 10 t at 0.5 m it is G_N M/(r c²) ≈ 1.5×10⁻²³, and a 1 mm height change gives gΔh/c² = 1.1×10⁻¹⁹ (v1–v7). Both are below 10⁻¹⁸, but both must be modelled. The v16 pitch's predicted amplitude of 5×10⁻¹⁹ does not follow from Eq. (9) for any stated L (Appendix A.13), so it should not be used as the pre-registered amplitude.

### 5.3 Protocol (proposed, building on the published v1–v7, v2 and v16 protocols)

1. **Pre-registration.** Before any data are taken, deposit the hypothesis, the predicted sign and amplitude, L, the modelled Δρ̄_L, σ_y, k, the integration time, the analysis code and the exclusion rules.
2. **Randomised AB blocks.** Alternate mass-NEAR (A) and mass-FAR (B) states, as in the published A–B–A cycle and the rail-mounted mass, in blocks whose order comes from a sealed random sequence generated by a third party. Include clock swaps, as in phase C of v2, to cancel clock-specific systematics.
3. **Double blinding.** The operator sees only block labels. The analyst receives frequency data with scrambled labels, which are revealed only after the analysis is frozen.
4. **Dummy control.** Make identical moves of a low-density dummy of equal outer dimensions, which reproduces vibration, thermal load and motor currents with negligible Δρ̄_L. A signal that also appears with the dummy is systematic.
5. **Environmental monitors.** Log the field at the clock with a fluxgate magnetometer, together with temperature, vibration, tilt and height. Any correlation of the frequency difference with these channels is treated as systematic unless it is pre-registered otherwise.
6. **Falsification criteria.** The law is bounded or rejected at this sensitivity if: (i) |Δf/f − GR model| < k σ_y, which yields bound (10); or (ii) the dummy reproduces the signal; or (iii) the signal tracks an environmental channel; or (iv) its sign or its scaling with Δρ̄_L, or with distance as L varies, contradicts Eq. (9). A positive result counts only after independent replication in another laboratory.

**Author-reported apparatus (not verifiable from the repository).** The author reports a laboratory set-up consisting of an 8-column solenoid, a granite chamber, 12/24 AB-swap blocks, and logging with a PicoScope and a Rigol oscilloscope. None of these is documented in the repository, the Zenodo records or the v17 package, and no data were supplied. They are therefore not verifiable from the published material and play no role in this manuscript.

## 6. Limitations

1. **No empirical advantage.** No TAL density term outperforms the RAR on SPARC in any reproducible test (§4.7–4.8). Version 17 tests baselines only.
2. **Undetermined constants.** C₀ and ρ_c were postulated, never derived, and ρ_c has taken three different values. L has taken at least five values. The multipliers M_ρ, M_B (v8.1) and κ (v8.2) have no derivation.
3. **GPS.** The GPS clock offset of about 38 μs/day is fully accounted for by GR and special relativity [5]. The author reports that applying TAL with fixed parameters to this offset gives a discrepancy of about 13,000×. The calculation is not published, so the figure is *author-reported, not independently verified*. If confirmed, it would falsify that parameterisation. Making a parameter density-dependent after the fact would add freedom, not a prediction.
4. **Dwarf galaxies.** The author's "240× worse than the RAR" result is *author-reported, not independently verified* (§4.6).
5. **Error model.** With a 20 km s⁻¹ floor active at every point, the test discriminates weakly, especially for dwarfs. The floor was chosen by the author, and the record does not show whether it was fixed before the outcomes were known.
6. **Baryonic inputs.** Υ = 1 is used for both disk and bulge, gas signs are discarded (298 points), and the runner replaces a first radius below 0.01 kpc by 0.16 kpc (this never triggers on the blind sample).
7. **Holdout.** The 30 holdout galaxies were "earlier explored", so they are not a pristine test set. The hash-split threshold is not published.
8. **Theoretical consistency.** The v1–v7 action is incomplete (field equation only to O(δf_L/δg)), and the v14–v16 multiplier has no action (§3.5). The redshift law of v14–v16 needs a metric before it can be compared with BAO data, as v16 itself notes.
9. **Documentation.** The v17 README refers to `sparse_full_runner.py` and `Rotmod_data.zip`, whereas the actual files are `sparc_full_runner.py` and `Rotmod_LTG.zip`. The runner's docstring describes an extracted `data/rotmod/` folder, whereas the code requires the zip. The v17 PDF also contains the sentence "γ_L = 10²³ conflicts with the lab" without units. The ORCID iD differs between versions.[^orcid]

## 7. Version history (v1–v17)

The v17 `CHANGELOG.md` lists only v17, v15, v13, v12.1 and v8.2. Table 5 reconstructs the full history from the Zenodo records (publication dates from the record metadata), the repository PDFs and the git history. The version labels are inconsistent in the sources themselves, and gaps are stated. Formulas for each entry are given in Appendix A.

**Table 5.** Version history.

| Label (as in source) | Record / file | Date | Main change |
|---|---|---|---|
| v1 | Zenodo 22983682: *publikacio*, *javitott* | 2026-09-27 (documents dated 2025-08-17 and 2026-09-26) | Axioms, Cauchy derivation, f_L = exp(−γ_Lρ̄_L), γ_L = G_N L²/c²; calibration γ_L0 = 1.27×10⁻¹³ m³ kg⁻¹ (L₀ ≈ 1.3×10⁷ m); the V/r ≈ 1 error corrected |
| v2 ("Verzio 2.0") | Zenodo 22984786 (zip; also 22984952, 22984968) | 2026-09-27 | 12-page corrected version; notation Gp ≡ γ_L; laboratory bound γ_lab < 10⁻²² m³ kg⁻¹; DE/EN/arXiv variants |
| v3 ("v3.0") | in the same zip: *Javitott v3 Sulyosabb*, *Publikacio 3oldal*, HTML | 2026-09-27 | Chameleon γ(ρ) = γ₀/(1+(ρ/ρ_c)ⁿ); Lagrangian density; chameleon action |
| V5 | Zenodo 23002701 | 2026-09-27 | Bilingual formalism; Sr clock protocol |
| 4-in-1 | Zenodo 23045757 | 2026-09-29 | Compilation of the four parts |
| "v1–v7" (PDF header: "v13 FINAL") | Zenodo 23066993 | 2026-09-30 | Full formalism: L₀ bounds, six theorems, action, field equation, Yb⁺/W protocol |
| v8.1 | Zenodo 23101489 | 2026-10-02 | Effective-density framework: mass-equivalent terms ρ_B, ρ_E, ρ_GR; multipliers M_ρ, M_B; ten galaxy densities ρ_tot = v_c²/(4πGr²) |
| v8.2 / v8.2.1 | Zenodo 23135605; repo v8.2.1.pdf | 2026-10-04 | Corrigendum R = (1+r)/(1+2r) (credited to J. Audet Palau); α ansatz; predicted evolution of R(a) with scale factor |
| — | Zenodo 23137284 (derivation.md) | 2026-10-04 | "150 SPARC" text; claimed χ² improvement 7.7×–16.8× (no data) |
| v8.3 | Zenodo 23138388; repo v8.3_.pdf | 2026-10-04 | Closed form ε = XG₃X/G₄φ − 2; f_L cancels in r; background ρ_m(a) |
| v9–v11 | **no separate record or PDF found** | — | The repo README header reads "v8.3 FINAL (v11)" |
| v12.1 | Zenodo 23148414 | 2026-10-05 | R(a) background correction, ρ_L → ρ_m(a) = ρ_m0 a⁻³ |
| v13 | Zenodo 23193655 | 2026-10-06 | α_B with full denominator H M*² |
| v14.0 | Zenodo 23217189 | 2026-10-07 | Postulate 3: C(ρ) = 1 + C₀e^(−ρ/ρ_c); F_TAL; redshift reinterpretation |
| v15 | Zenodo 23220318 (file is v14.1); GitHub README | 2026-10-07 | Constants C₀ = 5, ρ_c = 10⁻²³, L₀ = 1 kpc; invalid ODE derivation; unreproduced fit |
| v16.0 | Zenodo 23234825 | 2026-10-08 | Pair functional q, conservative potential U_r, photometric proxy, unified z(d); Pb–Sr laboratory pitch |
| v17 | Zenodo 23479940; repo v17_package | 2026-10-19 | Locked blind baryonic baseline audit, single error model, hash checks |

## 8. Reproduction guide

```
git clone https://github.com/unbekanntedavid-glitch/Time-Attraction.git
cd Time-Attraction/v17_package/github
mkdir -p data                        # blind_123.txt is already in data/
# download Rotmod_LTG.zip from the SPARC site (Lelli et al. 2016) and place it, still zipped, in data/
sha256sum data/Rotmod_LTG.zip data/blind_123.txt
#   0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588  data/Rotmod_LTG.zip
#   3637d2f5bda6b2b97e645c5fba5cf39d310a1826bee388286c2c74b9171aa283  data/blind_123.txt
python3 sparc_full_runner.py         # standard library only; writes output/
for f in output/*.csv; do cmp "$f" "run_output/$(basename "$f")" && echo "OK $f"; done
sha256sum output/*.csv run_output/*.csv   # each pair must match
```

The runner checks both hashes and the count of 2596 radii, and it exits with an error on any mismatch or missing galaxy. Its console output must report the six medians in Table 2. For this manuscript the procedure was repeated in a clean directory with Python 3.13.5, and all 12 CSV files matched byte-for-byte. The additional statistics (Tables 2–4) are produced by the two scripts listed in Appendix B.

## Data and code availability

Code, the blind list and reference outputs: [https://github.com/unbekanntedavid-glitch/Time-Attraction](https://github.com/unbekanntedavid-glitch/Time-Attraction) (`v17_package/`). Archived version: [https://doi.org/10.5281/zenodo.23479940](https://doi.org/10.5281/zenodo.23479940). SPARC data: Lelli, McGaugh & Schombert [1].

## Competing interests and funding

The author has previously declared "no conflict, no institutional funding" (Zenodo record 23045757, Part 4). No external funding supported this work.

## License

© 2026 David Szabolcsi. This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0), [https://creativecommons.org/licenses/by-nc/4.0/](https://creativecommons.org/licenses/by-nc/4.0/).

## References

1. Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass models for 175 disk galaxies with Spitzer photometry at 3.6 μm and accurate rotation curves. *The Astronomical Journal*, 152, 157. [https://doi.org/10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157)
2. McGaugh, S. S., Lelli, F., & Schombert, J. M. (2016). Radial acceleration relation in rotationally supported galaxies. *Physical Review Letters*, 117, 201101. [https://doi.org/10.1103/PhysRevLett.117.201101](https://doi.org/10.1103/PhysRevLett.117.201101)
3. Szabolcsi, D. (2026). *Time Attraction Law v17 – SPARC baryonic baselines locked blind audit*. Zenodo. [https://doi.org/10.5281/zenodo.23479940](https://doi.org/10.5281/zenodo.23479940) (concept DOI [https://doi.org/10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681)).
4. Szabolcsi, D. (2026). *Time-Attraction* (software repository). GitHub. [https://github.com/unbekanntedavid-glitch/Time-Attraction](https://github.com/unbekanntedavid-glitch/Time-Attraction)
5. Ashby, N. (2003). Relativity in the Global Positioning System. *Living Reviews in Relativity*, 6, 1. [https://doi.org/10.12942/lrr-2003-1](https://doi.org/10.12942/lrr-2003-1)
6. Bellini, E., & Sawicki, I. (2014). Maximal freedom at minimum cost: linear large-scale structure in general modifications of gravity. *Journal of Cosmology and Astroparticle Physics*, 2014(07), 050. [https://doi.org/10.1088/1475-7516/2014/07/050](https://doi.org/10.1088/1475-7516/2014/07/050) (cited in v8.2–v13 and v17).
7. Bothwell, T., et al. (2022). Resolving the gravitational redshift across a millimetre-scale atomic sample. *Nature*, 602, 420–424. [https://doi.org/10.1038/s41586-021-04349-7](https://doi.org/10.1038/s41586-021-04349-7) (cited in v17 and earlier versions).
8. Khoury, J., & Weltman, A. (2004). Chameleon fields: awaiting surprises for tests of gravity in space. *Physical Review Letters*, 93, 171104. [https://doi.org/10.1103/PhysRevLett.93.171104](https://doi.org/10.1103/PhysRevLett.93.171104) (cited in v2–v3).

Earlier Zenodo versions cited in Table 5 resolve as https://doi.org/10.5281/zenodo.NNNNNNNN with the record numbers given.

[^motto]: Motto used in the v16 and v17 README files.
[^vr]: The first manuscript derived T = Gρ/c² using V/r ≈ 1. The 2026-09-26 revision corrected this, noting that [V/r] = m², not a pure number.
[^errV]: The pre-v17 audits (§4.7) used the published errV only, which gives χ²_red values one to two orders of magnitude larger. Numbers from the two error models must not be compared directly.
[^orcid]: For example, 0009-0009-1901-4830 in v8.2, 0009-0000-6300-074X in v12–v14.1, and 0009-0005-6300-074X in v8.3, v15–v17 and this manuscript.

---

## Appendix A. Formula notebook, v1 → v17

*How to read this notebook.* Each page below corresponds to one version or record and keeps the order in which the author wrote things down: first the base formula, then multipliers and substitutions, then the derivation and the author's own remarks. Formulas are copied **as they appear in the source**, in the source's own ASCII or Unicode spelling (the original language is kept where it matters, with a translation in brackets). Margin remarks marked ✎ are editorial: dimension checks, status, and cross-references. A version for which no file exists is listed as such. The notebook covers the gravitational, cosmological and laboratory formulas of the published records.

### A.1 v1 — Zenodo 22983682 (record 2026-09-27; manuscript dated "2025. augusztus 17.")

*Definitions*

```
ρ̄_L(r) = 3/(4πL³) ∫_{|r′−r|<L} ρ(r′) d³r′                 (1)
G_N = 6.67430 × 10⁻¹¹ m³ kg⁻¹ s⁻²,   c = 299792458 m/s
γ_L = G_N L² / c² ,   [γ_L] = m³/kg  "idővonzási együttható L skálán"   (5)
f_L = dτ/dt  "időfaktor, dimenziótlan"                                  (6)
```
✎ The readings "L = 1 m laboratory" and "L ~ 10⁷ m planetary" are given in words.

*Axioms (Axiómák)*

```
1. f_L = f(ρ̄_L)    2. f_L(0) = 1    3. f_L monotonically decreasing
4. f(ρ1 + ρ2) = f(ρ1) f(ρ2), f continuous
```

*Derivation*

```
g(ρ) = −ln f(ρ)  ⇒  g(ρ1+ρ2) = g(ρ1) + g(ρ2),  g(0) = 0          (7)
Cauchy  ⇒  g(ρ) = γ_L ρ  ⇒  f_L(ρ̄_L) = exp(−γ_L ρ̄_L)              (8)
f_L ≈ 1 − γ_L ρ̄_L                                                  (9)
Δf_{A→B} = e^{−γ_L ρ_B} − e^{−γ_L ρ_A}                             (10)
Δf_{A→B} ≈ −γ_L(ρ_B − ρ_A),   Δf/f ≈ γ_L Δρ                         (11–12)
Δ(ns/s)_{A→B} ≈ γ_L (ρ_B − ρ_A) · 10⁹                               (13)
```
✎ Status: derived (form only). The sign convention changes between (11) and (12).

*Link to the potential*

```
M = 4πL³ρ̄_L/3,  Φ = −G_N M/L  ⇒  |Φ|/c² ≈ (4π/3)(G_N L²/c²) ρ̄_L ∼ γ_L ρ̄_L     (14)
```

*Extended model (first appearance of the energy-density terms)*

```
ρ_tot = ρ_b + ρ_DM + B²/(2µ0c²) + c²(R − 4Λ)/(8πG_N c²)              (15)
f(r,t) = exp[−γ_L ρ_tot(r,t)]                                        (16)
ρ_DM(r,t) = (1/(4πG_N r²)) d/dr [r v_c(r,t)²] − ρ_b(r,t)             (17)
```
✎ Eq. (17) infers dark matter from the rotation curve, i.e. the model includes ρ_DM rather than replacing it. Dimension check of the curvature term in (15): [R/G] = kg s² m⁻⁵, so (R − 4Λ)/(8πG_N) is not a density. The correct mass-equivalent density is c²(R − 4Λ)/(8πG_N) without the extra 1/c², which is the form v8.1 later uses (A.6).

*v1 corrected file (javitott, 2026-09-26): the substitutions that were fixed*

```
original:   Φ = G_N M/r = G_N ρV/r ,  "V/r ≈ 1"  ⇒  T = Gρ/c²      ✗ [V/r] = m²
corrected:  Φ ≈ (4π/3) G_N ρ r² ,   Φ/c² = (G_N r²/c²) ρ = γ_r ρ     (7)
calibration: Φ/c² ≈ G_N M/(R c²) ≈ 7.0 × 10⁻¹⁰ (Earth surface)
             γ_L0 = 7.0 × 10⁻¹⁰ / 5510 = 1.27 × 10⁻¹³ m³/kg           (8)
             L0 = √(γ_L0 c²/G_N) ≈ 1.3 × 10⁷ m  "Ez a Föld átmérője" [this is the Earth's diameter]   (9)
             γ_lab = G_N (1 m)²/c² ≈ 7.4 × 10⁻²⁸ m³/kg                   (10)
ρ_tot = ρ_b + ρ_DM + B²/(2µ0c²) + c²(R − 4Λ_cos)/(8πG_N c²) + ρ_other   (11)
"Nem G_p/ρ_tot, hanem exp(−γ_L ρ_tot)."  [not G_p/ρ_tot but exp(−γ_L ρ_tot)]
bound: γ_lab < 10⁻²² m³/kg (L = 1 m)                                   (14)
```
✎ γ_L0 is a one-point calibration, not a derivation. The same file states "a 1.27e−13 érték már ki van zárva laborban" [the value 1.27e−13 is already excluded in the laboratory]. Its dark-matter remark: γ_L ρ ~ 1 "would need L ~ 10²⁶ m … not physical".

### A.2 v2 — Zenodo 22984786 (zip "fwdmetaminden odsa.zip"; duplicates 22984952, 22984968), 2026-09-27

*"Javitott, merheto, 12 oldalas valtozat — Verzio 2.0" (2026.09.26).* This file records what the very first manuscript contained (formulas as quoted in v2):

```
A) exponential branch (original):  f(rho) = d tau/dt = exp(-G rho)
B) linear branch (later):          T = (G * rho)/c^2      "ahol G a Newton-fele"
old table: t_rel ∝ 1/rho  ("vakuum 1e30, viz 1")       ✎ removed in v2: a third model, inconsistent with exp
old calibration: G = (7e-17)/rho_Fold                 ✎ "a 7e-17 eredete nincs levezetve" [origin of 7e-17 not derived]
old dark-matter ansatz: v_obs^2 = (G_N M(r)/r) * exp(2 gamma rho_eff),  M_eff = M_barion * exp(2 gamma rho_eff)
```
✎ v2 removes the linear branch B (it reproduces GR's Φ/c²), the 1/ρ table and the e^{2γρ} dark-matter ansatz, and keeps only the exponential branch.

*v2 final collection*

```
d tau/dt (x) = exp[-gamma * rho_L(x)],  [gamma] = m3/kg,  bound gamma < 1e-22
weak field: d tau/dt ≈ 1 - gamma rho_L,   Delta f/f ≈ -gamma Delta rho_L
bound: Delta f/f < 1e-18 (lead, Delta rho ~ 1e4)  ⇒  gamma < 1e-18/1e4 = 1e-22 m3/kg
prediction: 10 t lead, L = 1.5 m:  Delta f/f ≈ 1.1e-18 * (gamma / 1e-22)
dark matter: exp(2 gamma rho_eff) ≈ 10  ⇒  rho_eff ≈ 1/gamma ≈ 1e22 kg/m3  ("neutroncsillag suruseg")
```

*"Hivatalos" (official) and DE/Gp variants: notation only.*

```
f_L = exp(−γ_L ρ̄_L),   γ_L ≡ G_p = G_N L²/c²    "G_p az idővonzás állandója"
L = 1 m:  γ_L^lab = 7.4 × 10⁻²⁸,  Δf/f ≈ 3.7 × 10⁻²⁴;  bound γ_L^lab < 10⁻²²
```

### A.3 v3 — same zip: "Az Idovonzas Torvenye – v3.0, Javitott, tudomanyosan sulyosabb valtozat"

```
1. d tau/dt (x) = exp[ -gamma(rho_L) * rho_L(x) ]
2. rho_L(x) = 3/(4 pi L3) * integral rho(x+x') d3x' , |x'| < L
3. gamma(rho) = gamma0 / (1 + (rho/rho_c)^n)        "rho_c kb 1e-18 kg/m3 galaktikus atlag, n kb 1-2"
4. L = -rho_L c2 exp(-gamma rho_L) - 0.5 k (nabla rho_L)2      "k merteke m5/(kg s2)"
5. weak field: L ≈ -rho_L c2 + gamma rho_L2 c2  ("taszito jellegu" [repulsive])
6. Delta f/f ≈ 1.1e-18 * (gamma0/1e-22) * screening factor   (10 t lead, L = 1.5 m)
7. gamma0 < 1e-22 m3/kg without screening; with screening gamma0 < 1e-12 "still allowed"
example: gamma0 = 1e-12, rho_c = 1e-15, n = 1:
         gamma_lab = 1e-12/(1+1e3/1e-15) = 1e-27 ;  gamma_gal = 1e-12
         exp(2 gamma rho) = exp(2e-12 * 1e-21) = 1 + 2e-33   "meg mindig semmi" [still nothing]
```
✎ The Lagrangian density (4) is postulated. The dimension given for k, m⁵/(kg s²), makes k(∇ρ)² have J m⁻³ ✓. No field equation is written out. ρ_c = 10⁻¹⁸ kg m⁻³ is the first of three ρ_c values.

*Early HTML ("gyönyörű publikációs anyag"): chameleon embedding*

```
f_cham(ρ) = exp[ −γ ρ / (1 + (ρ/ρ_c)ⁿ) ]          (figure: "ρ_c = 10³ kg/m³, n = 1")
S_φ = ∫ d⁴x √−g [ −½ (∂φ)² − V(φ) ] ,   V(φ) ≈ ½ m₀² φ²
S = ∫√−g [ −½(∂φ)² −½ m_eff² φ² ] + e^{βφ} L_m ,   φ = γρ
∇²φ = m_eff² φ + β ρ / M_Pl
m_eff²(ρ) = m₀² + α ρⁿ
γ = β/(M_Pl m_eff²)
v_eff²(r) = v_N²(r) · [1 + 2γ ρ(r)]^{−1} · f_cham^{−2}
```
✎ These are standard chameleon forms [8], postulated for TAL. β here is the chameleon matter coupling. ρ_c = 10³ kg m⁻³ is the second ρ_c value.

*"Publikacio 3oldal" (2026 Sept.):* the same law, with the bound restated: "Ha nincs jel 1e-18 szinten, gamma kisebb mint 2.5e-22 m3/kg" [if there is no signal at 1e-18, gamma < 2.5e-22]. It cites Ashby (2003) for GPS [5].

### A.4 V5 — Zenodo 23002701, and the 4-in-1 compilation — Zenodo 23045757 (2026-09-29)

```
f_L(ρ̄_L) = exp(−γ_L ρ̄_L),  γ_L = G_N L²/c²;   weak field f_L ≈ 1 − γ_L ρ̄_L,   Δf/f ≈ γ_L Δρ
Fig. 1: "L0 = 1.3 × 10⁷ m (γ = 1.27 × 10⁻¹³)"
sphere: |Φ|/c² = (4π/3) γ_L ρ̄_L
ρ_tot and ρ_DM as in A.1 (15)–(17)
```
*4-in-1, Part 4 ("Experimental and Computational Verification Proposal", V1.1, 2026-09-29): a momentum-density variant*

```
f_L = exp(-gamma_L * p_bar_L),   0 < f_L <= 1;    d tau_L = f_L * d tau_0
p_bar_L = (1/V) integral_V |rho * v| dV   [kg m^-2 s^-1]
gamma_L >= 0  [m^2 s kg^-1]      (abstract: "[m^3 s kg^-1]")
P1: Delta f / f = exp(-gamma_L * p_bar_L) - 1 ≈ -gamma_L * p_bar_L ;  "gamma_L ~10^-6 to 10^-9 (estimate from V5)"
P2: tau_L = tau_0 / f_L ,  Delta tau/tau_0 ≈ gamma_L * p_bar_L      (muons, tau_0 = 2.2 µs)
P3: dS/d tau_L = (1/f_L) dS/d tau_0  (Second Law preserved for gamma_L >= 0)
```
✎ This variant replaces mass density by momentum density. With p̄ in kg m⁻² s⁻¹, only m² s kg⁻¹ makes γ_L p̄ dimensionless, so the abstract's unit is inconsistent. The variant does not recur in later versions.

✎ Parts 1–3 of the compilation repeat the formulas above. The compilation also mentions "A GPS rendszer napi 38µs korrekciója" [the GPS 38 μs/day correction] as standard background.

### A.5 "v1–v7" complete formalism — Zenodo 23066993 (2026-09-30; PDF header "Version v13 FINAL"; DOI line "23045757 (v1-v6)")

*Definitions*

```
ρ̄_L(x) = 3/(4πL³) ∫_{|x′−x|≤L} ρ(x′) d³x′ = (K_L ∗ ρ)(x)                       (1)
K_L(r) = 3/(4πL³) Θ(L − r),  K̃_L(k) = 3(sin kL − kL cos kL)/(kL)³
γ_L = G_N L²/c²,  [γ_L ρ̄_L] = 1                                                (2)
f_L(ρ̄_L) = exp(−γ_L ρ̄_L),  dτ_L = f_L dt,  ρ̄_L = −ln f_L / γ_L               (7–8)
```

*Bounds on L₀ (postulate: L₀ fixed)*

```
lower: l_P = √(ħG/c³) = 1.616255 × 10⁻³⁵ m,  γ_P = 1.9 × 10⁻⁹⁷
upper: Δf/f < 10⁻¹⁸, Δρ̄ = 10⁴  ⇒  γ < 10⁻²²,  L = √(γc²/G) < 367 m             (3)
favoured: L0 = 10⁻⁴⁹ × 9 × 10¹⁶ / 6.67 × 10⁻¹¹ = 1.16 × 10⁻¹¹ m ≈ 1.2 × 10⁻¹¹ m     (4)
          ⇒ γ0 = 1.07 × 10⁻⁴⁹ m³/kg
```
✎ No radical is visible in the text extraction of Eqs. (3)–(4), but the numerical value corresponds to L₀ = √(γc²/G): √(10⁻⁴⁹ × 9×10¹⁶ / 6.67×10⁻¹¹) ≈ 1.16×10⁻¹¹ m ✓.

*Three heuristic derivations*

```
weak field:  |Φ_L| ≈ G_N ρ̄_L L²,  f_L ≈ 1 + Φ_L/c² ≈ 1 − γ_L ρ̄_L  ⇒ f_L = exp(−γ_L ρ̄_L)      (9–10)
max entropy: S = −Σ p ln p,  ⟨ρ⟩ = ρ̄_L  ⇒  p ∝ exp(−βρ),  γ = βV
Fisher:      I = ∫ (p′)²/p  minimisation  ⇒  exp(−γρ)
```
✎ None of the three fixes γ. The step from "≈ 1 − γρ" to "= exp(−γρ)" is an extrapolation.

*Theorems (Tétel 1–6)*

```
1 f_L > 0    2 df/dρ̄_L = −γ_L f_L < 0    3 0 < f_L ≤ 1    4 f_L(a+b) = f_L(a) f_L(b)
5 lim_{L→0} f_L = 1    6 (causality) dτ_L/dt > 0
```

*Action and field equation*

```
S = ∫ d⁴x √−g [ R/(16πG_N) + f_L(ρ̄_L) L_m ]                                  (11)
G_µν = 8πG_N f_L T_µν + O(δf_L/δg),   T^eff = f_L T                          (12)
H = (π φ̇ − L) f_L ,  stable if f_L > 0
```

*Experiment and astrophysics*

```
2 Yb⁺ clocks, tungsten cube 0.5³ m³, ρ = 19250 kg/m³, M = 2406 kg, Δρ̄ ≈ 10⁴;   Δf/f = γ_L Δρ̄   (13)
error budget RSS 1.1 × 10⁻¹⁸;  height 1 mm ⇒ gΔh/c² = 1.1 × 10⁻¹⁹
Φ_corr = −c² γ ρ̄,  ρ̄ = 10⁻²¹, γ_lab = 10⁻⁴⁹  ⇒  Φ_corr = −9 × 10⁻⁵⁴ J/kg  ⇒ 10⁻³⁴     (14)
neutron star: L = 10 km, γ = 7.41 × 10⁻²⁰, γρ = 0.0296, f = 0.9708 ⇒ 2.9 %          (15)
              L < 367 m, γ = 10⁻²², γρ = 4 × 10⁻⁵, f = 0.99996 ⇒ 0.004 %            (16)
TOV: dm/dr = 4πr² ρ f_L;   ρ_crit = 9 × 10⁻²⁷ ⇒ γρ = 9 × 10⁻⁷⁶
Popper: L = (Δf/f) · c²/(G Δρ)
"modell önmagában nem váltja ki sötét anyagot" [the model alone does not replace dark matter]
```
✎ From Δf/f = G L²Δρ/c², the Popper relation should read L = √[(Δf/f)·c²/(GΔρ)]; the radical is not visible in the text extraction. Theorem 6 is the causality condition invoked in §3.5.

### A.6 v8.1 — Zenodo 23101489 (2026-10-02): "Unified Effective Density Framework, Idokamra 1.0"

*Mass-equivalent density terms (as printed, ASCII)*

```
rho_B  = B^2/(2 mu0 c^2)
rho_E  = E/(c^2 V_eff)          [also: rho_E = epsilon0*E^2/(2*c2)]
rho_GR = (c^2/(8 pi G))*(R - 4 Lambda)
rho_eff = sum of the baryonic, dark-matter, magnetic, curvature and field terms,
          with rho_b, rho_DM weighted by M_rho and rho_B, rho_GR, rho_E weighted by M_B
```
✎ ρ_GR as printed is dimensionally correct: [c²R/G] = kg m⁻³ ✓. This repairs v1 Eq. (15). ρ_E = E/(c²V_eff) is a density ✓.

*Multipliers*

```
"M_rho ~1.2, M_B ~1.1 dimenziotlan kalibracios faktorok" [dimensionless calibration factors]
fractal: M_rho = (L/L0)^(Df-3)
```

*Substitutions and examples*

```
u_B = B^2/(2*mu0) [J/m3] ;  rho_B = u_B/c2 ;  B = 1.5 T → rho_B = 9.95e-12 kg/m3
E = 5e5 V/m → rho_E ~ 1.2e-17 kg/m3
galaxy: v_c^2/r = G*M/r2 → M = v_c^2*r/G ;  rho_tot(r) = v_c^2/(4 pi G r^2)
        e.g. Milky Way r = 8 kpc, v = 220 km/s → 9.47e-22 kg/m3 ; 10-galaxy average 2.47e-22 kg/m3
"Galaxis rho nem 1/r2 -> modell bukik" [if galaxy rho is not ∝ 1/r², the model fails]
```
✎ ρ_tot = v_c²/(4πGr²) is the density of a singular isothermal sphere, i.e. the total dynamical density including any dark matter. The quoted 9.95×10⁻¹² kg m⁻³ at 1.5 T is consistent with B²/(2µ₀c²) ≈ 4.4×10⁻¹² kg m⁻³ at 1 T, scaled by 1.5² ✓.

### A.7 v8.2 — Zenodo 23135605 and repository v8.2.1.pdf (2026-10-04)

```
R ≡ ΔΣ/Δµ = (α_B + α_M)/(α_B + 2 α_M),   r = α_M/α_B  ⇒  R(r) = (1 + r)/(1 + 2r)
"Not (r+1)/(r+2)":  r = 2 → 3/5 = 0.6 (correct) vs 3/4 = 0.75 (old)
inversion: r = (1 − R)/(2R − 1),   δr = δR/(2R − 1)²;   singular at R = 0.5
"2/3 ↔ α_M = α_B ≠ 0 ⇔ X G_3X − 2 G_4φ = 0"  (subclass G5 = 0, G4 = G4(φ))
model ansatz (Sec 4.2):
   α_B(a) = α_0 · f_L(ρ_L(a)) · S(φ)
   α_M(a) = α_0 · [f_L(ρ_L(a)) + (1 − f_L) · q(ρ_L)] · S(φ),   q → 1 as ρ_L → 0
low-density expansion:  R(a) = 2/3 + (1/9)(1 − r) + O((1 − r)²)
   r(a) = 1 − κ γ_L ρ_L(a) + …  ⇒  R(a) = 2/3 + (κ γ_L/9) ρ_L(a)
γ_L = G_N L²/c²  (v8.2.1 prints "γ_L = G * N_L^2 / c^2", with N_L called dimensionless)
```
✎ The expansion is correct (dR/dr|_{r=1} = −1/9). κ here is a dimensionless expansion coefficient, the only κ in the sources. As printed, a dimensionless N_L would give γ_L the unit m kg⁻¹, not the claimed m³ kg⁻¹ ✗. The dimensionally correct form is G L²/c² with L a length. The Horndeski relation is cited from Bellini & Sawicki [6] and Audet Palau (note not supplied; v16: "not verifiable here").

*derivation.md — Zenodo 23137284 (2026-10-04)*

```
α_M(ρ_L) = α_0 exp(−γ_L ρ_L);   α_B(ρ_L) = α_0 exp(−γ_L ρ_L)(1 + ε(ρ_L)),  ε → 0 as ρ_L → 0
ρ_L(a) = ρ_L0 a⁻³ ⇒ f_L(a) = exp(−γ_L ρ_L0 a⁻³);  toy r(a) = 1 + δ_r exp(−k a)
G_eff = G_N(1 + Δµ),  G_lens = G_N(1 + ΔΣ)
"χ² improvement with NEW formula: 7.7x (mock) to 16.8x (real SPARC)"
```
✎ No data or code were supplied for the 7.7×–16.8× claim, so it is author-reported and not verified. v8.3 itself states that R cannot be formed from SPARC.

### A.8 v8.3 — Zenodo 23138388 and repository v8.3_.pdf (2026-10-04)

```
M*² = 2 G4;   α_M = 2 φ̇ G4φ/(H M*²);   H M*² α_B = 2 φ̇ [X G3X − G4φ]
r = α_M/α_B = G4φ/(X G3X − G4φ);   1 + ε = 1/r  ⇒  ε = X G3X/G4φ − 2    "QED"
locus: ε = 0 ⇔ X G3X = 2 G4φ ⇔ r = 1 ⇔ R = 2/3
"f_L kiesik" [f_L drops out]:  r = α0 f_L / [α0 f_L (1+ε)] = 1/(1+ε)
dR/dr = −1/(1+2r)² < 0;   ε > 0 → r < 1 → R from above
background: ρ_L → ρ_m(a) = ρ_m0 a⁻³;  R(a) = (1 + r(a))/(1 + 2r(a)) binned in a
"SPARC … nincs per-galaxis lensing → ΔΣ nincs → R nem képezhető SPARC-ból" [no ΔΣ → R cannot be formed from SPARC]
```
✎ Derived algebraically within the cited Horndeski subclass. Because f_L cancels, the exponential law plays no role in this observable.

### A.9 v9, v10, v11 — no files

No Zenodo record, PDF or repository file is labelled v9, v10 or v11. The repository README at commit f6d18e9 is headed "Time-Attraction Law v8.3 FINAL (v11)", and v8.2.1 mentions a "Live Lab v9.0.0" web page. **Gap: no formulas can be attributed to v9–v11.**

### A.10 v12.1 — Zenodo 23148414 (2026-10-05) and v13 — Zenodo 23193655 (2026-10-06)

```
R(a) = (1 + r(a))/(1 + 2 r(a)),  dR/dr = −1/(1+2r)² < 0                       (3)
r(a) = (1 − R)/(2R − 1) = G4_φ/(X G3X − G4_φ) = 1/(1 + ε(a))                  (4)
ε(a) = X G3X/G4_φ − 2                                                        (5)
locus: ε = 0 ⇔ X G3X = 2 G4_φ ⇔ r = 1 ⇔ R = 2/3                               (6)
ρ_m(a) = ρ_m0 a⁻³,  a = 1/(1+z)
v13: α_B = 2 φ̇ [X G3X − G4_φ] / (H M²)   "(corrected, full denominator)"
v12 intro: "lab prediction Delta f/f > 1e-18 for Delta rho = 1e4 kg/m3"
```
✎ The R(a) and ε algebra is identical to v8.3. The background functions M_i(R(a)) used in the density sum are not specified.

### A.11 v14.0 — Zenodo 23217189 (2026-10-07)

```
Postulate 3: Effective gravitational attraction = Newton * C(rho)
C(rho) = 1 + C0 * exp(-rho / rho_c)
F_TAL = G m1 m2 / r^2 * C(rho)
z_TAL = integral_0^d (1 - X(l)) dl / lambda0 ,   X = 1/C   (definition of X as in v14.1)
"TAL achieves median chi2_red ~1.2 without dark halos"
```
✎ C(ρ) is explicitly a postulate. No numerical C₀ or ρ_c is given in v14.0 (C ~ 6 at ρ ~ 10⁻²⁴ implies C₀ ≈ 5). The χ²_red ~ 1.2 claim is unverified and was not reproduced (§3.4, §4.7).

### A.12 v15 — Zenodo 23220318 (attached file is "Time-Attraction_v14.1.pdf") and GitHub README (2026-10-07)

```
Postulate 1: a dimensionless variable X ∈ [0,1] obeys dX/dl = −X · rho_eff/(rho_c · L0),   with X = 1/C
Postulate 2: 1 + z = 1/X integrated
rho_local = (1/V_Jeans) ∫_{|x′−x|<0.3 kpc} rho_baryon dV
rho_integrated = (1/L0) ∫_0^d rho(l) dl,  L0 = 1 kpc = 3.086e19 m   "kg/m2 / m = kg/m3"
rho_eff = rho_local + rho_integrated
C(rho_eff) = 1 + C0 exp(−rho_eff/rho_c),  C0 = 5 dimensionless,  rho_c = 1e-23 kg/m3
C(0) = 6.0, C(1e-24) = 5.524, C(1e-23) = 2.839, C(5e-23) = 1.0337, C(1e-22) = 1.00023
dC/drho = −(C0/rho_c) exp(−rho/rho_c) < 0;   max boost √6 = 2.449
F_TAL = G m1 m2/r² · (1 + 5 exp(−rho_eff/1e-23))
v_bar² = Y_disk v_disk² + Y_bul v_bul² + v_gas²,  Y_disk = 0.5, Y_bul = 0.7;   rho ~ v_bar²/(G r)/hz,  hz = 0.2 kpc
v_TAL² = C(rho_eff) v_bar²;   chi2_red = chi2/(N − 2)
"Fit result: … C0=5.0 ±0.4, rho_c=(1.0 ±0.2)e-23 kg/m3, median chi2_red=1.18"
z_TAL = ∫_0^d (1 − 1/C(rho(l)))/lambda0 dl,  lambda0 = 1.3e26 m;   D_L = (1+z) ∫ dl · C(l)
H(a) = H0 √(Omega_m a⁻³ C(a)⁻¹),  Omega_m = 1.0;   "H0_local = √C · H0_global ≈ 73"
R = 1 + 0.03 exp(−((z − 0.65)/0.3)²)   (Euclid/DESI prediction)
README v15: "rho_eff = rho_local + rho_integrated with ∫ rho(r′)/|r−r′| dr′";  "r = rho_eff/rho_c";  "Locus R = 2/3"
```
✎ The "derivation" X = exp(−∫ρ dl/(ρ_c L₀)) = 1/C is invalid (§3.4). The fit is unreproduced; "v14.1-data.zip" does not exist. The README's ∫ρ(r′)/|r−r′| dr′ is a density only if read as a one-dimensional integral (with d³r′ it has the unit kg m⁻¹), and it differs from the 1/L₀ form of v14.1. Using r = ρ_eff/ρ_c conflicts with r = α_M/α_B.

### A.13 v16.0 — Zenodo 23234825 (2026-10-08)

```
q(x1,x2) = [rho_local(x1) + rho_local(x2)]/(2 rho_c) + (1/(L0 rho_c)) ∫_{x1}^{x2} rho_b(l) dl     (Eq 1)
rho_eff = rho_c · q,   C(r) = 1 + C0 exp(−q) = 1 + 5 exp(−q)                                  (Eq 2)
Fp,2 = −G m1 m2 C(r) n / r²,   A = G m1 m2                                                     (Eq 3)
U0 = −A C/r  ⇒  F0 − Fp = (A/r) ∇2 C   (≠ 0 unless ∇C = 0)
Ur(r) = −A ∫_r^∞ C(s)/s² ds,   −dUr/dr = −A C(r)/r² = Fp,r                                      (Eq 4)
rho_disk(R,0) = Ydisk · SBdisk(R) · M_sun/pc² / (2 hz)                                         (Eq 5)
Sigma_gas(r) = 1.33 M_HI exp(−r/Rg)/(2π Rg²),  Rg = 2 R_disk                                   (Eq 6)
rho_bulge proxy = 3 Ybulge Vbul² · 1e6 / (4π G r_m³)
z(d) = ∫_0^d [1 − 1/C(l)] dl / lambda0,   lambda0 = 1.3e26 m                                   (Eq 7–8)
D_L = (1+z) ∫_0^d C(l) dl                                                                       (Eq 9)
centred sphere: v_p² = GM C/r,  v0² = GM [C − rC′]/r,  κ_ep² = GM [C + rC′]/r³
constants: C0 = 5, rho_c = 1e-23, L0 = 1 kpc, hz = 0.2 kpc, Ydisk = 0.5, Ybulge = 0.7, lambda0 = 1.3e26 m
"Force multiplier ansatz remains unless derived from field dynamics."
NGC3198 r = 30.08 kpc: v_bar = 70.317 km/s (Ydisk 0.5), C_required = (146/70.317)² = 4.311
```
✎ The pair functional restores F₁ = −F₂ by construction. U_r makes the radial force conservative. The photometric density is "explicitly assumption". (Symbol κ_ep here denotes the epicyclic frequency, written κ in v16.)

*Laboratory pitch (same record)*

```
f_L = exp(−γ_L ρ̄_L),  γ_L = G_N L²/c²
10 t Pb (11.34 g/cm³) at 0.5 m vs 3.0 m from an Sr clock:  predicted Δf/f ~ 5 × 10⁻¹⁹
falsified "if Δf/f < 1 × 10⁻¹⁹ after 72 h averaging"
```
✎ The pitch does not state the L that produces 5×10⁻¹⁹. With L = 1 m, γ_L Δρ ≈ 7.4×10⁻²⁸ × 10⁴ ≈ 10⁻²³, so the stated amplitude is unexplained.

### A.14 v17 — Zenodo 23479940; repository v17_package (2026-10-19): current formulas

```
Axioms: f(0) = 1; f monotonically decreasing; f(ρ1 + ρ2) = f(ρ1) f(ρ2); f continuous
        → (Cauchy) f(ρ) = exp(−γ ρ)   "independent of the scale of γ_L and of V² = V_bar² · C"
γ_L = G_N L²/c²;   ρ_eff = baryonic (disk + bulge) density proxy
"The phenomenological application of the TAL factor is outside the scope of this audit."
split: SHA256("sparc-holdout-20261010:" + name);  123 blind galaxies, 2596 radii, 30 holdout
L0 = 2.36 kpc (median R_disk of the 123);  23 galaxies with T ≥ 10: median R_disk 0.83 kpc
σ = max(0.05 · Vobs, 20 km/s);  χ² = Σ[(Vpred − Vobs)/σ]²;  χ²_red = χ²/(N − 1)
V_bar = √(Vgas² + Vdisk² + Vbul²);  models: V_bar, 0.95 V_bar, 1.05 V_bar, √(Vgas² + Vdisk²), √(Vdisk² + Vbul²), |Vdisk|
"γL = 10²³ conflicts with the lab";  "scale-dependent γL nulls the galactic effect";
"none of the models is better than pure Vbar on the locked sample"
```
✎ The v17 formulas are the current, "academy-corrected" content. On the "none better" sentence, see §4.4.

### A.15 Synopsis: how the central quantities changed

| Quantity | v1–v3 | v1–v7 | v8.1–v13 | v14–v16 | v17 |
|---|---|---|---|---|---|
| Base law | f_L = exp(−γ_Lρ̄_L) | same + action | ρ_eff as a sum of mass-equivalent terms; R(r) | V² = C(ρ)V_bar² | f = exp(−γρ) (axiomatic); baselines only |
| L | 1 m; 1.3×10⁷ m | L₀ fixed, 1.2×10⁻¹¹ m | N_L "dimensionless" (v8.2) | L₀ = 1 kpc | 2.36 kpc (median R_disk) |
| ρ_c | 10⁻¹⁸ (v3); 10³ (HTML) | — | — | 10⁻²³ kg m⁻³ | not used |
| Multipliers | M_eff = M e^{2γρ} | — | M_ρ ≈ 1.2, M_B ≈ 1.1, (L/L₀)^{D_f−3}; κ | C₀ = 5 | none |
| Dark matter | "not replaced" | "not replaced" | ρ_DM ≈ 90 % of ρ_eff | "eliminated" (claim) | "not claimed" |

## Appendix B. Computations made for this manuscript

All additional numbers in Tables 2–4 and §4.1 come from the following scripts. They read only the v17 reference outputs (which are byte-identical to a fresh run), the SHA-256-verified SPARC Rotmod files, and the SPARC catalogue `SPARC_Lelli2016c.mrt`:

- `tal_v17_extra/extra_stats.py`: per-model distribution statistics; point-level checks (all 2596 σ values equal to the floor, maximum V_obs = 383 km s⁻¹); the RAR on the same points (Eq. 11; a₀ = 1.2×10⁻¹⁰ m s⁻²; the runner's parsing including its first-radius rule); paired Wilcoxon and sign tests (SciPy).
- `tal_v17_extra/dwarf_stats.py`: catalogue cross-match (Q, inclination, T, R_disk, V_flat); dwarf subsets; paired tests.
- Outputs: `extra_stats.json`, `dwarf_stats.json`, `rar_per_galaxy.csv`.

Checksums of the scripts and outputs are listed in `tal_v17_extra/SHA256SUMS`. Nothing in the git repository was modified.
