TAL COMPACT-STAR CONDITIONAL SENSITIVITY STUDY
Package version: prepared for Zenodo V21 (GitHub V19), 23 October 2026.
Licence: CC BY-NC 4.0. Rights held by David Szabolcsi.

Read note_HU_EN.md first: all results are conditional on density, geometry, EOS
and phenomenological coupling. Backlog task 50aa2db0 is unresolved.
Section 7 of the note is the V21 addendum (HU and EN): CSV-driven NICER
overlay, causal-onset densities and the turning-point stability heuristic.
Section 7.4 lists the fixes made after the Academia audit of this V21 package
(H1-H4 and endpoint/core-boundary tests) and the physics items that remain open.

REPRODUCTION
Run Python against these scripts by their extracted paths. They write outputs
next to themselves; no network or external EOS tables are needed.

1. python3 stellar_tal.py
2. python3 diagnostics.py

Equivalent split run (used for this V21 build). Run these five commands from the
extracted directory; the three --eos commands are independent and may be run in
separate shells at the same time, but --combine must wait for all three:
python3 stellar_tal.py --eos SLy
python3 stellar_tal.py --eos APR4
python3 stellar_tal.py --eos H4
python3 stellar_tal.py --combine
python3 diagnostics.py
Afterwards the *_SLy/_APR4/_H4 intermediate tables may be deleted.

Requires NumPy, SciPy, pandas, matplotlib. The versions of this run are recorded in
numerical_diagnostics.json (Python 3.13.5, NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.0,
matplotlib 3.11.2). OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 is advised.
Rerunning the unchanged earlier scripts in this environment reproduced the earlier
outputs numerically, not byte-for-byte (relative Mmax change <=1.7e-7, fixed-mass
radii <=1.4e-6 km); see note section 7.3.

stellar_tal.py options:
  --eos SLy|APR4|H4 : run one EOS and save suffixed intermediate tables.
  --combine        : combine all three EOS runs; requires those intermediates.
                     Also runs the causality/turning-point analysis, WD and plots.
  --wd-only        : recompute Sirius calculations only.
  --plots-only     : redraw M-R figure and radius-shift table from saved data.
  --causal-only    : recompute causal-onset and turning-point columns from saved data.
A default full run replaces base outputs. Run diagnostics.py afterwards to add
causality flags, 66 convergence comparisons, independent geometry checks and (V21
audit fixes) the causal-endpoint/core-boundary refinement tests
(causal_endpoint_convergence.csv) and the NICER symmetrized-error check.

NICER OVERLAY IS CSV-DRIVEN (V21)
plots() starts by reading observational_inputs.csv from the script directory and
draws every NICER ellipse, marker, label, colour and the colour caption from it;
no NICER numbers are hard-coded. To update a measurement, edit the CSV and run
  python3 stellar_tal.py --plots-only
Rows with used_in_plot=True are drawn. Ellipse semi-axes are
sqrt(-2 ln(1-p)) * R_sym_err_km and * M_sym_err_Msun for p = ellipse_joint_prob_1/2,
drawn as an uncorrelated symmetrized Gaussian (sym_err = mean of lower and upper error).
Before drawing, check_observational_inputs() verifies for every row that
M_sym_err_Msun and R_sym_err_km equal the mean of the lower and upper errors
(tolerance 1e-9); a mismatch raises ValueError naming the row, so the figure is not
drawn. PASS/FAIL is printed per source and quantity and stored in
numerical_diagnostics.json. If you update a measurement, update the asymmetric and
the symmetrized errors together.
plot_label is used in the "NICER (GR-derived)" legend (lower-left panel) and in the
colour caption below the figure (with plot_color_name); changing plot_label changes
the figure.
Figure styles: TAL curves are solid and GR curves dashed up to the causal endpoint
(the star whose centre is at the first c_s^2=c^2 density); the continuation up to
the sequence maximum is dotted and faded and is a potentially acausal benchmark of
the published fit. Open circle = TAL causal endpoint, open square = GR causal
endpoint; markers appear only where the sequence crosses the onset. The caption
lists which sequences are affected. (Before these audit fixes, the CSV-driven drawing
code was checked to be byte-identical to the earlier hard-coded figure; the present
figure differs intentionally because of the new styles, labels and markers.)
RA/Dec (J2000) and position epochs: ATNF Pulsar Catalogue v2.6.5,
https://www.atnf.csiro.au/research/pulsar/psrcat/ (Manchester et al. 2005, AJ 129, 1993),
queried 23 October 2026; cross-checked against SIMBAD (https://simbad.cds.unistra.fr).
They are metadata and are not used in the plot.

DATA
stellar_summary.csv : 33 peak summaries and fixed-mass radii, with sound-speed flags.
  V21 columns: log10_rho0c_at_<M> (central density of each fixed-mass star);
  causal_violation_flag, cs2_max_over_c2_in_Mmax_star, cs2_onset_log10_rho0_cgs,
  cs2_onset_log10_eps_over_c2_cgs, cs2_onset_percent_of_rho0c_at_Mmax,
  r_acausal_core_at_Mmax_km (areal radius of the c_s^2=c^2 core boundary, from
  the dense ODE solution at the threshold enthalpy), r_acausal_core_sampled_lower_at_Mmax_km
  (outermost sampled acausal point; a lower estimate only), M_max_causal_Msun, R_at_M_max_causal_km,
  log10_rho0c_at_M_max_causal, causal_at_1p4, causal_at_2p0;
  turning-point heuristic: stable_branch_log10_rho0c_min/max, stable_branch_M_min/max_Msun,
  dMdrho_sign_changes_scan_grid, dMdrho_sign_changes_all_points,
  dMdrho_flips_outside_peak_2e-4dex, dMdrho_flip_max_reversal_dM_Msun,
  stable_branch_at_1p4, stable_branch_at_2p0.
mr_sequences.csv    : 3113 equilibrium points; not all are stable. V21 columns:
  cs2_max_over_c2_all_steps (c_s^2=dP/d eps over every ODE step/output node/EOS
  boundary), acausal_point, r_acausal_core_km (core boundary, as above),
  r_acausal_core_sampled_lower_km (sampled lower estimate).
causal_endpoint_convergence.csv : standard vs finer (rtol 1e-11, ngrid 1601, nq 160)
  recomputation of all 33 causal endpoints and of the core boundary of all 14
  acausal Mmax stars, with c_s^2/c^2-1 at the threshold enthalpy.
causality_stability_log.txt : V21 per-EOS onset densities and per-sequence causal/
  turning-point log.
radius_shifts.csv   : 150 fixed-mass comparisons; missing means unavailable.
observational_inputs.csv : verified NICER marginal intervals, symmetrized errors,
  ellipse levels/alphas, labels/colours, RA/Dec, and plotting selection.
convergence.csv     : 66 standard/finer comparisons, peak and near-1.4 cases.
profiles_L10km.csv  : six representative local density/exponent profiles.
sirius_redshift.csv : conditional additive surface-clock shifts for two polytropes.
wd_summary.json     : adopted WD inputs, propagated residual and illustrative peaks.
numerical_diagnostics.json : versions, benchmarks, independent shell tests,
                            and profile-independent all-enclosed WD band.
geometry_tests.json : uniform-ball and EOS inversion tests.

Units are in column names. The EOS rho0 is rest-mass density, whereas TAL uses
rho=epsilon/c^2. The top-hat uses Euclidean coordinate volume in areal radius.
GR-derived NICER ellipses are independent symmetrized Gaussian joint-68/95%
approximations, not published posterior contours or a TAL likelihood.
Causality: APR4 is superluminal above rho0 = 1.4458e15 g/cm3 in GR and in all
modified branches; SLy above 2.0071e15 g/cm3 (only three mass-source maxima reach it);
H4 only above 1.0461e16 g/cm3 (never reached). All 1.4 and 2.0 Msun comparison
stars are causal. No causal completion was imposed; M_max_causal_Msun is the
endpoint of the published fit's causal range on the followed branch, not the
maximum mass of a causally completed EOS.
Stability: Mmax is an equilibrium-sequence maximum. The dM/drho_c>0 turning-point
criterion is a heuristic for the modified equations, not established radial
stability; a full radial-mode analysis is open.

DOCUMENTATION-ONLY CLARIFICATION (23 October 2026)
The speculative 182/183 point-count explanation in note section 7.4 was replaced
by a verified rule: count unique acausal_point=True rows whose stored central
density is <= the stored peak density of the same EOS/coupling/L sequence.
Compare CSV fields as exact decimal values, without additional rounding or an
inclusion tolerance. This gives 182 rows (173 APR4, 9 SLy); repeated appearances
across figure panels are not counted again. Only note_HU_EN.md, this README and
SHA256SUMS.txt changed in this documentation update. Code, numerical tables,
figures and scientific conclusions are unchanged. No new full run is claimed.

SHA256SUMS.txt records the uncompressed checksums of the supplied files.
