# SPDX-License-Identifier: CC-BY-NC-4.0
# Copyright (c) 2026 Szabolcsi David (ORCID 0009-0005-6300-074X). Commercial use only with the author's written permission.
"""Every quantitative claim of the article: id -> source file/key in the module archives, formula, format, unit."""

CLAIMS = []


def c(id, src, fmt, expr="a", unit="", desc=""):
    if not isinstance(src, list):
        src = [src]
    CLAIMS.append(dict(id=id, src=src, fmt=fmt, expr=expr, unit=unit, desc=desc))


REF = "reference/reference_results.json"
SET = "inputs/setup_declared.json"
BRJ = "outputs/bren_continuation_results.json"
BRV = "inputs/report_values.json"
GAP = "outputs/check_hust_free_gap.json"


def ck(*k):
    return ("clock", REF, "json", list(k))


def cs(*k):
    return ("clock", SET, "json", list(k))


def br(*k):
    return ("bren", BRJ, "json", list(k))


def bv(*k):
    return ("bren", BRV, "json", list(k))


def sc(f, *k):
    return ("subcm", f, "json", list(k))


PS = ["P1", "P2", "P3"]

# ---------------- branch, notation (section 2) ----------------
c("qAlpha", ck("parameter_sets", 0, "alpha"), ("sig", 6), desc="alpha = alpha_Y/(4 pi) for alpha_Y = 1 (q = 1/(4 pi))")
c("qCheck", ck("parameter_sets", 0, "alpha"), ("sig", 6), expr="1/(4*math.pi) if abs(a-1/(4*math.pi))<1e-15 else float('nan')",
  desc="q = 1/(4 pi) recomputed; equals the reference alpha of P1")
for i, P in enumerate(PS):
    c(f"L{P}", ck("parameter_sets", i, "L_m"), ("int",), expr="a/1000", unit="km", desc=f"{P}: L")
    c(f"aY{P}", ck("parameter_sets", i, "alpha_Y"), ("frac",), desc=f"{P}: alpha_Y")
    c(f"alpha{P}", ck("parameter_sets", i, "alpha"), ("sig", 4), desc=f"{P}: alpha = alpha_Y/(4 pi)")
    c(f"aprime{P}", ck("parameter_sets", i, "alpha_prime"), ("sig", 4), desc=f"{P}: alpha' = alpha/(1+alpha kappa)")
    c(f"Ginf{P}", ck("parameter_sets", i, "Ginf_over_Glab"), ("sig", 6), desc=f"{P}: G_inf/G_lab at s_cal")

# ---------------- calibration and setup (section 3) ----------------
c("scal", cs("calibration", "s_cal_m"), ("fix", 2), unit="m", desc="declared calibration separation s_cal")
c("scalLo", cs("calibration", "s_cal_variants_m", 0), ("fix", 2), unit="m", desc="s_cal sensitivity variant (low)")
c("scalHi", cs("calibration", "s_cal_variants_m", 1), ("fix", 2), unit="m", desc="s_cal sensitivity variant (high)")
c("Glab", cs("constants_SI", "G_lab_m3_kg_s2"), ("sci", 5), unit="m3 kg-1 s-2", desc="G_lab (CODATA 2018)")
c("rIn", cs("source_mass", "r_in_m"), ("fix", 2), unit="m", desc="granite shell inner radius")
c("rOut", cs("source_mass", "r_out_m"), ("fix", 2), unit="m", desc="granite shell outer radius")
c("hSrc", cs("source_mass", "height_m"), ("fix", 2), unit="m", desc="granite shell height")
c("rhoSrc", cs("source_mass", "rho_kg_m3"), ("int",), unit="kg m-3", desc="granite density")
c("Msrc", ck("masses", "S1", 0), ("fix", 1), unit="kg", desc="source mass (quadrature)")
c("AB", cs("clocks", "B_xyz_m", 0), ("fix", 2), unit="m", desc="A-B separation")
c("xSone", cs("source_positions", "S1_axis_xy_m", 0), ("fix", 2), unit="m", desc="S1 source axis x")
c("xStwo", cs("source_positions", "S2_axis_xy_m", 0), ("fix", 2), unit="m", desc="S2 source axis x")
c("dmin", cs("source_positions", "clock_to_surface_m"), ("fix", 2), unit="m", desc="clock to outer surface (minimum distance)")
c("DyGR", ck("GR", "Dy_GR"), ("sci", 4), desc="GR reference Delta y (S1-S2)")
c("DaGR", ck("GR", "Da_GR_A"), ("sci", 4), unit="m s-2", desc="GR reference Delta a at A (S1-S2)")
c("nPanel", ck("n_per_panel"), ("int",), desc="Gauss-Legendre nodes per panel")
c("quadLo", cs("lambda_scan", "quadrature_log10L_min"), ("int",), desc="log10 L_min of the quadrature grid (m)")
c("quadHi", cs("lambda_scan", "quadrature_log10L_max"), ("int",), desc="log10 L_max of the quadrature grid (m)")
c("boundLo", cs("lambda_scan", "bound_only_log10L_min"), ("int",), desc="log10 L_min of the bound-only grid (m)")
c("perDec", cs("lambda_scan", "quadrature_points_per_decade"), ("int",), desc="quadrature scan points per decade")
c("perDecB", cs("lambda_scan", "bound_only_points_per_decade"), ("int",), desc="bound-only scan points per decade")

# ---------------- results (section 5) ----------------
for i, P in enumerate(PS):
    c(f"dy{P}", ck("parameter_sets", i, "dy_cal"), ("sci", 4), desc=f"{P}: calibrated delta y (TAL-GR, S1-S2)")
    c(f"da{P}", ck("parameter_sets", i, "da_cal_A"), ("sci", 4), unit="m s-2", desc=f"{P}: calibrated delta a at A")
    if P == "P1":
        c(f"dyRaw{P}", ck("parameter_sets", i, "dy_raw"), ("sci", 3), desc=f"{P}: raw (uncalibrated) delta y, illustration")
        c(f"daRaw{P}", ck("parameter_sets", i, "da_raw_A"), ("sci", 3), unit="m s-2", desc=f"{P}: raw delta a, illustration")
    c(f"dyRel{P}", ck("parameter_sets", i, "dy_over_DyGR"), ("sci", 3), desc=f"{P}: delta y / Delta y_GR")
    c(f"daRel{P}", ck("parameter_sets", i, "da_over_DaGR"), ("sci", 3), desc=f"{P}: delta a / Delta a_GR")
    for tag, s in (("Lo", "0.05"), ("Hi", "0.2")):
        c(f"dy{tag}{P}", ck("parameter_sets", i, "s_cal_variants", s, "dy_cal"), ("sci", 4), desc=f"{P}: delta y at s_cal = {s} m")
        c(f"da{tag}{P}", ck("parameter_sets", i, "s_cal_variants", s, "da_cal_A"), ("sci", 4), unit="m s-2", desc=f"{P}: delta a at s_cal = {s} m")
    b = ("parameter_sets", i, "bounds")
    c(f"aLimLab{P}", ck(*b, "alpha_lim_lab_compilation"), ("sci", 3), desc=f"{P}: alpha_lim, lab compilation")
    c(f"RLab{P}", ck(*b, "R_lab_compilation"), ("sig", 4), desc=f"{P}: R, lab compilation")
    for v, vn in (("A", "A_spherical_raw"), ("B", "B_spherical_zero_sum"), ("C", "C_planar_raw")):
        c(f"aLimBren{v}{P}", ck(*b, "bren_B89_rho2000", vn), ("sci", 3), desc=f"{P}: alpha_lim BREN B89 rho 2000 variant {v}")
    c(f"RBrenMin{P}", ck(*b, "R_bren_min"), ("fix", 2), desc=f"{P}: smallest BREN R (variants A-C)")
    c(f"RBrenMax{P}", ck(*b, "R_bren_max"), ("fix", 2), desc=f"{P}: largest BREN R (variants A-C)")
    c(f"RCons{P}", ck(*b, "R_conservative"), ("fix", 2), desc=f"{P}: conservative R (smallest of lab and BREN)")
    r = ("parameter_sets", i, "resolution")
    c(f"SyOcxo{P}", ck(*r, "clock", "OCXO_8663_short_term|eps_rel_granite", "S"), ("sci", 2), desc=f"{P}: S_y, OCXO, eps granite")
    c(f"SySr{P}", ck(*r, "clock", "Sr_differential_Bothwell2022|eps_rel_precision", "S"), ("sci", 2), desc=f"{P}: S_y, Sr differential, eps 1e-4")
    c(f"SaG{P}", ck(*r, "accel", "torsion_balance_Schlamminger2008|eps_rel_granite", "S"), ("sci", 2), desc=f"{P}: S_a, torsion, eps granite")
    c(f"SaP{P}", ck(*r, "accel", "torsion_balance_Schlamminger2008|eps_rel_precision", "S"), ("sci", 2), desc=f"{P}: S_a, torsion, eps 1e-4")
    c(f"maxS{P}", ck(*r, "max_S"), ("sci", 3), desc=f"{P}: largest S over the declared scenarios")
    c(f"cat{P}", ck("decision", "all_sets", P), ("str",), desc=f"{P}: decision category")
    n = ("parameter_sets", i, "nonuniversal")
    c(f"dyKone{P}", ck(*n, "dy_extra_for_k1"), ("sci", 4), desc=f"{P}: addon delta y_extra for k_c = 1 (Appendix C)")
    c(f"kcNeed{P}", ck(*n, "k_c_needed_5sigma_Bothwell"), ("sci", 4), desc=f"{P}: addon k_c for S = 5 at the Sr differential level")
c("aLimKN", ck("parameter_sets", 2, "bounds", "alpha_lim_K90"), ("sig", 2), desc="K90 report limit at 1 km (1 sigma, correlated chi2), quoted")
c("RKN", ck("parameter_sets", 2, "bounds", "R_K90"), ("fix", 1), desc="P3: R against the quoted K90 limit")
c("aLimBNsp", ck("parameter_sets", 2, "bounds", "alpha_lim_B90_single_point"), ("sci", 3), desc="B90 conditional single-point 1 sigma limit at 1 km (context)")
c("RBNsp", ck("parameter_sets", 2, "bounds", "R_B90_single_point_conditional"), ("fix", 1), desc="P3: R against the B90 single-point value (context only)")
c("sigOcxo", cs("clock_uncertainty", "OCXO_8663_short_term", "value"), ("sci", 1), desc="OCXO 8663 short-term stability")
c("sigBeat", cs("clock_uncertainty", "OCXO_measured_beat_offset", "value"), ("sci", 3), desc="V20 measured mean beat offset (offset, not noise)")
c("sigSrSys", cs("clock_uncertainty", "Sr_systematic_Aeppli2024", "value"), ("sci", 2), desc="Sr single-clock systematic (Aeppli 2024)")
c("sigSrDiff", cs("clock_uncertainty", "Sr_differential_Bothwell2022", "value"), ("sci", 2), desc="Sr differential (Bothwell 2022)")
c("sigTors", cs("acceleration_uncertainty", "torsion_balance_Schlamminger2008", "value"), ("sci", 2), unit="m s-2", desc="torsion-balance differential acceleration (Schlamminger 2008)")
c("epsG", cs("model_uncertainty_scenarios", "eps_rel_granite"), ("sig", 3), desc="declared model-uncertainty proxy (granite densities)")
c("epsP", cs("model_uncertainty_scenarios", "eps_rel_precision"), ("sci", 1), desc="declared best-case model uncertainty")
c("Sthr", cs("resolution_threshold_sigma"), ("int",), desc="feasibility threshold S")
c("MinFactor", [cs("resolution_threshold_sigma"), ck("parameter_sets", 0, "resolution", "max_S"),
                ck("parameter_sets", 1, "resolution", "max_S"), ck("parameter_sets", 2, "resolution", "max_S")],
  ("sci", 4), expr="a/max(b,c,d)", desc="factor by which the calibrated signal falls short of S = 5 (best declared scenario, P1-P3)")
c("catRem", ck("decision", "alphaY1_allowed_remainder"), ("str",), desc="category of the bound-allowed alpha_Y = 1 remainder")
c("catBand", ck("decision", "alphaY1_transition_band_upper_edge"), ("str",), desc="category at the transition-band upper edge")
c("lamEdge", ck("edge", "lambda_edge_m"), ("fix", 2), expr="a*1e6", unit="um", desc="allowed edge (HUST 2020 |alpha| crossing), clock package")
c("logDyEdge", ck("edge", "log10_Dy_bound"), ("fix", 1), desc="log10 bound on |delta y| at the edge")
c("logDaEdge", ck("edge", "log10_Da_bound"), ("fix", 1), desc="log10 bound on |delta a| (m s-2) at the edge")
c("logDyPlus", ck("edge", "transition_band", "thesis_plus_alpha95", "log10_Dy_bound"), ("fix", 1), desc="log10 bound on |delta y| at 101.03 um")
c("logDaPlus", ck("edge", "transition_band", "thesis_plus_alpha95", "log10_Da_bound"), ("fix", 1), desc="log10 bound on |delta a| at 101.03 um")
c("Rfifty", ck("edge", "R_at_50um"), ("sig", 3), desc="R at 50 um (drawn |alpha| envelope)")
c("ImpFifty", ck("edge", "R_at_50um"), ("sig", 2), expr="1/a", desc="sensitivity gain needed at 50 um for R = 1")
c("lamFifty", ck("edge", "R_at_50um"), ("int",), expr="50.0 if a>0 else 0", unit="um", desc="wavelength of the R_at_50um key (label)")
c("xiTi", ck("parameter_sets", 0, "nonuniversal", "xi_Ti"), ("sci", 4), desc="Coulomb-fraction proxy Ti")
c("xiPt", ck("parameter_sets", 0, "nonuniversal", "xi_Pt"), ("sci", 4), desc="Coulomb-fraction proxy Pt")
c("etaNF", ck("parameter_sets", 0, "nonuniversal", "eta95"), ("sci", 4), desc="MICROSCOPE |eta| 95% proxy")
c("kcMax", ck("parameter_sets", 0, "nonuniversal", "k_c_max_long_range_proxy"), ("sci", 4), desc="long-range |k_c| proxy (not transferred to 1-10 km)")
c("nTests", ("clock", "RUN_LOG_r4.md", "re_count", r"^\s+PASS  T\d\d "), ("int",), desc="clock physics self-tests passed (r4 run log)")
c("nCtrl", ("clock", "RUN_LOG_r4.md", "re_count", r"^\s+PASS \(control failed as required\)"), ("int",), desc="clock broken controls that failed as required")
c("rtolClock", ("clock", "RUN_LOG_r4.md", "re", (r"rtol=1e-0?(\d+)\)", 1, "int")), ("int",), desc="clock reference tolerance exponent (rtol = 1e-8)")

# ---------------- BREN chain (section 6.1, Appendix A) ----------------
c("brAlpha", br("model_point", "alpha"), ("sig", 6), desc="BREN model point alpha")
c("brLam", br("model_point", "lambda_m"), ("int",), expr="a/1000", unit="km", desc="BREN model point lambda")
c("zNull", br("z0_m"), ("fix", 1), unit="m", desc="reference height z0 (K90)")
c("zKEN", bv("kasameyer1989", "tower_top_m"), ("int",), unit="m", desc="K89 top height (B89 point)")
c("zBN", br("bounds", "B90", "z_m"), ("fix", 2), unit="m", desc="B90 point height (K90 Fig. 1 top point)")
c("hKN", bv("kammeraad1990", "tower_height_m"), ("int",), unit="m", desc="K90 tower height")
c("Rtrunc", bv("kammeraad1990", "R_trunc_m"), ("int",), expr="a/1000", unit="km", desc="survey truncation radius")
c("nRings", br("base_config", "n_rings"), ("int",), desc="rings in the base layout (synthetic radii)")
c("nStations", br("base_config", "n_stations"), ("int",), desc="stations in the base layout")
c("nGrid", br("base_config", "info", "n"), ("int",), desc="DEM grid cells per side")
c("cellM", br("base_config", "info", "cell_m"), ("int",), unit="m", desc="DEM cell size")
c("nPrisms", br("base_config", "info", "n_prisms"), ("int",), desc="signed terrain prisms (base)")
c("aEarth", ("bren", "methodology_EN.md", "re", (r"a = (\d+) km", 1, "int")), ("int",), unit="km", desc="Earth radius a in the operator")
c("towerLat", bv("tower_location", "lat"), ("fix", 5), unit="deg N", desc="tower latitude")
c("towerLon", bv("tower_location", "lon"), ("fix", 5), expr="-a", unit="deg W", desc="tower longitude (west)")
for rho in (1800, 2000, 2200, 2670):
    c(f"rho{rho}", ("bren", BRJ, "jkey", (["lambda_scan", "A_spherical_raw", "10000"], f"rho_{rho}")), ("int",),
      expr="a" , unit="kg m-3", desc=f"density value {rho} (key rho_{rho})")
    CLAIMS[-1]["src"] = [("bren", BRJ, "re", (rf'"rho_({rho})"', 1, "int"))]
c("kRes", bv("kasameyer1989", "residual_top_uGal"), ("int",), unit="uGal", desc="K89 residual 93 uGal (sign not stated)")
c("kSig", bv("kasameyer1989", "residual_sigma_uGal"), ("int",), unit="uGal", desc="K89 uncertainty 95 uGal")
c("kSigma", ("bren", "inputs/v28_step27_reference.json", "json", ["k_sigma"]), ("fix", 2), desc="k = 1.96")
c("edgeBEN", br("bounds", "B89", "edge_uGal"), ("fix", 1), unit="uGal", desc="B89 edge -93 - 1.96*95")
c("yBN", ("bren", "inputs/kammeraad1990_fig1_method1_digitized.json", "json", ["Y_uGal", -1]), ("fix", 1), unit="uGal", desc="K90 Fig. 1 method 1 top point (digitized)")
c("sBN", ("bren", "inputs/kammeraad1990_fig1_method1_digitized.json", "json", ["sigma_uGal", -1]), ("fix", 1), unit="uGal", desc="K90 Fig. 1 top point sigma (digitized)")
c("edgeBN", br("bounds", "B90", "edge_uGal"), ("fix", 1), unit="uGal", desc="B90 edge")
c("nbias", br("newton_terrain_bias_check", "direct_minus_continued_uGal"), ("fix", 1), unit="uGal", desc="Newtonian terrain bias, direct - continued, top, rho 2200")
c("nbiasRep", br("newton_terrain_bias_check", "report_value_uGal"), ("int",), unit="uGal", desc="K90 terrain bias value")
c("nbiasRepSig", br("newton_terrain_bias_check", "report_sigma_uGal"), ("int",), unit="uGal", desc="K90 terrain bias uncertainty")
c("nbiasKEN", br("newton_terrain_bias_check", "kasameyer1989_rho1800_uGal_at_454"), ("fix", 1), unit="uGal", desc="same at 454 m, rho 1800")
c("nbiasKENrep", br("newton_terrain_bias_check", "kasameyer1989_report_value_uGal"), ("int",), unit="uGal", desc="K89 value (about 100 uGal)")
c("wsum", br("weights_top_minus_z0", "sum"), ("sig", 3), desc="raw weight sum with 240 km truncation")
TOP = "[r for r in a if r['z_m']==457.15][0]"
for key, nm, v in (("flat_uGal", "topFlat", "A_spherical_raw"), ("continuation_uGal", "topA", "A_spherical_raw"),
                   ("continuation_uGal", "topB", "B_spherical_zero_sum"), ("operator_only_uGal", "topOpA", "A_spherical_raw"),
                   ("terrain_part_uGal", "topTerr", "A_spherical_raw")):
    c(nm, br("tables_lambda10km", v), ("fix", 1), expr=TOP + f"['rho_2200']['{key}']", unit="uGal",
      desc=f"top 457.15 m, rho 2200, lambda 10 km: {key} ({v})")
c("topTerrPct", br("tables_lambda10km", "A_spherical_raw"), ("fix", 1),
  expr=f"100*{TOP}['rho_2200']['terrain_part_uGal']/abs({TOP}['rho_2200']['flat_uGal'])", unit="%",
  desc="terrain part relative to |flat| at the top")
VAR = (("A", "A_spherical_raw"), ("B", "B_spherical_zero_sum"), ("C", "C_planar_raw"), ("D", "D_spherical_raw_linear_bias"))
for rho in (2000, 2200, 2670):
    for bnd in ("B89", "B90"):
        t = "BEN" if bnd == "B89" else "BN"
        c(f"Rflat{t}{rho}", br("lambda_scan", "A_spherical_raw", "10000", f"rho_{rho}", bnd, "R_flat"), ("fix", 3),
          desc=f"R flat (z0 = 9.8 m), {bnd}, rho {rho}, lambda 10 km")
        for v, vn in VAR:
            c(f"R{v}{t}{rho}", br("lambda_scan", vn, "10000", f"rho_{rho}", bnd, "R_cont"), ("fix", 3),
              desc=f"R continuation variant {v}, {bnd}, rho {rho}, lambda 10 km")
c("RvNullTwo", br("v28_reproduction", "rho_2000", "R"), ("fix", 3), desc="V28 Step 27 R (z0 = 0), rho 2000, reproduced")
c("RvNullThree", br("v28_reproduction", "rho_2670", "R"), ("fix", 3), desc="V28 Step 27 R (z0 = 0), rho 2670, reproduced")
c("RvNullDiff", [br("v28_reproduction", "rho_2000", "max_rel_diff"), br("v28_reproduction", "rho_2670", "max_rel_diff")],
  ("int",), expr="max(a,b)", desc="largest relative difference of the V28 reproduction")
for lam in ("300", "1000", "3000", "10000", "50000"):
    c(f"lamS{lam}", ("bren", BRJ, "jkey", (["lambda_scan", "A_spherical_raw"], lam)), ("sig", 3), expr="a/1000", unit="km",
      desc=f"lambda scan point {lam} m")
    c(f"aFlat{lam}", br("lambda_scan", "A_spherical_raw", lam, "rho_2000", "B89", "alpha_lim_flat"), ("sig", 3),
      desc=f"alpha_lim flat, B89, rho 2000, lambda {lam} m")
    for v, vn in VAR[:2]:
        c(f"a{v}{lam}", br("lambda_scan", vn, lam, "rho_2000", "B89", "alpha_lim_cont"), ("sig", 3),
          desc=f"alpha_lim continuation {v}, B89, rho 2000, lambda {lam} m")
c("aKN", br("report_alpha_limit_check_1km", "report_1sigma_correlated_chi2"), ("sig", 2), desc="K90 limit at 1 km (1 sigma, correlated chi2), BREN input")
c("aBNsp", br("report_alpha_limit_check_1km", "this_work_single_top_point_1sigma_rho2200_contA"), ("sci", 3), desc="B90 single-point 1 sigma at 1 km, variant A, rho 2200")
CONV = br("convergence_lambda10km")
c("nConv", CONV, ("int",), expr="len(a)-1", desc="convergence configurations besides the base")
c("convBmax", CONV, ("fix", 2), expr="max(abs(v['contB_uGal_rho2200']-a['base']['contB_uGal_rho2200']) for v in a.values())",
  unit="uGal", desc="largest change of variant B over all configurations")
c("convAsixty", CONV, ("fix", 1), expr="a[[k for k in a if k.startswith('Rtrunc_')][0]]['contA_uGal_rho2200']", unit="uGal", desc="variant A at R_trunc 60 km")
c("convAfour", CONV, ("fix", 1), expr="a[[k for k in a if k.startswith('Rtrunc_')][-1]]['contA_uGal_rho2200']", unit="uGal", desc="variant A at R_trunc 480 km")
c("rtSixty", CONV, ("int",), expr="float([k for k in a if k.startswith('Rtrunc_')][0][7:-2])", unit="km", desc="smallest R_trunc of the convergence configurations (key name)")
c("rtFour", CONV, ("int",), expr="float([k for k in a if k.startswith('Rtrunc_')][-1][7:-2])", unit="km", desc="largest R_trunc of the convergence configurations (key name)")
c("convNmin", CONV, ("fix", 1), expr="min(v['newtonbias_uGal_rho2200'] for v in a.values())", unit="uGal", desc="smallest Newtonian bias over configurations")
c("convNmax", CONV, ("fix", 1), expr="max(v['newtonbias_uGal_rho2200'] for v in a.values())", unit="uGal", desc="largest Newtonian bias over configurations")
SS = ("bren", "outputs/check_station_surfaces.json", "json", ["configs"])
c("nSurf", SS, ("int",), expr="len(a)", desc="configurations in check_station_surfaces (base + 10)")
c("nSurfMis", SS, ("int",), expr="sum(v['mismatch_count'] for v in a.values())", desc="station-surface mismatches, all configurations")
SP = ("bren", "outputs/check_single_point_marginal.json", "json", ["rows"])
c("nCov", SP, ("int",), expr="len(a)", desc="covariance structures tried in check_single_point_marginal")
c("nCovSame", SP, ("int",), expr="sum(1 for r in a if r['unchanged'])", desc="structures with an unchanged single-point edge")
c("metricPct", ("bren", "outputs/check_projection_metric.json", "json", ["grids", "301", "max_abs_rel_change"]), ("fix", 2),
  expr="100*a", unit="%", desc="equirectangular east-west metric change, base grid")
ST = ("self_tests",)
c("stPrism", br(*ST, "prism_vs_brute_rel"), ("sci", 2), desc="analytic prism vs brute force")
c("stYnear", br(*ST, "yukawa_near_vs_brute_rel_lam1km"), ("sci", 2), desc="Yukawa near field vs brute force")
c("stHarm", br(*ST, "harmonic_point_mass_rel"), ("sci", 2), desc="harmonic point mass through the operator")

# ---------------- sub-cm chain (section 6.2, Appendix B) ----------------
c("scAlpha", sc("outputs/results.json", "alpha0"), ("sig", 6), desc="sub-cm alpha0")
c("CL", ("subcm", "outputs/results.json", "re", (r"(\d+)% CL", 1, "int")), ("int",), unit="%", desc="confidence level of the sub-cm bounds")
c("eHust", sc(GAP, "hust_2020_abs_um"), ("fix", 3), unit="um", desc="HUST 2020 |alpha| crossing")
c("eEW", sc(GAP, "eot_wash_2020_abs_um"), ("fix", 3), unit="um", desc="Eot-Wash 2020 |alpha| crossing (drawn curve)")
c("eAbs", sc(GAP, "thesis_abs_alpha95_um"), ("fix", 3), unit="um", desc="Lee thesis |alpha95| crossing")
c("ePlus", sc(GAP, "thesis_plus_alpha95_um"), ("fix", 3), unit="um", desc="Lee thesis +alpha95 crossing")
c("eGap", sc(GAP, "gap_width_um"), ("fix", 2), unit="um", desc="HUST-free gap width")
CR = "outputs/crossings.csv"
for nm, fig, cur in (("Hust", "Lee2020_fig5", "HUST_2020"), ("EW", "Lee2020_fig5", "Eot-Wash_2020"),
                     ("Abs", "Lee2020_thesis_table11_4", "abs_alpha95"), ("Plus", "Lee2020_thesis_table11_4", "plus_alpha95")):
    c(f"lo{nm}", ("subcm", CR, "csv", ({"figure": fig, "curve": cur}, "lambda_um_lo")), ("fix", 2), unit="um", desc=f"{cur} crossing, lower end")
    c(f"hi{nm}", ("subcm", CR, "csv", ({"figure": fig, "curve": cur}, "lambda_um_hi")), ("fix", 2), unit="um", desc=f"{cur} crossing, upper end")
c("eHustSix", ("subcm", CR, "csv", ({"figure": "Lee2020_fig5", "curve": "HUST_2016"}, "lambda_um")), ("fix", 2), unit="um", desc="HUST 2016 |alpha| crossing")
c("eAHN", ("subcm", CR, "csv", ({"figure": "AHN2003_fig5", "curve": "Eot-Wash"}, "lambda_um")), ("int",), expr="round(a)", unit="um",
  desc="AHN 2003 Fig. 5 envelope edge (Eot-Wash 2001)")
c("eFigFour", sc("outputs/results.json", "sources", "AHN2003_fig4_fill", "first_excluded_lambda_m"), ("fix", 2), expr="a*1e3", unit="mm",
  desc="AHN 2003 Fig. 4 excluded-fill left edge")
c("dchi", ("subcm", "inputs/frozen_digitization/Lee2020_thesis_table11_4.csv", "re", (r"Delta chi2 = (\d+\.\d+)", 1, "float")), ("fix", 2),
  desc="Delta chi2 of the +-alpha95 surface (thesis)")
c("nThesis", ("subcm", "inputs/frozen_digitization/Lee2020_thesis_table11_4.csv", "csv", ({}, "#count")), ("int",), desc="rows of thesis Table 11.4")
c("rowPlus", sc(GAP, "row_100um", "plus_alpha95"), ("fix", 4), desc="thesis row 100 um: +alpha95")
c("rowAbs", sc(GAP, "row_100um", "abs_alpha95"), ("fix", 4), desc="thesis row 100 um: |alpha95|")
c("rowLam", sc(GAP, "row_100um", "lambda_um"), ("int",), unit="um", desc="thesis row wavelength")
RT = "outputs/ratio_table.csv"
for lam in ("50", "100", "200", "1000"):
    c(f"rL{lam}", ("subcm", RT, "csv", ({"source": "Lee2020_thesis +alpha95", "lambda_um": lam}, "lambda_um")), ("int",), unit="um", desc=f"ratio-table wavelength {lam} um")
    c(f"rPlus{lam}", ("subcm", RT, "csv", ({"source": "Lee2020_thesis +alpha95", "lambda_um": lam}, "R")), ("sig", 3), desc=f"R (+alpha95, thesis) at {lam} um")
    c(f"rEnv{lam}", ("subcm", RT, "csv", ({"source": "envelope combined (|alpha|)", "lambda_um": lam}, "R")), ("sig", 3), desc=f"R (combined |alpha| envelope) at {lam} um")
VC = "outputs/validation_checks.csv"
c("nVal", ("subcm", VC, "csv", ({}, "#count")), ("int",), desc="sub-cm validation checks")
c("nValOk", ("subcm", VC, "csv", ({"result": "consistent"}, "#count")), ("int",), desc="sub-cm validation checks marked consistent")
c("valEW", ("subcm", VC, "csv", ({"check": "Lee2020 Fig5 Eot-Wash 2020 curve: lambda at |alpha|=1 [um]"}, "digitized")), ("fix", 2), unit="um",
  desc="Eot-Wash 2020 curve: lambda at |alpha| = 1, digitized")
c("valEWst", ("subcm", VC, "csv", ({"check": "Lee2020 Fig5 Eot-Wash 2020 curve: lambda at |alpha|=1 [um]"}, "stated")), ("fix", 1), unit="um",
  desc="same, stated in Lee 2020")

# ---------------- reproducibility (section 4.3) ----------------
LOGC = ("article", "logs/fresh_run_clock.txt")
c("pyVer", ("article", "logs/fresh_run_clock.txt", "re", (r"^python (\S+)", 1, "str")), ("str",), desc="Python version of the fresh runs")
c("npVer", ("article", "logs/fresh_run_clock.txt", "re", (r"^numpy (\S+)", 1, "str")), ("str",), desc="NumPy version")
c("spVer", ("article", "logs/fresh_run_clock.txt", "re", (r"^scipy (\S+)", 1, "str")), ("str",), desc="SciPy version")
c("mplVer", ("article", "logs/fresh_run_clock.txt", "re", (r"^matplotlib (\S+)", 1, "str")), ("str",), desc="matplotlib version")
c("pilVer", ("article", "logs/fresh_run_subcm.txt", "re", (r"^Pillow (\S+)", 1, "str")), ("str",), desc="Pillow version (sub-cm)")

# cross-module equalities checked by check_claims.py (value ids that must agree)
EQUAL = [("qAlpha", "brAlpha"), ("qAlpha", "scAlpha"), ("lamEdge", "eHust"), ("aLimKN", "aKN"), ("aLimBNsp", "aBNsp"),
         ("RBrenMinP1", "RABEN2000"), ("rho2000", "rho2000")]
