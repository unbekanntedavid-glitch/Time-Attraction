#!/usr/bin/env python3
"""
Step 26 re-audit controls (R1-R6, 5.1-5.3) for the TAL / DESI DR1 BAO package.

This script does NOT run the Block 1 scripts (run_step1_1*.py, run_step1_2_and_1_3.py,
run_step1_4_to_1_6.py). It imports tal_cosmology_common.py, reads the existing CSV/JSON/DESI
files, computes the re-audit controls and
  * writes step26_revision_controls.csv (one row per numeric control),
  * writes step26_r2_piecewise_transition.csv and step26_r5_sne_symmetric_offsets.csv,
  * applies the 5.3 / R2 relabels and the R5 symmetric columns to the existing CSVs,
  * synchronizes no_go_routes_comparison.csv with step2_nogo_inverse_ladder_rd_H0.csv,
  * merges a "step26_revision_controls" section into tal_hubble_desi_bao_summary.json,
    with NaN/inf -> null and json.dump(..., allow_nan=False), then validates strictly.
"""
import os
import json
import math
import shutil

import numpy as np
import pandas as pd
from scipy.integrate import quad

import tal_cosmology_common as tcc

OUT = tcc.OUT_DIR
JSON_PATH = os.path.join(OUT, "tal_hubble_desi_bao_summary.json")

OLD_LATE_ID = "LateTime_Rescaled_72p80_FixedEarly"
NEW_LATE_ID = "Inconsistent_DMoverS_FixedEarly_Illustration"
NEW_LATE_DESC = ("Inconsistent geometric illustration: D_M(z*)/s with unchanged early r_s, r_d "
                 "(no transition redshift; not a consistent late-time model)")
OLD_MT_ID = "Flat_H0_72p80_MatchThetaReint_Om0p2454"


def _nan_to_null(obj):
    if isinstance(obj, dict):
        return {k: _nan_to_null(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_nan_to_null(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return None if (math.isnan(v) or math.isinf(v)) else v
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def _strict_constant(name):
    raise ValueError(f"non-strict JSON constant found: {name}")


# ----------------------------------------------------------------------------- R1
def r1_normalization():
    beta = tcc.BETA_YUK
    k_old_over_k_new = (1.0 / (32.0 * math.pi * beta)) / (1.0 / (8.0 * math.pi * beta))
    # V26 (step23 eq. S_II): kinetic coefficient c^4/(16 pi G alpha_tilde_L), alpha_tilde_L = alpha_Y/(8 pi)
    alpha_tilde_L = tcc.ALPHA_Y / (8.0 * math.pi)
    k_v26_over_k_new = (1.0 / (16.0 * math.pi * alpha_tilde_L)) / (1.0 / (8.0 * math.pi * beta))
    beta_v26 = 2.0 * alpha_tilde_L * tcc.ALPHA_Y ** 2
    # V26 cmb_and_cosmology_audit.py "std" energy density: (phi/4pi)^2 c^2/(8 pi G L^2) vs action-consistent value
    om_phi_new = tcc.OMEGA_PHI_ATTR_0
    om_phi_v26_std = om_phi_new / (4.0 * math.pi) ** 3
    return {
        "kinetic_coefficient_old_draft": "c^4/(32 pi G beta) (Step 26 first draft text)",
        "kinetic_coefficient_new": "c^4/(8 pi G beta)",
        "K_old_over_K_new": k_old_over_k_new,
        "K_V26_over_K_new_at_alpha_Y_1": k_v26_over_k_new,
        "beta_Yuk_V26_at_alpha_Y_1": beta_v26,
        "beta_Yuk_V27": beta,
        "field_redefinition": "phi_old = 2 phi_new requires alpha_new = 2 alpha_old",
        "alpha_DEF_squared_new": tcc.ALPHA_Y ** 2 * beta,
        "alpha_DEF_squared_old_same_alpha": 4.0 * tcc.ALPHA_Y ** 2 * beta,
        "phi_attr_0_new": tcc.PHI_ATTR_0,
        "Omega_phi_attr_0_new": om_phi_new,
        "Omega_phi_attr_0_V26_std_convention": om_phi_v26_std,
        "V26_sources": [
            "scalar-tensor-yukawa-fit/step23_gw170817_and_lensing_en.md:19 (S_II kinetic c^4/(16 pi G alpha_tilde_L))",
            "scalar-tensor-yukawa-fit/step23_gw170817_and_lensing_en.md:29 (alpha_tilde_L = alpha_Y/(8 pi))",
            "scalar-tensor-yukawa-fit/README_en.md:50 (S_II with c^4/(2 G alpha_Y))",
            "scalar-tensor-yukawa-fit/gw170817_and_lensing_audit.py:104-113 (BETA_YUK = 2 alpha_tilde_L alpha_Y^2)",
            "scalar-tensor-yukawa-fit/step25_cmb_and_cosmology_en.md:37 ((Box - L^-2) phi = -(alpha_Y^2 G/c^4) eps_m)",
            "scalar-tensor-yukawa-fit/cmb_and_cosmology_audit.py:327,620 (phi_attr_unit = alpha_Y^2 gamma_L rho)",
            "scalar-tensor-yukawa-fit/cmb_and_cosmology_audit.py:330,623 (std (phi/4pi)^2 convention, 1/(4pi)^3 smaller)",
            "scalar-tensor-yukawa-fit/stellar_tal_v21.py:140-160 (top-hat fixed point x = L^2 (G/c^2) rho_bar)",
            "scalar-tensor-yukawa-fit/data/v21/profiles_L10km.csv (x ~ 0.0412 at r ~ 2.0-2.75 km, interior)",
        ],
        "compact_star_phi_0p0412": ("interior value (r ~ 2.0-2.75 km) of the V21 top-hat fixed-point profile, "
                                    "not a surface value and not a solution of the new Helmholtz-Yukawa action"),
        "factor_2_convention": "2 alpha^2 beta is a sensitivity control, not an equivalent convention",
    }


# ----------------------------------------------------------------------------- R2
def r2_piecewise(z_t_list=(2.5, 10.0, 100.0, 1000.0)):
    s = tcc.H0_LOCAL / tcc.H0_STEP25
    Hm4 = tcc.BG_STEP25_M4["Hz_fn"]
    zs = tcc.Z_STAR_PLANCK
    rs = tcc.RS_ZSTAR
    opts = dict(epsabs=1e-11, epsrel=1e-11, limit=500)
    dm_m4 = tcc.DM_of_z(zs, Hm4)
    rows = []
    for zt in z_t_list:
        low, _ = quad(lambda z: tcc.C_KM_S / (s * Hm4(z)), 0.0, zt, **opts)
        high, _ = quad(lambda z: tcc.C_KM_S / Hm4(z), zt, zs, **opts)
        dm = low + high
        # independent check: single quad of the piecewise H with a break point
        Hpw = lambda z, zt=zt: (s * Hm4(z)) if z < zt else Hm4(z)
        dm_chk, _ = quad(lambda z: tcc.C_KM_S / Hpw(z), 0.0, zs, points=[zt], **opts)
        th = 100.0 * rs / dm
        rows.append({
            "z_t": zt,
            "DM_zstar_Mpc": dm,
            "DM_zstar_check_Mpc": dm_chk,
            "theta_star_100": th,
            "pull_theta_vs_Planck_sigma": (th - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100,
            "H_jump_ratio_at_z_t": s,
            "BAO_bins_inside_rescaled_region": bool(zt > 2.33),
        })
    illus = 100.0 * rs / (dm_m4 / s)
    return {
        "s": s, "r_s_zstar_Mpc": rs, "z_star": zs, "DM_zstar_M4_Mpc": dm_m4,
        "theta_star_100_M4": 100.0 * rs / dm_m4,
        "theta_star_100_DMoverS_illustration": illus,
        "rows": rows,
        "reaudit_reference_table": {"2.5": [13428.2963, 1.075660], "10": [13156.6704, 1.097867],
                                    "100": [12920.1717, 1.117963], "1000": [12843.5875, 1.124630]},
        "note": ("Piecewise H = s H_M4 (z<z_t), H_M4 (z>=z_t) with a jump at z_t: a kinematic control, "
                 "not a dynamical model. D_M/s (100 theta* = 1.124750) has no transition redshift and is "
                 "reclassified as an inconsistent geometric illustration."),
    }


# ----------------------------------------------------------------------------- R4
def r4_joint_split():
    prof = tcc.compute_bao_profile_controls()
    om, h0, rd = prof["joint_best_params"]
    E_fn = lambda z: math.sqrt(om * ((1.0 + z) ** 3) + (1.0 - om))
    u = tcc.get_unit_bao_pred_for_E_fn(E_fn)
    d = (u / (h0 * rd)) - tcc.COBAYA_MEAN
    c2_bao = float(d @ tcc.COBAYA_INV_COV @ d)
    c2_rd = ((rd - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC) ** 2
    c2_h0 = ((h0 - tcc.H0_LOCAL) / tcc.H0_LOCAL_SIGMA) ** 2
    return {
        "best_fit_Omega_m": om, "best_fit_H0": h0, "best_fit_rd_Mpc": rd,
        "chi2_BAO": c2_bao, "chi2_rd_prior": c2_rd, "chi2_H0_prior": c2_h0,
        "chi2_joint": c2_bao + c2_rd + c2_h0,
        "chi2_joint_min_profile": prof["chi2_joint_min"],
        "Delta_chi2_vs_BAO_MLE": prof["dchi2_joint"],
        "sqrt_Delta_chi2": math.sqrt(prof["dchi2_joint"]),
        "rd_prior_sigma_Mpc": tcc.RD_PLANCK_ERR_MPC,
        "rd_error_label": "input Gaussian prior sigma (Planck r_d), not a fitted joint uncertainty",
    }


# ----------------------------------------------------------------------------- R5
def r5_symmetric(df):
    w = 1.0 / df["sigma_mu_rep"].values ** 2
    out = {}
    for key, col in (("M4", "delta_mu_M4_minus_ref"), ("lowOm_0p2448", "delta_mu_Om0p2448_minus_ref")):
        dmu = df[col].values
        off = float(np.sum(w * dmu) / np.sum(w))
        out[key] = {
            "best_offset_mag": off,
            "chi2_without_offset": float(np.sum(w * dmu ** 2)),
            "chi2_with_offset": float(np.sum(w * (dmu - off) ** 2)),
            "n_nodes": int(len(dmu)), "dof_with_offset": int(len(dmu) - 1),
        }
    out["lowest_node_z"] = float(df["z_bin"].iloc[0])
    out["lowest_node_delta_mu_M4_minus_ref"] = float(df["delta_mu_M4_minus_ref"].iloc[0])
    out["identity_5log10_H0ratio"] = 5.0 * math.log10(tcc.H0_LOCAL / tcc.H0_STEP25)
    out["note"] = ("Synthetic shape control: node and error changes do not imply a physical improvement; "
                   "the lowest-node value differs from the pure H0 identity because the reference has a different Omega_m.")
    return out, w


# ----------------------------------------------------------------------------- 5.2
def nu_fd_check(z_list=(0.0, 1.0, 10.0, 100.0, 1000.0, 3000.0)):
    """Fermi-Dirac energy integral for the single massive species vs the p = 1.83 interpolation,
    with both matched to the same relativistic (Omega_nu_rel_1) and non-relativistic (Omega_nu_0) asymptotes."""
    zeta3 = 1.2020569031595942
    F0 = 7.0 * math.pi ** 4 / 120.0
    nr = 1.5 * zeta3

    def F(y):
        v, _ = quad(lambda q: q * q * math.sqrt(q * q + y * y) / (math.exp(q) + 1.0), 0.0, 60.0 + y * 0.0,
                    epsabs=1e-14, epsrel=1e-12, limit=400)
        return v

    r = tcc.OMEGA_NU_0 / tcc.OMEGA_NU_REL_1_0
    y_per_a = r * F0 / nr       # large-y limit F/F0 -> y nr/F0 = r a
    bg = tcc.BG_STEP25_M4
    rows = []
    for z in z_list:
        a = 1.0 / (1.0 + z)
        exact = tcc.OMEGA_NU_REL_1_0 * (1.0 + z) ** 4 * F(y_per_a * a) / F0
        interp = tcc.neutrino_massive_density_E2(z)
        E2 = bg["E2_fn"](z)
        rel_nu = interp / exact - 1.0
        rel_H = 0.5 * (interp - exact) / E2
        rows.append({"z": z, "rho_nu_exact_E2": exact, "rho_nu_interp_E2": interp,
                     "rel_err_massive_nu": rel_nu, "rel_err_H_total": rel_H})
    # largest error over a dense grid
    zz = np.concatenate([np.linspace(0, 10, 101), np.logspace(1, 4, 121)])
    mx_nu, mx_h, z_nu, z_h = 0.0, 0.0, 0.0, 0.0
    for z in zz:
        a = 1.0 / (1.0 + z)
        ex = tcc.OMEGA_NU_REL_1_0 * (1.0 + z) ** 4 * F(y_per_a * a) / F0
        it = tcc.neutrino_massive_density_E2(z)
        e1 = abs(it / ex - 1.0); e2 = abs(0.5 * (it - ex) / bg["E2_fn"](z))
        if e1 > mx_nu: mx_nu, z_nu = e1, float(z)
        if e2 > mx_h: mx_h, z_h = e2, float(z)
    return {
        "interp_power_p": 1.83,
        "z0_normalization_drift_Omega": tcc.NEUTRINO_INTERP_Z0_OFFSET,
        "z0_normalization_drift_relative": tcc.NEUTRINO_INTERP_Z0_OFFSET / tcc.OMEGA_NU_0,
        "rows": rows,
        "max_abs_rel_err_massive_nu": mx_nu, "z_at_max_massive_nu": z_nu,
        "max_abs_rel_err_H_total": mx_h, "z_at_max_H_total": z_h,
        "note": ("1.3546e-7 is the z=0 normalization drift of the interpolation; the Fermi-Dirac "
                 "interpolation error is reported separately here."),
    }


# ----------------------------------------------------------------------------- review 2 (F01-F12)
V21_PROFILE_REL = os.path.join("..", "scalar-tensor-yukawa-fit", "data", "v21", "profiles_L10km.csv")
V21_PROFILE_SHA = "0e4dd1b2600dfbd33bb3a800b39608ed5342d50868cb1d50a4a9a30ab3f74f4f"


def review2_checks(summary, add):
    """Targeted checks for review items F01-F12 (all values recomputed or read from the regenerated outputs)."""
    import hashlib
    out = {}
    q = tcc.ALPHA_Y ** 2 * tcc.BETA_YUK
    # F01 / F02: G dictionary (text + one computed ratio)
    out["F01_frame_coupling"] = {
        "statement": "g_J = A^2 g_E, A = exp(alpha_Y phi): 1/G_T = A^-2/G_*, G_T = G_* A^2; phi < 0 decreases G_T",
        "separated_quantities": ["G_T (Jordan-frame tensor coupling)", "G_lab (operationally measured, ~G_*(1+q))",
                                 "2|phi| (bare magnitude proxy of section 6)"],
    }
    ratio = 1.0 / (1.0 + q)
    assert abs(ratio - 0.9262883) < 5e-8
    out["F02_G_dictionary"] = {
        "G_star": "G in the action (Einstein frame)",
        "G_T": "Jordan-frame tensor coupling, G_T = G_* A^2",
        "G_lab": "G_lab ~ G_*(1+q) for r << L, unscreened small bodies, q = alpha_Y^2 beta_Yuk",
        "G_N_in_growth_equation": "G_N := G_* in the perturbation equation",
        "q": q,
        "G_eff_over_G_lab_large_scale": ratio,
        "note": ("the 1e-41-level growth statement is a formal control at the same bare G; no f sigma8 exclusion "
                 "claim is made and no screening branch is introduced"),
    }
    add("F02", "G_eff_over_G_lab_1_over_1_plus_q", ratio, "", "q = alpha_Y^2 beta_Yuk = 1/(4 pi)")
    # F03: M3 own horizons (stored m_branches must match the review values)
    mb = summary["step2_m1_to_m4_and_inverse_ladder"]["m_branches"]
    ref3 = {"rs_zstar_m3_own_Mpc": 144.36622618, "rd_zdrag_m3_own_Mpc": 147.01159289,
            "theta100_m3": 0.41346227, "chi2_bao_m3_full": 2201.29927}
    for k, v in ref3.items():
        assert abs(mb[k] - v) < 1e-5 * max(1.0, abs(v)), (k, mb[k], v)
        add("F03", k, mb[k], "", f"review value {v}")
    out["F03_M3_own_horizons"] = {k: mb[k] for k in ref3} | {"formal_pull_theta_m3": mb["formal_pull_theta_m3"]}
    # F04: no exclusion wording in the status strings
    df2 = pd.read_csv(os.path.join(OUT, "step2_nogo_inverse_ladder_rd_H0.csv"))
    bad = [w for w in ("excluded by Planck", "Consistent with Planck", "sigma below", "SH0ES")
           for t in df2["status_summary"] if w in t]
    assert not bad, bad
    out["F04_status_labels"] = ["BAO fit", "deviation on illustrative prior", "pull vs fixed target"]
    # F05: assumed local H0 benchmark prior
    tp = summary["tal_fundamental_parameters"]
    if "H0_local_SH0ES_km_s_Mpc" in tp:
        tp["H0_local_benchmark_prior_km_s_Mpc"] = tp.pop("H0_local_SH0ES_km_s_Mpc")
    out["F05_local_H0"] = {"value": tcc.H0_LOCAL, "sigma": tcc.H0_LOCAL_SIGMA,
                           "label": "assumed local H0 benchmark prior",
                           "note": "not the Riess et al. 2022 (ApJL 934, L7) value 73.04 +/- 1.04; no SH0ES release is claimed"}
    # F06 / F08: text items. Their HU/EN wording (and the absence of the superseded wrong sentences) is checked
    # by the real text regression test test_step26_text_regression.py, run by run_all.py after the Markdown
    # generation (nonzero exit on failure; review T01). The two strings below are metadata only.
    out["F06_two_alpha2_beta_row"] = "separate double-strength sensitivity control; not the growth law of the published V26 Step25"
    out["F08_lowest_node_sentence"] = "removed; only the identity 5 log10(72.80/67.40) for identical E(z) is kept"
    # F07: canonical RSD maximum in the legend
    d14 = pd.read_csv(os.path.join(OUT, "step1_4_rsd_fsigma8_growth.csv"))
    mcan = float(d14["delta_fsigma8_var_can_z0norm"].max()); m2x = float(d14["delta_fsigma8_var_2x_z0norm"].max())
    assert abs(mcan / 1.666e-40 - 1.0) < 1e-3
    out["F07_rsd_legend"] = {"max_abs_delta_fsigma8_var_canonical": mcan, "max_abs_delta_fsigma8_var_2x_sensitivity": m2x}
    add("F07", "max_abs_delta_fsigma8_var_canonical", mcan, "", "shown in the RSD figure legend")
    # F09: SNe figure content
    d13 = pd.read_csv(os.path.join(OUT, "step1_3_pantheon_union3_sne_residuals.csv"))
    w = 1.0 / d13["sigma_mu_rep"].values ** 2
    chi_low_sh = float(np.sum((d13["delta_mu_Om0p2448_shifted_minus_ref"] / d13["sigma_mu_rep"]) ** 2))
    out["F09_sne_figure"] = {"N_nodes": int(len(d13)), "dof_after_one_profiled_offset": int(len(d13)) - 1,
                             "curves": ["M4 unshifted", "M4 profiled offset", "low-Om unshifted", "low-Om profiled offset"],
                             "chi2_low_Om_shifted": chi_low_sh,
                             "label": "profiled constant zero-point offset (not an observational calibration)"}
    add("F09", "chi2_low_Om_profiled_offset", chi_low_sh, "", "N=18, nu=17")
    # F10: V21 top-hat interior external benchmark
    d16 = pd.read_csv(os.path.join(OUT, "step1_6_local_yukawa_environments.csv"))
    r7 = d16[d16["env_id"] == "7_V21_TopHat_Interior_External_Benchmark"].iloc[0]
    assert pd.isna(r7["abs_phi_surface_R_gg_L"]) and abs(r7["abs_phi_bulk"] - 0.0412) < 1e-12
    f10 = {"env_id": r7["env_id"], "abs_phi_interior": 0.0412, "surface": None,
           "SLy_gradient": {"r_km": 2.0046409680935633, "x": 0.041239199439254325},
           "APR4_mass": {"r_km": 2.75417748649486, "x": 0.041233030079658145},
           "source_file": "scalar-tensor-yukawa-fit/data/v21/profiles_L10km.csv (V26, read-only)",
           "source_sha256": V21_PROFILE_SHA,
           "Delta_H0_bulk_km_s_Mpc": float(r7["Delta_H0_apparent_bulk_km_s_Mpc"]),
           "Delta_H0_status": "conditional illustration"}
    pth = os.path.normpath(os.path.join(OUT, V21_PROFILE_REL))
    if os.path.exists(pth):
        f10["source_sha256_verified_here"] = hashlib.sha256(open(pth, "rb").read()).hexdigest() == V21_PROFILE_SHA
    else:
        f10["source_sha256_verified_here"] = "source tree not present next to this package (value recorded from the V26 tree)"
    out["F10_v21_interior"] = f10
    # F11: PDF builder exit codes (documented; tested with PATH lacking pdflatex)
    out["F11_pdf_status"] = {"pdf_mode": "nonzero exit if pandoc/pdflatex missing, compile fails or PDF not newly written",
                             "tex_only_mode": "status TEX-ONLY, never counted as a PDF build"}
    # F12: two separated z_* shift diagnostics
    d12 = pd.read_csv(os.path.join(OUT, "step1_2_planck_cmb_rd_consistency.csv"))
    lo = d12[d12["model_id"] == "Flat_H0_72p80_Om0p244828"].iloc[0]
    assert abs(lo["theta_star_100_rs_only_shifted_fixed_DM"] - 1.04248582) < 5e-8, lo["theta_star_100_rs_only_shifted_fixed_DM"]
    assert abs(lo["theta_star_100_both_endpoints_shifted"] - 1.04250287) < 5e-8, lo["theta_star_100_both_endpoints_shifted"]
    out["F12_zstar_shift_diagnostics"] = {
        "rs_only_shifted_fixed_DM": "100 r_s(z_* + dz_s) / D_M(z_* = 1089.80)",
        "both_endpoints_shifted": "100 r_s(z_* + dz_s) / D_M(z_* + dz_s)",
        "low_Om_rs_only_shifted_fixed_DM": float(lo["theta_star_100_rs_only_shifted_fixed_DM"]),
        "low_Om_both_endpoints_shifted": float(lo["theta_star_100_both_endpoints_shifted"]),
        "low_Om_z_star_shifted": float(lo["z_star_shifted_HuSugiyama"]),
        "note": "Hu-Sugiyama fit shift diagnostics only; neither is a recombination solution",
    }
    add("F12", "low_Om_theta100_rs_only_shifted_fixed_DM", float(lo["theta_star_100_rs_only_shifted_fixed_DM"]))
    add("F12", "low_Om_theta100_both_endpoints_shifted", float(lo["theta_star_100_both_endpoints_shifted"]))
    print(f"F02: G_eff/G_lab = 1/(1+q) = {ratio:.7f}; F03: M3 r_s={mb['rs_zstar_m3_own_Mpc']:.8f}, r_d={mb['rd_zdrag_m3_own_Mpc']:.8f}, "
          f"100theta*={mb['theta100_m3']:.8f}, chi2_BAO={mb['chi2_bao_m3_full']:.5f}; F07 max_can={mcan:.4e}; "
          f"F12 low-Om rs-only={lo['theta_star_100_rs_only_shifted_fixed_DM']:.8f}, both={lo['theta_star_100_both_endpoints_shifted']:.8f}")
    return out


def review3_checks(summary, add):
    """Third re-audit (N01, N02, T01)."""
    out = {}
    # N01: rows without a boundary geometry carry no boundary value (null), bulk values unchanged (not halved)
    d16 = pd.read_csv(os.path.join(OUT, "step1_6_local_yukawa_environments.csv"))
    no_bnd = ["1_Cosmic_Background_z0", "2_Local_Cosmic_Void", "3_Milky_Way_Disk_Solar_Nbr"]
    rows = []
    for _, r in d16.iterrows():
        bulk_ref = (tcc.ALPHA_Y * tcc.GAMMA_L * r["rho_kg_m3"]) if not pd.isna(r["rho_kg_m3"]) else 0.0412
        assert abs(r["abs_phi_bulk"] / bulk_ref - 1.0) < 1e-14, r["env_id"]
        if r["env_id"] in no_bnd:
            assert pd.isna(r["abs_phi_surface_R_gg_L"]) and pd.isna(r["Delta_H0_apparent_surf_km_s_Mpc"]), r["env_id"]
            assert r["boundary_status"] == "n.a. (no boundary geometry)", r["env_id"]
        rows.append({"env_id": r["env_id"], "boundary_status": r["boundary_status"],
                     "abs_phi_bulk": float(r["abs_phi_bulk"]),
                     "abs_phi_boundary_R_gg_L": None if pd.isna(r["abs_phi_surface_R_gg_L"]) else float(r["abs_phi_surface_R_gg_L"])})
    out["N01_boundary_column"] = {
        "rule": "boundary value |phi|_bulk/2 only where a half-space boundary geometry is defined; otherwise null / n.a.; "
                "bulk values unchanged (not halved, not copied into the boundary column)",
        "rows_without_boundary_geometry": no_bnd, "rows": rows}
    add("N01", "rows_with_null_boundary_value", float(sum(1 for x in rows if x["abs_phi_boundary_R_gg_L"] is None)), "",
        "FLRW, void, Milky Way disk (no boundary geometry) + V21 interior benchmark")
    # N02: the two distinct, negligible scalar effects quoted in Section 1.7
    r4 = summary["step1_4_rsd_fsigma8_growth"]
    tp = summary["tal_fundamental_parameters"]
    out["N02_scalar_effects"] = {
        "Omega_phi_attr_z0_background_energy": tp["Omega_phi_attr_z0"],
        "delta_Geff_over_Gstar_z0_direct_Yukawa": r4["delta_Geff_over_GN_z0_canonical"],
        "max_abs_delta_fsigma8_variational_canonical": r4["max_abs_delta_fsigma8_variational_canonical"],
        "note": "background scalar-energy correction and the direct Yukawa perturbation of RSD are separate mechanisms"}
    add("N02", "delta_Geff_over_Gstar_z0", r4["delta_Geff_over_GN_z0_canonical"], "", "direct Yukawa perturbation (Section 5)")
    add("N02", "max_abs_delta_fsigma8_var_canonical", r4["max_abs_delta_fsigma8_variational_canonical"], "", "Section 5")
    out["T01_text_regression"] = {"test": "test_step26_text_regression.py", "invoked_by": "run_all.py",
                                  "on_failure": "nonzero exit"}
    print(f"N01: null boundary rows = {[x['env_id'] for x in rows if x['abs_phi_boundary_R_gg_L'] is None]}")
    return out


def main():
    print("=== Step 26 re-audit controls (no Block 1 re-run) ===")
    with open(JSON_PATH, "r", encoding="utf-8") as fh:
        summary = json.load(fh)

    rows = []
    def add(item, key, val, unit="", note=""):
        rows.append({"item": item, "quantity": key, "value": val, "unit": unit, "note": note})

    # R1
    r1 = r1_normalization()
    add("R1", "K_old_over_K_new", r1["K_old_over_K_new"], "", "c^4/(32 pi G beta) vs c^4/(8 pi G beta)")
    add("R1", "K_V26_over_K_new_at_alpha_Y_1", r1["K_V26_over_K_new_at_alpha_Y_1"], "", "V26 S_II kinetic term")
    add("R1", "alpha_DEF_squared_new", r1["alpha_DEF_squared_new"], "", "alpha^2 beta")
    add("R1", "alpha_DEF_squared_old_same_alpha", r1["alpha_DEF_squared_old_same_alpha"], "", "4 alpha^2 beta")
    add("R1", "Omega_phi_attr_0_new", r1["Omega_phi_attr_0_new"], "", "action-consistent")
    add("R1", "Omega_phi_attr_0_V26_std_convention", r1["Omega_phi_attr_0_V26_std_convention"], "", "V26 (phi/4pi)^2 convention")
    print(f"R1: K_old/K_new = {r1['K_old_over_K_new']:.4f}, K_V26/K_new = {r1['K_V26_over_K_new_at_alpha_Y_1']:.6f}")

    # R2
    r2 = r2_piecewise()
    pd.DataFrame(r2["rows"]).to_csv(os.path.join(OUT, "step26_r2_piecewise_transition.csv"), index=False)
    for r in r2["rows"]:
        ref = r2["reaudit_reference_table"][f"{r['z_t']:g}"]
        print(f"R2: z_t={r['z_t']:g}: D_M={r['DM_zstar_Mpc']:.4f} Mpc (check {r['DM_zstar_check_Mpc']:.4f}; re-audit {ref[0]}), "
              f"100theta*={r['theta_star_100']:.6f} (re-audit {ref[1]})")
        assert abs(r["DM_zstar_Mpc"] - ref[0]) < 2e-3 and abs(r["theta_star_100"] - ref[1]) < 2e-6
        add("R2", f"DM_zstar_z_t_{r['z_t']:g}", r["DM_zstar_Mpc"], "Mpc")
        add("R2", f"theta_star_100_z_t_{r['z_t']:g}", r["theta_star_100"], "")
    add("R2", "theta_star_100_DMoverS_illustration", r2["theta_star_100_DMoverS_illustration"], "", "inconsistent illustration")
    print(f"R2: D_M/s illustration 100theta* = {r2['theta_star_100_DMoverS_illustration']:.6f}")

    # R4
    r4 = r4_joint_split()
    for k in ("chi2_BAO", "chi2_rd_prior", "chi2_H0_prior", "chi2_joint"):
        add("R4", k, r4[k])
    print(f"R4: chi2_BAO={r4['chi2_BAO']:.8f}, chi2_rd={r4['chi2_rd_prior']:.8f}, chi2_H0={r4['chi2_H0_prior']:.8f}, "
          f"chi2_joint={r4['chi2_joint']:.8f}")
    assert abs(r4["chi2_joint"] - r4["chi2_joint_min_profile"]) < 1e-8

    # R5
    sne_path = os.path.join(OUT, "step1_3_pantheon_union3_sne_residuals.csv")
    df13 = pd.read_csv(sne_path)
    r5, _w = r5_symmetric(df13)
    off_low = r5["lowOm_0p2448"]["best_offset_mag"]
    df13["mu_Flat_72p8_Om0p2448_shifted"] = df13["mu_Flat_72p8_Om0p2448"] - off_low
    df13["delta_mu_Om0p2448_shifted_minus_ref"] = df13["mu_Flat_72p8_Om0p2448_shifted"] - df13["mu_synth_ref_72p8_Om0p334"]
    df13.to_csv(sne_path, index=False)
    r5rows = []
    for m in ("M4", "lowOm_0p2448"):
        r5rows.append({"model": m, **r5[m]})
        for k in ("best_offset_mag", "chi2_without_offset", "chi2_with_offset"):
            add("R5", f"{m}_{k}", r5[m][k])
    pd.DataFrame(r5rows).to_csv(os.path.join(OUT, "step26_r5_sne_symmetric_offsets.csv"), index=False)
    add("R5", "lowest_node_delta_mu_M4_minus_ref", r5["lowest_node_delta_mu_M4_minus_ref"], "mag")
    add("R5", "identity_5log10_H0ratio", r5["identity_5log10_H0ratio"], "mag")
    print(f"R5: M4 offset={r5['M4']['best_offset_mag']:.8f}, chi2 {r5['M4']['chi2_without_offset']:.5f} -> {r5['M4']['chi2_with_offset']:.6f}; "
          f"lowOm offset={off_low:.8f}, chi2 {r5['lowOm_0p2448']['chi2_without_offset']:.5f} -> {r5['lowOm_0p2448']['chi2_with_offset']:.6f}")

    # 5.2
    nu = nu_fd_check()
    for r in nu["rows"]:
        add("5.2", f"rel_err_massive_nu_z_{r['z']:g}", r["rel_err_massive_nu"])
        add("5.2", f"rel_err_H_total_z_{r['z']:g}", r["rel_err_H_total"])
    add("5.2", "max_abs_rel_err_massive_nu", nu["max_abs_rel_err_massive_nu"], "", f"z={nu['z_at_max_massive_nu']:.3g}")
    add("5.2", "max_abs_rel_err_H_total", nu["max_abs_rel_err_H_total"], "", f"z={nu['z_at_max_H_total']:.3g}")
    add("5.2", "z0_normalization_drift_Omega", nu["z0_normalization_drift_Omega"])
    r1000 = [r for r in nu["rows"] if r["z"] == 1000.0][0]
    print(f"5.2: z=1000 massive-nu rel.err = {r1000['rel_err_massive_nu']:.4e}, H rel.err = {r1000['rel_err_H_total']:.4e}; "
          f"max |nu| = {nu['max_abs_rel_err_massive_nu']:.4e} at z={nu['z_at_max_massive_nu']:.4g}, "
          f"max |H| = {nu['max_abs_rel_err_H_total']:.4e} at z={nu['z_at_max_H_total']:.4g}")

    # 5.3 + R2 relabel in the step1_2 CSV and JSON
    cmb_path = os.path.join(OUT, "step1_2_planck_cmb_rd_consistency.csv")
    df12 = pd.read_csv(cmb_path)
    mt = df12["model_id"].str.startswith("Flat_H0_72p80_MatchThetaReint")
    om_mt = float(df12.loc[mt, "Omega_m_0"].iloc[0]); wcb_mt = float(df12.loc[mt, "omega_cb_h2"].iloc[0])
    new_mt_id = "Flat_H0_72p80_MatchThetaReint_Om0p" + f"{om_mt:.4f}".split(".")[1]
    new_mt_desc = f"Flat H0=72.80 matching re-integrated 100*theta_* (Om={om_mt:.5f}, omega_cb={wcb_mt:.5f})"
    df12.loc[mt, "model_id"] = new_mt_id
    df12.loc[mt, "description"] = new_mt_desc
    lt = df12["model_id"].isin([OLD_LATE_ID, NEW_LATE_ID])
    df12.loc[lt, "model_id"] = NEW_LATE_ID
    df12.loc[lt, "description"] = NEW_LATE_DESC
    df12.to_csv(cmb_path, index=False)
    print(f"5.3: MatchTheta id -> {new_mt_id} (Om={om_mt:.17g}, omega_cb={wcb_mt:.10g})")
    add("5.3", "MatchTheta_Omega_m", om_mt); add("5.3", "MatchTheta_omega_cb_h2", wcb_mt)

    s12 = summary["step1_2_planck_cmb_and_illustrative_3d_prior"]
    s12["models"] = json.loads(df12.to_json(orient="records"))

    # R4 / step2 sync
    s2_path = os.path.join(OUT, "step2_nogo_inverse_ladder_rd_H0.csv")
    df2 = pd.read_csv(s2_path)
    shutil.copyfile(s2_path, os.path.join(OUT, "no_go_routes_comparison.csv"))
    summary["step2_m1_to_m4_and_inverse_ladder"]["inverse_ladder_scenarios"] = json.loads(df2.to_json(orient="records"))
    j3 = summary["step1_1_desi_dr1_bao"]["joint_3param_BAO_PlanckRd_LocalH0"]
    j3.update({"chi2_BAO": r4["chi2_BAO"], "chi2_rd_prior": r4["chi2_rd_prior"],
               "chi2_H0_prior": r4["chi2_H0_prior"], "chi2_joint": r4["chi2_joint"],
               "sqrt_Delta_chi2": r4["sqrt_Delta_chi2"], "rd_error_label": r4["rd_error_label"]})
    j3.pop("equivalent_gaussian_sigma", None)

    # R5 into step1_3 section
    s13 = summary["step1_3_synthetic_distance_modulus_shape_control"]
    s13["symmetric_zero_point_R5"] = r5

    # R6 metadata
    r6 = {
        "status": "assumed order-of-magnitude proxies; no BBN reaction network is run",
        "clock_proxy": ("bare magnitude proxy |delta G/G| ~ 2|phi| (review F01). With g_J = A^2 g_E, A = exp(alpha_Y phi), "
                        "the Jordan-frame Ricci term carries A^-2/G_* (A^-2 multiplies the Planck mass squared): "
                        "1/G_T = A^-2/G_*, i.e. G_T = G_* A^2 = G_* exp(2 alpha_Y phi), so phi < 0 decreases G_T. "
                        "The tensor coupling G_T, the operationally measured G_lab ~ G_*(1+q) and the bare 2|phi| proxy "
                        "are different quantities; sign and reference state are not propagated, only |2 phi| is used"),
        "sensitivity_Yp_per_dH_over_H": 0.08,
        "sensitivity_dlnDH_per_dH_over_H": 1.6,
        "sensitivity_reference": ("Kneller & Steigman 2004 (New J. Phys. 6, 117) fits: dYp ~ 0.16 (S-1), "
                                  "d ln(D/H) ~ 1.6 (S-1), S = H'/H; the code's 1.6 matches, its 0.08 is a factor 2 smaller"),
        "Yp_std": {"value": 0.24709, "source": "PRIMAT, Pitrou, Coc, Uzan & Vangioni 2018, Phys. Rep. 754, 1 (central value)"},
        "DH_e5_std": {"value": 2.5147, "status": "arithmetic proxy (review P02)",
                      "source": "linearized benchmark value without verified single-paper provenance "
                      "(PRIMAT 2018 gives 2.459); not a published BBN network prediction"},
        "Yp_obs": {"value": 0.245, "err": 0.003, "source": "PDG 2024 BBN review (Fields, Molaro & Sarkar), recommended value; "
                   "Aver et al. 2021 (JCAP 03, 027) give 0.2453 +/- 0.0034"},
        "DH_e5_obs": {"value": 2.527, "err": 0.030, "source": "Cooke, Pettini & Steidel 2018, ApJ 855, 102; the CSV field name "
                      "DH_e5_obs_PDG is a legacy label (PDG 2024 recommends 2.547 +/- 0.025)"},
    }
    # step1_5 is a list of rows in the master JSON; the R6 provenance lives in the revision section below.

    v51 = {
        "float64_branch": "geff_ratio = 1 + delta_Geff_canonical (source corrected; the stored float64 difference stays 0 at this precision)",
        "variational_headers": "absolute values |delta f sigma8|",
        "reaudit_independent_rescaled_check": {"max_abs_delta_fsigma8_can": 1.6663681521e-40,
                                               "max_rel_deviation": 5.1e-6,
                                               "source": "quoted from the re-audit; confirms 1.666e-40"},
    }

    # ------------------------------------------------------------------ review 2 (F01-F12)
    f_items = review2_checks(summary, add)
    n_items = review3_checks(summary, add)

    summary["step26_revision_controls"] = {
        "script": "run_step26_revision_controls.py",
        "R1_action_normalization": r1,
        "R2_piecewise_late_transition": r2,
        "R3_package": {"json": "strict (allow_nan=False, NaN -> null)",
                       "runner": "run_all.py", "requirements": "requirements.txt"},
        "R4_joint_chi2_split": r4,
        "R5_sne_symmetric_offsets": r5,
        "R6_bbn_proxy_provenance": r6,
        "item_5_1_growth": v51,
        "item_5_2_neutrino_fermi_dirac": nu,
        "item_5_3_relabel": {"old_id": OLD_MT_ID, "new_id": new_mt_id, "Omega_m": om_mt, "omega_cb_h2": wcb_mt,
                             "late_time_old_id": OLD_LATE_ID, "late_time_new_id": NEW_LATE_ID},
        "review2_F01_F12": f_items,
        "review3_N01_N02_T01": n_items,
    }

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "step26_revision_controls.csv"), index=False)

    summary = _nan_to_null(summary)
    with open(JSON_PATH, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False, allow_nan=False)
    with open(JSON_PATH, "r", encoding="utf-8") as fh:
        json.loads(fh.read(), parse_constant=_strict_constant)
    print(f"Saved strict JSON: {JSON_PATH}")
    print("Saved: step26_revision_controls.csv, step26_r2_piecewise_transition.csv, step26_r5_sne_symmetric_offsets.csv")


if __name__ == "__main__":
    main()
