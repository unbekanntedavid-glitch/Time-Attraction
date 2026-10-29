#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1.2 & Step 1.3:
- Step 1.2: Planck 2018 CMB Sound Horizon r_s(z_*), Drag Scale r_d(z_drag),
  Acoustic Scale 100*theta_*, Self-Consistent Re-Integration vs Late-Time
  Rescaling Controls, and Illustrative 3D Gaussian Prior on (omega_b, omega_cb, 100*theta_*).
- Step 1.3: Synthetic Distance-Modulus mu(z) Shape Control on an 18-Point
  Pantheon+-like Redshift Grid (z = 0.015 - 2.260) with a profiled constant zero-point
  offset (N = 18 nodes, nu = 17 dof after one profiled offset; not an observational calibration).
"""

import os
import math
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tal_cosmology_common as tcc

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_step1_2():
    bg_m4 = tcc.BG_STEP25_M4
    bg_pl = tcc.BG_PLANCK_67P36

    dm_zstar_m4 = tcc.DM_of_z(tcc.Z_STAR_PLANCK, bg_m4["Hz_fn"], 0.0, tcc.H0_STEP25)
    theta100_m4 = 100.0 * (tcc.RS_ZSTAR / dm_zstar_m4)

    zstar_eh_m4 = tcc.z_star_hu_sugiyama(tcc.OMEGA_B_H2, bg_m4["omega_m_h2"])
    zdrag_eh_m4 = tcc.z_drag_eisenstein_hu(tcc.OMEGA_B_H2, bg_m4["omega_m_h2"])

    # 1) Flat H0 = 72.80, Om = 0.244828 benchmark (Review Section 4.4)
    om_flat_2448 = 0.244828
    bg_flat_2448 = tcc.make_model_background(tcc.H0_LOCAL, om_flat_2448, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True)

    # 2) Flat H0 = 72.80 where re-integrated 100*theta_*_reint(fixed z_*) == theta100_m4 (1.041320)
    def _theta_reint_diff_flat(om_val):
        bg_tmp = tcc.make_model_background(tcc.H0_LOCAL, om_val, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True)
        rs_tmp = tcc.compute_sound_horizon_for_bg(tcc.Z_STAR_PLANCK, bg_tmp)
        dm_tmp = tcc.DM_of_z(tcc.Z_STAR_PLANCK, bg_tmp["Hz_fn"], 0.0, tcc.H0_LOCAL)
        return 100.0 * (rs_tmp / dm_tmp) - theta100_m4

    om_flat_728_reint_theta = brentq(_theta_reint_diff_flat, 0.15, 0.35)
    bg_flat_728_reint_theta = tcc.make_model_background(
        tcc.H0_LOCAL, om_flat_728_reint_theta, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True
    )

    # 3) Open H0 = 72.80 with fixed physical omega_cb = 0.14237 (Om = 0.26985) matching DM(z_*) = dm_zstar_m4
    def _dm_diff_open(ok_val):
        bg_tmp = tcc.make_model_background(tcc.H0_LOCAL, None, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, ok_val, exact_z0_norm=True)
        return tcc.DM_of_z(tcc.Z_STAR_PLANCK, bg_tmp["Hz_fn"], ok_val, tcc.H0_LOCAL) - dm_zstar_m4

    ok_open_728 = brentq(_dm_diff_open, -0.01, 0.10)
    bg_open_728 = tcc.make_model_background(tcc.H0_LOCAL, None, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, ok_open_728, exact_z0_norm=True)

    # 4) Standard thermal history at H0 = 72.80 with unchanged Omega_m = 0.314820 (so omega_cb = 0.16610)
    bg_same_om_728 = tcc.make_model_background(tcc.H0_LOCAL, tcc.OMEGA_M_0, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True)

    models_spec = [
        ("Step25_M4_TAL_Baseline", "Step 25 M4 / TAL baseline (H0=67.40, Om=0.31482, flat)", bg_m4, "self_consistent"),
        ("Planck2018_67p36_Ref", "Planck 2018 LCDM reference (H0=67.36, Om=0.31519, flat)", bg_pl, "self_consistent"),
        ("Rigid_AllZ_Rescaled_72p80", "All-z rigid rescaling H(z)->(72.80/67.40)H_M4(z) with unchanged c_s(z)", bg_m4, "all_z_rigid"),
        # Re-audit R2: this row is NOT a late-transition model (there is no transition redshift z_t);
        # it combines D_M(z*)/s with unchanged early r_s, r_d and is kept only as an inconsistent geometric illustration.
        # The piecewise late-transition control H = s*H_M4 (z < z_t), H_M4 (z >= z_t) is in run_step26_revision_controls.py.
        ("Inconsistent_DMoverS_FixedEarly_Illustration", "Inconsistent geometric illustration: D_M(z*)/s with unchanged early r_s, r_d (no transition redshift; not a consistent late-time model)", bg_m4, "late_time_only"),
        ("Flat_H0_72p80_SameOm0p3148", "Standard thermal history at H0=72.80 with unchanged Om=0.31482 (omega_cb=0.16610)", bg_same_om_728, "self_consistent"),
        ("Flat_H0_72p80_Om0p244828", "Flat H0=72.80 low-matter benchmark (Om=0.244828, omega_cb=0.12911)", bg_flat_2448, "self_consistent"),
        # Re-audit 5.3: id and description carry the actually solved Omega_m (about 0.2442), not a stale value.
        ("Flat_H0_72p80_MatchThetaReint_Om0p" + f"{om_flat_728_reint_theta:.4f}".split(".")[1],
         f"Flat H0=72.80 matching re-integrated 100*theta_* (Om={om_flat_728_reint_theta:.5f}, omega_cb={bg_flat_728_reint_theta['omega_cb_h2']:.5f})",
         bg_flat_728_reint_theta, "self_consistent"),
        ("Open_H0_72p80_FixedOmegaCB_OmK0p0098", "Open H0=72.80 with fixed omega_cb=0.14237 and Om_K=+0.00979", bg_open_728, "self_consistent"),
    ]

    rows = []
    scale_s = tcc.H0_LOCAL / tcc.H0_STEP25
    for m_id, desc, bg, mode in models_spec:
        if mode == "all_z_rigid":
            H0_v = tcc.H0_LOCAL
            om_m = tcc.OMEGA_M_0
            om_k = 0.0
            om_b_h2 = tcc.OMEGA_B_H2
            om_cb_h2 = tcc.OMEGA_CB_H2
            Hz_fn = lambda z, f=bg["Hz_fn"]: scale_s * f(z)
            rs_reint = tcc.RS_ZSTAR / scale_s
            rd_reint = tcc.RD_ZDRAG / scale_s
            rs_eh = rs_reint
            rd_eh = rd_reint
            dm_zstar = dm_zstar_m4 / scale_s
            zstar_shift = tcc.Z_STAR_PLANCK
            dm_zstar_shift = dm_zstar
        elif mode == "late_time_only":
            H0_v = tcc.H0_LOCAL
            om_m = tcc.OMEGA_M_0
            om_k = 0.0
            om_b_h2 = tcc.OMEGA_B_H2
            om_cb_h2 = tcc.OMEGA_CB_H2
            Hz_fn = lambda z, f=bg["Hz_fn"]: scale_s * f(z)
            rs_reint = tcc.RS_ZSTAR
            rd_reint = tcc.RD_ZDRAG
            rs_eh = tcc.RS_ZSTAR
            rd_eh = tcc.RD_ZDRAG
            dm_zstar = dm_zstar_m4 / scale_s
            zstar_shift = tcc.Z_STAR_PLANCK
            dm_zstar_shift = dm_zstar
        else:
            H0_v = bg["H0"]
            om_m = bg["Omega_m_0"]
            om_k = bg["Omega_k_0"]
            om_b_h2 = bg["omega_b_h2"]
            om_cb_h2 = bg["omega_cb_h2"]
            Hz_fn = bg["Hz_fn"]
            rs_reint = tcc.compute_sound_horizon_for_bg(tcc.Z_STAR_PLANCK, bg)
            rd_reint = tcc.compute_sound_horizon_for_bg(tcc.Z_DRAG_PLANCK, bg)
            dz_s = tcc.z_star_hu_sugiyama(om_b_h2, bg["omega_m_h2"]) - zstar_eh_m4
            dz_d = tcc.z_drag_eisenstein_hu(om_b_h2, bg["omega_m_h2"]) - zdrag_eh_m4
            rs_eh = tcc.compute_sound_horizon_for_bg(tcc.Z_STAR_PLANCK + dz_s, bg)
            rd_eh = tcc.compute_sound_horizon_for_bg(tcc.Z_DRAG_PLANCK + dz_d, bg)
            dm_zstar = tcc.DM_of_z(tcc.Z_STAR_PLANCK, Hz_fn, om_k, H0_v)
            # Review F12: D_M also evaluated at the shifted endpoint z_* + dz_s ("both endpoints shifted").
            zstar_shift = tcc.Z_STAR_PLANCK + dz_s
            dm_zstar_shift = tcc.DM_of_z(zstar_shift, Hz_fn, om_k, H0_v)

        theta100_reint = 100.0 * (rs_reint / dm_zstar)
        theta100_fixed_rs = 100.0 * (tcc.RS_ZSTAR / dm_zstar)
        # Review F12: two clearly separated Hu-Sugiyama shift diagnostics; neither is a recombination solution.
        #   rs_only_shifted_fixed_DM: r_s(z_* + dz_s) / D_M(z_* = 1089.80)  (only the r_s endpoint shifted)
        #   both_endpoints_shifted:   r_s(z_* + dz_s) / D_M(z_* + dz_s)     (r_s and D_M at the same shifted z_*)
        theta100_eh = 100.0 * (rs_eh / dm_zstar)
        theta100_eh_both = 100.0 * (rs_eh / dm_zstar_shift)

        pull_theta_reint = (theta100_reint - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100
        pull_theta_fixed_rs = (theta100_fixed_rs - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100
        pull_omega_cb = (om_cb_h2 - tcc.PLANCK_3D_MEAN[1]) / tcc.PLANCK_3D_SIGMA[1]
        pull_rd_planck = (rd_reint - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC

        chi2_3d_reint = tcc.planck_3d_illustrative_chi2(om_b_h2, om_cb_h2, theta100_reint)
        chi2_bao_reint_full = tcc.bao_chi2_cobaya_full(Hz_fn, rd_reint, om_k, H0_v)
        chi2_bao_reint_tab1 = tcc.bao_chi2_table1(Hz_fn, rd_reint, om_k, H0_v)
        chi2_bao_fixed_rd_full = tcc.bao_chi2_cobaya_full(Hz_fn, tcc.RD_ZDRAG, om_k, H0_v)

        rows.append({
            "model_id": m_id,
            "description": desc,
            "H0_km_s_Mpc": H0_v,
            "Omega_m_0": om_m,
            "Omega_K_0": om_k,
            "omega_b_h2": om_b_h2,
            "omega_cb_h2": om_cb_h2,
            "rs_zstar_reint_Mpc": rs_reint,
            "rd_zdrag_reint_Mpc": rd_reint,
            "rs_zstar_EH_shifted_Mpc": rs_eh,
            "rd_zdrag_EH_shifted_Mpc": rd_eh,
            "DM_zstar_Mpc": dm_zstar,
            "theta_star_100_reint": theta100_reint,
            "z_star_shifted_HuSugiyama": zstar_shift,
            "DM_zstar_shifted_Mpc": dm_zstar_shift,
            "theta_star_100_rs_only_shifted_fixed_DM": theta100_eh,
            "theta_star_100_both_endpoints_shifted": theta100_eh_both,
            "theta_star_100_fixed_early_rs": theta100_fixed_rs,
            "pull_theta_star_reint_sigma": pull_theta_reint,
            "pull_theta_star_fixed_rs_sigma": pull_theta_fixed_rs,
            "pull_omega_cb_planck_sigma": pull_omega_cb,
            "pull_rd_reint_planck_sigma": pull_rd_planck,
            "chi2_Planck_3D_illustrative_prior": chi2_3d_reint,
            "chi2_DESI_BAO_Cobaya_with_rd_reint": chi2_bao_reint_full,
            "chi2_DESI_BAO_Table1_with_rd_reint": chi2_bao_reint_tab1,
            "chi2_DESI_BAO_Cobaya_with_fixed_rd": chi2_bao_fixed_rd_full,
        })

    df_12 = pd.DataFrame(rows)
    csv_12 = os.path.join(OUT_DIR, "step1_2_planck_cmb_rd_consistency.csv")
    df_12.to_csv(csv_12, index=False)
    return df_12, csv_12

# -----------------------------------------------------------------------------
# Step 1.3: Synthetic Distance-Modulus mu(z) Shape Control (18 Redshift Nodes)
# -----------------------------------------------------------------------------
PANTHEON_GRID = [
    (0.015, 0.042),
    (0.023, 0.038),
    (0.035, 0.035),
    (0.050, 0.032),
    (0.075, 0.030),
    (0.100, 0.028),
    (0.150, 0.026),
    (0.200, 0.025),
    (0.250, 0.025),
    (0.300, 0.026),
    (0.400, 0.028),
    (0.500, 0.030),
    (0.650, 0.033),
    (0.800, 0.036),
    (1.000, 0.042),
    (1.300, 0.055),
    (1.700, 0.075),
    (2.260, 0.110),
]

def distance_modulus_mu(z, Hz_fn, omega_k=0.0, H0_for_curv=tcc.H0_STEP25):
    dm = tcc.DM_of_z(z, Hz_fn, omega_k, H0_for_curv)
    dl_mpc = (1.0 + z) * dm
    return 5.0 * math.log10(dl_mpc) + 25.0

def generate_sne_plots(df_13, delta_mu_zero, best_shift_mag, chi2_unshifted, chi2_shifted, lang="hu",
                       best_shift_low=0.0, chi2_low=0.0, chi2_low_shifted=0.0):
    # Review F09: both models get the same profiled constant zero-point offset; N = 18 nodes,
    # nu = 17 dof after one profiled offset. This is not an observational (Cepheid) calibration.
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11.5,
        "legend.fontsize": 9.0,
        "figure.dpi": 200,
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.0))
    z_arr = df_13["z_bin"].values

    if lang == "hu":
        lbl_synth = "Szintetikus Pantheon+-szeru kontrollgorbe (Om=0.334, H0=72.80)"
        lbl_m4_un = f"Step 25 M4 / TAL (H0=67.40, eltolas nelkul): Delta_mu(0) = +{delta_mu_zero:.4f} mag"
        lbl_m4_sh = (f"Step 25 M4 / TAL, profilozott konstans nullpont-eltolas -{best_shift_mag:.4f} mag:\n"
                     f"chi2 = {chi2_shifted:.2f} (N = 18 csomopont, nu = 17 szabadsagi fok)")
        lbl_245 = f"Alacsony Om=0.2448, H0=72.80 modell, eltolas nelkul: chi2 = {chi2_low:.2f} (N = 18)"
        lbl_245s = (f"Alacsony Om=0.2448, profilozott konstans nullpont-eltolas -{best_shift_low:.4f} mag:\n"
                    f"chi2 = {chi2_low_shifted:.2f} (N = 18 csomopont, nu = 17 szabadsagi fok)")
        t1 = "(a) Szintetikus tavolsagmodulus mu(z) alakkontroll (z = 0.015 - 2.260)"
        t2 = r"(b) Tavolsagmodulus-maradek $\Delta\mu(z) \equiv \mu_{\mathrm{modell}}(z) - \mu_{\mathrm{ref}}(z; 72.80)$"
        sup = ("V27 (Step 26) 1.3: Szintetikus SNe Ia tavolsagmodulus mu(z) alakkontroll (18 csomopont) "
               "profilozott konstans nullpont-eltolassal (nem megfigyelesi kalibracio)")
        fname = "step1_3_pantheon_sne_plot.png"
    else:
        lbl_synth = "Synthetic Pantheon+-like reference curve (Om=0.334, H0=72.80)"
        lbl_m4_un = f"Step 25 M4 / TAL (H0=67.40, unshifted): Delta_mu(0) = +{delta_mu_zero:.4f} mag"
        lbl_m4_sh = (f"Step 25 M4 / TAL, profiled constant zero-point offset -{best_shift_mag:.4f} mag:\n"
                     f"chi2 = {chi2_shifted:.2f} (N = 18 nodes, nu = 17 dof)")
        lbl_245 = f"Low Om=0.2448, H0=72.80 model, unshifted: chi2 = {chi2_low:.2f} (N = 18)"
        lbl_245s = (f"Low Om=0.2448, profiled constant zero-point offset -{best_shift_low:.4f} mag:\n"
                    f"chi2 = {chi2_low_shifted:.2f} (N = 18 nodes, nu = 17 dof)")
        t1 = "(a) Synthetic Distance-Modulus mu(z) Shape Control (z = 0.015 - 2.260)"
        t2 = r"(b) Distance-Modulus Residual $\Delta\mu(z) \equiv \mu_{\mathrm{model}}(z) - \mu_{\mathrm{ref}}(z; 72.80)$"
        sup = ("V27 (Step 26) 1.3: Synthetic SNe Ia Distance-Modulus mu(z) Shape Control (18 Redshift Nodes) "
               "with a Profiled Constant Zero-Point Offset (not an observational calibration)")
        fname = "step1_3_pantheon_sne_plot_en.png"

    ax1.errorbar(z_arr, df_13["mu_synth_ref_72p8_Om0p334"], yerr=df_13["sigma_mu_rep"],
                 fmt="o", color="black", capsize=3, ms=5, label=lbl_synth)
    ax1.plot(z_arr, df_13["mu_TAL_M4_67p4"], color="#1f4e78", lw=2.2, label=lbl_m4_un)
    ax1.plot(z_arr, df_13["mu_TAL_M4_shifted"], color="#2ca02c", lw=2.0, ls="--", label=lbl_m4_sh)
    ax1.set_xscale("log")
    ax1.set_xlabel("z (log skála / log scale)")
    ax1.set_ylabel(r"$\mu(z)$ [mag]")
    ax1.set_title(t1)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="lower right")

    ax2.axhline(0.0, color="black", lw=1.2, ls=":")
    ax2.errorbar(z_arr, np.zeros_like(z_arr), yerr=df_13["sigma_mu_rep"],
                 fmt="o", color="black", capsize=3, ms=4.5, alpha=0.65, label=lbl_synth)
    ax2.plot(z_arr, df_13["delta_mu_M4_minus_ref"], color="#1f4e78", lw=2.2, label=lbl_m4_un)
    ax2.plot(z_arr, df_13["delta_mu_M4_shifted_minus_ref"], color="#2ca02c", lw=2.2, ls="--", label=lbl_m4_sh)
    ax2.plot(z_arr, df_13["delta_mu_Om0p2448_minus_ref"], color="#c00000", lw=1.9, ls="-.", label=lbl_245)
    ax2.plot(z_arr, df_13["delta_mu_Om0p2448_shifted_minus_ref"], color="#ff7f0e", lw=1.9, ls=":", label=lbl_245s)
    ax2.set_xscale("log")
    ax2.set_xlabel("z (log skála / log scale)")
    ax2.set_ylabel(r"$\Delta\mu(z) = \mu_{\mathrm{model}}(z) - \mu_{\mathrm{ref}}(z)$ [mag]")
    ax2.set_title(t2)
    ax2.set_ylim(-0.12, 0.40)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper left", fontsize=7.6)

    fig.suptitle(sup, fontsize=11.8, fontweight="bold", y=1.02)
    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, fname)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path

def run_step1_3():
    bg_synth_ref = tcc.make_model_background(tcc.H0_LOCAL, 0.334, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True)
    bg_m4 = tcc.BG_STEP25_M4
    bg_low_om = tcc.make_model_background(tcc.H0_LOCAL, 0.244828, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True)

    # Pure identity for identical E(z): 5 log10(72.80/67.40) = +0.167357 mag. The actual lowest-node
    # difference of the tabulated models (different Omega_m) is 0.16782511 mag (re-audit R5).
    delta_mu_zero = 5.0 * math.log10(tcc.H0_LOCAL / tcc.H0_STEP25)   # +0.167357 mag

    rows = []
    for z, sig in PANTHEON_GRID:
        mu_ref = distance_modulus_mu(z, bg_synth_ref["Hz_fn"], 0.0, tcc.H0_LOCAL)
        mu_m4 = distance_modulus_mu(z, bg_m4["Hz_fn"], 0.0, tcc.H0_STEP25)
        mu_low = distance_modulus_mu(z, bg_low_om["Hz_fn"], 0.0, tcc.H0_LOCAL)
        rows.append({
            "z_bin": z,
            "sigma_mu_rep": sig,
            "mu_synth_ref_72p8_Om0p334": mu_ref,
            "mu_TAL_M4_67p4": mu_m4,
            "mu_Flat_72p8_Om0p2448": mu_low,
            "delta_mu_M4_minus_ref": mu_m4 - mu_ref,
            "delta_mu_ref_minus_M4": mu_ref - mu_m4,
            "delta_mu_Om0p2448_minus_ref": mu_low - mu_ref,
        })

    df_13 = pd.DataFrame(rows)
    weights = 1.0 / (df_13["sigma_mu_rep"].values ** 2)
    diff_unshifted = df_13["delta_mu_M4_minus_ref"].values
    best_shift_mag = float(np.sum(weights * diff_unshifted) / np.sum(weights))  # +0.175950 mag

    df_13["mu_TAL_M4_shifted"] = df_13["mu_TAL_M4_67p4"] - best_shift_mag
    df_13["delta_mu_M4_shifted_minus_ref"] = df_13["mu_TAL_M4_shifted"] - df_13["mu_synth_ref_72p8_Om0p334"]

    chi2_unshifted = float(np.sum(((df_13["delta_mu_M4_minus_ref"]) / df_13["sigma_mu_rep"]) ** 2))
    chi2_shifted = float(np.sum(((df_13["delta_mu_M4_shifted_minus_ref"]) / df_13["sigma_mu_rep"]) ** 2))
    chi2_low_om = float(np.sum(((df_13["delta_mu_Om0p2448_minus_ref"]) / df_13["sigma_mu_rep"]) ** 2))

    # Re-audit R5: symmetric zero-point treatment -- the low-Om control gets the same profiled constant offset.
    best_shift_low = float(np.sum(weights * df_13["delta_mu_Om0p2448_minus_ref"].values) / np.sum(weights))
    df_13["mu_Flat_72p8_Om0p2448_shifted"] = df_13["mu_Flat_72p8_Om0p2448"] - best_shift_low
    df_13["delta_mu_Om0p2448_shifted_minus_ref"] = df_13["mu_Flat_72p8_Om0p2448_shifted"] - df_13["mu_synth_ref_72p8_Om0p334"]

    csv_13 = os.path.join(OUT_DIR, "step1_3_pantheon_union3_sne_residuals.csv")
    df_13.to_csv(csv_13, index=False)

    chi2_low_shifted = float(np.sum(((df_13["delta_mu_Om0p2448_shifted_minus_ref"]) / df_13["sigma_mu_rep"]) ** 2))
    kw = dict(best_shift_low=best_shift_low, chi2_low=chi2_low_om, chi2_low_shifted=chi2_low_shifted)
    plot_hu = generate_sne_plots(df_13, delta_mu_zero, best_shift_mag, chi2_unshifted, chi2_shifted, lang="hu", **kw)
    plot_en = generate_sne_plots(df_13, delta_mu_zero, best_shift_mag, chi2_unshifted, chi2_shifted, lang="en", **kw)
    return df_13, csv_13, plot_hu, plot_en, delta_mu_zero, best_shift_mag, chi2_unshifted, chi2_shifted, chi2_low_om

def main():
    df_12, csv_12 = run_step1_2()
    df_13, csv_13, p_hu, p_en, dmu0, shift_best, c2_un, c2_sh, c2_low = run_step1_3()
    print("=== STEP 1.2 PLANCK CMB SOUND HORIZON & ILLUSTRATIVE 3D PRIOR SUMMARY ===")
    for _, r in df_12.iterrows():
        print(f"  {r['model_id']:38s}: rs_reint={r['rs_zstar_reint_Mpc']:.4f}, rd_reint={r['rd_zdrag_reint_Mpc']:.4f}, "
              f"100*theta*_reint={r['theta_star_100_reint']:.6f} ({r['pull_theta_star_reint_sigma']:+.2f}s), "
              f"pull_wcb={r['pull_omega_cb_planck_sigma']:+.2f}s, chi2_3D={r['chi2_Planck_3D_illustrative_prior']:.2f}, "
              f"chi2_BAO(rd_reint)={r['chi2_DESI_BAO_Cobaya_with_rd_reint']:.2f}")
    print("\n=== STEP 1.3 SYNTHETIC DISTANCE-MODULUS SHAPE CONTROL SUMMARY ===")
    print(f"  Zero-point offset Delta_mu(0) = mu_M4(67.40) - mu_ref(72.80) = +{dmu0:.6f} mag (reverse: -{dmu0:.6f} mag)")
    print(f"  Weighted best shift across z=0.015..2.260: +{shift_best:.6f} mag")
    print(f"  Synthetic shape control chi2 (unshifted M4 vs 72.80 ref): {c2_un:.4f} / 18 pts")
    print(f"  Synthetic shape control chi2 (zero-point shifted M4 vs Om=0.334 ref): {c2_sh:.4f} / 18 pts")
    print(f"  Synthetic shape control chi2 (Flat H0=72.80, Om=0.244828 unshifted): {c2_low:.4f} / 18 pts")
    print(f"Saved: {csv_12}, {csv_13}, {p_hu}, {p_en}")

if __name__ == "__main__":
    main()
