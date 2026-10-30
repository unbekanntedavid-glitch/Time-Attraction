#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2.1, Step 2.2 & Step 2.3:
- Step 2.1: M1-M4 Vacuum-Energy / Branch Comparison (including exact full-background
  root z_root = 0.295737 vs pure-matter 0.295878 for M2, and formal 1D theta_*
  pull for M3).
- Step 2.2 & 2.3: Analytical & Numerical No-Go Appendix and Inverse Distance
  Ladder (H0, r_d) Degeneracy Analysis (72.8 +/- 1.0 is an assumed local H0 benchmark prior,
  not a sourced SH0ES release; review F05), separating the fixed-P illustration
  (P = 9914.02 km/s -> r_d = 136.18 Mpc) from the profiled BAO + Planck r_d +
  local H0 statistics (P_best = 10193.67 +/- 126.47 km/s, r_d^req = 140.02 +/- 1.74 Mpc,
  H0 = 69.30 +/- 0.87 km/s/Mpc, joint Delta_chi2 = 6.97, sqrt(Delta_chi2) = 2.64;
  the joint chi2 is split into chi2_BAO, chi2_rd_prior and chi2_H0_prior (re-audit R4).
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

def analyze_m1_to_m4_branches():
    # M2 uncompensated flat: Omega_Lambda_eff = -0.685180, Omega_K = 0.0
    om_de_m2 = -tcc.OMEGA_LAMBDA_FLAT_0
    E2_m2_fn = lambda z: tcc.E2_general(z, tcc.OMEGA_CB_0, tcc.OMEGA_NU_0,
                                        tcc.OMEGA_R_MASSLESS_0, tcc.OMEGA_NU_REL_1_0,
                                        0.0, om_de_m2)
    E2_m2_z0 = E2_m2_fn(0.0)  # -0.37019956
    z_root_full = brentq(E2_m2_fn, 0.0, 1.0)  # 0.295737
    z_root_matter_approx = ((tcc.OMEGA_LAMBDA_FLAT_0 / tcc.OMEGA_M_0) ** (1.0 / 3.0)) - 1.0  # 0.295878

    # M3 uncompensated open: Omega_Lambda_eff = -0.685180, Omega_K = +1.370200 so E^2(0) == 1.0
    om_k_m3 = 1.0 - tcc.OMEGA_CB_0 - tcc.OMEGA_R_MASSLESS_0 - tcc.OMEGA_NU_INTERP_Z0 - om_de_m2  # +1.37019956
    E2_m3_fn = lambda z: tcc.E2_general(z, tcc.OMEGA_CB_0, tcc.OMEGA_NU_0,
                                        tcc.OMEGA_R_MASSLESS_0, tcc.OMEGA_NU_REL_1_0,
                                        om_k_m3, om_de_m2)
    Hz_m3_fn = lambda z: tcc.H0_STEP25 * math.sqrt(E2_m3_fn(z))
    dm_zstar_m3 = tcc.DM_of_z(tcc.Z_STAR_PLANCK, Hz_m3_fn, om_k_m3, tcc.H0_STEP25)
    # Review F03: the M3 control uses its OWN sound horizons, re-integrated with the M3 background
    # (same fixed z_*, z_drag; same omega_b, omega_gamma), not the flat-M4 RS_ZSTAR / RD_ZDRAG.
    bg_m3 = {"H0": tcc.H0_STEP25, "Omega_b_0": tcc.OMEGA_B_0, "Omega_gamma_0": tcc.OMEGA_GAMMA_0, "E2_fn": E2_m3_fn}
    rs_m3 = tcc.compute_sound_horizon_for_bg(tcc.Z_STAR_PLANCK, bg_m3)
    rd_m3 = tcc.compute_sound_horizon_for_bg(tcc.Z_DRAG_PLANCK, bg_m3)
    theta100_m3 = 100.0 * (rs_m3 / dm_zstar_m3)
    formal_pull_theta_m3 = (theta100_m3 - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100
    chi2_bao_m3_full = tcc.bao_chi2_cobaya_full(Hz_m3_fn, rd_m3, om_k_m3, tcc.H0_STEP25)
    chi2_bao_m3_tab1 = tcc.bao_chi2_table1(Hz_m3_fn, rd_m3, om_k_m3, tcc.H0_STEP25)
    # Former values with the frozen M4 horizons, kept only as a labelled comparison:
    theta100_m3_m4hor = 100.0 * (tcc.RS_ZSTAR / dm_zstar_m3)
    chi2_bao_m3_full_m4hor = tcc.bao_chi2_cobaya_full(Hz_m3_fn, tcc.RD_ZDRAG, om_k_m3, tcc.H0_STEP25)

    return {
        "E2_m2_z0": E2_m2_z0,
        "z_root_full": z_root_full,
        "z_root_matter_approx": z_root_matter_approx,
        "om_k_m3": om_k_m3,
        "om_de_m3": om_de_m2,
        "rs_zstar_m3_own_Mpc": rs_m3,
        "rd_zdrag_m3_own_Mpc": rd_m3,
        "horizon_definition_m3": "own r_s, r_d re-integrated with the M3 H(z) at fixed z_*=1089.80, z_drag=1059.94 (review F03)",
        "dm_zstar_m3": dm_zstar_m3,
        "theta100_m3": theta100_m3,
        "formal_pull_theta_m3": formal_pull_theta_m3,
        "chi2_bao_m3_full": chi2_bao_m3_full,
        "chi2_bao_m3_tab1": chi2_bao_m3_tab1,
        "theta100_m3_with_frozen_M4_horizon_superseded": theta100_m3_m4hor,
        "chi2_bao_m3_full_with_frozen_M4_horizon_superseded": chi2_bao_m3_full_m4hor,
    }

def joint_chi2_components(prof):
    """Re-audit R4: split the joint (Om, H0, rd) minimum into its BAO, r_d-prior and H0-prior parts."""
    om, h0, rd = prof["joint_best_params"]
    E_fn = lambda z: math.sqrt(om * ((1.0 + z) ** 3) + (1.0 - om))
    u = tcc.get_unit_bao_pred_for_E_fn(E_fn)
    d = (u / (h0 * rd)) - tcc.COBAYA_MEAN
    c2_bao = float(d @ tcc.COBAYA_INV_COV @ d)
    c2_rd = ((rd - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC) ** 2
    c2_h0 = ((h0 - tcc.H0_LOCAL) / tcc.H0_LOCAL_SIGMA) ** 2
    return {"chi2_BAO": c2_bao, "chi2_rd_prior": c2_rd, "chi2_H0_prior": c2_h0, "chi2_joint": c2_bao + c2_rd + c2_h0}

def build_step2_table(prof, m_branches):
    jc = joint_chi2_components(prof)
    # Also evaluate self-consistent Flat H0=72.80, Om=0.244828 and Open H0=72.80, Om_K=+0.009785
    bg_flat_2448 = tcc.make_model_background(tcc.H0_LOCAL, 0.244828, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True)
    rs_2448 = tcc.compute_sound_horizon_for_bg(tcc.Z_STAR_PLANCK, bg_flat_2448)
    rd_2448 = tcc.compute_sound_horizon_for_bg(tcc.Z_DRAG_PLANCK, bg_flat_2448)
    dm_zstar_2448 = tcc.DM_of_z(tcc.Z_STAR_PLANCK, bg_flat_2448["Hz_fn"], 0.0, tcc.H0_LOCAL)
    theta100_2448 = 100.0 * (rs_2448 / dm_zstar_2448)
    pull_theta_2448 = (theta100_2448 - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100
    pull_wcb_2448 = (bg_flat_2448["omega_cb_h2"] - tcc.PLANCK_3D_MEAN[1]) / tcc.PLANCK_3D_SIGMA[1]
    chi2_bao_2448_reint = tcc.bao_chi2_cobaya_full(bg_flat_2448["Hz_fn"], rd_2448, 0.0, tcc.H0_LOCAL)

    dm_zstar_m4 = tcc.DM_of_z(tcc.Z_STAR_PLANCK, tcc.BG_STEP25_M4["Hz_fn"], 0.0, tcc.H0_STEP25)
    theta100_m4 = 100.0 * (tcc.RS_ZSTAR / dm_zstar_m4)
    pull_theta_m4 = (theta100_m4 - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100
    pull_wcb_m4 = (tcc.OMEGA_CB_H2 - tcc.PLANCK_3D_MEAN[1]) / tcc.PLANCK_3D_SIGMA[1]

    def _dm_diff_open(ok_val):
        bg_tmp = tcc.make_model_background(tcc.H0_LOCAL, None, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, ok_val, exact_z0_norm=True)
        return tcc.DM_of_z(tcc.Z_STAR_PLANCK, bg_tmp["Hz_fn"], ok_val, tcc.H0_LOCAL) - dm_zstar_m4

    ok_open_728 = brentq(_dm_diff_open, -0.01, 0.10)
    bg_open_728 = tcc.make_model_background(tcc.H0_LOCAL, None, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, ok_open_728, exact_z0_norm=True)
    rs_open = tcc.compute_sound_horizon_for_bg(tcc.Z_STAR_PLANCK, bg_open_728)
    rd_open = tcc.compute_sound_horizon_for_bg(tcc.Z_DRAG_PLANCK, bg_open_728)
    dm_zstar_open = tcc.DM_of_z(tcc.Z_STAR_PLANCK, bg_open_728["Hz_fn"], ok_open_728, tcc.H0_LOCAL)
    theta100_open = 100.0 * (rs_open / dm_zstar_open)
    pull_theta_open = (theta100_open - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100
    chi2_bao_open_reint = tcc.bao_chi2_cobaya_full(bg_open_728["Hz_fn"], rd_open, ok_open_728, tcc.H0_LOCAL)

    rows = [
        {
            "scenario_id": "1_Step25_M4_TAL_Baseline_Fixed",
            "category": "Fixed_Baseline",
            "H0_km_s_Mpc": tcc.H0_STEP25,
            "Omega_m_0": tcc.OMEGA_M_0,
            "Omega_K_0": 0.0,
            "omega_cb_h2": tcc.OMEGA_CB_H2,
            "rd_used_Mpc": tcc.RD_ZDRAG,
            "rd_err_Mpc": 0.0,
            "P_H0_rd_km_s": prof["P_m4"],
            "P_H0_rd_err_km_s": 0.0,
            "chi2_DESI_BAO_Cobaya_full": prof["chi2_m4_full"],
            "chi2_DESI_BAO_Table1": prof["chi2_m4_tab1"],
            "theta_star_100_reint": theta100_m4,
            "pull_theta_star_reint_sigma": pull_theta_m4,
            "pull_omega_cb_planck_sigma": pull_wcb_m4,
            "pull_rd_vs_Planck1D_sigma": (tcc.RD_ZDRAG - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC,
            "pull_H0_vs_local_sigma": (tcc.H0_STEP25 - tcc.H0_LOCAL) / tcc.H0_LOCAL_SIGMA,
            "status_summary": "BAO fit: chi2_BAO(12x12)=20.79, p=0.054 (fixed M4 point); deviation on illustrative Planck prior: 100*theta*_reint pull +0.77 (1D prior); pull vs fixed target: (67.40-72.80)/1.00 = -5.40 against the assumed local H0 benchmark prior 72.8+/-1.0 (conditional pull, not a joint tension of two fitted measurements)",
        },
        {
            "scenario_id": "2_Profiled_M4_Shape_P_Free",
            "category": "Profiled_BAO_Control",
            "H0_km_s_Mpc": tcc.H0_LOCAL,
            "Omega_m_0": tcc.OMEGA_M_0,
            "Omega_K_0": 0.0,
            "omega_cb_h2": tcc.OMEGA_CB_H2,
            "rd_used_Mpc": prof["rd_req_m4_prof"],
            "rd_err_Mpc": prof["sigma_rd_m4_prof_fixedH0"],
            "P_H0_rd_km_s": prof["P_best_m4"],
            "P_H0_rd_err_km_s": prof["sigma_P_m4"],
            "chi2_DESI_BAO_Cobaya_full": prof["chi2_m4_prof_full"],
            "chi2_DESI_BAO_Table1": prof["chi2_m4_prof_tab1"],
            "theta_star_100_reint": theta100_m4,
            "pull_theta_star_reint_sigma": pull_theta_m4,
            "pull_omega_cb_planck_sigma": pull_wcb_m4,
            "pull_rd_vs_Planck1D_sigma": (prof["rd_req_m4_prof"] - tcc.RD_PLANCK_MPC) / math.sqrt(prof["sigma_rd_m4_prof_fixedH0"]**2 + tcc.RD_PLANCK_ERR_MPC**2),
            "pull_H0_vs_local_sigma": 0.0,
            "status_summary": "BAO fit: 1-param profile (fixed M4 Om=0.3148), P_best=10029.65+/-47.43 km/s, chi2_BAO=14.70; arithmetic consequence at the fixed target H0=72.80: rd=137.77+/-0.65 Mpc",
        },
        {
            "scenario_id": "3_Profiled_2Param_OmegaM_and_P_Free",
            "category": "Profiled_BAO_Control",
            "H0_km_s_Mpc": tcc.H0_LOCAL,
            "Omega_m_0": prof["om_mle_2p"],
            "Omega_K_0": 0.0,
            "omega_cb_h2": np.nan,
            "rd_used_Mpc": prof["rd_req_2p"],
            "rd_err_Mpc": prof["sigma_rd_2p_fixedH0"],
            "P_H0_rd_km_s": prof["P_best_2p"],
            "P_H0_rd_err_km_s": prof["sigma_P_2p"],
            "chi2_DESI_BAO_Cobaya_full": prof["chi2_2p_full"],
            "chi2_DESI_BAO_Table1": prof["chi2_2p_tab1"],
            "theta_star_100_reint": np.nan,
            "pull_theta_star_reint_sigma": np.nan,
            "pull_omega_cb_planck_sigma": np.nan,
            "pull_rd_vs_Planck1D_sigma": (prof["rd_req_2p"] - tcc.RD_PLANCK_MPC) / math.sqrt(prof["sigma_rd_2p_fixedH0"]**2 + tcc.RD_PLANCK_ERR_MPC**2),
            "pull_H0_vs_local_sigma": 0.0,
            "status_summary": "BAO fit: 2-param profile (Om=0.2939+/-0.0145, P=10193.67+/-126.47 km/s, chi2_BAO=12.74); at the fixed target H0=72.80: rd=140.02+/-1.74 Mpc; with the illustrative Planck rd prior: H0=69.30+/-0.87",
        },
        {
            "scenario_id": "4_Joint_3Param_BAO_PlanckRd_LocalH0",
            "category": "Joint_Profile_Test",
            "H0_km_s_Mpc": prof["joint_best_params"][1],
            "Omega_m_0": prof["joint_best_params"][0],
            "Omega_K_0": 0.0,
            "omega_cb_h2": np.nan,
            "rd_used_Mpc": prof["joint_best_params"][2],
            "rd_err_Mpc": tcc.RD_PLANCK_ERR_MPC,  # input Gaussian prior sigma on r_d, not a fitted joint uncertainty
            "P_H0_rd_km_s": prof["joint_best_params"][1] * prof["joint_best_params"][2],
            "P_H0_rd_err_km_s": np.nan,
            "chi2_DESI_BAO_Cobaya_full": jc["chi2_BAO"],
            "chi2_DESI_BAO_Table1": np.nan,
            "theta_star_100_reint": np.nan,
            "pull_theta_star_reint_sigma": np.nan,
            "pull_omega_cb_planck_sigma": np.nan,
            "pull_rd_vs_Planck1D_sigma": (prof["joint_best_params"][2] - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC,
            "pull_H0_vs_local_sigma": (prof["joint_best_params"][1] - tcc.H0_LOCAL) / tcc.H0_LOCAL_SIGMA,
            "status_summary": f"Joint 3-param (Om, H0, rd) control: DESI BAO + illustrative Planck rd prior + assumed local H0 benchmark prior 72.8+/-1.0: chi2_BAO={jc['chi2_BAO']:.4f}, chi2_rd={jc['chi2_rd_prior']:.4f}, chi2_H0={jc['chi2_H0_prior']:.4f}, chi2_joint={jc['chi2_joint']:.4f}; Delta_chi2 = {prof['dchi2_joint']:.2f} (sqrt(Delta_chi2) = {prof['sigma_joint_equiv']:.2f}); rd error column = input prior sigma",
        },
        {
            "scenario_id": "5_Rigid_LateTime_H0_72p80_FixedRd",
            "category": "NoGo_LateTime_Rescaling",
            "H0_km_s_Mpc": tcc.H0_LOCAL,
            "Omega_m_0": tcc.OMEGA_M_0,
            "Omega_K_0": 0.0,
            "omega_cb_h2": tcc.OMEGA_CB_H2,
            "rd_used_Mpc": tcc.RD_ZDRAG,
            "rd_err_Mpc": 0.0,
            "P_H0_rd_km_s": prof["P_rigid"],
            "P_H0_rd_err_km_s": 0.0,
            "chi2_DESI_BAO_Cobaya_full": prof["chi2_rig_full"],
            "chi2_DESI_BAO_Table1": prof["chi2_rig_tab1"],
            "theta_star_100_reint": 100.0 * (tcc.RS_ZSTAR / (dm_zstar_m4 / (tcc.H0_LOCAL / tcc.H0_STEP25))),
            "pull_theta_star_reint_sigma": (100.0 * (tcc.RS_ZSTAR / (dm_zstar_m4 / (tcc.H0_LOCAL / tcc.H0_STEP25))) - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100,
            "pull_omega_cb_planck_sigma": pull_wcb_m4,
            "pull_rd_vs_Planck1D_sigma": (tcc.RD_ZDRAG - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC,
            "pull_H0_vs_local_sigma": 0.0,
            "status_summary": "BAO fit: rigid all-z rescaling with unchanged rd gives a very poor DESI DR1 BAO fit (chi2_BAO=194.34, p=4.79e-35); its theta_star column uses D_M/s with unchanged r_s, an inconsistent geometric illustration (re-audit R2), not a late-transition model",
        },
        {
            "scenario_id": "6_Rigid_AllZ_H0_72p80_ScaledRd_FixedP",
            "category": "Fixed_P_Illustration",
            "H0_km_s_Mpc": tcc.H0_LOCAL,
            "Omega_m_0": tcc.OMEGA_M_0,
            "Omega_K_0": 0.0,
            "omega_cb_h2": tcc.OMEGA_CB_H2,
            "rd_used_Mpc": prof["rd_req_m4_fixed"],
            "rd_err_Mpc": 0.0,
            "P_H0_rd_km_s": prof["P_m4"],
            "P_H0_rd_err_km_s": 0.0,
            "chi2_DESI_BAO_Cobaya_full": prof["chi2_m4_full"],
            "chi2_DESI_BAO_Table1": prof["chi2_m4_tab1"],
            "theta_star_100_reint": theta100_m4,
            "pull_theta_star_reint_sigma": pull_theta_m4,
            "pull_omega_cb_planck_sigma": pull_wcb_m4,
            "pull_rd_vs_Planck1D_sigma": (prof["rd_req_m4_fixed"] - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC,
            "pull_H0_vs_local_sigma": 0.0,
            "status_summary": "Fixed-P illustration (P=9914.02 km/s -> rd=136.18 Mpc at the fixed target H0=72.80): conditional arithmetic offset -10.91 Mpc from the illustrative Planck rd value if early-universe physics is unmodified (no sigma label)",
        },
        {
            "scenario_id": "7_Flat_H0_72p80_Om0p244828_SelfConsistent",
            "category": "NoGo_LowMatter_Flat",
            "H0_km_s_Mpc": tcc.H0_LOCAL,
            "Omega_m_0": 0.244828,
            "Omega_K_0": 0.0,
            "omega_cb_h2": bg_flat_2448["omega_cb_h2"],
            "rd_used_Mpc": rd_2448,
            "rd_err_Mpc": 0.0,
            "P_H0_rd_km_s": tcc.H0_LOCAL * rd_2448,
            "P_H0_rd_err_km_s": 0.0,
            "chi2_DESI_BAO_Cobaya_full": chi2_bao_2448_reint,
            "chi2_DESI_BAO_Table1": tcc.bao_chi2_table1(bg_flat_2448["Hz_fn"], rd_2448, 0.0, tcc.H0_LOCAL),
            "theta_star_100_reint": theta100_2448,
            "pull_theta_star_reint_sigma": pull_theta_2448,
            "pull_omega_cb_planck_sigma": pull_wcb_2448,
            "pull_rd_vs_Planck1D_sigma": (rd_2448 - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC,
            "pull_H0_vs_local_sigma": 0.0,
            "status_summary": "BAO fit: chi2_BAO(12x12)=70.30 with re-integrated rd=150.65 Mpc; deviation on illustrative Planck prior: 100*theta*_reint=1.041745 (1D pull +2.18), omega_cb pull -10.98 on the assumed 1D prior (formal distances, not a Planck likelihood exclusion)",
        },
        {
            "scenario_id": "8_Open_H0_72p80_FixedOmegaCB_OmK0p0098",
            "category": "NoGo_Open_Curvature",
            "H0_km_s_Mpc": tcc.H0_LOCAL,
            "Omega_m_0": bg_open_728["Omega_m_0"],
            "Omega_K_0": ok_open_728,
            "omega_cb_h2": bg_open_728["omega_cb_h2"],
            "rd_used_Mpc": rd_open,
            "rd_err_Mpc": 0.0,
            "P_H0_rd_km_s": tcc.H0_LOCAL * rd_open,
            "P_H0_rd_err_km_s": 0.0,
            "chi2_DESI_BAO_Cobaya_full": chi2_bao_open_reint,
            "chi2_DESI_BAO_Table1": tcc.bao_chi2_table1(bg_open_728["Hz_fn"], rd_open, ok_open_728, tcc.H0_LOCAL),
            "theta_star_100_reint": theta100_open,
            "pull_theta_star_reint_sigma": pull_theta_open,
            "pull_omega_cb_planck_sigma": (bg_open_728["omega_cb_h2"] - tcc.PLANCK_3D_MEAN[1]) / tcc.PLANCK_3D_SIGMA[1],
            "pull_rd_vs_Planck1D_sigma": (rd_open - tcc.RD_PLANCK_MPC) / tcc.RD_PLANCK_ERR_MPC,
            "pull_H0_vs_local_sigma": 0.0,
            "status_summary": "BAO fit: chi2_BAO(12x12)=61.21 with re-integrated rd; deviation on illustrative Planck prior: 100*theta*_reint 1D pull +0.75 and omega_cb fixed at the M4 value 0.14237 (1D pull -0.02); Om_K=+0.0098 is not tested against any curvature likelihood in this package",
        },
    ]
    df_s2 = pd.DataFrame(rows)
    # Re-audit R4: explicit chi2 components for every row (BAO-only rows have zero prior terms)
    df_s2["chi2_BAO_only_Cobaya_full"] = df_s2["chi2_DESI_BAO_Cobaya_full"]
    df_s2["chi2_rd_prior"] = 0.0
    df_s2["chi2_H0_prior"] = 0.0
    df_s2["chi2_joint_total"] = df_s2["chi2_DESI_BAO_Cobaya_full"]
    j = df_s2["scenario_id"] == "4_Joint_3Param_BAO_PlanckRd_LocalH0"
    df_s2.loc[j, "chi2_rd_prior"] = jc["chi2_rd_prior"]
    df_s2.loc[j, "chi2_H0_prior"] = jc["chi2_H0_prior"]
    df_s2.loc[j, "chi2_joint_total"] = jc["chi2_joint"]
    df_s2["rd_err_kind"] = ["none", "fisher_fixed_H0", "fisher_fixed_H0", "input_prior_sigma", "none", "none", "none", "none"]
    csv_s2 = os.path.join(OUT_DIR, "step2_nogo_inverse_ladder_rd_H0.csv")
    df_s2.to_csv(csv_s2, index=False)
    return df_s2, csv_s2

def generate_inverse_ladder_plots(prof, lang="hu"):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11.2,
        "legend.fontsize": 8.6,
        "figure.dpi": 200,
    })

    h0_grid = np.linspace(63.5, 76.5, 300)
    P_2p = prof["P_best_2p"]
    sig_P_2p = prof["sigma_P_2p"]
    P_m4 = prof["P_m4"]

    rd_curve_2p = P_2p / h0_grid
    rd_2p_1s_lo = (P_2p - sig_P_2p) / h0_grid
    rd_2p_1s_hi = (P_2p + sig_P_2p) / h0_grid
    rd_2p_2s_lo = (P_2p - 2.0 * sig_P_2p) / h0_grid
    rd_2p_2s_hi = (P_2p + 2.0 * sig_P_2p) / h0_grid
    rd_curve_m4 = P_m4 / h0_grid

    fig, ax = plt.subplots(figsize=(9.4, 5.8))

    if lang == "hu":
        lbl_2p = (f"DESI DR1 BAO 2-param. profilozott hiperbola (rang-1 degeneracio):\n"
                  f"P_best = H0*r_d = {P_2p:.1f} +/- {sig_P_2p:.1f} km/s (Om_MLE = {prof['om_mle_2p']:.3f})")
        lbl_m4_h = f"Fix Step 25 M4 skalahiperbola: P_M4 = {P_m4:.1f} km/s (Om = 0.31482)"
        lbl_pl = f"Planck 2018 + TAL korai univerzum: r_d = {tcc.RD_PLANCK_MPC:.2f} +/- {tcc.RD_PLANCK_ERR_MPC:.2f} Mpc"
        lbl_sh = f"Feltetelezett lokalis H0 benchmark-prior: H0 = {tcc.H0_LOCAL:.2f} +/- {tcc.H0_LOCAL_SIGMA:.2f} km/s/Mpc"
        pt_m4 = "Step 25 M4 / TAL (67.40, 147.09 Mpc)"
        pt_bao_pl = f"DESI BAO + Planck r_d:\nH0 = {prof['H0_bao_planck_rd']:.2f} +/- {prof['sigma_H0_bao_planck_rd']:.2f}"
        pt_req_2p = f"Profilozott BAO H0=72.80-nal:\nr_d^req = {prof['rd_req_2p']:.2f} +/- {prof['sigma_rd_2p_fixedH0']:.2f} Mpc"
        pt_req_m4 = f"Fix P_M4 H0=72.80-nal:\nr_d = {prof['rd_req_m4_fixed']:.2f} Mpc"
        pt_rig = f"Merev H0=72.80, r_d=147.09:\nchi2_BAO = {prof['chi2_rig_full']:.1f}"
        note_txt = ("Megjegyzes: A BAO onmagaban a P = H0*r_d 1D szorzatot korlatozza (rang-1 a (H0, r_d) sikon).\n"
                    f"Kozos 3-param. (Om, H0, r_d) profil (DESI BAO + Planck r_d + feltetelezett lokalis H0 prior): "
                    f"Delta_chi2 = {prof['dchi2_joint']:.2f}, sqrt(Delta_chi2) = {prof['sigma_joint_equiv']:.2f}.")
        ttl = "V27 (Step 26) 2.3: Inverz távolságlétra a (H0, r_d) síkon — Profilozott BAO vs. Fix Step 25 M4"
        fname = "step2_3_inverse_ladder_rd_H0_exclusion.png"
    else:
        lbl_2p = (f"DESI DR1 BAO 2-param. profiled hyperbola (rank-1 degeneracy):\n"
                  f"P_best = H0*r_d = {P_2p:.1f} +/- {sig_P_2p:.1f} km/s (Om_MLE = {prof['om_mle_2p']:.3f})")
        lbl_m4_h = f"Fixed Step 25 M4 scale hyperbola: P_M4 = {P_m4:.1f} km/s (Om = 0.31482)"
        lbl_pl = f"Planck 2018 + TAL early universe: r_d = {tcc.RD_PLANCK_MPC:.2f} +/- {tcc.RD_PLANCK_ERR_MPC:.2f} Mpc"
        lbl_sh = f"Assumed local H0 benchmark prior: H0 = {tcc.H0_LOCAL:.2f} +/- {tcc.H0_LOCAL_SIGMA:.2f} km/s/Mpc"
        pt_m4 = "Step 25 M4 / TAL (67.40, 147.09 Mpc)"
        pt_bao_pl = f"DESI BAO + Planck r_d:\nH0 = {prof['H0_bao_planck_rd']:.2f} +/- {prof['sigma_H0_bao_planck_rd']:.2f}"
        pt_req_2p = f"Profiled BAO at H0=72.80:\nr_d^req = {prof['rd_req_2p']:.2f} +/- {prof['sigma_rd_2p_fixedH0']:.2f} Mpc"
        pt_req_m4 = f"Fixed P_M4 at H0=72.80:\nr_d = {prof['rd_req_m4_fixed']:.2f} Mpc"
        pt_rig = f"Rigid H0=72.80, r_d=147.09:\nchi2_BAO = {prof['chi2_rig_full']:.1f}"
        note_txt = ("Note: BAO alone constrains the 1D product P = H0*r_d (rank-1 in the (H0, r_d) plane).\n"
                    f"Joint 3-param. (Om, H0, r_d) profile (DESI BAO + Planck r_d + assumed local H0 prior): "
                    f"Delta_chi2 = {prof['dchi2_joint']:.2f}, sqrt(Delta_chi2) = {prof['sigma_joint_equiv']:.2f}.")
        ttl = "V27 (Step 26) 2.3: Inverse Distance Ladder in the (H0, r_d) Plane — Profiled BAO vs. Fixed Step 25 M4"
        fname = "step2_3_inverse_ladder_rd_H0_exclusion_en.png"

    # Plot BAO 1-sigma and 2-sigma bands around P_best_2p
    ax.fill_between(h0_grid, rd_2p_2s_lo, rd_2p_2s_hi, color="#2ca02c", alpha=0.14)
    ax.fill_between(h0_grid, rd_2p_1s_lo, rd_2p_1s_hi, color="#2ca02c", alpha=0.28, label=lbl_2p)
    ax.plot(h0_grid, rd_curve_2p, color="#1b7837", lw=2.3)
    ax.plot(h0_grid, rd_curve_m4, color="#1f4e78", lw=2.0, ls="--", label=lbl_m4_h)

    # Planck r_d horizontal band
    ax.axhspan(tcc.RD_PLANCK_MPC - 2.0 * tcc.RD_PLANCK_ERR_MPC,
               tcc.RD_PLANCK_MPC + 2.0 * tcc.RD_PLANCK_ERR_MPC, color="#1f4e78", alpha=0.15)
    ax.axhline(tcc.RD_PLANCK_MPC, color="#1f4e78", lw=1.6, ls=":", label=lbl_pl)

    # Local H0 vertical band
    ax.axvspan(tcc.H0_LOCAL - 2.0 * tcc.H0_LOCAL_SIGMA,
               tcc.H0_LOCAL + 2.0 * tcc.H0_LOCAL_SIGMA, color="#d95f02", alpha=0.12)
    ax.axvspan(tcc.H0_LOCAL - tcc.H0_LOCAL_SIGMA,
               tcc.H0_LOCAL + tcc.H0_LOCAL_SIGMA, color="#d95f02", alpha=0.24, label=lbl_sh)
    ax.axvline(tcc.H0_LOCAL, color="#d95f02", lw=1.6, ls="-.")

    # Key points
    ax.plot(tcc.H0_STEP25, tcc.RD_ZDRAG, "o", color="#1f4e78", ms=8, zorder=6)
    ax.annotate(pt_m4, xy=(tcc.H0_STEP25, tcc.RD_ZDRAG), xytext=(64.0, 149.2),
                arrowprops=dict(arrowstyle="->", color="#1f4e78", lw=1.1), fontsize=8.5, fontweight="bold", color="#1f4e78")

    ax.plot(prof["H0_bao_planck_rd"], tcc.RD_PLANCK_MPC, "D", color="#1b7837", ms=7.5, zorder=6)
    ax.annotate(pt_bao_pl, xy=(prof["H0_bao_planck_rd"], tcc.RD_PLANCK_MPC), xytext=(68.2, 151.3),
                arrowprops=dict(arrowstyle="->", color="#1b7837", lw=1.1), fontsize=8.4, fontweight="bold", color="#1b7837")

    ax.errorbar([tcc.H0_LOCAL], [prof["rd_req_2p"]], yerr=[prof["sigma_rd_2p_fixedH0"]],
                fmt="s", color="#1b7837", capsize=4, ms=7, zorder=6)
    ax.annotate(pt_req_2p, xy=(tcc.H0_LOCAL, prof["rd_req_2p"]), xytext=(69.3, 141.6),
                arrowprops=dict(arrowstyle="->", color="#1b7837", lw=1.1), fontsize=8.4, color="#1b7837")

    ax.plot(tcc.H0_LOCAL, prof["rd_req_m4_fixed"], "^", color="#1f4e78", ms=7.5, zorder=6)
    ax.annotate(pt_req_m4, xy=(tcc.H0_LOCAL, prof["rd_req_m4_fixed"]), xytext=(68.8, 134.3),
                arrowprops=dict(arrowstyle="->", color="#1f4e78", lw=1.1), fontsize=8.4, color="#1f4e78")

    ax.plot(tcc.H0_LOCAL, tcc.RD_ZDRAG, "X", color="#c00000", ms=9, zorder=6)
    ax.annotate(pt_rig, xy=(tcc.H0_LOCAL, tcc.RD_ZDRAG), xytext=(73.2, 149.0),
                arrowprops=dict(arrowstyle="->", color="#c00000", lw=1.1), fontsize=8.4, fontweight="bold", color="#c00000")

    ax.text(0.02, 0.03, note_txt, transform=ax.transAxes, fontsize=8.3,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8f9fa", edgecolor="#adb5bd", alpha=0.94))

    ax.set_xlim(63.8, 76.2)
    ax.set_ylim(132.5, 154.5)
    ax.set_xlabel(r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]")
    ax.set_ylabel(r"$r_d(z_{\mathrm{drag}})$ [Mpc]")
    ax.set_title(ttl)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, fname)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path

def main():
    prof = tcc.compute_bao_profile_controls()
    m_branches = analyze_m1_to_m4_branches()
    df_s2, csv_s2 = build_step2_table(prof, m_branches)
    p_hu = generate_inverse_ladder_plots(prof, lang="hu")
    p_en = generate_inverse_ladder_plots(prof, lang="en")

    print("=== STEP 2.1 M1-M4 BRANCH DIAGNOSTICS ===")
    print(f"  M2 uncompensated flat: E^2(0) = {m_branches['E2_m2_z0']:.6f} < 0 (unphysical under H0 > 0 z=0 normalization)")
    print(f"  M2 E^2(z) = 0 root: z_root (full rad+nu+matter) = {m_branches['z_root_full']:.6f}; (pure-matter approx) = {m_branches['z_root_matter_approx']:.6f}")
    print(f"  M3 uncompensated open: Omega_K = +{m_branches['om_k_m3']:.6f}, own r_s = {m_branches['rs_zstar_m3_own_Mpc']:.8f} Mpc, own r_d = {m_branches['rd_zdrag_m3_own_Mpc']:.8f} Mpc, 100*theta* = {m_branches['theta100_m3']:.8f} (formal 1D pull {m_branches['formal_pull_theta_m3']:.1f}), BAO chi2_Cobaya = {m_branches['chi2_bao_m3_full']:.5f}")

    print("\n=== STEP 2.2 & 2.3 INVERSE DISTANCE LADDER SUMMARY ===")
    for _, r in df_s2.iterrows():
        print(f"  {r['scenario_id']:42s}: H0={r['H0_km_s_Mpc']:.2f}, Om={r['Omega_m_0']:.4f}, rd={r['rd_used_Mpc']:.2f} Mpc, "
              f"P={r['P_H0_rd_km_s']:.1f} km/s, chi2_BAO={r['chi2_DESI_BAO_Cobaya_full']:.2f}")
    print(f"Saved: {csv_s2}, {p_hu}, {p_en}")

if __name__ == "__main__":
    main()
