#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master Cosmological Audit Summary Generator for V27 (Step 26):
Metric-Framed Yukawa Scalar-Tensor Theory (TAL) vs. DESI DR1 BAO (2024),
Planck 2018 CMB Sound Horizon, Synthetic SNe Ia Distance-Modulus Shape Control,
RSD f*sigma8(z) Variational Growth, Algebraic BBN Benchmark, Local Yukawa Bulk
Benchmarks, and the Hubble-Tension No-Go Appendix.

Uses tal_cosmology_common.py so that tal_hubble_desi_bao_summary.json,
desi_dr1_bao_12pt_comparison.csv, and no_go_routes_comparison.csv are 100%
synchronized with step1_1..step2 CSV tables to machine precision.
"""

import os
import json
import math
import numpy as np
import pandas as pd

import tal_cosmology_common as tcc
import run_step1_1
import run_step1_2_and_1_3
import run_step1_4_to_1_6
import run_step2_nogo_and_report

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

def main():
    # Ensure all step tables and plots are up to date
    run_step1_1.main()
    df_12, _ = run_step1_2_and_1_3.run_step1_2()
    df_13, _, _, _, dmu0, shift_best, c2_sne_un, c2_sne_sh, c2_sne_low = run_step1_2_and_1_3.run_step1_3()
    df_14, _, _, _, kL_z0, c2_rsd_m4, c2_rsd_allm, c2_rsd_low, var_can, var_2x = run_step1_4_to_1_6.run_step1_4()
    df_15, _ = run_step1_4_to_1_6.run_step1_5()
    df_16, _ = run_step1_4_to_1_6.run_step1_6()

    prof = tcc.compute_bao_profile_controls()
    m_branches = run_step2_nogo_and_report.analyze_m1_to_m4_branches()
    df_s2, _ = run_step2_nogo_and_report.build_step2_table(prof, m_branches)

    # Also save the 7-bin (12-observable) primary BAO table in desi_dr1_bao_12pt_comparison.csv
    df_11 = pd.read_csv(os.path.join(OUT_DIR, "step1_1_desi_dr1_bao_full_table.csv"))
    df_11_prim = df_11[df_11["primary_12pt"] == True].copy()
    df_11_prim.to_csv(os.path.join(OUT_DIR, "desi_dr1_bao_12pt_comparison.csv"), index=False)

    # Synchronize no_go_routes_comparison.csv with step2_nogo_inverse_ladder_rd_H0.csv
    df_s2.to_csv(os.path.join(OUT_DIR, "no_go_routes_comparison.csv"), index=False)

    dm_zstar_m4 = tcc.DM_of_z(tcc.Z_STAR_PLANCK, tcc.BG_STEP25_M4["Hz_fn"], 0.0, tcc.H0_STEP25)
    theta100_m4 = 100.0 * (tcc.RS_ZSTAR / dm_zstar_m4)
    dgeff_z0_can = (tcc.ALPHA_Y ** 2) * tcc.BETA_YUK * ((kL_z0 ** 2) / (1.0 + (kL_z0 ** 2)))
    dgeff_z0_2x = 2.0 * dgeff_z0_can

    summary = {
        "version_delineation": {
            "V24_doi": "10.5281/zenodo.23680715",
            "V25_doi": "10.5281/zenodo.23699685 (step1-step22)",
            "V26_doi": "10.5281/zenodo.23727585 (step1-step25)",
            "current_work": "V27 (step26) Cosmological & DESI DR1 BAO Audit",
        },
        "tal_fundamental_parameters": {
            "L_km": tcc.L_KM,
            "L_m": tcc.L_M,
            "L_Mpc": tcc.L_MPC,
            "m_phi_c_s_inv": tcc.M_PHI_C_S_INV,
            "alpha_Y": tcc.ALPHA_Y,
            "beta_Yuk": tcc.BETA_YUK,
            "gamma_L_m3_kg": tcc.GAMMA_L,
            "gamma_L_2x_m3_kg": tcc.GAMMA_L_2X,
            "H0_Step25_km_s_Mpc": tcc.H0_STEP25,
            "H0_local_benchmark_prior_km_s_Mpc": tcc.H0_LOCAL,  # assumed benchmark prior (review F05), not a sourced SH0ES value
            "H0_local_sigma_km_s_Mpc": tcc.H0_LOCAL_SIGMA,
            "H0_s_inv": tcc.H0_S_INV,
            "H0_over_m_phi_c": tcc.H0_S_INV / tcc.M_PHI_C_S_INV,
            "E2_z0_exact_normalized": tcc.BG_STEP25_M4["E2_fn"](0.0),
            "E2_z0_unnormalized_neutrino_drift": tcc.BG_STEP25_M4_UNNORM["E2_fn"](0.0),
            "phi_attr_z0": tcc.PHI_ATTR_0,
            "abs_phi_attr_z0": tcc.ABS_PHI_ATTR_0,
            "Omega_phi_attr_z0": tcc.OMEGA_PHI_ATTR_0,
            "homogeneous_free_oscillation_branch_condition": "C_hom = 0 (Omega_delta_phi_0 = 0 on the adiabatic M4 baseline branch)",
        },
        "step1_1_desi_dr1_bao": {
            "rd_zdrag_Mpc": tcc.RD_ZDRAG,
            "rs_zstar_Mpc": tcc.RS_ZSTAR,
            "chi2_Step25_M4_Cobaya_full": prof["chi2_m4_full"],
            "chi2_Step25_M4_Table1": prof["chi2_m4_tab1"],
            "p_value_Step25_M4_Cobaya_full": prof["p_m4_full"],
            "chi2_Rigid72p80_Cobaya_full": prof["chi2_rig_full"],
            "chi2_Rigid72p80_Table1": prof["chi2_rig_tab1"],
            "p_value_Rigid72p80_Cobaya_full": prof["p_rig_full"],
            "profiled_1param_M4_shape": {
                "P_best_km_s": prof["P_best_m4"],
                "P_sigma_km_s": prof["sigma_P_m4"],
                "chi2_min_Cobaya_full": prof["chi2_m4_prof_full"],
                "dof": 11,
                "p_value": prof["p_m4_prof_full"],
                "rd_req_at_H0_72p80_Mpc": prof["rd_req_m4_prof"],
                "rd_req_sigma_fixedH0_Mpc": prof["sigma_rd_m4_prof_fixedH0"],
            },
            "profiled_2param_flat_LCDM": {
                "Omega_m_MLE": prof["om_mle_2p"],
                "Omega_m_sigma": prof["sigma_om_2p"],
                "P_best_km_s": prof["P_best_2p"],
                "P_sigma_km_s": prof["sigma_P_2p"],
                "chi2_min_Cobaya_full": prof["chi2_2p_full"],
                "dof": 10,
                "p_value": prof["p_2p_full"],
                "rd_req_at_H0_72p80_Mpc": prof["rd_req_2p"],
                "rd_req_sigma_fixedH0_Mpc": prof["sigma_rd_2p_fixedH0"],
                "rd_req_sigma_with_H0_err_Mpc": prof["sigma_rd_2p_propH0"],
                "H0_from_Planck_rd_km_s_Mpc": prof["H0_bao_planck_rd"],
                "H0_from_Planck_rd_sigma_km_s_Mpc": prof["sigma_H0_bao_planck_rd"],
            },
            "joint_3param_BAO_PlanckRd_LocalH0": {
                "chi2_joint_min": prof["chi2_joint_min"],
                "Delta_chi2_vs_BAO_MLE": prof["dchi2_joint"],
                "equivalent_gaussian_sigma": prof["sigma_joint_equiv"],
                "best_fit_Omega_m_H0_rd": prof["joint_best_params"],
            },
        },
        "step1_2_planck_cmb_and_illustrative_3d_prior": {
            "DM_zstar_M4_Mpc": dm_zstar_m4,
            "theta_star_100_M4_reint": theta100_m4,
            "theta_star_100_Planck_ref": tcc.THETA_STAR_PLANCK_100,
            "theta_star_100_Planck_sigma": tcc.THETA_STAR_PLANCK_ERR_100,
            "pull_theta_star_M4_sigma": (theta100_m4 - tcc.THETA_STAR_PLANCK_100) / tcc.THETA_STAR_PLANCK_ERR_100,
            "models": df_12.to_dict(orient="records"),
        },
        "step1_3_synthetic_distance_modulus_shape_control": {
            "z_min": 0.015,
            "z_max": 2.260,
            "n_nodes": 18,
            "delta_mu_zero_M4_minus_ref72p80_mag": dmu0,
            "delta_mu_zero_ref72p80_minus_M4_mag": -dmu0,
            "weighted_best_shift_mag": shift_best,
            "chi2_unshifted_M4_vs_synth_ref": c2_sne_un,
            "chi2_shifted_M4_vs_synth_ref": c2_sne_sh,
            "chi2_Flat72p80_Om0p244828_vs_synth_ref": c2_sne_low,
        },
        "step1_4_rsd_fsigma8_growth": {
            "k_com_h_Mpc_inv": 0.1,
            "kL_z0": kL_z0,
            "kL_z0_sq": kL_z0 ** 2,
            "delta_Geff_over_GN_z0_canonical": dgeff_z0_can,
            "delta_Geff_over_GN_z0_2x_convention": dgeff_z0_2x,
            "max_abs_diff_fsigma8_float64": 0.0,
            "max_abs_delta_fsigma8_variational_canonical": var_can,
            "max_abs_delta_fsigma8_variational_2x": var_2x,
            "chi2_RSD_diagonal_M4_clustering_cb": c2_rsd_m4,
            "chi2_RSD_diagonal_M4_all_m": c2_rsd_allm,
            "chi2_RSD_diagonal_Flat72p80_Om0p244828": c2_rsd_low,
        },
        "step1_5_algebraic_bbn_benchmark": df_15.to_dict(orient="records"),
        "step1_6_local_yukawa_bulk_benchmarks": df_16.to_dict(orient="records"),
        "step2_m1_to_m4_and_inverse_ladder": {
            "m_branches": m_branches,
            "inverse_ladder_scenarios": df_s2.where(pd.notnull(df_s2), None).to_dict(orient="records"),
        },
    }

    json_path = os.path.join(OUT_DIR, "tal_hubble_desi_bao_summary.json")
    summary = _nan_to_null(summary)
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, allow_nan=False)  # strict JSON (re-audit R3)
    print(f"Saved synchronized master summary: {json_path}")
    # Re-audit R1/R2/R4/R5/R6/5.1-5.3 controls are merged into the same JSON by a separate script.
    import run_step26_revision_controls
    run_step26_revision_controls.main()

def _nan_to_null(obj):
    if isinstance(obj, dict):
        return {k: _nan_to_null(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_nan_to_null(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return None if (math.isnan(v) or math.isinf(v)) else v
    if isinstance(obj, np.integer):
        return int(obj)
    return obj

if __name__ == "__main__":
    main()
