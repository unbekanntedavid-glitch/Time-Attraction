#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1.1: DESI DR1 (2024, arXiv:2404.03002 & arXiv:2404.03000) BAO Audit
for the V26 (Step 25) Metric-Framed Yukawa Scalar-Tensor Theory (TAL).

Evaluates both the published DESI VI Table 1 summary (12 observables across 7
primary redshift bins + 3 single-tracer diagnostic rows) and the official
full-precision Cobaya 12-element data vector and 12x12 covariance matrix
(desi_2024_gaussian_bao_ALL_GCcomb_mean.txt / _cov.txt), including the 4
independent BAO profile controls (fixed M4, rigid H0=72.80, M4 shape with P
profiled, and 2-parameter flat LCDM with (Omega_m, P) profiled).
"""

import os
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tal_cosmology_common as tcc

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

PROVENANCE_NOTES = {
    "1_BGS":       "DESI VI (2404.03002) Table 1 primary isotropic BGS bin (1 observable in 12-pt fit)",
    "2_LRG1":      "DESI VI (2404.03002) Table 1 primary anisotropic LRG1 bin (2 observables in 12-pt fit)",
    "3_LRG2":      "DESI VI (2404.03002) Table 1 primary anisotropic LRG2 bin (2 observables in 12-pt fit)",
    "4_LRG3_ELG1": "DESI VI (2404.03002) Table 1 primary combined LRG3+ELG1 bin (2 observables in 12-pt fit)",
    "4a_LRG3_sep": "DESI III (2404.03000) Table 18 post-recon single-tracer LRG3 (DV/rd=19.82+/-0.18; excluded from 12-pt fit)",
    "4b_ELG1_955": "DESI III (2404.03000) Table 18 post-recon single-tracer ELG1 (DV/rd=20.01+/-0.41; excluded from 12-pt fit)",
    "4c_ELG1_845": "Chosen illustration (DV/rd=18.33+/-0.62) without a verified primary source; NOT a DESI measurement (not in DESI III Table 18 / VI Table 1); excluded from 12-pt fit (review P01)",
    "5_ELG2":      "DESI VI (2404.03002) Table 1 primary anisotropic ELG2 bin (2 observables in 12-pt fit)",
    "6_QSO":       "DESI VI (2404.03002) Table 1 primary isotropic QSO bin (1 observable in 12-pt fit)",
    "7_Lya_QSO":   "DESI VI (2404.03002) Table 1 primary anisotropic Lya QSO bin (2 observables in 12-pt fit)",
}

# Map primary bins to Cobaya 12-vector indices for exact per-bin full-precision chi2:
COBAYA_BIN_INDICES = {
    "1_BGS":       [0],
    "2_LRG1":      [1, 2],
    "3_LRG2":      [3, 4],
    "4_LRG3_ELG1": [5, 6],
    "5_ELG2":      [7, 8],
    "6_QSO":       [9],
    "7_Lya_QSO":   [10, 11],
}

def eval_bin_predictions(z, Hz_fn, rd_val):
    dm_rd = tcc.DM_of_z(z, Hz_fn, 0.0, tcc.H0_STEP25) / rd_val
    dh_rd = (tcc.C_KM_S / Hz_fn(z)) / rd_val
    dv_rd = (z * (dm_rd ** 2) * dh_rd) ** (1.0 / 3.0)
    return dm_rd, dh_rd, dv_rd

def eval_bin_chi2_table1(b, dm_rd, dh_rd, dv_rd):
    if b["mode"] == "iso":
        res_dv = (dv_rd - b["DV_rd"]) / b["DV_err"]
        return res_dv ** 2, np.nan, np.nan, res_dv
    res_dm = (dm_rd - b["DM_rd"]) / b["DM_err"]
    res_dh = (dh_rd - b["DH_rd"]) / b["DH_err"]
    r = b["r"]
    c2 = (res_dm ** 2 - 2.0 * r * res_dm * res_dh + res_dh ** 2) / (1.0 - r ** 2)
    res_dv = (dv_rd - b["DV_rd"]) / b["DV_err"] if "DV_rd" in b else np.nan
    return c2, res_dm, res_dh, res_dv

def eval_bin_chi2_cobaya_subblock(bin_id, dm_rd, dh_rd, dv_rd):
    if bin_id not in COBAYA_BIN_INDICES:
        return np.nan
    idx = COBAYA_BIN_INDICES[bin_id]
    sub_cov = tcc.COBAYA_COV[np.ix_(idx, idx)]
    sub_inv = np.linalg.inv(sub_cov)
    sub_mean = tcc.COBAYA_MEAN[idx]
    if len(idx) == 1:
        pred_v = np.array([dv_rd])
    else:
        pred_v = np.array([dm_rd, dh_rd])
    d = pred_v - sub_mean
    return float(d @ sub_inv @ d)

def generate_plots(df, prof, lang="hu"):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11.5,
        "legend.fontsize": 8.6,
        "figure.dpi": 200,
    })

    z_grid = np.linspace(0.05, 2.55, 300)
    Hz_m4 = tcc.BG_STEP25_M4["Hz_fn"]
    Hz_rig = lambda z: (tcc.H0_LOCAL / tcc.H0_STEP25) * Hz_m4(z)
    om_mle = prof["om_mle_2p"]
    P_2p = prof["P_best_2p"]
    Hz_prof2p = lambda z: (P_2p / tcc.RD_ZDRAG) * math.sqrt(om_mle * ((1.0 + z) ** 3) + (1.0 - om_mle))

    dm_m4, dh_m4, dv_m4 = [], [], []
    dm_rg, dh_rg, dv_rg = [], [], []
    dm_2p, dh_2p, dv_2p = [], [], []

    for zv in z_grid:
        a, b, c = eval_bin_predictions(zv, Hz_m4, tcc.RD_ZDRAG)
        dm_m4.append(a); dh_m4.append(b); dv_m4.append(c)
        a2, b2, c2 = eval_bin_predictions(zv, Hz_rig, tcc.RD_ZDRAG)
        dm_rg.append(a2); dh_rg.append(b2); dv_rg.append(c2)
        a3, b3, c3 = eval_bin_predictions(zv, Hz_prof2p, tcc.RD_ZDRAG)
        dm_2p.append(a3); dh_2p.append(b3); dv_2p.append(c3)

    dm_m4, dh_m4, dv_m4 = np.array(dm_m4), np.array(dh_m4), np.array(dv_m4)
    dm_rg, dh_rg, dv_rg = np.array(dm_rg), np.array(dh_rg), np.array(dv_rg)
    dm_2p, dh_2p, dv_2p = np.array(dm_2p), np.array(dh_2p), np.array(dv_2p)
    z_sqrt = np.sqrt(z_grid)

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.2))
    df_prim = df[df["primary_12pt"] == True]
    df_sep = df[df["primary_12pt"] == False]

    if lang == "hu":
        lbl_m4 = (f"Step 25 M4 / TAL (fix H0=67.40, rd=147.09 Mpc)\n"
                  f"chi2_Cobaya = {prof['chi2_m4_full']:.2f} (T1: {prof['chi2_m4_tab1']:.2f}), p = {prof['p_m4_full']:.3f}")
        lbl_2p = (f"Profilozott BAO kontroll (Om={om_mle:.3f}, P={P_2p:.1f} km/s)\n"
                  f"chi2_min = {prof['chi2_2p_full']:.2f} (dof=10, p = {prof['p_2p_full']:.3f})")
        lbl_rg = (f"Merev H0=72.80 + valtozatlan rd=147.09 Mpc\n"
                  f"chi2_Cobaya = {prof['chi2_rig_full']:.2f} (T1: {prof['chi2_rig_tab1']:.2f})")
        lbl_obs_p = "DESI DR1 (2024) fo bin (12 megfigyeles)"
        lbl_obs_s = "Kulonallo 1-traceres kontroll (LRG3 / ELG1);\nz=0.845: valasztott illusztracio, nem DESI-meres"
        t1 = "(a) Transzverzalis BAO: (D_M / r_d) / sqrt(z)"
        t2 = "(b) Radial BAO: z * (D_H / r_d)"
        t3 = "(c) Izotrop BAO: (D_V / r_d) / z^(2/3)"
        sup = ("V27 (Step 26) 1.1: DESI DR1 BAO (z = 0.1 - 4.2) osszevetes — "
               "Step 25 M4 vs. profilozott BAO kontroll vs. merev H0 = 72.80")
        fname = "step1_1_desi_dr1_bao_plot.png"
    else:
        lbl_m4 = (f"Step 25 M4 / TAL (fixed H0=67.40, rd=147.09 Mpc)\n"
                  f"chi2_Cobaya = {prof['chi2_m4_full']:.2f} (T1: {prof['chi2_m4_tab1']:.2f}), p = {prof['p_m4_full']:.3f}")
        lbl_2p = (f"Profiled BAO control (Om={om_mle:.3f}, P={P_2p:.1f} km/s)\n"
                  f"chi2_min = {prof['chi2_2p_full']:.2f} (dof=10, p = {prof['p_2p_full']:.3f})")
        lbl_rg = (f"Rigid H0=72.80 + unchanged rd=147.09 Mpc\n"
                  f"chi2_Cobaya = {prof['chi2_rig_full']:.2f} (T1: {prof['chi2_rig_tab1']:.2f})")
        lbl_obs_p = "DESI DR1 (2024) primary bins (12 observables)"
        lbl_obs_s = "Separate single-tracer controls (LRG3 / ELG1);\nz=0.845: chosen illustration, not a DESI measurement"
        t1 = "(a) Transverse BAO: (D_M / r_d) / sqrt(z)"
        t2 = "(b) Radial BAO: z * (D_H / r_d)"
        t3 = "(c) Isotropic BAO: (D_V / r_d) / z^(2/3)"
        sup = ("V27 (Step 26) 1.1: DESI DR1 BAO (z = 0.1 - 4.2) Comparison — "
               "Step 25 M4 vs. Profiled BAO Control vs. Rigid H0 = 72.80")
        fname = "step1_1_desi_dr1_bao_plot_en.png"

    # Panel 1: DM / rd / sqrt(z)
    ax = axes[0]
    ax.plot(z_grid, dm_m4 / z_sqrt, color="#1f4e78", lw=2.3, label=lbl_m4)
    ax.plot(z_grid, dm_2p / z_sqrt, color="#2ca02c", lw=1.9, ls="-.", label=lbl_2p)
    ax.plot(z_grid, dm_rg / z_sqrt, color="#c00000", lw=2.0, ls="--", label=lbl_rg)
    sub_dm = df_prim.dropna(subset=["DM_rd_obs"])
    ax.errorbar(sub_dm["z_eff"], sub_dm["DM_rd_obs"] / np.sqrt(sub_dm["z_eff"]),
                yerr=sub_dm["DM_rd_err"] / np.sqrt(sub_dm["z_eff"]),
                fmt="o", color="black", capsize=3.5, ms=5.5, label=lbl_obs_p)
    sub_dm_s = df_sep.dropna(subset=["DM_rd_obs"])
    ax.errorbar(sub_dm_s["z_eff"], sub_dm_s["DM_rd_obs"] / np.sqrt(sub_dm_s["z_eff"]),
                yerr=sub_dm_s["DM_rd_err"] / np.sqrt(sub_dm_s["z_eff"]),
                fmt="s", color="#666666", mfc="white", capsize=3.0, ms=5.0, label=lbl_obs_s)
    ax.set_xlabel("z_eff")
    ax.set_ylabel(r"$(D_M / r_d) / \sqrt{z}$")
    ax.set_title(t1)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")

    # Panel 2: z * DH / rd
    ax = axes[1]
    ax.plot(z_grid, z_grid * dh_m4, color="#1f4e78", lw=2.3, label="Step 25 M4 / TAL (67.40)")
    ax.plot(z_grid, z_grid * dh_2p, color="#2ca02c", lw=1.9, ls="-.", label=f"Profiled BAO (Om={om_mle:.3f})")
    ax.plot(z_grid, z_grid * dh_rg, color="#c00000", lw=2.0, ls="--", label="Rigid H0 = 72.80")
    sub_dh = df_prim.dropna(subset=["DH_rd_obs"])
    ax.errorbar(sub_dh["z_eff"], sub_dh["z_eff"] * sub_dh["DH_rd_obs"],
                yerr=sub_dh["z_eff"] * sub_dh["DH_rd_err"],
                fmt="o", color="black", capsize=3.5, ms=5.5, label="DESI DR1 D_H/r_d")
    sub_dh_s = df_sep.dropna(subset=["DH_rd_obs"])
    ax.errorbar(sub_dh_s["z_eff"], sub_dh_s["z_eff"] * sub_dh_s["DH_rd_obs"],
                yerr=sub_dh_s["z_eff"] * sub_dh_s["DH_rd_err"],
                fmt="s", color="#666666", mfc="white", capsize=3.0, ms=5.0, label="LRG3 (sep.)")
    ax.set_xlabel("z_eff")
    ax.set_ylabel(r"$z \cdot (D_H / r_d)$")
    ax.set_title(t2)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")

    # Panel 3: DV / rd / z^(2/3)
    ax = axes[2]
    ax.plot(z_grid, dv_m4 / (z_grid ** (2.0 / 3.0)), color="#1f4e78", lw=2.3, label="Step 25 M4 / TAL (67.40)")
    ax.plot(z_grid, dv_2p / (z_grid ** (2.0 / 3.0)), color="#2ca02c", lw=1.9, ls="-.", label=f"Profiled BAO (Om={om_mle:.3f})")
    ax.plot(z_grid, dv_rg / (z_grid ** (2.0 / 3.0)), color="#c00000", lw=2.0, ls="--", label="Rigid H0 = 72.80")
    sub_dv = df_prim.dropna(subset=["DV_rd_obs"])
    ax.errorbar(sub_dv["z_eff"], sub_dv["DV_rd_obs"] / (sub_dv["z_eff"] ** (2.0 / 3.0)),
                yerr=sub_dv["DV_rd_err"] / (sub_dv["z_eff"] ** (2.0 / 3.0)),
                fmt="o", color="black", capsize=3.5, ms=5.5, label="DESI DR1 D_V/r_d")
    sub_dv_s = df_sep.dropna(subset=["DV_rd_obs"])
    ax.errorbar(sub_dv_s["z_eff"], sub_dv_s["DV_rd_obs"] / (sub_dv_s["z_eff"] ** (2.0 / 3.0)),
                yerr=sub_dv_s["DV_rd_err"] / (sub_dv_s["z_eff"] ** (2.0 / 3.0)),
                fmt="s", color="#666666", mfc="white", capsize=3.0, ms=5.0, label="LRG3 / ELG1 (sep.)")
    ax.set_xlabel("z_eff")
    ax.set_ylabel(r"$(D_V / r_d) / z^{2/3}$")
    ax.set_title(t3)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    fig.suptitle(sup, fontsize=12.2, fontweight="bold", y=1.02)
    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, fname)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path

def main():
    Hz_m4 = tcc.BG_STEP25_M4["Hz_fn"]
    Hz_pl = tcc.BG_PLANCK_67P36["Hz_fn"]
    Hz_rig = lambda z: (tcc.H0_LOCAL / tcc.H0_STEP25) * Hz_m4(z)

    prof = tcc.compute_bao_profile_controls()
    om_mle = prof["om_mle_2p"]
    P_2p = prof["P_best_2p"]
    Hz_prof2p = lambda z: (P_2p / tcc.RD_ZDRAG) * math.sqrt(om_mle * ((1.0 + z) ** 3) + (1.0 - om_mle))

    rows = []
    for b in tcc.DESI_BINS:
        z = b["z_eff"]
        dm_m4, dh_m4, dv_m4 = eval_bin_predictions(z, Hz_m4, tcc.RD_ZDRAG)
        dm_pl, dh_pl, dv_pl = eval_bin_predictions(z, Hz_pl, tcc.RD_ZDRAG)
        dm_rg, dh_rg, dv_rg = eval_bin_predictions(z, Hz_rig, tcc.RD_ZDRAG)
        dm_2p, dh_2p, dv_2p = eval_bin_predictions(z, Hz_prof2p, tcc.RD_ZDRAG)

        c2_m4_t1, r_dm_m4, r_dh_m4, r_dv_m4 = eval_bin_chi2_table1(b, dm_m4, dh_m4, dv_m4)
        c2_pl_t1, _, _, _ = eval_bin_chi2_table1(b, dm_pl, dh_pl, dv_pl)
        c2_rg_t1, r_dm_rg, r_dh_rg, r_dv_rg = eval_bin_chi2_table1(b, dm_rg, dh_rg, dv_rg)
        c2_2p_t1, _, _, _ = eval_bin_chi2_table1(b, dm_2p, dh_2p, dv_2p)

        c2_m4_full = eval_bin_chi2_cobaya_subblock(b["bin_id"], dm_m4, dh_m4, dv_m4)
        c2_pl_full = eval_bin_chi2_cobaya_subblock(b["bin_id"], dm_pl, dh_pl, dv_pl)
        c2_rg_full = eval_bin_chi2_cobaya_subblock(b["bin_id"], dm_rg, dh_rg, dv_rg)
        c2_2p_full = eval_bin_chi2_cobaya_subblock(b["bin_id"], dm_2p, dh_2p, dv_2p)

        rho_m_z = tcc.OMEGA_M_0 * ((1.0 + z) ** 3) * tcc.RHO_CRIT_0_SI
        phi_attr_z = -tcc.ALPHA_Y * tcc.GAMMA_L * rho_m_z
        hz_s_inv = (Hz_m4(z) * 1000.0) / tcc.MPC_M
        delta_E2_rel = ((tcc.M_PHI_C_S_INV ** 2) * (phi_attr_z ** 2)) / (3.0 * tcc.BETA_YUK * (hz_s_inv ** 2))

        rows.append({
            "bin_id": b["bin_id"],
            "tracer": b["tracer"],
            "z_range": b["z_range"],
            "z_eff": z,
            "primary_12pt": b["primary_12pt"],
            "mode": b["mode"],
            "DM_rd_obs": b.get("DM_rd", np.nan),
            "DM_rd_err": b.get("DM_err", np.nan),
            "DH_rd_obs": b.get("DH_rd", np.nan),
            "DH_rd_err": b.get("DH_err", np.nan),
            "r_DM_DH": b.get("r", np.nan),
            "DV_rd_obs": b.get("DV_rd", np.nan),
            "DV_rd_err": b.get("DV_err", np.nan),
            "DM_rd_TAL_M4": dm_m4,
            "DH_rd_TAL_M4": dh_m4,
            "DV_rd_TAL_M4": dv_m4,
            "DM_rd_Profiled2P": dm_2p,
            "DH_rd_Profiled2P": dh_2p,
            "DV_rd_Profiled2P": dv_2p,
            "DM_rd_Rigid72p8": dm_rg,
            "DH_rd_Rigid72p8": dh_rg,
            "DV_rd_Rigid72p8": dv_rg,
            "pull_DM_TAL_M4_sigma": r_dm_m4,
            "pull_DH_TAL_M4_sigma": r_dh_m4,
            "pull_DV_TAL_M4_sigma": r_dv_m4,
            "pull_DM_Rigid72p8_sigma": r_dm_rg,
            "pull_DH_Rigid72p8_sigma": r_dh_rg,
            "pull_DV_Rigid72p8_sigma": r_dv_rg,
            "chi2_bin_TAL_M4_table1": c2_m4_t1,
            "chi2_bin_TAL_M4_cobaya_full": c2_m4_full,
            "chi2_bin_Planck67p36_table1": c2_pl_t1,
            "chi2_bin_Planck67p36_cobaya_full": c2_pl_full,
            "chi2_bin_Profiled2P_table1": c2_2p_t1,
            "chi2_bin_Profiled2P_cobaya_full": c2_2p_full,
            "chi2_bin_Rigid72p8_table1": c2_rg_t1,
            "chi2_bin_Rigid72p8_cobaya_full": c2_rg_full,
            "abs_phi_attr_z": abs(phi_attr_z),
            "delta_E2_over_E2_TAL": delta_E2_rel,
            "provenance_note": PROVENANCE_NOTES[b["bin_id"]],
        })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUT_DIR, "step1_1_desi_dr1_bao_full_table.csv")
    df.to_csv(csv_path, index=False)

    # Save dedicated BAO Profile Control Table (Review Section 4.6 & 6)
    prof_rows = [
        {
            "control_id": "1_Fixed_Step25_M4",
            "description": "Fixed Step 25 M4 / TAL baseline (no free BAO parameters; fixed H0=67.40, rd=147.0922 Mpc)",
            "Omega_m_0": tcc.OMEGA_M_0,
            "Omega_m_0_err": 0.0,
            "P_H0_rd_km_s": prof["P_m4"],
            "P_H0_rd_err_km_s": 0.0,
            "chi2_cobaya_full": prof["chi2_m4_full"],
            "chi2_table1_rounded": prof["chi2_m4_tab1"],
            "n_obs": 12,
            "n_fitted_params": 0,
            "dof": 12,
            "p_value_cobaya_full": prof["p_m4_full"],
            "implied_rd_at_H0_72p80_Mpc": prof["rd_req_m4_fixed"],
            "implied_rd_err_fixedH0_Mpc": 0.0,
        },
        {
            "control_id": "2_Rigid_H0_72p80_Fixed_rd",
            "description": "Fixed Step 25 M4 shape with rigid H0=72.80 km/s/Mpc and unchanged rd=147.0922 Mpc",
            "Omega_m_0": tcc.OMEGA_M_0,
            "Omega_m_0_err": 0.0,
            "P_H0_rd_km_s": prof["P_rigid"],
            "P_H0_rd_err_km_s": 0.0,
            "chi2_cobaya_full": prof["chi2_rig_full"],
            "chi2_table1_rounded": prof["chi2_rig_tab1"],
            "n_obs": 12,
            "n_fitted_params": 0,
            "dof": 12,
            "p_value_cobaya_full": prof["p_rig_full"],
            "implied_rd_at_H0_72p80_Mpc": tcc.RD_ZDRAG,
            "implied_rd_err_fixedH0_Mpc": 0.0,
        },
        {
            "control_id": "3_Fixed_M4_Shape_Profiled_P",
            "description": "Fixed Step 25 M4 shape (Omega_m=0.314820) with scale product P = H0*rd profiled",
            "Omega_m_0": tcc.OMEGA_M_0,
            "Omega_m_0_err": 0.0,
            "P_H0_rd_km_s": prof["P_best_m4"],
            "P_H0_rd_err_km_s": prof["sigma_P_m4"],
            "chi2_cobaya_full": prof["chi2_m4_prof_full"],
            "chi2_table1_rounded": prof["chi2_m4_prof_tab1"],
            "n_obs": 12,
            "n_fitted_params": 1,
            "dof": 11,
            "p_value_cobaya_full": prof["p_m4_prof_full"],
            "implied_rd_at_H0_72p80_Mpc": prof["rd_req_m4_prof"],
            "implied_rd_err_fixedH0_Mpc": prof["sigma_rd_m4_prof_fixedH0"],
        },
        {
            "control_id": "4_Profiled_2Param_OmegaM_and_P",
            "description": "Flat LCDM 2-parameter BAO shape control with both Omega_m and P = H0*rd profiled",
            "Omega_m_0": prof["om_mle_2p"],
            "Omega_m_0_err": prof["sigma_om_2p"],
            "P_H0_rd_km_s": prof["P_best_2p"],
            "P_H0_rd_err_km_s": prof["sigma_P_2p"],
            "chi2_cobaya_full": prof["chi2_2p_full"],
            "chi2_table1_rounded": prof["chi2_2p_tab1"],
            "n_obs": 12,
            "n_fitted_params": 2,
            "dof": 10,
            "p_value_cobaya_full": prof["p_2p_full"],
            "implied_rd_at_H0_72p80_Mpc": prof["rd_req_2p"],
            "implied_rd_err_fixedH0_Mpc": prof["sigma_rd_2p_fixedH0"],
        },
    ]
    df_prof = pd.DataFrame(prof_rows)
    prof_csv_path = os.path.join(OUT_DIR, "step1_1_desi_bao_profile_control.csv")
    df_prof.to_csv(prof_csv_path, index=False)

    plot_hu = generate_plots(df, prof, lang="hu")
    plot_en = generate_plots(df, prof, lang="en")

    print("=== STEP 1.1 DESI DR1 BAO SUMMARY ===")
    print(f"Sound horizon r_d(z_drag={tcc.Z_DRAG_PLANCK}) = {tcc.RD_ZDRAG:.4f} Mpc (Planck 2018: {tcc.RD_PLANCK_MPC} +/- {tcc.RD_PLANCK_ERR_MPC} Mpc)")
    print(f"1. Fixed Step 25 M4: chi2_Cobaya = {prof['chi2_m4_full']:.4f} (Table 1: {prof['chi2_m4_tab1']:.4f}), dof=12, p={prof['p_m4_full']:.4f}")
    print(f"2. Rigid H0=72.80:   chi2_Cobaya = {prof['chi2_rig_full']:.4f} (Table 1: {prof['chi2_rig_tab1']:.4f}), dof=12, p={prof['p_rig_full']:.4e}")
    print(f"3. M4 shape, P prof: P_best = {prof['P_best_m4']:.2f} +/- {prof['sigma_P_m4']:.2f} km/s, chi2_Cobaya = {prof['chi2_m4_prof_full']:.4f}, dof=11, p={prof['p_m4_prof_full']:.4f}")
    print(f"4. 2-param (Om, P):  Om_MLE = {prof['om_mle_2p']:.6f} +/- {prof['sigma_om_2p']:.5f}, P_best = {prof['P_best_2p']:.2f} +/- {prof['sigma_P_2p']:.2f} km/s, chi2_Cobaya = {prof['chi2_2p_full']:.4f}, dof=10, p={prof['p_2p_full']:.4f}")
    print(f"Saved: {csv_path}, {prof_csv_path}, {plot_hu}, {plot_en}")

if __name__ == "__main__":
    main()
