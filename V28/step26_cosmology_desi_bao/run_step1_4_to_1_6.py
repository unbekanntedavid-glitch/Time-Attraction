#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1.4, Step 1.5 & Step 1.6:
- Step 1.4: Redshift-Space Distortions (RSD) f*sigma8(z) Linear Growth Rate,
  First-Order Linear Variational ODE for the Yukawa Perturbation delta(f*sigma8),
  and Diagonal RSD Benchmark across 12 Measurements.
- Step 1.5: Algebraic BBN Order-of-Magnitude Benchmark (z ~ 10^8 - 10^9) for
  4He (Y_p) and D/H with Unified beta_Yuk = 1/(4*pi) Scalar Energy Density.
- Step 1.6: Constant-Density Bulk Benchmarks (phi_bulk = -gamma_L * rho),
  Finite-Sphere Surface Factor (R >> L), Full-Precision expm1 Clock Shifts,
  and the V21 Top-Hat Interior External Benchmark (|phi| = 0.0412, interior value, review F10).
"""

import os
import math
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tal_cosmology_common as tcc

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

SIGMA8_0_PLANCK = 0.8111

# 12 representative diagonal RSD f*sigma8(z) measurements (6dFGS, SDSS MGS, BOSS DR12, WiggleZ, eBOSS DR16, FastSound)
RSD_DATA = [
    ("6dFGS (Beutler+ 2012)",        0.067, 0.423, 0.055),
    ("SDSS MGS (Howlett+ 2015)",     0.150, 0.530, 0.160),
    ("BOSS DR12 (Alam+ 2017 z1)",    0.380, 0.497, 0.045),
    ("WiggleZ (Blake+ 2012 z1)",     0.440, 0.413, 0.080),
    ("BOSS DR12 (Alam+ 2017 z2)",    0.510, 0.459, 0.038),
    ("WiggleZ (Blake+ 2012 z2)",     0.600, 0.390, 0.063),
    ("BOSS DR12 (Alam+ 2017 z3)",    0.610, 0.436, 0.034),
    ("eBOSS DR16 LRG (Bautista+21)", 0.700, 0.473, 0.041),
    ("WiggleZ (Blake+ 2012 z3)",     0.730, 0.437, 0.072),
    ("eBOSS DR16 ELG (Tamone+ 20)",  0.850, 0.315, 0.095),
    ("FastSound (Okumura+ 2016)",    1.360, 0.482, 0.116),
    ("eBOSS DR16 QSO (Neveux+ 20)",  1.480, 0.462, 0.045),
]

def solve_growth_with_variational_yukawa(bg_dict, k_h_mpc=0.1, use_cb_only=True):
    """
    Solve the linear growth ODE in x = ln(a) from a_ini = 1e-3 to a = 1.0 using
    DOP853 with dense_output=True, rtol=1e-12, atol=1e-14, and simultaneously
    integrate the first-order linear variational ODE for delta_D(a) = D_TAL(a) - D_GR(a).
    """
    h = bg_dict["h"]
    om_source_0 = bg_dict["Omega_cb_0"] if use_cb_only else bg_dict["Omega_m_0"]
    E2_fn = bg_dict["E2_fn"]

    # Comoving wavenumber k in Mpc^-1 and m^-1:
    k_com_mpc_inv = k_h_mpc * h
    k_com_m_inv = k_com_mpc_inv / tcc.MPC_M
    kL_z0 = k_com_m_inv * tcc.L_M  # 2.184285e-20 at h=0.674, k=0.1 h/Mpc

    def _E2_and_dlnH_dlnA(a):
        z = (1.0 / a) - 1.0
        e2 = E2_fn(z)
        eps = 1e-6 * (1.0 + z)
        de2_dz = (E2_fn(z + eps) - E2_fn(z - eps)) / (2.0 * eps)
        dlnH_dlna = -0.5 * (1.0 + z) * de2_dz / e2
        om_src_a = om_source_0 * (a ** -3) / e2
        return e2, dlnH_dlna, om_src_a

    def _delta_Geff_over_GN_canonical(a):
        k_phys_L = kL_z0 / a
        k2 = k_phys_L * k_phys_L
        return (tcc.ALPHA_Y ** 2) * tcc.BETA_YUK * (k2 / (1.0 + k2))

    # 1) Direct float64 integration of GR and TAL (where 1 + delta_Geff/GN == 1.0 in IEEE-754 float64)
    def _ode_gr(x, y):
        a = math.exp(x)
        _, dlnH, om_a = _E2_and_dlnH_dlnA(a)
        D, dD = y
        ddD = -(2.0 + dlnH) * dD + 1.5 * om_a * D
        return [dD, ddD]

    def _ode_tal_float64(x, y):
        a = math.exp(x)
        _, dlnH, om_a = _E2_and_dlnH_dlnA(a)
        geff_ratio = 1.0 + _delta_Geff_over_GN_canonical(a)  # canonical normalization (re-audit 5.1)
        D, dD = y
        ddD = -(2.0 + dlnH) * dD + 1.5 * om_a * geff_ratio * D
        return [dD, ddD]

    # 2) Coupled 4-component variational ODE for (D_GR, dD_GR, delta_D_can, d_delta_D_can)
    def _ode_variational(x, y):
        a = math.exp(x)
        _, dlnH, om_a = _E2_and_dlnH_dlnA(a)
        dgeff_can = _delta_Geff_over_GN_canonical(a)
        D, dD, delD, d_delD = y
        ddD = -(2.0 + dlnH) * dD + 1.5 * om_a * D
        dd_delD = -(2.0 + dlnH) * d_delD + 1.5 * om_a * delD + 1.5 * om_a * dgeff_can * D
        return [dD, ddD, d_delD, dd_delD]

    a_ini = 1.0e-3
    x_ini = math.log(a_ini)
    sol_gr = solve_ivp(_ode_gr, [x_ini, 0.0], [a_ini, a_ini], method="DOP853",
                       rtol=1e-12, atol=1e-14, dense_output=True)
    sol_tal64 = solve_ivp(_ode_tal_float64, [x_ini, 0.0], [a_ini, a_ini], method="DOP853",
                          rtol=1e-12, atol=1e-14, dense_output=True)
    sol_var = solve_ivp(_ode_variational, [x_ini, 0.0], [a_ini, a_ini, 0.0, 0.0], method="DOP853",
                        rtol=1e-12, atol=1e-14, dense_output=True)

    D0_gr = float(sol_gr.sol(0.0)[0])
    D0_tal64 = float(sol_tal64.sol(0.0)[0])
    _, _, delD0_can_raw, _ = sol_var.sol(0.0)

    def eval_at_z(z):
        a = 1.0 / (1.0 + z)
        x = math.log(a)
        D_raw, dD_raw = sol_gr.sol(x)
        D64_raw, dD64_raw = sol_tal64.sol(x)
        _, _, delD_raw, d_delD_raw = sol_var.sol(x)

        D_norm = float(D_raw / D0_gr)
        f_gr = float(dD_raw / D_raw)
        fs8_gr = SIGMA8_0_PLANCK * f_gr * D_norm

        D64_norm = float(D64_raw / D0_tal64)
        f_tal64 = float(dD64_raw / D64_raw)
        fs8_tal64 = SIGMA8_0_PLANCK * f_tal64 * D64_norm
        diff_float64 = abs(fs8_tal64 - fs8_gr)

        # Variational perturbation of normalized growth rate f*sigma8(z) = sigma8_0 * d(D/D0)/d(ln a):
        delta_fs8_can_z0norm = SIGMA8_0_PLANCK * float((d_delD_raw - dD_raw * (delD0_can_raw / D0_gr)) / D0_gr)
        delta_fs8_2x_z0norm = 2.0 * delta_fs8_can_z0norm

        # Unnormalized (fixed a_ini amplitude) variational shift:
        delta_fs8_can_prim = SIGMA8_0_PLANCK * float(d_delD_raw / D0_gr)
        delta_fs8_2x_prim = 2.0 * delta_fs8_can_prim

        dgeff_can = _delta_Geff_over_GN_canonical(a)
        return {
            "D_norm": D_norm,
            "f_growth": f_gr,
            "fsigma8_GR": fs8_gr,
            "fsigma8_TAL_float64": fs8_tal64,
            "abs_diff_float64": diff_float64,
            "delta_Geff_over_GN_can": dgeff_can,
            "delta_Geff_over_GN_2x": 2.0 * dgeff_can,
            "delta_fsigma8_var_can_z0norm": abs(delta_fs8_can_z0norm),
            "delta_fsigma8_var_2x_z0norm": abs(delta_fs8_2x_z0norm),
            "delta_fsigma8_var_can_prim": abs(delta_fs8_can_prim),
            "delta_fsigma8_var_2x_prim": abs(delta_fs8_2x_prim),
        }

    return eval_at_z, kL_z0

def generate_rsd_plots(df_14, chi2_m4, chi2_low, max_var_can, lang="hu"):
    # Review F07: the legend shows the canonical (alpha_Y^2 beta_Yuk) maximum; the 2x value is only a
    # separate sensitivity control (CSV column delta_fsigma8_var_2x_z0norm).
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11.5,
        "legend.fontsize": 9.0,
        "figure.dpi": 200,
    })

    eval_m4, _ = solve_growth_with_variational_yukawa(tcc.BG_STEP25_M4, 0.1, True)
    bg_low = tcc.make_model_background(tcc.H0_LOCAL, 0.244828, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True)
    eval_low, _ = solve_growth_with_variational_yukawa(bg_low, 0.1, True)

    z_grid = np.linspace(0.0, 1.65, 200)
    fs8_m4_grid = [eval_m4(zv)["fsigma8_GR"] for zv in z_grid]
    fs8_low_grid = [eval_low(zv)["fsigma8_GR"] for zv in z_grid]

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    if lang == "hu":
        lbl_m4 = (f"Step 25 M4 / TAL (Om=0.3148, sigma8,0=0.8111)\n"
                  f"chi2_RSD,diag = {chi2_m4:.2f}/12; max|delta f*sigma8|_var,can = {max_var_can:.3e}")
        lbl_low = f"Alacsony Om=0.244828 (H0=72.80) kontroll: chi2_RSD,diag = {chi2_low:.2f}/12"
        lbl_obs = "RSD f*sigma8(z) meresek (6dFGS, SDSS, BOSS DR12, WiggleZ, eBOSS)"
        ttl = "V27 (Step 26) 1.4: Lineáris szerkezetnövekedés f*sigma8(z) — Diagonal RSD Benchmark"
        fname = "step1_4_rsd_fsigma8_plot.png"
    else:
        lbl_m4 = (f"Step 25 M4 / TAL (Om=0.3148, sigma8,0=0.8111)\n"
                  f"chi2_RSD,diag = {chi2_m4:.2f}/12; max|delta f*sigma8|_var,can = {max_var_can:.3e}")
        lbl_low = f"Low Om=0.244828 (H0=72.80) control: chi2_RSD,diag = {chi2_low:.2f}/12"
        lbl_obs = "RSD f*sigma8(z) measurements (6dFGS, SDSS, BOSS DR12, WiggleZ, eBOSS)"
        ttl = "V27 (Step 26) 1.4: Linear Structure Growth f*sigma8(z) — Diagonal RSD Benchmark"
        fname = "step1_4_rsd_fsigma8_plot_en.png"

    ax.plot(z_grid, fs8_m4_grid, color="#1f4e78", lw=2.4, label=lbl_m4)
    ax.plot(z_grid, fs8_low_grid, color="#c00000", lw=2.0, ls="--", label=lbl_low)
    ax.errorbar(df_14["z"], df_14["fsigma8_obs"], yerr=df_14["fsigma8_err"],
                fmt="o", color="black", capsize=3.5, ms=5.5, label=lbl_obs)
    ax.set_xlabel("z")
    ax.set_ylabel(r"$f\sigma_8(z)$")
    ax.set_title(ttl)
    ax.set_ylim(0.22, 0.68)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left")
    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, fname)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path

def run_step1_4():
    eval_m4_cb, kL_z0 = solve_growth_with_variational_yukawa(tcc.BG_STEP25_M4, 0.1, use_cb_only=True)
    eval_m4_allm, _ = solve_growth_with_variational_yukawa(tcc.BG_STEP25_M4, 0.1, use_cb_only=False)
    bg_low = tcc.make_model_background(tcc.H0_LOCAL, 0.244828, tcc.OMEGA_B_H2, tcc.OMEGA_NU_H2, 0.0, exact_z0_norm=True)
    eval_low_cb, _ = solve_growth_with_variational_yukawa(bg_low, 0.1, use_cb_only=True)

    rows = []
    for survey, z, obs, err in RSD_DATA:
        r_m4 = eval_m4_cb(z)
        r_m4_all = eval_m4_allm(z)
        r_low = eval_low_cb(z)
        pull_m4 = (r_m4["fsigma8_GR"] - obs) / err
        pull_low = (r_low["fsigma8_GR"] - obs) / err
        rows.append({
            "survey": survey,
            "z": z,
            "fsigma8_obs": obs,
            "fsigma8_err": err,
            "D_norm_M4": r_m4["D_norm"],
            "f_growth_M4": r_m4["f_growth"],
            "fsigma8_Step25_M4_TAL": r_m4["fsigma8_GR"],
            "fsigma8_LCDM_GR": r_m4["fsigma8_GR"],
            "fsigma8_Step25_M4_all_m": r_m4_all["fsigma8_GR"],
            "abs_diff_TAL_minus_GR_float64": r_m4["abs_diff_float64"],
            "delta_fsigma8_var_can_z0norm": r_m4["delta_fsigma8_var_can_z0norm"],
            "delta_fsigma8_var_2x_z0norm": r_m4["delta_fsigma8_var_2x_z0norm"],
            "delta_Geff_over_GN_can_at_k0p1h": r_m4["delta_Geff_over_GN_can"],
            "delta_Geff_over_GN_2x_at_k0p1h": r_m4["delta_Geff_over_GN_2x"],
            "pull_M4_sigma": pull_m4,
            "chi2_M4_point": pull_m4 ** 2,
            "fsigma8_Flat72p8_Om0p245": r_low["fsigma8_GR"],
            "pull_Flat72p8_Om0p245_sigma": pull_low,
            "chi2_Flat72p8_Om0p245_point": pull_low ** 2,
        })

    df_14 = pd.DataFrame(rows)
    csv_14 = os.path.join(OUT_DIR, "step1_4_rsd_fsigma8_growth.csv")
    df_14.to_csv(csv_14, index=False)

    chi2_m4 = float(df_14["chi2_M4_point"].sum())
    chi2_m4_allm = float(np.sum(((df_14["fsigma8_Step25_M4_all_m"] - df_14["fsigma8_obs"]) / df_14["fsigma8_err"]) ** 2))
    chi2_low = float(df_14["chi2_Flat72p8_Om0p245_point"].sum())
    max_var_can = float(df_14["delta_fsigma8_var_can_z0norm"].max())
    max_var_2x = float(df_14["delta_fsigma8_var_2x_z0norm"].max())

    p_hu = generate_rsd_plots(df_14, chi2_m4, chi2_low, max_var_can, lang="hu")
    p_en = generate_rsd_plots(df_14, chi2_m4, chi2_low, max_var_can, lang="en")
    return df_14, csv_14, p_hu, p_en, kL_z0, chi2_m4, chi2_m4_allm, chi2_low, max_var_can, max_var_2x

def run_step1_5():
    """
    Step 1.5: Algebraic BBN Order-of-Magnitude Benchmark (z = 1e8, 4e8, 1e9).
    Uses the unified action scalar energy density formula with beta_Yuk = 1/(4*pi)
    and standard linearized BBN sensitivity around Omega_b*h^2 = 0.02237 (eta_10 = 6.1245).
    """
    eta_10 = 273.78 * tcc.OMEGA_B_H2  # 6.124459
    Yp_std = 0.24709
    DH_e5_std = 2.5147
    Yp_obs, Yp_err = 0.2450, 0.0030
    DH_obs, DH_err = 2.5270, 0.0300

    z_epochs = [
        ("T_0p02MeV_end_BBN",      1.0e8),
        ("T_0p09MeV_D_bottleneck", 4.0e8),
        ("T_0p23MeV_np_freezeout", 1.0e9),
    ]

    om_r_total_early = tcc.OMEGA_GAMMA_0 * (1.0 + 0.22710731766 * tcc.N_EFF)
    rows = []
    for ep_name, z in z_epochs:
        zp1 = 1.0 + z
        E2_bbn = om_r_total_early * (zp1 ** 4) + tcc.OMEGA_M_0 * (zp1 ** 3)
        H_bbn_s_inv = tcc.H0_S_INV * math.sqrt(E2_bbn)
        rho_m_z = tcc.OMEGA_M_0 * (zp1 ** 3) * tcc.RHO_CRIT_0_SI
        phi_attr_z = -tcc.ALPHA_Y * tcc.GAMMA_L * rho_m_z
        abs_phi = abs(phi_attr_z)

        delta_H2_over_H2 = ((tcc.M_PHI_C_S_INV ** 2) * (phi_attr_z ** 2)) / (3.0 * tcc.BETA_YUK * (H_bbn_s_inv ** 2))
        delta_H_over_H = 0.5 * delta_H2_over_H2
        # Re-audit R6: ASSUMED order-of-magnitude proxies, not derived here.
        # (i) delta_G_over_G_clock = 2|phi| is a bare magnitude proxy (review F01). With g_J = A^2 g_E,
        #     A = exp(alpha_Y phi), the Jordan-frame Ricci term carries A^-2/G_* (A^-2 multiplies the Planck
        #     mass squared), so the tensor coupling is 1/G_T = A^-2/G_*, i.e. G_T = G_* A^2 = G_* exp(2 alpha_Y phi);
        #     phi < 0 therefore DECREASES G_T. G_T, the operationally measured G_lab ~ G_*(1+q) (q = alpha_Y^2 beta_Yuk,
        #     r << L, unscreened small bodies) and this bare 2|phi| proxy are three different quantities; the sign
        #     and the reference state (G measured today, |phi_0| ~ 2e-46) are not propagated -- only |2 phi| is used.
        # (ii) 0.08 and 1.6 are assumed linear sensitivities dYp ~ 0.08 dH/H and d ln(D/H) ~ 1.6 dH/H;
        #     Kneller & Steigman 2004 (New J. Phys. 6, 117) give dYp ~ 0.16 (S-1) and d ln(D/H) ~ 1.6 (S-1),
        #     S = H'/H, so 0.08 is an order-of-magnitude choice; no reaction network is run.
        # Standard values: Yp_std = 0.24709 (PRIMAT central value, Pitrou et al. 2018, Phys. Rep. 754, 1);
        # DH_e5_std = 2.5147 is an ARITHMETIC PROXY (review P02): a linearized benchmark value without verified
        # single-paper provenance; it is not a published BBN network prediction.
        # Observations: Yp = 0.245 +/- 0.003 (PDG 2024 BBN review, Fields, Molaro & Sarkar);
        # D/H = (2.527 +/- 0.030)e-5 (Cooke, Pettini & Steidel 2018, ApJ 855, 102); the field name
        # DH_e5_obs_PDG is a legacy label (PDG 2024 recommends 2.547 +/- 0.025).
        delta_G_over_G_clock = 2.0 * abs_phi

        dYp_tal = 0.08 * (delta_H_over_H + 0.5 * delta_G_over_G_clock)
        dDH_e5_tal = DH_e5_std * 1.6 * (delta_H_over_H + 0.5 * delta_G_over_G_clock)

        rows.append({
            "benchmark_type": "Algebraic_BBN_Order_of_Magnitude_Benchmark",
            "epoch_id": ep_name,
            "z_epoch": z,
            "H_BBN_s_inv": H_bbn_s_inv,
            "H_over_m_phi_c": H_bbn_s_inv / tcc.M_PHI_C_S_INV,
            "rho_m_kg_m3": rho_m_z,
            "abs_phi_attr": abs_phi,
            "delta_H_over_H_scalar_energy": delta_H_over_H,
            "delta_G_over_G_clock": delta_G_over_G_clock,
            "eta_10": eta_10,
            "Yp_TAL_pred": Yp_std + dYp_tal,
            "Yp_obs_PDG": Yp_obs,
            "Yp_obs_err": Yp_err,
            "pull_Yp_sigma": (Yp_std + dYp_tal - Yp_obs) / Yp_err,
            "DH_e5_TAL_pred": DH_e5_std + dDH_e5_tal,
            "DH_e5_obs_PDG": DH_obs,
            "DH_e5_obs_err": DH_err,
            "pull_DH_sigma": (DH_e5_std + dDH_e5_tal - DH_obs) / DH_err,
            "abs_shift_Yp_TAL": abs(dYp_tal),
            "abs_shift_DH_e5_TAL": abs(dDH_e5_tal),
        })

    df_15 = pd.DataFrame(rows)
    csv_15 = os.path.join(OUT_DIR, "step1_5_bbn_consistency.csv")
    df_15.to_csv(csv_15, index=False)
    return df_15, csv_15

def run_step1_6():
    """
    Step 1.6: Constant-Density Bulk Benchmarks (phi_bulk = -gamma_L * rho_bulk),
    Finite-Sphere Surface Factor (phi_surf ~ 0.5 * phi_bulk for R >> L),
    Full-Precision expm1 Clock Shifts, and the V21 Top-Hat Interior External Benchmark (review F10).
    """
    envs = [
        ("1_Cosmic_Background_z0",      "Homogeneous FLRW background at z=0 (bulk)",        tcc.RHO_M_0_SI, False, None, "Constant_Density_Bulk_Benchmark"),
        ("2_Local_Cosmic_Void",         "Cosmic void interior (0.2 * rho_m,0 bulk)",        0.2 * tcc.RHO_M_0_SI, False, None, "Constant_Density_Bulk_Benchmark"),
        ("3_Milky_Way_Disk_Solar_Nbr",  "Galactic disk ISM + stellar average (bulk)",       1.0e-20, False, None, "Constant_Density_Bulk_Benchmark"),
        ("4_Solar_Mean_Density_Bulk",   "Mean solar density constant-rho bulk benchmark",   1.408e3, True,  None, "Constant_Density_Bulk_Benchmark"),
        ("5_Earth_Mean_Density_Bulk",   "Mean terrestrial density constant-rho bulk",       5.514e3, True,  None, "Constant_Density_Bulk_Benchmark"),
        ("6_White_Dwarf_Core_Bulk",     "White dwarf interior constant-rho bulk benchmark", 1.0e9,   True,  None, "Constant_Density_Bulk_Benchmark"),
        # Review F10: interior value of the V21 top-hat fixed-point profile (read-only V26 tree,
        # scalar-tensor-yukawa-fit/data/v21/profiles_L10km.csv, SHA-256 0e4dd1b2...f74f4f):
        #   SLy,  mode "gradient", r = 2.0046409680935633 km, x = 0.041239199439254325 (line 237)
        #   APR4, mode "mass",     r = 2.75417748649486 km,   x = 0.041233030079658145 (line 452)
        # 0.0412 is the rounded common value; it is NOT a surface value (surface column = n.a.) and NOT a
        # solution of the present Helmholtz-Yukawa action. The Delta H0 column is a conditional illustration.
        ("7_V21_TopHat_Interior_External_Benchmark",
         "V21 top-hat interior external benchmark (interior |phi| = 0.0412; SLy gradient r=2.0046 km x=0.041239, "
         "APR4 mass r=2.7542 km x=0.041233; V26 data/v21/profiles_L10km.csv); not a surface value",
         float("nan"), False, 0.0412, "V21_TopHat_Interior_External_Benchmark"),
    ]

    rows = []
    for env_id, desc, rho_val, has_surf_half, bvp_override, row_type in envs:
        if bvp_override is not None:
            abs_phi_bulk = bvp_override
            abs_phi_surf = float("nan")  # review F10: no surface value exists for this interior benchmark
        else:
            abs_phi_bulk = tcc.ALPHA_Y * tcc.GAMMA_L * rho_val
            # Review N01: the R >> L half-space boundary value exists only where a boundary geometry is
            # defined (has_surf_half=True). For the homogeneous FLRW, void and Milky Way disk rows there is
            # no boundary geometry here: the boundary cell is n.a. (null), the bulk value is NOT copied into it
            # and NOT halved.
            abs_phi_surf = 0.5 * abs_phi_bulk if has_surf_half else float("nan")

        f_clock_minus_1 = float(np.expm1(-abs_phi_bulk))
        inv_f_clock_minus_1 = float(np.expm1(abs_phi_bulk))
        f_clock = 1.0 + f_clock_minus_1
        inv_f_clock = 1.0 + inv_f_clock_minus_1
        delta_H0_app_bulk = tcc.H0_STEP25 * inv_f_clock_minus_1
        delta_H0_app_surf = tcc.H0_STEP25 * float(np.expm1(abs_phi_surf)) if not math.isnan(abs_phi_surf) else float("nan")
        H0_app_bulk = tcc.H0_STEP25 + delta_H0_app_bulk

        rows.append({
            "env_id": env_id,
            "benchmark_category": row_type,
            "description": desc,
            "rho_kg_m3": rho_val,
            "abs_phi_bulk": abs_phi_bulk,
            "abs_phi_surface_R_gg_L": abs_phi_surf,
            "boundary_status": ("half-space boundary" if (bvp_override is None and has_surf_half)
                                else ("n.a. (interior benchmark, no surface value)" if bvp_override is not None
                                      else "n.a. (no boundary geometry)")),
            "f_clock_exp_phi": f_clock,
            "f_clock_minus_1_expm1": f_clock_minus_1,
            "inv_f_clock_exp_abs_phi": inv_f_clock,
            "inv_f_clock_minus_1_expm1": inv_f_clock_minus_1,
            "H0_apparent_bulk_km_s_Mpc": H0_app_bulk,
            "Delta_H0_apparent_bulk_km_s_Mpc": delta_H0_app_bulk,
            "Delta_H0_apparent_surf_km_s_Mpc": delta_H0_app_surf,
            "fraction_of_5p4_Hubble_gap_bulk": delta_H0_app_bulk / (tcc.H0_LOCAL - tcc.H0_STEP25),
            "provenance": ("V26 read-only tree scalar-tensor-yukawa-fit/data/v21/profiles_L10km.csv "
                           "(SHA-256 0e4dd1b2600dfbd33bb3a800b39608ed5342d50868cb1d50a4a9a30ab3f74f4f): "
                           "SLy gradient r=2.0046409680935633 km x=0.041239199439254325; "
                           "APR4 mass r=2.75417748649486 km x=0.041233030079658145; Delta H0 = conditional illustration"
                           if bvp_override is not None else "computed here: phi_bulk = alpha_Y gamma_L rho"),
        })

    df_16 = pd.DataFrame(rows)
    csv_16 = os.path.join(OUT_DIR, "step1_6_local_yukawa_environments.csv")
    df_16.to_csv(csv_16, index=False)
    return df_16, csv_16

def main():
    df_14, csv_14, p_hu, p_en, kL_z0, c2_m4, c2_m4_all, c2_low, m_can, m_2x = run_step1_4()
    df_15, csv_15 = run_step1_5()
    df_16, csv_16 = run_step1_6()

    dgeff_z0_can = (tcc.ALPHA_Y ** 2) * tcc.BETA_YUK * ((kL_z0 ** 2) / (1.0 + (kL_z0 ** 2)))
    dgeff_z0_2x = 2.0 * dgeff_z0_can
    print("=== STEP 1.4 RSD f*sigma8(z) LINEAR GROWTH & VARIATIONAL ODE SUMMARY ===")
    print(f"  At k = 0.1 h/Mpc = 0.0674 Mpc^-1, L = 10 km: k*L(z=0) = {kL_z0:.6e}, (k*L)^2 = {kL_z0**2:.6e}")
    print(f"  delta_Geff/GN (z=0): canonical = {dgeff_z0_can:.6e}, 2x convention = {dgeff_z0_2x:.6e}")
    print(f"  max|fsigma8_TAL - fsigma8_GR|_float64 = {df_14['abs_diff_TAL_minus_GR_float64'].max():.6e} (IEEE-754 underflow of 1 + 10^-40)")
    print(f"  max|delta(fsigma8)|_variational (z=0 norm): canonical = {m_can:.6e}; 2x convention = {m_2x:.6e}")
    print(f"  Diagonal RSD benchmark chi2 (Step 25 M4, clustering cb): {c2_m4:.4f} / 12 pts (all-m source: {c2_m4_all:.4f})")
    print(f"  Diagonal RSD benchmark chi2 (Flat H0=72.80, Om=0.244828): {c2_low:.4f} / 12 pts")

    print("\n=== STEP 1.5 ALGEBRAIC BBN ORDER-OF-MAGNITUDE BENCHMARK SUMMARY ===")
    for _, r in df_15.iterrows():
        print(f"  {r['epoch_id']:26s} (z={r['z_epoch']:.1e}): |phi|={r['abs_phi_attr']:.3e}, dH/H={r['delta_H_over_H_scalar_energy']:.3e}, "
              f"dG/G_clock={r['delta_G_over_G_clock']:.3e}, |dYp|={r['abs_shift_Yp_TAL']:.3e}, |d(D/H)e5|={r['abs_shift_DH_e5_TAL']:.3e}")

    print("\n=== STEP 1.6 LOCAL YUKAWA BULK BENCHMARKS & BVP REFERENCE SUMMARY ===")
    for _, r in df_16.iterrows():
        print(f"  {r['env_id']:29s}: |phi_bulk|={r['abs_phi_bulk']:.3e}, |phi_surf|={r['abs_phi_surface_R_gg_L']:.3e}, "
              f"Delta_H0_bulk={r['Delta_H0_apparent_bulk_km_s_Mpc']:.3e} km/s/Mpc")
    print(f"Saved: {csv_14}, {csv_15}, {csv_16}, {p_hu}, {p_en}")

if __name__ == "__main__":
    main()
