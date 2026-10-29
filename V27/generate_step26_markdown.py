#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V27 (Step 26): generate the HU and EN Markdown reports (and README.md) from the
single result object tal_hubble_desi_bao_summary.json and the Step 26 CSV tables.
Computed results are generated from the JSON/CSV sources; fixed input constants
(data values, priors, literature benchmarks) are stated in the code and cited in the text.
This script does not run any computation of Block 1 (run_step1_*); it only formats.
"""
import os, re, json, math, shutil, hashlib
import csv

D = os.path.dirname(os.path.abspath(__file__))
J = json.load(open(os.path.join(D, "tal_hubble_desi_bao_summary.json")))

def rcsv(name):
    with open(os.path.join(D, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

BAO = rcsv("step1_1_desi_dr1_bao_full_table.csv")
PROF = rcsv("step1_1_desi_bao_profile_control.csv")
CMB = rcsv("step1_2_planck_cmb_rd_consistency.csv")
SNE = rcsv("step1_3_pantheon_union3_sne_residuals.csv")
RSD = rcsv("step1_4_rsd_fsigma8_growth.csv")
BBN = rcsv("step1_5_bbn_consistency.csv")
LOC = rcsv("step1_6_local_yukawa_environments.csv")
NOGO = rcsv("no_go_routes_comparison.csv")

def f(x, n=3):
    return f"{float(x):.{n}f}"

def sci(x, n=3, signed=False):
    x = float(x)
    if x == 0.0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(abs(m), n) >= 10.0:
        e += 1; m = x / 10 ** e
    s = f"{m:+.{n}f}" if signed else f"{m:.{n}f}"
    return f"{s}\\times 10^{{{e}}}"

def num(x, n=3):
    """plain or scientific depending on magnitude (math-mode string)."""
    x = float(x)
    if x != 0 and (abs(x) < 1e-3 or abs(x) >= 1e5):
        return sci(x, n)
    return f(x, n)

TP = J["tal_fundamental_parameters"]
B1 = J["step1_1_desi_dr1_bao"]
P1 = B1["profiled_1param_M4_shape"]
P2 = B1["profiled_2param_flat_LCDM"]
P3 = B1["joint_3param_BAO_PlanckRd_LocalH0"]
C2 = J["step1_2_planck_cmb_and_illustrative_3d_prior"]
S3 = J["step1_3_synthetic_distance_modulus_shape_control"]
R4 = J["step1_4_rsd_fsigma8_growth"]
B5 = J["step1_5_algebraic_bbn_benchmark"]
L6 = J["step1_6_local_yukawa_bulk_benchmarks"]
N2 = J["step2_m1_to_m4_and_inverse_ladder"]
MB = N2["m_branches"]

beta = TP["beta_Yuk"]
H0P = 67.36
zstar, zdrag = 1089.80, 1059.94
m4 = [m for m in C2["models"] if m["model_id"] == "Step25_M4_TAL_Baseline"][0]
lowom = [m for m in C2["models"] if m["model_id"] == "Flat_H0_72p80_Om0p244828"][0]
rigid = [m for m in C2["models"] if m["model_id"] == "Rigid_AllZ_Rescaled_72p80"][0]

V = {
    "L_km": f(TP["L_km"], 0),
    "beta": f(beta, 6),
    "gammaL": sci(TP["gamma_L_m3_kg"], 5),
    "mphic": f(TP["m_phi_c_s_inv"], 4),
    "H0overm": sci(TP["H0_over_m_phi_c"], 3),
    "phi0": sci(TP["phi_attr_z0"], 4),
    "phistar": sci(TP["phi_attr_z0"] * (1 + zstar) ** 3, 3),
    "Omphi0": sci(TP["Omega_phi_attr_z0"], 3),
    "E2unnorm": f(TP["E2_z0_unnormalized_neutrino_drift"], 11),
    "E2offset": sci(TP["E2_z0_unnormalized_neutrino_drift"] - 1.0, 4),
    "rd": f(B1["rd_zdrag_Mpc"], 4), "rs": f(B1["rs_zstar_Mpc"], 4),
    "chi2M4": f(B1["chi2_Step25_M4_Cobaya_full"], 4), "chi2M4t": f(B1["chi2_Step25_M4_Table1"], 4),
    "pM4": f(B1["p_value_Step25_M4_Cobaya_full"], 4),
    "chi2R": f(B1["chi2_Rigid72p80_Cobaya_full"], 3), "chi2Rt": f(B1["chi2_Rigid72p80_Table1"], 3),
    "pR": sci(B1["p_value_Rigid72p80_Cobaya_full"], 2),
    "P1": f(P1["P_best_km_s"], 2), "P1s": f(P1["P_sigma_km_s"], 2), "chi2P1": f(P1["chi2_min_Cobaya_full"], 4),
    "pP1": f(P1["p_value"], 4), "rdP1": f(P1["rd_req_at_H0_72p80_Mpc"], 2), "rdP1s": f(P1["rd_req_sigma_fixedH0_Mpc"], 2),
    "Om2": f(P2["Omega_m_MLE"], 4), "Om2s": f(P2["Omega_m_sigma"], 4), "P2": f(P2["P_best_km_s"], 2),
    "P2s": f(P2["P_sigma_km_s"], 2), "chi2P2": f(P2["chi2_min_Cobaya_full"], 4), "pP2": f(P2["p_value"], 4),
    "rdP2": f(P2["rd_req_at_H0_72p80_Mpc"], 2), "rdP2s": f(P2["rd_req_sigma_fixedH0_Mpc"], 2),
    "rdP2sH": f(P2["rd_req_sigma_with_H0_err_Mpc"], 2), "H0Pl": f(P2["H0_from_Planck_rd_km_s_Mpc"], 2),
    "H0Pls": f(P2["H0_from_Planck_rd_sigma_km_s_Mpc"], 2),
    "chi2J": f(P3["chi2_joint_min"], 3), "dchi2J": f(P3["Delta_chi2_vs_BAO_MLE"], 3),
    "sigJ": f(P3["sqrt_Delta_chi2"], 2), "OmJ": f(P3["best_fit_Omega_m_H0_rd"][0], 4),
    "H0J": f(P3["best_fit_Omega_m_H0_rd"][1], 3), "rdJ": f(P3["best_fit_Omega_m_H0_rd"][2], 3),
    "PM4": f(N2["inverse_ladder_scenarios"][0]["P_H0_rd_km_s"], 2),
    "rdfixP": f(N2["inverse_ladder_scenarios"][5]["rd_used_Mpc"], 2),
    "th100M4": f(C2["theta_star_100_M4_reint"], 6), "th100ref": f(C2["theta_star_100_Planck_ref"], 5),
    "th100sig": f(C2["theta_star_100_Planck_sigma"], 5), "pullthM4": f(C2["pull_theta_star_M4_sigma"], 2),
    "chi23DM4": f(m4["chi2_Planck_3D_illustrative_prior"], 3),
    "rsLow": f(lowom["rs_zstar_reint_Mpc"], 4), "rdLow": f(lowom["rd_zdrag_reint_Mpc"], 4),
    "thLow": f(lowom["theta_star_100_reint"], 6), "thLowFix": f(lowom["theta_star_100_fixed_early_rs"], 6),
    "chi2Low": f(lowom["chi2_DESI_BAO_Cobaya_with_rd_reint"], 2), "chi2LowFix": f(lowom["chi2_DESI_BAO_Cobaya_with_fixed_rd"], 2),
    "pullwcbLow": f(lowom["pull_omega_cb_planck_sigma"], 2), "wcbLow": f(lowom["omega_cb_h2"], 5),
    "rsRig": f(rigid["rs_zstar_reint_Mpc"], 4), "rdRig": f(rigid["rd_zdrag_reint_Mpc"], 4),
    "thRig": f(rigid["theta_star_100_reint"], 6), "thRigFix": f(rigid["theta_star_100_fixed_early_rs"], 6),
    "zmin": f(S3["z_min"], 3), "zmax": f(S3["z_max"], 2), "nn": str(S3["n_nodes"]),
    "dmu": f(S3["delta_mu_zero_M4_minus_ref72p80_mag"], 6), "dmuexact": f(5 * math.log10(72.8 / 67.4), 6),
    "bshift": f(S3["weighted_best_shift_mag"], 5), "chi2Su": f(S3["chi2_unshifted_M4_vs_synth_ref"], 2),
    "chi2Ss": f(S3["chi2_shifted_M4_vs_synth_ref"], 4), "chi2Slow": f(S3["chi2_Flat72p80_Om0p244828_vs_synth_ref"], 2),
    "kL": sci(R4["kL_z0"], 5), "kL2": sci(R4["kL_z0_sq"], 4),
    "dGcan": sci(R4["delta_Geff_over_GN_z0_canonical"], 4), "dG2x": sci(R4["delta_Geff_over_GN_z0_2x_convention"], 4),
    "dfsCan": sci(R4["max_abs_delta_fsigma8_variational_canonical"], 3), "dfs2x": sci(R4["max_abs_delta_fsigma8_variational_2x"], 3),
    "chi2RSD": f(R4["chi2_RSD_diagonal_M4_clustering_cb"], 4), "chi2RSDm": f(R4["chi2_RSD_diagonal_M4_all_m"], 4),
    "chi2RSDlow": f(R4["chi2_RSD_diagonal_Flat72p80_Om0p244828"], 4), "kcom": f(R4["k_com_h_Mpc_inv"], 1),
    "eta10": f(B5[0]["eta_10"], 4), "Yp": f(B5[0]["Yp_TAL_pred"], 5), "YpO": f(B5[0]["Yp_obs_PDG"], 3),
    "YpE": f(B5[0]["Yp_obs_err"], 3), "pYp": f(B5[0]["pull_Yp_sigma"], 2, ), "DH": f(B5[0]["DH_e5_TAL_pred"], 4),
    "DHO": f(B5[0]["DH_e5_obs_PDG"], 3), "DHE": f(B5[0]["DH_e5_obs_err"], 3), "pDH": f(B5[0]["pull_DH_sigma"], 2),
    "HbbnMax": sci(max(b["H_over_m_phi_c"] for b in B5), 3), "phibbnMax": sci(max(b["abs_phi_attr"] for b in B5), 3),
    "dHbbnMax": sci(max(b["delta_H_over_H_scalar_energy"] for b in B5), 3),
    "dYpMax": sci(max(b["abs_shift_Yp_TAL"] for b in B5), 3),
    "E2M2": f(MB["E2_m2_z0"], 5), "zr_full": f(MB["z_root_full"], 6), "zr_mat": f(MB["z_root_matter_approx"], 6),
    "OmKM3": f(MB["om_k_m3"], 4), "DMM3": f(MB["dm_zstar_m3"], 1), "thM3": f(MB["theta100_m3"], 5),
    "pullM3": f(MB["formal_pull_theta_m3"], 1), "chi2M3": f(MB["chi2_bao_m3_full"], 1),
    "H0beta": f(H0P * (1 + beta), 4), "H0exp": f(H0P * math.exp(beta), 4),
    "ns_phi": f(L6[-1]["abs_phi_bulk"], 4), "ns_dH0": f(L6[-1]["Delta_H0_apparent_bulk_km_s_Mpc"], 3),
    "ns_frac": f(L6[-1]["fraction_of_5p4_Hubble_gap_bulk"], 3),
    "sun_dH0": sci(L6[3]["Delta_H0_apparent_bulk_km_s_Mpc"], 3),
}

RV = J["step26_revision_controls"]
RV1, RV2, RV4, RV5, RV6 = RV["R1_action_normalization"], RV["R2_piecewise_late_transition"], RV["R4_joint_chi2_split"], RV["R5_sne_symmetric_offsets"], RV["R6_bbn_proxy_provenance"]
NU = RV["item_5_2_neutrino_fermi_dirac"]
NU1000 = [r for r in NU["rows"] if r["z"] == 1000.0][0]
V.update({
    "R1Kratio": f(RV1["K_old_over_K_new"], 2), "R1Kv26": f(RV1["K_V26_over_K_new_at_alpha_Y_1"], 6),
    "R1betaV26": f(RV1["beta_Yuk_V26_at_alpha_Y_1"], 6), "R1aold": f(RV1["alpha_DEF_squared_old_same_alpha"], 6),
    "R1anew": f(RV1["alpha_DEF_squared_new"], 6), "R1Omv26": sci(RV1["Omega_phi_attr_0_V26_std_convention"], 3),
    "R2s": f(RV2["s"], 6), "R2illus": f(RV2["theta_star_100_DMoverS_illustration"], 6), "R2DMM4": f(RV2["DM_zstar_M4_Mpc"], 4),
    "R2rs": f(RV2["r_s_zstar_Mpc"], 6),
    "R4bao": f(RV4["chi2_BAO"], 8), "R4rd": f(RV4["chi2_rd_prior"], 8), "R4h0": f(RV4["chi2_H0_prior"], 8), "R4joint": f(RV4["chi2_joint"], 8),
    "R5offM4": f(RV5["M4"]["best_offset_mag"], 8), "R5offLow": f(RV5["lowOm_0p2448"]["best_offset_mag"], 8),
    "R5c0M4": f(RV5["M4"]["chi2_without_offset"], 5), "R5c1M4": f(RV5["M4"]["chi2_with_offset"], 6),
    "R5c0Low": f(RV5["lowOm_0p2448"]["chi2_without_offset"], 5), "R5c1Low": f(RV5["lowOm_0p2448"]["chi2_with_offset"], 6),
    "R5node": f(RV5["lowest_node_delta_mu_M4_minus_ref"], 8), "R5ident": f(RV5["identity_5log10_H0ratio"], 6),
    "NUdrift": sci(NU["z0_normalization_drift_Omega"], 4), "NUrel1000": sci(NU1000["rel_err_massive_nu"], 3, True),
    "NUH1000": sci(NU1000["rel_err_H_total"], 3, True), "NUmax": sci(NU["max_abs_rel_err_massive_nu"], 3),
    "NUmaxZ": f(NU["z_at_max_massive_nu"], 1), "NUHmax": sci(NU["max_abs_rel_err_H_total"], 3), "NUHmaxZ": f(NU["z_at_max_H_total"], 0),
    "MTom": f(RV["item_5_3_relabel"]["Omega_m"], 5), "MTwcb": f(RV["item_5_3_relabel"]["omega_cb_h2"], 5),
    "MTid": RV["item_5_3_relabel"]["new_id"].replace("_", "\\_"),
})
V["chi2J"] = f(RV4["chi2_joint"], 3)

def table_r2(lang):
    hdr = {"EN": [r"$z_t$", r"$D_M(z_*)$ [Mpc]", r"$100\theta_*$", r"pull vs $100\theta_*^{\rm Pl}$ (formal)"],
           "HU": [r"$z_t$", r"$D_M(z_*)$ [Mpc]", r"$100\theta_*$", r"pull a $100\theta_*^{\rm Pl}$-hoz (formális)"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|-----:|-----------:|---------:|--------:|"]
    for r in RV2["rows"]:
        out.append(f"| {r['z_t']:g} | {f(r['DM_zstar_Mpc'],4)} | {f(r['theta_star_100'],6)} | {f(r['pull_theta_vs_Planck_sigma'],1)} |")
    lab = {"EN": r"$D_M/s$ (no $z_t$; inconsistent illustration)", "HU": r"$D_M/s$ (nincs $z_t$; inkonzisztens szemléltetés)"}[lang]
    out.append(f"| {lab} | {f(RV2['DM_zstar_M4_Mpc']/RV2['s'],4)} | {f(RV2['theta_star_100_DMoverS_illustration'],6)} | — |")
    return "\n".join(out)

def table_r5(lang):
    hdr = {"EN": ["Model", "best offset [mag]", r"$\chi^2$ without offset", r"$\chi^2$ with offset (17 dof)"],
           "HU": ["Modell", u"legjobb eltolás [mag]", u"$\\chi^2$ eltolás nélkül", u"$\\chi^2$ eltolással (17 dof)"]}[lang]
    names = {"M4": ("M4 ($H_0=67.4$)", "M4 ($H_0=67.4$)"), "lowOm_0p2448": (r"Flat $H_0=72.8$, $\Omega_m=0.2448$", u"Sík $H_0=72.8$, $\\Omega_m=0.2448$")}
    i = 0 if lang == "EN" else 1
    out = ["| " + " | ".join(hdr) + " |", "|:--------------------|------:|------:|------:|"]
    for k in ("M4", "lowOm_0p2448"):
        r = RV5[k]
        out.append(f"| {names[k][i]} | {f(r['best_offset_mag'],8)} | {f(r['chi2_without_offset'],5)} | {f(r['chi2_with_offset'],6)} |")
    return "\n".join(out)

def table_nu(lang):
    hdr = {"EN": [r"$z$", r"$\rho_\nu^{\rm FD}$ [$E^2$ units]", r"$\rho_\nu^{p=1.83}$ [$E^2$ units]", "rel. error (massive)", r"rel. error ($H$)"],
           "HU": [r"$z$", u"$\\rho_\\nu^{\\rm FD}$ [$E^2$-egység]", u"$\\rho_\\nu^{p=1.83}$ [$E^2$-egység]", u"rel. hiba (tömeges)", u"rel. hiba ($H$)"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|-----:|-------:|-------:|------:|------:|"]
    for r in NU["rows"]:
        out.append(f"| {r['z']:g} | ${sci(r['rho_nu_exact_E2'],5)}$ | ${sci(r['rho_nu_interp_E2'],5)}$ | ${sci(r['rel_err_massive_nu'],3,True)}$ | ${sci(r['rel_err_H_total'],3,True)}$ |")
    return "\n".join(out)

TR = {"Lya QSO": r"Ly$\alpha$ QSO"}

def obs_cell(r):
    parts = []
    if r["DM_rd_obs"]:
        parts.append(f"$D_M$: ${f(r['DM_rd_obs'],2)}\\pm{f(r['DM_rd_err'],2)}$")
        parts.append(f"$D_H$: ${f(r['DH_rd_obs'],2)}\\pm{f(r['DH_rd_err'],2)}$ ($r={f(r['r_DM_DH'],3)}$)")
    else:
        parts.append(f"$D_V$: ${f(r['DV_rd_obs'],2)}\\pm{f(r['DV_rd_err'],2)}$")
    return "; ".join(parts)

def pred_cell(r, key):
    if r["DM_rd_obs"]:
        return f"${f(r['DM_rd_'+key],3)}$; ${f(r['DH_rd_'+key],3)}$"
    return f"${f(r['DV_rd_'+key],3)}$"

def table_bao(lang):
    hdr = {"EN": ["Tracer", r"$z_{\rm eff}$", r"DESI DR1 ($X/r_d$)", "M4 prediction", r"$\chi^2_{\rm bin}$ M4", r"$\chi^2_{\rm bin}$ 2P-prof.", r"$\chi^2_{\rm bin}$ rigid 72.8"],
           "HU": ["Tracer", r"$z_{\rm eff}$", r"DESI DR1 ($X/r_d$)", u"M4 jóslat", r"$\chi^2_{\rm bin}$ M4", r"$\chi^2_{\rm bin}$ 2P-prof.", r"$\chi^2_{\rm bin}$ merev 72.8"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|:-----------|-----:|:------------------------------------------|:--------------------|-----:|-----:|-----:|"]
    tot = [0.0, 0.0, 0.0]
    for r in BAO:
        prim = r["primary_12pt"] == "True"
        name = TR.get(r["tracer"], r["tracer"])
        if prim:
            c = [float(r["chi2_bin_TAL_M4_cobaya_full"]), float(r["chi2_bin_Profiled2P_cobaya_full"]), float(r["chi2_bin_Rigid72p8_cobaya_full"])]
            tot = [a + b for a, b in zip(tot, c)]
            cells = [f"**{name}**", f(r["z_eff"], 3), obs_cell(r), pred_cell(r, "TAL_M4"), f(c[0], 3), f(c[1], 3), f(c[2], 2)]
        else:
            c = [float(r["chi2_bin_TAL_M4_table1"]), float(r["chi2_bin_Profiled2P_table1"]), float(r["chi2_bin_Rigid72p8_table1"])]
            tag = {"EN": " (diagn.)", "HU": " (diagn.)"}[lang]
            if r["bin_id"] == "4c_ELG1_845":  # review P01
                name = "ELG1"
                tag = {"EN": " (chosen illustration, not a DESI measurement)", "HU": u" (választott illusztráció, nem DESI-mérés)"}[lang]
            cells = [f"*{name}{tag}*", f(r["z_eff"], 3), obs_cell(r), pred_cell(r, "TAL_M4"), f"*{f(c[0],3)}*", f"*{f(c[1],3)}*", f"*{f(c[2],2)}*"]
        out.append("| " + " | ".join(cells) + " |")
    lbl = {"EN": "**Sum, 12 primary observables (12×12 cov.)**", "HU": "**Összeg, 12 elsődleges megfigyelés (12×12 kov.)**"}[lang]
    out.append(f"| {lbl} | | | | **{f(tot[0],3)}** | **{f(tot[1],3)}** | **{f(tot[2],2)}** |")
    return "\n".join(out)

def table_prof(lang):
    names = {"EN": {"1_Fixed_Step25_M4": "Fixed M4 point (no free BAO parameter)",
                    "2_Rigid_H0_72p80_Fixed_rd": r"M4 shape, rigid $H_0=72.8$, unchanged $r_d$",
                    "3_Fixed_M4_Shape_Profiled_P": r"M4 shape, $P=H_0r_d$ profiled",
                    "4_Profiled_2Param_OmegaM_and_P": r"Flat $\Lambda$CDM shape, $\Omega_m$ and $P$ profiled"},
             "HU": {"1_Fixed_Step25_M4": u"Rögzített M4-pont (nincs szabad BAO-paraméter)",
                    "2_Rigid_H0_72p80_Fixed_rd": u"M4-alak, merev $H_0=72.8$, változatlan $r_d$",
                    "3_Fixed_M4_Shape_Profiled_P": u"M4-alak, $P=H_0r_d$ profilozva",
                    "4_Profiled_2Param_OmegaM_and_P": u"Sík $\\Lambda$CDM-alak, $\\Omega_m$ és $P$ profilozva"}}[lang]
    hdr = {"EN": ["Control", r"$\Omega_m$", r"$P=H_0r_d$ [km/s]", r"$\chi^2$ (12×12)", "dof", "$p$", r"$r_d$ needed at $H_0=72.8$ [Mpc]"],
           "HU": [u"Kontroll", r"$\Omega_m$", r"$P=H_0r_d$ [km/s]", r"$\chi^2$ (12×12)", "dof", "$p$", u"Szükséges $r_d$ $H_0=72.8$ mellett [Mpc]"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|:------------------------------------|-----:|------------:|------:|---:|------:|------------:|"]
    for r in PROF:
        om = f(r["Omega_m_0"], 4) + (f" $\\pm$ {f(r['Omega_m_0_err'],4)}" if float(r["Omega_m_0_err"]) > 0 else "")
        P = f(r["P_H0_rd_km_s"], 2) + (f" $\\pm$ {f(r['P_H0_rd_err_km_s'],2)}" if float(r["P_H0_rd_err_km_s"]) > 0 else "")
        p = float(r["p_value_cobaya_full"]); ps = f"${sci(p,2)}$" if p < 1e-3 else f(p, 4)
        rd = f(r["implied_rd_at_H0_72p80_Mpc"], 2) + (f" $\\pm$ {f(r['implied_rd_err_fixedH0_Mpc'],2)}" if float(r["implied_rd_err_fixedH0_Mpc"]) > 0 else "")
        if r["control_id"].startswith("2_"):
            rd = u"— (147.09 " + ("fixed)" if lang == "EN" else u"rögzített)")
        out.append(f"| {names[r['control_id']]} | {om} | {P} | {f(r['chi2_cobaya_full'],4)} | {r['dof']} | {ps} | {rd} |")
    return "\n".join(out)

CMB_NAMES = {
 "Step25_M4_TAL_Baseline": ("Step 25 M4 / TAL", "Step 25 M4 / TAL"),
 "Planck2018_67p36_Ref": ("Planck 2018 reference", "Planck 2018 referencia"),
 "Rigid_AllZ_Rescaled_72p80": (r"All-$z$ rigid $H\to sH$ (kinematic control)", r"Minden $z$-n merev $H\to sH$ (kinematikai kontroll)"),
 "Inconsistent_DMoverS_FixedEarly_Illustration": (r"*$D_M/s$ with fixed early $r_s,r_d$ (inconsistent geometric illustration, not a model)*", u"*$D_M/s$ rögzített korai $r_s,r_d$-vel (inkonzisztens geometriai szemléltetés, nem modell)*"),
 "Flat_H0_72p80_SameOm0p3148": (r"Flat $H_0=72.8$, same $\Omega_m$", u"Sík $H_0=72.8$, azonos $\\Omega_m$"),
 "Flat_H0_72p80_Om0p244828": (r"Flat $H_0=72.8$, $\Omega_m=0.244828$", u"Sík $H_0=72.8$, $\\Omega_m=0.244828$"),
 "Flat_H0_72p80_MatchThetaReint_Om0p2442": (r"Flat $H_0=72.8$, $\theta_*$-matched $\Omega_m\approx0.2442$", u"Sík $H_0=72.8$, $\\theta_*$-hez illesztett $\\Omega_m\\approx0.2442$"),
 "Open_H0_72p80_FixedOmegaCB_OmK0p0098": (r"Open $H_0=72.8$, fixed $\omega_{cb}$", u"Nyílt $H_0=72.8$, rögzített $\\omega_{cb}$"),
}

def table_cmb(lang):
    i = 0 if lang == "EN" else 1
    hdr = {"EN": ["Model", r"$H_0$", r"$\Omega_m$", r"$\Omega_K$", r"$\omega_{cb}$", r"$r_s(z_*)$ [Mpc]", r"$r_d(z_{\rm drag})$ [Mpc]", r"$100\theta_*$", r"$\chi^2_{\rm 3D,illustr.}$", r"$\chi^2_{\rm BAO}$ ($r_d$ re-int.)"],
           "HU": [u"Modell", r"$H_0$", r"$\Omega_m$", r"$\Omega_K$", r"$\omega_{cb}$", r"$r_s(z_*)$ [Mpc]", r"$r_d(z_{\rm drag})$ [Mpc]", r"$100\theta_*$", r"$\chi^2_{\rm 3D,szeml.}$", r"$\chi^2_{\rm BAO}$ (újraint. $r_d$)"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|:------------------------------|----:|------:|------:|------:|-------:|-------:|-------:|---------:|-------:|"]
    for r in CMB:
        x2 = float(r["chi2_Planck_3D_illustrative_prior"])
        out.append("| " + " | ".join([CMB_NAMES[r["model_id"]][i], f(r["H0_km_s_Mpc"], 2), f(r["Omega_m_0"], 4), f(r["Omega_K_0"], 4),
                    f(r["omega_cb_h2"], 5), f(r["rs_zstar_reint_Mpc"], 3), f(r["rd_zdrag_reint_Mpc"], 3), f(r["theta_star_100_reint"], 5),
                    (f(x2, 3) if x2 < 1000 else f(x2, 1)), f(r["chi2_DESI_BAO_Cobaya_with_rd_reint"], 2)]) + " |")
    return "\n".join(out)

def table_sne(lang):
    hdr = {"EN": [r"$z$ node", r"$\sigma_\mu$ [mag]", r"$\mu_{\rm ref}$ (72.8, 0.334)", r"$\mu_{\rm M4}$ (67.4)", r"$\Delta\mu_{\rm M4-ref}$", r"$\Delta\mu_{\rm M4,shift-ref}$", r"$\Delta\mu_{0.2448-{\rm ref}}$", r"$\Delta\mu_{0.2448,{\rm shift-ref}}$"],
           "HU": [r"$z$ csomópont", r"$\sigma_\mu$ [mag]", r"$\mu_{\rm ref}$ (72.8, 0.334)", r"$\mu_{\rm M4}$ (67.4)", r"$\Delta\mu_{\rm M4-ref}$", r"$\Delta\mu_{\rm M4,eltolt-ref}$", r"$\Delta\mu_{0.2448-{\rm ref}}$", r"$\Delta\mu_{0.2448,{\rm eltolt-ref}}$"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|------:|------:|----------:|----------:|---------:|---------:|---------:|---------:|"]
    for r in SNE:
        out.append("| " + " | ".join([f(r["z_bin"], 3), f(r["sigma_mu_rep"], 3), f(r["mu_synth_ref_72p8_Om0p334"], 4), f(r["mu_TAL_M4_67p4"], 4),
                    f"{float(r['delta_mu_M4_minus_ref']):+.4f}", f"{float(r['delta_mu_M4_shifted_minus_ref']):+.4f}", f"{float(r['delta_mu_Om0p2448_minus_ref']):+.4f}", f"{float(r['delta_mu_Om0p2448_shifted_minus_ref']):+.4f}"]) + " |")
    return "\n".join(out)

def table_rsd(lang):
    hdr = {"EN": ["Survey", r"$z$", r"$f\sigma_8$ obs.", r"$f\sigma_8$ M4", "pull M4", r"$f\sigma_8$ (72.8, 0.2448)", r"$\lvert\delta f\sigma_8\rvert$ variational", r"$\Delta G_{\rm eff}/G_N$"],
           "HU": [u"Felmérés", r"$z$", r"$f\sigma_8$ mért", r"$f\sigma_8$ M4", "pull M4", r"$f\sigma_8$ (72.8, 0.2448)", r"$\lvert\delta f\sigma_8\rvert$ variációs", r"$\Delta G_{\rm eff}/G_N$"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|:------------------------------|-----:|-----------:|-------:|------:|-------:|----------:|----------:|"]
    for r in RSD:
        out.append("| " + " | ".join([r["survey"], f(r["z"], 3), f"${f(r['fsigma8_obs'],3)}\\pm{f(r['fsigma8_err'],3)}$", f(r["fsigma8_Step25_M4_TAL"], 4),
                    f"{float(r['pull_M4_sigma']):+.2f}", f(r["fsigma8_Flat72p8_Om0p245"], 4), f"${sci(r['delta_fsigma8_var_can_z0norm'],2)}$",
                    f"${sci(r['delta_Geff_over_GN_can_at_k0p1h'],2)}$"]) + " |")
    return "\n".join(out)

def table_bbn(lang):
    lab = {"T_0p02MeV_end_BBN": ("end of BBN, $T\\approx 0.02$ MeV", u"BBN vége, $T\\approx 0.02$ MeV"),
           "T_0p09MeV_D_bottleneck": ("D bottleneck, $T\\approx 0.09$ MeV", u"D-szűk keresztmetszet, $T\\approx 0.09$ MeV"),
           "T_0p23MeV_np_freezeout": ("$z=10^9$, $T\\approx 0.23$ MeV", "$z=10^9$, $T\\approx 0.23$ MeV")}
    i = 0 if lang == "EN" else 1
    hdr = {"EN": ["Epoch", r"$z$", r"$H/(c/L)$", r"$\rho_m$ [kg/m$^3$]", r"$\lvert\bar\phi_{\rm attr}\rvert$", r"$\delta H/H$ (scalar energy)", r"$\delta G/G$ (clock)", r"$\lvert\delta Y_p\rvert$"],
           "HU": [u"Korszak", r"$z$", r"$H/(c/L)$", r"$\rho_m$ [kg/m$^3$]", r"$\lvert\bar\phi_{\rm attr}\rvert$", u"$\\delta H/H$ (skalárenergia)", u"$\\delta G/G$ (óra)", r"$\lvert\delta Y_p\rvert$"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|:------------------------------|------:|------:|------:|------:|------:|------:|------:|"]
    for r in BBN:
        out.append("| " + " | ".join([lab[r["epoch_id"]][i], f"${sci(r['z_epoch'],1)}$", f"${sci(r['H_over_m_phi_c'],3)}$", f"${sci(r['rho_m_kg_m3'],3)}$",
                   f"${sci(r['abs_phi_attr'],3)}$", f"${sci(r['delta_H_over_H_scalar_energy'],3)}$", f"${sci(r['delta_G_over_G_clock'],3)}$", f"${sci(r['abs_shift_Yp_TAL'],3)}$"]) + " |")
    return "\n".join(out)

def table_loc(lang):
    lab = {"1_Cosmic_Background_z0": ("Homogeneous FLRW, $z=0$", u"Homogén FLRW, $z=0$"),
           "2_Local_Cosmic_Void": (r"Cosmic void ($0.2\rho_{m,0}$)", u"Kozmikus üreg ($0.2\\rho_{m,0}$)"),
           "3_Milky_Way_Disk_Solar_Nbr": ("Milky Way disk (average)", u"Tejútrendszer-korong (átlag)"),
           "4_Solar_Mean_Density_Bulk": ("Mean solar density (bulk)", u"Nap átlagsűrűség (bulk)"),
           "5_Earth_Mean_Density_Bulk": ("Mean Earth density (bulk)", u"Föld átlagsűrűség (bulk)"),
           "6_White_Dwarf_Core_Bulk": ("White dwarf interior (bulk)", u"Fehér törpe belseje (bulk)"),
           "7_V21_TopHat_Interior_External_Benchmark": ("V21 top-hat interior external benchmark", u"V21 top-hat belső külső benchmark")}
    i = 0 if lang == "EN" else 1
    hdr = {"EN": ["Benchmark", r"$\rho$ [kg/m$^3$]", r"$\lvert\phi\rvert_{\rm bulk}=\alpha_Y\gamma_L\rho$", r"$\lvert\phi\rvert$ ($R\gg L$ edge)", r"${\rm expm1}(-\alpha_Y\lvert\phi\rvert)$", r"$\Delta H_0^{\rm app}$ [km/s/Mpc]", "fraction of 5.4 gap"],
           "HU": [u"Benchmark", r"$\rho$ [kg/m$^3$]", r"$\lvert\phi\rvert_{\rm bulk}=\alpha_Y\gamma_L\rho$", u"$\\lvert\\phi\\rvert$ ($R\\gg L$ perem)", r"${\rm expm1}(-\alpha_Y\lvert\phi\rvert)$", r"$\Delta H_0^{\rm app}$ [km/s/Mpc]", u"az 5.4-es rés hányada"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|:------------------------------|------:|------:|------:|------:|------:|------:|"]
    for r in LOC:
        def s(x):
            if x in ("", None) or (isinstance(x, str) and x.lower() == "nan"):
                return "n.a."
            x = float(x)
            if math.isnan(x):
                return "n.a."
            return f"${sci(x,3)}$" if (abs(x) < 1e-3 or abs(x) > 1e4) else f"${f(x,4)}$"
        out.append("| " + " | ".join([lab[r["env_id"]][i], s(r["rho_kg_m3"]), s(r["abs_phi_bulk"]), s(r["abs_phi_surface_R_gg_L"]),
                   s(r["f_clock_minus_1_expm1"]), s(r["Delta_H0_apparent_bulk_km_s_Mpc"]), s(r["fraction_of_5p4_Hubble_gap_bulk"])]) + " |")
    return "\n".join(out)

NOGO_NAMES = {
 "1_Step25_M4_TAL_Baseline_Fixed": ("1. Fixed Step 25 M4 baseline", u"1. Rögzített Step 25 M4 alapág"),
 "2_Profiled_M4_Shape_P_Free": (r"2. M4 shape, $P$ profiled", u"2. M4-alak, $P$ profilozva"),
 "3_Profiled_2Param_OmegaM_and_P_Free": (r"3. $\Omega_m$ and $P$ profiled", u"3. $\\Omega_m$ és $P$ profilozva"),
 "4_Joint_3Param_BAO_PlanckRd_LocalH0": (r"4. Joint BAO + Planck $r_d$ + assumed local $H_0$ prior", u"4. Közös BAO + Planck $r_d$ + feltételezett lokális $H_0$-prior"),
 "5_Rigid_LateTime_H0_72p80_FixedRd": (r"5. Rigid all-$z$ $H_0=72.8$, unchanged $r_d$", u"5. Merev, minden $z$-n $H_0=72.8$, változatlan $r_d$"),
 "6_Rigid_AllZ_H0_72p80_ScaledRd_FixedP": (r"6. Fixed-$P$ illustration", u"6. Rögzített-$P$ szemléltetés"),
 "7_Flat_H0_72p80_Om0p244828_SelfConsistent": (r"7. Flat 72.8, $\Omega_m=0.2448$, re-int. $r_d$", u"7. Sík 72.8, $\\Omega_m=0.2448$, újraint. $r_d$"),
 "8_Open_H0_72p80_FixedOmegaCB_OmK0p0098": (r"8. Open 72.8, fixed $\omega_{cb}$", u"8. Nyílt 72.8, rögzített $\\omega_{cb}$"),
}

def table_nogo(lang):
    i = 0 if lang == "EN" else 1
    hdr = {"EN": ["Scenario", r"$H_0$", r"$\Omega_m$", r"$\Omega_K$", r"$r_d$ [Mpc]", r"$P=H_0r_d$ [km/s]", r"$\chi^2_{\rm BAO}$ only (12×12)", r"$\Delta r_d$ vs Planck [Mpc]"],
           "HU": [u"Forgatókönyv", r"$H_0$", r"$\Omega_m$", r"$\Omega_K$", r"$r_d$ [Mpc]", r"$P=H_0r_d$ [km/s]", u"csak $\\chi^2_{\\rm BAO}$ (12×12)", u"$\\Delta r_d$ Planckhoz [Mpc]"]}[lang]
    out = ["| " + " | ".join(hdr) + " |", "|:------------------------------------|-----:|------:|------:|----------:|----------:|-------:|-------:|"]
    for r in NOGO:
        rd = f(r["rd_used_Mpc"], 2) + (f" $\\pm$ {f(r['rd_err_Mpc'],2)}" if r["rd_err_Mpc"] and float(r["rd_err_Mpc"]) > 0 else "")
        if r.get("rd_err_kind") == "input_prior_sigma":
            rd = f(r["rd_used_Mpc"], 2) + " (" + ("prior $\\sigma$ 0.26" if lang == "EN" else u"prior-$\\sigma$ 0.26") + ")"
        P = f(r["P_H0_rd_km_s"], 2) + (f" $\\pm$ {f(r['P_H0_rd_err_km_s'],2)}" if r["P_H0_rd_err_km_s"] and r["P_H0_rd_err_km_s"] != "nan" and float(r["P_H0_rd_err_km_s"]) > 0 else "")
        drd = f"{float(r['rd_used_Mpc']) - 147.09:+.2f}"
        out.append("| " + " | ".join([NOGO_NAMES[r["scenario_id"]][i], f(r["H0_km_s_Mpc"], 2), f(r["Omega_m_0"], 4), f(r["Omega_K_0"], 4), rd, P,
                   f(r["chi2_DESI_BAO_Cobaya_full"], 3), drd]) + " |")
    return "\n".join(out)

def filemeta(name):
    p = os.path.join(D, name)
    b = open(p, "rb").read()
    return f"`{name}`: {len(b):,} B, SHA-256 `{hashlib.sha256(b).hexdigest()}`"

DATA_FILES = ["desi_2024_gaussian_bao_ALL_GCcomb_mean.txt", "desi_2024_gaussian_bao_ALL_GCcomb_cov.txt",
              "tal_cosmology_common.py", "tal_hubble_desi_bao_summary.json"]

def subst(tpl, lang):
    tables = {"TABLE_BAO": table_bao(lang), "TABLE_PROF": table_prof(lang), "TABLE_CMB": table_cmb(lang), "TABLE_SNE": table_sne(lang),
              "TABLE_RSD": table_rsd(lang), "TABLE_BBN": table_bbn(lang), "TABLE_LOC": table_loc(lang), "TABLE_NOGO": table_nogo(lang),
              "DATA_FILES": "\n".join("- " + filemeta(n) for n in DATA_FILES),
              "TABLE_R2": table_r2(lang), "TABLE_R5": table_r5(lang), "TABLE_NU": table_nu(lang)}
    def rep(m):
        k = m.group(1)
        if k in tables: return tables[k]
        return V[k]
    out = re.sub(r"<<([A-Za-z0-9_]+)>>", rep, tpl)
    assert "<<" not in out, "unresolved placeholder"
    return out

# Review 2 (F01-F12) values, all read from the regenerated JSON
F2 = RV["review2_F01_F12"]
V.update({
    "Gratio": f(F2["F02_G_dictionary"]["G_eff_over_G_lab_large_scale"], 7),
    "qval": f(F2["F02_G_dictionary"]["q"], 7),
    "rsM3": f(MB["rs_zstar_m3_own_Mpc"], 8), "rdM3": f(MB["rd_zdrag_m3_own_Mpc"], 8),
    "thM3x": f(MB["theta100_m3"], 8), "chi2M3x": f(MB["chi2_bao_m3_full"], 5), "pullM3x": f(MB["formal_pull_theta_m3"], 2),
    "thM3old": f(MB["theta100_m3_with_frozen_M4_horizon_superseded"], 8),
    "chi2M3old": f(MB["chi2_bao_m3_full_with_frozen_M4_horizon_superseded"], 5),
    "thLowRsOnly": f(F2["F12_zstar_shift_diagnostics"]["low_Om_rs_only_shifted_fixed_DM"], 8),
    "thLowBoth": f(F2["F12_zstar_shift_diagnostics"]["low_Om_both_endpoints_shifted"], 8),
    "zsLowShift": f(F2["F12_zstar_shift_diagnostics"]["low_Om_z_star_shifted"], 4),
    "chi2LowSh": f(F2["F09_sne_figure"]["chi2_low_Om_shifted"], 6),
    "dfsCanF": sci(F2["F07_rsd_legend"]["max_abs_delta_fsigma8_var_canonical"], 4),
    "ns_dH0x": f(F2["F10_v21_interior"]["Delta_H0_bulk_km_s_Mpc"], 4),
    # Review N02: the two separate scalar effects quoted in Section 1.7, read from the JSON
    "N02dG": sci(R4["delta_Geff_over_GN_z0_canonical"], 1),
    "N02dfs": sci(R4["max_abs_delta_fsigma8_variational_canonical"], 2),
})

from step26_md_templates import TEMPLATE_EN, TEMPLATE_HU, README_TPL

def main():
    en = subst(TEMPLATE_EN, "EN"); hu = subst(TEMPLATE_HU, "HU")
    for name, txt in [("step26_desi_bao_planck_sne_rsd_bbn_nogo_en.md", en), ("step26_desi_bao_planck_sne_rsd_bbn_nogo_hu.md", hu),
                      ("v27_cosmology_desi_planck_sne_rsd_audit_en.md", en), ("v27_cosmology_desi_planck_sne_rsd_audit_hu.md", hu)]:
        open(os.path.join(D, name), "w", encoding="utf-8").write(txt)
        print(f"wrote {name}: {len(txt.encode('utf-8'))} B")
    rd = subst(README_TPL, "EN")
    open(os.path.join(D, "README.md"), "w", encoding="utf-8").write(rd)
    print(f"wrote README.md: {len(rd.encode('utf-8'))} B")

if __name__ == "__main__":
    main()
