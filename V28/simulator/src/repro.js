/* Reproduction table: JS recomputation vs package JSON values (shared by the page and verify_node.js). */
(function (root) {
  "use strict";
  root.TAL_reproRows = function (co) {
    co = co || {}; var T = root.TAL, D = root.TAL_DATA, R = D.ref;
    var c = co.ctl || T.baoControls(), f = co.fit || T.profile2p(), j = T.jointFit(), m3 = co.m3 || T.m3Branch(), cm = T.cmbForBackground(T.BG_M4);
    var bS = T.makeBackground(72.8, T.OMEGA_M_0), cS = T.cmbForBackground(bS), bL = T.makeBackground(72.8, 0.244828), cL = T.cmbForBackground(bL);
    var sn = T.sneControl(T.BG_M4), snL = T.sneControl(bL), g1 = T.growth(T.BG_M4, {}), g2 = T.growth(T.BG_M4, { cbOnly: false }), g3 = T.growth(bL, {});
    var r1 = T.rsdChi2(g1), bb = T.bbn(1e8, 10, 1), AS = 180 / Math.PI * 3600, s = R.sne.symmetric_zero_point_R5, J = R.joint, M = R.m3;
    var rows = [
      ["r_d(z_drag) M4 [Mpc]", T.RD_ZDRAG, R.rd_zdrag], ["r_s(z*) M4 [Mpc]", T.RS_ZSTAR, R.rs_zstar],
      ["χ² M4 (12×12)", c.chi2_m4, R.chi2_m4], ["χ² M4 (Table 1)", c.chi2_m4_tab1, R.chi2_m4_tab1],
      ["χ² rigid H0 = 72.8 (12×12)", c.chi2_rigid, R.chi2_rigid], ["χ² rigid H0 = 72.8 (Table 1)", c.chi2_rigid_tab1, R.chi2_rigid_tab1],
      ["1-param P [km/s]", c.P1, R.P1], ["1-param σ_P", c.sigmaP1, R.sigmaP1], ["1-param χ²", c.chi2_P1, R.chi2_P1], ["1-param r_d req", c.rdReq1, R.rdReq1],
      ["2-param Ω_m", f.om, R.om2], ["2-param σ(Ω_m)", f.sigmaOm, R.sigmaOm2], ["2-param P", f.P, R.P2], ["2-param σ_P", f.sigmaP, R.sigmaP2], ["2-param χ²", f.chi2, R.chi2_2p],
      ["r_d req @ 72.8", f.rdReq, R.rdReq2], ["σ r_d req (fixed H0)", f.sigmaRdFixedH0, R.sigmaRdReq2], ["σ r_d req (with H0)", f.sigmaRdWithH0, R.sigmaRdReq2H0],
      ["H0 from Planck r_d", f.H0fromPlanckRd, R.H0PlanckRd], ["σ H0 from Planck r_d", f.sigmaH0fromPlanckRd, R.sigmaH0PlanckRd],
      ["joint Ω_m", j.om, J.best_fit_Omega_m_H0_rd[0]], ["joint H0", j.h0, J.best_fit_Omega_m_H0_rd[1]], ["joint r_d", j.rd, J.best_fit_Omega_m_H0_rd[2]],
      ["joint χ²_BAO", j.parts.bao, J.chi2_BAO], ["joint χ²_rd", j.parts.rd, J.chi2_rd_prior], ["joint χ²_H0", j.parts.h0, J.chi2_H0_prior], ["joint χ²", j.chi2, J.chi2_joint],
      ["joint Δχ²", j.chi2 - f.chi2, J.Delta_chi2_vs_BAO_MLE],
      ["M3 r_s own", m3.rs, M.rs_zstar_m3_own_Mpc], ["M3 r_d own", m3.rd, M.rd_zdrag_m3_own_Mpc], ["M3 100θ*", m3.theta100, M.theta100_m3], ["M3 θ pull", m3.pullTheta, M.formal_pull_theta_m3],
      ["M3 χ²_BAO", m3.chi2_bao, M.chi2_bao_m3_full], ["M3 Ω_K", m3.om_k, M.om_k_m3], ["M3 Ω_DE", m3.om_de, M.om_de_m3], ["M2 E²(0)", m3.E2_m2_z0, M.E2_m2_z0], ["M2 root z", m3.z_root_m2, M.z_root_full],
      ["M4 100θ*", cm.theta100, R.theta_m4], ["M4 D_M(z*)", cm.dm, R.dm_zstar_m4], ["M4 Planck 3D χ²", cm.chi2_3d, R.planck3d_m4],
      ["Flat 72.8 same Ω_m r_d", cS.rd, R.flat728_sameOm.rd_zdrag_reint_Mpc], ["Flat 72.8 same Ω_m χ²_BAO", T.chi2BAO(bS.Hz, cS.rd), R.flat728_sameOm.chi2_DESI_BAO_Cobaya_with_rd_reint],
      ["Flat 72.8 Ω_m 0.2448 r_d", cL.rd, R.flat728_lowOm.rd_zdrag_reint_Mpc], ["Flat 72.8 Ω_m 0.2448 χ²_BAO", T.chi2BAO(bL.Hz, cL.rd), R.flat728_lowOm.chi2_DESI_BAO_Cobaya_with_rd_reint], ["Flat 72.8 Ω_m 0.2448 100θ*", cL.theta100, R.flat728_lowOm.theta_star_100_reint],
      ["SNe offset M4", sn.offset, s.M4.best_offset_mag], ["SNe χ² M4 unshifted", sn.chi2_unshifted, s.M4.chi2_without_offset], ["SNe χ² M4 shifted", sn.chi2_shifted, s.M4.chi2_with_offset],
      ["SNe offset low Ω_m", snL.offset, s.lowOm_0p2448.best_offset_mag], ["SNe χ² low Ω_m unshifted", snL.chi2_unshifted, s.lowOm_0p2448.chi2_without_offset], ["SNe χ² low Ω_m shifted", snL.chi2_shifted, s.lowOm_0p2448.chi2_with_offset],
      ["RSD χ² M4 (cb)", r1.chi2, R.rsd.chi2_RSD_diagonal_M4_clustering_cb], ["RSD χ² M4 (all m)", T.rsdChi2(g2).chi2, R.rsd.chi2_RSD_diagonal_M4_all_m], ["RSD χ² Flat 72.8 Ω_m 0.2448", T.rsdChi2(g3).chi2, R.rsd.chi2_RSD_diagonal_Flat72p80_Om0p244828],
      ["RSD max |δfσ8|", r1.maxDelta, R.rsd.max_abs_delta_fsigma8_variational_canonical], ["RSD kL(z=0)", g1.kL0, R.rsd.kL_z0],
      ["BBN H(z=10⁸) [s⁻¹]", bb.H, R.bbn[0].H_BBN_s_inv], ["BBN H/(m_φ c)", bb.H_over_mphi, R.bbn[0].H_over_m_phi_c], ["BBN |φ|", bb.absPhi, R.bbn[0].abs_phi_attr], ["BBN δH/H", bb.dHH, R.bbn[0].delta_H_over_H_scalar_energy], ["BBN Y_p", bb.Yp, R.bbn[0].Yp_TAL_pred],
      ["γ_L [m³/kg]", T.gammaL(10), R.gammaL], ["G_eff/G_lab (r ≫ L)", T.gDictionary(1, 0).Geff_over_Glab_large, R.G_ratio]
    ];
    R.solar.forEach(function (row) {
      var d = T.deflectionSchw(row.b_km, row.GM_over_c2_km);
      rows.push(["Sun α̂ exact, " + row.scale_label.replace("Solar_1Msun_", "").replace(/_/g, " ") + " [″]", d.alpha * AS, row.alpha_hat_exact_arcsec]);
      if (row.Cassini_roundtrip_log_delay_us !== null) rows.push(["Cassini round trip, " + row.scale_label.replace("Solar_1Msun_", "").replace(/_/g, " ") + " [μs]", T.cassiniDelay(row.b_km, row.GM_over_c2_km).rt, row.Cassini_roundtrip_log_delay_us]);
    });
    var r0 = T.deflectionSchw(R.solar[0].b_km, R.solar[0].GM_over_c2_km).r0;
    rows.push(["Sun Shapiro (10 r0) [μs]", T.shapiroSchw(r0, 10 * r0, R.solar[0].GM_over_c2_km), R.solar[0].Shapiro_two_sided_10r0_us]);
    var Mns = R.ns_ref.M_ADM_Msun * T.GM_SUN_KM;
    R.ns_lens.forEach(function (row) { rows.push(["NS α̂ Schw M_ADM, b = " + (+row.b_km).toPrecision(5) + " km [rad]", T.deflectionSchw(+row.b_km, Mns).alpha, +row.alpha_hat_Schw_MADM_rad]); });
    rows.push(["f_cut [Hz]", T.gwScalar(100, 10, 40).fcut, R.gw_meta.f_cut_Hz]);
    R.gw_rows.forEach(function (r) {
      var f = +r.frequency_Hz, g = T.gwScalar(f, 10, 40), lab = String(+f.toPrecision(6)) + " Hz";
      if (r.scalar_regime.indexOf("propagating") === 0) { rows.push(["scalar v_g/c @ " + lab, g.vg, +r.scalar_v_group_over_c]); rows.push(["scalar delay over 40 Mpc @ " + lab + " [s]", g.delay, +r.scalar_delay_at_40Mpc_s]); }
      else if (r.scalar_regime.indexOf("evanescent") === 0 && /LIGO|PSR_B|Sub_cutoff/.test(r.mode_label)) rows.push(["log10 scalar amplitude @ 40 Mpc, " + lab, g.log10amp, +r.log10_scalar_amp_at_40Mpc]);
    });
    var n3 = D.sparc.filter(function (g) { return g.n === "NGC3198"; })[0];
    rows.push(["NGC 3198 χ² (Υ = 0.5)", T.galaxyModel(n3, 0.5, 0.5, 10 / T.KPC_KM, 1).chi2_bar, 32639.3]);
    rows.push(["NGC 3198 χ² (Υ = 1.0)", T.galaxyModel(n3, 1.0, 1.0, 10 / T.KPC_KM, 1).chi2_bar, 12645.3]);
    var cq = T.galaxyModel(n3, 0.5, 0.5, 10 / T.KPC_KM, 1).pts.filter(function (p) { return p.R === 30.08; })[0];
    rows.push(["NGC 3198 C_req @ 30.08 kpc (Υ = 0.5)*", cq.Creq, R.macro.NGC3198_C_req_Upsilon_0p5]);
    R.local.forEach(function (e) { var dH = e.rho_kg_m3 === null ? T.H0_STEP25 * Math.expm1(e.abs_phi_bulk) : T.localBulk(e.rho_kg_m3, 10, 1).dH0; rows.push(["ΔH0 bulk, " + e.env_id.replace(/^\d+_/, "").replace(/_/g, " "), dH, e.Delta_H0_apparent_bulk_km_s_Mpc]); });
    return rows;
  };
})(typeof window !== "undefined" ? window : globalThis);
