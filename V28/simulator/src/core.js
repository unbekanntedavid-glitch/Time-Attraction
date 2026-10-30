/* Time Attraction Law V27 simulator: numerical core (port of tal_cosmology_common.py and Step 23/26 scripts). */
(function (root) {
  "use strict";
  var T = {};
  // ---------- constants (SI) ----------
  var G_SI = 6.67430e-11, C_SI = 299792458.0, C_KM_S = C_SI / 1000.0;
  var MPC_M = 3.085677581491367e22;
  T.G_SI = G_SI; T.C_SI = C_SI; T.C_KM_S = C_KM_S; T.MPC_M = MPC_M;
  T.MSUN_KG = 1.98840987e30;
  T.GM_SUN_KM = 1.476625038050125; // GM_sun/c^2 in km (V26 Step 23)
  T.R_SUN_KM = 695700.0; T.AU_KM = 149597870.7; T.PC_KM = 3.0856775814913673e13;
  T.KPC_KM = 1e3 * T.PC_KM; T.MPC_KM = 1e6 * T.PC_KM; T.YR_S = 365.25 * 86400.0;
  T.HBAR_C_EV_M = 1.973269804e-7;

  // ---------- TAL parameters ----------
  T.L_KM_DEF = 10.0; T.ALPHA_Y_DEF = 1.0; T.BETA = 1.0 / (4.0 * Math.PI);
  T.gammaL = function (L_km) { var L = L_km * 1000; return 4 * Math.PI * T.BETA * G_SI * L * L / (C_SI * C_SI); };
  T.qCoupling = function (alphaY) { return alphaY * alphaY * T.BETA; };

  // ---------- cosmology baseline ----------
  var H0_STEP25 = 67.40, H0_PLANCK = 67.36, H0_LOCAL = 72.80, H0_LOCAL_SIGMA = 1.0;
  var OMEGA_B_H2 = 0.02237, OMEGA_C_H2 = 0.12000, OMEGA_NU_H2 = 0.000645, OMEGA_GAMMA_H2 = 2.47282e-5, N_EFF = 3.046;
  var Z_STAR = 1089.80, Z_DRAG = 1059.94;
  var THETA100_PL = 1.04109, THETA100_PL_ERR = 0.00030, RD_PL = 147.09, RD_PL_ERR = 0.26;
  T.H0_STEP25 = H0_STEP25; T.H0_PLANCK = H0_PLANCK; T.H0_LOCAL = H0_LOCAL; T.H0_LOCAL_SIGMA = H0_LOCAL_SIGMA;
  T.OMEGA_B_H2 = OMEGA_B_H2; T.OMEGA_C_H2 = OMEGA_C_H2; T.OMEGA_NU_H2 = OMEGA_NU_H2;
  T.Z_STAR = Z_STAR; T.Z_DRAG = Z_DRAG; T.RD_PL = RD_PL; T.RD_PL_ERR = RD_PL_ERR;
  T.THETA100_PL = THETA100_PL; T.THETA100_PL_ERR = THETA100_PL_ERR;

  // ---------- Gauss-Legendre nodes ----------
  function glNodes(n) {
    var x = new Float64Array(n), w = new Float64Array(n), m = (n + 1) >> 1;
    for (var i = 0; i < m; i++) {
      var z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), z1, pp;
      do {
        var p1 = 1, p2 = 0;
        for (var j = 0; j < n; j++) { var p3 = p2; p2 = p1; p1 = ((2 * j + 1) * z * p2 - j * p3) / (j + 1); }
        pp = n * (z * p1 - p2) / (z * z - 1); z1 = z; z = z1 - p1 / pp;
      } while (Math.abs(z - z1) > 1e-15);
      x[i] = -z; x[n - 1 - i] = z; w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * pp * pp);
    }
    return { x: x, w: w };
  }
  T.glNodes = glNodes;
  var GL10 = glNodes(10), GL20 = glNodes(20), GL128 = glNodes(128);
  function glPanel(f, a, b, g) {
    var c = 0.5 * (a + b), h = 0.5 * (b - a), s = 0;
    for (var i = 0; i < g.x.length; i++) s += g.w[i] * f(c + h * g.x[i]);
    return s * h;
  }
  // adaptive Gauss-Legendre (10 vs 20 points, bisection)
  function integrate(f, a, b, tol, depth) {
    if (a === b) return 0;
    tol = tol || 1e-12; depth = depth || 0;
    var s1 = glPanel(f, a, b, GL10), s2 = glPanel(f, a, b, GL20);
    if (Math.abs(s2 - s1) <= tol * Math.max(Math.abs(s2), 1e-300) || depth > 40) return s2;
    var m = 0.5 * (a + b);
    return integrate(f, a, m, tol, depth + 1) + integrate(f, m, b, tol, depth + 1);
  }
  T.integrate = integrate;

  // ---------- background ----------
  function nuMassiveE2(z, om_nu0, om_nu_rel1) {
    if (om_nu0 <= 0) return om_nu_rel1 * Math.pow(1 + z, 4);
    var a = 1 / (1 + z), p = 1.83, ra = (om_nu0 / om_nu_rel1) * a;
    return om_nu_rel1 * Math.pow(1 + Math.pow(ra, p), 1 / p) * Math.pow(1 + z, 4);
  }
  T.nuMassiveE2 = nuMassiveE2;
  function E2general(z, p) {
    var zp1 = 1 + z, de;
    if (p.w0 === undefined || (p.w0 === -1 && !p.wa)) de = p.om_de;
    else de = p.om_de * Math.pow(zp1, 3 * (1 + p.w0 + p.wa)) * Math.exp(-3 * p.wa * z / zp1);
    return p.om_r * zp1 * zp1 * zp1 * zp1 + nuMassiveE2(z, p.om_nu, p.om_nu_r1) + p.om_cb * zp1 * zp1 * zp1 + p.om_k * zp1 * zp1 + de;
  }
  T.E2general = E2general;
  // make_model_background (exact z=0 normalisation by default)
  function makeBackground(H0, om_m, opts) {
    opts = opts || {};
    var om_b_h2 = opts.om_b_h2 !== undefined ? opts.om_b_h2 : OMEGA_B_H2;
    var om_nu_h2 = opts.om_nu_h2 !== undefined ? opts.om_nu_h2 : OMEGA_NU_H2;
    var om_k = opts.om_k || 0, exact = opts.exact !== false;
    var h = H0 / 100, h2 = h * h;
    var om_b = om_b_h2 / h2, om_nu = om_nu_h2 / h2, om_c, om_cb;
    if (om_m === null || om_m === undefined) { om_c = OMEGA_C_H2 / h2; om_cb = om_b + om_c; om_m = om_cb + om_nu; }
    else { om_cb = om_m - om_nu; om_c = om_cb - om_b; }
    var om_g = OMEGA_GAMMA_H2 / h2, om_nu_r1 = om_g * 0.22710731766 * (N_EFF / 3.0), om_r = om_g + 2 * om_nu_r1;
    var om_nu_z0 = nuMassiveE2(0, om_nu, om_nu_r1);
    var om_de = exact ? (1 - om_k - om_cb - om_r - om_nu_z0) : (1 - om_k - (om_cb + om_nu) - om_r);
    if (opts.om_de !== undefined) om_de = opts.om_de;
    var p = { om_r: om_r, om_nu: om_nu, om_nu_r1: om_nu_r1, om_cb: om_cb, om_k: om_k, om_de: om_de };
    var bg = { H0: H0, h: h, Omega_b_0: om_b, Omega_c_0: om_c, Omega_cb_0: om_cb, Omega_nu_0: om_nu, Omega_m_0: om_m,
      Omega_gamma_0: om_g, Omega_r_massless_0: om_r, Omega_nu_rel_1_0: om_nu_r1, Omega_k_0: om_k, Omega_de_0: om_de,
      omega_b_h2: om_b_h2, omega_cb_h2: om_cb * h2, omega_m_h2: om_m * h2, Omega_nu_interp_z0: om_nu_z0, params: p };
    bg.E2 = function (z) { return E2general(z, p); };
    bg.Hz = function (z) { return H0 * Math.sqrt(E2general(z, p)); };
    return bg;
  }
  T.makeBackground = makeBackground;

  function soundSpeed(z, om_b, om_g) { var rb = (3 * om_b) / (4 * om_g * (1 + z)); return 1 / Math.sqrt(3 * (1 + rb)); }
  // r_s(z_lower) = c/H0 * int_{z_lower}^{inf} c_s(z)/E(z) dz, evaluated as int_0^{a_l} c_s/(a^2 E) da
  function soundHorizon(zl, bg) {
    var al = 1 / (1 + zl), ob = bg.Omega_b_0, og = bg.Omega_gamma_0, E2 = bg.E2;
    var f = function (a) { if (a <= 0) a = 1e-300; var z = 1 / a - 1; return soundSpeed(z, ob, og) / (a * a * Math.sqrt(E2(z))); };
    // a -> 0 limit handled by geometric panels; integrand is finite there
    var s = 0, hi = al, lo;
    for (var k = 0; k < 40; k++) { lo = hi * 0.5; s += integrate(f, lo, hi, 1e-13); hi = lo; }
    s += integrate(f, 0, hi, 1e-13);
    return (C_KM_S / bg.H0) * s;
  }
  T.soundHorizon = soundHorizon;
  function chiOfZ(z, Hz) { return integrate(function (zp) { return C_KM_S / Hz(zp); }, 0, z, 1e-13); }
  function DMofZ(z, Hz, om_k, H0c) {
    var chi = chiOfZ(z, Hz);
    if (!om_k || Math.abs(om_k) < 1e-12) return chi;
    var dh0 = C_KM_S / (H0c || H0_STEP25), s = Math.sqrt(Math.abs(om_k)), x = s * chi / dh0;
    return om_k > 0 ? (dh0 / s) * Math.sinh(x) : (dh0 / s) * Math.sin(x);
  }
  T.DMofZ = DMofZ;
  T.zStarHuSugiyama = function (obh2, omh2) {
    var g1 = (0.0783 * Math.pow(obh2, -0.238)) / (1 + 39.5 * Math.pow(obh2, 0.763)), g2 = 0.560 / (1 + 21.1 * Math.pow(obh2, 1.81));
    return 1048.0 * (1 + 0.00124 * Math.pow(obh2, -0.738)) * (1 + g1 * Math.pow(omh2, g2));
  };
  T.zDragEH = function (obh2, omh2) {
    var b1 = 0.313 * Math.pow(omh2, -0.419) * (1 + 0.607 * Math.pow(omh2, 0.674)), b2 = 0.238 * Math.pow(omh2, 0.223);
    return 1291.0 * (Math.pow(omh2, 0.251) / (1 + 0.659 * Math.pow(omh2, 0.828))) * (1 + b1 * Math.pow(obh2, b2));
  };

  // ---------- baseline objects ----------
  var h674 = H0_STEP25 / 100;
  var OMEGA_M_0 = (OMEGA_B_H2 + OMEGA_C_H2 + OMEGA_NU_H2) / (h674 * h674);
  T.OMEGA_M_0 = OMEGA_M_0;
  T.BG_M4 = makeBackground(H0_STEP25, OMEGA_M_0);
  T.BG_PLANCK = makeBackground(H0_PLANCK, null);
  T.RS_ZSTAR = soundHorizon(Z_STAR, T.BG_M4);
  T.RD_ZDRAG = soundHorizon(Z_DRAG, T.BG_M4);
  var H0_S_INV = H0_STEP25 * 1000 / MPC_M;
  T.H0_S_INV = H0_S_INV;
  T.RHO_CRIT_0 = 3 * H0_S_INV * H0_S_INV / (8 * Math.PI * G_SI);
  T.OMEGA_LAMBDA_FLAT_0 = T.BG_M4.Omega_de_0;

  // ---------- linear algebra ----------
  function invert(Ain) {
    var n = Ain.length, A = Ain.map(function (r) { return r.slice(); }), I = [], i, j, k;
    for (i = 0; i < n; i++) { I.push([]); for (j = 0; j < n; j++) I[i].push(i === j ? 1 : 0); }
    for (i = 0; i < n; i++) {
      var piv = i; for (k = i + 1; k < n; k++) if (Math.abs(A[k][i]) > Math.abs(A[piv][i])) piv = k;
      var t = A[i]; A[i] = A[piv]; A[piv] = t; t = I[i]; I[i] = I[piv]; I[piv] = t;
      var d = A[i][i];
      for (j = 0; j < n; j++) { A[i][j] /= d; I[i][j] /= d; }
      for (k = 0; k < n; k++) if (k !== i) { var fct = A[k][i]; if (fct !== 0) for (j = 0; j < n; j++) { A[k][j] -= fct * A[i][j]; I[k][j] -= fct * I[i][j]; } }
    }
    return I;
  }
  T.invert = invert;
  function quadForm(a, M, b) { var s = 0; for (var i = 0; i < a.length; i++) { var r = 0; for (var j = 0; j < b.length; j++) r += M[i][j] * b[j]; s += a[i] * r; } return s; }
  T.quadForm = quadForm;

  // ---------- DESI DR1 BAO (12-vector + 12x12 covariance) ----------
  var D = root.TAL_DATA;
  T.DESI = D.desi;
  T.DESI_INVCOV = invert(D.desi.cov);
  T.DESI_SIG = D.desi.cov.map(function (r, i) { return Math.sqrt(r[i]); });

  function baoObs(z, Hz, rd, om_k, H0c) {
    var dm = DMofZ(z, Hz, om_k, H0c) / rd, dh = (C_KM_S / Hz(z)) / rd, dv = Math.pow(z * dm * dm * dh, 1 / 3);
    return { DM: dm, DH: dh, DV: dv };
  }
  T.baoObs = baoObs;
  function baoPred(Hz, rd, om_k, H0c) {
    var cache = {}, out = [];
    for (var i = 0; i < D.desi.z.length; i++) {
      var z = D.desi.z[i]; if (!cache[z]) cache[z] = baoObs(z, Hz, rd, om_k, H0c);
      var q = D.desi.q[i]; out.push(q === "DM_over_rs" ? cache[z].DM : (q === "DH_over_rs" ? cache[z].DH : cache[z].DV));
    }
    return out;
  }
  T.baoPred = baoPred;
  function chi2Vec(pred) { var d = pred.map(function (v, i) { return v - D.desi.mean[i]; }); return quadForm(d, T.DESI_INVCOV, d); }
  T.chi2Vec = chi2Vec;
  T.chi2BAO = function (Hz, rd, om_k, H0c) { return chi2Vec(baoPred(Hz, rd, om_k, H0c)); };
  T.chi2BAOTable1 = function (Hz, rd, om_k, H0c) {
    var c2 = 0;
    D.desi.table1.forEach(function (b) {
      if (!b.primary) return;
      var o = baoObs(b.z, Hz, rd, om_k, H0c);
      if (b.mode === "iso") c2 += Math.pow((o.DV - b.DV) / b.DVe, 2);
      else { var dM = (o.DM - b.DM) / b.DMe, dH = (o.DH - b.DH) / b.DHe; c2 += (dM * dM - 2 * b.r * dM * dH + dH * dH) / (1 - b.r * b.r); }
    });
    return c2;
  };
  // unit-scale predictions (P = H0 r_d = 1 km/s) for a dimensionless E(z)
  function unitPred(Efn) {
    var cache = {}, out = [];
    for (var i = 0; i < D.desi.z.length; i++) {
      var z = D.desi.z[i];
      if (!cache[z]) { var dm = integrate(function (zp) { return C_KM_S / Efn(zp); }, 0, z, 1e-13), dh = C_KM_S / Efn(z); cache[z] = [dm, dh, Math.pow(z * dm * dm * dh, 1 / 3)]; }
      var q = D.desi.q[i]; out.push(q === "DM_over_rs" ? cache[z][0] : (q === "DH_over_rs" ? cache[z][1] : cache[z][2]));
    }
    return out;
  }
  T.unitPred = unitPred;
  function profileP(u) {
    var W = T.DESI_INVCOV, y = D.desi.mean, uiu = quadForm(u, W, u), uiy = quadForm(u, W, y);
    var a = uiy / uiu, P = 1 / a, sa = 1 / Math.sqrt(uiu), sP = sa / (a * a);
    var c2 = chi2Vec(u.map(function (v) { return v * a; }));
    return { P: P, sigmaP: sP, chi2: c2 };
  }
  T.profileP = profileP;
  function flatE(om) { return function (z) { return Math.sqrt(om * Math.pow(1 + z, 3) + (1 - om)); }; }
  T.flatE = flatE;
  // Brent minimiser (1D)
  function brent(f, ax, bx, cx, tol) {
    var CG = 0.3819660112501051, a = Math.min(ax, cx), b = Math.max(ax, cx), x = bx, w = bx, v = bx, fx = f(x), fw = fx, fv = fx, d = 0, e = 0;
    for (var it = 0; it < 500; it++) {
      var xm = 0.5 * (a + b), tol1 = tol * Math.abs(x) + 1e-15, tol2 = 2 * tol1;
      if (Math.abs(x - xm) <= tol2 - 0.5 * (b - a)) break;
      if (Math.abs(e) > tol1) {
        var r = (x - w) * (fx - fv), q = (x - v) * (fx - fw), p = (x - v) * q - (x - w) * r; q = 2 * (q - r);
        if (q > 0) p = -p; q = Math.abs(q); var et = e; e = d;
        if (Math.abs(p) >= Math.abs(0.5 * q * et) || p <= q * (a - x) || p >= q * (b - x)) { e = x >= xm ? a - x : b - x; d = CG * e; }
        else { d = p / q; var u0 = x + d; if (u0 - a < tol2 || b - u0 < tol2) d = xm - x >= 0 ? tol1 : -tol1; }
      } else { e = x >= xm ? a - x : b - x; d = CG * e; }
      var u = Math.abs(d) >= tol1 ? x + d : x + (d >= 0 ? tol1 : -tol1), fu = f(u);
      if (fu <= fx) { if (u >= x) a = x; else b = x; v = w; fv = fw; w = x; fw = fx; x = u; fx = fu; }
      else { if (u < x) a = u; else b = u; if (fu <= fw || w === x) { v = w; fv = fw; w = u; fw = fu; } else if (fu <= fv || v === x || v === w) { v = u; fv = fu; } }
    }
    return { x: x, f: fx };
  }
  T.brent = brent;
  // Nelder-Mead (scipy-like defaults)
  function nelderMead(f, x0, opt) {
    opt = opt || {}; var n = x0.length, xatol = opt.xatol || 1e-8, fatol = opt.fatol || 1e-10, maxit = opt.maxiter || 200 * n * 10;
    var sim = [x0.slice()];
    for (var k = 0; k < n; k++) { var y = x0.slice(); y[k] = y[k] !== 0 ? 1.05 * y[k] : 0.00025; sim.push(y); }
    var fs = sim.map(f), it = 0;
    function order() { var idx = fs.map(function (_, i) { return i; }).sort(function (a, b) { return fs[a] - fs[b]; }); sim = idx.map(function (i) { return sim[i]; }); fs = idx.map(function (i) { return fs[i]; }); }
    order();
    while (it++ < maxit) {
      var mx = 0, mf = 0;
      for (var i = 1; i <= n; i++) { mf = Math.max(mf, Math.abs(fs[i] - fs[0])); for (var j = 0; j < n; j++) mx = Math.max(mx, Math.abs(sim[i][j] - sim[0][j])); }
      if (mx <= xatol && mf <= fatol) break;
      var xb = new Array(n).fill(0);
      for (i = 0; i < n; i++) for (j = 0; j < n; j++) xb[j] += sim[i][j] / n;
      var xr = xb.map(function (v, j) { return 2 * v - sim[n][j]; }), fr = f(xr);
      if (fr < fs[0]) {
        var xe = xb.map(function (v, j) { return 3 * v - 2 * sim[n][j]; }), fe = f(xe);
        if (fe < fr) { sim[n] = xe; fs[n] = fe; } else { sim[n] = xr; fs[n] = fr; }
      } else if (fr < fs[n - 1]) { sim[n] = xr; fs[n] = fr; }
      else {
        var shrink = false;
        if (fr < fs[n]) { var xc = xb.map(function (v, j) { return 1.5 * v - 0.5 * sim[n][j]; }), fc = f(xc); if (fc <= fr) { sim[n] = xc; fs[n] = fc; } else shrink = true; }
        else { var xcc = xb.map(function (v, j) { return 0.5 * v + 0.5 * sim[n][j]; }), fcc = f(xcc); if (fcc < fs[n]) { sim[n] = xcc; fs[n] = fcc; } else shrink = true; }
        if (shrink) for (i = 1; i <= n; i++) { sim[i] = sim[i].map(function (v, j) { return sim[0][j] + 0.5 * (v - sim[0][j]); }); fs[i] = f(sim[i]); }
      }
      order();
    }
    return { x: sim[0], f: fs[0], iterations: it };
  }
  T.nelderMead = nelderMead;

  // 2-parameter (Omega_m, P) flat profile with Fisher errors
  T.profile2p = function () {
    var cacheU = {};
    function uOf(om) { var k = om.toPrecision(15); if (!cacheU[k]) cacheU[k] = unitPred(flatE(om)); return cacheU[k]; }
    var r = brent(function (om) { return profileP(uOf(om)).chi2; }, 0.2, 0.294, 0.4, 1e-11);
    var om = r.x, pr = profileP(uOf(om));
    function vec(o, P) { return unitPred(flatE(o)).map(function (v) { return v / P; }); }
    var eo = 1e-5, eP = 1e-2, P = pr.P;
    var a1 = vec(om + eo, P), a2 = vec(om - eo, P), b1 = vec(om, P + eP), b2 = vec(om, P - eP);
    var J0 = a1.map(function (v, i) { return (v - a2[i]) / (2 * eo); }), J1 = b1.map(function (v, i) { return (v - b2[i]) / (2 * eP); });
    var W = T.DESI_INVCOV, F = [[quadForm(J0, W, J0), quadForm(J0, W, J1)], [quadForm(J1, W, J0), quadForm(J1, W, J1)]], Ci = invert(F);
    var sOm = Math.sqrt(Ci[0][0]), sP = Math.sqrt(Ci[1][1]);
    var rdReq = P / H0_LOCAL, sRd = sP / H0_LOCAL;
    return { om: om, sigmaOm: sOm, P: P, sigmaP: sP, chi2: pr.chi2, rdReq: rdReq, sigmaRdFixedH0: sRd,
      sigmaRdWithH0: rdReq * Math.sqrt(Math.pow(sP / P, 2) + Math.pow(H0_LOCAL_SIGMA / H0_LOCAL, 2)),
      H0fromPlanckRd: P / RD_PL, sigmaH0fromPlanckRd: (P / RD_PL) * Math.sqrt(Math.pow(sP / P, 2) + Math.pow(RD_PL_ERR / RD_PL, 2)) };
  };
  // joint (Omega_m, H0, r_d): DESI BAO + Planck r_d prior + assumed local H0 benchmark prior
  T.jointChi2Parts = function (om, h0, rd, h0Prior, h0Sig, rdPrior, rdSig) {
    var u = unitPred(flatE(om)), c2b = chi2Vec(u.map(function (v) { return v / (h0 * rd); }));
    return { bao: c2b, rd: Math.pow((rd - rdPrior) / rdSig, 2), h0: Math.pow((h0 - h0Prior) / h0Sig, 2) };
  };
  T.jointFit = function (h0Prior, h0Sig, rdPrior, rdSig) {
    h0Prior = h0Prior === undefined ? H0_LOCAL : h0Prior; h0Sig = h0Sig || H0_LOCAL_SIGMA;
    rdPrior = rdPrior === undefined ? RD_PL : rdPrior; rdSig = rdSig || RD_PL_ERR;
    var cacheU = {};
    function f(x) {
      var om = x[0], h0 = x[1], rd = x[2];
      if (om <= 0.05 || om >= 0.95 || h0 <= 40 || rd <= 100) return 1e9;
      var k = om.toPrecision(14); if (!cacheU[k]) cacheU[k] = unitPred(flatE(om));
      var P = h0 * rd, c2b = chi2Vec(cacheU[k].map(function (v) { return v / P; }));
      return c2b + Math.pow((rd - rdPrior) / rdSig, 2) + Math.pow((h0 - h0Prior) / h0Sig, 2);
    }
    var r = nelderMead(f, [0.294, 70.5, 146.9], { xatol: 1e-8, fatol: 1e-10 });
    r = nelderMead(f, r.x, { xatol: 1e-9, fatol: 1e-11 });
    var parts = T.jointChi2Parts(r.x[0], r.x[1], r.x[2], h0Prior, h0Sig, rdPrior, rdSig);
    return { om: r.x[0], h0: r.x[1], rd: r.x[2], chi2: r.f, parts: parts };
  };
  // profile controls (BAO)
  T.baoControls = function () {
    var m4 = T.BG_M4, RD = T.RD_ZDRAG;
    var c2m4 = T.chi2BAO(m4.Hz, RD), c2m4t = T.chi2BAOTable1(m4.Hz, RD);
    var Hrig = function (z) { return (H0_LOCAL / H0_STEP25) * m4.Hz(z); };
    var c2rig = T.chi2BAO(Hrig, RD), c2rigt = T.chi2BAOTable1(Hrig, RD);
    var u4 = unitPred(function (z) { return Math.sqrt(m4.E2(z)); }), p1 = profileP(u4);
    return { chi2_m4: c2m4, chi2_m4_tab1: c2m4t, chi2_rigid: c2rig, chi2_rigid_tab1: c2rigt,
      P_m4: H0_STEP25 * RD, P_rigid: H0_LOCAL * RD, P1: p1.P, sigmaP1: p1.sigmaP, chi2_P1: p1.chi2,
      rdReq1: p1.P / H0_LOCAL, sigmaRdReq1: p1.sigmaP / H0_LOCAL, rdReqFixedP: H0_STEP25 * RD / H0_LOCAL };
  };
  // M3 uncompensated open control with own horizons
  T.m3Branch = function () {
    var m4 = T.BG_M4, p = m4.params, om_de = -m4.Omega_de_0;
    var om_k = 1 - m4.Omega_cb_0 - m4.Omega_r_massless_0 - m4.Omega_nu_interp_z0 - om_de;
    var pp = { om_r: p.om_r, om_nu: p.om_nu, om_nu_r1: p.om_nu_r1, om_cb: p.om_cb, om_k: om_k, om_de: om_de };
    var E2 = function (z) { return E2general(z, pp); }, Hz = function (z) { return H0_STEP25 * Math.sqrt(E2(z)); };
    var bg = { H0: H0_STEP25, Omega_b_0: m4.Omega_b_0, Omega_gamma_0: m4.Omega_gamma_0, E2: E2, Hz: Hz, Omega_k_0: om_k };
    var rs = soundHorizon(Z_STAR, bg), rd = soundHorizon(Z_DRAG, bg), dm = DMofZ(Z_STAR, Hz, om_k, H0_STEP25);
    // M2 flat uncompensated: E^2(0) and root
    var pp2 = { om_r: p.om_r, om_nu: p.om_nu, om_nu_r1: p.om_nu_r1, om_cb: p.om_cb, om_k: 0, om_de: om_de };
    var E2m2 = function (z) { return E2general(z, pp2); }, lo = 0, hi = 1;
    for (var i = 0; i < 200; i++) { var mid = 0.5 * (lo + hi); if (E2m2(lo) * E2m2(mid) <= 0) hi = mid; else lo = mid; }
    var th = 100 * rs / dm;
    return { om_k: om_k, om_de: om_de, rs: rs, rd: rd, dm_zstar: dm, theta100: th, pullTheta: (th - THETA100_PL) / THETA100_PL_ERR,
      chi2_bao: T.chi2BAO(Hz, rd, om_k, H0_STEP25), Hz: Hz, bg: bg, E2_m2_z0: E2m2(0), z_root_m2: 0.5 * (lo + hi) };
  };
  // Planck 3D illustrative prior on (omega_b, omega_cb, 100 theta_*)
  var PL3_MEAN = [0.02237, 0.14240, 1.04109], PL3_SIG = [0.00015, 0.00121, 0.00030], PL3_CORR = [[1, -0.62, 0.45], [-0.62, 1, -0.41], [0.45, -0.41, 1]];
  var PL3_COV = PL3_CORR.map(function (r, i) { return r.map(function (c, j) { return c * PL3_SIG[i] * PL3_SIG[j]; }); });
  T.PL3_INV = invert(PL3_COV);
  T.planck3dChi2 = function (obh2, ocbh2, th100) { var d = [obh2 - PL3_MEAN[0], ocbh2 - PL3_MEAN[1], th100 - PL3_MEAN[2]]; return quadForm(d, T.PL3_INV, d); };
  T.cmbForBackground = function (bg, Hz) {
    Hz = Hz || bg.Hz;
    var rs = soundHorizon(Z_STAR, bg), rd = soundHorizon(Z_DRAG, bg), dm = DMofZ(Z_STAR, Hz, bg.Omega_k_0 || 0, bg.H0);
    var th = 100 * rs / dm;
    return { rs: rs, rd: rd, dm: dm, theta100: th, pull: (th - THETA100_PL) / THETA100_PL_ERR,
      chi2_3d: T.planck3dChi2(bg.omega_b_h2, bg.omega_cb_h2, th) };
  };

  // ---------- SNe Ia synthetic distance-modulus shape control ----------
  T.muOfZ = function (z, Hz, om_k, H0c) { return 5 * Math.log10((1 + z) * DMofZ(z, Hz, om_k, H0c)) + 25; };
  T.sneControl = function (bgModel, bgRef) {
    bgRef = bgRef || makeBackground(H0_LOCAL, 0.334);
    var g = D.sne_grid, wsum = 0, wd = 0, rows = [];
    g.forEach(function (n) {
      var mr = T.muOfZ(n[0], bgRef.Hz, 0, bgRef.H0), mm = T.muOfZ(n[0], bgModel.Hz, bgModel.Omega_k_0 || 0, bgModel.H0);
      var w = 1 / (n[1] * n[1]); wsum += w; wd += w * (mm - mr); rows.push({ z: n[0], sig: n[1], muRef: mr, muModel: mm, d: mm - mr });
    });
    var off = wd / wsum, c0 = 0, c1 = 0;
    rows.forEach(function (r) { c0 += Math.pow(r.d / r.sig, 2); c1 += Math.pow((r.d - off) / r.sig, 2); });
    return { rows: rows, offset: off, chi2_unshifted: c0, chi2_shifted: c1 };
  };

  // ---------- RSD growth with variational Yukawa perturbation ----------
  T.growth = function (bg, opts) {
    opts = opts || {}; var kh = opts.k_h || 0.1, cbOnly = opts.cbOnly !== false, alphaY = opts.alphaY === undefined ? 1 : opts.alphaY;
    var L_m = (opts.L_km === undefined ? 10 : opts.L_km) * 1000, sig8 = opts.sigma8 || 0.8111, N = opts.steps || 6000;
    var om0 = cbOnly ? bg.Omega_cb_0 : bg.Omega_m_0, kL0 = (kh * bg.h / MPC_M) * L_m, q = alphaY * alphaY * T.BETA;
    function coef(x) {
      var a = Math.exp(x), z = 1 / a - 1, e2 = bg.E2(z), eps = 1e-6 * (1 + z);
      var de2 = (bg.E2(z + eps) - bg.E2(z - eps)) / (2 * eps);
      var kk = kL0 / a, k2 = kk * kk;
      return [-0.5 * (1 + z) * de2 / e2, om0 * Math.pow(a, -3) / e2, q * (k2 / (1 + k2))];
    }
    function rhs(x, y) { var c = coef(x), d = 2 + c[0]; return [y[1], -d * y[1] + 1.5 * c[1] * y[0], y[3], -d * y[3] + 1.5 * c[1] * y[2] + 1.5 * c[1] * c[2] * y[0]]; }
    var x0 = Math.log(1e-3), h = -x0 / N, y = [1e-3, 1e-3, 0, 0], xs = [x0], ys = [y.slice()];
    for (var i = 0; i < N; i++) {
      var x = x0 + i * h, k1 = rhs(x, y);
      var y2 = y.map(function (v, j) { return v + 0.5 * h * k1[j]; }), k2 = rhs(x + 0.5 * h, y2);
      var y3 = y.map(function (v, j) { return v + 0.5 * h * k2[j]; }), k3 = rhs(x + 0.5 * h, y3);
      var y4 = y.map(function (v, j) { return v + h * k3[j]; }), k4 = rhs(x + h, y4);
      y = y.map(function (v, j) { return v + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]); });
      xs.push(x + h); ys.push(y.slice());
    }
    var D0 = y[0], delD0 = y[2];
    function at(z) {
      var x = Math.log(1 / (1 + z)), t = (x - x0) / h, i = Math.min(N - 1, Math.max(0, Math.floor(t))), f = t - i;
      // cubic Hermite interpolation using the ODE derivatives
      var ya = ys[i], yb = ys[i + 1], da = rhs(xs[i], ya), db = rhs(xs[i + 1], yb), out = [];
      for (var j = 0; j < 4; j++) {
        var h00 = 2 * f * f * f - 3 * f * f + 1, h10 = f * f * f - 2 * f * f + f, h01 = -2 * f * f * f + 3 * f * f, h11 = f * f * f - f * f;
        out.push(h00 * ya[j] + h10 * h * da[j] + h01 * yb[j] + h11 * h * db[j]);
      }
      var Dn = out[0] / D0, f_g = out[1] / out[0], fs8 = sig8 * f_g * Dn;
      var dfs8 = sig8 * (out[3] - out[1] * (delD0 / D0)) / D0;
      var kk = kL0 * (1 + z), k2 = kk * kk;
      return { D: Dn, f: f_g, fs8: fs8, dfs8: Math.abs(dfs8), dGeff: q * k2 / (1 + k2) };
    }
    return { at: at, kL0: kL0 };
  };
  T.rsdChi2 = function (gr, key) {
    var c = 0, rows = [], maxd = 0;
    D.rsd.forEach(function (r) { var g = gr.at(r[1]); c += Math.pow((g.fs8 - r[2]) / r[3], 2); maxd = Math.max(maxd, g.dfs8); rows.push({ name: r[0], z: r[1], obs: r[2], err: r[3], fs8: g.fs8, dfs8: g.dfs8 }); });
    return { chi2: c, rows: rows, maxDelta: maxd };
  };

  // ---------- algebraic BBN benchmark ----------
  T.bbn = function (z, L_km, alphaY) {
    var bg = T.BG_M4, zp1 = 1 + z, omr = bg.Omega_gamma_0 * (1 + 0.22710731766 * N_EFF);
    var E2 = omr * Math.pow(zp1, 4) + OMEGA_M_0 * Math.pow(zp1, 3), H = H0_S_INV * Math.sqrt(E2), mphi = C_SI / (L_km * 1000);
    var rho = OMEGA_M_0 * Math.pow(zp1, 3) * T.RHO_CRIT_0, phi = alphaY * T.gammaL(L_km) * rho;
    var dH = 0.5 * (mphi * mphi * phi * phi) / (3 * T.BETA * H * H), dG = 2 * phi;
    var Yp = 0.24709 + 0.08 * (dH + 0.5 * dG), DH = 2.5147 + 2.5147 * 1.6 * (dH + 0.5 * dG);
    return { H: H, H_over_mphi: H / mphi, rho: rho, absPhi: phi, dHH: dH, dGG: dG, Yp: Yp, DH: DH,
      pullYp: (Yp - 0.245) / 0.003, pullDH: (DH - 2.527) / 0.030 };
  };
  // local constant-density bulk benchmark: |phi| = alpha_Y gamma_L rho, apparent H0 = H0 e^{|phi|}
  T.localBulk = function (rho, L_km, alphaY) {
    var phi = alphaY * T.gammaL(L_km) * rho, H0a = H0_STEP25 * Math.exp(phi);
    return { absPhi: phi, fClock: Math.exp(-phi), H0app: H0a, dH0: H0_STEP25 * Math.expm1(phi), frac: H0_STEP25 * Math.expm1(phi) / (H0_LOCAL - H0_STEP25) };
  };

  // ---------- Yukawa force law and G dictionary ----------
  // Psi = -G_* M/r (1 + q e^{-r/L}),  q = alpha_Y^2 beta;  G_eff(r)/G_* = 1 + q (1 + r/L) e^{-r/L}
  T.geffRatio = function (r, L, alphaY) { var q = T.qCoupling(alphaY), x = r / L; return 1 + q * (1 + x) * Math.exp(-x); };
  T.potRatio = function (r, L, alphaY) { return 1 + T.qCoupling(alphaY) * Math.exp(-r / L); };
  T.phiPoint = function (M_kg, r_m, L_m, alphaY) { return -alphaY * T.BETA * (G_SI * M_kg / (C_SI * C_SI * r_m)) * Math.exp(-r_m / L_m); };
  T.gDictionary = function (alphaY, phi) {
    var q = T.qCoupling(alphaY), A = Math.exp(alphaY * phi);
    return { q: q, A: A, GT_over_Gstar: A * A, Glab_over_Gstar: 1 + q, Geff_over_Glab_large: 1 / (1 + q) };
  };
  // PPN slip and lensing/dynamics ratio (V26 Step 23 closed forms)
  T.gammaPPN = function (r, L, beta) { var x = r / L; if (x > 700) return 1 - 2 * beta * Math.exp(-x); var e = Math.exp(x); return (e - beta) / (e + beta); };
  T.log10OneMinusGamma = function (r, L, beta) { var x = r / L; if (x > 30) return Math.log10(2 * beta) - x / Math.LN10; return Math.log10(Math.abs(1 - T.gammaPPN(r, L, beta))); };
  T.lensDynRatio = function (r, L, beta) { var x = r / L; if (x > 700) return 1; var e = Math.exp(x); return e / (beta + e); };
  // exact vacuum Schwarzschild null-geodesic deflection (128-pt GL in u = sin(theta)), port of V26 Step 23
  T.deflectionSchw = function (b_km, M_km) {
    var nu = function (r) { return 0.5 * Math.log1p(-2 * M_km / r); };
    var r0;
    if (b_km <= 1e5) {
      var lo = Math.max(3 * M_km * 1.0000001, 1e-9), hi = b_km * 1.05, g = function (r) { return r * Math.exp(-nu(r)) - b_km; };
      for (var it = 0; it < 300; it++) { var mid = 0.5 * (lo + hi); if (g(lo) * g(mid) <= 0) hi = mid; else lo = mid; }
      r0 = 0.5 * (lo + hi);
    } else {
      var x = M_km / b_km; r0 = b_km * (1 - x - 1.5 * x * x - 4 * x * x * x);
      for (var k = 0; k < 3; k++) { var fv = r0 / Math.sqrt(1 - 2 * M_km / r0) - b_km, df = (1 - 3 * M_km / r0) / Math.pow(1 - 2 * M_km / r0, 1.5); r0 -= fv / df; }
    }
    var total = 0, hp = 0.5 * Math.PI;
    for (var i = 0; i < 128; i++) {
      var th = 0.5 * hp * (GL128.x[i] + 1), wt = 0.5 * hp * GL128.w[i], u = Math.sin(th), c = Math.cos(th), r = r0 / u;
      var Bm1 = -2 * M_km / r, B = 1 + Bm1, omu = c * c / (1 + u);
      var D0 = Math.log1p(-2 * M_km * omu / (r0 * (1 - 2 * M_km * u / r0)));
      var C = Math.expm1(D0) / (c * c), Y = Bm1 + B * C, opY = Math.max(1 + Y, 1e-30);
      total += wt * (-Y / (Math.sqrt(opY) * (1 + Math.sqrt(opY))));
    }
    return { alpha: 2 * total, r0: r0, alpha1PN: 4 * M_km / b_km };
  };
  T.shapiroSchw = function (r0, rEnd, M_km) {
    var wmax = Math.acosh(rEnd / r0), tot = 0;
    for (var i = 0; i < 128; i++) {
      var w = 0.5 * wmax * (GL128.x[i] + 1), wt = 0.5 * wmax * GL128.w[i], ch = Math.cosh(w), sh = Math.sinh(w), r = r0 * ch;
      var n = 0.5 * Math.log1p(-2 * M_km / r), A = Math.expm1(2 * n) - (2 * M_km / r) * Math.exp(2 * n);
      var Dl = Math.log1p((2 * M_km * sh * sh) / (r0 * ch * (ch + 1) * (1 - 2 * M_km / r0)));
      var Cf = -Math.expm1(Dl) / (sh * sh), X = A + Cf + A * Cf, opX = Math.max(1 + X, 1e-30);
      tot += wt * (r0 * ch * (-X / (Math.sqrt(opX) * (1 + Math.sqrt(opX)))));
    }
    return 2 * tot / C_KM_S * 1e6;
  };
  T.cassiniDelay = function (b_km, M_km) { var tM = M_km / C_KM_S * 1e6, ln = Math.log(4 * T.AU_KM * 8.43 * T.AU_KM / (b_km * b_km)); return { rt: 4 * tM * ln, rtm1: 4 * tM * (ln - 1) }; };
  // GW170817: tensor c_T = c; massive scalar dispersion omega^2 = c^2 k^2 + c^2/L^2
  T.gwScalar = function (f_hz, L_km, D_mpc) {
    var fcut = (C_KM_S / L_km) / (2 * Math.PI), ratio = f_hz / fcut, Dkm = D_mpc * T.MPC_KM, tl = Dkm / C_KM_S;
    if (ratio < 1) { var kap = (1 / L_km) * Math.sqrt(1 - ratio * ratio); return { fcut: fcut, ratio: ratio, regime: "evanescent", kappa: kap, lambda: 1 / kap, log10amp: -(Dkm * kap) / Math.LN10, vg: 0, delay: 0, cT: 1 }; }
    var vg = Math.sqrt(1 - 1 / (ratio * ratio));
    return { fcut: fcut, ratio: ratio, regime: "propagating", kappa: 0, vg: vg, vp: 1 / vg, log10amp: 0, delay: vg > 0 ? tl * (1 / vg - 1) : Infinity, cT: 1 };
  };

  // ---------- SPARC rotation curves ----------
  T.galaxyModel = function (g, Yd, Yb, L_kpc, alphaY) {
    var out = [], c2b = 0, c2m = 0, n = 0;
    g.d.forEach(function (p) {
      var R = p[0], V = p[1], e = p[2], vg = p[3], vd = p[4], vb = p[5];
      var vb2 = vg * Math.abs(vg) + Yd * vd * Math.abs(vd) + Yb * vb * Math.abs(vb);
      var vbar = Math.sign(vb2) * Math.sqrt(Math.abs(vb2));
      var boost = T.geffRatio(R, L_kpc, alphaY); // enclosed-mass approximation for a Yukawa term of range L
      var vm2 = vb2 * boost, vm = Math.sign(vm2) * Math.sqrt(Math.abs(vm2));
      c2b += Math.pow((V - vbar) / e, 2); c2m += Math.pow((V - vm) / e, 2); n++;
      out.push({ R: R, V: V, e: e, vbar: vbar, vmod: vm, Creq: (V * V) / vb2 });
    });
    return { pts: out, chi2_bar: c2b, chi2_mod: c2m, n: n };
  };
  T.bestUpsilon = function (g) {
    var r = brent(function (Y) { return T.galaxyModel(g, Y, 1.4 * Y, 1e-15, 0).chi2_bar; }, 0.05, 0.5, 2.0, 1e-8);
    return r.x;
  };

  // chi^2 survival function Q(k/2, x/2) (series / Lentz continued fraction)
  function lgamma(x) { var c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5], y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t); var s = 1.000000000190015; for (var j = 0; j < 6; j++) s += c[j] / ++y; return -t + Math.log(2.5066282746310005 * s / x); }
  T.chi2sf = function (x, k) {
    var a = k / 2, xx = x / 2; if (xx <= 0) return 1; var gl = lgamma(a);
    if (xx < a + 1) { var ap = a, sum = 1 / a, del = sum; for (var n = 0; n < 500; n++) { ap++; del *= xx / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-15) break; } return 1 - sum * Math.exp(-xx + a * Math.log(xx) - gl); }
    var b = xx + 1 - a, cc = 1e300, d = 1 / b, h = d;
    for (var i = 1; i < 500; i++) { var an = -i * (i - a); b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300; cc = b + an / cc; if (Math.abs(cc) < 1e-300) cc = 1e-300; d = 1 / d; var dl = d * cc; h *= dl; if (Math.abs(dl - 1) < 1e-15) break; }
    return Math.exp(-xx + a * Math.log(xx) - gl) * h;
  };
  T.normSigma = function (p) { // two-sided Gaussian equivalent of a p-value (bisection on erfc)
    if (!(p > 0)) return Infinity; var lo = 0, hi = 40;
    function erfc(x) { var t = 1 / (1 + 0.5 * x), y = t * Math.exp(-x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))); return y; }
    for (var i = 0; i < 200; i++) { var m = 0.5 * (lo + hi); if (erfc(m / Math.SQRT2) > p) lo = m; else hi = m; } return 0.5 * (lo + hi);
  };

  root.TAL = T;
})(typeof window !== "undefined" ? window : globalThis);
