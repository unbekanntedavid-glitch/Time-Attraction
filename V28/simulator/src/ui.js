/* Time Attraction Law V27 simulator: user interface. */
(function () {
  "use strict";
  var T = window.TAL, D = window.TAL_DATA, R = D.ref, I = window.TAL_I18N, P = window.TALPlot;
  var C = { acc: "#5cc8ff", acc2: "#ffb454", ok: "#5bd69a", vio: "#b58cff", bad: "#ff6b6b", txt: "#e8eefb", mute: "#8fa1c2" };
  var lang = "hu";
  try { var sv = window.localStorage.getItem("tal_lang"); if (sv === "en" || sv === "hu") lang = sv; } catch (e) { /* storage unavailable */ }
  function t(k) { return (I[lang] && I[lang][k] !== undefined) ? I[lang][k] : k; }
  function plain(html) { var d = document.createElement("div"); d.innerHTML = html; return d.textContent; }
  function $(id) { return document.getElementById(id); }
  function fmt(v, d) {
    d = d || 6; if (v === null || v === undefined || typeof v !== "number") return String(v);
    if (!isFinite(v)) return v > 0 ? "∞" : (v < 0 ? "−∞" : "NaN");
    var a = Math.abs(v); var s;
    if (a !== 0 && (a < 1e-4 || a >= 1e7)) { s = v.toExponential(d - 1).replace("e", "·10^"); s = s.replace(/\^\+?(-?)0*(\d+)/, function (m, sg, n) { return P.sup((sg ? "-" : "") + n); }).replace("^", ""); }
    else s = String(+v.toPrecision(d));
    return s.replace(/^-/, "−");
  }
  // ---------------- widgets ----------------
  var sliders = {};
  function slider(container, id, key, o, cb) {
    var wrap = document.createElement("div"); wrap.className = "sl";
    var lab = document.createElement("label"); lab.setAttribute("data-t", key); lab.innerHTML = t(key);
    var rg = document.createElement("input"); rg.type = "range"; rg.id = id;
    var nb = document.createElement("input"); nb.type = "number"; nb.step = "any";
    var s = { o: o, log: !!o.log, wrap: wrap, rg: rg, nb: nb, v: o.value };
    function setRange(o2) {
      s.o = o2; s.log = !!o2.log;
      rg.min = s.log ? Math.log10(o2.min) : o2.min; rg.max = s.log ? Math.log10(o2.max) : o2.max;
      rg.step = s.log ? (o2.lstep || 0.001) : (o2.step || "any");
    }
    setRange(o);
    s.setRange = setRange;
    s.set = function (v, silent) { s.v = v; rg.value = s.log ? Math.log10(v) : v; nb.value = +(+v).toPrecision(o.prec || 7); if (!silent && cb) cb(); };
    s.get = function () { return s.v; };
    s.disable = function (b) { rg.disabled = b; nb.disabled = b; wrap.classList.toggle("dis", b); };
    rg.addEventListener("input", function () { var x = +rg.value; s.v = s.log ? Math.pow(10, x) : x; nb.value = +s.v.toPrecision(o.prec || 7); if (cb) cb(); });
    nb.addEventListener("change", function () { var v = +nb.value; if (isFinite(v) && (!s.log || v > 0)) { s.v = v; rg.value = s.log ? Math.log10(v) : v; if (cb) cb(); } });
    wrap.appendChild(lab); wrap.appendChild(rg); wrap.appendChild(nb); container.appendChild(wrap);
    s.set(o.value, true); sliders[id] = s; return s;
  }
  function kv(el, rows) {
    el.innerHTML = rows.map(function (r) {
      if (r.h) return '<div class="hd">' + t(r.h) + "</div>";
      return "<div><span>" + t(r[0]) + "</span><b" + (r[2] ? ' class="' + r[2] + '"' : "") + ">" + r[1] + "</b></div>";
    }).join("");
  }
  function logspace(a, b, n) { var out = [], la = Math.log10(a), lb = Math.log10(b); for (var i = 0; i < n; i++) out.push(Math.pow(10, la + (lb - la) * i / (n - 1))); return out; }
  function linspace(a, b, n) { var out = []; for (var i = 0; i < n; i++) out.push(a + (b - a) * i / (n - 1)); return out; }

  // ---------------- Tab 1: Yukawa ----------------
  var BODIES = {
    sun: { M: T.MSUN_KG, R: T.R_SUN_KM, key: "body_sun" },
    earth: { M: 5.9722e24, R: 6371.0, key: "body_earth" },
    ns: { M: R.ns_ref.M_ADM_Msun * T.MSUN_KG, R: R.ns_ref.R_star_km, key: "body_ns" }
  };
  var yk = {};
  function initYukawa() {
    var c = $("yk-ctrl");
    yk.L = slider(c, "yk-L", "L_km", { min: 0.1, max: 1e5, value: 10, log: true }, updYukawa);
    yk.a = slider(c, "yk-a", "alphaY", { min: 0, max: 3, step: 0.001, value: 1 }, updYukawa);
    yk.x = slider(c, "yk-x", "r_over_L", { min: 0.01, max: 100, value: 1, log: true }, updYukawa);
    var row = document.createElement("label"); row.className = "row";
    row.innerHTML = '<span data-t="body">' + t("body") + '</span><select id="yk-body"><option value="ns" data-t="body_ns"></option><option value="sun" data-t="body_sun"></option><option value="earth" data-t="body_earth"></option></select>';
    c.appendChild(row); $("yk-body").addEventListener("change", updYukawa);
    var pn = document.createElement("p"); pn.className = "ptnote"; pn.id = "yk-point-note"; pn.setAttribute("data-t", "yk_point_note"); pn.innerHTML = t("yk_point_note"); c.appendChild(pn);
  }
  function updYukawa() {
    var L = yk.L.get(), a = yk.a.get(), x = yk.x.get(), r = x * L, b = BODIES[$("yk-body").value];
    var q = T.qCoupling(a), ge = T.geffRatio(r, L, a), pr = T.potRatio(r, L, a);
    var rr = Math.max(r, b.R), phi = T.phiPoint(b.M, rr * 1000, L * 1000, a), gd = T.gDictionary(a, phi);
    kv($("yk-out"), [
      ["o_beta", fmt(T.BETA, 8)], ["o_q", fmt(q, 8)], ["o_glab", fmt(1 + q, 9)], ["o_ginf", fmt(gd.Geff_over_Glab_large, 9), Math.abs(a - 1) < 1e-12 ? "ok" : ""],
      ["o_r", fmt(r, 6)], ["o_geffr", fmt(ge, 9)], ["o_geffr_lab", fmt(ge / (1 + q), 9)], ["o_pot", fmt(pr, 9)],
      ["o_phi", (phi !== 0 || a === 0 ? fmt(phi, 6) : "−10^" + (Math.log10(a * T.BETA * T.G_SI * b.M / (T.C_SI * T.C_SI * rr * 1000)) - (rr / L) / Math.LN10).toFixed(2)) + (r < b.R ? " (r = R)" : "")],
      ["o_A", fmt(gd.A, 12)], ["o_GT", fmt(gd.GT_over_Gstar, 12)],
      ["o_gammaL", fmt(T.gammaL(L), 7)], ["o_mphi", fmt(T.HBAR_C_EV_M / (L * 1000), 7)], ["o_fcut", fmt(T.C_KM_S / L / (2 * Math.PI), 7)]
    ]);
    var rs = logspace(0.01 * L, 100 * L, 320);
    P.plot($("yk-chart"), { xlog: true, xlabel: "r [km]", ylabel: "G_eff / G", series: [
      { label: "G_eff(r)/G*", data: rs.map(function (v) { return [v, T.geffRatio(v, L, a)]; }), color: C.acc },
      { label: "G_eff(r)/G_lab", data: rs.map(function (v) { return [v, T.geffRatio(v, L, a) / (1 + q)]; }), color: C.acc2 },
      { label: "Ψ/Ψ_N", data: rs.map(function (v) { return [v, T.potRatio(v, L, a)]; }), color: C.vio, dash: [6, 4], width: 1.6 },
      { data: [[r, ge]], type: "points", color: C.txt, size: 4.5 }
    ], vlines: [{ x: L, label: "L", color: C.mute }], hlines: [{ y: 1, color: C.mute }, { y: 1 / (1 + q), color: C.acc2, dash: [2, 4] }] });
    var rb = logspace(Math.max(b.R, 0.01 * L) , Math.max(b.R, L) * 200, 260), K = T.BETA * T.G_SI * b.M / (T.C_SI * T.C_SI);
    var ys = a > 0 ? rb.map(function (v) { return [v, Math.log10(2 * a * a * K / (v * 1000)) - (v / L) / Math.LN10]; }) : [];
    P.plot($("yk-chart2"), { xlog: true, height: 230, xlabel: "r [km]", ylabel: "log10 |A² − 1|", series: [{ label: plain(t(b.key)), data: ys, color: C.ok }],
      vlines: [{ x: L, label: "L", color: C.mute }] });
  }

  // ---------------- Tab 2: SPARC ----------------
  var sp = {}, LKPC_PKG = 10 / T.KPC_KM;
  function initSparc() {
    var sel = $("sp-gal"); D.sparc.forEach(function (g, i) { var o = document.createElement("option"); o.value = i; o.textContent = g.n; sel.appendChild(o); });
    sel.value = D.sparc.findIndex(function (g) { return g.n === "NGC3198"; });
    sel.addEventListener("change", updSparc);
    sp.Yd = slider($("sp-ctrl"), "sp-yd", "Yd", { min: 0.05, max: 2, step: 0.01, value: 0.5 }, updSparc);
    sp.Yb = slider($("sp-ctrl"), "sp-yb", "Yb", { min: 0.05, max: 2, step: 0.01, value: 0.5 }, updSparc);
    sp.L = slider($("sp-ctrl2"), "sp-L", "Lgal", { min: 0.1, max: 300, value: 10, log: true }, updSparc);
    sp.a = slider($("sp-ctrl2"), "sp-a", "alphaYgal", { min: 0, max: 3, step: 0.01, value: 1 }, updSparc);
    $("sp-explore").addEventListener("change", updSparc);
    $("sp-best").addEventListener("click", function () { var g = D.sparc[+$("sp-gal").value], y = T.bestUpsilon(g); sp.Yd.set(y, true); sp.Yb.set(1.4 * y, true); updSparc(); });
    $("sp-reset").addEventListener("click", function () { sp.Yd.set(0.5, true); sp.Yb.set(0.5, true); $("sp-explore").checked = false; updSparc(); });
  }
  function median(a) { var b = a.slice().sort(function (x, y) { return x - y; }), n = b.length; return n % 2 ? b[(n - 1) / 2] : 0.5 * (b[n / 2 - 1] + b[n / 2]); }
  function updSparc() {
    var ex = $("sp-explore").checked; sp.L.disable(!ex); sp.a.disable(!ex);
    var L = ex ? sp.L.get() : LKPC_PKG, a = ex ? sp.a.get() : 1, Yd = sp.Yd.get(), Yb = sp.Yb.get();
    var g = D.sparc[+$("sp-gal").value], m = T.galaxyModel(g, Yd, Yb, L, a), last = m.pts[m.pts.length - 1];
    kv($("sp-out"), [["s_n", m.n], ["s_D", fmt(g.D, 4)], ["s_rlast", fmt(last.R, 5)], ["s_c2b", fmt(m.chi2_bar, 7)], ["s_c2m", fmt(m.chi2_mod, 7)],
      ["s_c2n", fmt(m.chi2_mod / m.n, 5)], ["s_ratio", fmt(m.chi2_mod / m.chi2_bar, 7)], ["s_creq", fmt(last.Creq, 6)]]);
    var sq = function (Y, v) { return Math.sign(v) * Math.sqrt(Y) * Math.abs(v); };
    var ser = [
      { label: "V_obs", data: g.d.map(function (p) { return [p[0], p[1], p[2]]; }), type: "points", color: C.txt },
      { label: "V_gas", data: g.d.map(function (p) { return [p[0], p[3]]; }), color: C.ok, dash: [4, 3], width: 1.4 },
      { label: "√Υd·V_disk", data: g.d.map(function (p) { return [p[0], sq(Yd, p[4])]; }), color: C.vio, dash: [4, 3], width: 1.4 }
    ];
    if (g.d.some(function (p) { return p[5] !== 0; })) ser.push({ label: "√Υb·V_bul", data: g.d.map(function (p) { return [p[0], sq(Yb, p[5])]; }), color: C.bad, dash: [4, 3], width: 1.4 });
    ser.push({ label: "V_bar", data: m.pts.map(function (p) { return [p.R, p.vbar]; }), color: C.acc2, width: 2.4 });
    ser.push({ label: ex ? "V_mod (L_gal)" : "V_mod (L = 10 km)", data: m.pts.map(function (p) { return [p.R, p.vmod]; }), color: C.acc, dash: ex ? [] : [7, 5], width: 2 });
    P.plot($("sp-chart"), { xlabel: "R [kpc]", ylabel: "V [km/s]", title: g.n, series: ser, legend: "br" });
    var s2 = [{ label: "C_req = (V_obs/V_bar)²", data: m.pts.map(function (p) { return [p.R, p.Creq]; }), type: "points", color: C.acc2 }];
    if (ex) s2.push({ label: "G_eff(R)/G*", data: m.pts.map(function (p) { return [p.R, T.geffRatio(p.R, L, a)]; }), color: C.acc });
    P.plot($("sp-chart2"), { height: 230, xlabel: "R [kpc]", ylabel: "C", series: s2, hlines: [{ y: 1, color: C.mute }], legend: "tl" });
    var cb = [], cm = [], cr = [], sb = 0, sm = 0, np = 0;
    D.sparc.forEach(function (gg) { var mm = T.galaxyModel(gg, Yd, Yb, L, a); cb.push(mm.chi2_bar / mm.n); cm.push(mm.chi2_mod / mm.n); var lp = mm.pts[mm.pts.length - 1]; if (isFinite(lp.Creq)) cr.push(lp.Creq); sb += mm.chi2_bar; sm += mm.chi2_mod; np += mm.n; });
    kv($("sp-sample"), [["s_npts", np], ["s_med_c2n", fmt(median(cb), 5)], ["s_med_c2n_m", fmt(median(cm), 5)], ["s_med_creq", fmt(median(cr), 5)], ["s_sum", fmt(sb, 6) + " / " + fmt(sm, 6)]]);
  }

  // ---------------- Tab 3: GW & lensing ----------------
  var gw = {}, ln = {};
  function initGW() {
    gw.f = slider($("gw-ctrl"), "gw-f", "f_hz", { min: 1e-5, max: 1e5, value: 100, log: true }, updGW);
    gw.D = slider($("gw-ctrl"), "gw-D", "D_mpc", { min: 0.01, max: 3000, value: 40, log: true }, updGW);
    gw.L = slider($("gw-ctrl"), "gw-L", "L_km", { min: 0.1, max: 1e5, value: 10, log: true }, updGW);
    $("ln-body").addEventListener("change", function () { lensRange(true); updLens(); });
    ln.b = slider($("ln-ctrl"), "ln-b", "l_b", { min: 14, max: 1e4, value: 14.749, log: true }, updLens);
    ln.L = slider($("ln-ctrl"), "ln-L", "L_km", { min: 0.1, max: 1e5, value: 10, log: true }, updLens);
    $("ln-body").value = "ns"; lensRange(true);
  }
  function lensBody() {
    var ns = $("ln-body").value === "ns", M = ns ? R.ns_ref.M_ADM_Msun * T.GM_SUN_KM : T.GM_SUN_KM, Rb = ns ? R.ns_ref.R_star_km : T.R_SUN_KM;
    return { ns: ns, M: M, R: Rb, bmin: ns ? Rb / Math.sqrt(1 - 2 * M / Rb) : Rb };
  }
  function lensRange(reset) { var b = lensBody(); ln.b.setRange({ min: b.bmin, max: b.ns ? 1e4 : 1000 * b.R, log: true }); if (reset) ln.b.set(b.bmin, true); }
  function updGW() {
    var f = gw.f.get(), Dm = gw.D.get(), L = gw.L.get(), g = T.gwScalar(f, L, Dm), m = R.gw_meta;
    var inb = (0 >= m.GW170817_cT_lower_bound && 0 <= m.GW170817_cT_upper_bound);
    kv($("gw-out"), [["g_cT", "1"], ["g_dcT", "0"], ["g_bound", "[−3·10⁻¹⁵, +7·10⁻¹⁵] ✓ " + t("ok"), inb ? "ok" : "bad"], ["g_dt", "1.74"],
      ["g_fcut", fmt(g.fcut, 8)], ["g_mphi", fmt(T.HBAR_C_EV_M / (L * 1000), 6)], ["g_regime", t(g.regime === "evanescent" ? "reg_ev" : "reg_pr")],
      ["g_kappa", fmt(g.kappa, 6)], ["g_amp", fmt(g.log10amp, 8)], ["g_vg", fmt(g.vg, 10)], ["g_delay", fmt(g.delay, 6)]]);
    var fs = logspace(1e-5, 1e5, 600), fc = g.fcut;
    if (fc > 1e-5 && fc < 1e5) fs = fs.concat(logspace(fc * 1.0000001, fc * 1.5, 120)).sort(function (x, y) { return x - y; });
    P.plot($("gw-chart"), { xlog: true, xlabel: "f [Hz]", ylabel: "v / c", ylim: [-0.03, 1.08], series: [
      { label: "c_T / c (tensor)", data: [[1e-5, 1], [1e5, 1]], color: C.acc2, dash: [6, 4] },
      { label: "v_g / c (scalar)", data: fs.map(function (v) { return [v, T.gwScalar(v, L, Dm).vg]; }), color: C.acc },
      { data: [[f, g.vg]], type: "points", color: C.txt, size: 4.5 }
    ], vbands: [{ x0: 24, x1: 600, color: "rgba(91,214,154,0.12)" }], vlines: [{ x: fc, label: "f_cut", color: C.vio }], legend: "tl" });
  }
  function updLens() {
    var b = lensBody(), bk = Math.max(ln.b.get(), b.bmin), L = ln.L.get(), d = T.deflectionSchw(bk, b.M), AS = 180 / Math.PI * 3600;
    var g = T.gammaPPN(bk, L, T.BETA), lg = T.log10OneMinusGamma(bk, L, T.BETA), sh = T.shapiroSchw(d.r0, 10 * d.r0, b.M);
    var rows = [["l_b", fmt(bk, 9)], ["l_r0", fmt(d.r0, 9)],
      ["l_alpha", b.ns ? fmt(d.alpha, 9) + " rad" : fmt(d.alpha * AS, 10) + "″"], ["l_alpha1", b.ns ? fmt(d.alpha1PN, 9) + " rad" : fmt(d.alpha1PN * AS, 10) + "″"],
      ["l_gamma", fmt(g, 12)], ["l_log", fmt(lg, 8)], ["l_ratio", fmt(T.lensDynRatio(bk, L, T.BETA), 12)], ["l_shap", fmt(sh, 10)]];
    if (!b.ns) rows.push(["l_cass", fmt(T.cassiniDelay(bk, b.M).rt, 10)]);
    if (!b.ns) rows.push(["l_cassok", lg < Math.log10(2.3e-5) ? "✓" : "✗", lg < Math.log10(2.3e-5) ? "ok" : "bad"]);
    kv($("ln-out"), rows);
    var xmax = Math.max(1e3, 3 * bk / L), xs = logspace(0.1, xmax, 400);
    P.plot($("ln-chart"), { xlog: true, xlabel: "b / L", ylabel: "log10 |1 − γ_PPN|", series: [
      { label: "log10|1 − γ|", data: xs.map(function (x) { return [x, T.log10OneMinusGamma(x * L, L, T.BETA)]; }), color: C.acc },
      { data: [[bk / L, lg]], type: "points", color: C.acc2, size: 5 }
    ], hlines: [{ y: Math.log10(2.3e-5), color: C.bad }], legend: "tr" });
  }

  // ---------------- Tab 4: Cosmology ----------------
  var co = { fit: null, joint: null, ctl: null, m3: null };
  function initCosmo() {
    var c = $("co-ctrl");
    co.H0 = slider(c, "co-H0", "H0", { min: 55, max: 85, step: 0.01, value: 72.8 }, updCosmo);
    co.Om = slider(c, "co-Om", "Om", { min: 0.15, max: 0.5, step: 0.0001, value: T.OMEGA_M_0, prec: 9 }, updCosmo);
    co.rd = slider(c, "co-rd", "rd", { min: 120, max: 165, step: 0.01, value: T.RD_ZDRAG, prec: 9 }, updCosmo);
    co.P = slider(c, "co-P", "P", { min: 8000, max: 12500, step: 0.1, value: 72.8 * T.RD_ZDRAG, prec: 9 }, updCosmo);
    $("co-bg").addEventListener("change", updCosmo); $("co-rdmode").addEventListener("change", updCosmo);
    var pc = $("co-prior");
    co.h0p = slider(pc, "co-h0p", "h0prior", { min: 65, max: 80, step: 0.01, value: 72.8 }, null);
    co.h0s = slider(pc, "co-h0s", "h0sig", { min: 0.1, max: 5, step: 0.01, value: 1.0 }, null);
    co.rdp = slider(pc, "co-rdp", "rdprior", { min: 130, max: 160, step: 0.01, value: 147.09 }, null);
    co.rds = slider(pc, "co-rds", "rdsig", { min: 0.05, max: 5, step: 0.01, value: 0.26 }, null);
    $("co-refit").addEventListener("click", function () { runFits(); updCosmo(); });
    var pres = [["m4", "pre_m4"], ["rigid", "pre_rigid"], ["same", "pre_same"], ["low", "pre_low"], ["m3", "pre_m3"], ["p2", "pre_2p"], ["joint", "pre_joint"]];
    pres.forEach(function (p) { var b = document.createElement("button"); b.type = "button"; b.setAttribute("data-t", p[1]); b.innerHTML = t(p[1]); b.dataset.p = p[0]; b.addEventListener("click", function () { preset(p[0]); }); $("co-presets").appendChild(b); });
    runFits();
  }
  function runFits() {
    co.ctl = T.baoControls(); co.fit = T.profile2p(); co.m3 = T.m3Branch();
    co.joint = T.jointFit(co.h0p.get(), co.h0s.get(), co.rdp.get(), co.rds.get());
  }
  function preset(p) {
    var bg = $("co-bg"), rm = $("co-rdmode");
    if (p === "m4") { bg.value = "shape"; rm.value = "own"; co.H0.set(67.4, true); }
    if (p === "rigid") { bg.value = "shape"; rm.value = "rd"; co.H0.set(72.8, true); co.rd.set(T.RD_ZDRAG, true); }
    if (p === "same") { bg.value = "full"; rm.value = "own"; co.H0.set(72.8, true); co.Om.set(T.OMEGA_M_0, true); }
    if (p === "low") { bg.value = "full"; rm.value = "own"; co.H0.set(72.8, true); co.Om.set(0.244828, true); }
    if (p === "m3") { bg.value = "m3"; rm.value = "own"; }
    if (p === "p2") { bg.value = "flat2p"; rm.value = "P"; co.H0.set(72.8, true); co.Om.set(co.fit.om, true); co.P.set(co.fit.P, true); }
    if (p === "joint") { bg.value = "flat2p"; rm.value = "rd"; co.H0.set(co.joint.h0, true); co.Om.set(co.joint.om, true); co.rd.set(co.joint.rd, true); }
    co.preset = p; updCosmo();
  }
  var fullCache = {};
  function fullBG(H0, Om) { var k = H0.toPrecision(12) + "_" + Om.toPrecision(12); if (!fullCache[k]) { var bg = T.makeBackground(H0, Om); fullCache[k] = { bg: bg, cmb: T.cmbForBackground(bg) }; } return fullCache[k]; }
  function cosmoModel() {
    var mode = $("co-bg").value, H0 = co.H0.get(), Om = co.Om.get(), m = { mode: mode, om_k: 0, H0c: H0 };
    if (mode === "m3") { var b = co.m3; m.Hz = b.Hz; m.om_k = b.om_k; m.H0c = T.H0_STEP25; m.H0 = T.H0_STEP25; m.rs = b.rs; m.rdOwn = b.rd; m.theta = b.theta100; m.pl3 = T.planck3dChi2(T.BG_M4.omega_b_h2, T.BG_M4.omega_cb_h2, b.theta100); m.Om = T.BG_M4.Omega_m_0; }
    else if (mode === "shape") { var s = H0 / T.H0_STEP25; m.Hz = function (z) { return s * T.BG_M4.Hz(z); }; m.H0 = H0; m.rs = T.RS_ZSTAR / s; m.rdOwn = T.RD_ZDRAG / s; var cm = T.cmbForBackground(T.BG_M4); m.theta = cm.theta100; m.pl3 = cm.chi2_3d; m.Om = T.BG_M4.Omega_m_0; }
    else if (mode === "full") { var f = fullBG(H0, Om); m.Hz = f.bg.Hz; m.H0 = H0; m.rs = f.cmb.rs; m.rdOwn = f.cmb.rd; m.theta = f.cmb.theta100; m.pl3 = f.cmb.chi2_3d; m.Om = Om; }
    else { var E = T.flatE(Om), f2 = fullBG(H0, Om); m.Hz = function (z) { return H0 * E(z); }; m.H0 = H0; m.rs = f2.cmb.rs; m.rdOwn = f2.cmb.rd; m.theta = f2.cmb.theta100; m.pl3 = f2.cmb.chi2_3d; m.Om = Om; }
    var rm = $("co-rdmode").value;
    m.rd = rm === "own" ? m.rdOwn : (rm === "rd" ? co.rd.get() : co.P.get() / m.H0);
    return m;
  }
  function updCosmo() {
    var mode = $("co-bg").value, rm = $("co-rdmode").value, m = cosmoModel();
    co.H0.disable(mode === "m3"); co.Om.disable(mode === "m3" || mode === "shape");
    co.rd.disable(rm !== "rd"); co.P.disable(rm !== "P");
    if (rm !== "rd") co.rd.set(m.rd, true); if (rm !== "P") co.P.set(m.H0 * m.rd, true);
    if (mode === "m3") co.H0.set(T.H0_STEP25, true);
    var pred = T.baoPred(m.Hz, m.rd, m.om_k, m.H0c), c2 = T.chi2Vec(pred), c2t = T.chi2BAOTable1(m.Hz, m.rd, m.om_k, m.H0c), p = T.chi2sf(c2, 12);
    var thp = (m.theta - T.THETA100_PL) / T.THETA100_PL_ERR;
    kv($("co-out"), [["c_chi2", fmt(c2, 9), c2 < 25 ? "ok" : "bad"], ["c_chi2t", fmt(c2t, 9)], ["c_dof", "12"], ["c_p", fmt(p, 4)], ["c_sig", fmt(T.normSigma(p), 4)],
      ["c_rd", fmt(m.rd, 9)], ["c_P", fmt(m.H0 * m.rd, 9)], ["c_rs", fmt(m.rs, 9)], ["c_th", fmt(m.theta, 9)], ["c_thp", fmt(thp, 6), Math.abs(thp) < 3 ? "ok" : "bad"], ["c_pl3", fmt(m.pl3, 7)], ["c_ok", fmt(m.om_k, 7)]]);
    document.querySelectorAll("#co-presets button").forEach(function (b) { b.classList.toggle("on", b.dataset.p === co.preset); });
    co.preset = null;
    // BAO chart
    var zs = linspace(0.02, 2.5, 110), ref = { Hz: T.BG_M4.Hz, rd: T.RD_ZDRAG }, cur = [], m4 = [];
    zs.forEach(function (z) { cur.push(T.baoObs(z, m.Hz, m.rd, m.om_k, m.H0c)); m4.push(T.baoObs(z, ref.Hz, ref.rd, 0, T.H0_STEP25)); });
    var Q = { DM_over_rs: ["DM", C.acc], DH_over_rs: ["DH", C.acc2], DV_over_rs: ["DV", C.ok] }, ser = [];
    ["DM", "DH", "DV"].forEach(function (k, j) {
      var col = [C.acc, C.acc2, C.ok][j];
      ser.push({ label: k + "/r_d", data: zs.map(function (z, i) { return [z, cur[i][k]]; }), color: col });
      ser.push({ label: j === 0 ? "M4 (H0 = 67.4)" : null, data: zs.map(function (z, i) { return [z, m4[i][k]]; }), color: col, dash: [5, 4], width: 1.2 });
    });
    D.desi.q.forEach(function (q, i) { var qq = Q[q] || Q.DV_over_rs; ser.push({ data: [[D.desi.z[i], D.desi.mean[i], T.DESI_SIG[i]]], type: "points", color: qq[1], marker: "square" }); });
    P.plot($("co-bao"), { height: 380, xlabel: "z", ylabel: "D / r_d", series: ser, legend: "tl", xlim: [0, 2.5] });
    var pulls = pred.map(function (v, i) { return [i + 1, (D.desi.mean[i] - v) / T.DESI_SIG[i]]; });
    var ps = []; D.desi.q.forEach(function (q, i) { ps.push({ data: [pulls[i]], type: "points", color: (Q[q] || Q.DV_over_rs)[1], size: 5 }); });
    P.plot($("co-pull"), { height: 230, xlabel: "# (z = 0.30 … 2.33)", ylabel: "pull [σ]", xlim: [0.5, 12.5], series: ps, hlines: [{ y: 0, color: C.mute, dash: [] }, { y: 2, color: C.bad }, { y: -2, color: C.bad }] });
    updFits(m);
    // H(z)/(1+z)
    var hz = zs.map(function (z) { return [z, m.Hz(z) / (1 + z)]; }), hm = zs.map(function (z) { return [z, T.BG_M4.Hz(z) / (1 + z)]; }), pts = [];
    D.desi.q.forEach(function (q, i) { if (q === "DH_over_rs") { var z = D.desi.z[i], H = T.C_KM_S / (D.desi.mean[i] * m.rd); pts.push([z, H / (1 + z), H / (1 + z) * T.DESI_SIG[i] / D.desi.mean[i]]); } });
    P.plot($("co-hz"), { height: 230, xlabel: "z", ylabel: "H(z)/(1+z) [km/s/Mpc]", series: [{ label: "model", data: hz, color: C.acc }, { label: "M4", data: hm, color: C.acc2, dash: [5, 4] }, { label: "DESI DH (r_d model)", data: pts, type: "points", color: C.txt }], legend: "tr" });
  }
  function updFits(m) {
    var f = co.fit, j = co.joint, c = co.ctl, m3 = co.m3, h0p = co.h0p.get(), h0s = co.h0s.get();
    var rq = f.P / h0p, srq = f.sigmaP / h0p, srqH = rq * Math.sqrt(Math.pow(f.sigmaP / f.P, 2) + Math.pow(h0s / h0p, 2));
    kv($("co-fits"), [
      ["f_m4", fmt(c.chi2_m4, 9)], ["f_rigid", fmt(c.chi2_rigid, 9)],
      ["f_1p", "P = " + fmt(c.P1, 8) + " ± " + fmt(c.sigmaP1, 5) + ", χ² = " + fmt(c.chi2_P1, 7)],
      ["f_2p", "Ω_m = " + fmt(f.om, 7) + " ± " + fmt(f.sigmaOm, 4)], ["P", fmt(f.P, 9) + " ± " + fmt(f.sigmaP, 5)], ["χ²<sub>min</sub>", fmt(f.chi2, 9)],
      ["f_req", fmt(rq, 8) + " ± " + fmt(srq, 4)], ["f_reqH", "± " + fmt(srqH, 5)], ["f_h0pl", fmt(f.H0fromPlanckRd, 7) + " ± " + fmt(f.sigmaH0fromPlanckRd, 4)],
      ["f_joint", fmt(j.om, 6) + ", " + fmt(j.h0, 7) + ", " + fmt(j.rd, 8)],
      ["f_parts", fmt(j.parts.bao, 6) + " + " + fmt(j.parts.rd, 4) + " + " + fmt(j.parts.h0, 5) + " = " + fmt(j.chi2, 7)], ["f_dchi", fmt(j.chi2 - f.chi2, 7)],
      ["f_m3", fmt(m3.rs, 9) + " / " + fmt(m3.rd, 9)], ["f_m3b", fmt(m3.theta100, 7) + ", " + fmt(m3.pullTheta, 7) + ", " + fmt(m3.chi2_bao, 7)],
      ["f_m3k", fmt(m3.om_k, 8) + " / " + fmt(m3.om_de, 8)], ["f_m2", fmt(m3.E2_m2_z0, 8) + ", z = " + fmt(m3.z_root_m2, 8)]]);
    var Hs = linspace(55, 85, 120), Pf = f.P, sP = f.sigmaP;
    P.plot($("co-ladder"), { height: 380, xlabel: "H0 [km/s/Mpc]", ylabel: "r_d [Mpc]", xlim: [60, 80], ylim: [125, 160],
      hbands: [{ y0: co.rdp.get() - co.rds.get(), y1: co.rdp.get() + co.rds.get(), color: "rgba(255,180,84,0.28)" }],
      vbands: [{ x0: h0p - h0s, x1: h0p + h0s, color: "rgba(181,140,255,0.18)" }],
      series: [
        { label: "DESI BAO: r_d = P/H0 ± σ", type: "band", data: Hs.map(function (h) { return [h, (Pf - sP) / h, (Pf + sP) / h]; }), color: C.acc, alpha: 0.28 },
        { data: Hs.map(function (h) { return [h, Pf / h]; }), color: C.acc, width: 1.5 },
        { label: "Planck r_d = 147.09 ± 0.26", data: [[0, -1]], color: C.acc2, noRange: true },
        { label: "H0 prior " + fmt(h0p, 4) + " ± " + fmt(h0s, 3), data: [[0, -1]], color: C.vio, noRange: true },
        { label: "M4 (67.4, 147.09)", data: [[67.4, T.RD_ZDRAG]], type: "points", color: C.acc2, size: 5 },
        { label: "r_d req @ " + fmt(h0p, 4), data: [[h0p, rq, srq]], type: "points", color: C.ok, size: 5 },
        { label: "joint best fit", data: [[j.h0, j.rd]], type: "points", color: C.vio, size: 5.5 },
        { label: "current model", data: [[m.H0, m.rd]], type: "star", color: C.txt }
      ], legend: "tr" });
  }

  // ---------------- Tab 5: further probes ----------------
  var pr = {};
  function initProbes() {
    var c = $("pr-ctrl");
    pr.H0 = slider(c, "pr-H0", "H0", { min: 55, max: 85, step: 0.01, value: 67.4 }, updProbes);
    pr.Om = slider(c, "pr-Om", "Om", { min: 0.15, max: 0.5, step: 0.0001, value: T.OMEGA_M_0, prec: 9 }, updProbes);
    pr.L = slider(c, "pr-L", "L_km", { min: 0.1, max: 1e5, value: 10, log: true }, updProbes);
    pr.a = slider(c, "pr-a", "alphaY", { min: 0, max: 3, step: 0.001, value: 1 }, updProbes);
    pr.rho = slider(c, "pr-rho", "rho_custom", { min: 1e-28, max: 1e18, value: 1e3, log: true }, updProbes);
    [["m4", "pre_m4"], ["same", "pre_same"], ["low", "pre_low"]].forEach(function (p) {
      var b = document.createElement("button"); b.type = "button"; b.setAttribute("data-t", p[1]); b.innerHTML = t(p[1]);
      b.addEventListener("click", function () { pr.H0.set(p[0] === "m4" ? 67.4 : 72.8, true); pr.Om.set(p[0] === "low" ? 0.244828 : T.OMEGA_M_0, true); updProbes(); });
      $("pr-presets").appendChild(b);
    });
    $("pr-allm").addEventListener("change", updProbes);
  }
  function updProbes() {
    var H0 = pr.H0.get(), Om = pr.Om.get(), L = pr.L.get(), a = pr.a.get(), f = fullBG(H0, Om), bg = f.bg, cm = f.cmb;
    kv($("pr-cmb"), [["c_th", fmt(cm.theta100, 9)], ["c_thp", fmt(cm.pull, 6), Math.abs(cm.pull) < 3 ? "ok" : "bad"], ["c_rs", fmt(cm.rs, 9)], ["c_rd", fmt(cm.rd, 9)], ["c_pl3", fmt(cm.chi2_3d, 7)], ["ω<sub>cb</sub>", fmt(bg.omega_cb_h2, 7)]]);
    var s = T.sneControl(bg);
    P.plot($("pr-sne"), { xlabel: "z", ylabel: "μ_model − μ_ref [mag]", series: [
      { label: "unshifted", data: s.rows.map(function (r) { return [r.z, r.d, r.sig]; }), type: "points", color: C.acc2, hollow: true },
      { label: "shifted (offset " + fmt(s.offset, 4) + ")", data: s.rows.map(function (r) { return [r.z, r.d - s.offset, r.sig]; }), type: "points", color: C.acc }
    ], hlines: [{ y: 0, color: C.mute, dash: [] }], legend: "tr" });
    kv($("pr-sne-out"), [["n_off", fmt(s.offset, 9)], ["n_c0", fmt(s.chi2_unshifted, 9)], ["n_c1", fmt(s.chi2_shifted, 9)]]);
    var g = T.growth(bg, { L_km: L, alphaY: a, cbOnly: !$("pr-allm").checked }), r = T.rsdChi2(g), zz = linspace(0, 1.6, 120), q = T.qCoupling(a);
    P.plot($("pr-rsd"), { xlabel: "z", ylabel: "fσ8", series: [
      { label: "model", data: zz.map(function (z) { return [z, g.at(z).fs8]; }), color: C.acc },
      { label: "RSD data", data: D.rsd.map(function (d) { return [d[1], d[2], d[3]]; }), type: "points", color: C.txt }
    ], legend: "br" });
    kv($("pr-rsd-out"), [["r_chi2", fmt(r.chi2, 9)], ["r_max", fmt(r.maxDelta, 5)], ["r_kL", fmt(g.kL0, 6)], ["r_dG", fmt(q * g.kL0 * g.kL0 / (1 + g.kL0 * g.kL0), 5)]]);
    var rows = R.bbn.map(function (e) { var b = T.bbn(e.z_epoch, L, a); return "<tr><td>" + e.epoch_id.replace(/^T_(\d+)p(\d+)MeV_(.*)$/, function (m0, a1, a2, rest) { return "T = " + a1 + "." + a2 + " MeV (" + rest.replace(/_/g, " ") + ")"; }) + '</td><td class="n">' + fmt(e.z_epoch, 3) + '</td><td class="n">' + fmt(b.H, 6) + '</td><td class="n">' + fmt(b.H_over_mphi, 5) + '</td><td class="n">' + fmt(b.absPhi, 5) + '</td><td class="n">' + fmt(b.dHH, 5) + '</td><td class="n">' + fmt(b.Yp, 6) + " (" + fmt(b.pullYp, 3) + ')</td><td class="n">' + fmt(b.DH, 6) + " (" + fmt(b.pullDH, 3) + ")</td></tr>"; });
    $("pr-bbn").innerHTML = '<table class="tbl"><tr>' + ["b_epoch", "b_z", "b_H", "b_Hm", "b_phi", "b_dH", "b_Yp", "b_DH"].map(function (k) { return "<th>" + t(k) + "</th>"; }).join("") + "</tr>" + rows.join("") + "</table>";
    var lr = R.local.map(function (e) {
      var o = e.rho_kg_m3 === null ? null : T.localBulk(e.rho_kg_m3, L, a);
      var ph = o ? o.absPhi : e.abs_phi_bulk, dH = o ? o.dH0 : T.H0_STEP25 * Math.expm1(ph);
      return "<tr><td>" + e.env_id.replace(/^\d+_/, "").replace(/_/g, " ") + '</td><td class="n">' + (e.rho_kg_m3 === null ? "—" : fmt(e.rho_kg_m3, 4)) + '</td><td class="n">' + fmt(ph, 5) + '</td><td class="n">' + fmt(dH, 5) + '</td><td class="n">' + fmt(dH / (T.H0_LOCAL - T.H0_STEP25), 4) + "</td></tr>";
    });
    var oc = T.localBulk(pr.rho.get(), L, a);
    lr.push('<tr><td><b>' + t("rho_custom") + '</b></td><td class="n">' + fmt(pr.rho.get(), 4) + '</td><td class="n">' + fmt(oc.absPhi, 5) + '</td><td class="n">' + fmt(oc.dH0, 5) + '</td><td class="n">' + fmt(oc.frac, 4) + "</td></tr>");
    $("pr-local").innerHTML = '<table class="tbl"><tr>' + ["lc_env", "lc_rho", "lc_phi", "lc_dH0", "lc_frac"].map(function (k) { return "<th>" + t(k) + "</th>"; }).join("") + "</tr>" + lr.join("") + "</table>";
  }

  // ---------------- Tab 6: info + reproduction ----------------
  function reproRows() { return window.TAL_reproRows(co); }
  function relDev(a, b) { if (b === 0) return Math.abs(a); return Math.abs(a - b) / Math.abs(b); }
  function updInfo() {
    var L = 10;
    kv($("in-params"), [["L", "10 km"], ["α<sub>Y</sub>", "1"], ["β = 1/(4π)", fmt(T.BETA, 9)], ["q = α<sub>Y</sub>²β", fmt(T.qCoupling(1), 9)], ["γ<sub>L</sub> [m³/kg]", fmt(T.gammaL(L), 8)],
      ["m<sub>φ</sub>c² [eV]", fmt(T.HBAR_C_EV_M / (L * 1000), 6)], ["f<sub>cut</sub> [Hz]", fmt(R.gw_meta.f_cut_Hz, 8)], ["H₀ (M4)", "67.40"], ["ω<sub>b</sub>, ω<sub>c</sub>, ω<sub>ν</sub>", "0.02237, 0.12000, 0.000645"],
      ["N<sub>eff</sub>", "3.046"], ["z<sub>*</sub>, z<sub>drag</sub>", "1089.80, 1059.94"], ["Planck r<sub>d</sub>", "147.09 ± 0.26 Mpc"], ["Ω<sub>m</sub> (M4)", fmt(T.OMEGA_M_0, 9)]]);
    $("in-limits").innerHTML = ["lim1", "lim2", "lim3", "lim4", "lim5", "lim6"].map(function (k) { return "<li>" + t(k) + "</li>"; }).join("");
    if (!window.TAL_REPRO) window.TAL_REPRO = reproRows();
    var html = '<table class="tbl"><tr><th>' + t("rp_q") + "</th><th>" + t("rp_js") + "</th><th>" + t("rp_ref") + "</th><th>" + t("rp_rel") + "</th></tr>";
    window.TAL_REPRO.forEach(function (r) {
      var d = relDev(r[1], r[2]), cls = d < 1e-6 ? "ok" : (d < 1e-3 ? "" : "bad");
      html += "<tr><td>" + r[0] + '</td><td class="n">' + fmt(r[1], 12) + '</td><td class="n">' + fmt(r[2], 12) + '</td><td class="n ' + cls + '">' + (d === 0 ? "0" : d.toExponential(1)) + "</td></tr>";
    });
    $("in-repro").innerHTML = html + "</table>";
  }

  // ---------------- language, tabs ----------------
  var updaters = { yukawa: updYukawa, sparc: updSparc, gw: function () { updGW(); updLens(); }, cosmo: updCosmo, probes: updProbes, info: updInfo };
  var current = "cosmo";
  function applyLang() {
    document.documentElement.lang = lang;
    document.querySelectorAll("[data-t]").forEach(function (el) { el.innerHTML = t(el.getAttribute("data-t")); });
    var opt = { "co-bg": ["bg_shape", "bg_full", "bg_flat2p", "bg_m3"], "co-rdmode": ["rd_own", "rd_manual", "rd_P"], "ln-body": ["ln_sun", "ln_ns"] };
    Object.keys(opt).forEach(function (id) { var os = $(id).options; opt[id].forEach(function (k, i) { os[i].textContent = plain(t(k)); }); });
    document.querySelectorAll("#yk-body option").forEach(function (o) { o.textContent = plain(t(o.getAttribute("data-t"))); });
    $("lang-hu").classList.toggle("on", lang === "hu"); $("lang-en").classList.toggle("on", lang === "en");
  }
  function show(tab) {
    if (!updaters[tab]) tab = "cosmo"; current = tab;
    document.querySelectorAll(".tab").forEach(function (s) { s.classList.toggle("on", s.id === "tab-" + tab); });
    document.querySelectorAll("#tabs button").forEach(function (b) { b.classList.toggle("on", b.dataset.tab === tab); });
    updaters[tab]();
  }
  function parseHash() {
    var h = (location.hash || "").replace(/^#/, "").split("&"), o = {};
    h.forEach(function (p) { if (!p) return; var kvp = p.split("="); if (kvp.length === 1) o.tab = kvp[0]; else o[kvp[0]] = kvp[1]; });
    return o;
  }
  function setLang(l) { lang = l; try { window.localStorage.setItem("tal_lang", l); } catch (e) { /* storage unavailable */ } applyLang(); updaters[current](); }
  function boot() {
    var h = parseHash(); if (h.lang === "en" || h.lang === "hu") lang = h.lang;
    initYukawa(); initSparc(); initGW(); initCosmo(); initProbes();
    $("co-bg").value = "shape"; $("co-rdmode").value = "rd"; co.preset = "rigid";
    applyLang();
    document.querySelectorAll("#tabs button").forEach(function (b) { b.addEventListener("click", function () { var hh = parseHash(); location.hash = b.dataset.tab + (hh.lang ? "&lang=" + hh.lang : ""); }); });
    $("lang-hu").addEventListener("click", function () { setLang("hu"); }); $("lang-en").addEventListener("click", function () { setLang("en"); });
    window.addEventListener("hashchange", function () { var hh = parseHash(); if (hh.lang && hh.lang !== lang) { lang = hh.lang; applyLang(); } show(hh.tab || current); });
    var rt; window.addEventListener("resize", function () { clearTimeout(rt); rt = setTimeout(function () { updaters[current](); }, 120); });
    show(h.tab || "yukawa");
    document.body.setAttribute("data-ready", "1");
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})();
