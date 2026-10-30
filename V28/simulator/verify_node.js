#!/usr/bin/env node
/* Runs the computation scripts embedded in index.html under Node (vm) and compares them with the package JSON values. */
"use strict";
const fs = require("fs"), vm = require("vm"), path = require("path");
const html = fs.readFileSync(path.join(__dirname, "index.html"), "utf8");
const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m => m[1]);
const ctx = { console, Math, devicePixelRatio: 1 }; ctx.window = ctx; ctx.globalThis = ctx; vm.createContext(ctx);
// data, core, chart, i18n, repro (the last script is the DOM user interface and is not executed here)
scripts.slice(0, 5).forEach((s, i) => vm.runInContext(s, ctx, { filename: "inline_" + i + ".js" }));
const T = ctx.TAL, R = ctx.TAL_DATA.ref;
const t0 = Date.now(), rows = ctx.TAL_reproRows({});
// explicit checks of the UI computation paths used by the cosmology tab presets
const fit = T.profile2p(), j = T.jointFit();
const flat = (om, H0, rd) => { const E = T.flatE(om); return T.chi2BAO(z => H0 * E(z), rd, 0, H0); };
rows.push(["UI path: rigid preset χ² (shape ×72.8/67.4, r_d = 147.09)", T.chi2BAO(z => 72.8 / 67.4 * T.BG_M4.Hz(z), T.RD_ZDRAG), R.chi2_rigid]);
rows.push(["UI path: 2-param preset χ² (flat E, r_d = P/H0)", flat(fit.om, 72.8, fit.P / 72.8), R.chi2_2p]);
rows.push(["UI path: joint preset χ²_BAO", flat(j.om, j.h0, j.rd), R.joint.chi2_BAO]);
const m3 = T.m3Branch(); rows.push(["UI path: M3 preset χ²_BAO", T.chi2BAO(m3.Hz, m3.rd, m3.om_k, 67.4), R.m3.chi2_bao_m3_full]);
const ms = Date.now() - t0;
const out = rows.map(r => ({ quantity: r[0], js: r[1], package: r[2], rel_dev: r[2] === 0 ? Math.abs(r[1]) : Math.abs(r[1] - r[2]) / Math.abs(r[2]) }));
const tolStrict = 1e-6;
const loose = out.filter(r => r.rel_dev >= tolStrict);
fs.writeFileSync(path.join(__dirname, "verify_node_out.json"), JSON.stringify({ n_rows: out.length, max_rel_dev_strict_rows: Math.max(...out.filter(r => r.rel_dev < tolStrict).map(r => r.rel_dev)), rows_above_1e_6: loose, runtime_ms: ms, rows: out }, null, 1));
out.forEach(r => console.log((r.rel_dev < tolStrict ? "OK  " : "DIFF") + "  " + r.quantity.padEnd(62) + String(r.js).padStart(24) + String(r.package).padStart(24) + "  " + r.rel_dev.toExponential(2)));
console.log("\nrows:", out.length, " above 1e-6:", loose.length, " runtime:", ms, "ms");
