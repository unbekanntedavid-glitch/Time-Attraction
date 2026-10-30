/* Minimal canvas plotting helper (no external dependency). */
(function (root) {
  "use strict";
  function niceStep(range, n) {
    var raw = range / n, p = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / p;
    return (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p;
  }
  function fmtTick(v, log) {
    if (log) { var e = Math.round(Math.log10(v)); if (Math.abs(e) <= 2) return String(+v.toPrecision(3)); return "10" + sup(e); }
    var a = Math.abs(v); if (a !== 0 && (a < 1e-3 || a >= 1e5)) return v.toExponential(1);
    return String(+v.toPrecision(6));
  }
  var SUP = { "-": "\u207b", "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3", "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079" };
  function sup(n) { return String(n).split("").map(function (c) { return SUP[c] || c; }).join(""); }
  function plot(canvas, cfg) {
    var dpr = root.devicePixelRatio || 1, W = canvas.clientWidth || 600, H = cfg.height || canvas.clientHeight || 340;
    canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr); canvas.style.height = H + "px";
    var ctx = canvas.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    var th = cfg.theme || { bg: "#121a2b", grid: "rgba(160,180,220,0.13)", axis: "#8ea0c0", text: "#d9e2f2" };
    ctx.fillStyle = th.bg; ctx.fillRect(0, 0, W, H);
    var m = { l: 64, r: 16, t: cfg.title ? 30 : 14, b: 46 };
    var pw = W - m.l - m.r, ph = H - m.t - m.b;
    var xs = [], ys = [];
    (cfg.series || []).forEach(function (s) {
      if (s.noRange) return;
      (s.data || []).forEach(function (p) {
        if (!isFinite(p[0]) || !isFinite(p[1])) return; xs.push(p[0]);
        var e = p[2] || 0; ys.push(p[1] - e); ys.push(p[1] + e);
        if (s.type === "band") ys.push(p[2]);
      });
    });
    function rng(arr, lim, log) {
      if (lim) return lim.slice();
      var a = arr.filter(function (v) { return isFinite(v) && (!log || v > 0); });
      if (!a.length) return [0, 1];
      var lo = Math.min.apply(null, a), hi = Math.max.apply(null, a);
      if (log) return [lo / 1.3, hi * 1.3];
      if (lo === hi) { lo -= 1; hi += 1; }
      var pad = 0.06 * (hi - lo); return [lo - pad, hi + pad];
    }
    var xr = rng(xs, cfg.xlim, cfg.xlog), yr = rng(ys, cfg.ylim, cfg.ylog);
    function tx(v) { return cfg.xlog ? Math.log10(v) : v; }
    function ty(v) { return cfg.ylog ? Math.log10(v) : v; }
    var X0 = tx(xr[0]), X1 = tx(xr[1]), Y0 = ty(yr[0]), Y1 = ty(yr[1]);
    function px(v) { return m.l + (tx(v) - X0) / (X1 - X0) * pw; }
    function py(v) { return m.t + ph - (ty(v) - Y0) / (Y1 - Y0) * ph; }
    ctx.font = "12px system-ui, -apple-system, Segoe UI, sans-serif"; ctx.lineWidth = 1;
    function ticks(lo, hi, log) {
      var out = [];
      if (log) { var e0 = Math.ceil(Math.log10(lo) - 1e-9), e1 = Math.floor(Math.log10(hi) + 1e-9), st = Math.max(1, Math.ceil((e1 - e0 + 1) / 8)); for (var e = e0; e <= e1; e += st) out.push(Math.pow(10, e)); return out; }
      var s = niceStep(hi - lo, 6), v = Math.ceil(lo / s) * s; for (; v <= hi + 1e-12 * s; v += s) out.push(Math.abs(v) < s * 1e-9 ? 0 : v); return out;
    }
    ctx.strokeStyle = th.grid; ctx.fillStyle = th.axis; ctx.textAlign = "center"; ctx.textBaseline = "top";
    ticks(xr[0], xr[1], cfg.xlog).forEach(function (v) { var x = px(v); ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph); ctx.stroke(); ctx.fillText(fmtTick(v, cfg.xlog), x, m.t + ph + 6); });
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    ticks(yr[0], yr[1], cfg.ylog).forEach(function (v) { var y = py(v); ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(m.l + pw, y); ctx.stroke(); ctx.fillText(fmtTick(v, cfg.ylog), m.l - 6, y); });
    ctx.strokeStyle = th.axis; ctx.strokeRect(m.l, m.t, pw, ph);
    ctx.fillStyle = th.text; ctx.textAlign = "center"; ctx.textBaseline = "alphabetic";
    if (cfg.xlabel) ctx.fillText(cfg.xlabel, m.l + pw / 2, H - 8);
    if (cfg.title) { ctx.font = "600 13px system-ui, sans-serif"; ctx.fillText(cfg.title, m.l + pw / 2, 18); ctx.font = "12px system-ui, sans-serif"; }
    if (cfg.ylabel) { ctx.save(); ctx.translate(14, m.t + ph / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(cfg.ylabel, 0, 0); ctx.restore(); }
    ctx.save(); ctx.beginPath(); ctx.rect(m.l, m.t, pw, ph); ctx.clip();
    (cfg.vbands || []).forEach(function (b) { ctx.fillStyle = b.color; ctx.fillRect(px(b.x0), m.t, px(b.x1) - px(b.x0), ph); });
    (cfg.hbands || []).forEach(function (b) { ctx.fillStyle = b.color; ctx.fillRect(m.l, py(b.y1), pw, py(b.y0) - py(b.y1)); });
    (cfg.series || []).forEach(function (s) {
      var d = (s.data || []).filter(function (p) { return isFinite(p[0]) && isFinite(p[1]) && (!cfg.xlog || p[0] > 0) && (!cfg.ylog || p[1] > 0); });
      ctx.strokeStyle = s.color || "#5cc8ff"; ctx.fillStyle = s.color || "#5cc8ff"; ctx.lineWidth = s.width || 2;
      ctx.setLineDash(s.dash || []);
      if (s.type === "band") {
        ctx.globalAlpha = s.alpha || 0.25; ctx.beginPath();
        d.forEach(function (p, i) { var x = px(p[0]), y = py(p[1]); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); });
        for (var i = d.length - 1; i >= 0; i--) ctx.lineTo(px(d[i][0]), py(d[i][2]));
        ctx.closePath(); ctx.fill(); ctx.globalAlpha = 1;
      } else if (s.type === "points") {
        ctx.setLineDash([]); ctx.lineWidth = 1.3;
        d.forEach(function (p) {
          var x = px(p[0]), y = py(p[1]);
          if (p[2]) { ctx.beginPath(); ctx.moveTo(x, py(p[1] - p[2])); ctx.lineTo(x, py(p[1] + p[2])); ctx.stroke(); }
          ctx.beginPath(); if (s.marker === "square") ctx.rect(x - 3.5, y - 3.5, 7, 7); else ctx.arc(x, y, s.size || 3.6, 0, 2 * Math.PI);
          if (s.hollow) ctx.stroke(); else ctx.fill();
        });
      } else if (s.type === "star") {
        d.forEach(function (p) { var x = px(p[0]), y = py(p[1]); ctx.beginPath(); for (var k = 0; k < 10; k++) { var r = k % 2 ? 4 : 9, a = -Math.PI / 2 + k * Math.PI / 5; ctx.lineTo(x + r * Math.cos(a), y + r * Math.sin(a)); } ctx.closePath(); ctx.fill(); });
      } else {
        ctx.beginPath(); var started = false;
        d.forEach(function (p) { var x = px(p[0]), y = py(p[1]); if (started) ctx.lineTo(x, y); else { ctx.moveTo(x, y); started = true; } });
        ctx.stroke();
      }
      ctx.setLineDash([]);
    });
    (cfg.hlines || []).forEach(function (h) { ctx.strokeStyle = h.color || th.axis; ctx.setLineDash(h.dash || [5, 4]); ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(m.l, py(h.y)); ctx.lineTo(m.l + pw, py(h.y)); ctx.stroke(); ctx.setLineDash([]); });
    (cfg.vlines || []).forEach(function (h) { ctx.strokeStyle = h.color || th.axis; ctx.setLineDash(h.dash || [5, 4]); ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(px(h.x), m.t); ctx.lineTo(px(h.x), m.t + ph); ctx.stroke(); ctx.setLineDash([]); if (h.label) { ctx.fillStyle = h.color || th.axis; ctx.textAlign = "left"; ctx.fillText(h.label, px(h.x) + 4, m.t + 14); } });
    ctx.restore();
    var items = (cfg.series || []).filter(function (s) { return s.label; });
    if (items.length) {
      ctx.font = "11.5px system-ui, sans-serif"; var lw = 0;
      items.forEach(function (s) { lw = Math.max(lw, ctx.measureText(s.label).width); });
      var bw = lw + 38, bh = items.length * 17 + 8, pos = cfg.legend || "tr";
      var bx = pos.indexOf("l") >= 0 ? m.l + 8 : m.l + pw - bw - 8, by = pos.indexOf("b") >= 0 ? m.t + ph - bh - 8 : m.t + 8;
      ctx.fillStyle = "rgba(10,15,26,0.82)"; ctx.fillRect(bx, by, bw, bh); ctx.strokeStyle = th.grid; ctx.strokeRect(bx, by, bw, bh);
      items.forEach(function (s, i) {
        var y = by + 13 + i * 17; ctx.strokeStyle = s.color; ctx.fillStyle = s.color; ctx.lineWidth = 2; ctx.setLineDash(s.dash || []);
        if (s.type === "points" || s.type === "star") { ctx.beginPath(); ctx.arc(bx + 15, y, 3.5, 0, 2 * Math.PI); ctx.fill(); }
        else if (s.type === "band") { ctx.globalAlpha = 0.4; ctx.fillRect(bx + 6, y - 5, 18, 10); ctx.globalAlpha = 1; }
        else { ctx.beginPath(); ctx.moveTo(bx + 6, y); ctx.lineTo(bx + 24, y); ctx.stroke(); }
        ctx.setLineDash([]); ctx.fillStyle = th.text; ctx.textAlign = "left"; ctx.textBaseline = "middle"; ctx.fillText(s.label, bx + 30, y);
      });
    }
  }
  root.TALPlot = { plot: plot, sup: sup };
})(typeof window !== "undefined" ? window : globalThis);
