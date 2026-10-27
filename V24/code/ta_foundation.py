#!/usr/bin/env python3
"""
ta_foundation.py - "v27 foundation" mode (clean Python stand-in; NOT the user's app.tsk formula mode).

Separate from ta_pipeline.py (the photometric mode), which is kept unchanged.

REAL DATA ONLY
  Rotmod .dat files are read *directly from* /workspace/sparc-reference/Rotmod_LTG.zip after asserting
  its SHA256 == 0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588 (the hash quoted in the
  v16.0 PDF) and that it holds exactly 175 *_rotmod.dat files. The catalogue SPARC_Lelli2016c.mrt is
  asserted to have SHA256 5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243.
  No random/mock data is generated anywhere; the only RNG use is the seeded galaxy split (imported from
  ta_pipeline) and the seeded bootstrap of the paired test statistics.

BARYONS
  Vbar^2(r) = sign(Vgas) Vgas^2 + Ydisk Vdisk^2 + Ybul Vbul^2   (Ydisk configurable; default 0.5; Ybul 0.7)
  NGC3198 at 30.08 kpc: Vbar = 70.32 km/s for Ydisk 0.5 and 88.53 km/s for Ydisk 1.0.

DENSITY TEMPLATE (spherical, from the baryonic rotation template)
  rho(r) = 1/(4 pi G r^2) * d(r Vbar^2)/dr                [SI: r in m, V in m/s -> kg/m^3]
  Finite differences on the native SPARC radii (no resampling): numpy.gradient(r*Vbar^2, r, edge_order=1),
  i.e. 2nd-order non-uniform central difference at interior points
      f'(x_i) ~ [h_l^2 f_{i+1} - h_r^2 f_{i-1} + (h_r^2 - h_l^2) f_i] / (h_l h_r (h_l + h_r)),
  h_l = x_i - x_{i-1}, h_r = x_{i+1} - x_i, and 1st-order one-sided differences at the first/last radius.
  Edge points (first/last radius, one-sided) are FLAGGED 'edge' and excluded from C evaluation.

EFFECTIVE DENSITY (v16 Sec 3 structure, local + integrated)
  rho_int(r_i) = (1/L0) * Integral_0^{r_i} rho(s) ds, trapezoid on native radii using the SIGNED rho of all
  points where rho is finite (no clipping), rho held constant from r=0 to the innermost radius (v16 text).
  rho_eff = rho_local + rho_int.   L0 LOCKED (not fitted) = median catalogue Rdisk over the blind set.
  Variant 'local' uses rho_eff = rho_local only (no L0).
  C = 1 + C0 exp(-rho_eff/rho_c), C0=5, rho_c=1e-23 kg/m^3 (v16 constants); Vmod = sqrt(C) Vbar.

FLAGS (every flagged point listed; never clipped; excluded from C evaluation and from chi2 of ALL models,
so all paired comparisons use identical points):
  edge, rho_nan, rho<0, rho==0, rho_eff<=0, Vbar2<0, errV<=0
"""
import argparse, hashlib, io, json, math, os, sys, zipfile, datetime
import numpy as np
import ta_pipeline as tp

ZIP = "/workspace/sparc-reference/Rotmod_LTG.zip"
ZIP_SHA = "0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588"
MRT_SHA = "5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243"
G_SI = 6.67430e-11; KPC = tp.KPC
C0_V16, RHOC_V16, YBUL = 5.0, 1e-23, 0.7
A0_RAR = 1.2e-10
QUALITY = dict(Q_in=(1, 2), Inc_gt=30.0, min_points=5)   # v16 cuts (Q in {1,2}, Inc>30) + >=5 radii
OUT = os.path.join(tp.HERE, "results", "foundation")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def load_real_sparc():
    assert sha256(ZIP) == ZIP_SHA, "Rotmod_LTG.zip is not the verified SPARC file"
    assert sha256(tp.MRT) == MRT_SHA, "SPARC_Lelli2016c.mrt is not the verified SPARC file"
    gals = {}
    with zipfile.ZipFile(ZIP) as z:
        names = sorted(n for n in z.namelist() if n.endswith("_rotmod.dat"))
        assert len(names) == 175, f"expected 175 real SPARC rotmod files, got {len(names)}"
        for n in names:
            txt = z.read(n).decode()
            assert txt.startswith("# Distance"), f"{n}: not a SPARC rotmod header"
            a = np.loadtxt(io.StringIO(txt), comments="#", ndmin=2)
            assert a.shape[1] == 8 and a.shape[0] >= 1 and np.all(np.isfinite(a)), f"{n}: bad table"
            assert np.all(np.diff(a[:, 0]) > 0), f"{n}: radii not strictly increasing"
            gname = os.path.basename(n).replace("_rotmod.dat", "")
            gals[gname] = dict(zip(["R", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"], a.T))
    meta = tp.load_mrt()
    assert set(meta) == set(gals), "catalogue / rotmod name mismatch"
    for k in gals: gals[k]["meta"] = meta[k]
    return gals

def passes_quality(m, n):
    return m["Q"] in QUALITY["Q_in"] and m["Inc"] > QUALITY["Inc_gt"] and n >= QUALITY["min_points"]

def blind_set(gals, labels):
    return sorted(k for k, g in gals.items() if labels[k] in ("train", "test")
                  and passes_quality(g["meta"], len(g["R"])))

def compute(g, Ydisk, L0_kpc, variant):
    R = g["R"]; r_m = R * KPC
    Vbar2 = np.sign(g["Vgas"]) * g["Vgas"]**2 + Ydisk * g["Vdisk"]**2 + YBUL * g["Vbul"]**2   # (km/s)^2
    f = r_m * Vbar2 * 1e6                                                                       # m^3/s^2
    n = len(R)
    if n >= 2:
        dfdr = np.gradient(f, r_m, edge_order=1)
    else:
        dfdr = np.full(n, np.nan)
    rho = dfdr / (4 * np.pi * G_SI * r_m**2)                                                    # kg/m^3
    # integrated term (signed rho, trapezoid on native radii, constant to r=0)
    rho_int = np.full(n, np.nan)
    if variant == "integrated":
        fin = np.isfinite(rho)
        cum = np.zeros(n); acc = 0.0
        if fin[0]: acc = rho[0] * r_m[0]
        cum[0] = acc
        for i in range(1, n):
            if fin[i] and fin[i - 1]:
                acc += 0.5 * (rho[i] + rho[i - 1]) * (r_m[i] - r_m[i - 1])
            cum[i] = acc
        rho_int = cum / (L0_kpc * KPC)
        rho_eff = rho + rho_int
    else:
        rho_eff = rho.copy()
    flags = [[] for _ in range(n)]
    for i in range(n):
        if n < 3 or i == 0 or i == n - 1: flags[i].append("edge")
        if not np.isfinite(rho[i]): flags[i].append("rho_nan")
        elif rho[i] < 0: flags[i].append("rho<0")
        elif rho[i] == 0: flags[i].append("rho==0")
        if np.isfinite(rho_eff[i]) and rho_eff[i] <= 0 and "rho<0" not in flags[i]: flags[i].append("rho_eff<=0")
        if Vbar2[i] < 0: flags[i].append("Vbar2<0")
        if not g["errV"][i] > 0: flags[i].append("errV<=0")
    ok = np.array([len(x) == 0 for x in flags])
    C = np.full(n, np.nan); C[ok] = 1 + C0_V16 * np.exp(-rho_eff[ok] / RHOC_V16)
    Vbar = np.sign(Vbar2) * np.sqrt(np.abs(Vbar2))
    Vmod = np.full(n, np.nan); Vmod[ok] = np.sqrt(C[ok] * Vbar2[ok])
    gbar = np.where(Vbar2 > 0, Vbar2 * 1e6 / r_m, np.nan)
    with np.errstate(all="ignore"):
        Vrar = np.sqrt(gbar / (1 - np.exp(-np.sqrt(gbar / A0_RAR))) * r_m) / 1e3
    return dict(R=R, Vobs=g["Vobs"], errV=g["errV"], Vbar2=Vbar2, Vbar=Vbar, rho=rho, rho_int=rho_int,
                rho_eff=rho_eff, flags=flags, ok=ok, C=C, Vmod=Vmod, Vrar=Vrar)

def chi2_red(res, key):
    ok = res["ok"]
    if ok.sum() == 0: return None
    if key == "C1": V = np.sqrt(res["Vbar2"][ok])
    elif key == "model": V = res["Vmod"][ok]
    else: V = res["Vrar"][ok]
    return float(np.mean(((res["Vobs"][ok] - V) / res["errV"][ok])**2))

def paired_tests(d, seed=27, nboot=20000):
    from scipy.stats import wilcoxon, binomtest
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    nz = d[d != 0]
    w = wilcoxon(nz, alternative="greater", zero_method="wilcox") if len(nz) else None
    npos, nneg = int((nz > 0).sum()), int((nz < 0).sum())
    s = binomtest(npos, npos + nneg, 0.5, alternative="greater") if npos + nneg else None
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(nboot, len(d)))
    bm, bmed = d[idx].mean(1), np.median(d[idx], 1)
    return dict(n=len(d), n_nonzero=len(nz), mean=float(d.mean()), median=float(np.median(d)),
                wilcoxon_W=float(w.statistic) if w else None, wilcoxon_p_one_sided=float(w.pvalue) if w else None,
                sign_n_pos=npos, sign_n_neg=nneg, sign_p_one_sided=float(s.pvalue) if s else None,
                boot_mean_CI95=[float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5))],
                boot_median_CI95=[float(np.percentile(bmed, 2.5)), float(np.percentile(bmed, 97.5))],
                boot_n=nboot, boot_seed=seed)

def L0_from_blind(gals, blind):
    vals = [gals[k]["meta"]["Rdisk"] for k in blind]
    return float(np.median(vals)), dict(rule="median catalogue Rdisk (SPARC_Lelli2016c.mrt, [3.6] disk scale length, kpc) over blind set",
                                       blind_set_definition="train+revealed-test of ta_pipeline seed-27 split passing quality cuts",
                                       quality_cuts=QUALITY, N_blind=len(blind), user_blind_list_N=123,
                                       matches_user_N=(len(blind) == 123), mrt_sha256=MRT_SHA, locked=True, fitted=False)

def write_table(name, res, path_md, path_csv, Ydisk_pair=None):
    hdr = ["R_kpc", "Vobs", "errV", "Vbar_Y0.5", "Vbar_Y1.0", "rho_kgm3", "rho_int_kgm3", "rho_eff_kgm3", "flags", "C", "Vmod"]
    r05, r10 = Ydisk_pair
    with open(path_csv, "w") as fc, open(path_md, "w") as fm:
        fc.write(",".join(hdr) + "\n")
        fm.write(f"# {name} - foundation mode, density/C/Vmod computed with Ydisk={res['_Y']}, variant={res['_variant']}, L0={res['_L0']:.4g} kpc\n\n")
        fm.write("Real SPARC Rotmod_LTG (zip SHA256 verified). rho = d(r Vbar^2)/dr / (4 pi G r^2), numpy.gradient on native radii. "
                 "Flagged points: no C/Vmod (excluded). C = 1+5 exp(-rho_eff/1e-23). Velocities km/s.\n\n")
        fm.write("| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr) + "\n")
        for i in range(len(res["R"])):
            row = [res["R"][i], res["Vobs"][i], res["errV"][i], r05["Vbar"][i], r10["Vbar"][i], res["rho"][i],
                   res["rho_int"][i], res["rho_eff"][i], ";".join(res["flags"][i]) or "ok", res["C"][i], res["Vmod"][i]]
            fmt = [f"{x:.6g}" if isinstance(x, (float, np.floating)) else str(x) for x in row]
            fc.write(",".join(fmt) + "\n")
            b = "**" if abs(res["R"][i] - 30.08) < 1e-9 else ""
            fm.write("| " + " | ".join(f"{b}{x}{b}" for x in fmt) + " |\n")

def log(rec):
    rec["timestamp"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with open(tp.LOG, "a") as f: f.write(json.dumps(rec) + "\n")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["ngc3198", "run"])
    ap.add_argument("--ydisk", type=float, action="append", help="repeatable; default 0.5 and 1.0")
    ap.add_argument("--variant", choices=["integrated", "local"], default="integrated")
    a = ap.parse_args()
    Ys = a.ydisk or [0.5, 1.0]
    os.makedirs(OUT, exist_ok=True)
    gals = load_real_sparc()
    labels = tp.split(gals)                       # same seed-27 split as ta_pipeline
    blind = blind_set(gals, labels)
    L0, prov = L0_from_blind(gals, blind)
    print(f"blind set N={len(blind)} (user list N=123; match={len(blind)==123}); L0 = {L0:.4f} kpc (locked)")
    if a.cmd == "ngc3198":
        assert labels["NGC3198"] != "holdout"
        for Y in Ys:
            res = compute(gals["NGC3198"], Y, L0, a.variant); res.update(_Y=Y, _variant=a.variant, _L0=L0)
            pair = (compute(gals["NGC3198"], 0.5, L0, a.variant), compute(gals["NGC3198"], 1.0, L0, a.variant))
            tag = f"Y{Y}_{a.variant}"
            write_table("NGC3198", res, os.path.join(OUT, f"NGC3198_{tag}.md"), os.path.join(OUT, f"NGC3198_{tag}.csv"), pair)
            i = int(np.argmin(abs(res["R"] - 30.08)))
            print(json.dumps(dict(Ydisk=Y, R=res["R"][i], Vobs=res["Vobs"][i], errV=res["errV"][i], Vbar=res["Vbar"][i],
                                  rho=res["rho"][i], rho_int=res["rho_int"][i], rho_eff=res["rho_eff"][i],
                                  flags=res["flags"][i], C=res["C"][i], Vmod=res["Vmod"][i],
                                  n_flagged=int((~res["ok"]).sum()), n_points=len(res["R"]),
                                  flagged=[(float(res["R"][j]), res["flags"][j], float(res["rho"][j])) for j in range(len(res["R"])) if not res["ok"][j]])))
        return
    # run on all NON-HOLDOUT galaxies (holdout is never evaluated)
    visible = sorted(k for k in gals if labels[k] != "holdout")
    for Y in Ys:
        flagged_rows = []; per = {}
        for k in visible:
            res = compute(gals[k], Y, L0, a.variant)
            for j in range(len(res["R"])):
                if not res["ok"][j]:
                    flagged_rows.append((k, labels[k], float(res["R"][j]), ";".join(res["flags"][j]), float(res["rho"][j]),
                                         float(res["rho_eff"][j]) if np.isfinite(res["rho_eff"][j]) else float("nan")))
            per[k] = dict(split=labels[k], blind=k in blind, n=len(res["R"]), n_ok=int(res["ok"].sum()),
                          c1=chi2_red(res, "C1"), model=chi2_red(res, "model"), rar=chi2_red(res, "rar"),
                          n_rho_neg=sum("rho<0" in f for f in res["flags"]))
        tag = f"Y{Y}_{a.variant}"
        with open(os.path.join(OUT, f"flagged_points_{tag}.csv"), "w") as f:
            f.write("galaxy,split,R_kpc,flags,rho_kgm3,rho_eff_kgm3\n")
            for r in flagged_rows: f.write(",".join(str(x) if not isinstance(x, float) else f"{x:.6g}" for x in r) + "\n")
        with open(os.path.join(OUT, f"per_galaxy_{tag}.csv"), "w") as f:
            f.write("galaxy,split,blind,n,n_ok,n_rho_neg,chi2red_C1,chi2red_model,chi2red_RAR\n")
            for k, v in per.items():
                f.write(f"{k},{v['split']},{v['blind']},{v['n']},{v['n_ok']},{v['n_rho_neg']},{v['c1']},{v['model']},{v['rar']}\n")
        from collections import Counter
        fc = Counter(fl for r in flagged_rows for fl in r[3].split(";"))
        neg = [r for r in flagged_rows if "rho<0" in r[3]]
        near = [r for r in flagged_rows if abs(r[2] - 30.37) <= 0.5]
        near_neg = [r for r in neg if abs(r[2] - 30.37) <= 0.5]
        summary = dict(Ydisk=Y, variant=a.variant, L0_kpc=L0, n_gal_visible=len(visible),
                       n_points=sum(v["n"] for v in per.values()), n_points_used=sum(v["n_ok"] for v in per.values()),
                       n_points_flagged=len(flagged_rows), flag_counts=dict(fc),
                       n_rho_neg=len(neg), n_gal_with_rho_neg=len(set(r[0] for r in neg)),
                       near_30p37_within_0p5kpc_flagged=near, near_30p37_negative=near_neg)
        tests = {}
        for setname, keys in (("blind", [k for k in visible if per[k]["blind"]]), ("all_visible", visible)):
            ks = [k for k in keys if per[k]["model"] is not None]
            dC1 = [per[k]["c1"] - per[k]["model"] for k in ks]
            dRAR = [per[k]["rar"] - per[k]["model"] for k in ks]
            tests[setname] = dict(N_gal=len(ks),
                                  median_chi2red=dict(C1=float(np.median([per[k]["c1"] for k in ks])),
                                                      model=float(np.median([per[k]["model"] for k in ks])),
                                                      RAR=float(np.median([per[k]["rar"] for k in ks]))),
                                  model_vs_C1=paired_tests(dC1), model_vs_RAR=paired_tests(dRAR))
        summary["tests"] = tests
        json.dump(summary, open(os.path.join(OUT, f"summary_{tag}.json"), "w"), indent=2)
        log(dict(mode="foundation", expression="C = 1 + C0*exp(-rho_eff/rho_c), rho_eff = rho_local + (1/L0) int_0^r rho ds" if a.variant == "integrated"
                 else "C = 1 + C0*exp(-rho_local/rho_c)", density="rho = d(r Vbar^2)/dr / (4 pi G r^2), numpy.gradient native radii",
                 params=dict(C0=C0_V16, rho_c=RHOC_V16, Ydisk=Y, Ybul=YBUL, L0_kpc=(L0 if a.variant == "integrated" else None)),
                 L0_provenance=prov, data=dict(zip=ZIP, zip_sha256=ZIP_SHA, mrt_sha256=MRT_SHA, real_data_only=True),
                 split_seed=tp.SEED, holdout_revealed=False, holdout_evaluated=False,
                 chi2_red_definition="per galaxy mean of ((Vobs-V)/errV)^2 over unflagged points; same points for C1, model, RAR",
                 paired_metric="d_i = chi2_red(reference) - chi2_red(model); H1: d>0 (model better)",
                 stats=dict(summary, near_30p37_within_0p5kpc_flagged=len(near))))
        print(json.dumps({k: v for k, v in summary.items() if not k.startswith("near")}, indent=1))
        print("near 30.37 kpc (+/-0.5) flagged:", near); print("near 30.37 kpc negative:", near_neg)

if __name__ == "__main__":
    main()
