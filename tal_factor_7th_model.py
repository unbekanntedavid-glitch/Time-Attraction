#!/usr/bin/env python3
"""
TAL factor f_L = exp(-gamma_L * rho_bar) as a 7th model next to the six
baryonic baselines of the hash-locked v17 SPARC runner (sparc_full_runner.py).

Re-uses, unchanged, from the copied v17 runner: hash checks, Rotmod extraction,
blind-sample loading, parse_rotmod (incl. the R<0.01 -> 0.16 kpc fix),
sigma = max(0.05*Vobs, 20 km/s), chi2_red = chi2/(N-1), CSV layout.

Published inputs (see tal_factor_7th_model.md for citations):
  f_L(rho_bar_L) = exp(-gamma_L rho_bar_L),  gamma_L = G_N L^2 / c^2
  T_eff = f_L T  (G_mu_nu = 8 pi G_N f_L T_mu_nu),  TOV: dm/dr = 4 pi r^2 rho f_L
  L0 favoured = 1.2e-11 m; upper bound L < 367 m (gamma < 1e-22 m^3/kg)
  rho_disk(R,0) = Ydisk * SBdisk * Msun/pc^2 / (2 hz), hz = 0.2 kpc, Ydisk = 0.5, Ybulge = 0.7
Implementation choices are flagged in the report.
"""
import csv, math, statistics, sys, json, hashlib
from pathlib import Path
import importlib.util

SCRIPT_DIR = Path(__file__).resolve().parent

def _candidate_bases():
    """Script dir, then <script dir>/v17_package/github, then <cwd>/v17_package/github."""
    raw = [
        SCRIPT_DIR,
        SCRIPT_DIR / "v17_package" / "github",
        Path.cwd() / "v17_package" / "github",
    ]
    bases, seen = [], set()
    for base in raw:
        resolved = base.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        bases.append(resolved)
    return bases

def _resolve_inputs():
    bases = _candidate_bases()
    runner = None
    for base in bases:
        candidate = base / "sparc_full_runner.py"
        if candidate.is_file():
            runner = candidate
            break
    if runner is None:
        searched = "\n  ".join(str(b / "sparc_full_runner.py") for b in bases)
        sys.exit(f"ERROR: sparc_full_runner.py not found. Searched:\n  {searched}")
    data_root = None
    for base in bases:
        data = base / "data"
        if (data / "blind_123.txt").is_file() and (data / "Rotmod_LTG.zip").is_file():
            data_root = base
            break
    if data_root is None:
        searched = "\n  ".join(str(b / "data") for b in bases)
        sys.exit(
            "ERROR: data directory not found. sparc_full_runner.py expects "
            "data/blind_123.txt and data/Rotmod_LTG.zip. Searched:\n  " + searched
        )
    return runner, data_root

def _resolve_points_vbar(data_root, runner_path):
    """Baseline the gamma=0 check reads. The runner writes output/points_vbar.csv;
    the hash-locked copy in this repo is run_output/points_vbar.csv."""
    candidates, seen = [], set()
    for base in (data_root, runner_path.parent):
        for rel in (("output", "points_vbar.csv"), ("run_output", "points_vbar.csv")):
            path = base.joinpath(*rel)
            if path in seen:
                continue
            seen.add(path)
            candidates.append(path)
    for path in candidates:
        if path.is_file():
            return path
    searched = "\n  ".join(str(p) for p in candidates)
    sys.exit(f"ERROR: points_vbar.csv not found. Searched:\n  {searched}")

RUNNER_PATH, DATA_ROOT = _resolve_inputs()
spec = importlib.util.spec_from_file_location("runner", RUNNER_PATH)
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
# Runner reads blind_123.txt and Rotmod_LTG.zip from ROOT / "data".
R.ROOT = DATA_ROOT
R.DATA = DATA_ROOT / "data"
R.ROTMOD = R.DATA / "rotmod"
R.BLIND_FILE = R.DATA / "blind_123.txt"
R.OUTPUT = DATA_ROOT / "output"

G = 6.67430e-11; C = 299792458.0; C2 = C * C
MSUN = 1.98847e30; PC = 3.085677581e16; KPC = 1e3 * PC; MPC = 1e6 * PC
HZ = 0.2 * KPC; YDISK = 0.5; YBUL = 0.7

def gamma_of_L(L_m): return G * L_m**2 / C2
def L_of_gamma(g): return math.sqrt(g * C2 / G)

def parse_full(fp):
    """Same rows/filters as runner.parse_rotmod, plus SBdisk/SBbul columns."""
    base = R.parse_rotmod(fp)
    sb = []
    with open(fp) as f:
        for line in f:
            p = line.split()
            if len(p) < 6: continue
            try:
                float(p[0]); float(p[1]); float(p[3]); float(p[4]); float(p[5])
            except ValueError:
                continue
            sbd = float(p[6]) if len(p) > 6 else 0.0
            sbb = float(p[7]) if len(p) > 7 else 0.0
            sb.append((sbd, sbb))
    assert len(sb) == len(base), fp
    return [b + s for b, s in zip(base, sb)]

def rho_diskbulge(Rkpc, Vb, sbd, sbb):
    """Primary rho_bar proxy (disk midplane + bulge), kg/m^3."""
    rho_d = YDISK * sbd * MSUN / PC**2 / (2 * HZ)
    r = Rkpc * KPC
    rho_b = 3 * YBUL * (Vb * 1e3)**2 / (4 * math.pi * G * r**2)
    return rho_d + rho_b

def rho_enclosed(Rkpc, Vbar):
    """Alternative proxy: mean enclosed baryonic density 3 Vbar^2/(4 pi G R^2)."""
    r = Rkpc * KPC
    return 3 * (Vbar * 1e3)**2 / (4 * math.pi * G * r**2)

def fL_and_vpred(Vbar, gamma, rho, coupling):
    x = gamma * rho
    if coupling == "forward":            # V^2 = f_L Vbar^2
        f = math.exp(-x)
        return f, Vbar * math.sqrt(f)
    else:                                # V^2 = Vbar^2 / f_L
        f = math.exp(-x)
        try: return f, Vbar * math.exp(0.5 * x)
        except OverflowError: return f, float("inf")

def run(name, gamma, coupling, rhofn, blind, outdir, write=True):
    points, summary, missing = [], [], []
    for gal in sorted(blind):
        fp = R.ROTMOD / f"{gal}_rotmod.dat"
        if not fp.exists(): missing.append(gal); continue
        raw = parse_full(fp)
        if not raw: missing.append(gal); continue
        chi2 = 0.0; n = 0
        for Rk, Vobs, Vg, Vd, Vb, sbd, sbb in raw:
            Vbar = math.sqrt(Vg**2 + Vd**2 + Vb**2)
            rho = rhofn(Rk, Vbar, Vb, sbd, sbb)
            f, Vpred = fL_and_vpred(Vbar, gamma, rho, coupling)
            if Vobs > 0:
                s = R.sigma(Vobs); n += 1
                rr = (Vpred - Vobs) / s
                try: chi2 += rr**2          # identical arithmetic to the v17 runner
                except OverflowError: chi2 = float("inf")
                points.append({"model": name, "galaxy": gal, "R_kpc": Rk, "Vobs": Vobs,
                               "Vbar": Vbar, "Vpred": Vpred, "sigma": s,
                               "residual": rr,
                               "rho_bar_kg_m3": rho, "gammaL_rho": gamma * rho, "f_L": f})
        if n > 1:
            summary.append({"model": name, "galaxy": gal, "n_radii": n, "chi2_red": chi2 / (n - 1)})
    if missing: sys.exit(f"missing galaxies [{name}]: {missing}")
    if write:
        with open(outdir / f"points_{name}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(points[0].keys())); w.writeheader(); w.writerows(points)
        with open(outdir / f"summary_{name}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["model","galaxy","n_radii","chi2_red"]); w.writeheader(); w.writerows(summary)
    return points, summary

RHO_PRIMARY = lambda Rk, Vbar, Vb, sbd, sbb: rho_diskbulge(Rk, Vb, sbd, sbb)
RHO_ENC = lambda Rk, Vbar, Vb, sbd, sbb: rho_enclosed(Rk, Vbar)

def main():
    vbar_csv = _resolve_points_vbar(DATA_ROOT, RUNNER_PATH)
    out = SCRIPT_DIR / "output_7th"
    print(f"runner: {RUNNER_PATH}")
    print(f"data: {R.DATA}")
    print(f"blind_123.txt: {R.BLIND_FILE}")
    print(f"Rotmod_LTG.zip: {R.DATA / 'Rotmod_LTG.zip'}")
    print(f"points_vbar.csv: {vbar_csv}")
    print(f"output: {out}")
    R.verify_rotmod_archive()            # hash checks + extraction + 2596-point check
    blind = R.load_blind()
    out.mkdir(exist_ok=True)
    res = {}
    # sanity: gamma=0 must reproduce points_vbar.csv Vpred exactly
    pts0, _ = run("check_gamma0", 0.0, "forward", RHO_PRIMARY, blind, out, write=False)
    ref = list(csv.DictReader(open(vbar_csv)))
    assert len(ref) == len(pts0) == 2596, f"point count mismatch: ref={len(ref)} pts0={len(pts0)} (expected 2596)"
    same = all(str(p["Vpred"]) == r["Vpred"] and p["galaxy"] == r["galaxy"] for p, r in zip(pts0, ref)) and len(ref) == len(pts0)
    print(f"gamma=0 reproduces points_vbar.csv Vpred: {same} (N={len(pts0)})")
    if not same: sys.exit("ERROR: gamma=0 run does not reproduce points_vbar.csv")
    res["gamma0_check"] = same

    # rho statistics (primary and alternative proxies)
    rhos = [p["rho_bar_kg_m3"] for p in pts0]
    rho_enc = []
    for gal in sorted(blind):
        for Rk, Vobs, Vg, Vd, Vb, sbd, sbb in parse_full(R.ROTMOD / f"{gal}_rotmod.dat"):
            if Vobs > 0: rho_enc.append(rho_enclosed(Rk, math.sqrt(Vg**2+Vd**2+Vb**2)))
    def q(v, p):
        v = sorted(v); return v[min(len(v)-1, int(p*(len(v)-1)+0.5))]
    for lab, v in (("primary_disk_bulge", rhos), ("enclosed_baryonic", rho_enc)):
        pos = [x for x in v if x > 0]
        Lu = [C / math.sqrt(G * x) / MPC for x in pos]   # L with gamma_L rho = 1
        res[f"rho_{lab}"] = {"min": min(v), "median": statistics.median(v), "max": max(v),
                             "n_zero": sum(1 for x in v if x <= 0),
                             "L_unity_Mpc_median": statistics.median(Lu),
                             "L_unity_Mpc_p05": q(Lu, .05), "L_unity_Mpc_p95": q(Lu, .95)}
    print(json.dumps(res, indent=1))

    models = [
        # name, gamma, coupling, rho, label
        ("tal_fL", gamma_of_L(1.2e-11), "forward", RHO_PRIMARY, "PRIMARY: L0=1.2e-11 m (favoured), V^2=f_L Vbar^2"),
        ("tal_fL_L367m", gamma_of_L(367.0), "forward", RHO_PRIMARY, "upper bound L=367 m"),
        ("tal_fL_g1e-22", 1e-22, "forward", RHO_PRIMARY, "upper bound gamma=1e-22"),
        ("tal_fL_L0_2p36kpc", gamma_of_L(2.36 * KPC), "forward", RHO_PRIMARY, "v17 L0=2.36 kpc"),
        ("tal_fL_L0_2p36kpc_inv", gamma_of_L(2.36 * KPC), "inverse", RHO_PRIMARY, "v17 L0=2.36 kpc, inverse"),
        ("tal_fL_L64Mpc", gamma_of_L(64.2 * MPC), "forward", RHO_PRIMARY, "L=64.2 Mpc fwd"),
        ("tal_fL_L64Mpc_inv", gamma_of_L(64.2 * MPC), "inverse", RHO_PRIMARY, "L=64.2 Mpc inv"),
        ("tal_fL_L85Mpc", gamma_of_L(85.3 * MPC), "forward", RHO_PRIMARY, "L=85.3 Mpc fwd"),
        ("tal_fL_L85Mpc_inv", gamma_of_L(85.3 * MPC), "inverse", RHO_PRIMARY, "L=85.3 Mpc inv"),
        ("tal_fL_g1e23", 1e23, "forward", RHO_PRIMARY, "gamma=1e23 fwd"),
        ("tal_fL_g1e23_inv", 1e23, "inverse", RHO_PRIMARY, "gamma=1e23 inv"),
        ("tal_fL_enc_L85Mpc_inv", gamma_of_L(85.3 * MPC), "inverse", RHO_ENC, "enclosed rho, L=85.3 Mpc inv"),
        ("tal_fL_enc_L64Mpc_inv", gamma_of_L(64.2 * MPC), "inverse", RHO_ENC, "enclosed rho, L=64.2 Mpc inv"),
    ]
    # data-derived L_unity (median over sample) for both proxies
    for lab, rf, tag in (("primary_disk_bulge", RHO_PRIMARY, "Lunity"), ("enclosed_baryonic", RHO_ENC, "enc_Lunity")):
        Lu = res[f"rho_{lab}"]["L_unity_Mpc_median"]
        for cpl in ("forward", "inverse"):
            models.append((f"tal_fL_{tag}_{'inv' if cpl=='inverse' else 'fwd'}", gamma_of_L(Lu * MPC), cpl, rf,
                           f"data-derived median L_unity={Lu:.2f} Mpc ({lab}), {cpl}"))
    summaries = {}
    for name, g, cpl, rf, lab in models:
        pts, summ = run(name, g, cpl, rf, blind, out, write=True)
        dev = max(abs(p["f_L"] - 1.0) for p in pts)
        gr = max(p["gammaL_rho"] for p in pts)
        n_ident = sum(1 for p, r in zip(pts, ref) if str(p["Vpred"]) == r["Vpred"])
        med = statistics.median(s["chi2_red"] for s in summ)
        summaries[name] = {"label": lab, "gamma_m3_per_kg": g, "L_m": L_of_gamma(g), "coupling": cpl,
                           "max_abs_fL_minus_1": dev, "max_gammaL_rho": gr,
                           "points_with_Vpred_identical_to_vbar": n_ident, "n_points": len(pts),
                           "n_gal": len(summ), "median_chi2_red": med}
        print(f"{name:28s} gamma={g:.4e} max|f-1|={dev:.3e} max(g*rho)={gr:.3e} identical={n_ident}/{len(pts)} median={med:.4f}")
    res["models"] = summaries
    json.dump(res, open(out / "run_meta.json", "w"), indent=1)

    # fitted global L (NOT blind): grid scan, inverse and forward couplings, primary rho
    scan = []
    for cpl in ("inverse", "forward"):
        for k in range(0, 161):
            LMpc = 10 ** (-1 + k * 0.025)    # 0.1 Mpc .. 1e3 Mpc
            _, summ = run("scan", gamma_of_L(LMpc * MPC), cpl, RHO_PRIMARY, blind, out, write=False)
            v = [s["chi2_red"] for s in summ]
            scan.append({"coupling": cpl, "L_Mpc": LMpc, "median_chi2_red": statistics.median(v),
                         "mean_chi2_red": statistics.fmean(v)})
    with open(out / "fit_scan_L.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(scan[0].keys())); w.writeheader(); w.writerows(scan)
    for cpl in ("inverse", "forward"):
        best = min((s for s in scan if s["coupling"] == cpl), key=lambda s: s["median_chi2_red"])
        print("best fitted", best)
        name = f"tal_fL_FITTED_{cpl}"
        run(name, gamma_of_L(best["L_Mpc"] * MPC), cpl, RHO_PRIMARY, blind, out, write=True)
        res.setdefault("fitted", {})[cpl] = best
    json.dump(res, open(out / "run_meta.json", "w"), indent=1)

if __name__ == "__main__":
    main()
