#!/usr/bin/env python3
"""
ta_pipeline.py - Time-Attraction v27 WORK PIPELINE (clean Python stand-in)

IMPORTANT: This is NOT a replication of the user's app.tsk "formula mode"
TypeScript source (that source was never provided). It is an independent,
clean Python stand-in built from the v16.0 PDF definitions and public SPARC data.

Data : SPARC (Lelli, McGaugh & Schombert 2016) Rotmod_LTG + SPARC_Lelli2016c.mrt
Grid : locked to SPARC radii (no interpolation, no resampling).
Density (external, independent of rotation speeds; v16 Eq5/Eq6 photometric proxy):
   Sigma_disk = Ydisk * SBdisk                      [Msun/pc^2], Ydisk = 0.5
   Sigma_bul  = Ybul  * SBbul                       [Msun/pc^2], Ybul  = 0.7  (bulge treated with the same
                                                      thin-layer proxy; documented simplification)
   Sigma_gas  = 1.33 * M_HI * exp(-R/Rg) / (2 pi Rg^2), Rg = 2*Rdisk  (v16 Eq6; M_HI, Rdisk from mrt)
   rho        = (Sigma_disk + Sigma_bul + Sigma_gas) / (2 hz),  hz = 0.2 kpc      [kg/m^3]
   (local proxy only: v16's nonlocal rho_integrated term is NOT included; --gas none drops Eq6)
Baryons: Vbar^2 = sign(Vgas) Vgas^2 + Ydisk Vdisk^2 + Ybul Vbul^2   (signed gas)
Model : C(rho) user formula or v16 C = 1 + 5 exp(-rho/1e-23); Vmod = sqrt(C) * Vbar
        (accel mode: expression gives g_mod [m/s^2]; Vmod = sqrt(g_mod * R))
chi2  : per point ((Vobs - Vmod)/errV)^2; per galaxy chi2_g/N_g ("chi2_red", no dof subtraction).
Bad points (rho<=0, NaN, Vbar^2<0, errV<=0, non-finite model) are COUNTED, LISTED, and EXCLUDED
from chi2 - never clipped.
"""
import argparse, ast, json, math, os, sys, glob, datetime, hashlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROT_DIR = os.path.join(HERE, "data", "Rotmod_LTG")
MRT = os.path.join(HERE, "data", "SPARC_Lelli2016c.mrt")
LOG = os.path.join(HERE, "attempts_log.jsonl")

MSUN = 1.98847e30; PC = 3.085677581e16; KPC = 1e3 * PC
YDISK, YBUL, HZ_KPC = 0.5, 0.7, 0.2
SEED, FRAC = 27, (0.6, 0.2, 0.2)          # train / revealed-test / holdout

# ---------------------------------------------------------------- data
def load_mrt(path=MRT):
    """Whitespace-split parse (the data rows do not follow the byte ranges in the header exactly).
    Columns: Galaxy T D e_D f_D Inc e_Inc L36 e_L36 Reff SBeff Rdisk SBdisk MHI RHI Vflat e_Vflat Q [Ref]"""
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) < 18:
            continue
        try:
            rec = dict(T=int(p[1]), D=float(p[2]), Inc=float(p[5]), L36=float(p[7]), Rdisk=float(p[11]),
                       SBdisk0=float(p[12]), MHI=float(p[13]), RHI=float(p[14]), Vflat=float(p[15]),
                       e_Vflat=float(p[16]), Q=int(p[17]))
        except ValueError:
            continue
        out[p[0]] = rec
    return out

def load_galaxies(gas="v16eq6"):
    meta = load_mrt()
    gals = {}
    for f in sorted(glob.glob(os.path.join(ROT_DIR, "*_rotmod.dat"))):
        name = os.path.basename(f).replace("_rotmod.dat", "")
        a = np.loadtxt(f, comments="#")
        R, Vobs, errV, Vgas, Vdisk, Vbul, SBd, SBb = a.T
        m = meta[name]
        Sd = YDISK * SBd; Sb = YBUL * SBb
        if gas == "v16eq6" and m["MHI"] > 0 and m["Rdisk"] > 0:
            Rg = 2.0 * m["Rdisk"]                                   # kpc
            Sg = 1.33 * m["MHI"] * 1e9 * np.exp(-R / Rg) / (2 * np.pi * Rg**2) / 1e6   # Msun/pc^2
        else:
            Sg = np.zeros_like(R)
        Sigma = Sd + Sb + Sg                                          # Msun/pc^2
        rho = Sigma * MSUN / PC**2 / (2 * HZ_KPC * KPC)               # kg/m^3
        Vbar2 = np.sign(Vgas) * Vgas**2 + YDISK * Vdisk**2 + YBUL * Vbul**2
        gals[name] = dict(name=name, R=R, Vobs=Vobs, errV=errV, Vgas=Vgas, Vdisk=Vdisk, Vbul=Vbul,
                          SBdisk=SBd, SBbul=SBb, Sigma=Sigma, Sigma_gas=Sg, rho=rho, Vbar2=Vbar2, meta=m)
    return gals

def split(names, seed=SEED):
    names = sorted(names); rng = np.random.default_rng(seed); p = rng.permutation(len(names))
    n1 = int(round(FRAC[0] * len(names))); n2 = n1 + int(round(FRAC[1] * len(names)))
    lab = {}
    for i, k in enumerate(p):
        lab[names[k]] = "train" if i < n1 else ("test" if i < n2 else "holdout")
    return lab

# ---------------------------------------------------------------- safe expression
ALLOWED_FUNCS = {"exp": np.exp, "log": np.log, "log10": np.log10, "sqrt": np.sqrt, "abs": np.abs,
                 "tanh": np.tanh, "minimum": np.minimum, "maximum": np.maximum, "where": np.where,
                 "sin": np.sin, "cos": np.cos, "pow": np.power}
ALLOWED_NODES = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Call, ast.Name, ast.Load, ast.Constant,
                 ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd, ast.Mod,
                 ast.Compare, ast.Lt, ast.Gt, ast.LtE, ast.GtE, ast.IfExp)
VARS = ["rho", "R_kpc", "R_m", "Vbar", "Vbar2", "gbar", "Sigma", "rho_c", "C0", "pi", "G"]

def safe_compile(expr, extra_names=()):
    tree = ast.parse(expr, mode="eval")
    for node in ast.walk(tree):
        if not isinstance(node, ALLOWED_NODES):
            raise ValueError(f"disallowed syntax: {type(node).__name__}")
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in ALLOWED_FUNCS:
                raise ValueError(f"disallowed function: {ast.dump(node.func)}")
            if node.keywords:
                raise ValueError("keyword arguments not allowed")
        if isinstance(node, ast.Name) and node.id not in ALLOWED_FUNCS and node.id not in VARS \
                and node.id not in extra_names:
            raise ValueError(f"unknown name: {node.id}")
        if isinstance(node, ast.Constant) and not isinstance(node.value, (int, float)):
            raise ValueError("only numeric constants allowed")
    return compile(tree, "<formula>", "eval")

def evaluate(code, env):
    ns = {"__builtins__": {}}; ns.update(ALLOWED_FUNCS); ns.update(env)
    with np.errstate(all="ignore"):
        return eval(code, ns)

# ---------------------------------------------------------------- model + stats
def model_galaxy(g, mode, code, params):
    R = g["R"]; Vbar2 = g["Vbar2"]
    Vbar = np.sign(Vbar2) * np.sqrt(np.abs(Vbar2))
    env = dict(rho=g["rho"], R_kpc=R, R_m=R * KPC, Vbar=Vbar, Vbar2=Vbar2,
               gbar=Vbar2 * 1e6 / (R * KPC), Sigma=g["Sigma"], pi=math.pi, G=6.67430e-11,
               rho_c=1e-23, C0=5.0)
    env.update(params)
    val = np.broadcast_to(np.asarray(evaluate(code, env), dtype=float), R.shape)
    with np.errstate(all="ignore"):
        if mode == "C":
            Vmod = np.sqrt(val * Vbar2)
        else:  # accel
            Vmod = np.sqrt(val * R * KPC) / 1e3
    return val, Vmod

def flags_for(g, Vmod):
    bad = {}
    bad["rho<=0"] = ~(g["rho"] > 0)
    bad["rho_nan"] = ~np.isfinite(g["rho"])
    bad["Vbar2<0"] = g["Vbar2"] < 0
    bad["errV<=0"] = ~(g["errV"] > 0)
    bad["model_nonfinite"] = ~np.isfinite(Vmod)
    return bad

def run_attempt(gals, labels, mode, expr, params, reveal_holdout=False, note="", log=True, quiet=False):
    try:
        code = safe_compile(expr, extra_names=tuple(params))
    except (ValueError, SyntaxError) as e:
        rec = dict(timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(), mode=mode,
                   expression=expr, params=params, gas_model=GAS_MODEL, rejected=True, error=str(e), note=note)
        if log:
            with open(LOG, "a") as f: f.write(json.dumps(rec) + "\n")
        print("REJECTED:", e, file=sys.stderr); sys.exit(2)
    per = {}; flagged = []
    for name, g in gals.items():
        val, Vmod = model_galaxy(g, mode, code, params)
        bad = flags_for(g, Vmod)
        anybad = np.zeros_like(g["R"], dtype=bool)
        for k, b in bad.items():
            anybad |= b
            for i in np.where(b)[0]:
                flagged.append(dict(galaxy=name, R_kpc=float(g["R"][i]), flag=k,
                                    rho=float(g["rho"][i]), Vbar2=float(g["Vbar2"][i])))
        ok = ~anybad
        res = (g["Vobs"][ok] - Vmod[ok]) / g["errV"][ok]
        per[name] = dict(split=labels[name], n=int(ok.sum()), n_bad=int(anybad.sum()),
                         chi2=float(np.sum(res**2)),
                         sq=float(np.sum((g["Vobs"][ok] - Vmod[ok])**2)))
    stats = {}
    for s in ("all", "train", "test", "train+test", "holdout"):
        sel = [v for v in per.values() if s == "all" or v["split"] == s or (s == "train+test" and v["split"] in ("train", "test"))]
        sel_n = [v for v in sel if v["n"] > 0]
        N = sum(v["n"] for v in sel); X = sum(v["chi2"] for v in sel)
        st = dict(n_gal=len(sel), n_points=N, n_bad_points=sum(v["n_bad"] for v in sel),
                  chi2=X, chi2_per_point=X / N if N else None,
                  median_chi2_red=float(np.median([v["chi2"] / v["n"] for v in sel_n])) if sel_n else None,
                  rmse_kms=math.sqrt(sum(v["sq"] for v in sel) / N) if N else None)
        stats[s] = st
    visible = {k: v for k, v in stats.items() if k in ("train", "test", "train+test") or reveal_holdout}
    if not reveal_holdout:
        visible["holdout"] = "HIDDEN (use --reveal-holdout)"
        visible["all"] = "HIDDEN (includes holdout; use --reveal-holdout)"
    rec = dict(timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(), mode=mode,
               expression=expr, params=params, gas_model=GAS_MODEL, split_seed=SEED,
               split_fracs=FRAC, holdout_revealed=reveal_holdout, note=note,
               n_flag_entries=len(flagged),
               n_flagged_points_visible=sum(1 for x in set((d["galaxy"], d["R_kpc"]) for d in flagged) if labels[x[0]] != "holdout" or reveal_holdout), stats=visible)
    if log:
        with open(LOG, "a") as f:
            f.write(json.dumps(rec) + "\n")
    if not quiet:
        print(json.dumps(rec, indent=2))
    return rec, per, flagged, stats


# ---------------------------------------------------------------- fitting (train set only)
def fit_C0_rhoc(gals, labels, starts=((5.0, 1e-23), (1.0, 1e-22), (10.0, 1e-24), (3.0, 3e-23))):
    from scipy.optimize import minimize
    code = safe_compile("1 + C0*exp(-rho/rho_c)", extra_names=("C0", "rho_c"))
    train = [g for n, g in gals.items() if labels[n] == "train"]
    def chi2(theta):
        C0, rc = 10**theta[0], 10**theta[1]
        X = 0.0
        for g in train:
            _, Vmod = model_galaxy(g, "C", code, {"C0": C0, "rho_c": rc})
            bad = np.zeros_like(g["R"], dtype=bool)
            for b in flags_for(g, Vmod).values(): bad |= b
            ok = ~bad
            X += np.sum(((g["Vobs"][ok] - Vmod[ok]) / g["errV"][ok])**2)
        return X
    best = None; tried = []
    for C0, rc in starts:
        r = minimize(chi2, [math.log10(C0), math.log10(rc)], method="Nelder-Mead",
                     options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
        tried.append(dict(start=[C0, rc], C0=10**r.x[0], rho_c=10**r.x[1], chi2=float(r.fun), ok=bool(r.success)))
        if best is None or r.fun < best.fun: best = r
    # crude 1-sigma (delta chi2 = 1) from numerical Hessian in log space
    h = 1e-3; x0 = best.x; H = np.zeros((2, 2)); f0 = best.fun
    for i in range(2):
        for j in range(2):
            ei = np.eye(2)[i] * h; ej = np.eye(2)[j] * h
            H[i, j] = (chi2(x0 + ei + ej) - chi2(x0 + ei - ej) - chi2(x0 - ei + ej) + chi2(x0 - ei - ej)) / (4 * h * h)
    try:
        cov = np.linalg.inv(H / 2.0); err = np.sqrt(np.diag(cov)).tolist()
    except np.linalg.LinAlgError:
        err = [float("nan")] * 2
    return dict(C0=10**best.x[0], rho_c=10**best.x[1], log10_C0=best.x[0], log10_rho_c=best.x[1],
                sigma_log10=err, chi2_train=float(best.fun), starts=tried)

# ---------------------------------------------------------------- NGC 3198 table
def ngc3198_table(gals, out_csv, out_md):
    g = gals["NGC3198"]
    C = 1 + 5 * np.exp(-g["rho"] / 1e-23)
    Vb05 = np.sqrt(g["Vbar2"])
    Vb10 = np.sqrt(np.sign(g["Vgas"]) * g["Vgas"]**2 + 1.0 * g["Vdisk"]**2 + YBUL * g["Vbul"]**2)
    Vmod = np.sqrt(C) * Vb05
    rho_d = (YDISK * g["SBdisk"] + YBUL * g["SBbul"]) * MSUN / PC**2 / (2 * HZ_KPC * KPC)
    C_d = 1 + 5 * np.exp(-rho_d / 1e-23); Vmod_d = np.sqrt(C_d) * Vb05
    rows = []
    for i in range(len(g["R"])):
        rows.append([g["R"][i], g["Vobs"][i], g["errV"][i], g["Vgas"][i], g["Vdisk"][i], g["SBdisk"][i],
                     g["Sigma_gas"][i], g["rho"][i], C[i], Vb05[i], Vb10[i], Vmod[i],
                     (g["Vobs"][i] - Vmod[i]) / g["errV"][i], rho_d[i], C_d[i], Vmod_d[i]])
    hdr = ["R_kpc", "Vobs", "errV", "Vgas", "Vdisk", "SBdisk_Lpc2", "Sigma_gas_Msunpc2", "rho_kgm3", "C_v16",
           "Vbar_Y0.5", "Vbar_Y1.0", "Vmod_v16", "pull", "rho_diskonly", "C_diskonly", "Vmod_diskonly"]
    with open(out_csv, "w") as f:
        f.write(",".join(hdr) + "\n")
        for r in rows: f.write(",".join(f"{x:.6g}" for x in r) + "\n")
    with open(out_md, "w") as f:
        f.write("# NGC 3198 - SPARC Rotmod_LTG, all radii (locked grid)\n\n")
        f.write("Vbar(Y=0.5) = sqrt(sign(Vgas)Vgas^2 + 0.5 Vdisk^2 + 0.7 Vbul^2); Vbar(Y=1.0) same with Ydisk=1.0 "
                "(Vbul=0 for NGC3198). rho = v16 photometric proxy incl. v16 Eq6 gas. C_v16 = 1+5exp(-rho/1e-23). "
                "Vmod = sqrt(C_v16)*Vbar(Y=0.5). pull = (Vobs-Vmod)/errV. *_diskonly: same without the Eq6 gas term "
                "(rho=0 where SBdisk rounds to 0.00 -> C=6). Velocities km/s, rho kg/m^3. "
                "Python stand-in, not the app.tsk formula mode.\n\n")
        f.write("| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr) + "\n")
        for r in rows:
            mark = "**" if abs(r[0] - 30.08) < 1e-6 else ""
            f.write("| " + " | ".join(f"{mark}{x:.4g}{mark}" for x in r) + " |\n")
    return hdr, rows

# ---------------------------------------------------------------- dwarfs
def dwarf_list(gals):
    out = []
    for n, g in gals.items():
        m = g["meta"]
        why = []
        if n.upper().startswith("DDO"): why.append("DDO")
        if m["T"] >= 9: why.append(f"T={m['T']}")
        if 0 < m["Vflat"] < 60: why.append(f"Vflat={m['Vflat']}")
        if why:
            out.append(dict(galaxy=n, T=m["T"], Vflat=m["Vflat"], Q=m["Q"], n_points=len(g["R"]), why=why))
    return sorted(out, key=lambda d: d["galaxy"])

GAS_MODEL = "v16eq6"

def main():
    global GAS_MODEL
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gas", choices=["v16eq6", "none"], default="v16eq6")
    ap.add_argument("--reveal-holdout", action="store_true")
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("try", help="evaluate a formula on the locked SPARC grid and log it")
    t.add_argument("--mode", choices=["C", "accel"], default="C")
    t.add_argument("--expr", required=True, help="e.g. '1 + C0*exp(-rho/rho_c)'")
    t.add_argument("--param", action="append", default=[], help="name=value (repeatable)")
    t.add_argument("--note", default="")
    sub.add_parser("baselines", help="log C=1 and v16 C(rho) attempts")
    sub.add_parser("split", help="print the deterministic split")
    sub.add_parser("fit", help="fit C0, rho_c of 1+C0*exp(-rho/rho_c) on TRAIN only; logs an attempt")
    sub.add_parser("ngc3198", help="write results/NGC3198_table.{csv,md}")
    sub.add_parser("dwarfs", help="list SPARC dwarfs (DDO, T>=9, Vflat<60)")
    sub.add_parser("flags", help="write results/flagged_points.csv for the v16 baseline")
    a = ap.parse_args()
    GAS_MODEL = a.gas
    gals = load_galaxies(gas=a.gas); labels = split(gals)
    if a.cmd == "split":
        for s in ("train", "test"):
            print(s, sorted(k for k, v in labels.items() if v == s))
        print("holdout:", sum(v == "holdout" for v in labels.values()), "galaxies (names hidden unless --reveal-holdout)")
        if a.reveal_holdout:
            print("holdout", sorted(k for k, v in labels.items() if v == "holdout"))
    elif a.cmd == "baselines":
        run_attempt(gals, labels, "C", "1", {}, a.reveal_holdout, note="baseline: baryons only C=1")
        run_attempt(gals, labels, "C", "1 + C0*exp(-rho/rho_c)", {"C0": 5.0, "rho_c": 1e-23},
                    a.reveal_holdout, note="baseline: v16 C(rho)")
    elif a.cmd == "fit":
        r = fit_C0_rhoc(gals, labels)
        print(json.dumps(r, indent=2))
        run_attempt(gals, labels, "C", "1 + C0*exp(-rho/rho_c)", {"C0": r["C0"], "rho_c": r["rho_c"]},
                    a.reveal_holdout, note=f"fit on train only (Nelder-Mead, log-space); sigma_log10={r['sigma_log10']}; chi2_train={r['chi2_train']:.3f}")
    elif a.cmd == "ngc3198":
        os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
        hdr, rows = ngc3198_table(gals, os.path.join(HERE, "results", "NGC3198_table.csv"),
                                  os.path.join(HERE, "results", "NGC3198_table.md"))
        for r in rows:
            if abs(r[0] - 30.08) < 1e-6: print(dict(zip(hdr, r)))
    elif a.cmd == "dwarfs":
        for d in dwarf_list(gals): print(d)
    elif a.cmd == "flags":
        _, _, flagged, _ = run_attempt(gals, labels, "C", "1 + C0*exp(-rho/rho_c)", {}, a.reveal_holdout, log=False, quiet=True)
        os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
        with open(os.path.join(HERE, "results", f"flagged_points_gas-{a.gas}.csv"), "w") as f:
            f.write("galaxy,R_kpc,flag,rho_kgm3,Vbar2\n")
            for x in flagged:
                gname = x["galaxy"] if (labels[x["galaxy"]] != "holdout" or a.reveal_holdout) else "[holdout-hidden]"
                f.write(f"{gname},{x['R_kpc']},{x['flag']},{x['rho']:.6g},{x['Vbar2']:.6g}\n")
        print(len(flagged), "flag entries (a point can carry >1 flag)")
    elif a.cmd == "try":
        params = {}
        for p in a.param:
            k, v = p.split("="); params[k.strip()] = float(v)
        run_attempt(gals, labels, a.mode, a.expr, params, a.reveal_holdout, note=a.note)

if __name__ == "__main__":
    main()
