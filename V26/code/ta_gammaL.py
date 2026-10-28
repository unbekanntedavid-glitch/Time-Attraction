#!/usr/bin/env python3
"""
ta_gammaL.py - NEW script (2026-10-19). Clean Python stand-in, NOT the user's app.tsk.
Tests the ORIGINAL Time-Attraction form with gamma_L = G L^2 / c^2 on REAL SPARC data (non-holdout only).

Reuses (imports, does not modify) ta_foundation.py: verified real Rotmod_LTG.zip + mrt (SHA256 asserted),
seed-27 split, blind set (Q in {1,2}, Inc>30, >=5 radii; N=114), density rho = d(r Vbar^2)/dr/(4 pi G r^2)
("dinamikai rho"), flags (edge, rho<0, ... excluded identically for every model), paired tests.

MODELS (all with the same unflagged points):
  C1   : V^2 = Vbar^2                                         ("tiszta Vbar")
  smooth: f = 1/(1 + gamma_L rho)  ->  V^2 = Vbar^2 / f = Vbar^2 (1 + gamma_L rho)
  exp2 : V^2 = Vbar^2 exp(2 gamma_L rho)   (form quoted in Javitott_12oldal p10 l163; extra, labelled)
  RAR  : McGaugh+2016, a0 = 1.2e-10 m/s^2
  gamma_L(L) = G L^2 / c^2 with L = r (local radius) | Rdisk (catalogue, per galaxy) | L0 in {0.83,1.47,1.55,2.33} kpc
  Global-L scan: L on a log grid, chosen by min summed chi2 on TRAIN only, then evaluated on blind/all visible.
SCALE ANALYSIS: per unflagged point L_unity = c / sqrt(G rho) (gamma_L rho = 1) and
  L_req = c sqrt(((Vobs/Vbar)^2 - 1) / (G rho)) for Vobs > Vbar (scale that would make the smooth form fit exactly).
Holdout is never evaluated. Log: results/gamma_L/attempts_log_gammaL.jsonl (new file; attempts_log.jsonl untouched).
"""
import json, os, datetime
import numpy as np
import ta_pipeline as tp
import ta_foundation as tf

G = 6.67430e-11; c = 299792458.0; KPC = tf.KPC; MPC = 1e3 * KPC
OUT = os.path.join(tp.HERE, "results", "gamma_L"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "attempts_log_gammaL.jsonl")
L0S = {"0.83_user_supplied": 0.83, "1.47_median73dwarfs": 1.47, "1.55_median59nonholdoutdwarfs": 1.55, "2.33_locked": 2.33}

def log(rec):
    rec["timestamp"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with open(LOG, "a") as f: f.write(json.dumps(rec) + "\n")

def chi2(Vobs, err, V): return float(np.mean(((Vobs - V) / err) ** 2))

def main():
    gals = tf.load_real_sparc(); labels = tp.split(gals); blind = set(tf.blind_set(gals, labels))
    visible = sorted(k for k in gals if labels[k] != "holdout")
    assert len(blind) == 114 and len(visible) == 140
    summary = {}
    for Y in (0.5, 1.0):
        pts = {}
        for k in visible:
            res = tf.compute(gals[k], Y, 2.33, "local"); ok = res["ok"]
            if ok.sum() == 0: continue
            pts[k] = dict(r=res["R"][ok] * KPC, rho=res["rho"][ok], Vb2=res["Vbar2"][ok], Vobs=res["Vobs"][ok],
                          err=res["errV"][ok], Vrar=res["Vrar"][ok], Rd=gals[k]["meta"]["Rdisk"] * KPC,
                          Rlast=gals[k]["R"][-1] * KPC, split=labels[k])
        def model_V(p, form, Lspec):
            L = p["r"] if Lspec == "r" else (p["Rd"] if Lspec == "Rdisk" else Lspec * KPC)
            x = G * L ** 2 / c ** 2 * p["rho"]
            return np.sqrt(p["Vb2"] * ((1 + x) if form == "smooth" else np.exp(2 * x))), x
        def evaluate(form, Lspec):
            per = {}; xmax = 0.0; xmed = []
            for k, p in pts.items():
                V, x = model_V(p, form, Lspec); xmax = max(xmax, float(x.max())); xmed.append(float(np.median(x)))
                per[k] = (chi2(p["Vobs"], p["err"], np.sqrt(p["Vb2"])), chi2(p["Vobs"], p["err"], V), chi2(p["Vobs"], p["err"], p["Vrar"]))
            out = dict(max_gammaL_rho=xmax, median_of_galaxy_median_gammaL_rho=float(np.median(xmed)))
            for sname, keys in (("blind", [k for k in per if k in blind]), ("all_visible", list(per))):
                a = np.array([per[k] for k in keys])
                out[sname] = dict(N_gal=len(keys), median_chi2red=dict(C1=float(np.median(a[:, 0])), model=float(np.median(a[:, 1])), RAR=float(np.median(a[:, 2]))),
                                  model_vs_C1=tf.paired_tests(a[:, 0] - a[:, 1]), model_vs_RAR=tf.paired_tests(a[:, 2] - a[:, 1]))
            return out
        resY = {}
        for form in ("smooth", "exp2"):
            for lab, Lspec in [("L=r", "r"), ("L=Rdisk", "Rdisk")] + [(f"L0={k}", v) for k, v in L0S.items()]:
                resY[f"{form}|{lab}"] = evaluate(form, Lspec)
                log(dict(mode="gammaL", form=("V^2=Vbar^2(1+gamma_L rho) [f=1/(1+gamma_L rho)]" if form == "smooth" else "V^2=Vbar^2 exp(2 gamma_L rho)"),
                         gammaL="G L^2/c^2", L=lab, Ydisk=Y, density="rho=d(r Vbar^2)/dr/(4 pi G r^2) local", data_sha=dict(zip=tf.ZIP_SHA, mrt=tf.MRT_SHA),
                         split_seed=tp.SEED, holdout_evaluated=False, fitted=False, stats=resY[f"{form}|{lab}"]))
        # global-L scan on TRAIN only (smooth form)
        grid = np.logspace(0, 12, 241)   # kpc: 1 kpc .. 1e12 kpc
        tr = [k for k in pts if pts[k]["split"] == "train"]
        tot = []
        for Lk in grid:
            s = 0.0
            for k in tr:
                p = pts[k]; V, _ = model_V(p, "smooth", Lk); s += float(np.sum(((p["Vobs"] - V) / p["err"]) ** 2))
            tot.append(s)
        Lbest = float(grid[int(np.argmin(tot))])
        resY["smooth|L_globalfit_train"] = evaluate("smooth", Lbest); resY["smooth|L_globalfit_train"]["L_best_kpc"] = Lbest
        resY["smooth|L_globalfit_train"]["L_best_Mpc"] = Lbest / 1e3; resY["smooth|L_globalfit_train"]["grid"] = "logspace(0,12,241) kpc, train only"
        log(dict(mode="gammaL_globalL_fit", form="V^2=Vbar^2(1+G L^2 rho/c^2)", fit_set="train", L_best_kpc=Lbest, Ydisk=Y,
                 holdout_evaluated=False, fitted=True, stats=resY["smooth|L_globalfit_train"]))
        # scale analysis
        Lu, Lreq, Lreq_r, Lreq_Rlast, frac_need = [], [], [], [], []
        for k, p in pts.items():
            Lu += list(c / np.sqrt(G * p["rho"]))
            need = (p["Vobs"] ** 2 / p["Vb2"] - 1)
            m = need > 0; frac_need.append(m.mean())
            lr = c * np.sqrt(need[m] / (G * p["rho"][m]))
            Lreq += list(lr); Lreq_r += list(lr / p["r"][m]); Lreq_Rlast += list(lr / p["Rlast"])
        sc = dict(N_points=len(Lu), L_unity_median_Mpc=float(np.median(Lu) / MPC), L_unity_p05_p95_Mpc=[float(np.percentile(Lu, 5) / MPC), float(np.percentile(Lu, 95) / MPC)],
                  N_points_Vobs_gt_Vbar=len(Lreq), L_req_median_Mpc=float(np.median(Lreq) / MPC),
                  L_req_p05_p95_Mpc=[float(np.percentile(Lreq, 5) / MPC), float(np.percentile(Lreq, 95) / MPC)],
                  L_req_over_r_median=float(np.median(Lreq_r)), L_req_over_Rlast_median=float(np.median(Lreq_Rlast)),
                  Rlast_median_kpc=float(np.median([p["Rlast"] / KPC for p in pts.values()])))
        resY["scale_analysis"] = sc
        log(dict(mode="gammaL_scale_analysis", Ydisk=Y, holdout_evaluated=False, stats=sc))
        summary[f"Y{Y}"] = resY
    json.dump(summary, open(os.path.join(OUT, "summary_gammaL.json"), "w"), indent=1)
    rows = ["Ydisk,model,set,N,med_C1,med_model,med_RAR,vsC1_wilcoxon_p,vsC1_sign,vsC1_median_d,vsRAR_wilcoxon_p,vsRAR_sign,max_gammaL_rho"]
    for Yk, resY in summary.items():
        for mk, r in resY.items():
            if mk == "scale_analysis": continue
            for s in ("blind", "all_visible"):
                t = r[s]; a, b = t["model_vs_C1"], t["model_vs_RAR"]
                rows.append(f"{Yk[1:]},{mk},{s},{t['N_gal']},{t['median_chi2red']['C1']:.4g},{t['median_chi2red']['model']:.6g},{t['median_chi2red']['RAR']:.4g},"
                            f"{a['wilcoxon_p_one_sided']},{a['sign_n_pos']}/{a['sign_n_neg']},{a['median']:.3g},{b['wilcoxon_p_one_sided']},{b['sign_n_pos']}/{b['sign_n_neg']},{r['max_gammaL_rho']:.3g}")
    open(os.path.join(OUT, "summary_gammaL.csv"), "w").write("\n".join(rows) + "\n")
    print("\n".join(rows)); print(json.dumps({k: v["scale_analysis"] for k, v in summary.items()}, indent=1))
    for Yk, resY in summary.items(): print(Yk, "L_best_kpc", resY["smooth|L_globalfit_train"]["L_best_kpc"])

if __name__ == "__main__": main()
