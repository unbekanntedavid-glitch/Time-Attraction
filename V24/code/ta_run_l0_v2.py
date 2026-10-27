#!/usr/bin/env python3
"""
ta_run_l0_v2.py - L0 sensitivity run, v2: SUPERSEDES ta_run_l0.py.
Fix: additionally flags rho_eff = rho + rho_int <= 0 (as ta_foundation.py does). Without it, UGC01281
(negative inner Vbar^2 -> negative integral) got C up to 1.65e16 at L0=0.83 in ta_run_l0.py.
Log goes to results/ta_run_l0_v2/attempts_log_l0_v2.jsonl (new file) - not appended to attempts_log.jsonl.

(original header follows) ta_run_l0.py - L0 sensitivity run (new file; reuses ta_run.py functions read-only, no edits to existing files).
Clean Python stand-in; NOT the user's app.tsk formula mode.

L0 values (only the 'integrated' model depends on L0; C=1, local, RAR do not):
  0.83 kpc  USER-SUPPLIED: stated median Rdisk of David's 23 dwarfs (names NOT provided; unverified here)
  2.330 kpc locked value = median catalogue Rdisk (all 175 SPARC, and seed-27 blind set)
  1.470 kpc median catalogue Rdisk of the 73 SPARC dwarfs (DDO* or T>=9 or 0<Vflat<60)
Sets: NGC3198; SPARC dwarfs (non-holdout only); all non-holdout galaxies. Holdout (seed-27 split, 35 gal)
is never evaluated. Split re-implemented and asserted identical to ta_pipeline.split.
Per galaxy chi2_red = mean(((Vobs-V)/errV)^2) over unflagged points (edge, rho<0, Vbar2<0), identical
points for all models. Paired d_i = chi2_red(ref) - chi2_red(model); one-sided Wilcoxon signed-rank
(zeros dropped) and sign test, H1: model better (d>0).
"""
import json, os
import numpy as np
from scipy.stats import wilcoxon, binomtest
import ta_run as tr

L0S = [("0.83_user_supplied_23dwarfs", 0.83), ("2.330_locked_median_all175", 2.330), ("1.470_median_73_SPARC_dwarfs", 1.470)]
OUT = os.path.join(tr.HERE, "results", "ta_run_l0_v2")

def split(names, seed=27, frac=(0.6, 0.2, 0.2)):
    names = sorted(names); p = np.random.default_rng(seed).permutation(len(names))
    n1 = int(round(frac[0] * len(names))); n2 = n1 + int(round(frac[1] * len(names)))
    return {names[k]: ("train" if i < n1 else "test" if i < n2 else "holdout") for i, k in enumerate(p)}

def ptest(d):
    d = np.asarray(d, float); nz = d[d != 0]
    pos, neg = int((nz > 0).sum()), int((nz < 0).sum())
    w = wilcoxon(nz, alternative="greater").pvalue if len(nz) else None
    s = binomtest(pos, pos + neg, 0.5, alternative="greater").pvalue if pos + neg else None
    return dict(n=len(d), median_d=float(np.median(d)), wilcoxon_p=None if w is None else float(w),
                sign_pos=pos, sign_neg=neg, sign_p=None if s is None else float(s))

def main():
    os.makedirs(OUT, exist_ok=True)
    gals = tr.read_zip(); meta = tr.read_mrt(); dw = tr.dwarfs(meta)
    lab = split(gals)
    try:
        import ta_pipeline as tp; assert tp.split(gals) == lab; split_check = "identical to ta_pipeline.split"
    except Exception as e:
        split_check = f"could not verify: {e}"
    vis = sorted(k for k in gals if lab[k] != "holdout")
    dvis = [k for k in vis if k in dw]
    sub_med = float(np.median([meta[k]["Rdisk"] for k in dvis]))
    print(f"split check: {split_check}; non-holdout galaxies {len(vis)}, non-holdout SPARC dwarfs {len(dvis)} of {len(dw)}")
    print(f"median Rdisk: all 73 SPARC dwarfs 1.470 kpc; non-holdout dwarf subset ({len(dvis)}) {sub_med:.3f} kpc")
    sets = {"NGC3198": ["NGC3198"], "dwarfs_nonholdout": dvis, "all_nonholdout": vis}
    allres = dict(split_check=split_check, n_vis=len(vis), n_dwarf_vis=len(dvis), dwarf_vis=dvis,
                  median_Rdisk_dwarf_vis=sub_med, runs=[])
    rows = ["L0_label,L0_kpc,Ydisk,set,N_gal,n_rho_neg,n_gal_rho_neg,n_rho_eff_le0_extra,med_C1,med_local,med_integrated,med_RAR,"
            "local_vs_C1_wilcoxon_p,local_vs_C1_sign,local_vs_C1_sign_p,integ_vs_C1_wilcoxon_p,integ_vs_C1_sign,integ_vs_C1_sign_p,"
            "local_vs_RAR_wilcoxon_p,integ_vs_RAR_wilcoxon_p"]
    for lbl, L0 in L0S:
        for Y in (0.5, 1.0):
            per = {}
            for k in vis:
                R, Vobs, errV, Vbar, rho, rho_int, flags, ok, V = tr.models(gals[k], Y, L0)
                reff_bad = (rho + rho_int) <= 0
                for j in np.where(reff_bad)[0]:
                    if "rho<0" not in flags[j]: flags[j].append("rho_eff<=0")
                ok = ok & ~reff_bad
                neff = int(sum("rho_eff<=0" in f for f in flags))
                if ok.sum() == 0: continue
                c = {m: float(np.mean(((Vobs[ok] - v[ok]) / errV[ok])**2)) for m, v in V.items()}
                c["nneg"] = sum("rho<0" in f for f in flags); c["neff"] = neff; per[k] = c
            for sname, keys in sets.items():
                ks = [k for k in keys if k in per]
                med = {m: float(np.median([per[k][m] for k in ks])) for m in ("C=1", "v16_local", "v16_integrated", "RAR")}
                t = dict(local_vs_C1=ptest([per[k]["C=1"] - per[k]["v16_local"] for k in ks]),
                         integ_vs_C1=ptest([per[k]["C=1"] - per[k]["v16_integrated"] for k in ks]),
                         local_vs_RAR=ptest([per[k]["RAR"] - per[k]["v16_local"] for k in ks]),
                         integ_vs_RAR=ptest([per[k]["RAR"] - per[k]["v16_integrated"] for k in ks]))
                nneg = sum(per[k]["nneg"] for k in ks); ng = sum(per[k]["nneg"] > 0 for k in ks)
                neff_t = sum(per[k]["neff"] for k in ks); neff_g = [k for k in ks if per[k]["neff"] > 0]
                rec = dict(L0_label=lbl, L0_kpc=L0, Ydisk=Y, set=sname, N_gal=len(ks), n_rho_neg=nneg, n_gal_rho_neg=ng, n_rho_eff_le0_extra=neff_t, gal_rho_eff_le0=neff_g,
                           median_chi2_red=med, tests=t)
                if sname == "NGC3198": rec["NGC3198_chi2_red"] = per["NGC3198"]
                allres["runs"].append(rec)
                f = lambda x: "NA" if x is None else f"{x:.3g}"
                rows.append(f"{lbl},{L0},{Y},{sname},{len(ks)},{nneg},{ng},{neff_t},{med['C=1']:.4g},{med['v16_local']:.4g},"
                            f"{med['v16_integrated']:.4g},{med['RAR']:.4g},{f(t['local_vs_C1']['wilcoxon_p'])},"
                            f"{t['local_vs_C1']['sign_pos']}/{t['local_vs_C1']['sign_neg']},{f(t['local_vs_C1']['sign_p'])},"
                            f"{f(t['integ_vs_C1']['wilcoxon_p'])},{t['integ_vs_C1']['sign_pos']}/{t['integ_vs_C1']['sign_neg']},"
                            f"{f(t['integ_vs_C1']['sign_p'])},{f(t['local_vs_RAR']['wilcoxon_p'])},{f(t['integ_vs_RAR']['wilcoxon_p'])}")
    open(os.path.join(OUT, "summary.csv"), "w").write("\n".join(rows) + "\n")
    json.dump(allres, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    with open(os.path.join(OUT, "attempts_log_l0_v2.jsonl"), "a") as fh:
        import datetime
        for r in allres["runs"]:
            fh.write(json.dumps(dict(timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(), mode="ta_run_l0_v2",
                     expression="C=1 / v16 local / v16 integrated (rho_eff=rho+(1/L0)int rho ds) / RAR a0=1.2e-10",
                     params=dict(C0=5.0, rho_c=1e-23, Ydisk=r["Ydisk"], L0_kpc=r["L0_kpc"]),
                     L0_provenance=r["L0_label"] + (" (user-supplied; 23 dwarf names not provided)" if r["L0_kpc"] == 0.83 else ""),
                     holdout_revealed=False, holdout_evaluated=False, set=r["set"], stats=r)) + "\n")
    print("\n".join(rows))

if __name__ == "__main__":
    main()
