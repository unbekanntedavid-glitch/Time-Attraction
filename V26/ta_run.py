#!/usr/bin/env python3
"""
ta_run.py - simple standalone run (no dependency on ta_pipeline.py / ta_foundation.py).
Clean Python stand-in; NOT the user's app.tsk formula mode.

Inputs (real SPARC only, SHA256-asserted):
  /workspace/sparc-reference/Rotmod_LTG.zip   0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588
  data/SPARC_Lelli2016c.mrt                    5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243
Prints: zip contents summary, SPARC dwarf table (DDO*, T>=9, or 0<Vflat<60) with Rdisk,
        L0 = median Rdisk over (a) all 175, (b) dwarf set  [informational only],
        then runs NGC3198 with LOCKED L0 = 2.330 kpc:
          C=1, v16 local, v16 integrated, RAR(a0=1.2e-10), at Ydisk 0.5 and 1.0.
Density: rho = d(r Vbar^2)/dr / (4 pi G r^2) [kg/m^3], numpy.gradient(edge_order=1) on native radii;
         first/last radius flagged 'edge'; rho<0, Vbar^2<0 flagged; never clipped; flagged points excluded
         from chi2 and RMS of ALL models (identical point set).
Integrated: rho_eff = rho + (1/L0) * trapezoid integral of signed rho from 0 (constant below r1) to r.
chi2_red = mean(((Vobs-Vmod)/errV)^2) over unflagged points; RMS = sqrt(mean((Vmod-Vobs)^2)) [km/s].
"""
import hashlib, io, os, sys, zipfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ZIP = "/workspace/sparc-reference/Rotmod_LTG.zip"
MRT = os.path.join(HERE, "data", "SPARC_Lelli2016c.mrt")
ZIP_SHA = "0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588"
MRT_SHA = "5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243"
G = 6.67430e-11; KPC = 3.085677581e19
C0, RHO_C, YBUL, A0 = 5.0, 1e-23, 0.7, 1.2e-10
L0_LOCKED = 2.330
OUT = os.path.join(HERE, "results", "ta_run")

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

def read_zip():
    assert sha(ZIP) == ZIP_SHA, "zip hash mismatch"
    with zipfile.ZipFile(ZIP) as z:
        names = sorted(z.namelist())
        assert all(n.endswith("_rotmod.dat") for n in names) and len(names) == 175
        gals = {}
        for n in names:
            txt = z.read(n).decode(); assert txt.startswith("# Distance")
            gals[n.replace("_rotmod.dat", "")] = np.loadtxt(io.StringIO(txt), comments="#", ndmin=2)
    return gals

def read_mrt():
    assert sha(MRT) == MRT_SHA, "mrt hash mismatch"
    out = {}
    for line in open(MRT):
        p = line.split()
        if len(p) < 18: continue
        try:
            out[p[0]] = dict(T=int(p[1]), Inc=float(p[5]), Rdisk=float(p[11]), Vflat=float(p[15]), Q=int(p[17]))
        except ValueError:
            pass
    return out

def dwarfs(meta):
    d = {}
    for k, m in meta.items():
        why = [w for w, c in (("DDO", k.startswith("DDO")), ("T>=9", m["T"] >= 9), ("Vflat<60", 0 < m["Vflat"] < 60)) if c]
        if why: d[k] = why
    return d

def models(a, Ydisk, L0):
    R, Vobs, errV, Vgas, Vdisk, Vbul = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5]
    r = R * KPC
    Vb2 = np.sign(Vgas) * Vgas**2 + Ydisk * Vdisk**2 + YBUL * Vbul**2
    rho = np.gradient(r * Vb2 * 1e6, r, edge_order=1) / (4 * np.pi * G * r**2)
    acc = rho[0] * r[0]; cum = [acc]
    for i in range(1, len(R)):
        acc += 0.5 * (rho[i] + rho[i - 1]) * (r[i] - r[i - 1]); cum.append(acc)
    rho_int = np.array(cum) / (L0 * KPC)
    n = len(R); flags = []
    for i in range(n):
        f = []
        if i in (0, n - 1): f.append("edge")
        if rho[i] < 0: f.append("rho<0")
        if Vb2[i] < 0: f.append("Vbar2<0")
        flags.append(f)
    ok = np.array([not f for f in flags])
    with np.errstate(all="ignore"):
        gbar = Vb2 * 1e6 / r
        V = {"C=1": np.sqrt(Vb2),
             "v16_local": np.sqrt((1 + C0 * np.exp(-rho / RHO_C)) * Vb2),
             "v16_integrated": np.sqrt((1 + C0 * np.exp(-(rho + rho_int) / RHO_C)) * Vb2),
             "RAR": np.sqrt(gbar / (1 - np.exp(-np.sqrt(gbar / A0))) * r) / 1e3}
    return R, Vobs, errV, np.sqrt(np.abs(Vb2)), rho, rho_int, flags, ok, V

def main():
    os.makedirs(OUT, exist_ok=True)
    gals = read_zip(); meta = read_mrt()
    print(f"Rotmod_LTG.zip: {len(gals)} files, NGC3198 present: {'NGC3198' in gals}")
    assert set(gals) == set(meta)
    dw = dwarfs(meta)
    with open(os.path.join(OUT, "sparc_dwarfs_rdisk.csv"), "w") as f:
        f.write("galaxy,T,Vflat_kms,Q,Inc_deg,Rdisk_kpc,criteria\n")
        for k in sorted(dw):
            m = meta[k]; f.write(f"{k},{m['T']},{m['Vflat']},{m['Q']},{m['Inc']},{m['Rdisk']},{'+'.join(dw[k])}\n")
    print(f"SPARC dwarfs (DDO* or T>=9 or 0<Vflat<60): {len(dw)}; DDO: {sum(k.startswith('DDO') for k in dw)}")
    L0_all = float(np.median([m["Rdisk"] for m in meta.values()]))
    L0_dw = float(np.median([meta[k]["Rdisk"] for k in dw]))
    print(f"L0 (a) median Rdisk, all 175 unfiltered = {L0_all:.4f} kpc")
    print(f"L0 (b) median Rdisk, {len(dw)} SPARC dwarfs  = {L0_dw:.4f} kpc")
    print(f"L0 used (LOCKED, not fitted)              = {L0_LOCKED:.3f} kpc")
    lines = ["Ydisk,model,n_used,n_flagged,n_rho_neg,chi2_red,dV_30.08_kms,RMS_kms"]
    for Y in (0.5, 1.0):
        R, Vobs, errV, Vbar, rho, rho_int, flags, ok, V = models(gals["NGC3198"], Y, L0_LOCKED)
        i = int(np.where(np.isclose(R, 30.08))[0][0])
        nneg = sum("rho<0" in f for f in flags)
        print(f"\nNGC3198 Ydisk={Y}: Vbar(30.08)={Vbar[i]:.2f} km/s, rho(30.08)={rho[i]:.4e}, rho_int(30.08)={rho_int[i]:.4e} kg/m3; "
              f"flagged {int((~ok).sum())}/{len(R)}: " + ", ".join(f"{R[j]}:{'+'.join(flags[j])}({rho[j]:.3e})" for j in range(len(R)) if flags[j]))
        for k, v in V.items():
            c2 = float(np.mean(((Vobs[ok] - v[ok]) / errV[ok])**2)); rms = float(np.sqrt(np.mean((v[ok] - Vobs[ok])**2)))
            print(f"  {k:15s} chi2_red={c2:10.4f}  Vmod-Vobs @30.08 = {v[i]-Vobs[i]:+8.2f} km/s (Vmod={v[i]:.2f})  RMS={rms:7.2f} km/s")
            lines.append(f"{Y},{k},{int(ok.sum())},{int((~ok).sum())},{nneg},{c2:.6g},{v[i]-Vobs[i]:.6g},{rms:.6g}")
    open(os.path.join(OUT, "NGC3198_summary.csv"), "w").write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main()
