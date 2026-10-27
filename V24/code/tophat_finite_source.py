#!/usr/bin/env python3
"""
TAL mozgatható tömeg – órateszt: véges forrás top-hat integrál és érzékenység (STEP 1).

Publikált képletek (forrás: Time_Attraction.pdf, Zenodo 10.5281/zenodo.23066993):
  Def.1 Eq.(1), 3. old.:  rho_L(x) = 3/(4 pi L^3) * Int_{|x'-x|<=L} rho(x') d^3x'
  Def.2 Eq.(2), 3. old.:  gamma_L = G_N L^2 / c^2
  Eq.(13), 5. old.:       Delta f/f = gamma_L * Delta rho_L
  Def.3 Eq.(7), 4. old.:  f_L = exp(-gamma_L rho_L)  ->  pontos arány: f_near/f_far - 1 = exp(-gamma_L Delta rho_L) - 1
Egy homogén, rho_s sűrűségű forrás esetén Delta rho_L = rho_s * V_overlap(L) / (4/3 pi L^3),
ahol V_overlap az óra köré írt L sugarú gömb és a forrástest metszetének térfogata.
Az óra pontszerű (feltevés). A háttér (labor, Föld) a rögzített óránál mindkét állapotban
azonos, így a near/far különbségből kiesik (feltevés: a forrás által elhagyott tér vákuum/levegő
sűrűsége elhanyagolható – a levegőt nem modellezzük).
"""
import math, csv, json, os, sys
import numpy as np
from scipy.integrate import quad

G = 6.67430e-11          # Time_Attraction.pdf 3. old. (Def. 2 alatt)
C = 299792458.0          # Time_Attraction.pdf 3. old.
C2 = C * C
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CSV = os.path.join(OUT, "csv"); PLOTS = os.path.join(OUT, "plots")
os.makedirs(CSV, exist_ok=True); os.makedirs(PLOTS, exist_ok=True)

def gamma_L(L): return G * L * L / C2
def vball(L): return 4.0 / 3.0 * math.pi * L ** 3

# ---------- analitikus: gömb–gömb metszet ----------
def v_sphere_sphere(L, a, d):
    """L sugarú gömb (középpont: óra) és a sugarú gömb (középpont d távolságra) metszete."""
    if d >= L + a: return 0.0
    if d <= abs(L - a): return vball(min(L, a))
    return math.pi * (L + a - d) ** 2 * (d * d + 2 * d * a - 3 * a * a + 2 * d * L + 6 * a * L - 3 * L * L) / (12 * d)

def v_shell_centered(L, ri, ro):
    """Az órára centrált gömbhéj (ri..ro) és az L-gömb metszete."""
    return vball(min(L, ro)) - vball(min(L, ri))

# ---------- kvadratúra: gömb ∩ téglatest ----------
def _G(Y, Z, r):
    # {0<=u<=Y, 0<=v<=Z, u^2+v^2<=r^2} területe, Y,Z>=0
    if r <= 0 or Y <= 0 or Z <= 0: return 0.0
    Ym = min(Y, r)
    u0 = math.sqrt(max(r * r - Z * Z, 0.0))
    if u0 >= Ym: return Ym * Z
    F = lambda u: 0.5 * (u * math.sqrt(max(r * r - u * u, 0.0)) + r * r * math.asin(min(u / r, 1.0)))
    return Z * u0 + F(Ym) - F(u0)

def _S(Y, Z, r):
    return math.copysign(1.0, Y) * math.copysign(1.0, Z) * _G(abs(Y), abs(Z), r)

def rect_disk_area(y1, y2, z1, z2, r):
    return _S(y2, z2, r) - _S(y1, z2, r) - _S(y2, z1, r) + _S(y1, z1, r)

def v_ball_box(L, box):
    """Origó-centrált L sugarú gömb ∩ [x1,x2]x[y1,y2]x[z1,z2]; 1D kvadratúra x szerint,
    a keresztmetszet (téglalap ∩ kör) területe analitikus."""
    x1, x2, y1, y2, z1, z2 = box
    a, b = max(x1, -L), min(x2, L)
    if a >= b: return 0.0
    far = max(abs(x1), abs(x2)) ** 2 + max(abs(y1), abs(y2)) ** 2 + max(abs(z1), abs(z2)) ** 2
    if L * L >= far:   # a téglatest teljesen a gömbben: pontos érték
        return (x2 - x1) * (y2 - y1) * (z2 - z1)
    f = lambda x: rect_disk_area(y1, y2, z1, z2, math.sqrt(max(L * L - x * x, 0.0)))
    # törésponti helyek: ahol a kör sugara egyenlő a téglalap sarkai/élei távolságával
    pts = set()
    for yy in (0.0, y1, y2):
        for zz in (0.0, z1, z2):
            s = L * L - yy * yy - zz * zz
            if s > 0:
                for xx in (math.sqrt(s), -math.sqrt(s)):
                    if a < xx < b: pts.add(xx)
    val, err = quad(f, a, b, points=sorted(pts) or None, limit=400, epsabs=0, epsrel=1e-10)
    return val

# ---------- Monte Carlo ellenőrzés ----------
def mc_ball_fraction(L, inside_fn, n=2_000_000, seed=12345):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    r = L * rng.random(n) ** (1.0 / 3.0)
    p = v * r[:, None]
    k = inside_fn(p).sum()
    frac = k / n
    return frac * vball(L), math.sqrt(frac * (1 - frac) / n) * vball(L)

# ---------- forrásgeometriák a publikált protokollokból ----------
RHO_W, RHO_PB = 19250.0, 11340.0   # Time_Attraction.pdf 5. old.; v2 12oldal 7. old. / pitch: 11.34 g/cm3
M_W, A_W = 2406.0, 0.5             # Time_Attraction.pdf 5. old.: "wolfram kocka 0.5^3 m^3 ... M = 2406kg"
M_PB = 10000.0                     # 10 t – v2 (12oldal 9. old.; bővített 3. old.), v16 pitch 1. old.
V_PB = M_PB / RHO_PB
A_PB_CUBE = V_PB ** (1 / 3)        # SZCENÁRIÓ: alak nincs megadva a forrásban -> kocka
R_PB_SPH = (3 * V_PB / (4 * math.pi)) ** (1 / 3)  # SZCENÁRIÓ: gömb

def box_on_axis(d_center, wx, wy, wz):
    """Téglatest, középpontja az x tengelyen d_center távolságra, élei wx,wy,wz."""
    return (d_center - wx / 2, d_center + wx / 2, -wy / 2, wy / 2, -wz / 2, wz / 2)

def scenarios():
    S = []
    # v1–v7 W kocka: 0.5 m él, 0.5 m "táv" (középpont vagy felület – a forrás nem dönti el); "ki" = eltávolítva
    for lab, dc in (("center", 0.5), ("surface", 0.5 + A_W / 2)):
        S.append(dict(id=f"W_cube_v1v7_d0.5{lab}", src="Time_Attraction.pdf 5. old. (7. fejezet)",
                      rho=RHO_W, M=RHO_W * A_W ** 3, kind="box", near=box_on_axis(dc, A_W, A_W, A_W),
                      far=None, d_near_c=dc, d_far_c=math.inf,
                      note=f"0.5 m táv mint {'középpont' if lab=='center' else 'felület'}-óra távolság; távoli állapot = eltávolítva"))
    # 10 t ólom: közel/távol párok a forrásokból
    pairs = [("v2bov_0.5_2.0", 0.5, 2.0, "publikacio_bovitett.pdf 3. old.: r = 0.5–2 m"),
             ("v2resz_0.5_5.0", 0.5, 5.0, "Teljes_Reszletes.pdf 4. old.: 0.5 m (A) / 5 m (B)"),
             ("v16pitch_0.5_3.0", 0.5, 3.0, "Time_Attraction_LabxPitch.pdf 1. old.: NEAR 0.5 m / FAR 3.0 m")]
    for pid, dn, dfar, src in pairs:
        # kocka, középpont-olvasat (rés 0.5 - a/2 = 0.0205 m) és felület-olvasat
        for lab, off in (("center", 0.0), ("surface", A_PB_CUBE / 2)):
            S.append(dict(id=f"Pb10t_cube_{pid}_{lab}", src=src, rho=RHO_PB, M=M_PB, kind="box",
                          near=box_on_axis(dn + off, A_PB_CUBE, A_PB_CUBE, A_PB_CUBE),
                          far=box_on_axis(dfar + off, A_PB_CUBE, A_PB_CUBE, A_PB_CUBE),
                          d_near_c=dn + off, d_far_c=dfar + off,
                          note=f"SZCENÁRIÓ alak: kocka (él {A_PB_CUBE:.4f} m, a 10 t és 11340 kg/m3-ből); távolság mint {lab}"))
        # gömb, csak felület-olvasat (középpont-olvasatnál R=0.595 m > 0.5 m, az óra a testbe esne)
        S.append(dict(id=f"Pb10t_sphere_{pid}_surface", src=src, rho=RHO_PB, M=M_PB, kind="sphere",
                      a=R_PB_SPH, d_near_c=dn + R_PB_SPH, d_far_c=dfar + R_PB_SPH,
                      note=f"SZCENÁRIÓ alak: gömb (R {R_PB_SPH:.4f} m); távolság mint felület"))
    # torony-aspektus scan (a v2/v16 "torony" méretei nem találhatók): négyzetes alap w, magasság h=q*w
    for q in (2.0, 4.0):
        w = (V_PB / q) ** (1 / 3); h = q * w
        S.append(dict(id=f"Pb10t_tower_hw{int(q)}_v16pitch_surface", src="LabxPitch 1. old. távolságok; alak SZCENÁRIÓ",
                      rho=RHO_PB, M=M_PB, kind="box",
                      near=box_on_axis(0.5 + w / 2, w, w, h), far=box_on_axis(3.0 + w / 2, w, w, h),
                      d_near_c=0.5 + w / 2, d_far_c=3.0 + w / 2,
                      note=f"SZCENÁRIÓ: négyzetes alapú torony w={w:.4f} m, h={h:.4f} m (h/w={q}), óra a torony közepének magasságában"))
    # v2 gömbhéj a "ORA A" körül: L=1.5 m sugarú, 0.5 m vastag
    vol_shell = vball(1.5) - vball(1.0)
    S.append(dict(id="Pb_shell_v2_1.0-1.5_fullPb", src="Javitott_12oldal.pdf 9. old.", rho=RHO_PB,
                  M=RHO_PB * vol_shell, kind="shell", ri=1.0, ro=1.5, d_near_c=0.0, d_far_c=math.inf,
                  note=f"héj 1.0–1.5 m tömör ólomból: tömeg {RHO_PB*vol_shell/1000:.1f} t (≠ 10 t – a forrás belső ellentmondása)"))
    S.append(dict(id="Pb_shell_v2_1.0-1.5_10t_smeared", src="Javitott_12oldal.pdf 9. old.", rho=M_PB / vol_shell,
                  M=M_PB, kind="shell", ri=1.0, ro=1.5, d_near_c=0.0, d_far_c=math.inf,
                  note=f"10 t a 1.0–1.5 m héjba kenve: {M_PB/vol_shell:.1f} kg/m3 (nem ólomsűrűség)"))
    return S

def drho(sc, L):
    """Delta rho_L = rho_L(near) - rho_L(far) az óránál [kg/m3]."""
    if sc["kind"] == "box":
        vn = v_ball_box(L, sc["near"]); vf = v_ball_box(L, sc["far"]) if sc["far"] else 0.0
    elif sc["kind"] == "sphere":
        vn = v_sphere_sphere(L, sc["a"], sc["d_near_c"]); vf = v_sphere_sphere(L, sc["a"], sc["d_far_c"])
    elif sc["kind"] == "shell":
        vn = v_shell_centered(L, sc["ri"], sc["ro"]); vf = 0.0
    return sc["rho"] * (vn - vf) / vball(L)

# Óra-instabilitás szintek A FORRÁSOKBÓL
SIGMAS_SRC = {"1e-18 (Time_Attraction.pdf 3.,5. old.; v2 12oldal 9. old.)": 1e-18,
              "1.2e-18 Sr (LabxPitch 1. old.)": 1.2e-18}

def main():
    Ls = np.unique(np.concatenate([np.logspace(-3, 3, 241), [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 367.0]]))
    S = scenarios()
    rows, summ = [], []
    for sc in S:
        best = (0, None)
        for L in Ls:
            dr = drho(sc, L); g = gamma_L(L)
            y_lin = g * dr; y_exact = math.expm1(-g * dr)
            gmin1 = (1e-18 / dr) if dr > 0 else math.inf
            gmin2 = (2e-18 / dr) if dr > 0 else math.inf
            rows.append(dict(scenario=sc["id"], L_m=L, gamma_L_m3_per_kg=g, delta_rho_L_kg_m3=dr,
                             dff_lin_gammaL_drho=y_lin, dff_exact_expm1=y_exact,
                             gamma_min_k1_sigma1e18=gmin1, gamma_min_k2_sigma1e18=gmin2))
            if abs(y_lin) > abs(best[0]): best = (y_lin, L)
        dr1 = drho(sc, 1.0); dr367 = drho(sc, 367.0)
        drmax = max(abs(drho(sc, L)) for L in Ls)
        Ldrmax = Ls[int(np.argmax([abs(drho(sc, L)) for L in Ls]))]
        dn, dfr = sc["d_near_c"], sc["d_far_c"]
        if sc["kind"] == "shell":   # homogén héj potenciálja a középpontban: 2 pi G rho (ro^2 - ri^2) (pontos)
            newt = 2 * math.pi * G * sc["rho"] * (sc["ro"] ** 2 - sc["ri"] ** 2) / C2
        else:                       # pontszerű/gömb közelítés: G M (1/d_near - 1/d_far)/c^2 (gömbre pontos, kockára közelítő)
            newt = G * sc["M"] / C2 * (1 / dn - (0 if math.isinf(dfr) else 1 / dfr))
        summ.append(dict(scenario=sc["id"], source=sc["src"], note=sc["note"], rho_kg_m3=sc["rho"], M_kg=sc["M"],
                         d_near_center_m=dn, d_far_center_m=dfr,
                         drho_L1m=dr1, dff_L1m=gamma_L(1.0) * dr1, drho_L367m=dr367, dff_L367m=gamma_L(367.0) * dr367,
                         max_abs_drho=drmax, L_at_max_drho=Ldrmax, max_abs_dff=abs(best[0]), L_at_max_dff=best[1],
                         newton_GM_c2_1overr_diff=newt,
                         ratio_maxdff_to_1e18=abs(best[0]) / 1e-18, ratio_maxdff_to_newton=abs(best[0]) / newt))
    with open(os.path.join(CSV, "step1_tophat_scan.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(os.path.join(CSV, "step1_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)

    # ---- keresztellenőrzések ----
    checks = []
    a = A_PB_CUBE
    for L in (0.3, 0.6, 1.0, 1.5, 3.0):
        box = box_on_axis(0.5 + a / 2, a, a, a)
        q = v_ball_box(L, box)
        inside = lambda p, b=box: (p[:, 0] >= b[0]) & (p[:, 0] <= b[1]) & (np.abs(p[:, 1]) <= b[3]) & (np.abs(p[:, 2]) <= b[5])
        m, e = mc_ball_fraction(L, inside)
        checks.append(dict(test="ball∩cube kvadratúra vs MC", L_m=L, quad=q, mc=m, mc_1sigma=e, diff_sigma=(q - m) / e if e > 0 else 0))
    for L in (0.3, 0.8, 1.2, 2.0):
        R = R_PB_SPH; d = 0.5 + R
        an = v_sphere_sphere(L, R, d)
        inside = lambda p: ((p[:, 0] - d) ** 2 + p[:, 1] ** 2 + p[:, 2] ** 2) <= R * R
        m, e = mc_ball_fraction(L, inside)
        # gömb-ellenőrzés kvadratúrával is: forgásszimmetrikus integrál
        f = lambda x: math.pi * max(0.0, min(L * L - x * x, R * R - (x - d) ** 2)) if (L*L-x*x)>0 and (R*R-(x-d)**2)>0 else 0.0
        qq = quad(f, max(-L, d - R), min(L, d + R), limit=200)[0] if max(-L, d-R) < min(L, d+R) else 0.0
        checks.append(dict(test="gömb–gömb analitikus vs MC (és 1D kvadratúra)", L_m=L, quad=qq, mc=m, mc_1sigma=e,
                           analytic=an, diff_sigma=(an - m) / e if e > 0 else 0))
    # határesetek: L-gömb teljesen a kockában ill. kocka teljesen a gömbben
    box = (-1, 1, -1, 1, -1, 1)
    checks.append(dict(test="határeset: L=0.5 gömb 2 m-es kockában = 4/3 pi L^3", L_m=0.5, quad=v_ball_box(0.5, box), analytic=vball(0.5)))
    checks.append(dict(test="határeset: 2 m-es kocka L=2 gömbben = 8 m3", L_m=2.0, quad=v_ball_box(2.0, box), analytic=8.0))
    with open(os.path.join(CSV, "step1_crosschecks.csv"), "w", newline="") as fh:
        keys = sorted({k for c in checks for k in c}, key=lambda k: ["test","L_m","analytic","quad","mc","mc_1sigma","diff_sigma"].index(k))
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(checks)

    # ---- inverz kérdések a publikált állításokhoz ----
    inv = {}
    # (i) v16 pitch 5e-19: van-e L, ahol gamma_L*drho = 5e-19 (bármely alak-szcenárióra)?
    for sc in S:
        mx = max(gamma_L(L) * drho(sc, L) for L in np.logspace(-3, 3, 600))
        inv[sc["id"]] = {"max_dff_any_L": mx, "reaches_5e-19": bool(mx >= 5e-19), "reaches_1e-18": bool(mx >= 1e-18), "reaches_1e-19": bool(mx >= 1e-19)}
    json.dump(dict(A_PB_CUBE=A_PB_CUBE, R_PB_SPH=R_PB_SPH, V_PB=V_PB, inverse=inv), open(os.path.join(CSV, "step1_inverse_checks.json"), "w"), indent=1)
    return S, rows, summ, checks

if __name__ == "__main__":
    S, rows, summ, checks = main()
    for s in summ:
        print(f"{s['scenario']:42s} drho(1m)={s['drho_L1m']:.4g} dff(1m)={s['dff_L1m']:.3g} max|dff|={s['max_abs_dff']:.3g} @L={s['L_at_max_dff']:.3g} m  maxdrho={s['max_abs_drho']:.4g}@{s['L_at_max_drho']:.3g}  Newton={s['newton_GM_c2_1overr_diff']:.3g}")
    for c in checks: print(c)
