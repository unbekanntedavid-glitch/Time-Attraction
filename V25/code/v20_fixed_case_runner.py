# -*- coding: utf-8 -*-
"""EGYETLEN RÖGZÍTETT ESET (2026-10-21, akadémiai lektori javításokkal). V19: https://doi.org/10.5281/zenodo.23497128
[SZERZŐ] = a szerző (Szabolcsi D.) által megadott mért paraméter, 2026-10-21.
Rögzített eset: R = 1789 Ω, I = 14 mA/tekercs, 24 V BVS, párhuzamos, 8 tekercs, azonos áramirány, függőleges tengely;
tengelykör 0.45 m (a 300 mm-es tekercstest belső széle 0.30 m-en); menetátmérő 291.7 mm (mért, tolómérő, 5 leolvasás) a térmodellben;
óra középen, 300 mm-rel a tekercsek alsó síkja felett; referencia-OCXO 1.60 m-re, azonos magasságban, gránit nélkül.
Gránit: hengerhéj r = 0.09–0.22 m, h = 0.50 m, 2700 kg/m³. Kapszula: henger, belső r 0.05 m, külső r 0.08 m, h 0.24 m.
MINDEN Δf/f kimenet a DIFFERENCIÁLIS (fő − referencia) érték, előjellel: df_f_diff_signed = expm1(−γ·(ρ̄_fő − ρ̄_ref)) ≈ −γ·Δρ̄."""
import csv, json, math, os
import numpy as np
import fields as F, tophat_finite_source as T1
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); CSV = os.path.join(OUT, "csv"); PL = os.path.join(OUT, "plots")
G, C = T1.G, T1.C; C2 = C * C; MU0 = F.MU0; DAY = 86400.0; YEAR = 365.25 * DAY
def g4(x): return f"{x:.4g}"
def w(n, rows, h):
    with open(os.path.join(CSV, n), "w", newline="", encoding="utf-8") as f: c = csv.writer(f); c.writerow(h); c.writerows(rows)
S = {}
# ---------------- bemenetek
N, ell = 2000, 0.400; dw = 0.15e-3; A_w = math.pi * dw ** 2 / 4
R = 1789.0; I0 = 0.014; U_bat = 24.0; NC = 8; Rax = 0.45; hz = 0.200 - 0.300
RHO = 1.7241e-8; DENS = 8890.0          # IEC 60028
YIELD = 68.9e6; T_MELT = 1083.0; GAMMA_B = 1e-22
Dm = 0.2917                             # [SZERZŐ] mért menetátmérő (tolómérő, 0.02 mm felbontás, 5 × 291.7 mm, 2025-11-13)
a = Dm / 2; D_body = 0.300               # 300 mm = tekercstest pereme (külső méret)
L_wire = N * math.pi * Dm               # huzalhossz a mért geometriából
L_wire_R = R * A_w / RHO; D_R = L_wire_R / (N * math.pi)   # összevetés: R-ből (IEC ρ) származtatott hossz és átmérő
m_cu8 = 8 * DENS * A_w * L_wire
S.update(L_wire=L_wire, L_wire_R=L_wire_R, D_R=D_R, R_geom=RHO * L_wire / A_w, Dm=Dm, a=a, I_24V=U_bat / R, U_14mA=I0 * R, m_cu8=m_cu8)
w("derived_from_R1789.csv", [["huzalhossz/tekercs a mért menetátmérőből (N·π·D)", g4(L_wire), "m"], ["R20 a mért geometriából (ρ IEC)", g4(RHO * L_wire / A_w), "Ω"],
   ["huzalhossz R = 1789 Ω-ból (ρ IEC)", g4(L_wire_R), "m"], ["menetátmérő R-ből", f"{D_R:.5f}", "m"],
   ["24 V / 1789 Ω", g4(U_bat / R * 1e3), "mA"], ["14 mA × 1789 Ω", g4(I0 * R), "V"], ["Cu tömeg, 8 tekercs (IEC 8890 kg/m³)", f"{m_cu8:.4f}", "kg"]],
  ["mennyiség", "érték", "egység"])
# ---------------- geometria + tér
cs = np.array([[Rax * math.cos(2 * math.pi * i / NC), Rax * math.sin(2 * math.pi * i / NC), hz] for i in range(NC)])
SG = [1] * NC
pref = np.array([1.60, 0.0, 0.0])      # [SZERZŐ] 1 m a tekercstest külső felületétől (0.60 m) → 1.60 m
def Bn(p): return float(np.linalg.norm(F.coils_B(np.array([p], float), cs, a, ell, N, SG)[0]))
Bpa = Bn([0, 0, 0]); Bref_pa = Bn(pref); Bin_pa = MU0 * N / math.sqrt(ell ** 2 + Dm ** 2)
S.update(Bpa=Bpa, Bref_pa=Bref_pa, Bin_pa=Bin_pa, d_ref=1.60, d_ref_winding=float(min(np.hypot(pref[0] - c[0], pref[1] - c[1]) - a for c in cs)),
         gap_winding=float(2 * Rax * math.sin(math.pi / NC) - Dm), gap_body=float(2 * Rax * math.sin(math.pi / NC) - D_body))
# ---------------- differenciális térenergia: energy_integrator.py (GL-panelek, közös szögpontok, r = L/t távoli tartomány)
from energy_integrator import Integrator
from multiprocessing import Pool
Lg = np.geomspace(1e-2, 1e2, 81)
LCHK = (0.1, 0.3162, 0.631, 1.0, 2.0, 2.6, 3.162, 10.0, 31.62, 100.0)
VARIANTS = {"NHI24000 (alap)": dict(nhi=24000), "NHI6000": dict(nhi=6000), "NHI12000": dict(nhi=12000), "NHI12000, radiális panelek ×2": dict(nhi=12000, href=2)}
def run_variant(key):
    I = Integrator(cs, a, ell, N, pref, **VARIANTS[key])
    out = dict(chk=[I.dW(L) for L in LCHK], Wt=[I.w_tot(0), I.w_tot(1)], jump=[I.dW_inner(0.95), I.dW_A(0.95)])
    if key.startswith("NHI24000"):
        out["dW"] = [I.dW(L) for L in Lg]
        out["cont"] = [(L, I.dW_A(L), I.dW_B(L)) for L in (2.6, 2.7, 2.8, 2.9, 2.99, 3.0, 3.01)]
    return key, out
if __name__ == "__main__":
    with Pool(len(VARIANTS)) as pool: RES = dict(pool.map(run_variant, list(VARIANTS)))
    BASE = RES["NHI24000 (alap)"]; dWL = np.array(BASE["dW"])
    w("continuity_test.csv", [[g4(L), f"{A:.12e}", f"{B:.12e}", f"{(A - B) / B:.3e}", "folytonossági teszt" if L in (2.99, 3.0, 3.01) else "átfedési ellenőrzés 2.6–3 m"] for L, A, B in BASE["cont"]],
      ["L (m)", "ΔW A út: GL-panelek + r=R0/t farok (J/A²)", "ΔW B út: r=L/t, közvetlen külső különbség (J/A²)", "relatív eltérés (A−B)/B", "szerep"])
    keys = list(VARIANTS)
    conv = []
    for n_, L in enumerate(LCHK):
        vals = np.array([RES[k]["chk"][n_] for k in keys])
        conv.append([g4(L)] + [f"{v:.6g}" for v in vals] + [f"{vals.min():.6g}–{vals.max():.6g}", f"{vals.std(ddof=1):.3g}", f"{vals.std(ddof=1) / abs(vals.mean()):.3g}"])
    w("convergence_dW.csv", conv, ["L (m)"] + [f"ΔW(L) /I² (J/A²), {k}" for k in keys] + ["min-max tartomány (min–max)", "mintaszórás s (n−1, 4 futás) (J/A²)", "s / átlag"])
    S["W_total"] = {k: RES[k]["Wt"] for k in keys}; S["jump_LJ"] = {k: RES[k]["jump"] for k in keys}
    S["conv"] = {str(L): [RES[k]["chk"][n_] for k in keys] for n_, L in enumerate(LCHK)}
    S["cont"] = BASE["cont"]
    # független kereszt-ellenőrzés: közös térfogati Monte Carlo minták (a két indikátorfüggvény különbsége)
    rng = np.random.default_rng(7); mc = []
    for L in (0.3162, 0.631, 1.0, 2.0):
        Rs = L + np.linalg.norm(pref); n = 1_500_000
        u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1)[:, None]; X = u * Rs * rng.uniform(0, 1, n)[:, None] ** (1 / 3)
        dens = np.sum(F.coils_B(X, cs, a, ell, N * 1.0, SG) ** 2, -1) / (2 * MU0)
        ind = (np.linalg.norm(X, axis=1) < L).astype(float) - (np.linalg.norm(X - pref, axis=1) < L).astype(float)
        Vs = 4 / 3 * np.pi * Rs ** 3; val = Vs * np.mean(dens * ind); err = Vs * np.std(dens * ind) / math.sqrt(n)
        mc.append([L, val, err, float(np.interp(L, LCHK, BASE["chk"])) if L in LCHK else float(np.interp(L, Lg, dWL))])
    w("convergence_mc_shared_samples.csv", [[g4(L), g4(v), g4(e), g4(e / v), g4(s), g4((s - v) / e)] for L, v, e, s in mc],
      ["L (m)", "ΔW közös MC minták (J/A²)", "MC 1σ (J/A²)", "MC relatív hiba σ/ΔW", "ΔW GL-integrátor (alap)", "eltérés / σ"])
    S["mc"] = mc
    drho = dWL / C2 / T1.vball(Lg)
    df_g = np.expm1(-G * Lg ** 2 / C2 * drho * I0 ** 2); df_b = np.expm1(-GAMMA_B * drho * I0 ** 2)
    w("tal_differential_vs_L.csv", [[g4(L), g4(x), g4(y), g4(z), g4(u_), g4(v_)] for L, x, y, z, u_, v_ in zip(Lg, dWL, drho, df_g / I0 ** 2, df_g, df_b)],
      ["L (m)", "ΔW(L)/I² (J/A²)", "Δρ̄_L/I² fő−ref (kg/m³/A²)", "df_f_diff_signed γ_L / I² (1/A²)", "df_f_diff_signed γ_L, 14 mA", "df_f_diff_signed γ=1e-22, 14 mA"])
    ig = int(np.argmax(np.abs(df_g))); ib = int(np.argmax(np.abs(df_b)))
    kg = float(-G * Lg[ig] ** 2 / C2 * drho[ig]); kb = float(-GAMMA_B * drho[ib])
    S["TAL_diff"] = dict(k_g=kg, Lg=float(Lg[ig]), k_b=kb, Lb=float(Lg[ib]), neg_L=[float(L) for L, d in zip(Lg, dWL) if d < 0])
    # ---------------- statikus: réztekercsek (vékony héj r = Dm/2), közös minták
    rng = np.random.default_rng(1); M = 400000
    ph = rng.uniform(0, 2 * np.pi, M); zz = rng.uniform(-ell / 2, ell / 2, M); k_ = rng.integers(0, NC, M)
    P = np.stack([cs[k_, 0] + a * np.cos(ph), cs[k_, 1] + a * np.sin(ph), cs[k_, 2] + zz], -1)
    d0 = np.linalg.norm(P, axis=1); d1 = np.linalg.norm(P - pref, axis=1)
    cu = [[g4(L), g4(np.mean(d0 < L)), g4(np.mean(d1 < L)), g4(m_cu8 * (np.mean(d0 < L) - np.mean(d1 < L)) / T1.vball(L)),
           g4(-G * L * L / C2 * m_cu8 * (np.mean(d0 < L) - np.mean(d1 < L)) / T1.vball(L)), g4(-GAMMA_B * m_cu8 * (np.mean(d0 < L) - np.mean(d1 < L)) / T1.vball(L))] for L in Lg]
    w("copper_coils_static_differential.csv", cu, ["L (m)", "tömegtört a fő óra L-gömbjében", "u.a. referencia", "Δρ̄ fő−ref (kg/m³)", "df_f_diff_signed γ_L (statikus)", "df_f_diff_signed γ=1e-22 (statikus)"])
    cg = [float(r[4]) for r in cu]; cb = [float(r[5]) for r in cu]
    S["cu_static"] = dict(m=m_cu8, GR_main=float(-G * m_cu8 * np.mean(1 / d0) / C2), GR_ref=float(-G * m_cu8 * np.mean(1 / d1) / C2),
        ext_g=min(cg), ext_b=min(cb), dmin=float(d0.min()), dmax=float(d0.max()))
    # ---------------- statikus: gránit (csak a fő óra körül)
    RG_IN, RG_OUT, HG, RHO_G = 0.09, 0.22, 0.50, 2700.0
    MG = RHO_G * math.pi * (RG_OUT ** 2 - RG_IN ** 2) * HG
    Mg = 1000000; rr = np.sqrt(rng.uniform(RG_IN ** 2, RG_OUT ** 2, Mg)); pp = rng.uniform(0, 2 * np.pi, Mg); zg = rng.uniform(-HG / 2, HG / 2, Mg)
    Pg = np.stack([rr * np.cos(pp), rr * np.sin(pp), zg], -1); e0 = np.linalg.norm(Pg, axis=1); e1 = np.linalg.norm(Pg - pref, axis=1)
    gr = []
    for L in Lg:
        dr = MG * (np.mean(e0 < L) - np.mean(e1 < L)) / T1.vball(L)
        gr.append([g4(L), g4(np.mean(e0 < L)), g4(np.mean(e1 < L)), g4(dr), g4(-G * L * L / C2 * dr), g4(-GAMMA_B * dr)])
    w("granite_static_differential.csv", gr, ["L (m)", "gránit tömegtört a fő óra L-gömbjében", "u.a. referencia", "Δρ̄ fő−ref (kg/m³)", "df_f_diff_signed γ_L (statikus)", "df_f_diff_signed γ=1e-22 (statikus)"])
    gg = [float(r[4]) for r in gr]; gb = [float(r[5]) for r in gr]
    S["granite"] = dict(M=MG, V=MG / RHO_G, GR_main=float(-G * MG * np.mean(1 / e0) / C2), GR_ref=float(-G * MG * np.mean(1 / e1) / C2),
        ext_g=min(gg), Lg=float(Lg[int(np.argmin(gg))]), ext_b=min(gb), Lb=float(Lg[int(np.argmin(gb))]), dmin=float(e0.min()), dmax=float(e0.max()),
        gap_coil_body=float(Rax - D_body / 2 - RG_OUT), gap_coil_winding=float(Rax - a - RG_OUT))
    # ---------------- kapszula (mindkét OCXO-nál azonos → differenciális tag 0)
    CAP = dict(r_in=0.05, h=0.24, t_layers=0.03, r_out=0.08, gap_granite=RG_IN - 0.08)
    t_al = 1.2e-3; V_al = math.pi * ((CAP["r_in"] + t_al) ** 2 * (CAP["h"] + 2 * t_al) - CAP["r_in"] ** 2 * CAP["h"])
    CAP.update(V_al_housing=V_al, m_al_housing=2700.0 * V_al, m_cu_wire20m=DENS * math.pi * (0.75e-3) ** 2 * 20.0,
               r_wire_range=[CAP["r_in"] + t_al, CAP["r_out"]], turns_wire_range=[20.0 / (2 * math.pi * CAP["r_out"]), 20.0 / (2 * math.pi * (CAP["r_in"] + t_al))])
    S["capsule"] = CAP
    # ---------------- órák és áramtábla
    OFFSET = 1 / 3900 / 10e6                 # [SZERZŐ] átlagos beat-offset (1 h 5 min / 360°); nem zajküszöb
    SR = dict(s1=4.8e-17, floor=8.1e-19); Z_SR = -23.0e6 / 429.228004229873e12
    def tau(Sg, s1, k=2): return (2 * k * math.sqrt(2) * s1 / abs(Sg)) ** 2
    A_cmil = (dw / 0.0254e-3) ** 2; LOGT = math.log10(1 + (T_MELT - 20) / 254)
    def t_melt(I): return (A_cmil / I) ** 2 * LOGT / 33
    I_preece = 10244 * (dw / 0.0254) ** 1.5; I_cont = math.sqrt(800.0 * math.pi * Dm * ell / R)
    rows = []
    for lab, I in (("[SZERZŐ] mért 14 mA (RÖGZÍTETT ESET)", I0), ("[SZERZŐ] lista 1.5 A", 1.5), ("[SZERZŐ] lista 2 A", 2.0), ("[SZERZŐ] lista 3 A", 3.0)):
        V = I * R; P_ = I * I * R; Bc = Bpa * I; Br = Bref_pa * I; Bi = Bin_pa * I; Sg = kg * I * I; Sb = kb * I * I
        rows.append(dict(eset=lab, I=I, U=V, P1=P_, P8=8 * P_, Itot=8 * I,
            bat=("igen" if I < 0.02 else f"NEM ({V:.0f} V kellene)"), psu=("igen" if I < 0.02 else f"NEM ({V:.0f} V > 30 V; {8*I:.0f} A > 3 A)"),
            tmelt=t_melt(I), cont=I / I_cont, Bin=Bi, Bc=Bc, Br=Br, p=Bi ** 2 / (2 * MU0), pY=Bi ** 2 / (2 * MU0) / YIELD,
            Sg=Sg, dtg=Sg * DAY, dtg12=Sg * 12 * 3600, Sb=Sb, dtb=Sb * DAY, r_off=abs(Sg) / OFFSET, r_off_b=abs(Sb) / OFFSET,
            r_sr=abs(Sg) / SR["floor"], tau_sr=tau(Sg, SR["s1"]) / YEAR, tau_sr_b=tau(Sb, SR["s1"]) / YEAR,
            z_sr=Z_SR * Bc ** 2, k1=OFFSET / Bc, k2=OFFSET / Bc ** 2))
    H = list(rows[0].keys())
    w("fixed_case_and_feasibility.csv", [[r[k] if isinstance(r[k], str) else g4(r[k]) for k in H] for r in rows],
      ["eset", "I/tekercs (A)", "U = I·R (V)", "P/tekercs (W)", "P 8 tekercs (W)", "forrásáram (A)", "24 V BVS akku", "30 V/3 A táp",
       "Onderdonk olvadási idő (s)", "I / folyamatos határ", "B a tekercsben (T)", "B a fő órán (T)", "B a referencián (T)", "B²/2µ0 (Pa)", "nyomás / Cu folyáshatár",
       "df_f_diff_signed γ_L (max abs, L szerint)", "Δt_diff/nap γ_L (s)", "Δt_diff/12 h γ_L (s)", "df_f_diff_signed γ=1e-22", "Δt_diff/nap γ=1e-22 (s)",
       "abs(df_f_diff) / átlagos beat-offset 2.5641e-11 (γ_L)", "u.a. (γ=1e-22)", "abs(df_f_diff) / hivatkozott Sr rendszeres bizonytalanság 8.1e-19", "τ_AB Sr, fehér-FM képlet algebrai extrapolációja (év, γ_L)", "u.a. (év, γ=1e-22)",
       "Sr kvadratikus Zeeman Δf/f a fő órán", "főóra-egyenérték csatolási együttható k1 (ideális hatás = átlagos beat-offset) (1/T)", "főóra-egyenérték k2 (1/T²)"])
    S["rows"] = rows; S.update(I_preece=I_preece, I_cont=I_cont, OFFSET=OFFSET, SR=SR, Z_SR=Z_SR,
        I_needed_off_g=math.sqrt(OFFSET / abs(kg)), I_needed_off_b=math.sqrt(OFFSET / abs(kb)), Bz_sr_floor=math.sqrt(SR["floor"] / abs(Z_SR)),
        r_off_datasheet=abs(rows[0]["Sg"]) / 1e-11)
    json.dump(S, open(os.path.join(CSV, "fixed_case_numbers.json"), "w"), ensure_ascii=False, indent=1, default=float)
    # ---------------- ábra: differenciális görbék
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    pos = df_g < 0
    ax.loglog(Lg[pos], -df_g[pos], "-", color="C0", label="−df_f_diff (γ_L), fő − referencia")
    ax.loglog(Lg[pos], -df_b[pos], "-", color="C1", label="−df_f_diff (γ=1e-22, kereszt-ellenőrzés)")
    if (~pos).any(): ax.loglog(Lg[~pos], np.abs(df_g[~pos]), "x", color="C0", ms=4, label="abs, ahol ΔW < 0 (numerikus/geometriai)")
    for y, lab in ((OFFSET, "átlagos beat-offset 2.5641e-11 (nem zajküszöb)"), (8.1e-19, "hivatkozott Sr rendszeres bizonytalanság 8.1e-19")): ax.axhline(y, color="k", ls="--", lw=0.7); ax.text(1.2e-2, y * 2, lab, fontsize=7)
    ax.set_xlabel("L (m)"); ax.set_ylabel("abs(df_f_diff), 14 mA/tekercs"); ax.set_title("Differenciális TAL-jel (fő − referencia), 14 mA, menetátmérő 291.7 mm (mért)"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(PL, "fixed_signal_vs_L.png"), dpi=140)
    print(json.dumps(dict(TAL_diff=S["TAL_diff"], Bpa=Bpa, Bref=Bref_pa, mc=mc, W=S["W_total"])))
