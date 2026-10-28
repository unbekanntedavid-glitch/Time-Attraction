#!/usr/bin/env python3
r"""
True top-hat average density rho_bar_L around SPARC rotation-curve points,
and the TAL factor f_L = exp(-gamma_L rho_bar_L), gamma_L = G L^2 / c^2.

Purpose: Academia review request on the seventh-model note
(tal_factor_7th_model.md, V19, DOI 10.5281/zenodo.23497128).  That note
uses the midplane local density proxy at every L (incl. the fitted
L = 2.3714 Mpc).  Here rho_bar_L = 3/(4 pi L^3) * \int_{|x'-x|<=L} rho dV
is computed for a 3D baryonic model of every real SPARC galaxy.

Data (real only; no synthetic fallback - the script stops if data are missing):
  Rotmod_LTG.zip      SPARC (Lelli, McGaugh & Schombert 2016), SHA256 0a80cc90...f588
  SPARC_Lelli2016c.mrt SPARC Table 1 (Rdisk, MHI, RHI), SHA256 5aa0501f...8243
  blind_123.txt       hash-locked v17 sample list (only used as a flag)
  points_tal_fL_FITTED_forward.csv (repo, hetedik_modell_kimenet/) - used only
                      to verify that rho_midplane reproduces the note exactly.

Usage: python3 sparc_tophat.py --data DIR --repo REPO --out OUTDIR [--procs 8]
All densities in the output are SI (kg m^-3); radii and L in kpc.
"""
import argparse, csv, hashlib, io, json, math, sys, time, zipfile
from pathlib import Path
from multiprocessing import Pool
import numpy as np

# ---------------------------------------------------------------- constants
G = 6.67430e-11; C = 299792458.0; C2 = C * C
MSUN = 1.98847e30; PC = 3.085677581e16; KPC = 1e3 * PC; MPC = 1e6 * PC
RHO_UNIT = MSUN / KPC**3            # 1 Msun/kpc^3 in kg/m^3
G_KPC = G * MSUN / KPC / 1e6        # kpc (km/s)^2 / Msun
# seventh-model note constants (tal_factor_7th_model.py)
HZ_NOTE = 0.2                        # kpc
YDISK = 0.5; YBUL = 0.7
GAS_FACTOR = 1.33                    # He correction (SPARC convention)
# Planck 2018 (TT,TE,EE+lowE+lensing): H0=67.4, Omega_m=0.315, Omega_b h^2=0.0224
H0 = 67.4e3 / MPC
RHO_CRIT = 3 * H0**2 / (8 * math.pi * G)
OMEGA_M = 0.315; OMEGA_B = 0.0224 / 0.674**2
RHO_BG_M = OMEGA_M * RHO_CRIT; RHO_BG_B = OMEGA_B * RHO_CRIT

EXPECTED = {"Rotmod_LTG.zip": "0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588",
            "SPARC_Lelli2016c.mrt": "5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243",
            "blind_123.txt": "3637d2f5bda6b2b97e645c5fba5cf39d310a1826bee388286c2c74b9171aa283"}

L_FIT = 2.371373705661655 * 1e3      # kpc, run_meta.json forward fit
SPECIAL_L = {"fitted_2.3714Mpc": L_FIT, "v17_L0_2.36kpc": 2.36,
             "Lunity_27.64Mpc": 27.641074834585233e3, "L_64.2Mpc": 64.2e3, "L_85.3Mpc": 85.3e3}

def gamma_of_L(Lkpc): return G * (Lkpc * KPC)**2 / C2

def sha256(path):
    h = hashlib.sha256(); h.update(Path(path).read_bytes()); return h.hexdigest()

# ---------------------------------------------------------------- data input
def read_mrt(path):
    """Whitespace-split parse of SPARC Table 1 data rows (the byte offsets in the
    header do not match the file's actual spacing). Tokens: Galaxy T D e_D f_D Inc
    e_Inc L36 e_L36 Reff SBeff Rdisk SBdisk MHI RHI Vflat e_Vflat Q Ref."""
    gal = {}
    for line in Path(path).read_text().splitlines():
        p = line.split()
        if len(p) < 18: continue
        try:
            int(p[1]); vals = [float(x) for x in p[2:17]]; int(p[17])
        except ValueError:
            continue
        gal[p[0]] = dict(D=vals[0], L36=vals[5], Rdisk=vals[9], MHI=vals[11], RHI=vals[12])
    return gal

def parse_rotmod_text(txt):
    """Same rows as runner.parse_rotmod (>=6 columns, numeric) incl. R<0.01 -> 0.16 fix."""
    rows = []
    for line in txt.splitlines():
        p = line.split()
        if len(p) < 6: continue
        try:
            r = [float(p[0]), float(p[1]), float(p[3]), float(p[4]), float(p[5]),
                 float(p[6]) if len(p) > 6 else 0.0, float(p[7]) if len(p) > 7 else 0.0]
        except ValueError:
            continue
        rows.append(r)
    if rows and rows[0][0] < 0.01:
        rows[0][0] = 0.16
    return rows

# ---------------------------------------------------------------- quadrature
def gl(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1), 0.5 * w           # on [0,1]

# ---------------------------------------------------------------- galaxy model
class Galaxy:
    def __init__(self, name, rows, meta, hz_star=HZ_NOTE, hz_gas=HZ_NOTE, gas_scale=1.0):
        self.name = name
        a = np.array(rows, float)
        self.R = a[:, 0]; self.Vobs = a[:, 1]; self.Vgas = a[:, 2]; self.Vdisk = a[:, 3]
        self.Vbul = a[:, 4]; self.SBd = a[:, 5]; self.SBb = a[:, 6]
        self.meta = meta
        self.hs = hz_star; self.hg = hz_gas
        self.Rd = meta["Rdisk"] if meta["Rdisk"] > 0 else max(self.R[-1] / 4, 0.1)
        # stellar disk surface density, Msun/kpc^2
        self.Sig_d = YDISK * self.SBd * 1e6
        self.Rlast = self.R[-1]; self.Slast = self.Sig_d[-1]
        self.Rcut_d = self.Rlast + 15 * self.Rd if self.Slast > 0 else self.Rlast
        # gas disk: exponential, M = 1.33 MHI, Sigma_HI(RHI) = 1 Msun/pc^2
        MHI = meta["MHI"] * 1e9; RHI = meta["RHI"]
        self.Mg = GAS_FACTOR * MHI
        self.gas_rule = "none"
        self.Rg = 0.0
        if MHI > 0:
            if RHI > 0:
                # f(Rg) = MHI/(2 pi Rg^2) exp(-RHI/Rg) [Msun/kpc^2] = 1e6
                fmax = MHI / (2 * math.pi * (RHI / 2)**2) * math.exp(-2)
                if fmax >= 1e6:
                    lo, hi = RHI / 2, 1e4 * RHI
                    for _ in range(200):
                        mid = math.sqrt(lo * hi)
                        f = MHI / (2 * math.pi * mid**2) * math.exp(-RHI / mid)
                        lo, hi = (mid, hi) if f > 1e6 else (lo, mid)
                    self.Rg = math.sqrt(lo * hi); self.gas_rule = "Sigma_HI(RHI)=1 (outer root)"
                else:
                    self.Rg = RHI / 2; self.gas_rule = "no root; Rg=RHI/2"
            else:
                self.Rg = 2 * self.Rd; self.gas_rule = "RHI=0; Rg=2 Rdisk"
        self.Rg *= gas_scale
        self.Sig_g_RHI = (self.Mg / GAS_FACTOR / (2 * math.pi * self.Rg**2) * math.exp(-RHI / self.Rg) / 1e6) if (self.Mg > 0 and RHI > 0) else float('nan')
        self.Rcut_g = 25 * self.Rg if self.Mg > 0 else 0.0
        # bulge: spherical Hernquist fitted to the enclosed-mass curve implied by the
        # SPARC bulge velocities (SPARC computes V_bul for a spherical deprojection of SBbul):
        # M_enc(R_i) = Upsilon_bul * V_bul(R_i)^2 * R_i / G
        self.Mb = 0.0; self.ab = 0.0; self.bulge_fit = "none"; self.Mb_proj = 0.0
        m = (self.Vbul > 0) & (self.SBb > 0)
        if m.sum() >= 2:
            Rb = self.R[m]; Me = YBUL * self.Vbul[m]**2 * Rb / G_KPC
            best = None
            for ab in np.geomspace(1e-3, 50, 4000):
                prof = Rb**2 / (Rb + ab)**2
                off = np.mean(np.log(Me) - np.log(prof))
                res = np.sum((np.log(Me) - np.log(prof) - off)**2)
                if best is None or res < best[0]: best = (res, ab, math.exp(off))
            self.ab = best[1]; self.Mb = best[2]
            # projected-light cross-check: Upsilon_bul * int 2 pi R SBbul dR over tabulated radii
            Rp = np.concatenate([[0.0], self.R]); Sp = np.concatenate([[YBUL * self.SBb[0] * 1e6], YBUL * self.SBb * 1e6])
            self.Mb_proj = float(np.trapezoid(2 * math.pi * Rp * Sp, Rp))
            self.bulge_fit = (f"Hernquist fit to Ybul*Vbul^2*R/G: a={self.ab:.4g} kpc, M={self.Mb:.4g} Msun, "
                              f"rms_dex={math.sqrt(best[0]/m.sum())/math.log(10):.3f}, n={int(m.sum())}")
            if best[1] <= 1.0001e-3 or best[1] >= 49.99: self.bulge_fit += " (a at grid edge!)"
        self.Md = self.disk_total_mass()
        self.Mtot = self.Md + self.Mg + self.Mb
        # radius beyond which the ball centred at R0 surely contains everything
        self.Rmax_all = max(self.Rcut_d, self.Rcut_g, 0.0)

    @staticmethod
    def hern_proj(R, M, a):
        s = np.asarray(R, float) / a
        X = np.empty_like(s)
        lo = s < 1 - 1e-6; hi = s > 1 + 1e-6; mid = ~(lo | hi)
        X[lo] = np.arccosh(1 / s[lo]) / np.sqrt(1 - s[lo]**2)
        X[hi] = np.arccos(1 / s[hi]) / np.sqrt(s[hi]**2 - 1)
        X[mid] = 1.0
        out = np.empty_like(s)
        nm = ~mid
        out[nm] = M / (2 * math.pi * a**2 * (1 - s[nm]**2)**2) * ((2 + s[nm]**2) * X[nm] - 3)
        out[mid] = 2 * M / (15 * math.pi * a**2)
        return out

    def sigma_star(self, R):
        R = np.asarray(R, float)
        out = np.interp(R, self.R, self.Sig_d)          # R<R1: Sigma(R1) (flat core)
        beyond = R > self.Rlast
        out = np.where(beyond, self.Slast * np.exp(-(R - self.Rlast) / self.Rd), out)
        return np.where(R > self.Rcut_d, 0.0, out)

    def sigma_gas(self, R):
        R = np.asarray(R, float)
        if self.Mg <= 0: return np.zeros_like(R)
        out = self.Mg / (2 * math.pi * self.Rg**2) * np.exp(-R / self.Rg)
        return np.where(R > self.Rcut_g, 0.0, out)

    def rho_bulge(self, r):
        if self.Mb <= 0: return np.zeros_like(np.asarray(r, float))
        r = np.asarray(r, float)
        return self.Mb * self.ab / (2 * math.pi * r * (r + self.ab)**3)

    def disk_total_mass(self):
        brk = np.unique(np.concatenate([[0.0], self.R, [self.Rlast + k * self.Rd for k in range(1, 16)]]))
        brk = brk[brk <= self.Rcut_d + 1e-12]
        x, w = gl(16); M = 0.0
        for a, b in zip(brk[:-1], brk[1:]):
            R = a + (b - a) * x
            M += np.sum((b - a) * w * 2 * math.pi * R * self.sigma_star(R))
        return M

    def rho_local(self, R0):
        """Model density at (R0, z=0), Msun/kpc^3."""
        return (self.sigma_star(R0) / (2 * self.hs) + self.sigma_gas(R0) / (2 * self.hg)
                + (self.rho_bulge(R0) if self.Mb > 0 else 0.0))

    # ---- ball mass of an axisymmetric disk with exponential vertical profile
    def disk_ball_mass(self, sigma_fn, h, Rcut, R0, L, extra_brk, scale, nR=16, nT=32):
        Rmin = max(0.0, R0 - L); Rmax = min(R0 + L, Rcut)
        if Rmax <= Rmin: return 0.0
        brk = [Rmin, Rmax, abs(L - R0), R0] + list(extra_brk)
        brk = np.unique(np.clip(np.array(brk, float), Rmin, Rmax))
        dmax = max(min(L / 4, scale / 2), 1e-4)
        edges = []
        for a, b in zip(brk[:-1], brk[1:]):
            n = max(1, int(math.ceil((b - a) / dmax)))
            edges.append(np.linspace(a, b, n + 1)[:-1])
        edges = np.concatenate(edges + [[brk[-1]]])
        xR, wR = gl(nR)
        a = edges[:-1, None]; b = edges[1:, None]
        R = (a + (b - a) * xR).ravel(); wr = ((b - a) * wR).ravel()
        # vertical: u = 1 - exp(-z/h) on z>=0, factor 2 for both sides, rho dz = Sigma/2 du
        zmax2 = L**2 - (R - R0)**2
        zmax = np.sqrt(np.clip(zmax2, 0, None))
        umax = -np.expm1(-zmax / h)
        zf2 = L**2 - (R + R0)**2                      # full-ring threshold
        zfull = np.sqrt(np.clip(zf2, 0, None))
        ufull = np.where(zf2 > 0, -np.expm1(-zfull / h), 0.0)
        # integral over u in [ufull, umax] of dphi(u): substitution u = umax - (umax-ufull) t^2
        xT, wT = gl(nT)
        span = (umax - ufull)[:, None]
        u = umax[:, None] - span * xT[None, :]**2
        with np.errstate(divide="ignore"):
            z = -h * np.log1p(-np.clip(u, 0, 1.0))   # u==1 -> z=inf: outside ball, dphi=0 (density weight already in du)
        with np.errstate(divide="ignore", invalid="ignore"):
            cphi = (R[:, None]**2 + R0**2 + z**2 - L**2) / (2 * R[:, None] * R0)
        dphi = 2 * np.arccos(np.clip(cphi, -1, 1))
        Iu = 2 * math.pi * ufull + np.sum(dphi * 2 * span * xT[None, :] * wT[None, :], axis=1)
        # mass = int R dR Sigma(R) * (1/2 * 2) * Iu   [Sigma/2 per side, two sides]
        return float(np.sum(wr * R * sigma_fn(R) * Iu))

    def bulge_ball_mass(self, d, L, n=48):
        if self.Mb <= 0: return 0.0
        Menc = lambda r: self.Mb * r**2 / (r + self.ab)**2
        M = Menc(L - d) if L > d else 0.0
        lo, hi = abs(L - d), L + d
        # geometric sub-intervals (rho has a 1/r cusp)
        pts = np.unique(np.concatenate([[lo, hi], np.geomspace(max(lo, 1e-6 * hi), hi, 9)]))
        x, w = gl(n)
        for a, b in zip(pts[:-1], pts[1:]):
            r = a + (b - a) * x
            frac = (L**2 - (r - d)**2) / (4 * r * d)
            M += np.sum((b - a) * w * 4 * math.pi * r**2 * self.rho_bulge(r) * np.clip(frac, 0, 1))
        return float(M)

    def ball_mass(self, R0, L, nR=16, nT=32, nB=64, force_numeric=False):
        """Returns (M_star, M_gas, M_bulge) in Msun inside the ball."""
        if (not force_numeric) and L >= R0 + self.Rmax_all + 60 * max(self.hs, self.hg):
            # ball contains the whole (truncated) disk; bulge part analytic/numeric
            return self.Md, self.Mg if self.Mg > 0 else 0.0, self.bulge_ball_mass(R0, L, nB)
        ext = [r for r in self.R] + [self.Rlast + k * self.Rd for k in range(1, 16)]
        Ms = self.disk_ball_mass(self.sigma_star, self.hs, self.Rcut_d, R0, L, ext, self.Rd, nR, nT)
        Mg = 0.0
        if self.Mg > 0:
            extg = [k * self.Rg for k in range(1, 26)]
            Mg = self.disk_ball_mass(self.sigma_gas, self.hg, self.Rcut_g, R0, L, extg, self.Rg, nR, nT)
        return Ms, Mg, self.bulge_ball_mass(R0, L, nB)

def rho_mid_note(R_kpc, Vb, sbd):
    """Exact copy of rho_diskbulge() in tal_factor_7th_model.py (SI)."""
    rho_d = YDISK * sbd * MSUN / PC**2 / (2 * HZ_NOTE * KPC)
    r = R_kpc * KPC
    rho_b = 3 * YBUL * (Vb * 1e3)**2 / (4 * math.pi * G * r**2)
    return rho_d + rho_b

# ---------------------------------------------------------------- worker
GALS = {}
def init_worker(gal_inputs, hs, hg, gs=1.0):
    global GALS
    GALS = {n: Galaxy(n, rows, meta, hs, hg, gs) for n, (rows, meta) in gal_inputs.items()}

def L_grid():
    core = list(np.geomspace(0.1, 1e4, 64))            # 0.1 kpc .. 10 Mpc, 64 points
    rows = [(float(L), "grid") for L in core]
    for k, v in SPECIAL_L.items(): rows.append((v, k))
    rows.sort()
    return rows

def work(gname):
    g = GALS[gname]; out = []
    Ls = L_grid()
    for i in range(len(g.R)):
        if g.Vobs[i] <= 0: continue
        R0 = g.R[i]
        rmid = rho_mid_note(R0, g.Vbul[i], g.SBd[i])
        rloc = float(g.rho_local(R0)) * RHO_UNIT
        for L, flag in Ls:
            Ms, Mg, Mb = g.ball_mass(R0, L)
            rhoL = (Ms + Mg + Mb) / (4 / 3 * math.pi * L**3) * RHO_UNIT
            out.append((gname, i, R0, L, flag, rhoL, Ms, Mg, Mb, rmid, rloc))
    return out

# ---------------------------------------------------------------- checks
def numerical_checks(gals, rng):
    chk = {}
    # 1) analytic: uniform sheet Sigma0 with exponential vertical profile, ball fully inside in R
    class Sheet(Galaxy):
        def __init__(s): pass
    s = Sheet(); s.hs = s.hg = 0.2
    S0 = 1e8
    def M_exact(L, h):
        # int_{-L}^{L} pi (L^2 - z^2) S0/(2h) e^{-|z|/h} dz
        e = math.exp(-L / h)   # closed form of the integral above
        return math.pi * S0 * (L**2 * (1 - e) - (2 * h**2 - e * (L**2 + 2 * L * h + 2 * h**2)))
    rel = []
    for L in (0.01, 0.1, 0.2, 1.0, 5.0):
        R0 = 1000.0
        Mn = Galaxy.disk_ball_mass(s, lambda R: np.full_like(R, S0), 0.2, 1e9, R0, L, [], 1e3)
        rel.append(abs(Mn / M_exact(L, 0.2) - 1))
    chk["uniform_sheet_vs_analytic_max_rel_err"] = max(rel)
    # 2) node doubling on random (galaxy, R0, L) samples
    names = sorted(gals); diffs = []; worst = None
    Ls = [L for L, _ in L_grid()]
    for _ in range(400):
        g = gals[names[rng.integers(len(names))]]
        i = rng.integers(len(g.R)); L = Ls[rng.integers(len(Ls))]
        a = sum(g.ball_mass(g.R[i], L, force_numeric=True))
        b = sum(g.ball_mass(g.R[i], L, nR=32, nT=64, nB=128, force_numeric=True))
        d = abs(a / b - 1) if b > 0 else 0.0
        diffs.append(d)
        if worst is None or d > worst[0]: worst = (d, g.name, float(g.R[i]), L)
    chk["node_doubling_max_rel_diff"] = max(diffs)
    chk["node_doubling_median_rel_diff"] = float(np.median(diffs))
    chk["node_doubling_worst_case"] = worst
    # 3) shortcut (whole galaxy inside) vs forced numeric
    dd = []
    for n in names[:40]:
        g = gals[n]; R0 = g.R[len(g.R) // 2]; L = R0 + g.Rmax_all + 60 * g.hs + 1
        a = sum(g.ball_mass(R0, L)); b = sum(g.ball_mass(R0, L, nR=20, nT=48, force_numeric=True))
        dd.append(abs(a / b - 1))
    chk["full_containment_shortcut_vs_numeric_max_rel_diff"] = max(dd)
    # 4) L -> 0 limit: rho_bar_L -> local model density
    ll = []
    for n in names:
        g = gals[n]
        for i in (0, len(g.R) // 2):
            R0 = g.R[i]; L = 1e-4
            rhoL = sum(g.ball_mass(R0, L)) / (4 / 3 * math.pi * L**3)
            loc = float(g.rho_local(R0))
            if loc > 0: ll.append(abs(rhoL / loc - 1))
    chk["small_L_1e-4kpc_vs_local_density_max_rel_diff"] = max(ll)
    # 5) independent Monte Carlo (uniform points in the ball, 4e6 samples) for a few cases
    mc = []
    for n, frac, L in (("NGC3198", 0.5, 5.0), ("NGC7814", 0.2, 3.0), ("DDO154", 0.5, 2.0),
                       ("UGC02953", 0.3, 20.0), ("NGC2841", 0.1, 1.0)):
        if n not in gals: continue
        g = gals[n]; R0 = float(g.R[int(frac * (len(g.R) - 1))])
        N = 4_000_000
        v = rng.normal(size=(N, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
        p = v * (L * rng.random(N) ** (1 / 3))[:, None]
        x = p[:, 0] + R0; y = p[:, 1]; z = p[:, 2]
        Rc = np.hypot(x, y); r = np.sqrt(Rc**2 + z**2)
        rho = (g.sigma_star(Rc) / (2 * g.hs) * np.exp(-np.abs(z) / g.hs)
               + g.sigma_gas(Rc) / (2 * g.hg) * np.exp(-np.abs(z) / g.hg)
               + (g.rho_bulge(r) if g.Mb > 0 else 0.0))
        mean = rho.mean(); se = rho.std() / math.sqrt(N)
        quad = sum(g.ball_mass(R0, L)) / (4 / 3 * math.pi * L**3)
        mc.append(dict(galaxy=n, R0=R0, L=L, quad=quad, mc=float(mean), mc_se=float(se),
                       pull=float((quad - mean) / se)))
    chk["monte_carlo"] = mc
    return chk

# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--repo", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--hz_star", type=float, default=HZ_NOTE); ap.add_argument("--hz_gas", type=float, default=HZ_NOTE)
    ap.add_argument("--gas_scale", type=float, default=1.0); ap.add_argument("--tag", default=""); ap.add_argument("--galaxies", default="")
    ap.add_argument("--skip_checks", action="store_true")
    A = ap.parse_args()
    data = Path(A.data); repo = Path(A.repo); out = Path(A.out); out.mkdir(parents=True, exist_ok=True)
    files = {"Rotmod_LTG.zip": data / "Rotmod_LTG.zip", "SPARC_Lelli2016c.mrt": data / "SPARC_Lelli2016c.mrt",
             "blind_123.txt": repo / "v17_package/github/data/blind_123.txt"}
    for k, p in files.items():
        if not p.is_file(): sys.exit(f"ERROR: missing real input {p} (no synthetic fallback)")
        h = sha256(p)
        if h != EXPECTED[k]: sys.exit(f"ERROR: hash mismatch {k}: {h}")
        print(f"OK {k} {h}")
    mrt = read_mrt(files["SPARC_Lelli2016c.mrt"])
    blind = set(files["blind_123.txt"].read_text().split())
    gal_inputs = {}
    with zipfile.ZipFile(files["Rotmod_LTG.zip"]) as z:
        for nm in sorted(z.namelist()):
            if not nm.endswith("_rotmod.dat"): continue
            g = Path(nm).name[:-len("_rotmod.dat")]
            rows = parse_rotmod_text(z.read(nm).decode())
            if g not in mrt: sys.exit(f"ERROR: {g} not in SPARC table")
            gal_inputs[g] = (rows, mrt[g])
    if A.galaxies: gal_inputs = {k: v for k, v in gal_inputs.items() if k in A.galaxies.split(",")}
    print(f"galaxies: {len(gal_inputs)}; blind-123 present: {len(blind & set(gal_inputs))}")
    init_worker(gal_inputs, A.hz_star, A.hz_gas, A.gas_scale)
    tag = A.tag
    # ---- model summary per galaxy
    gm = []
    for n, g in GALS.items():
        gm.append(dict(galaxy=n, in_blind123=int(n in blind), Rdisk_kpc=g.Rd, R_last_kpc=g.Rlast,
                       M_star_disk_Msun=g.Md, M_gas_Msun=g.Mg, Rgas_kpc=g.Rg, gas_rule=g.gas_rule,
                       M_bulge_Msun=g.Mb, a_bulge_kpc=g.ab, bulge_fit=g.bulge_fit, M_bulge_projlight_Msun=g.Mb_proj, SigmaHI_at_RHI=g.Sig_g_RHI,
                       M_star_from_L36_Msun=(YDISK * g.meta["L36"] * 1e9)))
    t0 = time.time()
    with Pool(A.procs, initializer=init_worker, initargs=(gal_inputs, A.hz_star, A.hz_gas, A.gas_scale)) as pool:
        res = pool.map(work, sorted(GALS), chunksize=1)
    print(f"integration done in {time.time()-t0:.1f} s")
    rows = [r for rr in res for r in rr]
    # ---- verify rho_midplane against the note's own output (blind-123, 2596 points)
    ref = {}
    fp = repo / "hetedik_modell_kimenet/points_tal_fL_FITTED_forward.csv"
    for r in csv.DictReader(open(fp)):
        ref.setdefault(r["galaxy"], []).append((float(r["R_kpc"]), float(r["rho_bar_kg_m3"])))
    mine = {}
    seen = set()
    for r in rows:
        key = (r[0], r[1])
        if key in seen: continue
        seen.add(key); mine.setdefault(r[0], []).append((r[2], r[9]))
    nchk = 0; maxrel = 0.0
    for gname, lst in ref.items():
        if gname not in mine: continue
        assert len(lst) == len(mine[gname]), gname
        for (Ra, ra), (Rb, rb) in zip(lst, mine[gname]):
            assert abs(Ra - Rb) < 1e-12
            maxrel = max(maxrel, abs(ra - rb) / ra if ra > 0 else abs(rb)); nchk += 1
    print(f"rho_midplane vs note points_tal_fL_FITTED_forward.csv: {nchk} points, max rel diff {maxrel:.3e}")
    # ---- write full CSV
    hdr = ["galaxy", "in_blind123", "R_eval_kpc", "L_kpc", "L_flag", "rho_L", "rho_L_bg_baryon", "rho_L_bg_matter",
           "rho_midplane", "rho_local_model", "ratio_rhoL_over_rhomid", "one_minus_ratio_abs",
           "ratio_rhoLbgm_over_rhomid", "one_minus_ratio_bgm_abs", "gammaL_m3_per_kg",
           "gammaL_rhoL", "fL_TAL", "one_minus_fL_TAL_abs",
           "gammaL_rhoL_bg_matter", "one_minus_fL_TAL_bg_matter_abs",
           "gammaL_rhomid_note", "one_minus_fL_note_midplane_abs", "M_star_in_ball_Msun", "M_gas_in_ball_Msun", "M_bulge_in_ball_Msun"]
    full = out / f"sparc_tophat_full{tag}.csv"
    stats = []
    with open(full, "w", newline="") as f:
        f.write("# units: R_eval_kpc, L_kpc in kpc; all rho_* in kg m^-3 (SI); gammaL in m^3 kg^-1; gamma*rho dimensionless; masses in Msun\n")
        f.write(f"# rho_L = true top-hat mean of the galaxy model (no background); rho_L_bg_* add Planck-2018 Omega_b rho_crit={RHO_BG_B:.4e} or Omega_m rho_crit={RHO_BG_M:.4e}\n")
        f.write("# rho_midplane = seventh-model-note proxy (disk midplane + bulge enclosed-density, no gas), as in tal_factor_7th_model.py\n")
        f.write("# ratio_* = user-requested ratio rho_L/rho_midplane (NOT the TAL f_L); fL_TAL = exp(-gamma_L rho_L), gamma_L = G L^2/c^2; one_minus = -expm1(-x)\n")
        w = csv.writer(f); w.writerow(hdr); FL_COL = hdr.index("fL_TAL")
        for (gname, i, R0, L, flag, rhoL, Ms, Mg, Mb, rmid, rloc) in rows:
            gam = gamma_of_L(L)
            rb = rhoL + RHO_BG_B; rm = rhoL + RHO_BG_M
            ratio = rhoL / rmid if rmid > 0 else float("nan")
            ratm = rm / rmid if rmid > 0 else float("nan")
            x = gam * rhoL; xm = gam * rm; xn = gam * rmid
            rec = [gname, int(gname in blind), R0, L, flag, rhoL, rb, rm, rmid, rloc, ratio, abs(1 - ratio),
                   ratm, abs(1 - ratm), gam, x, math.exp(-x), -math.expm1(-x), xm, -math.expm1(-xm),
                   xn, -math.expm1(-xn), Ms, Mg, Mb]
            w.writerow([(f"{v:.17g}" if j == FL_COL else f"{v:.10g}") if isinstance(v, float) else v for j, v in enumerate(rec)])
            stats.append((gname, int(gname in blind), R0, L, flag, rhoL, rmid, ratio, abs(1 - ratio), x, -math.expm1(-x),
                          xm, -math.expm1(-xm), xn, -math.expm1(-xn), abs(1 - ratm)))
    print("wrote", full, len(rows), "rows")
    # ---- maxima
    def argmax(sel, col):
        best = None
        for s in stats:
            if not sel(s): continue
            v = s[col]
            if v != v: continue
            if best is None or v > best[col]: best = s
        return best
    def desc(s):
        if s is None: return None
        return dict(galaxy=s[0], in_blind123=s[1], R_kpc=s[2], L_kpc=s[3], L_flag=s[4], rho_L=s[5], rho_midplane=s[6],
                    ratio=s[7], one_minus_ratio_abs=s[8], gammaL_rhoL=s[9], one_minus_fL_TAL_abs=s[10],
                    gammaL_rhoL_bgm=s[11], one_minus_fL_TAL_bgm_abs=s[12], gammaL_rhomid_note=s[13],
                    one_minus_fL_note_abs=s[14])
    core = lambda s: s[3] <= 1e4 * (1 + 1e-12)
    summ = {"n_galaxies": len(GALS), "n_points": len(seen), "n_rows": len(rows),
            "n_points_rho_mid_zero": sum(1 for k in seen if False),
            "rho_bg_baryon": RHO_BG_B, "rho_bg_matter": RHO_BG_M, "rho_crit": RHO_CRIT,
            "rho_mid_verification": {"n_points": nchk, "max_rel_diff": maxrel}}
    zero_mid = sum(1 for s in stats if s[6] <= 0) // len(L_grid())
    summ["n_points_rho_mid_zero"] = zero_mid
    for lab, sel in (("all175_L_le_10Mpc", core), ("all175_incl_special_L_to_85Mpc", lambda s: True),
                     ("blind123_L_le_10Mpc", lambda s: core(s) and s[1] == 1)):
        summ[lab] = {"max_one_minus_ratio": desc(argmax(sel, 8)),
                     "max_one_minus_ratio_suppression_side(ratio<1)": desc(argmax(lambda s: sel(s) and s[7] < 1, 8)),
                     "max_one_minus_fL_TAL_galaxy_only": desc(argmax(sel, 10)),
                     "max_one_minus_fL_TAL_with_matter_bg": desc(argmax(sel, 12)),
                     "max_one_minus_fL_note_midplane_assumption": desc(argmax(sel, 14))}
    fit = [s for s in stats if s[4] == "fitted_2.3714Mpc"]
    summ["at_L_2.3714Mpc"] = {
        "true_galaxy_only_max_gammaL_rhoL": desc(max(fit, key=lambda s: s[9])),
        "true_with_matter_bg_max_gammaL_rhoL": desc(max(fit, key=lambda s: s[11])),
        "note_midplane_max_gammaL_rhomid": desc(max(fit, key=lambda s: s[13])),
        "median_true_gammaL_rhoL": float(np.median([s[9] for s in fit])),
        "median_true_gammaL_rhoL_bgm": float(np.median([s[11] for s in fit])),
        "median_note_gammaL_rhomid": float(np.median([s[13] for s in fit])),
        "gammaL": gamma_of_L(L_FIT), "gammaL_times_rho_bg_matter": gamma_of_L(L_FIT) * RHO_BG_M}
    # ---- refine the galaxy-only TAL maximum in L (fine scan around the top-5 galaxies' peak points)
    ref_rows = []
    pk = {}
    for t in stats:
        if t[10] > pk.get(t[0], (0, None))[0]: pk[t[0]] = (t[10], t[2])
    top = sorted(pk.items(), key=lambda kv: -kv[1][0])[:5]
    for gname, (_, R0) in top:
        g = GALS[gname]
        Lf = np.geomspace(0.1, 300.0, 1200)
        x = [gamma_of_L(L) * sum(g.ball_mass(R0, L)) / (4 / 3 * math.pi * L**3) * RHO_UNIT for L in Lf]
        k = int(np.argmax(x))
        ref_rows.append(dict(galaxy=gname, R_kpc=R0, L_peak_kpc=float(Lf[k]), gammaL_rhoL_peak=float(x[k]),
                             one_minus_fL_peak=float(-math.expm1(-x[k]))))
    summ["fine_L_refinement_top5"] = ref_rows
    print("fine L refinement:", ref_rows)
    # ---- per-galaxy max CSV
    pg = out / f"sparc_tophat_per_galaxy_max{tag}.csv"
    gmd = {d["galaxy"]: d for d in gm}
    with open(pg, "w", newline="") as f:
        f.write("# densities kg m^-3; L, R in kpc; masses Msun; core grid L<=10 Mpc unless column says 2.3714Mpc\n")
        cols = ["galaxy", "in_blind123", "n_points", "max_one_minus_ratio_abs", "R_at_max_ratio_kpc", "L_at_max_ratio_kpc",
                "max_one_minus_fL_TAL_abs", "R_at_max_fL_kpc", "L_at_max_fL_kpc", "rho_L_at_max_fL",
                "max_one_minus_fL_TAL_bg_matter_abs", "L_at_max_fL_bgm_kpc",
                "at2.37Mpc_max_one_minus_fL_TAL_true", "at2.37Mpc_max_one_minus_fL_TAL_true_bgm",
                "at2.37Mpc_max_one_minus_fL_note_midplane", "at2.37Mpc_max_gammaL_rhomid_note",
                "M_star_disk_Msun", "M_gas_Msun", "M_bulge_Msun", "M_bulge_projected_light_Msun", "M_star_from_0.5L36_Msun", "Rgas_kpc", "SigmaHI_model_at_RHI_Msun_pc2", "gas_rule", "bulge_fit"]
        w = csv.writer(f); w.writerow(cols)
        by = {}
        for s in stats: by.setdefault(s[0], []).append(s)
        for n in sorted(by):
            ss = [s for s in by[n] if core(s)]
            a = max((s for s in ss if s[8] == s[8]), key=lambda s: s[8], default=None)
            b = max(ss, key=lambda s: s[10]); c = max(ss, key=lambda s: s[12])
            ff = [s for s in by[n] if s[4] == "fitted_2.3714Mpc"]
            d = gmd[n]
            w.writerow([n, d["in_blind123"], len(ff),
                        f"{a[8]:.6g}" if a else "nan", a[2] if a else "", a[3] if a else "",
                        f"{b[10]:.6g}", b[2], f"{b[3]:.6g}", f"{b[5]:.6g}", f"{c[12]:.6g}", f"{c[3]:.6g}",
                        f"{max(s[10] for s in ff):.6g}", f"{max(s[12] for s in ff):.6g}",
                        f"{max(s[14] for s in ff):.6g}", f"{max(s[13] for s in ff):.6g}",
                        f"{d['M_star_disk_Msun']:.5g}", f"{d['M_gas_Msun']:.5g}", f"{d['M_bulge_Msun']:.5g}", f"{d['M_bulge_projlight_Msun']:.5g}",
                        f"{d['M_star_from_L36_Msun']:.5g}", f"{d['Rgas_kpc']:.4g}", f"{d['SigmaHI_at_RHI']:.4g}", d["gas_rule"], d["bulge_fit"]])
    print("wrote", pg)
    if not A.skip_checks:
        summ["numerical_checks"] = numerical_checks(GALS, np.random.default_rng(20261022))
    json.dump(summ, open(out / f"summary{tag}.json", "w"), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float))

if __name__ == "__main__":
    main()
