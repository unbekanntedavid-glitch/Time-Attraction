# -*- coding: utf-8 -*-
"""Véges szolenoid (áramlap) tere – Derby & Olbert, Am. J. Phys. 78, 229 (2010), doi:10.1119/1.3256157,
Bulirsch-féle általánosított másodfajú elliptikus integrál (cel). Ellenőrzés: köráram-összeg (2000 menet) és
a tengelyen érvényes zárt alak B = µ0 n I (ζ+/√(ζ+²+a²) − ζ−/√(ζ−²+a²))/2."""
import numpy as np
from scipy.special import ellipk, ellipe

MU0 = 4e-7 * np.pi

def cel(kc, p, c, s, iters=12):   # 12 iteráció: 40-hez képest <1e-15 eltérés (ellenőrizve)
    kc = np.abs(np.asarray(kc, float)); p = np.asarray(p, float) * np.ones_like(kc)
    c = np.asarray(c, float) * np.ones_like(kc); s = np.asarray(s, float) * np.ones_like(kc)
    # csak p > 0 ág (itt mindig p>0, kivéve ρ=a, amit a hívó elkerül)
    k = kc.copy(); em = np.ones_like(k)
    pp = np.sqrt(p); ss = s / pp; cc = c.copy()
    f = cc; cc = cc + ss / pp; g = k / pp; ss = 2 * (ss + f * g); pp = g + pp; g = em.copy(); em = k + em; kk = k.copy()
    for _ in range(iters):
        k = 2 * np.sqrt(kk); kk = k * em; f = cc; cc = cc + ss / pp; g = kk / pp
        ss = 2 * (ss + f * g); pp = g + pp; g = em.copy(); em = k + em
    return (np.pi / 2) * (ss + cc * em) / (em * (em + pp))

def solenoid_B(rho, z, a, ell, NI):
    """Saját hengerkoordinátás tér (Bρ, Bz) egy a sugarú, ell hosszú, origó-középpontú áramlapra; NI = menet·áram."""
    rho = np.asarray(rho, float); z = np.asarray(z, float)
    rho = np.where(np.abs(rho - a) < 1e-9 * a, a * (1 + 1e-9), rho)
    b = ell / 2; n = NI / ell; B0 = MU0 * n / np.pi
    out_r = 0.0; out_z = 0.0
    for sgn, zeta in ((+1, z + b), (-1, z - b)):
        den = np.sqrt(zeta ** 2 + (rho + a) ** 2)
        alpha = a / den; beta = zeta / den; gam = (a - rho) / (a + rho)
        kk = np.sqrt((zeta ** 2 + (a - rho) ** 2) / (zeta ** 2 + (a + rho) ** 2))
        out_r = out_r + sgn * alpha * cel(kk, 1.0, 1.0, -1.0)
        out_z = out_z + sgn * beta * cel(kk, gam ** 2, 1.0, gam)
    return B0 * out_r, B0 * a / (a + rho) * out_z

def loop_B(rho, z, a, I):
    """Egyetlen köráram tere (klasszikus elliptikus alak)."""
    rho = np.asarray(rho, float) + 1e-15; z = np.asarray(z, float)
    m = 4 * a * rho / ((a + rho) ** 2 + z ** 2)
    K, E = ellipk(m), ellipe(m); C = MU0 * I / (2 * np.pi)
    d = np.sqrt((a + rho) ** 2 + z ** 2); q = (a - rho) ** 2 + z ** 2
    Bz = C / d * (K + (a * a - rho ** 2 - z ** 2) / q * E)
    Br = C * z / (rho * d) * (-K + (a * a + rho ** 2 + z ** 2) / q * E)
    return Br, Bz

def coils_B(P, centers, a, ell, NI, senses):
    """Függőleges tengelyű tekercsek összege. P: (...,3) pontok; centers: (n,3) tekercsközéppontok."""
    P = np.asarray(P, float); B = np.zeros_like(P)
    for (cx, cy, cz), s in zip(centers, senses):
        dx, dy, dz = P[..., 0] - cx, P[..., 1] - cy, P[..., 2] - cz
        r = np.hypot(dx, dy)
        Br, Bz = solenoid_B(r, dz, a, ell, NI * s)
        with np.errstate(invalid="ignore", divide="ignore"):
            ux = np.where(r > 0, dx / np.where(r > 0, r, 1), 0.0); uy = np.where(r > 0, dy / np.where(r > 0, r, 1), 0.0)
        B[..., 0] += Br * ux; B[..., 1] += Br * uy; B[..., 2] += Bz
    return B

def nagaoka(a, ell):
    """Nagaoka-együttható (áramlap): k² = 4a²/(4a²+ℓ²)."""
    k2 = 4 * a * a / (4 * a * a + ell * ell); k = np.sqrt(k2); kp = np.sqrt(1 - k2)
    K, E = ellipk(k2), ellipe(k2)
    return 4 / (3 * np.pi * kp) * ((kp * kp / k2) * (K - E) + E - k)
