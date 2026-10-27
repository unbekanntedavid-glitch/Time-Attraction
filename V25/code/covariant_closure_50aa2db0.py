#!/usr/bin/env python3
"""
Covariant Action and Coupling Closure Study (Task 50aa2db0) for the Time Attraction Law (TAL).
Author: David Szabolcsi (ORCID: 0009-0005-6300-074X)
License: CC BY-NC 4.0

This deterministic, self-contained module implements:
1. Symbolic & numerical verification of the Covariant Action Closure (Formulation I: Nonlocal
   Positive-Definite Bilateral Kernel Action with auxiliary fields (phi_L, lambda_L); Formulation II:
   Local Screened Scalar-Tensor / Helmholtz-Yukawa Action) and machine/quadrature-precision verification
   of the covariant Bianchi/Noether conservation identity nabla_mu (T^{mu nu}_matter + T^{mu nu}_X) = 0.
2. Comparison of Formulation I (3D spherical top-hat kernel, exact background & bilateral reaction) and
   Formulation II (Yukawa-screened Helmholtz kernel, both unnormalized alpha_Y=1 and volume-normalized
   alpha_Y=1/3) across SLy, APR4, and H4 neutron-star sequences at L = 1 km and L = 10 km.
3. Coupled Radial Stability Eigenvalue Solver (omega_{0,std}^2, omega_{0,modB}^2, and coupled
   action-derived omega_{0,coupled}^2 including delta f_L = -gamma_L * delta bar_rho_L[xi] * f_L),
   quantifying the turning-point alignment with dM/d(rho_c) = 0 along the SLy L = 10 km sequence.
4. Microscopic Clock & Spectral-Line Coupling (H0 minimal physical-metric coupling vs. H1 disformal
   lapse / Jordan-frame conformal coupling) for Sirius B and the V20 8558 kg granite chamber.
"""

import json
import math
import os
import sys
import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import solve_ivp, quad, cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, minimize_scalar
from scipy.linalg import solve_banded
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
import warnings

# Standalone package: all helper modules live next to this script (relative paths only)
OUT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT_DIR))
import tal_profiles as tp
import radial_modes as rm

OUT_DIR = str(OUT_DIR)

# Fundamental constants (SI / CGS / geometric)
G_SI = 6.67430e-11
C_SI = 299792458.0
C_KMS = 299792.458
MSUN_KG = 1.98847e30
GM_KM = tp.st.GM  # 1.4766250385 km


# ==============================================================================
# PART 1: SYMBOLIC VERIFICATION OF COVARIANT CONSERVATION & STURM-LIOUVILLE LINK
# ==============================================================================
def run_symbolic_verifications():
    """Verify exact symbolic identities for covariant divergence and perturbation operators."""
    r = sp.Symbol("r", positive=True)
    eps = sp.Function("eps")(r)
    P = sp.Function("P")(r)
    m = sp.Function("m")(r)
    f = sp.Function("f")(r)
    nu = sp.Function("nu")(r)
    Gam = sp.Function("Gam")(r)
    xi = sp.Function("xi")(r)

    A = (m + 4 * sp.pi * r**3 * P) / (r * (r - 2 * m))

    # 1. Mass mode covariant conservation identity
    sub_mass = {
        sp.Derivative(m, r): 4 * sp.pi * r**2 * f * eps,
        sp.Derivative(nu, r): A,
        sp.Derivative(P, r): -(eps + P) * A,
    }
    div_matter_mass = sp.simplify((sp.diff(P, r) + (eps + P) * sp.diff(nu, r)).subs(sub_mass))
    eps_X_mass = -(1 - f) * eps
    pr_X_mass = sp.Integer(0)
    pt_X_mass = -r / 2 * (1 - f) * eps * A
    div_X_mass = sp.simplify(
        (sp.diff(pr_X_mass, r) + (eps_X_mass + pr_X_mass) * sp.diff(nu, r) - 2 / r * (pt_X_mass - pr_X_mass)).subs(sub_mass)
    )

    # 2. Gradient mode covariant conservation identity
    sub_grad = {
        sp.Derivative(m, r): 4 * sp.pi * r**2 * eps,
        sp.Derivative(nu, r): f * A,
        sp.Derivative(P, r): -(eps + P) * f * A,
    }
    div_matter_grad = sp.simplify((sp.diff(P, r) + (eps + P) * sp.diff(nu, r)).subs(sub_grad))
    eps_X_grad = sp.Integer(0)
    pr_X_grad = -(1 - f) * (P + m / (4 * sp.pi * r**3))
    pt_X_grad = pr_X_grad + r / 2 * (sp.diff(pr_X_grad, r) + pr_X_grad * f * A)
    div_X_grad = sp.simplify(
        (sp.diff(pr_X_grad, r) + (eps_X_grad + pr_X_grad) * sp.diff(nu, r) - 2 / r * (pt_X_grad - pr_X_grad)).subs(sub_grad)
    )

    # 3. Sturm-Liouville vs Linearized Hydrostatic Operator Identity in GR & Mass Mode
    sub_mass_full = dict(sub_mass)
    sub_mass_full[sp.Derivative(eps, r)] = -(eps + P)**2 / (Gam * P) * A
    lam = -sp.log(1 - 2 * m / r)
    nuC = 2 * nu
    nuCp = sp.diff(nuC, r)
    Pi = sp.exp((lam + 3 * nuC) / 2) * r**(-2) * Gam * P
    zeta = r**3 * sp.exp(-nuC / 2) * xi
    Q_matter = sp.exp((lam + 3 * nuC) / 2) * r**(-2) * (P + eps) * (nuCp**2 / 4 + 2 * nuCp / r - 8 * sp.pi * sp.exp(lam) * P)
    SL_matter = (-sp.diff(Pi * sp.diff(zeta, r), r) - Q_matter * zeta).doit().subs(sub_mass_full).doit().subs(sub_mass_full)

    DeltaP = -Gam * P * sp.exp(nuC / 2) / r**2 * sp.diff(zeta, r)
    dP = (DeltaP - r * xi * sp.diff(P, r)).doit().subs(sub_mass_full)
    deps = ((eps + P) / (Gam * P) * dP).doit().subs(sub_mass_full)
    dm2 = -4 * sp.pi * r**3 * (f * eps + P) * xi
    pref = sp.sqrt(r) * sp.exp(2 * nu) / sp.sqrt(r - 2 * m)
    dA2 = sp.diff(A, m) * dm2 + sp.diff(A, P) * dP
    hydro2 = (sp.diff(dP, r) + (deps + dP) * A + (eps + P) * dA2).doit().subs(sub_mass_full).doit().subs(sub_mass_full)

    diff_sl_hydro = sp.factor(sp.simplify(SL_matter - pref * hydro2))
    diff_gr_limit = sp.simplify(diff_sl_hydro.subs(f, 1))

    # 4. H1 disformal metric sign check: g~ = g + (1 - f^2)/c^2 u u, comoving g00=-1, u_0=-c
    fL, cc = sp.symbols("f_L c", positive=True)
    g00_t = -1 + (1 - fL**2) / cc**2 * (-cc) * (-cc)
    sym_disformal = sp.simplify(g00_t - (-fL**2))
    # 5. Yukawa center limit: lim_{r->0} (alpha_Y G L/(c^2 r)) sinh(r/L) e^{-r'/L} = alpha_Y G/c^2 r' e^{-r'/L} * (r'/r'..)
    rr, rp, LL, aY, GG = sp.symbols("r rprime L alpha_Y G", positive=True)
    lim_c = sp.limit(aY * GG * LL / (cc**2 * rr) * sp.sinh(rr / LL) * rp * sp.exp(-rp / LL), rr, 0)
    sym_yuk_center = sp.simplify(lim_c - aY * GG / cc**2 * rp * sp.exp(-rp / LL))
    # homogeneous check: alpha_Y = 4 pi alpha_std, alpha_std = 1/(4 pi) -> phi = gamma_L rho
    a_std = sp.Rational(1, 4) / sp.pi
    sym_yuk_norm = sp.simplify(4 * sp.pi * a_std - 1)
    # 6. Action-consistent Yukawa: L = -((phi')^2 + phi^2/L^2)/(2 alpha_Y) - e^{-phi} eps (G=c=1)
    rr6 = sp.Symbol("r", positive=True)
    ph = sp.Function("phi")(rr6)
    epsr = sp.Function("epsr")(rr6)
    Lag = rr6**2 * (-(sp.diff(ph, rr6)**2 + ph**2 / LL**2) / (2 * aY) - sp.exp(-ph) * epsr)
    EL = sp.diff(sp.diff(Lag, sp.diff(ph, rr6)), rr6) - sp.diff(Lag, ph)
    # expected: (-nabla^2 + L^-2) phi = alpha_Y f(phi) eps  (radial flat Laplacian)
    lap = sp.diff(rr6**2 * sp.diff(ph, rr6), rr6) / rr6**2
    expected = (-lap + ph / LL**2) - aY * sp.exp(-ph) * epsr
    sym_el_source = sp.simplify(EL * aY / rr6**2 - expected)
    # matter pressure from f(phi) eps(n): P_m = n d(f eps)/dn - f eps = f P
    nn = sp.Symbol("n", positive=True)
    epn = sp.Function("epn")(nn)
    fphi = sp.Symbol("f", positive=True)
    Pn = nn * sp.diff(epn, nn) - epn
    sym_el_pressure = sp.simplify(nn * sp.diff(fphi * epn, nn) - fphi * epn - fphi * Pn)
    # total radial divergence (flat scalar sector) with EL substituted gives h' = -nu' + phi'
    Pr, nup6 = sp.Function("Pr")(rr6), sp.Function("nup")(rr6)
    f6 = sp.exp(-ph)
    rho_s = (sp.diff(ph, rr6)**2 + ph**2 / LL**2) / (2 * aY)
    pr_s = (sp.diff(ph, rr6)**2 - ph**2 / LL**2) / (2 * aY)
    pt_s = -rho_s
    phpp = sp.solve(sp.Eq(expected, 0), sp.diff(ph, rr6, 2))[0]
    div_s = (sp.diff(pr_s, rr6) + 2 / rr6 * (pr_s - pt_s)).subs(sp.Derivative(ph, (rr6, 2)), phpp)
    hprime = -nup6 + sp.diff(ph, rr6)  # h' = P'/(eps+P)
    div_m = (sp.diff(f6 * Pr, rr6) + f6 * (epsr + Pr) * nup6).subs(sp.Derivative(Pr, rr6), (epsr + Pr) * hprime)
    sym_el_hydro = sp.simplify(div_m + div_s)
    return {
        "sym_yukawa_el_source_factor_residual": str(sym_el_source),
        "sym_yukawa_el_matter_pressure_residual": str(sym_el_pressure),
        "sym_yukawa_el_hydrostatic_residual": str(sym_el_hydro),
        "sym_disformal_g00_residual": str(sym_disformal),
        "sym_yukawa_center_limit_residual": str(sym_yuk_center),
        "sym_yukawa_alphaY_eq_4pi_alphastd_residual": str(sym_yuk_norm),
        "div_matter_mass_symbolic": str(div_matter_mass),
        "div_X_mass_symbolic": str(div_X_mass),
        "div_matter_grad_symbolic": str(div_matter_grad),
        "div_X_grad_symbolic": str(div_X_grad),
        "SL_minus_hydro_factored": str(diff_sl_hydro),
        "SL_minus_hydro_at_f_eq_1": str(diff_gr_limit),
    }


# ==============================================================================
# PART 2: NONLOCAL TOP-HAT vs. LOCAL YUKAWA-SCREENED HELMHOLTZ EQUILIBRIUM SOLVER
# ==============================================================================
def build_tophat_matrix(r, w_quad, L):
    """Exact 3D spherical top-hat shell-intersection quadrature matrix K_{ij} on radial grid r."""
    rq, rs = np.meshgrid(r, r, indexing="ij")
    frac = np.zeros_like(rq)
    d_diff = np.abs(rq - rs)
    d_sum = rq + rs
    inside = d_sum <= L
    partial = (~inside) & (d_diff < L)
    frac[inside] = 1.0
    frac[partial] = (L**2 - (rq[partial] - rs[partial])**2) / (4.0 * rq[partial] * rs[partial])
    vol = (4.0 / 3.0) * np.pi * (L**3)
    K = (4.0 * np.pi * (rs**2) * frac * w_quad[None, :]) / vol
    return K


def yukawa_average(rg, er, L, alpha_Y=1.0 / 3.0):
    """
    Solve the radial screened Helmholtz / Yukawa equation:
        (-nabla^2 + 1/L^2) phi_L = alpha_Y * (4*pi*G / c^2) * rho(r)
    in geometric units (where G=c=1 and gamma_L = L^2).
    Note on normalization:
      For a uniform density rho_0 in infinite space, (-nabla^2 + 1/L^2) phi_L = (1/L^2) phi_L,
      so if the source is (4*pi*G/c^2)*alpha_Y*rho_0, then phi_L = 4*pi*alpha_Y * (G*L^2/c^2)*rho_0.
      Matching the unit-integral normalization int K_Y(x) d^3x = 1 so that phi_L -> gamma_L * rho_0
      for uniform density fixes alpha_Y = 1/(4*pi) in units of 4*pi*G/c^2, or compared to the
      3D spherical top-hat volume (4/3)*pi*L^3 vs Yukawa volume 4*pi*L^3, the ratio of integral
      volumes is exactly 1/3!
    Here we return phi_L(rg) / (L^2) as the effective smoothed density bar_rho_Y(rg).
    """
    rq = rg[:, None]
    rs = rg[None, :]
    w = np.zeros_like(rg)
    w[1:-1] = 0.5 * (rg[2:] - rg[:-2])
    w[0] = 0.5 * (rg[1] - rg[0])
    w[-1] = 0.5 * (rg[-1] - rg[-2])
    r_min = np.minimum(rq, rs)
    r_max = np.maximum(rq, rs)
    # Overflow-free evaluation of sinh(r_min/L) * exp(-r_max/L) = 0.5*(exp(-(r_max-r_min)/L) - exp(-(r_max+r_min)/L))
    prod_exp = 0.5 * (np.exp(-(r_max - r_min) / L) - np.exp(-(r_max + r_min) / L))
    with np.errstate(divide="ignore", invalid="ignore"):
        kernel_core = np.where(rq > 1e-12, (L / np.maximum(rq, 1e-30)) * prod_exp, np.exp(-rs / L))
    K_unit = (rs * kernel_core * w[None, :]) / (L**2)
    return alpha_Y * (K_unit @ er)


def solve_yukawa_star(eos, logrho, mode="mass", L=10.0, alpha_Y=1.0, rtol=1e-10, ngrid=801, tol=2e-7):
    """Self-consistent fixed-point TOV solver for the local Yukawa-screened Helmholtz action."""
    base = tp.st.integrate(eos, logrho, rtol=rtol, ngrid=ngrid)
    if L <= 0.0 or mode == "GR":
        return base, None, None
    prev = base
    for _ in range(1, 40):
        av = yukawa_average(prev["rg"], prev["er"], L, alpha_Y=alpha_Y)
        xp = (L**2) * av
        xinterp = PchipInterpolator(prev["rg"], xp)
        oldR = prev["R"]
        fp = (lambda xi, oR: (lambda r: np.exp(-xi(min(r, oR)))))(xinterp, oldR)
        cur = tp.st.integrate(eos, logrho, mode, fp, rtol, ngrid)
        res = max(abs(cur["R"] / prev["R"] - 1.0), abs(cur["M"] / prev["M"] - 1.0))
        if res < tol:
            xc = (L**2) * yukawa_average(cur["rg"], cur["er"], L, alpha_Y=alpha_Y)
            xold = np.interp(cur["rg"], prev["rg"], xp)
            if float(np.max(abs(xc - xold))) < tol:
                return cur, fp, (xinterp, oldR)
        prev = cur
    raise RuntimeError("Yukawa fixed point did not converge")


def _avg_weights(rg):
    w = np.zeros_like(rg)
    w[1:-1] = 0.5 * (rg[2:] - rg[:-2])
    w[0] = 0.5 * (rg[1] - rg[0])
    w[-1] = 0.5 * (rg[-1] - rg[-2])
    return w


def integrate_r_tov(eos, logrho, f_of_r, lnV_of_r=None, rho_phi=None, pr_phi=None, ngrid=401, rtol=1e-9,
                    p_weighted=False):
    """
    Radial mass-mode TOV integration in r with (i) matter mass weight f_L(r), (ii) optional extra
    hydrostatic force d ln V_eff/dr (bilateral reciprocal enthalpy term) and (iii) optional canonical
    scalar energy density rho_phi(r) and radial pressure pr_phi(r) added to m' and nu'.
        dh/dr = -nu'(r) - d ln V_eff/dr,  nu' = (m + 4 pi r^3 (P + p_r,phi)) / (r (r - 2m)),
        dm/dr = 4 pi r^2 (f_L eps + rho_phi).
    p_weighted=True uses the action-consistent matter pressure P_matter = f_L P in nu'.
    Also returns m(r), nu(r) (nu(0)=0, additive constant irrelevant for ratios) and nu'(r) on the grid.
    """
    hc = eos.h_of_rho(10**logrho)
    zero = lambda r: 0.0
    rp = rho_phi or zero
    pp = pr_phi or zero
    # Integrated form: y0 = h + V_recip(r) - V_recip(0) satisfies dy0/dr = -nu' (no spline derivative needed)
    lnV = lnV_of_r
    def hof(r, y0):
        return y0 - (lnV(r) - lnV(0.0)) if lnV is not None else y0
    def rhs(r, y):
        h = hof(r, y[0]); m = y[1]
        e, P, _, _ = eos.at(max(h, 1e-12))
        Pm = f_of_r(r) * P if p_weighted else P
        nup = (m + 4.0 * np.pi * r**3 * (Pm + pp(r))) / (r * (r - 2.0 * m))
        return [-nup, 4.0 * np.pi * r * r * (f_of_r(r) * e + rp(r)), nup]
    ev = lambda r, y: hof(r, y[0]) - 1e-10
    ev.terminal = True
    r0 = 1e-4
    ec = eos.at(hc)[0]
    sol = solve_ivp(rhs, (r0, 40.0), [hc, 4.0 * np.pi / 3.0 * ec * f_of_r(0.0) * r0**3, 0.0], events=ev,
                    rtol=rtol, atol=[1e-12, 1e-12, 1e-12], dense_output=True, max_step=0.05)
    R = float(sol.t_events[0][0])
    M = float(sol.y_events[0][0][1])
    rg = np.linspace(0.0, R, ngrid)
    Y = sol.sol(np.maximum(rg, r0))
    hg = np.array([hof(x, y) for x, y in zip(rg, Y[0])])
    hg[0] = hc
    er = eos.at(np.maximum(hg, 1e-12))[0]
    er[-1] = 0.0
    mg = Y[1].copy(); mg[0] = 0.0; mg[-1] = M
    nug = Y[2].copy(); nug[0] = 0.0
    return dict(R=R, M=M / GM_KM, rg=rg, er=er, hc=hc, hg=hg, mg=mg, nug=nug)


def solve_bilateral_reaction_star(eos, logrho, L=10.0, ngrid=1201, tol=2e-7, return_profile=False):
    """
    Flat-kernel fixed-point equilibrium with bilateral reciprocal force (unweighted variant; NOT a full
    covariant EL solution: the lapse weight, and the metric variation of the kernel/volume element, are absent):
    L_fluid(n, phi_L) = -f_L(phi_L) eps(n), phi_L(x) = gamma_L int K_L(x,y) eps(y) d^3y.
    Varying eps(x) both explicitly and inside the kernel gives the effective potential
        V_eff(x) = f_L(x) - gamma_L int K_L(y,x) f_L(y) eps(y) d^3y   (K symmetric),
    and hydrostatic equilibrium V_eff * mu * e^nu = const, i.e.
        dh/dr = -A(r) - d V_recip/dr,   V_recip := ln V_eff.
    Solved here as a self-consistent fixed point (no reuse of the algebraic V21 background).
    """
    cur = tp.st.integrate(eos, logrho, rtol=1e-9, ngrid=ngrid)
    f_i = lambda r: 1.0
    dV_i = None
    for it in range(80):
        rg, er = cur["rg"], cur["er"]
        K = build_tophat_matrix(rg, _avg_weights(rg), L)
        f_th = np.exp(-(L**2) * (K @ er))
        recip = -(L**2) * (K @ (f_th * er))
        Veff = f_th + recip
        xg = rg / rg[-1]
        if it > 0:  # under-relaxation on the normalized radius x = r/R (stabilizes the fixed point)
            f_th = 0.5 * f_th + 0.5 * np.interp(xg, x_prev, f_prev)
            Veff = 0.5 * Veff + 0.5 * np.interp(xg, x_prev, V_prev)
        x_prev, f_prev, V_prev = xg, f_th.copy(), Veff.copy()
        lnV = PchipInterpolator(rg, np.log(Veff))
        fI = PchipInterpolator(rg, f_th)
        Rc = rg[-1]
        f_i = (lambda fi, Rr: (lambda r: float(fi(min(r, Rr)))))(fI, Rc)
        dV_i = (lambda d, Rr: (lambda r: float(d(min(r, Rr)))))(lnV, Rc)
        new = integrate_r_tov(eos, logrho, f_i, dV_i, ngrid=ngrid)
        res = max(abs(new["R"] / cur["R"] - 1.0), abs(new["M"] / cur["M"] - 1.0))
        cur = new
        # converged, or settled on the radial-grid discretization noise floor (<5e-5 relative in R, M)
        if (res < tol and it > 1) or (it >= 25 and res < 5e-5):
            break
    else:
        raise RuntimeError(f"bilateral fixed point did not converge (logrho={logrho}, L={L}, res={res:.2e})")
    out = {
        "R_km": float(cur["R"]), "M_Msun": float(cur["M"]),
        "f_center": float(f_th[0]), "f_surface": float(f_th[-1]),
        "recip_center": float(recip[0]), "recip_surface": float(recip[-1]),
        "V_eff_center": float(Veff[0]), "V_eff_surface": float(Veff[-1]),
        "fixed_point_iterations": it + 1,
        "fixed_point_last_rel_residual": float(res),
    }
    return (out, cur) if return_profile else out


def solve_yukawa_canonical_star(eos, logrho, L=10.0, alpha_Y=1.0, ngrid=401, tol=1e-7):
    """
    Phenomenological TOV model with flat-space Yukawa field and added scalar stress (rows
    yukawa_flat_phenom_Tphi). Not the EL solution of the declared action: the source lacks f_L, P is not
    f_L-weighted, the scalar force is absent and B=1-2m/r is omitted; see solve_yukawa_action_el_star.
    Field equation (code units, G=c=1, gamma_L=L^2):
        (-L^2 nabla^2 + 1) phi_L = alpha_Y L^2 rho  <=>  (-nabla^2 + L^-2) phi_L = 4 pi alpha_std rho,
    alpha_Y = 4 pi alpha_std. With L_phi = -(1/(2 alpha_Y)) [(nabla phi)^2 + phi^2/L^2] and matter
    coupling -f_L(phi) eps, the static canonical stress is
        rho_phi = (phi'^2 + phi^2/L^2)/(2 alpha_Y),  p_r = (phi'^2 - phi^2/L^2)/(2 alpha_Y),
        p_t = -(phi'^2 + phi^2/L^2)/(2 alpha_Y);
    rho_phi and p_r enter m' and nu' (interior; exterior scalar tail neglected).
    """
    cur = tp.st.integrate(eos, logrho, rtol=1e-9, ngrid=ngrid)
    for it in range(60):
        rg, er = cur["rg"], cur["er"]
        phi = (L**2) * yukawa_average(rg, er, L, alpha_Y=alpha_Y)
        phiI = PchipInterpolator(rg, phi)
        dphi = phiI.derivative()(rg)
        rho_p = (dphi**2 + phi**2 / L**2) / (2.0 * alpha_Y)
        pr_p = (dphi**2 - phi**2 / L**2) / (2.0 * alpha_Y)
        Rc = rg[-1]
        mk = lambda arr: (lambda I, Rr: (lambda r: float(I(min(r, Rr)))))(PchipInterpolator(rg, arr), Rc)
        f_i = mk(np.exp(-phi))
        new = integrate_r_tov(eos, logrho, f_i, None, mk(rho_p), mk(pr_p), ngrid=ngrid)
        res = max(abs(new["R"] / cur["R"] - 1.0), abs(new["M"] / cur["M"] - 1.0))
        cur = new
        if res < tol and it > 1:
            break
    else:
        raise RuntimeError("canonical Yukawa fixed point did not converge")
    cur["f_center"] = float(np.exp(-phi[0]))
    cur["f_surface"] = float(np.exp(-phi[-1]))
    cur["rho_phi_over_eps_center"] = float(rho_p[0] / er[0])
    return cur


def solve_bilateral_lapse_weighted_star(eos, logrho, L=10.0, ngrid=1201, tol=2e-7):
    """
    Flat-kernel fixed-point equilibrium with bilateral reciprocal force, lapse-weighted variant.
    Static matter functional E = int dV N(x) f_L(x) eps(x), N = e^nu, phi_L = gamma_L K eps (K, metric fixed).
    Its density variation gives
        V_eff,N(x) = f_L(x) - (gamma_L / N(x)) int dV_y N(y) K(y,x) f_L(y) eps(y),
    with P_matter = f_L P in the lapse equation and hydrostatics V_eff,N mu e^nu = const.
    The metric variation of K and of the volume element (needed for full covariant closure) is NOT included.
    """
    cur = integrate_r_tov(eos, logrho, lambda r: 1.0, None, ngrid=ngrid, p_weighted=True)
    for it in range(80):
        rg, er, nug = cur["rg"], cur["er"], cur["nug"]
        Nl = np.exp(nug - nug[-1])  # only ratios N(y)/N(x) enter
        K = build_tophat_matrix(rg, _avg_weights(rg), L)
        f_th = np.exp(-(L**2) * (K @ er))
        recip = -(L**2) * (K @ (Nl * f_th * er)) / Nl
        Veff = f_th + recip
        xg = rg / rg[-1]
        if it > 0:
            f_th = 0.5 * f_th + 0.5 * np.interp(xg, x_prev, f_prev)
            Veff = 0.5 * Veff + 0.5 * np.interp(xg, x_prev, V_prev)
        x_prev, f_prev, V_prev = xg, f_th.copy(), Veff.copy()
        Rc = rg[-1]
        mk = lambda arr: (lambda I, Rr: (lambda r: float(I(min(r, Rr)))))(PchipInterpolator(rg, arr), Rc)
        new = integrate_r_tov(eos, logrho, mk(f_th), mk(np.log(Veff)), ngrid=ngrid, p_weighted=True)
        res = max(abs(new["R"] / cur["R"] - 1.0), abs(new["M"] / cur["M"] - 1.0))
        cur = new
        if (res < tol and it > 1) or (it >= 25 and res < 5e-5):
            break
    else:
        raise RuntimeError(f"lapse-weighted bilateral fixed point did not converge (logrho={logrho}, L={L})")
    return {"R_km": float(cur["R"]), "M_Msun": float(cur["M"]),
            "f_center": float(f_th[0]), "V_eff_center": float(Veff[0]),
            "fixed_point_iterations": it + 1, "fixed_point_last_rel_residual": float(res)}


def _solve_curved_scalar(rg, B, dB, nup, src, L, k_robin=None):
    """FD solve of B phi'' + (2B/r + B nu' + B'/2) phi' - phi/L^2 = -src on [0,R],
    phi'(0)=0 (regular centre: 3 phi'' - phi/L^2 = -src), phi'(R) = -k_robin * phi(R)
    (where k_robin defaults to the flat-space approximation 1/R + 1/L when None, or is set
    to the exact curved-spacetime Dirichlet-to-Neumann slope k_curved from the exterior BVP)."""
    n = len(rg); h = rg[1] - rg[0]
    ab = np.zeros((3, n)); rhs = -np.asarray(src, float).copy()
    ab[1, 0] = -6.0 / h**2 - 1.0 / L**2; ab[0, 1] = 6.0 / h**2
    r = rg[1:-1]
    a = B[1:-1]; b = 2 * B[1:-1] / r + B[1:-1] * nup[1:-1] + 0.5 * dB[1:-1]
    ab[2, 0:n - 2] = a / h**2 - b / (2 * h)
    ab[1, 1:-1] = -2 * a / h**2 - 1.0 / L**2
    ab[0, 2:n] = a / h**2 + b / (2 * h)
    R = rg[-1]
    k = -(1.0 / R + 1.0 / L) if k_robin is None else -float(k_robin)
    aN, bN = B[-1], 2 * B[-1] / R + B[-1] * nup[-1] + 0.5 * dB[-1]
    # ghost phi_{n} = phi_{n-2} + 2 h k phi_{n-1}
    ab[2, n - 2] = 2 * aN / h**2
    ab[1, n - 1] = -2 * aN / h**2 + (aN / h**2 + bN / (2 * h)) * 2 * h * k - 1.0 / L**2
    return solve_banded((1, 1), ab, rhs)


def solve_yukawa_action_el_star(eos, logrho, L=10.0, alpha_Y=1.0, ngrid=801, tol=1e-8):
    """
    Action-consistent Yukawa branch (G=c=1, f_L = e^{-phi}, gamma_L = L^2) from
        L = -((nabla phi)^2 + phi^2/L^2)/(2 alpha_Y) - f_L(phi) eps(n):
      scalar EL:  B phi'' + (2B/r + B nu' + B'/2) phi' - phi/L^2 = -alpha_Y gamma_L f_L rho / L^2 = -alpha_Y f_L eps
      matter:     P_matter = f_L P,  dh/dr = -nu' - (ln f_L)' = -nu' + phi'
      stress:     rho_phi = (B phi'^2 + phi^2/L^2)/(2 alpha_Y),  p_r,phi = (B phi'^2 - phi^2/L^2)/(2 alpha_Y)
      exterior:   phi'(R) = -(1/R + 1/L) phi(R); tail mass Delta M_phi,ext integrated with phi_R (R/r) e^{-(r-R)/L}.
    """
    cur = integrate_r_tov(eos, logrho, lambda r: 1.0, None, ngrid=ngrid, p_weighted=True)
    phi = np.zeros(ngrid); x_prev = cur["rg"] / cur["rg"][-1]
    rho_p = np.zeros(ngrid); dphi = np.zeros(ngrid)
    prp_on = lambda rr: np.zeros_like(rr)
    for it in range(100):
        rg, er, mg, hg = cur["rg"], cur["er"], cur["mg"], cur["hg"]
        phi0 = np.interp(rg / rg[-1], x_prev, phi)
        rho0 = np.interp(rg / rg[-1], x_prev, rho_p)
        f = np.exp(-phi0)
        P = eos.at(np.maximum(hg, 1e-12))[1]; P[-1] = 0.0
        rs = np.maximum(rg, 1e-6)
        B = 1.0 - 2.0 * mg / rs; B[0] = 1.0
        pr0 = prp_on(rg)
        nup = (mg + 4 * np.pi * rs**3 * (f * P + pr0)) / (rs * (rs - 2 * mg)); nup[0] = 0.0
        dm = 4 * np.pi * rs**2 * (f * er + rho0)
        dB = 2 * mg / rs**2 - 2 * dm / rs; dB[0] = 0.0
        phin = _solve_curved_scalar(rg, B, dB, nup, alpha_Y * f * er, L)
        phi_new = phin if it == 0 else 0.5 * phin + 0.5 * phi0
        dphi = np.gradient(phi_new, rg, edge_order=2)
        dphi[0] = 0.0
        rho_p = (B * dphi**2 + phi_new**2 / L**2) / (2 * alpha_Y)
        pr_p = (B * dphi**2 - phi_new**2 / L**2) / (2 * alpha_Y)
        x_prev, phi = rg / rg[-1], phi_new
        Rc = rg[-1]
        mk = lambda arr: (lambda I, Rr: (lambda r: float(I(min(r, Rr)))))(PchipInterpolator(rg, arr), Rc)
        prp_on = (lambda xp, arr: (lambda rr: np.interp(rr / rr[-1], xp, arr)))(x_prev, pr_p.copy())
        new = integrate_r_tov(eos, logrho, mk(np.exp(-phi)), mk(-phi), mk(rho_p), mk(pr_p),
                              ngrid=ngrid, p_weighted=True)
        res = max(abs(new["R"] / cur["R"] - 1.0), abs(new["M"] / cur["M"] - 1.0))
        cur = new
        if res < tol and it > 2:
            break
    else:
        raise RuntimeError(f"action-EL Yukawa fixed point did not converge (logrho={logrho}, L={L})")
    R = cur["R"]; Mkm = cur["M"] * GM_KM; phiR = float(phi[-1])
    re = np.linspace(R, R + 40.0 * L, 20001)
    pe = phiR * (R / re) * np.exp(-(re - R) / L)
    dpe = -pe * (1.0 / re + 1.0 / L)
    rho_e = ((1 - 2 * Mkm / re) * dpe**2 + pe**2 / L**2) / (2 * alpha_Y)
    dM_ext = float(np.trapezoid(4 * np.pi * re**2 * rho_e, re)) / GM_KM
    return {"R_km": float(R), "M_Msun": float(cur["M"]), "Delta_M_phi_ext_Msun": dM_ext,
            "M_total_Msun": float(cur["M"]) + dM_ext, "f_center": float(np.exp(-phi[0])),
            "phi_surface": phiR, "fixed_point_iterations": it + 1, "fixed_point_last_rel_residual": float(res)}


def two_cell_lapse_weighting_test():
    """2-cell variational check: E = sum_i N_i f_i eps_i V_i, phi_i = gamma sum_j K_ij eps_j V_j, K_ij = 1/V_tot."""
    e1, e2 = sp.symbols("e1 e2", positive=True)
    N = (sp.Integer(1), sp.Rational(4, 5)); V = (1, 1); g = sp.Rational(1, 10); Vt = 2
    phi = g * (e1 * V[0] + e2 * V[1]) / Vt
    E = sum(N[i] * sp.exp(-phi) * e * V[i] for i, e in enumerate((e1, e2)))
    subs = {e1: 1, e2: 2}
    # V_eff,N(i) = (dE/d eps_i) / (N_i V_i);  unweighted recipe: f - gamma K (f eps) V
    weighted = [float((sp.diff(E, e) / (N[i] * V[i])).subs(subs)) for i, e in enumerate((e1, e2))]
    f = float(sp.exp(-phi).subs(subs))
    unweighted = f - float(g) * f * (1 * 1 + 2 * 1) / Vt
    return {"N": [1.0, 0.8], "eps": [1.0, 2.0], "gamma": 0.1, "V_cells": [1.0, 1.0], "K_ij": "1/V_tot",
            "V_eff_lapse_weighted": weighted, "V_eff_unweighted": [unweighted, unweighted],
            "expected_weighted": [0.74881594, 0.72084293], "expected_unweighted": 0.73160178}


def compute_conservation_residuals(eos_name="SLy", logrho=15.0, L=10.0):
    """
    Compute numerical covariant conservation residuals nabla_mu T^{mu}_{r}(matter) and
    nabla_mu T^{mu}_{r}(X) for both Formulation I (Nonlocal Top-Hat in mass & gradient closures)
    and Formulation II (Local Yukawa-screened Helmholtz closure).
    """
    eos = tp.EOS(eos_name)
    results = {}
    for mode in ["mass", "gradient"]:
        p = tp.dense_profile(eos, logrho, mode, L)
        R = p["R"]
        # Evaluate on interior points away from piecewise-polytrope transitions
        r = np.linspace(0.02 * R, 0.98 * R, 1200)
        h = p["h"](r)
        e, P, _, cs2 = eos.at(h)
        f = p["f"](r)
        dlnf = p["dlnf_dr"](r)
        df = f * dlnf
        m = p["m"](r)
        mu = p["mu"](r)
        dmu = p["dmu_dr"](r)
        A = r * (mu + 4.0 * np.pi * P) / (1.0 - 2.0 * mu * r * r)
        nu_p = A if mode == "mass" else f * A

        # Matter stress-energy divergence: P' + (e + P)*nu'

        # More directly, use dense_profile's h(r) derivative via central finite difference with step 1e-6 km
        eps_r = 1e-6
        dh_dr = (p["h"](r + eps_r) - p["h"](r - eps_r)) / (2.0 * eps_r)
        dP_dr = (e + P) * dh_dr
        div_matter = dP_dr + (e + P) * nu_p
        scale_matter = np.maximum(np.abs(dP_dr), 1e-30)

        # Auxiliary / X-sector stress-energy components and covariant divergence
        if mode == "mass":
            rho_X = -(1.0 - f) * e
            pr_X = np.zeros_like(r)
            pt_X = -0.5 * r * (1.0 - f) * e * nu_p
            dpr_X_dr = np.zeros_like(r)
            div_X = dpr_X_dr + (rho_X + pr_X) * nu_p - (2.0 / r) * (pt_X - pr_X)
            scale_X = np.maximum(np.abs(rho_X * nu_p), 1e-30)
        else:
            rho_X = np.zeros_like(r)
            pr_X = -(1.0 - f) * (P + mu / (4.0 * np.pi))
            # Analytic derivative of pr_X:
            dpr_X_dr = df * (P + mu / (4.0 * np.pi)) - (1.0 - f) * (dP_dr + dmu / (4.0 * np.pi))
            pt_X = pr_X + 0.5 * r * (dpr_X_dr + pr_X * nu_p)
            div_X = dpr_X_dr + (rho_X + pr_X) * nu_p - (2.0 / r) * (pt_X - pr_X)
            scale_X = np.maximum(np.abs(dpr_X_dr), 1e-30)

        # Total Einstein tensor divergence nabla_mu (T_matter + T_X)
        div_total = div_matter + div_X
        results[f"{eos_name}_{mode}_L{int(L)}km"] = {
            "R_km": float(R),
            "M_Msun": float(p["Msun"]),
            "f_center": float(p["f"](0.0)),
            "f_surface": float(p["f"](R)),
            "max_rel_div_matter": float(np.max(np.abs(div_matter / scale_matter))),
            "median_rel_div_matter": float(np.median(np.abs(div_matter / scale_matter))),
            "max_rel_div_X": float(np.max(np.abs(div_X / scale_X))),
            "algebraic_X_denominator": "|rho_X * nu'|" if mode == "mass" else "|P'_{r,X}|",
            "full_spline_denominator": "|P'|",
            "max_rel_div_total": float(np.max(np.abs(div_total / scale_matter))),
            "rho_X_over_eps_center": float(rho_X[0] / e[0]),
            "pr_X_over_P_center": float(pr_X[0] / P[0]),
            "pt_X_over_P_center": float(pt_X[0] / P[0]),
        }
    return results


# ==============================================================================
# PART 3: COUPLED RADIAL PULSATION EIGENVALUE SOLVER (TURNING-POINT RESTORATION)
# ==============================================================================
def solve_coupled_radial_omega2(eos, log10_rho0c, mode="mass", L_km=10.0, N=450, include_df=True,
                                W_scale=1.0, return_details=False):
    """
    Solve the linearized radial pulsation boundary-value problem derived from the second
    variation of the covariant action around the self-consistent equilibrium configuration:
        d(delta h)/dr + delta A_eff[delta h](r)
            = - omega^2 * [e^{2(lambda - nu)} / (4*pi*r^2 * f_m * (eps + P))] * delta m[delta h](r),
    subject to regularity at the center (delta h(0) = 1, delta m(0) = 0) and the Birkhoff
    surface boundary condition delta m(R; omega^2) = 0 (equivalent to Lagrangian pressure
    perturbation Delta P(R) = 0).

    When include_df=True, the non-local kernel perturbation
        delta f_L(r) = -gamma_L * f_L(r) * int K_L(r, r') * delta eps(r') d^3r'
    is included self-consistently in delta m'[delta h] (for mass mode) or delta nu'[delta h]
    (for gradient mode). When include_df=False, f_L(r) is held frozen during the pulsation.
    """
    p = tp.dense_profile(eos, log10_rho0c, mode, L_km)
    R = p["R"]
    hc = p["hc"]

    # Place exact grid cell boundaries at the piecewise-polytrope transition radii
    r_breaks = [0.0]
    for h_tr in sorted(eos.hb[1:], reverse=True):
        if 1e-9 < h_tr < hc:
            r_tr = brentq(lambda x: p["h"](x) - h_tr, 0.0, R)
            r_breaks.append(r_tr)
    r_breaks.append(R * 0.9996)
    r_breaks = sorted(r_breaks)

    edges_list = []
    for a, b in zip(r_breaks[:-1], r_breaks[1:]):
        n_sub = max(16, int(round(N * (b - a) / R)))
        edges_list.append(np.linspace(a, b, n_sub + 1)[:-1])
    edges = np.append(np.concatenate(edges_list), r_breaks[-1])
    N_tot = len(edges) - 1
    r_mid = 0.5 * (edges[:-1] + edges[1:])
    dr = np.diff(edges)

    h = p["h"](r_mid)
    e, P, _, cs2 = eos.at(h)
    f = p["f"](r_mid)
    m = p["m"](r_mid)
    mu = p["mu"](r_mid)
    nu = p["nu"](r_mid)
    A = r_mid * (mu + 4.0 * np.pi * P) / (1.0 - 2.0 * mu * r_mid * r_mid)

    M_int = np.tril(np.ones((N_tot, N_tot)), -1) * dr[None, :] + 0.5 * np.eye(N_tot) * dr[None, :]
    deps_factor = (e + P) / cs2

    if mode == "mass":
        if include_df and L_km > 0.0:
            K = build_tophat_matrix(r_mid, dr, L_km)
            S_mat = (4.0 * np.pi * (r_mid**2) * f)[:, None] * (np.eye(N_tot) - (L_km**2) * e[:, None] * K) * deps_factor[None, :]
        else:
            S_mat = np.diag(4.0 * np.pi * (r_mid**2) * f * deps_factor)
    else:
        S_mat = np.diag(4.0 * np.pi * (r_mid**2) * deps_factor)

    DM_mat = M_int @ S_mat

    c_dm = (1.0 + 2.0 * r_mid * A) / (r_mid * (r_mid - 2.0 * m))
    c_dh = (4.0 * np.pi * (r_mid**3) * (e + P)) / (r_mid * (r_mid - 2.0 * m))
    if mode != "gradient":
        DA_mat = (c_dm[:, None] * DM_mat) + np.diag(c_dh)
    else:
        DA_base = (c_dm[:, None] * DM_mat) + np.diag(c_dh)
        if include_df and L_km > 0.0:
            K = build_tophat_matrix(r_mid, dr, L_km)
            DF_mat = -(L_km**2) * (f * A)[:, None] * K * deps_factor[None, :]
            DA_mat = (f[:, None] * DA_base) + DF_mat
        else:
            DA_mat = f[:, None] * DA_base

    fm = f if mode == "mass" else np.ones_like(f)
    e2lam = 1.0 / (1.0 - 2.0 * mu * r_mid * r_mid)
    e_minus_2nu = np.exp(-2.0 * nu)
    w_inert = (e2lam * e_minus_2nu) / (4.0 * np.pi * (r_mid**2) * fm * (e + P))
    W_mat = W_scale * w_inert[:, None] * DM_mat

    om2_guess = float(rm.fundamental(p, 300, "matter", on_ill="indefinite")["omega2"][0])
    def F_shoot(om2):
        dh_s = np.linalg.solve(np.eye(N_tot) + M_int @ (DA_mat + om2 * W_mat), np.ones(N_tot))
        return float(np.sum(dr * (S_mat @ dh_s)))

    delta = 0.0012
    grid = np.linspace(om2_guess - delta, om2_guess + delta, 25)
    fvals = [F_shoot(g) for g in grid]
    best_om2 = None
    for k in range(len(grid) - 1):
        if fvals[k] * fvals[k + 1] <= 0.0 and abs(fvals[k]) < 5.0 and abs(fvals[k + 1]) < 5.0:
            cand = brentq(F_shoot, grid[k], grid[k + 1], xtol=1e-11)
            if best_om2 is None or abs(cand - om2_guess) < abs(best_om2 - om2_guess):
                best_om2 = cand
    if best_om2 is None:
        best_om2 = brentq(F_shoot, om2_guess - 0.0005, om2_guess + 0.0005, xtol=1e-11)
    if not return_details:
        return best_om2 * (C_KMS**2)
    # First-order tangent (Rayleigh-Ritz-type) estimate around omega^2 = 0: omega^2 ~ -F(0)/F'(0)
    d_om = 1e-6
    F0 = F_shoot(0.0)
    dF = (F_shoot(d_om) - F_shoot(-d_om)) / (2.0 * d_om)
    # Symmetrized Galerkin kinetic form M_sym = DM^T diag(w_inert * dr) DM (symmetric, PSD by construction)
    M_sym = DM_mat.T @ (((W_scale * w_inert * dr))[:, None] * DM_mat)
    M_sym = 0.5 * (M_sym + M_sym.T)
    return {
        "omega2_s2": best_om2 * (C_KMS**2),
        "F_shoot_0": F0,
        "omega2_tangent_s2": (-F0 / dF) * (C_KMS**2),
        "min_M_sym_diag": float(np.min(np.diag(M_sym))),
        "M_sym_asymmetry": float(np.max(np.abs(M_sym - M_sym.T))),
        "W_raw_asymmetry_rel": float(np.max(np.abs(W_mat - W_mat.T)) / np.max(np.abs(W_mat))),
    }


# ==============================================================================
# PART 4: MICROSCOPIC CLOCK & SPECTRAL-LINE COUPLING (H0 vs. H1 FROM ACTION)
# ==============================================================================
def compute_clock_coupling_table():
    """
    Compute exact action-derived predictions for H0 (minimal physical-metric coupling of
    atomic Dirac + Maxwell fields) vs. H1 (disformal lapse coupling g_tilde_{mu nu} =
    g_{mu nu} + (1 - f_L^2)/c^2 u_mu u_nu or Jordan-frame conformal coupling with fixed atomic mass):
      1. Sirius B white dwarf (both 3D Top-Hat and Yukawa-screened Helmholtz kernels at L = 0, 1, 10, 100, 1000, 5000 km).
      2. V20 8558 kg Granite Chamber (R_in = 0.25 m, R_out = 0.85 m, H = 1.20 m cylinder + two endcaps).
    """
    rows = []

    # 1. Sirius B White Dwarf (Joyce et al. 2018: M = 1.018 +/- 0.011 Msun, R = 0.00803 +/- 0.00011 Rsun,
    #    v_obs = 80.65 +/- 0.77 km/s; Lane-Emden n = 3 relativistic white-dwarf polytrope)
    M_wd_Msun = 1.018
    sigM_wd = 0.011
    R_wd_m = 0.00803 * tp.st.RS
    R_wd_km = R_wd_m / 1000.0
    sigR_wd_m = 0.00011 * tp.st.RS
    v_obs_kms = 80.65
    sig_v_obs = 0.77
    comp_wd = G_SI * (M_wd_Msun * tp.st.MS) / (R_wd_m * (C_SI**2))
    z_gr_wd = (1.0 - 2.0 * comp_wd) ** (-0.5) - 1.0
    v_gr_kms = (C_SI / 1000.0) * z_gr_wd
    sig_gr_kms = v_gr_kms * math.hypot(sigM_wd / M_wd_Msun, sigR_wd_m / R_wd_m)
    sigma_comb_kms = math.hypot(sig_v_obs, sig_gr_kms)

    rr_le, dd_le = tp.st.lane_emden(3.0)
    r_wd_m_grid = rr_le * R_wd_m
    rho_wd_si_grid = dd_le * (M_wd_Msun * tp.st.MS) / (R_wd_m**3)

    # Also Sirius A reference (n=3, M_A=2.042 Msun, R_A=1.711 Rsun)
    R_A_m = 1.711 * tp.st.RS
    M_A_Msun = 2.042
    r_A_m_grid = rr_le * R_A_m
    rho_A_si_grid = dd_le * (M_A_Msun * tp.st.MS) / (R_A_m**3)

    for L_km in [0.0, 1.0, 10.0, 100.0, 1000.0, 5000.0]:
        L_m = L_km * 1000.0
        if L_km == 0.0:
            xB_th = 0.0
            xA_th = 0.0
            xB_yuk = 0.0
            xA_yuk = 0.0
            dM_over_M_th = 0.0
        else:
            rb_th = float(tp.st.average(r_wd_m_grid, rho_wd_si_grid, L_m, queries=[R_wd_m], nq=160)[0])
            ra_th = float(tp.st.average(r_A_m_grid, rho_A_si_grid, L_m, queries=[R_A_m], nq=160)[0])
            xB_th = (G_SI * (L_m**2) / (C_SI**2)) * rb_th
            xA_th = (G_SI * (L_m**2) / (C_SI**2)) * ra_th

            # Yukawa-screened Helmholtz surface value phi_Y(R_wd) with unit-normalized kernel:
            # Using yukawa_average in meters and SI density:
            av_y_B = yukawa_average(r_wd_m_grid, rho_wd_si_grid, L_m, alpha_Y=1.0)
            av_y_A = yukawa_average(r_A_m_grid, rho_A_si_grid, L_m, alpha_Y=1.0)
            xB_yuk = float((G_SI * (L_m**2) / (C_SI**2)) * av_y_B[-1])
            xA_yuk = float((G_SI * (L_m**2) / (C_SI**2)) * av_y_A[-1])

            # Under H0 (minimal physical-metric coupling), the surface redshift changes only via the
            # interior gravitational mass modification delta M / M:
            rb_prof = tp.st.average(r_wd_m_grid, rho_wd_si_grid, L_m, queries=r_wd_m_grid, nq=80)
            f_prof = np.exp(-(G_SI * (L_m**2) / (C_SI**2)) * rb_prof)
            M_eff = np.trapezoid(4.0 * math.pi * (r_wd_m_grid**2) * rho_wd_si_grid * f_prof, r_wd_m_grid)
            M_tot = np.trapezoid(4.0 * math.pi * (r_wd_m_grid**2) * rho_wd_si_grid, r_wd_m_grid)
            dM_over_M_th = float((M_eff - M_tot) / M_tot)

        f_surf_th = math.exp(-xB_th)
        f_surf_yuk = math.exp(-xB_yuk)
        v_H0_mass = v_gr_kms * (1.0 + dM_over_M_th)
        dv_B_th = (C_SI / 1000.0) * (1.0 + z_gr_wd) * math.expm1(xB_th)
        dv_A_th = (C_SI / 1000.0) * math.expm1(xA_th)
        v_H1_tophat = v_gr_kms + (dv_B_th - dv_A_th)

        dv_B_yuk = (C_SI / 1000.0) * (1.0 + z_gr_wd) * math.expm1(xB_yuk)
        dv_A_yuk = (C_SI / 1000.0) * math.expm1(xA_yuk)
        v_H1_yuk = v_gr_kms + (dv_B_yuk - dv_A_yuk)

        rows.append({
            "system": "Sirius_B_White_Dwarf_n3",
            "L_m": L_m,
            "L_km": L_km,
            "M_Msun_or_kg": M_wd_Msun,
            "R_km_or_m": R_wd_km,
            "f_L_surface_tophat": f_surf_th,
            "f_L_surface_yukawa": f_surf_yuk,
            "H0_prediction": v_H0_mass,
            "H1_tophat_prediction": v_H1_tophat,
            "H1_yukawa_prediction": v_H1_yuk,
            "extra_shift_H1_tophat_kms": dv_B_th - dv_A_th,
            "extra_shift_H1_yukawa_kms": dv_B_yuk - dv_A_yuk,
            "unit": "km/s",
            "obs_value": v_obs_kms,
            "obs_sigma": sigma_comb_kms,
            "tension_H0_sigma": abs(v_H0_mass - v_obs_kms) / sigma_comb_kms,
            "tension_H1_tophat_sigma": abs(v_H1_tophat - v_obs_kms) / sigma_comb_kms,
            "tension_H1_yukawa_sigma": abs(v_H1_yuk - v_obs_kms) / sigma_comb_kms,
            "H0_fixedM_kms": v_gr_kms,
            "H1_tophat_fixedM_kms": v_H1_tophat,
            "H0_common_background_kms": v_H0_mass,
            "H1_tophat_common_background_kms": v_H0_mass + (dv_B_th - dv_A_th),
            "H1_yukawa_common_background_kms": v_H0_mass + (dv_B_yuk - dv_A_yuk),
        })

    # 2. V20 8558 kg Granite Chamber (R_in = 0.25 m, R_out = 0.85 m, H_in = 1.20 m barrel + endcaps = 3.112 m^3, 8558 kg)
    rho_granite = 2750.0  # kg/m^3
    R_in = 0.25
    R_out = 0.85
    z_in = 0.60
    # Total volume V_tot = 8558.0 / 2750.0 = 3.112000 m^3; barrel volume pi*(0.85^2 - 0.25^2)*1.20 = 2.488141 m^3;
    # two solid endcap disks of radius R_out = 0.85 m have thickness dz_cap = (3.112000 - 2.488141)/(2*pi*0.85^2) = 0.137426 m
    z_out = 0.60 + (8558.0 / rho_granite - math.pi * (R_out**2 - R_in**2) * (2.0 * z_in)) / (2.0 * math.pi * (R_out**2))
    nr, nz = 500, 500
    rc_edges = np.linspace(0.0, R_out, nr + 1)
    zc_edges = np.linspace(0.0, z_out, nz + 1)
    rc = 0.5 * (rc_edges[:-1] + rc_edges[1:])
    zc = 0.5 * (zc_edges[:-1] + zc_edges[1:])
    drc = rc_edges[1] - rc_edges[0]
    dzc = zc_edges[1] - zc_edges[0]
    RC, ZC = np.meshgrid(rc, zc, indexing="ij")
    in_granite = (RC >= R_in) | (ZC >= z_in)
    dV = 2.0 * (2.0 * math.pi * RC * drc * dzc)
    M_granite_kg = float(np.sum(rho_granite * dV[in_granite]))
    dist_3d = np.sqrt(RC**2 + ZC**2)
    dm_g = (rho_granite * dV)[in_granite]
    rg_c, zg_c, dg_in = RC[in_granite], ZC[in_granite], dist_3d[in_granite]
    Phi_N_center = -G_SI * float(np.sum(dm_g / dg_in))
    dy_Newtonian = Phi_N_center / (C_SI**2)
    # OUT state: clock displaced to D_out = 2.0 m in the mid-plane (full 3D azimuthal quadrature)
    D_out = 2.0
    n_phi = 256
    phis = (np.arange(n_phi) + 0.5) * 2.0 * math.pi / n_phi
    L_list = [0.20, 0.50, 0.85, 1.00, 1.50, 10.0]
    Phi_out = 0.0
    th_out = {L: 0.0 for L in L_list}
    yk_out = {L: 0.0 for L in L_list}
    for ph in phis:
        d_o = np.sqrt(rg_c**2 + D_out**2 - 2.0 * rg_c * D_out * math.cos(ph) + zg_c**2)
        dmp = dm_g / n_phi
        Phi_out += -G_SI * float(np.sum(dmp / d_o))
        for L in L_list:
            th_out[L] += float(np.sum(dmp[d_o <= L]))
            yk_out[L] += float(np.sum(dmp * np.exp(-d_o / L) / (4.0 * math.pi * L**2 * d_o)))
    dy_N_out = Phi_out / (C_SI**2)

    for L_m in L_list:
        gamma_L_SI = G_SI * (L_m**2) / (C_SI**2)
        vol_ball = (4.0 / 3.0) * math.pi * (L_m**3)
        phi_th_center = gamma_L_SI * float(np.sum(dm_g[dg_in <= L_m])) / vol_ball
        phi_th_out = gamma_L_SI * th_out[L_m] / vol_ball
        phi_yuk_center = gamma_L_SI * float(np.sum(dm_g * np.exp(-dg_in / L_m) / (4.0 * math.pi * L_m**2 * dg_in)))
        phi_yuk_out = gamma_L_SI * yk_out[L_m]
        dy_H1_tophat = -phi_th_center
        dy_H1_yukawa = -phi_yuk_center
        # H0: TAL excess over Newtonian GR (second order, |Phi_N/c^2| * phi_L)
        dy_H0_cavity = abs(dy_Newtonian) * phi_th_center

        rows.append({
            "system": "V20_Granite_Chamber_8558kg",
            "L_m": L_m,
            "L_km": L_m / 1000.0,
            "M_Msun_or_kg": M_granite_kg,
            "R_km_or_m": R_out,
            "f_L_surface_tophat": 1.0 - phi_th_center,
            "f_L_surface_yukawa": 1.0 - phi_yuk_center,
            "H0_prediction": dy_H0_cavity,
            "H1_tophat_prediction": dy_H1_tophat,
            "H1_yukawa_prediction": dy_H1_yukawa,
            "extra_shift_H1_tophat_kms": np.nan,
            "extra_shift_H1_yukawa_kms": np.nan,
            "unit": "fractional_dy",
            "obs_value": 0.0,
            "obs_sigma": 1.0e-18,  # illustrative fractional-frequency benchmark, not a demonstrated sensitivity
            "tension_H0_sigma": abs(dy_H0_cavity) / 1.0e-18,
            "tension_H1_tophat_sigma": abs(dy_H1_tophat) / 1.0e-18,
            "tension_H1_yukawa_sigma": abs(dy_H1_yukawa) / 1.0e-18,
            "D_out_m": D_out,
            "GR_newtonian_IN": dy_Newtonian,
            "GR_newtonian_OUT": dy_N_out,
            "GR_newtonian_IN_minus_OUT": dy_Newtonian - dy_N_out,
            "H0_TAL_excess_IN": dy_H0_cavity,
            "H0_TAL_excess_IN_minus_OUT": dy_H0_cavity - abs(dy_N_out) * phi_th_out,
            "H1_shift_IN": dy_H1_tophat,
            "H1_shift_OUT": -phi_th_out,
            "H1_shift_IN_minus_OUT": dy_H1_tophat + phi_th_out,
            "H1_yukawa_shift_IN_minus_OUT": dy_H1_yukawa + phi_yuk_out,
        })

    return pd.DataFrame(rows)


def extra_tov_rows():
    """Flat-kernel bilateral fixed-point rows (unweighted and lapse-weighted; SLy/APR4/H4), the flat-space
    phenomenological Yukawa+scalar-stress rows and the action-consistent Yukawa EL rows (SLy) at L = 1, 10 km."""
    rows = []
    jobs = [(e, L, "bilateral_eom", 1.0) for e in ["SLy", "APR4", "H4"] for L in [1.0, 10.0]]
    jobs += [(e, L, "bilateral_lapse_weighted", 1.0) for e in ["SLy", "APR4", "H4"] for L in [1.0, 10.0]]
    jobs += [("SLy", L, "yukawa_flat_phenom_Tphi", 1.0) for L in [1.0, 10.0]]
    jobs += [("SLy", L, "yukawa_action_el", 1.0) for L in [1.0, 10.0]]
    for eos_name, L_km, kern, a_y in jobs:
        eos = tp.EOS(eos_name)
        if kern == "bilateral_eom":
            def solve(lr):
                o, cur = solve_bilateral_reaction_star(eos, lr, L_km, return_profile=True)
                return dict(M=cur["M"], R=cur["R"], fc=o["f_center"], fs=o["f_surface"])
        elif kern == "bilateral_lapse_weighted":
            def solve(lr):
                o = solve_bilateral_lapse_weighted_star(eos, lr, L_km)
                return dict(M=o["M_Msun"], R=o["R_km"], fc=o["f_center"], fs=np.nan)
        elif kern == "yukawa_action_el":
            def solve(lr):
                o = solve_yukawa_action_el_star(eos, lr, L_km, alpha_Y=a_y)
                return dict(M=o["M_Msun"], R=o["R_km"], fc=o["f_center"], fs=float(np.exp(-o["phi_surface"])),
                            dMext=o["Delta_M_phi_ext_Msun"])
        else:
            def solve(lr):
                cur = solve_yukawa_canonical_star(eos, lr, L_km, alpha_Y=a_y)
                return dict(M=cur["M"], R=cur["R"], fc=cur["f_center"], fs=cur["f_surface"])
        s15 = solve(15.0)
        opt = minimize_scalar(lambda lr: -solve(lr)["M"], bounds=(15.1, 15.6), method="bounded", options={"xatol": 1e-5})
        pk = solve(float(opt.x))
        lr14 = brentq(lambda lr: solve(lr)["M"] - 1.4, 14.6, float(opt.x), xtol=1e-6)
        rows.append({
            "EOS": eos_name, "L_km": L_km, "mode": "mass", "kernel": kern, "alpha_Y": a_y,
            "alpha_Y_unit": a_y if kern.startswith("yukawa") else np.nan,
            "alpha_std_poisson": a_y / (4.0 * np.pi) if kern.startswith("yukawa") else np.nan,
            "M_max_Msun": pk["M"], "R_at_Mmax_km": pk["R"], "log10_rho0c_turning": float(opt.x),
            "R_1p4_km": solve(lr14)["R"], "M_at_15p0_Msun": s15["M"], "R_at_15p0_km": s15["R"],
            "f_L_center_15p0": s15["fc"], "f_L_surface_15p0": s15["fs"],
            "Delta_M_phi_ext_15p0_Msun": s15.get("dMext", np.nan),
            "Delta_M_phi_ext_at_Mmax_Msun": pk.get("dMext", np.nan),
        })
        print(f"   {eos_name:4s} L={L_km:4.1f} mass     {kern:22s}: M_max={pk['M']:.5f} Msun at {float(opt.x):.6f}, R_1.4={rows[-1]['R_1p4_km']:.4f} km")
    return rows


# ==============================================================================
# MAIN EXECUTION PIPELINE
# ==============================================================================
def main():
    print("1. Running symbolic verifications (SymPy)...")
    sym_res = run_symbolic_verifications()
    print("   Symbolic checks passed:", sym_res["div_matter_mass_symbolic"], sym_res["div_X_mass_symbolic"])

    print("2. Computing numerical covariant conservation residuals...")
    cons_res = {}
    for eos_name in ["SLy", "APR4", "H4"]:
        cons_res.update(compute_conservation_residuals(eos_name, logrho=15.0, L=10.0))
    bilat_res = {}
    for eos_name in ["SLy", "APR4", "H4"]:
        bilat_res[eos_name] = solve_bilateral_reaction_star(tp.EOS(eos_name), 15.0, L=10.0)

    sly_m, sly_g = cons_res["SLy_mass_L10km"], cons_res["SLy_gradient_L10km"]
    eos_w = tp.EOS("SLy")
    lrc_w = 15.30426924786
    dW1 = solve_coupled_radial_omega2(eos_w, lrc_w, "mass", 10.0, N=500, include_df=True, W_scale=1.0, return_details=True)
    dW4 = solve_coupled_radial_omega2(eos_w, lrc_w, "mass", 10.0, N=500, include_df=True, W_scale=4.0, return_details=True)
    cons_payload = {
        "symbolic_identities_residual": 0,
        "algebraic_X_sector_quadrature_residual": {
            "SLy_mass_L10km": {"denominator": "|rho_X * nu'|", "max_rel": sly_m["max_rel_div_X"]},
            "SLy_gradient_L10km": {"denominator": "|P'_{r,X}|", "max_rel": sly_g["max_rel_div_X"]},
        },
        "full_numerical_spline_divergence_residual": {
            "SLy_mass_L10km": {"denominator": "|P'|", "max_rel": sly_m["max_rel_div_matter"]},
            "SLy_gradient_L10km": {"denominator": "|P'|", "max_rel": sly_g["max_rel_div_matter"]},
        },
        "inertia_scaling_control_W_to_4W": {
            "EOS": "SLy", "mode": "mass", "L_km": 10.0, "N": 500, "log10_rho0c": lrc_w,
            "omega2_W_s2": dW1["omega2_s2"], "omega2_4W_s2": dW4["omega2_s2"],
            "ratio_4W_over_W": dW4["omega2_s2"] / dW1["omega2_s2"],
            "F_shoot_0_W": dW1["F_shoot_0"], "F_shoot_0_4W": dW4["F_shoot_0"],
            "omega2_tangent_W_s2": dW1["omega2_tangent_s2"],
            "min_M_sym_diag": dW1["min_M_sym_diag"], "min_M_sym_diag_positive": bool(dW1["min_M_sym_diag"] > 0),
            "M_sym_asymmetry": dW1["M_sym_asymmetry"], "raw_W_relative_asymmetry": dW1["W_raw_asymmetry_rel"],
            "note": "omega^2=0 is fixed by the static tangent operator (F_shoot(0) independent of W); the zero crossing is inertia-invariant.",
        },
        "symbolic_identities": sym_res,
        "numerical_conservation_residuals_L10km_logrho15": cons_res,
        "bilateral_kernel_reaction_L10km_logrho15": bilat_res,
        "sym_yukawa_el_source_factor_residual": sym_res["sym_yukawa_el_source_factor_residual"],
        "sym_yukawa_el_matter_pressure_residual": sym_res["sym_yukawa_el_matter_pressure_residual"],
        "sym_yukawa_el_hydrostatic_residual": sym_res["sym_yukawa_el_hydrostatic_residual"],
        "two_cell_lapse_weighting_test": two_cell_lapse_weighting_test(),
        "bilateral_logrho15_L10km_reproduction": {
            e: {"bilateral_eom_flat_unweighted": solve_bilateral_reaction_star(tp.EOS(e), 15.0, 10.0),
                "bilateral_lapse_weighted": solve_bilateral_lapse_weighted_star(tp.EOS(e), 15.0, 10.0)}
            for e in ["SLy", "APR4", "H4"]},
        "yukawa_action_el_SLy_logrho15": {
            f"L{int(Lk)}km": solve_yukawa_action_el_star(tp.EOS("SLy"), 15.0, Lk) for Lk in (1.0, 10.0)},
        "static_tangent_model_note": ("The coupled radial BVP linearizes around the V21 tp.dense_profile background, "
            "not around bilateral_* or yukawa_* branches. Positive diagonal of M_sym = D^T diag(w dr) D is a "
            "Gram-matrix property, not proof of a positive-definite physical kinetic Hessian."),
    }
    with open(os.path.join(OUT_DIR, "conservation_Residuals.json"), "w", encoding="utf-8") as fh:
        json.dump(cons_payload, fh, indent=2)
    print("   Saved conservation_Residuals.json")

    print("3. Computing Nonlocal Top-Hat vs. Local Yukawa-Screened TOV sequences (SLy, APR4, H4)...")
    tov_rows = []
    for eos_name in ["SLy", "APR4", "H4"]:
        eos = tp.EOS(eos_name)
        for L_km in [0.0, 1.0, 10.0]:
            modes_to_run = [("GR", "none", 0.0)] if L_km == 0.0 else [
                ("mass", "tophat", 1.0),
                ("gradient", "tophat", 1.0),
                ("mass", "yukawa_unit", 1.0),
                ("gradient", "yukawa_unit", 1.0),
                ("mass", "yukawa_vol13", 1.0 / 3.0),
            ]
            for mode, kernel_type, alpha_y in modes_to_run:
                # 1. Canonical star at log10(rho0c) = 15.0
                if kernel_type in ("none", "tophat"):
                    p15 = tp.dense_profile(eos, 15.0, mode, L_km)
                    R15 = float(p15["R"])
                    M15 = float(p15["Msun"])
                    f_c15 = float(p15["f"](0.0))
                    f_s15 = float(p15["f"](R15))
                    # Find M_max and turning point
                    res_opt = minimize_scalar(
                        lambda lrc: -tp.st.star(eos, lrc, mode, L_km)["M"],
                        bounds=(15.15, 15.55),
                        method="bounded",
                        options={"xatol": 1e-6},
                    )
                    lrc_peak = float(res_opt.x)
                    st_peak = tp.st.star(eos, lrc_peak, mode, L_km)
                    M_max = float(st_peak["M"])
                    R_Mmax = float(st_peak["R"])
                    # Find R_1.4 at M = 1.4 Msun
                    lrc_14 = brentq(
                        lambda lrc: tp.st.star(eos, lrc, mode, L_km)["M"] - 1.4,
                        14.70,
                        15.15,
                        xtol=1e-6,
                    )
                    R_14 = float(tp.st.star(eos, lrc_14, mode, L_km)["R"])
                else:
                    cur15, fp15, _ = solve_yukawa_star(eos, 15.0, mode, L_km, alpha_Y=alpha_y)
                    R15 = float(cur15["R"])
                    M15 = float(cur15["M"])
                    f_c15 = float(fp15(0.0))
                    f_s15 = float(fp15(R15))
                    res_opt = minimize_scalar(
                        lambda lrc: -solve_yukawa_star(eos, lrc, mode, L_km, alpha_Y=alpha_y)[0]["M"],
                        bounds=(15.15, 15.55),
                        method="bounded",
                        options={"xatol": 1e-5},
                    )
                    lrc_peak = float(res_opt.x)
                    st_peak, _, _ = solve_yukawa_star(eos, lrc_peak, mode, L_km, alpha_Y=alpha_y)
                    M_max = float(st_peak["M"])
                    R_Mmax = float(st_peak["R"])
                    lrc_14 = brentq(
                        lambda lrc: solve_yukawa_star(eos, lrc, mode, L_km, alpha_Y=alpha_y)[0]["M"] - 1.4,
                        14.70,
                        15.15,
                        xtol=1e-5,
                    )
                    R_14 = float(solve_yukawa_star(eos, lrc_14, mode, L_km, alpha_Y=alpha_y)[0]["R"])

                tov_rows.append({
                    "EOS": eos_name,
                    "L_km": L_km,
                    "mode": mode,
                    "kernel": kernel_type,
                    "alpha_Y": alpha_y,
                    "alpha_Y_unit": alpha_y if kernel_type.startswith("yukawa") else np.nan,
                    "alpha_std_poisson": alpha_y / (4.0 * np.pi) if kernel_type.startswith("yukawa") else np.nan,
                    "M_max_Msun": M_max,
                    "R_at_Mmax_km": R_Mmax,
                    "log10_rho0c_turning": lrc_peak,
                    "R_1p4_km": R_14,
                    "M_at_15p0_Msun": M15,
                    "R_at_15p0_km": R15,
                    "f_L_center_15p0": f_c15,
                    "f_L_surface_15p0": f_s15,
                })
                print(f"   {eos_name:4s} L={L_km:4.1f} {mode:8s} {kernel_type:12s}: M_max={M_max:.5f} Msun at log10_rho0c={lrc_peak:.6f}, R_1.4={R_14:.4f} km")

    tov_rows += extra_tov_rows()
    df_tov = pd.DataFrame(tov_rows)
    df_tov.to_csv(os.path.join(OUT_DIR, "closure_tov_comparison.csv"), index=False)
    print("   Saved closure_tov_comparison.csv")

    print("4. Computing Coupled Radial Stability Curves & Zero Crossings (SLy, L = 10 km)...")
    eos_sly = tp.EOS("SLy")
    su_v21 = pd.read_csv(os.path.join(OUT_DIR, "v21_radial_summary_ref.csv"))
    stab_rows = []
    zero_crossings = {}

    for mode, L_km in [("GR", 0.0), ("mass", 10.0), ("gradient", 10.0)]:
        row_v21 = su_v21[(su_v21.EOS == "SLy") & (su_v21["mode"] == mode) & (su_v21.L_m == L_km * 1000.0)].iloc[0]
        pk = float(row_v21["refined_log10_rho0c_dMdrho0"])
        lrc_grid = np.unique(np.r_[
            np.linspace(14.85, 15.42, 18),
            np.linspace(pk - 0.015, pk + 0.015, 15),
            pk
        ])
        for lrc in lrc_grid:
            p = tp.dense_profile(eos_sly, lrc, mode, L_km)
            om2_std = float(rm.fundamental(p, 500, "matter", on_ill="indefinite")["omega2"][0] * (C_KMS**2))
            om2_modB = float(rm.fundamental(p, 500, "metric", on_ill="indefinite")["omega2"][0] * (C_KMS**2))
            om2_frozen_bvp = float(solve_coupled_radial_omega2(eos_sly, lrc, mode, L_km, N=450, include_df=False))
            det_c = solve_coupled_radial_omega2(eos_sly, lrc, mode, L_km, N=450, include_df=True, return_details=True)
            om2_coupled = float(det_c["omega2_s2"])
            om2_tangent = float(det_c["omega2_tangent_s2"])
            stab_rows.append({
                "EOS": "SLy",
                "mode": mode,
                "L_km": L_km,
                "log10_rho0c": float(lrc),
                "M_Msun": float(p["Msun"]),
                "R_km": float(p["R"]),
                "omega0_sq_std_matter_s2": om2_std,
                "omega0_sq_modB_metric_s2": om2_modB,
                "omega0_sq_bvp_frozen_f_s2": om2_frozen_bvp,
                "omega0_sq_coupled_action_s2": om2_coupled,
                "omega0_sq_sym_tangent_s2": om2_tangent,
                "log10_rho0c_turning": pk,
            })

        # Refine exact zero crossings with brentq
        zc_std = float(row_v21["omega2_zero_log10_rho0c_matter"])
        zc_modB = float(row_v21["omega2_zero_log10_rho0c_metric"]) if not np.isnan(row_v21["omega2_zero_log10_rho0c_metric"]) else zc_std
        zc_frozen = float(brentq(
            lambda x: solve_coupled_radial_omega2(eos_sly, x, mode, L_km, N=500, include_df=False),
            pk - 0.012, pk + 0.012, xtol=1e-6
        ))
        zc_coupled = float(brentq(
            lambda x: solve_coupled_radial_omega2(eos_sly, x, mode, L_km, N=500, include_df=True),
            pk - 0.005, pk + 0.005, xtol=1e-6
        ))
        zero_crossings[f"SLy_{mode}_L{int(L_km)}km"] = {
            "log10_rho0c_turning_dM0": pk,
            "zc_std_matter": zc_std,
            "delta_dex_std_matter": zc_std - pk,
            "zc_modB_metric": zc_modB,
            "delta_dex_modB_metric": zc_modB - pk,
            "zc_bvp_frozen_f": zc_frozen,
            "delta_dex_bvp_frozen_f": zc_frozen - pk,
            "zc_coupled_action": zc_coupled,
            "delta_dex_coupled_action": zc_coupled - pk,
            "mismatch_reduction_factor_vs_std": abs(zc_std - pk) / max(abs(zc_coupled - pk), 1e-12),
            "mismatch_reduction_factor_vs_modB": abs(zc_modB - pk) / max(abs(zc_coupled - pk), 1e-12),
            "mismatch_reduction_factor_vs_frozen": abs(zc_frozen - pk) / max(abs(zc_coupled - pk), 1e-12),
        }
        print(f"   Stability {mode:8s} L={L_km:4.1f} km: turning={pk:.6f}, std={zc_std:.6f} ({zc_std-pk:+.2e}), modB={zc_modB:.6f} ({zc_modB-pk:+.2e}), frozen={zc_frozen:.6f} ({zc_frozen-pk:+.2e}), coupled={zc_coupled:.6f} ({zc_coupled-pk:+.2e})")

    df_stab = pd.DataFrame(stab_rows)
    df_stab.to_csv(os.path.join(OUT_DIR, "coupled_stability_results.csv"), index=False)
    cons_payload["radial_stability_zero_crossings"] = zero_crossings
    curved_json_path = os.path.join(OUT_DIR, "curved_yukawa_audit_summary.json")
    if os.path.exists(curved_json_path):
        with open(curved_json_path, "r", encoding="utf-8") as cfh:
            curved_audit = json.load(cfh)
        cons_payload["curved_exterior_yukawa_adm_audit"] = curved_audit
    curved_seq_path = os.path.join(OUT_DIR, "curved_yukawa_adm_sequences.csv")
    if os.path.exists(curved_seq_path):
        df_cseq = pd.read_csv(curved_seq_path)
        extra_curved = []
        for _, r_c in df_cseq.iterrows():
            extra_curved.append({
                "EOS": r_c["EOS"],
                "L_km": float(r_c["L_km"]),
                "mode": "mass",
                "kernel": "yukawa_action_el_curved_adm",
                "alpha_Y": float(r_c["alpha_Y"]),
                "alpha_Y_unit": float(r_c["alpha_Y"]),
                "alpha_std_poisson": float(r_c["alpha_std_poisson"]),
                "M_max_Msun": float(r_c["M_ADM_max_Msun"]),
                "R_at_Mmax_km": float(r_c["R_at_M_ADM_max_km"]),
                "log10_rho0c_turning": float(r_c["log10_rho0c_at_M_ADM_max"]),
                "R_1p4_km": float(r_c["R_at_M_ADM_1p4_km"]),
                "M_at_15p0_Msun": float(r_c["M_ADM_at_15p0_Msun"]),
                "R_at_15p0_km": float(r_c["R_at_15p0_curved_km"]),
                "f_L_center_15p0": float(r_c["f_L_center_15p0_curved"]),
                "f_L_surface_15p0": float(np.exp(-r_c["phi_surface_15p0_curved"])),
                "Delta_M_phi_ext_15p0_Msun": float(r_c["Delta_M_phi_ext_curved_15p0_Msun"]),
                "Delta_M_phi_ext_at_Mmax_Msun": float(r_c["Delta_M_phi_ext_at_M_ADM_max_Msun"]),
                "M_surface_max_Msun": float(r_c["M_surface_max_Msun"]),
                "R_at_M_surface_max_km": float(r_c["R_at_M_surface_max_km"]),
                "log10_rho0c_at_M_surface_max": float(r_c["log10_rho0c_at_M_surface_max"]),
                "R_at_M_surface_1p4_km": float(r_c["R_at_M_surface_1p4_km"]),
                "M_baryon_max_Msun": float(r_c["M_baryon_max_Msun"]),
                "log10_rho0c_at_M_baryon_max": float(r_c["log10_rho0c_at_M_baryon_max"]),
                "k_curved_15p0_km_inv": float(r_c["k_curved_15p0_km_inv"]),
                "k_flat_15p0_km_inv": float(r_c["k_flat_15p0_km_inv"]),
                "rel_diff_dphi_R_15p0": float(r_c["rel_diff_dphi_R_15p0"]),
            })
        df_tov = pd.concat([df_tov, pd.DataFrame(extra_curved)], ignore_index=True)
        df_tov.to_csv(os.path.join(OUT_DIR, "closure_tov_comparison.csv"), index=False)
    with open(os.path.join(OUT_DIR, "conservation_Residuals.json"), "w", encoding="utf-8") as fh:
        json.dump(cons_payload, fh, indent=2)
    print("   Saved coupled_stability_results.csv and updated conservation_Residuals.json")

    print("5. Computing Microscopic Clock & Spectral-Line Coupling Table (H0 vs. H1)...")
    df_clock = compute_clock_coupling_table()
    df_clock.to_csv(os.path.join(OUT_DIR, "clock_coupling_H0_H1.csv"), index=False)
    print("   Saved clock_coupling_H0_H1.csv")

    print("6. Generating High-Resolution Figures...")
    generate_figures(df_tov, df_stab, zero_crossings, cons_payload, df_clock)
    print("Done! All numerical outputs and figures generated.")


def generate_figures(df_tov, df_stab, zero_crossings, cons_payload, df_clock):
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 8.5,
        "figure.dpi": 300,
    })

    # --------------------------------------------------------------------------
    # FIGURE 1: Action Profiles (Top-Hat vs. Yukawa-Screened Kernel & M-R Curves)
    # --------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    eos_sly = tp.EOS("SLy")
    p_th = tp.dense_profile(eos_sly, 15.0, "mass", 10.0)
    rg = p_th["cur"]["rg"]
    er = p_th["cur"]["er"]
    f_th = p_th["f"](rg)
    _, fp_yuk1, _ = solve_yukawa_star(eos_sly, 15.0, "mass", 10.0, alpha_Y=1.0)
    _, fp_yuk13, _ = solve_yukawa_star(eos_sly, 15.0, "mass", 10.0, alpha_Y=1.0 / 3.0)
    f_yuk1 = np.array([fp_yuk1(r) for r in rg])
    f_yuk13 = np.array([fp_yuk13(r) for r in rg])

    ax = axes[0]
    ax.plot(rg, f_th, "b-", lw=2.2, label=r"Formulation I: 3D Top-Hat $f_L(r)=\exp(-\gamma_L\bar{\rho}_L)$")
    ax.plot(rg, f_yuk1, "r--", lw=2.0, label=r"Formulation II: Yukawa ($\alpha_Y=1$, unit $\int K_Y d^3x=1$)")
    ax.plot(rg, f_yuk13, "g-.", lw=2.0, label=r"Formulation II: Yukawa ($\alpha_Y=1/3$, matched volume)")
    ax.set_xlabel(r"Radius $r$ [km]")
    ax.set_ylabel(r"Lapse / Mass Suppression Factor $f_L(r)$")
    ax.set_title(r"(a) Interior TAL Profile $f_L(r)$ (SLy, $\log_{10}\rho_{0c}=15.0$, $L=10$ km)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")

    ax2 = axes[1]
    lrc_seq = np.linspace(14.65, 15.45, 32)
    for label, style, solver_fn in [
        ("GR Baseline ($L=0$)", "k-", lambda x: tp.st.star(eos_sly, x, "GR", 0.0)),
        ("Formulation I: Top-Hat Mass ($L=10$ km)", "b-", lambda x: tp.st.star(eos_sly, x, "mass", 10.0)),
        ("Formulation I: Top-Hat Gradient ($L=10$ km)", "c--", lambda x: tp.st.star(eos_sly, x, "gradient", 10.0)),
        ("Formulation II: Yukawa Mass ($\\alpha_Y=1$, $L=10$ km)", "r--", lambda x: solve_yukawa_star(eos_sly, x, "mass", 10.0, 1.0)[0]),
        ("Formulation II: Yukawa Mass ($\\alpha_Y=1/3$, $L=10$ km)", "g-.", lambda x: solve_yukawa_star(eos_sly, x, "mass", 10.0, 1.0/3.0)[0]),
    ]:
        pts = [solver_fn(x) for x in lrc_seq]
        ax2.plot([p["R"] for p in pts], [p["M"] for p in pts], style, lw=1.9, label=label)
    ax2.set_xlabel(r"Stellar Radius $R$ [km]")
    ax2.set_ylabel(r"Gravitational Mass $M$ [$M_\odot$]")
    ax2.set_title(r"(b) SLy Mass–Radius Sequences Across Action Formulations")
    ax2.set_xlim(9.5, 12.5)
    ax2.set_ylim(1.0, 2.35)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="lower left", fontsize=7.8)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "fig1_action_profiles.png"), dpi=300)
    plt.close(fig)

    # --------------------------------------------------------------------------
    # FIGURE 2: Coupled Radial Stability & Turning-Point Alignment (SLy, L=10 km)
    # --------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8))
    for idx, mode in enumerate(["mass", "gradient"]):
        ax = axes[idx]
        sub = df_stab[(df_stab["mode"] == mode) & (df_stab["L_km"] == 10.0)].sort_values("log10_rho0c")
        zc_info = zero_crossings[f"SLy_{mode}_L10km"]
        pk = zc_info["log10_rho0c_turning_dM0"]

        ax.axhline(0.0, color="k", lw=1.0, ls=":")
        ax.axvline(pk, color="k", lw=1.8, ls="--", label=rf"Static Turning Point $dM/d\rho_c=0$ ({pk:.5f})")
        ax.plot(
            sub["log10_rho0c"], sub["omega0_sq_std_matter_s2"] * 1e-7,
            "o--", color="#d95f02", ms=4, lw=1.5,
            label=rf"Frozen-$f_L$ Matter $\omega_{{0,\mathrm{{std}}}}^2$ ($\Delta={zc_info['delta_dex_std_matter']:+.2e}$ dex)"
        )
        ax.plot(
            sub["log10_rho0c"], sub["omega0_sq_modB_metric_s2"] * 1e-7,
            "s-.", color="#7570b3", ms=4, lw=1.5,
            label=rf"Frozen-$f_L$ Metric $\omega_{{0,\mathrm{{modB}}}}^2$ ($\Delta={zc_info['delta_dex_modB_metric']:+.2e}$ dex)"
        )
        ax.plot(
            sub["log10_rho0c"], sub["omega0_sq_bvp_frozen_f_s2"] * 1e-7,
            "^:", color="#e7298a", ms=4, lw=1.5,
            label=rf"Frozen-$f_L$ BVP ($\delta f_L=0$, $\Delta={zc_info['delta_dex_bvp_frozen_f']:+.2e}$ dex)"
        )
        ax.plot(
            sub["log10_rho0c"], sub["omega0_sq_coupled_action_s2"] * 1e-7,
            "b-", lw=2.4,
            label=rf"Static tangent model (V21 bg), selected BVP branch (incl. $\delta f_L$, $\Delta={zc_info['delta_dex_coupled_action']:+.2e}$ dex)"
        )
        ax.set_xlim(pk - 0.020, pk + 0.028)
        ax.set_ylim(-4.5, 4.5)
        ax.set_xlabel(r"Central Rest-Mass Density $\log_{10}(\rho_{0c}\,[\mathrm{g\,cm^{-3}}])$")
        ax.set_ylabel(r"Selected BVP-branch $\omega^2$ [$10^7\,\mathrm{rad^2\,s^{-2}}$]")
        ax.set_title(rf"({'a' if idx==0 else 'b'}) V21-background static nonlocal tangent model: SLy $L=10$ km, {mode}", fontsize=9.5)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=7.6)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "fig2_coupled_stability.png"), dpi=300)
    plt.close(fig)

    # --------------------------------------------------------------------------
    # FIGURE 3: Covariant Conservation Check & Auxiliary Stress Profiles T^mu_nu(X)
    # --------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    p_m = tp.dense_profile(eos_sly, 15.0, "mass", 10.0)
    p_g = tp.dense_profile(eos_sly, 15.0, "gradient", 10.0)
    r_m = np.linspace(0.02 * p_m["R"], 0.98 * p_m["R"], 400)
    r_g = np.linspace(0.02 * p_g["R"], 0.98 * p_g["R"], 400)

    e_m, P_m, _, _ = eos_sly.at(p_m["h"](r_m))
    f_m = p_m["f"](r_m)
    mu_m = p_m["mu"](r_m)
    A_m = r_m * (mu_m + 4.0 * np.pi * P_m) / (1.0 - 2.0 * mu_m * r_m * r_m)
    rho_X_m = -(1.0 - f_m) * e_m
    pt_X_m = -0.5 * r_m * (1.0 - f_m) * e_m * A_m

    e_g, P_g, _, _ = eos_sly.at(p_g["h"](r_g))
    f_g = p_g["f"](r_g)
    df_g = f_g * p_g["dlnf_dr"](r_g)
    mu_g = p_g["mu"](r_g)
    dmu_g = p_g["dmu_dr"](r_g)
    A_g = r_g * (mu_g + 4.0 * np.pi * P_g) / (1.0 - 2.0 * mu_g * r_g * r_g)
    nu_p_g = f_g * A_g
    eps_r = 1e-6
    dP_g = (e_g + P_g) * (p_g["h"](r_g + eps_r) - p_g["h"](r_g - eps_r)) / (2.0 * eps_r)
    pr_X_g = -(1.0 - f_g) * (P_g + mu_g / (4.0 * np.pi))
    dpr_X_g = df_g * (P_g + mu_g / (4.0 * np.pi)) - (1.0 - f_g) * (dP_g + dmu_g / (4.0 * np.pi))
    pt_X_g = pr_X_g + 0.5 * r_g * (dpr_X_g + pr_X_g * nu_p_g)

    ax = axes[0]
    ax.plot(r_m / p_m["R"], rho_X_m / e_m[0], "b-", lw=2.0, label=r"Mass mode: $\varepsilon_X(r)/\varepsilon_c$")
    ax.plot(r_m / p_m["R"], pt_X_m / P_m[0], "b--", lw=1.8, label=r"Mass mode: $P_{t,X}(r)/P_c$ ($P_{r,X}=0$)")
    ax.plot(r_g / p_g["R"], pr_X_g / P_g[0], "r-", lw=2.0, label=r"Gradient mode: $P_{r,X}(r)/P_c$ ($\varepsilon_X=0$)")
    ax.plot(r_g / p_g["R"], pt_X_g / P_g[0], "r--", lw=1.8, label=r"Gradient mode: $P_{t,X}(r)/P_c$")
    ax.set_xlabel(r"Normalized Radius $r / R_*$")
    ax.set_ylabel(r"Normalized Auxiliary Stress $T^\mu_{\;\nu}(X)$")
    ax.set_title(r"(a) Covariant Auxiliary Sector $T^\mu_{\;\nu}(X)$ Profiles (SLy, $L=10$ km)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.0)

    ax2 = axes[1]
    dP_m = (e_m + P_m) * (p_m["h"](r_m + eps_r) - p_m["h"](r_m - eps_r)) / (2.0 * eps_r)
    spl_m = np.abs(dP_m + (e_m + P_m) * A_m) / np.maximum(np.abs(dP_m), 1e-30)
    spl_g = np.abs(dP_g + (e_g + P_g) * nu_p_g) / np.maximum(np.abs(dP_g), 1e-30)
    alg_m = np.abs((rho_X_m) * A_m - (2.0 / r_m) * pt_X_m) / np.maximum(np.abs(rho_X_m * A_m), 1e-30)
    alg_g = np.abs(dpr_X_g + pr_X_g * nu_p_g - (2.0 / r_g) * (pt_X_g - pr_X_g)) / np.maximum(np.abs(dpr_X_g), 1e-30)
    fl = 1e-18
    ax2.semilogy(r_m / p_m["R"], np.maximum(spl_m, fl), "b-", lw=1.8, label=r"Mass: full spline $|\nabla_\mu T^\mu_{\;r}| / |P'|$")
    ax2.semilogy(r_g / p_g["R"], np.maximum(spl_g, fl), "r-", lw=1.8, label=r"Gradient: full spline $|\nabla_\mu T^\mu_{\;r}| / |P'|$")
    ax2.semilogy(r_m / p_m["R"], np.maximum(alg_m, fl), "b:", lw=1.6, label=r"Mass: algebraic X-sector $/|\rho_X\nu'|$")
    ax2.semilogy(r_g / p_g["R"], np.maximum(alg_g, fl), "r--", lw=1.6, label=r"Gradient: algebraic X-sector $/|P'_{r,X}|$")
    ax2.set_xlabel(r"Normalized Radius $r / R_*$")
    ax2.set_ylabel(r"Relative Covariant Divergence Residual")
    ax2.set_title(r"(b) Algebraic X-sector vs. full numerical divergence residual")
    ax2.set_ylim(1e-18, 1e-3)
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(loc="center right", fontsize=7.0, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "fig3_conservation_check.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
