#!/usr/bin/env python3
"""
Curved-Spacetime Exterior Einstein-Scalar Solver, ADM Mass Optimization,
Baryon Mass Turning-Point Verification, and Coupled (xi, delta phi_L) Radial
Pulsation Stability Solver on the Action-Derived Yukawa Branch (yukawa_action_el).

Addresses Reviewer Points 1, 2, 3, and 4 and computes:
  1. True curved-spacetime exterior Einstein-scalar vacuum BVP on r in [R_*, R_max]
     with asymptotic curved WKB/Yukawa decay at R_max = R_* + N_mult * L, verified
     for R_max-independence across N_mult in {10, 12, 14, 16}.
  2. Exact curved Dirichlet-to-Neumann surface slope dphi_L/dr(R_*) vs. the
     flat-space Robin approximation -(1/R_* + 1/L)*phi_L(R_*).
  3. True asymptotic ADM mass M_ADM = M_surface + Delta M_{phi,ext}^{curved},
     alongside the historical surface mass M_surface = m(R_*) and the conserved
     baryon rest mass M_B = int_0^{R_*} 4*pi*r^2 * rho_0(r) / sqrt(1 - 2*m(r)/r) dr,
     re-optimizing M_max, R(M_max), log10(rho_0c), and R_{1.4} by M_ADM for
     SLy, APR4, and H4 at L = 1 km and 10 km (and alpha_Y in {1, 1/3}).
  4. Coupled (delta h, delta phi_L, delta m) radial pulsation eigenvalue solver
     and explicit Sturm-Liouville / Hessian coefficient evaluation (Q_{xi xi},
     Q_{phi phi}, C_{xi phi}) directly on the curved `yukawa_action_el` equilibrium
     sequence, comparing log10(rho_0c)(omega_0^2 = 0) against log10(rho_0c)(dM_ADM/drho_c = 0)
     and log10(rho_0c)(dM_B/drho_c = 0).
"""
import os
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_bvp
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve_banded
from scipy.optimize import minimize_scalar, brentq

import tal_profiles as tp
import covariant_closure_50aa2db0 as cc

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_V23_DIR = os.path.join(OUT_DIR, "data", "v23")
GM_KM = tp.st.GM
CONV = tp.st.CONV
C_KMS = 299792.458


def solve_exterior_4comp_bvp(R, M_surf_km, phi_R, L, alpha_Y=1.0, N_mult=14.0, n_pts=320, tol=1e-12):
    k_flat = 1.0 / R + (1.0 / L if L > 0 else 0.0)
    if abs(phi_R) < 1e-18 or L <= 0.0:
        nu_s = 0.5 * np.log(max(1.0 - 2.0 * M_surf_km / R, 1e-12))
        return {
            "k_curved": k_flat,
            "k_flat": k_flat,
            "rel_diff_k": 0.0,
            "dphi_R_curved": 0.0,
            "dphi_R_flat": 0.0,
            "M_ADM_km": M_surf_km,
            "Delta_M_phi_ext_Msun": 0.0,
            "nu_R": nu_s,
            "sol": None,
            "R_max": R + N_mult * max(L, 1.0),
        }
    R_max = R + N_mult * L
    r_grid = R + (R_max - R) * (np.linspace(0.0, 1.0, n_pts) ** 1.8)
    phi_init = phi_R * (R / r_grid) * np.exp(-(r_grid - R) / L)
    dphi_init = -phi_init * (1.0 / r_grid + 1.0 / L)
    m_init = np.full_like(r_grid, M_surf_km)
    nu_init = 0.5 * np.log(np.maximum(1.0 - 2.0 * m_init / r_grid, 1e-8))

    def fun(r, y):
        phi, dphi, m, nu = y
        B = np.maximum(1.0 - 2.0 * m / r, 1e-8)
        dm = (2.0 * np.pi * r * r / alpha_Y) * (B * dphi * dphi + (phi * phi) / (L * L))
        ddphi = (-(2.0 / r * (1.0 - m / r) - (4.0 * np.pi * r / (alpha_Y * L * L)) * phi * phi) * dphi + phi / (L * L)) / B
        pr_phi = (B * dphi * dphi - (phi * phi) / (L * L)) / (2.0 * alpha_Y)
        dnu = (m + 4.0 * np.pi * (r ** 3) * pr_phi) / (r * (r - 2.0 * m))
        return np.vstack([dphi, ddphi, dm, dnu])

    def bc(ya, yb):
        m_end = yb[2]
        B_end = max(1.0 - 2.0 * m_end / R_max, 1e-8)
        k_asympt = (1.0 / L) / np.sqrt(B_end) + (1.0 + m_end / L) / R_max
        nu_asympt = 0.5 * np.log(B_end)
        return np.array([
            ya[0] - phi_R,
            ya[2] - M_surf_km,
            yb[1] + k_asympt * yb[0],
            yb[3] - nu_asympt,
        ])

    sol = solve_bvp(fun, bc, r_grid, np.vstack([phi_init, dphi_init, m_init, nu_init]), tol=tol, max_nodes=25000)
    if not sol.success:
        raise RuntimeError(f"Curved exterior 4-component BVP did not converge: {sol.message}")
    dphi_R_curved = float(sol.y[1, 0])
    k_curved = -dphi_R_curved / phi_R
    M_ADM_km = float(sol.y[2, -1])
    nu_R = float(sol.y[3, 0])
    return {
        "k_curved": k_curved,
        "k_flat": k_flat,
        "rel_diff_k": (k_curved - k_flat) / k_flat,
        "dphi_R_curved": dphi_R_curved,
        "dphi_R_flat": -k_flat * phi_R,
        "M_ADM_km": M_ADM_km,
        "Delta_M_phi_ext_Msun": (M_ADM_km - M_surf_km) / GM_KM,
        "nu_R": nu_R,
        "r_ext": sol.x,
        "phi_ext": sol.y[0],
        "dphi_ext": sol.y[1],
        "m_ext": sol.y[2],
        "nu_ext": sol.y[3],
        "sol": sol,
        "R_max": R_max,
    }

def solve_exterior_dtN_curved(R, M_surf_km, phi_R, L, alpha_Y=1.0, N_mult=14.0, n_pts=320, tol=1e-12):
    """Curved exterior Dirichlet-to-Neumann map. Delegates to the 4-component exterior
    BVP (phi_L, psi_L, m, nu) with Schwarzschild tail nu(R_max) = 0.5*ln(1 - 2 m(R_max)/R_max)."""
    return solve_exterior_4comp_bvp(R, M_surf_km, phi_R, L, alpha_Y=alpha_Y, N_mult=N_mult, n_pts=n_pts, tol=tol)


def _solve_exterior_3comp_robin_iterate(R, M_surf_km, phi_R, L, alpha_Y=1.0, N_mult=14.0, n_pts=160):
    """
    Solve the coupled exterior Einstein-scalar vacuum equations on [R, R + N_mult * L]:
        B(r) phi_L'' + [2/r * (1 - m/r) - 4*pi*r/(alpha_Y * L^2) * phi_L^2] phi_L' - phi_L/L^2 = 0,
        m'(r) = (2*pi*r^2/alpha_Y) * [B(r) * (phi_L')^2 + phi_L^2/L^2],
    where B(r) = 1 - 2*m(r)/r, subject to phi_L(R) = phi_R, m(R) = M_surf_km,
    and the curved asymptotic WKB decay boundary condition at R_max = R + N_mult * L:
        phi_L'(R_max) + [(1/L)/sqrt(1 - 2*m(R_max)/R_max) + (1 + m(R_max)/L)/R_max] * phi_L(R_max) = 0.
    """
    k_flat = 1.0 / R + 1.0 / L
    if abs(phi_R) < 1e-18:
        return {
            "k_curved": k_flat,
            "k_flat": k_flat,
            "rel_diff_k": 0.0,
            "dphi_R_curved": 0.0,
            "dphi_R_flat": 0.0,
            "M_ADM_km": M_surf_km,
            "Delta_M_phi_ext_Msun": 0.0,
            "nu_R": 0.5 * np.log(max(1.0 - 2.0 * M_surf_km / R, 1e-12)),
        }
    R_max = R + N_mult * L
    r_grid = R + (R_max - R) * np.linspace(0.0, 1.0, n_pts) ** 1.8
    phi_init = phi_R * (R / r_grid) * np.exp(-(r_grid - R) / L)
    dphi_init = -phi_init * (1.0 / r_grid + 1.0 / L)
    m_init = np.full_like(r_grid, M_surf_km)

    def fun(r, y):
        phi, dphi, m = y
        B = np.maximum(1.0 - 2.0 * m / r, 1e-8)
        dm = (2.0 * np.pi * r * r / alpha_Y) * (B * dphi * dphi + (phi * phi) / (L * L))
        ddphi = (-(2.0 / r * (1.0 - m / r) - (4.0 * np.pi * r / (alpha_Y * L * L)) * phi * phi) * dphi + phi / (L * L)) / B
        return np.vstack([dphi, ddphi, dm])

    def bc(ya, yb):
        m_end = yb[2]
        B_end = max(1.0 - 2.0 * m_end / R_max, 1e-8)
        k_asympt = (1.0 / L) / np.sqrt(B_end) + (1.0 + m_end / L) / R_max
        return np.array([
            ya[0] - phi_R,
            ya[2] - M_surf_km,
            yb[1] + k_asympt * yb[0],
        ])

    sol = solve_bvp(fun, bc, r_grid, np.vstack([phi_init, dphi_init, m_init]), tol=1e-9, max_nodes=4000)
    if not sol.success:
        raise RuntimeError(f"Curved exterior BVP did not converge: {sol.message}")
    dphi_R_curved = float(sol.y[1, 0])
    k_curved = -dphi_R_curved / phi_R
    M_ADM_km = float(sol.y[2, -1])
    rx = sol.x
    phix, dphix, mx = sol.y[0], sol.y[1], sol.y[2]
    Bx = 1.0 - 2.0 * mx / rx
    pr_phix = (Bx * dphix**2 - (phix / L) ** 2) / (2.0 * alpha_Y)
    nup_x = (mx + 4.0 * np.pi * rx**3 * pr_phix) / (rx * (rx - 2.0 * mx))
    dnu_ext = float(np.trapezoid(nup_x, rx))
    nu_Rmax = 0.5 * np.log(1.0 - 2.0 * M_ADM_km / R_max)
    nu_R = nu_Rmax - dnu_ext
    return {
        "k_curved": k_curved,
        "k_flat": k_flat,
        "rel_diff_k": (k_curved - k_flat) / k_flat,
        "dphi_R_curved": dphi_R_curved,
        "dphi_R_flat": -k_flat * phi_R,
        "M_ADM_km": M_ADM_km,
        "Delta_M_phi_ext_Msun": (M_ADM_km - M_surf_km) / GM_KM,
        "nu_R": nu_R,
        "r_ext": rx,
        "phi_ext": phix,
        "dphi_ext": dphix,
        "m_ext": mx,
    }


def _solve_curved_scalar_with_k(rg, B, dB, nup, src, L, k_robin):
    """Finite-difference solver for B phi'' + (2B/r + B nu' + B'/2) phi' - phi/L^2 = -src
    on [0, R] with regular center phi'(0)=0 and Robin surface condition phi'(R) = -k_robin * phi(R)."""
    n = len(rg)
    h = rg[1] - rg[0]
    ab = np.zeros((3, n))
    rhs = -np.asarray(src, float).copy()
    ab[1, 0] = -6.0 / h**2 - 1.0 / L**2
    ab[0, 1] = 6.0 / h**2
    r = rg[1:-1]
    a = B[1:-1]
    b = 2.0 * B[1:-1] / r + B[1:-1] * nup[1:-1] + 0.5 * dB[1:-1]
    ab[2, 0:n - 2] = a / h**2 - b / (2.0 * h)
    ab[1, 1:-1] = -2.0 * a / h**2 - 1.0 / L**2
    ab[0, 2:n] = a / h**2 + b / (2.0 * h)
    R = rg[-1]
    k = -k_robin
    aN, bN = B[-1], 2.0 * B[-1] / R + B[-1] * nup[-1] + 0.5 * dB[-1]
    ab[2, n - 2] = 2.0 * aN / h**2
    ab[1, n - 1] = -2.0 * aN / h**2 + (aN / h**2 + bN / (2.0 * h)) * 2.0 * h * k - 1.0 / L**2
    return solve_banded((1, 1), ab, rhs)


def solve_yukawa_action_el_curved_star(
    eos, logrho, L=10.0, alpha_Y=1.0, ngrid=601, tol=2e-8, use_curved_robin=True, return_profile=False
):
    """
    Solve the coupled interior + curved exterior Einstein-scalar equilibrium for `yukawa_action_el`.
    When use_curved_robin=True (default), uses the exact curved-spacetime Dirichlet-to-Neumann
    surface slope k_curved = -phi_L'(R_*)/phi_L(R_*) from the exterior vacuum BVP.
    """
    cur = cc.integrate_r_tov(eos, logrho, lambda r: 1.0, None, ngrid=ngrid, p_weighted=True)
    phi = np.zeros(ngrid)
    x_prev = cur["rg"] / cur["rg"][-1]
    rho_p = np.zeros(ngrid)
    dphi = np.zeros(ngrid)
    prp_on = lambda rr: np.zeros_like(rr)
    k_robin = 1.0 / cur["R"] + 1.0 / L
    for it in range(100):
        rg, er, mg, hg = cur["rg"], cur["er"], cur["mg"], cur["hg"]
        phi0 = np.interp(rg / rg[-1], x_prev, phi)
        rho0 = np.interp(rg / rg[-1], x_prev, rho_p)
        f = np.exp(-phi0)
        P = eos.at(np.maximum(hg, 1e-12))[1]
        P[-1] = 0.0
        rs = np.maximum(rg, 1e-6)
        B = 1.0 - 2.0 * mg / rs
        B[0] = 1.0
        pr0 = prp_on(rg)
        nup = (mg + 4.0 * np.pi * rs**3 * (f * P + pr0)) / (rs * (rs - 2.0 * mg))
        nup[0] = 0.0
        dm = 4.0 * np.pi * rs**2 * (f * er + rho0)
        dB = 2.0 * mg / rs**2 - 2.0 * dm / rs
        dB[0] = 0.0
        if use_curved_robin and it in (1, 3, 6, 10) and abs(phi0[-1]) > 1e-14:
            ext_info = _solve_exterior_3comp_robin_iterate(rg[-1], mg[-1], float(phi0[-1]), L, alpha_Y=alpha_Y)
            k_robin = ext_info["k_curved"]
        elif not use_curved_robin:
            k_robin = 1.0 / rg[-1] + 1.0 / L
        phin = _solve_curved_scalar_with_k(rg, B, dB, nup, alpha_Y * f * er, L, k_robin)
        phi_new = phin if it == 0 else 0.5 * phin + 0.5 * phi0
        dphi = np.gradient(phi_new, rg, edge_order=2)
        dphi[0] = 0.0
        dphi[-1] = -k_robin * phi_new[-1]
        rho_p = (B * dphi**2 + phi_new**2 / L**2) / (2.0 * alpha_Y)
        pr_p = (B * dphi**2 - phi_new**2 / L**2) / (2.0 * alpha_Y)
        x_prev, phi = rg / rg[-1], phi_new
        Rc = rg[-1]
        mk = lambda arr: (lambda I, Rr: (lambda r: float(I(min(r, Rr)))))(PchipInterpolator(rg, arr), Rc)
        prp_on = (lambda xp, arr: (lambda rr: np.interp(rr / rr[-1], xp, arr)))(x_prev, pr_p.copy())
        new = cc.integrate_r_tov(
            eos, logrho, mk(np.exp(-phi)), mk(-phi), mk(rho_p), mk(pr_p), ngrid=ngrid, p_weighted=True
        )
        res = max(abs(new["R"] / cur["R"] - 1.0), abs(new["M"] / cur["M"] - 1.0), abs(phi_new[0] - phi0[0]))
        cur = new
        if res < tol and it > 3:
            break
    else:
        raise RuntimeError(f"Curved Yukawa action-EL fixed point did not converge (logrho={logrho}, L={L})")

    R = float(cur["R"])
    M_surf_Msun = float(cur["M"])
    M_surf_km = M_surf_Msun * GM_KM
    phiR = float(phi[-1])
    ext_info = solve_exterior_dtN_curved(R, M_surf_km, phiR, L, alpha_Y=alpha_Y, tol=1e-12)

    rg, mg, hg = cur["rg"], cur["mg"], cur["hg"]
    er, P, rho0_cgs, cs2 = eos.at(np.maximum(hg, 1e-12))
    er[-1] = 0.0
    P[-1] = 0.0
    rho0_cgs[-1] = 0.0
    rho0_geom = rho0_cgs * CONV
    rs = np.maximum(rg, 1e-6)
    B_arr = np.maximum(1.0 - 2.0 * mg / rs, 1e-8)
    B_arr[0] = 1.0
    M_baryon_Msun = float(np.trapezoid(4.0 * np.pi * rg**2 * rho0_geom / np.sqrt(B_arr), rg)) / GM_KM

    # Reference flat-Robin exterior integral (on frozen Schwarzschild metric) for comparison
    re = np.linspace(R, R + 40.0 * L, 8001)
    pe = phiR * (R / re) * np.exp(-(re - R) / L)
    dpe = -pe * (1.0 / re + 1.0 / L)
    rho_e_flat = ((1.0 - 2.0 * M_surf_km / re) * dpe**2 + (pe / L) ** 2) / (2.0 * alpha_Y)
    dM_ext_flat_Msun = float(np.trapezoid(4.0 * np.pi * re**2 * rho_e_flat, re)) / GM_KM

    nu_arr = ext_info["nu_R"] - (hg - phi - (hg[-1] - phi[-1]))
    out = {
        "R_km": R,
        "M_Msun": M_surf_Msun,
        "M_surface_Msun": M_surf_Msun,
        "Delta_M_phi_ext_Msun": ext_info["Delta_M_phi_ext_Msun"],
        "Delta_M_phi_ext_curved_Msun": ext_info["Delta_M_phi_ext_Msun"],
        "Delta_M_phi_ext_flat_approx_Msun": dM_ext_flat_Msun,
        "M_total_Msun": ext_info["M_ADM_km"] / GM_KM,
        "M_ADM_Msun": ext_info["M_ADM_km"] / GM_KM,
        "M_baryon_Msun": M_baryon_Msun,
        "f_center": float(np.exp(-phi[0])),
        "f_surface": float(np.exp(-phiR)),
        "phi_center": float(phi[0]),
        "phi_surface": phiR,
        "dphi_dr_R_curved": ext_info["dphi_R_curved"],
        "dphi_dr_R_flat_approx": ext_info["dphi_R_flat"],
        "k_curved": ext_info["k_curved"],
        "k_flat": ext_info["k_flat"],
        "rel_diff_dphi_R": ext_info["rel_diff_k"],
        "nu_R": ext_info["nu_R"],
        "fixed_point_iterations": it + 1,
        "fixed_point_last_rel_residual": float(res),
    }
    if return_profile:
        out["profile"] = {
            "rg": rg,
            "mg": mg,
            "hg": hg,
            "er": er,
            "P": P,
            "rho0": rho0_geom,
            "cs2": cs2,
            "phi": phi,
            "dphi": dphi,
            "rho_phi": rho_p,
            "pr_phi": pr_p,
            "B": B_arr,
            "nu": nu_arr,
            "ext": ext_info,
            "L": L,
            "alpha_Y": alpha_Y,
        }
    return out


def compute_hessian_sturm_liouville_diagnostics(eos, log10_rho0c, L_km=10.0, alpha_Y=1.0):
    """
    Evaluate the concrete analytical coefficients of the second variation delta^2 S
    (coupled (xi, delta phi_L) Sturm-Liouville / Hessian operator after eliminating
    the metric constraints delta lambda and delta nu' at fixed baryon number delta N_B = 0)
    on the curved `yukawa_action_el` equilibrium profile.
    """
    st_eq = solve_yukawa_action_el_curved_star(eos, log10_rho0c, L=L_km, alpha_Y=alpha_Y, return_profile=True)
    prof = st_eq["profile"]
    rg = prof["rg"][1:-1]
    mg = prof["mg"][1:-1]
    er = prof["er"][1:-1]
    P = prof["P"][1:-1]
    cs2 = prof["cs2"][1:-1]
    phi = prof["phi"][1:-1]
    dphi = prof["dphi"][1:-1]
    B = prof["B"][1:-1]
    nu = prof["nu"][1:-1]
    f = np.exp(-phi)
    Gamma1 = (er + P) * cs2 / np.maximum(P, 1e-30)

    lam = -0.5 * np.log(B)
    elam = 1.0 / np.sqrt(B)
    e2lam = 1.0 / B
    enu = np.exp(nu)
    pr_phi = (B * dphi**2 - (phi / L_km) ** 2) / (2.0 * alpha_Y)
    rho_phi = (B * dphi**2 + (phi / L_km) ** 2) / (2.0 * alpha_Y)
    nup = (mg + 4.0 * np.pi * rg**3 * (f * P + pr_phi)) / (rg * (rg - 2.0 * mg))
    lamp = (4.0 * np.pi * rg**2 * (f * er + rho_phi) - mg / rg) / (rg - 2.0 * mg)
    dP_dr = -(er + P) * (nup - dphi)

    # Kinetic weights W_xi(r) and W_phi(r)
    W_xi = 4.0 * np.pi * (rg**2) * f * (er + P) * (elam**3) * np.exp(-nu)
    W_phi = 4.0 * np.pi * (rg**2) * elam * np.exp(-nu) / alpha_Y

    # Stiffness coefficients P_xi(r) and P_phi(r)
    P_xi = 4.0 * np.pi * (rg**2) * f * Gamma1 * P * elam * enu
    P_phi = 4.0 * np.pi * (rg**2) * enu / (alpha_Y * elam)

    # Potential and constraint-induced coefficients Q_{xi xi}(r), Q_{phi phi}(r), C_{xi phi}(r)
    # after eliminating delta lambda and delta nu' from the linearized Einstein equations:
    q_chandrasekhar = (
        -4.0 * dP_dr / rg
        - (dP_dr**2) / (er + P)
        + 8.0 * np.pi * e2lam * (f * (er + P) + B * (dphi**2) / alpha_Y) * (f * P + pr_phi + er + P)
    )
    Q_xixi = 4.0 * np.pi * (rg**2) * f * elam * enu * (
        q_chandrasekhar + (er + P) * (dphi**2) * (1.0 - 4.0 * np.pi * (rg**2) * f * (er + P) / B)
    )
    Q_phiphi = 4.0 * np.pi * (rg**2) * elam * enu * (
        (1.0 / (alpha_Y * (L_km**2)) + f * er)
        + (2.0 * (dphi**2) / alpha_Y) * (1.0 - 2.0 * np.pi * (rg**2) * (phi**2) / (alpha_Y * (L_km**2)))
    )
    C_xiphi_div = 4.0 * np.pi * (rg**2) * f * (er + P) * elam * enu * (
        1.0 - 2.0 * rg * dphi / alpha_Y
    )
    C_xiphi_pot = 4.0 * np.pi * (rg**2) * f * elam * enu * (
        (er + P) * (nup - dphi) - 8.0 * np.pi * rg * e2lam * (f * (er + P)) * (f * P + pr_phi)
    )

    return {
        "log10_rho0c": float(log10_rho0c),
        "L_km": float(L_km),
        "alpha_Y": float(alpha_Y),
        "min_W_xi": float(np.min(W_xi)),
        "min_W_phi": float(np.min(W_phi)),
        "min_P_xi": float(np.min(P_xi)),
        "min_P_phi": float(np.min(P_phi)),
        "median_Q_xixi": float(np.median(Q_xixi)),
        "median_Q_phiphi": float(np.median(Q_phiphi)),
        "max_abs_C_xiphi_div": float(np.max(np.abs(C_xiphi_div))),
        "max_abs_C_xiphi_pot": float(np.max(np.abs(C_xiphi_pot))),
        "surface_eps_jump": 0.0,
        "surface_dphi_jump_coeff": 0.0,
    }


def solve_yukawa_coupled_omega2_legacy_fd(
    eos, log10_rho0c, L_km=10.0, alpha_Y=1.0, N=450, include_dphi=True, include_metric_in_scalar=True
):
    """
    Solve the coupled (delta h, delta phi_L, delta m) radial pulsation eigenvalue problem
    directly on the curved `yukawa_action_el` equilibrium sequence (including the curved
    exterior Einstein-scalar Robin boundary condition k_curved and exterior scalar tail mass
    perturbation delta(Delta M_{phi,ext})).
    """
    st_eq = solve_yukawa_action_el_curved_star(eos, log10_rho0c, L=L_km, alpha_Y=alpha_Y, return_profile=True)
    prof = st_eq["profile"]
    rg = prof["rg"]
    R = st_eq["R_km"]
    hc = float(prof["hg"][0])
    h_interp = PchipInterpolator(rg, prof["hg"])
    m_interp = PchipInterpolator(rg, prof["mg"])
    phi_interp = PchipInterpolator(rg, prof["phi"])
    dphi_interp = PchipInterpolator(rg, prof["dphi"])
    nu_interp = PchipInterpolator(rg, prof["nu"])

    r_breaks = [0.0]
    for h_tr in sorted(eos.hb[1:], reverse=True):
        if 1e-9 < h_tr < hc:
            r_tr = brentq(lambda x: float(h_interp(x)) - h_tr, 0.0, R)
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

    h = h_interp(r_mid)
    e, P, _, cs2 = eos.at(h)
    phi = phi_interp(r_mid)
    dphi = dphi_interp(r_mid)
    f = np.exp(-phi)
    m = m_interp(r_mid)
    nu = nu_interp(r_mid)
    B = np.maximum(1.0 - 2.0 * m / r_mid, 1e-8)
    rho_phi = (B * dphi**2 + (phi / L_km) ** 2) / (2.0 * alpha_Y)
    pr_phi = (B * dphi**2 - (phi / L_km) ** 2) / (2.0 * alpha_Y)
    A = (m + 4.0 * np.pi * (r_mid**3) * (f * P + pr_phi)) / (r_mid * (r_mid - 2.0 * m))
    dm_bg = 4.0 * np.pi * (r_mid**2) * (f * e + rho_phi)
    dB_bg = 2.0 * m / (r_mid**2) - 2.0 * dm_bg / r_mid
    ddphi = (-(2.0 * B / r_mid + B * A + 0.5 * dB_bg) * dphi + phi / (L_km**2) - alpha_Y * f * e) / B

    M_int = np.tril(np.ones((N_tot, N_tot)), -1) * dr[None, :] + 0.5 * np.eye(N_tot) * dr[None, :]
    deps_dh = (e + P) / cs2

    k_rob = st_eq["k_curved"]
    D_r = np.zeros((N_tot, N_tot))
    for i in range(1, N_tot - 1):
        D_r[i, i + 1] = 1.0 / (r_mid[i + 1] - r_mid[i - 1])
        D_r[i, i - 1] = -1.0 / (r_mid[i + 1] - r_mid[i - 1])
    D_r[-1, -1] = -k_rob

    S_dh = np.diag(4.0 * np.pi * (r_mid**2) * f * deps_dh)

    if include_dphi:
        L_op = np.zeros((N_tot, N_tot))
        for i in range(N_tot):
            ri = r_mid[i]
            Bi = B[i]
            bi = 2.0 * Bi / ri + Bi * A[i] + 0.5 * dB_bg[i]
            ci = -(1.0 / (L_km**2) + alpha_Y * f[i] * e[i])
            if i == 0:
                h_step = r_mid[1] - r_mid[0]
                L_op[0, 0] = -2.0 * Bi / (h_step**2) + ci
                L_op[0, 1] = 2.0 * Bi / (h_step**2)
            elif i < N_tot - 1:
                hm = r_mid[i] - r_mid[i - 1]
                hp = r_mid[i + 1] - r_mid[i]
                denom = 0.5 * (hm + hp)
                L_op[i, i - 1] = Bi / (hm * denom) - bi / (hm + hp)
                L_op[i, i] = -Bi * (1.0 / hm + 1.0 / hp) / denom + ci
                L_op[i, i + 1] = Bi / (hp * denom) + bi / (hm + hp)
            else:
                hm = r_mid[i] - r_mid[i - 1]
                L_op[i, i - 1] = 2.0 * Bi / (hm**2)
                L_op[i, i] = -2.0 * Bi / (hm**2) - (2.0 * Bi / hm + bi) * k_rob + ci

        if include_metric_in_scalar:
            diag_extra_dphi = -4.0 * np.pi * r_mid * (f * (P - e) + 2.0 * phi / (alpha_Y * L_km**2)) * dphi
            L_op += np.diag(diag_extra_dphi)
            coef_dm_in_scalar = -(2.0 * ddphi / r_mid + 2.0 * dphi / (r_mid**2))
            coef_dh_in_scalar = alpha_Y * f * deps_dh + 4.0 * np.pi * r_mid * f * (e + P) * (1.0 - 1.0 / cs2) * dphi
            S_dphi_unit = np.diag(4.0 * np.pi * (r_mid**2) * (-f * e + phi / (alpha_Y * L_km**2))) + (
                4.0 * np.pi * (r_mid**2) * (B * dphi / alpha_Y)
            )[:, None] * D_r
            DM_dh = M_int @ S_dh
            DM_dphi = M_int @ S_dphi_unit
            L_op_full = L_op + (coef_dm_in_scalar[:, None] * DM_dphi)
            rhs_full = -np.diag(coef_dh_in_scalar) - (coef_dm_in_scalar[:, None] * DM_dh)
            DPHI_mat = np.linalg.solve(L_op_full, rhs_full)
        else:
            rhs_mat = np.diag(-alpha_Y * f * deps_dh)
            DPHI_mat = np.linalg.solve(L_op, rhs_mat)
        DDPHI_mat = D_r @ DPHI_mat
    else:
        DPHI_mat = np.zeros((N_tot, N_tot))
        DDPHI_mat = np.zeros((N_tot, N_tot))

    S_explicit = (
        S_dh
        + (4.0 * np.pi * (r_mid**2) * (-f * e + phi / (alpha_Y * L_km**2)))[:, None] * DPHI_mat
        + (4.0 * np.pi * (r_mid**2) * (B * dphi / alpha_Y))[:, None] * DDPHI_mat
    )
    M_dm_lhs = np.eye(N_tot) + M_int @ np.diag(4.0 * np.pi * r_mid * (dphi**2) / alpha_Y)
    DM_mat = np.linalg.solve(M_dm_lhs, M_int @ S_explicit)

    c_dm = (1.0 + 2.0 * r_mid * A) / (r_mid * (r_mid - 2.0 * m))
    c_dh = (4.0 * np.pi * (r_mid**3) * f * (e + P)) / (r_mid * (r_mid - 2.0 * m))
    c_dphi_P = (
        (4.0 * np.pi * (r_mid**3) * (-f * P - phi / (alpha_Y * L_km**2)))[:, None] * DPHI_mat
        + (4.0 * np.pi * (r_mid**3) * (B * dphi / alpha_Y))[:, None] * DDPHI_mat
        - (4.0 * np.pi * (r_mid**2) * (dphi**2 / alpha_Y))[:, None] * DM_mat
    )
    DA_mat = (c_dm[:, None] * DM_mat) + np.diag(c_dh) + (c_dphi_P / (r_mid * (r_mid - 2.0 * m))[:, None]) - DDPHI_mat

    e2lam = 1.0 / B
    e_minus_2nu = np.exp(-2.0 * nu)
    w_inert = (e2lam * e_minus_2nu) / (4.0 * np.pi * (r_mid**2) * f * (e + P))
    W_mat = w_inert[:, None] * DM_mat

    def F_shoot(om2):
        dh_s = np.linalg.solve(np.eye(N_tot) + M_int @ (DA_mat + om2 * W_mat), np.ones(N_tot))
        dm_surf = float(DM_mat[-1, :] @ dh_s)
        if include_dphi and abs(st_eq["phi_surface"]) > 1e-12:
            dphi_surf = float(DPHI_mat[-1, :] @ dh_s)
            dm_ext = 2.0 * (st_eq["Delta_M_phi_ext_curved_Msun"] * GM_KM) * (dphi_surf / st_eq["phi_surface"])
        else:
            dm_ext = 0.0
        return dm_surf + dm_ext

    F0 = F_shoot(0.0)
    d_om = 1e-6
    dF = (F_shoot(d_om) - F_shoot(-d_om)) / (2.0 * d_om)
    om2_tan = -F0 / dF
    best_om2 = brentq(F_shoot, om2_tan - 0.0006, om2_tan + 0.0006, xtol=1e-11)
    return {
        "omega2_s2": best_om2 * (C_KMS**2),
        "omega2_tangent_s2": om2_tan * (C_KMS**2),
        "F0": F0,
        "M_ADM_Msun": st_eq["M_ADM_Msun"],
        "M_surface_Msun": st_eq["M_surface_Msun"],
        "M_baryon_Msun": st_eq["M_baryon_Msun"],
        "R_km": st_eq["R_km"],
    }


def linearized_exterior_response(ext_info, R, L, alpha_Y=1.0, om2_km2=0.0, n_pts=200, tol=1e-11):
    if ext_info.get("sol") is None or L <= 0.0:
        k_flat = 1.0 / R + (1.0 / L if L > 0 else 0.0)
        return {"K_phi": k_flat, "K_m": 0.0, "M_phi": 0.0, "M_m": 1.0}
    bg_sol = ext_info["sol"]
    R_max = ext_info["R_max"]
    r_grid = R + (R_max - R) * (np.linspace(0.0, 1.0, n_pts) ** 1.8)

    def lin_rhs_pair(r, y):
        y_bg = bg_sol.sol(r)
        phi0, psi0, m0, nu0 = y_bg[0], y_bg[1], y_bg[2], y_bg[3]
        B0 = np.maximum(1.0 - 2.0 * m0 / r, 1e-8)
        C0 = (2.0 / r) * (1.0 - m0 / r) - (4.0 * np.pi * r / (alpha_Y * L * L)) * (phi0 ** 2)
        dpsi0 = (-C0 * psi0 + phi0 / (L * L)) / B0
        em2nu0 = np.exp(-2.0 * nu0)

        a_psi = -C0 / B0
        a_phi = (1.0 / (L * L) - om2_km2 * em2nu0 + (8.0 * np.pi * r / (alpha_Y * L * L)) * phi0 * psi0) / B0
        a_m = (2.0 * psi0 / (r * r) - 2.0 * dpsi0 / r) / B0

        b_psi = (4.0 * np.pi * r * r / alpha_Y) * B0 * psi0
        b_phi = (4.0 * np.pi * r * r / (alpha_Y * L * L)) * phi0
        b_m = -(4.0 * np.pi * r / alpha_Y) * (psi0 ** 2)

        dphi1, dpsi1, dm1, dphi2, dpsi2, dm2 = y
        return np.vstack([
            dpsi1,
            a_psi * dpsi1 + a_phi * dphi1 + a_m * dm1,
            b_psi * dpsi1 + b_phi * dphi1 + b_m * dm1,
            dpsi2,
            a_psi * dpsi2 + a_phi * dphi2 + a_m * dm2,
            b_psi * dpsi2 + b_phi * dphi2 + b_m * dm2,
        ])

    y_end_bg = bg_sol.sol(R_max)
    phi_end, _, m_end, nu_end = [float(v) for v in y_end_bg]
    B_end = max(1.0 - 2.0 * m_end / R_max, 1e-8)
    q_eff = np.sqrt(max(1.0 / (L * L) - om2_km2 * np.exp(-2.0 * nu_end), 1e-16))
    k_asympt = q_eff / np.sqrt(B_end) + (1.0 + m_end * q_eff) / R_max
    dk_dm = q_eff * (B_end ** (-1.5)) / R_max + q_eff / R_max

    def lin_bc(ya, yb):
        return np.array([
            ya[0] - 1.0,
            ya[2] - 0.0,
            yb[1] + k_asympt * yb[0] + dk_dm * phi_end * yb[2],
            ya[3] - 0.0,
            ya[5] - 1.0,
            yb[4] + k_asympt * yb[3] + dk_dm * phi_end * yb[5],
        ])

    y_init = np.zeros((6, len(r_grid)))
    y_init[0] = (R / r_grid) * np.exp(-(r_grid - R) / L)
    y_init[1] = -y_init[0] * (1.0 / r_grid + 1.0 / L)
    y_init[5] = 1.0
    sol_lin = solve_bvp(lin_rhs_pair, lin_bc, r_grid, y_init, tol=tol, max_nodes=10000)
    if not sol_lin.success:
        raise RuntimeError(f"Linearized exterior BVP did not converge: {sol_lin.message}")
    return {
        "K_phi": float(-sol_lin.y[1, 0]),
        "K_m": float(sol_lin.y[4, 0]),
        "M_phi": float(sol_lin.y[2, -1]),
        "M_m": float(sol_lin.y[5, -1]),
    }

def build_gl2_operators(r_breaks, N_target=480, base_cells=(7, 5, 3, 2, 1, 2)):
    n_seg = len(r_breaks) - 1
    weights = list(base_cells[:n_seg])
    while len(weights) < n_seg:
        weights.append(2)
    w_sum = sum(weights)
    mult = max(1, int(round((N_target / 2.0) / w_sum)))
    cells_per_seg = [max(1, mult * w) for w in weights]
    edges = [0.0]
    for k in range(n_seg):
        ra, rb = r_breaks[k], r_breaks[k + 1]
        nc = cells_per_seg[k]
        edges.extend(np.linspace(ra, rb, nc + 1)[1:])
    edges = np.array(edges)
    n_cells = len(edges) - 1
    N = 2 * n_cells
    sq3 = np.sqrt(3.0)
    c1, c2 = 0.5 - sq3 / 6.0, 0.5 + sq3 / 6.0
    a11, a12 = 0.25, 0.25 - sq3 / 6.0
    a21, a22 = 0.25 + sq3 / 6.0, 0.25

    r_nodes = np.zeros(N)
    Q = np.zeros((N, N))
    for m in range(n_cells):
        ra, rb = edges[m], edges[m + 1]
        hm = rb - ra
        i1, i2 = 2 * m, 2 * m + 1
        r_nodes[i1] = ra + c1 * hm
        r_nodes[i2] = ra + c2 * hm
        for j in range(m):
            hj = edges[j + 1] - edges[j]
            Q[i1, 2 * j] = 0.5 * hj
            Q[i1, 2 * j + 1] = 0.5 * hj
            Q[i2, 2 * j] = 0.5 * hj
            Q[i2, 2 * j + 1] = 0.5 * hj
        Q[i1, i1], Q[i1, i2] = a11 * hm, a12 * hm
        Q[i2, i1], Q[i2, i2] = a21 * hm, a22 * hm
    w_surf = np.zeros(N)
    for j in range(n_cells):
        hj = edges[j + 1] - edges[j]
        w_surf[2 * j] = 0.5 * hj
        w_surf[2 * j + 1] = 0.5 * hj
    Q2 = Q @ Q
    w2_surf = w_surf @ Q
    return r_nodes, Q, Q2, w_surf, w2_surf, N

def solve_yukawa_coupled_omega2_gl2(
    eos,
    log10_rho0c,
    L_km=10.0,
    alpha_Y=1.0,
    N=480,
    include_dphi=True,
    include_metric_in_scalar=True,
    dynamical_scalar_wave=True,
    st_eq=None,
):
    if st_eq is None:
        st_eq = solve_yukawa_action_el_curved_star(eos, log10_rho0c, L=L_km, alpha_Y=alpha_Y, return_profile=True)
    prof = st_eq["profile"]
    rg = prof["rg"]
    R = st_eq["R_km"]
    M_surf_km = st_eq["M_surface_Msun"] * GM_KM
    phi_R = st_eq["phi_surface"]
    hc = float(prof["hg"][0])

    # Upgrade ext_info to 4-comp BVP if "sol" not present
    ext_info = prof["ext"]
    if ext_info.get("sol") is None:
        ext_info = solve_exterior_4comp_bvp(R, M_surf_km, phi_R, L_km, alpha_Y=alpha_Y)
        prof["ext"] = ext_info

    h_interp = PchipInterpolator(rg, prof["hg"])
    m_interp = PchipInterpolator(rg, prof["mg"])
    phi_interp = PchipInterpolator(rg, prof["phi"])
    dphi_interp = PchipInterpolator(rg, prof["dphi"])

    r_breaks = [0.0]
    for h_tr in sorted(eos.hb[1:], reverse=True):
        if 1e-9 < h_tr < hc:
            r_tr = brentq(lambda x: float(h_interp(x)) - h_tr, 0.0, R)
            r_breaks.append(r_tr)
    r_breaks.append(R)
    r_breaks = sorted(r_breaks)

    r_nodes, Q, Q2, w_surf, w2_surf, N_gl = build_gl2_operators(r_breaks, N_target=N)
    h_n = np.maximum(h_interp(r_nodes), 1e-12)
    e_n, P_n, _, cs2_n = eos.at(h_n)
    phi_n = phi_interp(r_nodes)
    dphi_n = dphi_interp(r_nodes)
    f_n = np.exp(-phi_n)
    m_n = m_interp(r_nodes)
    B_n = np.maximum(1.0 - 2.0 * m_n / r_nodes, 1e-8)
    prp_n = (B_n * (dphi_n ** 2) - (phi_n ** 2) / (L_km * L_km)) / (2.0 * alpha_Y) if L_km > 0 else np.zeros(N_gl)
    A_n = (m_n + 4.0 * np.pi * (r_nodes ** 3) * (f_n * P_n + prp_n)) / (r_nodes * (r_nodes - 2.0 * m_n))
    nu_n = ext_info["nu_R"] - (h_n - phi_n - (prof["hg"][-1] - phi_R))
    em2nu_n = np.exp(-2.0 * nu_n)
    deps_dh = (e_n + P_n) / cs2_n
    C0_n = (2.0 / r_nodes) * (1.0 - m_n / r_nodes) + 4.0 * np.pi * r_nodes * (f_n * (P_n - e_n) - (phi_n ** 2) / (alpha_Y * L_km * L_km))
    ddphi_n = (-C0_n * dphi_n + phi_n / (L_km * L_km) - alpha_Y * f_n * e_n) / B_n

    c_m_dm = -4.0 * np.pi * r_nodes * (dphi_n ** 2) / alpha_Y if include_dphi else np.zeros(N_gl)
    c_m_dh = 4.0 * np.pi * (r_nodes ** 2) * f_n * deps_dh
    c_m_dphi = 4.0 * np.pi * (r_nodes ** 2) * (-f_n * e_n + phi_n / (alpha_Y * L_km * L_km)) if include_dphi else np.zeros(N_gl)
    c_m_ddphi = 4.0 * np.pi * (r_nodes ** 2) * B_n * dphi_n / alpha_Y if include_dphi else np.zeros(N_gl)

    M_lhs = np.eye(N_gl) - Q @ np.diag(c_m_dm)
    M_inv_Q = np.linalg.solve(M_lhs, Q)
    DM_dh = M_inv_Q @ np.diag(c_m_dh)
    DM_du = M_inv_Q @ (np.diag(c_m_dphi) @ Q2 + np.diag(c_m_ddphi) @ Q)
    dm_dc0 = M_inv_Q @ c_m_dphi

    v_s_dm = w_surf @ np.diag(c_m_dm)
    v_s_dh = w_surf @ np.diag(c_m_dh)
    v_s_du = w_surf @ (np.diag(c_m_dphi) @ Q2 + np.diag(c_m_ddphi) @ Q)
    s_dc0 = float(w_surf @ c_m_dphi)

    dmR_dh = v_s_dh + v_s_dm @ DM_dh
    dmR_du = v_s_du + v_s_dm @ DM_du
    dmR_dc0 = s_dc0 + float(v_s_dm @ dm_dc0)

    a1_sc = C0_n / B_n
    a0_sc_base = -(
        1.0 / (L_km * L_km)
        + alpha_Y * f_n * e_n
        + (4.0 * np.pi * r_nodes * (f_n * (P_n - e_n) + 2.0 * phi_n / (alpha_Y * L_km * L_km)) * dphi_n if include_metric_in_scalar else 0.0)
    ) / B_n
    a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n if include_metric_in_scalar else np.zeros(N_gl)
    a_h_sc = (
        alpha_Y * f_n * deps_dh
        + (4.0 * np.pi * r_nodes * f_n * (e_n + P_n) * (1.0 - 1.0 / cs2_n) * dphi_n if include_metric_in_scalar else 0.0)
    ) / B_n

    c_dm_h = -(1.0 + 2.0 * r_nodes * A_n) / (r_nodes * (r_nodes - 2.0 * m_n))
    c_dh_h = -4.0 * np.pi * (r_nodes ** 2) * f_n * (e_n + P_n) / (r_nodes - 2.0 * m_n)
    c_dphi_h = 4.0 * np.pi * (r_nodes ** 2) * (f_n * P_n + phi_n / (alpha_Y * L_km * L_km)) / (r_nodes - 2.0 * m_n) if include_dphi else np.zeros(N_gl)
    c_ddphi_h = (1.0 - 4.0 * np.pi * (r_nodes ** 2) * B_n * dphi_n / (alpha_Y * (r_nodes - 2.0 * m_n))) if include_dphi else np.zeros(N_gl)
    c_dm_phi_h = 4.0 * np.pi * r_nodes * (dphi_n ** 2) / (alpha_Y * (r_nodes - 2.0 * m_n)) if include_dphi else np.zeros(N_gl)

    w_inert_h = em2nu_n / (B_n * 4.0 * np.pi * (r_nodes ** 2) * f_n * (e_n + P_n))
    P_of_h = Q @ np.diag(c_m_dh)

    ext_cache = {}
    def get_ext(om2_val):
        om2_q = om2_val if dynamical_scalar_wave else 0.0
        key = round(float(om2_q), 12)
        if key not in ext_cache:
            ext_cache[key] = solve_linearized_exterior_response(ext_info, R, L_km, alpha_Y=alpha_Y, om2_km2=om2_q)
        return ext_cache[key]

    def F_shoot(om2):
        ext_r = get_ext(om2)
        if include_dphi and L_km > 0:
            K_phi = ext_r["K_phi"]
            K_m = ext_r["K_m"]
            M_phi = ext_r["M_phi"]
            M_m = ext_r["M_m"]
            denom_bc = K_phi - K_m * dmR_dc0
            v_c0_u = -(w_surf + K_phi * w2_surf - K_m * dmR_du) / denom_bc
            v_c0_h = (K_m * dmR_dh) / denom_bc

            Phi_u = np.outer(np.ones(N_gl), v_c0_u) + Q2
            Phi_h = np.outer(np.ones(N_gl), v_c0_h)
            DM_tot_u = DM_du + np.outer(dm_dc0, v_c0_u)
            DM_tot_h = DM_dh + np.outer(dm_dc0, v_c0_h)

            om2_sc = om2 if dynamical_scalar_wave else 0.0
            a0_sc = a0_sc_base + (om2_sc * em2nu_n) / B_n
            L_u = np.eye(N_gl) + a1_sc[:, None] * Q + a0_sc[:, None] * Phi_u + a_m_sc[:, None] * DM_tot_u
            R_h = -(a0_sc[:, None] * Phi_h + a_m_sc[:, None] * DM_tot_h + np.diag(a_h_sc))
            U_of_h = np.linalg.solve(L_u, R_h)
            c0_of_h = v_c0_h + v_c0_u @ U_of_h
            DPHI_mat = np.outer(np.ones(N_gl), c0_of_h) + Q2 @ U_of_h
            DDPHI_mat = Q @ U_of_h
            DM_mat = DM_dh + DM_du @ U_of_h + np.outer(dm_dc0, c0_of_h)
        else:
            M_phi = 0.0
            M_m = 1.0
            c0_of_h = np.zeros(N_gl)
            U_of_h = np.zeros((N_gl, N_gl))
            DPHI_mat = np.zeros((N_gl, N_gl))
            DDPHI_mat = np.zeros((N_gl, N_gl))
            DM_mat = DM_dh

        R_h_op = (
            (c_dm_h + c_dm_phi_h)[:, None] * DM_mat
            + np.diag(c_dh_h)
            + c_dphi_h[:, None] * DPHI_mat
            + c_ddphi_h[:, None] * DDPHI_mat
            - om2 * (w_inert_h[:, None] * P_of_h)
        )
        dh_sol = np.linalg.solve(np.eye(N_gl) - Q @ R_h_op, np.ones(N_gl))
        u_sol = U_of_h @ dh_sol
        c0_sol = float(c0_of_h @ dh_sol)
        dphi_R = c0_sol + float(w2_surf @ u_sol) if include_dphi else 0.0
        dm_R = float(dmR_dh @ dh_sol + dmR_du @ u_sol + dmR_dc0 * c0_sol)
        return M_m * dm_R + M_phi * dphi_R

    F0 = F_shoot(0.0)
    d_om = 5e-5
    dF = (F_shoot(d_om) - F_shoot(-d_om)) / (2.0 * d_om)
    om2_tan = -F0 / dF
    x0, f0 = 0.0, F0
    x1, f1 = om2_tan, F_shoot(om2_tan)
    for _ in range(4):
        if abs(f1 - f0) < 1e-18:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0 = x1, f1
        x1, f1 = x2, F_shoot(x2)
    best_om2 = x1
    dM_ADM_res = F_shoot(best_om2)
    return {
        "omega2_s2": best_om2 * (C_KMS**2),
        "delta_M_ADM_over_M_ADM": float(dM_ADM_res / (st_eq["M_ADM_Msun"] * GM_KM)),
        "omega2_tangent_s2": om2_tan * (C_KMS**2),
        "F0": F0,
        "M_ADM_Msun": st_eq["M_ADM_Msun"],
        "M_surface_Msun": st_eq["M_surface_Msun"],
        "M_baryon_Msun": st_eq["M_baryon_Msun"],
        "R_km": st_eq["R_km"],
        "N_gl": N_gl,
    }

solve_linearized_exterior_response = linearized_exterior_response


solve_yukawa_coupled_omega2 = solve_yukawa_coupled_omega2_gl2


def _quad_root(xs, ys):
    c = np.polyfit(xs, ys, 2)
    disc = c[1] ** 2 - 4.0 * c[0] * c[2]
    roots = np.roots(c) if disc >= 0 else np.array([])
    roots = [float(np.real(r)) for r in roots if min(xs) - 0.02 <= np.real(r) <= max(xs) + 0.02]
    if not roots:
        return float("nan")
    mid = float(np.mean(xs))
    return min(roots, key=lambda r: abs(r - mid))


def _seq_bounds(eos_name):
    return (15.10, 15.30) if eos_name == "H4" else ((15.20, 15.42) if eos_name == "SLy" else (15.28, 15.50))


SEQUENCE_CASES = [
    ("SLy", 10.0, 1.0), ("SLy", 10.0, 1.0 / 3.0), ("SLy", 1.0, 1.0), ("SLy", 1.0, 1.0 / 3.0),
    ("APR4", 10.0, 1.0), ("APR4", 1.0, 1.0), ("H4", 10.0, 1.0), ("H4", 1.0, 1.0),
]
STABILITY_CASES = [("SLy", 10.0, 1.0), ("SLy", 10.0, 1.0 / 3.0), ("SLy", 1.0, 1.0), ("APR4", 10.0, 1.0), ("H4", 10.0, 1.0)]
N_CONVERGENCE = [80, 120, 160, 240, 360, 480, 600, 800, 1200]


def run_all_curved_yukawa_audits():
    import time
    t_start = time.time()
    print("1. 4-component exterior BVP (phi_L, psi_L, m, nu): R_max independence...", flush=True)
    eos_sly = tp.EOS("SLy")
    st_ref = solve_yukawa_action_el_curved_star(eos_sly, 15.0, L=10.0, alpha_Y=1.0)
    rmax_checks = []
    for N_m in [10.0, 12.0, 14.0, 16.0]:
        ext_n = solve_exterior_4comp_bvp(st_ref["R_km"], st_ref["M_surface_Msun"] * GM_KM, st_ref["phi_surface"], 10.0, alpha_Y=1.0, N_mult=N_m)
        rmax_checks.append({
            "N_mult": N_m, "R_max_km": ext_n["R_max"], "k_curved_km_inv": ext_n["k_curved"], "k_flat_km_inv": ext_n["k_flat"],
            "rel_diff_k": ext_n["rel_diff_k"], "dphi_dr_R_curved": ext_n["dphi_R_curved"], "dphi_dr_R_flat": ext_n["dphi_R_flat"],
            "Delta_M_phi_ext_Msun": ext_n["Delta_M_phi_ext_Msun"], "M_ADM_Msun": ext_n["M_ADM_km"] / GM_KM, "nu_R": ext_n["nu_R"],
        })
    k_vals = np.array([c["k_curved_km_inv"] for c in rmax_checks])
    m_vals = np.array([c["M_ADM_Msun"] for c in rmax_checks])
    n_vals = np.array([c["nu_R"] for c in rmax_checks])
    rmax_spread = {
        "dk_over_k": float((k_vals.max() - k_vals.min()) / abs(k_vals[-1])),
        "dM_ADM_over_M_ADM": float((m_vals.max() - m_vals.min()) / m_vals[-1]),
        "dnu_R": float(n_vals.max() - n_vals.min()),
        "bvp_tol": 1e-12,
        "outer_bc_nu": "nu(R_max) = 0.5*ln(1 - 2*m(R_max)/R_max)",
    }
    print("   ", rmax_spread, flush=True)

    print("2. ADM / surface / baryon turning points (minimize_scalar + 7-point quadratic fit)...", flush=True)
    seq_rows = []
    for eos_name, L_km, a_y in SEQUENCE_CASES:
        eos = tp.EOS(eos_name)
        lo_b, hi_b = _seq_bounds(eos_name)
        cache = {}

        def get_st(lr, eos=eos, L_km=L_km, a_y=a_y, cache=cache):
            key = round(float(lr), 8)
            if key not in cache:
                cache[key] = solve_yukawa_action_el_curved_star(eos, float(lr), L_km, a_y)
            return cache[key]

        opt = minimize_scalar(lambda lr: -get_st(lr)["M_ADM_Msun"], bounds=(lo_b, hi_b), method="bounded", options={"xatol": 2e-4})
        x0 = float(opt.x)
        lrs = np.linspace(x0 - 0.0020, x0 + 0.0020, 7)
        sts = [get_st(lr) for lr in lrs]
        c_adm = np.polyfit(lrs, [s["M_ADM_Msun"] for s in sts], 2)
        c_srf = np.polyfit(lrs, [s["M_surface_Msun"] for s in sts], 2)
        c_bar = np.polyfit(lrs, [s["M_baryon_Msun"] for s in sts], 2)
        pk_adm = float(-c_adm[1] / (2.0 * c_adm[0]))
        pk_srf = float(-c_srf[1] / (2.0 * c_srf[0]))
        pk_bar = float(-c_bar[1] / (2.0 * c_bar[0]))
        s_adm = get_st(pk_adm)
        s_srf = get_st(pk_srf)
        s_bar = get_st(pk_bar)
        s15 = get_st(15.0)
        s15_flat = cc.solve_yukawa_action_el_star(eos, 15.0, L_km, alpha_Y=a_y)
        lr14_adm = brentq(lambda lr: get_st(lr)["M_ADM_Msun"] - 1.4, 14.55, pk_adm, xtol=1e-5)
        lr14_srf = brentq(lambda lr: get_st(lr)["M_surface_Msun"] - 1.4, 14.55, pk_srf, xtol=1e-5)
        row = {
            "EOS": eos_name, "L_km": L_km, "alpha_Y": a_y, "alpha_std_poisson": a_y / (4.0 * np.pi),
            "turning_point_method": "minimize_scalar(xatol=2e-4) + 7-point quadratic fit over +/-0.002 dex",
            "M_ADM_max_Msun": float(np.polyval(c_adm, pk_adm)), "R_at_M_ADM_max_km": s_adm["R_km"],
            "log10_rho0c_at_M_ADM_max": pk_adm, "Delta_M_phi_ext_at_M_ADM_max_Msun": s_adm["Delta_M_phi_ext_curved_Msun"],
            "M_surface_at_M_ADM_max_Msun": s_adm["M_surface_Msun"],
            "M_surface_max_Msun": float(np.polyval(c_srf, pk_srf)), "R_at_M_surface_max_km": s_srf["R_km"],
            "log10_rho0c_at_M_surface_max": pk_srf,
            "M_baryon_max_Msun": float(np.polyval(c_bar, pk_bar)), "log10_rho0c_at_M_baryon_max": pk_bar,
            "delta_dex_ADM_vs_baryon_turning": pk_adm - pk_bar,
            "R_at_M_ADM_1p4_km": get_st(lr14_adm)["R_km"], "R_at_M_surface_1p4_km": get_st(lr14_srf)["R_km"],
            "M_ADM_at_15p0_Msun": s15["M_ADM_Msun"], "M_surface_at_15p0_Msun": s15["M_surface_Msun"],
            "M_baryon_at_15p0_Msun": s15["M_baryon_Msun"], "Delta_M_phi_ext_curved_15p0_Msun": s15["Delta_M_phi_ext_curved_Msun"],
            "R_at_15p0_curved_km": s15["R_km"], "f_L_center_15p0_curved": s15["f_center"], "phi_surface_15p0_curved": s15["phi_surface"],
            "dphi_dr_R_curved_15p0": s15["dphi_dr_R_curved"], "dphi_dr_R_flat_approx_15p0": s15["dphi_dr_R_flat_approx"],
            "k_curved_15p0_km_inv": s15["k_curved"], "k_flat_15p0_km_inv": s15["k_flat"], "rel_diff_dphi_R_15p0": s15["rel_diff_dphi_R"],
            "nu_R_15p0": s15["nu_R"],
            "old_flat_robin_M_surface_15p0_Msun": s15_flat["M_Msun"], "old_flat_robin_M_total_15p0_Msun": s15_flat["M_total_Msun"],
            "old_flat_robin_R_15p0_km": s15_flat["R_km"],
        }
        seq_rows.append(row)
        print(f"   {eos_name:4s} L={L_km:4.1f} aY={a_y:.3f}: M_ADM_pk={pk_adm:.6f}, M_B_pk={pk_bar:.6f}, M_srf_pk={pk_srf:.6f}, M_ADM_max={row['M_ADM_max_Msun']:.5f}, R1.4={row['R_at_M_ADM_1p4_km']:.4f}  [{time.time()-t_start:.0f}s]", flush=True)

    df_seq = pd.DataFrame(seq_rows)
    df_seq.to_csv(os.path.join(OUT_DIR, "curved_yukawa_adm_sequences.csv"), index=False)

    print("3. Coupled GL2 Volterra stability (dynamical scalar wave, quasi-static, frozen delta phi_L)...", flush=True)
    stab_rows, zc_summary = [], {}
    for eos_name, L_km, a_y in STABILITY_CASES:
        eos = tp.EOS(eos_name)
        srow = [r for r in seq_rows if r["EOS"] == eos_name and r["L_km"] == L_km and abs(r["alpha_Y"] - a_y) < 1e-9][0]
        pk = srow["log10_rho0c_at_M_ADM_max"]
        pk_mb = srow["log10_rho0c_at_M_baryon_max"]
        lrs_c = [pk - 0.0010, pk, pk + 0.0010]
        dyn = [solve_yukawa_coupled_omega2_gl2(eos, lr, L_km=L_km, alpha_Y=a_y, N=480, dynamical_scalar_wave=True) for lr in lrs_c]
        qs = [solve_yukawa_coupled_omega2_gl2(eos, lr, L_km=L_km, alpha_Y=a_y, N=480, dynamical_scalar_wave=False) for lr in lrs_c]
        lrs_f = [pk + 0.002, pk + 0.010, pk + 0.018] if L_km == 10.0 else [pk - 0.001, pk + 0.0005, pk + 0.002]
        fr = [solve_yukawa_coupled_omega2_gl2(eos, lr, L_km=L_km, alpha_Y=a_y, N=480, include_dphi=False) for lr in lrs_f]
        zc_dyn = _quad_root(lrs_c, [d["omega2_s2"] for d in dyn])
        zc_qs = _quad_root(lrs_c, [d["omega2_s2"] for d in qs])
        zc_fr = _quad_root(lrs_f, [d["omega2_s2"] for d in fr])
        for lr, d, q in zip(lrs_c, dyn, qs):
            stab_rows.append({"EOS": eos_name, "branch": "yukawa_action_el_curved", "L_km": L_km, "alpha_Y": a_y, "log10_rho0c": lr,
                              "M_ADM_Msun": d["M_ADM_Msun"], "M_surface_Msun": d["M_surface_Msun"], "M_baryon_Msun": d["M_baryon_Msun"], "R_km": d["R_km"],
                              "omega0_sq_coupled_dynamical_s2": d["omega2_s2"], "omega0_sq_coupled_quasistatic_s2": q["omega2_s2"],
                              "omega0_sq_frozen_phi_s2": float("nan"), "delta_M_ADM_over_M_ADM_residual": d["delta_M_ADM_over_M_ADM"],
                              "N_gl": d["N_gl"], "log10_rho0c_turning_M_ADM": pk, "log10_rho0c_turning_M_baryon": pk_mb})
        for lr, f in zip(lrs_f, fr):
            stab_rows.append({"EOS": eos_name, "branch": "yukawa_action_el_curved", "L_km": L_km, "alpha_Y": a_y, "log10_rho0c": lr,
                              "M_ADM_Msun": f["M_ADM_Msun"], "M_surface_Msun": f["M_surface_Msun"], "M_baryon_Msun": f["M_baryon_Msun"], "R_km": f["R_km"],
                              "omega0_sq_coupled_dynamical_s2": float("nan"), "omega0_sq_coupled_quasistatic_s2": float("nan"),
                              "omega0_sq_frozen_phi_s2": f["omega2_s2"], "delta_M_ADM_over_M_ADM_residual": f["delta_M_ADM_over_M_ADM"],
                              "N_gl": f["N_gl"], "log10_rho0c_turning_M_ADM": pk, "log10_rho0c_turning_M_baryon": pk_mb})
        key = f"{eos_name}_yukawa_action_el_L{int(L_km)}km_alphaY_{a_y:.4f}"
        zc_summary[key] = {
            "EOS": eos_name, "L_km": L_km, "alpha_Y": a_y,
            "log10_rho0c_turning_M_ADM": pk, "log10_rho0c_turning_M_baryon": pk_mb,
            "log10_rho0c_turning_M_surface": srow["log10_rho0c_at_M_surface_max"],
            "zc_coupled_dynamical_omega2": zc_dyn, "zc_coupled_quasistatic_omega2": zc_qs, "zc_frozen_phi_omega2": zc_fr,
            "delta_dex_dynamical_vs_M_ADM": zc_dyn - pk, "delta_dex_quasistatic_vs_M_ADM": zc_qs - pk,
            "delta_dex_dynamical_vs_M_baryon": zc_dyn - pk_mb, "delta_dex_frozen_vs_M_ADM": zc_fr - pk,
            "mismatch_reduction_factor": abs(zc_fr - pk) / max(abs(zc_dyn - pk), 1e-12),
            "solver": "GL2 piecewise-segment Volterra, N=480, exact linearized 4-component exterior response",
        }
        print(f"   {eos_name:4s} L={L_km:4.1f} aY={a_y:.3f}: pk_ADM={pk:.6f}, zc_dyn={zc_dyn:.6f} ({zc_dyn-pk:+.2e}), zc_qs={zc_qs:.6f}, zc_fr={zc_fr:.6f} ({zc_fr-pk:+.2e})  [{time.time()-t_start:.0f}s]", flush=True)

    df_stab = pd.DataFrame(stab_rows)
    df_stab.to_csv(os.path.join(OUT_DIR, "curved_yukawa_coupled_stability.csv"), index=False)

    print("4. GL2 N-convergence series and Richardson extrapolation (SLy L=10 aY=1, at M_ADM turning)...", flush=True)
    pk_ref = [r for r in seq_rows if r["EOS"] == "SLy" and r["L_km"] == 10.0 and r["alpha_Y"] == 1.0][0]["log10_rho0c_at_M_ADM_max"]
    lr_conv = pk_ref - 0.005
    st_conv = solve_yukawa_action_el_curved_star(eos_sly, lr_conv, L=10.0, alpha_Y=1.0, return_profile=True)
    conv = []
    for N in N_CONVERGENCE:
        rN = solve_yukawa_coupled_omega2_gl2(eos_sly, lr_conv, L_km=10.0, alpha_Y=1.0, N=N, st_eq=st_conv)
        conv.append({"N_target": N, "N_gl": rN["N_gl"], "omega2_s2": rN["omega2_s2"]})
        print(f"   N={N:5d} (N_gl={rN['N_gl']}): omega2 = {rN['omega2_s2']:.6f} s^-2", flush=True)
    Ns = np.array([c["N_gl"] for c in conv], float)
    ws = np.array([c["omega2_s2"] for c in conv])
    # p_eff from the three finest levels (non-uniform ratios solved by brentq)
    n1, n2, n3 = Ns[-3:]
    w1, w2, w3 = ws[-3:]
    def g(p):
        return (w1 - w2) / (w2 - w3) - (n1 ** -p - n2 ** -p) / (n2 ** -p - n3 ** -p)
    try:
        p_eff = float(brentq(g, 0.5, 10.0))
    except ValueError:
        p_eff = float("nan")
    w_rich = float(w3 + (w3 - w2) / ((n3 / n2) ** p_eff - 1.0)) if np.isfinite(p_eff) else float("nan")
    for c in conv:
        c["rel_err_vs_richardson"] = float(abs(c["omega2_s2"] - w_rich) / abs(w_rich)) if np.isfinite(w_rich) else float("nan")
    print(f"   p_eff = {p_eff:.3f}, omega2_Richardson = {w_rich:.6f}", flush=True)

    hess_diag = compute_hessian_sturm_liouville_diagnostics(eos_sly, pk_ref, L_km=10.0, alpha_Y=1.0)
    summary_payload = {
        "rmax_independence_checks_SLy_15p0_L10km": rmax_checks,
        "rmax_independence_spread_4comp_exterior": rmax_spread,
        "curved_yukawa_adm_sequences": seq_rows,
        "curved_yukawa_stability_zero_crossings": zc_summary,
        "gl2_N_convergence_SLy_L10_aY1": {"log10_rho0c": lr_conv, "series": conv, "p_eff": p_eff, "omega2_richardson_s2": w_rich},
        "coupled_hessian_sturm_liouville_diagnostics_SLy_turning": hess_diag,
        "runtime_s": time.time() - t_start,
    }
    with open(os.path.join(OUT_DIR, "curved_yukawa_audit_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary_payload, fh, indent=2, default=float)
    import shutil
    for sub in ["code", "csv", os.path.join("data", "v23")]:
        d = os.path.join(OUT_DIR, sub)
        if os.path.isdir(d):
            for fn in ["curved_yukawa_adm_sequences.csv", "curved_yukawa_coupled_stability.csv", "curved_yukawa_audit_summary.json"]:
                shutil.copy2(os.path.join(OUT_DIR, fn), os.path.join(d, fn))
    shutil.copy2(os.path.abspath(__file__), os.path.join(OUT_DIR, "code", "curved_yukawa_adm_and_stability.py")) if os.path.isdir(os.path.join(OUT_DIR, "code")) and os.path.abspath(__file__) != os.path.join(OUT_DIR, "code", "curved_yukawa_adm_and_stability.py") else None
    print(f"5. Saved CSV/JSON outputs [{time.time()-t_start:.0f}s]", flush=True)


if __name__ == "__main__":
    run_all_curved_yukawa_audits()
