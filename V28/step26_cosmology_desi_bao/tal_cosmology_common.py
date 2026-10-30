#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified Cosmological and Scalar-Tensor Background Module for V27 (Step 26).
Used by all Step 26 scripts to ensure 100% identical background evolution,
exact E^2(0) = 1.0 normalization, self-consistent sound horizons, and unified
SI / natural unit scalar-tensor normalization.
"""

import os
import sys
import math
import hashlib

# Safe fallback if invoked with a system Python where user site-packages is not loaded
try:
    import scipy  # noqa: F401
except ImportError:
    for _cand in [
        os.path.expanduser("~/.local/lib/python3.13/site-packages"),
        "/workspace/scratch/venv_openitems/lib/python3.13/site-packages",
    ]:
        if os.path.isdir(_cand) and _cand not in sys.path:
            sys.path.append(_cand)

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize, brentq
from scipy.stats import chi2 as chi2_dist

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# -----------------------------------------------------------------------------
# 1. Fundamental Physical Constants and V26 / Step 25 TAL Parameters (SI Units)
# -----------------------------------------------------------------------------
G_SI = 6.67430e-11
C_SI = 299792458.0
C_KM_S = C_SI / 1000.0
MPC_M = 3.085677581491367e22
GYR_S = 365.25 * 86400.0 * 1.0e9

L_KM = 10.0
L_M = L_KM * 1000.0
L_MPC = L_M / MPC_M                          # 3.2407792894443654e-19 Mpc
M_PHI_C_S_INV = C_SI / L_M                   # 29979.2458 s^-1

ALPHA_Y = 1.0
BETA_YUK = 1.0 / (4.0 * math.pi)             # 0.07957747154594767

# Unified Action & Sign Convention (Review Section 4.1 & 5.6):
# S = \int d^4x \sqrt{-g} [ (c^4 / 16\pi G) R
#     - (c^4 / 4\pi G \beta_Yuk) (1/2 g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi + \phi^2 / (2 L^2))
#     - \rho_{vac,eff} c^2 ] + S_m[\Psi_m, A^2(\phi) g_{\mu\nu}],
# where A(\phi) = exp(+\alpha_Y \phi), so g_{\mu\nu}^{(J)} = exp(+2 \alpha_Y \phi) g_{\mu\nu}^{(E)}.
# With (-+++) signature, T^{\mu\nu} = (2/\sqrt{-g}) \delta S_m / \delta g_{\mu\nu} and
# T_m = g_{\mu\nu} T^{\mu\nu} = -\rho_m c^2 + 3 p_m < 0.
# Variation w.r.t. \phi gives the covariant wave equation in m^-2:
#   \Box \phi - L^-2 \phi = -(4\pi G / c^4) \alpha_Y \beta_Yuk T_m
#                         = +(4\pi G / c^2) \alpha_Y \beta_Yuk (\rho_m - 3 p_m / c^2)   [m^-2]
# - In the static point-mass limit (\Box\phi = \nabla^2\phi):
#   \nabla^2\phi - L^-2\phi = +(4\pi G / c^2) \alpha_Y \beta_Yuk M \delta^{(3)}(r)
#   ==> \phi(r) = -\alpha_Y \beta_Yuk (G M / (c^2 r)) exp(-r/L) < 0,
#   and local proper time is d\tau_J = A(\phi) d\tau_E = exp(\alpha_Y \phi) d\tau_E = exp(-\alpha_Y |\phi|) d\tau_E.
# - In homogeneous FLRW with (-+++) signature (\Box\phi = -c^-2 (\ddot{\phi} + 3 H \dot{\phi})),
#   multiplying by -c^2 gives the cosmic-time Klein-Gordon equation in s^-2:
#   \ddot{\phi} + 3 H \dot{\phi} + (c^2 / L^2) \phi = -4\pi G \alpha_Y \beta_Yuk (\rho_m - 3 p_m / c^2)   [s^-2]
# - Since 4\pi \beta_Yuk = 1, the quasi-static adiabatic attractor (\ddot{\phi}, 3H\dot{\phi} -> 0) is:
#   \bar{\phi}_{attr}(z) = -\alpha_Y \gamma_L (\rho_m(z) - 3 p_m(z)/c^2) < 0,
#   with \gamma_L = 4\pi \beta_Yuk G L^2 / c^2 = G L^2 / c^2 = 7.426160270019091e-20 m^3/kg.
GAMMA_L = 4.0 * math.pi * BETA_YUK * G_SI * (L_M ** 2) / (C_SI ** 2)   # 7.426160270019091e-20 m^3/kg
GAMMA_L_2X = 2.0 * GAMMA_L                                             # 1.4852320540038182e-19 m^3/kg

# -----------------------------------------------------------------------------
# 2. Cosmological Baseline Parameters (Step 25 M4 and Planck 2018)
# -----------------------------------------------------------------------------
H0_STEP25 = 67.40
H0_PLANCK = 67.36
# Review F05: 72.80 +/- 1.00 is an ASSUMED LOCAL H0 BENCHMARK PRIOR (a chosen target value).
# No specific SH0ES publication or data release is claimed for it (e.g. Riess et al. 2022,
# ApJL 934, L7 quotes 73.04 +/- 1.04, which is a different number).
H0_LOCAL = 72.80
H0_LOCAL_SIGMA = 1.00

H674 = H0_STEP25 / 100.0
H0_S_INV = (H0_STEP25 * 1000.0) / MPC_M
RHO_CRIT_0_SI = 3.0 * (H0_S_INV ** 2) / (8.0 * math.pi * G_SI)

OMEGA_B_H2 = 0.02237
OMEGA_C_H2 = 0.12000
OMEGA_CB_H2 = OMEGA_B_H2 + OMEGA_C_H2        # 0.14237
OMEGA_NU_H2 = 0.000645
OMEGA_M_H2 = OMEGA_CB_H2 + OMEGA_NU_H2       # 0.143015
OMEGA_GAMMA_H2 = 2.47282e-5
N_EFF = 3.046

OMEGA_B_0 = OMEGA_B_H2 / (H674 ** 2)
OMEGA_C_0 = OMEGA_C_H2 / (H674 ** 2)
OMEGA_CB_0 = OMEGA_B_0 + OMEGA_C_0
OMEGA_NU_0 = OMEGA_NU_H2 / (H674 ** 2)
OMEGA_M_0 = OMEGA_CB_0 + OMEGA_NU_0          # 0.314819598482

OMEGA_GAMMA_0 = OMEGA_GAMMA_H2 / (H674 ** 2)
OMEGA_NU_REL_1_0 = OMEGA_GAMMA_0 * 0.22710731766 * (N_EFF / 3.0)
OMEGA_R_MASSLESS_0 = OMEGA_GAMMA_0 + 2.0 * OMEGA_NU_REL_1_0

def neutrino_massive_density_E2(z, omega_nu_0=OMEGA_NU_0, omega_nu_rel_1=OMEGA_NU_REL_1_0):
    """Smooth massive-neutrino transition from relativistic (z >> z_nr) to non-relativistic (z -> 0)."""
    if omega_nu_0 <= 0.0:
        return omega_nu_rel_1 * ((1.0 + z) ** 4)
    a = 1.0 / (1.0 + z)
    p = 1.83
    ratio_a = (omega_nu_0 / omega_nu_rel_1) * a
    f_trans = (1.0 + (ratio_a ** p)) ** (1.0 / p)
    return omega_nu_rel_1 * f_trans * ((1.0 + z) ** 4)

# Exact z=0 massive neutrino contribution and exact flatness closure E^2(0) == 1.0:
OMEGA_NU_INTERP_Z0 = neutrino_massive_density_E2(0.0, OMEGA_NU_0, OMEGA_NU_REL_1_0)
NEUTRINO_INTERP_Z0_OFFSET = OMEGA_NU_INTERP_Z0 - OMEGA_NU_0   # 1.354616e-7
OMEGA_LAMBDA_FLAT_UNNORM_0 = 1.0 - OMEGA_M_0 - OMEGA_R_MASSLESS_0
OMEGA_LAMBDA_FLAT_0 = 1.0 - OMEGA_CB_0 - OMEGA_R_MASSLESS_0 - OMEGA_NU_INTERP_Z0

# Quasi-static adiabatic scalar attractor and energy density at z=0:
RHO_M_0_SI = OMEGA_M_0 * RHO_CRIT_0_SI
PHI_ATTR_0 = -ALPHA_Y * GAMMA_L * RHO_M_0_SI                  # -1.994897e-46
ABS_PHI_ATTR_0 = abs(PHI_ATTR_0)                              # +1.994897e-46
# From \rho_{\phi,mass} = (c^2 / (8\pi G \beta_Yuk)) [ (\dot{\phi}/c)^2 + (\phi/L)^2 ]
# and \rho_{crit,0} = 3 H_0^2 / (8\pi G):
# \Omega_{\phi,attr}(0) = [ (c/L)^2 \phi_{attr}(0)^2 ] / (3 \beta_Yuk H_0^2)
OMEGA_PHI_ATTR_0 = ((M_PHI_C_S_INV ** 2) * (PHI_ATTR_0 ** 2)) / (3.0 * BETA_YUK * (H0_S_INV ** 2))

Z_STAR_PLANCK = 1089.80
Z_DRAG_PLANCK = 1059.94
THETA_STAR_PLANCK_100 = 1.04109
THETA_STAR_PLANCK_ERR_100 = 0.00030
RD_PLANCK_MPC = 147.09
RD_PLANCK_ERR_MPC = 0.26
RS_PLANCK_MPC = 144.43
RS_PLANCK_ERR_MPC = 0.26

# Illustrative 3D Gaussian prior on (omega_b, omega_cb, 100*theta_*)
# Marginal means and sigmas from Planck 2018 VI Table 2 (base_plikHM_TTTEEE_lowl_lowE_lensing);
# local correlation matrix is an illustrative Gaussian approximation (Review Section 4.5).
PLANCK_3D_MEAN = np.array([0.02237, 0.14240, 1.04109])
PLANCK_3D_SIGMA = np.array([0.00015, 0.00121, 0.00030])
PLANCK_3D_CORR = np.array([
    [ 1.00, -0.62,  0.45],
    [-0.62,  1.00, -0.41],
    [ 0.45, -0.41,  1.00],
])
PLANCK_3D_COV = np.outer(PLANCK_3D_SIGMA, PLANCK_3D_SIGMA) * PLANCK_3D_CORR
PLANCK_3D_INV_COV = np.linalg.inv(PLANCK_3D_COV)

def planck_3d_illustrative_chi2(om_b_h2, om_cb_h2, theta100):
    vec = np.array([om_b_h2, om_cb_h2, theta100])
    diff = vec - PLANCK_3D_MEAN
    return float(diff @ PLANCK_3D_INV_COV @ diff)

def E2_general(z, omega_cb=OMEGA_CB_0, omega_nu_0=OMEGA_NU_0,
               omega_r_massless=OMEGA_R_MASSLESS_0, omega_nu_rel_1=OMEGA_NU_REL_1_0,
               omega_k=0.0, omega_de=OMEGA_LAMBDA_FLAT_0, w0=-1.0, wa=0.0):
    zp1 = 1.0 + z
    if w0 == -1.0 and wa == 0.0:
        de_term = omega_de
    else:
        de_term = omega_de * (zp1 ** (3.0 * (1.0 + w0 + wa))) * math.exp(-3.0 * wa * z / zp1)
    return (
        omega_r_massless * (zp1 ** 4)
        + neutrino_massive_density_E2(z, omega_nu_0, omega_nu_rel_1)
        + omega_cb * (zp1 ** 3)
        + omega_k * (zp1 ** 2)
        + de_term
    )

def make_model_background(H0_val=H0_STEP25, om_m_val=None, om_b_h2=OMEGA_B_H2,
                          om_nu_h2=OMEGA_NU_H2, omega_k=0.0, exact_z0_norm=True):
    """Build self-consistent background dictionary and H(z) callable for any (H0, Omega_m, omega_b, Omega_K)."""
    h = H0_val / 100.0
    om_b = om_b_h2 / (h ** 2)
    om_nu = om_nu_h2 / (h ** 2)
    if om_m_val is None:
        om_c = OMEGA_C_H2 / (h ** 2)
        om_cb = om_b + om_c
        om_m = om_cb + om_nu
    else:
        om_m = om_m_val
        om_cb = om_m - om_nu
        om_c = om_cb - om_b
    om_cb_h2 = om_cb * (h ** 2)
    om_g = OMEGA_GAMMA_H2 / (h ** 2)
    om_nu_r1 = om_g * 0.22710731766 * (N_EFF / 3.0)
    om_r_m = om_g + 2.0 * om_nu_r1
    om_nu_z0 = neutrino_massive_density_E2(0.0, om_nu, om_nu_r1)
    if exact_z0_norm:
        om_de = 1.0 - omega_k - om_cb - om_r_m - om_nu_z0
    else:
        om_de = 1.0 - omega_k - (om_cb + om_nu) - om_r_m

    def E2_fn(z):
        return E2_general(z, om_cb, om_nu, om_r_m, om_nu_r1, omega_k, om_de)

    def Hz_fn(z):
        return H0_val * math.sqrt(E2_fn(z))

    return {
        "H0": H0_val,
        "h": h,
        "Omega_b_0": om_b,
        "Omega_c_0": om_c,
        "Omega_cb_0": om_cb,
        "Omega_nu_0": om_nu,
        "Omega_m_0": om_m,
        "Omega_gamma_0": om_g,
        "Omega_r_massless_0": om_r_m,
        "Omega_nu_rel_1_0": om_nu_r1,
        "Omega_k_0": omega_k,
        "Omega_de_0": om_de,
        "omega_b_h2": om_b_h2,
        "omega_cb_h2": om_cb_h2,
        "omega_m_h2": om_m * (h ** 2),
        "E2_fn": E2_fn,
        "Hz_fn": Hz_fn,
    }

def sound_speed_over_c(z, omega_b=OMEGA_B_0, omega_gamma=OMEGA_GAMMA_0):
    r_b = (3.0 * omega_b) / (4.0 * omega_gamma * (1.0 + z))
    return 1.0 / math.sqrt(3.0 * (1.0 + r_b))

def z_star_hu_sugiyama(om_b_h2=OMEGA_B_H2, om_m_h2=OMEGA_M_H2):
    """Hu & Sugiyama (1996, ApJ 471, 542) fitting formula for recombination redshift z_*."""
    g1 = (0.0783 * (om_b_h2 ** -0.238)) / (1.0 + 39.5 * (om_b_h2 ** 0.763))
    g2 = 0.560 / (1.0 + 21.1 * (om_b_h2 ** 1.81))
    return 1048.0 * (1.0 + 0.00124 * (om_b_h2 ** -0.738)) * (1.0 + g1 * (om_m_h2 ** g2))

def z_drag_eisenstein_hu(om_b_h2=OMEGA_B_H2, om_m_h2=OMEGA_M_H2):
    """Eisenstein & Hu (1998, ApJ 496, 605) Eq. (4) fitting formula for baryon-drag redshift z_drag."""
    b1 = 0.313 * (om_m_h2 ** -0.419) * (1.0 + 0.607 * (om_m_h2 ** 0.674))
    b2 = 0.238 * (om_m_h2 ** 0.223)
    return 1291.0 * ((om_m_h2 ** 0.251) / (1.0 + 0.659 * (om_m_h2 ** 0.828))) * (1.0 + b1 * (om_b_h2 ** b2))

def compute_sound_horizon_for_bg(z_lower, bg_dict):
    c_over_H0 = C_KM_S / bg_dict["H0"]
    om_b = bg_dict["Omega_b_0"]
    om_g = bg_dict["Omega_gamma_0"]
    E2_fn = bg_dict["E2_fn"]
    val, _ = quad(
        lambda z: sound_speed_over_c(z, om_b, om_g) / math.sqrt(E2_fn(z)),
        z_lower, np.inf, epsabs=1e-11, epsrel=1e-11, limit=300
    )
    return c_over_H0 * val

def DM_of_z(z, Hz_fn, omega_k=0.0, H0_for_curv=H0_STEP25):
    chi, _ = quad(lambda zp: C_KM_S / Hz_fn(zp), 0.0, z, epsabs=1e-11, epsrel=1e-11, limit=250)
    if abs(omega_k) < 1e-12:
        return chi
    dh0 = C_KM_S / H0_for_curv
    sqrt_ok = math.sqrt(abs(omega_k))
    x = sqrt_ok * (chi / dh0)
    if omega_k > 0.0:
        return (dh0 / sqrt_ok) * math.sinh(x)
    return (dh0 / sqrt_ok) * math.sin(x)

# Baseline Step 25 M4 and Planck 2018 background objects
BG_STEP25_M4 = make_model_background(H0_STEP25, OMEGA_M_0, OMEGA_B_H2, OMEGA_NU_H2, 0.0, exact_z0_norm=True)
BG_STEP25_M4_UNNORM = make_model_background(H0_STEP25, OMEGA_M_0, OMEGA_B_H2, OMEGA_NU_H2, 0.0, exact_z0_norm=False)
BG_PLANCK_67P36 = make_model_background(H0_PLANCK, None, OMEGA_B_H2, OMEGA_NU_H2, 0.0, exact_z0_norm=True)

RS_ZSTAR = compute_sound_horizon_for_bg(Z_STAR_PLANCK, BG_STEP25_M4)
RD_ZDRAG = compute_sound_horizon_for_bg(Z_DRAG_PLANCK, BG_STEP25_M4)

# -----------------------------------------------------------------------------
# 3. DESI DR1 BAO Data (Published Table 1 + Full-Precision Cobaya 12-Vector & Covariance)
# -----------------------------------------------------------------------------
# Note on single-tracer rows (Review Section 6 & 5.5):
#   - 4a_LRG3_sep (z_eff=0.922): DESI III (arXiv:2404.03000) Table 18 post-recon LRG3, DV/rd = 19.82 +/- 0.18
#   - 4b_ELG1_955 (z_eff=0.955): DESI III Table 18 post-recon ELG1, DV/rd = 20.01 +/- 0.41 (excluded from 12-pt primary fit)
#   - 4c_ELG1_845 (z_eff=0.845): CHOSEN ILLUSTRATION without a verified primary source (review P01); it is NOT a DESI
#     measurement (not in DESI III Table 18 or DESI VI Table 1) and is excluded from the 12-pt primary fit
DESI_BINS = [
    {"bin_id": "1_BGS",       "tracer": "BGS",             "z_range": "0.1-0.4",   "z_eff": 0.295, "primary_12pt": True,  "mode": "iso",   "DV_rd": 7.93,  "DV_err": 0.15},
    {"bin_id": "2_LRG1",      "tracer": "LRG1",            "z_range": "0.4-0.6",   "z_eff": 0.510, "primary_12pt": True,  "mode": "aniso", "DM_rd": 13.62, "DM_err": 0.25, "DH_rd": 20.98, "DH_err": 0.61, "r": -0.445, "DV_rd": 12.56, "DV_err": 0.15},
    {"bin_id": "3_LRG2",      "tracer": "LRG2",            "z_range": "0.6-0.8",   "z_eff": 0.706, "primary_12pt": True,  "mode": "aniso", "DM_rd": 16.85, "DM_err": 0.32, "DH_rd": 20.08, "DH_err": 0.60, "r": -0.420, "DV_rd": 15.90, "DV_err": 0.20},
    {"bin_id": "4_LRG3_ELG1", "tracer": "LRG3+ELG1",       "z_range": "0.8-1.1",   "z_eff": 0.930, "primary_12pt": True,  "mode": "aniso", "DM_rd": 21.71, "DM_err": 0.28, "DH_rd": 17.88, "DH_err": 0.35, "r": -0.389, "DV_rd": 19.86, "DV_err": 0.17},
    {"bin_id": "4a_LRG3_sep", "tracer": "LRG3 (sep.)",     "z_range": "0.8-1.1",   "z_eff": 0.922, "primary_12pt": False, "mode": "aniso", "DM_rd": 21.81, "DM_err": 0.31, "DH_rd": 17.83, "DH_err": 0.38, "r": -0.393, "DV_rd": 19.82, "DV_err": 0.18},
    {"bin_id": "4b_ELG1_955", "tracer": "ELG1 (z=0.955)",  "z_range": "0.8-1.1",   "z_eff": 0.955, "primary_12pt": False, "mode": "iso",   "DV_rd": 20.01, "DV_err": 0.41},
    {"bin_id": "4c_ELG1_845", "tracer": "ELG1 (z=0.845, illustration)", "z_range": "0.8-1.1",   "z_eff": 0.845, "primary_12pt": False, "mode": "iso",   "DV_rd": 18.33, "DV_err": 0.62},
    {"bin_id": "5_ELG2",      "tracer": "ELG2",            "z_range": "1.1-1.6",   "z_eff": 1.317, "primary_12pt": True,  "mode": "aniso", "DM_rd": 27.79, "DM_err": 0.69, "DH_rd": 13.82, "DH_err": 0.42, "r": -0.444, "DV_rd": 24.13, "DV_err": 0.36},
    {"bin_id": "6_QSO",       "tracer": "QSO",             "z_range": "0.8-2.1",   "z_eff": 1.491, "primary_12pt": True,  "mode": "iso",   "DV_rd": 26.07, "DV_err": 0.67},
    {"bin_id": "7_Lya_QSO",   "tracer": "Lya QSO",         "z_range": "1.77-4.16", "z_eff": 2.330, "primary_12pt": True,  "mode": "aniso", "DM_rd": 39.71, "DM_err": 0.94, "DH_rd": 8.52,  "DH_err": 0.17, "r": -0.477},
]

def load_cobaya_desi_full():
    mean_path = os.path.join(OUT_DIR, "desi_2024_gaussian_bao_ALL_GCcomb_mean.txt")
    cov_path = os.path.join(OUT_DIR, "desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
    z_list, val_list, q_list = [], [], []
    with open(mean_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            z_list.append(float(parts[0]))
            val_list.append(float(parts[1]))
            q_list.append(parts[2])
    cov = np.loadtxt(cov_path)
    inv_cov = np.linalg.inv(cov)
    return np.array(z_list), np.array(val_list), q_list, cov, inv_cov

COBAYA_Z, COBAYA_MEAN, COBAYA_Q, COBAYA_COV, COBAYA_INV_COV = load_cobaya_desi_full()

def bao_chi2_table1(Hz_fn, rd_val=RD_ZDRAG, omega_k=0.0, H0_for_curv=H0_STEP25):
    c2 = 0.0
    for b in DESI_BINS:
        if not b["primary_12pt"]:
            continue
        z = b["z_eff"]
        DM_rd = DM_of_z(z, Hz_fn, omega_k, H0_for_curv) / rd_val
        DH_rd = (C_KM_S / Hz_fn(z)) / rd_val
        DV_rd = (z * (DM_rd ** 2) * DH_rd) ** (1.0 / 3.0)
        if b["mode"] == "iso":
            c2 += ((DV_rd - b["DV_rd"]) / b["DV_err"]) ** 2
        else:
            dM = (DM_rd - b["DM_rd"]) / b["DM_err"]
            dH = (DH_rd - b["DH_rd"]) / b["DH_err"]
            r = b["r"]
            c2 += (dM ** 2 - 2.0 * r * dM * dH + dH ** 2) / (1.0 - r ** 2)
    return float(c2)

def bao_chi2_cobaya_full(Hz_fn, rd_val=RD_ZDRAG, omega_k=0.0, H0_for_curv=H0_STEP25):
    pred = np.zeros(len(COBAYA_Z))
    cache = {}
    for i, (zv, qname) in enumerate(zip(COBAYA_Z, COBAYA_Q)):
        if zv not in cache:
            dm = DM_of_z(zv, Hz_fn, omega_k, H0_for_curv) / rd_val
            dh = (C_KM_S / Hz_fn(zv)) / rd_val
            dv = (zv * (dm ** 2) * dh) ** (1.0 / 3.0)
            cache[zv] = (dm, dh, dv)
        dm, dh, dv = cache[zv]
        if qname == "DM_over_rs":
            pred[i] = dm
        elif qname == "DH_over_rs":
            pred[i] = dh
        elif qname == "DV_over_rs":
            pred[i] = dv
    diff = pred - COBAYA_MEAN
    return float(diff @ COBAYA_INV_COV @ diff)

def get_unit_bao_pred_for_E_fn(E_fn):
    """Predict BAO observables at unit scale P = H0 * r_d = 1 km/s."""
    pred_u = np.zeros(len(COBAYA_Z))
    cache = {}
    for i, (zv, qname) in enumerate(zip(COBAYA_Z, COBAYA_Q)):
        if zv not in cache:
            dm_u, _ = quad(lambda zp: C_KM_S / E_fn(zp), 0.0, zv, epsabs=1e-12, epsrel=1e-12, limit=300)
            dh_u = C_KM_S / E_fn(zv)
            dv_u = (zv * (dm_u ** 2) * dh_u) ** (1.0 / 3.0)
            cache[zv] = (dm_u, dh_u, dv_u)
        dm_u, dh_u, dv_u = cache[zv]
        if qname == "DM_over_rs":
            pred_u[i] = dm_u
        elif qname == "DH_over_rs":
            pred_u[i] = dh_u
        elif qname == "DV_over_rs":
            pred_u[i] = dv_u
    return pred_u

def compute_bao_profile_controls():
    """
    Compute the 4 independent BAO profile controls (Review Section 4.6 & 6)
    plus inverse distance ladder derived quantities.
    """
    # 1. Fixed Step 25 M4 baseline
    P_m4 = H0_STEP25 * RD_ZDRAG
    chi2_m4_full = bao_chi2_cobaya_full(BG_STEP25_M4["Hz_fn"], RD_ZDRAG)
    chi2_m4_tab1 = bao_chi2_table1(BG_STEP25_M4["Hz_fn"], RD_ZDRAG)
    p_m4_full = float(chi2_dist.sf(chi2_m4_full, 12))

    # 2. Fixed Step 25 M4 shape, rigid H0 = 72.80 + unchanged rd = 147.092243 Mpc
    P_rigid = H0_LOCAL * RD_ZDRAG
    Hz_rigid = lambda z: (H0_LOCAL / H0_STEP25) * BG_STEP25_M4["Hz_fn"](z)
    chi2_rig_full = bao_chi2_cobaya_full(Hz_rigid, RD_ZDRAG)
    chi2_rig_tab1 = bao_chi2_table1(Hz_rigid, RD_ZDRAG)
    p_rig_full = float(chi2_dist.sf(chi2_rig_full, 12))

    # 3. Fixed Step 25 M4 shape (Om = 0.314820), P = H0 * rd profiled
    E_m4 = lambda z: math.sqrt(BG_STEP25_M4["E2_fn"](z))
    u_m4 = get_unit_bao_pred_for_E_fn(E_m4)
    u_inv_u = float(u_m4 @ COBAYA_INV_COV @ u_m4)
    u_inv_y = float(u_m4 @ COBAYA_INV_COV @ COBAYA_MEAN)
    alpha_best_m4 = u_inv_y / u_inv_u
    sigma_alpha_m4 = 1.0 / math.sqrt(u_inv_u)
    P_best_m4 = 1.0 / alpha_best_m4
    sigma_P_m4 = sigma_alpha_m4 / (alpha_best_m4 ** 2)
    diff_m4_prof = (u_m4 / P_best_m4) - COBAYA_MEAN
    chi2_m4_prof_full = float(diff_m4_prof @ COBAYA_INV_COV @ diff_m4_prof)
    Hz_m4_prof = lambda z: (P_best_m4 / RD_ZDRAG) * math.sqrt(BG_STEP25_M4["E2_fn"](z))
    chi2_m4_prof_tab1 = bao_chi2_table1(Hz_m4_prof, RD_ZDRAG)
    p_m4_prof_full = float(chi2_dist.sf(chi2_m4_prof_full, 11))

    # 4. Flat LCDM 2-parameter shape control (Om and P profiled)
    def _chi2_flat_2p(om):
        if om <= 0.05 or om >= 0.95:
            return 1e9, 0.0, None
        E_fn = lambda z: math.sqrt(om * ((1.0 + z) ** 3) + (1.0 - om))
        u = get_unit_bao_pred_for_E_fn(E_fn)
        uiu = float(u @ COBAYA_INV_COV @ u)
        uiy = float(u @ COBAYA_INV_COV @ COBAYA_MEAN)
        a_best = uiy / uiu
        p_best = 1.0 / a_best
        d = (u * a_best) - COBAYA_MEAN
        c2 = float(d @ COBAYA_INV_COV @ d)
        return c2, p_best, u

    res_2p = minimize(lambda x: _chi2_flat_2p(x[0])[0], [0.294], method="Nelder-Mead",
                      options={"xatol": 1e-10, "fatol": 1e-12})
    om_mle = float(res_2p.x[0])
    chi2_2p_full, P_best_2p, _ = _chi2_flat_2p(om_mle)
    Hz_2p = lambda z: (P_best_2p / RD_ZDRAG) * math.sqrt(om_mle * ((1.0 + z) ** 3) + (1.0 - om_mle))
    chi2_2p_tab1 = bao_chi2_table1(Hz_2p, RD_ZDRAG)
    p_2p_full = float(chi2_dist.sf(chi2_2p_full, 10))

    # 2x2 Fisher matrix at (om_mle, P_best_2p)
    def _vec_pred_2p(om, P):
        E_fn = lambda z: math.sqrt(om * ((1.0 + z) ** 3) + (1.0 - om))
        return get_unit_bao_pred_for_E_fn(E_fn) / P

    eps_om, eps_P = 1e-5, 1e-2
    dv_dom = (_vec_pred_2p(om_mle + eps_om, P_best_2p) - _vec_pred_2p(om_mle - eps_om, P_best_2p)) / (2.0 * eps_om)
    dv_dP = (_vec_pred_2p(om_mle, P_best_2p + eps_P) - _vec_pred_2p(om_mle, P_best_2p - eps_P)) / (2.0 * eps_P)
    J = np.vstack([dv_dom, dv_dP])
    F_2p = J @ COBAYA_INV_COV @ J.T
    cov_2p = np.linalg.inv(F_2p)
    sigma_om_2p = math.sqrt(cov_2p[0, 0])
    sigma_P_2p = math.sqrt(cov_2p[1, 1])

    # Derived inverse-ladder quantities:
    rd_req_m4_fixed = P_m4 / H0_LOCAL                          # 136.1816 Mpc
    rd_req_m4_prof = P_best_m4 / H0_LOCAL                      # 137.7699 Mpc
    sigma_rd_m4_prof_fixedH0 = sigma_P_m4 / H0_LOCAL           # 0.6515 Mpc
    rd_req_2p = P_best_2p / H0_LOCAL                           # 140.0229 Mpc
    sigma_rd_2p_fixedH0 = sigma_P_2p / H0_LOCAL                # 1.7372 Mpc
    sigma_rd_2p_propH0 = rd_req_2p * math.sqrt((sigma_P_2p / P_best_2p) ** 2 + (H0_LOCAL_SIGMA / H0_LOCAL) ** 2)

    H0_bao_planck_rd = P_best_2p / RD_PLANCK_MPC               # 69.3023 km/s/Mpc
    sigma_H0_bao_planck_rd = H0_bao_planck_rd * math.sqrt(
        (sigma_P_2p / P_best_2p) ** 2 + (RD_PLANCK_ERR_MPC / RD_PLANCK_MPC) ** 2
    )

    # Joint 3-parameter (Om, H0, rd) profile against DESI BAO + Planck rd prior + local H0 prior
    def _joint_chi2_3p(x):
        om, h0, rd = x
        if om <= 0.05 or om >= 0.95 or h0 <= 40.0 or rd <= 100.0:
            return 1e9
        E_fn = lambda z: math.sqrt(om * ((1.0 + z) ** 3) + (1.0 - om))
        u = get_unit_bao_pred_for_E_fn(E_fn)
        d = (u / (h0 * rd)) - COBAYA_MEAN
        c2_bao = float(d @ COBAYA_INV_COV @ d)
        c2_rd = ((rd - RD_PLANCK_MPC) / RD_PLANCK_ERR_MPC) ** 2
        c2_h0 = ((h0 - H0_LOCAL) / H0_LOCAL_SIGMA) ** 2
        return c2_bao + c2_rd + c2_h0

    res_joint = minimize(_joint_chi2_3p, [0.294, 70.5, 146.9], method="Nelder-Mead",
                         options={"xatol": 1e-8, "fatol": 1e-10})
    chi2_joint_min = float(res_joint.fun)
    dchi2_joint = chi2_joint_min - chi2_2p_full
    sigma_joint_equiv = math.sqrt(dchi2_joint)

    return {
        "P_m4": P_m4,
        "chi2_m4_full": chi2_m4_full,
        "chi2_m4_tab1": chi2_m4_tab1,
        "p_m4_full": p_m4_full,
        "P_rigid": P_rigid,
        "chi2_rig_full": chi2_rig_full,
        "chi2_rig_tab1": chi2_rig_tab1,
        "p_rig_full": p_rig_full,
        "P_best_m4": P_best_m4,
        "sigma_P_m4": sigma_P_m4,
        "chi2_m4_prof_full": chi2_m4_prof_full,
        "chi2_m4_prof_tab1": chi2_m4_prof_tab1,
        "p_m4_prof_full": p_m4_prof_full,
        "om_mle_2p": om_mle,
        "sigma_om_2p": sigma_om_2p,
        "P_best_2p": P_best_2p,
        "sigma_P_2p": sigma_P_2p,
        "chi2_2p_full": chi2_2p_full,
        "chi2_2p_tab1": chi2_2p_tab1,
        "p_2p_full": p_2p_full,
        "rd_req_m4_fixed": rd_req_m4_fixed,
        "rd_req_m4_prof": rd_req_m4_prof,
        "sigma_rd_m4_prof_fixedH0": sigma_rd_m4_prof_fixedH0,
        "rd_req_2p": rd_req_2p,
        "sigma_rd_2p_fixedH0": sigma_rd_2p_fixedH0,
        "sigma_rd_2p_propH0": sigma_rd_2p_propH0,
        "H0_bao_planck_rd": H0_bao_planck_rd,
        "sigma_H0_bao_planck_rd": sigma_H0_bao_planck_rd,
        "chi2_joint_min": chi2_joint_min,
        "dchi2_joint": dchi2_joint,
        "sigma_joint_equiv": sigma_joint_equiv,
        "joint_best_params": [float(v) for v in res_joint.x],
    }

if __name__ == "__main__":
    print(f"E2_M4(0) exact = {BG_STEP25_M4['E2_fn'](0.0):.16f}")
    print(f"E2_M4_unnorm(0) = {BG_STEP25_M4_UNNORM['E2_fn'](0.0):.16f}")
    print(f"GAMMA_L = {GAMMA_L:.6e} m^3/kg, PHI_ATTR_0 = {PHI_ATTR_0:.6e}, OMEGA_PHI_ATTR_0 = {OMEGA_PHI_ATTR_0:.6e}")
    print(f"RS_ZSTAR = {RS_ZSTAR:.6f} Mpc, RD_ZDRAG = {RD_ZDRAG:.6f} Mpc")
    prof = compute_bao_profile_controls()
    print(f"BAO chi2 M4 (Table 1) = {prof['chi2_m4_tab1']:.6f}, (Cobaya full) = {prof['chi2_m4_full']:.6f}")
    print(f"Profiled M4 shape: P_best = {prof['P_best_m4']:.2f} +/- {prof['sigma_P_m4']:.2f}, chi2_full = {prof['chi2_m4_prof_full']:.4f}")
    print(f"Profiled 2-param: Om_MLE = {prof['om_mle_2p']:.6f} +/- {prof['sigma_om_2p']:.5f}, P_best = {prof['P_best_2p']:.2f} +/- {prof['sigma_P_2p']:.2f}, chi2_full = {prof['chi2_2p_full']:.4f}")
    print(f"Joint (Om, H0, rd) control with the assumed local H0 benchmark prior: Delta_chi2 = {prof['dchi2_joint']:.4f} (sqrt(Delta_chi2) = {prof['sigma_joint_equiv']:.4f}; not a likelihood significance)")
