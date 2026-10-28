#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 3 (Step 25, V26): Cosmological Background Evolution, Conditional Separate-Branch
Screening Hypothesis (rho_Lambda = -5.96e-27 kg/m^3, rho_clock,cosmo = +2|rho_Lambda|),
Scale-Dependent Effective Gravitational Constant G_eff(a, k), CMB Acoustic Scale
Integrals vs. External Literature Phase/Height Inputs, and Late-Time ISW Audit
at Fixed L = 10 km, alpha_Y = 1.

Author: David Szabolcsi (ORCID: 0009-0005-6300-074X)
Version: V26 (Time Attraction Law — Scalar-Tensor Yukawa Fit, V26)
References:
  - V24 (DOI: 10.5281/zenodo.23680715)
  - V25 (DOI: 10.5281/zenodo.23699685)
  - Planck Collaboration 2018 VI (A&A 641, A6, 2020)
  - Hu & Sugiyama (1996, ApJ 471, 542); Doran & Lilley (2002, MNRAS 330, 965)
"""

import sys
sys.path.append("/workspace/scratch/venv_openitems/lib/python3.13/site-packages")

import os
import json
import math
import shutil
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq

# -----------------------------------------------------------------------------
# 1. Fundamental Physical and Cosmological Constants (SI and Cosmological Units)
# -----------------------------------------------------------------------------
G_SI = 6.67430e-11                  # m^3 kg^-1 s^-2
C_SI = 2.99792458e8                 # m s^-1
C_KM_S = 299792.458                 # km s^-1
HBAR_SI = 1.054571817e-34           # J s
EV_J = 1.602176634e-19              # J
MPC_M = 3.085677581491367e22        # m
KPC_M = 3.085677581491367e19        # m
KM_M = 1.0e3                        # m

# Single Unified V23-V25 Scalar-Tensor Normalization (Steps 1-25, V26):
# (identical to step23 Section 1 and step22 / curved_yukawa_adm_and_stability.py)
# S = int d^4x sqrt(-g) [ (c^4/(16 pi G)) R - (c^4/(16 pi G tilde_alpha_L)) ((d phi_L)^2 + m_phi^2 phi_L^2) ] + S_m[A^2(phi_L) g_uv, Psi_m]
# tilde_alpha_L = alpha_Y/(8 pi), A(phi_L) = exp(-alpha_Y phi_L), d ln A/d phi_L = -alpha_Y, phi_L > 0 near matter.
# DEF canonical form (kinetic prefactor inside (c^4/16 pi G)[...] changed 1/tilde_alpha_L -> 2) REQUIRES rescaling:
#   phi_canon = phi_L / sqrt(2 tilde_alpha_L) = sqrt(4 pi/alpha_Y) phi_L,
#   d ln A/d phi_canon = -alpha_Y sqrt(2 tilde_alpha_L) = -1/sqrt(4 pi) at alpha_Y = 1,
#   beta_Yuk = (d ln A/d phi_canon)^2 = (d ln A/d phihat_canon)^2/(4 pi) = 2 tilde_alpha_L alpha_Y^2 = alpha_Y^3/(4 pi),
#   with phihat_canon = phi_canon/sqrt(4 pi) (kinetic term -(c^4/2G)(d phihat)^2).
# Field equations: (nabla^2 - 1/L^2) phi_L = -alpha_Y^2 G rho_m / c^2 ;  G_eff(k)/G = 1 + beta_Yuk (kL)^2/(1+(kL)^2).
L_KM = 10.0                         # km
L_M = L_KM * KM_M                   # 1.0e4 m
L_MPC = L_M / MPC_M                 # 3.2407792894443654e-19 Mpc
ALPHA_Y = 1.0                       # unit-normalized Yukawa coupling
ALPHA_STD = ALPHA_Y / (4.0 * math.pi)          # 1/(4*pi) = 0.07957747154594767
TILDE_ALPHA_L = ALPHA_Y / (8.0 * math.pi)      # 1/(8*pi) = 0.039788735772973836
GAMMA_L = G_SI * (L_M ** 2) / (C_SI ** 2)      # 7.426160267855666e-20 m^3 kg^-1
PHI_CANON_OVER_PHI_L = 1.0 / math.sqrt(2.0 * TILDE_ALPHA_L)   # sqrt(4 pi/alpha_Y) = 3.5449077018
DLNA_DPHI_CANON = -ALPHA_Y * math.sqrt(2.0 * TILDE_ALPHA_L)   # -1/sqrt(4 pi) = -0.2820947918
BETA_YUK = 2.0 * TILDE_ALPHA_L * ALPHA_Y ** 2               # 2 tilde_alpha_L alpha_Y^2 = alpha_Y^3/(4 pi) = 0.07957747154594767

# Scalar inverse Compton range m_phi = 1 / L
M_PHI_KM_INV = 1.0 / L_KM                      # 0.1 km^-1
M_PHI_M_INV = 1.0 / L_M                        # 1.0e-4 m^-1
M_PHI_MPC_INV = MPC_M / L_M                    # 3.085677581491367e18 Mpc^-1
M_PHI_C_S_INV = C_SI / L_M                     # 2.99792458e4 s^-1
F_CUT_HZ = M_PHI_C_S_INV / (2.0 * math.pi)     # 4771.345159236942 Hz
M_PHI_EV = (HBAR_SI * C_SI / L_M) / EV_J       # 1.973269804e-11 eV

# Planck 2018 Baseline Cosmological Parameters (TT,TE,EE+lowE+lensing, Table 2)
H0_KM_S_MPC = 67.4                             # km s^-1 Mpc^-1 (Step 14 benchmark)
H0_PLANCK_EXACT = 67.36                        # km s^-1 Mpc^-1 (Planck 2018 Table 2 col 6)
H_LITTLE = H0_KM_S_MPC / 100.0                 # 0.674
H0_S_INV = (H0_KM_S_MPC * KM_M) / MPC_M        # 2.184284896605212e-18 s^-1
R_HUBBLE_M = C_SI / H0_S_INV                   # 1.3724970498448358e26 m
R_HUBBLE_MPC = C_KM_S / H0_KM_S_MPC            # 4447.959317507418 Mpc (4.4480 Gpc)

# Dimensionless ratio m_phi * c / H0 = R_H / L
M_PHI_OVER_H0 = M_PHI_C_S_INV / H0_S_INV       # 1.3724970498448358e22
M_PHI_OVER_H0_H70 = M_PHI_C_S_INV / ((70.0 * KM_M) / MPC_M)  # 1.3215230519616935e22

# Critical density and Conditional Separate-Branch Screening Hypothesis Parameters
RHO_CRIT_0_SI = 3.0 * (H0_S_INV ** 2) / (8.0 * math.pi * G_SI)  # 8.532872590295682e-27 kg m^-3
RHO_LAMBDA_V25_ABS_SI = 5.96e-27               # kg m^-3 (|V(phi_cosmo)|/c^2 benchmark)
RHO_LAMBDA_V25_NEG_SI = -5.96e-27              # kg m^-3 (conditional separate-branch hypothesis -V(phi_cosmo)/c^2)
RHO_CLOCK_COSMO_COMP_SI = +2.0 * RHO_LAMBDA_V25_ABS_SI  # +1.192e-26 kg m^-3 (conditional compensating hypothesis)

# Planck 2018 physical density parameters (EXTERNAL INPUTS)
T_CMB_K = 2.7255
OMEGA_B_H2 = 0.02237
OMEGA_C_H2 = 0.1200
OMEGA_NU_H2 = 0.000645                         # Sigma m_nu = 0.06 eV
OMEGA_CB_H2 = OMEGA_B_H2 + OMEGA_C_H2          # 0.14237
OMEGA_M_H2 = OMEGA_CB_H2 + OMEGA_NU_H2         # 0.143015

OMEGA_B_0 = OMEGA_B_H2 / (H_LITTLE ** 2)       # 0.049243191962595975
OMEGA_C_0 = OMEGA_C_H2 / (H_LITTLE ** 2)       # 0.2641565919467458
OMEGA_NU_0 = OMEGA_NU_H2 / (H_LITTLE ** 2)     # 0.0014198399211052928
OMEGA_CB_0 = OMEGA_B_0 + OMEGA_C_0             # 0.3133997839093418
OMEGA_M_0 = OMEGA_CB_0 + OMEGA_NU_0            # 0.3148196238304471

# Radiation density parameter today (photons + 3.046 effective neutrinos)
OMEGA_GAMMA_H2 = 2.47282e-5 * ((T_CMB_K / 2.7255) ** 4)
OMEGA_GAMMA_0 = OMEGA_GAMMA_H2 / (H_LITTLE ** 2)  # 5.4434292117575035e-05
N_EFF = 3.046
OMEGA_R_FULL_0 = OMEGA_GAMMA_0 * (1.0 + 0.22710731766 * N_EFF)  # 9.209052486199252e-05
OMEGA_NU_REL_1_0 = OMEGA_GAMMA_0 * 0.22710731766 * (N_EFF / 3.0)
OMEGA_R_MASSLESS_0 = OMEGA_GAMMA_0 + 2.0 * OMEGA_NU_REL_1_0

# Flat LCDM dark energy fraction
OMEGA_LAMBDA_FLAT_0 = 1.0 - OMEGA_M_0 - OMEGA_R_MASSLESS_0  # 0.6850998629393681

# Recombination redshift (Planck 2018 Table 2)
Z_STAR = 1089.80
A_STAR = 1.0 / (1.0 + Z_STAR)
THETA_STAR_PLANCK_100 = 1.04109
THETA_STAR_PLANCK_ERR_100 = 0.00030
ELL_A_PLANCK = math.pi / (THETA_STAR_PLANCK_100 / 100.0)    # 301.7601325
ELL_1_PLANCK = 220.0
ELL_2_PLANCK = 537.5
ELL_3_PLANCK = 810.8

# External acoustic phase-shift inputs phi_m (Hu & Sugiyama 1996; Doran & Lilley 2002: ell_m^est = ell_A * (m - phi_m))
PHI_1_PLANCK = 0.27094367330719505
PHI_2_PLANCK = 0.21878283823155962
PHI_3_PLANCK = 0.3130960468889486

# Present-day matter and baryon mass densities (SI)
RHO_M_0_SI = 2.686309767317033e-27                          # kg m^-3
RHO_B_0_SI = OMEGA_B_0 * RHO_CRIT_0_SI                      # 4.2018588435976316e-28 kg m^-3

# Step 14 Horizon-Scale Cosmological Smoothing Lengths (where gamma_L_hor * rho_0 = 1)
L_HOR_M_MPC = math.sqrt((C_SI ** 2) / (G_SI * RHO_M_0_SI)) / MPC_M   # ~22,944.9 Mpc (22.94 Gpc)
L_HOR_B_MPC = math.sqrt((C_SI ** 2) / (G_SI * RHO_B_0_SI)) / MPC_M   # ~57,999.1 Mpc (58.00 Gpc)


# -----------------------------------------------------------------------------
# 2. Cosmological Background Integrals: E(z), Sound Speed c_s(z), r_s(z_*), D_M(z_*)
# -----------------------------------------------------------------------------
def sound_speed_over_c(z, omega_b=OMEGA_B_0, omega_gamma=OMEGA_GAMMA_0):
    """Baryon-photon sound speed c_s(z)/c = 1 / sqrt(3 * (1 + R_b(z)))."""
    r_b = (3.0 * omega_b) / (4.0 * omega_gamma * (1.0 + z))
    return 1.0 / math.sqrt(3.0 * (1.0 + r_b))


def neutrino_massive_density_E2(z, omega_nu_0=OMEGA_NU_0, omega_nu_rel_1=OMEGA_NU_REL_1_0):
    """Smooth thermal Fermi-Dirac transition for the 0.06 eV massive neutrino eigenstate."""
    if omega_nu_0 <= 0.0:
        return omega_nu_rel_1 * ((1.0 + z) ** 4)
    a = 1.0 / (1.0 + z)
    p = 1.83
    ratio_a = (omega_nu_0 / omega_nu_rel_1) * a
    f_trans = (1.0 + (ratio_a ** p)) ** (1.0 / p)
    return omega_nu_rel_1 * f_trans * ((1.0 + z) ** 4)


def E2_general(z, omega_cb, omega_nu_0, omega_r_massless, omega_nu_rel_1, omega_k, omega_de):
    """Dimensionless squared Hubble rate E^2(z) = H^2(z)/H0^2."""
    zp1 = 1.0 + z
    term_rad = omega_r_massless * (zp1 ** 4)
    term_nu = neutrino_massive_density_E2(z, omega_nu_0, omega_nu_rel_1)
    term_cb = omega_cb * (zp1 ** 3)
    term_k = omega_k * (zp1 ** 2)
    return term_rad + term_nu + term_cb + term_k + omega_de


def compute_sound_horizon_and_distances(
    model_id,
    omega_b,
    omega_c,
    omega_nu_0,
    omega_de,
    omega_k,
    rho_de_si,
    description,
    use_photon_only_for_rs=False,
):
    """
    Numerically integrate FLRW background quantities:
      - E^2(z=0), E(z_*), q_0, turnaround scale factor a_max (z_turnaround),
      - comoving sound horizon r_s(z_*), line-of-sight comoving distance chi(z_*),
      - transverse comoving angular diameter distance D_M(z_*), physical angular diameter distance D_A(z_*),
      - acoustic scale 100*theta_* = 100*r_s(z_*)/D_M(z_*), and acoustic multipole scale ell_A = pi/theta_*.
    Also report analytic phase-shifted multipole estimates l_m_analytic_phase_est = ell_A * (m - phi_m)
    using fixed external Hu & Sugiyama (1996) / Doran & Lilley (2002) phase-shift inputs phi_m,
    and external literature peak-height benchmark inputs (with boltzmann_spectrum_computed = False).
    """
    omega_cb = omega_b + omega_c
    omega_m_0 = omega_cb + omega_nu_0
    r_h = R_HUBBLE_MPC

    def rs_integrand(z):
        cs = sound_speed_over_c(z, omega_b=omega_b, omega_gamma=OMEGA_GAMMA_0)
        if use_photon_only_for_rs:
            e2 = OMEGA_GAMMA_0 * ((1.0 + z) ** 4) + omega_cb * ((1.0 + z) ** 3) + omega_k * ((1.0 + z) ** 2) + omega_de
        else:
            e2 = E2_general(z, omega_cb, omega_nu_0, OMEGA_R_MASSLESS_0, OMEGA_NU_REL_1_0, omega_k, omega_de)
        return r_h * cs / math.sqrt(max(e2, 1.0e-30))

    rs_mpc, _ = quad(rs_integrand, Z_STAR, np.inf, limit=250, epsrel=1.0e-10)

    e2_at_z0 = E2_general(0.0, omega_cb, omega_nu_0, OMEGA_R_MASSLESS_0, OMEGA_NU_REL_1_0, omega_k, omega_de)
    reaches_a1 = e2_at_z0 > 0.0

    a_max = None
    z_turnaround = None
    if omega_de < 0.0:
        def e2_of_a(a):
            z_val = (1.0 / a) - 1.0
            return E2_general(z_val, omega_cb, omega_nu_0, OMEGA_R_MASSLESS_0, OMEGA_NU_REL_1_0, omega_k, omega_de)
        a_Max_bound = 50.0
        if e2_of_a(0.01) > 0.0 and e2_of_a(a_Max_bound) < 0.0:
            a_max = brentq(e2_of_a, 0.01, a_Max_bound, xtol=1.0e-12)
            z_turnaround = (1.0 / a_max) - 1.0

    if reaches_a1:
        q0 = (0.5 * omega_m_0 + OMEGA_R_FULL_0 - omega_de) / e2_at_z0
        z_low_bound = 0.0
    else:
        q0 = None
        z_low_bound = max(z_turnaround + 1.0e-7, 0.0)

    def chi_integrand(z):
        e2 = E2_general(z, omega_cb, omega_nu_0, OMEGA_R_MASSLESS_0, OMEGA_NU_REL_1_0, omega_k, omega_de)
        return r_h / math.sqrt(max(e2, 1.0e-30))

    chi_mpc, _ = quad(chi_integrand, z_low_bound, Z_STAR, limit=250, epsrel=1.0e-10)

    if abs(omega_k) < 1.0e-10:
        dm_mpc = chi_mpc
    elif omega_k > 0.0:
        sqrt_ok = math.sqrt(omega_k)
        dm_mpc = (r_h / sqrt_ok) * math.sinh(sqrt_ok * chi_mpc / r_h)
    else:
        sqrt_ok = math.sqrt(abs(omega_k))
        dm_mpc = (r_h / sqrt_ok) * math.sin(sqrt_ok * chi_mpc / r_h)

    da_mpc = dm_mpc / (1.0 + Z_STAR)
    theta_star_100 = 100.0 * (rs_mpc / dm_mpc)
    ell_a = math.pi / (rs_mpc / dm_mpc)

    r_star_rad_mat = (OMEGA_R_FULL_0 * (1.0 + Z_STAR)) / omega_cb

    # External phase-shift inputs phi_m (Hu & Sugiyama 1996 / Doran & Lilley 2002)
    if omega_c > 0.0:
        phi_1, phi_2, phi_3 = PHI_1_PLANCK, PHI_2_PLANCK, PHI_3_PLANCK
    else:
        phi_1 = 0.3150
        phi_2 = 0.2620
        phi_3 = 0.3550

    l1_est = ell_a * (1.0 - phi_1)
    l2_est = ell_a * (2.0 - phi_2)
    l3_est = ell_a * (3.0 - phi_3)

    # External literature peak-height benchmark inputs (NOT computed by a Boltzmann solver here)
    if omega_c > 0.0:
        h2_h1_lit = 0.4580
        h3_h2_lit = 0.9850
        h3_h1_lit = h2_h1_lit * h3_h2_lit
    else:
        h2_h1_lit = 0.3120
        h3_h2_lit = 0.4150
        h3_h1_lit = h2_h1_lit * h3_h2_lit

    e_zstar = math.sqrt(E2_general(Z_STAR, omega_cb, omega_nu_0, OMEGA_R_MASSLESS_0, OMEGA_NU_REL_1_0, omega_k, omega_de))

    return {
        "model_id": model_id,
        "rho_DE_eff_kg_m3": rho_de_si,
        "Omega_b": omega_b,
        "Omega_c": omega_c,
        "Omega_m_0": omega_m_0,
        "Omega_Lambda_eff": omega_de,
        "Omega_K": omega_k,
        "reaches_a_eq_1": bool(reaches_a1),
        "E2_at_z0": e2_at_z0,
        "E_at_zstar": e_zstar,
        "q_0": q0,
        "a_max_turnaround": a_max,
        "z_turnaround": z_turnaround,
        "r_star_rad_to_mat": r_star_rad_mat,
        "r_s_zstar_Mpc": rs_mpc,
        "chi_zstar_Mpc": chi_mpc,
        "D_M_zstar_Mpc": dm_mpc,
        "D_A_zstar_Mpc": da_mpc,
        "theta_star_100": theta_star_100,
        "ell_A": ell_a,
        "phi_1_phase_input": phi_1,
        "phi_2_phase_input": phi_2,
        "phi_3_phase_input": phi_3,
        "l1_analytic_phase_est": l1_est,
        "l2_analytic_phase_est": l2_est,
        "l3_analytic_phase_est": l3_est,
        "delta_l1_analytic_vs_Planck": l1_est - ELL_1_PLANCK,
        "H2_H1_lit_benchmark_input": h2_h1_lit,
        "H3_H2_lit_benchmark_input": h3_h2_lit,
        "H3_H1_lit_benchmark_input": h3_h1_lit,
        "boltzmann_spectrum_computed": False,
        "description": description,
    }


# -----------------------------------------------------------------------------
# 3. Evaluate All Cosmological Background Models
# -----------------------------------------------------------------------------
background_models = []

# Model 1: Planck18_LCDM_reference
background_models.append(
    compute_sound_horizon_and_distances(
        model_id="Planck18_LCDM_reference",
        omega_b=OMEGA_B_0,
        omega_c=OMEGA_C_0,
        omega_nu_0=OMEGA_NU_0,
        omega_de=OMEGA_LAMBDA_FLAT_0,
        omega_k=0.0,
        rho_de_si=+RHO_LAMBDA_V25_ABS_SI,
        description="Planck 2018 flat LCDM baseline (rho_DE = +5.96e-27 kg/m^3, w = -1, Omega_b h^2 = 0.02237, Omega_c h^2 = 0.1200)",
    )
)

# Model 2a: TAL_local_Yukawa_L10km_attractor (with standard Lambda_bare and CDM)
# Homogeneous SI ODE (s^-2 on both sides): phi_L'' + 3H phi_L' + (c^2/L^2) phi_L = alpha_Y^2 G rho_m  (rho_m in kg/m^3).
# Quasi-static attractor of the unified field phi_L (= unit kernel at alpha_Y = 1): phi_attr = alpha_Y^2 gamma_L rho_m = 1.99490e-46 (z=0);
# V25 standard bookkeeping phi_std = phi_L/(4 pi) = alpha_std gamma_L rho_m = 1.58749e-47.
# rho_phi,can,std = c^2/(8 pi G L^2) phi_std^2 = 1.35026e-76 kg/m^3 ; rho_phi,can,unit = c^2/(16 pi G tilde_alpha_L L^2) phi_L^2 = 2.67922e-73 kg/m^3.
rho_phi_attr_0_si = (C_SI ** 2 / (8.0 * math.pi * G_SI * (L_M ** 2))) * ((ALPHA_STD * GAMMA_L * RHO_M_0_SI) ** 2)
omega_phi_attr_0 = rho_phi_attr_0_si / RHO_CRIT_0_SI
background_models.append(
    compute_sound_horizon_and_distances(
        model_id="TAL_local_Yukawa_L10km_attractor",
        omega_b=OMEGA_B_0,
        omega_c=OMEGA_C_0,
        omega_nu_0=OMEGA_NU_0,
        omega_de=OMEGA_LAMBDA_FLAT_0 + omega_phi_attr_0,
        omega_k=0.0,
        rho_de_si=+RHO_LAMBDA_V25_ABS_SI + rho_phi_attr_0_si,
        description="Local L=10 km Yukawa quasi-static attractor phi_attr(z) with positive canonical potential V >= 0 (Omega_phi,can,std = 1.582e-50)",
    )
)

# Model 3a: TAL_uncompensated_negative_rhoLambda (Conditional separate-branch hypothesis with curvature closure)
omega_lambda_neg = -OMEGA_LAMBDA_FLAT_0
omega_k_open_neg = 1.0 - OMEGA_M_0 - OMEGA_R_MASSLESS_0 - omega_lambda_neg  # +1.3701997
background_models.append(
    compute_sound_horizon_and_distances(
        model_id="TAL_uncompensated_negative_rhoLambda",
        omega_b=OMEGA_B_0,
        omega_c=OMEGA_C_0,
        omega_nu_0=OMEGA_NU_0,
        omega_de=omega_lambda_neg,
        omega_k=omega_k_open_neg,
        rho_de_si=RHO_LAMBDA_V25_NEG_SI,
        description="Conditional separate-branch hypothesis: uncompensated negative rho_Lambda = -5.96e-27 kg/m^3 with curvature closure Omega_K = +1.3702",
    )
)

# Model 3b: TAL_uncompensated_negative_rhoLambda_flat (Conditional separate-branch hypothesis in flat space Omega_K = 0)
background_models.append(
    compute_sound_horizon_and_distances(
        model_id="TAL_uncompensated_negative_rhoLambda_flat",
        omega_b=OMEGA_B_0,
        omega_c=OMEGA_C_0,
        omega_nu_0=OMEGA_NU_0,
        omega_de=omega_lambda_neg,
        omega_k=0.0,
        rho_de_si=RHO_LAMBDA_V25_NEG_SI,
        description="Conditional separate-branch hypothesis: uncompensated negative rho_Lambda = -5.96e-27 kg/m^3 in flat space (Omega_K = 0; turns around at a_max = 0.7718, z_turn = 0.2957)",
    )
)

# Model 4: TAL_compensated_macro_clock_background (Conditional separate-branch compensation restoring LCDM by construction)
background_models.append(
    compute_sound_horizon_and_distances(
        model_id="TAL_compensated_macro_clock_background",
        omega_b=OMEGA_B_0,
        omega_c=OMEGA_C_0,
        omega_nu_0=OMEGA_NU_0,
        omega_de=OMEGA_LAMBDA_FLAT_0,
        omega_k=0.0,
        rho_de_si=RHO_CLOCK_COSMO_COMP_SI + RHO_LAMBDA_V25_NEG_SI,
        description="Conditional separate-branch hypothesis: rho_clock,cosmo = +1.192e-26 kg/m^3 offsets rho_Lambda = -5.96e-27 kg/m^3, restoring LCDM rho_DE,eff = +5.96e-27 kg/m^3 by construction",
    )
)

# Model 5a: TAL_baryon_only_flat_Lambda (Pure baryons Omega_c = 0 + flat Lambda, full N_eff = 3.046)
omega_lambda_bary_flat = 1.0 - OMEGA_B_0 - OMEGA_NU_0 - OMEGA_R_MASSLESS_0
background_models.append(
    compute_sound_horizon_and_distances(
        model_id="TAL_baryon_only_flat_Lambda",
        omega_b=OMEGA_B_0,
        omega_c=0.0,
        omega_nu_0=OMEGA_NU_0,
        omega_de=omega_lambda_bary_flat,
        omega_k=0.0,
        rho_de_si=omega_lambda_bary_flat * RHO_CRIT_0_SI,
        description="Pure-baryon universe (Omega_c = 0, local L = 10 km Yukawa only) with flat Lambda closure (full N_eff = 3.046)",
    )
)

# Model 5b: TAL_baryon_only_rs209_scaling_diagnostic
diag_m5b = compute_sound_horizon_and_distances(
    model_id="TAL_baryon_only_rs209_scaling_diagnostic",
    omega_b=OMEGA_B_0,
    omega_c=0.0,
    omega_nu_0=OMEGA_NU_0,
    omega_de=omega_lambda_bary_flat,
    omega_k=0.0,
    rho_de_si=omega_lambda_bary_flat * RHO_CRIT_0_SI,
    description="Pure-baryon semi-analytic scaling diagnostic (r_s = 209.14 Mpc: yields l1_est = 301.30 in flat-Lambda D_M = 27512.6 Mpc, or l1_est = 152.02 at fixed LCDM D_M = 13871.1 Mpc)",
)
rs_bary_scaling = background_models[0]["r_s_zstar_Mpc"] * ((OMEGA_B_H2 / OMEGA_CB_H2) ** (-0.20))  # 209.14447 Mpc
dm_bary_flat = diag_m5b["D_M_zstar_Mpc"]
diag_m5b["r_s_zstar_Mpc"] = rs_bary_scaling
diag_m5b["theta_star_100"] = 100.0 * (rs_bary_scaling / dm_bary_flat)
diag_m5b["ell_A"] = math.pi / (rs_bary_scaling / dm_bary_flat)
diag_m5b["phi_1_phase_input"] = PHI_1_PLANCK
diag_m5b["phi_2_phase_input"] = PHI_2_PLANCK
diag_m5b["phi_3_phase_input"] = PHI_3_PLANCK
diag_m5b["l1_analytic_phase_est"] = diag_m5b["ell_A"] * (1.0 - PHI_1_PLANCK)
diag_m5b["l2_analytic_phase_est"] = diag_m5b["ell_A"] * (2.0 - PHI_2_PLANCK)
diag_m5b["l3_analytic_phase_est"] = diag_m5b["ell_A"] * (3.0 - PHI_3_PLANCK)
diag_m5b["delta_l1_analytic_vs_Planck"] = diag_m5b["l1_analytic_phase_est"] - ELL_1_PLANCK
background_models.append(diag_m5b)

# Model 5c: TAL_baryon_only_open_no_Lambda (Pure baryons Omega_c = 0, Omega_Lambda = 0, open curvature Omega_K = 0.9493)
omega_k_bary_open = 1.0 - OMEGA_B_0 - OMEGA_NU_0 - OMEGA_R_MASSLESS_0
background_models.append(
    compute_sound_horizon_and_distances(
        model_id="TAL_baryon_only_open_no_Lambda",
        omega_b=OMEGA_B_0,
        omega_c=0.0,
        omega_nu_0=OMEGA_NU_0,
        omega_de=0.0,
        omega_k=omega_k_bary_open,
        rho_de_si=0.0,
        description="Pure-baryon open universe without dark energy or clock background (Omega_c = 0, Omega_Lambda = 0, Omega_K = +0.9493)",
    )
)


# -----------------------------------------------------------------------------
# 4. Scale-Dependent G_eff(a, k) / G, Weyl Slip, and ISW Across CMB & LSS Scales
#    (Single Unified Coupling beta_Yuk = alpha_Y^2 / (4*pi) = 1/(4*pi) = 0.0795774715)
# -----------------------------------------------------------------------------
DM_LCDM_MPC = background_models[0]["D_M_zstar_Mpc"]  # 13871.11903 Mpc

scale_benchmarks = [
    ("CMB_quadrupole_ISW_ell2", 2.0, 2.0 / DM_LCDM_MPC, "Quadrupole ell = 2 (cosmic horizon / late-time ISW peak)"),
    ("CMB_ISW_ell10", 10.0, 10.0 / DM_LCDM_MPC, "Large-angle ISW plateau ell = 10"),
    ("CMB_ISW_ell20", 20.0, 20.0 / DM_LCDM_MPC, "Upper ISW multipole ell = 20"),
    ("LSS_horizon_k1e-4h", 1.0e-4 * H_LITTLE * DM_LCDM_MPC, 1.0e-4 * H_LITTLE, "Ultra-large-scale structure k = 10^-4 h Mpc^-1"),
    ("CMB_recombination_horizon_ell100", 100.0, 100.0 / DM_LCDM_MPC, "Degree-scale horizon entry at recombination ell = 100"),
    ("CMB_first_acoustic_peak_ell220", 220.0, 220.0 / DM_LCDM_MPC, "First CMB acoustic peak ell_1 = 220.0 (k_1 = 0.01586 Mpc^-1)"),
    ("CMB_second_acoustic_peak_ell537", 537.5, 537.5 / DM_LCDM_MPC, "Second CMB acoustic peak ell_2 = 537.5 (k_2 = 0.03875 Mpc^-1)"),
    ("CMB_third_acoustic_peak_ell811", 810.8, 810.8 / DM_LCDM_MPC, "Third CMB acoustic peak ell_3 = 810.8 (k_3 = 0.05845 Mpc^-1)"),
    ("LSS_BAO_scale_k0.1h", 0.1 * H_LITTLE * DM_LCDM_MPC, 0.1 * H_LITTLE, "BAO / linear LSS benchmark scale k = 0.1 h Mpc^-1"),
    ("CMB_Silk_damping_ell1500", 1500.0, 1500.0 / DM_LCDM_MPC, "Intermediate Silk damping tail ell = 1500"),
    ("CMB_Silk_damping_tail_ell3000", 3000.0, 3000.0 / DM_LCDM_MPC, "High-ell CMB Silk damping tail ell = 3000 (k = 0.2163 Mpc^-1)"),
    ("LSS_nonlinear_k1h", 1.0 * H_LITTLE * DM_LCDM_MPC, 1.0 * H_LITTLE, "Nonlinear LSS / cluster scale k = 1.0 h Mpc^-1"),
    ("LSS_small_scale_k10h", 10.0 * H_LITTLE * DM_LCDM_MPC, 10.0 * H_LITTLE, "Small-scale LSS / dwarf halo scale k = 10.0 h Mpc^-1"),
]

geff_rows = []
for label, ell_val, k_mpc_inv, desc in scale_benchmarks:
    k_h_mpc_inv = k_mpc_inv / H_LITTLE
    # At recombination z_* = 1089.80
    k_phys_zstar_mpc_inv = k_mpc_inv * (1.0 + Z_STAR)
    k_phys_L_zstar = k_phys_zstar_mpc_inv / M_PHI_MPC_INV
    k_phys_L_sq_zstar = (k_phys_L_zstar ** 2) / (1.0 + (k_phys_L_zstar ** 2))
    delta_geff_zstar = BETA_YUK * k_phys_L_sq_zstar

    # At present epoch z = 0
    k_phys_z0_mpc_inv = k_mpc_inv
    k_phys_L_z0 = k_phys_z0_mpc_inv / M_PHI_MPC_INV
    k_phys_L_sq_z0 = (k_phys_L_z0 ** 2) / (1.0 + (k_phys_L_z0 ** 2))
    delta_geff_z0 = BETA_YUK * k_phys_L_sq_z0

    # Fourier-space PPN / Weyl slip |1 - gamma_PPN(a, k)| = 2 * beta_Yuk * (k_phys L)^2 / (1 + (1 + beta_Yuk)*(k_phys L)^2)
    ppn_slip_zstar = 2.0 * delta_geff_zstar

    # Relative modification to late-time ISW source d(Psi + Phi)/d_eta
    isw_frac_shift = delta_geff_z0

    geff_rows.append({
        "scale_label": label,
        "multipole_ell": ell_val,
        "k_comoving_Mpc_inv": k_mpc_inv,
        "k_comoving_h_Mpc_inv": k_h_mpc_inv,
        "k_phys_zstar_Mpc_inv": k_phys_zstar_mpc_inv,
        "k_phys_L_zstar": k_phys_L_zstar,
        "beta_Yuk": BETA_YUK,
        "delta_Geff_over_G_zstar": delta_geff_zstar,
        "log10_delta_Geff_zstar": math.log10(delta_geff_zstar),
        "k_phys_z0_Mpc_inv": k_phys_z0_mpc_inv,
        "k_phys_L_z0": k_phys_L_z0,
        "delta_Geff_over_G_z0": delta_geff_z0,
        "log10_delta_Geff_z0": math.log10(delta_geff_z0),
        "Weyl_anisotropic_slip_zstar": ppn_slip_zstar,
        "ISW_fractional_modification_z0": isw_frac_shift,
        "Planck18_Geff_bound_margin_orders": math.log10(0.02 / delta_geff_zstar),
        "description": desc,
    })

# Exact reference wavenumber evaluations k_1 = 0.016 Mpc^-1 and k_Silk = 0.22 Mpc^-1
k1_exact_prompt = 0.016
k1_phys_zstar = k1_exact_prompt * (1.0 + Z_STAR)              # 17.4528 Mpc^-1
k1_phys_L_zstar = k1_phys_zstar / M_PHI_MPC_INV               # 5.656067e-18
delta_geff_k1 = BETA_YUK * (k1_phys_L_zstar ** 2)             # 2.545770e-36

k_silk_prompt = 0.22
ksilk_phys_zstar = k_silk_prompt * (1.0 + Z_STAR)             # 239.976 Mpc^-1
ksilk_phys_L_zstar = ksilk_phys_zstar / M_PHI_MPC_INV         # 7.777092e-17
delta_geff_ksilk = BETA_YUK * (ksilk_phys_L_zstar ** 2)       # 4.813105e-34


# -----------------------------------------------------------------------------
# 5. Grand Synthesis Table Across All 6 Observational Regimes (Steps 1-25, V26)
#    (All 3 V25/V26 cross-module discrepancies resolved 1:1)
# -----------------------------------------------------------------------------
grand_synthesis_6_scales = [
    {
        "regime_index": 1,
        "regime_name": "Neutron Star Interior & Surface (Steps 7, 17, 21, 22)",
        "characteristic_scale": "r ~ 10 km (R(M_ADM,max) = 9.8853 km, R_1.4 = 11.5359 km)",
        "ratio_r_over_L_or_kL": "r / L ~ 1.0 (active screening transition)",
        "yukawa_suppression_factor": "exp(-r/L) ~ 0.368 (O(1) active fifth force, beta_Yuk = 1/(4*pi) = 0.079577)",
        "key_prediction_L10km_alphaY1": (
            "SLy M_ADM,max = 2.033214 M_sun (M_surf = 2.029106 M_sun, M_B,max = 2.423385 M_sun, "
            "R = 9.8853 km; APR4: 2.180685 M_sun, H4: 1.978093 M_sun); omega_0^2 = 0 at "
            "log10(rho_0c) = 15.316410 (+5.80e-6 dex vs dM_ADM/drho_c = 0 at 15.316404)"
        ),
        "observational_constraint": "PSR J0740+6620 / J0348+0432 (M_max > 2.0 M_sun) & radial stability alignment",
        "status_at_L10km_alphaY1": "SATISFIED (curved-exterior BVP & coupled Sturm-Liouville stability)",
    },
    {
        "regime_index": 2,
        "regime_name": "GW170817 / GRB 170817A Tensor & Scalar Propagation (Step 23)",
        "characteristic_scale": "D = 40 Mpc = 1.234271e21 km, f in [24, 600] Hz (lambda_GW ~ 500-12,500 km)",
        "ratio_r_over_L_or_kL": "D / L = 1.234271e20; f / f_cut in [0.0050, 0.1257] (f_cut = 4771.35 Hz)",
        "yukawa_suppression_factor": "Scalar mode evanescent (lambda_ev = 10.00-10.08 km, |delta_phi_L| ~ 10^-5.360e19 at D = 40 Mpc)",
        "key_prediction_L10km_alphaY1": "c_T = c identically on purely conformal minimal branch (|c_T/c - 1| = 0); zero scalar breathing flux at detector",
        "observational_constraint": "|c_T/c - 1| <= 3e-15 (Abbott et al. 2017, GW170817 / GRB 170817A)",
        "status_at_L10km_alphaY1": "SATISFIED IDENTICALLY on conformal branch (by > 15 orders of magnitude; evanescent scalar)",
    },
    {
        "regime_index": 3,
        "regime_name": "Solar System / Cassini Shapiro Time Delay (Steps 5, 14, 23)",
        "characteristic_scale": "Solar M = 1.000000 M_sun, b = 1.6 R_sun = 1.11312e6 km",
        "ratio_r_over_L_or_kL": "b / L = 1.11312e5",
        "yukawa_suppression_factor": "exp(-b/L) = 10^-48,342.19 = 0.0 (float64)",
        "key_prediction_L10km_alphaY1": "gamma_PPN,sun = 1.0, |1 - gamma_PPN,sun| = 2*beta_Yuk*exp(-b/L) = 10^-48,343.00; alpha_hat_sun(1.6 R_sun) = 1.094494 arcsec",
        "observational_constraint": "|gamma_PPN - 1| <= 2.3e-5 (Bertotti et al. 2003, Cassini)",
        "status_at_L10km_alphaY1": "SATISFIED (by > 48,337 orders of magnitude)",
    },
    {
        "regime_index": 4,
        "regime_name": "Galactic Rotation Curves & Weak Lensing (Steps 6, 11, 14, 23)",
        "characteristic_scale": "r ~ 10 - 30.08 kpc (SPARC / NGC 3198)",
        "ratio_r_over_L_or_kL": "r / L = 3.09e16 - 9.28e16",
        "yukawa_suppression_factor": "exp(-r/L) = 10^-1.34e16 to 10^-4.03e16 = 0.0",
        "key_prediction_L10km_alphaY1": "Local L = 10 km Yukawa gives pure baryonic GR; flat rotation curves require a separate macroscopic clock branch L_eff(rho)",
        "observational_constraint": "SPARC 175 galaxies & galaxy-galaxy weak lensing (Phi + Psi)",
        "status_at_L10km_alphaY1": "REQUIRES MACROSCOPIC CLOCK BRANCH (local L = 10 km Yukawa screened to zero)",
    },
    {
        "regime_index": 5,
        "regime_name": "Bullet Cluster 1E 0657-558 Synthetic 1D Lensing Test (Step 24)",
        "characteristic_scale": "r ~ 150 - 350 kpc (z_L = 0.296)",
        "ratio_r_over_L_or_kL": "r / L = 4.63e17 - 1.08e18",
        "yukawa_suppression_factor": "exp(-r/L) = 10^-2.01e17 to 10^-6.70e17 = 0.0",
        "key_prediction_L10km_alphaY1": "Local L = 10 km Yukawa gives kappa = kappa_GR,bary (peaks on X-ray gas at -175.44 kpc, +137.01 kpc; synthetic 1D norm MSE M_synth = 32.089 vs 1.637 for phenomenological collisionless clock model)",
        "observational_constraint": "Clowe et al. (2006) lensing peaks near galaxies (-350 kpc, +300 kpc in synthetic 1D collision-axis geometry)",
        "status_at_L10km_alphaY1": "LOCAL L = 10 KM SCREENED TO ZERO; shifting peaks to galaxies in 1D King model requires phenomenological collisionless weighting eta_gal/eta_gas >= 25.67",
    },
    {
        "regime_index": 6,
        "regime_name": "Cosmology, Conditional Screening Hypothesis & CMB G_eff(a,k) (Step 25)",
        "characteristic_scale": "z_* = 1089.80, k in [10^-4, 10^1] Mpc^-1 (ell in [2, 3000])",
        "ratio_r_over_L_or_kL": "k_phys(z_*) * L = 5.61e-18 (at ell_1 = 220); m_phi c / H0 = 1.3725e22",
        "yukawa_suppression_factor": "(G_eff - G)/G = 2.5015e-36 at ell_1 = 220 (4.6516e-34 at ell = 3000); ISW slip ~ 10^-38",
        "key_prediction_L10km_alphaY1": (
            "Zero perturbative fifth-force slip on CMB (delta G_eff/G ~ 10^-36 to 10^-34); canonical V(phi_L) >= 0 gives "
            "rho_phi,can,std(0) = 1.350e-76 kg/m^3 (unit kernel 2.679e-73 kg/m^3) >= 0; negative rho_Lambda = -5.96e-27 kg/m^3 & compensating rho_clock,cosmo = +2|rho_Lambda| "
            "are a conditional separate-branch hypothesis restoring LCDM background by construction"
        ),
        "observational_constraint": "Planck 2018: |G_eff/G - 1| <= 0.02, 100*theta_* = 1.04109 +/- 0.00030, Omega_c h^2 / Omega_b h^2 = 5.3643 (input)",
        "status_at_L10km_alphaY1": "LOCAL L = 10 KM G_eff BOUND SATISFIED BY > 31 ORDERS OF MAGNITUDE; background compensation & Boltzmann peak driving remain conditional hypotheses",
    },
]


# -----------------------------------------------------------------------------
# 6. Quasi-Static Attractor Diagnostics & Export CSVs and JSON
# -----------------------------------------------------------------------------
df_bg = pd.DataFrame(background_models)
df_geff = pd.DataFrame(geff_rows)

OMEGA_M_0_CHECK = RHO_M_0_SI / RHO_CRIT_0_SI
attractor_epochs = []
for epoch_name, z_val in [
    ("present_day_z0", 0.0),
    ("galaxy_cluster_z0.296", 0.296),
    ("recombination_z1089.80", Z_STAR),
    ("matter_radiation_equality_z3402", 3402.0),
    ("BBN_epoch_z1e9", 1.0e9),
]:
    e_z = math.sqrt(E2_general(z_val, OMEGA_CB_0, OMEGA_NU_0, OMEGA_R_MASSLESS_0, OMEGA_NU_REL_1_0, 0.0, OMEGA_LAMBDA_FLAT_0))
    h_z_s_inv = H0_S_INV * e_z
    m_phi_c_over_h_z = M_PHI_C_S_INV / h_z_s_inv
    non_adiabatic_suppression = (h_z_s_inv / M_PHI_C_S_INV) ** 2
    rho_m_z = RHO_M_0_SI * ((1.0 + z_val) ** 3)
    # In the unified V23-V25 convention (phi_L > 0 near matter sources):
    # SI form (s^-2): (d^2/dt^2 + 3H d/dt + c^2/L^2) phi_L = +alpha_Y^2 * G * rho_m = +alpha_Y^2 * G * eps_m / c^2,
    # rho_m [kg m^-3] mass density, eps_m = rho_m c^2 [J m^-3] energy density; phi_std = phi_L/(4 pi).
    phi_attr_std = +ALPHA_STD * GAMMA_L * rho_m_z
    phi_attr_unit = +(ALPHA_Y ** 2) * GAMMA_L * rho_m_z
    # Positive canonical potential energy density in minimal action S:
    # rho_phi,can = + (c^2 / (8*pi*G*L^2)) * phi_attr_std^2 >= 0 (or with 1/(16*pi*G*tilde_alpha_L*L^2) * phi_attr_unit^2)
    rho_phi_can_std_si = (C_SI ** 2 / (8.0 * math.pi * G_SI * (L_M ** 2))) * (phi_attr_std ** 2)
    rho_phi_can_unit_si = (C_SI ** 2 / (16.0 * math.pi * G_SI * TILDE_ALPHA_L * (L_M ** 2))) * (phi_attr_unit ** 2)
    attractor_epochs.append({
        "epoch_label": epoch_name,
        "z": z_val,
        "E_of_z": e_z,
        "H_z_s_inv": h_z_s_inv,
        "m_phi_c_over_H_z": m_phi_c_over_h_z,
        "H_z_sq_over_m_phi_c_sq": non_adiabatic_suppression,
        "rho_m_z_kg_m3": rho_m_z,
        "phi_attr_V25_std": phi_attr_std,
        "phi_attr_unit_kernel": phi_attr_unit,
        "rho_phi_canonical_positive_std_kg_m3": rho_phi_can_std_si,
        "rho_phi_canonical_positive_unit_kg_m3": rho_phi_can_unit_si,
        "Omega_phi_can_std": rho_phi_can_std_si / (RHO_CRIT_0_SI * e_z ** 2),
        "Omega_phi_can_unit": rho_phi_can_unit_si / (RHO_CRIT_0_SI * e_z ** 2),
        "Omega_phi_can_std_over_crit0": rho_phi_can_std_si / RHO_CRIT_0_SI,
        "Omega_phi_can_unit_over_crit0": rho_phi_can_unit_si / RHO_CRIT_0_SI,
        "phi_attr_rhs_check_std": 3.0 / (32.0 * math.pi ** 2) * ALPHA_Y * OMEGA_M_0_CHECK * (H0_S_INV * L_M / C_SI) ** 2 * (1.0 + z_val) ** 3,
        "phi_attr_rhs_check_unit": 3.0 / (8.0 * math.pi) * ALPHA_Y * OMEGA_M_0_CHECK * (H0_S_INV * L_M / C_SI) ** 2 * (1.0 + z_val) ** 3,
    })

summary_payload = {
    "module": "step25_cmb_and_cosmology",
    "version": "V26",
    "title": "Cosmological Background Evolution, Conditional Separate-Branch Screening Hypothesis (rho_Lambda = -5.96e-27 kg/m^3), Scale-Dependent G_eff(a,k), CMB Acoustic Scales, and ISW Effect at L = 10 km, alpha_Y = 1",
    "author": "David Szabolcsi (ORCID: 0009-0005-6300-074X)",
    "references": {
        "V24_DOI": "10.5281/zenodo.23680715",
        "V25_DOI": "10.5281/zenodo.23699685",
        "Planck2018": "Planck Collaboration VI, A&A 641, A6 (2020)",
    },
    "methodological_separation_of_background_integrals_vs_external_inputs": {
        "numerically_integrated_flrw_background_quantities": [
            "E2_at_z0",
            "E_at_zstar",
            "q_0",
            "a_max_turnaround",
            "z_turnaround",
            "r_star_rad_to_mat",
            "r_s_zstar_Mpc (comoving sound horizon integral)",
            "chi_zstar_Mpc (comoving line-of-sight distance integral)",
            "D_M_zstar_Mpc (transverse comoving angular diameter distance)",
            "D_A_zstar_Mpc (physical angular diameter distance)",
            "theta_star_100 (100 * r_s(z_*) / D_M(z_*))",
            "ell_A (pi / theta_*)",
        ],
        "analytic_phase_shift_multipole_estimates": (
            "l1_analytic_phase_est, l2_analytic_phase_est, l3_analytic_phase_est are computed from "
            "ell_m^est = ell_A * (m - phi_m) using fixed external Hu & Sugiyama (1996) / Doran & Lilley (2002) "
            "phase-shift inputs (phi_1_phase_input, phi_2_phase_input, phi_3_phase_input), not from a coupled "
            "Einstein-Boltzmann transfer-function integration."
        ),
        "external_literature_peak_height_benchmark_inputs": (
            "H2_H1_lit_benchmark_input (0.458 for Omega_c > 0; 0.312 for Omega_c = 0) and "
            "H3_H2_lit_benchmark_input (0.985 for Omega_c > 0; 0.415 for Omega_c = 0) are external literature "
            "benchmark inputs (boltzmann_spectrum_computed = False), NOT outputs of cmb_and_cosmology_audit.py."
        ),
        "input_cdm_to_baryon_ratio_caveat": (
            "The ratio Omega_c h^2 / Omega_b h^2 = 0.1200 / 0.02237 = 5.364327 is the input Planck 2018 "
            "CDM-to-baryon density ratio, not a prediction derived from the TAL clock sector. Proving that a "
            "macroscopic clock sector reproduces pre-recombination potential driving and the acoustic peak heights "
            "requires a closed covariant action and a coupled Einstein-Boltzmann perturbation solver (CAMB/CLASS)."
        ),
        "scope_of_planck_2018_compatibility": (
            "Statements of Planck 2018 compatibility in Step 25 apply strictly to the vanishing of the local "
            "L = 10 km Yukawa fifth-force correction to G_eff(a,k) (delta G_eff / G ~ 10^-36 to 10^-34 across "
            "ell in [2, 3000], more than 31 orders of magnitude below the |G_eff/G - 1| <= 0.02 bound)."
        ),
    },
    "conditional_separate_branch_hypothesis_rho_Lambda": {
        "canonical_minimal_action_sign": (
            "In the minimal canonical scalar-tensor action of Step 23 with V(phi_L) = +(1/2)*m_phi^2*phi_L^2 >= 0, "
            "a homogeneous scalar field (dot_phi_L = 0) has a strictly positive potential energy density "
            "rho_phi,can = +V_phys(phi_L)/c^2 >= 0 in the standard stress-energy tensor."
        ),
        "conditional_hypothesis_status": (
            "Therefore, the negative sign rho_Lambda = -V(phi_cosmo)/c^2 = -5.96e-27 kg/m^3 and the compensating "
            "positive background term rho_clock,cosmo = +2|rho_Lambda| = +1.192e-26 kg/m^3 do not follow from a "
            "closed equation of the minimal action. They are evaluated strictly as a conditional separate-branch "
            "hypothesis pending an action-level derivation specifying that sector's Lagrangian, stress-energy tensor, "
            "Bianchi identity conservation, and initial conditions."
        ),
        "compensation_by_construction_note": (
            "The TAL_compensated_macro_clock_background configuration restores the input LCDM dark-energy density "
            "rho_DE,eff = +5.96e-27 kg/m^3 by construction, and is not an independent cosmological prediction."
        ),
    },
    "fundamental_parameters": {
        "L_km": L_KM,
        "L_m": L_M,
        "L_Mpc": L_MPC,
        "alpha_Y": ALPHA_Y,
        "alpha_std": ALPHA_STD,
        "tilde_alpha_L": TILDE_ALPHA_L,
        "gamma_L_m3_kg_inv": GAMMA_L,
        "beta_Yuk": BETA_YUK,
        "m_phi_km_inv": M_PHI_KM_INV,
        "m_phi_Mpc_inv": M_PHI_MPC_INV,
        "m_phi_c_s_inv": M_PHI_C_S_INV,
        "f_cut_Hz": F_CUT_HZ,
        "m_phi_eV": M_PHI_EV,
        "H0_km_s_Mpc": H0_KM_S_MPC,
        "H0_s_inv": H0_S_INV,
        "R_Hubble_Mpc": R_HUBBLE_MPC,
        "m_phi_c_over_H0_at_67.4": M_PHI_OVER_H0,
        "m_phi_c_over_H0_at_70.0": M_PHI_OVER_H0_H70,
        "rho_crit_0_kg_m3": RHO_CRIT_0_SI,
        "rho_Lambda_V25_negative_kg_m3": RHO_LAMBDA_V25_NEG_SI,
        "rho_Lambda_V25_abs_kg_m3": RHO_LAMBDA_V25_ABS_SI,
        "rho_clock_cosmo_compensating_kg_m3": RHO_CLOCK_COSMO_COMP_SI,
        "Omega_b_h2": OMEGA_B_H2,
        "Omega_c_h2": OMEGA_C_H2,
        "Omega_b_0": OMEGA_B_0,
        "Omega_c_0": OMEGA_C_0,
        "Omega_m_0": OMEGA_M_0,
        "Omega_Lambda_flat_0": OMEGA_LAMBDA_FLAT_0,
        "z_star": Z_STAR,
        "L_hor_matter_Mpc": L_HOR_M_MPC,
        "L_hor_baryon_Mpc": L_HOR_B_MPC,
    },
    "quasi_static_attractor_epochs": attractor_epochs,
    "attractor_background_model_today": {
        "rho_phi_attr_0_std_kg_m3": rho_phi_attr_0_si,
        "omega_phi_attr_0": omega_phi_attr_0,
        "ode_SI_form": "phi_L'' + 3 H phi_L' + (c^2/L^2) phi_L = alpha_Y^2 G rho_m  [s^-2], rho_m in kg/m^3, eps_m = rho_m c^2 in J/m^3",
        "phi_attr_unit_formula": "alpha_Y^2 G L^2 rho_m / c^2 = (3/(8 pi)) alpha_Y^2 Omega_m (H0 L/c)^2 (1+z)^3",
        "phi_attr_std_formula": "phi_L/(4 pi) = (alpha_Y/(4 pi)) G L^2 rho_m / c^2 = (3/(32 pi^2)) alpha_Y Omega_m (H0 L/c)^2 (1+z)^3",
        "H0_L_over_c": H0_S_INV * L_M / C_SI,
        "Omega_m_0": OMEGA_M_0_CHECK,
        "rho_crit_0_kg_m3": RHO_CRIT_0_SI,
    },
    "normalization": {
        "tilde_alpha_L": TILDE_ALPHA_L,
        "phi_canon_over_phi_L": PHI_CANON_OVER_PHI_L,
        "dlnA_dphi_L": -ALPHA_Y,
        "dlnA_dphi_canon": DLNA_DPHI_CANON,
        "beta_Yuk": BETA_YUK,
        "beta_Yuk_formula": "(d ln A/d phi_canon)^2 = 2 tilde_alpha_L alpha_Y^2 = alpha_Y^3/(4 pi)",
        "Geff_formula": "G_eff/G = 1 + beta_Yuk (kL/a)^2/(1+(kL/a)^2)",
    },
    "cosmological_background_models": background_models,
    "cmb_geff_benchmarks": {
        "prompt_k1_0.016_Mpc_inv": {
            "k_comoving_Mpc_inv": k1_exact_prompt,
            "k_phys_zstar_Mpc_inv": k1_phys_zstar,
            "k_phys_L_zstar": k1_phys_L_zstar,
            "beta_Yuk": BETA_YUK,
            "delta_Geff_over_G_zstar": delta_geff_k1,
        },
        "prompt_kSilk_0.22_Mpc_inv": {
            "k_comoving_Mpc_inv": k_silk_prompt,
            "k_phys_zstar_Mpc_inv": ksilk_phys_zstar,
            "k_phys_L_zstar": ksilk_phys_L_zstar,
            "beta_Yuk": BETA_YUK,
            "delta_Geff_over_G_zstar": delta_geff_ksilk,
        },
        "scale_grid": geff_rows,
    },
    "cmb_acoustic_peak_inputs_and_ISW": {
        "Planck18_external_benchmark_inputs": {
            "theta_star_100": THETA_STAR_PLANCK_100,
            "theta_star_err_100": THETA_STAR_PLANCK_ERR_100,
            "ell_A": ELL_A_PLANCK,
            "ell_1": ELL_1_PLANCK,
            "ell_2": ELL_2_PLANCK,
            "ell_3": ELL_3_PLANCK,
            "H2_H1_lit_benchmark_input": 0.4580,
            "H3_H2_lit_benchmark_input": 0.9850,
            "boltzmann_spectrum_computed": False,
        },
        "input_Planck18_CDM_to_baryon_ratio_at_zstar": OMEGA_C_H2 / OMEGA_B_H2,
        "Weyl_sum_identity_Jordan_vs_Einstein": "tilde_Psi + tilde_Phi = Psi_E + Phi_E (exact at linear order in conformal Jordan frame)",
        "ISW_anisotropic_slip_at_ell_le_20": geff_rows[2]["Weyl_anisotropic_slip_zstar"],
        "ISW_growth_modification_at_ell_le_20": geff_rows[2]["ISW_fractional_modification_z0"],
    },
    "grand_synthesis_6_scales": grand_synthesis_6_scales,
}

script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(script_dir) if os.path.basename(script_dir) == "code" else script_dir
out_csv_bg = os.path.join(base_dir, "cosmological_background_models.csv")
out_csv_geff = os.path.join(base_dir, "cmb_geff_and_acoustic_scales.csv")
out_json = os.path.join(base_dir, "cmb_and_cosmology_summary.json")

df_bg.to_csv(out_csv_bg, index=False, float_format="%.10g")
df_geff.to_csv(out_csv_geff, index=False, float_format="%.10g")
with open(out_json, "w", encoding="utf-8") as f:
    json.dump(summary_payload, f, indent=2, ensure_ascii=False, allow_nan=False)

script_self = os.path.join(base_dir, "cmb_and_cosmology_audit.py")
shutil.copy2(script_self, os.path.join(base_dir, "code", "cmb_and_cosmology_audit.py"))

for fname in [
    "cosmological_background_models.csv",
    "cmb_geff_and_acoustic_scales.csv",
    "cmb_and_cosmology_summary.json",
]:
    src = os.path.join(base_dir, fname)
    for subdir in ["code", "csv", os.path.join("data", "v23")]:
        dst = os.path.join(base_dir, subdir, fname)
        shutil.copy2(src, dst)

print("Module 3 (step25, V26) cosmology and CMB audit executed and synced.")
