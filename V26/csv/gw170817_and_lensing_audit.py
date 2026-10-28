#!/usr/bin/env python3
"""
Step 23 (V26) — Module 1: GW170817 Tensor Wave Speed (c_T = 1), Massive Scalar Mode
Dispersion, and Gravitational Lensing (Weyl Potential Phi + Psi, Scale-Dependent
PPN Slip gamma_PPN(r), Curved-Spacetime Null-Geodesic Deflection Angle alpha_hat(b),
Shapiro Time Delay, and Weak-Lensing Convergence kappa) at L = 10 km, alpha_Y = 1.

Author: David Szabolcsi (ORCID: 0009-0005-6300-074X)
License: CC BY-NC 4.0

Computes and verifies:
  1. Symbolic (SymPy) and numerical verification of the Horndeski mapping in both
     the Einstein frame and the Jordan frame:
       - Einstein frame: G_4 = c^4/(16*pi*G) = const, G_{4,X} = 0, G_5 = 0 -> c_T = c.
       - Jordan frame (single V26 convention: g^(J)_{mu nu} = A^2(phi_L) g_{mu nu}
         = exp(-2 alpha_Y phi_L) g_{mu nu}, A(phi_L) = exp(-alpha_Y phi_L), phi_L > 0 near matter,
         beta_Yuk = (d ln A/d phi_canon)^2 = 2*tilde_alpha_L*alpha_Y^2 = alpha_Y^3/(4 pi)): G_4_tilde(phi_L) = c^4/(16*pi*G) * A^{-2}(phi_L),
         partial G_4_tilde / partial X = 0, G_5_tilde = 0 -> c_T^2 = c^2 identically,
         (c_T - c)/c = 0 across all scales and parameters, satisfying the
         GW170817 / GRB 170817A bound -3e-15 <= (c_T - c)/c <= +7e-15.
         Scope: proven on the purely conformal minimal Einstein-Hilbert + canonical scalar +
         universal Maxwell branch (G_{4,X} = 0, G_5 = 0); not automatically on a
         disformal/refractive branch with D(phi_L, X) != 0.
  2. Massive scalar breathing mode dispersion (omega^2 = c^2 k^2 + c^2 / L^2) for
     L = 10 km (m_phi = 0.1 km^{-1}, f_cut = c / (2*pi*L) = 4771.345159 Hz,
     m_phi * hbar * c = 1.9732698e-11 eV):
       - Evanescent decay length lambda_ev(f) = L / sqrt(1 - (f/f_cut)^2) in the
         sub-cutoff binary pulsar and LIGO GW170817 inspiral band (f < f_cut).
       - Phase velocity v_p(f), group velocity v_g(f), and arrival time delay
         Delta t_scalar(f) over D = 40 Mpc for super-cutoff modes (f > f_cut).
  3. Exact weak-field potentials in Einstein and Jordan frames, conformal
     cancellation in the linear Weyl lensing potential Phi_lens = (Psi_tilde + Phi_tilde)/2,
     and scale-dependent PPN slip gamma_PPN(r) and lensing-to-dynamical ratio
     Phi_lens(r) / Psi_tilde(r) from r = 10 km out to r = 1 Mpc.
  4. Exact curved-spacetime null-geodesic light deflection angle alpha_hat(b) and
     Shapiro time delay Delta t_Shapiro(b) integrated on the V25 4-component
     exterior Yukawa BVP solution (nu(r), m(r), phi_L(r), phi_L'(r)) for the
     reference SLy neutron star (log10(rho_0c) = 15.0, L = 10 km, alpha_Y = 1,
     R_* = 11.46572423 km, M_surf = 1.50766122 M_sun, M_ADM = 1.50999814 M_sun),
     compared against vacuum Schwarzschild(M_surf), vacuum Schwarzschild(M_ADM),
     and the disformal refractive branch across b in [14.749 km, 1 Mpc].
  5. Weak-lensing surface mass density and convergence kappa(b) across neutron-star,
     Solar-System (Cassini), galactic (SPARC / NGC 3198), and galaxy-cluster
     (Bullet Cluster 1E 0657-558) scales at L = 10 km, alpha_Y = 1.
"""

import os
import sys
import json
import shutil
import math

for extra_site in [
    "/workspace/scratch/venv_openitems/lib/python3.13/site-packages",
    "/tmp/tv/lib/python3.13/site-packages",
]:
    if os.path.isdir(extra_site) and extra_site not in sys.path:
        sys.path.append(extra_site)

import numpy as np
import pandas as pd
import sympy as sp
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(script_dir) if os.path.basename(script_dir) == "code" else script_dir
OUT_DIR = base_dir
for _p in (script_dir, base_dir):
    if _p not in sys.path:
        sys.path.insert(0, _p)
if OUT_DIR not in sys.path:
    sys.path.insert(0, OUT_DIR)

import tal_profiles as tp
import covariant_closure_50aa2db0 as cc
import curved_yukawa_adm_and_stability as cy

CODE_DIR = os.path.join(OUT_DIR, "code")
CSV_DIR = os.path.join(OUT_DIR, "csv")
DATA_V23_DIR = os.path.join(OUT_DIR, "data", "v23")

# Fundamental physical constants (CODATA / IAU exact consistent with step1-step22)
C_MS = 299792458.0                  # m / s
C_KMS = 299792.458                  # km / s
G_SI = 6.67430e-11                  # m^3 kg^{-1} s^{-2}
HBAR_EV_S = 6.582119569e-16         # eV * s
HBAR_C_EV_M = 1.973269804e-7        # eV * m
GM_SUN_KM = float(tp.st.GM)         # 1.4766250385 km
M_SUN_KG = 1.98847e30               # kg
R_SUN_KM = 695700.0                 # km
AU_KM = 149597870.7                 # km
PC_KM = 3.0856775814913673e13       # km
KPC_KM = 1.0e3 * PC_KM              # km
MPC_KM = 1.0e6 * PC_KM              # km
YR_S = 365.25 * 86400.0             # Julian year in seconds

# Baseline parameters for Module 1
L_KM = 10.0
L_M = L_KM * 1000.0
ALPHA_Y = 1.0
ALPHA_STD = ALPHA_Y / (4.0 * math.pi)
# Single action S[phi_L] (Step 22): kinetic prefactor c^4/(16 pi G tilde_alpha_L), tilde_alpha_L = alpha_Y/(8 pi),
# matter coupling A(phi_L) = exp(-alpha_Y phi_L), d ln A/d phi_L = -alpha_Y.
# DEF canonical field (kinetic prefactor 1/tilde_alpha_L -> 2): phi_canon = phi_L / sqrt(2 tilde_alpha_L),
# d ln A/d phi_canon = -alpha_Y sqrt(2 tilde_alpha_L); beta_Yuk = (d ln A/d phi_canon)^2 = 2 tilde_alpha_L alpha_Y^2.
ALPHA_TILDE_L = ALPHA_Y / (8.0 * math.pi)
PHI_CANON_OVER_PHI_L = 1.0 / math.sqrt(2.0 * ALPHA_TILDE_L)          # sqrt(4 pi/alpha_Y) = 3.5449077018
DLNA_DPHI_CANON = -ALPHA_Y * math.sqrt(2.0 * ALPHA_TILDE_L)          # alpha_DEF = d ln A/d phi_canon = -1/sqrt(4 pi) = -0.2820947918
ALPHA_DEF = DLNA_DPHI_CANON                                          # First-order DEF matter coupling derivative d ln A/d phi_canon
BETA_DEF_SECOND_DERIV = 0.0                                          # Second-order DEF derivative d^2 ln A/d phi_canon^2 = 0 for A = exp(-alpha_Y phi_L)
BETA_YUK = 2.0 * ALPHA_TILDE_L * (ALPHA_Y ** 2)                      # Fifth-force strength alpha_DEF^2 = 2 tilde_alpha_L alpha_Y^2 = alpha_Y^3/(4 pi) = 0.0795774715
BETA_DEF = abs(DLNA_DPHI_CANON)                                      # |alpha_DEF| = |d ln A/d phi_canon| (retained for backward key compatibility)
OMEGA_BD = 1.0 / (2.0 * BETA_YUK) - 1.5                              # Canonical single-action BD parameter 2 pi/alpha_Y^3 - 3/2 = 4.7831853072
OMEGA_BD_LEGACY_STEP5 = (1.0 - 3.0 * ALPHA_STD) / (2.0 * ALPHA_STD)  # Legacy step5 4*pi-convention BD parameter 4 pi/alpha_Y - 3/2 = 11.0663706144
GAMMA_L_SI = G_SI * (L_M ** 2) / (C_MS ** 2)

GL_NODES_128, GL_WEIGHTS_128 = np.polynomial.legendre.leggauss(128)


def run_symbolic_horndeski_and_lensing_checks():
    """Symbolic verification of Horndeski c_T = c and weak-field conformal cancellation."""
    phi, X, alpha_y, c_sym, G_sym, r_sym, L_sym, M_sym = sp.symbols(
        "phi X alpha_y c G r L M", positive=True, real=True
    )
    psi_E, phi_E = sp.symbols("psi_E phi_E", real=True)

    # 1. Einstein-frame Horndeski G_4 and c_T^2
    G4_E = c_sym**4 / (16 * sp.pi * G_sym)
    G4_E_X = sp.diff(G4_E, X)
    cT2_E = sp.simplify(G4_E / (G4_E - 2 * X * G4_E_X) * c_sym**2)

    # 2. Jordan-frame Horndeski G_4_tilde(phi) for the single V26 convention A(phi) = exp(-alpha_y * phi)
    A_conf_tal = sp.exp(-alpha_y * phi)
    G4_J_tal = G4_E * A_conf_tal**(-2)
    G4_J_tal_X = sp.diff(G4_J_tal, X)
    cT2_J_tal = sp.simplify(G4_J_tal / (G4_J_tal - 2 * X * G4_J_tal_X) * c_sym**2)

    # 3. Conformal cancellation in weak-field lensing potential Phi_lens = (Psi_tilde + Phi_tilde)/2
    phi_pert = sp.Symbol("phi_pert", real=True)
    Psi_tilde = psi_E - alpha_y * c_sym**2 * phi_pert
    Phi_tilde = phi_E + alpha_y * c_sym**2 * phi_pert
    Phi_lens_sum = sp.simplify(Psi_tilde + Phi_tilde)
    conformal_cancellation_residual = sp.simplify(Phi_lens_sum - (psi_E + phi_E))

    # 4. PPN slip gamma_PPN(r) = Phi_tilde(r) / Psi_tilde(r) and lensing ratio (1 + gamma_PPN)/2
    beta_yuk = sp.Symbol("beta_yuk", positive=True)
    exp_factor = sp.exp(-r_sym / L_sym)
    Psi_yuk = -(G_sym * M_sym / r_sym) * (1 + beta_yuk * exp_factor)
    Phi_yuk = -(G_sym * M_sym / r_sym) * (1 - beta_yuk * exp_factor)
    gamma_ppn_sym = sp.simplify(Phi_yuk / Psi_yuk)
    Phi_lens_yuk = sp.simplify(sp.Rational(1, 2) * (Psi_yuk + Phi_yuk))
    lensing_to_dyn_ratio = sp.simplify(Phi_lens_yuk / Psi_yuk)
    half_one_plus_gamma = sp.simplify((1 + gamma_ppn_sym) / 2)
    ratio_identity_residual = sp.simplify(lensing_to_dyn_ratio - half_one_plus_gamma)

    return {
        "Einstein_frame_G4_X": str(G4_E_X),
        "Einstein_frame_cT2_minus_c2": int(sp.simplify(cT2_E - c_sym**2)),
        "convention": "A(phi_L) = exp(-alpha_Y phi_L), phi_L > 0 near matter, g_J = exp(-2 alpha_Y phi_L) g, phi_canon = phi_L/sqrt(2 tilde_alpha_L), beta_Yuk = (d ln A/d phi_canon)^2 = 2 tilde_alpha_L alpha_Y^2 = alpha_Y^3/(4 pi)",
        "Jordan_frame_exp_minus_alphaY_phi_G4_X": str(G4_J_tal_X),
        "Jordan_frame_exp_minus_alphaY_phi_cT2_minus_c2": int(sp.simplify(cT2_J_tal - c_sym**2)),
        "cT_scope": "proven on conformal minimal EH + canonical scalar + universal Maxwell branch (G4_X = 0, G5 = 0); not automatic for disformal/refractive branch D(phi_L, X) != 0",
        "weak_field_Psi_tilde_plus_Phi_tilde": str(Phi_lens_sum),
        "conformal_cancellation_residual": int(conformal_cancellation_residual),
        "gamma_PPN_formula": str(gamma_ppn_sym),
        "Phi_lens_formula": str(Phi_lens_yuk),
        "lensing_to_dynamical_ratio_formula": str(lensing_to_dyn_ratio),
        "lensing_ratio_minus_half_one_plus_gamma_residual": int(ratio_identity_residual),
    }


def compute_gw170817_dispersion_table():
    """
    Compute GW170817 tensor speed c_T = c and massive scalar mode dispersion
    (omega^2 = c^2 k^2 + c^2/L^2) across binary-pulsar, LIGO inspiral/merger,
    and post-merger kHz frequencies for L = 10 km, alpha_Y = 1.
    """
    m_phi_km_inv = 1.0 / L_KM
    m_phi_m_inv = 1.0 / L_M
    omega_cut_rad_s = C_KMS / L_KM
    f_cut_hz = omega_cut_rad_s / (2.0 * math.pi)
    m_phi_ev = HBAR_C_EV_M / L_M

    D_gw170817_mpc = 40.0
    D_gw170817_km = D_gw170817_mpc * MPC_KM
    t_light_s = D_gw170817_km / C_KMS
    t_light_yr = t_light_s / YR_S

    freq_rows = [
        ("PSR_B1913p16_orbital", 3.58299e-5, "Binary pulsar orbital frequency (1/P_b)"),
        ("PSR_J0737m3039_orbital", 1.13199e-4, "Double pulsar orbital frequency (1/P_b)"),
        ("LISA_mHz_band", 1.0e-2, "Space-based mHz gravitational-wave band"),
        ("GW170817_LIGO_entry_24Hz", 24.0, "GW170817 LIGO band entry frequency"),
        ("GW170817_LIGO_100Hz", 100.0, "GW170817 LIGO maximum sensitivity band"),
        ("GW170817_LIGO_300Hz", 300.0, "GW170817 late inspiral band"),
        ("GW170817_LIGO_600Hz", 600.0, "GW170817 high-frequency inspiral band"),
        ("NS_post_merger_1500Hz", 1500.0, "Neutron-star tidal / bar-mode kHz band"),
        ("NS_f_mode_2500Hz", 2500.0, "Neutron-star fundamental quadrupole f-mode"),
        ("NS_high_3000Hz", 3000.0, "Upper LIGO post-merger search band"),
        ("Sub_cutoff_4500Hz", 4500.0, "Just below Yukawa scalar mass cutoff"),
        ("Yukawa_cutoff_f_cut", f_cut_hz, "Exact scalar cutoff frequency f_cut = c/(2*pi*L)"),
        ("Super_cutoff_4800Hz", 4800.0, "Just above Yukawa scalar cutoff (propagating)"),
        ("Super_cutoff_5000Hz", 5000.0, "Super-cutoff 5 kHz scalar breathing mode"),
        ("Super_cutoff_7500Hz", 7500.0, "Super-cutoff 7.5 kHz scalar breathing mode"),
        ("Super_cutoff_10000Hz", 10000.0, "Super-cutoff 10 kHz scalar breathing mode"),
    ]

    rows = []
    for label, f_hz, desc in freq_rows:
        omega = 2.0 * math.pi * f_hz
        ratio_f = f_hz / f_cut_hz
        c_T_over_c = 1.0
        frac_cT_minus_c = 0.0

        if abs(ratio_f - 1.0) < 1e-14:
            regime = "cutoff_threshold"
            k_real_km_inv = 0.0
            kappa_ev_km_inv = 0.0
            lambda_ev_km = float("inf")
            log10_suppression_40Mpc = 0.0
            v_group_over_c = 0.0
            v_phase_over_c = float("inf")
            scalar_delay_40Mpc_s = float("inf")
            scalar_delay_40Mpc_yr = float("inf")
        elif ratio_f < 1.0:
            regime = "evanescent_Yukawa_screened"
            k_real_km_inv = 0.0
            kappa_ev_km_inv = m_phi_km_inv * math.sqrt(1.0 - ratio_f**2)
            lambda_ev_km = 1.0 / kappa_ev_km_inv
            log10_suppression_40Mpc = -(D_gw170817_km / lambda_ev_km) / math.log(10.0)
            v_group_over_c = 0.0
            v_phase_over_c = 0.0
            scalar_delay_40Mpc_s = 0.0
            scalar_delay_40Mpc_yr = 0.0
        else:
            regime = "propagating_dispersive"
            k_real_km_inv = m_phi_km_inv * math.sqrt(ratio_f**2 - 1.0)
            kappa_ev_km_inv = 0.0
            lambda_ev_km = 0.0
            log10_suppression_40Mpc = 0.0
            v_group_over_c = math.sqrt(1.0 - (1.0 / ratio_f) ** 2)
            v_phase_over_c = 1.0 / v_group_over_c
            scalar_delay_40Mpc_s = t_light_s * (1.0 / v_group_over_c - 1.0)
            scalar_delay_40Mpc_yr = scalar_delay_40Mpc_s / YR_S

        rows.append({
            "mode_label": label,
            "frequency_Hz": f_hz,
            "omega_rad_s": omega,
            "f_over_f_cut": ratio_f,
            "scalar_regime": regime,
            "c_T_over_c_Einstein": c_T_over_c,
            "c_T_over_c_Jordan": c_T_over_c,
            "frac_c_T_minus_c": frac_cT_minus_c,
            "k_real_km_inv": k_real_km_inv,
            "kappa_evanescent_km_inv": kappa_ev_km_inv,
            "lambda_evanescent_km": lambda_ev_km,
            "log10_scalar_amp_at_40Mpc": log10_suppression_40Mpc,
            "scalar_v_group_over_c": v_group_over_c,
            "scalar_v_phase_over_c": v_phase_over_c,
            "scalar_delay_at_40Mpc_s": scalar_delay_40Mpc_s,
            "scalar_delay_at_40Mpc_yr": scalar_delay_40Mpc_yr,
            "description": desc,
        })

    df = pd.DataFrame(rows)
    meta = {
        "L_km": L_KM,
        "alpha_Y": ALPHA_Y,
        "m_phi_km_inv": m_phi_km_inv,
        "m_phi_m_inv": m_phi_m_inv,
        "omega_cut_rad_s": omega_cut_rad_s,
        "f_cut_Hz": f_cut_hz,
        "m_phi_eV": m_phi_ev,
        "D_GW170817_Mpc": D_gw170817_mpc,
        "D_GW170817_km": D_gw170817_km,
        "light_travel_time_40Mpc_s": t_light_s,
        "light_travel_time_40Mpc_yr": t_light_yr,
        "GW170817_observed_dt_s": 1.74,
        "GW170817_cT_lower_bound": -3.0e-15,
        "GW170817_cT_upper_bound": +7.0e-15,
        "TAL_predicted_frac_cT_minus_c_Einstein": 0.0,
        "TAL_predicted_frac_cT_minus_c_Jordan": 0.0,
    }
    return df, meta


def build_v25_reference_ns_exterior():
    """
    Solve the exact V25 reference SLy neutron star (log10(rho_0c) = 15.0, L = 10 km, alpha_Y = 1)
    and return the smooth BVP polynomial evaluator for m(r), nu(r), phi_L(r), phi_L'(r) on r >= R_*.
    """
    eos = tp.EOS("SLy")
    star = cy.solve_yukawa_action_el_curved_star(
        eos, 15.0, L=L_KM, alpha_Y=ALPHA_Y, ngrid=601, return_profile=True
    )
    R_star = float(star["R_km"])
    M_surf_Msun = float(star["M_surface_Msun"])
    M_surf_km = M_surf_Msun * GM_SUN_KM
    phi_R = float(star["phi_surface"])

    ext = star["profile"]["ext"]
    bvp_sol = ext["sol"]
    M_ADM_km = float(ext["M_ADM_km"])
    M_ADM_Msun = M_ADM_km / GM_SUN_KM
    R_max = float(ext["R_max"])

    phi_Rmax = float(ext["phi_ext"][-1])
    B_Rmax = 1.0 - 2.0 * M_ADM_km / R_max
    k_asympt = (1.0 / L_KM) / math.sqrt(B_Rmax)

    def eval_fields(r):
        r = float(r)
        if r <= R_max:
            y_v = bvp_sol.sol(r)
            phi_v = float(y_v[0])
            dphi_v = float(y_v[1])
            m_v = float(y_v[2])
            nu_v = float(y_v[3])
        else:
            m_v = M_ADM_km
            nu_v = 0.5 * math.log1p(-2.0 * M_ADM_km / r)
            dr = r - R_max
            if dr / L_KM > 700.0:
                phi_v = 0.0
                dphi_v = 0.0
            else:
                phi_v = phi_Rmax * (R_max / r) * math.exp(-k_asympt * dr)
                dphi_v = -phi_v * (1.0 / r + k_asympt)
        return m_v, nu_v, phi_v, dphi_v

    prof = star["profile"]
    rg = np.asarray(prof["rg"], dtype=float)
    er = np.asarray(prof["er"], dtype=float)
    phi_int = np.asarray(prof["phi"], dtype=float)
    src_int = np.exp(-phi_int) * er
    integrand_meff = 4.0 * math.pi * (rg ** 2) * np.where(
        rg > 1e-8, (L_KM / np.maximum(rg, 1e-8)) * np.sinh(rg / L_KM), 1.0
    ) * src_int
    M_eff_flat_km = float(np.trapezoid(integrand_meff, rg))
    M_eff_flat_Msun = M_eff_flat_km / GM_SUN_KM

    return {
        "star": star,
        "ext": ext,
        "R_star_km": R_star,
        "R_max_km": R_max,
        "M_surf_Msun": M_surf_Msun,
        "M_surf_km": M_surf_km,
        "M_ADM_Msun": M_ADM_Msun,
        "M_ADM_km": M_ADM_km,
        "Delta_M_phi_ext_Msun": float(ext["Delta_M_phi_ext_Msun"]),
        "Delta_M_phi_ext_km": M_ADM_km - M_surf_km,
        "phi_R": phi_R,
        "dphi_R_curved": float(ext["dphi_R_curved"]),
        "k_curved_km_inv": float(ext["k_curved"]),
        "nu_R_exact": float(ext["nu_R"]),
        "nu_R_Schw_Msurf": 0.5 * math.log1p(-2.0 * M_surf_km / R_star),
        "nu_R_Schw_MADM": 0.5 * math.log1p(-2.0 * M_ADM_km / R_star),
        "M_eff_flat_Msun": M_eff_flat_Msun,
        "eval_fields": eval_fields,
    }


def build_branch_ii_conf_trace_ns_exterior():
    """
    Solve the self-consistent Branch II-Conf (universal conformal Jordan coupling
    S_m[A^2(phi_L) g_{mu nu}, Psi_m], A(phi_L) = exp(-alpha_Y phi_L), sqrt(-g_tilde) = exp(-4 alpha_Y phi_L) sqrt(-g))
    interior + 4-component curved-exterior BVP for the same SLy neutron star
    (log10(rho_0c) = 15.0, L = 10 km, alpha_Y = 1):
      - Einstein-frame matter density & pressure:
          rho_m^(E,II-Conf) = exp(-4 alpha_Y phi_L) * eps(h),
          P_m^(E,II-Conf)   = exp(-4 alpha_Y phi_L) * P(h),
      - Scalar trace source (T_tilde = -eps + 3 P):
          B phi_L'' + (2B/r + B nu' + B'/2) phi_L' - phi_L/L^2 = -alpha_Y^2 exp(-4 alpha_Y phi_L) (eps - 3 P),
      - Jordan-frame hydrostatic Euler-TOV balance:
          dh/dr = dP/(eps + P) = -nu'(r) + alpha_Y phi_L'(r),
      - Exterior vacuum (r > R_*, eps = P = 0):
          Identical 4-component curved Einstein-scalar ODE system as Branch II-Fluid,
          evaluated at the Branch II-Conf surface boundary values (R_*, M_surf, phi_L(R_*)).
    """
    eos = tp.EOS("SLy")
    ngrid = 601

    def integrate_conf_tov(logrho, phi_fn, rho_phi_fn, pr_phi_fn, rtol=1e-9):
        hc = eos.h_of_rho(10.0 ** logrho)
        phi0_val = float(phi_fn(0.0))

        def hof(r, y0):
            return y0 + ALPHA_Y * (float(phi_fn(r)) - phi0_val)

        def rhs(r, y):
            h = hof(r, y[0])
            m = y[1]
            e, P, _, _ = eos.at(max(h, 1e-12))
            a4 = math.exp(-4.0 * ALPHA_Y * float(phi_fn(r)))
            nup = (m + 4.0 * math.pi * (r ** 3) * (a4 * P + float(pr_phi_fn(r)))) / (r * (r - 2.0 * m))
            dm = 4.0 * math.pi * (r ** 2) * (a4 * e + float(rho_phi_fn(r)))
            return [-nup, dm, nup]

        ev = lambda r, y: hof(r, y[0]) - 1e-10
        ev.terminal = True
        r0 = 1e-4
        ec = eos.at(hc)[0]
        a4_0 = math.exp(-4.0 * ALPHA_Y * phi0_val)
        sol = solve_ivp(
            rhs,
            (r0, 40.0),
            [hc, (4.0 * math.pi / 3.0) * (a4_0 * ec + float(rho_phi_fn(0.0))) * (r0 ** 3), 0.0],
            events=ev,
            rtol=rtol,
            atol=[1e-12, 1e-12, 1e-12],
            dense_output=True,
            max_step=0.05,
        )
        R = float(sol.t_events[0][0])
        M_km = float(sol.y_events[0][0][1])
        rg = np.linspace(0.0, R, ngrid)
        Y = sol.sol(np.maximum(rg, r0))
        hg = np.array([hof(r, y0) for r, y0 in zip(rg, Y[0])], dtype=float)
        hg[0] = hc
        hg[-1] = 0.0
        mg = Y[1].copy()
        mg[0] = 0.0
        er, Pr, _, _ = eos.at(np.maximum(hg, 1e-12))
        er[-1] = 0.0
        Pr[-1] = 0.0
        return {"R": R, "M_Msun": M_km / GM_SUN_KM, "rg": rg, "mg": mg, "hg": hg, "er": er, "Pr": Pr}

    zero_fn = lambda r: 0.0
    cur = integrate_conf_tov(15.0, zero_fn, zero_fn, zero_fn)
    phi = np.zeros(ngrid)
    x_prev = cur["rg"] / cur["rg"][-1]
    rho_p = np.zeros(ngrid)
    pr_p = np.zeros(ngrid)
    k_robin = 1.0 / cur["R"] + 1.0 / L_KM
    for it in range(100):
        rg, er, Pr, mg, hg = cur["rg"], cur["er"], cur["Pr"], cur["mg"], cur["hg"]
        phi0 = np.interp(rg / rg[-1], x_prev, phi)
        rho0 = np.interp(rg / rg[-1], x_prev, rho_p)
        pr0 = np.interp(rg / rg[-1], x_prev, pr_p)
        a4 = np.exp(-4.0 * ALPHA_Y * phi0)
        rs = np.maximum(rg, 1e-6)
        B = 1.0 - 2.0 * mg / rs
        B[0] = 1.0
        nup = (mg + 4.0 * math.pi * (rs ** 3) * (a4 * Pr + pr0)) / (rs * (rs - 2.0 * mg))
        nup[0] = 0.0
        dm = 4.0 * math.pi * (rs ** 2) * (a4 * er + rho0)
        dB = 2.0 * mg / (rs ** 2) - 2.0 * dm / rs
        dB[0] = 0.0
        if it in (1, 3, 6, 10) and abs(phi0[-1]) > 1e-14:
            ext_it = cy._solve_exterior_3comp_robin_iterate(rg[-1], mg[-1], float(phi0[-1]), L_KM, alpha_Y=ALPHA_Y)
            k_robin = ext_it["k_curved"]
        src_conf = (ALPHA_Y ** 2) * a4 * (er - 3.0 * Pr)
        phin = cy._solve_curved_scalar_with_k(rg, B, dB, nup, src_conf, L_KM, k_robin)
        phi_new = phin if it == 0 else 0.5 * phin + 0.5 * phi0
        dphi = np.gradient(phi_new, rg, edge_order=2)
        dphi[0] = 0.0
        dphi[-1] = -k_robin * phi_new[-1]
        rho_p = (B * (dphi ** 2) + (phi_new ** 2) / (L_KM ** 2)) / (2.0 * ALPHA_Y)
        pr_p = (B * (dphi ** 2) - (phi_new ** 2) / (L_KM ** 2)) / (2.0 * ALPHA_Y)
        x_prev, phi = rg / rg[-1], phi_new
        Rc = rg[-1]
        mk = lambda arr: (lambda I, Rr: (lambda r: float(I(min(r, Rr)))))(PchipInterpolator(rg, arr), Rc)
        new = integrate_conf_tov(15.0, mk(phi), mk(rho_p), mk(pr_p))
        res = max(abs(new["R"] / cur["R"] - 1.0), abs(new["M_Msun"] / cur["M_Msun"] - 1.0), abs(phi_new[0] - phi0[0]))
        cur = new
        if res < 2e-8 and it > 3:
            break

    # Evaluate interior scalar ODE residual and boundary derivatives on converged grid
    rg_f, er_f, Pr_f, mg_f, hg_f = cur["rg"], cur["er"], cur["Pr"], cur["mg"], cur["hg"]
    phi_f = np.interp(rg_f / rg_f[-1], x_prev, phi)
    a4_f = np.exp(-4.0 * ALPHA_Y * phi_f)
    rs_f = np.maximum(rg_f, 1e-6)
    B_f = 1.0 - 2.0 * mg_f / rs_f
    B_f[0] = 1.0
    dphi_f = np.gradient(phi_f, rg_f, edge_order=2)
    dphi_f[0] = 0.0
    dphi_f[-1] = -k_robin * phi_f[-1]
    rho_pf = (B_f * (dphi_f ** 2) + (phi_f ** 2) / (L_KM ** 2)) / (2.0 * ALPHA_Y)
    pr_pf = (B_f * (dphi_f ** 2) - (phi_f ** 2) / (L_KM ** 2)) / (2.0 * ALPHA_Y)
    nup_f = (mg_f + 4.0 * math.pi * (rs_f ** 3) * (a4_f * Pr_f + pr_pf)) / (rs_f * (rs_f - 2.0 * mg_f))
    nup_f[0] = 0.0
    dm_f = 4.0 * math.pi * (rs_f ** 2) * (a4_f * er_f + rho_pf)
    dB_f = 2.0 * mg_f / (rs_f ** 2) - 2.0 * dm_f / rs_f
    dB_f[0] = 0.0
    src_conf_f = (ALPHA_Y ** 2) * a4_f * (er_f - 3.0 * Pr_f)
    phi_banded_f = cy._solve_curved_scalar_with_k(rg_f, B_f, dB_f, nup_f, src_conf_f, L_KM, k_robin)
    h_step = float(rg_f[1] - rg_f[0])
    b_coef_f = 2.0 * B_f[1:-1] / rg_f[1:-1] + B_f[1:-1] * nup_f[1:-1] + 0.5 * dB_f[1:-1]
    ddphi_iter = (phi_f[2:] - 2.0 * phi_f[1:-1] + phi_f[:-2]) / (h_step ** 2)
    dphi_iter_c = (phi_f[2:] - phi_f[:-2]) / (2.0 * h_step)
    ode_res_iter = B_f[1:-1] * ddphi_iter + b_coef_f * dphi_iter_c - phi_f[1:-1] / (L_KM ** 2) + src_conf_f[1:-1]
    ddphi_band = (phi_banded_f[2:] - 2.0 * phi_banded_f[1:-1] + phi_banded_f[:-2]) / (h_step ** 2)
    dphi_band_c = (phi_banded_f[2:] - phi_banded_f[:-2]) / (2.0 * h_step)
    ode_res_band = B_f[1:-1] * ddphi_band + b_coef_f * dphi_band_c - phi_banded_f[1:-1] / (L_KM ** 2) + src_conf_f[1:-1]

    aN_f = float(B_f[-1])
    bN_f = float(2.0 * B_f[-1] / rg_f[-1] + B_f[-1] * nup_f[-1] + 0.5 * dB_f[-1])
    dphi_int_robin = float(-k_robin * phi[-1])
    dphi_int_ghost = float(
        (-(2.0 * aN_f / (h_step ** 2)) * phi[-2] + (2.0 * aN_f / (h_step ** 2) + 1.0 / (L_KM ** 2)) * phi[-1])
        / (2.0 * h_step * (aN_f / (h_step ** 2) + bN_f / (2.0 * h_step)))
    )
    dphi_int_fd2 = float((3.0 * phi[-1] - 4.0 * phi[-2] + phi[-3]) / (2.0 * h_step))

    R_star = float(cur["R"])
    M_surf_Msun = float(cur["M_Msun"])
    M_surf_km = M_surf_Msun * GM_SUN_KM
    phi_R = float(phi[-1])
    ext = cy.solve_exterior_dtN_curved(R_star, M_surf_km, phi_R, L_KM, alpha_Y=ALPHA_Y, tol=1e-12)
    bvp_sol = ext["sol"]
    M_ADM_km = float(ext["M_ADM_km"])
    M_ADM_Msun = M_ADM_km / GM_SUN_KM
    R_max = float(ext["R_max"])
    phi_Rmax = float(ext["phi_ext"][-1])
    B_Rmax = 1.0 - 2.0 * M_ADM_km / R_max
    k_asympt = (1.0 / L_KM) / math.sqrt(B_Rmax)

    def eval_fields(r):
        r = float(r)
        if r <= R_max:
            y_v = bvp_sol.sol(r)
            return float(y_v[2]), float(y_v[3]), float(y_v[0]), float(y_v[1])
        else:
            dr = r - R_max
            phi_v = 0.0 if dr / L_KM > 700.0 else phi_Rmax * (R_max / r) * math.exp(-k_asympt * dr)
            dphi_v = 0.0 if dr / L_KM > 700.0 else -phi_v * (1.0 / r + k_asympt)
            return M_ADM_km, 0.5 * math.log1p(-2.0 * M_ADM_km / r), phi_v, dphi_v

    return {
        "R_star_km": R_star,
        "R_max_km": R_max,
        "M_surf_Msun": M_surf_Msun,
        "M_surf_km": M_surf_km,
        "M_ADM_Msun": M_ADM_Msun,
        "M_ADM_km": M_ADM_km,
        "Delta_M_phi_ext_Msun": float(ext["Delta_M_phi_ext_Msun"]),
        "Delta_M_phi_ext_rel_Msurf": float(ext["Delta_M_phi_ext_Msun"]) / M_surf_Msun,
        "phi_0": float(phi[0]),
        "phi_R": phi_R,
        "dphi_R_int_robin": dphi_int_robin,
        "dphi_R_int_ghost": dphi_int_ghost,
        "dphi_R_int_fd2": dphi_int_fd2,
        "dphi_R_curved": float(ext["dphi_R_curved"]),
        "k_robin_int_km_inv": float(k_robin),
        "k_curved_km_inv": float(ext["k_curved"]),
        "abs_diff_k_int_vs_k_curved": abs(float(k_robin) - float(ext["k_curved"])),
        "rel_diff_k_int_vs_k_curved": abs(float(k_robin) - float(ext["k_curved"])) / float(ext["k_curved"]),
        "abs_diff_dphi_int_robin_vs_ext_curved": abs(dphi_int_robin - float(ext["dphi_R_curved"])),
        "rel_diff_dphi_int_robin_vs_ext_curved": abs(dphi_int_robin - float(ext["dphi_R_curved"])) / abs(float(ext["dphi_R_curved"])),
        "abs_diff_dphi_int_fd2_vs_ext_curved": abs(dphi_int_fd2 - float(ext["dphi_R_curved"])),
        "max_interior_scalar_ode_residual_iter_km_inv2": float(np.max(np.abs(ode_res_iter))),
        "max_interior_scalar_ode_residual_banded_km_inv2": float(np.max(np.abs(ode_res_band))),
        "max_exterior_4comp_bvp_rms_residual": float(np.max(bvp_sol.rms_residuals)),
        "nu_R_exact": float(ext["nu_R"]),
        "nu_R_Schw_Msurf": 0.5 * math.log1p(-2.0 * M_surf_km / R_star),
        "nu_R_Schw_MADM": 0.5 * math.log1p(-2.0 * M_ADM_km / R_star),
        "center_3P_over_eps": float(3.0 * cur["Pr"][0] / cur["er"][0]),
        "center_trace_factor_1_minus_3P_over_eps": float(1.0 - 3.0 * cur["Pr"][0] / cur["er"][0]),
        "eval_fields": eval_fields,
    }


def compute_exact_deflection_angle(b_km, eval_fields, R_star_km, R_max_km, M_asymp_km, branch="conformal"):
    """
    Compute the exact curved-spacetime null-geodesic deflection angle alpha_hat(b) [rad]
    using 128-point Gauss-Legendre quadrature in u = sin(theta) with cancellation-free
    algebraic formulation from b = 14.749 km out to b = 1 Mpc.
    """
    def nu_eff_fn(r):
        _, nu_v, phi_v, _ = eval_fields(r)
        return nu_v - phi_v if branch == "disformal" else nu_v

    def b_of_r0(r0):
        return r0 * math.exp(-nu_eff_fn(r0))

    if b_km <= 1.0e5:
        r0_low = R_star_km
        b_low = b_of_r0(r0_low)
        if b_km < b_low:
            raise ValueError(f"Impact parameter b={b_km} km is inside stellar surface b_min={b_low:.6f} km")
        r0 = float(brentq(lambda r: b_of_r0(r) - b_km, r0_low, b_km * 1.05, xtol=1e-14, rtol=1e-14))
    else:
        x = M_asymp_km / b_km
        r0 = b_km * (1.0 - x - 1.5 * (x ** 2) - 4.0 * (x ** 3))
        for _ in range(3):
            f = r0 / math.sqrt(1.0 - 2.0 * M_asymp_km / r0) - b_km
            df = (1.0 - 3.0 * M_asymp_km / r0) / ((1.0 - 2.0 * M_asymp_km / r0) ** 1.5)
            r0 -= f / df

    _, nu0, phi0, _ = eval_fields(r0)
    nu_eff_0 = nu0 - phi0 if branch == "disformal" else nu0

    half_pi = 0.5 * math.pi
    thetas = 0.5 * half_pi * (GL_NODES_128 + 1.0)
    w_scaled = 0.5 * half_pi * GL_WEIGHTS_128

    total = 0.0
    for th, wt in zip(thetas, w_scaled):
        u = math.sin(th)
        c_th = math.cos(th)
        r = r0 / u
        m_v, nu_v, phi_v, _ = eval_fields(r)
        B_minus_1 = -2.0 * m_v / r
        B_v = 1.0 + B_minus_1

        if r0 > R_max_km and phi0 == 0.0 and phi_v == 0.0:
            one_minus_u = (c_th * c_th) / (1.0 + u)
            Delta0 = math.log1p(-2.0 * M_asymp_km * one_minus_u / (r0 * (1.0 - 2.0 * M_asymp_km * u / r0)))
        else:
            nu_eff_v = nu_v - phi_v if branch == "disformal" else nu_v
            Delta0 = 2.0 * (nu_eff_0 - nu_eff_v)

        expm1_Delta0 = math.expm1(Delta0)
        C_fac = expm1_Delta0 / (c_th * c_th)
        Y = B_minus_1 + B_v * C_fac
        one_plus_Y = max(1.0 + Y, 1e-30)
        integrand = -Y / (math.sqrt(one_plus_Y) * (1.0 + math.sqrt(one_plus_Y)))
        total += wt * integrand

    return float(2.0 * total), float(r0)


def compute_curved_shapiro_delay_two_sided(r0_km, r_end_km, eval_fields, R_max_km, M_asymp_km, branch="conformal"):
    """
    Compute the curved-spacetime two-sided grazing Shapiro time delay [us] from
    -r_end to +r_end with closest approach radius r0_km using cancellation-free
    128-point Gauss-Legendre quadrature in r = r0 * cosh(w).
    """
    if r_end_km <= r0_km:
        return 0.0

    _, nu0, phi0, _ = eval_fields(r0_km)
    nu_eff_0 = nu0 - phi0 if branch == "disformal" else nu0

    w_max = math.acosh(r_end_km / r0_km)
    ws = 0.5 * w_max * (GL_NODES_128 + 1.0)
    w_scaled = 0.5 * w_max * GL_WEIGHTS_128

    total_km = 0.0
    for w, wt in zip(ws, w_scaled):
        ch = math.cosh(w)
        sh = math.sinh(w)
        r = r0_km * ch
        m_v, nu_v, phi_v, _ = eval_fields(r)
        nu_eff_v = nu_v - phi_v if branch == "disformal" else nu_v

        expm1_2nu = math.expm1(2.0 * nu_eff_v)
        A_fac = expm1_2nu - (2.0 * m_v / r) * math.exp(2.0 * nu_eff_v)

        if r0_km > R_max_km and phi0 == 0.0 and phi_v == 0.0:
            Delta = math.log1p(
                (2.0 * M_asymp_km * sh * sh)
                / (r0_km * ch * (ch + 1.0) * (1.0 - 2.0 * M_asymp_km / r0_km))
            )
        else:
            Delta = 2.0 * (nu_eff_v - nu_eff_0)

        expm1_Delta = math.expm1(Delta)
        C_fac = -expm1_Delta / (sh * sh)
        X = A_fac + C_fac + A_fac * C_fac
        one_plus_X = max(1.0 + X, 1e-30)
        delta_integrand = -X / (math.sqrt(one_plus_X) * (1.0 + math.sqrt(one_plus_X)))
        total_km += wt * (r0_km * ch * delta_integrand)

    return float((2.0 * total_km / C_KMS) * 1.0e6)


def compute_lensing_and_ppn_table(ns_ref):
    """
    Compute the exact curved-spacetime deflection angle alpha_hat(b), Shapiro delay,
    scale-dependent PPN slip gamma_PPN(r), lensing-to-dynamical ratio Phi_lens/Psi_tilde,
    and weak-lensing convergence ratio across all scales from the neutron-star surface
    (r = 11.47 km) out to galaxy clusters (1 Mpc).
    """
    R_star = ns_ref["R_star_km"]
    R_max = ns_ref["R_max_km"]
    M_surf_km = ns_ref["M_surf_km"]
    M_ADM_km = ns_ref["M_ADM_km"]
    eval_yukawa = ns_ref["eval_fields"]

    def eval_schw_surf(r):
        r = float(r)
        return M_surf_km, 0.5 * math.log1p(-2.0 * M_surf_km / r), 0.0, 0.0

    def eval_schw_adm(r):
        r = float(r)
        return M_ADM_km, 0.5 * math.log1p(-2.0 * M_ADM_km / r), 0.0, 0.0

    b_surf_yukawa = R_star * math.exp(-ns_ref["nu_R_exact"])
    b_surf_disf = R_star * math.exp(-(ns_ref["nu_R_exact"] - ns_ref["phi_R"]))
    b_grazing = max(b_surf_yukawa, b_surf_disf) * 1.0001

    scale_rows_spec = [
        ("NS_surface_grazing", b_grazing, "Grazing ray at NS surface r_0 ~ R_* = 11.466 km (b = 14.749 km)"),
        ("NS_near_halo_15km", 15.0, "Strong-field inner Yukawa cloud (b = 15 km)"),
        ("NS_Yukawa_cloud_20km", 20.0, "Inner Yukawa scalar cloud (b = 20 km = 2 L)"),
        ("NS_Yukawa_cloud_30km", 30.0, "Mid Yukawa scalar cloud (b = 30 km = 3 L)"),
        ("NS_Yukawa_cloud_50km", 50.0, "Outer Yukawa scalar cloud (b = 50 km = 5 L)"),
        ("NS_Yukawa_tail_100km", 100.0, "Asymptotic Yukawa tail (b = 100 km = 10 L)"),
        ("NS_Yukawa_tail_200km", 200.0, "Screened transition zone (b = 200 km = 20 L)"),
        ("NS_far_1000km", 1000.0, "Far zone around neutron star (b = 1000 km = 100 L)"),
        ("Sirius_B_radius_5844km", 5844.0, "White dwarf surface scale (b = 5844 km)"),
        ("NS_far_field_b_6p957e5_km", R_SUN_KM, "1.50999814 M_sun NS far-field benchmark impact parameter (b = 6.957e5 km)"),
        ("NS_far_field_b_1p113e6_km", 1.6 * R_SUN_KM, "1.50999814 M_sun NS far-field benchmark impact parameter (b = 1.11312e6 km)"),
        ("NS_far_field_b_1AU", AU_KM, "1.50999814 M_sun NS far-field benchmark impact parameter (b = 1 AU = 1.496e8 km)"),
        ("Outer_Solar_10_AU", 10.0 * AU_KM, "Saturn / outer Solar System scale (b = 1.496e9 km)"),
        ("Subparsec_1_pc", PC_KM, "Stellar / Einstein-ring parsec scale (b = 1 pc = 3.086e13 km)"),
        ("Galactic_1_kpc", KPC_KM, "Inner galactic bulge/disk scale (b = 1 kpc = 3.086e16 km)"),
        ("NGC3198_30p08_kpc", 30.08 * KPC_KM, "SPARC NGC 3198 outer disk control radius (b = 30.08 kpc)"),
        ("Bullet_Cluster_250_kpc", 250.0 * KPC_KM, "Bullet Cluster 1E 0657-558 subcluster offset (b = 250 kpc)"),
        ("Cluster_virial_1_Mpc", MPC_KM, "Galaxy cluster virial radius scale (b = 1 Mpc = 3.086e19 km)"),
    ]

    rows = []
    for label, b_km, desc in scale_rows_spec:
        r_over_L = b_km / L_KM
        if r_over_L > 700.0:
            exp_minus_r_L = 0.0
            log10_exp_minus_r_L = -r_over_L / math.log(10.0)
        else:
            exp_minus_r_L = math.exp(-r_over_L)
            log10_exp_minus_r_L = -r_over_L / math.log(10.0)

        # Single V26 normalization: beta_Yuk = alpha_Y^2 / (4*pi) = 1/(4*pi)
        gamma_ppn_tal = (1.0 - BETA_YUK * exp_minus_r_L) / (1.0 + BETA_YUK * exp_minus_r_L)
        abs_one_minus_gamma_tal = 2.0 * BETA_YUK * exp_minus_r_L / (1.0 + BETA_YUK * exp_minus_r_L)
        phi_lens_over_psi_tal = 1.0 / (1.0 + BETA_YUK * exp_minus_r_L)
        log10_abs_1_minus_gamma = math.log10(2.0 * BETA_YUK) + log10_exp_minus_r_L

        alpha_yuk_conf_rad, r0_yuk = compute_exact_deflection_angle(
            b_km, eval_yukawa, R_star, R_max, M_ADM_km, branch="conformal"
        )
        alpha_yuk_disf_rad, _ = compute_exact_deflection_angle(
            b_km, eval_yukawa, R_star, R_max, M_ADM_km, branch="disformal"
        )
        alpha_schw_adm_rad, _ = compute_exact_deflection_angle(
            b_km, eval_schw_adm, R_star, 0.0, M_ADM_km, branch="conformal"
        )
        alpha_schw_surf_rad, _ = compute_exact_deflection_angle(
            b_km, eval_schw_surf, R_star, 0.0, M_surf_km, branch="conformal"
        )

        rel_diff_vs_schw_adm = (alpha_yuk_conf_rad - alpha_schw_adm_rad) / alpha_schw_adm_rad
        rel_diff_vs_schw_surf = (alpha_yuk_conf_rad - alpha_schw_surf_rad) / alpha_schw_surf_rad
        rel_diff_disf_vs_conf = (alpha_yuk_disf_rad - alpha_yuk_conf_rad) / alpha_yuk_conf_rad

        m_at_r0_km, _, phi_at_r0, _ = eval_yukawa(r0_yuk)
        m_at_r0_Msun = m_at_r0_km / GM_SUN_KM

        r_end_shapiro_km = max(1.0e6, 10.0 * r0_yuk)
        dt_shap_yuk_conf_us = compute_curved_shapiro_delay_two_sided(
            r0_yuk, r_end_shapiro_km, eval_yukawa, R_max, M_ADM_km, branch="conformal"
        )
        dt_shap_yuk_disf_us = compute_curved_shapiro_delay_two_sided(
            r0_yuk, r_end_shapiro_km, eval_yukawa, R_max, M_ADM_km, branch="disformal"
        )
        dt_shap_schw_adm_us = compute_curved_shapiro_delay_two_sided(
            r0_yuk, r_end_shapiro_km, eval_schw_adm, 0.0, M_ADM_km, branch="conformal"
        )
        dt_shap_schw_surf_us = compute_curved_shapiro_delay_two_sided(
            r0_yuk, r_end_shapiro_km, eval_schw_surf, 0.0, M_surf_km, branch="conformal"
        )

        rows.append({
            "source_body": "NS_SLy_1p51Msun",
            "scale_label": label,
            "b_km": b_km,
            "r0_closest_km": r0_yuk,
            "b_over_L": r_over_L,
            "exp_minus_b_over_L": exp_minus_r_L,
            "log10_exp_minus_b_over_L": log10_exp_minus_r_L,
            "m_enclosed_at_r0_Msun": m_at_r0_Msun,
            "phi_L_at_r0": phi_at_r0,
            "beta_Yuk": BETA_YUK,
            "gamma_PPN": gamma_ppn_tal,
            "abs_1_minus_gamma_PPN": abs_one_minus_gamma_tal,
            "log10_abs_1_minus_gamma_PPN": log10_abs_1_minus_gamma,
            "Phi_lens_over_Psi_dyn": phi_lens_over_psi_tal,
            "alpha_hat_Yukawa_conf_rad": alpha_yuk_conf_rad,
            "alpha_hat_Yukawa_disf_rad": alpha_yuk_disf_rad,
            "alpha_hat_Schw_MADM_rad": alpha_schw_adm_rad,
            "alpha_hat_Schw_Msurf_rad": alpha_schw_surf_rad,
            "rel_diff_Yukawa_conf_vs_Schw_MADM": rel_diff_vs_schw_adm,
            "rel_diff_Yukawa_conf_vs_Schw_Msurf": rel_diff_vs_schw_surf,
            "rel_diff_Yukawa_disf_vs_conf": rel_diff_disf_vs_conf,
            "Shapiro_delay_Yukawa_conf_us": dt_shap_yuk_conf_us,
            "Shapiro_delay_Yukawa_disf_us": dt_shap_yuk_disf_us,
            "Shapiro_delay_Schw_MADM_us": dt_shap_schw_adm_us,
            "Shapiro_delay_Schw_Msurf_us": dt_shap_schw_surf_us,
            "description": desc,
        })

    return pd.DataFrame(rows)


def compute_solar_cassini_1msun_table():
    """Dedicated 1.000000 M_sun Solar-Cassini benchmark (vacuum Schwarzschild, GM_sun/c^2 = 1.4766250385 km)."""
    M_km = 1.0 * GM_SUN_KM
    r_E = AU_KM
    r_S = 8.43 * AU_KM

    def eval_sun(r):
        r = float(r)
        return M_km, 0.5 * math.log1p(-2.0 * M_km / r), 0.0, 0.0

    spec = [
        ("Solar_1Msun_limb_1_Rsun", R_SUN_KM, True),
        ("Solar_1Msun_Cassini_1p6_Rsun", 1.6 * R_SUN_KM, True),
        ("Solar_1Msun_10_Rsun", 10.0 * R_SUN_KM, False),
        ("Solar_1Msun_Earth_1_AU", AU_KM, False),
    ]
    rad2as = 180.0 / math.pi * 3600.0
    rows = []
    for label, b_km, cassini in spec:
        a_ex, r0 = compute_exact_deflection_angle(b_km, eval_sun, 0.0, 0.0, M_km)
        a_1pn = 4.0 * M_km / b_km
        dt10 = compute_curved_shapiro_delay_two_sided(r0, 10.0 * r0, eval_sun, 0.0, M_km)
        tM_us = M_km / C_KMS * 1.0e6
        if cassini:
            lnarg = 4.0 * r_E * r_S / (b_km ** 2)
            dt_rt = 4.0 * tM_us * math.log(lnarg)
            dt_rt_m1 = 4.0 * tM_us * (math.log(lnarg) - 1.0)
        else:
            dt_rt = float("nan")
            dt_rt_m1 = float("nan")
        log10_g = math.log10(2.0 * BETA_YUK) - (b_km / L_KM) / math.log(10.0)
        rows.append({
            "source_body": "Sun_1p00Msun",
            "scale_label": label,
            "M_Msun": 1.0,
            "GM_over_c2_km": M_km,
            "b_km": b_km,
            "b_over_R_sun": b_km / R_SUN_KM,
            "r0_closest_km": r0,
            "alpha_hat_exact_arcsec": a_ex * rad2as,
            "alpha_hat_1PN_arcsec": a_1pn * rad2as,
            "Shapiro_two_sided_10r0_us": dt10,
            "Cassini_roundtrip_log_delay_us": dt_rt,
            "Cassini_roundtrip_log_delay_minus1_us": dt_rt_m1,
            "beta_Yuk": BETA_YUK,
            "log10_abs_1_minus_gamma_PPN": log10_g,
        })
    return pd.DataFrame(rows)


def compute_macroscopic_lensing_and_bullet_cluster_bridge():
    """
    Quantify weak-lensing convergence kappa(b) = Sigma_lens(b)/Sigma_crit on
    galactic (NGC 3198, R = 30.08 kpc) and cluster (Bullet Cluster 1E 0657-558,
    b = 250 kpc) scales for:
      (i)  The pure local Yukawa mode at L = 10 km, alpha_Y = 1,
      (ii) The macroeconomic / disformal / chameleon extensions required for
           non-baryonic lensing enhancement on kpc-Mpc scales.
    """
    R_ngc_kpc = 30.08
    V_obs_kms = 146.0
    V_bar_u05_kms = 70.32
    V_bar_u10_kms = 88.53
    C_req_u05 = (V_obs_kms / V_bar_u05_kms) ** 2
    C_req_u10 = (V_obs_kms / V_bar_u10_kms) ** 2

    rho_gal_kg_m3 = 1.0e-24
    phi_L_local_10km_gal = ALPHA_Y * GAMMA_L_SI * rho_gal_kg_m3
    Phi_N_over_c2_gal = (V_bar_u05_kms / C_KMS) ** 2

    Sigma_crit_Msun_kpc2 = 3.10e9
    Sigma_bary_gal_peak_Msun_kpc2 = 0.08e9
    Sigma_obs_lensing_peak_Msun_kpc2 = 0.72e9
    kappa_bary_gal_peak = Sigma_bary_gal_peak_Msun_kpc2 / Sigma_crit_Msun_kpc2
    kappa_obs_peak = Sigma_obs_lensing_peak_Msun_kpc2 / Sigma_crit_Msun_kpc2
    required_kappa_boost_bullet = kappa_obs_peak / kappa_bary_gal_peak

    rho_icm_kg_m3 = 1.0e-22
    phi_L_local_10km_cluster = ALPHA_Y * GAMMA_L_SI * rho_icm_kg_m3

    return {
        "NGC3198_R_kpc": R_ngc_kpc,
        "NGC3198_b_over_L10km": (R_ngc_kpc * KPC_KM) / L_KM,
        "NGC3198_C_req_Upsilon_0p5": C_req_u05,
        "NGC3198_C_req_Upsilon_1p0": C_req_u10,
        "NGC3198_phi_L_local_L10km": phi_L_local_10km_gal,
        "NGC3198_Phi_N_over_c2": Phi_N_over_c2_gal,
        "NGC3198_ratio_phi_L10km_to_Phi_N": phi_L_local_10km_gal / Phi_N_over_c2_gal,
        "NGC3198_kappa_lens_over_kappa_bary_local_L10km": 1.0,
        "Bullet_Cluster_offset_kpc": 250.0,
        "Bullet_Cluster_b_over_L10km": (250.0 * KPC_KM) / L_KM,
        "Bullet_Cluster_Sigma_crit_Msun_kpc2": Sigma_crit_Msun_kpc2,
        "Bullet_Cluster_kappa_bary_gal_subcluster": kappa_bary_gal_peak,
        "Bullet_Cluster_kappa_obs_subcluster": kappa_obs_peak,
        "Bullet_Cluster_required_kappa_ratio": required_kappa_boost_bullet,
        "Bullet_Cluster_phi_L_local_L10km_in_ICM": phi_L_local_10km_cluster,
        "Bullet_Cluster_kappa_lens_over_kappa_bary_local_L10km": 1.0,
        "Bullet_Cluster_conformal_branch_kappa_over_kappa_bary": 1.0,
        "required_disformal_or_Weyl_potential_condition": (
            "Phi_lens_eff = Phi_N - (1/2) c^2 phi_macro with L_eff(rho) ~ 10-300 kpc; "
            "purely conformal coupling A^2(phi) g_{mu nu} cancels at linear order in (Psi_tilde + Phi_tilde)/2, "
            "and single-scale L = 10 km Yukawa is exponentially screened (exp(-b/L) = 0) at kpc-Mpc scales."
        ),
    }


def compute_branch_ii_conf_bvp_convergence_and_mesh_refinement():
    """
    Targeted numerical verification of the Step 23 Branch II-Conf SLy log10(rho_0c) = 15.0
    trace-sourced interior + 4-component curved-exterior BVP:
      1. Converged interior scalar ODE residual max_{r in (0, R_*)} |res_ODE| and exterior
         4-component BVP RMS collocation residual.
      2. Outer fixed-point iteration convergence history and comparison between the
         4-step Robin update schedule (it in {1, 3, 6, 10}) and every-iteration Robin updates
         (both 3-component and 4-component exterior DtN maps).
      3. Mesh refinement study across N = 601, N = 2000, and N = 4000 radial grid points
         (with high-precision DOP853 TOV integration rtol_ivp=1e-11, atol=1e-13, max_step=0.02 km,
         tol_iter=1e-10, and every-step curved Robin updates), testing empirical second-order
         O(h^2) convergence of R_*, M_ADM, phi_L(R_*), and the independent one-sided
         second-order finite-difference surface derivative phi'_{L,int,FD2}(R_*).
    """
    eos = tp.EOS("SLy")

    def run_variant(ngrid, update_every=True, use_4comp_robin=False, tol_iter=1e-10, rtol_ivp=1e-11, use_dop853=True):
        def integrate_conf_tov(logrho, phi_fn, rho_phi_fn, pr_phi_fn):
            hc = eos.h_of_rho(10.0 ** logrho)
            phi0_val = float(phi_fn(0.0))

            def hof(r, y0):
                return y0 + ALPHA_Y * (float(phi_fn(r)) - phi0_val)

            def rhs(r, y):
                h = hof(r, y[0])
                m = y[1]
                e, P, _, _ = eos.at(max(h, 1e-12))
                a4 = math.exp(-4.0 * ALPHA_Y * float(phi_fn(r)))
                nup = (m + 4.0 * math.pi * (r ** 3) * (a4 * P + float(pr_phi_fn(r)))) / (r * (r - 2.0 * m))
                dm = 4.0 * math.pi * (r ** 2) * (a4 * e + float(rho_phi_fn(r)))
                return [-nup, dm, nup]

            ev = lambda r, y: hof(r, y[0]) - 1e-10
            ev.terminal = True
            r0 = 1e-4
            ec = eos.at(hc)[0]
            a4_0 = math.exp(-4.0 * ALPHA_Y * phi0_val)
            ivp_kwargs = {
                "events": ev,
                "rtol": rtol_ivp,
                "atol": [1e-13, 1e-13, 1e-13] if use_dop853 else [1e-12, 1e-12, 1e-12],
                "dense_output": True,
                "max_step": 0.02 if use_dop853 else 0.05,
            }
            if use_dop853:
                ivp_kwargs["method"] = "DOP853"
            sol = solve_ivp(
                rhs,
                (r0, 40.0),
                [hc, (4.0 * math.pi / 3.0) * (a4_0 * ec + float(rho_phi_fn(0.0))) * (r0 ** 3), 0.0],
                **ivp_kwargs,
            )
            R = float(sol.t_events[0][0])
            M_km = float(sol.y_events[0][0][1])
            rg = np.linspace(0.0, R, ngrid)
            Y = sol.sol(np.maximum(rg, r0))
            hg = np.array([hof(r, y0) for r, y0 in zip(rg, Y[0])], dtype=float)
            hg[0] = hc
            hg[-1] = 0.0
            mg = Y[1].copy()
            mg[0] = 0.0
            er, Pr, _, _ = eos.at(np.maximum(hg, 1e-12))
            er[-1] = 0.0
            Pr[-1] = 0.0
            return {"R": R, "M_Msun": M_km / GM_SUN_KM, "rg": rg, "mg": mg, "hg": hg, "er": er, "Pr": Pr}

        zero_fn = lambda r: 0.0
        cur = integrate_conf_tov(15.0, zero_fn, zero_fn, zero_fn)
        phi = np.zeros(ngrid)
        x_prev = cur["rg"] / cur["rg"][-1]
        rho_p = np.zeros(ngrid)
        pr_p = np.zeros(ngrid)
        k_robin = 1.0 / cur["R"] + 1.0 / L_KM
        history = []
        for it in range(100):
            rg, er, Pr, mg, hg = cur["rg"], cur["er"], cur["Pr"], cur["mg"], cur["hg"]
            phi0 = np.interp(rg / rg[-1], x_prev, phi)
            rho0 = np.interp(rg / rg[-1], x_prev, rho_p)
            pr0 = np.interp(rg / rg[-1], x_prev, pr_p)
            a4 = np.exp(-4.0 * ALPHA_Y * phi0)
            rs = np.maximum(rg, 1e-6)
            B = 1.0 - 2.0 * mg / rs
            B[0] = 1.0
            nup = (mg + 4.0 * math.pi * (rs ** 3) * (a4 * Pr + pr0)) / (rs * (rs - 2.0 * mg))
            nup[0] = 0.0
            dm = 4.0 * math.pi * (rs ** 2) * (a4 * er + rho0)
            dB = 2.0 * mg / (rs ** 2) - 2.0 * dm / rs
            dB[0] = 0.0
            if ((update_every and it >= 1) or (not update_every and it in (1, 3, 6, 10))) and abs(phi0[-1]) > 1e-14:
                if use_4comp_robin:
                    ext_it = cy.solve_exterior_dtN_curved(rg[-1], mg[-1], float(phi0[-1]), L_KM, alpha_Y=ALPHA_Y, tol=1e-12)
                else:
                    ext_it = cy._solve_exterior_3comp_robin_iterate(rg[-1], mg[-1], float(phi0[-1]), L_KM, alpha_Y=ALPHA_Y)
                k_robin = ext_it["k_curved"]
            src_conf = (ALPHA_Y ** 2) * a4 * (er - 3.0 * Pr)
            phin = cy._solve_curved_scalar_with_k(rg, B, dB, nup, src_conf, L_KM, k_robin)
            phi_new = phin if it == 0 else 0.5 * phin + 0.5 * phi0
            dphi = np.gradient(phi_new, rg, edge_order=2)
            dphi[0] = 0.0
            dphi[-1] = -k_robin * phi_new[-1]
            rho_p = (B * (dphi ** 2) + (phi_new ** 2) / (L_KM ** 2)) / (2.0 * ALPHA_Y)
            pr_p = (B * (dphi ** 2) - (phi_new ** 2) / (L_KM ** 2)) / (2.0 * ALPHA_Y)
            x_prev, phi = rg / rg[-1], phi_new
            Rc = rg[-1]
            mk = lambda arr: (lambda I, Rr: (lambda r: float(I(min(r, Rr)))))(PchipInterpolator(rg, arr), Rc)
            new = integrate_conf_tov(15.0, mk(phi), mk(rho_p), mk(pr_p))
            res = max(abs(new["R"] / cur["R"] - 1.0), abs(new["M_Msun"] / cur["M_Msun"] - 1.0), abs(phi_new[0] - phi0[0]))
            history.append({
                "iter": int(it),
                "step_residual": float(res),
                "max_unrelaxed_phi_step": float(np.max(np.abs(phin - phi0))),
                "R_star_km": float(new["R"]),
                "M_surf_Msun": float(new["M_Msun"]),
                "phi_0": float(phi_new[0]),
                "phi_R": float(phi_new[-1]),
                "k_robin_km_inv": float(k_robin),
            })
            cur = new
            if res < tol_iter and it > 3:
                break

        rg, er, Pr, mg, hg = cur["rg"], cur["er"], cur["Pr"], cur["mg"], cur["hg"]
        phi_f = np.interp(rg / rg[-1], x_prev, phi)
        a4 = np.exp(-4.0 * ALPHA_Y * phi_f)
        rs = np.maximum(rg, 1e-6)
        B = 1.0 - 2.0 * mg / rs
        B[0] = 1.0
        dphi_f = np.gradient(phi_f, rg, edge_order=2)
        dphi_f[0] = 0.0
        dphi_f[-1] = -k_robin * phi_f[-1]
        rho_pf = (B * (dphi_f ** 2) + (phi_f ** 2) / (L_KM ** 2)) / (2.0 * ALPHA_Y)
        pr_pf = (B * (dphi_f ** 2) - (phi_f ** 2) / (L_KM ** 2)) / (2.0 * ALPHA_Y)
        nup = (mg + 4.0 * math.pi * (rs ** 3) * (a4 * Pr + pr_pf)) / (rs * (rs - 2.0 * mg))
        nup[0] = 0.0
        dm = 4.0 * math.pi * (rs ** 2) * (a4 * er + rho_pf)
        dB = 2.0 * mg / (rs ** 2) - 2.0 * dm / rs
        dB[0] = 0.0
        src_conf = (ALPHA_Y ** 2) * a4 * (er - 3.0 * Pr)
        phi_exact_grid = cy._solve_curved_scalar_with_k(rg, B, dB, nup, src_conf, L_KM, k_robin)
        h_step = float(rg[1] - rg[0])
        b_coef = 2.0 * B[1:-1] / rg[1:-1] + B[1:-1] * nup[1:-1] + 0.5 * dB[1:-1]
        ddphi_fd = (phi_f[2:] - 2.0 * phi_f[1:-1] + phi_f[:-2]) / (h_step ** 2)
        dphi_fd = (phi_f[2:] - phi_f[:-2]) / (2.0 * h_step)
        ode_res_iter = B[1:-1] * ddphi_fd + b_coef * dphi_fd - phi_f[1:-1] / (L_KM ** 2) + src_conf[1:-1]
        ddphi_ex = (phi_exact_grid[2:] - 2.0 * phi_exact_grid[1:-1] + phi_exact_grid[:-2]) / (h_step ** 2)
        dphi_ex = (phi_exact_grid[2:] - phi_exact_grid[:-2]) / (2.0 * h_step)
        ode_res_ex = B[1:-1] * ddphi_ex + b_coef * dphi_ex - phi_exact_grid[1:-1] / (L_KM ** 2) + src_conf[1:-1]

        aN = float(B[-1])
        bN = float(2.0 * B[-1] / rg[-1] + B[-1] * nup[-1] + 0.5 * dB[-1])
        dphi_robin = float(-k_robin * phi[-1])
        dphi_ghost = float(
            (-(2.0 * aN / (h_step ** 2)) * phi[-2] + (2.0 * aN / (h_step ** 2) + 1.0 / (L_KM ** 2)) * phi[-1])
            / (2.0 * h_step * (aN / (h_step ** 2) + bN / (2.0 * h_step)))
        )
        dphi_fd2 = float((3.0 * phi[-1] - 4.0 * phi[-2] + phi[-3]) / (2.0 * h_step))
        R_star = float(cur["R"])
        M_surf_Msun = float(cur["M_Msun"])
        M_surf_km = M_surf_Msun * GM_SUN_KM
        phi_R = float(phi[-1])
        ext = cy.solve_exterior_dtN_curved(R_star, M_surf_km, phi_R, L_KM, alpha_Y=ALPHA_Y, tol=1e-12)
        bvp_sol = ext["sol"]
        M_ADM_km = float(ext["M_ADM_km"])
        R_max = float(ext["R_max"])
        phi_Rmax = float(ext["phi_ext"][-1])
        B_Rmax = 1.0 - 2.0 * M_ADM_km / R_max
        k_asympt = (1.0 / L_KM) / math.sqrt(B_Rmax)

        def eval_fields(r):
            r = float(r)
            if r <= R_max:
                y_v = bvp_sol.sol(r)
                return float(y_v[2]), float(y_v[3]), float(y_v[0]), float(y_v[1])
            else:
                dr = r - R_max
                phi_v = 0.0 if dr / L_KM > 700.0 else phi_Rmax * (R_max / r) * math.exp(-k_asympt * dr)
                dphi_v = 0.0 if dr / L_KM > 700.0 else -phi_v * (1.0 / r + k_asympt)
                return M_ADM_km, 0.5 * math.log1p(-2.0 * M_ADM_km / r), phi_v, dphi_v

        b_graz = max(R_star * math.exp(-float(ext["nu_R"])), R_star * math.exp(-(float(ext["nu_R"]) - phi_R))) * 1.0001
        alpha_graz, r0_graz = compute_exact_deflection_angle(b_graz, eval_fields, R_star, R_max, M_ADM_km, branch="conformal")
        return {
            "N_grid": int(ngrid),
            "h_step_km": h_step,
            "iterations": len(history),
            "R_star_km": R_star,
            "M_surf_Msun": M_surf_Msun,
            "M_ADM_Msun": M_ADM_km / GM_SUN_KM,
            "phi_0": float(phi[0]),
            "phi_R": phi_R,
            "dphi_int_robin_km_inv": dphi_robin,
            "dphi_int_ghost_km_inv": dphi_ghost,
            "dphi_int_fd2_km_inv": dphi_fd2,
            "dphi_ext_curved_km_inv": float(ext["dphi_R_curved"]),
            "abs_diff_dphi_robin_vs_ext_km_inv": abs(dphi_robin - float(ext["dphi_R_curved"])),
            "rel_diff_dphi_robin_vs_ext": abs(dphi_robin - float(ext["dphi_R_curved"])) / abs(float(ext["dphi_R_curved"])),
            "abs_diff_dphi_ghost_vs_ext_km_inv": abs(dphi_ghost - float(ext["dphi_R_curved"])),
            "abs_diff_dphi_fd2_vs_ext_km_inv": abs(dphi_fd2 - float(ext["dphi_R_curved"])),
            "k_int_km_inv": float(k_robin),
            "k_ext_km_inv": float(ext["k_curved"]),
            "abs_diff_k_int_vs_ext_km_inv": abs(float(k_robin) - float(ext["k_curved"])),
            "rel_diff_k_int_vs_ext": abs(float(k_robin) - float(ext["k_curved"])) / float(ext["k_curved"]),
            "max_interior_scalar_ode_residual_iter_km_inv2": float(np.max(np.abs(ode_res_iter))),
            "max_interior_scalar_ode_residual_banded_km_inv2": float(np.max(np.abs(ode_res_ex))),
            "max_exterior_4comp_bvp_rms_residual": float(np.max(bvp_sol.rms_residuals)),
            "b_graz_km": b_graz,
            "r0_graz_km": r0_graz,
            "alpha_hat_Yukawa_conf_graz_rad": alpha_graz,
            "history": history,
        }

    baseline_601 = run_variant(601, update_every=False, use_4comp_robin=False, tol_iter=2e-8, rtol_ivp=1e-9, use_dop853=False)
    every_step_4comp_601 = run_variant(601, update_every=True, use_4comp_robin=True, tol_iter=1e-11, rtol_ivp=1e-11, use_dop853=False)
    mesh_601 = run_variant(601, update_every=True, use_4comp_robin=False, tol_iter=1e-10, rtol_ivp=1e-11, use_dop853=True)
    mesh_2000 = run_variant(2000, update_every=True, use_4comp_robin=False, tol_iter=1e-10, rtol_ivp=1e-11, use_dop853=True)
    mesh_4000 = run_variant(4000, update_every=True, use_4comp_robin=False, tol_iter=1e-10, rtol_ivp=1e-11, use_dop853=True)

    # Strip full history from mesh_2000 and mesh_4000 to keep JSON concise while retaining full history for baseline_601 and every_step_4comp_601
    mesh_601_summary = {k: v for k, v in mesh_601.items() if k != "history"}
    mesh_2000_summary = {k: v for k, v in mesh_2000.items() if k != "history"}
    mesh_4000_summary = {k: v for k, v in mesh_4000.items() if k != "history"}

    return {
        "baseline_N601_4step_robin_schedule": baseline_601,
        "every_step_4comp_robin_N601_tol1e11": every_step_4comp_601,
        "mesh_refinement_DOP853_every_step_robin": {
            "N_601": mesh_601_summary,
            "N_2000": mesh_2000_summary,
            "N_4000": mesh_4000_summary,
            "convergence_diffs_N2000_vs_N4000": {
                "abs_diff_R_star_km": abs(mesh_4000["R_star_km"] - mesh_2000["R_star_km"]),
                "rel_diff_R_star": abs(mesh_4000["R_star_km"] - mesh_2000["R_star_km"]) / mesh_4000["R_star_km"],
                "abs_diff_M_surf_Msun": abs(mesh_4000["M_surf_Msun"] - mesh_2000["M_surf_Msun"]),
                "rel_diff_M_surf": abs(mesh_4000["M_surf_Msun"] - mesh_2000["M_surf_Msun"]) / mesh_4000["M_surf_Msun"],
                "abs_diff_M_ADM_Msun": abs(mesh_4000["M_ADM_Msun"] - mesh_2000["M_ADM_Msun"]),
                "rel_diff_M_ADM": abs(mesh_4000["M_ADM_Msun"] - mesh_2000["M_ADM_Msun"]) / mesh_4000["M_ADM_Msun"],
                "abs_diff_phi_0": abs(mesh_4000["phi_0"] - mesh_2000["phi_0"]),
                "rel_diff_phi_0": abs(mesh_4000["phi_0"] - mesh_2000["phi_0"]) / mesh_4000["phi_0"],
                "abs_diff_phi_R": abs(mesh_4000["phi_R"] - mesh_2000["phi_R"]),
                "rel_diff_phi_R": abs(mesh_4000["phi_R"] - mesh_2000["phi_R"]) / mesh_4000["phi_R"],
                "abs_diff_dphi_ext_curved_km_inv": abs(mesh_4000["dphi_ext_curved_km_inv"] - mesh_2000["dphi_ext_curved_km_inv"]),
                "rel_diff_dphi_ext_curved": abs(mesh_4000["dphi_ext_curved_km_inv"] - mesh_2000["dphi_ext_curved_km_inv"]) / abs(mesh_4000["dphi_ext_curved_km_inv"]),
                "abs_diff_alpha_hat_graz_rad": abs(mesh_4000["alpha_hat_Yukawa_conf_graz_rad"] - mesh_2000["alpha_hat_Yukawa_conf_graz_rad"]),
                "rel_diff_alpha_hat_graz": abs(mesh_4000["alpha_hat_Yukawa_conf_graz_rad"] - mesh_2000["alpha_hat_Yukawa_conf_graz_rad"]) / mesh_4000["alpha_hat_Yukawa_conf_graz_rad"],
                "fd2_boundary_slope_error_ratio_N2000_over_N4000": mesh_2000["abs_diff_dphi_fd2_vs_ext_km_inv"] / mesh_4000["abs_diff_dphi_fd2_vs_ext_km_inv"],
                "theoretical_second_order_ratio_N2000_over_N4000": (mesh_2000["h_step_km"] / mesh_4000["h_step_km"]) ** 2,
            },
        },
    }


def compute_branch_fluid_vs_conf_interior_comparison(ns_fluid, ns_conf):
    """
    Compare the Step 22 Branch II-Fluid background (L_m = -exp(-phi_L) eps(n), scalar source
    -alpha_Y exp(-phi_L) eps) against the self-consistent Step 23 Branch II-Conf trace-sourced
    background (S_m[exp(-2 alpha_Y phi_L) g, Psi], scalar source -alpha_Y^2 exp(-4 alpha_Y phi_L)(eps - 3P))
    across interior BVP quantities and curved-spacetime deflection angles.
    """
    def eval_branch_lensing_sample(ns_obj):
        R_s = ns_obj["R_star_km"]
        R_m = ns_obj["R_max_km"]
        M_s_km = ns_obj["M_surf_km"]
        M_a_km = ns_obj["M_ADM_km"]
        ev_y = ns_obj["eval_fields"]
        ev_adm = lambda r: (M_a_km, 0.5 * math.log1p(-2.0 * M_a_km / float(r)), 0.0, 0.0)
        ev_srf = lambda r: (M_s_km, 0.5 * math.log1p(-2.0 * M_s_km / float(r)), 0.0, 0.0)
        b_graz = max(
            R_s * math.exp(-ns_obj["nu_R_exact"]),
            R_s * math.exp(-(ns_obj["nu_R_exact"] - ns_obj["phi_R"])),
        ) * 1.0001
        b_shared_conf_graz = max(
            ns_conf["R_star_km"] * math.exp(-ns_conf["nu_R_exact"]),
            ns_conf["R_star_km"] * math.exp(-(ns_conf["nu_R_exact"] - ns_conf["phi_R"])),
        ) * 1.0001
        sample_b = [
            ("surface_grazing", b_graz),
            ("b_15p064914km_shared_conf_grazing", b_shared_conf_graz),
            ("b_20km", 20.0),
            ("b_50km", 50.0),
            ("b_200km", 200.0),
            ("b_1Rsun", R_SUN_KM),
        ]
        out = {}
        for lbl, b_km in sample_b:
            a_conf, r0 = compute_exact_deflection_angle(b_km, ev_y, R_s, R_m, M_a_km, branch="conformal")
            a_adm, _ = compute_exact_deflection_angle(b_km, ev_adm, R_s, 0.0, M_a_km, branch="conformal")
            a_srf, _ = compute_exact_deflection_angle(b_km, ev_srf, R_s, 0.0, M_s_km, branch="conformal")
            out[lbl] = {
                "b_km": b_km,
                "r0_km": r0,
                "alpha_hat_Yukawa_conf_rad": a_conf,
                "alpha_hat_Schw_MADM_rad": a_adm,
                "alpha_hat_Schw_Msurf_rad": a_srf,
                "rel_diff_vs_Schw_MADM": (a_conf - a_adm) / a_adm,
                "rel_diff_vs_Schw_Msurf": (a_conf - a_srf) / a_srf,
            }
        return out

    fluid_sample = eval_branch_lensing_sample(ns_fluid)
    conf_sample = eval_branch_lensing_sample(ns_conf)
    phi0_fluid = float(ns_fluid["star"]["profile"]["phi"][0])
    return {
        "EOS": "SLy",
        "log10_rho0c": 15.0,
        "L_km": L_KM,
        "alpha_Y": ALPHA_Y,
        "exterior_vacuum_ode_identical": True,
        "interior_weak_field_dust_limit_identical_at_alphaY_eq_1": True,
        "interior_relativistic_bvp_identical": False,
        "reason_for_interior_difference": (
            "Branch II-Fluid (Step 22) uses L_m = -exp(-phi_L) eps(n), giving Einstein-frame matter density "
            "exp(-phi_L) eps, pressure exp(-phi_L) P, scalar source -alpha_Y exp(-phi_L) eps, and Euler-TOV "
            "h' = -nu' + phi_L'. Branch II-Conf (Step 23) uses S_m[exp(-2 alpha_Y phi_L) g, Psi], giving "
            "Einstein-frame matter density exp(-4 alpha_Y phi_L) eps, pressure exp(-4 alpha_Y phi_L) P, "
            "scalar trace source -alpha_Y^2 exp(-4 alpha_Y phi_L)(eps - 3P), and Euler-TOV h' = -nu' + alpha_Y phi_L'. "
            f"At the center of the SLy log10(rho_0c)=15.0 star, 3P_c/eps_c = {ns_conf['center_3P_over_eps']:.6f}, so 1 - 3P_c/eps_c = {ns_conf['center_trace_factor_1_minus_3P_over_eps']:.6f}."
        ),
        "Branch_II_Fluid_Step22": {
            "matter_action": "L_m = -exp(-phi_L) eps(n)",
            "scalar_source_rhs": "-(alpha_Y G / c^2) exp(-phi_L) eps",
            "euler_tov": "h' = -nu' + phi_L'",
            "R_star_km": ns_fluid["R_star_km"],
            "M_surf_Msun": ns_fluid["M_surf_Msun"],
            "M_ADM_Msun": ns_fluid["M_ADM_Msun"],
            "Delta_M_phi_ext_Msun": ns_fluid["Delta_M_phi_ext_Msun"],
            "Delta_M_phi_ext_rel_Msurf": ns_fluid["Delta_M_phi_ext_Msun"] / ns_fluid["M_surf_Msun"],
            "phi_L_center": phi0_fluid,
            "phi_L_R": ns_fluid["phi_R"],
            "dphi_L_dr_R_curved_km_inv": ns_fluid["dphi_R_curved"],
            "k_curved_km_inv": ns_fluid["k_curved_km_inv"],
            "nu_R_exact": ns_fluid["nu_R_exact"],
            "lensing_samples": fluid_sample,
        },
        "Branch_II_Conf_Step23_trace": {
            "matter_action": "S_m[exp(-2 alpha_Y phi_L) g_{mu nu}, Psi_m]",
            "scalar_source_rhs": "-(alpha_Y^2 G / c^2) exp(-4 alpha_Y phi_L) (eps - 3 P)",
            "euler_tov": "h' = -nu' + alpha_Y phi_L'",
            "center_3P_over_eps": ns_conf["center_3P_over_eps"],
            "center_trace_factor_1_minus_3P_over_eps": ns_conf["center_trace_factor_1_minus_3P_over_eps"],
            "R_star_km": ns_conf["R_star_km"],
            "M_surf_Msun": ns_conf["M_surf_Msun"],
            "M_ADM_Msun": ns_conf["M_ADM_Msun"],
            "Delta_M_phi_ext_Msun": ns_conf["Delta_M_phi_ext_Msun"],
            "Delta_M_phi_ext_rel_Msurf": ns_conf["Delta_M_phi_ext_rel_Msurf"],
            "phi_L_center": ns_conf["phi_0"],
            "phi_L_R": ns_conf["phi_R"],
            "dphi_L_dr_R_int_robin_km_inv": ns_conf["dphi_R_int_robin"],
            "dphi_L_dr_R_int_ghost_km_inv": ns_conf["dphi_R_int_ghost"],
            "dphi_L_dr_R_int_fd2_km_inv": ns_conf["dphi_R_int_fd2"],
            "dphi_L_dr_R_curved_km_inv": ns_conf["dphi_R_curved"],
            "k_robin_int_km_inv": ns_conf["k_robin_int_km_inv"],
            "k_curved_km_inv": ns_conf["k_curved_km_inv"],
            "abs_diff_k_int_vs_k_curved": ns_conf["abs_diff_k_int_vs_k_curved"],
            "rel_diff_k_int_vs_k_curved": ns_conf["rel_diff_k_int_vs_k_curved"],
            "abs_diff_dphi_int_robin_vs_ext_curved": ns_conf["abs_diff_dphi_int_robin_vs_ext_curved"],
            "rel_diff_dphi_int_robin_vs_ext_curved": ns_conf["rel_diff_dphi_int_robin_vs_ext_curved"],
            "abs_diff_dphi_int_fd2_vs_ext_curved": ns_conf["abs_diff_dphi_int_fd2_vs_ext_curved"],
            "max_interior_scalar_ode_residual_iter_km_inv2": ns_conf["max_interior_scalar_ode_residual_iter_km_inv2"],
            "max_interior_scalar_ode_residual_banded_km_inv2": ns_conf["max_interior_scalar_ode_residual_banded_km_inv2"],
            "max_exterior_4comp_bvp_rms_residual": ns_conf["max_exterior_4comp_bvp_rms_residual"],
            "nu_R_exact": ns_conf["nu_R_exact"],
            "lensing_samples": conf_sample,
            "bvp_convergence_and_mesh_refinement": compute_branch_ii_conf_bvp_convergence_and_mesh_refinement(),
        },
    }


def main():
    sym_checks = run_symbolic_horndeski_and_lensing_checks()
    df_gw, gw_meta = compute_gw170817_dispersion_table()
    ns_ref = build_v25_reference_ns_exterior()
    ns_conf_ref = build_branch_ii_conf_trace_ns_exterior()
    branch_comparison = compute_branch_fluid_vs_conf_interior_comparison(ns_ref, ns_conf_ref)
    df_lens_ns = compute_lensing_and_ppn_table(ns_ref)
    df_sun = compute_solar_cassini_1msun_table()
    df_lens = pd.concat([df_lens_ns, df_sun], ignore_index=True, sort=False)
    macro_bridge = compute_macroscopic_lensing_and_bullet_cluster_bridge()

    csv_gw_path = os.path.join(OUT_DIR, "gw170817_dispersion_and_cT.csv")
    csv_lens_path = os.path.join(OUT_DIR, "lensing_deflection_and_ppn_slip.csv")
    json_summary_path = os.path.join(OUT_DIR, "gw170817_and_lensing_summary.json")

    df_gw.to_csv(csv_gw_path, index=False)
    df_lens.to_csv(csv_lens_path, index=False)

    summary_payload = {
        "module": "Step 23 (V26) — Module 1: GW170817 (c_T = 1) and Gravitational Lensing at L = 10 km, alpha_Y = 1",
        "author": "David Szabolcsi (ORCID: 0009-0005-6300-074X)",
        "baseline_parameters": {
            "L_km": L_KM,
            "L_m": L_M,
            "alpha_Y": ALPHA_Y,
            "alpha_std": ALPHA_STD,
            "beta_Yuk": BETA_YUK,
            "convention": "A(phi_L) = exp(-alpha_Y phi_L), phi_L > 0, g_J = exp(-2 alpha_Y phi_L) g",
            "tilde_alpha_L": ALPHA_TILDE_L,
            "alpha_DEF_dlnA_dphi_canon": ALPHA_DEF,
            "beta_DEF_d2lnA_dphi_canon2": BETA_DEF_SECOND_DERIV,
            "beta_DEF": BETA_DEF,
            "omega_BD": OMEGA_BD,
            "omega_BD_legacy_step5_4pi_conv": OMEGA_BD_LEGACY_STEP5,
            "phi_canon_over_phi_L": PHI_CANON_OVER_PHI_L,
            "dlnA_dphi_canon": DLNA_DPHI_CANON,
            "dlnA_dphi_L": -ALPHA_Y,
            "gamma_L_m3_kg_inv": GAMMA_L_SI,
        },
        "symbolic_horndeski_and_lensing_checks": sym_checks,
        "gw170817_and_scalar_dispersion_meta": gw_meta,
        "branch_ii_fluid_vs_branch_ii_conf_interior_bvp_comparison": branch_comparison,
        "v25_reference_neutron_star_SLy_15p0": {
            "background_branch_label": "Branch_II_Fluid_Step22_diagnostic_background",
            "EOS": "SLy",
            "log10_rho0c": 15.0,
            "R_star_km": ns_ref["R_star_km"],
            "M_surf_Msun": ns_ref["M_surf_Msun"],
            "M_ADM_Msun": ns_ref["M_ADM_Msun"],
            "Delta_M_phi_ext_Msun": ns_ref["Delta_M_phi_ext_Msun"],
            "M_eff_flat_Msun": ns_ref["M_eff_flat_Msun"],
            "phi_L_R": ns_ref["phi_R"],
            "dphi_L_dr_R_curved_km_inv": ns_ref["dphi_R_curved"],
            "k_curved_km_inv": ns_ref["k_curved_km_inv"],
            "nu_R_exact": ns_ref["nu_R_exact"],
            "nu_R_Schw_Msurf": ns_ref["nu_R_Schw_Msurf"],
            "nu_R_Schw_MADM": ns_ref["nu_R_Schw_MADM"],
        },
        "gw170817_dispersion_rows": df_gw.to_dict(orient="records"),
        "lensing_and_ppn_rows": df_lens_ns.to_dict(orient="records"),
        "solar_cassini_1Msun_benchmark": {
            "M_Msun": 1.0, "GM_sun_over_c2_km": GM_SUN_KM, "R_sun_km": R_SUN_KM,
            "r_E_AU": 1.0, "r_S_AU": 8.43,
            "rows": df_sun.to_dict(orient="records"),
        },
        "macroscopic_lensing_and_bullet_cluster_bridge": macro_bridge,
    }

    with open(json_summary_path, "w", encoding="utf-8") as f:
        json.dump(summary_payload, f, indent=2, allow_nan=True)

    for d in [CODE_DIR, CSV_DIR, DATA_V23_DIR]:
        os.makedirs(d, exist_ok=True)
    shutil.copy2(__file__, os.path.join(CODE_DIR, os.path.basename(__file__)))
    for p in [csv_gw_path, csv_lens_path, json_summary_path]:
        for target_dir in [CODE_DIR, CSV_DIR, DATA_V23_DIR]:
            shutil.copy2(p, os.path.join(target_dir, os.path.basename(p)))

    bvp_ver = branch_comparison["Branch_II_Conf_Step23_trace"]["bvp_convergence_and_mesh_refinement"]
    m2000 = bvp_ver["mesh_refinement_DOP853_every_step_robin"]["N_2000"]
    m4000 = bvp_ver["mesh_refinement_DOP853_every_step_robin"]["N_4000"]
    mdiff = bvp_ver["mesh_refinement_DOP853_every_step_robin"]["convergence_diffs_N2000_vs_N4000"]
    ev4 = bvp_ver["every_step_4comp_robin_N601_tol1e11"]

    print("=== MODULE 1 (STEP 23) AUDIT SUMMARY ===")
    print(f"f_cut = {gw_meta['f_cut_Hz']:.6f} Hz, omega_cut = {gw_meta['omega_cut_rad_s']:.4f} rad/s, m_phi = {gw_meta['m_phi_eV']:.7e} eV")
    print(f"V25 SLy 15.0 NS (Branch II-Fluid): R_* = {ns_ref['R_star_km']:.8f} km, M_surf = {ns_ref['M_surf_Msun']:.8f} M_sun, M_ADM = {ns_ref['M_ADM_Msun']:.8f} M_sun")
    print(f"Surface nu_exact = {ns_ref['nu_R_exact']:.10f}, nu_Schw(M_surf) = {ns_ref['nu_R_Schw_Msurf']:.10f}, nu_Schw(M_ADM) = {ns_ref['nu_R_Schw_MADM']:.10f}")
    print("--- Step 23 Branch II-Conf SLy 15.0 BVP Verification ---")
    print(f"Branch II-Conf (N=601 baseline): R_* = {ns_conf_ref['R_star_km']:.8f} km, M_surf = {ns_conf_ref['M_surf_Msun']:.8f} M_sun, M_ADM = {ns_conf_ref['M_ADM_Msun']:.8f} M_sun, phi_L(0) = {ns_conf_ref['phi_0']:.8f}, phi_L(R_*) = {ns_conf_ref['phi_R']:.8f}")
    print(f"Interior ODE residual max|res_ODE|: iter = {ns_conf_ref['max_interior_scalar_ode_residual_iter_km_inv2']:.3e} km^-2, banded = {ns_conf_ref['max_interior_scalar_ode_residual_banded_km_inv2']:.3e} km^-2; Exterior 4-comp BVP rms_res = {ns_conf_ref['max_exterior_4comp_bvp_rms_residual']:.3e}")
    print(f"Surface slope matching (N=601, 4-step Robin): dphi_int_robin = {ns_conf_ref['dphi_R_int_robin']:.12f} km^-1, dphi_ext_curved = {ns_conf_ref['dphi_R_curved']:.12f} km^-1, |diff| = {ns_conf_ref['abs_diff_dphi_int_robin_vs_ext_curved']:.3e} km^-1 (|Delta k| = {ns_conf_ref['abs_diff_k_int_vs_k_curved']:.3e} km^-1)")
    print(f"Surface slope matching (N=601, every-step 4-comp Robin, tol=1e-11): dphi_int_robin = {ev4['dphi_int_robin_km_inv']:.12f} km^-1, dphi_ext_curved = {ev4['dphi_ext_curved_km_inv']:.12f} km^-1, |diff| = {ev4['abs_diff_dphi_robin_vs_ext_km_inv']:.3e} km^-1 (|Delta k| = {ev4['abs_diff_k_int_vs_ext_km_inv']:.3e} km^-1), max|res_ODE| = {ev4['max_interior_scalar_ode_residual_iter_km_inv2']:.3e} km^-2")
    print(f"Mesh refinement N=2000 vs N=4000: |Delta R_*| = {mdiff['abs_diff_R_star_km']:.3e} km (rel {mdiff['rel_diff_R_star']:.3e}), |Delta M_ADM| = {mdiff['abs_diff_M_ADM_Msun']:.3e} M_sun (rel {mdiff['rel_diff_M_ADM']:.3e}), |Delta phi_L(R_*)| = {mdiff['abs_diff_phi_R']:.3e} (rel {mdiff['rel_diff_phi_R']:.3e}), FD2 slope error ratio = {mdiff['fd2_boundary_slope_error_ratio_N2000_over_N4000']:.4f} (theory {mdiff['theoretical_second_order_ratio_N2000_over_N4000']:.4f})")


if __name__ == "__main__":
    main()
