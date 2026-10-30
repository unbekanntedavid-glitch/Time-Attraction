#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 27 (part of V28, Time Attraction Law - Scalar-Tensor Yukawa Fit): measurable local consequence of the FIXED Yukawa branch of the
TAL scalar-tensor Yukawa fit (L = 10 km, alpha_Y = 1, beta_Yuk = 1/(4 pi), C = 0
attractor), tested against primary geophysical measurements near lambda = 10 km.

Convention used everywhere (standard Yukawa parametrization):
    V(r) = - G_inf m1 m2 / r * (1 + alpha * exp(-r/lambda))
G_inf is the long-range constant (G_* of the package); a laboratory experiment
at r << lambda measures G_lab = G_inf (1 + alpha).

Model point (from the package, see MODEL_SOURCES):
    alpha = q = alpha_Y^2 * beta_Yuk = 1/(4 pi) = 0.0795775, lambda = L = 10 km,
    positive (attractive), universal (all matter couples to A^2(phi) g), no screening.

The script only uses numbers stated in the cited primary papers (or quoted from
them in a cited secondary paper) plus explicitly labelled assumptions.
Ratio to bound: R = alpha_model / alpha_lim(10 km). R > 1 means the model point
lies outside the allowed region under the same convention and assumptions; it is
NOT a significance in sigma.

REVIEWER GUIDE: residual -> alpha chain (see step27_rederivation_method_en.md)
  Ocean (Zumberge et al. 1991 data, secondary-source values from March et al. 2019):
    1. yukawa_slab_field_infinite_layers()     volume integral of the Yukawa field of horizontal layers
    2. ocean_model_residual_slope()             model residual Delta g_m(z) and its fitted slope, same
                                                definition as the measured residual (March eq. 73)
    3. ocean_alpha_lim_from_residual()          alpha where the model slope reaches the edge of the
                                                measured slope interval (k_sigma * quoted uncertainty);
                                                exact inversion of s(alpha) = S alpha/(1+alpha)
    4. march_eq81_F_difference()                cross-check against March eq. (81) as printed
  Towers (Thomas 1989, Romaides 1994/1997, BREN report OSTI 6190822):
    5. tower_model_residual_uGal()              half-space Yukawa residual obs - continued prediction
    6. tower_alpha_lim_from_residual()          alpha where the model residual reaches mu - k_sigma*sigma
                                                (exact: alpha_lim = E/(K - E), G_inf = G_lab/(1+alpha))
    back_substitution_tests()                   each limit point must return the allowed residual edge
  Erie tower (Cruz et al. 1991):
    7. cruz_conditional_scaling_factor()        CONDITIONAL transfer 1 km -> 10 km (assumption-laden)
"""
import json, math, os, hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "step27_results.json")

# ---------------------------------------------------------------- model point
ALPHA_Y = 1.0
BETA_YUK = 1.0 / (4.0 * math.pi)
Q = ALPHA_Y ** 2 * BETA_YUK                 # 0.0795774715
LAMBDA_M = 10.0e3
G_LAB = 6.67430e-11                         # CODATA 2018 value used by the package
G_INF = G_LAB / (1.0 + Q)                   # long-range G_* AT THE MODEL POINT q only (G_lab fixed by calibration)
TWO_PI_G = 2.0 * math.pi * G_INF            # used only where alpha = q

def g_inf_of(alpha):
    """Calibration: the measured G_lab is held fixed, so for ANY trial alpha the long-range
    constant is G_inf(alpha) = G_lab/(1+alpha) (review r1, S01)."""
    return G_LAB / (1.0 + alpha)
Z95 = 1.959963984540054

MODEL_SOURCES = [
    {"file": "hubble_desi_bao_test/tal_cosmology_common.py", "lines": "43-49",
     "content": "L_KM = 10.0; ALPHA_Y = 1.0; BETA_YUK = 1/(4 pi)"},
    {"file": "hubble_desi_bao_test/tal_cosmology_common.py", "lines": "50-64",
     "content": "action with kinetic term c^4/(4 pi G beta_Yuk), matter coupled to A^2(phi) g, "
                "A = exp(alpha_Y phi); point source phi = -alpha_Y beta_Yuk G M exp(-r/L)/(c^2 r). "
                "Test-mass potential c^2 alpha_Y phi = -alpha_Y^2 beta_Yuk G M exp(-r/L)/r, i.e. alpha = q"},
    {"file": "hubble_desi_bao_test/step26_md_templates.py", "lines": "103",
     "content": "G_lab: for r << L and unscreened small bodies Phi = -G_* M (1 + q e^{-r/L})/r, "
                "G_lab ~ G_*(1+q), q = alpha_Y^2 beta_Yuk"},
    {"file": "hubble_desi_bao_test/step26_md_templates.py", "lines": "106",
     "content": "... no screening branch is introduced."},
    {"file": "hubble_desi_bao_test/run_step26_revision_controls.py", "lines": "244-253",
     "content": "G_eff/G_lab = 1/(1+q) = 0.9262883; 'no screening branch is introduced'"},
    {"file": "scalar-tensor-yukawa-fit/step5_coupling_dictionary_ppn_screening_en.md", "lines": "54, 73, 101",
     "content": "alpha_std = alpha_Y/(4 pi) = 2 beta^2; G_eff(r) = G[1 + alpha_std (1 + r/L) e^{-r/L}]; "
                "r << L: G_eff -> G(1 + alpha_std). ('screened' at line 100 only means exp(-r/L) decay)"},
    {"file": "scalar-tensor-yukawa-fit/step11_chameleon_screening_multiscale_eft_synthesis_en.md", "lines": "18",
     "content": "chameleon/symmetron screening is a standalone nonlinear extension hypothesis beyond the "
                "quadratic constant-alpha_Y Yukawa action, i.e. not part of the fixed branch"},
]

# ------------------------------------------------------------------ helpers
def yukawa_slab_field_infinite_layers(z, layers, lam, alpha, G=None):
    """Downward Yukawa acceleration at depth z (m, positive down) from horizontal
    layers [(top, bottom, rho)], bottom may be inf. Infinite-slab volume integral:
    g_Y = 2 pi G_inf alpha * int rho(z') sign(z'-z) exp(-|z'-z|/lam) dz'."""
    s = 0.0
    for a, b, rho in layers:
        if z <= a:
            tb = 0.0 if math.isinf(b) else math.exp(-(b - z) / lam)
            s += rho * lam * (math.exp(-(a - z) / lam) - tb)
        elif (not math.isinf(b)) and z >= b:
            s -= rho * lam * (math.exp(-(z - b) / lam) - math.exp(-(z - a) / lam))
        else:
            below = rho * lam * (1.0 - (0.0 if math.isinf(b) else math.exp(-(b - z) / lam)))
            above = rho * lam * (1.0 - math.exp(-(z - a) / lam))
            s += below - above
    G = g_inf_of(alpha) if G is None else G
    return 2.0 * math.pi * G * alpha * s

slab_yukawa_g = yukawa_slab_field_infinite_layers   # backward-compatible alias

# ======================================================= 1. ocean experiment
# Zumberge et al. 1991 PRL 67, 3051: G = (6.677 +/- 0.013)e-11 (abstract).
# Inputs quoted from March et al. 2019 PRD 100, 042002 (arXiv:1904.12789), which
# re-analysed the Zumberge data: eq.(46) layers, eq.(74) slope, eq.(79) densities,
# profiles z_N in [500, 4800] m, G = G_N/(1+alpha Phi) with Phi ~ 1 for lambda >= 10 m (eq. 65).
OCEAN = {
    "z_w_m": 5000.0, "z_c_minus_z_w_m": 7000.0,
    "rho_w_mean": 1037.05, "rho_c": 2860.0, "rho_m": 3400.0,
    "z_profile_m": [500.0, 4800.0],
    "slope_obs_mGal_per_km": -0.060, "slope_sigma_mGal_per_km": 0.172,
}

def ocean_model_residual_slope(alpha, lam):
    """STEP 2 (ocean). Model residual with the experiment's own definition
    (March et al. 2019, eq. 73; model form eqs. 63-64):
        Delta g(z) = [g(z) - g(z1)] - beta_N (z - z1),
    beta_N = reference gradient computed with the LAB constant G_N (= G_LAB here).
    In the Yukawa world the true long-range constant is G_inf = G_N/(1+alpha) (March eq. 65,
    Phi ~ 1 for lambda >= 10 m), so
        Delta g_m(z) = -4 pi (G_inf - G_N) rho_w (z - z1)                       (Bouguer-type term)
                       + 2 pi G_inf alpha [F(z) - F(z1)]                         (Yukawa, volume integral)
    The free-air part (G M_earth from orbits) is identical in data and reference and cancels.
    The slope is a least-squares line over z in [500, 4800] m (the measured profile range).
    Returns (slope in mGal/km, z grid, residual in m/s^2)."""
    o = OCEAN
    zc = o["z_w_m"] + o["z_c_minus_z_w_m"]
    layers = [(0.0, o["z_w_m"], o["rho_w_mean"]), (o["z_w_m"], zc, o["rho_c"]), (zc, math.inf, o["rho_m"])]
    g_inf = g_inf_of(alpha)
    z = np.linspace(o["z_profile_m"][0], o["z_profile_m"][1], 431)
    # Newtonian Bouguer part difference: -4 pi (G_inf - G_lab) rho_w z ; free-air term identical (GM from orbits)
    newton_diff = -4.0 * math.pi * (g_inf - G_LAB) * o["rho_w_mean"] * (z - z[0])
    yuk = np.array([slab_yukawa_g(zz, layers, lam, alpha, G=g_inf) for zz in z])
    yuk = yuk - yuk[0]
    resid = newton_diff + yuk                         # m/s^2
    slope = np.polyfit(z, resid, 1)[0]                # s^-2
    return slope * 1e5 * 1e3, z, resid                # mGal/km

def ocean_sensitivity_S(lam):
    """alpha -> 0 sensitivity S (mGal/km per unit alpha), with G_lab fixed. Both terms of the model
    residual scale exactly as alpha/(1+alpha), so S = s(alpha) (1+alpha)/alpha for any alpha."""
    a = 1e-6
    s1, _, _ = ocean_residual_slope(a, lam)
    return s1 * (1.0 + a) / a

def ocean_alpha_lim_from_residual(lam, k_sigma=Z95, linear=False):
    """STEP 3 (ocean). With G_lab fixed, s(alpha) = S alpha/(1+alpha) exactly, so the edge
    edge = slope_obs + k_sigma * sigma_quoted (S > 0 for alpha > 0) gives the EXACT inversion
        alpha_lim = edge / (S - edge)            (primary, review r1 S01)
    and the linear approximation alpha_lim ~ edge / S (reported separately, linear=True).
    k_sigma = 1.96 with sigma read as 1 sigma Gaussian: upper end of the central 95% interval,
    i.e. a one-sided 97.5% upper limit on alpha > 0 (ASSUMPTION; the meaning of +/-0.172 in the
    original paper was not verified). k_sigma = 1 reproduces March et al.'s interval, eqs. (75)-(76)."""
    sens = ocean_sensitivity_S(lam)                   # mGal/km per unit alpha (alpha -> 0)
    o = OCEAN
    lo = o["slope_obs_mGal_per_km"] - k_sigma * o["slope_sigma_mGal_per_km"]
    hi = o["slope_obs_mGal_per_km"] + k_sigma * o["slope_sigma_mGal_per_km"]
    edge = hi if sens > 0 else lo
    if linear:
        return edge / sens, sens
    return edge / (sens - edge), sens

ocean_residual_slope = ocean_model_residual_slope     # backward-compatible aliases
ocean_alpha_lim = ocean_alpha_lim_from_residual

def march_eq81_F_difference(lam):
    """STEP 4 (cross-check). March et al. 2019 eq. (81) as printed (p. 12):
        F(z_w) - F(0) = lam [ rho_w (e^{-z_w/lam} - 1) - (rho_m - rho_c)(1 - e^{z_w/lam})(e^{-z_w/lam} + e^{-z_c/lam}) ]
    compared with our own layer integral F(z_w) - F(0) (same layers, same lam)."""
    o = OCEAN
    zw = o["z_w_m"]; zc = zw + o["z_c_minus_z_w_m"]
    rw, rc, rm = o["rho_w_mean"], o["rho_c"], o["rho_m"]
    march = lam * (rw * (math.exp(-zw / lam) - 1.0) - (rm - rc) * (1.0 - math.exp(zw / lam)) * (math.exp(-zw / lam) + math.exp(-zc / lam)))
    layers = [(0.0, zw, rw), (zw, zc, rc), (zc, math.inf, rm)]
    ours = (yukawa_slab_field_infinite_layers(zw, layers, lam, 1.0, G=1.0) - yukawa_slab_field_infinite_layers(0.0, layers, lam, 1.0, G=1.0)) / (2.0 * math.pi)
    # independent brute-force numerical integral of rho(z') sign(z'-z) exp(-|z'-z|/lam)
    zp = np.linspace(0.0, 40.0 * lam + zc, 2000001)
    rho = np.where(zp < zw, rw, np.where(zp < zc, rc, rm))
    def F(z):
        f = rho * np.sign(zp - z) * np.exp(-np.abs(zp - z) / lam)
        return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(zp)))
    brute = F(zw) - F(0.0)
    # sensitivity (mGal/km per unit alpha, small alpha) if March eq. 81 as printed were used instead
    bouguer = 4.0 * math.pi * G_LAB * rw * 1e8                 # mGal/km per unit alpha (small alpha)
    # alpha -> 0 sensitivities: G_inf -> G_lab (review r1 S01); chord values over 0..z_w
    yuk_ours = 2.0 * math.pi * G_LAB * ours / zw * 1e8
    yuk_march = 2.0 * math.pi * G_LAB * march / zw * 1e8
    return {"F_diff_ours_kg_m2": ours, "F_diff_bruteforce_kg_m2": brute, "F_diff_march_eq81_as_printed_kg_m2": march,
            "difference_ours_minus_march_kg_m2": ours - march,
            "analytic_difference_formula": "lam (1 - e^{-z_w/lam}) (2 rho_c - rho_w - rho_m)",
            "analytic_difference_value": lam * (1.0 - math.exp(-zw / lam)) * (2 * rc - rw - rm),
            "endpoint_sensitivity_bouguer_term_mGal_per_km": bouguer,
            "endpoint_sensitivity_yukawa_ours_mGal_per_km": yuk_ours,
            "endpoint_sensitivity_yukawa_march81_mGal_per_km": yuk_march}

# ======================================================= 2. tower experiments
# Half-space idealisation: surface data are harmonically continued; the Yukawa part of
# the surface field (2 pi G alpha rho lam, ~constant over the survey) is continued as a
# constant, while the true Yukawa field decays as exp(-h/lam). Residual (obs - pred):
#   Delta g(h) = -2 pi G_inf(alpha) alpha rho lam (1 - exp(-h/lam)) = -K alpha/(1+alpha)
#   with G_inf(alpha) = G_lab/(1+alpha) and K = 2 pi G_lab rho lam (1 - exp(-h/lam)) 1e8 [uGal]
# G normalisation: continuation uses measured surface gravity, so no lab G enters the
# leading term (lab G only in near-tower terrain corrections, neglected here).
RHO_CASES = [2000.0, 2670.0]                         # ASSUMED near-surface densities (not from the papers)
TOWERS = [
    {"id": "thomas1989", "h_m": 465.0, "mu_uGal": -60.0, "sigma_uGal": 95.0,
     "source": "Thomas et al. 1989 PRL 63, 1902 (abstract): (-60 +/- 95) 1e-8 m/s^2 at tower top"},
    {"id": "romaides1994_493m", "h_m": 493.0, "mu_uGal": -33.0, "sigma_uGal": 30.0,
     "source": "Romaides et al. 1994 PRD 50, 3608 (abstract): largest discrepancy -33 +/- 30 uGal at 493 m"},
    {"id": "romaides1997_568m", "h_m": 568.0, "mu_uGal": 32.0, "sigma_uGal": 32.0,
     "source": "Romaides et al. 1997 PRD 55, 4532 (abstract): 32 +/- 32 uGal at 568 m"},
]

def tower_model_residual_uGal(alpha, lam, rho, h):
    """STEP 5 (towers). Residual = observed - harmonically continued prediction at height h.
    Half-space of density rho (ASSUMED), flat surface: the Yukawa field at the surface,
    2 pi G_inf alpha rho lam, is continued as a constant, while the true field above the
    surface is 2 pi G_inf alpha rho lam e^{-h/lam}. Hence
        Delta g(h) = -2 pi G_inf(alpha) alpha rho lam (1 - e^{-h/lam})   [uGal],
        G_inf(alpha) = G_lab/(1+alpha)  (calibration, consistent for every alpha; alpha = 1 -> G_N/2)
    The Newtonian part is continued exactly."""
    return -2.0 * math.pi * g_inf_of(alpha) * alpha * rho * lam * (1.0 - math.exp(-h / lam)) * 1e8

def tower_K_uGal(lam, rho, h):
    """K = 2 pi G_lab rho lam (1 - e^{-h/lam}) 1e8 [uGal]; Delta g(alpha) = -K alpha/(1+alpha)."""
    return 2.0 * math.pi * G_LAB * rho * lam * (1.0 - math.exp(-h / lam)) * 1e8

def tower_alpha_lim_from_residual(lam, rho, h, mu, sig, k_sigma=Z95):
    """STEP 6 (towers). Exact root of Delta g(alpha) = mu - k_sigma*sig with G_inf(alpha) = G_lab/(1+alpha):
        E = -(mu - k_sigma*sig) > 0,   alpha_lim = E / (K - E).
    With sig read as 1 sigma Gaussian and k_sigma = 1.96: one-sided 97.5% upper limit on alpha > 0 (ASSUMPTION)."""
    K = tower_K_uGal(lam, rho, h)
    E = -(mu - k_sigma * sig)                                 # most negative allowed value, sign-flipped
    if E <= 0:
        return 0.0
    if E >= K:
        return math.inf
    return E / (K - E)

tower_residual_uGal = tower_model_residual_uGal       # backward-compatible aliases
tower_alpha_lim = tower_alpha_lim_from_residual

# BREN tower, LLNL report OSTI 6190822, p. "upward continuation" section: "Our predictions are within
# 93 uGal of observations at all heights on the tower, and the uncertainties at each height are as large
# as 95 uGal ... single standard deviation uncertainty from sampling errors". The SIGN of the 93 uGal is
# not stated, so we take the conservative case mu = -93 uGal (same sign as the alpha > 0 signal).
BREN = {"id": "bren_tower_osti6190822", "h_m": 454.0, "mu_uGal": -93.0, "sigma_uGal": 95.0}

def cruz_conditional_scaling_factor(h=300.0, lam_from=1000.0, lam_to=LAMBDA_M):
    """STEP 7 (Erie). CONDITIONAL. Ratio of the half-space residual response lam(1 - e^{-h/lam})
    at lam_from and lam_to for a top-of-tower residual at h. Assumes: (i) the published bound is set
    by the top-of-tower residual, (ii) flat half-space of uniform density, (iii) constant continued
    Yukawa surface field, (iv) the same G normalization at both lambdas. The actual response function
    of the Erie measurement (source geometry, survey radius, G normalization) is NOT known to us
    (full text not obtained), so the result is not a transferable bound."""
    S = lambda lam: lam * (1.0 - math.exp(-h / lam))
    return S(lam_from) / S(lam_to)

# ======================================================= 3. summary figure (context)
with open(os.path.join(HERE, "reference_inputs", "adelberger2003_fig4_digitized.json")) as f:
    DIG = json.load(f)
dig_lg = np.array([p[0] for p in DIG]); dig_a = np.array([p[1] for p in DIG])
dig_10km = float(10 ** np.interp(4.0, dig_lg, np.log10(dig_a)))

def back_substitution_tests(a_oc, a_oc_1s, eq81, bren, towers_per):
    """Review r1 (S01): every limit point must return the allowed residual edge.
    Ocean: model slope s(alpha_lim) recomputed with the full slab model (polyfit) vs edge.
    Towers/BREN: Delta g(alpha_lim) from tower_model_residual_uGal vs mu - k sigma.
    Pass if the relative deviation is below 1e-9."""
    tol = 1e-9
    out = []
    o = OCEAN
    for name, a, k in (("ocean_exact_k1.96", a_oc, Z95), ("ocean_exact_k1_march_interval", a_oc_1s, 1.0)):
        edge = o["slope_obs_mGal_per_km"] + k * o["slope_sigma_mGal_per_km"]
        got = ocean_model_residual_slope(a, LAMBDA_M)[0]
        rel = abs(got - edge) / abs(edge)
        out.append({"bound": name, "alpha_lim": a, "edge": edge, "model_at_alpha_lim": got, "unit": "mGal/km",
                    "rel_dev": rel, "pass": bool(rel < tol)})
    # eq. (81) scenario: chord model s81(alpha) = S81 alpha/(1+alpha)
    S81 = eq81["S81"]; edge = o["slope_obs_mGal_per_km"] + Z95 * o["slope_sigma_mGal_per_km"]
    got = S81 * eq81["alpha_lim_exact"] / (1.0 + eq81["alpha_lim_exact"])
    rel = abs(got - edge) / abs(edge)
    out.append({"bound": "ocean_march_eq81_scenario_exact", "alpha_lim": eq81["alpha_lim_exact"], "edge": edge,
                "model_at_alpha_lim": got, "unit": "mGal/km (endpoint-chord model)", "rel_dev": rel, "pass": bool(rel < tol)})
    for tid, (t, per) in towers_per.items():
        for key, v in per.items():
            rho = float(key.split("_")[1])
            edge = t["mu_uGal"] - Z95 * t["sigma_uGal"]
            got = tower_model_residual_uGal(v["alpha_lim_10km_rederived_95"], LAMBDA_M, rho, t["h_m"])
            rel = abs(got - edge) / abs(edge)
            out.append({"bound": f"{tid}/{key}", "alpha_lim": v["alpha_lim_10km_rederived_95"], "edge": edge,
                        "model_at_alpha_lim": got, "unit": "uGal", "rel_dev": rel, "pass": bool(rel < tol)})
    return {"tolerance_rel": tol, "all_pass": all(x["pass"] for x in out), "tests": out}

def r1(x, n=2):
    return float(f"{x:.{n}g}")

def main():
    res = {"step": "Step 27 (part of V28) - fixed Yukawa branch local gravity test",
           "series": "Time Attraction Law - Scalar-Tensor Yukawa Fit, V28 (V28 = V27 Step 26 cosmology package + interactive simulator + this Step 27)",
           "doi_V28": None,
           "prior_dois": {"V24": "10.5281/zenodo.23680715", "V25": "10.5281/zenodo.23699685", "V26": "10.5281/zenodo.23727585", "V27": "10.5281/zenodo.23761574"},
           "convention": "V = -G_inf m1 m2/r (1 + alpha exp(-r/lambda)); G_lab = G_inf (1+alpha) for r << lambda",
           "model": {"alpha_Y": ALPHA_Y, "beta_Yuk": BETA_YUK, "alpha_eff_q": Q, "lambda_m": LAMBDA_M,
                     "sign": "positive (attractive)", "universality": "universal: all matter couples to A^2(phi) g",
                     "screening": "none in the fixed branch (package: 'no screening branch is introduced')",
                     "G_eff_over_G_lab": 1.0 / (1.0 + Q), "sources": MODEL_SOURCES},
           "ratio_to_bound_definition": ("R = alpha_model / alpha_lim(lambda=10 km). Interpretable only in the "
               "convention above and under the stated assumptions; R > 1: model point outside the published/re-derived "
               "allowed region. R is not a significance and is never quoted in sigma."),
           "records": []}

    # --- ocean
    s_model, zz, resid = ocean_residual_slope(Q, LAMBDA_M)
    a_lim_oc, sens_oc = ocean_alpha_lim(LAMBDA_M)                              # exact inversion (primary)
    a_lim_oc_lin, _ = ocean_alpha_lim_from_residual(LAMBDA_M, linear=True)    # linear approximation edge/S
    a_lim_oc_1s, _ = ocean_alpha_lim_from_residual(LAMBDA_M, k_sigma=1.0)   # March eqs. (75)-(76) interval, exact
    xchk = march_eq81_F_difference(LAMBDA_M)
    o = OCEAN
    edge95 = o["slope_obs_mGal_per_km"] + Z95 * o["slope_sigma_mGal_per_km"]
    S81 = xchk["endpoint_sensitivity_bouguer_term_mGal_per_km"] + xchk["endpoint_sensitivity_yukawa_march81_mGal_per_km"]
    S_chord_ours = xchk["endpoint_sensitivity_bouguer_term_mGal_per_km"] + xchk["endpoint_sensitivity_yukawa_ours_mGal_per_km"]
    eq81 = {"S81": S81, "alpha_lim_exact": edge95 / (S81 - edge95), "alpha_lim_linear": edge95 / S81}
    pull_like = (s_model - o["slope_obs_mGal_per_km"]) / o["slope_sigma_mGal_per_km"]
    lam_scan = np.logspace(0, 5, 41)
    oc_scan = [[float(l), float(ocean_alpha_lim(l)[0])] for l in lam_scan]
    res["records"].append({
        "id": "zumberge1991_ocean",
        "citation": "M. A. Zumberge, J. A. Hildebrand, J. M. Stevenson, R. L. Parker, A. D. Chave, M. E. Ander, F. N. Spiess, "
                    "Phys. Rev. Lett. 67, 3051 (1991), doi:10.1103/PhysRevLett.67.3051; inputs as quoted by "
                    "R. March et al., Phys. Rev. D 100, 042002 (2019), arXiv:1904.12789, eqs. (46),(65),(73)-(76),(79)",
        "observable": "slope of gravity residual vs depth in the water column, 500-4800 m (G_ocean vs G_lab)",
        "G_handling": "residual defined with lab G_N; long-range G = G_N/(1+alpha Phi), Phi ~ 1 for lambda >= 10 m",
        "geometry": "volume-integrated infinite slabs (water 0-5 km, crust 5-12 km, mantle half-space)",
        "assumptions": "alpha > 0 universal Yukawa, no screening; 'extra force' of the NMC model of March et al. dropped",
        "CL_published": "abstract: 2 per mille fractional G uncertainty; slope -0.060 +/- 0.172 mGal/km treated as 1 sigma (assumption)",
        "published_bound": "alpha < 0.002 for 1 m <= lambda <= 5000 m (Zumberge et al., as quoted by March et al. 2019); "
                           "no published value at lambda = 10 km",
        "alpha_lim_10km_published": None,
        "alpha_lim_10km_rederived_95": a_lim_oc,
        "alpha_lim_10km_rederived_95_method": "exact inversion alpha_lim = edge/(S - edge), s(alpha) = S alpha/(1+alpha) with G_lab fixed",
        "alpha_lim_10km_rederived_95_linear_approx": a_lim_oc_lin,
        "ratio_to_bound_R_linear_approx": Q / a_lim_oc_lin,
        "alpha_lim_source": "re-derived in this work (idealized slab model above; edge of the quoted error x1.96 under an ASSUMED "
                            "Gaussian 1-sigma reading of +/-0.172, i.e. one-sided 97.5% under that assumption); not a published number",
        "input_data_status": "secondary-source data (not directly verified primary data): Zumberge et al. 1991 values as quoted by March et al. 2019; the Zumberge paper itself was not obtained",
        "provenance": "re-derived in this work from secondary-source data (unverified re-derivation)",
        "alpha_lim_10km_rederived_march_interval_1x": a_lim_oc_1s,
        "alpha_lim_variant_note": "k_sigma = 1 reproduces the interval March et al. use in eqs. (75)-(76) (-0.232..+0.112 mGal/km); exact inversion",
        "open_questions": "original uncertainty, weighting of the measured profiles and the model dependence of pressure-derived depth "
                          "are NOT reconstructed; the analytic layer check only verifies the idealized response",
        "cross_check_march_eq81": xchk,
        "R_status": "CONDITIONAL: denominator alpha_lim ~ 0.0023 is an unverified re-derivation",
        "reading_uncertainty": "model-dependent: layer densities/thicknesses of March et al. eq.(46); topography ignored",
        "model_point": {"alpha": Q, "lambda_m": LAMBDA_M},
        "model_observable": {"computed": True, "slope_pred_mGal_per_km": s_model,
                             "slope_obs_mGal_per_km": o["slope_obs_mGal_per_km"],
                             "slope_sigma_mGal_per_km": o["slope_sigma_mGal_per_km"],
                             "difference_in_quoted_sigma_units_info_only": pull_like,
                             "sensitivity_mGal_per_km_per_unit_alpha": sens_oc,
                             "apparent_G_ocean_over_G_lab_minus_1": -s_model * 1e-8 / (4 * math.pi * G_LAB * o["rho_w_mean"])},
        "ratio_to_bound_R": Q / a_lim_oc,
        "alpha_lim_scan_rederived": oc_scan,
    })

    # --- towers
    jek = {"id": "jekeli1990_tower",
           "citation": "C. Jekeli, D. H. Eckhardt, A. J. Romaides, Phys. Rev. Lett. 64, 1204 (1990), doi:10.1103/PhysRevLett.64.1204",
           "observable": "gravity on a ~600 m tower (WTVD, N. Carolina) vs upward-continued surface gravity",
           "G_handling": "upward continuation of measured surface gravity (no lab-G calibration of the leading term)",
           "geometry": "per abstract: refined analysis with detailed topography (not inspected here)",
           "assumptions": "Yukawa-potential coupling; sign / lambda range / CL not stated in the abstract",
           "CL_published": "not stated in abstract",
           "published_bound": "'constrains the Yukawa-potential coupling constant alpha to be less than 0.001' (abstract text)",
           "alpha_lim_10km_published": None,
           "abstract_alpha_value": 1e-3,
           "lambda_status": "unknown",
           "alpha_lim_source": "abstract text value; lambda range of validity unknown (full text paywalled)",
           "provenance": "primary abstract text, lambda unspecified",
           "applied": False,
           "role": "text-only reference, NOT applied: not plotted, not in the 10 km R summary (review r1, S02)",
           "reading_uncertainty": "text value; applicability at 10 km unknown",
           "model_point": {"alpha": Q, "lambda_m": LAMBDA_M},
           "model_observable": {"computed": False, "note": "tower residual for this tower not computed (no published residual in abstract)"},
           "ratio_to_bound_R": None,
           "R_status": "not applied: lambda range, sign and CL of the abstract value unknown"}
    res["records"].append(jek)
    towers_per = {}
    for t in TOWERS:
        per = {}
        for rho in RHO_CASES:
            al = tower_alpha_lim(LAMBDA_M, rho, t["h_m"], t["mu_uGal"], t["sigma_uGal"])
            per[f"rho_{int(rho)}"] = {"alpha_lim_10km_rederived_95": al,
                                      "model_residual_uGal": tower_residual_uGal(Q, LAMBDA_M, rho, t["h_m"]),
                                      "ratio_to_bound_R": Q / al}
        towers_per[t["id"]] = (t, per)
        res["records"].append({
            "id": t["id"], "citation": t["source"],
            "observable": f"tower-minus-continuation gravity residual at h = {t['h_m']:.0f} m",
            "G_handling": "harmonic upward continuation of measured surface gravity; leading term G-independent",
            "geometry": "idealized half-space Yukawa, volume-integrated: Delta g = -2 pi G_inf(alpha) alpha rho lambda (1 - e^{-h/lambda}), "
                        "G_inf(alpha) = G_lab/(1+alpha); terrain ignored",
            "open_questions": "actual surface sampling, terrain, depth-dependent density and the joint Yukawa response of the Poisson "
                              "continuation are NOT computed",
            "assumptions": "alpha > 0 universal, no screening; near-surface density ASSUMED 2000 or 2670 kg/m^3; "
                           "published +/- treated as 1 sigma",
            "CL_published": "abstract quotes value +/- uncertainty only",
            "published_bound": "abstract gives residual only, no alpha at 10 km",
            "alpha_lim_10km_published": None,
            "alpha_lim_source": "re-derived in this work (exact inversion E/(K-E); edge mu - 1.96 sigma under an ASSUMED Gaussian 1-sigma reading); not a published number",
            "input_data_status": "primary-paper abstract text (full text not obtained)",
            "provenance": "re-derived in this work from primary abstract values (unverified re-derivation)",
            "R_status": "conditional: denominator is an unverified re-derivation (assumed density, half-space, sigma read as 1 sigma)",
            "reading_uncertainty": "dominated by assumed density (factor 1.34) and the half-space idealisation",
            "model_point": {"alpha": Q, "lambda_m": LAMBDA_M},
            "model_observable": {"computed": True, "residual_obs_uGal": t["mu_uGal"], "residual_sigma_uGal": t["sigma_uGal"],
                                 "by_density": {k: {"model_residual_uGal": v["model_residual_uGal"]} for k, v in per.items()}},
            "by_density": per,
        })

    res["records"].append({
        "id": "adelberger2003_fig4_context",
        "citation": "E. G. Adelberger, B. R. Heckel, A. E. Nelson, Annu. Rev. Nucl. Part. Sci. 53, 77 (2003), arXiv:hep-ph/0307284, "
                    "Fig. 4 (based on Fischbach & Talmadge 1999, Fig. 2.13); compilation, NOT a primary measurement",
        "observable": "compilation of geophysical / Earth-LAGEOS constraints",
        "CL_published": "95%",
        "alpha_lim_10km_published": r1(dig_10km, 1),
        "alpha_lim_source": "digitized from the figure in this work (pixel trace of the shaded-region edge)",
        "provenance": "digitized compilation (not primary)",
        "reading_uncertainty": "about +/-0.1 dex (factor ~1.3); one significant digit",
        "model_point": {"alpha": Q, "lambda_m": LAMBDA_M},
        "model_observable": {"computed": False},
        "ratio_to_bound_R": r1(Q / dig_10km, 1),
        "R_status": "context only: digitized compilation, not used for the decision",
    })


    # --- Erie tower, detailed paper: primary numeric bound at lambda = 1000 m
    cruz_scale = cruz_conditional_scaling_factor(300.0, 1000.0, LAMBDA_M)
    res["records"].append({
        "id": "cruz1991_erie_tower",
        "citation": "J. Y. Cruz, J. C. Harrison, C. C. Speake, T. M. Niebauer, M. P. McHugh, P. T. Keyser, J. E. Faller, J. Makinen, "
                    "R. B. Beruff, J. Geophys. Res. 96(B12), 20073 (1991), doi:10.1029/91JB01756 (PRL version: Speake et al., "
                    "PRL 65, 1967 (1990), doi:10.1103/PhysRevLett.65.1967)",
        "observable": "gravity at 8 heights on a 300 m tower vs Newtonian prediction from surface gravity; rms difference 10 uGal",
        "G_handling": "prediction from measured surface gravity (no lab-G calibration of the leading term)",
        "geometry": "per abstract: prediction from surface gravity values (full analysis not inspected)",
        "assumptions": "bound quoted on |alpha| (sign-independent); universality implied by ISL test; CL not stated in abstract",
        "CL_published": "not stated in abstract",
        "published_bound": "'the absolute value of alpha must be less than 0.001 for lambda = 1000 m' (abstract text, Crossref record)",
        "alpha_lim_published_value": 1e-3, "alpha_lim_published_lambda_m": 1000.0,
        "alpha_lim_10km_published": None,
        "alpha_lim_10km_conditional_estimate": 1e-3 * cruz_scale,
        "alpha_lim_10km_status": ("CONDITIONAL ESTIMATE, not a bound. A 1 km limit cannot be transferred to 10 km using the "
                                  "Yukawa exponential alone: it needs the full response function of the tower measurement "
                                  "difference, including source geometry and G normalization, which we do not have (full text: HTTP 403)."),
        "conditional_scaling_assumptions": ["the published 1 km limit is set by the top-of-tower residual at h = 300 m",
                                            "flat half-space of uniform density around the tower",
                                            "the Yukawa part of the surface field is continued as a constant",
                                            "the same G normalization at both lambdas",
                                            "factor = S(1 km)/S(10 km), S(lam) = lam (1 - e^{-h/lam}) = " + format(cruz_scale, ".4f")],
        "alpha_lim_source": "primary abstract text at lambda = 1 km; no 10 km value published",
        "input_data_status": "primary-paper abstract text (Crossref record); full text not obtained",
        "provenance": "primary numeric at 1 km (abstract text); 10 km only as a conditional estimate",
        "reading_uncertainty": "text value (1 significant digit); 10 km transfer not verifiable",
        "model_point": {"alpha": Q, "lambda_m": LAMBDA_M},
        "model_observable": {"computed": False, "note": "per-height residuals not available in the abstract"},
        "ratio_to_bound_R_at_1km_value": Q / 1e-3,
        "ratio_to_bound_R_conditional_10km": Q / (1e-3 * cruz_scale),
        "ratio_to_bound_R": None,
        "R_status": "conditional estimate only; not used for the decision",
    })
    # --- BREN tower, LLNL conference report (OSTI 6190822): preliminary residual
    bren = {}
    for rho in RHO_CASES:
        al = tower_alpha_lim_from_residual(LAMBDA_M, rho, BREN["h_m"], BREN["mu_uGal"], BREN["sigma_uGal"])
        bren[f"rho_{int(rho)}"] = {"alpha_lim_10km_rederived_95": al,
                                   "model_residual_uGal": tower_model_residual_uGal(Q, LAMBDA_M, rho, BREN["h_m"]),
                                   "ratio_to_bound_R": Q / al}
    towers_per[BREN["id"]] = (BREN, bren)
    # illustrative only: bounded +/-100 uGal systematic (report p. 12) taken at its worst-case edge
    bren_nuis = {f"rho_{int(r)}": tower_alpha_lim_from_residual(LAMBDA_M, r, BREN["h_m"], BREN["mu_uGal"] - 100.0, BREN["sigma_uGal"]) for r in RHO_CASES}
    res["records"].append({
        "id": "bren_tower_osti6190822",
        "citation": "J. Thomas, P. Kasameyer et al., 'A test of Newton's Law of Gravity using the BREN Tower, Nevada', LLNL "
                    "conference report, OSTI 6190822 (https://www.osti.gov/biblio/6190822); preliminary version of PRL 63, 1902",
        "observable": ("report text: predictions within 93 uGal of observations at all heights, uncertainties up to 95 uGal "
                       "(1 sigma, sampling); sign not stated -> conservative mu = -93 uGal, sigma = 95 uGal at h = 454 m"),
        "G_handling": "Laplace upward continuation of surface data", "geometry": "half-space Yukawa (ours), terrain ignored",
        "assumptions": "as for the other towers; density assumed", "CL_published": "+/- only",
        "published_bound": "no alpha value in report", "alpha_lim_10km_published": None,
        "alpha_lim_source": "re-derived in this work (one-sided 97.5% under Gaussian 1 sigma)",
        "input_data_status": "primary document (LLNL report, preliminary results), obtained from OSTI",
        "provenance": "re-derived in this work from a primary preliminary report (unverified re-derivation)",
        "R_status": "conditional: denominator is an unverified re-derivation",
        "reading_uncertainty": "assumed density, half-space idealisation",
        "systematic_note": ("the report's ~100 uGal spherical-approximation systematic needs a separate nuisance treatment; it cannot be "
                            "added automatically as an independent Gaussian sigma to the 95 uGal"),
        "illustrative_bounded_nuisance_100uGal": {"alpha_edge": bren_nuis,
            "note": "ILLUSTRATIVE ONLY: a bounded +/-100 uGal nuisance taken at its worst-case edge; not a new data fit and not an empirical bound"},
        "model_point": {"alpha": Q, "lambda_m": LAMBDA_M},
        "model_observable": {"computed": True, "by_density": {k: {"model_residual_uGal": v["model_residual_uGal"]} for k, v in bren.items()}},
        "by_density": bren,
    })
    # --- Fischbach & Talmadge 1996 compilation (context, alpha < 0 panel)
    res["records"].append({
        "id": "fischbach_talmadge1996_fig1_context",
        "citation": "E. Fischbach, C. Talmadge, arXiv:hep-ph/9606249 (1996), Fig. 1 (2 sigma); compilation, NOT primary",
        "observable": "compilation (Tower, Earth-LAGEOS segments around 10 km)",
        "CL_published": "2 sigma",
        "sign_note": "Fig. 1 vertical axis is -alpha, i.e. the plotted constraint is for alpha < 0 (repulsive); NOT directly "
                     "applicable to the model's alpha > 0",
        "alpha_lim_10km_published": None,
        "alpha_lim_10km_digitized": 4e-4,
        "alpha_lim_source": "digitized in this work from the rendered figure (pixel trace); 3e-4 to 5e-4 depending on axis calibration",
        "provenance": "digitized compilation (not primary), opposite sign",
        "reading_uncertainty": "about +/-0.15 dex; one significant digit",
        "model_point": {"alpha": Q, "lambda_m": LAMBDA_M},
        "model_observable": {"computed": False},
        "ratio_to_bound_R": None,
    })
    res["primary_search_log"] = [
        "Jekeli et al. PRL 64, 1204 (1990): only abstract accessible (APS paywall); 'alpha < 0.001', lambda range and CL not stated",
        "Speake et al. PRL 65, 1967 (1990) / Cruz et al. JGR 96, 20073 (1991): abstract gives |alpha| < 0.001 at lambda = 1000 m; Wiley PDF 403",
        "Romaides et al. PRD 50, 3608 (1994) and PRD 55, 4532 (1997): abstracts give residuals only; full text paywalled",
        "Thomas et al. PRL 63, 1902 (1989): abstract residual only; LLNL OSTI reports 6190822 and 7173840 opened: residuals, no alpha-lambda numbers",
        "Zumberge et al. PRL 67, 3051 (1991): abstract (G = 6.677 +/- 0.013); per March et al. 2019 the published bound alpha < 0.002 covers 1-5000 m only",
        "March et al. PRD 100, 042002 (2019), arXiv:1904.12789: bounds for the NMC model (with extra force), not a pure-Yukawa table at 10 km",
        "Eckhardt et al. PRL 60, 2567 (1988): original WTVD claim, superseded by Jekeli et al. 1990",
        "Fischbach & Talmadge hep-ph/9606249 Fig. 1: alpha < 0 panel only",
        "Adelberger et al. 2003 Fig. 4: |alpha| compilation, digitized",
        "Not reached (paywalled / no open copy found in this pass): Stacey et al. RMP 59, 157 (1987); Moore et al. PRD 38, 1023 (1988); "
        "Ander et al. PRL 62, 985 (1989); Achilli et al. 1997; Bartlett & Tew PRL 63, 1531 (1989); Rapp 1987; Fischbach & Talmadge 1999 book tables",
    ]
    res["decisive_bound"] = {
        "record": "zumberge1991_ocean",
        "alpha_lim_10km": a_lim_oc, "R": Q / a_lim_oc,
        "alpha_lim_method": "exact inversion edge/(S - edge)",
        "alpha_lim_10km_linear_approx": a_lim_oc_lin, "R_linear_approx": Q / a_lim_oc_lin,
        "R_rounded_denominator_0p0023": Q / 0.0023,
        "R_name": "ratio to bound (határhoz viszonyított arány); not a significance, never in sigma",
        "R_status": "CONDITIONAL: denominator alpha_lim ~ 0.0023 is an unverified re-derivation",
        "provenance": "re-derived in this work from secondary-source data (not directly verified primary data) (Zumberge et al. 1991 values as quoted by March et al. 2019); "
                      "no primary source found that states a numeric alpha limit at exactly 10 km",
        "variant_march_interval_1x": {"alpha_lim_10km": a_lim_oc_1s, "R": Q / a_lim_oc_1s},
        "scenario_march_eq81_as_printed": {"alpha_lim_10km": eq81["alpha_lim_exact"], "R": Q / eq81["alpha_lim_exact"],
                                           "alpha_lim_linear_approx": eq81["alpha_lim_linear"],
                                           "note": "separate scenario, not a statistical error band; the eq. (81) discrepancy is not resolvable by a "
                                                   "sign or G-convention change; the authors' intent/correction was not verified"},
        "scenarios_note": "R values across error readings, layer formulas and assumed densities are scenarios, not one combined confidence interval",
        "not_applied": "Cruz et al. 1991 |alpha| < 0.001 is stated at 1 km only (10 km value only a conditional estimate); Jekeli et al. 1990 "
                       "alpha < 0.001 has an unknown lambda range; neither enters the 10 km decision",
    }
    res["missing_primary_sources"] = [
        {"source": "Zumberge et al., PRL 67, 3051 (1991), doi:10.1103/PhysRevLett.67.3051",
         "status": "not obtained; ocean inputs taken from March et al. 2019 (secondary source)",
         "needed_for": "residual definition, meaning of +/-0.172 mGal/km, corrections and their uncertainties"},
        {"source": "Cruz et al., JGR 96, 20073 (1991), doi:10.1029/91JB01756", "status": "abstract only; full text not obtained (publisher PDF returned HTTP 403)",
         "needed_for": "response function of the Erie measurement, CL, lambda dependence of the bound"},
        {"source": "Jekeli, Eckhardt, Romaides, PRL 64, 1204 (1990)", "status": "abstract only (paywalled)",
         "needed_for": "lambda range, sign and CL of alpha < 0.001 (full text)"},
        {"source": "Thomas et al., PRL 63, 1902 (1989)", "status": "abstract only; preliminary LLNL report OSTI 6190822 obtained",
         "needed_for": "published per-height residuals, sign, CL, terrain treatment"},
        {"source": "Romaides et al., PRD 50, 3608 (1994) and PRD 55, 4532 (1997)", "status": "abstracts only (paywalled)",
         "needed_for": "residual definitions, densities, any alpha-lambda curve"},
        {"source": "Stacey et al. RMP 59, 157 (1987); Moore et al. PRD 38, 1023 (1988); Ander et al. PRL 62, 985 (1989); "
                   "Bartlett & Tew PRL 63, 1531 (1989); Achilli et al. 1997; Rapp 1987; Fischbach & Talmadge 1999 (book)",
         "status": "not reached (no open copy found)", "needed_for": "possible primary numeric bounds near 10 km"},
    ]
    res["review_trace"] = {
        "ocean": {"slope_obs_mGal_per_km": OCEAN["slope_obs_mGal_per_km"], "sigma_quoted_mGal_per_km": OCEAN["slope_sigma_mGal_per_km"],
                  "sensitivity_mGal_per_km_per_unit_alpha": sens_oc,
                  "of_which_bouguer_G_renormalization_endpoint": xchk["endpoint_sensitivity_bouguer_term_mGal_per_km"],
                  "of_which_yukawa_volume_integral_endpoint": xchk["endpoint_sensitivity_yukawa_ours_mGal_per_km"],
                  "edge_95_two_sided_mGal_per_km": OCEAN["slope_obs_mGal_per_km"] + Z95 * OCEAN["slope_sigma_mGal_per_km"],
                  "alpha_lim": a_lim_oc, "alpha_lim_linear_approx": a_lim_oc_lin,
                  "S_chord_0_5000m_ours_mGal_per_km": S_chord_ours, "edge_march_1x_mGal_per_km": OCEAN["slope_obs_mGal_per_km"] + OCEAN["slope_sigma_mGal_per_km"],
                  "alpha_lim_march_1x": a_lim_oc_1s, "model_slope_at_q_mGal_per_km": s_model,
                  "check_slope_q_vs_S_q_over_1_plus_q": s_model / (sens_oc * Q / (1.0 + Q)),
                  "variant_march_eq81_as_printed_endpoint_chord": {
                      "sensitivity_mGal_per_km_per_unit_alpha": xchk["endpoint_sensitivity_bouguer_term_mGal_per_km"] + xchk["endpoint_sensitivity_yukawa_march81_mGal_per_km"],
                      "alpha_lim": eq81["alpha_lim_exact"], "alpha_lim_linear_approx": eq81["alpha_lim_linear"],
                      "R": Q / eq81["alpha_lim_exact"],
                      "note": "separate scenario: if March eq. (81) as printed were the layer integral; our direct and brute-force integrals "
                              "disagree with it; not resolvable by a sign or G-convention change; author intent not verified"},
                  "endpoint_chord_vs_fitted_slope_note": "endpoint sensitivities are chord values over 0-5000 m; the fitted slope uses 500-4800 m"},
        "bren": {"h_m": BREN["h_m"], "mu_uGal": BREN["mu_uGal"], "sigma_uGal": BREN["sigma_uGal"],
                 "edge_uGal": BREN["mu_uGal"] - Z95 * BREN["sigma_uGal"],
                 "K_uGal_alpha_to_0_sensitivity": {f"rho_{int(r)}": tower_K_uGal(LAMBDA_M, r, BREN["h_m"]) for r in RHO_CASES},
                 "residual_at_alpha_1_uGal_G_N_over_2": {f"rho_{int(r)}": tower_model_residual_uGal(1.0, LAMBDA_M, r, BREN["h_m"]) for r in RHO_CASES},
                 "K_all_towers_uGal": {t["id"]: {f"rho_{int(r)}": tower_K_uGal(LAMBDA_M, r, t["h_m"]) for r in RHO_CASES} for t in TOWERS + [BREN]},
                 "alpha_lim": {k: v["alpha_lim_10km_rederived_95"] for k, v in bren.items()},
                 "variant_if_sign_were_positive_mu_plus93": {f"rho_{int(r)}": tower_alpha_lim_from_residual(LAMBDA_M, r, BREN["h_m"], 93.0, 95.0) for r in RHO_CASES},
                 "illustrative_bounded_nuisance_100uGal": bren_nuis},
    }

    # verdict
    res["back_substitution_tests"] = back_substitution_tests(a_lim_oc, a_lim_oc_1s, eq81, bren, towers_per)
    Rs = [Q / a_lim_oc] + [v["ratio_to_bound_R"] for r in res["records"] if "by_density" in r for v in r["by_density"].values()]
    res["verdict"] = {
        "label": "Step 27 (part of V28)",
        "scope": "fixed, unscreened, universal local Yukawa branch (alpha = 0.0795775, lambda = 10 km) only; not all TAL variants",
        "min_R": min(Rs), "max_R": max(Rs),
        "not_a_reopening_of": "Step 26 (V26/V27 kinetic normalization) is not re-opened; S01 concerns only holding G_lab fixed in a new parameter-limit search",
        "status": ("strong tension only in this sense: in the idealized layer/half-space models, under the assumed error reading, a large "
                   "numerical discrepancy results. Not an empirical exclusion with a confidence level."),
        "status_hu": ("erős feszültség csak ebben az értelemben: az idealizált réteg-/féltérmodellekben, a feltételezett hibaértelmezés mellett "
                      "nagy numerikus eltérés adódik. Nem konfidenciaszinttel rendelkező empirikus kizárás."),
        "R_name": "ratio to bound (határhoz viszonyított arány)",
        "R_status": "CONDITIONAL: denominator alpha_lim ~ 0.0023 is an unverified re-derivation",
        "R_decisive": Q / a_lim_oc, "R_decisive_rounded_denominator_0p0023": Q / 0.0023,
        "statement_en": ("Our own re-derivation indicates strong tension with the fixed, unscreened Yukawa branch; confirming an "
                         "exclusion requires verification of the measurement response and the uncertainty treatment."),
        "statement_hu": ("A saját újraszármaztatás erős feszültséget jelez a rögzített, árnyékolatlan Yukawa-ággal; a kizárás "
                         "megerősítéséhez a mérési válasz és a bizonytalanságkezelés ellenőrzése szükséges."),
        "R_range_note": ("min_R/max_R span the re-derived ocean and tower records (Jekeli text value not applied; Erie only conditional); "
                         "all CONDITIONAL. The range collects scenarios (error readings, layer formulas, assumed densities), it is NOT a "
                         "combined confidence interval, and R is never a sigma significance"),
        "data_status": "Zumberge inputs are secondary-source data via March et al. 2019 (not directly verified primary data)",
        "missing_primary_sources_short": ["Zumberge et al. 1991 (not obtained)", "Cruz et al. 1991 (full text)", "Jekeli et al. 1990 (full text)"],
        "untested_candidates_to_retest": ["screening of the local Yukawa for Earth-scale/laboratory sources (not in the package; a new branch)",
                                          "alpha <~ 1e-3 at lambda = 10 km (i.e. alpha_Y^2 beta_Yuk reduced by >~ 80x; a new branch)",
                                          "a different lambda (e.g. lambda < 1 cm) with unchanged alpha: an UNTESTED candidate to be re-tested; "
                                          "Step 27 did not examine whether that range is allowed, so it is NOT a verified escape range"]}
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=2, allow_nan=False)

    # ------------------------------------------------ figure 1: alpha-lambda plane
    for lang in ("en", "hu"):
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        msk = dig_lg <= 5.5
        ax.plot(10 ** dig_lg[msk], dig_a[msk], color="0.45", lw=1.8,
                label=("Adelberger+2003 Fig.4 (digitized, approx., 95%)" if lang == "en"
                       else "Adelberger+2003 4. ábra (digitalizált, közelítő, 95%)"))
        ax.plot([1, 5000], [2e-3, 2e-3], color="tab:blue", lw=2.2,
                label=("Zumberge+1991 ocean: alpha<0.002, 1-5000 m (quoted)" if lang == "en"
                       else "Zumberge+1991 óceán: alpha<0,002, 1-5000 m (idézett)"))
        sc = np.array(oc_scan)
        ax.plot(sc[:, 0], sc[:, 1], color="tab:blue", ls="--", lw=1.4,
                label=("ocean, own bound, CONDITIONAL (secondary-source data; idealized layer model; quoted error x1.96, assumed Gaussian 1σ)" if lang == "en"
                       else "óceán, saját korlát, FELTÉTELES (másodlagos forrás adatai; idealizált rétegmodell; idézett hiba ×1.96, feltételezett Gauss 1σ)"))
        for t, mk in zip(TOWERS, ["s", "^", "v"]):
            vals = [tower_alpha_lim(LAMBDA_M, rho, t["h_m"], t["mu_uGal"], t["sigma_uGal"]) for rho in RHO_CASES]
            ax.errorbar([LAMBDA_M], [math.sqrt(vals[0] * vals[1])], yerr=[[math.sqrt(vals[0]*vals[1]) - min(vals)], [max(vals) - math.sqrt(vals[0]*vals[1])]],
                        fmt=mk, color="tab:green", ms=6, capsize=3,
                        label=(f"{t['id']} re-derived, conditional (idealized half-space, rho 2000-2670)" if lang == "en"
                               else f"{t['id']} újraszámolva, feltételes (idealizált féltér, rho 2000-2670)"))
        ax.plot([1000.0], [1e-3], "P", color="tab:purple", ms=8,
                label=("Cruz+1991 Erie tower: |alpha|<0.001 at 1 km (abstract; 10 km transfer only conditional)" if lang == "en"
                       else "Cruz+1991 Erie-torony: |alpha|<0,001 1 km-nél (absztrakt; 10 km-re csak feltételesen)"))
        ax.plot([LAMBDA_M], [Q], "*", color="tab:red", ms=16,
                label=("TAL fixed branch: alpha=0.0796, lambda=10 km" if lang == "en"
                       else "TAL rögzített ág: alpha=0,0796, lambda=10 km"))
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(1e-1, 1e8); ax.set_ylim(1e-6, 1)
        ax.set_xlabel("lambda [m]"); ax.set_ylabel("|alpha|")
        ax.set_title("Step 27 (part of V28): Yukawa alpha-lambda, fixed branch vs bounds" if lang == "en"
                     else "27. lépés (a V28 része): Yukawa alpha-lambda, rögzített ág és korlátok")
        ax.legend(fontsize=6, loc="lower left"); ax.grid(alpha=0.3, which="both")
        fig.tight_layout(); fig.savefig(os.path.join(HERE, f"step27_alpha_lambda_{lang}.png"), dpi=160); plt.close(fig)

    # ------------------------------------------------ figure 2: observables
    for lang in ("en", "hu"):
        fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
        ax = axs[0]
        ax.plot(zz / 1e3, resid * 1e5, color="tab:red",
                label=("model residual (idealized layer/half-space model)" if lang == "en"
                       else "modell-reziduum (idealizált réteg-/féltérmodell)"))
        zline = np.array([0.5, 4.8])
        o = OCEAN
        ax.fill_between(zline, (o["slope_obs_mGal_per_km"] - Z95 * o["slope_sigma_mGal_per_km"]) * (zline - 0.5),
                        (o["slope_obs_mGal_per_km"] + Z95 * o["slope_sigma_mGal_per_km"]) * (zline - 0.5),
                        color="tab:blue", alpha=0.3, label=("quoted error ×1.96; assumed Gaussian 1σ" if lang == "en" else "idézett hiba ×1.96; feltételezett Gauss 1σ"))
        ax.set_xlabel("depth [km]" if lang == "en" else "mélység [km]")
        ax.set_ylabel("Delta g residual [mGal]" if lang == "en" else "Delta g reziduum [mGal]")
        ax.set_title("Ocean (Zumberge+1991 setup)" if lang == "en" else "Óceán (Zumberge+1991 elrendezés)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        ax = axs[1]
        hh = np.linspace(0, 600, 121)
        for rho, ls in zip(RHO_CASES, ["-", "--"]):
            ax.plot(hh, [tower_residual_uGal(Q, LAMBDA_M, rho, h) / 1e3 for h in hh], ls, color="tab:red",
                    label=(f"model, rho={int(rho)} (idealized half-space)" if lang == "en" else f"modell, rho={int(rho)} (idealizált féltér)"))
        for t in TOWERS:
            ax.errorbar([t["h_m"]], [t["mu_uGal"] / 1e3], yerr=[Z95 * t["sigma_uGal"] / 1e3], fmt="o", color="k", capsize=3)
        ax.set_xlabel("height [m]" if lang == "en" else "magasság [m]")
        ax.set_ylabel("residual [mGal]" if lang == "en" else "reziduum [mGal]")
        ax.set_title("Towers (quoted error ×1.96; assumed Gaussian 1σ)" if lang == "en" else "Tornyok (idézett hiba ×1.96; feltételezett Gauss 1σ)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        fig.suptitle("Step 27 (part of V28)" if lang == "en" else "27. lépés (a V28 része)", fontsize=10)
        fig.tight_layout(); fig.savefig(os.path.join(HERE, f"step27_observables_{lang}.png"), dpi=160); plt.close(fig)

    print(json.dumps({"q": Q, "ocean_slope_pred": s_model, "ocean_S": sens_oc, "ocean_alpha_lim": a_lim_oc, "R_ocean": Q / a_lim_oc,
                      "ocean_alpha_lim_linear": a_lim_oc_lin, "R_ocean_linear": Q / a_lim_oc_lin,
                      "ocean_1x": a_lim_oc_1s, "R_1x": Q / a_lim_oc_1s, "eq81": eq81, "R_eq81": Q / eq81["alpha_lim_exact"],
                      "march_endpoint": {k: v for k, v in xchk.items() if k.startswith("endpoint") or k.startswith("F_diff")},
                      "S_chord_ours": S_chord_ours, "bren_nuis": bren_nuis,
                      "backsub_all_pass": res["back_substitution_tests"]["all_pass"],
                      "backsub_max_rel": max(x["rel_dev"] for x in res["back_substitution_tests"]["tests"]),
                      "minR": min(Rs), "maxR": max(Rs),
                      "dig_10km": dig_10km,
                      "towers": {r["id"]: r["by_density"] for r in res["records"] if "by_density" in r}}, indent=1))

if __name__ == "__main__":
    main()
