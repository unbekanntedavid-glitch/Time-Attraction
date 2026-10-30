#!/usr/bin/env python3
"""Build data.js (embedded inputs + package reference numbers) from the V27 and V26 trees (read-only).

Input paths are portable (see ../README_DATA_en.md):
  V27 (Step 26) tree: --v27 DIR, or env TAL_V27_DIR; default: ../../step26_cosmology_desi_bao
                      (inside an extracted TAL_V28/ folder).
  V26 tree:           --v26 DIR, or env TAL_V26_DIR; required, not part of the V28 package
                      (V26, doi:10.5281/zenodo.23727585, or the GitHub V26/ folder).
  SPARC rotmod files: --sparc DIR, or env TAL_SPARC_DIR; default: <V26>/data/Rotmod_LTG.
  Output:             --out FILE; default: data.js next to this script.
"""
import os, sys, json, glob, csv, argparse
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--v27", default=os.environ.get("TAL_V27_DIR", os.path.join(HERE, "..", "..", "step26_cosmology_desi_bao")))
ap.add_argument("--v26", default=os.environ.get("TAL_V26_DIR"))
ap.add_argument("--sparc", default=os.environ.get("TAL_SPARC_DIR"))
ap.add_argument("--out", default=os.path.join(HERE, "data.js"))
args = ap.parse_args()
V27 = os.path.abspath(args.v27)
if not args.v26:
    sys.exit("build_data.py: the V26 input tree is required (--v26 DIR or TAL_V26_DIR); see ../README_DATA_en.md")
V26 = os.path.abspath(args.v26)
SPARC = os.path.abspath(args.sparc) if args.sparc else os.path.join(V26, "data", "Rotmod_LTG")
need = [os.path.join(V27, n) for n in ("tal_cosmology_common.py", "run_step1_2_and_1_3.py", "run_step1_4_to_1_6.py", "tal_hubble_desi_bao_summary.json")]
need += [os.path.join(V26, n) for n in ("gw170817_and_lensing_summary.json", "lensing_deflection_and_ppn_slip.csv", "gw170817_dispersion_and_cT.csv")]
miss = [p for p in need if not os.path.isfile(p)]
if not glob.glob(os.path.join(SPARC, "*_rotmod.dat")): miss.append(os.path.join(SPARC, "*_rotmod.dat"))
if miss:
    sys.exit("build_data.py: missing input files:\n  " + "\n  ".join(miss))
sys.path.insert(0, V27)
import tal_cosmology_common as tcc  # noqa: E402
import importlib.util  # noqa: E402

def load_mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

s13 = load_mod("s13", os.path.join(V27, "run_step1_2_and_1_3.py"))
s14 = load_mod("s14", os.path.join(V27, "run_step1_4_to_1_6.py"))

desi = {"z": [float(v) for v in tcc.COBAYA_Z], "mean": [float(v) for v in tcc.COBAYA_MEAN], "q": list(tcc.COBAYA_Q),
        "cov": [[float(c) for c in row] for row in tcc.COBAYA_COV], "table1": []}
for b in tcc.DESI_BINS:
    desi["table1"].append({"id": b["bin_id"], "tracer": b["tracer"], "z": b["z_eff"], "primary": b["primary_12pt"], "mode": b["mode"],
                           "DM": b.get("DM_rd"), "DMe": b.get("DM_err"), "DH": b.get("DH_rd"), "DHe": b.get("DH_err"), "r": b.get("r"),
                           "DV": b.get("DV_rd"), "DVe": b.get("DV_err")})

summ = json.load(open(os.path.join(V27, "tal_hubble_desi_bao_summary.json")))
s11 = summ["step1_1_desi_dr1_bao"]; mb = summ["step2_m1_to_m4_and_inverse_ladder"]["m_branches"]
s12 = summ["step1_2_planck_cmb_and_illustrative_3d_prior"]; s13j = summ["step1_3_synthetic_distance_modulus_shape_control"]
s14j = summ["step1_4_rsd_fsigma8_growth"]; bbn = summ["step1_5_algebraic_bbn_benchmark"]
rc = summ["step26_revision_controls"]["review2_F01_F12"]
models = {m["model_id"]: m for m in s12["models"]}
ref = {
    "rd_zdrag": s11["rd_zdrag_Mpc"], "rs_zstar": s11["rs_zstar_Mpc"],
    "chi2_m4": s11["chi2_Step25_M4_Cobaya_full"], "chi2_m4_tab1": s11["chi2_Step25_M4_Table1"],
    "chi2_rigid": s11["chi2_Rigid72p80_Cobaya_full"], "chi2_rigid_tab1": s11["chi2_Rigid72p80_Table1"],
    "P1": s11["profiled_1param_M4_shape"]["P_best_km_s"], "sigmaP1": s11["profiled_1param_M4_shape"]["P_sigma_km_s"],
    "chi2_P1": s11["profiled_1param_M4_shape"]["chi2_min_Cobaya_full"], "rdReq1": s11["profiled_1param_M4_shape"]["rd_req_at_H0_72p80_Mpc"],
    "om2": s11["profiled_2param_flat_LCDM"]["Omega_m_MLE"], "sigmaOm2": s11["profiled_2param_flat_LCDM"]["Omega_m_sigma"],
    "P2": s11["profiled_2param_flat_LCDM"]["P_best_km_s"], "sigmaP2": s11["profiled_2param_flat_LCDM"]["P_sigma_km_s"],
    "chi2_2p": s11["profiled_2param_flat_LCDM"]["chi2_min_Cobaya_full"],
    "rdReq2": s11["profiled_2param_flat_LCDM"]["rd_req_at_H0_72p80_Mpc"], "sigmaRdReq2": s11["profiled_2param_flat_LCDM"]["rd_req_sigma_fixedH0_Mpc"],
    "sigmaRdReq2H0": s11["profiled_2param_flat_LCDM"]["rd_req_sigma_with_H0_err_Mpc"],
    "H0PlanckRd": s11["profiled_2param_flat_LCDM"]["H0_from_Planck_rd_km_s_Mpc"], "sigmaH0PlanckRd": s11["profiled_2param_flat_LCDM"]["H0_from_Planck_rd_sigma_km_s_Mpc"],
    "joint": s11["joint_3param_BAO_PlanckRd_LocalH0"],
    "m3": mb, "theta_m4": s12["theta_star_100_M4_reint"], "dm_zstar_m4": s12["DM_zstar_M4_Mpc"],
    "planck3d_m4": models["Step25_M4_TAL_Baseline"]["chi2_Planck_3D_illustrative_prior"],
    "flat728_sameOm": models["Flat_H0_72p80_SameOm0p3148"], "flat728_lowOm": models["Flat_H0_72p80_Om0p244828"],
    "sne": s13j, "rsd": {k: s14j[k] for k in s14j if not isinstance(s14j[k], (list, dict))}, "bbn": bbn,
    "G_ratio": rc["F02_G_dictionary"]["G_eff_over_G_lab_large_scale"], "local": summ["step1_6_local_yukawa_bulk_benchmarks"],
    "gammaL": summ["tal_fundamental_parameters"]["gamma_L_m3_kg"],
}
# V26 Step 23 reference values (GW170817 and lensing)
gw = json.load(open(os.path.join(V26, "gw170817_and_lensing_summary.json")))
ref["gw_meta"] = gw["gw170817_and_scalar_dispersion_meta"]
ref["solar"] = gw["solar_cassini_1Msun_benchmark"]["rows"]
ref["ns_ref"] = gw["v25_reference_neutron_star_SLy_15p0"]
lens_rows = []
with open(os.path.join(V26, "lensing_deflection_and_ppn_slip.csv")) as f:
    for r in csv.DictReader(f):
        if r["source_body"].startswith("NS"):
            lens_rows.append({k: r[k] for k in ("scale_label", "b_km", "gamma_PPN", "Phi_lens_over_Psi_dyn", "alpha_hat_Yukawa_conf_rad", "alpha_hat_Schw_MADM_rad", "Shapiro_delay_Schw_MADM_us")})
ref["ns_lens"] = lens_rows
gw_rows = []
with open(os.path.join(V26, "gw170817_dispersion_and_cT.csv")) as f:
    for r in csv.DictReader(f):
        gw_rows.append({k: r[k] for k in ("mode_label", "frequency_Hz", "scalar_regime", "scalar_v_group_over_c", "scalar_delay_at_40Mpc_s", "log10_scalar_amp_at_40Mpc")})
ref["gw_rows"] = gw_rows
ref["macro"] = gw["macroscopic_lensing_and_bullet_cluster_bridge"]

# SPARC rotation curves (175 galaxies)
sparc = []
for fn in sorted(glob.glob(os.path.join(SPARC, "*_rotmod.dat"))):
    name = os.path.basename(fn).replace("_rotmod.dat", ""); dist = None; pts = []
    for line in open(fn):
        if line.startswith("#"):
            if "Distance" in line: dist = float(line.split("=")[1].split()[0])
            continue
        p = line.split()
        if len(p) >= 6: pts.append([float(x) for x in p[:6]])
    sparc.append({"n": name, "D": dist, "d": pts})

import math
def clean(o):
    if isinstance(o, float): return o if math.isfinite(o) else None
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list): return [clean(v) for v in o]
    return o
data = {"desi": desi, "sne_grid": [list(x) for x in s13.PANTHEON_GRID], "rsd": [list(x) for x in s14.RSD_DATA], "sparc": sparc, "ref": ref}
txt = "window.TAL_DATA = " + json.dumps(clean(data), separators=(",", ":"), allow_nan=False) + ";\n"
txt = txt.replace("window.TAL_DATA", "(typeof window !== 'undefined' ? window : globalThis).TAL_DATA", 1)
open(args.out, "w").write(txt)
print("data.js", len(txt), "bytes;", len(sparc), "galaxies")
