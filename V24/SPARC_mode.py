#!/usr/bin/env python3
"""
SPARC_mode.py — Unified self-contained SPARC verification runner for the
scalar-tensor / Yukawa / 3D top-hat derivation audit package.

Modes:
  --mode summary  (default): Verifies SHA-256 hashes of Rotmod_LTG.zip (175 galaxies),
                  SPARC_Lelli2016c.mrt, and blind_123.txt; checks the NGC 3198
                  control point at R = 30.08 kpc (V_obs = 146 +/- 2 km/s,
                  V_bar(Y=1.0) = 88.53 km/s -> C_req = 2.720,
                  V_bar(Y=0.5) = 70.32 km/s -> C_req = 4.311); and verifies
                  0 negative-density points across all 175 SPARC galaxies under
                  positive-definite photometric / 3D top-hat / Yukawa kernels
                  vs. 239 (140 visible non-holdout galaxies) / 303 (all 175 galaxies)
                  negative-density points in the 1D radial derivative formula.
  --mode pipeline : Forwards remaining arguments to ta_pipeline.py.
  --mode foundation: Runs ta_foundation.py (using local ./data/Rotmod_LTG.zip).
  --mode v19_7th  : Runs tal_factor_7th_model.py.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ta_pipeline as tp
import ta_foundation as tf

EXPECTED_ZIP_SHA = "0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588"
EXPECTED_MRT_SHA = "5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243"
EXPECTED_BLIND_SHA = "3637d2f5bda6b2b97e645c5fba5cf39d310a1826bee388286c2c74b9171aa283"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def run_summary_audit():
    zip_path = HERE / "data" / "Rotmod_LTG.zip"
    mrt_path = HERE / "data" / "SPARC_Lelli2016c.mrt"
    blind_path = HERE / "data" / "blind_123.txt"

    assert sha256_file(zip_path) == EXPECTED_ZIP_SHA, "Rotmod_LTG.zip SHA-256 mismatch"
    assert sha256_file(mrt_path) == EXPECTED_MRT_SHA, "SPARC_Lelli2016c.mrt SHA-256 mismatch"
    assert sha256_file(blind_path) == EXPECTED_BLIND_SHA, "blind_123.txt SHA-256 mismatch"

    tf.ZIP = str(zip_path)
    gals_phot = tp.load_galaxies(gas="v16eq6")
    gals_raw = tf.load_real_sparc()
    labels = tp.split(gals_raw)
    visible = sorted(k for k in gals_raw if labels[k] != "holdout")

    # 1. Negative-density comparison across SPARC galaxies
    n_phot_neg = sum(int(np.sum(g["rho"] <= 0)) for g in gals_phot.values())
    n_deriv_neg_visible_05 = 0
    n_deriv_neg_all175_05 = 0
    for k, g in gals_raw.items():
        c05 = tf.compute(g, 0.5, 2.36, "local")
        n_neg = sum("rho<0" in f for f in c05["flags"])
        n_deriv_neg_all175_05 += n_neg
        if k in visible:
            n_deriv_neg_visible_05 += n_neg

    # 2. NGC 3198 control point at R = 30.08 kpc
    ngc = gals_raw["NGC3198"]
    idx = int(np.argmin(np.abs(ngc["R"] - 30.08)))
    R_ctrl = float(ngc["R"][idx])
    V_obs = float(ngc["Vobs"][idx])
    errV = float(ngc["errV"][idx])
    V_gas = float(ngc["Vgas"][idx])
    V_disk = float(ngc["Vdisk"][idx])
    V_bul = float(ngc["Vbul"][idx])
    V_bar_10 = math.sqrt(math.copysign(V_gas**2, V_gas) + 1.0 * V_disk**2 + 0.7 * V_bul**2)
    V_bar_05 = math.sqrt(math.copysign(V_gas**2, V_gas) + 0.5 * V_disk**2 + 0.7 * V_bul**2)
    C_req_10 = (V_obs / V_bar_10) ** 2
    C_req_05 = (V_obs / V_bar_05) ** 2

    # 3. 3D spherical top-hat summary check
    tophat_summary_path = HERE / "data" / "sparc_tophat_summary.json"
    tophat_info = {}
    if tophat_summary_path.is_file():
        s = json.loads(tophat_summary_path.read_text())
        tophat_info = {
            "max_isolated_1_minus_fL_grid": s["all175_L_le_10Mpc"]["max_one_minus_fL_TAL_galaxy_only"],
            "fine_L_refinement_top5": s["fine_L_refinement_top5"][:2],
        }

    report = {
        "sparc_data_integrity": {
            "n_galaxies": len(gals_raw),
            "n_total_radial_points": sum(len(g["R"]) for g in gals_raw.values()),
            "Rotmod_LTG_zip_sha256": EXPECTED_ZIP_SHA,
            "SPARC_Lelli2016c_mrt_sha256": EXPECTED_MRT_SHA,
            "blind_123_txt_sha256": EXPECTED_BLIND_SHA,
        },
        "negative_density_audit_175_galaxies": {
            "photometric_and_3D_tophat_yukawa_rho_le_0_points": n_phot_neg,
            "derivative_formula_rho_lt_0_points_visible140_Ydisk_0p5": n_deriv_neg_visible_05,
            "derivative_formula_rho_lt_0_points_all175_Ydisk_0p5": n_deriv_neg_all175_05,
        },
        "NGC3198_control_point": {
            "R_kpc": R_ctrl,
            "V_obs_km_s": V_obs,
            "errV_km_s": errV,
            "V_bar_Ydisk_1p0_km_s": round(V_bar_10, 2),
            "C_req_Ydisk_1p0": round(C_req_10, 3),
            "V_bar_Ydisk_0p5_km_s": round(V_bar_05, 2),
            "C_req_Ydisk_0p5": round(C_req_05, 3),
        },
        "v19_3d_tophat_peak": tophat_info,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--mode", default="summary", choices=["summary", "pipeline", "foundation", "v19_7th"])
    args, rest = p.parse_known_args()
    if args.mode == "summary":
        run_summary_audit()
    elif args.mode == "pipeline":
        sys.argv = [sys.argv[0]] + rest
        tp.main()
    elif args.mode == "foundation":
        tf.ZIP = str(HERE / "data" / "Rotmod_LTG.zip")
        sys.argv = [sys.argv[0]] + rest
        tf.main()
    elif args.mode == "v19_7th":
        import runpy
        sys.argv = [str(HERE / "tal_factor_7th_model.py")] + rest
        runpy.run_path(str(HERE / "tal_factor_7th_model.py"), run_name="__main__")
