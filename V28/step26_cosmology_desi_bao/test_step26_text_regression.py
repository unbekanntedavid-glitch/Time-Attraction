#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 26 (V27) text regression test (review T01).

Checks the generated HU/EN Markdown reports (and their v27_* copies), the local-benchmark CSV and the
RSD figure source for the text items F06, F07, F08, N01 and N02:
  F06  correct sentence on the 2 alpha_Y^2 beta_Yuk row (separate double-strength sensitivity control,
       not the growth law of the published V26 Step25); no V26 factor-2 growth-law attribution
  F08  identity 5 log10(72.8/67.4) for identical E(z), node difference labelled a model difference;
       0.167357 mag never attributed to the lowest node
  N01  no bulk value in the "R >> L edge" column of the FLRW, void and Milky Way disk rows (n.a. / null)
  N02  Section 1.7 sentence separating Omega_phi ~ 1e-47 from the direct RSD Yukawa term (numbers from JSON)
  F07  canonical maximum string in the RSD figure legend source
Exit status 0 if every check passes, 1 otherwise. Usage: python3 test_step26_text_regression.py [directory]
"""
import os
import re
import sys
import csv
import json
import math

D = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
FAIL = []


def check(cond, msg):
    print(("  PASS " if cond else "  FAIL ") + msg)
    if not cond:
        FAIL.append(msg)


def read(name):
    with open(os.path.join(D, name), encoding="utf-8") as fh:
        return fh.read()


def sci(x, n):
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(abs(m), n) >= 10.0:
        e += 1
        m = x / 10 ** e
    return f"{m:.{n}f}\\times 10^{{{e}}}"


J = json.loads(read("tal_hubble_desi_bao_summary.json"))
R4 = J["step1_4_rsd_fsigma8_growth"]
N02_DG = sci(R4["delta_Geff_over_GN_z0_canonical"], 1)
N02_DFS = sci(R4["max_abs_delta_fsigma8_variational_canonical"], 2)
DMU = f"{5 * math.log10(72.8 / 67.4):.6f}"

REQUIRED = {
    "EN": [
        ("F06", r"The $2\alpha_Y^2\beta_{\rm Yuk}$ row is a separate, double-strength sensitivity control; "
                r"it is not the growth law of the published V26 Step25 (Section 5)."),
        ("F08", r"For identical $E(z)$ geometry the offset would be the constant identity $5\log_{10}(72.8/67.4)=+" + DMU),
        ("F08", r"is a model difference (the reference has $\Omega_m=0.334$), not this identity."),
        ("N02", r"(the background scalar-energy correction is $\Omega_\phi\sim10^{-47}$; the direct Yukawa perturbation "
                r"of RSD is a separate negligible term quantified in Section 5: $\delta G_{\rm eff}/G_*\approx" + N02_DG +
                r"$, $\delta f\sigma_8\approx" + N02_DFS + r"$)"),
    ],
    "HU": [
        ("F06", r"A $2\alpha_Y^2\beta_{\rm Yuk}$ sor külön, kétszeres erősségű érzékenységi kontroll; "
                r"nem a publikált V26 Step25 növekedési törvénye (5. fejezet)."),
        ("F08", r"Azonos $E(z)$ geometria esetén az eltolás a konstans $5\log_{10}(72.8/67.4)=+" + DMU),
        ("F08", r"modellkülönbség (a referencia $\Omega_m=0.334$), nem ez az azonosság."),
        ("N02", r"(a háttér skalárenergia-korrekciója $\Omega_\phi\sim10^{-47}$; az RSD közvetlen Yukawa-perturbációja "
                r"külön, az 5. fejezetben számszerűsített elhanyagolható tag: $\delta G_{\rm eff}/G_*\approx" + N02_DG +
                r"$, $\delta f\sigma_8\approx" + N02_DFS + r"$)"),
    ],
}

# superseded wrong sentences (verbatim fragments of earlier texts) and generic patterns
FORBIDDEN_FRAGMENTS = [
    ("F06", "factor-2 convention of the earlier V26 growth equation"),
    ("F06", "korábbi V26 növekedési egyenlet kétszeres konvenció"),
    ("F08", "at the lowest node is $+0.167357"),
    ("F08", "a legalsó csomóponton $+0.167357"),
]
FORBIDDEN_REGEX = [
    ("F08", re.compile(r"(lowest[- ]node|first node|legalsó csomópont|első csomópont)[^.\n]{0,80}0\.167357")),
    ("F08", re.compile(r"0\.167357[^.\n]{0,80}(lowest[- ]node|first node|legalsó csomópont|első csomópont)")),
]
NEG_EN = "not the growth law of the published V26 Step25"
NEG_HU = "nem a publikált V26 Step25 növekedési törvénye"


def sentences(txt):
    return re.split(r"(?<=[.;])\s+|\n", txt)


def table_rows(txt, labels):
    out = {}
    for line in txt.splitlines():
        if not line.startswith("| "):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells and cells[0] in labels:
            out[cells[0]] = cells
    return out


N01_LABELS = {"EN": ["Homogeneous FLRW, $z=0$", r"Cosmic void ($0.2\rho_{m,0}$)", "Milky Way disk (average)"],
              "HU": ["Homogén FLRW, $z=0$", r"Kozmikus üreg ($0.2\rho_{m,0}$)", "Tejútrendszer-korong (átlag)"]}
N01_IDS = ["1_Cosmic_Background_z0", "2_Local_Cosmic_Void", "3_Milky_Way_Disk_Solar_Nbr"]
HALF_LABELS = {"EN": ["Mean solar density (bulk)", "Mean Earth density (bulk)", "White dwarf interior (bulk)"],
               "HU": ["Nap átlagsűrűség (bulk)", "Föld átlagsűrűség (bulk)", "Fehér törpe belseje (bulk)"]}


def check_md(name, lang):
    print(f"[{name}]")
    txt = read(name)
    for item, s in REQUIRED[lang]:
        check(s in txt, f"{item} required sentence present ({lang})")
    for item, s in FORBIDDEN_FRAGMENTS:
        check(s not in txt, f"{item} wrong sentence absent: '{s[:50]}'")
    for item, rx in FORBIDDEN_REGEX:
        m = rx.search(txt)
        check(m is None, f"{item} 0.167357 not attributed to a node" + (f" (found: '{m.group(0)[:70]}')" if m else ""))
    # F06 generic: every sentence linking the 2 alpha^2 beta row to a V26 growth law must carry the negation
    bad = [s for s in sentences(txt)
           if "V26" in s and re.search(r"growth (law|equation)|növekedési (törvény|egyenlet)", s)
           and r"2\alpha" in s and NEG_EN not in s and NEG_HU not in s]
    check(not bad, "F06 no V26 factor-2 growth-law attribution" + (f" (found: '{bad[0][:90]}')" if bad else ""))
    # N01: section 7 table
    rows = table_rows(txt, N01_LABELS[lang] + HALF_LABELS[lang])
    for lab in N01_LABELS[lang]:
        r = rows.get(lab)
        check(r is not None and len(r) >= 4 and r[3] == "n.a." and r[2] != "n.a.",
              f"N01 row '{lab}': edge column n.a., bulk column present" + (f" (edge cell: {r[3]})" if r and len(r) >= 4 else ""))
    for lab in HALF_LABELS[lang]:
        r = rows.get(lab)
        check(r is not None and len(r) >= 4 and r[3] != "n.a.", f"N01 row '{lab}': half-space edge value present")


def check_csv():
    print("[step1_6_local_yukawa_environments.csv]")
    with open(os.path.join(D, "step1_6_local_yukawa_environments.csv"), newline="", encoding="utf-8") as fh:
        rows = {r["env_id"]: r for r in csv.DictReader(fh)}
    for i in N01_IDS:
        r = rows[i]
        check(r["abs_phi_surface_R_gg_L"] == "" and r["Delta_H0_apparent_surf_km_s_Mpc"] == ""
              and r["abs_phi_bulk"] != "" and r.get("boundary_status") == "n.a. (no boundary geometry)",
              f"N01 CSV {i}: edge and edge-Delta-H0 empty (null), bulk kept, boundary_status n.a.")
    jrows = {r["env_id"]: r for r in J["step1_6_local_yukawa_bulk_benchmarks"]}
    for i in N01_IDS:
        r = jrows[i]
        check(r["abs_phi_surface_R_gg_L"] is None and r["Delta_H0_apparent_surf_km_s_Mpc"] is None
              and r.get("boundary_status") == "n.a. (no boundary geometry)", f"N01 JSON {i}: null edge values, status n.a.")


def check_rsd_source():
    print("[run_step1_4_to_1_6.py RSD legend source]")
    src = read("run_step1_4_to_1_6.py")
    canon = "max|delta f*sigma8|_var,can = {max_var_can:.3e}"
    check(src.count(canon) >= 2, "F07 canonical legend string in HU and EN label")
    check("_var,2x" not in src.split("def generate_rsd_plots", 1)[1].split("\ndef ", 1)[0],
          "F07 no 2x value in the legend")
    check(len(re.findall(r"= generate_rsd_plots\(df_14, chi2_m4, chi2_low, max_var_can, lang=", src)) == 2
          and 'max_var_can = float(df_14["delta_fsigma8_var_can_z0norm"].max())' in src,
          "F07 legend receives the canonical column maximum")


def main():
    print(f"=== Step 26 text regression test (T01) in {D} ===")
    for name, lang in [("step26_desi_bao_planck_sne_rsd_bbn_nogo_en.md", "EN"), ("step26_desi_bao_planck_sne_rsd_bbn_nogo_hu.md", "HU"),
                       ("v27_cosmology_desi_planck_sne_rsd_audit_en.md", "EN"), ("v27_cosmology_desi_planck_sne_rsd_audit_hu.md", "HU")]:
        check_md(name, lang)
    check_csv()
    check_rsd_source()
    if FAIL:
        print(f"text regression: FAILED ({len(FAIL)} checks)")
        sys.exit(1)
    print("text regression: all checks passed")


if __name__ == "__main__":
    main()
