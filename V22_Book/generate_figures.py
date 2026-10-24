import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 9.5,
    "axes.titlesize": 10.0,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0,
})

# Load SPARC data
with open("/workspace/artifacts/sparc_tophat/sparc_tophat_per_galaxy_max.csv", "r", encoding="utf-8") as f:
    pg_rows = [r for r in csv.DictReader(filter(lambda l: not l.startswith("#"), f))]

with open("/workspace/artifacts/academia_review_files/fitted_L_scores.csv", "r", encoding="utf-8") as f:
    sc_rows = list(csv.DictReader(f))

L_b123 = [float(r["L_at_max_fL_kpc"]) for r in pg_rows if int(r["in_blind123"]) == 1]
f_b123 = [float(r["max_one_minus_fL_TAL_abs"]) for r in pg_rows if int(r["in_blind123"]) == 1]
L_oth  = [float(r["L_at_max_fL_kpc"]) for r in pg_rows if int(r["in_blind123"]) == 0]
f_oth  = [float(r["max_one_minus_fL_TAL_abs"]) for r in pg_rows if int(r["in_blind123"]) == 0]

top_row = max(pg_rows, key=lambda r: float(r["max_one_minus_fL_TAL_abs"]))
L_top = float(top_row["L_at_max_fL_kpc"])
f_top = float(top_row["max_one_minus_fL_TAL_abs"])

vbar_vals = [float(r["vbar"]) for r in sc_rows]
true_vals = [float(r["true"]) for r in sc_rows]
prox_vals = [float(r["proxy"]) for r in sc_rows]

# Load TOV M-R sequences
with open("/workspace/artifacts/ns_tal_V21_2026-10-23/mr_sequences.csv", "r", encoding="utf-8") as f:
    mr_rows = list(csv.DictReader(f))

# Load radial stability grid
with open("/workspace/artifacts/ns_tal_open_items_2026-10-23/outputs/item2_radial_grid.csv", "r", encoding="utf-8") as f:
    rad_rows = list(csv.DictReader(f))

# Load 8-coil lab differential vs L
with open("/workspace/artifacts/final_report/csv/tal_differential_vs_L.csv", "r", encoding="utf-8") as f:
    lab_rows = list(csv.DictReader(f))

L_lab = [float(r["L (m)"]) for r in lab_rows]
dy_gL = [abs(float(r["df_f_diff_signed γ_L, 14 mA"])) for r in lab_rows]
dy_fx = [abs(float(r["df_f_diff_signed γ=1e-22, 14 mA"])) for r in lab_rows]

LANGS = {
    "en": {
        "sparc_oth": "Other SPARC galaxies (n=52)",
        "sparc_b123": "Blind-123 subsample (n=123)",
        "sparc_max": "Max: UGC03546 (2.28e-7 grid / 2.303e-7 fine scan)",
        "sparc_xl1": "Smoothing scale L at galaxy maximum [kpc]",
        "sparc_yl1": "Maximum galaxy-attributable |1 - f_L|",
        "sparc_t1": "(a) 3D Top-Hat Peak |1 - f_L| (sparc_tophat.zip, N=175)",
        "sparc_prox": "1D midplane proxy (N=175 med.=2.672; Blind-123=2.708)",
        "sparc_vbar": "Baryonic baseline vbar (N=175 med.=2.898; Blind-123=3.143)",
        "sparc_true": "3D spherical top-hat (N=175 med.=2.898; Blind-123=3.143)",
        "sparc_xl2": "Reduced chi-squared per galaxy (fitted_L_scores.csv)",
        "sparc_yl2": "Number of galaxies (out of 175)",
        "sparc_t2": "(b) Reduced Chi-Squared at L = 2.3714 Mpc",
        "tov_t1": "(a) Conditional Mass-Source Mode (SLy, APR4, H4)",
        "tov_t2": "(b) Conditional Gradient-Scaling Mode (SLy, APR4, H4)",
        "tov_xl": "Areal radius R [km]",
        "tov_yl": "Gravitational mass M [M_sun]",
        "tov_gr": "GR (f=1)",
        "tov_ac": "acausal (c_s > c)",
        "stab_t1": "(a) Frozen-f 'matter' Operator (SLy, APR4, H4)",
        "stab_t2": "(b) Frozen-f 'metric' Operator (x = indefinite W<=0)",
        "stab_xl": "log10(central rest-mass density rho_0c [g/cm^3])",
        "stab_yl": "Signed fundamental frequency sign(omega_0^2) f_0 [kHz]",
        "stab_ac": "APR4 acausal (c_s>c)",
        "lab_gL": r"$|\Delta f/f_{\mathrm{diff}}|$, $\gamma_L = G L^2 / c^2$ (14 mA/coil)",
        "lab_fx": r"$|\Delta f/f_{\mathrm{diff}}|$, $\gamma = 10^{-22}\,\mathrm{m^3\,kg^{-1}}$ (reference)",
        "lab_pk": r"Sampled 25-pt grid max: $-1.2017\times 10^{-48}$ at $L = 0.631\,\mathrm{m}$",
        "lab_xl": "Smoothing radius L on 25-point logarithmic grid [m]",
        "lab_yl": "Differential frequency shift |Delta f / f|",
        "lab_t": "8-Coil Laboratory Experiment: Conditional Differential Shift vs. L",
    },
    "hu": {
        "sparc_oth": "További SPARC-galaxisok (n=52)",
        "sparc_b123": "Vak 123-as részhalmaz (n=123)",
        "sparc_max": "Max.: UGC03546 (2,28e-7 rács / 2,303e-7 finomított)",
        "sparc_xl1": "Simítási skála L a galaxis-maximumnál [kpc]",
        "sparc_yl1": "Maximális galaktikus |1 - f_L|",
        "sparc_t1": "(a) 3D top-hat maximum |1 - f_L| (sparc_tophat.zip, N=175)",
        "sparc_prox": "1D középsík-proxy (N=175 med.=2,672; Vak-123=2,708)",
        "sparc_vbar": "Barionikus alap (vbar, N=175 med.=2,898; Vak-123=3,143)",
        "sparc_true": "Valódi 3D top-hat (N=175 med.=2,898; Vak-123=3,143)",
        "sparc_xl2": "Redukált khi-négyzet galaxisonként (fitted_L_scores.csv)",
        "sparc_yl2": "Galaxisok száma (175-ből)",
        "sparc_t2": "(b) Redukált khi-négyzet L = 2,3714 Mpc-nél",
        "tov_t1": "(a) Feltételes tömegforrás-mód (SLy, APR4, H4)",
        "tov_t2": "(b) Feltételes gradiens-mód (SLy, APR4, H4)",
        "tov_xl": "Areális sugár R [km]",
        "tov_yl": "Gravitációs tömeg M [M_nap]",
        "tov_gr": "GR (f=1)",
        "tov_ac": "akauzális (c_s > c)",
        "stab_t1": "(a) Befagyasztott-f 'matter' operátor (SLy, APR4, H4)",
        "stab_t2": "(b) Befagyasztott-f 'metric' operátor (x = indefinit W<=0)",
        "stab_xl": "log10(központi nyugalmi sűrűség rho_0c [g/cm^3])",
        "stab_yl": "Előjeles alapfrekvencia sign(omega_0^2) f_0 [kHz]",
        "stab_ac": "APR4 akauzális (c_s>c)",
        "lab_gL": r"$|\Delta f/f_{\mathrm{diff}}|$, $\gamma_L = G L^2 / c^2$ (14 mA/tekercs)",
        "lab_fx": r"$|\Delta f/f_{\mathrm{diff}}|$, $\gamma = 10^{-22}\,\mathrm{m^3\,kg^{-1}}$ (referencia)",
        "lab_pk": r"Maximum a 25 pontos L-rácson: $-1{,}2017\times 10^{-48}$ ($L = 0{,}631\,\mathrm{m}$)",
        "lab_xl": "Simítási sugár L a 25 pontos logaritmikus rácson [m]",
        "lab_yl": "Differenciális frekvenciaeltolódás |Delta f / f|",
        "lab_t": "8 tekercses laborkísérlet: feltételes differenciális jel L függvényében",
    },
    "de": {
        "sparc_oth": "Weitere SPARC-Galaxien (n=52)",
        "sparc_b123": "Blinde 123er-Teilstichprobe (n=123)",
        "sparc_max": "Max.: UGC03546 (2,28e-7 Gitter / 2,303e-7 verfeinert)",
        "sparc_xl1": "Glättungsradius L beim Galaxien-Maximum [kpc]",
        "sparc_yl1": "Maximales galaxieneigenes |1 - f_L|",
        "sparc_t1": "(a) 3D-Top-Hat-Maximum |1 - f_L| (sparc_tophat.zip, N=175)",
        "sparc_prox": "1D-Mittelebenen-Proxy (N=175 Med.=2,672; Blind-123=2,708)",
        "sparc_vbar": "Baryonische Basis (vbar, N=175 Med.=2,898; Blind-123=3,143)",
        "sparc_true": "Echtes 3D-Kugelmittel (N=175 Med.=2,898; Blind-123=3,143)",
        "sparc_xl2": "Reduziertes Chi-Quadrat pro Galaxie (fitted_L_scores.csv)",
        "sparc_yl2": "Anzahl der Galaxien (von 175)",
        "sparc_t2": "(b) Reduziertes Chi-Quadrat bei L = 2,3714 Mpc",
        "tov_t1": "(a) Bedingter Massenquellen-Modus (SLy, APR4, H4)",
        "tov_t2": "(b) Bedingter Gradienten-Modus (SLy, APR4, H4)",
        "tov_xl": "Arealer Radius R [km]",
        "tov_yl": "Gravitationsmasse M [M_Sonne]",
        "tov_gr": "ART (f=1)",
        "tov_ac": "akausal (c_s > c)",
        "stab_t1": "(a) Eingefrorener-f 'matter'-Operator (SLy, APR4, H4)",
        "stab_t2": "(b) Eingefrorener-f 'metric'-Operator (x = indefinit W<=0)",
        "stab_xl": "log10(zentrale Ruhemasseedichte rho_0c [g/cm^3])",
        "stab_yl": "Vorzeichenbehaftete Grundfrequenz sign(omega_0^2) f_0 [kHz]",
        "stab_ac": "APR4 akausal (c_s>c)",
        "lab_gL": r"$|\Delta f/f_{\mathrm{diff}}|$, $\gamma_L = G L^2 / c^2$ (14 mA/Spule)",
        "lab_fx": r"$|\Delta f/f_{\mathrm{diff}}|$, $\gamma = 10^{-22}\,\mathrm{m^3\,kg^{-1}}$ (Referenz)",
        "lab_pk": r"Maximum auf 25-Punkte-L-Gitter: $-1{,}2017\times 10^{-48}$ ($L = 0{,}631\,\mathrm{m}$)",
        "lab_xl": "Glättungsradius L auf dem logarithmischen 25-Punkte-Gitter [m]",
        "lab_yl": "Differenzielle Frequenzverschiebung |Delta f / f|",
        "lab_t": "8-Spulen-Laborexperiment: Bedingte differenzielle Verschiebung vs. L",
    }
}

colors_eos = {"SLy": "#1f77b4", "APR4": "#d62728", "H4": "#2ca02c"}

for lang, T in LANGS.items():
    LANG_UP = lang.upper()
    # 1. SPARC figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 3.2), dpi=220)
    ax1.scatter(L_oth, f_oth, s=20, color="#888888", alpha=0.75, label=T["sparc_oth"])
    ax1.scatter(L_b123, f_b123, s=22, color="#1f77b4", alpha=0.85, label=T["sparc_b123"])
    ax1.scatter([L_top], [f_top], s=65, color="#d62728", zorder=5, label=T["sparc_max"])
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(T["sparc_xl1"])
    ax1.set_ylabel(T["sparc_yl1"])
    ax1.set_title(T["sparc_t1"], fontweight="bold")
    ax1.grid(True, which="both", linestyle=":", alpha=0.5)
    ax1.legend(loc="lower right", fontsize=7.5, framealpha=0.92)

    bins = np.logspace(-1.5, 2.2, 32)
    ax2.hist(prox_vals, bins=bins, color="#87ceeb", alpha=0.55, label=T["sparc_prox"])
    ax2.hist(vbar_vals, bins=bins, histtype="step", color="black", linewidth=1.6, label=T["sparc_vbar"])
    ax2.hist(true_vals, bins=bins, histtype="step", color="#d62728", linestyle="--", linewidth=1.6, label=T["sparc_true"])
    ax2.set_xscale("log")
    ax2.set_xlabel(T["sparc_xl2"])
    ax2.set_ylabel(T["sparc_yl2"])
    ax2.set_title(T["sparc_t2"], fontweight="bold")
    ax2.grid(True, which="both", linestyle=":", alpha=0.5)
    ax2.legend(loc="upper right", fontsize=7.3, framealpha=0.92)

    plt.tight_layout()
    fig.savefig(f"/workspace/artifacts/tal_booklet_2026-10-24/figures/sparc_tophat_{LANG_UP}.pdf")
    fig.savefig(f"/workspace/artifacts/tal_booklet_2026-10-24/figures/fig1_sparc_{lang}.png", dpi=220)
    plt.close(fig)

    # 2. TOV M-R sequences figure (SLy, APR4, H4 across mass and gradient modes, showing causal & acausal parts + NICER bands)
    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(10.2, 3.5), dpi=220, sharey=True)
    for mode_key, ax_curr, title_str in [
        ("mass", ax2a, T["tov_t1"]),
        ("gradient", ax2b, T["tov_t2"])
    ]:
        # NICER 1-sigma bands (J0740, J0437, J0030)
        ax_curr.axhspan(2.073 - 0.069, 2.073 + 0.069, color="#ff9896", alpha=0.22, label="PSR J0740+6620" if mode_key == "mass" else None)
        ax_curr.axhspan(1.418 - 0.037, 1.418 + 0.037, color="#aec7e8", alpha=0.25, label="PSR J0437-4715" if mode_key == "mass" else None)
        for eos in ["SLy", "APR4", "H4"]:
            gr_caus = [r for r in mr_rows if r["EOS"] == eos and r["mode"] == "GR" and r["acausal_point"] == "False"]
            gr_acaus = [r for r in mr_rows if r["EOS"] == eos and r["mode"] == "GR" and r["acausal_point"] == "True"]
            ax_curr.plot([float(r["R_km"]) for r in gr_caus], [float(r["M_Msun"]) for r in gr_caus],
                         linestyle="--", color=colors_eos[eos], linewidth=1.6, label=f"{eos} {T['tov_gr']}")
            if gr_acaus:
                ax_curr.plot([float(r["R_km"]) for r in gr_acaus], [float(r["M_Msun"]) for r in gr_acaus],
                             linestyle=":", color=colors_eos[eos], linewidth=1.2, alpha=0.65)
            tal_caus = [r for r in mr_rows if r["EOS"] == eos and r["mode"] == mode_key and abs(float(r["L_m"]) - 10000.0) < 1.0 and r["acausal_point"] == "False"]
            tal_acaus = [r for r in mr_rows if r["EOS"] == eos and r["mode"] == mode_key and abs(float(r["L_m"]) - 10000.0) < 1.0 and r["acausal_point"] == "True"]
            ax_curr.plot([float(r["R_km"]) for r in tal_caus], [float(r["M_Msun"]) for r in tal_caus],
                         linestyle="-", color=colors_eos[eos], linewidth=1.8, label=f"{eos} (L=10 km)")
            if tal_acaus:
                ax_curr.plot([float(r["R_km"]) for r in tal_acaus], [float(r["M_Msun"]) for r in tal_acaus],
                             linestyle=":", color=colors_eos[eos], linewidth=1.4, alpha=0.75,
                             label=f"{eos} {T['tov_ac']}" if eos == "APR4" else None)
        ax_curr.set_xlabel(T["tov_xl"])
        ax_curr.set_xlim(9.2, 15.2)
        ax_curr.set_ylim(0.5, 2.45)
        ax_curr.set_title(title_str, fontweight="bold")
        ax_curr.grid(True, linestyle=":", alpha=0.5)
        ax_curr.legend(loc="lower right", fontsize=7.2, ncol=2, framealpha=0.92)
    ax2a.set_ylabel(T["tov_yl"])
    plt.tight_layout()
    fig2.savefig(f"/workspace/artifacts/tal_booklet_2026-10-24/figures/fig2_tov_{lang}.png", dpi=220)
    plt.close(fig2)

    # 3. Radial stability figure (2-panel: matter vs metric across SLy, APR4, H4 and mass/gradient modes at L=10 km)
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(10.2, 3.45), dpi=220, sharey=True)
    for var_key, ax_curr, title_str in [
        ("matter", ax3a, T["stab_t1"]),
        ("metric", ax3b, T["stab_t2"])
    ]:
        ax_curr.axvspan(15.1601, 15.65, color="#ffbb78", alpha=0.18, label=T["stab_ac"])
        for eos in ["SLy", "APR4", "H4"]:
            gr_r = [r for r in rad_rows if r["EOS"] == eos and r["mode"] == "GR" and r["variant"] == var_key]
            ax_curr.plot([float(r["log10_rho0c"]) for r in gr_r], [float(r["f0_kHz"]) for r in gr_r],
                         linestyle="--", color=colors_eos[eos], linewidth=1.4, label=f"{eos} GR")
            m_r = [r for r in rad_rows if r["EOS"] == eos and r["mode"] == "mass" and r["variant"] == var_key and abs(float(r["L_m"]) - 10000.0) < 1.0]
            ax_curr.plot([float(r["log10_rho0c"]) for r in m_r], [float(r["f0_kHz"]) for r in m_r],
                         linestyle="-", color=colors_eos[eos], linewidth=1.6, label=f"{eos} mass (10 km)")
            g_r = [r for r in rad_rows if r["EOS"] == eos and r["mode"] == "gradient" and r["variant"] == var_key and abs(float(r["L_m"]) - 10000.0) < 1.0]
            xs_g = [float(r["log10_rho0c"]) for r in g_r]
            ys_g = [float(r["f0_kHz"]) for r in g_r]
            ax_curr.plot(xs_g, ys_g, linestyle="-.", color=colors_eos[eos], linewidth=1.4,
                         marker="x" if var_key == "metric" else None, markersize=3.5, markevery=4,
                         label=f"{eos} grad (10 km)")
        ax_curr.axhline(0.0, color="black", linestyle=":", linewidth=1.1)
        ax_curr.set_xlim(14.55, 15.60)
        ax_curr.set_ylim(-1.6, 4.2)
        ax_curr.set_xlabel(T["stab_xl"])
        ax_curr.set_title(title_str, fontweight="bold")
        ax_curr.grid(True, linestyle=":", alpha=0.5)
        ax_curr.legend(loc="lower left", fontsize=6.8, ncol=3, framealpha=0.92)
    ax3a.set_ylabel(T["stab_yl"])
    plt.tight_layout()
    fig3.savefig(f"/workspace/artifacts/tal_booklet_2026-10-24/figures/fig3_stability_{lang}.png", dpi=220)
    plt.close(fig3)

    # 4. 8-coil laboratory signal vs L figure
    fig4, ax4 = plt.subplots(figsize=(8.4, 3.1), dpi=220)
    ax4.plot(L_lab, dy_gL, color="#12284C", linewidth=2.1, marker="o", markersize=3.8, label=T["lab_gL"])
    ax4.plot(L_lab, dy_fx, color="#777777", linestyle="--", linewidth=1.5, label=T["lab_fx"])
    ax4.scatter([0.630957], [1.2017108e-48], color="#d62728", s=60, zorder=5, label=T["lab_pk"])
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.set_xlabel(T["lab_xl"])
    ax4.set_ylabel(T["lab_yl"])
    ax4.set_title(T["lab_t"], fontweight="bold")
    ax4.grid(True, which="both", linestyle=":", alpha=0.5)
    ax4.legend(loc="lower left", fontsize=8.0, framealpha=0.92)
    plt.tight_layout()
    fig4.savefig(f"/workspace/artifacts/tal_booklet_2026-10-24/figures/fig4_lab_{lang}.png", dpi=220)
    fig4.savefig(f"/workspace/artifacts/tal_booklet_2026-10-24/figures/fig1_signal_vs_L_{LANG_UP}.pdf")
    plt.close(fig4)

print("Generated all 12 PNG figures and PDF figures for EN, HU, DE successfully.")
