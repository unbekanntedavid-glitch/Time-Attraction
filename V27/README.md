# Time Attraction Law — Scalar-Tensor Yukawa Fit (V27): Step 26 Cosmological Audit / 26. lépés kozmológiai audit

**Author / Szerző:** Dávid Szabolcsi (`ORCID: 0009-0005-6300-074X`)  
**Date / Dátum:** 2026-10-28 (revised after the Academia Co-scientist review, the re-audit and the second review F01–F12 and the third re-audit N01, N02, T01 / a bírálat, az újrabírálat, az F01–F12 második bírálat és a harmadik újraaudit N01, N02, T01 után átdolgozva: 2026-10-29)

**Preceding DOIs / Előzmény DOI-k:**

- V24: `10.5281/zenodo.23680715` (strictly V24 / szigorúan V24)
- V25: `10.5281/zenodo.23699685` (strictly V25, `step1`–`step22`)
- V26: `10.5281/zenodo.23727585` (strictly V26, `step1`–`step25`)

## Reports / Jelentések

- EN: `step26_desi_bao_planck_sne_rsd_bbn_nogo_en.md` / `.pdf` (copy: `v27_cosmology_desi_planck_sne_rsd_audit_en.md`)
- HU: `step26_desi_bao_planck_sne_rsd_bbn_nogo_hu.md` / `.pdf` (copy: `v27_cosmology_desi_planck_sne_rsd_audit_hu.md`)

Computed results are generated from the JSON/CSV sources (`tal_hubble_desi_bao_summary.json`, Step 26 CSVs) by `generate_step26_markdown.py`; fixed input constants are stated in the code. / A számítási eredményeket a JSON/CSV-forrásokból generáljuk (`tal_hubble_desi_bao_summary.json`, a 26. lépés CSV-i) a `generate_step26_markdown.py` segítségével; a rögzített bemeneti állandók a kódban szerepelnek.

## Key results / Fő eredmények (branch: adiabatic attractor, $C=0$)

- Unified normalization: $\gamma_L=GL^2/c^2=7.42616\times 10^{-20}$ m^3/kg, $\bar\phi_{\rm attr}(0)=-1.9949\times 10^{-46}$, $\Omega_{\phi,\rm attr}(0)=3.140\times 10^{-47}$, $E^2(0)=1$ (DE closure).
- DESI DR1 BAO (12×12): M4 $\chi^2=20.7852$ ($p=0.0536$); rigid 72.8 $\chi^2=194.342$; 2-param profile $\Omega_m=0.2939$, $P=10193.67$ km/s, $\chi^2=12.7405$; $r_d$ needed at 72.8: $140.02\pm1.74$ Mpc; $H_0$ with Planck $r_d$: $69.30\pm0.87$; joint $\chi^2=15.67431585+0.06234684+3.97660606=19.71326876$ (BAO + $r_d$ prior + $H_0$ prior), $\Delta\chi^2=6.973$.
- Planck acoustic scale (re-integrated $r_s$, $r_d$): M4 $100\theta_*=1.041320$.
- Late transition (piecewise $H=sH_{\rm M4}$ for $z<z_t$): see Section 3 table; the $D_M/s$ value $100\theta_*=1.124750$ is an inconsistent illustration.
- Synthetic SNe shape control (not a catalogue likelihood), symmetric offsets: M4 $0.17595029$ mag, $\chi^2$ $522.00704\to0.916614$; low-$\Omega_m$ $0.04159799$ mag, $\chi^2$ $51.90332\to22.777527$.
- RSD diagonal benchmark $\chi^2=10.1233$; variational $\max|\delta f\sigma_8|=1.666\times 10^{-40}$.
- Algebraic BBN benchmark ($10^5\,{\rm D/H}^{\rm std}=2.5147$ is an arithmetic proxy); bulk local Yukawa benchmarks with `expm1`; the V21 top-hat interior external benchmark $|\phi|=0.0412$ is an external input.
- M3 control with its own horizons: $r_s=144.36622618$ Mpc, $r_d=147.01159289$ Mpc, $100\theta_*=0.41346227$, $\chi^2_{\rm BAO}=2201.29927$.
- $72.8\pm1.0$ km/s/Mpc is an assumed local $H_0$ benchmark prior (not a sourced SH0ES value). / A $72.8\pm1.0$ km/s/Mpc feltételezett lokális $H_0$ benchmark-prior.
- Second review F01–F12: Revision log 2 (Section 12 of the reports); open limitations P01–P09: Section 9.1; third re-audit N01, N02, T01: Revision log 3 (Section 13); text regression test: `test_step26_text_regression.py` (run by `run_all.py`). / Második bírálat: Revíziós napló 2 (12. fejezet); nyitott korlátok: 9.1 fejezet; harmadik újraaudit: Revíziós napló 3 (13. fejezet); szöveges regressziós teszt: `test_step26_text_regression.py`.

## Build / Fordítás

Dependencies / Függőségek: `pip install -r requirements.txt` (numpy, scipy, pandas, matplotlib); PDF: pandoc + pdflatex.

Recommended / Ajánlott: `python3 run_all.py` (full run including Block 1, a few seconds / az 1. blokkal együtt, néhány másodperc) or `python3 run_all.py --skip-block1` — checks the SHA-256 of the DESI input files, then runs in order / ellenőrzi a DESI-bemenetek SHA-256 értékét, majd sorban futtatja:

1. `run_step2_nogo_and_report.py` (M1–M4 branches, inverse ladder CSV and figures)
2. `run_step26_revision_controls.py` (re-audit controls R1–R6, 5.1–5.3; merges into the strict JSON)
3. `generate_step26_markdown.py` (Markdown from JSON/CSV)
   then `test_step26_text_regression.py` (HU/EN text regression test, review T01; nonzero exit on failure / szöveges regressziós teszt, hiba esetén nem nulla kilépési kód)
4. `compile_step26_clean_pdfs.py` (only with `--pdf`; pandoc + pdflatex; with `--pdf` a missing pandoc/pdflatex or a PDF that is not newly written gives a nonzero exit code; `--tex-only` writes only the `.tex` files with the separate status TEX-ONLY / `--pdf` esetén hiányzó pandoc/pdflatex vagy nem újonnan írt PDF nem nulla kilépési kódot ad; a `--tex-only` csak `.tex` fájlt ír, külön TEX-ONLY státusszal)

Without `--skip-block1`, `run_all.py` first runs `hubble_desi_bao_tal_audit.py` (Block 1: `run_step1_1.py`, `run_step1_2_and_1_3.py`, `run_step1_4_to_1_6.py`, then step 2 and the re-audit controls). / `--skip-block1` nélkül előbb a `hubble_desi_bao_tal_audit.py` fut (1. blokk, majd a 2. lépés és az újrabírálati kontrollok).

Input data / Bemeneti adatok:

- `desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`: 376 B, SHA-256 `dd2873a0b88459a491af3c0c0307ba059f62df9211d5b976760f310565a1be68`
- `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`: 2,170 B, SHA-256 `bbafa9074b51cf1a45e0d10e4f37db8c0e80a5d1d1788857abb7fc49fb21abcc`
- `tal_cosmology_common.py`: 24,847 B, SHA-256 `54c5636b35631af314f501060e39adf723a46cada78db15526b0591a87c2c109`
- `tal_hubble_desi_bao_summary.json`: 52,790 B, SHA-256 `94ba04e4c09fc8a7d6a94eed04eecf03c9f794444a3a9f646a58f80aa1e3ca53`

## Academia Co-scientist review (round 4, 2026-10-29)

The `academia_review/` folder contains the independent reviewer's closing report on the package with ZIP SHA-256 `e08c8dfcf5454c318762a9d154ad219ad2fc0b4907b465c6dafc8f3638c58150`, reproduced verbatim:

- `step26_academia_coscientist_review4_lezarasi_jelentes_hu.md`: closing report (HU). N01, N02 and T01 are closed, the F01–F12 closures remain valid, and no new required change is listed.
- `step26_academia_coscientist_review4_ellenorzesi_jegyzokonyv.json`: the reviewer's verification record, with input hashes, the run log, the numerical comparison, independent controls and the 16 mutation tests.

These two files are the reviewer's text and are not edited by the author.

License / Licenc: CC BY-NC 4.0 International
