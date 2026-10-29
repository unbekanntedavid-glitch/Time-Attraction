# Changelog

## Time Attraction Law — Scalar-Tensor Yukawa Fit, V27 — 2026-10-29

DOI: [10.5281/zenodo.23761574](https://doi.org/10.5281/zenodo.23761574)

| Version | DOI |
|---|---|
| V27 | [10.5281/zenodo.23761574](https://doi.org/10.5281/zenodo.23761574) |
| V26 | [10.5281/zenodo.23727585](https://doi.org/10.5281/zenodo.23727585) |
| V25 | [10.5281/zenodo.23699685](https://doi.org/10.5281/zenodo.23699685) |
| V24 | [10.5281/zenodo.23680715](https://doi.org/10.5281/zenodo.23680715) |

Zenodo V27: [10.5281/zenodo.23761574](https://doi.org/10.5281/zenodo.23761574). Earlier V26 reference: [10.5281/zenodo.23727585](https://doi.org/10.5281/zenodo.23727585) *(refers strictly to V26, `step1`–`step25`)*. Earlier V25 reference: [10.5281/zenodo.23699685](https://doi.org/10.5281/zenodo.23699685) *(refers strictly to V25)*. Earlier V24 reference: [10.5281/zenodo.23680715](https://doi.org/10.5281/zenodo.23680715) *(refers strictly to V24)*.

`step26` cosmology audit: Hubble tension vs. DESI DR1 BAO for the TAL scalar-tensor Yukawa model ($L = 10\;\text{km}$, $\alpha_Y = 1$). It uses the 12×12 DESI covariance with profiled $(\Omega_m, P)$ controls, model-specific re-integrated $r_s$/$r_d$, an illustrative Planck 3D prior, a synthetic SNe shape control, variational RSD, an algebraic BBN proxy, local Yukawa benchmarks, and no-go routes for a $\beta$ multiplier ($H_0$ 67.4 → 72.8). With $H_0 = 72.8$ held rigid, $\chi^2_{\text{BAO}} = 194.342$; the profiled fit gives $12.7405$. The required $r_d$ at $H_0 = 72.8$ is $140.02 \pm 1.74$ Mpc, and the joint $\chi^2$ is $15.6743 + 0.0623 + 3.9766 = 19.7133$. License: CC BY-NC 4.0.

Review fixes (four Academia Co-scientist rounds):
- **R1–R6:**
  - V26 normalization cross-check ($K_{V26}/K_{\text{new}} = 1$);
  - true piecewise late transition ($100\theta_*$ = 1.075660, 1.097867, 1.117963, 1.124630 for $z_t$ = 2.5, 10, 100, 1000);
  - $\chi^2$ split (R4);
  - symmetric SNe zero-point (M4 $\chi^2 = 0.9166$);
  - Fermi–Dirac neutrino check;
  - strict JSON with no NaN;
  - a self-contained `run_all.py`.
- **F01–F12:**
  - Jordan-frame coupling $G_T = G_* A^2$;
  - $G$ dictionary $G_{\text{lab}} \approx G_*(1+q)$, ratio 0.9262883;
  - M3 uses its own horizons ($r_s = 144.36622618$, $r_d = 147.01159289$ Mpc, $\chi^2_{\text{BAO}} = 2201.29927$);
  - 72.8±1.0 relabelled as an assumed local $H_0$ benchmark prior;
  - V21 profile values checked;
  - both $\theta_*$ variants reported (1.04248582, 1.04250287);
  - `--pdf` exits with code 2 when pdflatex is missing;
  - label, figure and export fixes.
- **N01:** the boundary cell is n.a. (JSON/CSV `null`) for the homogeneous FLRW, void and Milky Way disk rows, with a new `boundary_status` field.
- **N02:** the $\Omega_\phi \sim 10^{-47}$ background energy is separated from the direct RSD Yukawa term $\delta f\sigma_8 \approx 1.67\times10^{-40}$.
- **T01:** HU/EN text regression test `test_step26_text_regression.py` (77/77 checks, mutation-tested).

- Added [`V27/`](V27/): HU/EN `step26` Markdown chapters and PDFs (13 pages each), Python sources (`tal_cosmology_common.py`, `run_all.py`, `run_step1_*.py`, `run_step2_nogo_and_report.py`, `run_step26_revision_controls.py`), DESI DR1 inputs, CSV/PNG/JSON outputs, `academia_review/` (round-4 report and verification log, verbatim), and `SHA256SUMS.txt` (46 verified entries).

## Time Attraction Law — Scalar-Tensor Yukawa Fit, V26 — 2026-10-28

DOI: [10.5281/zenodo.23727585](https://doi.org/10.5281/zenodo.23727585)

| Version | DOI |
|---|---|
| V26 | [10.5281/zenodo.23727585](https://doi.org/10.5281/zenodo.23727585) |
| V25 | [10.5281/zenodo.23699685](https://doi.org/10.5281/zenodo.23699685) |
| V24 | [10.5281/zenodo.23680715](https://doi.org/10.5281/zenodo.23680715) |

Zenodo V26: [10.5281/zenodo.23727585](https://doi.org/10.5281/zenodo.23727585). Earlier V25 reference: [10.5281/zenodo.23699685](https://doi.org/10.5281/zenodo.23699685) *(refers strictly to V25)*. Earlier V24 reference: [10.5281/zenodo.23680715](https://doi.org/10.5281/zenodo.23680715) *(refers strictly to V24)*.

Extends the scalar-tensor, conformal/disformal, and Helmholtz–Yukawa analysis (`step1`–`step22`) with three bilingual HU/EN chapters (`step23`–`step25`), standalone Python solvers (`gw170817_and_lensing_audit.py`, `bullet_cluster_lensing_audit.py`, `cmb_and_cosmology_audit.py`), 9 CSV/JSON output tables, and 27-page HU/EN module PDFs (`step23_25_gw170817_lensing_bullet_cmb_hu.pdf`, `step23_25_gw170817_lensing_bullet_cmb_en.pdf`):
- **`step23` (GW170817 & Gravitational Lensing):** Horndeski tensor speed $G_4 = M_{\text{Pl}}^2/2 = \text{const.},\; G_5 = 0 \implies c_T = 1$ exactly ($m_g = 0$), scalar mode Compton wavelength $L = 10\;\text{km}$ ($m_\phi c^2 = 1.973\times 10^{-11}\;\text{eV}$), conformal photon null geodesics ($\gamma_{\text{PPN}} = 1$, $\Phi_{\text{lens}} = \Phi_N$) vs. disformal H1 clock coupling ($\tilde{g}_{\mu\nu}^{\text{clock}} = g_{\mu\nu} + c^{-2}(1 - e^{-2\phi_L})u_\mu u_\nu$), and 4-component exterior Yukawa BVP comparison between `Branch II-Fluid` ($\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$, source $-(\alpha_Y G/c^2)e^{-\phi_L}\varepsilon$, $M_{\text{ADM}} = 1.509998\,M_\odot$) and covariant trace-sourced `Branch II-Conf` ($S_m[e^{-2\alpha_Y\phi_L}g_{\mu\nu},\Psi_m]$, trace source $-(\alpha_Y^2 G/c^2)e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$, $M_{\text{ADM}} = 1.562035\,M_\odot$) on the $L = 10\;\text{km}$, $\alpha_Y = 1$, SLy $\log_{10}\rho_{0c} = 15.0$ reference star (`DOP853` `max_step = 0.02 km`, `tol_iter = 1e-10` mesh refinement $N = 601 \to 2000 \to 4000$: enforced Robin derivative agreement $\lvert\phi'_{L,\text{int,Robin}}(R_*) - \phi'_{L,\text{ext,curved}}(R_*)\rvert = 3.937\times 10^{-15}\;\text{km}^{-1}$, independent one-sided FD2 derivative difference $\lvert\phi'_{L,\text{int,FD2}}(R_*) - \phi'_{L,\text{ext,curved}}(R_*)\rvert = 4.7228\times 10^{-10}\;\text{km}^{-1}$ at $N = 4000$ converging as $\mathcal{O}(h^2)$ with empirical ratio $4.0332$).
- **`step24` (Bullet Cluster 1E 0657-558):** 1D convergence ($\kappa$) and peak-offset analysis across Model A ($L = 10\;\text{km}$, gas + stars: $\Delta x_{\text{peak}} = 0.00\;\text{kpc}$), Model B ($L = 62\;\text{kpc}$, unscreened galactic Yukawa kernel on gas + stars: $\Delta x_{\text{peak}} = 0.00\;\text{kpc}$, peak stays on X-ray gas at $+150\;\text{kpc}$), and Model C ($L = 10\;\text{km}$ + collisionless component: $\Delta x_{\text{peak}} = 169.39\;\text{kpc}$ vs. observed $170\;\text{kpc}$).
- **`step25` (FLRW Background & CMB Acoustic Scales):** Scale-dependent effective gravitational coupling $\Delta G_{\text{eff}}(k, z_*)/G_N = 2.55\times 10^{-36}$ at $\ell_1 = 220$ and $4.65\times 10^{-34}$ at $\ell = 3000$ ($z_* = 1089.80$, $L = 10\;\text{km}$, $\Omega_{\phi,\text{can}}(z_*) = 7.45\times 10^{-42}$), and sound-horizon shift $\lvert\Delta r_s/r_s\rvert < 10^{-15}$ without a collisionless cold component vs. $+31.94\%$ in a baryon-only cosmology. License: CC BY-NC 4.0.

- Added [`V26/`](V26/): bilingual HU/EN chapters `step1`–`step25`, `README.md`, `README_en.md`, `MANIFEST.txt`, `MANIFEST_en.txt`, `step23_25_gw170817_lensing_bullet_cmb_hu.pdf` (27 pages), `step23_25_gw170817_lensing_bullet_cmb_en.pdf` (27 pages), `V26_english_audit_package.pdf` / `V25_english_audit_package.pdf` (167 pages), `V26_english_translations.zip` / `V25_english_translations.zip` (27 files), `scalar-tensor-yukawa-fit_audit_package.zip` (393 files), `v26_hianypotlas_szamitas.md`, `v26_hianypotlas_ellenorzese_2026-10-28.md`, `SHA256SUMS_step23_25.txt` (20 verified entries), and `SHA256SUMS.txt` (392 verified entries).

## Time Attraction Law — Scalar-Tensor Yukawa Fit, V25 — 2026-10-27

DOI: [10.5281/zenodo.23699685](https://doi.org/10.5281/zenodo.23699685)

| Version | DOI |
|---|---|
| V25 | [10.5281/zenodo.23699685](https://doi.org/10.5281/zenodo.23699685) |
| V24 | [10.5281/zenodo.23680715](https://doi.org/10.5281/zenodo.23680715) |

Earlier V24 base reference: [10.5281/zenodo.23680715](https://doi.org/10.5281/zenodo.23680715) *(refers strictly to V24; DOI [10.5281/zenodo.23699685](https://doi.org/10.5281/zenodo.23699685) refers to V25. The V25 sign fixes, diagnostics, and conditional stability closure belong to V25)*.

Scalar-tensor, conformal/disformal, and Helmholtz–Yukawa (`yukawa_action_el`, Branch II-Fluid) analysis across galactic (175 SPARC galaxies, NGC 3198 control), laboratory (8558 kg granite chamber, 8-coil solenoid), white-dwarf (Sirius B conditional H0/H1), and relativistic neutron-star regimes (`step1`–`step22` in Hungarian and English). V25 incorporates the 4-component curved exterior Einstein–scalar vacuum BVP ($y=[\phi_L,\phi_L',m,\nu]$, $\nu_R(15.0)=-0.2463066872$), the exterior metric-coupling sign fix (`a_m = (2*psi0/r^2 + 2*dpsi0/r)/B0`) and the exact left-hand side sign `scalar_m_sign = -1.0` in the coupled GL2 radial pulsation operator ($L_u u = R_h\,\delta h$), multi-EOS central finite-difference Jacobian verification (`full_jacobian_multi_eos_and_points`), $s=r^2$ PCHIP + adaptive Gauss–Kronrod quadrature baryonic rest mass ($M_{B,\max}^{\text{quad}} = 2.42339\,M_\odot$ for SLy $L=10\;\text{km}$, $\alpha_Y=1$) vs. trapezoidal and single-pass ODE baryonic mass turning-point diagnostics (`baryonic_mass_and_first_law_audit`), and eigenfunction vs. static-tangent comparisons (`eigenmode_xi_nodes_and_static_tangent_audit`: dynamical $\omega_0^2=0$ root at $15.316410$ vs. ADM peak $15.316404$, $+5.80\times 10^{-6}\;\text{dex}$; interior relative $L^2$ differences $\delta h=5.04\times 10^{-4}$, $\delta m=1.77\times 10^{-4}$, $\delta\phi_L=2.27\times 10^{-4}$, $\delta\phi_L'=3.09\times 10^{-4}$ and $3.54\%$ free-surface displacement difference, with the stability closing statement formulated conditionally pending independent verification of the coupled variational derivation and boundary-condition matching). License: CC BY-NC 4.0.

- Added [`V25/`](V25/): bilingual HU/EN chapters `step1`–`step22`, `README.md`, `README_en.md`, `MANIFEST.txt`, `MANIFEST_en.txt`, `V25_english_audit_package.pdf` (109 pages), `V25_english_translations.zip`, `scalar-tensor-yukawa-fit_audit_package.zip` (341 files), `SHA256SUMS.txt` (341 verified entries), Python solvers, and CSV/JSON outputs.


## Time Attraction Law — Scalar-Tensor Yukawa Fit, V24 — 2026-10-27

Zenodo V24: [10.5281/zenodo.23680715](https://doi.org/10.5281/zenodo.23680715)

This is the V24 version of the Time Attraction Law (TAL) research, an independent theoretical physics program developing a scalar-tensor/dilaton framework for gravity. It documents derivation steps 1–22: dilaton definitions from V19–V23, Einstein↔Jordan/Clock conformal and disformal transformations, the Helmholtz–Yukawa action and the massive dilaton Klein–Gordon equation, covariant conservation and modified hydrostatic equilibrium, the TAL↔scalar-tensor/Brans–Dicke coupling dictionary with PPN parameters, and SPARC galactic consistency analysis (175 galaxies: 0 negative-density points with the 3D top-hat/Yukawa approach vs 239/303 with the derivative formula; NGC 3198 control: V_obs=146 km/s, V_bar=88.53 km/s, C_req=2.720). All Python programs, input and output data files, CSV/JSON tables, and SHA-256 checksums are included for full reproducibility. License: CC BY-NC 4.0.

- Added [`V24/`](V24/): Markdown chapters `step1`–`step22`, Python programs, input and output data files, CSV/JSON tables, `README.md`, `MANIFEST.txt`, `SHA256SUMS.txt` (313 verified entries), and `scalar-tensor-yukawa-fit_audit_package.zip`.
- English translation files: `step1_*_en.md`–`step22_*_en.md`, `README_en.md`, `MANIFEST_en.txt`.

## V23 (Zenodo V23) — 2026-10-25

Zenodo V23: [10.5281/zenodo.23639084](https://doi.org/10.5281/zenodo.23639084)  
Previous version: V22 Book Edition (Zenodo V22, [10.5281/zenodo.23604423](https://doi.org/10.5281/zenodo.23604423)).

- Added [`V23/`](V23/): Task `50aa2db0` bilingual (EN/HU) 4-page technical reports (`tal_covariant_closure_50aa2db0_EN_2026-10-24.pdf`, `tal_covariant_closure_50aa2db0_HU_2026-10-24.pdf` and `.tex` sources), standalone archive (`tal_covariant_closure_50aa2db0_2026-10-24.zip`, SHA-256 `807082095cffb35433fbd720f5e04c0e4770723c4111fc3f3ce9289c1f7f4508`), Python solvers (`covariant_closure_50aa2db0.py` and bundled V21 modules), 4 CSV/JSON data tables, 3 figures, and `SHA256SUMS.txt` (17 verified entries).
- Updated `README.md`, `CITATION.cff`, and `CHANGELOG.md` with DOI `10.5281/zenodo.23639084`.

## V22 Book Edition (Zenodo V22) — 2026-10-24

Zenodo V22 Book Edition: [10.5281/zenodo.23604423](https://doi.org/10.5281/zenodo.23604423)  
Previous version: Zenodo V21 (= GitHub V19/V21, [10.5281/zenodo.23602501](https://doi.org/10.5281/zenodo.23602501)).

- Added [`V22_Book/`](V22_Book/): 15-page, 12-chapter (+ 2-page Round-3 technical reproduction audit appendix) research monograph synthesizing Versions 1–21 of the Time Attraction Law in English (`tal_book_EN_2026-10-24.pdf` & `.tex`, primary, byte-identical to Zenodo record `23604423`), Hungarian (`tal_book_HU_2026-10-24.pdf` & `.tex`, byte-identical to Zenodo record `23604423`), and German (`tal_book_DE_2026-10-24.pdf` & `.tex`), together with `generate_figures.py`, 12 figures in `figures/`, and `SHA256SUMS.txt` (21 verified entries).
- Updated `README.md` and `CITATION.cff` with DOI `10.5281/zenodo.23604423`.

## V18 (GitHub) — 2026-10-21

GitHub V18 = Zenodo V20 (10.5281/zenodo.23540741)

Előző verzió: V19 – 7th model (10.5281/zenodo.23497128).
A Zenodo V20 rekord négy fájlja változatlanul a `V18/zenodo_v20/` mappában van. A `V18/final/` a V20 jelentés-PDF-hez
képest a revision3 review két szövegjavítását (és a review szerinti két pontosítást) tartalmazza; a kód és a CSV/JSON fájlok
bájtra azonosak a V20 csomagéval.

### Változások a V19 óta (a Zenodo V20 jelentésben és a `V18/final`-ban egyaránt)
- **Integrálfolytonosság:** L ≥ 2.6 m-re közvetlen külső különbségi integrál r = L/t helyettesítéssel (32 pontos
  Gauss–Legendre), levágás és radiális interpoláció nélkül; a korábbi 0.01268 J/A²-es mesterséges eltolás megszűnt.
  Folytonossági teszt L = 2.99 / 3.00 / 3.01 m-nél: relatív eltérés 0 / −7.4·10⁻¹⁶ / −1.5·10⁻¹⁵; a 2.6–2.9 m-es
  átfedésben ≤ 2.9·10⁻¹⁵ (`csv/continuity_test.csv`).
- **Háttértér-kereszttag lehatárolása:** a §4 polaritásmondata kiegészítve a 2B_bg·B_coil kereszttaggal; a közölt
  −1.20×10⁻⁴⁸ érték az izolált páros tagra vonatkozik, nem minden lehetséges egyszerű ON/OFF mérés teljes előrejelzése.
- **Biztonsági mondat cseréje:** az Onderdonk-becslés csak idealizált adiabatikus olvadási idő; ebből biztonságos
  impulzushossz nem állapítható meg, a számított nagyáramú esetek nem üzemeltetési javaslatok.
- **291.7 mm mért menetátmérő, mérési jegyzőkönyvvel (1.2 pont):** 2025-11-13, tolómérő (0.02 mm felbontás), mérte
  Szabolcsi Dávid, 5 egymást követő leolvasás, mind 291.7 mm (átlag 291.7 mm, mintaszórás 0); a 300 mm a tekercstest
  pereme (külső méret). A térmodell a mért 291.7 mm-rel számol (B a fő órán 11.63 µT).
- **Konvergencia-közlés:** a szórás definíciója (a 4 futás mintaszórása, n − 1 nevezővel), teljes tartomány és
  mintaszórás táblázatban; a négy MC relatív hiba külön: 2.17 % / 0.99 % / 1.14 % / 4.04 % (L = 0.316 / 0.631 / 1.0 / 2.0 m).
- **„Főóra-egyenérték” megnevezés** a táblázatban (2.20·10⁻⁶ /T; 0.189 /T²).

### A revision3 review pontjai, a `V18/final`-ban alkalmazva (a Zenodo V20 PDF-ben még nem)
- **Távoli „0 szórás” megjegyzés** (konvergenciatáblázat, L = 2.6–100 m): „A négy futás ebben a tartományban azonos nlo
  és nt beállítást használ. A nulla mintaszórás ezért nem távoli hibabecslés; a távoli konvergenciát nlo és nt külön
  változtatásával kell vizsgálni.”
- **Hőkontroll átfogalmazása** (§4, 5. pont): a „Így a be/ki különbség nem tartalmaz hőmérsékleti lépést.” mondat helyett:
  „Az azonos összteljesítményű, kis maradékterű kontroll célja a kapcsolással korrelált termikus változás csökkentése.
  A be- és kikapcsolt állapotok hőmérsékleti azonosságát és a kontroll maradékterét méréssel ellenőrizni kell; a maradék
  eltérést a frekvenciabecslés bizonytalanságába kell foglalni.”
- **L = 0.95 m illesztési hiba elismerése:** a V20 PDF már közölte a teljes energiák 1.4·10⁻³ J/A²-es (0.06 %) eltérését
  és az ebből adódó ugrást; a `V18/final` a pontos értékeket is megadja: belső forma 2.17738947, külső forma
  2.17597607 J/A², illesztési hiba 0.00141340 J/A² (a belső érték 0.0649 %-a).
- **§5 első pontja:** a −1.202·10⁻⁴⁸ az „izolált, polaritásra páros tekercs-B² tag” számított differenciális TAL-jele.

### V20 kézirat (angol, Foundations of Physics)

- Hozzáadva: `V18/manuscript/tal_v20_fop.pdf` és `V18/manuscript/tal_v20_fop_manuscript.zip` (LaTeX-forrás, ábrák, szkript).
- A kézirat tartalmazza a revision3 review két szövegcseréjét (a távoli nulla mintaszórásról szóló megjegyzést és a hőkontroll-mondatot), és hivatkozik a `revision3_checks.json`-ra.
- A Zenodo-feltöltést a szerző maga végzi; ez a commit a Zenodóhoz nem nyúl.

### Ellenőrzések
- `V18/revision3_checks.json`: a szerző pipeline-ja generálta 2026-10-21-én (`V18/make_revision3_checks.py`), nem a bíráló
  eredeti fájlja (a bíráló `check_revision3.py` szkriptje nem áll rendelkezésre). Tartalma: 18/18 SHA-256 egyezés a
  `V18/final/SHA256SUMS.txt` alapján; B a fő és a referenciaórán független köráram-összeggel (8 × 2000 menet) vs. a csomag
  térfüggvénye (relatív eltérés −4.3·10⁻⁸, ill. −7.6·10⁻⁹; 11.6329 µT, ill. 0.43054 µT); célzott főcsúcs-újrafuttatás
  L = 0.630957 m-nél (−1.2017·10⁻⁴⁸, egyezik a JSON-nal); külső integrálok L = 2.6 / 3.0 / 10 / 100 m-nél, nlo és nt
  együttes finomításával (6000 → 12000, 32 → 64; relatív változás −9.1·10⁻⁷ … −3.0·10⁻⁴); az L = 0.95 m-es maradék
  számtana (0.00141340 J/A², 0.0649 %, egyezik a teljes energiák különbségével).

### Repó-szintű változások
- `LICENSE`: a Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) hivatalos jogi szövege.
- `CITATION.cff`: érvényes CFF formátum, licenc `CC-BY-NC-4.0`, DOI 10.5281/zenodo.23540741 (Zenodo V20 / GitHub V18),
  javított GitHub URL (`unbekanntedavid-glitch`).
- Nem része a csomagnak: NIST SP 1065 PDF.
