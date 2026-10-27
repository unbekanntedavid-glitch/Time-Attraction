# Skalár–tenzor és Yukawa/diszform hatáselv-illesztés (V25; V23 $\leftrightarrow$ Jordan–Brans–Dicke / Dilaton keretrendszer) — 1–22. lépéses levezetési és reprodukciós csomag
# Scalar-Tensor and Yukawa/Disformal Action Mapping (V23 $\leftrightarrow$ Jordan–Brans–Dicke / Dilaton Framework) — 22-Step Derivation and Audit Reproduction Package

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-27  
**Munkakönyvtár / Working Directory:** `scalar-tensor-yukawa-fit`  
**Önálló bírálati archívum / Self-Contained Audit Archive:** `scalar-tensor-yukawa-fit_audit_package.zip`  
**Licenc / License:** CC BY-NC 4.0 — Kereskedelmi hasznosítás kizárólag a szerző előzetes írásbeli engedélyével.  
**Publikációs státusz és verzióelhatárolás (V24 vs. V25):** A korábbi Zenodo-publikáció (**V24**, DOI: `10.5281/zenodo.23680715`) kizárólag a korábbi V24-es kiadásra mutat. A `curved_yukawa_adm_and_stability.py` előjeljavításai (a `linearized_exterior_response` 593. sorában az `a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0` külső metrikus csatolási tag, valamint a `solve_yukawa_coupled_omega2_gl2` bal oldali $L_u u = R_h\,\delta h$ operátorában az egzakt `scalar_m_sign = -1.0`, azaz `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n` alapértelmezetté tétele), a `curved_yukawa_audit_summary.json` három új diagnosztikája (`full_jacobian_multi_eos_and_points`, `baryonic_mass_and_first_law_audit`, `eigenmode_xi_nodes_and_static_tangent_audit`), továbbá az elméleti ellenőrzés lezárása kizárólag az új **V25-ös** verzióhoz tartozik — a korábbi V24-es DOI-ra soha nem szabad úgy hivatkozni, mintha e javítások és lezárások már abban szerepeltek volna.

**V25 diagnosztikai összefoglaló (a `curved_yukawa_audit_summary.json` és `curved_yukawa_coupled_stability.csv` pontos adatai):**
- Deriváltellenőrzés (SLy, $\log_{10}\rho_{0c}=15{,}0$, $L=10$ km): `exterior_dphi_prime_dm_R` $=-1{,}683651\times10^{-4}\;\text{km}^{-2}$ (rel. eltérés $3{,}98\times10^{-10}$); `interior_r6km_a_m_sc` ($\partial\phi_L''/\partial m$) $=-1{,}728334\times10^{-4}\;\text{km}^{-3}$ (rel. eltérés $6{,}08\times10^{-9}$).
- Teljes 4×4 Jacobi: 4 konfiguráció (SLy és APR4), 108 ellenőrzés, legnagyobb relatív hiba $2{,}74\times10^{-5}$.
- Barionikus tömeg / első főtétel: az ADM és az adaptív Gauss–Kronrod barionikus fordulópont ($M_B^{\text{quad}}$) közötti `delta_dex_ADM_vs_baryon_turning` $-4{,}58\times10^{-6}$ … $+8{,}43\times10^{-6}$ dex; nem-konform többlet $8{,}23\times10^{-5}\,M_\odot$/dex.
- Fő stabilitási zérushely és sajátmódus-audit (`scalar_m_sign = -1.0` egzakt LHS előjellel): a `curved_yukawa_coupled_stability.csv` fő dinamikus $\omega_0^2=0$ zérushelye SLy $L=10\;\text{km},\;\alpha_Y=1$ esetén **`15.316410`** ($+5{,}80\times10^{-6}$ dex az ADM-csúcs `15.316404` és $+1{,}33\times10^{-6}$ dex az $M_B^{\text{quad}}$-csúcs `15.316409` felett; mind az 5 sorozaton $\le 6{,}9\times 10^{-6}$ dex eltérés, csomópontmentes $n=0$ alapmódus és $(0{,}6\text{–}5{,}0)\times 10^{-4}$ relatív $L^2$-egyezés a statikus érintővel; a pluszjeles `scalar_m_sign = +1.0` diagnosztikai gyök `15.319129`, $+2{,}725\times 10^{-3}$ dex — részletek: `step22`, 4. és 7. szakasz).

**Kapcsolódó hivatalos archív DOI-k / Official Archival DOIs:**  
- **V14:** `10.5281/zenodo.23217189`  
- **V17 (SPARC 155 / 123 vakminta):** `10.5281/zenodo.19666220`  
- **V19 (3D gömbi top-hat SPARC, 175 galaxis):** `10.5281/zenodo.23497128`  
- **V20 (Kétórás gránitkamra és 8 tekercses mágneses teszt):** `10.5281/zenodo.23540741`  
- **V21 (Kompakt csillagok, TOV és Szíriusz B audit):** `10.5281/zenodo.23602501`  
- **V22 Book Edition (Háromnyelvű monográfia, HU/EN/DE):** `10.5281/zenodo.23604423`  
- **V23 (Kovariáns hatáslezárás és csatolt radiális stabilitás, `50aa2db0`):** `10.5281/zenodo.23639084`  

---

## 1. Kutatási kérdés, elméleti keretrendszer és modellág-taxonómia (Research Question, Theoretical Framework & Branch Taxonomy)

Ez a bírálati csomag a **Time Attraction Law (Idővonzástörvény, TAL)** V19–V23-as kézirataiban és monográfia-kiadásában (`10.5281/zenodo.23604423`, `10.5281/zenodo.23639084`) kidolgozott kovariáns hatásjelölteknek — az **I. megfogalmazásnak** (nemlokális, pozitív definit bilaterális top-hat mag-hatás segédmezőkkel) és a **II. megfogalmazásnak** (lokális, árnyékolt Helmholtz–Yukawa skalármezős hatáselv) — az általános **skalár–tenzor gravitációs elméletekkel** (Jordan–Brans–Dicke-család, alacsony energiás húrelméleti **dilaton**-csatolás és Bekenstein-típusú **diszform** órametrikák) való egzakt, 22 lépéses matematikai megfeleltetését és számszerű ellenőrzését tartalmazza.

### 1.1. Alapegyenletek és kanonikus normálás (Core Equations & Canonical Normalization)
1. **Exponenciális időfaktor és skalár-/dilatonmező (`step1`):**
   $$f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(\mathbf{x})\bigr) \;\equiv\; e^{-\phi_L(\mathbf{x})}, \qquad \phi_L(\mathbf{x}) \;=\; \gamma_L\,\bar{\rho}_L(\mathbf{x}) \;\ge\; 0,$$
   ahol a kanonikus csatolási együttható a gravitációs térfogati kapacitásból ($C = \frac{4\pi G}{c^2} = 9{,}33306\times 10^{-27}\;\text{m/kg}$) származik:
   $$\gamma_L \;=\; \frac{G L^2}{c^2} \;=\; \frac{C}{4\pi}\,L^2, \qquad m_\phi \;=\; \frac{\hbar}{L c}.$$
2. **Konform Jordan-metrika vs. közeghez kötött diszform órametrika (`step2`, `step5`, `step8`, `step14`):**
   $$\tilde{g}_{\mu\nu}^{\text{conf}} \;=\; A^2(\phi_L)\,g_{\mu\nu} \;=\; f_L^2\,g_{\mu\nu} \qquad \text{vs.} \qquad \tilde{g}_{\mu\nu}^{\text{disf}} \;=\; g_{\mu\nu} + \frac{1 - f_L^2}{c^2}\,u_\mu u_\nu \;\Longrightarrow\; \tilde{g}_{00}^{\text{disf}} \;=\; f_L^2\,g_{00}\;(\text{együttmozgó keretben}).$$
   Együttmozgó ($u^\mu = (c/\sqrt{-g_{00}},\mathbf{0})$, $u_\mu = (-\sqrt{-g_{00}}\,c,\mathbf{0}) \implies u_0 u_0/c^2 = -g_{00}$) barionos közegben és atomórákban mindkét metrika ugyanazt a sajátidő-faktort ($d\tilde{\tau} = f_L\,d\tau_g$) adja. Ultrarelativisztikus fénysugarakra ($g_{\mu\nu}k^\mu k^\nu = 0$) és PPN-tesztekre azonban a két ág élesen elválik:
   - Az univerzális **konform Jordan-ágban (`Branch II-Conf`)** $\sqrt{-\tilde{g}^{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$, a Maxwell-hatás konform-invariáns, és érvényes a standard Brans–Dicke / PPN szótár:
     $$\alpha_{\text{std}} \;=\; 2\beta^2 \;=\; 2\tilde{\alpha}_L \;=\; \frac{\alpha_Y}{4\pi}, \qquad \omega_{\text{BD}} \;=\; \frac{1 - 3\alpha_{\text{std}}}{2\alpha_{\text{std}}} \;=\; \frac{4\pi}{\alpha_Y} - \frac{3}{2}, \qquad \gamma_{\text{PPN}}^{\text{conf}}(r) \;=\; \frac{1 - \alpha_{\text{std}}(1+r/L)e^{-r/L}}{1 + \alpha_{\text{std}}(1+r/L)e^{-r/L}}.$$
   - A közeghez kötött **diszform ágban (`Branch II-Disf`)** $\sqrt{-\tilde{g}^{\text{disf}}} = f_L\sqrt{-g}$, a mágneses energia ($\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$) közvetlenül forrásozza $\phi_L$-t (`step15`), a fényt pedig $n_t(s) = f_L(s)^{-1} = e^{+\phi_L(s)}$ akromatikus törésmutató téríti el (`step9`, `step14`), így erre az ágra a konform Brans–Dicke $\gamma_{\text{PPN}}^{\text{conf}}$ képlet nem alkalmazható.
3. **V23 lokális Helmholtz–Yukawa hatáselv, egzakt normálás és görbült külső peremérték-feladat (`step3`, `step5`, `step7`, `step22`):**
   A skalár kinetikus és tömegtag SI-egységekben vett pontos együtthatója ($\alpha_Y \equiv 8\pi\tilde{\alpha}_L = 4\pi\alpha_{\text{std}}$ definícióval):
   $$\boxed{\;S_{\text{II}}[g_{\mu\nu},\phi_L,\Psi] \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}\left(R \;-\; \frac{1}{\tilde{\alpha}_L}\Bigl(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{\phi_L^2}{L^2}\Bigr)\right) \;+\; e^{-\phi_L}\,\mathcal{L}_{\text{matter}}\right] \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}R \;-\; \frac{c^4}{2 G\,\alpha_Y}\Bigl(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{\phi_L^2}{L^2}\Bigr) \;+\; e^{-\phi_L}\,\mathcal{L}_{\text{matter}}\right],\;}$$
   amely geometrizált ($G = c = 1$) egységekben a V23 kód $\mathcal{L}_\phi = -\frac{1}{2\alpha_Y}\bigl((\nabla\phi_L)^2 + \phi_L^2/L^2\bigr)$ alakját és az egzakt görbült téridőbeli Euler–Lagrange-egyenletet (`yukawa_action_el`, $B(r) = 1 - 2m(r)/r$) adja:
   $$B\,\phi_L'' \;+\; \left[\frac{2B}{r} + B\nu' + \frac{B'}{2}\right]\phi_L' \;-\; \frac{\phi_L}{L^2} \;=\; -\frac{\alpha_Y G}{c^2}\,e^{-\phi_L}\,\varepsilon(r), \qquad \rho_\phi(r) \;=\; \frac{c^2}{2 G\,\alpha_Y}\left[B\bigl(\phi_L'\bigr)^2 + \frac{\phi_L^2}{L^2}\right].$$
   A csillag felszínén kívüli görbült vákuumban ($r \ge R_*$, ahol $\varepsilon = P = 0$, $B\nu' + B'/2 = \frac{2m}{r^2} - \frac{4\pi r}{\alpha_Y L^2}\phi_L^2$) a csatolt külső Einstein–skalár egyenletrendszer (`step22`, `curved_yukawa_adm_and_stability.py`):
   $$\left(1 - \frac{2m(r)}{r}\right)\phi_L'' \;+\; \left[\frac{2}{r}\left(1 - \frac{m(r)}{r}\right) - \frac{4\pi r}{\alpha_Y L^2}\phi_L^2\right]\phi_L' \;-\; \frac{\phi_L}{L^2} \;=\; 0, \qquad \frac{dm(r)}{dr} \;=\; 4\pi r^2\rho_\phi(r),$$
   amely egyértelműen meghatározza a valódi görbült Dirichlet–Neumann-meredekséget $\phi_{L,\text{curved}}'(R_*) = -k_{\text{curved}}\phi_L(R_*)$ (SLy $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $\alpha_Y=1$ esetén $k_{\text{curved}} = 0{,}233542\;\text{km}^{-1}$, ami $+24{,}74\%$-kal meredekebb a lapos téri $k_{\text{flat}} = 1/R_* + 1/L = 0{,}187216\;\text{km}^{-1}$ Robin-közelítésnél) és az aszimptotikus **ADM-tömeget**:
   $$M_{\text{ADM}} \;\equiv\; \lim_{r\to\infty} m(r) \;=\; M_{\text{surf}} \;+\; \Delta M_{\phi,\text{ext}}^{\text{curved}}, \qquad M_{\text{surf}} \;\equiv\; m(R_*).$$

### 1.2. A modellágak egzakt taxonómiája (Explicit Branch Taxonomy Table)

A csomagban vizsgált öt elméleti ág feltevéseit, momentumait és érvényességi körét az alábbi táblázat különíti el (`step2`, `step5`, `step10`, `step11`, `step21`, `step22`):

| Modellág / Branch | Hatáselv / metrikus definíció | Mag-momentumok ($\bar{\rho}_L(\mathbf{x})$ sorfejtése) | Belső hidrosztatikai egyensúly | Külső ($r > R_*$) vákuum és tömeg | Érvényességi kör és lépések |
|---|---|---|---|---|---|
| **Branch I-A** (V21 nemlokális 3D top-hat: `tophat_mass`, `tophat_gradient`) | Fenomenologikus 3D gömbi top-hat átlag + algebrai Bianchi X-szektor | $\left(1 + \frac{L^2}{10}\nabla^2 + \mathcal{O}(L^4)\right)\rho$ | `mass`: $h' = -\nu'$; `gradient`: $h' = -\nu' = -f_L A_{\text{GR}}$ | Kompakt tartó $[R_*, R_*+L]$; $m(r) = M_{\text{surf}}$ konstans $r \ge R_*$-ra | V21 769 modell, V23 Eq. (8)–(9) statikus érintő BVP (`step7`, `step16`–`step18`, `step21`) |
| **Branch I-B** (V23 I. megfogalmazás: `bilateral_lapse_weighted`) | Segédmezős ($\phi_L,\lambda_L$) bilaterális mag-hatás: $S_m = -\int\sqrt{-g}\,e^{-\phi_L}\varepsilon$ | $\left(1 + \frac{L^2}{10}\nabla^2 + \mathcal{O}(L^4)\right)\rho$ | $h' = -\nu' - (\ln f_{\text{recip}})'$ (reciprok lapse-súlyozott erő) | Kompakt tartó $[R_*, R_*+L]$ | V23 reciprok kétoldali potenciál $V_{\text{eff},N}(r)$ és Synge-határvariáció (`step10`, `step13`) |
| **Branch II-Fluid** (V23 II. megfogalmazás: `yukawa_action_el` / `yukawa_action_el_curved_adm`) | Lokális Helmholtz–Yukawa hatás folyadéksűrűség-csatolással: $\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$ | $(1 - L^2\nabla^2)^{-1}\rho = \left(1 + L^2\nabla^2 + \mathcal{O}(L^4)\right)\rho$ | $h'(r) = -\nu'(r) + \phi_L'(r)$ (`sym_yukawa_el_hydrostatic = 0`) | Görbült Einstein–skalár farok $r \in [R_*,\infty)$; $M_{\text{ADM}} = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$ | Görbült külső ADM-tömeg, tömeg–sugár görbék és csatolt $(\delta h,\delta\phi_L,\delta m)$ stabilitás (`step3`, `step4`, `step7`, `step12`, `step22`) |
| **Branch II-Conf** (Univerzális konform Jordan-ág) | $\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$, $\sqrt{-\tilde{g}^{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$, $S_m[\tilde{g}^{\text{conf}},\Psi]$ | $(1 - L^2\nabla^2)^{-1}\rho_{\text{trace}}$ | $h' = -\nu' + \phi_L'$ (nyommentes EM nem forrásozza $\phi_L$-t) | Yukawa skalárfarok; konform-invariáns null-geodetikusok | Standard Brans–Dicke / PPN szótár ($\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$, $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2}$, $\gamma_{\text{PPN}}^{\text{conf}}(r)$; `step1`, `step2`, `step5`, `step14`) |
| **Branch II-Disf** (Közeghez kötött tisztán időbeli diszform ág) | $\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$, $\sqrt{-\tilde{g}^{\text{disf}}} = f_L\sqrt{-g}$ | Top-hat vagy Yukawa mag az együttmozgó energiasűrűségre | Együttmozgó órákra $d\tilde{\tau} = f_L d\tau_g$; $u_\mu$-variáció extra kényszert ad | Mágneses energia ($\rho_B = B^2/(2\mu_0 c^2)$) forrásozza $\phi_L$-t; $n_t = e^{+\phi_L}$ törésmutató | V20 8 tekercses mágneses teszt, V20 gránitkamra H1, Szíriusz B H1, VHE gamma-diszperzió (`step2`, `step8`, `step9`, `step14`, `step15`, `step19`, `step20`) |
| **Branch II-Cham** (Nemlineáris kaméleon/szimmetron kiterjesztés) | Hipotetikus sűrűségfüggő $\beta_{\text{eff}}(\rho)$ és önhatási potenciál $V_0(\phi_L)$ | Sűrűségfüggő $L_{\text{eff}}(\rho) = [\partial^2 V_{\text{eff}}/\partial\phi_L^2]^{-1/2}$ | Árnyékolt (vékony héj) nagy sűrűségnél, felerősített kis sűrűségnél | Környezetfüggő Yukawa-hatótáv | Külön elméleti kiterjesztési hipotézis a konstans $\alpha_Y$-ú V23 hatáson túl (`step6`, `step11`) |

*(Megjegyzés a top-hat $\leftrightarrow$ Yukawa megfeleltetéshez: a 3D gömbi top-hat mag másodrendű Taylor-momentuma $\frac{1}{10}L^2\nabla^2$, míg a Yukawa Green-függvényé $L^2\nabla^2$, ezért az $\alpha_Y = 1$ választás a nulladrendű homogén [$\nabla^2\rho \to 0$] határesetet illeszti; a gradiens-érzékeny másodrendű momentum illesztéséhez $L_{\text{Yukawa}} = L_{\text{top-hat}}/\sqrt{10}$ vagy $\alpha_Y$ átskálázása szükséges — lásd `step10`.)*

---

## 2. A 22 levezetési lépés áttekintése (Step-by-Step Summary of Steps 1–22)

Az alábbi táblázat felsorolja a `scalar-tensor-yukawa-fit` mappában található **22 levezetési dokumentumot** (`step1`–`step22`), azok pontos fájlnevét, méretét, SHA-256 ellenőrzőösszegét és egymondatos tartalmi összefoglalóját:

| Lépés / Step | Fájlnév / File Name | Méret / Size (B) | SHA-256 (első 16 hex / prefix) | Egymondatos tartalmi összefoglaló (HU / EN One-Line Gist) |
|---:|---|---:|---|---|
| **1** | `step1_dilaton_form.md` | 11 20511 282ecd894436341cbdc6d12c0fd49bb3a50
| **2** | `step2_conformal_disformal.md` | 18 40218 4257a3761a4adcca6fd5a9fe423508e3329
| **3** | `step3_yukawa_helmholtz_eom.md` | 12 64712 67077fa5ad7646bad0b9ca162a9375eb71c
| **4** | `step4_covariant_conservation_hydrostatic.md` | 12 35612 37952b1cab653da99d68f5446c737f3c39f
| **5** | `step5_coupling_dictionary_ppn_screening.md` | 14 94914 972d9c9d69ff93e662533ce437c1eec32ba
| **6** | `step6_sparc_galactic_ngc3198_consistency.md` | 14 61514 6928fbdd2157d0a4aa91d3c17156dbf5ed4
| **7** | `step7_neutron_star_tov_stability.md` | 21 21520 579d656ac9249a495ae97e00c320ea3efa9
| **8** | `step8_microscopic_clock_coupling_siriusb_lab_v20.md` | 16 30616 329f4d790f47b5ab8d48e3645451604bbc8
| **9** | `step9_vhe_photon_propagation_hess_dispersion.md` | 12 16412 187af5b8d0ab66ff02c9e735e9bfd15e9bc
| **10** | `step10_bilateral_kernel_reciprocity_eft_synthesis.md` | 13 42013 4433cfe4aaa90e66f3f4e336d901146cfea
| **11** | `step11_chameleon_screening_multiscale_eft_synthesis.md` | 13 29713 32069c908b620c756263d7262216b373998
| **12** | `step12_second_variation_coupled_pulsation_hessian.md` | 19 54519 43466e0fdb1964263f1c98a4463a55c1f07
| **13** | `step13_curved_window_metric_variation_synge_world_function.md` | 11 59511 618261d16eadbf89cf4ec8dcd74bd038ec2
| **14** | `step14_gravitational_lensing_null_geodesics_cosmological_background.md` | 16 32416 347d377321aba8ddfff05725c4639528c5f
| **15** | `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md` | 13 34013 38675332aacb1edc9f3ff2ed9a3ea6f3080
| **16** | `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md` | 14 17714 20054fd2d65267c46330af30fe0290528d5
| **17** | `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md` | 14 38314 406fce1edf9f49c8d62b94fc6996982b636
| **18** | `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance.md` | 21 17821 1356ba1c66a31368362da9604ff0b15e0ad
| **19** | `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy.md` | 18 88218 905bf6eea32f1a730c04aeec2bc59136cb8
| **20** | `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature.md` | 17 36717 390bedb9d8511f3c287f373f8084f4257e2
| **21** | `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse.md` | 18 43118 01129eb46950ed70c5794ff920e1b578faf
| **22** | `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md` | 31 63331 00542082fc8ccaaceb120e350e8972f3109

### 2.1. Angol nyelvű fordítási fájlok (`step1_*_en.md`–`step22_*_en.md`, `README_en.md`, `MANIFEST_en.txt` — English Translation Files)

| Lépés / Step | Angol fájlnév / English File Name | Méret / Size (B) | SHA-256 (első 16 hex / prefix) | Tartalmi összefoglaló / Description |
|---:|---|---:|---|---|
| **1 (EN)** | `step1_dilaton_form_en.md` | 10 89010 890 83a2e8e88808373983a2e8e888083739
| **2 (EN)** | `step2_conformal_disformal_en.md` | 18 01418 014 abd9db8d104e3bcfabd9db8d104e3bcf
| **3 (EN)** | `step3_yukawa_helmholtz_eom_en.md` | 12 25812 258 8ca1db6f110e26178ca1db6f110e2617
| **4 (EN)** | `step4_covariant_conservation_hydrostatic_en.md` | 11 98611 986 1deaba07182e7cfa1deaba07182e7cfa
| **5 (EN)** | `step5_coupling_dictionary_ppn_screening_en.md` | 14 70114 701 e57d506cdd3ba02be57d506cdd3ba02b
| **6 (EN)** | `step6_sparc_galactic_ngc3198_consistency_en.md` | 13 87013 870 6f75c1f0ccba91de6f75c1f0ccba91de
| **7 (EN)** | `step7_neutron_star_tov_stability_en.md` | 20 54620 546 e219a8cfaf4c81a2e219a8cfaf4c81a2
| **8 (EN)** | `step8_microscopic_clock_coupling_siriusb_lab_v20_en.md` | 15 62015 620 b6caf38ba27692eab6caf38ba27692ea
| **9 (EN)** | `step9_vhe_photon_propagation_hess_dispersion_en.md` | 11 45411 454 dabd6f01ef17584bdabd6f01ef17584b
| **10 (EN)** | `step10_bilateral_kernel_reciprocity_eft_synthesis_en.md` | 12 89012 890 748620c6254121e1748620c6254121e1
| **11 (EN)** | `step11_chameleon_screening_multiscale_eft_synthesis_en.md` | 12 50812 508 375040d0e3dc34c7375040d0e3dc34c7
| **12 (EN)** | `step12_second_variation_coupled_pulsation_hessian_en.md` | 19 01019 010 798f21c32075f3a7798f21c32075f3a7
| **13 (EN)** | `step13_curved_window_metric_variation_synge_world_function_en.md` | 11 23311 233 dddf81914b76e505dddf81914b76e505
| **14 (EN)** | `step14_gravitational_lensing_null_geodesics_cosmological_background_en.md` | 15 51515 515 43051f7ff852c02543051f7ff852c025
| **15 (EN)** | `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20_en.md` | 12 94712 947 8b537673b566d4cb8b537673b566d4cb
| **16 (EN)** | `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress_en.md` | 13 53013 530 f3886101a664595bf3886101a664595b
| **17 (EN)** | `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness_en.md` | 13 76513 765 672dc1adfa28d46f672dc1adfa28d46f
| **18 (EN)** | `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance_en.md` | 20 83220 832 ad7267dc00729932ad7267dc00729932
| **19 (EN)** | `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy_en.md` | 18 06718 067 9fab3349db68588b9fab3349db68588b
| **20 (EN)** | `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature_en.md` | 16 82016 820 800ff0bb69664079800ff0bb69664079
| **21 (EN)** | `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse_en.md` | 17 54217 542 106b9ee27d97f4a4106b9ee27d97f4a4
| **22 (EN)** | `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability_en.md` | 36 89631 863dc84080449d2fcbac89fccc6bcfe0ffe
| **README (EN)** | `README_en.md` | 52 90652 426f7aaf312e782aa9116cdfaa6c46c4051
| **MANIFEST (EN)** | `MANIFEST_en.txt` | — | — | Angol nyelvű V25 fájljegyzék a magyar és az angol fájlok pontos méretével és SHA-256 ellenőrzőösszegével. |

---

## 3. Főbb számszerű eredmények összefoglalása (Summary of Key Numerical Results)

| Fizikai tartomány / Physical Domain | Mennyiség / Quantity | Egzakt érték / Exact Numerical Value | Forráslépés és szkript / Source Step & Script |
|---|---|---|---|
| **Skalár–tenzor / Yukawa szótár (`Branch II-Conf`)** | Kanonikus Yukawa-erősség ($\alpha_Y = 1$) és Brans–Dicke $\omega_{\text{BD}}$ | $\alpha_{\text{std}} = 2\beta^2 = 2\tilde{\alpha}_L = \frac{\alpha_Y}{4\pi} = \frac{1}{4\pi} \approx 0{,}0795775$; $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2} \approx 11{,}06637 \approx 11{,}07$ ($\alpha_Y=1$) ill. $36{,}1991$ ($\alpha_Y=1/3$) | `step5`, `covariant_closure_50aa2db0.py` |
| **SPARC 175 galaxis ($\rho < 0$ pozitivitási teszt)** | Negatív sűrűségű pontok száma (3391 sugárponton) | **0 pont** (fotometriai + 3D top-hat + 3D Yukawa magok pozitív definitségi ellenőrzése) vs. **239 pont** (140 nem-holdout galaxis, $\Upsilon_d=0{,}5$) / **301–303 pont** (175 galaxis, gömbi derivált-képlet) | `step6`, `SPARC_mode.py`, `sparc_tophat.py` |
| **NGC 3198 negatív kontroll ($R = 30{,}08\;\text{kpc}$)** | Megfigyelt $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$ vs. izolált konstans csatolású $C_{\text{pred}}^{\text{isol}}$ | $V_{\text{obs}} = 146 \pm 2\;\text{km/s}$; $\Upsilon_d=1{,}0$: $V_{\text{bar}} = 88{,}53\;\text{km/s}$, **$C_{\text{req}} = 2{,}720$**, de az izolált konstans csatolású modellre **$C_{\text{pred}}^{\text{isol}} \le 1 + \alpha_{\text{std}} = 1{,}0796$** ($V_{\text{pred}}^{\text{isol}} = 88{,}53\text{–}91{,}99\;\text{km/s}$; $\Upsilon_d=0{,}5$: $V_{\text{bar}}=70{,}32\;\text{km/s}$, $C_{\text{req}}=4{,}311$, $V_{\text{pred}}^{\text{isol}} = 70{,}32\text{–}73{,}06\;\text{km/s}$) | `step6`, `SPARC_mode.py`, `data/NGC3198_table.csv` |
| **Kozmológiai horizontskálák ($\Omega_b$ vs. $\Omega_m$)** | $L_{\text{hor}} = \frac{c}{\sqrt{G\rho_{\text{bg}}}} = \sqrt{\frac{8\pi}{3\Omega}}\,\frac{c}{H_0}$ ($H_0 = 67{,}4\;\text{km/s/Mpc}$) | **Barionos háttér ($\Omega_b = 0{,}0493$, $\rho_{b,0} = 4{,}2067\times 10^{-28}\;\text{kg/m}^3$):** $\sqrt{8\pi/(3\Omega_b)} = 13{,}031 \implies L_{\text{hor}}^{(b)} = 1{,}7892\times 10^{27}\;\text{m} \approx \mathbf{57{,}98\;\text{Gpc}}$;<br>**Anyagi háttér ($\Omega_m = 0{,}315$, $\rho_{m,0} = 2{,}686\times 10^{-27}\;\text{kg/m}^3$):** $\sqrt{8\pi/(3\Omega_m)} = 5{,}160 \implies L_{\text{hor}}^{(m)} = 7{,}081\times 10^{26}\;\text{m} \approx \mathbf{22{,}94\;\text{Gpc}}$ | `step14`, `sparc_tophat.py` |
| **Bianchi / Noether megmaradás** | $\max_r |r\cdot\nabla_\mu(T^{\mu}_{\;r,\text{matter}} + T^{\mu}_{\;r,X})|/\varepsilon_c$ ($s=r^2$ PCHIP) | Szimbolikus: `0.0`; folytonos: $1{,}82\times 10^{-11}$ (`tophat_mass`), $1{,}59\times 10^{-11}$ (`tophat_gradient`), $3{,}65\times 10^{-10}$ (`yukawa_action_el`) | `step4`, `step16`, `step21`, `covariant_closure_50aa2db0.py` |
| **Görbült külső Yukawa BVP és ADM-tömeg (`Branch II-Fluid`)** | $R_{\max}$-konvergencia ($[R_*, R_*+N_{\text{mult}}L]$, $N_{\text{mult}}\in\{10,12,14,16\}$), $k_{\text{curved}}$ vs. $k_{\text{flat}}$ és $M_{\text{ADM},\max}$ | $R_{\max}$-szórás: $\Delta k_{\text{curved}}/k_{\text{curved}} = 6{,}17\times 10^{-12}$, $\Delta M_{\text{ADM}}/M_{\text{ADM}} = 5{,}54\times 10^{-13}$; SLy $15{,}0$, $L=10\;\text{km}$, $\alpha_Y=1$: $k_{\text{curved}} = 0{,}233542\;\text{km}^{-1}$ vs. $k_{\text{flat}} = 0{,}187216\;\text{km}^{-1}$ (**$+24{,}74\%$**); $M_{\text{ADM},\max}^{\text{SLy}} = \mathbf{2{,}03321\,M_\odot}$ ($M_{\text{surf},\max} = 2{,}02911\,M_\odot$, $M_{B,\max}^{\text{quad}} = \mathbf{2{,}42339\,M_\odot}$, $\Delta M_{\phi,\text{ext}}^{\text{curved}} = 0{,}00411\,M_\odot$), $R_{1{,}4}^{\text{ADM}} = \mathbf{11{,}5359\;\text{km}}$ ($M_{\text{ADM},\max}^{\text{APR4}} = 2{,}18068\,M_\odot$, $M_{\text{ADM},\max}^{\text{H4}} = 1{,}97809\,M_\odot$) | `step7`, `step22`, `curved_yukawa_adm_and_stability.py` |
| **Csatolt $(\delta h,\delta\phi_L,\delta m)$ stabilitás a görbült `yukawa_action_el` ágon** | $\omega_0^2 = 0$ zérushely vs. $\frac{dM_{\text{ADM}}}{d\rho_c} = 0$ és $\frac{dM_B}{d\rho_c} = 0$ (`scalar_m_sign = -1.0`) | SLy $L=10\;\text{km}$, $\alpha_Y=1$: $\log_{10}\rho_{0c}(M_{\text{ADM},\max}) = 15{,}316404$, $\log_{10}\rho_{0c}(M_{B,\max}^{\text{quad}}) = 15{,}316409$, csatolt $\omega_0^2=0$: **$15{,}316410$** ($\mathbf{+5{,}80\times 10^{-6}\;\text{dex}}$, **$2097\times$-es** csökkenés a befagyasztott-$\phi_L$ $+1{,}22\times 10^{-2}\;\text{dex}$ eltolódásához képest); $\alpha_Y=1/3$: **$15{,}307170$** (**$+4{,}20\times 10^{-6}\;\text{dex}$**, **$995\times$** csökkenés); $L=1\;\text{km}$: **$15{,}303080$** (**$-6{,}86\times 10^{-6}\;\text{dex}$**, **$165\times$** csökkenés); APR4 ($L=10\;\text{km}$): **$15{,}293025$** (**$+3{,}57\times 10^{-6}\;\text{dex}$**, **$2392\times$** csökkenés); H4 ($L=10\;\text{km}$): **$15{,}227783$** (**$-1{,}34\times 10^{-6}\;\text{dex}$**, **$17593\times$** csökkenés) | `step12`, `step22`, `curved_yukawa_adm_and_stability.py` |
| **V21 nemlokális 3D top-hat érintő BVP (`Branch I-A`)** | SLy $L = 10\;\text{km}$ `mass` / `gradient` $\omega_0^2 = 0$ vs. $\partial M/\partial\rho_c = 0$ és $W\to 4W$ skálázás | `mass` módban: $15{,}29931$ vs. $15{,}29924$ ($+6{,}97\times 10^{-5}\;\text{dex}$, $58{,}3\times$ javulás); `gradient` módban: $15{,}30113$ vs. $15{,}30106$ ($+7{,}39\times 10^{-5}\;\text{dex}$, $59{,}3\times$ javulás); tehetetlenség-skálázás: egzaktul $\omega^2_{4W}/\omega^2_W = 0{,}2500$ | `step7`, `step17`, `step18`, `covariant_closure_50aa2db0.py` |
| **Szíriusz B H0 vs. H1 vöröseltolódás** | $v_{\text{obs}} = 80{,}65\pm 0{,}77\;\text{km/s}$ ($M=1{,}018\,M_\odot$, $R=5846{,}6\;\text{km}$) | H0: $|\Delta v|\le 1{,}27\;\text{m/s}$ (nincs korlát); feltételes H1 $2\sigma$ ($\Delta\chi^2=4$) kizárja: $L\in[4{,}08\times 10^6, 3{,}37\times 10^7]\;\text{m}$ ($n=1{,}5$, $\max\Delta\chi^2=52{,}3$ / $72{,}0$) | `step8`, `step19`, `covariant_closure_50aa2db0.py` |
| **Gránitkamra ($8558\;\text{kg}$) IN–OUT teszt** | Feltételes H1 differenciális órajel $\Delta y_{\text{H1}}(\text{IN}-\text{OUT})$ | $L=1{,}00\;\text{m}$-nél: Top-hat $-1{,}436\times 10^{-24}$ ($\equiv 0$ ha $L\ge 2{,}944\;\text{m}$); Yukawa $-3{,}781\times 10^{-25}$ ($L=10\;\text{m}$-nél $-4{,}818\times 10^{-25}$) | `step8`, `step19`, `covariant_closure_50aa2db0.py` |
| **8 tekercses szolenoid ($N=2000$, $14\;\text{mA}$)** | Mágneses csúcsjel ($L_p=0{,}630957\;\text{m}$), $r=L/t$ kvadratúra és $L=0{,}95\;\text{m}$ | $\Delta W(L_p)=1{,}96114452\;\text{J/A}^2$, $k_g = -6{,}13117767\times 10^{-45}\;\text{A}^{-2}$, $\Delta y(L_p) = -1{,}2017108\times 10^{-48} \sim -10^{-48}$; A vs. B út ($L=3\;\text{m}$): $7{,}44\times 10^{-16}$; $L=0{,}95\;\text{m}$ azonosság: $8{,}88\times 10^{-16}\;\text{J/A}^2$ | `step15`, `step20`, `make_revision3_checks.py` |

---


### 3.1. V25 javítások és diagnosztikák: 4-komponensű külső BVP, egzakt külső linearizáció, dinamikus skalárhullám és GL2 N-konvergencia

*(Forrás: `curved_yukawa_adm_and_stability.py`, `curved_yukawa_audit_summary.json`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_adm_sequences.csv`.)*

- **Külső lapse (5. hiba):** a külső $\nu$ negyedik ODE-komponensként együtt integrálódik $(\phi_L,\psi_L,m,\nu)$, `solve_bvp` `tol=1e-12`, végfeltétel $\nu(R_{\max})=\tfrac12\ln(1-2m(R_{\max})/R_{\max})$. SLy, $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $N_{\text{mult}}\in\{10,12,14,16\}$ szórások: $\delta k/k=6{,}17\times 10^{-12}$, $\delta M_{\text{ADM}}/M_{\text{ADM}}=5{,}54\times 10^{-13}$, $\delta\nu_R=1{,}63\times 10^{-14}$.

| $N_{\text{mult}}$ | $R_{\max}$ (km) | $k_{\text{curved}}$ (km$^{-1}$) | $M_{\text{ADM}}$ ($M_\odot$) | $\nu_R$ |
|---|---|---|---|---|
| `10` | `111.4657` | `0.233542329280` | `1.509998138706` | `-0.246306687228` |
| `12` | `131.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `14` | `151.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `16` | `171.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |

- **Külső linearizáció és véges differenciás deriválttesztek (SLy, $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $\alpha_Y=1$):** a `linearized_exterior_response` (593. sor: `a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0`) és a belső skalároperátor bal oldali ($L_u u = R_h\,\delta h$) `a_m_sc` együtthatója (`scalar_m_sign = -1.0`: `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n`) centrális véges differenciákkal összevetve:
  - Külső válasz $\partial\phi_L'(R_*)/\partial m(R_*)$: linearizált BVP **`-1.68365059e-04`** $\text{km}^{-2}$ vs. 4-komponensű külső BVP centrális véges differenciája **`-1.68365059e-04`** $\text{km}^{-2}$ (relatív eltérés: **$3{,}98\times 10^{-10}$**).
  - Belső $r = 6{,}0\;\text{km}$ pontbeli skalár-tömeg derivált ($\partial\phi_L''/\partial m$): analitikus **`-1.72833429e-04`** $\text{km}^{-3}$ vs. belső skalár RHS centrális véges differenciája **`-1.72833428e-04`** $\text{km}^{-3}$ (relatív eltérés: **$6{,}08\times 10^{-9}$**).
- **Dinamikus skalárblokk:** a skalárblokkban szerepel a $-\omega^2 e^{-2\nu}\delta\phi_L$ tag, és a külső választ $\omega$-függő BVP adja (`dynamical_scalar_wave=True`); megoldó: GL2 darabonkénti szegmenses Volterra-diszkretizáció, $N=480$, `scalar_m_sign = -1.0` egzakt LHS előjellel. A kvázisztatikus és a dinamikus gyök a közölt 6 tizedesjegyen egyezik.
- **Fordulópontok, bariontömeg-definíció és zérushelyek:** `minimize_scalar` + 7 pontos kvadratikus illesztés $\pm0{,}002$ dex-en; a zérushelyek a gyök körüli $\pm0{,}0005$ dex-es 3 pontos kvadratikus illesztésből. A táblázatokban szereplő $M_B$ az $s=r^2$ PCHIP-interpolált profilokon végzett adaptív Gauss–Kronrod kvadratúrával (`epsrel=1e-12`) definiált $M_B^{\text{quad}}$ (SLy $L=10\;\text{km}$, $\alpha_Y=1$: $M_{B,\max}^{\text{quad}} = 2{,}42338503\,M_\odot$, csúcs $15{,}316409$, szemben az $N_{\text{grid}}=2401$ trapézszabály $M_B^{\text{trapz}}=2{,}42338667\,M_\odot$ / $15{,}316418$ és az egyszeri menetű 5. ODE-komponens $M_B^{\text{ODE}}=2{,}42203137\,M_\odot$ / $15{,}316522$ értékével).

| Sorozat | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dinamikus $\omega_0^2=0$ | $\Delta_{\text{ADM}}$ dinamikus (dex) | $\Delta_{M_B}$ dinamikus (dex) | Kvázisztatikus $\omega_0^2=0$ | Befagyasztott $\delta\phi_L=0$ | $\Delta_{\text{ADM}}$ befagyasztott (dex) |
|---|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316409` | **`15.316410`** | **`+5.80e-06`** | `+1.33e-06` | `15.316410` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.307169` | **`15.307170`** | **`+4.20e-06`** | `+6.41e-07` | `15.307170` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303082` | **`15.303080`** | **`-6.86e-06`** | `-2.61e-06` | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.293025` | **`15.293025`** | **`+3.57e-06`** | `-1.40e-07` | `15.293025` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227787` | **`15.227783`** | **`-1.34e-06`** | `-4.05e-06` | `15.227783` | `15.251384` | `+2.36e-02` |

**GL2 konvergencia az $\omega_0^2=0$ zérushelyen ($N$ és $R_{\max}$ szerint) és $\log_{10}\rho_{0c}=15{,}311404$-nél** (SLy, $L=10\;\text{km}$, $\alpha_Y=1$): a gyökön $p_{\text{eff,root}}=4{,}248$, Richardson-gyök $\log_{10}\rho_{0c,\text{Rich}}(\omega_0^2=0)=`15.31640981`$ ($N=480$-nál a hiba $6{,}72\times 10^{-8}$ dex; $R_{\max}\in[R_*+10L,R_*+20L]$ közötti szórás $2{,}42\times 10^{-13}$ dex); $\log_{10}\rho_{0c}=15{,}311404$-nél $p_{\text{eff}}=5{,}223$, Richardson-érték $\omega_0^2=6369417{,}611\;\text{s}^{-2}$ ($N=480$-nál a relatív eltérés $1{,}29\times 10^{-5}$).

| $N$ | $\log_{10}\rho_{0c}(\omega_0^2=0)$ | $|\Delta\log_{10}\rho_{0c}|$ a Richardson-gyöktől (dex) | $\Delta_{\text{ADM}}$ (dex) | $\Delta_{M_B}$ (dex) | $\omega_0^2(15{,}311404)$ ($\text{s}^{-2}$) | Rel. eltérés ($\omega_0^2$) |
|---|---|---|---|---|---|---|
| `80` | `15.31645808` | `4.83e-05` | `+5.40e-05` | `+4.95e-05` | `6432300.503336` | `9.87e-03` |
| `120` | `15.31642223` | `1.24e-05` | `+1.82e-05` | `+1.37e-05` | `6385729.603235` | `2.56e-03` |
| `160` | `15.31641426` | `4.46e-06` | `+1.02e-05` | `+5.72e-06` | `6375276.486537` | `9.20e-04` |
| `240` | `15.31641079` | `9.80e-07` | `+6.71e-06` | `+2.24e-06` | `6370703.105233` | `2.02e-04` |
| `360` | `15.31641001` | `2.06e-07` | `+5.94e-06` | `+1.47e-06` | `6369680.451299` | `4.13e-05` |
| `480` | `15.31640987` | `6.72e-08` | `+5.80e-06` | `+1.33e-06` | `6369499.539632` | `1.29e-05` |
| `600` | `15.31640983` | `2.56e-08` | `+5.76e-06` | `+1.29e-06` | `6369446.532399` | `4.54e-06` |
| `800` | `15.31640982` | `7.55e-09` | `+5.74e-06` | `+1.27e-06` | `6369424.048161` | `1.01e-06` |
| `1200` | `15.31640981` | `1.35e-09` | `+5.73e-06` | `+1.26e-06` | `6369418.385568` | `1.22e-07` |

**Értelmezési megjegyzés (az elméleti lezárás határa és feltételes stabilitási záró állítás):** Az egzakt bal oldali előjelű (`scalar_m_sign = -1.0`) GL2-operátor gyöke $N$ és $R_{\max}$ szerint a közölt pontossággal ($<10^{-7}\;\text{dex}$) konvergál (`15.31640981`), és mind az 5 sorozaton $\le 6{,}9\times 10^{-6}\;\text{dex}$-en belül egybeesik az ADM-tömegmaximummal (SLy $L=10\;\text{km}$, $\alpha_Y=1$ esetén $15{,}316410$ vs. $15{,}316404$, $+5{,}80\times 10^{-6}\;\text{dex}$; a `scalar_m_sign = +1.0` diagnosztikai összevetést a `step22` 7.4. szakasza dokumentálja). SLy $L=10\;\text{km}$, $\alpha_Y=1$-re a sajátfüggvény és a statikus érintő relatív $L^2$-eltérései $\delta h = 5{,}04\times 10^{-4}$, $\delta m = 1{,}77\times 10^{-4}$, $\delta\phi_L = 2{,}27\times 10^{-4}$, $\delta\phi_L' = 3{,}09\times 10^{-4}$, míg a szabadfelszíni elmozdulások relatív különbsége $3{,}54\%$ (`rel_diff_xi_R_exact_minus_vs_static_tangent = 0.035406`). A mínuszjeles gyök és a rácskonvergencia komoly előrelépés, de a $3{,}54\%$-os szabadfelszíni eltérés és az ODE vs. kvadratúra bariontömeg-fordulópontok különbsége miatt önmagában nem független bizonyítás a teljes radiális alapspektrumra és a szabadfelszín-variációra; ezért a stabilitási záró állítás addig feltételesen értendő, amíg a variációs és peremfeltételi ellenőrzés függetlenül le nem zárul.

---

## 4. Felhasznált helyi kéziratok (Manuscripts Used — `manuscripts/`)

Minden hivatkozott V19–V23 LaTeX (`.tex`) forrásfájlt hiánytalanul elhelyeztünk a `manuscripts/` almappában:

1. **`manuscripts/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`** és **`manuscripts/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex`**  
   - *Verzió és DOI:* V23 (`50aa2db0`), Zenodo DOI `10.5281/zenodo.23639084`.  
   - *Szerepe a levezetésben:* A két kovariáns hatásjelölt (I. bilaterális mag-hatás és II. lokális Helmholtz–Yukawa skalár–tenzor hatás), a Bianchi-maradékok, a csatolt radiális stabilitási BVP (`eq:factored_diff`, `eq:coupled_bvp`), valamint a H0/H1 óracsatolás (Szíriusz B és $8558\;\text{kg}$-os gránitkamra) elsődleges forrása.
2. **`manuscripts/hu_p1.tex`**, **`manuscripts/hu_p2.tex`**, **`manuscripts/hu_p3.tex`**, **`manuscripts/en_part1.tex`**, **`manuscripts/en_part2.tex`**  
   - *Verzió és DOI:* V22 Book Edition fejezetforrások, Zenodo DOI `10.5281/zenodo.23604423`.  
   - *Szerepe a levezetésben:* A Cauchy-levezetés (2. fejezet), a 3D top-hat és SPARC 175 galaxis elemzés (4–5. fejezet), a V17 H.E.S.S./MAGIC gamma-időzítési audit (6. fejezet), a 7 szegmenses politróp TOV és $s=r^2$ regularizáció (7. fejezet), a csatolt radiális pulzáció (8. fejezet), a Szíriusz B vöröseltolódás (9. fejezet), a $8558\;\text{kg}$-os gránitkamra és a 8 tekercses szolenoid kísérlet (10. fejezet), valamint a numerikus Függelék szó szerinti idézeteinek forrása.
3. **`manuscripts/tal_book_HU_2026-10-24.tex`**, **`manuscripts/tal_book_EN_2026-10-24.tex`**, **`manuscripts/tal_book_DE_2026-10-24.tex`**  
   - *Verzió és DOI:* V22 Book Edition egybefüggő magyar, angol és német LaTeX monográfia-fájlok, Zenodo DOI `10.5281/zenodo.23604423`.
4. **`manuscripts/sparc_tophat_record_2026-10-22_v5.tex`**  
   - *Verzió és DOI:* V19 SPARC 3D gömbi top-hat kézirat, Zenodo DOI `10.5281/zenodo.23497128`.
5. **`manuscripts/tal_v20_fop.tex`**  
   - *Verzió és DOI:* V20 laboratóriumi kétórás és 8 tekercses mágneses kézirat (*Foundations of Physics* formátum), Zenodo DOI `10.5281/zenodo.23540741`.
6. **`manuscripts/tal_v21_fop.tex`**  
   - *Verzió és DOI:* V21 kompakt csillagászati TOV és Szíriusz B kézirat, Zenodo DOI `10.5281/zenodo.23602501`.

---

## 5. Felhasznált Python-programok és SPARC-adatfájlok (Python Programs & SPARC Data Files)

### 5.1. Python-programok (a gyökérmappában és a `code/` almappában)
1. **`curved_yukawa_adm_and_stability.py`** — A `yukawa_action_el` (**Branch II-Fluid**) csatolt belső + görbült külső ($r \in [R_*, R_*+12L]$) Einstein–skalár peremérték-feladatát, az $R_{\max}$-függetlenséget, a görbült vs. lapos Robin-meredekséget, az $M_{\text{ADM}}$-re és $M_B$-re újraszámolt tömeg–sugár görbéket (`curved_yukawa_adm_sequences.csv`) és a csatolt $(\delta h, \delta\phi_L, \delta m)$ dinamikus radiális pulzációs sajátértékeket (`curved_yukawa_coupled_stability.csv`, `curved_yukawa_audit_summary.json`) előállító modul.
2. **`covariant_closure_50aa2db0.py`** — A V23 (`10.5281/zenodo.23639084`) fő számítási és szimbolikus ellenőrző modulja (SymPy Bianchi-azonosságok, I. és II. megfogalmazás TOV-összehasonlítása — beleértve a `yukawa_action_el_curved_adm` sorokat —, csatolt `coupled_bvp` / `sym_galerkin` radiális stabilitás, $W\to 4W$ skálázási invariancia, Szíriusz B és $8558\;\text{kg}$-os gránitkamra H0/H1 óracsatolás).
3. **`radial_modes.py`** — A Chandrasekhar-féle Sturm–Liouville radiális pulzációs operátor (`std`, `modA`, `modB` és `fluid` vs. `metric` tehetetlenségi súly) végeselemes / sávos mátrixú sajátérték-megoldója.
4. **`tal_profiles.py`** — Az $s = r^2$ páros paritású centrális PCHIP-interpolációt (`Audit Fix O1`), a 3D top-hat héj-integrált és a `gradient` módú $r\in[R_*, R_*+L]$ külső vákuum-halo lapse-korrekciót ($\Delta\nu_{\text{ext}}$) megvalósító profilmodul.
5. **`stellar_tal_v21.py`** (és `code/stellar_tal.py`, `code/diagnostics.py`) — A Read et al. (2009) 7 szegmenses politróp állapotegyenleteket (SLy, APR4, H4), a log-entalpia ($h$) TOV-integrátort (769 modell, 33 kauzális végpont) és a fehér törpe (Szíriusz B) integrátort tartalmazó V21/V23 modul.
6. **`make_revision3_checks.py`** (valamint `code/energy_integrator.py`, `code/fields.py`, `code/tophat_finite_source.py`, `code/v20_fixed_case_runner.py`, `code/build_pdf.py`) — A V20 (`10.5281/zenodo.23540741`) 8 tekercses szolenoid elliptikus Biot–Savart ellenőrzését, az $L_p = 0{,}630957\;\text{m}$-es csúcs-újraszámolást, az $r=L/t$ kompaktifikált külső integrált (A és B út) és az $L=0{,}95\;\text{m}$-es belső/külső átkapcsolási azonosságot futtató szkriptek.
7. **`sparkto_b_raw.py`** és **`sparc_tophat.py`** — A 175 SPARC-galaxis 3D tengelyszimmetrikus barionos korong+gáz+dudor modelljének 3D gömbi top-hat átlagsűrűségét ($\bar{\rho}_L(\mathbf{x})$) és $f_L = \exp(-\gamma_L\bar{\rho}_L)$ faktorát kiszámító V19-es program.
8. **`SPARC_mode.py`**, **`ta_pipeline.py`**, **`ta_foundation.py`**, **`sparc_full_runner.py`**, **`tal_factor_7th_model.py`**, **`ta_gammaL.py`**, **`ta_run.py`**, **`ta_run_l0.py`**, **`ta_run_l0_v2.py`** — Az SPARC 175 / 155 / 123 galaxisos fotometriai és derivált-alapú (`foundation`) sűrűség-összehasonlítást (0 vs. 239 / 301–303 negatív sűrűségű pont), az NGC 3198 ($R=30{,}08\;\text{kpc}$) kontrollpontot és a V17/V19 7. modell forgásigörbe-illesztéseket futtató programok.

### 5.2. SPARC és V19–V23 adatfájlok (`data/`, `csv/`, `run_output/`)
- **`data/Rotmod_LTG.zip`** (SHA-256: `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`) és **`data/Rotmod_LTG/*_rotmod.dat`** — Mind a 175 SPARC-galaxis eredeti Lelli, McGaugh & Schombert (2016) forgásigörbe-adatfájlja (köztük `data/Rotmod_LTG/NGC3198_rotmod.dat`).
- **`data/SPARC_Lelli2016c.mrt`** (SHA-256: `5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243`) — Az SPARC 1. táblázata (távolságok, inklinációk, $L_{3.6}$, $R_{\text{disk}}$, $M_{\text{HI}}$, $R_{\text{HI}}$, minőségi zászlók).
- **`data/blind_123.txt`** (és `v17_package/github/data/blind_123.txt`, SHA-256: `3637d2f5bda6b2b97e645c5fba5cf39d310a1826bee388286c2c74b9171aa283`) — A V17/V19 hash-zárt 123 galaxisos vakminta listája.
- **`data/NGC3198_table.csv`** és **`data/NGC3198_table.md`** — Az NGC 3198 összes (43) radiális pontjának részletes táblázata (kiemelve az $R=30{,}08\;\text{kpc}$ kontrollpontot).
- **`data/sparc_tophat_per_galaxy_max.csv`** és **`data/sparc_tophat_summary.json`** — A 175 SPARC-galaxisra kiszámított V19 3D gömbi top-hat csúcsértékek és Monte Carlo keresztellenőrzések.
- **`data/v19_7th_model/`** és **`run_output/`** — A V17/V19 barionos alapvonalak (`points_vbar.csv`, `summary_vbar.csv`) és a 7. modell (`points_tal_fL_FITTED_forward.csv`, `fit_scan_L.csv`, `comparison.json`, `run_meta.json`) referencia-kimenetei.
- **`csv/`** és **`data/v23/`** — A V20 laboratóriumi gránitkamra és 8 tekercses szolenoid referenciafájljai, valamint a V23 (`50aa2db0`) és a 22. lépésben készült görbült külső ADM-tömeg és csatolt stabilitás eredményfájljai (`curved_yukawa_adm_sequences.csv`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_audit_summary.json`, `conservation_Residuals.json`, `closure_tov_comparison.csv`, `coupled_stability_results.csv`, `clock_coupling_H0_H1.csv`).
- **`data/v21/`** és **`observational_inputs.csv`** — A V21/V23 769 csillagmodell (`mr_sequences.csv`, `stellar_summary.csv`, `causal_endpoint_convergence.csv`, `profiles_L10km.csv`, `radius_shifts.csv`, `sirius_redshift.csv`, `numerical_diagnostics.json`, `wd_summary.json`) és a NICER/LIGO/pulsar megfigyelési bemenetek.

---

## 6. Pontos reprodukciós utasítások (Exact Reproduction Instructions)

### 6.1. Szoftverkövetelmények (Software Dependencies)
- **Operációs rendszer:** Linux / macOS / POSIX kompatibilis környezet  
- **Python verzió:** `Python >= 3.10` (tesztelve Python 3.11–3.13 alatt; a lokális munkaállomáson a `/workspace/scratch/venv_openitems/bin/python3` virtuális környezet minden tudományos csomagot előre telepítve tartalmaz)  
- **Szükséges Python-csomagok:** `numpy`, `scipy`, `pandas`, `sympy`, `matplotlib`  
  Telepítés (ha szükséges):
  ```bash
  pip install numpy scipy pandas sympy matplotlib
  ```

### 6.2. A programok futtatása és a várt kimenetek (Running Each Program & Expected Outputs)
Minden alábbi parancs közvetlenül a kicsomagolt `scalar-tensor-yukawa-fit` munkakönyvtárból futtatható hálózati kapcsolat nélkül:

1. **Görbült külső `yukawa_action_el` Einstein–skalár BVP, ADM-tömeg ($M_{\text{ADM}}$), tömeg–sugár görbék és csatolt $(\delta h, \delta\phi_L, \delta m)$ radiális stabilitás (`curved_yukawa_adm_and_stability.py`):**
   ```bash
   python3 curved_yukawa_adm_and_stability.py
   ```
   - *Várt kimenet (~2–3 perc):*  
     - Létrehozza a `curved_yukawa_adm_sequences.csv`, `curved_yukawa_coupled_stability.csv` és `curved_yukawa_audit_summary.json` fájlokat (valamint frissíti a `data/v23/` almappát).  
     - Ellenőrzi az $R_{\max}$-függetlenséget ($\Delta k_{\text{curved}}/k_{\text{curved}} = 6{,}17\times 10^{-12}$, $\Delta M_{\text{ADM}}/M_{\text{ADM}} = 5{,}54\times 10^{-13}$), a görbült vs. lapos Robin-meredekséget ($+24{,}74\%$ SLy $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $\alpha_Y=1$ esetén), az ADM-tömegmaximumokat ($M_{\text{ADM},\max}^{\text{SLy},10\text{km}} = 2{,}03321\,M_\odot$, $R_{1{,}4}^{\text{ADM}} = 11{,}5359\;\text{km}$) a véges differenciás deriváltteszteket ($\partial\phi_L'(R_*)/\partial m(R_*)=-1{,}68365\times 10^{-4}\;\text{km}^{-2}$, belső $r=6\;\text{km}$ `a_m_sc` $=-1{,}72833\times 10^{-4}\;\text{km}^{-3}$) és a csatolt stabilitási zérushelyeket (`scalar_m_sign = -1.0`: $\log_{10}\rho_{0c}(\omega_0^2=0)=15{,}316410$, $\Delta_{\text{ADM}} = +5{,}80\times 10^{-6}\;\text{dex}$, $L=1\;\text{km}$-nél $-6{,}86\times 10^{-6}\;\text{dex}$).

2. **SPARC 175 galaxis ($\rho < 0$ audit), NGC 3198 ($R = 30{,}08\;\text{kpc}$) kontrollpont és V19 3D top-hat csúcsok ellenőrzése (`SPARC_mode.py`):**
   ```bash
   python3 SPARC_mode.py --mode summary
   ```
   - *Várt kimenet (JSON a standard kimeneten, ~0,2 mp):*  
     - `"photometric_and_3D_tophat_yukawa_rho_le_0_points": 0` vs. `"derivative_formula_rho_lt_0_points_visible140_Ydisk_0p5": 239` (és `301` mind a 175 galaxison).  
     - `"NGC3198_control_point"`: `"R_kpc": 30.08`, `"V_obs_km_s": 146.0`, `"V_bar_Ydisk_1p0_km_s": 88.53`, `"C_req_Ydisk_1p0": 2.72`, `"V_bar_Ydisk_0p5_km_s": 70.32`, `"C_req_Ydisk_0p5": 4.311`.

3. **SPARC 3D gömbi top-hat átlagsűrűség és $f_L$ faktor kiszámítása (`sparkto_b_raw.py` / `sparc_tophat.py`):**
   ```bash
   python3 sparkto_b_raw.py --galaxies NGC3198,UGC02487,UGC03546 --skip_checks --out ./output_tophat_check
   ```

4. **V20 nyolctekercses mágnesestér-kísérlet, elliptikus Biot–Savart, $r=L/t$ külső kvadratúra és $L=0{,}95\;\text{m}$-es azonosság (`make_revision3_checks.py`):**
   ```bash
   python3 make_revision3_checks.py . ./csv/revision3_checks_reproduced.json
   ```

5. **V23 kovariáns hatáslezárás, Bianchi-maradékok, csatolt TOV stabilitási BVP, Szíriusz B és $8558\;\text{kg}$-os gránitkamra (`covariant_closure_50aa2db0.py`):**
   ```bash
   python3 covariant_closure_50aa2db0.py
   ```

6. **1–22. levezetési lépések és az összes csomagfájl SHA-256 integritásának ellenőrzése:**
   ```bash
   sha256sum -c SHA256SUMS.txt
   ```
