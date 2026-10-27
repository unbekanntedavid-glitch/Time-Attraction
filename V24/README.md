# Skalár–tenzor és Yukawa/diszform hatáselv-illesztés (V23 $\leftrightarrow$ Jordan–Brans–Dicke / Dilaton keretrendszer) — 1–22. lépéses levezetési és reprodukciós csomag
# Scalar-Tensor and Yukawa/Disformal Action Mapping (V23 $\leftrightarrow$ Jordan–Brans–Dicke / Dilaton Framework) — 22-Step Derivation and Audit Reproduction Package

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Munkakönyvtár / Working Directory:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit`  
**Önálló bírálati archívum / Self-Contained Audit Archive:** `scalar-tensor-yukawa-fit_audit_package.zip`  
**Licenc / License:** CC BY-NC 4.0 — Kereskedelmi hasznosítás kizárólag a szerző előzetes írásbeli engedélyével.  
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

Az alábbi táblázat felsorolja a `/workspace/ta-v27-work/scalar-tensor-yukawa-fit` mappában található **22 levezetési dokumentumot** (`step1`–`step22`), azok pontos fájlnevét, méretét, SHA-256 ellenőrzőösszegét és egymondatos tartalmi összefoglalóját:

| Lépés / Step | Fájlnév / File Name | Méret / Size (B) | SHA-256 (első 16 hex / prefix) | Egymondatos tartalmi összefoglaló (HU / EN One-Line Gist) |
|---:|---|---:|---|---|
| **1** | `step1_dilaton_form.md` | 11 282 | `6d12c0fd49bb3a50` | A Cauchy-additivitásból ($f(\rho_1+\rho_2)=f(\rho_1)f(\rho_2)$) és a $C = 4\pi G/c^2$ térfogati kapacitásból levezeti az $f_L = \exp(-\gamma_L\bar{\rho}_L) = e^{-\phi_L}$ alak és a húrelméleti dilaton-csatolás ($A(\phi)=e^{\beta\phi}$) matematikai ekvivalenciáját. / Derives the exact equivalence between the Cauchy-derived exponential factor $f_L = \exp(-\gamma_L\bar{\rho}_L) = e^{-\phi_L}$ and the string-theory dilaton matter coupling $A(\phi)=e^{\beta\phi}$. |
| **2** | `step2_conformal_disformal.md` | 18 425 | `5a9fe423508e3329` | Levezeti a konform ($\tilde{g}_{\mu\nu}^{\text{conf}} = f_L^2 g_{\mu\nu}$, $\sqrt{-\tilde{g}}=f_L^4\sqrt{-g}$) és a közeghez kötött diszform ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$, $\sqrt{-\tilde{g}}=f_L\sqrt{-g}$) transzformációkat és a modellágak taxonómiáját. / Derives the Jacobian, inverse metric, comoving-clock equivalence, and explicit branch taxonomy of the conformal and disformal metric frames. |
| **3** | `step3_yukawa_helmholtz_eom.md` | 12 670 | `9ca162a9375eb71c` | A V23 II. megfogalmazásából variációs elvvel levezeti a görbült téridőbeli tömeges Klein–Gordon / Helmholtz–Yukawa-egyenletet (`yukawa_action_el`), a Yukawa Green-függvényt és a lapos téri Robin-közelítést. / Variationally derives the curved-spacetime massive Klein–Gordon / Helmholtz–Yukawa field equation, Green function, and flat-space Robin approximation at $r=R_*$. |
| **4** | `step4_covariant_conservation_hydrostatic.md` | 12 379 | `8f5446c737f3c39f` | Diffeomorfizmus-invarianciából (Noether/Bianchi) bebizonyítja a $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$ impulzuscsere-tételt ($\nabla_\mu T^{\mu\nu}_{\text{matter}} = \mathcal{L}_{\text{matter}}\nabla^\nu f_L$) és a módosított Euler–TOV egyensúlyt ($h' = -\nu' + \phi_L'$). / Proves exact covariant stress-energy conservation $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$ and the modified hydrostatic Euler–TOV balance $h' = -\nu' + \phi_L'$. |
| **5** | `step5_coupling_dictionary_ppn_screening.md` | 14 972 | `33ce437c1eec32ba` | Felállítja az egzakt $-\frac{c^4}{2G\alpha_Y}$ normálású skalár–tenzor szótárat ($\alpha_{\text{std}} = 2\beta^2 = \frac{\alpha_Y}{4\pi} = \frac{1}{4\pi} \approx 0{,}07958$, $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2} \approx 11{,}07$ ill. $36{,}20$), és szétválasztja a konform (`Branch II-Conf`) és diszform (`Branch II-Disf`) PPN-viselkedést. / Establishes the exact $-\frac{c^4}{2G\alpha_Y}$ normalized coupling dictionary ($\alpha_{\text{std}} = 1/(4\pi)$, $\omega_{\text{BD}} \approx 11.07$) and separates conformal vs. disformal PPN behavior. |
| **6** | `step6_sparc_galactic_ngc3198_consistency.md` | 14 692 | `1d3c17156dbf5ed4` | Igazolja a 3D magok pozitív definitségét (**0 negatív sűrűségű pont** a 175 SPARC-galaxison vs. 239 a derivált-képlettel), és negatív kontrollként szembeállítja az NGC 3198 ($R = 30{,}08\;\text{kpc}$) megfigyelt $C_{\text{req}} = 2{,}720$ / $4{,}311$ igényét az izolált konstans csatolású modell $C_{\text{pred}}^{\text{isol}} \le 1{,}0796$ ($V_{\text{pred}}^{\text{isol}} = 88{,}53\text{–}91{,}99\;\text{km/s}$) előrejelzésével. / Verifies 0 negative-density points on 175 SPARC galaxies and contrasts NGC 3198's required $C_{\text{req}}=2.720$ against the isolated constant-coupling prediction $C_{\text{pred}}^{\text{isol}} \le 1.0796$ as a negative control. |
| **7** | `step7_neutron_star_tov_stability.md` | 20 067 | `ce336cc89fa20b2b` | Levezeti a V21/V23 görbült téridőbeli TOV-egyenleteket, összehasonlítja a lapos Robin-közelítést ($M_{\text{surf}}$) a görbült külső ADM-tömeggel ($M_{\text{ADM}}$), és szétválasztja a V21 nemlokális érintő BVP-t (`Branch I-A`, $|\Delta\log_{10}\rho_{0c}| = 7\times 10^{-5}\;\text{dex}$) a `yukawa_action_el` csatolt stabilitásától (`Branch II-Fluid`). / Derives the curved-spacetime TOV equations, compares flat-Robin surface mass vs. curved-exterior ADM mass, and separates Branch I-A vs. Branch II-Fluid stability. |
| **8** | `step8_microscopic_clock_coupling_siriusb_lab_v20.md` | 16 329 | `8e3645451604bbc8` | A Dirac–Maxwell hatásból levezeti a H0 (minimális) és feltételes H1 ($\tilde{g}_{00}=f_L^2 g_{00}$) atomóra-frekvenciaeltolódást, a Szíriusz B $2\sigma$-s H1 kizárási sávját ($[4{,}08\times 10^6, 3{,}37\times 10^7]\;\text{m}$), a $8558\;\text{kg}$-os gránitkamra ($\sim 10^{-24}$) és a 8 tekercses szolenoid ($\sim 10^{-48}$) jeleit. / Derives H0 vs. conditional H1 clock coupling from the Dirac–Maxwell action and quantifies Sirius B redshift exclusions, the $8558\;\text{kg}$ granite chamber ($\sim 10^{-24}$), and the 8-coil solenoid ($\sim 10^{-48}$). |
| **9** | `step9_vhe_photon_propagation_hess_dispersion.md` | 12 187 | `9e735e9bfd15e9bc` | Levezeti a diszform metrikabeli $n_t(s)=f_L(s)^{-1}$ akromatikus effektív törésmutatót, a Shapiro-típusú időkésést és a V17 többforrású H.E.S.S./MAGIC gamma-időzítési korlátokat ($\tau_1 = -6{,}1\pm 15{,}3\;\text{s/TeV}$, $|\Delta c_\gamma/c| \lesssim 10^{-21}$). / Derives the achromatic refractive index $n_t(s)=f_L(s)^{-1}$ in the disformal metric and the multi-source H.E.S.S./MAGIC VHE gamma-ray timing audit. |
| **10** | `step10_bilateral_kernel_reciprocity_eft_synthesis.md` | 13 443 | `4e336d901146cfea` | Levezeti a V23 I. megfogalmazásának segédmezős ($\phi_L,\lambda_L$) bilaterális mag-hatását, a lapse-súlyozott reciprocitási potenciált ($V_{\text{eff},N}(r)$), és tisztázza a nulladrendű ($\nabla^2\rho\to 0$) vs. másodrendű ($\frac{1}{10}L^2\nabla^2$ vs. $L^2\nabla^2$) momentum-megfeleltetést. / Derives Formulation I's auxiliary-field bilateral kernel action, reciprocal potential $V_{\text{eff},N}(r)$, and zeroth- vs. second-moment correspondence to Formulation II. |
| **11** | `step11_chameleon_screening_multiscale_eft_synthesis.md` | 13 320 | `3d7262216b373998` | Levezeti a környezetfüggő kaméleon/szimmetron effektív potenciált ($V_{\text{eff}}(\phi_L;\rho)$) és a sűrűségfüggő korrelációs hosszt ($L_{\text{eff}}(\rho)$), kifejezetten külön nemlineáris kiterjesztési hipotézisként (`Branch II-Cham`) kezelve azt a konstans $\alpha_Y$-ú V23 hatáson túl. / Formulates environment-dependent chameleon/symmetron screening (`Branch II-Cham`) as an explicit nonlinear extension hypothesis beyond the constant-$\alpha_Y$ V23 action. |
| **12** | `step12_second_variation_coupled_pulsation_hessian.md` | 18 743 | `1c2f9401c56e140b` | A hatás második variációjából ($\delta^2 S$) és a $(\delta\lambda,\delta m,\delta\nu')$ metrikus kényszerek kiküszöböléséből explicit alakban levezeti a $\mathcal{Q}_{\xi\xi}, \mathcal{Q}_{\phi\phi}, \mathcal{C}_{\xi\phi}, \mathcal{C}_{\phi\xi}$ operátor-együtthatókat, a mozgó felszíni illesztést és a görbült `yukawa_action_el` sorozat numerikus stabilitását ($2{,}04\times 10^{-5}\;\text{dex}$). / Explicitly derives $\mathcal{Q}_{\xi\xi}, \mathcal{Q}_{\phi\phi}, \mathcal{C}_{\xi\phi}, \mathcal{C}_{\phi\xi}$ after metric constraint elimination, moving-surface matching, and numerical stability on the curved `yukawa_action_el` sequence. |
| **13** | `step13_curved_window_metric_variation_synge_world_function.md` | 11 618 | `ec8dcd74bd038ec2` | A Synge-féle világfüggvénnyel ($\sigma(x,y) = \frac{1}{2}s_{\text{geo}}^2$) levezeti a görbült téridőbeli bilaterális top-hat ablak metrikus variációját ($\delta[\sqrt{h}K_L]/\delta g^{\mu\nu}$) és a $\delta_D(L-s_{\text{geo}})$ Dirac-delta határfeszültség-járulékot. / Derives the exact metric variation $\delta[\sqrt{h}K_L]/\delta g^{\mu\nu}$ of the curved-spacetime top-hat window via Synge's world function and the $\delta_D(L-s_{\text{geo}})$ surface stress. |
| **14** | `step14_gravitational_lensing_null_geodesics_cosmological_background.md` | 16 347 | `05725c4639528c5f` | Levezeti a konform vs. diszform gravitációs fényelhajlást, az FLRW kozmológiai háttérfaktor kiesését $\nabla\ln f_L$-ből, és egzaktul kiszámítja mind a barionos ($\Omega_b = 0{,}0493 \implies \sqrt{8\pi/(3\Omega_b)} = 13{,}031$, $L_{\text{hor}}^{(b)} \approx 57{,}98\;\text{Gpc}$), mind az összesített anyagi ($\Omega_m = 0{,}315 \implies \sqrt{8\pi/(3\Omega_m)} = 5{,}160$, $L_{\text{hor}}^{(m)} \approx 22{,}94\;\text{Gpc}$) horizontskálát. / Derives conformal vs. disformal light deflection, FLRW cancellation in $\nabla\ln f_L$, and exact horizon scales for both $\Omega_b=0.0493$ ($57.98\;\text{Gpc}$) and $\Omega_m=0.315$ ($22.94\;\text{Gpc}$). |
| **15** | `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md` | 13 386 | `ff2ed9a3ea6f3080` | Megmutatja, hogy bár a klasszikus Maxwell-tenzor nyommentes ($T^\mu_{\;\mu,\text{EM}}=0$), a közeghez kötött diszform csatolás ($u_\mu u_\nu T^{\mu\nu}_{\text{EM}} = u_{\text{EM}}$) közvetlenül a $\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$ mágneses energiasűrűséggel forrásozza $\phi_L$-t, és bizonyítja a V20 polaritás-páros kereszttag-kiejtését. / Proves that while $T^\mu_{\;\mu,\text{EM}}=0$ suppresses conformal trace coupling, disformal coupling $u_\mu u_\nu T^{\mu\nu}_{\text{EM}}$ sources $\phi_L$ via $\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$ with exact polarity-even cross-term cancellation. |
| **16** | `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md` | 14 200 | `0af30fe0290528d5` | A Bianchi-azonosságból ($\nabla_\mu G^{\mu\nu}\equiv 0$) levezeti a V21/V23 algebrai X-szektorbeli feszültségrekonstrukciót (`tophat_mass` és `tophat_gradient`, $\max|p_{tX}/P_c|=0{,}140$), és összeveti a dinamikus Yukawa skalár-feszültségtenzorral ($T^\mu_{\;\nu}(\phi)$). / Derives the V21/V23 algebraic X-sector anisotropic stress reconstruction from the Bianchi identity ($\max|p_{tX}/P_c|=0.140$) and compares it with the dynamical Yukawa stress tensor. |
| **17** | `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md` | 14 406 | `b94fc6996982b636` | Egzaktul levezeti a V23 (8) faktorizált Chandrasekhar–hidrosztatikai különbség-azonosságát (`eq:factored_diff`) és bebizonyítja, hogy `gradient` módban a `metric` tehetetlenségi súly ($\mathcal{W}_{\text{metric}} \propto f_L A_{\text{GR}} - A_{\text{GR}} \le 0$) $L \ge 367\;\text{m}$-nél nem-pozitívvá válik. / Analytically derives the V23 Eq. (8) factorized Chandrasekhar–hydrostatic difference identity and proves the non-positivity of the `gradient`-mode `metric` inertia weight for $L \ge 367\;\text{m}$. |
| **18** | `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance.md` | 20 862 | `0ebe52095efe2a94` | Levezeti a V23 (9) nemlokális statikus érintőoperátort (`eq:coupled_bvp` a V21 `Branch I-A` háttéren), a 3D gömbhéj-metszeti $K_{ij}$ magot ($1{,}214\times 10^{-6}$ pontosság), a Rayleigh–Ritz/Galerkin szimmetrizálást és a $W \to \kappa W \implies \omega^2_{4W}/\omega^2_W = 0{,}2500$ skálázási invarianciát. / Derives the V23 Eq. (9) nonlocal static tangent BVP on the V21 background, the 3D spherical-shell kernel $K_{ij}$, and the exact inertia-scaling invariance $\omega^2_{4W}/\omega^2_W = 0.2500$. |
| **19** | `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy.md` | 18 905 | `4aeec2bc59136cb8` | Levezeti a V20 $8558\;\text{kg}$-os gránitkamra 3D IN–OUT differenciális kvadratúráját a három geometriai $L$-tartományban ($\Delta y_{\text{H1}}^{\text{top-hat}}(1\;\text{m}) = -1{,}436\times 10^{-24}$, $\equiv 0$ ha $L\ge 2{,}944\;\text{m}$; Yukawa: $-3{,}781\times 10^{-25}$), valamint a Szíriusz B–Szíriusz A `common_background` vs. `fixedM` vöröseltolódást. / Derives the 3D IN–OUT differential quadrature of the $8558\;\text{kg}$ granite chamber across three geometric $L$-regimes and the Sirius B–Sirius A `common_background` vs. `fixedM` redshift proxy. |
| **20** | `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature.md` | 17 390 | `f373f8084f4257e2` | Levezeti a V20 8 tekercses ($N=2000$, $D_{\text{eff}}=291{,}7\;\text{mm}$, $I_0=14\;\text{mA}$) elliptikus Biot–Savart térszámítását, az $L_p = 0{,}630957\;\text{m}$-es csúcsjelet ($\Delta y = -1{,}2017\times 10^{-48}$), az $r=L/t$ kompaktifikált külső kvadratúrát ($7{,}44\times 10^{-16}$ A–B egyezés) és az $L=0{,}95\;\text{m}$-es átkapcsolási azonosságot. / Derives the 8-coil elliptic Biot–Savart field, peak signal $\Delta y(L_p) = -1.2017\times 10^{-48}$, $r=L/t$ compactified exterior quadrature ($7.44\times 10^{-16}$ agreement), and the $L=0.95\;\text{m}$ switching identity. |
| **21** | `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse.md` | 18 011 | `94ff920e1b578faf` | Levezeti a Read et al. (2009) 7 szegmenses politróp log-entalpia ($h$) TOV-transzformációt, az $s=r^2$ páros paritású centrális PCHIP-regularizációt (`Audit Fix O1`, $10^{-11}$ Bianchi-reziduum), a `gradient` módú $r\in[R_*,R_*+L]$ külső vákuum-halo lapse-integrált ($\Delta\nu_{\text{ext}}$), és szembeállítja a V21 $h'=-\nu'$ egyensúlyt a `yukawa_action_el` $h'=-\nu'+\phi_L'$ egyenletével. / Derives the 7-segment piecewise-polytrope log-enthalpy TOV transformation, $s=r^2$ even-parity PCHIP regularization, exterior vacuum-halo lapse $\Delta\nu_{\text{ext}}$, and contrasts V21 $h'=-\nu'$ against Yukawa $h'=-\nu'+\phi_L'$. |
| **22** | `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md` | 16 654 | `dce145dee013eb92` | Egzaktul megoldja a `yukawa_action_el` ág görbült külső Einstein–skalár vákuumrendszerét ($r \in [R_*,\infty)$), meghatározza az aszimptotikus ADM-tömeget ($M_{\text{ADM}} = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$), újraszámolja a tömeg–sugár görbéket (SLy, APR4, H4; $L=1, 10\;\text{km}$; $\alpha_Y=1, 1/3$) és a csatolt $(\delta h, \delta\phi_L, \delta m)$ radiális pulzációs stabilitást ($2{,}04\times 10^{-5}\;\text{dex}$ egyezés). / Solves the curved exterior Einstein–scalar vacuum BVP of `yukawa_action_el`, determines $M_{\text{ADM}}$, recomputes all mass–radius curves, and solves the coupled $(\delta h,\delta\phi_L,\delta m)$ radial stability problem. |

### 2.1. Angol nyelvű fordítási fájlok (`step1_*_en.md`–`step20_*_en.md`, `README_en.md`, `MANIFEST_en.txt` — English Translation Files)

| Lépés / Step | Angol fájlnév / English File Name | Méret / Size (B) | SHA-256 (első 16 hex / prefix) | Tartalmi összefoglaló / Description |
|---:|---|---:|---|---|
| **1 (EN)** | `step1_dilaton_form_en.md` | 10 967 | `2c9a146a645f68da` | English translation of `step1_dilaton_form.md`: Cauchy-derived exponential factor $f_L = \exp(-\gamma_L\bar{\rho}_L) = e^{-\phi_L}$ and string-theory dilaton matter coupling $A(\phi)=e^{\beta\phi}$. |
| **2 (EN)** | `step2_conformal_disformal_en.md` | 18 037 | `11e46f3c0e249a66` | English translation of `step2_conformal_disformal.md`: conformal and disformal metric frames, Jacobian, inverse metric, comoving-clock equivalence, and explicit branch taxonomy. |
| **3 (EN)** | `step3_yukawa_helmholtz_eom_en.md` | 12 281 | `39d4c08f70704587` | English translation of `step3_yukawa_helmholtz_eom.md`: curved-spacetime massive Klein–Gordon / Helmholtz–Yukawa field equation, Green function, and flat-space Robin approximation. |
| **4 (EN)** | `step4_covariant_conservation_hydrostatic_en.md` | 12 009 | `3f3ac1ebe1b6a984` | English translation of `step4_covariant_conservation_hydrostatic.md`: covariant stress-energy conservation $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$ and modified Euler–TOV balance. |
| **5 (EN)** | `step5_coupling_dictionary_ppn_screening_en.md` | 14 724 | `1c281eb0af68813f` | English translation of `step5_coupling_dictionary_ppn_screening.md`: exact $-\frac{c^4}{2G\alpha_Y}$ normalized coupling dictionary ($\alpha_{\text{std}} = 1/(4\pi)$, $\omega_{\text{BD}} \approx 11.07$) and PPN behavior. |
| **6 (EN)** | `step6_sparc_galactic_ngc3198_consistency_en.md` | 13 947 | `ea83e05ff331e3ef` | English translation of `step6_sparc_galactic_ngc3198_consistency.md`: 0 negative-density points on 175 SPARC galaxies and NGC 3198 ($R=30.08\;\text{kpc}$) negative control. |
| **7 (EN)** | `step7_neutron_star_tov_stability_en.md` | 19 506 | `d0be25723a09fdcb` | English translation of `step7_neutron_star_tov_stability.md`: curved-spacetime TOV equations, flat-Robin surface mass vs. curved-exterior ADM mass, and Branch I-A vs. II-Fluid stability. |
| **8 (EN)** | `step8_microscopic_clock_coupling_siriusb_lab_v20_en.md` | 15 643 | `3332edd07261c63f` | English translation of `step8_microscopic_clock_coupling_siriusb_lab_v20.md`: H0 vs. conditional H1 clock coupling, Sirius B redshift, $8558\;\text{kg}$ granite chamber, and 8-coil solenoid. |
| **9 (EN)** | `step9_vhe_photon_propagation_hess_dispersion_en.md` | 11 477 | `003941bcc65c727f` | English translation of `step9_vhe_photon_propagation_hess_dispersion.md`: achromatic refractive index $n_t(s)=f_L(s)^{-1}$ in the disformal metric and H.E.S.S./MAGIC VHE timing audit. |
| **10 (EN)** | `step10_bilateral_kernel_reciprocity_eft_synthesis_en.md` | 12 913 | `399f6fce0c03dc32` | English translation of `step10_bilateral_kernel_reciprocity_eft_synthesis.md`: Formulation I auxiliary-field bilateral kernel action, reciprocal potential $V_{\text{eff},N}(r)$, and moment matching. |
| **11 (EN)** | `step11_chameleon_screening_multiscale_eft_synthesis_en.md` | 12 531 | `4ad3ed4131afd2be` | English translation of `step11_chameleon_screening_multiscale_eft_synthesis.md`: environment-dependent chameleon/symmetron screening (`Branch II-Cham`) nonlinear extension hypothesis. |
| **12 (EN)** | `step12_second_variation_coupled_pulsation_hessian_en.md` | 18 198 | `919d54ad8e3455a7` | English translation of `step12_second_variation_coupled_pulsation_hessian.md`: second variation $\delta^2 S$, explicit $(\xi,\delta\phi_L)$ pulsation operator coefficients, and positive-definite kinetic Hessian. |
| **13 (EN)** | `step13_curved_window_metric_variation_synge_world_function_en.md` | 11 256 | `e2ae56ad1b197b06` | English translation of `step13_curved_window_metric_variation_synge_world_function.md`: metric variation $\delta[\sqrt{h}K_L]/\delta g^{\mu\nu}$ via Synge's world function and $\delta_D(L-d_\Sigma)$ boundary stress. |
| **14 (EN)** | `step14_gravitational_lensing_null_geodesics_cosmological_background_en.md` | 15 538 | `632373a149b21224` | English translation of `step14_gravitational_lensing_null_geodesics_cosmological_background.md`: conformal vs. disformal light deflection, NICER ray-tracing, and $\Omega_b$ ($57.98\;\text{Gpc}$) vs. $\Omega_m$ ($22.94\;\text{Gpc}$) horizon scales. |
| **15 (EN)** | `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20_en.md` | 12 993 | `d038c9b3de65558a` | English translation of `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md`: electromagnetic trace-freeness vs. disformal magnetic source $\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$ and polarity-even cross-term cancellation. |
| **16 (EN)** | `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress_en.md` | 13 553 | `24622e0f5d74bbe4` | English translation of `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md`: V21/V23 algebraic X-sector anisotropic stress reconstruction from the Bianchi identity vs. dynamical Yukawa stress tensor. |
| **17 (EN)** | `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness_en.md` | 13 788 | `9a04ce99cd38bbe2` | English translation of `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md`: V23 Eq. (8) factorized Chandrasekhar–hydrostatic difference identity and `gradient`-mode `metric` inertia-weight indefiniteness. |
| **18 (EN)** | `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance_en.md` | 20 532 | `d04a05e7c0507687` | English translation of `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance.md`: V23 Eq. (9) nonlocal static tangent BVP, 3D spherical-shell kernel $K_{ij}$, and $W\to\kappa W$ inertia-scaling invariance. |
| **19 (EN)** | `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy_en.md` | 18 090 | `0b4953fec2eb9fce` | English translation of `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy.md`: $8558\;\text{kg}$ granite-chamber 3D IN–OUT quadrature across three geometric $L$-regimes and Sirius B–Sirius A common proxy. |
| **20 (EN)** | `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature_en.md` | 16 843 | `ec747ce9b2db70b9` | English translation of `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature.md`: V20 8-coil elliptic Biot–Savart loop integral, $r=L/t$ compactified exterior quadrature, and $L=0.95\;\text{m}$ switching identity. |
| **README (EN)** | `README_en.md` | 45 954 | `70374049f360b95f` | Angol nyelvű áttekintő és reprodukciós dokumentáció (`README.md` angol megfelelője). / Full English edition of `README.md`. |
| **MANIFEST (EN)** | `MANIFEST_en.txt` | — | — | Angol nyelvű fájljegyzék a magyar és a 22 angol fájl pontos méretével és SHA-256 hash-ével. / English file manifest with exact sizes and SHA-256 hashes. |

---

## 3. Főbb számszerű eredmények összefoglalása (Summary of Key Numerical Results)

| Fizikai tartomány / Physical Domain | Mennyiség / Quantity | Egzakt érték / Exact Numerical Value | Forráslépés és szkript / Source Step & Script |
|---|---|---|---|
| **Skalár–tenzor / Yukawa szótár (`Branch II-Conf`)** | Kanonikus Yukawa-erősség ($\alpha_Y = 1$) és Brans–Dicke $\omega_{\text{BD}}$ | $\alpha_{\text{std}} = 2\beta^2 = 2\tilde{\alpha}_L = \frac{\alpha_Y}{4\pi} = \frac{1}{4\pi} \approx 0{,}0795775$; $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2} \approx 11{,}06637 \approx 11{,}07$ ($\alpha_Y=1$) ill. $36{,}1991$ ($\alpha_Y=1/3$) | `step5`, `covariant_closure_50aa2db0.py` |
| **SPARC 175 galaxis ($\rho < 0$ pozitivitási teszt)** | Negatív sűrűségű pontok száma (3391 sugárponton) | **0 pont** (fotometriai + 3D top-hat + 3D Yukawa magok pozitív definitségi ellenőrzése) vs. **239 pont** (140 nem-holdout galaxis, $\Upsilon_d=0{,}5$) / **301–303 pont** (175 galaxis, gömbi derivált-képlet) | `step6`, `SPARC_mode.py`, `sparc_tophat.py` |
| **NGC 3198 negatív kontroll ($R = 30{,}08\;\text{kpc}$)** | Megfigyelt $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$ vs. izolált konstans csatolású $C_{\text{pred}}^{\text{isol}}$ | $V_{\text{obs}} = 146 \pm 2\;\text{km/s}$; $\Upsilon_d=1{,}0$: $V_{\text{bar}} = 88{,}53\;\text{km/s}$, **$C_{\text{req}} = 2{,}720$**, de az izolált konstans csatolású modellre **$C_{\text{pred}}^{\text{isol}} \le 1 + \alpha_{\text{std}} = 1{,}0796$** ($V_{\text{pred}}^{\text{isol}} = 88{,}53\text{–}91{,}99\;\text{km/s}$; $\Upsilon_d=0{,}5$: $V_{\text{bar}}=70{,}32\;\text{km/s}$, $C_{\text{req}}=4{,}311$, $V_{\text{pred}}^{\text{isol}} = 70{,}32\text{–}73{,}06\;\text{km/s}$) | `step6`, `SPARC_mode.py`, `data/NGC3198_table.csv` |
| **Kozmológiai horizontskálák ($\Omega_b$ vs. $\Omega_m$)** | $L_{\text{hor}} = \frac{c}{\sqrt{G\rho_{\text{bg}}}} = \sqrt{\frac{8\pi}{3\Omega}}\,\frac{c}{H_0}$ ($H_0 = 67{,}4\;\text{km/s/Mpc}$) | **Barionos háttér ($\Omega_b = 0{,}0493$, $\rho_{b,0} = 4{,}2067\times 10^{-28}\;\text{kg/m}^3$):** $\sqrt{8\pi/(3\Omega_b)} = 13{,}031 \implies L_{\text{hor}}^{(b)} = 1{,}7892\times 10^{27}\;\text{m} \approx \mathbf{57{,}98\;\text{Gpc}}$;<br>**Anyagi háttér ($\Omega_m = 0{,}315$, $\rho_{m,0} = 2{,}686\times 10^{-27}\;\text{kg/m}^3$):** $\sqrt{8\pi/(3\Omega_m)} = 5{,}160 \implies L_{\text{hor}}^{(m)} = 7{,}081\times 10^{26}\;\text{m} \approx \mathbf{22{,}94\;\text{Gpc}}$ | `step14`, `sparc_tophat.py` |
| **Bianchi / Noether megmaradás** | $\max_r |r\cdot\nabla_\mu(T^{\mu}_{\;r,\text{matter}} + T^{\mu}_{\;r,X})|/\varepsilon_c$ ($s=r^2$ PCHIP) | Szimbolikus: `0.0`; folytonos: $1{,}82\times 10^{-11}$ (`tophat_mass`), $1{,}59\times 10^{-11}$ (`tophat_gradient`), $3{,}65\times 10^{-10}$ (`yukawa_action_el`) | `step4`, `step16`, `step21`, `covariant_closure_50aa2db0.py` |
| **Görbült külső Yukawa BVP és ADM-tömeg (`Branch II-Fluid`)** | $R_{\max}$-konvergencia ($[R_*, R_*+N_{\text{mult}}L]$, $N_{\text{mult}}\in\{10,12,14,16\}$), $k_{\text{curved}}$ vs. $k_{\text{flat}}$ és $M_{\text{ADM},\max}$ | $R_{\max}$-szórás: $\Delta k_{\text{curved}}/k_{\text{curved}} = 6{,}17\times 10^{-12}$, $\Delta M_{\text{ADM}}/M_{\text{ADM}} = 5{,}54\times 10^{-13}$; SLy $15{,}0$, $L=10\;\text{km}$, $\alpha_Y=1$: $k_{\text{curved}} = 0{,}233542\;\text{km}^{-1}$ vs. $k_{\text{flat}} = 0{,}187216\;\text{km}^{-1}$ (**$+24{,}74\%$**); $M_{\text{ADM},\max}^{\text{SLy}} = \mathbf{2{,}03321\,M_\odot}$ ($M_{\text{surf},\max} = 2{,}02911\,M_\odot$, $\Delta M_{\phi,\text{ext}}^{\text{curved}} = 0{,}00411\,M_\odot$), $R_{1{,}4}^{\text{ADM}} = \mathbf{11{,}5359\;\text{km}}$ ($M_{\text{ADM},\max}^{\text{APR4}} = 2{,}18068\,M_\odot$, $M_{\text{ADM},\max}^{\text{H4}} = 1{,}97809\,M_\odot$) | `step7`, `step22`, `curved_yukawa_adm_and_stability.py` |
| **Csatolt $(\delta h,\delta\phi_L,\delta m)$ stabilitás a görbült `yukawa_action_el` ágon** | $\omega_0^2 = 0$ zérushely vs. $\frac{dM_{\text{ADM}}}{d\rho_c} = 0$ és $\frac{dM_B}{d\rho_c} = 0$ | SLy $L=10\;\text{km}$, $\alpha_Y=1$: $\log_{10}\rho_{0c}(M_{\text{ADM},\max}) = 15{,}316404$, $\log_{10}\rho_{0c}(M_{B,\max}) = 15{,}316388$, csatolt $\omega_0^2=0$: **$15{,}316425$** ($\mathbf{+2{,}04\times 10^{-5}\;\text{dex}}$, **$595\times$-ös** javulás a befagyasztott-$\phi_L$ $+1{,}22\times 10^{-2}\;\text{dex}$ hibájához képest); $\alpha_Y=1/3$: **$+5{,}84\times 10^{-6}\;\text{dex}$** (**$715\times$** javulás); $L=1\;\text{km}$: **$-6{,}86\times 10^{-6}\;\text{dex}$** (**$165\times$** javulás); APR4 ($L=10\;\text{km}$): **$+2{,}07\times 10^{-5}\;\text{dex}$** (**$413\times$** javulás); H4 ($L=10\;\text{km}$): **$+5{,}71\times 10^{-6}\;\text{dex}$** (**$4137\times$** javulás) | `step12`, `step22`, `curved_yukawa_adm_and_stability.py` |
| **V21 nemlokális 3D top-hat érintő BVP (`Branch I-A`)** | SLy $L = 10\;\text{km}$ `mass` / `gradient` $\omega_0^2 = 0$ vs. $\partial M/\partial\rho_c = 0$ és $W\to 4W$ skálázás | `mass` módban: $15{,}29931$ vs. $15{,}29924$ ($+6{,}97\times 10^{-5}\;\text{dex}$, $58{,}3\times$ javulás); `gradient` módban: $15{,}30113$ vs. $15{,}30106$ ($+7{,}39\times 10^{-5}\;\text{dex}$, $59{,}3\times$ javulás); tehetetlenség-skálázás: egzaktul $\omega^2_{4W}/\omega^2_W = 0{,}2500$ | `step7`, `step17`, `step18`, `covariant_closure_50aa2db0.py` |
| **Szíriusz B H0 vs. H1 vöröseltolódás** | $v_{\text{obs}} = 80{,}65\pm 0{,}77\;\text{km/s}$ ($M=1{,}018\,M_\odot$, $R=5846{,}6\;\text{km}$) | H0: $|\Delta v|\le 1{,}27\;\text{m/s}$ (nincs korlát); feltételes H1 $2\sigma$ ($\Delta\chi^2=4$) kizárja: $L\in[4{,}08\times 10^6, 3{,}37\times 10^7]\;\text{m}$ ($n=1{,}5$, $\max\Delta\chi^2=52{,}3$ / $72{,}0$) | `step8`, `step19`, `covariant_closure_50aa2db0.py` |
| **Gránitkamra ($8558\;\text{kg}$) IN–OUT teszt** | Feltételes H1 differenciális órajel $\Delta y_{\text{H1}}(\text{IN}-\text{OUT})$ | $L=1{,}00\;\text{m}$-nél: Top-hat $-1{,}436\times 10^{-24}$ ($\equiv 0$ ha $L\ge 2{,}944\;\text{m}$); Yukawa $-3{,}781\times 10^{-25}$ ($L=10\;\text{m}$-nél $-4{,}818\times 10^{-25}$) | `step8`, `step19`, `covariant_closure_50aa2db0.py` |
| **8 tekercses szolenoid ($N=2000$, $14\;\text{mA}$)** | Mágneses csúcsjel ($L_p=0{,}630957\;\text{m}$), $r=L/t$ kvadratúra és $L=0{,}95\;\text{m}$ | $\Delta W(L_p)=1{,}96114452\;\text{J/A}^2$, $k_g = -6{,}13117767\times 10^{-45}\;\text{A}^{-2}$, $\Delta y(L_p) = -1{,}2017108\times 10^{-48} \sim -10^{-48}$; A vs. B út ($L=3\;\text{m}$): $7{,}44\times 10^{-16}$; $L=0{,}95\;\text{m}$ azonosság: $8{,}88\times 10^{-16}\;\text{J/A}^2$ | `step15`, `step20`, `make_revision3_checks.py` |

---


### 3.1. V27 re-review javítások: 4-komponensű külső BVP, egzakt külső linearizáció, dinamikus skalárhullám és GL2 N-konvergencia

*(Forrás: `curved_yukawa_adm_and_stability.py`, `curved_yukawa_audit_summary.json`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_adm_sequences.csv`.)*

- **Külső lapse (5. hiba):** a külső $\nu$ negyedik ODE-komponensként együtt integrálódik $(\phi_L,\psi_L,m,\nu)$, `solve_bvp` `tol=1e-12`, végfeltétel $\nu(R_{\max})=\tfrac12\ln(1-2m(R_{\max})/R_{\max})$. SLy, $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $N_{\text{mult}}\in\{10,12,14,16\}$ szórások: $\delta k/k=6{,}17\times 10^{-12}$, $\delta M_{\text{ADM}}/M_{\text{ADM}}=5{,}54\times 10^{-13}$, $\delta\nu_R=1{,}63\times 10^{-14}$.

| $N_{\text{mult}}$ | $R_{\max}$ (km) | $k_{\text{curved}}$ (km$^{-1}$) | $M_{\text{ADM}}$ ($M_\odot$) | $\nu_R$ |
|---|---|---|---|---|
| `10` | `111.4657` | `0.233542329280` | `1.509998138706` | `-0.246306687228` |
| `12` | `131.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `14` | `151.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `16` | `171.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |

- **Külső linearizáció (2. hiba):** a `linearized_exterior_response` a $(\delta\phi_L(R_*),\delta m(R_*),\omega^2)$ szerinti egzakt elsőrendű variációs BVP-ből adja $\delta\phi_L'(R_*)$-t, $\delta m(R_*)$-t és $\delta M_{\text{ADM}}=\delta m(\infty)$-t; a sajátérték-feltétel $\delta M_{\text{ADM}}=0$.
- **Dinamikus skalárblokk (1. és 3. hiba):** a skalárblokkban szerepel a $-\omega^2 e^{-2\nu}\delta\phi_L$ tag, és a külső választ $\omega$-függő BVP adja (`dynamical_scalar_wave=True`); megoldó: GL2 darabonkénti szegmenses Volterra-diszkretizáció, $N=480$. A kvázisztatikus és a dinamikus gyök a közölt 6 tizedesjegyen egyezik.
- **Fordulópontok:** `minimize_scalar` + 7 pontos kvadratikus illesztés $\pm0{,}002$ dex-en; a zérushelyek 3 pontos kvadratikus illesztésből.

| Sorozat | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dinamikus $\omega_0^2=0$ | $\Delta$ dinamikus (dex) | Kvázisztatikus $\omega_0^2=0$ | Befagyasztott $\delta\phi_L=0$ | $\Delta$ befagyasztott (dex) |
|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316388` | **`15.316425`** | **`+2.04e-05`** | `15.316425` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.306886` | **`15.307172`** | **`+5.84e-06`** | `15.307172` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303330` | **`15.303080`** | **`-6.86e-06`** | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.292989` | **`15.293042`** | **`+2.07e-05`** | `15.293042` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227981` | **`15.227790`** | **`+5.71e-06`** | `15.227790` | `15.251385` | `+2.36e-02` |

**GL2 $N$-konvergencia** (SLy, $L=10\;\text{km}$, $\alpha_Y=1$, $\log_{10}\rho_{0c}=15{,}311404$): $p_{\text{eff}}=5{,}23$ (a három legfinomabb szintből), Richardson-érték $\omega_0^2=6387684{,}936\;\text{s}^{-2}$; $N=480$-nál a relatív hiba $\approx1{,}28\times 10^{-5}$. A sorozat monoton, előjelváltás nélkül konvergál (a korábbi FD-megoldó $N$-függő előjelváltása nem jelentkezik).

| $N$ | $\omega_0^2$ ($\text{s}^{-2}$) | relatív eltérés a Richardson-értéktől |
|---|---|---|
| `80` | `6450570.762819` | `9.84e-03` |
| `120` | `6403998.296759` | `2.55e-03` |
| `160` | `6393544.377742` | `9.17e-04` |
| `240` | `6388970.545464` | `2.01e-04` |
| `360` | `6387947.756886` | `4.11e-05` |
| `480` | `6387766.842088` | `1.28e-05` |
| `600` | `6387713.828704` | `4.52e-06` |
| `800` | `6387691.354875` | `1.00e-06` |
| `1200` | `6387685.705914` | `1.21e-07` |

**Értelmezési korlát:** az $\omega_0^2=0$ és az $M_{\text{ADM}}$-maximum egyezése ugyanazon egyenletrendszer belső konzisztenciavizsgálata; a $\sim10^{-5}$ dex eltérés a 3 pontos gyökillesztés és a 7 pontos fordulópont-illesztés felbontásán belül értendő, nem önálló fizikai teszt.

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
Minden alábbi parancs közvetlenül a kicsomagolt `/workspace/ta-v27-work/scalar-tensor-yukawa-fit` munkakönyvtárból futtatható hálózati kapcsolat nélkül:

1. **Görbült külső `yukawa_action_el` Einstein–skalár BVP, ADM-tömeg ($M_{\text{ADM}}$), tömeg–sugár görbék és csatolt $(\delta h, \delta\phi_L, \delta m)$ radiális stabilitás (`curved_yukawa_adm_and_stability.py`):**
   ```bash
   python3 curved_yukawa_adm_and_stability.py
   ```
   - *Várt kimenet (~2–3 perc):*  
     - Létrehozza a `curved_yukawa_adm_sequences.csv`, `curved_yukawa_coupled_stability.csv` és `curved_yukawa_audit_summary.json` fájlokat (valamint frissíti a `data/v23/` almappát).  
     - Ellenőrzi az $R_{\max}$-függetlenséget ($\Delta k_{\text{curved}}/k_{\text{curved}} = 6{,}17\times 10^{-12}$, $\Delta M_{\text{ADM}}/M_{\text{ADM}} = 5{,}54\times 10^{-13}$), a görbült vs. lapos Robin-meredekséget ($+24{,}74\%$ SLy $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $\alpha_Y=1$ esetén), az ADM-tömegmaximumokat ($M_{\text{ADM},\max}^{\text{SLy},10\text{km}} = 2{,}03321\,M_\odot$, $R_{1{,}4}^{\text{ADM}} = 11{,}5359\;\text{km}$) és a csatolt stabilitási zérushelyeket ($\Delta\log_{10}\rho_{0c} = 2{,}04\times 10^{-5}\;\text{dex}$, $595\times$-ös javulás).

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
