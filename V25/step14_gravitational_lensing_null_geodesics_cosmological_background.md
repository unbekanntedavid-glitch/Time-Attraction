# 14. lépés: Konform vs. diszform gravitációs lencsézés (null-geodetikus fényelhajlási szög, $\hat{\alpha}_{\text{def}}(b)$), NICER sugárkövetés és kozmológiai háttérbe ágyazás ($L_{\text{hor}}^{(m)} \approx 22{,}94\;\text{Gpc}$ vs. $L_{\text{hor}}^{(b)} \approx 57{,}98\;\text{Gpc}$)
# Step 14: Conformal vs. Disformal Gravitational Lensing (Null-Geodesic Light Deflection Angle $\hat{\alpha}_{\text{def}}(b)$), NICER Ray-Tracing, and Cosmological Background Embedding ($L_{\text{hor}}^{(m)} \approx 22.94\;\text{Gpc}$ vs. $L_{\text{hor}}^{(b)} \approx 57.98\;\text{Gpc}$)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step14_gravitational_lensing_null_geodesics_cosmological_background.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step13_curved_window_metric_variation_synge_world_function.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (2. fejezet, 199. sor), `hu_p2.tex` (4. fejezet, 31. sor; 5. fejezet, 52–55. sor), és `hu_p3.tex` (9. fejezet, 59–60. sor; 11. fejezet, 104. sor; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (3. szakasz, 161–186. sor; DOI: `10.5281/zenodo.23639084`)

---

## 1. Szó szerinti idézetek a V22–V23 kéziratokból

Az 5. lépésben (`step5_coupling_dictionary_ppn_screening.md`) levezettük a gyenge téri PPN $\gamma_{\text{PPN}}$ paramétert, a 9. lépésben (`step9_vhe_photon_propagation_hess_dispersion.md`) pedig az egyenes látóirány menti VHE gamma-futásiidő-diszperziót. Ebben a lépésben három, egymással szorosan összefüggő kérdést vezetünk le: **(A)** milyen eltérést okoz a tisztán időbeli diszform órametrika ($\tilde{g}_{00} = f_L^2 g_{00}$, $\tilde{g}_{ij} = g_{ij}$, **Branch II-Disf**) a **fényelhajlási szögben (gravitációs lencsézésben, $\hat{\alpha}_{\text{def}}(b)$)** a konform $\tilde{g}_{\mu\nu}^{\text{conf}} = f_L^2 g_{\mu\nu}$ (**Branch II-Conf**) esethez képest; **(B)** hogyan jelenik meg a **kozmológiai háttérbe ágyazás** ($\gamma_L\rho_{\text{bg}} = 1$) az exponenciális faktorban, külön-külön kiszámolva a barionos ($\Omega_b = 0{,}0493 \implies L_{\text{hor}}^{(b)} \approx 57{,}98\;\text{Gpc}$) és az összesített kozmikus anyagsűrűséghez ($\Omega_m \approx 0{,}315 \implies L_{\text{hor}}^{(m)} \approx 22{,}94\;\text{Gpc}$) tartozó horizontskálát; és **(C)** mit jelent ez a **NICER pulzár-sugárkövetési (ray-tracing) likelihoodokra** (PSR J0030+0451, PSR J0740+6620, PSR J0437-4715).

1. **Az időbeli törésmutató és a kozmológiai háttérbe ágyazás a V22 monográfiában (`hu_p1.tex` 199. sor és `hu_p2.tex` 54. sor):**
   - `hu_p1.tex`, 199. sor:
     > *„Az időáramlást úgy értelmezem, mint amely effektív helyi hordozósűrűséggel rendelkezik (analóg módon egy $n_t(\mathbf{x}) = f_L(\mathbf{x})^{-1} = \exp(+\gamma_L \bar{\rho}_L(\mathbf{x}))$ törésmutatójú optikai közeggel vagy egy összenyomható időfolyadékkal).”*
   - `hu_p2.tex`, 5. fejezet, 54. sor:
     > *„**Kozmológiai háttérbe ágyazás ($L \approx 22{,}94\,\text{Gpc}$):** Ahhoz, hogy a $\bar{\rho}_L$ megaparszekes skálán elérje a $\gamma_L\bar{\rho}_L \sim 1$ értéket, a $B(\mathbf{x},L)$ gömbnek magában kell foglalnia a környező intergalaktikus és kozmológiai $\rho_{\text{bg}} \approx \Omega_b\rho_{\text{crit}}$ háttérsűrűséget. Egyenletes kozmikus barionsűrűségnél a $\gamma_L\bar{\rho}_{\text{bg}} = 1$ feltétel a formális $L \approx 22{,}94\,\text{Gpc}$ ($7{,}08\times 10^{26}\,\text{m}$) horizontskálán áll elő.”*

2. **A NICER pulzár-likelihoodok a V22 monográfiában (`hu_p3.tex`, 59–60. és 104. sor):**
   - `hu_p3.tex`, 9. fejezet, 60. sor:
     > *„A három NICER-pulzárra (PSR J0030+0451, PSR J0740+6620, PSR J0437-4715) az oksági végponti ($M_{\text{max}}(\text{causal})$) tömegnormálással kiszámított együttes log-likelihood megmutatja, hogy $L = 1\,\text{km}$-nél a GR-hez viszonyított eltolódás legfeljebb **$4{,}8\times 10^{-4}\,\text{dex}$** ($\Delta\ln\mathcal{L} \approx 1{,}1\times 10^{-3}$). Ennek következtében **a jelenlegi NICER tömeg–sugár adatok egyik hipotézis alatt sem képesek statisztikailag megkülönböztetni az $L$ skálát a GR-től**.”*

---

## 2. Konform vs. tisztán időbeli diszform gravitációs lencsézés ($\hat{\alpha}_{\text{def}}(b)$)

Tekintsünk egy sztatikus, gyenge gravitációs teret, amelyben az Einstein-keretbeli gravitációs metrika ($g_{\mu\nu}$) vonaleleme izotróp koordinátákban:
$$ds^2 \;=\; g_{\mu\nu}dx^\mu dx^\nu \;=\; -\left(1 + \frac{2\Phi_N(\mathbf{x})}{c^2}\right)c^2 dt^2 \;+\; \left(1 - \frac{2\Psi_N(\mathbf{x})}{c^2}\right)d\mathbf{x}^2, \qquad \Phi_N(\mathbf{x}) = \Psi_N(\mathbf{x}) = -\frac{G M(<\!r)}{r}.$$

### 2.1. Fényelhajlás tisztán konform Jordan-csatolás esetén (Branch II-Conf: $\tilde{g}_{\mu\nu}^{\text{conf}} = f_L^2 g_{\mu\nu}$)
Ha a Maxwell-mező a konform módon átméretezett $\tilde{g}_{\mu\nu}^{\text{conf}} = f_L(\mathbf{x})^2 g_{\mu\nu} = e^{-2\phi_L(\mathbf{x})}g_{\mu\nu}$ metrikához csatolódik, akkor $D=4$ téridő-dimenzióban a $\sqrt{-\tilde{g}^{\text{conf}}}\tilde{g}_{\text{conf}}^{\mu\alpha}\tilde{g}_{\text{conf}}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} = \sqrt{-g}g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}$ konform-invariancia miatt a fénysugarak null-geodetikusai ($\tilde{g}_{\mu\nu}^{\text{conf}}dx^\mu dx^\nu = 0 \iff g_{\mu\nu}dx^\mu dx^\nu = 0$) **függetlenek a $\phi_L(\mathbf{x})$ skalármezőtől**. A Fermat-elv ($\delta\int n_{\text{conf}}(\mathbf{x})\,dl = 0$) szerinti effektív optikai törésmutató:
$$n_{\text{conf}}(\mathbf{x}) \;=\; \sqrt{\frac{g_{rr}}{-g_{00}}} \;\simeq\; 1 \;-\; \frac{\Phi_N(\mathbf{x}) + \Psi_N(\mathbf{x})}{c^2} \;=\; 1 \;-\; \frac{2\Phi_N(\mathbf{x})}{c^2}.$$
Következésképpen egy $b$ ütközési paraméterrel érkező fénysugár elhajlási szöge ($\hat{\alpha}_{\text{def}}(b) = -\int_{-\infty}^{+\infty}\nabla_\perp\ln n\,dz$) tisztán konform csatolásnál **kizárólag a barionos Weyl-potenciált ($\Phi_{\text{lens}}^{\text{conf}} = \frac{1}{2}(\Phi_N + \Psi_N) = \Phi_N$) érzi**:
$$\boxed{\;\hat{\alpha}_{\text{def}}^{\text{conf}}(b) \;=\; \frac{2}{c^2}\int_{-\infty}^{+\infty}\nabla_\perp\Phi_N\bigl(\sqrt{b^2+z^2}\bigr)\,dz \;=\; \frac{4 G M_{\text{bary}}(<\!b)}{c^2 b}.\;}$$
Vagyis tisztán konform csatolásban a skalármező ($\phi_L$) megnöveli ugyan a nemrelativisztikus csillagok és gázok körpályás sebességét ($V_{\text{rot}}^2(r) = V_{\text{bar}}^2(r) - c^2 r\phi_L'(r)$), **de nem járul hozzá a gravitációs lencsézéshez**!

### 2.2. Fényelhajlás a V23 tisztán időbeli diszform órametrikájában (Branch II-Disf: $\tilde{g}_{00} = f_L^2 g_{00}$, $\tilde{g}_{ij} = g_{ij}$)
Ezzel szemben, ha az elektromágneses hullámok is a V23-as diszform fizikai metrikán ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$) terjednek, akkor a null-kúp egyenlete ($d\tilde{s}_{\text{disf}}^2 = 0$):
$$0 \;=\; \tilde{g}_{00}^{\text{disf}}\,c^2 dt^2 \;+\; \tilde{g}_{ij}^{\text{disf}}\,dx^i dx^j \;=\; -f_L(\mathbf{x})^2\left(1 + \frac{2\Phi_N(\mathbf{x})}{c^2}\right)c^2 dt^2 \;+\; \left(1 - \frac{2\Psi_N(\mathbf{x})}{c^2}\right)d\mathbf{x}^2.$$
Ebből a Fermat-elvhez tartozó effektív optikai törésmutató (`hu_p1.tex`, 199. sor: $n_t(\mathbf{x}) = f_L(\mathbf{x})^{-1} = e^{+\phi_L(\mathbf{x})}$):
$$n_{\text{disf}}(\mathbf{x}) \;=\; \frac{c}{v_{\text{fázis}}(\mathbf{x})} \;=\; f_L(\mathbf{x})^{-1}\sqrt{\frac{1 - 2\Psi_N/c^2}{1 + 2\Phi_N/c^2}} \;\simeq\; 1 \;-\; \frac{\Phi_N(\mathbf{x}) + \Psi_N(\mathbf{x})}{c^2} \;+\; \phi_L(\mathbf{x}).$$
Behelyettesítve a $\nabla_\perp\ln n_{\text{disf}} \simeq -\frac{1}{c^2}\nabla_\perp(\Phi_N + \Psi_N - c^2\phi_L)$ transzverzális gradienst a Fermat-integrálba, a diszform lencsézési potenciál és fényelhajlási szög:
$$\boxed{\;\Phi_{\text{lens}}^{\text{disf}}(\mathbf{x}) \;=\; \Phi_N(\mathbf{x}) \;-\; \frac{1}{2}c^2\,\phi_L(\mathbf{x}), \qquad \hat{\alpha}_{\text{def}}^{\text{disf}}(b) \;=\; \frac{4 G M_{\text{bary}}(<\!b)}{c^2 b} \;-\; \int_{-\infty}^{+\infty}\frac{b}{r}\,\frac{d\phi_L(r)}{dr}\,dz.\;}$$
Mivel a vonzó időpotenciálra $\frac{d\phi_L}{dr} < 0$ (a sűrűség kifelé csökken), a második tag **pozitív (befelé irányuló) többlet-fényelhajlást** ad: a diszform időtörésmutató ($n_t = e^{+\phi_L}$) a galaxisok és galaxishalmazok körül nemcsak a forgási görbéket, hanem a gravitációs lencsézést is felerősíti!

---

## 3. Kozmológiai háttérbe ágyazás: egzakt aritmetika $\Omega_b = 0{,}0493$ ($57{,}98\;\text{Gpc}$) és $\Omega_m = 0{,}315$ ($22{,}94\;\text{Gpc}$) esetén

Bontsuk fel a kozmikus tömeg-energia sűrűséget egy térben homogén FLRW háttérsűrűségre ($\rho_{\text{bg}}(t) = \Omega_{\text{bg}}\,\rho_{\text{crit},0}$) és a lokális túlsűrűségre ($\delta\rho(\mathbf{x})$):
$$\rho(\mathbf{x}) \;=\; \rho_{\text{bg}} \;+\; \delta\rho(\mathbf{x}) \quad\implies\quad \bar{\rho}_L(\mathbf{x}) \;=\; \frac{3}{4\pi L^3}\int_{|\mathbf{y}-\mathbf{x}|\le L}\bigl[\rho_{\text{bg}} + \delta\rho(\mathbf{y})\bigr]d^3y \;=\; \rho_{\text{bg}} \;+\; \delta\bar{\rho}_L(\mathbf{x}).$$

### 3.1. Általános képlet és a két kozmológiai sűrűségparaméter ($\Omega_b$ vs. $\Omega_m$) pontos szétválasztása
A Planck 2018 (TT,TE,EE+lowE+lensing) paraméterekkel:
$$H_0 \;=\; 67{,}4\;\text{km\,s}^{-1}\text{Mpc}^{-1} \;=\; 2{,}184285\times 10^{-18}\;\text{s}^{-1}, \qquad R_H \;\equiv\; \frac{c}{H_0} \;=\; 1{,}37250\times 10^{26}\;\text{m} \;=\; 4{,}4480\;\text{Gpc},$$
$$\rho_{\text{crit},0} \;=\; \frac{3 H_0^2}{8\pi G} \;=\; 8{,}53287\times 10^{-27}\;\text{kg/m}^3.$$
Az egységnyi háttérkitevő ($\gamma_L\rho_{\text{bg}} = 1$, ahol $\gamma_L = \frac{G L_{\text{hor}}^2}{c^2}$) egzakt általános megoldása bármely $\Omega_{\text{bg}} \equiv \rho_{\text{bg}}/\rho_{\text{crit},0}$ sűrűséghányadra:
$$\gamma_L\,\rho_{\text{bg}} \;=\; \frac{G L_{\text{hor}}^2}{c^2}\,\Omega_{\text{bg}}\,\frac{3 H_0^2}{8\pi G} \;=\; \frac{3\,\Omega_{\text{bg}}}{8\pi}\left(\frac{L_{\text{hor}}}{c/H_0}\right)^2 \;=\; 1 \quad\Longleftrightarrow\quad \boxed{\;L_{\text{hor}}(\Omega_{\text{bg}}) \;=\; \sqrt{\frac{8\pi}{3\,\Omega_{\text{bg}}}}\,\frac{c}{H_0} \;=\; \frac{c}{\sqrt{G\,\rho_{\text{bg}}}}.\;}$$
Figyeljük meg, hogy a Newton-féle $G$ állandó pontosan kiesik a $\gamma_L\rho_{\text{bg}}$ szorzatból. A V22 `hu_p2.tex` 54. sorában idézett számértékek pontos tisztázásához elengedhetetlen szétválasztani a tisztán barionos ($\Omega_b$) és az összesített kozmikus anyag ($\Omega_m$) esetet:

| Háttérsűrűség-modell | Sűrűségparaméter $\Omega_{\text{bg}}$ | Háttérsűrűség $\rho_{\text{bg}} = \Omega_{\text{bg}}\rho_{\text{crit},0}$ ($\text{kg/m}^3$) | Dimenziótlan szorzó $\sqrt{\frac{8\pi}{3\Omega_{\text{bg}}}}$ | Horizontskála $L_{\text{hor}}$ ($\text{m}$) | Horizontskála $L_{\text{hor}}$ ($\text{Gpc}$) | Megjegyzés / Forrás |
|---|---|---|---|---|---|---|
| **Tisztán barionos kozmikus háttér ($\Omega_b$)** | $\Omega_b = 0{,}0493$ | $4{,}2067\times 10^{-28}$ | **$13{,}0306$** | **$1{,}7885\times 10^{27}$** | **$57{,}96\text{--}57{,}98\;\text{Gpc}$** | Egzakt $\gamma_L\rho_{b,0} = 1$ megoldás $\Omega_b = 0{,}0493$-ra |
| **Kozmikus anyagsűrűség-háttér ($\Omega_m$)** | $\Omega_m = 0{,}3147\text{--}0{,}3150$ | $2{,}6855\text{--}2{,}6879\times 10^{-27}$ | **$5{,}1583\text{--}5{,}1605$** | **$7{,}080\text{--}7{,}083\times 10^{26}$** | **$22{,}94\text{--}22{,}95\;\text{Gpc}$** | A V22 `hu_p2.tex` 54. sorában szereplő $7{,}08\times 10^{26}\;\text{m} \approx 22{,}94\;\text{Gpc}$ érték ($\Omega_m \approx 0{,}315$-nek felel meg!) |

*(Megjegyzés a V22 `hu_p2.tex` 54. sorához: a kézirat szövege szóban $\Omega_b\rho_{\text{crit}}$-et említ, de a közölt $7{,}08\times 10^{26}\;\text{m} = 22{,}94\;\text{Gpc}$ számérték aritmetikailag az $\Omega_m \approx 0{,}315$, $\rho_{\text{bg}} = \frac{1}{\gamma_L(L_{\text{hor}})} = 2{,}686\times 10^{-27}\;\text{kg/m}^3$ kozmikus anyagsűrűségből és a $\sqrt{8\pi/(3\Omega_m)} = 5{,}160$ szorzóból adódik, míg $\Omega_b = 0{,}0493$-ra $\sqrt{8\pi/(3\Omega_b)} = 13{,}031$ és $L_{\text{hor}}^{(b)} = 1{,}789\times 10^{27}\;\text{m} \approx 57{,}98\;\text{Gpc}$.)*

### 3.2. A homogén háttérfaktor leválasztása a lokális dinamikai és lencsézési gradiensekről
A Cauchy-féle multiplikatív tulajdonság (`hu_p1.tex`, Eq. 2: $f(\rho_1+\rho_2) = f(\rho_1)f(\rho_2)$) miatt az exponenciális időfaktor szorzatként bomlik fel:
$$f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L[\rho_{\text{bg}} + \delta\bar{\rho}_L(\mathbf{x})]\bigr) \;=\; \underbrace{e^{-\gamma_L\rho_{\text{bg}}}}_{f_{L,\text{bg}}\;(\text{globális állandó})}\,\cdot\,\underbrace{e^{-\gamma_L\,\delta\bar{\rho}_L(\mathbf{x})}}_{f_{L,\text{loc}}(\mathbf{x})}.$$
Mivel $f_{L,\text{bg}}$ térben állandó ($\nabla f_{L,\text{bg}} = 0$), a globális időkoordináta $d\tilde{t} = f_{L,\text{bg}}\,dt$ átméretezésével a homogén háttérfaktor beolvad a kozmológiai háttérórába, míg a lokális dinamikai gyorsulásban ($\mathbf{a}_{\text{eff}} = -\nabla\Phi_N + c^2\nabla\phi_L$) és a diszform fényelhajlásban ($\hat{\alpha}_{\text{def}}^{\text{disf}}$) szereplő logaritmikus gradiens **független a homogén $\rho_{\text{bg}}$ háttértől**:
$$\nabla\ln f_L(\mathbf{x}) \;=\; -\nabla\phi_L(\mathbf{x}) \;=\; -\gamma_L\,\nabla\delta\bar{\rho}_L(\mathbf{x}).$$

---

## 4. Következmény a NICER pulzár-sugárkövetésre (PSR J0030+0451, PSR J0740+6620, PSR J0437-4715)

Egy $R_*$ sugarú neutroncsillag felszínéről ($r \ge R_*$) induló röntgenfotonok sugárkövetésénél (*ray-tracing*, Beloborodov / Pechenick–Ftaclas–Cohen-formula) a külső vákuumban ($\rho(r > R_*) = 0$):
- **3D top-hat mag esetén ($L \le 1\;\text{km} \ll R_*$):** A csillag felszínén kívül ($r > R_* + L$) a top-hat gömb nem metsz bele a csillagba, így $\bar{\rho}_L(r > R_* + L) = 0 \implies f_L = 1$, a külső fényelhajlás egzaktul Schwarzschild-geodetikusokat követ, és $L = 1\;\text{km}$-nél a tömeg–sugár görbe relatív eltolódása mindössze $4{,}8\times 10^{-4}\;\text{dex}$ ($\Delta\ln\mathcal{L} \approx 1{,}1\times 10^{-3}$, V22 `hu_p3.tex` 60. sor).
- **Yukawa-mag esetén (`yukawa_action_el`):** A külső skalárprofil a görbült külső téridőben $r \to \infty$-re exponenciálisan lecseng ($\phi_L(r) \propto r^{-(1+M_{\text{ADM}}/L)}e^{-r/L}$), a csillag aszimptotikus gravitációs tömege pedig a 7. és 22. lépésben (`step7_neutron_star_tov_stability.md`, `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`) levezetett $M_{\text{ADM}} = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$ ADM-tömegbe olvad be.

---

## 5. Ellenőrzések (Checks)

1. **Vákuum- és GR-határeset ($\phi_L \to 0 \implies f_L \to 1$):**
   - Ha $\phi_L(\mathbf{x}) \to 0$, akkor $\frac{d\phi_L}{dr} \to 0$, így mind a konform ($\hat{\alpha}_{\text{def}}^{\text{conf}}$), mind a diszform ($\hat{\alpha}_{\text{def}}^{\text{disf}}$) fényelhajlási szög visszaadja az Einstein-féle $\hat{\alpha}_{\text{GR}}(b) = \frac{4GM}{c^2 b}$ értéket $\checkmark$.
2. **A $L_{\text{hor}}^{(m)} \approx 22{,}94\;\text{Gpc}$ és $L_{\text{hor}}^{(b)} \approx 57{,}98\;\text{Gpc}$ számszerű ellenőrzése:**
   - $\Omega_m = 0{,}3147 \implies \rho_{m,0} = 2{,}6855\times 10^{-27}\;\text{kg/m}^3 \implies \sqrt{8\pi/(3\Omega_m)} = 5{,}1605 \implies L_{\text{hor}}^{(m)} = 7{,}083\times 10^{26}\;\text{m} = 22{,}95\;\text{Gpc}$ (és $\rho_{\text{bg}} = 2{,}686\times 10^{-27}\;\text{kg/m}^3$-ra $7{,}08\times 10^{26}\;\text{m} = 22{,}94\;\text{Gpc}$) $\checkmark$.
   - $\Omega_b = 0{,}0493 \implies \rho_{b,0} = 4{,}2067\times 10^{-28}\;\text{kg/m}^3 \implies \sqrt{8\pi/(3\Omega_b)} = 13{,}0306 \implies L_{\text{hor}}^{(b)} = 1{,}7885\times 10^{27}\;\text{m} = 57{,}96\text{--}57{,}98\;\text{Gpc}$ $\checkmark$.
