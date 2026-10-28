# 17. lépés: A V23 (8) faktorizált Chandrasekhar–hidrosztatikai különbség-azonosság (`eq:factored_diff`) és a `gradient` módbeli `metric` tehetetlenségi súly ($\mathcal{W}_{\text{metric}}(r) \le 0$, $L \ge 367\;\text{m}$) egzakt analitikus levezetése
# Step 17: Exact Analytical Derivation of the V23 Eq. (8) Factorized Chandrasekhar–Hydrostatic Difference Identity (`eq:factored_diff`) and the `gradient`-Mode `metric` Inertia-Weight Indefiniteness ($\mathcal{W}_{\text{metric}}(r) \le 0$, $L \ge 367\;\text{m}$)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (2. szakasz, 121–126. sor, Eq. 8 és 2. táblázat lábjegyzete, 152. sor; DOI: `10.5281/zenodo.23639084`)  
- V23 ellenőrző kódok: `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/covariant_closure_50aa2db0.py` (100–120. sor) és `radial_modes.py` (20–39. sor)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (8. fejezet 2. pont, 409–419. sor; DOI: `10.5281/zenodo.23604423`)

---

## 1. Szó szerinti idézetek a V22 és V23 kéziratokból

A 16. lépésben (`step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md`) levezettük a V23 (4)–(5) egyenletekben szereplő háttérbeli $T^\mu_{\;\nu}(X)$ feszültségtenzort. A V23 kézirat 2. szakasza (121–126. sor, Eq. 8 és 152. sor) és a V22 monográfia 8. fejezete (409–418. sor) megmutatja, hogyan okozza ugyanez a háttér-módosítás a két befagyasztott-$f_L$ ($\delta f_L = 0$) Sturm–Liouville-operátor (`matter` és `metric`) eltérését a $dM/d\rho_c = 0$ statikus fordulóponttól:

1. **A V23 (8) faktorizált Chandrasekhar–hidrosztatikai különbség-azonosság (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 2. szakasz, 121–125. sor, Eq. 8):**
   > *„Annak megértéséhez, hogy miért tér el mindkét befagyasztott $f_L$-ű operátor a $dM/d\rho_c = 0$ statikus tömegmaximumtól, először `SymPy`-ban meghatároztam a befagyasztott $f_L$-ű Chandrasekhar-féle $\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi]$ operátor és a linearizált $\mathcal{H}_{\text{lin}}[\xi] = (\delta P)' + (\delta\varepsilon+\delta P)\nu' + (\varepsilon+P)\delta\nu'\big|_{\delta f_L=0}$ hidrosztatikai operátor egzakt különbségét tömeg-módban:”*
   $$\begin{equation}
     \mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} = \frac{4\pi r^{5/2}(f_L - 1)\Gamma P \varepsilon e^{2\nu}}{(r - 2m)^{5/2}}\Big[(4\pi r^3 P - 3r + 7m)\xi - r(r-2m)\xi'\Big]. \tag{V23 Eq. 8}
   \end{equation}$$
   > *„A (8) kifejezés $f_L = 1$ (GR) esetén azonosságszerűen eltűnik, mivel a Chandrasekhar-féle parciális integrálás a GR-beli $\lambda' + 2\nu' = 8\pi r e^\lambda(\varepsilon+P)$ Einstein-azonosságot helyettesíti be. Ha azonban $f_L \ne 1$, akkor $\lambda' + 2\nu' = 8\pi r e^\lambda(f_L\varepsilon + P)$, ami előállítja a (8)-beli $(f_L-1)$ maradéktagot, RÁADÁSUL elhagyja a $\delta f_L[\delta\varepsilon,\xi]$ nemlokális perturbációt.”*

2. **A `metric` operátor $\mathcal{W}(r) \le 0$ tehetetlenségi súly-indefinitsége $L \ge 367\;\text{m}$ felett (`tal_book_HU_2026-10-24.tex`, 417. sor és V23 2. táblázat lábjegyzete, 152. sor):**
   - `tal_book_HU_2026-10-24.tex`, 417. sor:
     > *„**Közelítő `metric` operátor:** $L \ge 367\,\text{m}$ felett 12 sorozaton (609 modellen) a $\mathcal{W}(r)$ tehetetlenségi súly indefinitté válik ($\mathcal{W}\le 0$), és a kiválasztott pozitív normájú ágon negatív $\omega_0^2$ eredményeket ad; a fizikai stabilitás megállapítása nyitott (`50aa2db0` feladat).”*
   - `tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 152. sor:
     > *„$^*$Gradiens-módban $L=10\,\text{km}$-nél a közelítő befagyasztott $f_L$-ű `metric` operátor súlya a felszín közelében indefinitté válik ($\mathcal{W}(r)\le 0$); a (9) statikus érintőmodell pozitív Gram-súlyt használ (ez nem a kinetikus pozitivitás bizonyítéka).”*

---

## 2. A V23 (8) faktorizált különbség-azonosság (`eq:factored_diff`) lépésről lépésre történő levezetése

### 2.1. A két operátor definíciója a V23 `covariant_closure_50aa2db0.py` (101–120. sor) konvenciójában
A `covariant_closure_50aa2db0.py` (103–118. sor) és a `radial_modes.py` (3–7. sor) jelöléseiben ($G = c = 1$, $g_{tt} = -e^{\nu_C}$ ahol $\nu_C \equiv 2\nu$, és $e^\lambda \equiv \bigl(1 - \frac{2m}{r}\bigr)^{-1} = \frac{r}{r - 2m}$, így $e^{\lambda/2} = \sqrt{\frac{r}{r-2m}}$):
- A `mass` módbeli háttéregyenletek (`sub_mass_full`):
  $$m'(r) \;=\; 4\pi r^2 f_L(r)\,\varepsilon(r), \qquad \nu'(r) \;=\; A(r) \;\equiv\; \frac{m + 4\pi r^3 P}{r(r - 2m)}, \qquad P'(r) \;=\; -(\varepsilon + P)A(r), \qquad \varepsilon'(r) \;=\; -\frac{(\varepsilon + P)^2}{\Gamma P}A(r).$$
- A **befagyasztott-$f_L$ Chandrasekhar-féle `matter` Sturm–Liouville-operátor** ($\zeta(r) \equiv r^3 e^{-\nu(r)}\xi(r)$, ahol $\xi(r) = \Delta r/r$):
  $$\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] \;\equiv\; -\frac{d}{dr}\!\left(\Pi(r)\frac{d\zeta}{dr}\right) \;-\; Q_{\text{matter}}(r)\,\zeta(r),$$
  ahol az együtthatók (`covariant_closure_50aa2db0.py`, 106–108. sor):
  $$\Pi(r) \;=\; \sqrt{\frac{r}{r-2m}}\,e^{3\nu} r^{-2}\,\Gamma P, \qquad Q_{\text{matter}}(r) \;=\; \sqrt{\frac{r}{r-2m}}\,e^{3\nu} r^{-2}\,(\varepsilon + P)\left(\nu'^2 + \frac{4\nu'}{r} - \frac{8\pi r\,P}{r - 2m}\right).$$
- A **linearizált hidrosztatikai operátor** $\delta f_L = 0$ mellett (`covariant_closure_50aa2db0.py`, 111–117. sor):
  $$\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} \;\equiv\; \frac{d(\delta P)}{dr} \;+\; (\delta\varepsilon + \delta P)\,A(r) \;+\; (\varepsilon + P)\,\delta A(r)\Big|_{\delta f_L=0},$$
  ahol az Euler-féle nyomás-, energiasűrűség- és tömegperturbáció:
  $$\Delta P \;=\; -\Gamma P\,e^\nu r^{-2}\zeta'(r), \qquad \delta P \;=\; \Delta P - r\xi P'(r), \qquad \delta\varepsilon \;=\; \frac{\varepsilon + P}{\Gamma P}\delta P, \qquad \delta m\Big|_{\delta f_L=0} \;=\; -4\pi r^3\bigl(f_L\varepsilon + P\bigr)\xi,$$
  a geometria-faktor linearizált megváltozása pedig:
  $$\delta A\Big|_{\delta f_L=0} \;=\; \frac{\partial A}{\partial m}\,\delta m\Big|_{\delta f_L=0} \;+\; \frac{\partial A}{\partial P}\,\delta P \;=\; \frac{r + 8\pi r^3 P}{r(r-2m)^2}\,\delta m\Big|_{\delta f_L=0} \;+\; \frac{4\pi r^2}{r - 2m}\,\delta P.$$

### 2.2. Miért keletkezik az $(f_L - 1)$-gyel arányos maradék a parciális integrálásnál?
Fejtsük ki a $\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}e^{2\nu}\mathcal{H}_{\text{lin}}[\xi]\big|_{\delta f_L=0}$ különbséget!
1. Mivel $\delta P = -e^{-\nu}\sqrt{\frac{r-2m}{r}}\,\Pi(r)\zeta'(r) + r(\varepsilon + P)A\xi$, a $-\frac{d}{dr}(\Pi\zeta')$ vezető másodrendű derivált tag ($\propto \xi''(r)$) **azonosan kiesik** a $-\sqrt{\frac{r}{r-2m}}e^{2\nu}\frac{d(\delta P)}{dr}$ vezető tagjával!
2. A fennmaradó elsőrendű ($\propto \xi'(r)$) és nulladrendű ($\propto \xi(r)$) tagok összevonásakor a kulcslépés a $\sqrt{\frac{r}{r-2m}}e^{2\nu}$ metrikus szorzó radiális deriváltja:
   $$\frac{d}{dr}\!\left(\ln\left[\sqrt{\frac{r}{r-2m}}\,e^{2\nu}\right]\right) \;=\; \frac{1}{2}\lambda'(r) \;+\; 2\nu'(r) \;=\; \frac{r m'(r) - m(r)}{r(r-2m)} \;+\; \frac{2(m + 4\pi r^3 P)}{r(r-2m)} \;=\; \frac{m + 4\pi r^3\bigl[f_L(r)\varepsilon + 2P\bigr]}{r(r-2m)}.$$
   Általános relativitáselméletben ($f_L = 1$) a számlálóban $f_L\varepsilon = \varepsilon$ áll, és ugyanígy a $\delta m_{\text{GR}} = -4\pi r^3(\varepsilon + P)\xi$ tömegperturbációban is $(\varepsilon + P)$ szerepel, így a két kifejezés különbsége azonosan eltűnik.
3. Ha azonban $f_L(r) \neq 1$, akkor mind a háttérbeli $m'(r) = 4\pi r^2 f_L\varepsilon$ deriváltból, mind a $\delta m\big|_{\delta f_L=0} = -4\pi r^3(f_L\varepsilon + P)\xi$ perturbációból egy $(f_L - 1)\varepsilon$-nal arányos eltolódás marad vissza:
   - A $\xi'(r)$-es tag együtthatója (amely a $\Delta P = -\Gamma P(3\xi - r\nu'\xi + r\xi')$ kifejezés $r\xi'$ tagjának és a $\frac{1}{2}\lambda'\big|_{f_L} - \frac{1}{2}\lambda'\big|_{f_L=1} = \frac{4\pi r^2(f_L - 1)\varepsilon}{r - 2m}$ szorzónak a szorzatából adódik):
     $$\text{Coeff}(\xi') \;=\; -\sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\Gamma P\,r\cdot\frac{4\pi r^2(f_L - 1)\varepsilon}{r - 2m}\cdot\frac{r}{r-2m} \;=\; -\frac{4\pi r^{7/2}(f_L - 1)\Gamma P\varepsilon e^{2\nu}}{(r - 2m)^{3/2}}.$$
   - A $\xi(r)$-es tag együtthatója a $(3 - r\nu')\Delta P$ tényezőből és a $\frac{\partial A}{\partial m}\bigl(\delta m|_{f_L} - \delta m|_{f_L=1}\bigr)$ tagból összegezve pontosan $\frac{4\pi r^{5/2}(f_L - 1)\Gamma P\varepsilon e^{2\nu}}{(r-2m)^{5/2}}(4\pi r^3 P - 3r + 7m)$.
   Kiemelve a közös $\frac{4\pi r^{5/2}(f_L - 1)\Gamma P\varepsilon e^{2\nu}}{(r-2m)^{5/2}}$ faktort, egzaktul megkapjuk a **V23 (8) egyenletet**:
   $$\boxed{\;\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} \;=\; \frac{4\pi r^{5/2}(f_L - 1)\Gamma P \varepsilon e^{2\nu}}{(r - 2m)^{5/2}}\Big[(4\pi r^3 P - 3r + 7m)\xi - r(r-2m)\xi'\Big].\;}$$

---

## 3. Miért válik a `metric` operátor $\mathcal{W}_{\text{metric}}(r)$ súlya negatívvá a felszín közelében ($L \ge 367\;\text{m}$, `gradient` mód)?

A `radial_modes.py` (11–12. és 30–38. sor) definíciója szerint a `metric` operátorváltozat a Sturm–Liouville-féle $Q(r)$ és $\mathcal{W}(r)$ együtthatókban az anyagi $(\varepsilon, P)$ mennyiségeket a módosított háttér effektív Einstein-tenzor forrásaival ($\varepsilon_{\text{eff}} = \varepsilon + \varepsilon_X$, $P_{r,\text{eff}} = P + P_{r,X}$) helyettesíti:
$$\mathcal{W}_{\text{metric}}(r) \;=\; e^{(3\lambda - \nu_C)/2}\,r^4\,\bigl[\varepsilon_{\text{eff}}(r) \;+\; P_{r,\text{eff}}(r)\bigr].$$

1. **A `mass` módban (`radial_modes.py`, 31–32. sor):**
   A 16. lépésben levezetett V23 Eq. (4) szerint $\varepsilon_{\text{eff}} = \varepsilon + \varepsilon_X = f_L(r)\varepsilon(r) > 0$ és $P_{r,\text{eff}} = P(r) > 0$, így:
   $$\varepsilon_{\text{eff}}^{\text{mass}}(r) + P_{r,\text{eff}}^{\text{mass}}(r) \;=\; f_L(r)\,\varepsilon(r) \;+\; P(r) \;>\; 0 \qquad \forall\,r \in [0, R_*),$$
   tehát `mass` módban a tehetetlenségi súly mindenhol pozitív marad.

2. **A `gradient` módban (`radial_modes.py`, 33–34. sor):**
   A 16. lépésben levezetett V23 Eq. (5) szerint $\varepsilon_{\text{eff}} = \varepsilon(r)$, míg $P_{r,\text{eff}}(r) = P(r) + P_{r,X}(r) = f_L(r)P(r) - \bigl(1 - f_L(r)\bigr)\frac{m(r)}{4\pi r^3}$. Ezért az effektív inerciális összeg:
   $$\boxed{\;\varepsilon_{\text{eff}}^{\text{grad}}(r) + P_{r,\text{eff}}^{\text{grad}}(r) \;=\; \varepsilon(r) \;+\; f_L(r)\,P(r) \;-\; \bigl(1 - f_L(r)\bigr)\,\frac{m(r)}{4\pi r^3}.\;}$$
   Vizsgáljuk meg ennek a kifejezésnek a viselkedését a csillag felszíne ($r \to R_*$) közelében!
   - A felszínen az anyagi nyomás és energiasűrűség eltűnik: $\lim_{r\to R_*} P(r) = 0$ és $\lim_{r\to R_*} \varepsilon(r) = 0$.
   - Ezzel szemben az átlagos bezárt sűrűség ($\frac{m(R_*)}{4\pi R_*^3} = \frac{M}{4\pi R_*^3} > 0$) szigorúan pozitív, és bármely $L > 0$ nemlokális simítási sugárra a felszíni gömb belenyúlik a csillag belsejébe, így $\bar{\rho}_L(R_*) > 0 \implies 1 - f_L(R_*) > 0$!
   - Következésképpen a felszín közelében ($r \to R_*$) a kifejezés határértéke **szigorúan negatív**:
     $$\boxed{\;\lim_{r\to R_*}\Bigl[\varepsilon_{\text{eff}}^{\text{grad}}(r) + P_{r,\text{eff}}^{\text{grad}}(r)\Bigr] \;=\; -\bigl(1 - f_L(R_*)\bigr)\,\frac{M}{4\pi R_*^3} \;<\; 0 \qquad (\forall\,L > 0)!\;}$$
   Amikor $L \ge 367\;\text{m}$, az $1 - f_L(R_*)$ faktor és a kéregbeli sűrűség-lefutás elég széles ahhoz, hogy az $N = 1200$ elemű végeselemes kvadratúra-rács (`radial_modes.py`, 45–54. sor) legkülső Gauss–Legendre pontjai beleessenek az $\varepsilon(r) + f_L(r)P(r) < (1 - f_L(r))\frac{m(r)}{4\pi r^3}$ kéregzónába, ahol $\mathcal{W}_{\text{metric}}(r) \le 0$ (12 sorozaton, 609 modellen az `item2_radial_grid.csv`-ben). Ez bizonyítja analitikusan, hogy a `metric` súly-indefinitség nem numerikus hiba, hanem annak a következménye, hogy a `metric` ansatz a statikus háttér-feszültséget ($P_{r,X} < 0$) tehetetlenségi tömegként helyettesítette be, míg a 12. lépésben (`step12_second_variation_coupled_pulsation_hessian.md`) a kovariáns hatás $\delta^2 S$ második variációjából levezetett valódi kinetikus súlyok ($\mathcal{W}_\xi \propto f_L(\varepsilon+P) > 0$ és $\mathcal{W}_\phi \propto \alpha_Y^{-1} > 0$) minden $r \in (0, R_*)$-en és minden $L > 0$-ra pozitív definitek!

---

## 4. Ellenőrzések (Checks)

1. **$f_L \to 1$ (GR) határérték a V23 (8) egyenletben:**
   - Ha $f_L(r) \equiv 1$, akkor az $(f_L - 1)$ szorzó miatt a V23 (8) jobb oldala azonosan nulla (`diff_gr_limit = 0` a `covariant_closure_50aa2db0.py` 120. sorában), így $\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] \equiv \sqrt{\frac{r}{r-2m}}e^{2\nu}\mathcal{H}_{\text{lin}}[\xi]$, és az $\omega_0^2 = 0$ sajátérték pontosan egybeesik a $dM/d\rho_c = 0$ statikus fordulóponttal $\checkmark$.
2. **A `gradient` módbeli $\mathcal{W}_{\text{metric}}(r)$ előjelváltási sugara ($r_{\text{sign}}$):**
   - Az előjelváltás pontosan ott következik be a kéregben, ahol $\varepsilon(r_{\text{sign}}) + f_L(r_{\text{sign}})P(r_{\text{sign}}) = \bigl(1 - f_L(r_{\text{sign}})\bigr)\frac{m(r_{\text{sign}})}{4\pi r_{\text{sign}}^3}$; amikor $L \to 0$, $1 - f_L(R_*) \to 0$, így $r_{\text{sign}} \to R_*$ (kiszorul a kvadratúra-pontok közül), míg $L \ge 367\;\text{m}$ esetén $r_{\text{sign}} < r_{\text{quad,max}} < R_*$ $\checkmark$.
