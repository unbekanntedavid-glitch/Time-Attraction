# 6. lépés: Galaktikus forgási görbék (155 / 175 SPARC-galaxis), a pozitív definit 3D top-hat / Yukawa mag-pozitivitás és az NGC 3198 ($R = 30{,}08\text{ kpc}$) negatív kontrollvizsgálata
# Step 6: Galactic Rotation Curves (155 / 175 SPARC Galaxies), Positive-Definite 3D Top-Hat / Yukawa Kernel Positivity, and the NGC 3198 ($R = 30.08\text{ kpc}$) Negative Control Check

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step6_sparc_galactic_ngc3198_consistency.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step5_coupling_dictionary_ppn_screening.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (3. fejezet, 256–265. sor, Eq. 6) és `hu_p2.tex` (4–5. fejezet, 1–61. sor, Eq. 7–8; DOI: `10.5281/zenodo.23604423`)  
- V17 és V19 Zenodo rekordok: DOI `10.5281/zenodo.19666220` (V17 SPARC 155 / 123 vakminta) és DOI `10.5281/zenodo.23497128` (V19 3D gömbi top-hat, 175 SPARC-galaxis)  
- Helyi SPARC és NGC 3198 munkadokumentumok: `docs/README_draft.md` (6–24. sor), `negative_density.md` (1–12. és 55–63. sor), `foundation_mode.md` (13–15. és 50–63. sor)

---

## 1. Szó szerinti idézetek a V22 monográfiából és a helyi SPARC-dokumentumokból

Az 1–5. lépésben felépítettük a V23-as kovariáns skalár–tenzor / Yukawa keretrendszert és ellenőriztük a gyenge téri / Naprendszerbeli (PPN) határértéket. Ebben a lépésben vizsgáljuk meg a **galaktikus skálán a 3D mag matematikai pozitivitását (175 SPARC-galaxison a negatív sűrűségű pontok kizárását)** és az **NGC 3198 ($R = 30{,}08\text{ kpc}$) negatív kontrollpontot**, szigorúan szétválasztva a megfigyelés által megkövetelt erősítést ($C_{\text{req}}$) és az izolált, rögzített $\alpha_{\text{std}} = 1/(4\pi)$ csatolású modell által ténylegesen jósolt erősítést ($C_{\text{pred}}^{\text{isol}}$).

1. **A 3D gömbi top-hat definíciója és a negatív sűrűség ($\rho < 0$) kizárása (`hu_p2.tex`, 12–21. sor, Eq. 7–8):**
   > *„Pontonkénti deriváltak vagy lineáris közelítések helyett az időbeli csillapítás forrását tetszőleges $\mathbf{x} \in \mathbb{R}^3$ pontban a fizikai $\rho(\mathbf{y}) = T_{00}(\mathbf{y})/c^2 \ge 0$ tömeg-energia sűrűség $L > 0$ sugarú $B(\mathbf{x},L)$ gömbön belüli, súlyozatlan háromdimenziós top-hat átlagaként definiálom:”*
   $$\begin{equation}
     \bar{\rho}_L(\mathbf{x}) \;\equiv\; \frac{3}{4\pi L^3}\int_{|\mathbf{y}-\mathbf{x}|\le L} \rho(\mathbf{y})\,d^3y \;=\; \frac{M_{\text{enc}}(\mathbf{x};L)}{\frac{4}{3}\pi L^3}, \qquad
     f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(\mathbf{x})\bigr). \tag{V22 Eq. 7}
   \end{equation}$$
   > *„Ennek következtében az exponenciális csillapítási faktor kielégíti az egzakt analitikus korlátokat: $0 < f_L(\mathbf{x}) \le 1 \quad \forall\,\mathbf{x}\in\mathbb{R}^3,\;\forall\,L > 0$. Sem negatív sűrűség, sem előjelváltás, sem az $f_L$-ből eredő metrikus horizont-átfordulás nem fordulhat elő.”*

2. **Az SPARC 175 / 155 / 123 galaxisos minta és a 3D tengelyszimmetrikus korongmodell (`hu_p2.tex`, 47–55. sor):**
   > *„A SPARC főkatalógusa **175 késői típusú és törpegalaxist** tartalmaz $3{,}6\,\mu\text{m}$-es felületi fotometriával és felbontott $\text{H\,\textsc{i}}/\text{H}\alpha$ forgási görbékkel. A V17-es feldolgozásom során 20 súlyos kinematikai aszimmetriát vagy inklinációs hibajelzést ($i < 30^\circ$ vagy $Q=3$) mutató rendszer kiszűrése után állt elő a **155 galaxisból álló feldolgozott forgásigörbe-minta**. A szigorú adathigiéniai vágások (legalább 5 radiális mérési pont és hiánytalan $V_{\text{gas}}, V_{\text{disk}}, V_{\text{bulge}}$ barionos dekompozíció) kijelölték a **123 galaxisos vakteszt-mintát**, amely **2596 radiális sebességmérést** foglal magában.”*
   > *„A 19. verzióban a vékonykorong fősíkbeli $\rho_{\text{mid}}(r) = \Sigma_{\text{bary}}(r)/(2h_z)$ közelítést felváltottam az egzakt 3D tengelyszimmetrikus $\rho_{\text{bary}}(R,z) = \frac{\Sigma_{\text{bary}}(R)}{2h_z}\exp(-|z|/h_z)$ tömegeloszlással ($h_z = 0{,}196\,R_d^{0{,}633}$ skálamagasság mellett) [...]. Mind a 175 SPARC-galaxison és az $L \in [10^{16}, 10^{23}]\,\text{m}$ tartományban a maximális izolált galaktikus csillapítás $\max |1 - f_L| = 2{,}284\times 10^{-7}$ a logaritmikus rácson (és **$2{,}300\times 10^{-7}$** folytonos optimalizálással $L \approx 23{,}09\,\text{kpc}$-nél a nagy tömegű **UGC02487** spirálgalaxisban).”*

3. **Az NGC 3198 ($R = 30{,}08\text{ kpc}$) sebességértékeinek pontos tisztázása (`docs/README_draft.md`, 6–19. sor és `foundation_mode.md`, 13–14. sor):**
   > *„At $R = 30.08\text{ kpc}$ (`NGC3198_rotmod.dat`, $D = 13.8\text{ Mpc}$):*
   > - **`Vobs` = `146 ± 2 km/s`**: *Observed rotation speed. '145' in earlier text is consistent with this within $1\sigma$. Neighbouring points: 148 at 28.16 kpc, 147 at 32.14 kpc.*
   > - **`Vbar` ($\Upsilon_{\text{disk}} = 0.5, \Upsilon_{\text{bul}} = 0.7$, signed gas) = `70.32 km/s`**: *Baryonic speed with the documented v16 mass-to-light ratio.*
   > - **`Vbar` ($\Upsilon_{\text{disk}} = 1.0$) = `88.53 km/s`**: *Baryonic speed computed from the raw SPARC columns $\sqrt{V_{\text{gas}}^2 + V_{\text{disk}}^2}$. **'88.5' is a baryonic value with $\Upsilon_{\text{disk}} = 1$, not an observed speed.***
   > *The required boost is $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$, which is $4.311$ for $\Upsilon = 0.5$ and $2.720$ for $\Upsilon = 1.0$.”*

---

## 2. A galaktikus körpályasebesség levezetése és a 3D mag pozitivitási ellenőrzése

### 2.1. Körpályasebesség a módosított Euler-egyenletből ($h' = -\nu' + \phi_L'$)
A 4. és 5. lépésben (`step4_covariant_conservation_hydrostatic.md`, `step5_coupling_dictionary_ppn_screening.md`) levezettük, hogy a kovariáns megmaradási törvényből ($\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$) egy galaktikus korong fősíkjában ($z = 0$) keringő barionos elemre ható radiális effektív gyorsulás:
$$a_{r,\text{eff}}(r) \;=\; -\frac{d\Phi_N(r)}{dr} \;+\; c^2\,\frac{d\phi_L(r)}{dr} \;=\; -\frac{V_{\text{bar}}^2(r)}{r} \;-\; c^2\,\frac{d\ln f_L(r)}{dr}.$$
Az egyensúlyi körpályasebesség ($V_{\text{rot}}^2(r) = -r\,a_{r,\text{eff}}(r)$) így két módon írható fel:
1. **Additív skalár-gradiens (V19–V23 kanonikus skalár–tenzor / Yukawa) alak:**
   $$\boxed{\;V_{\text{rot}}^2(r) \;=\; V_{\text{bar}}^2(r) \;-\; c^2\,r\,\frac{d\phi_L(r)}{dr} \;=\; V_{\text{bar}}^2(r) \;+\; c^2\,r\,\gamma_L\,\left|\frac{d\bar{\rho}_L(r)}{dr}\right|\;}$$
   (mivel $\frac{d\bar{\rho}_L}{dr} < 0 \implies \frac{d\phi_L}{dr} < 0$, a skalárgradiens befelé mutató többletvonzást ad).
2. **Multiplikatív $C(\rho)$ erőszorzó (V14–V16 fenomenológiai) alak (`hu_p1.tex`, Eq. 6):**
   $$V_{\text{rot}}^2(r) \;=\; C(\rho(r))\,V_{\text{bar}}^2(r), \qquad C(\rho) \;=\; 1 + C_0\exp\!\left(-\frac{\rho}{\rho_c}\right).$$
   A két felírás közötti pontos algebrai kapcsolat bármely $r$ sugáron:
   $$C_{\text{eff}}(r) \;\equiv\; \frac{V_{\text{rot}}^2(r)}{V_{\text{bar}}^2(r)} \;=\; 1 \;-\; \frac{c^2\,r}{V_{\text{bar}}^2(r)}\,\frac{d\phi_L(r)}{dr}.$$

### 2.2. Mit bizonyít a „0 negatív sűrűségű pont 175 SPARC-galaxison” eredmény, és mit nem?
A bírálati egyértelműség érdekében szigorúan rögzítjük a 175 SPARC-galaxison végzett V19 (`sparc_tophat.py`) teszt pontos jelentését:
- **A korábbi gömbi derivált-képlet hibája (`negative_density.md`, 11. sor; `foundation_mode.md`, 16–18. és 55–63. sor):**
  Ha a sűrűséget a forgási görbe numerikus deriváltjából ($\rho_{\text{deriv}}(r) = \frac{1}{4\pi G r^2}\frac{d(r V_{\text{bar}}^2)}{dr}$) becsüljük, a nem-monoton gáz- és korongprofilok miatt a 140 nem-holdout SPARC-galaxisban **239** ($\Upsilon_{\text{disk}}=0{,}5$) illetve **303** ($\Upsilon_{\text{disk}}=1{,}0$) belső ponton $\rho_{\text{deriv}}(r) < 0$ lép fel (ebből 22 pont esik a $29{,}5\text{--}31{,}5\;\text{kpc}$ sávba, a $30{,}37\;\text{kpc}$-hez legközelebb az **NGC 4013** $R = 30{,}31\;\text{kpc}$-nél $\rho = -1{,}56\times 10^{-24}\;\text{kg/m}^3$-mal).
- **A V19 3D top-hat és a V23 Yukawa matematikai pozitivitása:**
  Mind a 3D gömbi top-hat mag ($K_L(\mathbf{x},\mathbf{y}) = \frac{\Theta(L - |\mathbf{x}-\mathbf{y}|)}{\frac{4}{3}\pi L^3} \ge 0$), mind a 3. lépésben levezetett 3D Yukawa Green-függvény ($G_L(\mathbf{x}-\mathbf{y}) = \frac{e^{-|\mathbf{x}-\mathbf{y}|/L}}{4\pi L^2|\mathbf{x}-\mathbf{y}|} > 0$) **szigorúan pozitív maggal integrálja a fotometriai felületi fényességből származó $\rho_{\text{bary}}(R,z) \ge 0$ sűrűséget**. Így minden galaxis minden sugárpontján:
  $$\bar{\rho}_L(r) > 0, \qquad \phi_L^{\text{Yuk}}(r) > 0, \qquad 0 < f_L(r) \le 1 \quad \checkmark$$
- **Kritikus módszertani megkülönböztetés:** Ez az eredmény a **3D integrálmag matematikai pozitivitásának igazolása** a 175 SPARC-galaxis pozitív fotometriai tömegmodelljén (`n_negative_density_points = 0`), **nem pedig** a konstans $\alpha_Y = 1$ ($\alpha_{\text{std}} = 1/(4\pi)$) Yukawa-hatás forgásigörbe-illesztési sikere! Ahogy az alábbi 3. szakaszban számszerűen kimutatjuk, az izolált, rögzített $\alpha_{\text{std}} = 1/(4\pi)$ csatolású modell önmagában nem reprodukálja a lapos galaktikus forgási görbéket.

---

## 3. Egzakt számszerű negatív kontroll: Izolált rögzített csatolású ($\alpha_{\text{std}} = 1/(4\pi)$) modell vs. NGC 3198 ($R = 30{,}08\text{ kpc}$)

A bírálati észrevételeknek megfelelően az alábbi táblázatban külön oszlopokban tüntetjük fel az `NGC3198_rotmod.dat` ($R = 30{,}08\;\text{kpc}$) kontrollpont megfigyelt sebességét ($V_{\text{obs}}$), barionos sebességét ($V_{\text{bar}}$), a megfigyeléshez **szükséges** erősítést ($C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$), valamint az izolált, rögzített $\gamma_L = G L^2/c^2$ ($\alpha_Y = 1$, $\alpha_{\text{std}} = 1/(4\pi) \approx 0{,}07958$) modell által **ténylegesen jósolt** erősítést ($C_{\text{pred}}^{\text{isol}}$) és forgási sebességet ($V_{\text{pred}}^{\text{isol}} = \sqrt{C_{\text{pred}}^{\text{isol}}}\,V_{\text{bar}}$):

| $\Upsilon_{\text{disk}},\,\Upsilon_{\text{bul}}$ választás | Megfigyelt $V_{\text{obs}}$ ($\text{km/s}$) | Barionos $V_{\text{bar}}$ ($\text{km/s}$) | Szükséges erősítés $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2$ | Jósolt erősítés $C_{\text{pred}}^{\text{isol}}$ (izolált 3D top-hat, $\gamma_L = GL^2/c^2$) | Jósolt erősítés $C_{\text{pred}}^{\text{isol,max}}$ (pontszerű Yukawa felső korlát, $1+\alpha_{\text{std}}$) | Jósolt sebesség $V_{\text{pred}}^{\text{isol}}$ (izolált $\alpha_{\text{std}} = 1/(4\pi)$, $\text{km/s}$) | Megfigyelt mínusz jósolt eltérés $\Delta V = V_{\text{obs}} - V_{\text{pred}}^{\text{isol}}$ |
|---|---|---|---|---|---|---|---|
| **$\Upsilon_{\text{disk}} = 1{,}0,\;\Upsilon_{\text{bul}} = 0$** (nyers SPARC oszlopok) | $146{,}0 \pm 2{,}0$ | $88{,}53$ | **$2{,}720$** | $1 + \mathcal{O}(10^{-7}) \approx 1{,}0000002$ | $\le 1 + \frac{1}{4\pi} = \mathbf{1{,}07958}$ | **$88{,}53\text{--}91{,}99$** | **$-54{,}0\text{ \dots }-57{,}5\;\text{km/s}$** (nem írja le a lapos görbét!) |
| **$\Upsilon_{\text{disk}} = 0{,}5,\;\Upsilon_{\text{bul}} = 0{,}7$** (V16/V17 tömeg–fény arány) | $146{,}0 \pm 2{,}0$ | $70{,}32$ | **$4{,}311$** | $1 + \mathcal{O}(10^{-7}) \approx 1{,}0000002$ | $\le 1 + \frac{1}{4\pi} = \mathbf{1{,}07958}$ | **$70{,}32\text{--}73{,}06$** | **$-72{,}9\text{ \dots }-75{,}7\;\text{km/s}$** (nem írja le a lapos görbét!) |

### Fizikai következtetés (Negatív kontroll):
1. Egy izolált galaxisban a gravitációs potenciál mélysége $\frac{G M_{\text{gal}}}{c^2 R} \sim \frac{V_{\text{bar}}^2}{c^2} \sim 10^{-7}$, és a kanonikus $\gamma_L = \frac{G L^2}{c^2}$ csatolással számolt 3D top-hat csillapítás mind a 175 SPARC-galaxison $\max |1 - f_L| \le 2{,}30\times 10^{-7}$ (`sparc_tophat_summary.json`).
2. A rögzített $\alpha_Y = 1$ ($\alpha_{\text{std}} = \frac{1}{4\pi} \approx 0{,}07958$) Yukawa ötödik erő még az optimális $r \sim L$ átmeneti zónában is legfeljebb $1 + \alpha_{\text{std}} \approx 1{,}0796$-szorosára növeli $V_{\text{bar}}^2$-et (ami $\Upsilon_{\text{disk}}=1{,}0$-nál $91{,}99\;\text{km/s}$-ot, $\Upsilon_{\text{disk}}=0{,}5$-nél $73{,}06\;\text{km/s}$-ot ad a megfigyelt $146\;\text{km/s}$ helyett).
3. Következésképpen az NGC 3198 ($R = 30{,}08\;\text{kpc}$) kontrollpont **egzakt negatív kontrollként** bizonyítja, hogy az izolált, konstans $\alpha_{\text{std}} = 1/(4\pi)$ csatolású V23 Yukawa-modell önmagában **nem** magyarázza meg a galaktikus lapos forgási görbéket; ehhez vagy a V16/V17 fenomenologikus $C(\rho) = 1 + C_0 e^{-\rho/\rho_c}$ fősíkbeli erősítésre ($C_0 \approx 3{,}34\text{--}5{,}09$), vagy egy a konstans-$\alpha_Y$ V23 hatáson túlmutató, sűrűségfüggő csatolású ($\beta_{\text{eff}}(\rho)$, $\alpha_{\text{eff}}^{\text{gal}} = 2\beta_{\text{eff}}^2 \approx 1{,}72\text{--}3{,}31$) kaméleon/szimmetron kiterjesztésre (lásd a 11. lépést, **Branch II-Cham**) van szükség.

---

## 4. Ellenőrzések (Checks)

1. **NGC 3198 $R = 30{,}08\text{ kpc}$ számszerű azonosságok (`C_{\text{req}}$ vs. `C_{\text{pred}}^{\text{isol}}$):**
   - $V_{\text{bar}}(\Upsilon_{\text{disk}}=1{,}0) = 88{,}53\;\text{km/s} \implies C_{\text{req}} = (146 / 88{,}53)^2 = 2{,}7199 \approx 2{,}720$, míg $C_{\text{pred}}^{\text{isol}} \le 1 + \frac{1}{4\pi} = 1{,}07958 \implies V_{\text{pred}}^{\text{isol}} \le 91{,}99\;\text{km/s} < 146\;\text{km/s}$ $\checkmark$.
   - $V_{\text{bar}}(\Upsilon_{\text{disk}}=0{,}5) = 70{,}32\;\text{km/s} \implies C_{\text{req}} = (146 / 70{,}32)^2 = 4{,}3107 \approx 4{,}311$, míg $C_{\text{pred}}^{\text{isol}} \le 1{,}07958 \implies V_{\text{pred}}^{\text{isol}} \le 73{,}06\;\text{km/s} < 146\;\text{km/s}$ $\checkmark$.
2. **Pozitív definitség ellenőrzése az összes (175) SPARC-galaxison:**
   - Fotometriai 3D korongmodell ($\rho_{\text{bary}}(R,z) \ge 0$) + 3D gömbi top-hat vagy 3D Yukawa Green-függvény: $\rho \le 0$ pontok száma **0** mind a 175 galaxisban (`negative_density.md`, 5–7. sor és `hu_p2.tex`, 17–21. sor) $\checkmark$.
