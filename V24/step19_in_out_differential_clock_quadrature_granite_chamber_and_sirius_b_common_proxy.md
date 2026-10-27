# 19. lépés: A V20 $8558\;\text{kg}$-os gránitkamra 3D IN–OUT differenciális kvadratúrája (három geometriai $L$-tartomány, Top-Hat vs. Yukawa) és a Szíriusz B–Szíriusz A közös effektívtömeg-proxy (`common_background` vs. `fixedM`) egzakt levezetése
# Step 19: Exact Derivation of the V20 $8558\;\text{kg}$ Granite-Chamber 3D IN–OUT Differential Quadrature (Three Geometric $L$-Regimes, Top-Hat vs. Yukawa) and the Sirius B–Sirius A Common Effective-Mass Proxy (`common_background` vs. `fixedM`)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (3. szakasz és 3. táblázat + lábjegyzet, 161–189. sor; DOI: `10.5281/zenodo.23639084`)  
- V23 ellenőrző kód és adatsor: `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/covariant_closure_50aa2db0.py` (4. rész, 700–885. sor) és `clock_coupling_H0_H1.csv` (1–13. sor)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (9. fejezet, 426–455. sor és 10. fejezet 1. pont, 463–468. sor; DOI: `10.5281/zenodo.23604423`)

---

## 1. Szó szerinti idézetek a V22 és V23 kéziratokból

A 8. lépésben (`step8_microscopic_clock_coupling_siriusb_lab_v20.md`) levezettük a Dirac–Maxwell atomóra-hatást ($S_{\text{atomic}}[\tilde{g}_{\mu\nu}, \psi_e, A_\mu]$) H0 ($\tilde{g}_{\mu\nu} = g_{\mu\nu}$) és feltételes H1 ($\tilde{g}_{00} = f_L^2 g_{00}$) alatt. A V23 kézirat 3. szakaszának 3. táblázata (170–189. sor) és a `covariant_closure_50aa2db0.py` 4. része (700–885. sor) azonban két további egzakt matematikai szerkezetet tartalmaz, amelyeket eddig még nem vezettünk le lépésről lépésre:

1. **A V23 3. táblázat lábjegyzete a gránitkamra IN–OUT differenciális jeléről és a Szíriusz B közös effektívtömeg-proxyról (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 188. sor):**
   > *„V20 sorok: H0 oszlop = a newtoni GR feletti TAL-többlet IN ($\mathbf{r}=\mathbf{0}$) állásban; H1 oszlopok = IN-eltolódás (IN$-$OUT differenciális jel, $D_{\text{out}}=2{,}0\,\text{m}$). Newtoni GR: $\Phi_N(\text{IN})/c^2 = -9{,}269\times 10^{-24}$, IN$-$OUT $= -6{,}094\times 10^{-24}$; a $\sigma_y=10^{-18}$ szemléltető referenciaérték. Szíriusz B közös effektívtömeg-proxyn (rögzített $n=3$ Lane–Emden-profil a H0 effektívtömeg-eltolással, nem újonnan megoldott egyensúly) $L=10^3\,\text{km}$-nél: H0 $80{,}69839$, H1 top-hat $80{,}69878\,\text{km/s}$.”*

2. **A $8558\;\text{kg}$-os gránitkamra geometriai és tömegmérlege (`tal_book_HU_2026-10-24.tex`, 10. fejezet 1. pont, 466. sor):**
   > *„**Pontos geometria- és tömegmérleg:** Az $R_{\text{in}} = 0{,}25\,\text{m}$, $R_{\text{out}} = 0{,}85\,\text{m}$, $H = 1{,}20\,\text{m}$ méretű nyitott hengerpalást térfogata $V_{\text{cyl}} = \pi(0{,}85^2 - 0{,}25^2)\times 1{,}20 = 2{,}48814\,\text{m}^3$, tömege **$M_{\text{cyl}} = 6842{,}39\,\text{kg}$**; a névleges $V_{\text{tot}} = 3{,}1102\,\text{m}^3$ (**$M_{\text{tot}} = 8558\,\text{kg}$**) értékhez további $0{,}6220\,\text{m}^3$ ($1710{,}6\,\text{kg}$) záró/tartó gránitelem (vagy $R_{\text{out}} = 0{,}942\,\text{m}$ külső sugár) tartozik.”*

---

## 2. A V20 $8558\;\text{kg}$-os gránitkamra 3D IN–OUT kvadratúrája és a három geometriai $L$-tartomány levezetése

### 2.1. A hengerpalást + zárókorongok ($\Omega_{\text{granite}}$) tartománya és az off-axis távolságfüggvény
A `covariant_closure_50aa2db0.py` (802–822. sor) a $\rho_{\text{granite}} = 2750\;\text{kg/m}^3$ sűrűségű, $M_{\text{tot}} = 8558\;\text{kg}$ tömegű ($V_{\text{tot}} = M_{\text{tot}}/\rho_{\text{granite}} = 8558/2750 = 3{,}112000\;\text{m}^3$) gránitkamrát hengerkoordinátákban ($\varrho \in [0, R_{\text{out}}]$, $z \in [-z_{\text{out}}, +z_{\text{out}}]$, $\varphi \in [0, 2\pi)$) a $z \leftrightarrow -z$ tükörszimmetria kihasználásával a következő tartományként építi fel:
- **Nyitott hengerpalást (barrel):** $R_{\text{in}} = 0{,}25\;\text{m} \le \varrho \le R_{\text{out}} = 0{,}85\;\text{m}$, $|z| \le z_{\text{in}} = 0{,}60\;\text{m}$ ($H_{\text{in}} = 2 z_{\text{in}} = 1{,}20\;\text{m}$), amelynek térfogata és tömege:
  $$V_{\text{cyl}} \;=\; \pi\bigl(R_{\text{out}}^2 - R_{\text{in}}^2\bigr)\,(2 z_{\text{in}}) \;=\; 2{,}488141\;\text{m}^3, \qquad M_{\text{cyl}} \;=\; \rho_{\text{granite}}\,V_{\text{cyl}} \;=\; 6842{,}39\;\text{kg}.$$
- **Két tömör zárókorong (endcaps):** $0 \le \varrho \le R_{\text{out}} = 0{,}85\;\text{m}$, $z_{\text{in}} \le |z| \le z_{\text{out}}$, ahol a zárókorongok vastagsága $\Delta z_{\text{cap}} = z_{\text{out}} - z_{\text{in}}$ (`covariant_closure_50aa2db0.py`, 808. sor):
  $$\Delta z_{\text{cap}} \;=\; \frac{V_{\text{tot}} - V_{\text{cyl}}}{2\pi R_{\text{out}}^2} \;=\; \frac{3{,}112000 - 2{,}488141}{2\pi\times 0{,}85^2} \;=\; 0{,}137426\;\text{m} \quad\Longrightarrow\quad z_{\text{out}} \;=\; 0{,}60 + 0{,}137426 \;=\; 0{,}737426\;\text{m}.$$
  Így a gránit által elfoglalt $\Omega_{\text{granite}}$ tartomány karakterisztikus feltétele a fél-meridiánsíkon ($\varrho \in [0, R_{\text{out}}]$, $z \in [0, z_{\text{out}}]$):
  $$\Omega_{\text{granite}} \;=\; \bigl\{(\varrho, z) \in [0, R_{\text{out}}]\times [0, z_{\text{out}}] \;:\; (\varrho \ge R_{\text{in}}) \;\lor\; (z \ge z_{\text{in}})\bigr\}, \qquad dV \;=\; 2\cdot\bigl(2\pi\varrho\,d\varrho\,dz\bigr),$$
  ahol a $2$-es szorzó a $\pm z$ félterek összegzéséből származik ($N_\varrho \times N_z = 500\times 500$ cellán $M_{\text{granite,num}} = 8558{,}23\;\text{kg}$).

Egy $(\varrho, \varphi, \pm z) \in \Omega_{\text{granite}}$ forráspont euklideszi távolsága a két óraállástól (`covariant_closure_50aa2db0.py`, 820. és 834. sor):
1. **IN állás ($\mathbf{r}_{\text{in}} = (0,0,0)$, a kamra geometriai középpontja):** az azimutszögtől függetlenül:
   $$d_{\text{in}}(\varrho, z) \;=\; \sqrt{\varrho^2 + z^2}.$$
2. **OUT állás ($\mathbf{r}_{\text{out}} = (D_{\text{out}}, 0, 0)$, $D_{\text{out}} = 2{,}0\;\text{m}$ a fősíkban):** a koszinusztétel alapján:
   $$d_{\text{out}}(\varrho, z, \varphi) \;=\; \sqrt{\varrho^2 + D_{\text{out}}^2 - 2\varrho D_{\text{out}}\cos\varphi + z^2},$$
   amelyet a kód $N_\varphi = 256$ egyenletes azimutszög-ponton ($\varphi_k = \frac{2\pi(k + 1/2)}{256}$, $dm_p = dm_g / N_\varphi$) integrál. Ebből a newtoni gravitációs potenciál és a standard GR vöröseltolódás-különbség (`clock_coupling_H0_H1.csv`, 8–13. sor):
   $$\frac{\Phi_N(\text{IN})}{c^2} \;=\; -\frac{G}{c^2}\int_{\Omega_{\text{granite}}}\frac{dm_g}{d_{\text{in}}} \;=\; -9{,}268639\times 10^{-24}, \qquad \frac{\Phi_N(\text{OUT})}{c^2} \;=\; -\frac{G}{c^2}\int_{\Omega_{\text{granite}}}\frac{dm_p}{d_{\text{out}}} \;=\; -3{,}174686\times 10^{-24},$$
   $$\Delta y_{\text{GR}}(\text{IN}-\text{OUT}) \;=\; \frac{\Phi_N(\text{IN}) - \Phi_N(\text{OUT})}{c^2} \;=\; -6{,}093953\times 10^{-24}.$$

### 2.2. A négy geometriai küszöbtávolság és a Top-Hat vs. Yukawa IN–OUT viselkedés egzakt bizonyítása
A kamra $\Omega_{\text{granite}}$ tartományának az IN és OUT pontoktól mért szélső euklideszi távolságai:
- **IN állásban ($\mathbf{r}_{\text{in}} = \mathbf{0}$):**
  $$d_{\min}^{\text{IN}} \;=\; \min_{(\varrho,z)\in\Omega_{\text{granite}}}\sqrt{\varrho^2 + z^2} \;=\; \min(R_{\text{in}}, z_{\text{in}}) \;=\; \min(0{,}25,\,0{,}60) \;=\; 0{,}25\;\text{m},$$
  $$d_{\max}^{\text{IN}} \;=\; \max_{(\varrho,z)\in\Omega_{\text{granite}}}\sqrt{\varrho^2 + z^2} \;=\; \sqrt{R_{\text{out}}^2 + z_{\text{out}}^2} \;=\; \sqrt{0{,}85^2 + 0{,}737426^2} \;=\; 1{,}12530\;\text{m}.$$
- **OUT állásban ($D_{\text{out}} = 2{,}00\;\text{m}$):**
  $$d_{\min}^{\text{OUT}} \;=\; \min_{\Omega_{\text{granite}}} d_{\text{out}} \;=\; D_{\text{out}} - R_{\text{out}} \;=\; 2{,}00 - 0{,}85 \;=\; 1{,}15\;\text{m} \quad (\varphi = 0,\; \varrho = R_{\text{out}},\; z = 0),$$
  $$d_{\max}^{\text{OUT}} \;=\; \max_{\Omega_{\text{granite}}} d_{\text{out}} \;=\; \sqrt{(D_{\text{out}} + R_{\text{out}})^2 + z_{\text{out}}^2} \;=\; \sqrt{2{,}85^2 + 0{,}737426^2} \;=\; 2{,}94385\;\text{m} \quad (\varphi = \pi,\; \varrho = R_{\text{out}},\; z = z_{\text{out}}).$$

Ezekből a küszöbtávolságokból ($0{,}25\;\text{m} < 1{,}1253\;\text{m} < 1{,}15\;\text{m} < 2{,}9439\;\text{m}$) közvetlenül adódik a 3D gömbi top-hat ($\phi_L^{\text{top-hat}}(\mathbf{r}) = \frac{3 G}{4\pi c^2 L}\int_{d \le L} dm$) és a lokális Yukawa-mag ($\phi_L^{\text{Yuk}}(\mathbf{r}) = \frac{G}{4\pi c^2}\int \frac{e^{-d/L}}{d}\,dm$ $\alpha_Y = 1$ mellett) közötti éles fizikai különbség (`clock_coupling_H0_H1.csv`, 8–13. sor):

1. **Üreg-tartomány ($L < d_{\min}^{\text{IN}} = 0{,}25\;\text{m}$, pl. $L = 0{,}20\;\text{m}$):**
   A $B(\mathbf{0}, L)$ és $B(\mathbf{r}_{\text{out}}, L)$ gömbök egyaránt a gránitfalon kívül (az üreg belsejében, illetve a kamrán kívül) maradnak, így $\phi_L^{\text{top-hat}}(\text{IN}) = \phi_L^{\text{top-hat}}(\text{OUT}) = 0$, tehát $\Delta y_{\text{H1}}^{\text{top-hat}} = 0$. Ezzel szemben a Yukawa-féle $e^{-d/L}/d$ Green-függvény az üregen át is érzékeli a falat: $\Delta y_{\text{H1}}^{\text{Yuk}}(\text{IN}) = -3{,}690\times 10^{-26}$ és $\Delta y_{\text{H1}}^{\text{Yuk}}(\text{IN}-\text{OUT}) = -3{,}683\times 10^{-26}$.
2. **Tisztán belső top-hat tartomány ($d_{\min}^{\text{IN}} < L < d_{\min}^{\text{OUT}} = 1{,}15\;\text{m}$, pl. $L = 0{,}50,\;0{,}85,\;1{,}00\;\text{m}$ a V23 3. táblázatban):**
   Mivel $L < 1{,}15\;\text{m}$, az OUT állás $B(\mathbf{r}_{\text{out}}, L)$ gömbje **nem metszi a gránitkamrát** ($d_{\text{out}} \ge 1{,}15\;\text{m} > L$), így:
   $$\boxed{\;\phi_L^{\text{top-hat}}(\text{OUT}) \;\equiv\; 0 \quad\Longrightarrow\quad \Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}-\text{OUT}) \;\equiv\; \Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}) \qquad \forall\,L < 1{,}15\;\text{m}!\;}$$
   Ezért áll a V23 3. táblázatában (183–185. sor) pontosan ugyanaz a szám a zárójel előtt és a zárójelben a Top-Hat oszlopban: $-3{,}318\times 10^{-25}\;(-3{,}318\times 10^{-25})$ ($L=0{,}50\;\text{m}$), $-1{,}303\times 10^{-24}\;(-1{,}303\times 10^{-24})$ ($L=0{,}85\;\text{m}$) és $-1{,}436\times 10^{-24}\;(-1{,}436\times 10^{-24})$ ($L=1{,}00\;\text{m}$). Ezzel szemben a Yukawa-mag exponenciális farka $D_{\text{out}} = 2{,}0\;\text{m}$-nél sem zérus ($\phi_L^{\text{Yuk}}(\text{OUT}) = +3{,}735\times 10^{-26}$ $L = 1{,}00\;\text{m}$-nél), ezért a Yukawa differenciális jel kisebb az abszolút IN-jelnél: $-3{,}781\times 10^{-25}\;(-3{,}407\times 10^{-25})$ $L = 1{,}00\;\text{m}$-nél!
3. **Részleges külső átfedés ($d_{\min}^{\text{OUT}} < L < d_{\max}^{\text{OUT}}$, pl. $L = 1{,}50\;\text{m}$):**
   Itt az IN gömb már az egész kamrát befogja ($L = 1{,}50\;\text{m} > d_{\max}^{\text{IN}} = 1{,}1253\;\text{m}$), míg az OUT gömb részben metszi a kamrát: $\Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}) = -1{,}0115\times 10^{-24}$, $\Delta y_{\text{H1}}^{\text{top-hat}}(\text{OUT}) = -1{,}0064\times 10^{-25}$, így $\Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}-\text{OUT}) = -9{,}1086\times 10^{-25}$.
4. **Közös bezárási tartomány ($L \ge d_{\max}^{\text{OUT}} = 2{,}94385\;\text{m}$, pl. $L = 10{,}0\;\text{m}$):**
   Amikor $L \ge 2{,}944\;\text{m}$, mind a $B(\mathbf{0}, L)$, mind a $B(\mathbf{r}_{\text{out}}, L)$ gömb hiánytalanul tartalmazza az egész $M_{\text{tot}} = 8558\;\text{kg}$ gránitkamrát! Ezért a 3D gömbi top-hat átlag mindkét állásban pontosan $\bar{\rho}_L^{\text{top-hat}}(\text{IN}) = \bar{\rho}_L^{\text{top-hat}}(\text{OUT}) = \frac{3 M_{\text{tot}}}{4\pi L^3}$, így a **Top-Hat IN$-$OUT differenciális jel azonosan eltűnik**:
   $$\boxed{\;\Delta y_{\text{H1}}^{\text{top-hat}}(\text{IN}-\text{OUT})\Big|_{L \ge 2{,}944\;\text{m}} \;=\; 0\quad(\text{numerikusan } -5{,}28\times 10^{-40}\text{ kerekítési zaj } L=10\;\text{m}\text{-nél}),\;}$$
   miközben a Yukawa-mag a $1/d$ súlyozás miatt $L = 10{,}0\;\text{m}$-nél is határozott differenciális jelet ad: $\Delta y_{\text{H1}}^{\text{Yuk}}(\text{IN}-\text{OUT}) = -4{,}818\times 10^{-25}$!
5. **A másodrendű H0 TAL-többlet az üregben (`dy_H0_cavity`):**
   H0 alatt nincs közvetlen $f_L$ óraszorzó, így a newtoni GR feletti többlet a kamra effektív gravitáló tömegének $1 - f_L \simeq \phi_L$ módosításából származik (`covariant_closure_50aa2db0.py`, 852. sor):
   $$\Delta y_{\text{H0}}^{\text{excess}}(\text{IN}) \;=\; \left|\frac{\Phi_N(\text{IN})}{c^2}\right|\,\phi_L^{\text{top-hat}}(\text{IN}) \;=\; 9{,}268639\times 10^{-24}\times 1{,}436124\times 10^{-24} \;=\; +1{,}331091\times 10^{-47} \quad (L = 1{,}00\;\text{m}).$$

---

## 3. A Szíriusz B–Szíriusz A differenciális vöröseltolódás és a `common_background` vs. `fixedM` proxy levezetése

A `covariant_closure_50aa2db0.py` (710–799. sor) a Szíriusz B fehér törpét ($M_B = 1{,}018 \pm 0{,}011\,M_\odot$, $R_B = 0{,}00803 \pm 0{,}00011\,R_\odot = 5586{,}471\;\text{km}$) és a Szíriusz A fősorozati referenciacsillagot ($M_A = 2{,}042\,M_\odot$, $R_A = 1{,}711\,R_\odot = 1{,}19034\times 10^6\;\text{km}$) $n = 3$ Lane–Emden politrópként modellezi:
1. **A GR alapvonal és a kombinált bizonytalanság (`covariant_closure_50aa2db0.py`, 719–723. sor):**
   A felszíni kompaktság $\mathcal{C}_B = \frac{GM_B}{R_B c^2} = 2{,}690878\times 10^{-4}$, amiből az egzakt Schwarzschild-vöröseltolódás és sebesség-egyenérték:
   $$z_{\text{GR},B} \;=\; \bigl(1 - 2\mathcal{C}_B\bigr)^{-1/2} - 1 \;=\; 2{,}6918805\times 10^{-4}, \qquad v_{\text{GR}} \;=\; c\,z_{\text{GR},B} \;=\; 80{,}700550\;\text{km/s}.$$
   A paraméter-bizonytalanság $\sigma_{\text{GR}} = v_{\text{GR}}\sqrt{(\sigma_M/M_B)^2 + (\sigma_R/R_B)^2} = 1{,}4080\;\text{km/s}$, amely a megfigyelési $\sigma_{\text{obs}} = 0{,}77\;\text{km/s}$ hibával kvadratikusan összegezve $\sigma_{\text{comb}} = \sqrt{\sigma_{\text{obs}}^2 + \sigma_{\text{GR}}^2} = 1{,}604806\;\text{km/s}$-ot ad (`clock_coupling_H0_H1.csv`, 2. sor).
2. **A Szíriusz B–Szíriusz A differenciális H1 óraeltolódás (`covariant_closure_50aa2db0.py`, 767–773. sor):**
   Mivel az asztrometriai gravitációs vöröseltolódást a Szíriusz B és a Szíriusz A színképvonalainak különbségéből mérik, és H1 alatt a felszíni vöröseltolódás $1 + z_{\text{H1}} = (1 + z_{\text{metric}})\,f_L(R_*)^{-1} = (1 + z_{\text{metric}})\,e^{+x_*}$ (ahol $x_* \equiv \gamma_L\bar{\rho}_L(R_*)$), a közvetlen H1 többletsebesség a két csillag között:
   $$\boxed{\;\Delta v_{\text{H1,extra}}(L) \;=\; \Delta v_B(L) \;-\; \Delta v_A(L) \;=\; c\,(1 + z_{\text{GR},B})\bigl(e^{x_B(L)} - 1\bigr) \;-\; c\,\bigl(e^{x_A(L)} - 1\bigr).\;}$$
   - $L = 1000\;\text{km} = 10^6\;\text{m}$-nél: $\Delta v_{\text{H1,extra}}^{\text{top-hat}} = +0{,}000392\;\text{km/s}$ és $\Delta v_{\text{H1,extra}}^{\text{Yuk}} = +0{,}044215\;\text{km/s}$.
   - $L = 5000\;\text{km} = 5\times 10^6\;\text{m}$-nél: $\Delta v_{\text{H1,extra}}^{\text{top-hat}} = +4{,}688557\;\text{km/s}$ ($v_{\text{H1}}^{\text{top-hat}} = 85{,}38911\;\text{km/s}$, $2{,}95\sigma$) és $\Delta v_{\text{H1,extra}}^{\text{Yuk}} = +2{,}150821\;\text{km/s}$ ($v_{\text{H1}}^{\text{Yuk}} = 82{,}85137\;\text{km/s}$, $1{,}37\sigma$).
3. **A rögzített tömegű (`fixedM`) és a közös effektívtömeg-proxy (`common_background`) közötti különbség (`covariant_closure_50aa2db0.py`, 758–766. és 794–798. sor; V23 188. sor):**
   - Ha a csillag $M_B = 1{,}018\,M_\odot$ gravitáló tömegét a pályaasztrometriából rögzítettnek tekintjük (`fixedM`), akkor a külső Schwarzschild-metrika változatlan, így $v_{\text{H0,fixedM}} \equiv v_{\text{GR}} = 80{,}70055\;\text{km/s}$, és $v_{\text{H1,fixedM}}(L) = v_{\text{GR}} + \Delta v_{\text{H1,extra}}(L)$ ($80{,}70094\;\text{km/s}$ $L = 10^3\;\text{km}$-nél).
   - Ha viszont rögzített $n = 3$ Lane–Emden barionos $\rho_B(r)$ profil mellett számoljuk ki a `mass` módbeli effektív tömegcsökkenést (`common_background` proxy):
     $$\frac{\Delta M_{\text{eff}}(L)}{M_B} \;=\; \frac{\int_0^{R_B} 4\pi r^2\rho_B(r)\,\bigl[f_L(r) - 1\bigr]\,dr}{\int_0^{R_B} 4\pi r^2\rho_B(r)\,dr} \;<\; 0,$$
     akkor a H0 előrejelzés $v_{\text{H0,mass}}(L) = v_{\text{GR}}\bigl(1 + \frac{\Delta M_{\text{eff}}(L)}{M_B}\bigr)$ ($80{,}69839\;\text{km/s}$ $L = 10^3\;\text{km}$-nél és $80{,}69487\;\text{km/s}$ $L = 5\times 10^3\;\text{km}$-nél), míg az ugyanerre a közös effektívtömeg-proxyra épülő H1 érték:
     $$\boxed{\;v_{\text{H1,common}}(L) \;=\; v_{\text{H0,mass}}(L) \;+\; \Delta v_{\text{H1,extra}}(L) \;=\; 80{,}698385 + 0{,}000392 \;=\; 80{,}69878\;\text{km/s}\quad(L = 10^3\;\text{km}\text{-nél, Top-Hat}).\;}$$

---

## 4. Ellenőrzések (Checks)

1. **A Top-Hat és a Yukawa $\phi_L(\text{IN})$ aránya nagy $L$-nél ($L = 10{,}0\;\text{m} \gg d_{\max}^{\text{IN}} = 1{,}125\;\text{m}$):**
   - Amikor $L \gg d_{\max}^{\text{IN}}$, a 3D gömbi top-hat átlagban az egész $M_{\text{tot}}$ tömeg benne van a $\frac{4}{3}\pi L^3$ gömbben, így $\phi_L^{\text{top-hat}}(\text{IN}) = \frac{G L^2}{c^2}\frac{3 M_{\text{tot}}}{4\pi L^3} = \frac{3 G M_{\text{tot}}}{4\pi c^2 L} = 1{,}517259\times 10^{-25}$ ($L = 10\;\text{m}$-nél), míg a Yukawa-integrálban $e^{-d_{\text{in}}/L} \approx 1 - d_{\text{in}}/L$, így:
     $$\phi_L^{\text{Yuk}}(\text{IN})\Big|_{L \gg d_{\max}^{\text{IN}}} \;=\; \frac{G}{4\pi c^2}\int_{\Omega_{\text{granite}}}\frac{e^{-d_{\text{in}}/L}}{d_{\text{in}}}\,dm_g \;\approx\; \frac{1}{4\pi}\left|\frac{\Phi_N(\text{IN})}{c^2}\right| - \frac{G M_{\text{tot}}}{4\pi c^2 L} \;=\; 7{,}3757\times 10^{-25} - 0{,}5058\times 10^{-25} \;\approx\; 6{,}87\times 10^{-25},$$
     amely kiválóan egyezik a `clock_coupling_H0_H1.csv` 13. sorában szereplő egzakt kvadratúra-értékkel ($6{,}888137\times 10^{-25}$) $\checkmark$.
2. **A Szíriusz B `common_background` és `fixedM` különbség ellenőrzése $L = 10^3\;\text{km}$-nél:**
   - $v_{\text{H1,fixedM}}^{\text{top-hat}} - v_{\text{H1,common}}^{\text{top-hat}} = 80{,}70094147 - 80{,}69877704 = +0{,}00216443\;\text{km/s} \equiv v_{\text{GR}} - v_{\text{H0,mass}}$ ($80{,}70054975 - 80{,}69838532 = +0{,}00216443\;\text{km/s}$) gépi pontossággal $\checkmark$.
