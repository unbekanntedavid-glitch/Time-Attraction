# 20. lépés: A V20 nyolctekercses mágnesestér-kísérlet elliptikus Biot–Savart hurokintegrálja, az $r = L/t$ kompaktifikált külső kvadratúra ($3{,}7\times 10^{-16}$ / $7{,}4\times 10^{-16}$ egyezés) és az $L = 0{,}95\;\text{m}$-es átkapcsolási azonosság egzakt levezetése
# Step 20: Exact Derivation of the V20 Eight-Coil Elliptic Biot–Savart Loop Integral, the $r = L/t$ Compactified Exterior Quadrature ($3.7\times 10^{-16}$ / $7.4\times 10^{-16}$ Agreement), and the $L = 0.95\;\text{m}$ Switching Identity

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step19_in_out_differential_clock_quadrature_granite_chamber_and_sirius_b_common_proxy.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (10. fejezet 2. pont, 470–485. sor és Függelék, 557–618. sor; DOI: `10.5281/zenodo.23604423`)  
- V20 (`10.5281/zenodo.23540741`) ellenőrző kód és adatsor: `/workspace/ta-v18/V20/make_revision3_checks.py` (1–77. sor) és `/workspace/ta-v18/V20/revision3_checks.json` (123–197. sor)  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (28. sor; DOI: `10.5281/zenodo.23639084`)

---

## 1. Szó szerinti idézetek a V22 és V20 kéziratokból

A 15. lépésben (`step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md`) levezettük, hogy a nyommentes Maxwell-tenzor miért kelt $\rho_B(\mathbf{x}) = \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0 c^2}$ forrást a V23 diszform metrikában, a 19. lépésben pedig a V20 „A” ágának ($8558\;\text{kg}$-os gránitkamra) 3D IN–OUT kvadratúráját zártuk le. A V22 monográfia 10. fejezetének 2. pontja (470–485. sor) és Függeléke (557–618. sor), valamint a V20 `make_revision3_checks.py` (25–75. sor) tartalmazza a V20 „B” ágának (a nyolctekercses mágnesestér-energetikai kísérletnek) egzakt elliptikus integrálos térszámítását, a két külső kvadratúra-utat (A és B út) és az $L = 0{,}95\;\text{m}$-es belső/külső átkapcsolási azonosságot:

1. **A 8 tekercses elrendezés paraméterei és a $L = 0{,}630957\;\text{m}$-es rácsmaximum (`tal_book_HU_2026-10-24.tex`, 10. fejezet 2. pont, 473–478. sor):**
   > *„**Auditált 8 tekercses hardverparaméterek:** 8 darab egyforma szolenoid, egyenként $l_{\text{coil}} = 400\,\text{mm}$ hosszúsággal, külső szerkezeti **$\varnothing 300\,\text{mm}$ peremhatár-átmérővel**, effektív **$\varnothing 291{,}7\,\text{mm}$ réz menetátmérővel** ($R_{\text{coil}} = 0{,}14585\,\text{m}$), $N = 2000$ menet $0{,}15\,\text{mm}$-es zománcozott rézhuzalból, tekercsenként $R_{\text{DC}} = 1789\,\Omega$ egyenáramú ellenállással, $I = 14{,}0\,\text{mA}$ árammal meghajtva ($0{,}351\,\text{W}$/tekercs, összesen $2{,}805\,\text{W}$).”*
   > *„**Maximum a mintavételezett 25 pontos $L$-rácson ($L = 0{,}630957\,\text{m}$):** A 8 tekercs mágneses energiasűrűségét a 25 pontos logaritmikus $L \in [0{,}1, 3{,}0]\,\text{m}$ rácson kiértékelve a feltételes polaritás-páros $\Delta y(L) = \frac{1}{2}[\Delta y(+I) + \Delta y(-I)]$ eltolódás a rácson **$L = 0{,}630957\,\text{m}$**-nél veszi fel diszkrét maximumát (nem folytonos szélsőérték-illesztés):*
   $$\begin{equation}
     \Delta y_{\text{rácsmax}}(L = 0{,}630957\,\text{m}) \;=\; -1{,}2017108\times 10^{-48}. \tag{V22 Eq. 14}
   \end{equation}$$
   > *Itt elkülönítem: (i) a két determinisztikus integráló rutin $3{,}7\times 10^{-16}$-os egyezését $L=3{,}0\,\text{m}$-nél, (ii) a Sobol-kvadratúra $-0{,}0722\%$-os ellenőrzését, (iii) az $L=0{,}95\,\text{m}$-es $1{,}23\times 10^{-3}$ illesztési maradékot és (iv) a valós mérési érzékenységet.”*

2. **Az $r = L/t$ kompaktifikált külső integrál (A és B út) és az $L = 0{,}95\;\text{m}$-es belső/külső váltás (`tal_book_HU_2026-10-24.tex`, Függelék, 581. és 615. sor):**
   - 581. sor:
     > *„Az új A út panelintegrált és végtelenbe menő farkat, a B út közvetlen r=L/t helyettesítést használ. A 2,6--3,01 m-es tartományban a célzott újrafuttatás a két utat körülbelül **$10^{-15}$ relatív pontossággal** azonosnak adja.”*
   - 615. sor:
     > *„**B. Az új, kisebb illesztési hiba közölt és valós.** A jelentés már feltünteti, hogy L=0,95 m-nél a belső/külső váltás nem pontosan folytonos. Az archivált értékek **2,17738947 és 2,17597607 J/A$^2$**, különbségük **0,00141340 J/A$^2$**, azaz a belső érték **0,0649\%-a**. Ez a összenergia-becslések eltérésével algebrailag egyezik; nem a régi 3 m-es nagy hiba fennmaradása.”*

---

## 2. A köráram-hurkok egzakt elliptikus Biot–Savart tere és a 8 tekercses szuperpozíció

### 2.1. Egyetlen $a = D_{\text{menet}}/2$ sugarú köráram-hurok off-axis mágneses tere (`make_revision3_checks.py`, 30–35. sor)
Tekintsünk egy $a = D_{\text{menet}}/2 = 0{,}2917 / 2 = 0{,}14585\;\text{m}$ sugarú, az $xy$-síkban fekvő köráram-hurkot, amelyben $I$ áram folyik! Hengerkoordinátákban ($\varrho = \sqrt{x^2+y^2}$, $z$) az azimutális vektorpotenciál $A_\varphi(\varrho, z)$ a Biot–Savart-törvényből:
$$A_\varphi(\varrho, z) \;=\; \frac{\mu_0 I\,a}{4\pi}\int_0^{2\pi}\frac{\cos\phi'\,d\phi'}{\sqrt{a^2 + \varrho^2 + z^2 - 2 a\varrho\cos\phi'}}.$$
A $\phi' = \pi - 2\theta$ helyettesítéssel ($\cos\phi' = 2\sin^2\theta - 1$) és az elliptikus paraméter (`make_revision3_checks.py`, 31–32. sor):
$$m(\varrho, z; a) \;\equiv\; \frac{4 a\varrho}{(a + \varrho)^2 + z^2} \;\in\; [0, 1), \qquad d(\varrho, z; a) \;\equiv\; \sqrt{(a + \varrho)^2 + z^2}, \qquad q(\varrho, z; a) \;\equiv\; (a - \varrho)^2 + z^2$$
bevezetésével a $\mathbf{B} = \nabla\times(A_\varphi\hat{\boldsymbol{\varphi}}) = -\frac{\partial A_\varphi}{\partial z}\hat{\boldsymbol{\varrho}} + \frac{1}{\varrho}\frac{\partial(\varrho A_\varphi)}{\partial\varrho}\hat{\mathbf{z}}$ komponensek zárt alakban kifejezhetők az első- és másodfajú elliptikus integrálokkal, $K(m) = \int_0^{\pi/2}(1 - m\sin^2\theta)^{-1/2}d\theta$ és $E(m) = \int_0^{\pi/2}(1 - m\sin^2\theta)^{1/2}d\theta$ (`make_revision3_checks.py`, 33–34. sor):
$$\boxed{\;\begin{aligned}
  B_z(\varrho, z; a, I) &\;=\; \frac{\mu_0 I}{2\pi\,d(\varrho,z;a)}\left[K(m) \;+\; \frac{a^2 - \varrho^2 - z^2}{q(\varrho,z;a)}\,E(m)\right], \\[6pt]
  B_\varrho(\varrho, z; a, I) &\;=\; \frac{\mu_0 I\,z}{2\pi\,\varrho\,d(\varrho,z;a)}\left[-K(m) \;+\; \frac{a^2 + \varrho^2 + z^2}{q(\varrho,z;a)}\,E(m)\right] \qquad (\varrho > 0,\;\; B_\varrho(0,z) = 0).
\end{aligned}\;}$$

### 2.2. A $N_c = 8 \times N = 2000$ hurok szuperpozíciója és a Derby–Olbert áramlemez-határérték (`make_revision3_checks.py`, 25–51. sor)
A V20 elrendezésben a $N_c = 8$ szolenoid középpontjai az $R_{\text{ax}} = 0{,}45\;\text{m}$ sugarú körön, $z_c = 0{,}200 - 0{,}300 = -0{,}100\;\text{m}$ magasságban helyezkednek el:
$$\mathbf{c}_k \;=\; \left(R_{\text{ax}}\cos\frac{2\pi k}{8},\; R_{\text{ax}}\sin\frac{2\pi k}{8},\; z_c\right), \qquad k = 0, 1, \dots, 7.$$
Mindegyik tekercs $\ell = 0{,}400\;\text{m}$ hosszúságú és $N = 2000$ menetet tartalmaz a $z_n = \bigl(\frac{n + 1/2}{N} - \frac{1}{2}\bigr)\ell$ ($n = 0,\dots,N-1$) középhelyeken. Tetszőleges $\mathbf{P} = (x,y,z)$ pontban az eredő mágneses tér:
$$\mathbf{B}_{\text{loops}}(\mathbf{P}; I_0) \;=\; \sum_{k=0}^{7}\sum_{n=0}^{N-1}\left[B_\varrho(\varrho_k, z - z_c - z_n; a, I_0)\,\frac{(x - c_{k,x})\hat{\mathbf{x}} + (y - c_{k,y})\hat{\mathbf{y}}}{\varrho_k} \;+\; B_z(\varrho_k, z - z_c - z_n; a, I_0)\,\hat{\mathbf{z}}\right],$$
ahol $\varrho_k = \sqrt{(x - c_{k,x})^2 + (y - c_{k,y})^2}$. Amikor $N \to \infty$ rögzített $N I_0 / \ell$ felületi áramsűrűség mellett, a középérték-összeg $\frac{1}{N}\sum_{n=0}^{N-1} f(z_n) = \frac{1}{\ell}\int_{-\ell/2}^{+\ell/2} f(z')\,dz' + \mathcal{O}(N^{-2})$ másodrendű pontossággal tart a folytonos Derby–Olbert hengeres áramlemez-integrálhoz (`fields.coils_B`). Ennek megfelelően $I_0 = 14{,}0\;\text{mA}$-nél a két független számítás (`revision3_checks.json`, 124–136. sor):
- **Fő órán ($\mathbf{P}_{\text{main}} = (0,0,0)$):** $|\mathbf{B}_{\text{loops}}| = 11{,}63290030\;\mu\text{T}$, $|\mathbf{B}_{\text{sheet}}| = 11{,}63289980\;\mu\text{T}$ (relatív eltérés $-4{,}28\times 10^{-8}$).
- **Referencia-órán ($\mathbf{P}_{\text{ref}} = (1{,}60, 0, 0)\;\text{m}$):** $|\mathbf{B}_{\text{loops}}| = 0{,}430539231\;\mu\text{T}$, $|\mathbf{B}_{\text{sheet}}| = 0{,}430539228\;\mu\text{T}$ (relatív eltérés $-7{,}56\times 10^{-9}$).

---

## 3. A differenciális mágneses energia ($\Delta W(L)$), az $r = L/t$ kompaktifikált külső kvadratúra és az $L = 0{,}95\;\text{m}$-es azonosság

### 3.1. A polaritás-páros relatív frekvenciaeltolódás kifejezése $\Delta W(L)$-lel (`make_revision3_checks.py`, 53–58. sor)
Definíció szerint az egységnyi áramra ($I = 1\;\text{A}$) normált mágneses energia egy $\mathbf{P}$ középpontú, $L$ sugarú $B(\mathbf{P}, L)$ gömbben:
$$W_{\text{in}}(\mathbf{P}; L) \;\equiv\; \frac{1}{2\mu_0 I_0^2}\int_{|\mathbf{x}-\mathbf{P}|\le L}|\mathbf{B}(\mathbf{x}; I_0)|^2\,d^3x \qquad [\text{J/A}^2],$$
a fő óra ($\mathbf{P}_{\text{main}} = \mathbf{0}$) és a referencia-óra ($\mathbf{P}_{\text{ref}}$) közötti differenciális bezárt energia pedig $\Delta W(L) \equiv W_{\text{in}}(\mathbf{P}_{\text{main}}; L) - W_{\text{in}}(\mathbf{P}_{\text{ref}}; L)$. Mivel a 3D top-hat átlagolt mágneses tömegsűrűség-különbség $\Delta\bar{\rho}_{B,L} = \frac{I_0^2\,\Delta W(L)}{c^2\,V_L}$ (ahol $V_L = \frac{4}{3}\pi L^3$), a feltételes H1 óraeltolódási érzékenységi együttható ($k_g(L)$, $[\text{A}^{-2}]$) és a $I_0 = 14\;\text{mA}$-es jel (`make_revision3_checks.py`, 55–57. sor; `revision3_checks.json`, 139–145. sor):
$$\boxed{\;k_g(L) \;=\; -\frac{G L^2}{c^2}\frac{\Delta W(L)}{c^2\,\frac{4}{3}\pi L^3} \;=\; -\frac{3 G\,\Delta W(L)}{4\pi c^4 L}, \qquad \Delta y(L) \;=\; k_g(L)\,I_0^2.\;}$$
A 25 pontos logaritmikus rácson a csúcspontban ($L_p = 10^{-0{,}2}\;\text{m} = 0{,}6309573444801932\;\text{m}$):
$$\Delta W(L_p) \;=\; 1{,}9611445207162488\;\text{J/A}^2 \quad\Longrightarrow\quad k_g(L_p) \;=\; -6{,}131177672606294\times 10^{-45}\;\text{A}^{-2},$$
$$\Delta y(L_p) \;=\; k_g(L_p)\times (0{,}014\;\text{A})^2 \;=\; -1{,}2017108238308338\times 10^{-48}.$$

### 3.2. Miért előnyös a külső komplementer integrál nagy $L$-nél, és hogyan működik az $r = L/t$ kompaktifikáció?
Amikor $L > L_{\text{jump}} = 0{,}95\;\text{m}$ olyan nagy, hogy mind a fő óra $B(\mathbf{P}_{\text{main}}, L)$ gömbje, mind a referencia-óra $B(\mathbf{P}_{\text{ref}}, L)$ gömbje tartalmazza az intenzív tekercsközel-zónát, a két belső integrál közvetlen kivonása ($W_{\text{in,main}}(L) - W_{\text{in,ref}}(L)$) két szinte egyenlő szám különbsége lenne! Euklideszi térben azonban az egész $\mathbb{R}^3$-ra integrált mágneses energia **eltolás-invariáns**:
$$W_{\text{tot}} \;\equiv\; \frac{1}{2\mu_0 I_0^2}\int_{\mathbb{R}^3}|\mathbf{B}(\mathbf{x}; I_0)|^2\,d^3x \;=\; W_{\text{in}}(\mathbf{P}_{\text{main}}; L) + W_{\text{ext}}(\mathbf{P}_{\text{main}}; L) \;=\; W_{\text{in}}(\mathbf{P}_{\text{ref}}; L) + W_{\text{ext}}(\mathbf{P}_{\text{ref}}; L),$$
ahol $W_{\text{ext}}(\mathbf{P}; L) \equiv \frac{1}{2\mu_0 I_0^2}\int_{|\mathbf{x}-\mathbf{P}| \ge L}|\mathbf{B}(\mathbf{x}; I_0)|^2\,d^3x$ a gömbön kívüli (komplementer) energia. Mivel az egzakt összenergia kiesik ($W_{\text{tot,main}}^{\text{exact}} - W_{\text{tot,ref}}^{\text{exact}} \equiv 0$), a belső energiakülönbség pontosan egyenlő a sima, szingularitásmentes külső integrálok negált különbségével:
$$\boxed{\;\Delta W_{\text{out}}(L) \;\equiv\; -\Bigl[W_{\text{ext}}(\mathbf{P}_{\text{main}}; L) \;-\; W_{\text{ext}}(\mathbf{P}_{\text{ref}}; L)\Bigr] \;=\; \frac{1}{2\mu_0 I_0^2}\int_{r \ge L} r^2\,dr\int_{S^2} d\Omega\,\Bigl[|\mathbf{B}(\mathbf{P}_{\text{ref}} + r\hat{\mathbf{n}})|^2 - |\mathbf{B}(\mathbf{P}_{\text{main}} + r\hat{\mathbf{n}})|^2\Bigr].\;}$$

A V20 `energy_integrator.py` két független úton számítja ki ezt a $[L, \infty)$ félegyenesre vett külső integrált (`make_revision3_checks.py`, 62–68. sor; `tal_book_HU_2026-10-24.tex`, 581. sor):
- **„B” út (`dW_B`): Egzakt $r = L / t$ kompaktifikáció véges levágás nélkül.**
  Bevezetve a dimenziótlan $t \equiv L / r \in (0, 1]$ változót, a differenciál $dr = -\frac{L}{t^2}\,dt$, a térfogati mérték pedig $r^2\,dr = -\frac{L^3}{t^4}\,dt$. Így a végtelen $[L, \infty)$ tartomány egzaktul a kompakt $(0, 1]$ intervallumra képeződik le:
  $$\boxed{\;\Delta W_{\text{B}}(L) \;=\; \frac{L^3}{2\mu_0 I_0^2}\int_0^1 \frac{dt}{t^4}\int_{S^2} d\Omega\,\Bigl[\bigl|\mathbf{B}(\mathbf{P}_{\text{ref}} + \tfrac{L}{t}\hat{\mathbf{n}})\bigr|^2 \;-\; \bigl|\mathbf{B}(\mathbf{P}_{\text{main}} + \tfrac{L}{t}\hat{\mathbf{n}})\bigr|^2\Bigr].\;}$$
  Vizsgáljuk meg az integráns viselkedését a $t \to 0^+$ ($r = L/t \to \infty$) határesetben! Mivel a 8 zárt áramhurok össztöltése nulla, a távoli mágneses tér vezető rendben dipólustér: $|\mathbf{B}(r\hat{\mathbf{n}})| = \mathcal{O}(r^{-3}) = \mathcal{O}(t^3/L^3)$, vagyis az energiasűrűség $|\mathbf{B}|^2 = \mathcal{O}(t^6/L^6)$, a két eltolt gömb különbsége pedig $|\mathbf{B}(\mathbf{P}_{\text{ref}} + \frac{L}{t}\hat{\mathbf{n}})|^2 - |\mathbf{B}(\mathbf{P}_{\text{main}} + \frac{L}{t}\hat{\mathbf{n}})|^2 = \mathcal{O}(t^7/L^7)$. A $\frac{L^3}{t^4}$ Jacobi-faktorral megszorozva az integráns $t \to 0^+$-nál **reguláris polinomként ($\mathcal{O}(t^2)$–$\mathcal{O}(t^3)$)** tűnik el! Ezért az $n_t = 32$ (vagy $64$) pontos Gauss–Legendre kvadratúra a $(0,1]$ intervallumon gépi pontosságot ér el, és egyúttal bebizonyítja a távoli térbeli skálázást:
  $$\Delta W(L) \;\propto\; L^{-3}\;\text{(kvadrupólus/eltolási szimmetria esetén } L^{-5}\text{)} \quad\Longrightarrow\quad |\Delta y(L)| \;\propto\; \frac{\Delta W(L)}{L} \;\propto\; L^{-4}\text{--}L^{-6}.$$
- **„A” út (`dW_A`): Panelintegrál + analitikus végtelen multipólus-farok.**
  $L = 3{,}0\;\text{m}$-nél az „A” út (`dW_A = 0,0023313149428584538 J/A^2`) és a „B” út (`dW_B = 0,0023313149428584555 J/A^2`) relatív eltérése $-7{,}44\times 10^{-16}$ (`revision3_checks.json`, 164–165. sor; a két determinisztikus útvonal $3{,}7\times 10^{-16}$-os egyezése `tal_book_HU_2026-10-24.tex`, 478. sor).

### 3.3. Az $L_{\text{jump}} = 0{,}95\;\text{m}$-es belső/külső átkapcsolási ugrás egzakt algebrai azonossága (`make_revision3_checks.py`, 70–74. sor)
Ha a belső és a külső kvadratúrát véges $n_{\text{hi}} = 24000$ szögpontszámmal kiértékelve összeadjuk, a két koordináta-középpont körüli numerikus összenergia (`revision3_checks.json`, 186–195. sor):
$$W_{\text{tot,main}}^{\text{num}} \;=\; W_{\text{in,main}}^{\text{num}}(L) + W_{\text{ext,main}}^{\text{num}}(L) \;=\; 2{,}237861258956951\;\text{J/A}^2,$$
$$W_{\text{tot,ref}}^{\text{num}} \;=\; W_{\text{in,ref}}^{\text{num}}(L) + W_{\text{ext,ref}}^{\text{num}}(L) \;=\; 2{,}2364478629240048\;\text{J/A}^2.$$
Kivonva a két egyenletet egymásból, minden $L$-re (így az $L_{\text{jump}} = 0{,}95\;\text{m}$ átkapcsolási pontban is) fennáll az egzakt algebrai azonosság:
$$\boxed{\;\Delta W_{\text{in}}^{\text{num}}(0{,}95\;\text{m}) \;-\; \Delta W_{\text{out}}^{\text{num}}(0{,}95\;\text{m}) \;\equiv\; W_{\text{tot,main}}^{\text{num}} \;-\; W_{\text{tot,ref}}^{\text{num}} \;=\; 0{,}001413396032946\;\text{J/A}^2\;}$$
($2{,}17738947 - 2{,}17597607 = 0{,}00141340\;\text{J/A}^2$, a belső érték $0{,}06491\%$-a; az algebrai különbség-ellenőrzés maradéka `identity_check_abs` $= 8{,}88\times 10^{-16}\;\text{J/A}^2$).

---

## 4. Ellenőrzések (Checks)

1. **Az egyenáramú tekercsellenállás konzisztenciája a mért $D_{\text{menet}} = 291{,}7\;\text{mm}$ átmérővel (`tal_book_HU_2026-10-24.tex`, 577. sor):**
   - $N = 2000$ menet $d_{\text{wire}} = 0{,}15\;\text{mm}$ átmérőjű rézhuzalból ($A_{\text{wire}} = \frac{\pi}{4}d_{\text{wire}}^2 = 1{,}767146\times 10^{-8}\;\text{m}^2$, $\rho_{\text{Cu}}(20^\circ\text{C}) = 1{,}7241\times 10^{-8}\;\Omega\,\text{m}$) és $D_{\text{menet}} = 0{,}2917\;\text{m}$ középátmérővel: $L_{\text{wire}} = \pi D_{\text{menet}} N = 1832{,}805\;\text{m}$, amiből $R_{\text{DC,calc}} = 1788{,}160\;\Omega$, amely mindössze $-0{,}0470\%$-kal tér el a névleges $1789\;\Omega$-tól $\checkmark$.
2. **A csúcsjel ellenőrzése az (1) képlettel $L_p = 0{,}63095734448\;\text{m}$-nél:**
   - $k_g(L_p) = -\frac{3\times 6{,}67430\times 10^{-11}\times 1{,}961144520716}{4\pi\times (299792458)^4\times 0{,}63095734448} = -6{,}1311776726\times 10^{-45}\;\text{A}^{-2}$, és $\Delta y(L_p) = k_g(L_p)\times 0{,}014^2 = -1{,}2017108238\times 10^{-48}$, pontosan reprodukálva a V20/V22 fő eredményét $\checkmark$.
