# 7. lépés: Kompakt objektumok, görbült téridőbeli skalár–tenzor / Yukawa TOV-egyenletek (lapos Robin-közelítés vs. görbült külső ADM-tömeg) és a statikus nemlokális érintőoperátor stabilitási konzisztenciája (V21 / V23)
# Step 7: Compact Objects, Curved-Spacetime Scalar-Tensor / Yukawa TOV Equations (Flat-Robin Approximation vs. Curved Exterior ADM Mass), and Static Nonlocal Tangent Stability Consistency (V21 / V23)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step7_neutron_star_tov_stability.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step6_sparc_galactic_ngc3198_consistency.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p2.tex` (7. fejezet, 102–124. sor, Eq. 10) és `hu_p3.tex` (8. fejezet, 1–28. sor, Eq. 11–12; DOI: `10.5281/zenodo.23604423`)  
- V21 Zenodo rekord: DOI `10.5281/zenodo.23602501` (Neutroncsillag feltételes TOV érzékenységi és radiális pulzációs vizsgálat)  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.1–1.2. és 2. szakasz, 37–159. sor, Eq. 4–9, 1–2. táblázat; DOI: `10.5281/zenodo.23639084`)  
- Görbült külső `yukawa_action_el` ADM-tömeg és csatolt stabilitás: `curved_yukawa_adm_and_stability.py` és `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`

---

## 1. Szó szerinti idézetek a V22 monográfiából és a V23 (`50aa2db0`) kéziratból

Ebben a lépésben a **kompakt objektumok és a neutroncsillag TOV-stabilitás (V21 / V23)** egyenleteit vezetjük le, szigorúan szétválasztva:
1. a **Branch I-A** (V21 nemlokális 3D top-hat `tophat_mass` és `tophat_gradient`) háttérre felírt V23 Eq. (9) nemlokális statikus érintőoperátort, valamint
2. a **Branch II-Fluid** (V23 lokális `yukawa_action_el`) csatolt Einstein–skalár–folyadék rendszert, amelynél külön összehasonlítjuk a V23 kézirat 1. táblázatában alkalmazott **lapos téri Robin-közelítést és felszíni tömeget** ($M_{\text{surf}} = m(R_*)$) a **csatolt görbült külső vákuummegoldásból ($r \in [R_*, \infty)$) meghatározott egzakt Dirichlet–Neumann-meredekséggel és aszimptotikus ADM-tömeggel** ($M_{\text{ADM}} = \lim_{r\to\infty} m(r) = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$, részletes levezetés a 22. lépésben, `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`).

1. **A V21-beli két feltételes TOV-mód, az $1{,}214\times 10^{-6}$ mag-pontosság, a 769 csillagmodell (benne a 33 oksági végpont) és az $e^\nu F$ módszertani kikötése (`hu_p2.tex`, 7. fejezet, 107–124. sor):**
   > *„Auditált kódunkban az analitikus héjmag és a független 3D kvadratúra közötti mért maximális relatív eltérés az összes profilon **$1{,}214\times 10^{-6}$**. Ezután két feltételes Tolman–Oppenheimer–Volkoff (TOV) lezárást integráltam:”*
   $$\begin{equation}
     \textbf{(A) \texttt{mass} mód:}\quad \frac{dm_{\text{eff}}}{dr} = 4\pi r^2 \rho(r)\,f_L(r), \qquad
     \textbf{(B) \texttt{gradient} mód:}\quad \frac{d\Phi_{\text{eff}}}{dr} = f_L(r)\,\frac{G\bigl(m(r) + 4\pi r^3 P(r)/c^2\bigr)}{c^2 r\bigl(r - 2Gm(r)/c^2\bigr)}, \tag{V22 Eq. 10}
   \end{equation}$$
   > *„33 egyensúlyi sorozatot integráltam három nukleáris állapotegyenletre — **SLy** (nukleonikus), **APR4** (variációs nukleonikus; az oksági $c_s = c$ hangsebességhatárt **$\rho_0 \approx 1{,}4458\times 10^{15}\,\text{g/cm}^3$** nyugalmi tömegsűrűségnél éri el) és **H4** (merev hiperonikus) — 11 logaritmikusan elosztott $L \in [1\,\text{m}, 100\,\text{km}]$ simítási skálán. A V21 archív alaptáblázata **182 egyensúlyi modellt** tartalmaz, míg a sűrű diagnosztikai stabilitási rács **összesen 769 csillagmodellből** áll (beleértve mind a **33 egzakt oksági végpontot**, ahol $c_s/c = 1$).”*
   > *„A TOV-egyenletek értékelésekor kifejezetten rögzítem: **a módosított $e^\nu F$ (vagy $\Phi_{\text{eff}}'$) hidrosztatikai faktor megjelenése a feltételes nyomásgradiens-egyenletben egy algebrai egyensúlyi feltevés, és önmagában nem vezeti le, hogyan csatolódnak az atomórák, az atommag-átmenetek vagy a felszíni fotonok vöröseltolódásai az $f_L$ faktorhoz.**”*

2. **Az algebrai X-szektorbeli feszültségrekonstrukció és a befagyasztott-$f_L$ Chandrasekhar-operátorok korlátai (`hu_p3.tex`, 8. fejezet, 8–22. sor és `tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 48–52. sor):**
   > *„A 769 auditált neutroncsillag-modellen a relatív tangenciális feszültség korlátos marad, maximuma **$\max |p_{tX}/P_c| = 0{,}140$** (a centrális nyomás $14{,}0\%$-a). Hangsúlyozom, hogy ez az X-szektorbeli rekonstrukció a háttérbeli Bianchi-azonosság algebrai könyvelési lezárása, még nem hatáselvből levezetett dinamikus mezőegyenlet.”*
   $$\begin{align}
     \textbf{Tömeg-mód } (m' = 4\pi r^2 f_L\varepsilon,\; \nu' = A_{\text{geo}}): &\quad \varepsilon_X = -(1-f_L)\varepsilon,\quad P_{r,X} = 0,\quad P_{t,X} = -\frac{r}{2}(1-f_L)\varepsilon\nu', \tag{V23 Eq. 4}\\
     \textbf{Gradiens-mód } (m' = 4\pi r^2\varepsilon,\; \nu' = f_L A_{\text{geo}}): &\quad \varepsilon_X = 0,\quad P_{r,X} = -(1-f_L)\!\left(P + \frac{m}{4\pi r^3}\right),\; P_{t,X} = P_{r,X} + \frac{r}{2}\bigl(P'_{r,X} + P_{r,X}\nu'\bigr). \tag{V23 Eq. 5}
   \end{align}$$

3. **A V23 görbült téridőbeli Yukawa hatás-EL (`yukawa_action_el`) egyenletei (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 65. sor):**
   > *„A `yukawa_action_el` sorok az $\mathcal{L}=-[(\nabla\phi)^2+\phi^2/L^2]/(2\alpha_Y)-f_L\varepsilon$ hatásból következnek ($G=c=1$; `SymPy`: `sym_yukawa_el_source_factor` és `_hydrostatic` maradék $=0$): $B\phi''+(2B/r+B\nu'+B'/2)\phi'-\phi/L^2=-\alpha_Y f_L\varepsilon$, $B=1-2m/r$, $P_{\text{matter}}=f_LP$, $h'=-\nu'+\phi'$, $\rho_\phi=(B\phi'^2+\phi^2/L^2)/(2\alpha_Y)$, külső illesztés $\phi'(R)=-(1/R+1/L)\phi(R)$; a $\Delta M_{\phi,\text{ext}}$ külső farok külön szerepel (SLy, $L=10\text{ km}$: $3{,}4\times 10^{-3}\,M_\odot$ $\log_{10}\rho_{0c}=15$-nél, $9{,}0\times 10^{-3}\,M_\odot$ $M_{\max}$-nál).”*

4. **A befagyasztott-$f_L$ Sturm–Liouville maradék és a nemlokális statikus érintőoperátor (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 2. szakasz, 122–135. sor, Eq. 8–9):**
   $$\begin{equation}
     \mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} = \frac{4\pi r^{5/2}(f_L - 1)\Gamma P \varepsilon e^{2\nu}}{(r - 2m)^{5/2}}\Big[(4\pi r^3 P - 3r + 7m)\xi - r(r-2m)\xi'\Big]. \tag{V23 Eq. 8}
   \end{equation}$$
   $$\begin{equation}
     \delta f_L(r_i) = -\gamma_L f_L(r_i)\sum_j K_{ij}\,\frac{\varepsilon(r_j)+P(r_j)}{c_s^2(r_j)}\,\delta h(r_j), \qquad
     \begin{cases}
       \dfrac{d(\delta h)}{dr} = -\delta\nu'\!\left[\delta h, \delta m, \delta f_L[\delta h]\right] + \omega_0^2\,e^{\lambda-2\nu}\xi, \\[4pt]
       \dfrac{d(\delta m)}{dr} = 4\pi r^2\!\left[f_{\text{mass}}\frac{\varepsilon+P}{c_s^2}\delta h + \varepsilon\,\delta f_{\text{mass}}[\delta h]\right].
     \end{cases} \tag{V23 Eq. 9}
   \end{equation}$$

---

## 2. A görbült téridőbeli skalár–tenzor / Yukawa TOV-rendszer (`yukawa_action_el`) levezetése

### 2.1. A görbült téridőbeli d'Alembert-operátor és a Yukawa mezőegyenlet
Sztatikus, gömbszimmetrikus téridőben ($G = c = 1$ geometrizált egységekben):
$$ds^2 \;=\; -e^{2\nu(r)}dt^2 + B(r)^{-1}dr^2 + r^2(d\theta^2 + \sin^2\theta\,d\varphi^2), \qquad B(r) \;\equiv\; e^{-2\lambda(r)} \;=\; 1 - \frac{2m(r)}{r},$$
a metrika determinánsa $\sqrt{-g} = e^{\nu(r)} B(r)^{-1/2} r^2 \sin\theta$. Egy csak $r$-től függő $\phi_L(r)$ skalármezőre a kovariáns Laplace–Beltrami / d'Alembert-operátor:
$$\Box_g \phi_L \;=\; \frac{1}{\sqrt{-g}}\partial_r\!\left(\sqrt{-g}\,g^{rr}\,\partial_r\phi_L\right) \;=\; \frac{B^{1/2}}{e^\nu r^2}\frac{d}{dr}\!\left(e^\nu B^{1/2} r^2 \phi_L'\right) \;=\; B\phi_L'' + \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L'.$$
Az $\mathcal{L}_{\phi+m} = -\frac{1}{2\alpha_Y}\left[g^{\mu\nu}\nabla_\mu\phi_L\nabla_\nu\phi_L + \frac{\phi_L^2}{L^2}\right] - f_L(\phi_L)\varepsilon(n)$ ($f_L = e^{-\phi_L}$, $\alpha_Y = 4\pi\alpha_{\text{std}} = 8\pi\tilde{\alpha}_L = 2\beta^2\cdot 4\pi$) Lagrangiánus $\phi_L$ szerinti Euler–Lagrange-egyenlete ($\frac{1}{\sqrt{-g}}\partial_\mu(\sqrt{-g}\frac{\partial\mathcal{L}}{\partial(\partial_\mu\phi_L)}) - \frac{\partial\mathcal{L}}{\partial\phi_L} = 0$, figyelembe véve, hogy $\frac{\partial(-e^{-\phi_L}\varepsilon)}{\partial\phi_L} = +f_L\varepsilon$):
$$\boxed{\;B\phi_L'' \;+\; \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L' \;-\; \frac{\phi_L}{L^2} \;=\; -\alpha_Y\,f_L(\phi_L)\,\varepsilon(r).\;}$$
Ez pontosan a V23 1.2. szakaszában (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 65. sor) szereplő `yukawa_action_el` egyenlet, ahol a forrástag a dilatonszerű $f_L = e^{-\phi_L}$ önkonzisztens súlyfaktorral szorzódik (`sym_yukawa_el_source_factor = 0`).

### 2.2. A csatolt Einstein–Yukawa–Euler (TOV) egyenletrendszer, a görbült külső peremfeltétel és az ADM-tömeg ($M_{\text{ADM}}$ vs. $M_{\text{surf}}$)
A 4. lépésben (`step4_covariant_conservation_hydrostatic.md`) levezetett $T^\mu_{\;\nu}(\text{matter}) = \text{diag}(-f_L\varepsilon, f_L P, f_L P, f_L P)$ és $T^\mu_{\;\nu}(\phi)$ tenzorokból a $(tt)$ és $(rr)$ Einstein-egyenletek, valamint a $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$ megmaradási törvény adják a csatolt, hatáselvből zárt TOV-rendszert:
$$\begin{aligned}
  \rho_\phi(r) &\;=\; -T^t_{\;t}(\phi) \;=\; \frac{1}{2\alpha_Y}\!\left(B(r)\bigl[\phi_L'(r)\bigr]^2 + \frac{\phi_L(r)^2}{L^2}\right), \qquad
  p_{r,\phi}(r) \;=\; T^r_{\;r}(\phi) \;=\; \frac{1}{2\alpha_Y}\!\left(B(r)\bigl[\phi_L'(r)\bigr]^2 - \frac{\phi_L(r)^2}{L^2}\right), \\[4pt]
  \frac{dm(r)}{dr} &\;=\; 4\pi r^2\bigl[f_L(r)\varepsilon(r) + \rho_\phi(r)\bigr], \\[4pt]
  \frac{d\nu(r)}{dr} &\;=\; \frac{m(r) + 4\pi r^3\bigl[f_L(r)P(r) + p_{r,\phi}(r)\bigr]}{r\bigl(r - 2m(r)\bigr)}, \\[4pt]
  \frac{dh(r)}{dr} &\;=\; \frac{P'(r)}{\varepsilon(r) + P(r)} \;=\; -\frac{d\nu(r)}{dr} \;+\; \frac{d\phi_L(r)}{dr}.
\end{aligned}$$

**Kritikus módszertani tisztázás (Point 1 és Point 2 — részletes levezetés a 22. lépésben):**
1. **A lapos téri Robin-közelítés és a valódi görbült külső Dirichlet–Neumann-meredekség különbsége:**
   - A V23 kézirat (`covariant_closure_50aa2db0.py` eredeti `yukawa_action_el` sora) a csillag $r = R_*$ felszínén a lapos téri Helmholtz-megoldás ($\phi_L^{\text{flat}}(r) = \phi_L(R_*)\frac{R_*}{r}e^{-(r-R_*)/L}$) deriváltját használta peremfeltételként:
     $$\phi_{L,\text{flat}}'(R_*) \;=\; -k_{\text{flat}}\,\phi_L(R_*), \qquad k_{\text{flat}} \;=\; \frac{1}{R_*} + \frac{1}{L}.$$
   - A valóságban azonban a csillag felszínén kívüli vákuumban ($r \ge R_*$, ahol $\varepsilon = P = 0$) a téridő erősen görbült ($B(R_*) = 1 - 2M_{\text{surf}}/R_* \approx 0{,}50\text{--}0{,}62$). A külső vákuumban az egzakt azonosság $B\nu' + B'/2 = \frac{2m(r)}{r^2} - \frac{4\pi r}{\alpha_Y L^2}\phi_L(r)^2$, így a külső csatolt Einstein–skalár egyenletrendszer:
     $$\left(1 - \frac{2m(r)}{r}\right)\phi_L''(r) + \left[\frac{2}{r}\left(1 - \frac{m(r)}{r}\right) - \frac{4\pi r}{\alpha_Y L^2}\phi_L(r)^2\right]\phi_L'(r) - \frac{\phi_L(r)}{L^2} \;=\; 0, \qquad m'(r) \;=\; 4\pi r^2\rho_\phi(r).$$
     Mivel $B(R_*) < 1$, a $\phi_L/L^2$ tömegtag effektív lokális csillapítási rátája a felszínnél $\sim B(R_*)^{-1/2}/L > 1/L$. Ahogy a 22. lépésben (`curved_yukawa_adm_and_stability.py`) bebizonyítjuk, SLy $\log_{10}\rho_{0c} = 15{,}0$, $L = 10\;\text{km}$, $\alpha_Y = 1$ esetén a valódi görbült Dirichlet–Neumann-meredekség $k_{\text{curved}} = -\phi_{L,\text{curved}}'(R_*)/\phi_L(R_*) = 0{,}233542\;\text{km}^{-1}$, ami **$+24{,}74\%$-kal nagyobb**, mint a lapos Robin-közelítés $k_{\text{flat}} = 0{,}187216\;\text{km}^{-1}$ értéke!
2. **Felszíni tömeg ($M_{\text{surf}}$) vs. aszimptotikus ADM-tömeg ($M_{\text{ADM}}$):**
   - A V23 kézirat 1. táblázatában a `yukawa_action_el` sor $M_{\max}$ és $R_{1{,}4}$ oszlopai a $m(R_*) = M_{\text{surf}}$ felszíni tömegre voltak optimalizálva, míg a $\Delta M_{\phi,\text{ext}}$ külső skalártömeg külön oszlopban szerepelt.
   - Egy távoli megfigyelő (és a bardeen–carter–hawking-i / Harrison–Wheeler-i fordulóponti tétel) számára a fizikai gravitációs tömeg az aszimptotikus **ADM-tömeg**:
     $$\boxed{\;M_{\text{ADM}} \;\equiv\; \lim_{r \to \infty} m(r) \;=\; M_{\text{surf}} \;+\; \Delta M_{\phi,\text{ext}}^{\text{curved}}, \qquad \Delta M_{\phi,\text{ext}}^{\text{curved}} \;=\; \int_{R_*}^\infty 4\pi r^2 \rho_\phi(r)\,dr.\;}$$
   - A görbült külső megoldással újraszámolt `yukawa_action_el_curved_adm` sorozatokon (SLy, $L = 10\;\text{km}$, $\alpha_Y = 1$):
     - $M_{\text{ADM},\max} = \mathbf{2{,}03321\,M_\odot}$ ($\log_{10}\rho_{0c}(M_{\text{ADM},\max}) = \mathbf{15{,}316404}$, $R(M_{\text{ADM},\max}) = 9{,}8853\;\text{km}$, $\Delta M_{\phi,\text{ext}}^{\text{curved}} = 0{,}00411\,M_\odot$, $M_{\text{surf}} = 2{,}02911\,M_\odot$),
     - a konzervált barionos nyugalmi tömeg ($M_B = \int_0^{R_*} 4\pi r^2 B(r)^{-1/2}\rho_0(r)\,dr$) maximuma $\log_{10}\rho_{0c}(M_{B,\max}^{\text{quad}}) = \mathbf{15{,}316409}$ ($|\Delta\log_{10}\rho_{0c}| = 4{,}47\times 10^{-6}\;\text{dex}$ egyezés az ADM-tömeg fordulópontjával!),
     - $R(M_{\text{ADM}} = 1{,}4\,M_\odot) = \mathbf{11{,}5359\;\text{km}}$ (míg $R(M_{\text{surf}} = 1{,}4\,M_\odot) = 11{,}5348\;\text{km}$),
     - és a `yukawa_action_el` háttéren megoldott csatolt $(\delta h, \delta\phi_L)$ dinamikus radiális pulzációs sajátérték zérushelye (`scalar_m_sign = -1.0` fő operátorral) $\log_{10}\rho_{0c}(\omega_0^2 = 0) = \mathbf{15{,}316410}$ ($\mathbf{+5{,}80\times 10^{-6}\;\text{dex}}$ eltéréssel a $dM_{\text{ADM}}/d\rho_c = 0$ fordulóponttól és $+1{,}33\times 10^{-6}\;\text{dex}$-szel a $dM_B/d\rho_c = 0$ fordulóponttól; $L=1\;\text{km}$-nél $\mathbf{15{,}303080}$, $-6{,}86\times 10^{-6}\;\text{dex}$; lásd a 12. és 22. lépést).

---

## 3. A V21 nemlokális 3D top-hat háttér (**Branch I-A**) stabilitási vizsgálata: befagyasztott-$f_L$ maradék vs. $\delta f_L[\delta h]$ érintőoperátor

*(Megjegyzés a modellágak szétválasztásához — Point 4: Az alábbi 3.1–3.2. szakaszok kifejezetten a V21 nemlokális 3D gömbi top-hat `tophat_mass` és `tophat_gradient` [**Branch I-A**] hátterekre vonatkoznak, míg a `yukawa_action_el` [**Branch II-Fluid**] csatolt stabilitását a 12. és 22. lépés tartalmazza.)*

### 3.1. A befagyasztott-$f_L$ Chandrasekhar-operátor algebrai maradéktagja (V23 Eq. 8)
Általános relativitáselméletben ($f_L \equiv 1$) a sztatikus egyensúlyi sorozat mentén vett $\partial/\partial\rho_c$ derivált egy nulla sajátfrekvenciájú ($\omega_0^2 = 0$) radiális próbaperturbációt állít elő a $dM/d\rho_c = 0$ extremális pontban. A klasszikus Chandrasekhar-féle másodrendű Sturm–Liouville-operátor ($\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi]$) levezetése azonban parciális integrálás közben behelyettesíti a GR-beli
$$\lambda' + \nu' \;=\; 4\pi r e^{2\lambda} (\varepsilon + P) \qquad (\text{csak } f_L \equiv 1 \text{ esetén!})$$
háttérazonosságot. Ha a háttér a V21 `mass` mód szerint módosul, akkor a Chandrasekhar-operátor és a $\delta f_L = 0$ mellett vett linearizált hidrosztatikai operátor ($\mathcal{H}_{\text{lin}}[\xi]\big|_{\delta f_L=0}$) különbsége nem tűnik el, hanem pontosan a V23 (8) egyenletében megadott $(f_L - 1)$-gyel arányos maradékot adja:
$$\mathcal{L}_{\text{SL}}^{\text{matter}}[\xi] - \sqrt{\frac{r}{r-2m}}\,e^{2\nu}\,\mathcal{H}_{\text{lin}}[\xi]\Big|_{\delta f_L=0} \;=\; \frac{4\pi r^{5/2}(f_L - 1)\Gamma P \varepsilon e^{2\nu}}{(r - 2m)^{5/2}}\Big[(4\pi r^3 P - 3r + 7m)\xi - r(r-2m)\xi'\Big].$$

### 3.2. A nemlokális $\delta f_L[\delta h]$ érintőperturbáció beépítése a V21 háttéren (V23 Eq. 9 és 2. táblázat)
Amikor a központi sűrűséget $\rho_{0c} \to \rho_{0c} + \delta\rho_{0c}$-vel megváltoztatjuk a V21 3D top-hat háttéren, az energiasűrűség Euler-perturbációja $\delta\varepsilon(r) = \frac{\varepsilon(r)+P(r)}{c_s^2(r)}\,\delta h(r)$, ami a nemlokális magon keresztül megváltoztatja az $f_L$ faktort is:
$$\delta f_L(r_i) \;=\; -\gamma_L\,f_L(r_i)\,\sum_j K_{ij}\,\frac{\varepsilon(r_j)+P(r_j)}{c_s^2(r_j)}\,\delta h(r_j).$$
Ezt a $\delta f_L[\delta h]$ operátort visszaírva a V21 háttér körüli linearizált $(\delta h, \delta m)$ peremérték-feladatba (V23 Eq. 9), $\omega_0^2 = 0$-nál a rendszer egzaktul megegyezik a V21 statikus egyensúlyi egyenletek $\log_{10}\rho_{0c}$ szerinti érintőoperátorával. Ennek eredményeként az SLy ($L = 10\text{ km}$) V21 sorozaton a $\omega_0^2 = 0$ zérushely és a $dM/d\rho_c = 0$ statikus fordulópont eltérése (V23 2. táblázat):
- **GR alapvonal ($L = 0\text{ km}$):** $\Delta\log_{10}\rho_{0c} = +7{,}12\times 10^{-5}\;\text{dex}$ (diszkretizációs alapérték),
- **V21 SLy `mass` mód ($L = 10\text{ km}$, fordulópont $\log_{10}\rho_{0c} = 15{,}29924$):** a befagyasztott-$f_L$ BVP $-4{,}06\times 10^{-3}\;\text{dex}$ (és a `matter` SL $-1{,}34\times 10^{-3}\;\text{dex}$) eltéréséről **$\mathbf{+6{,}97\times 10^{-5}\;\text{dex}}$**-re csökken (**$58{,}3\times$-os**, illetve **$19{,}2\times$-es** javulás),
- **V21 SLy `gradient` mód ($L = 10\text{ km}$, fordulópont $\log_{10}\rho_{0c} = 15{,}30106$):** a befagyasztott-$f_L$ BVP $-4{,}39\times 10^{-3}\;\text{dex}$ (és a `matter` SL $-3{,}57\times 10^{-3}\;\text{dex}$) eltéréséről **$\mathbf{+7{,}39\times 10^{-5}\;\text{dex}}$**-re csökken (**$59{,}3\times$-os**, illetve **$48{,}3\times$-os** javulás).

---

## 4. Ellenőrzések (Checks)

1. **$f_L \to 1$ ($L \to 0$ vagy $\alpha_Y \to 0$) GR-határeset:**
   - $\phi_L \to 0 \implies \rho_\phi \to 0,\; p_{r,\phi} \to 0,\; f_L \to 1$, a `yukawa_action_el` TOV-rendszer visszaadja a standard GR TOV-egyenleteket ($M_{\max}^{\text{SLy}} = 2{,}04841\,M_\odot$, $M_{\max}^{\text{APR4}} = 2{,}18800\,M_\odot$, $M_{\max}^{\text{H4}} = 2{,}01385\,M_\odot$), a V23 Eq. (8) maradéktag pedig $(f_L - 1) = 0$ miatt azonosan eltűnik $\checkmark$.
2. **Görbült külső `yukawa_action_el` ADM-tömeg és csatolt stabilitási fordulópont egyezése (lásd `step22`):**
   - SLy $L = 10\;\text{km}$, $\alpha_Y = 1$: $\log_{10}\rho_{0c}(dM_{\text{ADM}}/d\rho_c = 0) = 15{,}316404$, $\log_{10}\rho_{0c}(dM_B/d\rho_c = 0) = 15{,}316409$, csatolt $(\delta h, \delta\phi_L)$ zérushely $\log_{10}\rho_{0c}(\omega_0^2 = 0) = 15{,}316410$ ($\Delta_{\text{ADM}} = +5{,}80\times 10^{-6}\;\text{dex}$, $2097\times$-es csökkenés a befagyasztott-$\phi_L$ $+1{,}22\times 10^{-2}\;\text{dex}$ eltolódásához képest; $L=1\;\text{km}$-nél $15{,}303080$ vs. $15{,}303086$, $-6{,}86\times 10^{-6}\;\text{dex}$, $165\times$-ös csökkenés) $\checkmark$.
3. **Tehetetlenség-skálázási invariancia a V21 zérushelyen ($\omega_0^2 = 0$):**
   - A $W \to 4W$ skálázás $\log_{10}\rho_{0c} = 15{,}30427$-nél pontosan $\omega^2_{4W}/\omega^2_W = 0{,}2500$ arányt ad ($\omega^2_W = 4{,}897\times 10^6\;\text{s}^{-2}$), míg az $\omega_0^2 = 0$ zérushelyen $F(0)$ invariáns, igazolva, hogy a V21-beli $+6{,}97\times 10^{-5}\;\text{dex}$ és $+7{,}39\times 10^{-5}\;\text{dex}$ egyezés a statikus nemlokális érintőoperátor egzakt konzisztenciájának következménye $\checkmark$.

> **V25 frissítés (görbült `yukawa_action_el` ág):** 4-komponensű külső BVP ($\delta\nu_R=1{,}63\times 10^{-14}$), 108-pontos több-EOS / több-sugárpont $4\times 4$ Jacobi-ellenőrzés ($\max\text{rel.\ hiba} = 2{,}74\times 10^{-5}$; SLy $15{,}0$, $L=10\;\text{km}$: $\partial\phi_L'(R_*)/\partial m(R_*)=-1{,}68365\times 10^{-4}\;\text{km}^{-2}$, rel. eltérés $3{,}98\times 10^{-10}$; belső $r=6\;\text{km}$ `a_m_sc` $=-1{,}72833\times 10^{-4}\;\text{km}^{-3}$, rel. eltérés $6{,}08\times 10^{-9}$), háromféle barionos tömeg ($M_B^{\text{quad}}, M_B^{\text{ODE}}, M_B^{\text{trapz}}$) és első főtétel-audit, valamint csomópontmentes ($n=0$) alapmódus és statikus érintőprojekció ($|\langle\boldsymbol{\Psi}_0,\boldsymbol{\Psi}_{\text{tan}}\rangle|\in[0{,}9963,\;0{,}99999993]$). SLy $L=10$ km, $\alpha_Y=1$: $M_{\text{ADM}}$-fordulópont $15{,}316404$, $M_B^{\text{quad}}$-fordulópont $15{,}316409$ ($M_B^{\text{trapz}}$: $15{,}316389$, $M_B^{\text{ODE}}$: $15{,}316522$), fő dinamikus $\omega_0^2=0$ gyök (`scalar_m_sign=-1.0`) $15{,}316410$ ($+5{,}80\times 10^{-6}$ dex; pluszjeles `scalar_m_sign=+1.0` diagnosztikával $15{,}319129$, $+2{,}72\times 10^{-3}$ dex). Részletes táblák: `step22`, 4. és 7. szakasz.
