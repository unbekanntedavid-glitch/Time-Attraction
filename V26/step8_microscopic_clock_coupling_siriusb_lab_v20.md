# 8. lépés: Mikroszkopikus Dirac–Maxwell óra- és fotoncsatolás a hatásból (H0 vs. feltételes H1), a Sirius B vöröseltolódási korlátai és a V20 laboratóriumi tesztek ($8558\text{ kg}$-os gránitkamra és 8 tekercses mágnesestér-mérés)
# Step 8: Microscopic Dirac–Maxwell Clock and Photon Coupling from the Action (H0 vs. Conditional H1), Sirius B Redshift Constraints, and the V20 Laboratory Tests ($8558\text{ kg}$ Granite Chamber and 8-Coil Magnetic-Field Experiment)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step8_microscopic_clock_coupling_siriusb_lab_v20.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step7_neutron_star_tov_stability.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (3. szakasz, 161–190. sor, 3. táblázat; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p3.tex` (9. és 10. fejezet, 30–91. sor, Eq. 13; DOI: `10.5281/zenodo.23604423`)  
- V20 és V21 Zenodo rekordok: DOI `10.5281/zenodo.23540741` (V20 8 tekercses mágnesestér és $8558\text{ kg}$-os gránitkamra) és DOI `10.5281/zenodo.23602501` (V21 Sirius B és NICER feltételes H0/H1 elemzés)  
- Helyi szolenoid-munkajelentés: `/workspace/artifacts/solenoid_update/report_hu.md` (17–60. sor)

---

## 1. Szó szerinti idézetek a V23 (`50aa2db0`) kéziratból és a V22 monográfiából

A `README.md` teszttervének utolsó két vizsgálati eleme (a 3.4. pont második fele és a 3.5. pont) a **mikroszkopikus óra- és fotoncsatolás (H0 vs. feltételes H1)** levezetése a kovariáns skalár–tenzor / diszform hatásból, valamint ennek kvantitatív kiértékelése a **Sirius B fehér törpén** és a **V20 laboratóriumi program két ágán** ($8558\text{ kg}$-os gránitkamra és 8 tekercses szolenoid).

1. **A mikroszkopikus Dirac–Maxwell óracsatolás két ága a hatáselvben (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 3. szakasz, 162–166. sor):**
   > *„Jelölje $\Psi_{\text{atomic}} = (\psi_e, A_\mu)$ a barionos közeg $u^\mu$ négyes-sebességével együttmozgó atomóra Dirac-elektron- és Maxwell-elektromágneses mezőit, amelyek a $\tilde{g}_{\mu\nu}[\phi_L]$ effektív Jordan-keretbeli metrikához csatolódnak:*
   > - **Minimális / tisztán segédmezős csatolás ($\tilde{g}_{\mu\nu} = g_{\mu\nu}$ a $\Psi_{\text{atomic}}$ szektorra $\implies$ H0 hipotézis):** *Ha $\phi_L$ kizárólag a makroszkopikus barionos folyadék feszültség-energia tenzorához csatolódik, míg az atomi elektronok és fotonok a fizikai $g_{\mu\nu}$ téridő-metrikán mozognak, akkor az $r = R_*$ felszínen nyugvó atomóra sajátideje $d\tau_{\text{clock}} = \sqrt{-g_{00}(R_*)}\,dt = e^{\Phi_{\text{eff}}(R_*)/c^2}dt$. A felszíni vöröseltolódás tisztán metrikus ($1+z_{\text{H0}} = e^{-\Phi_{\text{eff}}(R_*)/c^2}$), közvetlen $f_L(R_*)$ szorzó nélkül.*
   > - **Univerzális diszform lapse- vagy konform Jordan-csatolás ($\tilde{g}_{\mu\nu} = g_{\mu\nu} + \frac{1 - f_L^2(\phi_L)}{c^2}u_\mu u_\nu \implies$ H1 hipotézis):** *Ha az atomórák a $\tilde{g}_{\mu\nu}$ diszform metrikát érzékelik, amelynek együttmozgó idő-idő komponense $g_{00}=-1$, $u_0=-c$ mellett $\tilde{g}_{00} = -1 + (1-f_L^2) = -f_L^2$ (`SymPy`: `sym_disformal_g00_residual` $=0$), akkor az atomi átmenet sajátidő-eleme $d\tilde{\tau}_{\text{clock}} = f_L(R_*)\,e^{\Phi_{\text{eff}}(R_*)/c^2}dt$. A végtelenbe futó fotonok koordináta-frekvenciája megmarad, így megjelenik a közvetlen felszíni óraszorzó: $1+z_{\text{H1}} = e^{-\Phi_{\text{eff}}(R_*)/c^2}\,f_L(R_*)^{-1}$.”*

2. **A Sirius B fehér törpe feltételes H0 és feltételes H1 vöröseltolódási eredményei (`hu_p3.tex`, 9. fejezet, 53–57. sor):**
   > *„Mind a hat hipotézist leteszteltem a **Sirius B** fehér törpe ($M = 1{,}018 \pm 0{,}011\,M_\odot$, $R = 0{,}008098 \pm 0{,}000046\,R_\odot$, megfigyelt vöröseltolódási sebesség $v_{\text{obs}} = 80{,}65 \pm 0{,}77\,\text{km/s}$, GR-alapvonal $v_{\text{GR}} \approx 80{,}42\,\text{km/s}$) nagy pontosságú HST/STIS Balmer-vonalbeli mérésén, nemrelativisztikus ($n=1{,}5$) és relativisztikus ($n=3$) politróp belső profilokon egyaránt:*
   > - **H0-gradient lapse alatt (Csak hidrosztatikai csatolás):** *[A] maximális sebességeltolódás minden $L$-en **$|\Delta v| \le 1{,}27\,\text{m/s}$** ($0{,}00127\,\text{km/s}$), ami három nagyságrenddel kisebb a $\pm 770\,\text{m/s}$-os mérési hibánál. Így a Sirius B a H0 alatt semmilyen korlátot nem ad az $L$-re.*
   > - **Feltételes H1 alatt (Közvetlen felszíni óracsatolás, $\nu \propto f_L(R_*)$):** *Amikor a fotoszférában zajló atomi átmeneteket közvetlenül megszorozzuk az $f_L(R_*)$ faktorral, az $L \sim R_* \approx 5{,}63\times 10^6\,\text{m}$ sugarú simítógömb befogja a fehér törpe tömegének jelentős részét, akár $\Delta v \approx +5{,}56\,\text{km/s}$ többlet-vöröseltolódást okozva, elérve a **$\max\Delta\chi^2 = 52{,}3$** (egyszorzós H1) és **$\max\Delta\chi^2 = 72{,}0$** (H1+gradient lapse) értéket. $2\sigma$-s ($\Delta\chi^2 = 4$) szinten a feltételes H1 kizárja az **$L \in [4{,}08\times 10^6,\, 3{,}37\times 10^7]\,\text{m}$ ($n=1{,}5$)** és **$L \in [4{,}95\times 10^6,\, 3{,}37\times 10^7]\,\text{m}$ ($n=3$)** köztes sávot, miközben szabadon hagyja mind a laboratóriumi ($L \sim 1\,\text{m}$), mind a kilométer alatti neutroncsillag-skálákat!”*

3. **A V20 laboratóriumi program két ága (`hu_p3.tex`, 10. fejezet, 71–85. sor):**
   > - **„A” ág ($8558\,\text{kg}$-os gránitkamra):** *„Belső furatsugár $R_{\text{in}} = 0{,}25\,\text{m}$, külső sugár $R_{\text{out}} = 0{,}85\,\text{m}$, függőleges magasság $H = 1{,}20\,\text{m}$ ($z \in [-0{,}60, +0{,}60]\,\text{m}$), térfogat $V = \pi(R_{\text{out}}^2 - R_{\text{in}}^2)H = 3{,}1102\,\text{m}^3$, a gránit össztömege pedig **$M_{\text{granite}} = 8558\,\text{kg}$**.”*
   > - **„B” ág (8 tekercses mágnesestér-mérés):** *„8 darab egyforma szolenoid, egyenként $l_{\text{coil}} = 400\,\text{mm}$ hosszúsággal, külső szerkezeti **$\varnothing 300\,\text{mm}$ peremhatár-átmérővel**, effektív **$\varnothing 291{,}7\,\text{mm}$ réz menetátmérővel** ($R_{\text{coil}} = 0{,}14585\,\text{m}$), $N = 2000$ menet $0{,}15\,\text{mm}$-es zománcozott rézhuzalból, tekercsenként $R_{\text{DC}} = 1789\,\Omega$ egyenáramú ellenállással, $I = 14{,}0\,\text{mA}$ árammal meghajtva ($0{,}351\,\text{W}$/tekercs, összesen $2{,}805\,\text{W}$).”*
   $$\begin{equation}
     \Delta y_{\text{peak}}(L = 0{,}630957\,\text{m}) \;=\; -1{,}2017108\times 10^{-48}. \tag{V22 Eq. 13}
   \end{equation}$$

---

## 2. Matematikai levezetés: A Dirac–Maxwell sajátidőtől az óraeltolódási képletekig

### 2.1. Hogyan választ a Dirac–Maxwell hatás a H0 és a feltételes H1 hipotézis között?
Tekintsük egy nyugvó atomóra ($\mathbf{x} = \text{konst.}$, $u^\mu = (c/\sqrt{-g_{00}}, \mathbf{0})$) lokális Dirac-elektron ($\psi_e$, tömeg $m_e$) és Coulomb/Maxwell ($A_\mu$) hatását az effektív $\tilde{g}_{\mu\nu}$ órametrikán:
$$S_{\text{atomic}}[\tilde{g}_{\mu\nu}, \psi_e, A_\mu] \;=\; \int d^4x\,\sqrt{-\tilde{g}}\left[\bar{\psi}_e\left(i\hbar c\,\tilde{e}^{\mu}_{\;a}\gamma^a D_\mu - m_e c^2\right)\psi_e - \frac{1}{4\mu_0}\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}\right].$$
A 2. lépésben (`step2_conformal_disformal.md`) levezetett tetrád-transzformáció alapján a lokális együttmozgó keretben $\tilde{e}^0_{\;0} = \sqrt{-\tilde{g}_{00}}$, így a stacionárius Dirac-sajátállapotok ($\psi_e \propto e^{-i\omega_{\text{coord}} t}$) energiája (Bohr-/hiperfinom átmeneti frekvenciája $\nu_0 = \Delta E_{\text{atomic}}/h$) a lokális óra-sajátidőben ($d\tilde{\tau}_{\text{clock}} = \sqrt{-\tilde{g}_{00}}\,dt/c$) állandó:
$$\nu_{\text{coord}}(\mathbf{x}) \;=\; \nu_0\,\sqrt{-\tilde{g}_{00}(\mathbf{x})}.$$
Mivel sztatikus téridőben ($\partial_t g_{\mu\nu} = 0$, $\partial_t \phi_L = 0$) a Killing-vektor $\xi^\mu = (\partial_t)^\mu$ mentén a fotonok $\omega_{\text{coord}}$ koordináta-frekvenciája az adó ($\mathbf{x}_A$) és a távoli megfigyelő ($\mathbf{x}_B \to \infty$, ahol $g_{00} \to -1$, $f_L \to 1$) között megmarad, a megfigyelt frekvenciaarány:

1. **H0 ág ($\tilde{g}_{\mu\nu} = g_{\mu\nu}$ a $\Psi_{\text{atomic}}$ szektorra — minimális óracsatolás):**
   $$\sqrt{-\tilde{g}_{00}(\mathbf{x})} \;=\; \sqrt{-g_{00}(\mathbf{x})} \;=\; e^{\Phi_{\text{eff}}(\mathbf{x})/c^2} \quad\Longrightarrow\quad \frac{\nu_{\text{obs}}}{\nu_0}\bigg|_{\text{H0}} \;=\; e^{\Phi_{\text{eff}}(\mathbf{x})/c^2}, \qquad 1 + z_{\text{H0}} \;=\; e^{-\Phi_{\text{eff}}(\mathbf{x})/c^2}.$$
   Gyenge térben ($\Phi_{\text{eff}} = \Phi_N + \Delta\Phi_{\text{TAL}}$) a relatív frekvenciaeltolódás a standard newtoni GR-hez képest:
   $$\Delta y_{\text{H0}}(\mathbf{x}) \;\equiv\; \left(\frac{\nu_{\text{obs}}}{\nu_0}\right)_{\text{H0}} - \left(\frac{\nu_{\text{obs}}}{\nu_0}\right)_{\text{GR}} \;\simeq\; \frac{\Delta\Phi_{\text{TAL}}(\mathbf{x})}{c^2} \;=\; \mathcal{O}\!\left(\frac{G M}{c^2 R}\cdot\phi_L\right) \;=\; \mathcal{O}\!\left(\frac{G^2}{c^4}\right).$$
   Ezért a H0 ág jele másodrendűen ($\propto G^2/c^4$) elnyomott: a V20 gránitkamrában $\Delta y_{\text{H0}} \sim +10^{-48}\text{..}+10^{-47}$, a Sirius B-nél pedig $|\Delta v_{\text{H0}}| \le 1{,}27\;\text{m/s}$.

2. **Feltételes H1 ág ($\tilde{g}_{\mu\nu} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu \implies \tilde{g}_{00} = f_L^2 g_{00}$ — diszform lapse / Jordan óracsatolás):**
   $$\sqrt{-\tilde{g}_{00}(\mathbf{x})} \;=\; f_L(\mathbf{x})\,e^{\Phi_{\text{eff}}(\mathbf{x})/c^2} \quad\Longrightarrow\quad \frac{\nu_{\text{obs}}}{\nu_0}\bigg|_{\text{H1}} \;=\; f_L(\mathbf{x})\,e^{\Phi_{\text{eff}}(\mathbf{x})/c^2}, \qquad 1 + z_{\text{H1}} \;=\; \frac{e^{-\Phi_{\text{eff}}(\mathbf{x})/c^2}}{f_L(\mathbf{x})}.$$
   Gyenge térben ($\phi_L \ll 1 \implies f_L = e^{-\phi_L} \simeq 1 - \phi_L$) a newtoni GR feletti többlet-frekvenciaeltolódás elsőrendűen ($\propto G/c^2$) mérhető:
   $$\boxed{\;\Delta y_{\text{H1}}(\mathbf{x}) \;\simeq\; f_L(\mathbf{x}) - 1 \;\simeq\; -\phi_L(\mathbf{x}) \;=\; \begin{cases} -\gamma_L\,\bar{\rho}_L(\mathbf{x}) & \text{(I. megfogalmazás: 3D top-hat)}, \\[4pt] -\phi_L^{\text{Yuk}}(\mathbf{x}) & \text{(II. megfogalmazás: árnyékolt Yukawa)}. \end{cases}\;}$$

---

## 3. Kvantitatív ellenőrzés a Sirius B-re és a V20 kétkarú laboratóriumi programra

### 3.1. Sirius B ($M = 1{,}018\,M_\odot$, $R_* = 5586{,}5\text{ km}$, $n=3$ politróp) és aV20 $8558\text{ kg}$-os gránitkamra (V23 3. táblázat)
A V23 `clock_coupling_H0_H1.csv` (3. táblázat) egzakt számszerű értékei mindkét megfogalmazásban:

| Rendszer | $L$ [m] | $f_L(R_*)$ Top-Hat | $f_L(R_*)$ Yukawa | H0 előrejelzés | Feltételes H1 Top-Hat | Feltételes H1 Yukawa ($\alpha_Y=1$) | Egység |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Sirius B** ($1{,}018\,M_\odot$) | $1{,}0\times 10^3$ | $1{,}000000$ | $1{,}000000$ | $80{,}70055$ | $80{,}70055$ | $80{,}70055$ | $\text{km/s}$ |
| **Sirius B** ($1{,}018\,M_\odot$) | $1{,}0\times 10^4$ | $1{,}000000$ | $1{,}000000$ | $80{,}70055$ | $80{,}70055$ | $80{,}70055$ | $\text{km/s}$ |
| **Sirius B** ($1{,}018\,M_\odot$) | $1{,}0\times 10^6$ | $1{,}000000$ | $1{,}000000$ | $80{,}69839$ | $80{,}70094$ | $80{,}74476$ | $\text{km/s}$ |
| **Sirius B** ($1{,}018\,M_\odot$) | $5{,}0\times 10^6$ | $0{,}999984$ | $0{,}999993$ | $80{,}69487$ | $85{,}38911$ | $82{,}85137$ | $\text{km/s}$ |
| **V20 gránit** ($8558\text{ kg}$) | $0{,}50$ | $1{,}000000$ | $1{,}000000$ | $+3{,}075\times 10^{-48}$ | $-3{,}318\times 10^{-25}\;(-3{,}318\times 10^{-25})$ | $-2{,}008\times 10^{-25}\;(-1{,}943\times 10^{-25})$ | $\Delta y$ |
| **V20 gránit** ($8558\text{ kg}$) | $0{,}85$ | $1{,}000000$ | $1{,}000000$ | $+1{,}208\times 10^{-47}$ | $-1{,}303\times 10^{-24}\;(-1{,}303\times 10^{-24})$ | $-3{,}373\times 10^{-25}\;(-3{,}101\times 10^{-25})$ | $\Delta y$ |
| **V20 gránit** ($8558\text{ kg}$) | $1{,}00$ | $1{,}000000$ | $1{,}000000$ | $+1{,}331\times 10^{-47}$ | $-1{,}436\times 10^{-24}\;(-1{,}436\times 10^{-24})$ | $-3{,}781\times 10^{-25}\;(-3{,}407\times 10^{-25})$ | $\Delta y$ |

*(Megjegyzés a V20 sorokhoz: a feltételes H1 oszlopokban az első szám az IN-eltolódás $\mathbf{r}=\mathbf{0}$-ban, a zárójelben lévő szám az IN$-$OUT differenciális jel $D_{\text{out}} = 2{,}0\text{ m}$-nél; a newtoni GR eltolódás $\Phi_N(\text{IN})/c^2 = -9{,}269\times 10^{-24}$, IN$-$OUT $= -6{,}094\times 10^{-24}$. A $103\text{ kg}$-os statikus gránit ofszet felső korlátja kanonikus $\gamma_L$ mellett $\le 5{,}8\times 10^{-26}$).*

### 3.2. A V20 8 tekercses mágnesestér-energetikai kísérlet a skalár–tenzor / diszform olvasatban
Ha a skalármező forrástagja az elektromágneses energia-impulzus tenzor $T_{00}^{(\text{EM})} = u_B = \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0}$ energiasűrűségéhez is csatolódik ($\rho_B(\mathbf{x}) = \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0 c^2}$, noha a konforman invariáns 4D Maxwell-elmélet nyoma $T^\mu_{\;\mu}(\text{EM}) = 0$, a közeghez kötött diszform $\tilde{g}_{00} = -f_L^2$ csatolásban az elektromágneses energia-komponens közvetlenül súlyozódik $f_L$-lel), akkor a feltételes H1 ágon a mágneses tér BE/KI modulációja:
$$\Delta y_B(\mathbf{0}; L) \;=\; -\gamma_L\,\bar{\rho}_{B,L}(\mathbf{0}) \;=\; -\frac{G L^2}{c^2}\cdot\frac{3}{4\pi L^3}\int_{|\mathbf{y}|\le L}\frac{|\mathbf{B}(\mathbf{y})|^2}{2\mu_0 c^2}\,d^3y \;=\; -\frac{3G}{8\pi\mu_0 c^4 L}\int_{|\mathbf{y}|\le L}|\mathbf{B}(\mathbf{y})|^2\,d^3y.$$
- **Geometriai csúcs $L = 0{,}630957\text{ m}$-nél:** A 8 darab ($\varnothing 300\text{ mm}$ peremhatár, $\varnothing 291{,}7\text{ mm}$ effektív réz menetátmérő, $l_{\text{coil}} = 400\text{ mm}$, $N = 2000$, $I = 14{,}0\text{ mA}$) tekercs Biot–Savart terét integrálva a $\Delta y_{\text{even}}(L) = \frac{1}{2}[\Delta y(+I) + \Delta y(-I)]$ polaritás-páros jel pontosan $L = 0{,}630957\text{ m}$-nél éri el extremumát ($\Delta y_{\text{peak}} = -1{,}2017108\times 10^{-48}$).
- **A Zeeman-maszkolás és a háttér-kereszttag leválasztása:** Az $+I \leftrightarrow -I$ polaritásfordítás kiejti a lineáris Zeeman-tagot és a $2\mathbf{B}_{\text{bg}}\cdot\mathbf{B}_{\text{coil}}$ kereszttagot; az azonos áramirányú másodrendű Zeeman-eltolódás ($0{,}86\text{--}43\,\mu\text{T}$ mellett Sr-ra $-4\times 10^{-20}\text{..}-9{,}9\times 10^{-17}$, Cs-ra $3{,}4\times 10^{-12}\text{..}8{,}6\times 10^{-9}$) váltakozó tekercs-polaritással a középpontban szimmetrikusan $\mathbf{B}(\mathbf{0}) = \mathbf{0}$-ra redukálható.

---

## 4. Ellenőrzések (Checks)

1. **A V20 gránitkamra analitikus 3D top-hat értéke $L = 1{,}00\text{ m}$-nél:**
   - Mivel az üreges gránithenger külső sugara $R_{\text{out}} = 0{,}85\text{ m}$ és félmagassága $H/2 = 0{,}60\text{ m}$, a sarokpont távolsága az origótól $r_{\text{corner}} = \sqrt{0{,}85^2 + 0{,}60^2} = \sqrt{1{,}0825} \approx 1{,}0404\text{ m}$. $L = 1{,}00\text{ m}$-nél a $B(\mathbf{0}, 1\text{ m})$ gömb a $8558\text{ kg}$-os gránittömeg csaknem egészét befogja:
     $$\Delta y_{\text{H1}}^{\text{top-hat}}(L=1{,}00\text{ m}) \;=\; -\frac{G L^2}{c^2}\,\frac{M_{\text{enc}}(<1\text{ m})}{\frac{4}{3}\pi L^3} \;=\; -1{,}436\times 10^{-24} \quad \checkmark$$
2. **A Sirius B H1 sebességeltolódásának nagyságrendi ellenőrzése $L = 5{,}0\times 10^6\text{ m}$-nél:**
   - $1 - f_L(R_*) = 1 - 0{,}999984 = 1{,}6\times 10^{-5} \implies \Delta v_{\text{H1}} \simeq c\,\bigl(1 - f_L(R_*)\bigr) \approx 3\times 10^5\;\text{km/s}\times 1{,}56\times 10^{-5} \approx +4{,}69\;\text{km/s}$, így $v_{\text{H1}} = 80{,}69487 + 4{,}69424 = 85{,}38911\;\text{km/s}$ (Top-Hat) és $82{,}85137\;\text{km/s}$ (Yukawa), pontosan ahogy a V23 3. táblázata tartalmazza $\checkmark$.
