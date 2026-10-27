# 5. lépés: Csatolási állandók szótára, gyenge téri Yukawa ötödik erő és konform vs. diszform PPN fényelhajlás
# Step 5: Dictionary of Coupling Constants, Weak-Field Yukawa Fifth Force, and Conformal vs. Disformal PPN Light Deflection

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step5_coupling_dictionary_ppn_screening.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md`, `step2_conformal_disformal.md`, `step3_yukawa_helmholtz_eom.md`, `step4_covariant_conservation_hydrostatic.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` magyar kézirat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.2. szakasz, Eq. 6–7, 56–65. sor; 3. szakasz, 161–167. sor; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (3. fejezet, 245–255. sor, Eq. 4–5) és `hu_p2.tex` (4. fejezet, 23–38. sor; DOI: `10.5281/zenodo.23604423`)

---

## 1. Szó szerinti idézetek a V22 és V23 kéziratokból

A 2. lépésben (`step2_conformal_disformal.md`) rögzítettük a V19–V23 modellágak taxonómiáját, és megmutattuk, hogy különbséget kell tenni az univerzális konform Jordan-metrika (**Branch II-Conf**), az együttmozgó időirányú diszform metrika (**Branch II-Disf**) és az effektív folyadéksűrűség-csatolás (**Branch II-Fluid**) között.

1. **A $C$ csatolási állandó, $\gamma_L$ és a kompaktsági skála a V22 monográfiában (`hu_p1.tex`, 247–254. sor, Eq. 4–5):**
   $$\begin{equation}
     C \;\equiv\; \frac{4\pi G}{c^2} \;=\; 9{,}33306\times 10^{-27}\;\text{m}\,\text{kg}^{-1}, \qquad
     \gamma_L \;=\; \frac{G\,L^2}{c^2} \;=\; \frac{C\,L^2}{4\pi} \;=\; 7{,}425915\times 10^{-28}\,L^2\;\text{m}^3\,\text{kg}^{-1}. \tag{V22 Eq. 4–5}
   \end{equation}$$
   > *„Amikor egy $M$ ponttömeget $L$ sugarú gömbön átlagolunk, átlagsűrűsége $\bar{\rho}_L = \frac{3M}{4\pi L^3}$, így a kitevő $\gamma_L \bar{\rho}_L = \frac{3GM}{4\pi c^2 L} = \frac{3\,C\,M}{(4\pi)^2 L}$, ami a $3/(4\pi)$ geometriai szorzóig visszaadja a dimenziótlan $GM/(c^2 L)$ kompaktsági skálát.”*

2. **A Naprendszerbeli árnyékolás a V22 4. fejezetében (`hu_p2.tex`, 32. és 38. sor):**
   > *„Szubgalaktikus $L \gtrsim 1\,\text{pc}$ mellett $\gamma_L\bar{\rho}_L^\odot = \frac{3GM_\odot}{4\pi c^2 L} \lesssim 10^{-14}$; a bolygópályák érintetlenek maradnak! [...] amikor $L$ meghaladja a bolygópályák méretét, a Nap $M_\odot$ tömege felhígul a $\frac{4}{3}\pi L^3$ térfogaton, így az $f_L$ bolygópályák menti gradiensei a Cassini- és efemerisz-korlátok alá nyomódnak.”*

3. **A V23 Yukawa-normálási állandók (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 65. sor):**
   > *„Az egységre normált mag $\alpha_Y = 1$-nek, azaz $\alpha_{\text{std}} = 1/(4\pi)\approx 0{,}0796$-nak felel meg, ami visszaadja a homogén top-hat határesetet ($\phi_L = \gamma_L\rho_{\{src\}}$); az $\alpha_Y = 1/3$ ($\alpha_{\text{std}} = 1/(12\pi)$) a térfogatra illesztett változat.”*

---

## 2. Gyenge téri effektív potenciál, Yukawa ötödik erő és az effektív gravitációs állandó ($G_{\text{eff}}(r)$)

A 4. lépésben (`step4_covariant_conservation_hydrostatic.md`) levezettük, hogy a kovariáns megmaradási törvényből ($\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{(\phi)}) = 0$) a módosított hidrosztatikai / Euler-egyenlet:
$$h'(r) \;=\; -\nu'(r) \;+\; \phi_L'(r), \qquad \mathbf{a}_{\text{eff}}(\mathbf{r}) \;=\; -\nabla\Phi_{\text{eff}}(\mathbf{r}) \;=\; -\nabla\Bigl(\Phi_N(\mathbf{r}) \;-\; c^2\,\phi_L(\mathbf{r})\Bigr).$$

### 2.1. Külső Yukawa-potenciál kompakt forrás körül ($r > R_*$)
A 3. lépésben (`step3_yukawa_helmholtz_eom.md`, V23 Eq. 7) felírt lapos téri Green-függvényes megoldás egy $M$ tömegű, $R_*$ sugarú gömbszimmetrikus test külsejében ($r > R_*$, gyenge téri közelítésben):
$$\phi_L^{\text{Yuk}}(r) \;=\; \alpha_Y\,\frac{G}{c^2}\,\frac{L\,e^{-r/L}}{r}\int_0^{R_*} r'\,\sinh\!\left(\frac{r'}{L}\right)\rho_{\text{src}}(r')\,dr' \;=\; \frac{\alpha_{\text{std}}\,G M_{\text{eff}}(L)}{c^2\,r}\,e^{-r/L},$$
ahol $\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$, és pontszerű határesetben ($R_* \ll L$, ahol $\sinh(r'/L) \simeq r'/L$):
$$M_{\text{eff}}(L) \;=\; 4\pi\int_0^{R_*} r'^2\left[\frac{L}{r'}\sinh\!\left(\frac{r'}{L}\right)\right]\rho_{\text{src}}(r')\,dr' \;\xrightarrow{R_*\ll L}\; \int_0^{R_*} 4\pi r'^2\rho_{\text{src}}(r')\,dr' \;=\; M.$$
Behelyettesítve a $\Phi_N(r) = -\frac{GM}{r}$ newtoni potenciál mellé az effektív potenciálba:
$$\Phi_{\text{eff}}(r) \;=\; \Phi_N(r) \;-\; c^2\,\phi_L^{\text{Yuk}}(r) \;=\; -\frac{G M}{r}\Bigl(1 \;+\; \alpha_{\text{std}}\,e^{-r/L}\Bigr).$$
Differenciálva $r$ szerint megkapjuk a távolságfüggő effektív gravitációs állandót:
$$\boxed{\;a_{\text{eff}}(r) \;=\; -\frac{G_{\text{eff}}(r)\,M}{r^2}, \qquad G_{\text{eff}}(r) \;=\; G\left[1 \;+\; \alpha_{\text{std}}\left(1 + \frac{r}{L}\right)e^{-r/L}\right].\;}$$

### 2.2. Kanonikus skalár–tenzor normálás és az $\alpha_{\text{std}} = 2\beta^2$ kapcsolat
Írjuk át a V23 Eq. (8) skalár-hatást:
$$S_\phi \;=\; -\int d^4x\,\sqrt{-g}\,\frac{c^4}{16\pi G}\frac{1}{\tilde{\alpha}_L}\left(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{\phi_L^2}{L^2}\right) \;=\; -\int d^4x\,\sqrt{-g}\,\frac{c^4}{2 G\,\alpha_Y}\left(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{\phi_L^2}{L^2}\right)$$
(ahol $\alpha_Y = 8\pi\tilde{\alpha}_L = 4\pi\alpha_{\text{std}}$, tehát $\tilde{\alpha}_L = \frac{\alpha_{\text{std}}}{2} = \frac{\alpha_Y}{8\pi}$, és $G=c=1$ egységekben a skalár Lagrangiánus pontosan $-\frac{1}{2\alpha_Y}[(\nabla\phi_L)^2 + \phi_L^2/L^2]$) a Damour–Esposito-Farèse (1992) / Damour–Polyakov (1994) kanonikus Einstein-keretbeli $\varphi$ skalármezőre, amelynek kinetikus tagja $-\frac{c^4}{4\pi G}\,\frac{1}{2}(\nabla\varphi)^2 = -\frac{c^4}{8\pi G}(\nabla\varphi)^2$:
$$\frac{c^4}{2 G\,\alpha_Y}(\nabla\phi_L)^2 \;=\; \frac{c^4}{8\pi G\,\alpha_{\text{std}}}(\nabla\phi_L)^2 \;\equiv\; \frac{c^4}{8\pi G}\,2(\nabla\varphi)^2 \quad\implies\quad \phi_L \;=\; \sqrt{2\alpha_{\text{std}}}\,\varphi \;=\; 2\beta\,\varphi.$$
Mivel az anyagcsatolási tényező $f_L = e^{-\phi_L} = e^{-2\beta\varphi}$, a dimenziótlan skalár–anyag csatolási állandó ($\beta$) és a V23-as $\alpha_{\text{std}}, \alpha_Y, \tilde{\alpha}_L$ paraméterek közötti egzakt kapcsolat:
$$\boxed{\;\alpha_{\text{std}} \;=\; 2\beta^2 \;=\; 2\tilde{\alpha}_L \;=\; \frac{\alpha_Y}{4\pi}, \qquad \beta \;=\; \sqrt{\frac{\alpha_{\text{std}}}{2}} \;=\; \sqrt{\frac{\alpha_Y}{8\pi}}.\;}$$
- Az **egységre normált top-hat megfelelőre** ($\alpha_Y = 1$):
  $$\alpha_{\text{std}} \;=\; \frac{1}{4\pi} \;\approx\; 0{,}079577, \qquad \beta \;=\; \frac{1}{\sqrt{8\pi}} \;\approx\; 0{,}199471, \qquad 2\omega_{\text{BD}} + 3 \;=\; \frac{1}{\beta^2} \;=\; 8\pi \implies \omega_{\text{BD}} \;=\; 4\pi - \frac{3}{2} \;\approx\; 11{,}066.$$
- A **térfogatra illesztett változatra** ($\alpha_Y = 1/3$):
  $$\alpha_{\text{std}} \;=\; \frac{1}{12\pi} \;\approx\; 0{,}026526, \qquad \beta \;=\; \frac{1}{\sqrt{24\pi}} \;\approx\; 0{,}115165, \qquad \omega_{\text{BD}} \;=\; 12\pi - \frac{3}{2} \;\approx\; 36{,}199.$$

---

## 3. A három csatolási ág (Branch II-Conf, Branch II-Fluid, Branch II-Disf) szigorú szétválasztása a PPN és fényelhajlási szektorban

A bírálati észrevételeknek megfelelően rögzítjük, hogy a $\tilde{g}_{\mu\nu}^{\text{conf}}$ konform, a $\tilde{g}_{\mu\nu}^{\text{disf}}$ diszform és az $S_m = -\int d^4x\sqrt{-g}\,e^{-\phi_L}\varepsilon(n)$ effektív folyadék-csatolás **három különböző fizikai feltevés** az elektromágneses/nullgeodetikus szektorra nézve, ezért PPN-viselkedésük is eltérő:

1. **Branch II-Conf — Univerzális konform Jordan-csatolás ($\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$, $\sqrt{-\tilde{g}_{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$):**
   - Ebben az ágban minden anyagi és elektromágneses mező univerzálisan a $\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$ metrikához csatolódik.
   - Mivel a Maxwell-hatás $D=4$ dimenzióban konform-invariáns ($\sqrt{-\tilde{g}_{\text{conf}}}\,\tilde{g}_{\text{conf}}^{\mu\alpha}\tilde{g}_{\text{conf}}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} = \sqrt{-g}\,g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}$), a fotonok null-geodetikusai ($\tilde{g}_{\mu\nu}^{\text{conf}}dx^\mu dx^\nu = 0 \iff g_{\mu\nu}dx^\mu dx^\nu = 0$) kizárólag az Einstein-keretbeli $g_{\mu\nu}$ metrikát érzékelik, és az elektromágneses tér nem forrása a $\phi_L$ skalárnak ($T^\mu_{\;\mu}(\text{EM}) = 0$).
   - **Kizárólag erre a Branch II-Conf ágra** érvényes a standard skalár–tenzor / Brans–Dicke PPN-paraméter (Will, *Living Rev. Relativ.* **17**, 4, 2014):
     $$\boxed{\;\gamma_{\text{PPN}}^{\text{conf}}(r) \;=\; \frac{1 - \alpha_{\text{std}}\,e^{-r/L}}{1 + \alpha_{\text{std}}\,e^{-r/L}},\qquad \gamma_{\text{PPN}}^{\text{conf}}(r \ll L) \;=\; \frac{2\omega_{\text{BD}} + 1}{2\omega_{\text{BD}} + 3} \;=\; \frac{4\pi - 1}{4\pi + 1} \;\approx\; 0{,}8526\;(\alpha_Y=1),\qquad \gamma_{\text{PPN}}^{\text{conf}}(r \gg L) \;\to\; 1.\;}$$

2. **Branch II-Fluid + $H_0$ — Effektív barionos folyadék-csatolás minimális foton-/atomi metrikával ($S_m = -\int d^4x\sqrt{-g}\,e^{-\phi_L}\varepsilon(n)$, $\Psi_{\text{atomic/EM}}$ a $g_{\mu\nu}$-höz csatolódik):**
   - A makroszkopikus barionos folyadékra hat az ötödik erő ($h' = -\nu' + \phi_L'$), így két makroszkopikus test közötti effektív gravitációs állandó $G_{\text{eff}}(r) = G\bigl[1 + \alpha_{\text{std}}(1 + r/L)e^{-r/L}\bigr]$.
   - Azonban $r \gg L$ tartományban (például ha $L$ szubmilliméteres vagy csillagászati skálán $r_{\text{AU}} \gg L$) a Yukawa-tag exponenciálisan lecseng ($e^{-r/L} \to 0$), a külső téridő aszimptotikusan Schwarzschild-geometria ($M_{\text{ADM}}$ tömeggel), és $\gamma_{\text{PPN}}(r \gg L) = 1$, $\beta_{\text{PPN}}(r \gg L) = 1$.

3. **Branch II-Disf + $H_1$ — Együttmozgó időirányú diszform csatolás ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$, $\sqrt{-\tilde{g}_{\text{disf}}} = f_L\sqrt{-g}$):**
   - Ebben az ágban a térfogatelem $\sqrt{-\tilde{g}_{\text{disf}}} = e^{-\phi_L}\sqrt{-g}$ pontosan megegyezik a Branch II-Fluid $e^{-\phi_L}$ súlyfaktorával, de a négyessebesség $u^\mu$ explicite belép a metrikába (ami vagy együttmozgó khronon/éter rögzítést igényel, vagy szabad $u^\mu$-variáció esetén további kényszeregyenletet ad).
   - Mivel a diszform transzformáció nem konform, a Maxwell-szektorban megtöri a konform invarianciát: egyrészt a statikus mágneses tér $\rho_B = \frac{|\mathbf{B}|^2}{2\mu_0 c^2}$ effektív tömegsűrűséggel forrása lesz $\phi_L$-nek (15. és 20. lépés), másrészt a fénysugarak a $\tilde{g}_{\mu\nu}^{\text{disf}}k^\mu k^\nu = 0$ diszform fénykúpon terjednek, $n_{\text{refr}}(r) = f_L(r)^{-1} = e^{+\phi_L(r)}$ effektív törésmutatóval (14. lépés).
   - Következésképpen a Branch II-Disf ág fényelhajlását **nem** a konform $\gamma_{\text{PPN}}^{\text{conf}} = (1 - \alpha_{\text{std}}e^{-r/L})/(1 + \alpha_{\text{std}}e^{-r/L})$ Brans–Dicke-képlet írja le, hanem a 14. lépésben levezetett $\hat{\alpha}_{\text{disf}}(b) = \frac{4GM}{c^2 b} - \int_{-\infty}^\infty \frac{\partial_\perp\phi_L}{c^2}dz$ refraktív elhajlási szög.

---

## 4. Összefoglaló csatolási szótár (TAL $\leftrightarrow$ Skalár–tenzor / Dilaton)

| Fizikai mennyiség | Time Attraction Law (V19–V23) | Skalár–tenzor / Dilaton megfelelő | Egzakt átszámítás és ág-specifikus megjegyzés |
| :--- | :--- | :--- | :--- |
| **Dimenziótlan skalármező** | $\phi_L(\mathbf{x})$ | $2\beta\varphi(\mathbf{x})$ (vagy $2\phi_{\text{dil}}(\mathbf{x})$) | $f_L(\mathbf{x}) = e^{-\phi_L} \leftrightarrow e^{-2\beta\varphi} \leftrightarrow e^{-2\phi_{\text{dil}}}$ |
| **Alapvető gravitációs–időbeli állandó** | $C \equiv \frac{4\pi G}{c^2} = 9{,}33306\times 10^{-27}\;\text{m/kg}$ | $8\pi G / c^4$ ($\times c^2/2$) | $\frac{G}{c^2} = \frac{C}{4\pi} = 7{,}425915\times 10^{-28}\;\text{m/kg}$ |
| **Skálafüggő sűrűség-csatolás** | $\gamma_L = \frac{G L^2}{c^2} = \frac{C L^2}{4\pi}$ ($\text{m}^3/\text{kg}$) | $\frac{G}{c^2\mu_\phi^2}$ | $[\gamma_L] = \text{m}^3\,\text{kg}^{-1} = [\text{sűrűség}]^{-1}$ |
| **Simítási sugár / Hatótávolság** | $L$ (3D top-hat / Yukawa-skála) | Compton-hullámhossz $\lambda_\phi = \mu_\phi^{-1}$ | $m_\phi = \frac{\hbar}{L c}$ (skalár-/dilaton-tömeg) |
| **Hatás-normálási állandók** | $\alpha_Y = 4\pi\alpha_{\text{std}} = 8\pi\tilde{\alpha}_L$ | $8\pi\beta^2 = \frac{4\pi}{2\omega_{\text{BD}} + 3}$ | Hatásbeli skalár előfaktor: $-\frac{c^4}{16\pi G}\frac{1}{\tilde{\alpha}_L} = -\frac{c^4}{2 G\alpha_Y}$ ($G=c=1$-ben $-\frac{1}{2\alpha_Y}$) |
| **Konform ág (Branch II-Conf)** | $\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$ | Jordan–Brans–Dicke metrika | $\sqrt{-\tilde{g}_{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$; $\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$, $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2}$ |
| **Diszform ág (Branch II-Disf)** | $\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$ | Bekenstein-féle diszform metrika | $\tilde{g}_{00} = f_L^2 g_{00},\; \tilde{g}_{ij} = g_{ij},\; \sqrt{-\tilde{g}_{\text{disf}}} = f_L\sqrt{-g}$; mágneses $\rho_B$ forrás és $n=e^{\phi_L}$ törésmutató |

---

## 5. Ellenőrzések (Checks)

1. **Yukawa-határértékek ($r \ll L$ és $r \gg L$):**
   - $r \gg L$ (erősen árnyékolt tartomány, $e^{-r/L} \to 0$): $G_{\text{eff}}(r) \to G$ és $\gamma_{\text{PPN}}^{\text{conf}}(r) \to 1$, tehát a modell minden $r \gg L$ távolságon visszaadja az általános relativitáselméletet $\checkmark$.
   - $r \ll L$ (árnyékolatlan kis távolság, $e^{-r/L} \to 1$): $G_{\text{eff}}(r) \to G(1 + \alpha_{\text{std}}) = G(1 + 2\beta^2)$, ami pontosan megegyezik a standard skalár–tenzor gravitáció kéttest-kölcsönhatási képletével (Damour & Esposito-Farèse 1992, Eq. 2.12) $\checkmark$.
2. **A 3D top-hat Naprendszerbeli hígulási korlátja (V22 `hu_p2.tex`, 32. sor):**
   A Nap tömegére ($M_\odot = 1{,}98847\times 10^{30}\;\text{kg}$) $L = 1\;\text{pc} = 3{,}0857\times 10^{16}\;\text{m}$ mellett:
   $$\gamma_L\bar{\rho}_L^\odot \;=\; \frac{3 G M_\odot}{4\pi c^2 L} \;=\; \frac{3 \times 1{,}4766\times 10^3\;\text{m}}{4\pi \times 3{,}0857\times 10^{16}\;\text{m}} \;\approx\; 1{,}14\times 10^{-14} \;\ll\; 10^{-5} \quad (\text{Cassini PPN-korlát}) \quad \checkmark$$
