# 9. lépés: Nagyenergiájú fotonterjedés a diszform/konform skalár–tenzor metrikában, az $n_t(s) = f_L(s)^{-1}$ effektív törésmutató és a többforrású V17 H.E.S.S./MAGIC gamma-időzítési audit
# Step 9: High-Energy Photon Propagation in the Disformal/Conformal Scalar-Tensor Metric, the Effective Refractive Index $n_t(s) = f_L(s)^{-1}$, and the Multi-Source V17 H.E.S.S./MAGIC Gamma-Ray Timing Audit

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step9_vhe_photon_propagation_hess_dispersion.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step8_microscopic_clock_coupling_siriusb_lab_v20.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p2.tex` (6. fejezet, 63–99. sor, Eq. 9 és az 5 forrásos H.E.S.S./MAGIC kontrolltáblázat; DOI: `10.5281/zenodo.23604423`)  
- V17 Zenodo rekord: DOI `10.5281/zenodo.19666220` (SPARC 155 galaxis és többforrású H.E.S.S. VHE gamma-időzítési audit)  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1. és 3. szakasz, Eq. 1 és 161–166. sor; DOI: `10.5281/zenodo.23639084`)

---

## 1. Szó szerinti idézetek a V22 monográfia 6. fejezetéből (`hu_p2.tex`)

A 8. lépésben (`step8_microscopic_clock_coupling_siriusb_lab_v20.md`) levezettük az együttmozgó atomóra ($\Psi_{\text{atomic}} = (\psi_e, A_\mu)$) stacionárius sajátidő- és vöröseltolódás-csatolását a $\tilde{g}_{\mu\nu}$ effektív metrikán. A következő logikai lépés a Maxwell-fotonmező ($F_{\mu\nu}$) **asztronomikus és kozmológiai bázistávolságokon ($s \in [0, D]$) történő terjedésének** levezetése a skalár–tenzor / diszform hatáselvből, összevetve a V22 monográfia 6. fejezetében dokumentált többforrású V17 H.E.S.S./MAGIC gamma-időzítési audittal:

1. **Az időbeli hordozó helyi törésmutatója ($n_t(s) = f_L(s)^{-1}$), az akromatikus késleltetés és a diszperziós paraméterezés (`hu_p2.tex`, 6. fejezet, 68–74. sor, Eq. 9):**
   > *„Amikor egy VHE gamma-foton kozmikus voidokon, filamenteken, gazdagalaxis-halókon és galaktikus csillagközi struktúrákon halad át egy $s \in [0, D]$ látóirány mentén, az időbeli hordozó helyi törésmutatója $n_t(s) = f_L(s)^{-1} = \exp(+\gamma_L \bar{\rho}_L(s))$. Tiszta metrikus (geometriai) határesetben, ahol minden fotonenergia azonosan csatolódik az effektív $f_L(s)$ lapse-faktorhoz, a beérkezési idő akromatikus Shapiro-jellegű többletkésést szenved: $\Delta t_{\text{geom}} = c^{-1}\int_0^D \bigl[\exp(\gamma_L\bar{\rho}_L(s)) - 1\bigr]\,ds$. Ha azonban a nagy $E$ energiájú foton-hullámcsomagok véges hullámhosszú diszperziót szenvednek az $L$ skálájú $\nabla\bar{\rho}_L$ térbeli gradienseken szóródva, az $E_1$ és $E_2$ energiacsatornák közötti energiafüggő időkésés paraméterezett alakja:”*
   $$\begin{equation}
     \Delta t(E_1, E_2) \;=\; \tau_n\,\frac{E_2^n - E_1^n}{E_{\text{ref}}^n}\int_0^D \gamma_L\,\bar{\rho}_L(s)\,(1+z(s))^n\,ds, \qquad n \in \{1, 2\}, \tag{V22 Eq. 9}
   \end{equation}$$
   > *„ahol $\tau_n$ a diszperziós csatolási együttható, $D(z)$ pedig a vöröseltolódással súlyozott együttmozgó úthossz.”*

2. **Az ötforrásos V17 H.E.S.S. / MAGIC kontrolltáblázat (`hu_p2.tex`, 6. fejezet, 78–91. sor):**

   | VHE-forrás | Vöröseltolódás / Táv. | Energiasáv ($\text{TeV}$) | Megfigyelt késés / Korlát | V17 audit-eredmény |
   | :--- | :--- | :--- | :--- | :--- |
   | **PKS 2155-304 (H.E.S.S.)** | $z = 0{,}116$ | $0{,}2\text{--}4{,}0\,\text{TeV}$ | $\Delta t \approx +28 \pm 12\,\text{s}$ | Marginális $\sim 2\sigma$ késés; a belső SSC-emissziós profil dominálja. |
   | **Mrk 501 (1997 / H.E.S.S.)** | $z = 0{,}034$ | $0{,}15\text{--}10\,\text{TeV}$ | $\Delta t \approx +240 \pm 60\,\text{s}$ | $3{,}4\times$ kisebb $z$-nél nagyobb késés sérti a monoton távolságskálázást! |
   | **GRB 190114C (MAGIC)** | $z = 0{,}4245$ | $0{,}2\text{--}1{,}0\,\text{TeV}$ | $|\Delta t / \Delta E| \le 19\,\text{s/TeV}$ | Szigorú kozmológiai korlát, amely kizárja a nagy lineáris vákuumdiszperziót. |
   | **Crab pulzár (H.E.S.S./MAGIC)** | $d \approx 2{,}0\,\text{kpc}$ | $0{,}1\text{--}1{,}5\,\text{TeV}$ | $|\Delta t| < 0{,}4\,\text{ms}$ | Szoros galaktikus korongbeli korlát sűrű ISM-en át ($\bar{n}_H \sim 1\,\text{cm}^{-3}$). |
   | **Vela pulzár (H.E.S.S. II)** | $d \approx 0{,}287\,\text{kpc}$ | $0{,}02\text{--}0{,}1\,\text{TeV}$ | $|\Delta t| < 0{,}2\,\text{ms}$ | Szub-ezredmásodperces fázisilleszkedés a sűrű Vela SNR-burkon át. |

3. **A V17 többforrású audit két fő következtetése (`hu_p2.tex`, 6. fejezet, 95–98. sor):**
   > - **A naiv távolságarányos diszperzió elvetése:** *„Mivel az Mrk 501 ($z=0{,}034$) nagyobb látszólagos kitörési késést mutat, mint a PKS 2155-304 ($z=0{,}116$), miközben a GRB 190114C ($z=0{,}4245$) és a Crab/Vela pulzárok szoros null-korlátokat adnak, a PKS 2155-304 $28\,\text{s}$-os késésének egyenletes intergalaktikus diszperziós állandóval való magyarázata megcáfolódik. Először a blazár belső szinkrotron-ön-Compton (SSC) gyorsítási időskáláit kell modellezni.”*
   > - **Egyedi környezeti túlsűrűség-előrejelzés a $\Lambda$CDM-mel szemben:** *„A Lorentz-invariancia-sértő (LIV) modellekkel ellentétben, amelyek pusztán vöröseltolódás-függő $\Delta t \propto \int_0^z \frac{(1+z')\,dz'}{H(z')}$ késést jósolnak, a Time Attraction Law azt jósolja, hogy bármely valódi nem-geometriai időzítési maradéknak a barionos túlsűrűség-integrállal ($\int \bar{\rho}_L(s)\,ds$) kell korrelálnia a gazdag galaxishalmazokon (például Coma vagy Perseus) áthaladó látóirányok mentén, az azonos vöröseltolódású kozmikus voidokhoz képest.”*

---

## 2. Matematikai levezetés: Fénykúp, effektív törésmutató és EFT-diszperzió a skalár–tenzor keretben

### 2.1. Konform invariancia vs. diszform fénykúp: Az $n_t(s) = f_L(s)^{-1}$ törésmutató levezetése
Vizsgáljuk meg a szabad Maxwell-hatást ($S_{\text{EM}} = -\frac{1}{4\mu_0}\int d^4x\sqrt{-\tilde{g}}\,\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}$) a 2. lépésben (`step2_conformal_disformal.md`) levezetett két metrikus transzformáció alatt:

1. **Tiszta konform transzformáció ($\tilde{g}_{\mu\nu}^{\text{conf}} = f_L^2 g_{\mu\nu}$):**
   $D = 4$ téridő-dimenzióban $\sqrt{-\tilde{g}^{\text{conf}}} = f_L^4\sqrt{-g}$ és $\tilde{g}_{\text{conf}}^{\mu\alpha}\tilde{g}_{\text{conf}}^{\nu\beta} = f_L^{-4}g^{\mu\alpha}g^{\nu\beta}$, így a $f_L^4 f_L^{-4} = 1$ tényező kiesik:
   $$S_{\text{EM}}\!\left[f_L^2 g_{\mu\nu}, A_\mu\right] \;=\; S_{\text{EM}}\!\left[g_{\mu\nu}, A_\mu\right].$$
   Emiatt egy tisztán konform skalár–tenzor csatolásban a fénykúp ($\tilde{g}_{\mu\nu}^{\text{conf}}k^\mu k^\nu = 0 \iff g_{\mu\nu}k^\mu k^\nu = 0$) azonos a háttérmetrika fénykúpjával.

2. **A V23 közeghez kötött diszform órametrikája ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$):**
   Ha az elektromágneses hullámok a barionos közeg $u^\mu$ négyessebességével definiált diszform $\tilde{g}_{\mu\nu}^{\text{disf}}$ metrikához csatolódnak, a lokális nyugalmi keretben ($g_{\mu\nu} \simeq \eta_{\mu\nu}$, $u^\mu = (c, \mathbf{0})$) a vonalelem:
   $$d\tilde{s}_{\text{disf}}^2 \;=\; -f_L(s)^2 c^2 dt^2 + ds^2, \qquad ds^2 = dx^2 + dy^2 + dz^2.$$
   A geometriai optikai (eikonál / null-geodetikus) határesetben $d\tilde{s}_{\text{disf}}^2 = 0$, amiből a foton koordináta-fénysebessége és a közeg effektív időbeli törésmutatója:
   $$v_{\text{eff}}(s) \;=\; \frac{ds}{dt} \;=\; c\,f_L(s) \;=\; c\,\exp\!\bigl(-\phi_L(s)\bigr) \quad\Longrightarrow\quad \boxed{\;n_t(s) \;\equiv\; \frac{c}{v_{\text{eff}}(s)} \;=\; f_L(s)^{-1} \;=\; \exp\!\bigl(+\phi_L(s)\bigr) \;=\; \exp\!\bigl(+\gamma_L\,\bar{\rho}_L(s)\bigr).\;}$$
   Az $s \in [0, D]$ útvonalon való integrálás közvetlenül kiadja a V22 `hu_p2.tex` 69. sorában szereplő akromatikus Shapiro-jellegű futásiidő-többletet:
   $$\Delta t_{\text{geom}} \;=\; \int_0^D \left(\frac{1}{v_{\text{eff}}(s)} - \frac{1}{c}\right)ds \;=\; \frac{1}{c}\int_0^D \Bigl[\exp\!\bigl(\gamma_L\bar{\rho}_L(s)\bigr) - 1\Bigr]\,ds \;\simeq\; \frac{1}{c}\int_0^D \phi_L(s)\,ds \qquad (\phi_L \ll 1).$$
   Mivel $\Delta t_{\text{geom}}$ nem függ a foton $E = \hbar\omega$ energiájától, azonos időben kibocsátott $E_1$ és $E_2$ energiájú fotonok között a tiszta eikonál-határesetben $\Delta t_{\text{geom}}(E_1, E_2) = 0$.

### 2.2. Véges hullámhosszú skalár–foton EFT-operátorok és a V22 Eq. (9) diszperziós képlet levezetése
Amikor a nagyenergiájú (VHE, $E \sim 0{,}1\text{--}10\;\text{TeV}$) gamma-fotonok hullámhossza kölcsönhatásba lép a skalármező $\phi_L = \gamma_L\bar{\rho}_L$ (vagy $\phi_L^{\text{Yuk}}$) gradienseivel, az effektív térelméleti (EFT) kifejtésben a vezető dimenzió-$4+n$ ($n \in \{1,2\}$) skalár–Maxwell operátor a közeg $u^\mu$ négyessebességével:
$$\Delta\mathcal{L}_{\text{disp}}^{(n)} \;=\; -\frac{c\,\tau_n}{4\mu_0}\,\frac{\phi_L(x)}{E_{\text{ref}}^n}\,F_{\mu\nu}\,\bigl(i\hbar\,u^\alpha\nabla_\alpha\bigr)^n F^{\mu\nu}.$$
Síkbeli hullámközelítésben ($A_\mu \propto e^{i k_\alpha x^\alpha}$, ahol a lokális együttmozgó fotonenergia $E_{\text{loc}}(s) = -\hbar u^\alpha k_\alpha = E_{\text{obs}}\,(1 + z(s))$) a módosított diszperziós reláció:
$$\omega^2 \;=\; c^2 k^2\left[1 - \frac{2c\,\tau_n}{n+1}\,\phi_L(s)\,\left(\frac{E_{\text{obs}}(1+z(s))}{E_{\text{ref}}}\right)^n\right].$$
A csoportsebesség ($v_g = \partial\omega/\partial k$) reciprokát integrálva két különböző megfigyelt energiacsatorna ($E_1 < E_2$) között pontosan a V22 (9) egyenletét kapjuk:
$$\boxed{\;\Delta t(E_1, E_2) \;=\; t_{\text{arr}}(E_2) - t_{\text{arr}}(E_1) \;=\; \tau_n\,\frac{E_2^n - E_1^n}{E_{\text{ref}}^n}\int_0^D \phi_L(s)\,\bigl(1 + z(s)\bigr)^n\,ds \;=\; \tau_n\,\frac{E_2^n - E_1^n}{E_{\text{ref}}^n}\int_0^D \gamma_L\,\bar{\rho}_L(s)\,\bigl(1 + z(s)\bigr)^n\,ds.\;}$$

---

## 3. Ellenőrzések és fizikai következmények (Checks)

1. **Miért zárja ki az 5 forrásos V17 audit az univerzális vákuumdiszperziót ($\phi_L = \text{konst.}$), és hogyan védi meg a lokális $\phi_L(\mathbf{x})$ csatolást?**
   - Ha a diszperziót sűrűségfüggetlen vákuum-LIV okozná ($\phi_L(s) \equiv \text{konst.}$), akkor $\Delta t / \Delta E$ szigorúan monoton nőne a $z$ vöröseltolódással. Ezzel szemben az **Mrk 501** ($z = 0{,}034$, $\Delta t \approx +240 \pm 60\;\text{s}$) $3{,}4\times$ kisebb vöröseltolódásnál $\sim 8{,}6\times$ nagyobb késést mutat, mint a **PKS 2155-304** ($z = 0{,}116$, $\Delta t \approx +28 \pm 12\;\text{s}$), míg a $z = 0{,}4245$-nél lévő **GRB 190114C** $|\Delta t/\Delta E| \le 19\;\text{s/TeV}$ és a galaktikus **Crab** ($|\Delta t| < 0{,}4\;\text{ms}$) / **Vela** ($|\Delta t| < 0{,}2\;\text{ms}$) szoros felső korlátot ad $\checkmark$.
2. **A skalár–tenzor / Yukawa TAL megkülönböztető előrejelzése a LIV-modellekkel szemben:**
   - Mivel a skalármező $\phi_L(s)$ a látóirány menti barionos tömegsűrűséggel arányos ($\phi_L = \gamma_L\bar{\rho}_L$ vagy $\phi_L^{\text{Yuk}}$), azonos $z$ vöröseltolódású források esetén a látóirány menti integrál $\mathcal{I}_L \equiv \int_0^D \phi_L(s)(1+z(s))^n\,ds$ arányos a látóirányba eső galaxishalmazok (pl. Coma, Perseus) barionos túlsűrűségével a kozmikus voidokhoz képest $\checkmark$.
3. **Vákuum-határeset ($\bar{\rho}_L \to 0 \implies \phi_L \to 0$):**
   - $\phi_L(s) \to 0 \implies n_t(s) = e^0 = 1$, $\Delta t_{\text{geom}} = 0$ és $\Delta t(E_1, E_2) = 0$: a fényterjedés egzaktul visszaadja a speciális/általános relativitáselméleti vákuum-fénykúpot $\checkmark$.
