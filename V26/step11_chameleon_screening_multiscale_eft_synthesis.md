# 11. lépés: Környezetfüggő skalár–tenzor árnyékolási kiterjesztési hipotézis (Branch II-Cham: $V_{\text{eff}}(\phi_L; \rho)$, kaméleon/szimmetron $\beta_{\text{eff}}(\rho)$), az $L_{\text{eff}}(\rho)$ korrelációs hossz és a többskálás EFT-szintézis
# Step 11: Environment-Dependent Scalar-Tensor Screening Extension Hypothesis (Branch II-Cham: $V_{\text{eff}}(\phi_L; \rho)$, Chameleon/Symmetron $\beta_{\text{eff}}(\rho)$), Correlation Length $L_{\text{eff}}(\rho)$, and Multi-Scale EFT Synthesis

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step11_chameleon_screening_multiscale_eft_synthesis.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step10_bilateral_kernel_reciprocity_eft_synthesis.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (1. és 3. fejezet, 171–180. és 256–265. sor, Eq. 6), `hu_p2.tex` (4–5. fejezet, 24–61. sor) és `hu_p3.tex` (11–12. fejezet, 99–120. sor; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1–3. szakasz; DOI: `10.5281/zenodo.23639084`)

---

## 1. Szó szerinti idézetek és a Branch II-Cham kiterjesztés módszertani elhatárolása

Az 1–10. lépésben rögzített, állandó $L$ korrelációs/Yukawa-hossz és állandó $\alpha_Y \in \{1, 1/3\}$ ($\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$) csatolás mellett vezettük le a V23 I. (3D gömbi top-hat, **Branch I-A/I-B**) és II. (lokális kvadratikus Helmholtz–Yukawa skalár–tenzor, **Branch II-Fluid / II-Conf / II-Disf**) megfogalmazás egyenleteit. Ahogy a 6. lépésben (`step6_sparc_galactic_ngc3198_consistency.md`) az NGC 3198 ($R = 30{,}08\;\text{kpc}$) negatív kontrollvizsgálatával kimutattuk, az állandó $\alpha_{\text{std}} = 1/(4\pi) \approx 0{,}07958$ csatolású izolált V23 modell legfeljebb $C_{\text{pred}}^{\text{isol}} \le 1 + \alpha_{\text{std}} = 1{,}0796$-szoros sebességnégyzet-erősítést ad, szemben a megfigyelt $C_{\text{req}} = 2{,}720\text{--}4{,}311$ értékkel.

**Módszertani elhatárolás (Point 5):** Az ebben a lépésben vizsgált sűrűségfüggő $\beta_{\text{eff}}(\rho)$ és $V(\phi_L)$ kaméleon/szimmetron mechanizmus (**Branch II-Cham** a 2. lépés taxonómiájában) **önálló nemlineáris továbbfejlesztési hipotézis**, amely túlmutat a V23 (`50aa2db0`) kvadratikus, konstans $\alpha_Y$-ú Yukawa-hatásán, és a V22 monográfia 11. fejezetének **(d)** nyitott feladatára (`hu_p3.tex`, 104. sor: *„az $L$ fizikai skálaválasztási mechanizmusával”*) ad lehetséges 4D effektív térelméleti (EFT) keretet.

1. **A V5–V16 kaméleon- és Horndeski-előzmények, valamint a $C(\rho)$ erőszorzó (`hu_p1.tex`, 172–174. és 256–265. sor, Eq. 6):**
   > *„V5 – V7: Gyenge tér és kaméleon-vizsgálatok — Skalár–tenzor és kaméleon-jellegű árnyékolási analógiákat vizsgáltam, hogy a helyi sűrűséggradienseket összekapcsoljam az alacsony sűrűségű térbeli effektív gyorsulással.”*
   > *„A 14–16. verziókban (DOI `10.5281/zenodo.23217189`) a galaktikus forgási görbékre vonatkozó első fenomenológiai vizsgálataim során az effektív gravitációs erősítést az alábbi sűrűségfüggő erőszorzóval paramétereztem:”*
   $$\begin{equation}
     C(\rho) \;=\; 1 + C_0\,\exp\!\left(-\frac{\rho}{\rho_c}\right), \qquad V_{\text{rot}}^2(r) \;=\; C(\rho(r))\,V_{\text{bar}}^2(r), \tag{V22 Eq. 6}
   \end{equation}$$
   > *„ahol $\rho_c$ a kritikus átmeneti sűrűség, a kezdeti munkákban pedig $C_0 = 5$-öt állítottam be (így $\rho \ll \rho_c$ esetén $C(\rho) \to 6$, nagy sűrűségnél pedig $\rho \gg \rho_c \implies C(\rho) \to 1$, visszaadva a Naprendszerbeli newtoni határesetet).”*

2. **A V22 monográfia 11. fejezetének (d) nyitott elméleti feladata (`hu_p3.tex`, 104. sor):**
   > *„**(d) Önkonzisztens, többcsatornás megfigyelési kiértékelés:** Az (a)–(c) pontokban levezetett hatáselv és anyagcsatolás birtokában a galaktikus forgási görbéket (a kozmikus $\rho_{\text{bg}}$ háttérbe ágyazással együtt), a halmaz-túlsűrűségeken áthaladó VHE gamma-időzítést, a Sirius B-t és a NICER-poszteriorokat együttesen, **az $L$ fizikai skálaválasztási mechanizmusával** kell illeszteni.”*

---

## 2. Matematikai levezetés: A környezetfüggő effektív potenciál $V_{\text{eff}}(\phi_L; \rho)$ és az $L_{\text{eff}}(\rho)$ skálaválasztás

### 2.1. Általános önkölcsönható skalár–tenzor hatás és az effektív potenciál
Általánosítsuk a V23 II. megfogalmazásának (Eq. 8) kvadratikus tömegtagját ($\frac{1}{2L_0^2}\phi_L^2$) egy általános $V(\phi_L)$ önkölcsönhatási potenciálra a 4 dimenziós skalár–tenzor hatásban:
$$S \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}R \;-\; \frac{\kappa_\phi}{2}g^{\mu\nu}\nabla_\mu\phi_L\nabla_\nu\phi_L \;-\; \kappa_\phi\,V(\phi_L) \;+\; f_L(\phi_L)\,\mathcal{L}_m\right], \qquad \kappa_\phi \;=\; \frac{c^4}{4\pi G\alpha_{\text{std}}} \;=\; \frac{c^4}{G\alpha_Y}, \quad f_L(\phi_L) = e^{-\phi_L}.$$
*(Megjegyzés a normáláshoz: $\kappa_\phi/2 = \frac{c^4}{2 G\alpha_Y} = \frac{c^4}{16\pi G}\frac{1}{\tilde{\alpha}_L}$, pontosan megegyezve a 3. és 7. lépésben használt V23 Eq. 8 normálással!)* Nemrelativisztikus barionos közegben ($\mathcal{L}_m = -\rho c^2$) a $\phi_L$ mezőre vonatkozó Euler–Lagrange-egyenlet:
$$\Box_g\phi_L \;=\; \frac{\partial V_{\text{eff}}(\phi_L; \rho)}{\partial\phi_L}, \qquad \boxed{\;V_{\text{eff}}(\phi_L; \rho) \;\equiv\; V(\phi_L) \;+\; \frac{\rho c^2}{\kappa_\phi}\,e^{-\phi_L} \;=\; V(\phi_L) \;+\; \frac{\alpha_Y G}{c^2}\,\rho\,e^{-\phi_L}\;}$$
(ahol $\alpha_Y = 4\pi\alpha_{\text{std}} = 8\pi\tilde{\alpha}_L$).

### 2.2. Egyensúlyi minimum $\phi_{\min}(\rho)$ és a sűrűségfüggő Yukawa-hossz $L_{\text{eff}}(\rho)$
Egy lokálisan kvázi-homogén $\rho$ sűrűségű tartományban a skalármező a $\frac{\partial V_{\text{eff}}}{\partial\phi_L}\big|_{\phi_{\min}} = 0$ minimumot veszi fel, a kis $\delta\phi_L = \phi_L - \phi_{\min}(\rho)$ perturbációk pedig a lokális görbület által meghatározott **környezetfüggő Yukawa-hosszal** ($L_{\text{eff}}(\rho)$) terjednek:
$$\frac{\partial V_{\text{eff}}}{\partial\phi_L}\bigg|_{\phi_{\min}} \;=\; V'\bigl(\phi_{\min}(\rho)\bigr) \;-\; \frac{\alpha_Y G}{c^2}\,\rho\,e^{-\phi_{\min}(\rho)} \;=\; 0,$$
$$\boxed{\;L_{\text{eff}}^{-2}(\rho) \;\equiv\; \mu_{\text{eff}}^2(\rho) \;\equiv\; \frac{\partial^2 V_{\text{eff}}(\phi_L; \rho)}{\partial\phi_L^2}\bigg|_{\phi_L = \phi_{\min}(\rho)} \;=\; V''\bigl(\phi_{\min}(\rho)\bigr) \;+\; \frac{\alpha_Y G}{c^2}\,\rho\,e^{-\phi_{\min}(\rho)}.\;}$$
Vizsgáljuk meg ennek a képletnek a két alapvető fizikai esetét:

1. **A V23 kvadratikus (állandó csupasz tömegű) határesete (**Branch II-Fluid**, $V(\phi_L) = \frac{\phi_L^2}{2L_0^2}$):**
   Itt $V'(\phi_L) = \phi_L / L_0^2$ és $V''(\phi_L) = L_0^{-2}$. A minimumfeltétel egzaktul a Lambert-$W$ függvénnyel írható fel ($\phi_{\min}e^{\phi_{\min}} = \alpha_Y\gamma_{L_0}\rho$, így $\phi_{\min}(\rho) \simeq \alpha_Y\gamma_{L_0}\rho$ kis $\phi_L$-re), a lokális effektív Yukawa-hossz pedig már a puszta $e^{-\phi_L}\rho c^2$ csatolás miatt is **csökken a sűrűség növekedésével**:
   $$L_{\text{eff}}^{-2}(\rho) \;=\; L_0^{-2}\Bigl(1 + \phi_{\min}(\rho)\Bigr) \;\simeq\; L_0^{-2}\Bigl(1 + \alpha_Y\,\gamma_{L_0}\,\rho\Bigr) \quad\Longrightarrow\quad L_{\text{eff}}(\rho) \;=\; \frac{L_0}{\sqrt{1 + \phi_{\min}(\rho)}}.$$

2. **Nemlineáris kaméleon / szimmetron kiterjesztési hipotézis (**Branch II-Cham**, `hu_p1.tex` V5–V7):**
   Ha a V23 kvadratikus modellen túlmutatóan a csupasz potenciál $V(\phi_L) = \Lambda_0^4\,\exp(\sigma\phi_L)$ vagy $V(\phi_L) = M_0^4\phi_L^{-n}$ ($n > 0$), és a csatolási függvény általános nemlineáris $f_{\text{eff}}(\phi_L) = \exp[-A(\phi_L)]$ alakú, amelynek effektív meredeksége a lokális minimumban $\beta_{\text{eff}}(\rho) = \beta_0\exp\!\left(-\frac{\rho}{2\rho_c}\right)$, akkor az effektív gravitációs szorzó ($G_{\text{eff}}/G = 1 + 2\beta_{\text{eff}}^2(\rho)$) reprodukálja a V14–V16 fenomenológiai $C(\rho)$ erőszorzót (`hu_p1.tex`, Eq. 6):
   $$\boxed{\;\frac{G_{\text{eff}}(\rho)}{G} \;=\; 1 + 2\beta_{\text{eff}}^2(\rho) \;=\; 1 + 2\beta_0^2\exp\!\left(-\frac{\rho}{\rho_c}\right) \;\equiv\; C(\rho), \qquad C_0 \;\equiv\; 2\beta_0^2 \;=\; \alpha_{\text{std},0}.\;}$$
   - **Galaktikus külső tartomány ($\rho \ll \rho_c \sim 10^{-23}\text{--}10^{-21}\;\text{kg/m}^3$):** ebben a nemlineáris kiterjesztési hipotézisben az árnyékolás feloldódik, $2\beta_{\text{eff}}^2(\rho) \to C_0$ (az NGC 3198 $R = 30{,}08\text{ kpc}$ pontján $C_{\text{req}} - 1 = 1{,}720$ [$\Upsilon=1{,}0$] illetve $3{,}311$ [$\Upsilon=0{,}5$]).
   - **Naprendszerbeli és földi nagy sűrűségű tartomány ($\rho \gg \rho_c$):** az exponenciális faktor elnyomja a hosszú hatótávolságú skalár-erőszorzót ($C(\rho) - 1 = C_0 e^{-\rho/\rho_c} \to 0$), és csak a rövid $L$ skálájú kanonikus $\gamma_L = G L^2/c^2$ ($\alpha_{\text{std}} = 1/(4\pi)$) V23 módus marad meg.

---

## 3. Többskálás EFT-áttekintés a V17–V23 vizsgált tartományokon (ág-megjelöléssel)

Az alábbi táblázat minden vizsgált fizikai rendszerhez egyértelműen megadja, hogy melyik modellágon (**Branch I-A**, **Branch II-Fluid**, **Branch II-Disf**, **Branch II-Cham**) értelmezett az adott eredmény:

| Fizikai skála és rendszer | Modellág (`step2` taxonómia) | Jellemző sűrűség $\rho$ [$\text{kg/m}^3$] | Vizsgált $L$ skála | $\phi_L = \gamma_L\bar{\rho}_L$ (vagy $\phi_L^{\text{Yuk}}$) | Fő megfigyelhető mennyiség és V17–V23 bírálati eredmény |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **V20 8 tekercses mágneses tér** (`step8`, `step20`) | **Branch II-Disf / I-A** ($H_1$) | $\rho_B = \frac{B^2}{2\mu_0 c^2} \sim 10^{-22}$ | $L_{\text{peak}} = 0{,}630957\text{ m}$ | $1{,}2017108\times 10^{-48}$ | Polaritás-páros $\Delta y_{\text{peak}} = -1{,}2017108\times 10^{-48}$; Zeeman + kereszttag kiejtve. |
| **V20 $8558\text{ kg}$ gránitkamra** (`step8`, `step19`) | **Branch I-A / II-Fluid** ($H_0$ vs. $H_1$) | $2{,}75\times 10^3$ | $L = 0{,}50\text{--}1{,}00\text{ m}$ | $1{,}436\times 10^{-24}$ ($1\text{ m}$) | $H_0$: $+1{,}331\times 10^{-47}$; feltételes $H_1$: $-1{,}436\times 10^{-24}$ (top-hat), $-3{,}781\times 10^{-25}$ (Yukawa). |
| **Szíriusz B fehér törpe** (`step8`, `step19`) | **Branch I-A / II-Fluid** ($H_0$ vs. $H_1$) | $\bar{\rho} \sim 3\times 10^9$ | $L = 1\text{--}5000\text{ km}$ | $1{,}6\times 10^{-5}$ ($5000\text{ km}$) | $H_0$: $|\Delta v|\le 1{,}27\text{ m/s}$ (kompatibilis); $H_1$: $L \ge 1\text{ km}$-re kizárva ($L \lesssim 10^6\text{ m}$ lokális magra). |
| **Neutroncsillagok** (`step7`, `step18`, `step22`) | **Branch I-A** vs. **Branch II-Fluid** | $10^{17}\text{--}1{,}45\times 10^{18}$ | $L = 1\text{--}10\text{ km}$ | $0{,}016\text{--}0{,}050$ ($10\text{ km}$) | Görbült külső `yukawa_action_el` $M_{\text{ADM}}$ és csatolt $(\delta h,\delta\phi_L)$ stabilitás (`step22`); V21 nemlokális érintő BVP (`step18`). |
| **Naprendszer (PPN / Cassini)** (`step5`) | **Branch II-Conf** vs. **II-Disf / I-A** | $\bar{\rho}_L^\odot \propto L^{-3}$ | $L = 1\text{ CsE}\text{--}1\text{ pc}$ | $1{,}14\times 10^{-14}$ ($1\text{ pc}$) | 3D top-hat hígulás: $|1-f_L^\odot| \approx 1{,}14\times 10^{-14} \ll 10^{-5}$; $r \gg L$-re $\gamma_{\text{PPN}} \to 1$. |
| **SPARC (175/155) és NGC 3198** (`step6`) | **Branch I-A / II-Fluid** (negatív kontroll) vs. **II-Cham** | $10^{-24}\text{--}10^{-20}$ | $L \sim 1\text{--}23{,}09\text{ kpc}$ | $\le 2{,}300\times 10^{-7}$ (izolált) | 0 db $\rho \le 0$ pont (mag-pozitivitás); NGC 3198 ($30{,}08\text{ kpc}$): $C_{\text{req}} = 2{,}720 / 4{,}311$, míg konstans $\alpha_{\text{std}}$: $C_{\text{pred}}^{\text{isol}} \le 1{,}0796$. |
| **V17 H.E.S.S./MAGIC gamma-időzítés** (`step9`) | **Branch II-Disf** ($H_1$ EM) | $\rho_{\text{IGM/ISM}} \sim 10^{-27}\text{--}10^{-21}$ | $z = 0{,}034\text{--}0{,}4245$ | $\int_0^D \phi_L(s)\,ds$ | Naiv $z$-arányos LIV kizárva; halmaz-túlsűrűséggel ($\int\bar{\rho}_L\,ds$) korreláló előrejelzés. |

---

## 4. Ellenőrzések (Checks)

1. **A $C(\rho) = 1 + 2\beta_{\text{eff}}^2(\rho)$ megfeleltetés határértékei (Branch II-Cham):**
   - $\rho \to 0$ (mély galaktikus külső tartomány): $2\beta_{\text{eff}}^2(0) = 2\beta_0^2 = C_0 \implies C(0) = 1 + C_0$ $\checkmark$.
   - $\rho \gg \rho_c$ (csillagbelső, Naprendszer, laboratórium): $2\beta_{\text{eff}}^2(\rho) = C_0 e^{-\rho/\rho_c} \to 0 \implies C(\rho) \to 1$ $\checkmark$.
2. **Az $L_{\text{eff}}(\rho) = L_0 / \sqrt{1 + \phi_{\min}(\rho)}$ kvadratikus Yukawa-határeset ellenőrzése (Branch II-Fluid):**
   - Kis sűrűségnél ($\phi_{\min}(\rho) \to 0$) $L_{\text{eff}}(\rho) \to L_0$, visszaadva a V23 II. megfogalmazásának állandó $L_0$ Yukawa-hosszát $\checkmark$.
