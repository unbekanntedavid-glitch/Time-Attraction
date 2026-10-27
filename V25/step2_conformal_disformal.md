# 2. lépés: A konform, a diszform és az effektív folyadék-csatolású Jordan-metrikák ($\tilde{g}_{\mu\nu}$), valamint a V19–V23 modellágak taxonómiája
# Step 2: Conformal, Disformal, and Effective Fluid-Coupled Jordan Metrics ($\tilde{g}_{\mu\nu}$) and the Taxonomy of V19–V23 Model Branches

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step2_conformal_disformal.md`  
**Előzmény / Previous Step:** `step1_dilaton_form.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` magyar kézirat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1., 1.1., 1.2. és 3. szakasz, Eq. 1–12; DOI: `10.5281/zenodo.23639084`)  
- V23 (`50aa2db0` angol kézirat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex`  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` és `hu_p2.tex` (DOI: `10.5281/zenodo.23604423`)  
- V23 SymPy ellenőrző kód: `covariant_closure_50aa2db0.py` (105–138. sor)

---

## 1. Szó szerinti idézetek a V23 kéziratból (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`)

### 1.1. A kovariáns hatás felírása (V23 EN/HU 1. szakasz, Eq. (1), 33–35. sor)
A V23 kézirat 1. szakasza (*„Kovariáns hatáselv: két explicit megfogalmazás”*) a $(-+++)$ metrikus előjel-konvencióval az alábbi kovariáns hatást rögzíti:
$$S\bigl[g_{\mu\nu},\,\phi_L,\,\Psi_{\text{matter}}\bigr] \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}\,R(g) \;+\; \mathcal{L}_{\text{matter}}\!\bigl(\Psi_{\text{matter}},\,\tilde{g}_{\mu\nu}(g,\phi_L)\bigr) \;+\; \mathcal{L}_{\phi}\!\bigl(g_{\mu\nu},\,\phi_L\bigr)\right]. \tag{V23 Eq. 1}$$

### 1.2. A mikroszkopikus $H_0$ vs. $H_1$ csatolás és a két effektív Jordan-metrika szó szerinti definíciója (V23 3. szakasz, Eq. (11)–(12), 95–101. sor)
A V23 kézirat 3. szakasza (*„Mikroszkopikus óra- és színképvonal-csatolás a hatásból ($H_0$ vs. $H_1$)”*) szó szerint így definiálja a két csatolási ágat és a két Jordan-metrikát:
> *„Az (1) kovariáns hatásban az atomi órák és az elektromágneses színképvonalak mikroszkopikus csatolását egyértelműen meghatározza az, hogy az $\mathcal{L}_{\text{atom}}(\psi_e, A_\mu)$ Dirac + Maxwell Lagrangiánus melyik metrikával húzza össze az indexeket:*
> - **Minimális fizikai metrikás csatolás ($H_0$ — alapértelmezett kovariáns lezárás):** *Az $\mathcal{L}_{\text{atom}}(\psi_e, A_\mu; g_{\mu\nu})$ közvetlenül a fizikai $g_{\mu\nu}$ téridő-metrikához csatolódik, miközben a makroszkopikus barionos energia–impulzus forrást az $f_L$ súlyozza az I. vagy II. megfogalmazáson keresztül. Ekkor egy helyhez kötött atom sajátideje $d\tau_{H_0} = \sqrt{-g_{00}}\,dt = e^{\nu(r)}dt$, a lokális átmeneti frekvenciák $\omega_{\text{atom}} = m_e c^2 \alpha_{\text{EM}}^2/(2\hbar)$ téridő-állandók, és az atomi órák pontosan az $f_L$-lel módosított Einstein-egyenletekből származó metrikus lapse-függvényt mérik. A makroszkopikus óra-anomália kizárólag a módosított tömegből eredő $\Delta\Phi_N$ gravitációs potenciálváltozáson keresztül jelenik meg:*
>   $$\left.\frac{\Delta\nu}{\nu}\right|_{H_0} \;=\; \frac{\Delta\Phi_{\text{metric}}}{c^2}, \qquad \Delta z_{H_0}(\text{Szíriusz B}) \;=\; z_{\text{GR}}\!\left(\frac{\Delta M_B}{M_B}\right)_{\text{TAL}}. \tag{V23 Eq. 11}$$
> - **Diszform / konform sajátidő-csatolás ($H_1$ — közvetlen óracsatolási hipotézis):** *Az $\mathcal{L}_{\text{atom}}(\psi_e, A_\mu; \tilde{g}_{\mu\nu})$ az effektív*
>   $$\tilde{g}_{\mu\nu}^{\text{disf}} \;=\; g_{\mu\nu} \;+\; \frac{1 - f_L(\phi_L)^2}{c^2}\,u_\mu u_\nu \qquad\text{vagy}\qquad \tilde{g}_{\mu\nu}^{\text{conf}} \;=\; f_L(\phi_L)^2\,g_{\mu\nu} \tag{V23 Eq. 12}$$
>   *(rögzített $m_e$-vel) metrikához csatolódik, így a sajátidő-elem $d\tilde{\tau}_{H_1} = \sqrt{-\tilde{g}_{00}}\,dt = f_L(\mathbf{x})\,e^{\nu(\mathbf{x})}dt$. Ez közvetlen nem-metrikus frekvenciaeltolódást indukál:*
>   $$\left.\frac{\Delta\nu}{\nu}\right|_{H_1}(\mathbf{x}) \;=\; \frac{\Delta\Phi_{\text{metric}}(\mathbf{x})}{c^2} \;+\; \bigl(f_L(\mathbf{x}) - 1\bigr) \;\simeq\; \frac{\Delta\Phi_{\text{metric}}(\mathbf{x})}{c^2} \;-\; \phi_L(\mathbf{x}).\text{”}$$

---

## 2. A konform ($\tilde{g}_{\mu\nu}^{\text{conf}}$) transzformáció levezetése és a Brans–Dicke-kapcsolat

Az 1. lépésben (`step1_dilaton_form.md`) igazoltuk, hogy $f_L(\phi_L) = e^{-\phi_L}$. Behelyettesítve ezt a V23 Eq. (12) második egyenletébe:
$$\boxed{\;\tilde{g}_{\mu\nu}^{\text{conf}}(x) \;=\; f_L(\phi_L(x))^2\,g_{\mu\nu}(x) \;=\; e^{-2\phi_L(x)}\,g_{\mu\nu}(x),\qquad \Omega(\phi_L) \;\equiv\; f_L(\phi_L) \;=\; e^{-\phi_L}.\;}$$

### 2.1. Inverz metrika, determináns és sajátidő-elem
1. **Inverz konform metrika:** Mivel $\tilde{g}_{\mu\alpha}^{\text{conf}}\tilde{g}_{\text{conf}}^{\alpha\nu} = \delta_\mu^\nu$,
   $$\tilde{g}_{\text{conf}}^{\mu\nu} \;=\; f_L(\phi_L)^{-2}\,g^{\mu\nu} \;=\; e^{+2\phi_L}\,g^{\mu\nu}.$$
2. **Négyes térfogatelem (determináns):** $D=4$ téridő-dimenzióban $\det(e^{-2\phi_L}g_{\mu\nu}) = (e^{-2\phi_L})^4 \det(g_{\mu\nu})$, tehát:
   $$\sqrt{-\tilde{g}_{\text{conf}}} \;=\; f_L(\phi_L)^4\,\sqrt{-g} \;=\; e^{-4\phi_L}\,\sqrt{-g}.$$
3. **Statikus megfigyelő sajátidő-eleme ($dx^i = 0$):** Statikus gömbszimmetrikus fizikai metrikában ($g_{00} = -e^{2\nu(r)}$):
   $$d\tilde{\tau}_{\text{conf}}^2 \;=\; -\frac{1}{c^2}\tilde{g}_{00}^{\text{conf}}\,dt^2 \;=\; -f_L^2\frac{g_{00}}{c^2}\,dt^2 \;=\; f_L(r)^2\,e^{2\nu(r)}\,dt^2 \quad\implies\quad \boxed{\;d\tilde{\tau}_{\text{conf}} \;=\; f_L(r)\,e^{\nu(r)}\,dt \;=\; e^{\nu(r) - \phi_L(r)}\,dt.\;}$$

### 2.2. A Ricci-skalár konform transzformációja (Jordan-keret $\leftrightarrow$ Einstein-keret)
Standard $D=4$ konform geometriából ismert (pl. Wald 1984, Appendix D; Fujii & Maeda 2003), hogy ha $g_{\mu\nu} = \Omega^{-2}\tilde{g}_{\mu\nu}^{\text{conf}} = e^{+2\phi_L}\tilde{g}_{\mu\nu}^{\text{conf}}$ (azaz $\tilde{g}_{\mu\nu}^{\text{conf}} = \Omega^2 g_{\mu\nu}$, ahol $\Omega = e^{-\phi_L}$), akkor a két Ricci-skalár kapcsolata:
$$R(g) \;=\; \Omega^2\left[\tilde{R}(\tilde{g}_{\text{conf}}) \;+\; 6\,\tilde{\Box}\ln\Omega \;-\; 6\,\tilde{g}_{\text{conf}}^{\mu\nu}\,\partial_\mu\ln\Omega\,\partial_\nu\ln\Omega\right].$$
Ennek megfelelően az Einstein–Hilbert-tag átírása Jordan-keretbe ($\sqrt{-g} = \Omega^{-4}\sqrt{-\tilde{g}_{\text{conf}}}$):
$$\int d^4x\,\sqrt{-g}\,\frac{c^4}{16\pi G}\,R(g) \;=\; \int d^4x\,\sqrt{-\tilde{g}_{\text{conf}}}\,\frac{c^4}{16\pi G}\left[\Omega^{-2}\tilde{R}(\tilde{g}_{\text{conf}}) \;-\; 6\,\Omega^{-4}\tilde{g}_{\text{conf}}^{\mu\nu}\partial_\mu\Omega\,\partial_\nu\Omega\right],$$
ahol a $\nabla_\mu(\sqrt{-g}\,6\Omega\nabla^\mu\Omega)$ határdivergencia-tagot parciális integrálással elhagytuk. Bevezetve a **Brans–Dicke skalármezőt**:
$$\boxed{\;\Phi_{\text{BD}}(x) \;\equiv\; \Omega(x)^{-2} \;=\; f_L(x)^{-2} \;=\; e^{+2\phi_L(x)},\;}$$
amelyre $\frac{\partial_\mu\Phi_{\text{BD}}}{\Phi_{\text{BD}}} = +2\partial_\mu\phi_L$, a V23 II. megfogalmazásának lokális skalár-kinetikus tagja ($\mathcal{L}_\phi = -\frac{c^4}{8\pi G\tilde{\alpha}_L}\bigl[\frac{1}{2}g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + \frac{1}{2L^2}\phi_L^2\bigr]$, V23 Eq. 8) Jordan-keretben:
$$-\frac{1}{2}g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L \;\longrightarrow\; -\frac{\omega_{\text{BD}}}{\Phi_{\text{BD}}}\,\tilde{g}_{\text{conf}}^{\mu\nu}\partial_\mu\Phi_{\text{BD}}\,\partial_\nu\Phi_{\text{BD}}.$$
Az együtthatók összevetése adja a standard skalár–tenzor azonosságot:
$$2\omega_{\text{BD}} + 3 \;=\; \frac{1}{\tilde{\alpha}_L} \;=\; \frac{1}{2\beta^2} \quad\Longleftrightarrow\quad \boxed{\;\omega_{\text{BD}} \;=\; \frac{1}{2\tilde{\alpha}_L} - \frac{3}{2} \;=\; \frac{1}{4\beta^2} - \frac{3}{2}.\;}$$
*(A $\tilde{\alpha}_L$, $\alpha_Y = 8\pi\tilde{\alpha}_L$, $\alpha_{\text{std}} = 2\tilde{\alpha}_L = \frac{\alpha_Y}{4\pi}$ és $\omega_{\text{BD}}$ pontos számértékeit és a Cassini PPN-korlátot az 5. lépésben elemezzük részletesen.)*

---

## 3. A diszform ($\tilde{g}_{\mu\nu}^{\text{disf}}$) csatolás és az algebrai projektorszétbontás

A V23 Eq. (12) első egyenlete egy tisztán időirányú **Bekenstein-típusú diszform metrikát** (Bekenstein 1993, *Phys. Rev. D* 48, 3641) definiál a közeg $u^\mu$ négyessebességével ($g_{\mu\nu}u^\mu u^\nu = -c^2$, $u_\mu = g_{\mu\nu}u^\nu$):
$$\boxed{\;\tilde{g}_{\mu\nu}^{\text{disf}} \;=\; g_{\mu\nu} \;+\; \frac{1 - f_L(\phi_L)^2}{c^2}\,u_\mu u_\nu \;=\; g_{\mu\nu} \;+\; \frac{1 - e^{-2\phi_L}}{c^2}\,u_\mu u_\nu.\;}$$

### 3.1. Egzakt projektorszétbontás (SymPy-ellenőrzött azonosság, `covariant_closure_50aa2db0.py` 108–125. sor)
Vezessük be a téridő-metrika $u^\mu$-re merőleges (térbeli) projektorát:
$$h_{\mu\nu} \;\equiv\; g_{\mu\nu} + \frac{1}{c^2}u_\mu u_\nu \qquad\left(h_{\mu\nu}u^\nu = 0,\quad h^\mu_{\;\alpha}h^\alpha_{\;\nu} = h^\mu_{\;\nu},\quad h^\mu_{\;\mu} = 3\right).$$
Ekkor $g_{\mu\nu} = h_{\mu\nu} - \frac{1}{c^2}u_\mu u_\nu$, így a diszform metrika szeparált projektor-alakja:
$$\boxed{\;\tilde{g}_{\mu\nu}^{\text{disf}} \;=\; h_{\mu\nu} \;-\; f_L(\phi_L)^2\,\frac{u_\mu u_\nu}{c^2} \;=\; \left(g_{\mu\nu} + \frac{u_\mu u_\nu}{c^2}\right) \;-\; e^{-2\phi_L}\,\frac{u_\mu u_\nu}{c^2}.\;}$$
Ebből azonnal látszik a fizikai jelentés:
1. **A térbeli geometriát a diszform csatolás érintetlenül hagyja:** minden $u^\mu$-re merőleges térbeli $\delta x^\mu$ elmozdulásra ($u_\mu \delta x^\mu = 0$) $\tilde{g}_{\mu\nu}^{\text{disf}}\delta x^\mu\delta x^\nu = g_{\mu\nu}\delta x^\mu\delta x^\nu$.
2. **Az időirányú (lapse) komponenst pontosan $f_L(\phi_L) = e^{-\phi_L}$ faktorral skálázza:** egy a közeggel együttmozgó óra ($\delta x^\mu = \frac{u^\mu}{c}d\tau$) sajátideje:
   $$d\tilde{\tau}_{\text{disf}}^2 \;=\; -\frac{1}{c^2}\tilde{g}_{\mu\nu}^{\text{disf}}\,dx^\mu dx^\nu \;=\; f_L(\phi_L)^2\left(-\frac{1}{c^2}g_{\mu\nu}dx^\mu dx^\nu\right) \;=\; f_L(\phi_L)^2\,d\tau^2.$$
3. **Inverz diszform metrika és determináns:** Mivel $h^{\mu\alpha}u_\alpha = 0$ és $\bigl(-\frac{u_\mu u_\alpha}{c^2}\bigr)\bigl(-\frac{u^\alpha u^\nu}{c^2}\bigr) = -\frac{u_\mu u^\nu}{c^2}$, az inverz metrika és a térfogatelem egzaktul:
   $$\boxed{\;\tilde{g}_{\text{disf}}^{\mu\nu} \;=\; h^{\mu\nu} \;-\; f_L(\phi_L)^{-2}\,\frac{u^\mu u^\nu}{c^2} \;=\; g^{\mu\nu} \;-\; \frac{f_L(\phi_L)^{-2} - 1}{c^2}\,u^\mu u^\nu,\qquad \sqrt{-\tilde{g}_{\text{disf}}} \;=\; f_L(\phi_L)\,\sqrt{-g} \;=\; e^{-\phi_L}\sqrt{-g}.\;}$$

### 3.2. Általános statikus gömbszimmetrikus háttér és a `SymPy` ellenőrzés
Statikus gömbszimmetrikus fizikai háttérmetrikában:
$$ds^2 \;=\; g_{\mu\nu}dx^\mu dx^\nu \;=\; -e^{2\nu(r)}c^2 dt^2 + e^{2\lambda(r)}dr^2 + r^2 d\Omega^2, \qquad g_{00} = -e^{2\nu(r)},$$
egy nyugvó (együttmozgó) megfigyelő négyessebessége $u^\mu = \bigl(c/\sqrt{-g_{00}},\,0,\,0,\,0\bigr) = \bigl(c\,e^{-\nu(r)},\,0,\,0,\,0\bigr)$, kovariáns alakja pedig $u_\mu = g_{\mu\nu}u^\nu = \bigl(-\sqrt{-g_{00}}\,c,\,0,\,0,\,0\bigr) = \bigl(-c\,e^{\nu(r)},\,0,\,0,\,0\bigr)$. Így:
$$\frac{u_0 u_0}{c^2} \;=\; -g_{00} \;=\; +e^{2\nu(r)},$$
és ennek behelyettesítésével a diszform metrika $(0,0)$-komponense:
$$\tilde{g}_{00}^{\text{disf}} \;=\; g_{00} \;+\; \bigl(1 - f_L(r)^2\bigr)\frac{u_0 u_0}{c^2} \;=\; g_{00} \;+\; \bigl(1 - f_L(r)^2\bigr)\bigl(-g_{00}\bigr) \;=\; f_L(r)^2\,g_{00} \;=\; -f_L(r)^2\,e^{2\nu(r)}.$$
A `covariant_closure_50aa2db0.py` (108–125. sor) a lokalizált ortonormált tetrádban ($g_{\hat{a}\hat{b}} = \eta_{\hat{a}\hat{b}} = \text{diag}(-1,1,1,1)$, $u^{\hat{a}} = (c,0,0,0)$, $u_{\hat{a}} = (-c,0,0,0)$) ugyanezt az azonosságot ellenőrzi: $\tilde{g}_{\hat{0}\hat{0}}^{\text{disf}} = -1 + (1 - f_L^2) = -f_L^2$, így az általános koordinátás alak $\tilde{g}_{00}^{\text{disf}} = e_{\;0}^{\hat{0}} e_{\;0}^{\hat{0}} \tilde{g}_{\hat{0}\hat{0}}^{\text{disf}} = e^{2\nu}(-f_L^2) = f_L^2 g_{00}$.

---

## 4. A V19–V23 modellágak explicit taxonómiája (Branch Taxonomy Table)

A bírálati átláthatóság érdekében az alábbi táblázat egyértelműen szétválasztja a V19–V23 kéziratokban és a jelen levezetési láncban szereplő modellágakat, azok hatáselvét/definícióját, térfogatelemét, hidrosztatikai egyenletét, óra-/fénycsatolását és érvényességi körét:

| Modellág kódja | Megnevezés és kézirat-forrás | Skalár / mag-definíció és 2. momentum | Anyagi csatolás / térfogatelem | Belső hidrosztatikai egyensúly | Elektromágneses / óra / PPN viselkedés |
|---|---|---|---|---|---|
| **Branch I-A** (`tophat_mass`, `tophat_gradient`) | V21 fenomenologikus 3D gömbi top-hat TOV (`stellar_tal_v21.py`, V23 Eq. 4–5) | Nemlokális 3D gömbi ablak: $\phi_L = \gamma_L(K_L*\varepsilon)$, $K_L = \frac{3}{4\pi L^3}\Theta(L-r)$, $1 + \frac{L^2}{10}\nabla^2 + \dots$ | Fenomenologikus $f_L$ szorzó $m'$-ben (`mass`) vagy $\nu'$-ben (`gradient`), algebrai $T^\mu_{\;\nu}(X)$ lezárással | $h'(r) = -\nu'(r)$ (`mass`: $\nu'=A_{\text{geo}}$; `gradient`: $\nu'=f_L A_{\text{geo}}$) | $H_0$: minimális fizikai metrika; $H_1$: közvetlen $f_L$ óraeltolódás (Szíriusz B kizárja $L\ge 1\;\text{km}$-re) |
| **Branch I-B** (`bilateral_eom`, `bilateral_lapse_weighted`) | V23 I. megfogalmazás: bilaterális nemlokális mag-hatás (V23 Eq. 2–3, 6–7) | Szimmetrikus kétpont-mag $K_L(x,y) = K_L(y,x)$, $1 + \frac{L^2}{10}\nabla^2 + \dots$ | $S_{\text{I},m} = -\int d^4x\sqrt{-g}\,e^{-\phi_L}\varepsilon$ (fix metrikájú euklideszi mag-közelítésben) | $V_{\text{eff}}(r) = f_L(r) - \gamma_L\int K_L f_L\varepsilon$ reciprocitási erővel | Makroszkopikus hidrosztatikai tesztág; a metrika-variáció a magban nincs lezárva |
| **Branch II-Fluid** (`yukawa_action_el`) | V23 II. megfogalmazás: lokális Yukawa-árnyékolt Einstein–skalár–folyadék hatás (V23 Eq. 8–10) | Lokális Helmholtz–Yukawa PDE: $(-\Box_g + L^{-2})\phi_L = \frac{\alpha_Y G}{c^2}e^{-\phi_L}\varepsilon$, $1 + L^2\nabla^2 + \dots$ | $S_m = -\int d^4x\sqrt{-g}\,e^{-\phi_L}\varepsilon(n)$, $P_{\text{matter}} = e^{-\phi_L}P$, $\rho_\phi = \frac{B(\phi_L')^2 + \phi_L^2/L^2}{2\alpha_Y}$ | $h'(r) = -\nu'(r) + \phi_L'(r)$ (ötödik erő + görbült külső skalárfarok $M_{\text{ADM}}$) | $H_0$ alatt az atomi órák a fizikai $g_{\mu\nu}$-höz csatolódnak; $M_{\text{ADM}}$ és $M_B$ fordulópontja egybeesik |
| **Branch II-Conf** ($\tilde{g}_{\mu\nu}^{\text{conf}} = e^{-2\phi_L}g_{\mu\nu}$) | Univerzális konform Jordan–Brans–Dicke / dilaton metrika (V23 Eq. 12b, `step1`, `step5`) | Ugyanaz a lokális Yukawa PDE; $\Phi_{\text{BD}} = e^{+2\phi_L}$ | $\sqrt{-\tilde{g}_{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$; konform invariáns Maxwell-hatás ($\tilde{F}_{\mu\nu}\tilde{F}^{\mu\nu}\sqrt{-\tilde{g}} = F_{\mu\nu}F^{\mu\nu}\sqrt{-g}$) | Nemrelativisztikus határesetben ($\varepsilon \simeq \rho_0 c^2 \gg 3P$) megegyezik a $\beta = \sqrt{\alpha_{\text{std}}/2}$ ötödik erővel | Erהre az ágra érvényes a standard konform Brans–Dicke / PPN szótár: $\alpha_{\text{std}} = \frac{\alpha_Y}{4\pi}$, $\omega_{\text{BD}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2} \approx 11{,}07$ ($\alpha_Y=1$), fényelhajlás konform invariáns |
| **Branch II-Disf** ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$) | Együttmozgó időirányú diszform óra- és EM-metrika (V23 Eq. 12a, V20/V22 szolenoid, `step14`–`step15`) | Ugyanaz a lokális Yukawa PDE, mágneses $\rho_B = \frac{|\mathbf{B}|^2}{2\mu_0 c^2}$ forrással kiegészítve | $\sqrt{-\tilde{g}_{\text{disf}}} = f_L\sqrt{-g} = e^{-\phi_L}\sqrt{-g}$ (egybeesik a Branch II-Fluid $e^{-\phi_L}$ súlyával!) | Ha $u^\mu$-t együttmozgó khronon/éter mező rögzíti, az energia-súly $e^{-\phi_L}$; szabad $u^\mu$-variáció extra kényszert ad | Törésmutató $n_{\text{refr}} = e^{+\phi_L}$ (refraktív fényelhajlás), mágneses szolenoid-forrás ($\rho_B$); **nem** a konform Brans–Dicke PPN-képlet írja le |
| **Branch II-Cham** ($\beta_{\text{eff}}(\rho)$ kiterjesztés) | Sűrűségfüggő kaméleon / szimmetron nemlineáris EFT-hipotézis (`step11`) | Sűrűségfüggő effektív tömeg $m_{\text{eff}}^2(\rho)$ és csatolás $\beta_{\text{eff}}(\rho)$ | Nemlineáris $V(\phi) + A(\phi)\rho$ kiterjesztés a kvadratikus konstans-$\alpha_Y$ V23 hatáson túl | Környezetfüggő ötödik erő (erős árnyékolás laborban/NS-ben, erősítés galaktikus ritka közegben) | **Önálló továbbfejlesztési hipotézis**, nem a konstans $\alpha_Y$-ú V23 `yukawa_action_el` egyenletek következménye |

---

## 5. Ellenőrzések (Checks)

1. **Határérték-ellenőrzés ($\phi_L \to 0$, azaz $f_L \to 1$):**
   - Ha $\phi_L = 0$ ($f_L = 1$), akkor $\tilde{g}_{\mu\nu}^{\text{conf}} = g_{\mu\nu}$ és $\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + 0 = g_{\mu\nu}$, tehát mindkét Jordan-metrika visszatér a fizikai $g_{\mu\nu}$ Einstein-metrikához, és $d\tilde{\tau}_{H_1} = d\tau_{H_0}$ $\checkmark$.
2. **Sajátidő-egyezés nyugvó órákra ($u^\mu$ mentén):**
   - Egy a barionos közeggel együttmozgó órára ($dx^\mu \parallel u^\mu$) a konform és a diszform metrika **azonos** sajátidő-elemet és frekvenciaeltolódást ad:
     $$d\tilde{\tau}_{\text{conf}} \;=\; d\tilde{\tau}_{\text{disf}} \;=\; f_L(x)\,d\tau \;=\; e^{-\phi_L(x)}\,\sqrt{-g_{00}}\,dt,$$
     amit a `covariant_closure_50aa2db0.py` (`sym_res["proper_time_ratio_H1_over_H0"] = "f_L"`, `sym_res["weak_field_H1_clock_shift"] = "Phi_N/c^2 - phi_L"`) szimbolikusan is igazol $\checkmark$.
3. **Különbség a térbeli szektorban, a térfogatelemben és a fényterjedésben:**
   - Míg $\sqrt{-\tilde{g}_{\text{conf}}} = e^{-4\phi_L}\sqrt{-g}$ és a konform metrika a nullgeodetikusokat (fénysugarakat) változatlanul hagyja ($\tilde{g}_{\mu\nu}^{\text{conf}}k^\mu k^\nu = 0 \iff g_{\mu\nu}k^\mu k^\nu = 0$), addig a diszform metrikára $\sqrt{-\tilde{g}_{\text{disf}}} = e^{-\phi_L}\sqrt{-g}$, és a fényre effektív törésmutatót ($n = f_L^{-1} = e^{+\phi_L}$) indukál, miközben a térbeli távolságokat nem skálázza $\checkmark$.
