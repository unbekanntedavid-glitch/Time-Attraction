# 4. lépés: Kovariáns energia–impulzus-megmaradás, skalár–anyag impulzuscsere és a módosított hidrosztatikai egyensúly ($h' = -\nu' + \phi_L'$)
# Step 4: Covariant Stress-Energy Conservation, Scalar-Matter Momentum Exchange, and Modified Hydrostatic Equilibrium ($h' = -\nu' + \phi_L'$)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step4_covariant_conservation_hydrostatic.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md`, `step2_conformal_disformal.md`, `step3_yukawa_helmholtz_eom.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` magyar kézirat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.1. és 1.2. alfejezet, 47–65. sor, Eq. 4–5 és Eq. 6–7; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (2. fejezet, 205–207. sor; DOI: `10.5281/zenodo.23604423`)  
- V23 SymPy ellenőrző kód: `covariant_closure_50aa2db0.py` (72–98. és 143–160. sor) és `conservation_Residuals.json`

---

## 1. Szó szerinti idézetek a V22 és V23 kéziratokból

A 3. lépésben (`step3_yukawa_helmholtz_eom.md`) levezettük a $\phi_L$ skalár-/dilatonmező Euler–Lagrange-egyenletét. A következő szükségszerű elméleti lépés a $g^{\mu\nu}$ metrika szerinti variáció és a **Bianchi-azonosságból ($\nabla_\mu G^{\mu\nu} \equiv 0$) következő kovariáns energia–impulzus-megmaradás** ($\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_X) = 0$), valamint az ebből adódó **módosított hidrosztatikai egyensúly (ötödik erő / skalár-gradiens erő)** levezetése.

1. **A fizikai probléma megfogalmazása a V22 monográfiában (`hu_p1.tex`, 206. sor):**
   > *„Ugyanakkor az $f_L(\mathbf{x})$ külső szorzóként való beillesztése a Tolman–Oppenheimer–Volkoff-egyenletekbe vagy az órafrekvencia-képletekbe $\nabla_\mu T^{\mu\nu}_{\text{matter}} \neq 0$ maradékot hagy maga után, ha nem egyensúlyozza ki egy dinamikus mezőszektor. Ezért a dinamikai ‚miért’-et az `50aa2db0` feladatként fogalmaztam meg: a nemlokális mag olyan kovariáns skalár–tenzor vagy nemlokális $S = \int d^4x\sqrt{-g}\,\bigl[\frac{c^4}{16\pi G}R + \mathcal{L}_{\text{matter}}(\psi, \tilde{g}_{\mu\nu}) + \mathcal{L}_\phi(\phi,\nabla\phi)\bigr]$ hatáselvbe emelése, amelynek Euler–Lagrange-egyenlete kvázisztatikus határesetben visszaadja a $\phi(\mathbf{x}) \to \gamma_L\bar{\rho}_L(\mathbf{x})$ profilt, miközben garantálja a $\nabla_\mu(T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_X) = 0$ megmaradási törvényt.”*

2. **Az I. megfogalmazásbeli algebrai $X$-szektor rekonstrukció a V23-ban (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 47–52. sor, Eq. 4–5):**
   > *„Statikus gömbszimmetrikus belső térben, ha a $T^\mu_{\;\nu}(X) = \text{diag}(-\varepsilon_X, P_{r,X}, P_{t,X}, P_{t,X})$ segédszektor veszi fel az ablak-feszültséget a $G^\mu_{\;\nu} = \frac{8\pi G}{c^4}(T^\mu_{\;\nu,\text{matter}} + T^\mu_{\;\nu,X})$ egyenlet jobb oldalán, míg a folyadékelemek kielégítik a $\nabla_\mu T^{\mu\nu}_{\text{matter}} = 0$ egyenletet ($P' = -(\varepsilon+P)\nu'$ geometrizált egységekben, $G=c=1$, $\nu = \Phi_{\text{eff}}/c^2$, $g_{tt}=-e^{2\nu}=-e^{2\Phi_{\text{eff}}/c^2}$), akkor a $\nabla_\mu G^\mu_{\;r} \equiv 0$ Bianchi-azonosság mindkét V21 módban egyértelműen rögzíti $T^\mu_{\;\nu}(X)$-et:”*
   $$\begin{aligned}
     \textbf{Tömeg-mód } (m' = 4\pi r^2 f_L\varepsilon,\; \nu' = A_{\text{geo}}): &\quad \varepsilon_X = -(1-f_L)\varepsilon,\quad P_{r,X} = 0,\quad P_{t,X} = -\frac{r}{2}(1-f_L)\varepsilon\nu', \tag{V23 Eq. 4}\\
     \textbf{Gradiens-mód } (m' = 4\pi r^2\varepsilon,\; \nu' = f_L A_{\text{geo}}): &\quad \varepsilon_X = 0,\quad P_{r,X} = -(1-f_L)\!\left(P + \frac{m}{4\pi r^3}\right),\; P_{t,X} = P_{r,X} + \frac{r}{2}\bigl(P'_{r,X} + P_{r,X}\nu'\bigr), \tag{V23 Eq. 5}
   \end{aligned}$$
   > *„ahol $A_{\text{geo}}(r) = (m + 4\pi r^3 P)/[r(r-2m)]$. A `SymPy`-ban elvégzett szimbolikus ellenőrzés (`conservation_Residuals.json`) igazolja, hogy $\nabla_\mu T^\mu_{\;r}(\text{matter}) \equiv 0$ és $\nabla_\mu T^\mu_{\;r}(X) \equiv 0$ azonosságszerűen fennáll.”*

3. **A II. megfogalmazásbeli lokális skalár–tenzor (Yukawa) feszültség–energia tenzor és a hatással konzisztens hidrosztatikai egyensúly (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 65. sor):**
   > *„a $\delta S_\phi/\delta g^{\mu\nu}$ variáció pedig a lokális kovariáns $T_{\mu\nu}^{(\phi)} = \frac{c^4}{8\pi G\tilde{\alpha}_L}\bigl[\nabla_\mu\phi_L\nabla_\nu\phi_L - \frac{1}{2}g_{\mu\nu}\bigl((\nabla\phi_L)^2 + L^{-2}\phi_L^2\bigr)\bigr]$ feszültség-energia tenzort szolgáltatja. [...] A `yukawa_action_el` sorok az $\mathcal{L}=-[(\nabla\phi)^2+\phi^2/L^2]/(2\alpha_Y)-f_L\varepsilon$ hatásból következnek ($G=c=1$; `SymPy`: `sym_yukawa_el_source_factor` és `_hydrostatic` maradék $=0$): $B\phi''+(2B/r+B\nu'+B'/2)\phi'-\phi/L^2=-\alpha_Y f_L\varepsilon$, $B=1-2m/r$, $P_{\text{matter}}=f_LP$, $h'=-\nu'+\phi'$, $\rho_\phi=(B\phi'^2+\phi^2/L^2)/(2\alpha_Y)$, külső illesztés $\phi'(R)=-(1/R+1/L)\phi(R)$.”*

---

## 2. A skalár–anyag impulzuscsere és a módosított hidrosztatikai egyensúly levezetése

### 2.1. Az effektív anyagnyomás ($P_{\text{matter}} = f_L P$) levezetése
A V23 `yukawa_action_el` ágában (és általános Einstein-keretbeli skalár–tenzor elméletben, ahol az anyagi Lagrangiánus $\mathcal{L}_{\text{fluid}}(n,\phi_L) = -f_L(\phi_L)\,\varepsilon(n) = -e^{-\phi_L}\varepsilon(n)$, $n$ pedig a barion-részecskeszámsűrűség) a termodinamikai nyomást a barionszám-megmaradás melletti Legendre-transzformáció definiálja:
$$P_{\text{matter}} \;=\; n\,\frac{\partial\bigl(f_L(\phi_L)\,\varepsilon(n)\bigr)}{\partial n} \;-\; f_L(\phi_L)\,\varepsilon(n) \;=\; f_L(\phi_L)\Bigl(n\,\frac{d\varepsilon}{dn} - \varepsilon(n)\Bigr) \;=\; f_L(\phi_L)\,P(n).$$
Így az $f_L$-lel súlyozott ideális folyadék feszültség–energia tenzora:
$$T^{\mu\nu}_{\text{matter}} \;=\; f_L(\phi_L)\,\Bigl[\bigl(\varepsilon + P\bigr)\frac{u^\mu u^\nu}{c^2} + P\,g^{\mu\nu}\Bigr] \quad\implies\quad \varepsilon_{\text{matter}} = f_L\varepsilon, \qquad P_{\text{matter}} = f_L P.$$

### 2.2. A skalármező feszültség–energia tenzorának ($\nabla^\mu T_{\mu\nu}^{(\phi)}$) divergenciája
Geometrizált egységekben ($G = c = 1$, ahol $\alpha_Y = 8\pi\tilde{\alpha}_L$, így $\frac{c^4}{8\pi G\tilde{\alpha}_L} = \frac{1}{\alpha_Y}$) a skalármező feszültség–energia tenzora:
$$T_{\mu\nu}^{(\phi)} \;=\; \frac{1}{\alpha_Y}\left[\nabla_\mu\phi_L\nabla_\nu\phi_L - \frac{1}{2}g_{\mu\nu}\left(g^{\alpha\beta}\nabla_\alpha\phi_L\nabla_\beta\phi_L + \frac{1}{L^2}\phi_L^2\right)\right].$$
Számítsuk ki a kovariáns divergenciáját ($\nabla^\mu T_{\mu\nu}^{(\phi)}$):
$$\nabla^\mu T_{\mu\nu}^{(\phi)} \;=\; \frac{1}{\alpha_Y}\left[\left(\Box_g\phi_L\right)\nabla_\nu\phi_L + \left(\nabla^\mu\phi_L\right)\nabla_\mu\nabla_\nu\phi_L - \left(\nabla^\alpha\phi_L\right)\nabla_\nu\nabla_\alpha\phi_L - \frac{1}{L^2}\phi_L\nabla_\nu\phi_L\right].$$
Mivel skalármezőre a második kovariáns deriváltak felcserélhetők ($\nabla_\mu\nabla_\nu\phi_L = \nabla_\nu\nabla_\mu\phi_L$), a második és harmadik tag kiesik:
$$\nabla^\mu T_{\mu\nu}^{(\phi)} \;=\; \frac{1}{\alpha_Y}\left(\Box_g\phi_L - \frac{1}{L^2}\phi_L\right)\nabla_\nu\phi_L.$$
Most helyettesítsük be a 3. lépésben (`step3_yukawa_helmholtz_eom.md`) levezetett egzakt `yukawa_action_el` Euler–Lagrange-mezőegyenletet:
$$\left(-\Box_g + \frac{1}{L^2}\right)\phi_L \;=\; \alpha_Y\,f_L(\phi_L)\,\varepsilon \quad\Longleftrightarrow\quad \Box_g\phi_L - \frac{1}{L^2}\phi_L \;=\; -\alpha_Y\,f_L\,\varepsilon.$$
Ezzel a skalárszektor divergenciája zárt alakban:
$$\boxed{\;\nabla^\mu T_{\mu\nu}^{(\phi)} \;=\; -f_L(\phi_L)\,\varepsilon\,\nabla_\nu\phi_L \;=\; \varepsilon\,\nabla_\nu f_L(\phi_L).\;}$$

### 2.3. Az együttes megmaradás ($\nabla^\mu(T_{\mu\nu}^{\text{matter}} + T_{\mu\nu}^{(\phi)}) = 0$) és a módosított Euler-/TOV-egyenlet
Az Einstein-egyenletek ($G_{\mu\nu} = 8\pi(T_{\mu\nu}^{\text{matter}} + T_{\mu\nu}^{(\phi)})$) bal oldalára a kontraktált Bianchi-azonosság ($\nabla^\mu G_{\mu\nu} \equiv 0$) miatt a jobb oldal együttes divergenciája azonosan eltűnik:
$$\nabla^\mu T_{\mu\nu}^{\text{matter}} \;+\; \nabla^\mu T_{\mu\nu}^{(\phi)} \;=\; 0 \quad\implies\quad \nabla^\mu T_{\mu\nu}^{\text{matter}} \;=\; +f_L(\phi_L)\,\varepsilon\,\nabla_\nu\phi_L.$$
Statikus gömbszimmetrikus téridőben ($ds^2 = -e^{2\nu(r)}dt^2 + B(r)^{-1}dr^2 + r^2 d\Omega^2$, $B(r) = 1 - 2m(r)/r$) az anyagi tenzor radiális ($\nu = r$) kovariáns divergenciája:
$$\nabla_\mu T^\mu_{\;r}(\text{matter}) \;=\; \frac{d}{dr}\bigl(f_L P\bigr) \;+\; f_L(\varepsilon + P)\,\nu' \;=\; f_L\,P' \;-\; f_L\,P\,\phi_L' \;+\; f_L(\varepsilon + P)\,\nu',$$
ahol felhasználtuk, hogy $f_L' = \frac{d}{dr}(e^{-\phi_L}) = -f_L\,\phi_L'$. Összeadva a skalárszektor $\nabla_\mu T^\mu_{\;r}(\phi) = -f_L\,\varepsilon\,\phi_L'$ divergenciájával:
$$\nabla_\mu T^\mu_{\;r}(\text{matter}) + \nabla_\mu T^\mu_{\;r}(\phi) \;=\; f_L\Bigl[P' + (\varepsilon + P)\nu' - (\varepsilon + P)\phi_L'\Bigr] \;=\; 0.$$
Mivel $f_L = e^{-\phi_L} > 0$, eloszthatunk $f_L(\varepsilon + P)$-vel, és bevezetve a fajlagos entalpiát ($dh = \frac{dP}{\varepsilon + P}$):

$$\boxed{\;h'(r) \;\equiv\; \frac{P'(r)}{\varepsilon(r) + P(r)} \;=\; -\nu'(r) \;+\; \phi_L'(r) \;=\; -\nu'(r) \;-\; \frac{d\ln f_L(r)}{dr}.\;}$$

---

## 3. Fizikai értelmezés a skalár–tenzor / dilaton keretrendszerben

A fenti $h' = -\nu' + \phi_L'$ összefüggés a skalár–tenzor gravitáció (Damour & Esposito-Farèse 1992; Fujii & Maeda 2003) és a Time Attraction Law közötti közvetlen dinamikai híd:
1. **Az effektív potenciál szétválása:** A barionos elemre ható eredő egyensúlyi erő két részből áll:
   - a $-\nu'(r) = -\frac{1}{c^2}\frac{d\Phi_{\text{eff}}}{dr}$ szokásos metrikus gravitációs vonzásból, és
   - a $+\phi_L'(r) = -\frac{d\ln f_L}{dr}$ **skalár-/időkoherencia-gradiens erőből** (a skalár–tenzor irodalomban *fifth force* / dilaton-gradiens erő).
2. **Gyenge téri / galaktikus határérték:**
   Nemrelativisztikus közegben ($P \ll \varepsilon \simeq \rho c^2$) egy próbatest vagy folyadékelem effektív gyorsulása:
   $$\mathbf{a}_{\text{eff}} \;=\; -\nabla\Phi_N \;+\; c^2\,\nabla\phi_L \;=\; -\nabla\left(\Phi_N - c^2\phi_L\right) \;=\; -\nabla\left(\Phi_N + c^2\ln f_L\right).$$
   Mivel egy tömegkoncentráció körül a sűrűség befelé nő ($\frac{d\bar{\rho}_L}{dr} < 0$), ezért $\phi_L(r)$ befelé nő ($\phi_L'(r) < 0$, illetve $\frac{d\ln f_L}{dr} > 0$), tehát a $+\phi_L'(r)$ tag a hidrosztatikai egyensúlyban és a pályadinamikában egyaránt a sűrűséggradiens által hajtott járulékot ad.

---

## 4. Ellenőrzések (Checks)

1. **SymPy szimbolikus azonosságok a V23 csomagban (`covariant_closure_50aa2db0.py`, 72–98. és 150–160. sor; `conservation_Residuals.json`):**
   - `sym_yukawa_el_matter_pressure_residual = 0`: igazolja a $P_{\text{matter}} = n\frac{\partial(f_L\varepsilon)}{\partial n} - f_L\varepsilon = f_L P$ azonosságot $\checkmark$.
   - `sym_yukawa_el_hydrostatic_residual = 0`: igazolja, hogy a $\nabla_\mu T^\mu_{\;r}(\text{matter}) + \nabla_\mu T^\mu_{\;r}(\phi) = 0$ egyenletből a Yukawa EL-egyenlet behelyettesítésével egzaktul $h' = -\nu' + \phi_L'$ következik $\checkmark$.
   - `div_matter_mass_symbolic = 0`, `div_X_mass_symbolic = 0`, `div_matter_grad_symbolic = 0`, `div_X_grad_symbolic = 0`: igazolja az I. megfogalmazásbeli V23 Eq. (4)–(5) algebrai $X$-szektor megmaradását $\checkmark$.
2. **Általános relativitáselméleti (GR) határérték ($\phi_L \to 0$, $f_L \to 1$):**
   Amikor $\phi_L \equiv 0$ ($\phi_L' = 0$, $f_L = 1$), a skalár feszültség–energia tenzor eltűnik ($T_{\mu\nu}^{(\phi)} = 0$), az effektív nyomás $P_{\text{matter}} \to P$, a hidrosztatikai egyensúly pedig visszatér a standard GR Tolman–Oppenheimer–Volkoff-egyenlethez ($h' = -\nu'$, azaz $P' = -(\varepsilon + P)\nu'$) $\checkmark$.
