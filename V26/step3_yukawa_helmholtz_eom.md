# 3. lépés: A V23 lokális Helmholtz–Yukawa-hatáselv és mezőegyenlet leképezése a tömeges skalár–tenzor / dilaton Klein–Gordon-egyenletre
# Step 3: Mapping the V23 Local Screened Helmholtz–Yukawa Action and Field Equation to the Massive Scalar-Tensor / Dilaton Klein–Gordon Equation

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step3_yukawa_helmholtz_eom.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md`, `step2_conformal_disformal.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` magyar kézirat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.2. alfejezet, 56–65. sor, Eq. 6–7; DOI: `10.5281/zenodo.23639084`)  
- V23 (`50aa2db0` angol kézirat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` (Sec. 1.2, 56–65. sor, Eq. 6–7)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (200–206. és 245–254. sor; DOI: `10.5281/zenodo.23604423`)  
- V23 SymPy ellenőrző kód: `covariant_closure_50aa2db0.py` (126–149. sor) és `conservation_Residuals.json`

---

## 1. Szó szerinti idézetek a V23 (`50aa2db0`) és V22 kéziratokból

Az 1. lépésben (`step1_dilaton_form.md`) nyitott kérdésként rögzítettük, hogyan válik a nemlokális $\phi_L(\mathbf{x}) = \gamma_L\bar{\rho}_L(\mathbf{x})$ kitevő önálló kinetikus taggal rendelkező, lokális dinamikus skalármezővé. A V23 kézirat **1.2. alfejezete** (*II. megfogalmazás: lokális árnyékolt skalár-tenzor (Helmholtz–Yukawa) hatás*, 56–65. sor, **Eq. 6–7**) ezt szó szerint a következőképpen definiálja:

- **A skalármező Lagrangiánusa és mezőegyenlete (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 57–60. sor, Eq. 6):**
  > *„A II. megfogalmazásban a nemlokális integrálmagot lokális kanonikus árnyékolt skalárhatással helyettesítem:”*
  $$\begin{equation} \mathcal{L}_\phi = -\frac{c^4}{8\pi G \tilde{\alpha}_L}\!\left[\frac{1}{2}g^{\mu\nu}\nabla_\mu\phi_L\nabla_\nu\phi_L + \frac{1}{2L^2}\phi_L^2\right], \quad \left(-\Box_g + \frac{1}{L^2}\right)\phi_L = \frac{4\pi G}{c^2}\alpha_{\text{std}}\,\rho_{\text{src}} \iff (-L^2\nabla^2+1)\phi_L = \alpha_Y\gamma_L\rho_{\text{src}}. \tag{V23 Eq. 6} \end{equation}$$

- **A statikus gömbszimmetrikus Green-függvényes megoldás és a skalár feszültség–energia tenzor (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 61–65. sor, Eq. 7):**
  > *„Statikus gömbi rendszerben a $-\frac{1}{r^2}\frac{d}{dr}(r^2\phi_L') + L^{-2}\phi_L = \frac{4\pi G}{c^2}\alpha_{\text{std}}\rho_{\text{src}}$ ($\alpha_Y = 4\pi\alpha_{\text{std}}$) radiális Helmholtz-operátor egzakt, pozitív-definit Green-függvényes megoldása”*
  $$\begin{equation} \phi_L^{\text{Yuk}}(r) = \alpha_Y\,\frac{G}{c^2}\int_0^{R_*} dr'\,\frac{r' L}{r}\,\sinh\!\frac{r_<}{L}\,e^{-r_>/L}\rho_{\text{src}}(r'), \qquad \phi_L^{\text{Yuk}}(0) = \alpha_Y\,\frac{G}{c^2}\int_0^{R_*} r' e^{-r'/L}\rho_{\text{src}}(r')\,dr', \tag{V23 Eq. 7} \end{equation}$$
  > *„a $\delta S_\phi/\delta g^{\mu\nu}$ variáció pedig a lokális kovariáns $T_{\mu\nu}^{(\phi)} = \frac{c^4}{8\pi G\tilde{\alpha}_L}\bigl[\nabla_\mu\phi_L\nabla_\nu\phi_L - \frac{1}{2}g_{\mu\nu}\bigl((\nabla\phi_L)^2 + L^{-2}\phi_L^2\bigr)\bigr]$ feszültség-energia tenzort szolgáltatja. (`SymPy`: `sym_yukawa_center_limit_residual` $=0$, nincs többlet $L$-faktor). Az egységre normált mag $\alpha_Y = 1$-nek, azaz $\alpha_{\text{std}} = 1/(4\pi)\approx 0{,}0796$-nak felel meg, ami visszaadja a homogén top-hat határesetet ($\phi_L = \gamma_L\rho_{\text{src}}$); az $\alpha_Y = 1/3$ ($\alpha_{\text{std}} = 1/(12\pi)$) a térfogatra illesztett változat.”*

- **A hatással konzisztens (`yukawa_action_el`) Euler–Lagrange-ág (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 65. sor):**
  > *„A `yukawa_action_el` sorok az $\mathcal{L}=-[(\nabla\phi)^2+\phi^2/L^2]/(2\alpha_Y)-f_L\varepsilon$ hatásból következnek ($G=c=1$; `SymPy`: `sym_yukawa_el_source_factor` és `_hydrostatic` maradék $=0$): $B\phi''+(2B/r+B\nu'+B'/2)\phi'-\phi/L^2=-\alpha_Y f_L\varepsilon$, $B=1-2m/r$, $P_{\text{matter}}=f_LP$, $h'=-\nu'+\phi'$, $\rho_\phi=(B\phi'^2+\phi^2/L^2)/(2\alpha_Y)$, külső illesztés $\phi'(R)=-(1/R+1/L)\phi(R)$.”*

---

## 2. Az Euler–Lagrange-mezőegyenlet levezetése és a tömeges dilaton / skalár–tenzor megfeleltetés

### 2.1. Variáció $\phi_L$ szerint (az anyagszektor $f_L(\phi_L) = e^{-\phi_L}$ csatolásával)
Tekintsük a V23 II. megfogalmazásának lokális hatását, ahol az anyagi energiasűrűség az 1. és 2. lépésben levezetett $f_L(\phi_L) = e^{-\phi_L}$ faktorral csatolódik ($\mathcal{L}_m = -f_L(\phi_L)\,\rho c^2 = -e^{-\phi_L}\varepsilon$):
$$S_{\text{II}}[g_{\mu\nu},\phi_L,\Psi] \;=\; \int d^4x\sqrt{-g}\left[\frac{c^4}{16\pi G}R \;-\; \frac{c^4}{8\pi G\tilde{\alpha}_L}\left(\frac{1}{2}g^{\mu\nu}\nabla_\mu\phi_L\nabla_\nu\phi_L + \frac{1}{2L^2}\phi_L^2\right) \;-\; e^{-\phi_L}\,\rho c^2\right].$$
Variáljuk a hatást $\phi_L$ szerint ($\delta S_{\text{II}} / \delta\phi_L = 0$):
$$\frac{1}{\sqrt{-g}}\frac{\delta S_{\text{II}}}{\delta\phi_L} \;=\; \frac{c^4}{8\pi G\tilde{\alpha}_L}\left(\Box_g\phi_L - \frac{1}{L^2}\phi_L\right) \;+\; e^{-\phi_L}\,\rho c^2 \;=\; 0,$$
ahol $\Box_g\phi_L = \frac{1}{\sqrt{-g}}\partial_\mu\bigl(\sqrt{-g}\,g^{\mu\nu}\partial_\nu\phi_L\bigr)$. Átrendezve és megszorozva $-\frac{8\pi G\tilde{\alpha}_L}{c^4}$-nel:
$$\left(-\Box_g + \frac{1}{L^2}\right)\phi_L \;=\; \frac{8\pi G}{c^2}\,\tilde{\alpha}_L\,\rho_{\text{src}}, \qquad \rho_{\text{src}} \;=\; \begin{cases} \rho & \text{(gyenge téri / lineáris forrás, } f_L \simeq 1\text{)}, \\[3pt] f_L(\phi_L)\,\rho = e^{-\phi_L}\rho & \text{(egzakt } \texttt{yukawa\_action\_el} \text{ forrás)}. \end{cases}$$

### 2.2. Az állandók azonosítása ($\tilde{\alpha}_L$, $\alpha_{\text{std}}$, $\alpha_Y$, $\gamma_L$ és $C$)
A V23 Eq. (6)-tal való tagonkénti egybevetésből:
1. **A Poisson- és a Yukawa-normálás kapcsolata:**
   $$\frac{8\pi G}{c^2}\,\tilde{\alpha}_L \;=\; \frac{4\pi G}{c^2}\,\alpha_{\text{std}} \;=\; C\,\alpha_{\text{std}} \quad\implies\quad \tilde{\alpha}_L \;=\; \frac{\alpha_{\text{std}}}{2},$$
   ahol $C = \frac{4\pi G}{c^2}$ a V22 monográfia (`hu_p1.tex`, Eq. 4) alapvető $C$ csatolási állandója.
2. **A $\gamma_L = \frac{G L^2}{c^2}$ faktor kiemelése ($-L^2\Box_g + 1$ alak):**
   Ha a mezőegyenlet mindkét oldalát megszorozzuk $L^2$-tel (kvázisztatikus határesetben $\Box_g \to \nabla^2$):
   $$(-L^2\nabla^2 + 1)\phi_L \;=\; 4\pi\alpha_{\text{std}}\,\left(\frac{G L^2}{c^2}\right)\rho_{\text{src}} \;=\; \alpha_Y\,\gamma_L\,\rho_{\text{src}}, \qquad \boxed{\;\alpha_Y \;\equiv\; 4\pi\,\alpha_{\text{std}} \;=\; 8\pi\,\tilde{\alpha}_L.\;}$$
   Vagy ekvivalens felírásban (ahogy a `README.md` 1.1. pontja írja):
   $$(\nabla^2 - L^{-2})\phi_L \;=\; -\frac{4\pi G}{c^2}\alpha_{\text{std}}\,\rho_{\text{src}} \qquad \left(\text{illetve }\phi_L \to -\phi_L\text{ előjel-konvencióban } +\frac{\alpha_{\text{std}}G}{c^2}\rho_{\text{src}}\text{ a } 4\pi\text{ faktor beolvasztásával}\right).$$

### 2.3. Megfeleltetés a tömeges skalár–tenzor / dilaton Klein–Gordon-egyenletnek
A húrelméleti dilaton vagy általános skalár–tenzor gravitációban (pl. Damour & Polyakov 1994; Fujii & Maeda 2003; Khoury & Weltman 2004) egy $m_\phi$ tömegű, kvadratikus önkölcsönhatási potenciállal ($V(\varphi) = \frac{1}{2}\mu_\phi^2\varphi^2$, ahol $\mu_\phi = \frac{m_\phi c}{\hbar}$) rendelkező kanonikus skalármező mezőegyenlete:
$$\left(-\Box_g + \mu_\phi^2\right)\phi_L \;=\; \frac{4\pi G}{c^2}\,\alpha_{\text{std}}\,\rho_{\text{src}}.$$
Ebből a TAL simítási skálája ($L$) és a skalár-/dilatonmező effektív Compton-hullámhossza, illetve tömege közötti egzakt szótári megfeleltetés:
$$\boxed{\;\mu_\phi \;\equiv\; \frac{1}{L}, \qquad m_\phi \;=\; \frac{\hbar}{L\,c}\;}$$
- **Fizikai jelentés:** Amit a V19/V21 nemlokális megfogalmazása egy $L$ sugarú 3D gömbi top-hat ablakolásként ír le, az a V23 II. megfogalmazásában pontosan egy **$\lambda_{\text{Compton}} = L$ hatótávolságú (azaz $m_\phi = \hbar/(Lc)$ tömegű) Yukawa-árnyékolt skalár-/dilatonmező** cserefolyamatának felel meg.

---

## 3. Hogyan adja vissza a Yukawa-megoldás a 3D top-hat $\phi_L = \gamma_L\bar{\rho}_L$ törvényt?

1. **A 3D Yukawa Green-függvény:**
   A $(-L^2\nabla^2 + 1)G_L(\mathbf{x}-\mathbf{y}) = \delta^{(3)}(\mathbf{x}-\mathbf{y})$ egyenlet alapmegoldása a háromdimenziós euklideszi térben:
   $$G_L(\mathbf{x}-\mathbf{y}) \;=\; \frac{\exp\!\left(-\frac{|\mathbf{x}-\mathbf{y}|}{L}\right)}{4\pi L^2\,|\mathbf{x}-\mathbf{y}|}, \qquad \int_{\mathbb{R}^3} G_L(\mathbf{r})\,d^3r \;=\; \frac{4\pi}{4\pi L^2}\int_0^\infty r\,e^{-r/L}\,dr \;=\; \frac{L^2}{L^2} \;=\; 1.$$
   Így a V23 Eq. (6) egyenlet megoldása tetszőleges $\rho_{\text{src}}(\mathbf{y})$ eloszlásra egy **pozitív definit Yukawa-simítás**:
   $$\phi_L^{\text{Yuk}}(\mathbf{x}) \;=\; \alpha_Y\,\gamma_L\int_{\mathbb{R}^3} G_L(\mathbf{x}-\mathbf{y})\,\rho_{\text{src}}(\mathbf{y})\,d^3y \;=\; \alpha_Y\,\frac{G}{c^2}\int_{\mathbb{R}^3}\frac{e^{-|\mathbf{x}-\mathbf{y}|/L}}{4\pi\,|\mathbf{x}-\mathbf{y}|}\,\rho_{\text{src}}(\mathbf{y})\,d^3y.$$
2. **Homogén / lassan változó sűrűség határesete ($\rho_{\text{src}} = \text{const.}$):**
   Mivel $\int_{\mathbb{R}^3} G_L(\mathbf{r})\,d^3r = 1$, homogén közegben:
   $$\phi_L^{\text{Yuk}} \;=\; \alpha_Y\,\gamma_L\,\rho_{\text{src}}.$$
   Választva az **egységre normált magot** ($\alpha_Y = 1 \iff \alpha_{\text{std}} = \frac{1}{4\pi} \approx 0{,}079577$):
   $$\phi_L^{\text{Yuk}} \;=\; \gamma_L\,\rho_{\text{src}} \quad\implies\quad f_L \;=\; e^{-\phi_L^{\text{Yuk}}} \;=\; \exp(-\gamma_L\rho_{\text{src}}),$$
   ami pontosan visszaadja az 1. lépésben (`step1_dilaton_form.md`) vizsgált TAL exponenciális törvényt!
3. **Ponttömeg / kompakt forrás határesete ($R_* \ll r$):**
   Egy $M = \int \rho_{\text{src}}\,d^3y$ tömegű kompakt forráson kívül ($r > R_*$) a V23 Eq. (7) megoldás klasszikus ötödik erő (*fifth-force*) Yukawa-potenciállá egyszerűsödik:
   $$\phi_L^{\text{Yuk}}(r) \;\propto\; \alpha_Y\,\frac{G M}{c^2 r}\,e^{-r/L}.$$

---

## 4. Ellenőrzések (Checks)

### 4.1. Dimenzióanalízis a V23 Eq. (6) és Eq. (7) egyenletekre
- **Bal oldal a V23 Eq. (6)-ban:** $\phi_L$ dimenziótlan ($[\phi_L] = 1$), a d'Alembert-/Laplace-operátor dimenziója $[\Box_g] = [\nabla^2] = \text{m}^{-2}$, és $[L^{-2}] = \text{m}^{-2}$, tehát:
  $$\left[\left(-\Box_g + \frac{1}{L^2}\right)\phi_L\right] \;=\; \text{m}^{-2}.$$
- **Jobb oldal a V23 Eq. (6)-ban:** $[G/c^2] = \text{m}\,\text{kg}^{-1}$, $[\alpha_{\text{std}}] = 1$, $[\rho_{\text{src}}] = \text{kg}\,\text{m}^{-3}$, tehát:
  $$\left[\frac{4\pi G}{c^2}\alpha_{\text{std}}\,\rho_{\text{src}}\right] \;=\; \left(\text{m}\,\text{kg}^{-1}\right)\cdot \left(\text{kg}\,\text{m}^{-3}\right) \;=\; \text{m}^{-2} \quad \checkmark$$
- **A V23 Eq. (7) középponti integrálja ($\phi_L^{\text{Yuk}}(0)$):**
  $$\left[\alpha_Y\,\frac{G}{c^2}\int_0^{R_*} r'\,e^{-r'/L}\rho_{\text{src}}(r')\,dr'\right] \;=\; 1 \cdot \left(\text{m}\,\text{kg}^{-1}\right)\cdot \text{m}\cdot 1\cdot \left(\text{kg}\,\text{m}^{-3}\right)\cdot \text{m} \;=\; \text{m}^3\,\text{kg}^{-1}\cdot \text{kg}\,\text{m}^{-3} \;=\; 1 \quad (\text{dimenziótlan}) \quad \checkmark$$

### 4.2. SymPy-ellenőrzések egyezése a V23 csomaggal (`conservation_Residuals.json`)
A V23 `covariant_closure_50aa2db0.py` (126–148. sor) három szimbolikus ellenőrzést futtat erre a szektorra, amelyek mindegyike egzakt $0$ maradékot ad:
1. `sym_yukawa_center_limit_residual = 0`: $\lim_{r\to 0} \frac{L}{r}\sinh(r/L) = 1$, így $\phi_L^{\text{Yuk}}(0)$ nem tartalmaz többlet $L$-szorzót $\checkmark$.
2. `sym_yukawa_alphaY_eq_4pi_alphastd_residual = 0`: $\alpha_{\text{std}} = \frac{1}{4\pi} \implies \alpha_Y = 4\pi\alpha_{\text{std}} = 1 \implies \phi_L = \gamma_L\rho$ homogén közegben $\checkmark$.
3. `sym_yukawa_el_source_factor_residual = 0` és `sym_yukawa_el_matter_pressure_residual = 0`: az $\mathcal{L} \supset -e^{-\phi_L}\varepsilon(n)$ anyagtag variációja a $-\alpha_Y f_L\varepsilon$ forrástagot, a $P_{\text{matter}} = f_L P$ effektív nyomást és a $h' = -\nu' + \phi_L'$ módosított hidrosztatikai egyensúlyt adja $\checkmark$.
