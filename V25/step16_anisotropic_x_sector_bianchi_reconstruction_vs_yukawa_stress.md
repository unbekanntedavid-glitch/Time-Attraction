# 16. lépés: A V21/V23 algebrai X-szektorbeli anizotróp feszültségrekonstrukció (tömeg- és gradiens-mód, V23 Eq. 4–5) levezetése a Bianchi-azonosságból és összevetése a dinamikus Yukawa skalár-feszültségtenzorral ($T^\mu_{\;\nu}(\phi)$)
# Step 16: Derivation of the V21/V23 Algebraic X-Sector Anisotropic Stress Reconstruction (Mass and Gradient Modes, V23 Eqs. 4–5) from the Bianchi Identity and Comparison with the Dynamical Yukawa Scalar Stress Tensor ($T^\mu_{\;\nu}(\phi)$)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step16_anisotropic_x_sector_bianchi_reconstruction_vs_yukawa_stress.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.1. szakasz, 47–53. sor, Eq. 4–5 és 1(a–b). ábra; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (8. fejezet 1. pont, 400–408. sor; DOI: `10.5281/zenodo.23604423`) és `hu_p3.tex` (8. fejezet 1. pont, 6–12. sor)

---

## 1. Szó szerinti idézetek a V22 és V23 kéziratokból

A 4. és 7. lépésben (`step4_covariant_conservation_hydrostatic.md`, `step7_neutron_star_tov_stability.md`) levezettük a II. megfogalmazás (`yukawa_action_el`) lokális kovariáns feszültség-energia tenzorát ($T_{\mu\nu}^{(\phi)}$), még nem vezettük le azonban lépésről lépésre a V23 kézirat **(4)-es és (5)-ös egyenleteit** (`eq:mass_X` és `eq:grad_X`), a V22 monográfia 8. fejezetében szereplő nemnulla centrális határértéket ($p_{tX}(0) = p_{rX}(0) = -(1-f_c)(\varepsilon_c/3 + P_c)$), valamint a háromszintű kovariáns divergencia-reziduál hierarchiát (`conservation_Residuals.json`):

1. **A V23 (4)–(5) egyenletei és a háromszintű Bianchi-maradék (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 1.1. szakasz, 47–53. sor):**
   > *„Statikus gömbszimmetrikus belső térben, ha a $T^\mu_{\;\nu}(X) = \text{diag}(-\varepsilon_X, P_{r,X}, P_{t,X}, P_{t,X})$ segédszektor veszi fel az ablak-feszültséget a $G^\mu_{\;\nu} = \frac{8\pi G}{c^4}(T^\mu_{\;\nu,\text{matter}} + T^\mu_{\;\nu,X})$ egyenlet jobb oldalán, míg a folyadékelemek kielégítik a $\nabla_\mu T^{\mu\nu}_{\text{matter}} = 0$ egyenletet ($P' = -(\varepsilon+P)\nu'$ geometrizált egységekben, $G=c=1$, $\nu = \Phi_{\text{eff}}/c^2$, $g_{tt}=-e^{2\nu}=-e^{2\Phi_{\text{eff}}/c^2}$), akkor a $\nabla_\mu G^\mu_{\;r} \equiv 0$ Bianchi-azonosság mindkét V21 módban egyértelműen rögzíti $T^\mu_{\;\nu}(X)$-et:”*
   $$\begin{align}
     \textbf{Tömeg-mód } (m' = 4\pi r^2 f_L\varepsilon,\; \nu' = A_{\text{geo}}): &\quad \varepsilon_X = -(1-f_L)\varepsilon,\quad P_{r,X} = 0,\quad P_{t,X} = -\frac{r}{2}(1-f_L)\varepsilon\nu', \tag{V23 Eq. 4}\\
     \textbf{Gradiens-mód } (m' = 4\pi r^2\varepsilon,\; \nu' = f_L A_{\text{geo}}): &\quad \varepsilon_X = 0,\quad P_{r,X} = -(1-f_L)\!\left(P + \frac{m}{4\pi r^3}\right),\; P_{t,X} = P_{r,X} + \frac{r}{2}\bigl(P'_{r,X} + P_{r,X}\nu'\bigr), \tag{V23 Eq. 5}
   \end{align}$$
   > *„ahol $A_{\text{geo}}(r) = (m + 4\pi r^3 P)/[r(r-2m)]$. A `SymPy`-ban elvégzett szimbolikus ellenőrzés (`conservation_Residuals.json`) igazolja, hogy $\nabla_\mu T^\mu_{\;r}(\text{matter}) \equiv 0$ és $\nabla_\mu T^\mu_{\;r}(X) \equiv 0$ azonosságszerűen fennáll. Az 1. táblázatban és az 1(b). ábrán a (4)–(5) összefüggések kiértékelése a numerikus TOV-rácsokon ($L=10\,\text{km}$, $\log_{10}\rho_{0c}=15{,}0$) három különálló szintet ad (`conservation_Residuals.json`): (i) szimbolikus azonosságok, maradék pontosan $0$; (ii) algebrai X-szektor kvadratúra-maradék $\max|\nabla_\mu T^\mu_{\;r}(X)|/|\rho_X\nu'| = 4{,}44\times 10^{-16}$ (mass) és $/|P'_{r,X}| = 4{,}39\times 10^{-14}$ (gradient); (iii) a numerikusan (spline-nal) differenciált háttér divergenciája $\max|\nabla_\mu T^\mu_{\;r}|/|P'| = 2{,}70\times 10^{-5}$ (mass) és $2{,}36\times 10^{-5}$ (gradient), amelyet $h(r)$ véges differenciás deriváltja határoz meg, nem a hatás.”*

2. **A gradiensmód középponti határértéke a V22 monográfiában (`tal_book_HU_2026-10-24.tex`, 8. fejezet 1. pont, 404–408. sor):**
   > *„A két beillesztési mód háttérbeli X-szektora eltérő:*
   > $$\textbf{Gradiensmód:}\quad \rho_X = 0,\quad p_{rX} = -(1-f)\!\left(P + \frac{m}{4\pi r^3}\right),\quad p_{tX} = p_{rX} + \frac{r}{2}\bigl(p'_{rX} + p_{rX}\nu'\bigr),$$
   > *amelynek középponti határértéke $p_{tX}(0) = p_{rX}(0) = -(1-f_c)(\varepsilon_c/3 + P_c)$, általában nem nulla. A **tömegforrásmódban** ezzel szemben $\rho_X = -(1-f)\varepsilon$, $p_{rX} = 0$ és $p_{tX} = -\frac{r}{2}(1-f)\varepsilon\nu'$ ($p_{tX}(0)=0$). A 36 konfigurációs stresszdiagnosztikai táblában szereplő $0{,}140$-es érték az adott konfiguráció háttérbeli könyvelési diagnosztikája, amely nem helyettesíti a kovariáns hatáselvet (`50aa2db0` feladat).”*

---

## 2. A V23 (4) és (5) egyenletek egzakt levezetése az Einstein-egyenletekből és a Bianchi-azonosságból

Sztatikus, gömbszimmetrikus téridőben ($ds^2 = -e^{2\nu(r)}dt^2 + e^{2\lambda(r)}dr^2 + r^2 d\Omega^2$, geometrizált $G = c = 1$ egységekben, ahol $e^{-2\lambda(r)} \equiv 1 - \frac{2m(r)}{r}$) az Einstein-tenzor nemnulla vegyes komponensei:
$$G^t_{\;t} \;=\; -\frac{2m'(r)}{r^2}, \qquad G^r_{\;r} \;=\; -\frac{2m(r)}{r^3} \;+\; \frac{2}{r}\left(1 - \frac{2m(r)}{r}\right)\nu'(r).$$
Írjuk fel az Einstein-egyenleteket a barionos ideális folyadék ($T^\mu_{\;\nu,\text{matter}} = \text{diag}(-\varepsilon, P, P, P)$) és az anizotróp $X$-szektor ($T^\mu_{\;\nu,X} = \text{diag}(-\varepsilon_X, P_{r,X}, P_{t,X}, P_{t,X})$) összegére:
$$G^t_{\;t} \;=\; -8\pi\bigl(\varepsilon + \varepsilon_X\bigr), \qquad G^r_{\;r} \;=\; 8\pi\bigl(P + P_{r,X}\bigr), \qquad G^\theta_{\;\theta} \;=\; G^\phi_{\;\phi} \;=\; 8\pi\bigl(P + P_{t,X}\bigr).$$
Ha a folyadék kielégíti a $P'(r) = -(\varepsilon + P)\nu'(r)$ egyenletet, akkor $\nabla_\mu T^\mu_{\;r}(\text{matter}) = P' + (\varepsilon + P)\nu' \equiv 0$. Mivel a kontrakciós Bianchi-azonosság miatt $\nabla_\mu G^\mu_{\;r} \equiv 0$, az $X$-szektornak önmagában is ki kell elégítenie a radiális kovariáns megmaradási egyenletet:
$$\boxed{\;\nabla_\mu T^\mu_{\;r}(X) \;=\; P'_{r,X}(r) \;+\; \bigl[\varepsilon_X(r) + P_{r,X}(r)\bigr]\nu'(r) \;+\; \frac{2}{r}\bigl[P_{r,X}(r) - P_{t,X}(r)\bigr] \;=\; 0.\;}$$

### 2.1. Levezetés a V21/V23 `mass` (tömegforrás-) módban (V23 Eq. 4)
A `mass` módban a tömegfüggvény és a lapse-gradiens definíciója:
$$m'(r) \;=\; 4\pi r^2 f_L(r)\,\varepsilon(r), \qquad \nu'(r) \;=\; A_{\text{geo}}(r) \;\equiv\; \frac{m(r) + 4\pi r^3 P(r)}{r\bigl(r - 2m(r)\bigr)}.$$
1. **A $(tt)$ Einstein-egyenletből ($\varepsilon_X$):**
   $$-\frac{2m'}{r^2} = -8\pi f_L\varepsilon = -8\pi(\varepsilon + \varepsilon_X) \quad\implies\quad \varepsilon_X(r) \;=\; -\bigl(1 - f_L(r)\bigr)\,\varepsilon(r).$$
2. **Az $(rr)$ Einstein-egyenletből ($P_{r,X}$):**
   Mivel $\nu' = A_{\text{geo}}$, a $G^r_{\;r}$ kifejezése pontosan $8\pi P(r)$-et ad, tehát:
   $$8\pi P(r) = 8\pi\bigl(P(r) + P_{r,X}(r)\bigr) \quad\implies\quad P_{r,X}(r) \;\equiv\; 0.$$
3. **A radiális Bianchi-azonosságból ($P_{t,X}$):**
   Behelyettesítve a $P_{r,X} = 0$ és $P'_{r,X} = 0$ értékeket a $\nabla_\mu T^\mu_{\;r}(X) = 0$ egyenletbe:
   $$\varepsilon_X(r)\,\nu'(r) \;-\; \frac{2}{r}P_{t,X}(r) \;=\; 0 \quad\implies\quad \boxed{\;P_{t,X}(r) \;=\; \frac{r}{2}\,\varepsilon_X(r)\,\nu'(r) \;=\; -\frac{r}{2}\bigl(1 - f_L(r)\bigr)\,\varepsilon(r)\,\nu'(r),\;}$$
   amely pontosan a **V23 (4) egyenlet** (`eq:mass_X`), és mivel $\nu'(r) = \mathcal{O}(r)$ az origóban, $P_{t,X}(r) = \mathcal{O}(r^2) \xrightarrow{r\to 0} 0$.

### 2.2. Levezetés a V21/V23 `gradient` (gradiens-) módban (V23 Eq. 5) és a centrális izotrópia bizonyítása
A `gradient` módban a tömegfüggvény és a lapse-gradiens definíciója:
$$m'(r) \;=\; 4\pi r^2\,\varepsilon(r), \qquad \nu'(r) \;=\; f_L(r)\,A_{\text{geo}}(r) \;=\; f_L(r)\,\frac{m(r) + 4\pi r^3 P(r)}{r\bigl(r - 2m(r)\bigr)}.$$
1. **A $(tt)$ Einstein-egyenletből ($\varepsilon_X$):**
   $$-\frac{2m'}{r^2} = -8\pi\varepsilon = -8\pi(\varepsilon + \varepsilon_X) \quad\implies\quad \varepsilon_X(r) \;\equiv\; 0.$$
2. **Az $(rr)$ Einstein-egyenletből ($P_{r,X}$):**
   Behelyettesítve a $\nu'(r) = f_L(r)\frac{m + 4\pi r^3 P}{r(r-2m)}$ összefüggést a $G^r_{\;r} = 8\pi(P + P_{r,X})$ egyenletbe:
   $$P(r) + P_{r,X}(r) \;=\; \frac{G^r_{\;r}}{8\pi} \;=\; -\frac{m(r)}{4\pi r^3} \;+\; f_L(r)\left(\frac{m(r)}{4\pi r^3} + P(r)\right) \quad\implies\quad \boxed{\;P_{r,X}(r) \;=\; -\bigl(1 - f_L(r)\bigr)\left(P(r) + \frac{m(r)}{4\pi r^3}\right).\;}$$
3. **A radiális Bianchi-azonosságból ($P_{t,X}$):**
   Mivel $\varepsilon_X = 0$, a $\nabla_\mu T^\mu_{\;r}(X) = P'_{r,X} + P_{r,X}\nu' + \frac{2}{r}(P_{r,X} - P_{t,X}) = 0$ egyenletből kifejezve $P_{t,X}(r)$-et:
   $$\boxed{\;P_{t,X}(r) \;=\; P_{r,X}(r) \;+\; \frac{r}{2}\Bigl[P'_{r,X}(r) \;+\; P_{r,X}(r)\,\nu'(r)\Bigr],\;}$$
   ami pontosan a **V23 (5) egyenlet** (`eq:grad_X`).
4. **Középponti határérték és izotrópia ($r \to 0$):**
   A csillag középpontja körül ($r \to 0$) $m(r) = \frac{4\pi}{3}\varepsilon_c r^3 + \mathcal{O}(r^5)$, így $\lim_{r\to 0}\frac{m(r)}{4\pi r^3} = \frac{\varepsilon_c}{3}$. Mivel az origóban minden reguláris skalárprofil deriváltja nulla ($P'(0) = 0$, $f_L'(0) = 0$, $\nu'(0) = 0$), ezért $P'_{r,X}(r) = \mathcal{O}(r)$, tehát az $\frac{r}{2}[P'_{r,X} + P_{r,X}\nu'] = \mathcal{O}(r^2)$ anizotrópia-tag eltűnik az origóban, és:
   $$\boxed{\;\lim_{r\to 0} P_{t,X}(r) \;=\; \lim_{r\to 0} P_{r,X}(r) \;=\; -\bigl(1 - f_c\bigr)\left(\frac{\varepsilon_c}{3} + P_c\right) \;\neq\; 0, \qquad \lim_{r\to 0}\bigl[P_{t,X}(r) - P_{r,X}(r)\bigr] \;=\; 0.\;}$$
   Ezért nem lép fel geometriai kúpszingularitás az origóban: bár a centrális nyomás nem nulla, az $X$-szektor feszültsége az origóban **izotróp** ($P_{r,X}(0) = P_{t,X}(0)$)!

---

## 3. Összevetés a II. megfogalmazás (`yukawa_action_el`) dinamikus skalár-feszültségtenzorával

A V23 II. megfogalmazásában (`yukawa_action_el`, `step4` és `step7`) nem algebrai rekonstrukcióval határozzuk meg az $X$-szektort, hanem közvetlenül a kovariáns hatásból ($S = \int d^4x\sqrt{-g}\bigl[\frac{R}{16\pi} - \frac{(\nabla\phi_L)^2 + \phi_L^2/L^2}{2\alpha_Y} - f_L(\phi_L)\varepsilon(n)\bigr]$). Ekkor a standard GR-beli $(\varepsilon, P)$ forráshoz képest megjelenő effektív $X$-szektor két dinamikai részből áll össze:
$$T^\mu_{\;\nu,X,\text{dyn}} \;=\; \bigl(T^\mu_{\;\nu,\text{matter}} - T^\mu_{\;\nu,\text{GR}}\bigr) \;+\; T^\mu_{\;\nu}(\phi),$$
azaz komponensenként ($B(r) = 1 - 2m(r)/r$):
$$\begin{aligned}
  \varepsilon_{X,\text{dyn}}(r) &\;=\; -\bigl(1 - f_L(r)\bigr)\,\varepsilon(r) \;+\; \frac{1}{2\alpha_Y}\left(B(r)\bigl[\phi_L'(r)\bigr]^2 + \frac{\phi_L(r)^2}{L^2}\right), \\
  P_{r,X,\text{dyn}}(r) &\;=\; -\bigl(1 - f_L(r)\bigr)\,P(r) \;+\; \frac{1}{2\alpha_Y}\left(B(r)\bigl[\phi_L'(r)\bigr]^2 - \frac{\phi_L(r)^2}{L^2}\right), \\
  P_{t,X,\text{dyn}}(r) &\;=\; -\bigl(1 - f_L(r)\bigr)\,P(r) \;-\; \frac{1}{2\alpha_Y}\left(B(r)\bigl[\phi_L'(r)\bigr]^2 + \frac{\phi_L(r)^2}{L^2}\right).
\end{aligned}$$
Ebből három mély szerkezeti kapcsolat látható az I. és II. megfogalmazás között:
1. **A tömeg-módbeli $-(1-f_L)\varepsilon$ tag eredete:** A dinamikus $\varepsilon_{X,\text{dyn}}(r)$ első tagja egzaktul a V21/V23 `mass` mód $\varepsilon_X = -(1-f_L)\varepsilon$ felhígulási tagja, amelyet a skalármező pozitív $\rho_\phi = \frac{1}{2\alpha_Y}(B[\phi_L']^2 + \phi_L^2/L^2)$ saját-energiasűrűsége egészít ki.
2. **Az anizotrópia ($\Pi_\phi \equiv P_{r,X,\text{dyn}} - P_{t,X,\text{dyn}}$) tisztán gradiens-eredete:**
   $$\boxed{\;\Pi_\phi(r) \;\equiv\; P_{r,X,\text{dyn}}(r) - P_{t,X,\text{dyn}}(r) \;=\; p_{r,\phi}(r) - p_{t,\phi}(r) \;=\; \frac{B(r)}{\alpha_Y}\bigl[\phi_L'(r)\bigr]^2 \;\ge\; 0.\;}$$
   Mivel $\phi_L'(r) = \mathcal{O}(r)$ az origó közelében, a dinamikus anizotrópia is $\Pi_\phi(r) = \mathcal{O}(r^2) \xrightarrow{r\to 0} 0$, tehát az origóban a nyomás egzaktul izotróp: $P_{r,X,\text{dyn}}(0) = P_{t,X,\text{dyn}}(0) = -(1-f_c)P_c - \frac{\phi_L(0)^2}{2\alpha_Y L^2} < 0$.
3. **Miért váltja ki a $\phi_L$ Euler–Lagrange-egyenlete a külön Bianchi-rekonstrukciót?**
   Míg az algebrai V21 módokban $P_{t,X}(r)$-et utólag, a $\nabla_\mu T^\mu_{\;r}(X) = 0$ egyenletből kellett kiszámítani, addig a `yukawa_action_el` rendszerben a 4. lépésben (`step4_covariant_conservation_hydrostatic.md`) bizonyított $\nabla_\mu\bigl(T^\mu_{\;r,\text{matter}} + T^\mu_{\;r}(\phi)\bigr) \equiv 0$ Noether-azonosság a $\phi_L$ mezőegyenlet következtében automatikusan teljesül!

---

## 4. Ellenőrzések (Checks)

1. **Szimbolikus Bianchi-azonosság ellenőrzése mindkét V23 algebrai módra:**
   - **Tömeg-mód:** $P'_{r,X} + (\varepsilon_X + P_{r,X})\nu' + \frac{2}{r}(P_{r,X} - P_{t,X}) = 0 + \varepsilon_X\nu' + \frac{2}{r}\bigl(0 - \frac{r}{2}\varepsilon_X\nu'\bigr) = \varepsilon_X\nu' - \varepsilon_X\nu' \equiv 0$ $\checkmark$.
   - **Gradiens-mód:** $P'_{r,X} + (0 + P_{r,X})\nu' + \frac{2}{r}\Bigl(P_{r,X} - \bigl[P_{r,X} + \frac{r}{2}(P'_{r,X} + P_{r,X}\nu')\bigr]\Bigr) = P'_{r,X} + P_{r,X}\nu' - (P'_{r,X} + P_{r,X}\nu') \equiv 0$ $\checkmark$.
2. **$f_L \to 1$ ($\phi_L \to 0$) GR-határérték:**
   - Mind a `mass` módban ($\varepsilon_X, P_{r,X}, P_{t,X} \to 0$), mind a `gradient` módban ($\varepsilon_X, P_{r,X}, P_{t,X} \to 0$), mind a dinamikus `yukawa_action_el` esetben ($\varepsilon_{X,\text{dyn}}, P_{r,X,\text{dyn}}, P_{t,X,\text{dyn}} \to 0$) az $X$-szektor azonosan eltűnik, és visszakapjuk az izotróp GR TOV-egyenleteket $\checkmark$.
