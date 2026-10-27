# 1. lépés: A Time Attraction Law $f_L = \exp(-\gamma_L \bar{\rho}_L)$ faktorának és a húrelméleti dilaton-csatolásnak a matematikai megfeleltetése
# Step 1: Mathematical Equivalence of the Form of the TAL Factor $f_L = \exp(-\gamma_L \bar{\rho}_L)$ and the String-Theory Dilaton Coupling

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step1_dilaton_form.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0`): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_EN_2026-10-24.tex` & `tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` & `tal_book_HU_2026-10-24.tex` (DOI: `10.5281/zenodo.23604423`)  
- V19 (7. modell): 3D gömbi top-hat törvény (DOI: `10.5281/zenodo.23497128`)  
- Helyi munkajegyzet: `docs/C_derivation_status.md` (15. sor: `Time_Attraction.pdf (v1–v7)`)

---

## 1. A Time Attraction Law (TAL) szó szerinti definíciói a helyi kéziratokból

A helyi V19–V23 kéziratokban a pozitív definit exponenciális időáramlási / csillapítási faktor $f_L(\mathbf{x})$, a 3D gömbi top-hat átlagsűrűség $\bar{\rho}_L(\mathbf{x})$, valamint a $\gamma_L$ skálafüggő csatolási együttható (és a $C$ csatolási állandó) szó szerint az alábbi alakban szerepel:

### 1.1. Az exponenciális törvény és az állandók definíciója (V22 `hu_p1.tex`, 141–144. és 245–253. sor)
- **Exponenciális időáramlási faktor:**
  $$f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(\mathbf{x})\bigr), \qquad 0 < f_L(\mathbf{x}) \le 1$$
- **3D gömbi top-hat átlagsűrűség ($L$ simítási sugárral, $V_L = \frac{4}{3}\pi L^3$ térfogattal):**
  $$\bar{\rho}_L(\mathbf{x}) \;=\; \frac{3}{4\pi L^3}\int_{|\mathbf{y}-\mathbf{x}|\le L}\rho(\mathbf{y})\,d^3y \;\ge\; 0$$
- **Alapvető gravitációs–időbeli állandók (SI egységekben, CODATA $G = 6{,}67430\times 10^{-11}\;\text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}$, $c = 2{,}99792458\times 10^{8}\;\text{m}\,\text{s}^{-1}$):**
  $$C \;\equiv\; \frac{4\pi G}{c^2} \;=\; 9{,}33306\times 10^{-27}\;\text{m}\,\text{kg}^{-1}, \qquad \frac{G}{c^2} \;=\; \frac{C}{4\pi} \;=\; 7{,}425915\times 10^{-28}\;\text{m}\,\text{kg}^{-1}$$
  $$\gamma_L \;=\; \frac{G\,L^2}{c^2} \;=\; \frac{C\,L^2}{4\pi} \;=\; 7{,}425915\times 10^{-28}\,L^2\;\text{m}^3\,\text{kg}^{-1}$$

### 1.2. A kovariáns segédskalár és az anyag-Lagrangiánus csatolása a V23-ban (`tal_covariant_closure_50aa2db0_EN_2026-10-24.tex`, 32–46. és 54. sor; `C_derivation_status.md`, 15. sor)
- A korai `Time_Attraction.pdf` (v1–v7, idézve a `docs/C_derivation_status.md` 15. sorában) és a V23-as `50aa2db0` kézirat (1. szakasz, Eq. 1–3 és 54. sor) az anyaghatást az $f_L$ exponenciális faktorral súlyozza:
  $$S \;=\; \int d^4x\sqrt{-g}\left[\frac{c^4}{16\pi G}R + \mathcal{L}_X(g_{\mu\nu}, \phi_L, \lambda_L) + \mathcal{L}_{\text{matter}}\right], \qquad \mathcal{L}_{\text{fluid}}(n,\phi_L) \;=\; -f_L(\phi_L)\,\varepsilon(n)$$
  illetve az 1–7. verziós felírásban:
  $$S \;=\; \int d^4x\sqrt{-g}\left[\frac{c^4}{16\pi G}R + f_L\,\mathcal{L}_m\right]$$
- A V23 I. megfogalmazásában (Eq. 2–3) a dimenziótlan $\phi_L(x)$ skalárt a $\lambda_L(x)$ Lagrange-multiplikátor köti a $\Sigma_u$ térbeli hiperfelületen vett 3D gömbi top-hat átlaghoz ($K_L(x,y) = \delta(\tau_u(x)-\tau_u(y))\,\Theta(L - d_\Sigma(x,y))/V_L(x)$):
  $$\frac{\delta S}{\delta \lambda_L(x)} = 0 \quad\implies\quad \phi_L(x) \;=\; \gamma_L\,\bar{\rho}_L(x), \qquad f_L(x) \;=\; e^{-\phi_L(x)} \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(x)\bigr).$$

---

## 2. A húrelmélet alacsony energiás effektív hatása és a dilaton-csatolás konvenciói

A húrelmélet (boszonikus és szuperhúrelmélet) alacsony energiás, 4 dimenzióra redukált effektív térelméletében a gravitációs szektor szükségszerűen tartalmaz egy skalármezőt, a **dilatont ($\phi$)** (lásd pl. Polchinski: *String Theory*, Vol. I, Sec. 3.7; Gasperini & Veneziano, *Phys. Rept.* **373**, 1, 2003; Damour & Polyakov, *Nucl. Phys. B* **423**, 532, 1994).

1. **Húr-keret (*String frame*):**
   Fa-szinten (*tree-level*, a húrcsatolási állandó $g_s = e^{\phi}$ négyzetével, $g_s^{-2} = e^{-2\phi}$ súlyozva) a gravitációs és a zárt/nyílt húr szektor effektív hatása:
   $$S_{\text{string}} \;=\; \int d^4x\sqrt{-g_{(s)}}\,e^{-2\phi}\left[\frac{c^4}{16\pi G_0}\Big(R_{(s)} + 4\,g_{(s)}^{\mu\nu}\partial_\mu\phi\,\partial_\nu\phi\Big) - \dots\right] + \int d^4x\sqrt{-g_{(s)}}\,e^{-\eta\phi}\,\mathcal{L}_m\!\left(g_{(s)\mu\nu},\Psi\right),$$
   ahol:
   - zárt húr / fa-szintű szektorban $\eta = 2$, vagyis az anyag-Lagrangiánust közvetlenül az **$e^{-2\phi}\,\mathcal{L}_m$** exponenciális faktor szorozza;
   - nyílt húr / D-brán (Born–Infeld) szektorban $\eta = 1$, vagyis a csatolási faktor **$e^{-\phi}\,\mathcal{L}_m$**.

2. **Einstein-keret (*Einstein frame*) — kanonikus $R$-tag melletti exponenciális anyagcsatolás:**
   Ahhoz, hogy a görbületi skalár $R$ előtti faktor (miként a TAL V23 Eq. (1)-ben: $\frac{c^4}{16\pi G}R$) kanonikus, $\phi$-től független alakot öltsön, elvégezzük a $g_{(s)\mu\nu} = e^{2\phi} g_{\mu\nu}$ konform átméretezést (Damour & Polyakov 1994; Fujii & Maeda: *The Scalar-Tensor Theory of Gravitation*, Cambridge Univ. Press, 2003). Ekkor a görbületi tag tisztán Einstein–Hilbert-alakú lesz, míg az anyagszektorban egy exponenciális dilaton-szorzó jelenik meg:
   $$S_{\text{Einstein}} \;=\; \int d^4x\sqrt{-g}\left[\frac{c^4}{16\pi G}\Big(R - 2\,g^{\mu\nu}\partial_\mu\phi\,\partial_\nu\phi\Big) \;+\; e^{-2\beta_{\text{dil}}\phi}\,\mathcal{L}_m(\Psi, g_{\mu\nu})\right],$$
   (ahol a konvenciótól és az anyagszektortól függően az exponenciális szorzó szokásos jelölése $e^{-2\phi}$, $e^{-\phi}$ vagy általános skalár–tenzor jelöléssel $e^{-2\beta_{\text{dil}}\phi}$).

---

## 3. A formális matematikai megfeleltetés (Step 1 Mapping)

Összevetve a TAL anyagcsatolási tagját ($f_L(\mathbf{x})\,\mathcal{L}_m$, ahol $f_L(\mathbf{x}) = e^{-\phi_L(\mathbf{x})} = \exp(-\gamma_L\bar{\rho}_L(\mathbf{x}))$) a húrelméleti dilaton exponenciális anyagcsatolásával ($e^{-2\phi}\,\mathcal{L}_m$ a $2\phi$-konvencióban, illetve $e^{-\phi}\,\mathcal{L}_m$ az $1\phi$-konvencióban), a két kifejezés **matematikai függvényalakja pontosan megegyezik**:

$$\boxed{\;f_L(\mathbf{x}) \;=\; \exp\!\bigl(-\gamma_L\,\bar{\rho}_L(\mathbf{x})\bigr) \;=\; e^{-\phi_L(\mathbf{x})} \quad\Longleftrightarrow\quad e^{-2\phi(\mathbf{x})} \;\;\text{(vagy } e^{-\phi(\mathbf{x})}\text{)}\;}$$

azaz a dimenziótlan kitevők közvetlen megfeleltetése:

$$\boxed{\;\phi_L(\mathbf{x}) \;\equiv\; \gamma_L\,\bar{\rho}_L(\mathbf{x}) \;=\; \frac{G L^2}{c^2}\,\bar{\rho}_L(\mathbf{x}) \quad\Longleftrightarrow\quad 2\phi(\mathbf{x}) \;\;\text{(a } e^{-2\phi} \text{ konvencióban), illetve } \phi(\mathbf{x}) \;\;\text{(az } e^{-\phi} \text{ konvencióban).}\;}$$

### Rögzített nyitott kérdés a későbbi lépésekre (nem ebben a lépésben oldjuk meg)
- **A hasonlóság:** Mindkét elméletben egy dimenziótlan skalármennyiség negatív exponenciális függvénye ($e^{-\phi_L}$ illetve $e^{-2\phi}$) szorozza az anyag-Lagrangiánust.
- **A különbség (nyitott kérdés a 2. és további lépésekre):** A Time Attraction Law alapdefiníciójában (és a V23 I. megfogalmazásában) a $\phi_L(\mathbf{x}) = \gamma_L \bar{\rho}_L(\mathbf{x})$ kitevőt algebrailag / nemlokális integrálon keresztül közvetlenül az $L$ skálán átlagolt $\bar{\rho}_L$ tömeg-energia sűrűség határozza meg, míg a húrelméletben a $\phi$ dilaton egy önálló dinamikus mező saját kinetikus taggal ($\propto (\partial\phi)^2$). Ennek a dinamikus kinetikus tagnak és a V23 II. megfogalmazásában szereplő lokális Helmholtz–Yukawa-operátornak $(-L^2\Box_g + 1)\phi_L = \alpha_Y \gamma_L \rho_{\text{src}}$ az összekapcsolását a későbbi lépésekben vizsgáljuk, itt nem nyúlunk hozzá.

---

## 4. Ellenőrzések (Checks)

### 4.1. Dimenzióanalízis ($\gamma_L$ inverz sűrűség dimenziója)
Ellenőrizzük, hogy a $\gamma_L \bar{\rho}_L$ kitevő dimenziótlan-e SI alapegységekben ($\text{m}, \text{kg}, \text{s}$):
- A Newton-féle gravitációs állandó dimenziója:
  $$[G] \;=\; \text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}$$
- A fénysebesség négyzetének dimenziója:
  $$[c^2] \;=\; \text{m}^2\,\text{s}^{-2}$$
- A tömeg–hossz arány ($G/c^2$, illetve $C = 4\pi G/c^2$) dimenziója:
  $$\left[\frac{G}{c^2}\right] \;=\; \frac{\text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}}{\text{m}^2\,\text{s}^{-2}} \;=\; \text{m}\,\text{kg}^{-1}$$
- A simítási skála négyzetének dimenziója:
  $$[L^2] \;=\; \text{m}^2$$
- Következésképpen a $\gamma_L = \frac{G L^2}{c^2}$ csatolási együttható dimenziója:
  $$[\gamma_L] \;=\; \left[\frac{G}{c^2}\right]\cdot [L^2] \;=\; \left(\text{m}\,\text{kg}^{-1}\right)\cdot \text{m}^2 \;=\; \text{m}^3\,\text{kg}^{-1} \;=\; \frac{1}{\text{kg}\,\text{m}^{-3}} \;=\; \frac{1}{[\text{sűrűség}]} \quad \checkmark$$
- Mivel a 3D térfogati tömegsűrűség dimenziója $[\bar{\rho}_L] = \text{kg}\,\text{m}^{-3}$, a szorzat:
  $$[\gamma_L\,\bar{\rho}_L] \;=\; \left(\text{m}^3\,\text{kg}^{-1}\right)\cdot \left(\text{kg}\,\text{m}^{-3}\right) \;=\; 1 \quad (\text{dimenziótlan}) \quad \checkmark$$

*(Megjegyzés a kiírási ellenőrzéshez: a $[G/c^2]$ hányados önmagában $\text{m}\,\text{kg}^{-1}$ [lineáris tömegsűrűség reciproka], és az $L^2$ [$\text{m}^2$] faktorral való szorzás teszi $\gamma_L = G L^2/c^2$ dimenzióját $\text{m}^3\,\text{kg}^{-1}$-vé, ami pontosan a 3D térfogati sűrűség [$\text{kg}\,\text{m}^{-3}$] reciproka, összhangban a V22 `hu_p1.tex` Eq. (3) és Eq. (5) egyenleteivel.)*

### 4.2. Határértékek és előjelek ellenőrzése
1. **Vákuum-határérték ($\bar{\rho}_L \to 0$):**
   Amikor az $L$ sugarú környezetben az átlagolt tömeg-energia sűrűség eltűnik ($\bar{\rho}_L \to 0$), a kitevő $\phi_L = \gamma_L \bar{\rho}_L \to 0$, így:
   $$\lim_{\bar{\rho}_L \to 0} f_L(\bar{\rho}_L) \;=\; \exp(0) \;=\; 1 \quad \checkmark$$
   Ez pontosan visszaadja a módosítatlan standard anyagcsatolást ($\mathcal{L}_{\text{matter}} \to 1\cdot \mathcal{L}_m$), ami a dilaton-képben a $\phi \to 0$ ($e^{-2\phi} \to 1$) aszimptotikus vákuumnak felel meg.
2. **Pozitív sűrűségű tartomány ($\bar{\rho}_L > 0$):**
   Mivel $G > 0$, $c > 0$, $L > 0 \implies \gamma_L > 0$, és a 3D gömbi top-hat átlag nemnegatív energiafeltétel ($\rho(\mathbf{y}) \ge 0$) mellett szigorúan pozitív ($\bar{\rho}_L > 0$), ezért $\gamma_L \bar{\rho}_L > 0$, tehát:
   $$0 \;<\; f_L(\bar{\rho}_L) \;=\; \exp(-\gamma_L \bar{\rho}_L) \;<\; 1 \qquad (\forall\,\bar{\rho}_L > 0) \quad \checkmark$$
   és a függvény szigorúan monoton csökkenő:
   $$\frac{df_L}{d\bar{\rho}_L} \;=\; -\gamma_L\,\exp(-\gamma_L \bar{\rho}_L) \;<\; 0 \quad \checkmark$$
