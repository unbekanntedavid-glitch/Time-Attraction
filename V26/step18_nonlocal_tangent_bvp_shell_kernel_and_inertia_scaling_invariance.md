# 18. lépés: A V23 (9) nemlokális statikus érintőoperátor (`eq:coupled_bvp`), a 3D gömbhéj-metszeti $K_{ij}$ kvadratúra-mag és a $W \to \kappa W$ tehetetlenség-skálázási invariancia egzakt levezetése
# Step 18: Exact Derivation of the V23 Eq. (9) Nonlocal Static Tangent Operator (`eq:coupled_bvp`), the 3D Spherical-Shell Intersection Kernel $K_{ij}$, and the $W \to \kappa W$ Inertia-Scaling Invariance

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step18_nonlocal_tangent_bvp_shell_kernel_and_inertia_scaling_invariance.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (2. szakasz, 127–153. sor, Eq. 9 és 2. táblázat; DOI: `10.5281/zenodo.23639084`)  
- V23 ellenőrző kód és eredmények: `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/covariant_closure_50aa2db0.py` (179–191. és 583–695. sor) és `coupled_stability_results.csv`  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/tal_book_HU_2026-10-24.tex` (7. fejezet 1. pont, 374. sor és 8. fejezet 2. pont, 409–419. sor; DOI: `10.5281/zenodo.23604423`)

---

## 1. Szó szerinti idézetek a V22 és V23 kéziratokból

A 17. lépésben (`step17_chandrasekhar_factored_difference_eq8_and_metric_weight_indefiniteness.md`) levezettük a V23 (8) egyenletet (`eq:factored_diff`) és a `gradient` módbeli `metric` operátor $\mathcal{W}_{\text{metric}}(r) \le 0$ felszíni előjelváltását. A V23 kézirat 2. szakaszában közvetlenül ezután következik a **V23 (9) egyenlet (`eq:coupled_bvp`)**, a szakaszonkénti politróp $r_{\text{trans}}$ átmeneti sugarakhoz illesztett nemlokális $K_{ij}$ mag-perturbáció, a $W \to 4W$ tehetetlenség-skálázási invariancia, az $\omega^2 \simeq -F(0)/F'(0)$ elsőrendű érintőbecslés és az $\mathbf{M}_{\text{sym}} = \mathbf{D}^T\mathrm{diag}(w\,\Delta r)\mathbf{D}$ szimmetrizált Galerkin-mátrix:

1. **A V23 (9) csatolt statikus nemlokális érintőmodell és kvázi-dinamikus BVP-teszt (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 2. szakasz, 127–136. sor, Eq. 9):**
   > *„A csatolt rendszert (a V21 `dense_profile` háttér körül, nem az új ágak körül) a $\delta h(r) = \delta P/(\varepsilon+P)$ entalpiaperturbációra és a $\delta m(r)$ tömegperturbációra linearizálom a szakaszonkénti politróp egzakt $r_{\text{trans}}$ átmeneti sugaraihoz illesztett rácson, beépítve a nemlokális mag-perturbációs operátort:”*
   $$\begin{equation}
     \delta f_L(r_i) = -\gamma_L f_L(r_i)\sum_j K_{ij}\,\frac{\varepsilon(r_j)+P(r_j)}{c_s^2(r_j)}\,\delta h(r_j), \qquad
     \begin{cases}
       \dfrac{d(\delta h)}{dr} = -\delta\nu'\!\left[\delta h, \delta m, \delta f_L[\delta h]\right] + \omega_0^2\,e^{\lambda-2\nu}\xi, \\[4pt]
       \dfrac{d(\delta m)}{dr} = 4\pi r^2\!\left[f_{\text{mass}}\frac{\varepsilon+P}{c_s^2}\delta h + \varepsilon\,\delta f_{\text{mass}}[\delta h]\right].
     \end{cases} \tag{V23 Eq. 9}
   \end{equation}$$
   > *„$\omega_0^2 = 0$-nál a (9) rendszer pontosan a statikus egyensúlyi egyenletek $\log_{10}\rho_{0c}$ szerinti érintőoperátora. Mint a 2. táblázat és a 3. ábra mutatja, a (9) csatolt sajátérték-feladat megoldása az SLy $L=10\,\text{km}$ sorozatok mentén $\sim 7\times 10^{-5}\,\text{dex}$-re ad zérushelyet a $dM/d\rho_c=0$ ponttól (a tiszta GR diszkretizációs alapérték $7{,}12\times 10^{-5}\,\text{dex}$). Mivel $\omega_0^2=0$-nál a feladat $\mathcal{H}_{\text{lin}}[\delta h]=0$-ra redukálódik, ez az egyezés a statikus nemlokális linearizációt igazolja, nem a dinamikai tehetetlenséget: a $W\to 4W$ skálázás $\log_{10}\rho_{0c}=15{,}30427$-nél $\omega^2_{4W}/\omega^2_W = 0{,}2500$ arányt ad ($\omega^2_W=4{,}897\times 10^{6}\,\text{s}^{-2}$), azonos $F(0)=0{,}023252$ lövőfüggvény mellett, így a zérushely tehetetlenség-invariáns. A nyers BVP $W$ mátrix nem szimmetrikus; a szimmetrizált Galerkin-alak $\mathbf{M}_{\text{sym}} = \mathbf{D}^{T}\mathrm{diag}(w\,\Delta r)\mathbf{D}$ szimmetrikus, $\min\mathrm{diag}\,\mathbf{M}_{\text{sym}} = 5{,}68\times 10^{-9}>0$; ez Gram-mátrix-tulajdonság, nem a fizikai kinetikus Hesse-mátrix pozitív definitségének bizonyítéka. Az `omega0_sq_sym_tangent_s2` oszlop az elsőrendű érintő-becslést ($\omega^2\simeq -F(0)/F'(0)$) adja.”*

2. **A 3D gömbhéj-metszeti integrál (`tal_book_HU_2026-10-24.tex`, 7. fejezet 1. pont, 374. sor):**
   > *„Egy $R_*$ sugarú gömbszimmetrikus csillag belsejében az $r$ radiális koordinátánál vett egzakt 3D top-hat $\bar{\rho}_L(r)$ átlag analitikusan számítható a $B(\mathbf{x},L)$ gömb és az $r' \in [|r-L|, r+L]$ gömbhéjak metszetén, amiből $f_L(r) = \exp(-\gamma_L\bar{\rho}_L(r))$.”*

---

## 2. A 3D gömbhéj-metszeti $K_{ij}$ kvadratúra-mag és az $r_{\text{trans}}$ átmeneti rács egzakt levezetése

### 2.1. A 3D gömbi top-hat gömbhéj-metszeti térszög-aránya (`build_tophat_matrix`)
Legyen $\mathbf{x}_i = (0,0,r_i)$ a lekérdezési pont ($r_i > 0$), és tekintsünk egy $r_j > 0$ sugarú, $\Delta r_j$ vastagságú koncentrikus gömbhéjat a csillag középpontja körül! A gömbhéj egy $\mathbf{y} = (r_j\sin\theta\cos\varphi,\,r_j\sin\theta\sin\varphi,\,r_j\cos\theta)$ pontjának euklideszi távolságnégyzete $\mathbf{x}_i$-től a koszinusztétel szerint:
$$|\mathbf{y} - \mathbf{x}_i|^2 \;=\; r_i^2 \;+\; r_j^2 \;-\; 2\,r_i\,r_j\cos\theta.$$
A $|\mathbf{y} - \mathbf{x}_i| \le L$ feltétel ekvivalens azzal, hogy:
$$\cos\theta \;\ge\; \mu_{\min}(r_i, r_j; L) \;\equiv\; \frac{r_i^2 + r_j^2 - L^2}{2\,r_i\,r_j}.$$
Az $r_j$ sugarú gömbhéjnak a $B(\mathbf{x}_i, L)$ gömbbe eső felszínaránya ($\mathrm{frac}(r_i, r_j; L) = \frac{1}{4\pi}\int_{|\mathbf{y}-\mathbf{x}_i|\le L}\sin\theta\,d\theta\,d\varphi = \frac{1}{2}\int_{\max(-1, \mu_{\min})}^{1} d(\cos\theta)$) három geometriai tartományra bomlik:
1. **Egészében bezárt gömbhéj ($r_i + r_j \le L$):** ekkor $r_i^2 + r_j^2 - L^2 \le -2 r_i r_j \implies \mu_{\min} \le -1$, így az egész gömbhéj $B(\mathbf{x}_i, L)$ belsejébe esik: $\mathrm{frac}(r_i, r_j; L) = 1$.
2. **Részlegesen metszett gömbhéj ($|r_i - r_j| < L < r_i + r_j$):** ekkor $-1 < \mu_{\min} < +1$, így a gömbsüveg térszög-aránya:
   $$\mathrm{frac}(r_i, r_j; L) \;=\; \frac{1 - \mu_{\min}}{2} \;=\; \frac{2 r_i r_j - (r_i^2 + r_j^2 - L^2)}{4\,r_i\,r_j} \;=\; \frac{L^2 - (r_i - r_j)^2}{4\,r_i\,r_j}.$$
3. **Kívül eső gömbhéj ($|r_i - r_j| \ge L$):** ekkor $\mu_{\min} \ge +1$, így a metszet üres: $\mathrm{frac}(r_i, r_j; L) = 0$.

Mivel a $j$-edik gömbhéj térfogata $4\pi r_j^2\,\Delta r_j$, a $B(\mathbf{x}_i, L)$ gömb térfogata pedig $V_L = \frac{4}{3}\pi L^3$, a diszkrét 3D top-hat magmátrix eleme pontosan (`covariant_closure_50aa2db0.py`, 188–190. sor):
$$\boxed{\;K_{ij} \;=\; \frac{4\pi r_j^2\,\mathrm{frac}(r_i, r_j; L)\,\Delta r_j}{\frac{4}{3}\pi L^3} \;=\; \frac{3\,r_j^2\,\Delta r_j}{L^3} \begin{cases}
  1, & r_i + r_j \le L, \\[4pt]
  \dfrac{L^2 - (r_i - r_j)^2}{4\,r_i\,r_j}, & |r_i - r_j| < L < r_i + r_j, \\[6pt]
  0, & |r_i - r_j| \ge L.
\end{cases}\;}$$
Következésképpen a $\bar{\rho}_L(r_i) = \sum_j K_{ij}\,\varepsilon(r_j)$ nemlokális átlag és az $f_L(r_i) = \exp(-\gamma_L\bar{\rho}_L(r_i))$ faktor első variációja a $\delta\varepsilon(r_j) = \frac{d\varepsilon}{dh}\big|_{r_j}\delta h(r_j) = \frac{\varepsilon(r_j)+P(r_j)}{c_s^2(r_j)}\delta h(r_j)$ entalpiaperturbáció szerint pontosan a V23 (9) egyenlet első fele:
$$\delta f_L(r_i) \;=\; -\gamma_L\,f_L(r_i)\sum_{j=1}^{N_{\text{tot}}} K_{ij}\,\frac{\varepsilon(r_j) + P(r_j)}{c_s^2(r_j)}\,\delta h(r_j).$$

### 2.2. Miért szükséges a cellahatárok $r_{\text{trans}}$ átmeneti sugarakhoz való illesztése?
A szakaszonkénti politróp (piecewise-polytrope) állapotegyenletekben (SLy, APR4, H4) a nyomás $P(h)$ és az energiasűrűség $\varepsilon(h)$ folytonos a $h_{\text{tr}} \in \text{EOS}.h_b$ politróp-határokon, de az adiabatikus index $\Gamma_k$ és így a hangsebesség-négyzet $c_s^2(h) = \frac{dP}{d\varepsilon}$ **szakadásos ugrást** szenved a $h(r_{\text{trans}}) = h_{\text{tr}}$ átmeneti sugaraknál! Ha egy kvadratúra-cella átlépne egy $r_{\text{trans}}$ sugáron, akkor a $\frac{\varepsilon+P}{c_s^2}\delta h$ integráns ugrása $\mathcal{O}(\Delta r)$ elsőrendű kvadratúra-hibát okozna, elrontva a $dM/d\rho_c = 0$ fordulópont-illeszkedést. Ezért a `solve_coupled_radial_omega2` (`covariant_closure_50aa2db0.py`, 604–616. sor) Brent-módszerrel (`brentq(lambda x: p["h"](x) - h_tr, 0.0, R)`) egzaktul megkeresi az összes $0 < r_{\text{trans}} < R_*$ átmeneti sugarat, és ezeket rögzített cellahatárként építi be a rácsba.

---

## 3. A V23 (9) csatolt BVP mátrixalakja, a $W \to \kappa W$ skálázási invariancia és a Galerkin-szimmetrizálás

### 3.1. Az operátormátrixok ($\mathbf{S}$, $\mathbf{D}_M$, $\mathbf{D}_A$, $\mathbf{W}$) és a Fredholm-féle lövőfüggvény $F(\omega^2)$
Geometriai egységekben ($G = c = 1$, ahol $\gamma_L = L^2$) a V23 (9) rendszer két elsőrendű egyenlete a cellaközéppontok $\delta\mathbf{h} = (\delta h(r_1),\dots,\delta h(r_{N_{\text{tot}}}))^T$ vektorára (`covariant_closure_50aa2db0.py`, 629–660. sor):
1. **A linearizált tömegforrás-mátrix ($\mathbf{S}$) és az integrált tömegperturbáció ($\delta\mathbf{m} = \mathbf{D}_M\,\delta\mathbf{h}$):**
   $$\frac{d(\delta m)}{dr} \;=\; (\mathbf{S}\,\delta\mathbf{h})_i, \qquad \mathbf{S} \;=\; \begin{cases}
     \mathrm{diag}\!\bigl(4\pi r^2 f_L\bigr)\,\bigl[\mathbf{I} - L^2\,\mathrm{diag}(\varepsilon)\,\mathbf{K}\bigr]\,\mathrm{diag}\!\left(\dfrac{\varepsilon+P}{c_s^2}\right), & \text{(\texttt{mass} mód)}, \\[8pt]
     \mathrm{diag}\!\left(4\pi r^2 \dfrac{\varepsilon+P}{c_s^2}\right), & \text{(\texttt{gradient} mód)}.
   \end{cases}$$
   A trapézszabályú alsó háromszög-integrálmátrixszal $\bigl(\mathbf{M}_{\text{int}}\bigr)_{ij} = \Delta r_j\,\Theta(i - j) + \frac{1}{2}\Delta r_i\,\delta_{ij}$ a tömegperturbáció-operátor $\mathbf{D}_M \equiv \mathbf{M}_{\text{int}}\mathbf{S}$, vagyis $\delta\mathbf{m} = \mathbf{D}_M\,\delta\mathbf{h}$ (automatikusan kielégítve a $\delta m(0) = 0$ centrális peremfeltételt).
2. **A linearizált hidrosztatikai gradiens-mátrix ($\mathbf{D}_A$):**
   Mivel $A(r) = \frac{m + 4\pi r^3 P}{r(r-2m)}$, a parciális deriváltak $c_{dm}(r) \equiv \frac{\partial A}{\partial m} = \frac{1 + 2 r A}{r(r-2m)}$ és $c_{dh}(r) \equiv (\varepsilon+P)\frac{\partial A}{\partial P} = \frac{4\pi r^3(\varepsilon+P)}{r(r-2m)}$, így az alap-perturbáció $\mathbf{D}_{A,\text{base}} = \mathrm{diag}(c_{dm})\,\mathbf{D}_M + \mathrm{diag}(c_{dh})$, és:
   $$\delta\boldsymbol{\nu}' \;=\; \mathbf{D}_A\,\delta\mathbf{h}, \qquad \mathbf{D}_A \;=\; \begin{cases}
     \mathbf{D}_{A,\text{base}}, & \text{(\texttt{mass} mód)}, \\[6pt]
     \mathrm{diag}(f_L)\,\mathbf{D}_{A,\text{base}} \;-\; L^2\,\mathrm{diag}(f_L A)\,\mathbf{K}\,\mathrm{diag}\!\left(\dfrac{\varepsilon+P}{c_s^2}\right), & \text{(\texttt{gradient} mód)}.
   \end{cases}$$
3. **Az adiabatikus tömegmegmaradási kapcsolat és a tehetetlenségi mátrix ($\mathbf{W}$):**
   Adiabatikus radiális elmozdulásnál a Lagrange-féle tömegmegmaradás $\Delta m = 0 \implies \delta m(r) = -4\pi r^2 f_{\text{mass}}(\varepsilon + P)\,r\xi(r)$ ($f_{\text{mass}} = f_L$ a `mass` módban és $f_{\text{mass}} = 1$ a `gradient` módban), amiből az elmozdulás kifejezhető a tömegperturbációval:
   $$-\omega^2 e^{2(\lambda-\nu)} r\xi(r) \;=\; +\omega^2\,w_{\text{inert}}(r)\,\delta m(r), \qquad w_{\text{inert}}(r) \;\equiv\; \frac{e^{2(\lambda-\nu)}}{4\pi r^2 f_{\text{mass}}(r)\bigl(\varepsilon(r) + P(r)\bigr)} \;>\; 0.$$
   Ezért a tehetetlenségi mátrix $\mathbf{W}(\kappa) \equiv \kappa\,\mathrm{diag}(w_{\text{inert}})\,\mathbf{D}_M$ (ahol $\kappa \equiv W_{\text{scale}} > 0$ a skálázási faktor), és a V23 (9) entalpiaegyenlet $\delta h(0) = 1$ centrális normálással integrálva a következő diszkrét Fredholm-rendszert adja (`covariant_closure_50aa2db0.py`, 664–665. sor):
   $$\boxed{\;\Bigl[\mathbf{I} \;+\; \mathbf{M}_{\text{int}}\bigl(\mathbf{D}_A \;+\; \omega^2\,\mathbf{W}(\kappa)\bigr)\Bigr]\,\delta\mathbf{h}(\omega^2; \kappa) \;=\; \mathbf{1}, \qquad F(\omega^2; \kappa) \;\equiv\; \delta m(R_*; \omega^2, \kappa) \;=\; \boldsymbol{\Delta r}^T \mathbf{S}\,\delta\mathbf{h}(\omega^2; \kappa) \;=\; 0.\;}$$
   A felszíni $F(\omega^2; \kappa) = \delta m(R_*; \omega^2, \kappa) = 0$ Birkhoff-peremfeltétel a $\delta m(R_*) = -4\pi R_*^3 f_{\text{mass}}(\varepsilon + P)\xi(R_*)$ és $\Delta P(R_*) = 0$ szabad felszíni peremfeltétellel egyenértékű.

### 3.2. A $W \to \kappa W$ tehetetlenség-skálázási invariancia és az elsőrendű érintőtétel bizonyítása
A fenti Fredholm-egyenletből három egzakt algebrai következmény adódik:
1. **Egzakt skálázási azonosság ($\forall\,\kappa > 0$):**
   Mivel $\mathbf{W}(\kappa) = \kappa\,\mathbf{W}(1)$, a $\mathbf{D}_A + \omega^2\mathbf{W}(\kappa) = \mathbf{D}_A + (\kappa\omega^2)\mathbf{W}(1)$ kombináció kizárólag a $\tilde{\omega}^2 \equiv \kappa\,\omega^2$ szorzaton keresztül függ az $\omega^2$ sajátértéktől és a $\kappa$ tehetetlenségi skálától:
   $$\delta\mathbf{h}(\omega^2; \kappa) \;\equiv\; \delta\mathbf{h}(\kappa\omega^2; 1) \quad\Longrightarrow\quad F(\omega^2; \kappa) \;\equiv\; F(\kappa\omega^2; 1).$$
   Következésképpen, ha $\omega^2_W$ a sajátérték $\kappa = 1$ mellett ($F(\omega^2_W; 1) = 0$), akkor $\kappa$-szoros tehetetlenségi súly mellett a pontos sajátérték:
   $$\boxed{\;\omega^2_{\kappa W} \;=\; \frac{\omega^2_W}{\kappa}, \qquad \frac{\omega^2_{4W}}{\omega^2_W} \;=\; \frac{1}{4} \;=\; 0{,}2500\quad\text{egzaktul!}\;}$$
   Speciálisan $\log_{10}\rho_{0c} = 15{,}30427$-nél (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 135. sor) $F(0; \kappa) \equiv F(0; 1) = 0{,}023252$ ($\kappa$-tól függetlenül), $\omega^2_W = 4{,}897\times 10^6\;\text{s}^{-2}$, és $\omega^2_{4W}/\omega^2_W = 0{,}2500$.
2. **A statikus fordulópont ($\omega_0^2 = 0$) tehetetlenség-invarianciája:**
   $\omega^2 = 0$-nál a $\omega^2\mathbf{W}(\kappa)$ tag kiesik, így $(\mathbf{I} + \mathbf{M}_{\text{int}}\mathbf{D}_A)\delta\mathbf{h}(0) = \mathbf{1}$. Mivel $\mathbf{D}_A$ és $\mathbf{S}$ a nemlokális statikus TOV-egyenletek egzakt Frechet-deriváltja a $\delta h(0)$ centrális entalpiaperturbáció szerint, a $\delta\mathbf{h}(0)$ megoldás arányos a $\frac{\partial h(r; \rho_{0c})}{\partial\log_{10}\rho_{0c}}$ statikus sorozat-deriválttal, a lövőfüggvény pedig:
   $$F(0) \;=\; \delta m(R_*; 0) \;=\; \left(\frac{\partial h_c}{\partial\log_{10}\rho_{0c}}\right)^{-1}\frac{dM}{d\log_{10}\rho_{0c}}.$$
   Ezért $F(0) = 0 \iff \frac{dM}{d\rho_c} = 0$, ami megmagyarázza, miért csökken a fordulópont-eltérés $\Delta\log_{10}\rho_{0c} = +6{,}97\times 10^{-5}\;\text{dex}$-re (`mass`, $19{,}2\times$/$58{,}3\times$ javulás) és $+7{,}39\times 10^{-5}\;\text{dex}$-re (`gradient`, $48{,}3\times$/$59{,}3\times$ javulás), megegyezve a tiszta GR diszkretizációs alapvonallal ($+7{,}12\times 10^{-5}\;\text{dex}$, V23 2. táblázat).
3. **Az elsőrendű érintő-becslés (`omega0_sq_sym_tangent_s2`) és a Galerkin-féle Gram-mátrix ($\mathbf{M}_{\text{sym}}$):**
   Az $F(\omega^2) = 0$ gyök $\omega^2 = 0$ körüli elsőrendű Taylor-sorfejtéséből $F(0) + F'(0)\,\omega^2 + \mathcal{O}(\omega^4) = 0$, vagyis (`covariant_closure_50aa2db0.py`, 680–690. sor):
   $$\omega^2_{\text{tangent}} \;\simeq\; -\frac{F(0)}{F'(0)}.$$
   Noha a nyers BVP-beli $\mathbf{W} = \mathrm{diag}(w_{\text{inert}})\mathbf{D}_M$ mátrix nem szimmetrikus (mivel $\mathbf{D}_M = \mathbf{M}_{\text{int}}\mathbf{S}$ alsó háromszög-integrátorral szorzott), a $\delta\mathbf{m} = \mathbf{D}_M\delta\mathbf{h}$ tömegperturbációkkal felírt kvadratikus kinetikus alak szimmetrizált Galerkin-mátrixa (`covariant_closure_50aa2db0.py`, 685–686. sor):
   $$\boxed{\;\mathbf{M}_{\text{sym}} \;\equiv\; \mathbf{D}_M^T\,\mathrm{diag}\!\bigl(w_{\text{inert}}\,\Delta r\bigr)\,\mathbf{D}_M \;=\; \mathbf{M}_{\text{sym}}^T.\;}$$
   Mivel $w_{\text{inert}}(r_i)\,\Delta r_i > 0$ minden $i = 1,\dots,N_{\text{tot}}$ cellában (mind `mass`, mind `gradient` módban, mivel a nevezőben az anyagi $\varepsilon + P > 0$ áll!), bármely $\mathbf{v} \neq \mathbf{0}$ vektorra $\mathbf{v}^T\mathbf{M}_{\text{sym}}\mathbf{v} = \sum_{i=1}^{N_{\text{tot}}} w_{\text{inert}}(r_i)\Delta r_i\,(\mathbf{D}_M\mathbf{v})_i^2 \ge 0$, és főátlójának minimuma $\min\mathrm{diag}\,\mathbf{M}_{\text{sym}} = 5{,}68\times 10^{-9} > 0$. Ahogy a V23 kézirat (135. és 152. sor) hangsúlyozza, ez konstrukció szerinti **Gram-mátrix-tulajdonság** ($\mathbf{B}^T\mathbf{B}$ alak, ahol $\mathbf{B} = \mathrm{diag}(\sqrt{w_{\text{inert}}\Delta r})\mathbf{D}_M$), míg a valódi dinamikai kinetikus Hesse-mátrix pozitív definitségét a 12. lépésben (`step12_second_variation_coupled_pulsation_hessian.md`) a hatás $\delta^2 S$ második variációjából vezettük le.

---

## 4. Ellenőrzések (Checks)

1. **A 3D gömbhéj-metszeti mag $K_{ij}$ normálása ($\sum_j K_{ij} = 1$):**
   - Ha a csillag sugara $R_* > r_i + L$ és a sűrűség homogén ($\varepsilon(r_j) \equiv 1$), akkor folytonos határesetben:
     $$\sum_j K_{ij} \;\longrightarrow\; \frac{3}{L^3}\left[\int_0^{L - r_i} r'^2\,dr' \;+\; \int_{|L - r_i|}^{L + r_i} r'^2\,\frac{L^2 - (r_i - r')^2}{4\,r_i\,r'}\,dr'\right] \;=\; 1 \qquad \forall\,r_i \ge 0 \;\checkmark.$$
2. **Az $\omega^2_{\text{tangent}} \simeq -F(0)/F'(0)$ és $\omega^2_{\text{BVP}}$ egyezése a fordulópontnál (`coupled_stability_results.csv`):**
   - A `coupled_stability_results.csv` 23. sorában (SLy GR, $\log_{10}\rho_{0c} = 15{,}302445$) $\omega^2_{\text{BVP}} = 1{,}04011\times 10^5\;\text{s}^{-2}$ és $\omega^2_{\text{tangent}} = 1{,}03992\times 10^5\;\text{s}^{-2}$ ($0{,}018\%$ eltérés); az 56. sorban (SLy `mass`, $L = 10\;\text{km}$, $\log_{10}\rho_{0c} = 15{,}308269$) pedig $\omega^2_{\text{BVP}} = 9{,}84022\times 10^4\;\text{s}^{-2}$ és $\omega^2_{\text{tangent}} = 9{,}83845\times 10^4\;\text{s}^{-2}$ ($0{,}018\%$ eltérés), igazolva az elsőrendű érintő-linearizáció egzaktságát az $\omega_0^2 = 0$ fordulópont környezetében $\checkmark$.

---

## 5. A V21 nemlokális 3D top-hat érintőmodell (`step18`) és a V23 lokális görbült `yukawa_action_el` csatolt stabilitás (`step12`, `step22`) explícit szétválasztása

A modellágak egyértelmű megkülönböztetése érdekében rögzítjük a két stabilitási számítás közötti egzakt hatókört:

| Tulajdonság | **Branch I-A**: V21 nemlokális 3D top-hat (`step18`, `covariant_closure_50aa2db0.py`) | **Branch II-Fluid**: V23 görbült külső `yukawa_action_el` (`step12`, `step22`, `curved_yukawa_adm_and_stability.py`) |
|---|---|---|
| **Egyensúlyi háttér** | V21 `tophat_mass` és `tophat_gradient` (`tp.dense_profile`) | Csatolt belső + görbült külső ($r \in [0, R_*+12L]$) `yukawa_action_el` |
| **Perturbációs változók** | $(\delta h, \delta m)$ a nemlokális $\delta f_L = -\gamma_L f_L K_{ij}\delta\varepsilon_j$ héjmaggal | Csatolt $(\delta h, \delta\phi_L, \delta m)$ a belső térben és a görbült külső vákuumfarokban |
| **Tömeg-fordulópont (SLy, $L=10\;\text{km}$)** | $\log_{10}\rho_{0c}(M_{\max}) = 15{,}29924$ (`mass`), $15{,}30106$ (`gradient`) | $\log_{10}\rho_{0c}(M_{\text{ADM},\max}) = 15{,}316404$ ($\alpha_Y=1$), $15{,}307166$ ($\alpha_Y=1/3$) |
| **$\omega_0^2 = 0$ zérushely** | $15{,}29931$ (`mass`), $15{,}30113$ (`gradient`) | $15{,}316410$ ($\alpha_Y=1$), $15{,}307170$ ($\alpha_Y=1/3$) |
| **Fordulópont-eltérés $\Delta\log_{10}\rho_{0c}$** | $+6{,}97\times 10^{-5}\;\text{dex}$ (`mass`), $+7{,}39\times 10^{-5}\;\text{dex}$ (`gradient`) | $+5{,}80\times 10^{-6}\;\text{dex}$ ($\alpha_Y=1$, $2097\times$ csökkenés), $+4{,}20\times 10^{-6}\;\text{dex}$ ($\alpha_Y=1/3$, $995\times$ csökkenés; $L=1\;\text{km}$: $-6{,}86\times 10^{-6}\;\text{dex}$, $165\times$) |

> **V25 frissítés (görbült `yukawa_action_el` ág):** 4-komponensű külső BVP ($\delta\nu_R=1{,}63\times 10^{-14}$), véges differenciákkal igazolt külső ($\partial\phi_L'(R_*)/\partial m(R_*)=-1{,}68365\times 10^{-4}\;\text{km}^{-2}$) és belső ($r=6\;\text{km}$ $\partial\phi_L''/\partial m=-1{,}72833\times 10^{-4}\;\text{km}^{-3}$, `scalar_m_sign = -1.0` fő operátor) deriváltak, dinamikus skalárhullám, GL2 gyökkonvergencia ($p_{\text{eff,root}}=4{,}248$). SLy $L=10$ km, $\alpha_Y=1$: $M_{\text{ADM}}$-fordulópont $15{,}316404$, fő dinamikus $\omega_0^2=0$ gyök $15{,}316410$ ($+5{,}80\times 10^{-6}$ dex). Részletes táblák: `step22`, 4. és 7. szakasz.
