# 10. lépés: A nemlokális bilaterális mag-hatás (I. megfogalmazás), a lapse-súlyozott reciprocitási potenciál ($V_{\text{eff},N}$), valamint az I. és II. (Yukawa skalár–tenzor) megfogalmazás közötti momentum-/EFT-híd
# Step 10: The Nonlocal Bilateral Kernel Action (Formulation I), Lapse-Weighted Reciprocity Potential ($V_{\text{eff},N}$), and the Moment/EFT Bridge Between Formulations I and II (Yukawa Scalar-Tensor)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step10_bilateral_kernel_reciprocity_eft_synthesis.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step9_vhe_photon_propagation_hess_dispersion.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.1. és 1.2. szakasz, 37–66. sor, Eq. 2–3, 1. táblázat; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p3.tex` (11. és 12. fejezet, 94–120. sor; DOI: `10.5281/zenodo.23604423`)

---

## 1. Szó szerinti idézetek a V23 (`50aa2db0`) kéziratból és a V22 monográfiából

Az 1–9. lépésben végigvezettük a II. megfogalmazást (a lokális árnyékolt Helmholtz–Yukawa skalár–tenzor hatást) és annak megfigyelési és laboratóriumi következményeit. Ebben a lépésben levezetjük az **I. megfogalmazást (nemlokális bilaterális mag-hatás $\phi_L, \lambda_L$ segédmezőkkel, súlyozatlan `bilateral_eom` és lapse-súlyozott `bilateral_lapse_weighted` reciprok egyensúly)**, valamint az **I. (3D gömbi top-hat, Branch I-A/I-B) és a II. (Yukawa skalár–tenzor, Branch II-Fluid, $\alpha_Y = 1$ vs. $\alpha_Y = 1/3$) megfogalmazást összekötő egzakt momentum-kifejtést**, szigorúan rögzítve, hogy a két mag másodrendű gradiens-momentuma ($\frac{1}{10}L^2\nabla^2$ vs. $L^2\nabla^2$) eltér egymástól.

1. **A V23 I. megfogalmazása: Nemlokális pozitív-definit mag-hatás és Lagrange-multiplikátor (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 1.1. szakasz, 38–46. sor, Eq. 2–3):**
   > *„Az I. megfogalmazásban a $\phi_L(x)$ segédskalárt egy $\lambda_L(x)$ Lagrange-multiplikátor mező kapcsolja a folyadék $u^\mu$ négyes-sebességére ortogonális $\Sigma_u$ hiperfelületen ($h_{\mu\nu} = g_{\mu\nu} + u_\mu u_\nu/c^2$) vett térbeli top-hat átlaghoz:”*
   $$\begin{equation}
     S_{\text{I},X}[\phi_L, \lambda_L, g_{\mu\nu}] = -\frac{c^4}{8\pi G \gamma_L}\int d^4x\sqrt{-g}\,\lambda_L(x)\!\left[\phi_L(x) - \gamma_L \int_{\Sigma_u} d^4y\sqrt{-g(y)}\,K_L(x,y)\,\rho_{\text{src}}(y)\right], \tag{V23 Eq. 2}
   \end{equation}$$
   > *„ahol $K_L(x,y) = \delta(\tau_u(x)-\tau_u(y))\,\Theta(L - d_\Sigma(x,y))/V_L(x)$ az $L$ sugarú, $V_L = \frac{4}{3}\pi L^3$ saját-térfogatú pozitív-definit 3D gömbi mag. A $\lambda_L(x)$, $\phi_L(x)$ és $\rho_{\text{src}}(x)$ szerinti variáció az alábbi Euler–Lagrange-lezárást adja:”*
   $$\begin{equation}
     \frac{\delta S}{\delta \lambda_L(x)} = 0 \implies \phi_L(x) = \gamma_L\bar{\rho}_L(x),\quad f_L(x) = e^{-\phi_L(x)}; \qquad
     \frac{\delta S}{\delta \phi_L(x)} = 0 \implies \lambda_L(x) = \frac{8\pi G\gamma_L}{c^4}\frac{\partial \mathcal{L}_{\text{eff}}}{\partial \phi_L}(x). \tag{V23 Eq. 3}
   \end{equation}$$

2. **A bilaterális mag-reciprocitás ($x \leftrightarrow y$), a lapse-súlyozott $V_{\text{eff},N}$ és a kétcellás variációs próba (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 54. sor):**
   > *„Ha $\lambda_L$-t a folyadékba olvasztjuk $\mathcal{L}_{\text{fluid}}(n,\phi_L) = -f_L(\phi_L)\,\varepsilon(n)$ alakban, és $\varepsilon$-t explicit módon és a magon belül is variáljuk, akkor $V_{\text{eff}}(x) = f_L(x) - \gamma_L\!\int K_L(y,x) f_L(y)\varepsilon(y)\,d^3y$ és a módosított egyensúly $dh/dr = -A(r) - dV_{\text{recip}}/dr$, $V_{\text{recip}} = \ln V_{\text{eff}}$. Ezt önkonzisztens fixpontként oldom meg (`bilateral_eom` sorok az 1. táblázatban); eltér az algebrailag lezárt V21 mass-módtól. Az $N=e^\nu$ lapse-szal súlyozott $\int dV\,N f_L\varepsilon$ variációja (rögzített mag és metrika mellett) $V_{\text{eff},N}(x)=f_L-\frac{\gamma_L}{N(x)}\!\int dV_y\,N(y)K(y,x)f_L\varepsilon$, amelyet a lapse-egyenletben $P_{\text{matter}}=f_LP$ mellett oldok meg (`bilateral_lapse_weighted` sorok; kétcellás próba $N=(1;0{,}8)$, $\varepsilon=(1;2)$, $\gamma=0{,}1$: $(0{,}74882;0{,}72084)$, súlyozatlanul $0{,}73160$). A kovariáns lezáráshoz még a mag és a térfogatelem metrikus variációja szükséges; egyik sor sem a (2) hatás EL-megoldása.”*

3. **Az $\alpha_Y = 1$ (egységre normált) és $\alpha_Y = 1/3$ (térfogatra illesztett) Yukawa-normálás (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 65. sor):**
   > *„Az egységre normált mag $\alpha_Y = 1$-nek, azaz $\alpha_{\text{std}} = 1/(4\pi)\approx 0{,}0796$-nak felel meg, ami visszaadja a homogén top-hat határesetet ($\phi_L = \gamma_L\rho_{\text{src}}$); az $\alpha_Y = 1/3$ ($\alpha_{\text{std}} = 1/(12\pi)$) a térfogatra illesztett változat.”*

---

## 2. Matematikai levezetés

### 2.1. A bilaterális mag-reciprocitás ($x \leftrightarrow y$) és a lapse-súlyozott $V_{\text{eff},N}(x)$ levezetése
Tekintsük a sztatikus energiafunkcionált, amelyben az anyag-Lagrangiánus $\mathcal{L}_{\text{fluid}}(n, \phi_L) = -f_L(\phi_L)\,\varepsilon(n)$ ($f_L(x) = \exp[-\gamma_L\int dV_y\,K_L(x,y)\varepsilon(y)]$) az $N(x) = e^{\nu(x)}$ lapse-faktorral súlyozva szerepel (rögzített háttérmetrika és rögzített $K_L(x,y)$ mag mellett):
$$E_{\text{matter}}[\varepsilon] \;=\; \int dV_x\,N(x)\,f_L(x)\,\varepsilon(x), \qquad \phi_L(x) \;=\; \gamma_L\int dV_y\,K_L(x,y)\,\varepsilon(y).$$
Variáljuk $E_{\text{matter}}[\varepsilon]$-t az $\varepsilon(x)$ energiasűrűség szerint! Mivel $\varepsilon$ két helyen jelenik meg — egyszer explicit szorzóként az $x$ pontban, egyszer pedig az $f_L(y) = \exp[-\gamma_L\int dV_x\,K_L(y,x)\varepsilon(x)]$ kitevőjében minden olyan $y$ pontban, amelynek simítógömbje tartalmazza $x$-et — a funkcionális derivált a láncszabály és az $x \leftrightarrow y$ néma integrációs változók felcserélése (*bilaterális mag-reciprocitás*) alapján:
$$\delta E_{\text{matter}} \;=\; \int dV_x\,N(x)\,f_L(x)\,\delta\varepsilon(x) \;+\; \int dV_y\,N(y)\,\varepsilon(y)\,\frac{\partial f_L(y)}{\partial\phi_L(y)}\,\gamma_L\!\int dV_x\,K_L(y,x)\,\delta\varepsilon(x).$$
Mivel $\frac{\partial f_L(y)}{\partial\phi_L(y)} = -f_L(y)$, a második integrálban felcserélve a $dV_y$ és $dV_x$ integrálási sorrendet, majd kiemelve $N(x)\,\delta\varepsilon(x)$-et:
$$\delta E_{\text{matter}} \;=\; \int dV_x\,N(x)\,\underbrace{\left[f_L(x) \;-\; \frac{\gamma_L}{N(x)}\int dV_y\,N(y)\,K_L(y,x)\,f_L(y)\,\varepsilon(y)\right]}_{\equiv\; V_{\text{eff},N}(x)}\delta\varepsilon(x).$$
Ezzel egzaktul levezettük a V23 1.1. szakaszában (54. sor) szereplő két bilaterális potenciált:
1. **Lapse-súlyozott bilaterális potenciál (`bilateral_lapse_weighted`):**
   $$\boxed{\;V_{\text{eff},N}(x) \;=\; f_L(x) \;-\; \frac{\gamma_L}{N(x)}\int dV_y\,N(y)\,K_L(y,x)\,f_L(y)\,\varepsilon(y).\;}$$
2. **Súlyozatlan bilaterális potenciál (`bilateral_eom`, $N(y)/N(x) \to 1$):**
   $$\boxed{\;V_{\text{eff}}(x) \;=\; f_L(x) \;-\; \gamma_L\int dV_y\,K_L(y,x)\,f_L(y)\,\varepsilon(y).\;}$$
A hidrosztatikai Euler-egyenletben ($\delta\varepsilon/\delta n = (\varepsilon + P)/n = \mu e^h$) a kémiai potenciál állandóságának feltétele ($N(r)\,V_{\text{eff}}(r)\,e^{h(r)} = \text{konst.}$) logaritmikus deriválás után közvetlenül adja a módosított egyensúlyi egyenletet:
$$\frac{dh}{dr} \;=\; -\frac{d\nu}{dr} \;-\; \frac{d V_{\text{recip}}}{dr}, \qquad V_{\text{recip}}(r) \;\equiv\; \ln V_{\text{eff}}(r).$$
Mivel $\gamma_L\int K_L f_L\varepsilon\,dV_y > 0$, a bilaterális reciprok tag **csökkenti** $V_{\text{eff}}$-et a csillag sűrű magjában ($V_{\text{eff}}(0) < f_L(0)$), így befelé mutató többlet-összehúzó erőt ad, ami megmagyarázza, miért ad a `bilateral_eom` és a `bilateral_lapse_weighted` ág a GR-nél kisebb maximális tömeget ($M_{\max}^{\text{SLy}} = 1{,}93214\,M_\odot$ és $1{,}90217\,M_\odot$ $L = 10\text{ km}$-nél, szemben a V21 `mass` mód $2{,}09822\,M_\odot$ értékével; V23 1. táblázat).

### 2.2. Egzakt momentum-kifejtési összevetés: miben egyezik meg és miben különbözik a 3D top-hat (Branch I) és a Yukawa (Branch II) megfogalmazás?
Fejtsük Taylor-sorba a $\rho_{\text{src}}(\mathbf{x} + \mathbf{u})$ forrássűrűséget az $\mathbf{x}$ pont körül a 3D gömbi top-hat átlagban ($|\mathbf{u}| \le L$):
$$\bar{\rho}_L(\mathbf{x}) \;=\; \frac{3}{4\pi L^3}\int_{|\mathbf{u}|\le L}\!\left[\rho_{\text{src}}(\mathbf{x}) + u_i\partial_i\rho_{\text{src}}(\mathbf{x}) + \frac{1}{2}u_i u_j\partial_i\partial_j\rho_{\text{src}}(\mathbf{x}) + \dots\right]d^3u.$$
Gömbiszimmetria miatt $\int u_i\,d^3u = 0$ és $\int u_i u_j\,d^3u = \frac{1}{3}\delta_{ij}\int_0^L 4\pi u^4\,du = \frac{4\pi}{15}L^5\,\delta_{ij}$, így a 3D top-hat operátor gradiens-kifejtése:
$$\boxed{\;\phi_L^{\text{top-hat}}(\mathbf{x}) \;=\; \gamma_L\,\bar{\rho}_L(\mathbf{x}) \;=\; \gamma_L\left(1 \;+\; \frac{L^2}{10}\nabla^2 \;+\; \mathcal{O}(L^4\nabla^4)\right)\rho_{\text{src}}(\mathbf{x}).\;}$$
Másrészt a II. megfogalmazás Helmholtz–Yukawa egyenletének ($(-L^2\nabla^2 + 1)\phi_L^{\text{Yuk}} = \alpha_Y\gamma_L\rho_{\text{src}}$) formális inverze kis gradienseknél:
$$\boxed{\;\phi_L^{\text{Yuk}}(\mathbf{x}) \;=\; \alpha_Y\,\gamma_L\,\bigl(1 - L^2\nabla^2\bigr)^{-1}\rho_{\text{src}}(\mathbf{x}) \;=\; \alpha_Y\,\gamma_L\left(1 \;+\; L^2\nabla^2 \;+\; \mathcal{O}(L^4\nabla^4)\right)\rho_{\text{src}}(\mathbf{x}).\;}$$

**Kritikus elméleti megkülönböztetés (Point 5):**
1. **Zeroth-moment (homogén $\nabla^2\rho_{\text{src}} \to 0$) határérték:**
   Az $\alpha_Y = 1$ ($\alpha_{\text{std}} = \frac{1}{4\pi} \approx 0{,}07958$) normálás mellett a 3D top-hat és a Yukawa Green-függvény **kizárólag a nulladik momentumban** ($\int K_L(\mathbf{u})\,d^3u = \int G_L(\mathbf{u})\,d^3u = 1$, azaz homogén közegben: $\phi_L^{\text{Yuk}} = \gamma_L\rho_{\text{src}} = \phi_L^{\text{top-hat}}$) egyezik meg!
2. **Eltérés a másodrendű gradiens-momentumban ($\frac{1}{10}L^2\nabla^2$ vs. $L^2\nabla^2$) és a kompakt tartóban:**
   Inhomogén közegben ($\nabla^2\rho_{\text{src}} \neq 0$) a két modell **nem** azonos effektív térelmélet:
   - rögzített $L$ és $\alpha_Y = 1$ mellett a másodrendű $\nabla^2\rho_{\text{src}}$ együttható 10-szeres faktorral különbözik ($\frac{1}{10}L^2$ a 3D top-hat magra, szemben az $1\cdot L^2$ értékkel a Yukawa Green-függvényre; a másodrendű momentum egyezéséhez $\ell_Y = L/\sqrt{10}$-re lenne szükség);
   - a 3D top-hat mag kompakt tartójú ($r \le L$, így a csillagon kívüli halo pontosan $r = R_* + L$-nél véget ér), míg a Yukawa Green-függvény végtelen tartományon exponenciálisan lecsengő ($\propto e^{-r/L}/r$) görbült külső skalárfarkat és $\Delta M_{\phi,\text{ext}}^{\text{curved}}$ ADM-tömegjárulékot hoz létre;
   - az $\alpha_Y = \frac{1}{3}$ ($\alpha_{\text{std}} = \frac{1}{12\pi} \approx 0{,}02653$) változat a $\frac{4}{3}\pi L^3$ gömbtérfogat és a $4\pi L^3$ Yukawa-normálási térfogat arányát ($\frac{4\pi/3}{4\pi} = \frac{1}{3}$) illeszti a csillagcentrumban ($r \to 0$).

---

## 3. Ellenőrzések (Checks)

1. **A V23-beli kétcellás variációs próba (`bilateral_lapse_weighted` vs. `bilateral_eom`) egzakt reprodukálása:**
   - Legyen két egyenlő térfogatú cella ($\Delta V_1 = \Delta V_2 = 1/2$, szimmetrikus átfedő mag $K_{ij} = 1$, így $\bar{\varepsilon} = \frac{1}{2}(\varepsilon_1 + \varepsilon_2)$), ahol a V23 54. sora szerint $N = (1{,}0;\; 0{,}8)$, $\varepsilon = (1{,}0;\; 2{,}0)$ és $\gamma_L = 0{,}1$:
     - Átlagsűrűség: $\bar{\varepsilon} = \frac{1 + 2}{2} = 1{,}5 \implies \phi_L = \gamma_L\bar{\varepsilon} = 0{,}15 \implies f_L(1) = f_L(2) = e^{-0{,}15} = 0{,}86070798$.
     - **Súlyozatlan (`bilateral_eom`):**
       $$V_{\text{eff}}(1) = V_{\text{eff}}(2) = f_L\bigl(1 - \gamma_L\bar{\varepsilon}\bigr) = e^{-0{,}15}\times (1 - 0{,}15) = 0{,}85\,e^{-0{,}15} = \mathbf{0{,}73160178} \approx \mathbf{0{,}73160} \quad \checkmark$$
     - **Lapse-súlyozott (`bilateral_lapse_weighted`):**
       A lapse-súlyozott forrásátlag $\overline{N\varepsilon} = \frac{1}{2}(N_1\varepsilon_1 + N_2\varepsilon_2) = \frac{1}{2}(1{,}0\times 1{,}0 + 0{,}8\times 2{,}0) = \frac{2{,}6}{2} = 1{,}3$:
       $$V_{\text{eff},N}(1) = f_L\left(1 - \frac{\gamma_L\,\overline{N\varepsilon}}{N_1}\right) = e^{-0{,}15}\left(1 - \frac{0{,}13}{1{,}0}\right) = 0{,}87\,e^{-0{,}15} = \mathbf{0{,}74881594} \approx \mathbf{0{,}74882} \quad \checkmark$$
       $$V_{\text{eff},N}(2) = f_L\left(1 - \frac{\gamma_L\,\overline{N\varepsilon}}{N_2}\right) = e^{-0{,}15}\left(1 - \frac{0{,}13}{0{,}8}\right) = 0{,}8375\,e^{-0{,}15} = \mathbf{0{,}72084293} \approx \mathbf{0{,}72084} \quad \checkmark$$
   - Mindhárom szám jegyről jegyre megegyezik a V23 kézirat (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 54. sor) szövegével.
