# 13. lépés: A görbült téridőbeli bilaterális top-hat ablak és térfogatelem metrikus variációja ($\delta[\sqrt{h}\,K_L(x,y)]/\delta g^{\mu\nu}$), a Synge-féle világfüggvény és a felszíni $\delta_D(L - d_\Sigma)$ határfeszültség
# Step 13: Metric Variation of the Curved-Spacetime Bilateral Top-Hat Window and Volume Element ($\delta[\sqrt{h}\,K_L(x,y)]/\delta g^{\mu\nu}$), Synge's World Function, and the $\delta_D(L - d_\Sigma)$ Boundary Surface Stress

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step13_curved_window_metric_variation_synge_world_function.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step12_second_variation_coupled_pulsation_hessian.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (Absztrakt, 28. sor és 1.1. szakasz, 38–55. sor, Eq. 2–3; DOI: `10.5281/zenodo.23639084`)  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p2.tex` (4. és 7. fejezet) és `hu_p3.tex` (8. és 11. fejezet; DOI: `10.5281/zenodo.23604423`)

---

## 1. Szó szerinti idézetek a V23 (`50aa2db0`) kéziratból

A 10. lépésben (`step10_bilateral_kernel_reciprocity_eft_synthesis.md`) levezettük a V23 I. megfogalmazásának súlyozatlan (`bilateral_eom`) és lapse-súlyozott (`bilateral_lapse_weighted`) reciprok egyensúlyát **rögzített mag és metrika mellett**. A V23 kézirat absztraktja (28. sor) és 1.1. szakasza (38–54. sor) azonban kifejezetten rögzíti, hogy az I. megfogalmazás (Eq. 2) kovariáns Euler–Lagrange-lezárásához még egyetlen nyitott matematikai lépés tartozik: **a görbült téri $K_L(x,y)$ ablakfüggvény és a térfogatelem $g_{\mu\nu}$ metrika szerinti variációja**:

1. **A V23 I. megfogalmazásának kovariáns hatása és a $\Sigma_u$ hiperfelületi mag (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 1.1. szakasz, 38–46. sor, Eq. 2–3):**
   > *„Az I. megfogalmazásban a $\phi_L(x)$ segédskalárt egy $\lambda_L(x)$ Lagrange-multiplikátor mező kapcsolja a folyadék $u^\mu$ négyes-sebességére ortogonális $\Sigma_u$ hiperfelületen ($h_{\mu\nu} = g_{\mu\nu} + u_\mu u_\nu/c^2$) vett térbeli top-hat átlaghoz:”*
   $$\begin{equation}
     S_{\text{I},X}[\phi_L, \lambda_L, g_{\mu\nu}] = -\frac{c^4}{8\pi G \gamma_L}\int d^4x\sqrt{-g}\,\lambda_L(x)\!\left[\phi_L(x) - \gamma_L \int_{\Sigma_u} d^4y\sqrt{-g(y)}\,K_L(x,y)\,\rho_{\text{src}}(y)\right], \tag{V23 Eq. 2}
   \end{equation}$$
   > *„ahol $K_L(x,y) = \delta(\tau_u(x)-\tau_u(y))\,\Theta(L - d_\Sigma(x,y))/V_L(x)$ az $L$ sugarú, $V_L = \frac{4}{3}\pi L^3$ saját-térfogatú pozitív-definit 3D gömbi mag.”*

2. **A V23 kéziratban nyitottként megjelölt görbült téri ablak- és térfogatelem-variáció (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 28. és 54. sor):**
   > - **Absztrakt (28. sor):** *„[E]z nem dinamikai bizonyíték a **görbült téri ablak-variációra**, amely nyitott kérdés marad.”*
   > - **1.1. szakasz (54. sor):** *„Az $N=e^\nu$ lapse-szal súlyozott $\int dV\,N f_L\varepsilon$ variációja (rögzített mag és metrika mellett) $V_{\text{eff},N}(x)=f_L-\frac{\gamma_L}{N(x)}\!\int dV_y\,N(y)K(y,x)f_L\varepsilon$, amelyet a lapse-egyenletben $P_{\text{matter}}=f_LP$ mellett oldok meg (`bilateral_lapse_weighted` sorok; kétcellás próba $N=(1;0{,}8)$, $\varepsilon=(1;2)$, $\gamma=0{,}1$: $(0{,}74882;0{,}72084)$, súlyozatlanul $0{,}73160$). **A kovariáns lezáráshoz még a mag és a térfogatelem metrikus variációja szükséges; egyik sor sem a (2) hatás EL-megoldása.**”*

---

## 2. Matematikai levezetés: A görbült téri térfogatelem és a geodetikus $d_\Sigma(x,y)$ ablak metrikus variációja

### 2.1. A $\Sigma_u$ hiperfelületi integrálmérték és a saját-gömbtérfogat $V_L(x)$ variációja
Sztatikus téridőben ($ds^2 = -N(x)^2 c^2 dt^2 + h_{ij}(x)dx^i dx^j$, ahol $N(x) = e^{\nu(x)}$ a lapse-függvény, $h_{ij}(x)$ pedig a $\Sigma_u$ térbeli hiperfelület indukált 3-metrikája, $\sqrt{-g(x)} = N(x)\sqrt{h(x)}$) a $\tau_u(x) = \tau_u(y)$ egyidejű hiperfelületre leszűkített mag-integrál:
$$\bar{\rho}_L(x) \;=\; \int_{\Sigma_u} d^3y\,\sqrt{h(y)}\,\mathcal{K}_L(x,y)\,\varepsilon(y), \qquad \mathcal{K}_L(x,y) \;\equiv\; \frac{\Theta\!\bigl(L - d_\Sigma(x,y)\bigr)}{V_L(x)},$$
ahol $d_\Sigma(x,y)$ az $x$ és $y$ pontokat $\Sigma_u$-n összekötő minimális térbeli geodetikus ívhossz, a $B_\Sigma(x,L) = \{y \in \Sigma_u : d_\Sigma(x,y) \le L\}$ geodetikus gömb saját-térfogata pedig:
$$V_L(x) \;=\; \int_{\Sigma_u} d^3y\,\sqrt{h(y)}\,\Theta\!\bigl(L - d_\Sigma(x,y)\bigr).$$
A V23 (3) egyenlete alapján a Lagrange-multiplikátor a megoldáson ($\phi_L = \gamma_L\bar{\rho}_L$, $\mathcal{L}_{\text{eff}} = -f_L(\phi_L)\varepsilon$):
$$\lambda_L(x) \;=\; \frac{8\pi G\gamma_L}{c^4}\,\frac{\partial\bigl(-e^{-\phi_L(x)}\varepsilon(x)\bigr)}{\partial\phi_L(x)} \;=\; \frac{8\pi G\gamma_L}{c^4}\,f_L(x)\,\varepsilon(x).$$
Amikor a $S_{\text{I}} = \int dt\int d^3x\,N(x)\sqrt{h(x)}\,\bigl[-f_L(x)\varepsilon(x)\bigr]$ hatást a térbeli $h_{ij}(z)$ metrika szerint variáljuk (rögzített $N(x)$ és $\varepsilon(x)$ mellett), a $\delta\phi_L(x) = \gamma_L\,\delta\bar{\rho}_L(x)$ variációban három különálló metrikus járulék jelenik meg:
$$\delta\bar{\rho}_L(x) \;=\; \int_{\Sigma_u} d^3y\,\delta\Bigl[\sqrt{h(y)}\,\mathcal{K}_L(x,y)\Bigr]\varepsilon(y) \;=\; \delta\bar{\rho}_L^{(\text{vol})}(x) \;+\; \delta\bar{\rho}_L^{(\text{norm})}(x) \;+\; \delta\bar{\rho}_L^{(\text{bdry})}(x).$$

1. **A $\sqrt{h(y)}$ térfogatelem közvetlen variációja ($\delta\bar{\rho}_L^{(\text{vol})}$):**
   Mivel $\delta\sqrt{h(y)} = \frac{1}{2}\sqrt{h(y)}\,h^{ij}(y)\,\delta h_{ij}(y)$:
   $$\delta\bar{\rho}_L^{(\text{vol})}(x) \;=\; \frac{1}{2}\int_{\Sigma_u} d^3y\,\sqrt{h(y)}\,\mathcal{K}_L(x,y)\,\varepsilon(y)\,h^{ij}(y)\,\delta h_{ij}(y).$$
2. **A $V_L(x)$ normálási térfogat variációja ($\delta\bar{\rho}_L^{(\text{norm})}$):**
   Mivel $\delta\bigl(V_L(x)^{-1}\bigr) = -V_L(x)^{-2}\,\delta V_L(x)$, a térfogatelem-variáció a normálásból a helyi $\bar{\rho}_L(x)$ átlaggal súlyozott levonást ad:
   $$\delta\bar{\rho}_L^{(\text{norm})}\big|_{\text{vol}}(x) \;=\; -\bar{\rho}_L(x)\,\frac{1}{2}\int_{\Sigma_u} d^3y\,\sqrt{h(y)}\,\mathcal{K}_L(x,y)\,h^{ij}(y)\,\delta h_{ij}(y).$$
   A két térfogatelem-tag összege önmagában is figyelemre méltó: **homogén közegben ($\varepsilon(y) = \bar{\rho}_L(x)$) a két térfogatelem-variáció ($\delta\bar{\rho}_L^{(\text{vol})} + \delta\bar{\rho}_L^{(\text{norm})}\big|_{\text{vol}}$) egzaktul kiejti egymást!**

### 2.2. A geodetikus távolság $\delta d_\Sigma(x,y)$ variációja (Synge-féle világfüggvény) és a $\delta_D(L - d_\Sigma)$ határfeszültség
Jelölje $\gamma_{xy}: [0, d_\Sigma(x,y)] \to \Sigma_u$, $l \mapsto z^k(l)$ az $x$ és $y$ pontokat összekötő, $l$ saját-ívhosszal paraméterezett minimális térbeli geodetikus görbét, amelynek egységnyi érintővektora $\hat{t}^k(l) = \frac{dz^k}{dl}$ ($h_{ij}\hat{t}^i\hat{t}^j = 1$). A Synge-féle térbeli világfüggvény ($\sigma_\Sigma(x,y) = \frac{1}{2}d_\Sigma(x,y)^2$) első variációja rögzített végpontok ($x, y$) mellett a geodetikus egyenlet miatt nem függ a geodetikus pálya $\delta z^k(l)$ eltolódásától, hanem kizárólag a metrika $\delta h_{ij}$ megváltozásának a geodetikus menti vonalintegrálja:
$$\boxed{\;\delta d_\Sigma(x,y) \;=\; \frac{1}{2}\int_0^{d_\Sigma(x,y)} dl\,\hat{t}^i(l)\,\hat{t}^j(l)\,\delta h_{ij}\bigl(\gamma_{xy}(l)\bigr) \;=\; \frac{1}{2}\int_{\Sigma_u} d^3z\,\delta^{(3)}\!\bigl(z - \gamma_{xy}(l)\bigr)\,\hat{t}^i(l)\hat{t}^j(l)\,\delta h_{ij}(z)\,dl.\;}$$
Mivel a Heaviside-féle lépcsőfüggvény deriváltja a Dirac-delta ($\frac{\partial}{\partial d_\Sigma}\Theta(L - d_\Sigma) = -\delta_D(L - d_\Sigma)$), az ablakhatár elmozdulása egy **gömbfelszíni ($d_\Sigma(x,y) = L$) határfeszültséget** állít elő:
$$\boxed{\;\delta\bar{\rho}_L^{(\text{bdry})}(x) \;=\; -\frac{1}{2 V_L(x)}\oint_{d_\Sigma(x,y)=L} dA_y\,\bigl[\varepsilon(y) - \bar{\rho}_L(x)\bigr]\int_0^L dl\,\hat{t}^i(l)\hat{t}^j(l)\,\delta h_{ij}\bigl(\gamma_{xy}(l)\bigr),\;}$$
ahol $dA_y$ a $d_\Sigma(x,y) = L$ geodetikus gömb saját-felszíneleme, az $\varepsilon(y) - \bar{\rho}_L(x)$ tényező pedig a $\Theta(L - d_\Sigma)$ számlálóbeli és a $V_L(x)$ nevezőbeli határvariációjának egyesítéséből származik.

### 2.3. Kiértékelés gömbszimmetrikus neutroncsillag-belsőben ($h_{rr} = e^{2\lambda(r)} = B(r)^{-1}$) és kapcsolat a II. (Yukawa) megfogalmazással
Gömbszimmetrikus metrikában ($ds_\Sigma^2 = e^{2\lambda(r)}dr^2 + r^2 d\Omega^2$, ahol $e^{-2\lambda(r)} = B(r) = 1 - 2m(r)/r$) a $\delta h_{rr}(r) = 2 e^{2\lambda(r)}\delta\lambda(r)$ radiális metrikavariáció az $(rr)$ Einstein-egyenletben ($G^r_{\;r} = 8\pi T^r_{\;r,\text{tot}}$) az alábbi nemlokális radiális ablak-nyomást ($P_{r,\text{window}}(r)$) adja hozzá a `bilateral_lapse_weighted` egyenletekhez:
$$P_{r,\text{window}}(r) \;=\; -\frac{1}{4\pi r^2 N(r) e^{\lambda(r)}}\,\frac{\delta S_{\text{I},X}}{\delta\lambda(r)}.$$
Ez a kifejezés megvilágítja a V23 kézirat két alapvető megállapítását:
1. **Miért tér el a síkbeli maggal számolt `bilateral_lapse_weighted` a (2) hatás EL-megoldásától?**
   Mert a síkbeli mag ($d_{\text{flat}}(x,y) = |\mathbf{x} - \mathbf{y}|$) elhagyja az $e^{\lambda(r)} = (1 - 2m/r)^{-1/2}$ görbült térfogatelem-tényezőt és a $\delta d_\Sigma(x,y)$ geodetikus határ-vonalintegrált.
2. **Miért küszöböli ki a II. (Helmholtz–Yukawa) megfogalmazás (`yukawa_action_el`) a $\delta_D(L - d_\Sigma)$ felszíni disztribúciós határfeszültséget?**
   Amikor az éles $\Theta(L - d_\Sigma)$ top-hat ablakot a sima 3D Helmholtz–Yukawa rezolvenssel ($(-\Box_\Sigma + L^{-2})^{-1}$) helyettesítjük, a geodetikus vonalintegrálokból származó nemlokális határfeszültség egzaktul átalakul a 7. lépésben (`step7_neutron_star_tov_stability.md`) levezetett **lokális, folytonosan differenciálható skalár-nyomássá**:
   $$p_{r,\phi}(r) \;=\; \frac{1}{2\alpha_Y}\left(e^{-2\lambda(r)}\bigl[\phi_L'(r)\bigr]^2 - \frac{\phi_L(r)^2}{L^2}\right),$$
   amelyben a $g^{rr} = e^{-2\lambda(r)} = B(r)$ tényező variációja ($\delta g^{rr} = -2e^{-2\lambda}\delta\lambda$) automatikusan és algebrailag tartalmazza a görbült téri távolságmérték összes metrikus visszahatását!

---

## 3. Ellenőrzések (Checks)

1. **Homogén sűrűségű határeset ($\varepsilon(y) \equiv \varepsilon_0 = \text{konst.}$):**
   - Ha $\varepsilon(y) = \varepsilon_0$ állandó a $d_\Sigma(x,y) \le L$ gömbön belül, akkor $\bar{\rho}_L(x) = \varepsilon_0$, így mind a térfogatelem-variáció ($\varepsilon(y) - \bar{\rho}_L(x) = 0$), mind a $\delta_D(L - d_\Sigma)$ felszíni határfeszültség ($\varepsilon(y)\big|_{d_\Sigma=L} - \bar{\rho}_L(x) = 0$) **azonosan eltűnik**: $\delta\bar{\rho}_L(x) \equiv 0$ $\checkmark$.
2. **Sík téridőbeli ($\lambda(r) \to 0$, $N(r) \to 1$) határeset:**
   - $h_{ij} \to \delta_{ij} \implies d_\Sigma(x,y) \to |\mathbf{x} - \mathbf{y}|$, $V_L(x) \to \frac{4}{3}\pi L^3$, és a variációs egyenlet visszaadja a 10. lépésben levezetett `bilateral_eom` / `bilateral_lapse_weighted` kifejezéseket $\checkmark$.
