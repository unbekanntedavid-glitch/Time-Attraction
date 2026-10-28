# 12. lépés: A skalár–tenzor / Yukawa hatás második variációja ($\delta^2 S$), a csatolt $(\xi, \delta\phi_L)$ radiális pulzációs operátor explicit együtthatói és a pozitív definit fizikai kinetikus Hesse-mátrix levezetése
# Step 12: Second Variation ($\delta^2 S$) of the Scalar-Tensor / Yukawa Action, Explicit Coefficients of the Coupled $(\xi, \delta\phi_L)$ Radial Pulsation Operator, and Derivation of the Positive-Definite Physical Kinetic Hessian

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step12_second_variation_coupled_pulsation_hessian.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step11_chameleon_screening_multiscale_eft_synthesis.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p3.tex` (8. fejezet, 13–23. sor, Eq. 12, és 11. fejezet **(b)** feladat, 102. sor; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.2. szakasz, 65. sor és 2. szakasz, 120–153. sor, Eq. 8–9, 2. táblázat; DOI: `10.5281/zenodo.23639084`)  
- V21 Zenodo rekord: DOI `10.5281/zenodo.23602501`  
- Numerikus implementáció és ellenőrzés a `yukawa_action_el` sorozaton: `curved_yukawa_adm_and_stability.py`, `curved_yukawa_coupled_stability.csv` és `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`

---

## 1. Szó szerinti idézetek a V22 monográfiából és a V23 (`50aa2db0`) kéziratból

A 7. lépésben (`step7_neutron_star_tov_stability.md`) levezettük a görbült téridőbeli `yukawa_action_el` TOV-háttéregyenleteket, valamint a V21-háttér körüli **statikus ($\omega_0^2 = 0$) nemlokális érintőmodellt** ($\delta f_L[\delta h]$). Mind a V22 monográfia 11. fejezetének **(b)** nyitott feladata (`hu_p3.tex`, 102. sor), mind a V23 kézirat 2. szakasza (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 134–136. sor) kifejezetten megkülönbözteti a statikus érintőmodellt (és annak Gram-mátrixát) a **valódi, a hatás $\delta^2 S$ második variációjából levezetett csatolt dinamikai pulzációs rendszertől és fizikai kinetikus Hesse-mátrixtól**:

1. **A V22 monográfia 11. fejezetének (b) elméleti feladata (`hu_p3.tex`, 102. sor):**
   > *„**(b) A csatolásból levezetett valódi módosított radiális stabilitásvizsgálat:** Amint az (a) feladat megadja a $\delta\phi$ és $\delta T^{\mu\nu}_X$ perturbációk dinamikus mezőegyenleteit, a Chandrasekhar-féle radiális pulzációs operátort variációs úton, a hatás $\delta^2 S$ második variációjából kell újra levezetni, kiküszöbölve az $L \ge 367\,\text{m}$ felett megjelölt nem-variációs metrikus operátor-patológiát és feloldva a befagyasztott-$f$ ($\delta f_L = 0$) közelítést.”*

2. **A V21 `gradient`/`metric` operátor $L \ge 367\text{ m}$ feletti nem-variációs jellege (`hu_p3.tex`, 8. fejezet, 21. sor):**
   > *„**`gradient`/`metric` operátor ($L \ge 367\,\text{m}$ felett nem-variációsként megjelölve):** A $\nu' = 2\Phi_{\text{eff}}'$ metrikus lapse-deriváltat helyettesíti a $\mathcal{Q}(r)$ standard GR görbületi tagjaiba anélkül, hogy figyelembe venné az X-szektorbeli $\delta p_{tX}$ tangenciális feszültség perturbációját. Auditom nyíltan megjelöli ezt az operátort mint **nem pozitív definitet és nem-variációst $L \ge 367\,\text{m}$ felett 12 sorozaton (609 modellen)**, ahol a kompenzálatlan maradéktagok már az emelkedő $\partial M/\partial\rho_c > 0$ ágon is nemfizikai $\omega_0^2 < 0$ sajátértékeket okoznak.”*

3. **A V23 statikus érintőmodell tehetetlenség-invarianciája és a Gram-mátrix vs. fizikai kinetikus Hesse-mátrix különbsége (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 2. szakasz, 134–136. sor):**
   > *„Mivel $\omega_0^2=0$-nál a feladat $\mathcal{H}_{\text{lin}}[\delta h]=0$-ra redukálódik, ez az egyezés a statikus nemlokális linearizációt igazolja, nem a dinamikai tehetetlenséget: a $W\to 4W$ skálázás $\log_{10}\rho_{0c}=15{,}30427$-nél $\omega^2_{4W}/\omega^2_W = 0{,}2500$ arányt ad ($\omega^2_W=4{,}897\times 10^{6}\,\text{s}^{-2}$), azonos $F(0)=0{,}023252$ lövőfüggvény mellett, így a zérushely tehetetlenség-invariáns. A nyers BVP $W$ mátrix nem szimmetrikus; a szimmetrizált Galerkin-alak $\mathbf{M}_{\text{sym}} = \mathbf{D}^{T}\mathrm{diag}(w\,\Delta r)\mathbf{D}$ szimmetrikus, $\min\mathrm{diag}\,\mathbf{M}_{\text{sym}} = 5{,}68\times 10^{-9}>0$; ez Gram-mátrix-tulajdonság, nem a fizikai kinetikus Hesse-mátrix pozitív definitségének bizonyítéka. Az `omega0_sq_sym_tangent_s2` oszlop az elsőrendű érintő-becslést ($\omega^2\simeq -F(0)/F'(0)$) adja.”*

---

## 2. Matematikai levezetés: A $\delta^2 S$ második variáció, a metrikus kényszerek kiküszöbölése és az explicit együtthatók

### 2.1. Időfüggő gömbszimmetrikus perturbációk és a $(\delta\lambda, \delta m, \delta\nu')$ metrikus kényszerek
Tekintsük a V23 II. megfogalmazásának (`yukawa_action_el`, **Branch II-Fluid**, $G = c = 1$) kovariáns hatását időfüggő, gömbszimmetrikus téridőben ($ds^2 = -e^{2\nu(r,t)}dt^2 + e^{2\lambda(r,t)}dr^2 + r^2 d\Omega^2$, ahol a sztatikus háttéren $B_0(r) = e^{-2\lambda_0(r)} = 1 - 2m_0(r)/r$):
$$S[g_{\mu\nu}, \phi_L, \Psi_{\text{fluid}}] \;=\; \int dt\int_0^\infty dr\,4\pi r^2 e^{\nu+\lambda}\left[\frac{R}{16\pi} - \frac{1}{2\alpha_Y}\left(-e^{-2\nu}(\dot{\phi}_L)^2 + e^{-2\lambda}(\phi_L')^2 + \frac{\phi_L^2}{L^2}\right) - e^{-\phi_L}\varepsilon(n)\right].$$
Vezessük be a sztatikus egyensúlyi háttér ($\nu_0(r), \lambda_0(r), \phi_{L,0}(r), \varepsilon_0(r), P_0(r)$) körüli kis radiális Lagrange-féle folyadékelmozdulást $\xi(r,t)$ (amelyre az entalpiaperturbáció $\delta h(r,t) = \frac{\delta P(r,t)}{\varepsilon_0(r) + P_0(r)}$) és az Euler-féle skalármező-perturbációt $\delta\phi_L(r,t)$:
$$r_{\text{fluid}}(r,t) \;=\; r + \xi(r,t), \qquad \phi_L(r,t) \;=\; \phi_{L,0}(r) + \delta\phi_L(r,t).$$

1. **A $(tr)$ és $(tt)$ Einstein-egyenletekből adódó metrikus kényszerek ($\delta\lambda$ és $\delta m$):**
   A gömbszimmetrikus $(tr)$ Einstein-egyenlet ($\delta G_t^{\;r} = 8\pi\,\delta T_t^{\;r}$) idő szerinti integrálásából a belső térben ($0 < r < R_*$):
   $$\frac{\partial(\delta\lambda)}{\partial t} \;=\; -4\pi r e^{2\lambda_0}\left[f_{L,0}(\varepsilon_0 + P_0)\,\dot{\xi} - \frac{1}{\alpha_Y}\phi_{L,0}'\,\delta\dot{\phi}_L\right],$$
   míg a $(tt)$ Einstein-egyenletből ($\delta m'(r,t) = 4\pi r^2\,\delta(-T_t^{\;t})$, ahol $\delta B = -2\delta m / r = -2 B_0\,\delta\lambda$) rögzített barionszám ($\delta N_B = 0$) mellett a tömegperturbáció parciális deriváltja az $(\delta h, \delta\phi_L, \delta\phi_L', \delta m)$ változókban:
   $$\boxed{\;\frac{\partial(\delta m)}{\partial r} \;=\; 4\pi r^2\left[f_{L,0}\frac{\varepsilon_0+P_0}{c_s^2}\,\delta h - f_{L,0}\varepsilon_0\,\delta\phi_L + \frac{1}{\alpha_Y}\left(B_0\phi_{L,0}'\,\delta\phi_L' + \frac{\phi_{L,0}}{L^2}\delta\phi_L - \frac{(\phi_{L,0}')^2}{r}\delta m\right)\right].\;}$$
   Ekvivalens módon, a Lagrange-féle elmozdulással $\xi(r)$ felírva (a barionszám-megmaradási kényszer $\delta n = -\frac{1}{r^2 e^{\lambda_0}}(r^2 e^{\lambda_0} n_0\xi)' - n_0\delta\lambda$ behelyettesítésével és parciális integrálással) a tömegperturbáció:
   $$\delta m(r) \;=\; -4\pi r^2 f_{L,0}(\varepsilon_0 + P_0)\,\xi(r) \;+\; \int_0^r 4\pi \tilde{r}^2\left[\frac{B_0\phi_{L,0}'}{\alpha_Y}\delta\phi_L' + \left(\frac{\phi_{L,0}}{\alpha_Y L^2} - f_{L,0}\varepsilon_0\right)\delta\phi_L + f_{L,0}(\varepsilon_0+P_0)\phi_{L,0}'\,\xi - \frac{(\phi_{L,0}')^2}{\alpha_Y\tilde{r}}\delta m\right]d\tilde{r}.$$
2. **Az $(rr)$ Einstein-egyenletből adódó lapse-gradiens kényszer ($\delta\nu'$):**
   Mivel $\nu_0'(r) = \frac{m_0 + 4\pi r^3 P_{\text{tot},0}}{r(r-2m_0)}$, ahol $P_{\text{tot},0} = f_{L,0}P_0 + p_{r,\phi,0}$ és $p_{r,\phi,0} = \frac{1}{2\alpha_Y}[B_0(\phi_{L,0}')^2 - \phi_{L,0}^2/L^2]$, az összesített nyomásperturbáció:
   $$\delta P_{\text{tot}} \;=\; f_{L,0}(\varepsilon_0 + P_0)\,\delta h - f_{L,0}P_0\,\delta\phi_L + \frac{1}{\alpha_Y}\left(B_0\phi_{L,0}'\,\delta\phi_L' - \frac{\phi_{L,0}}{L^2}\delta\phi_L - \frac{(\phi_{L,0}')^2}{r}\delta m\right),$$
   amelyből a lapse-gradiens perturbációja:
   $$\boxed{\;\delta\nu'(r) \;=\; \frac{1 + 2 r \nu_0'(r)}{r\bigl(r - 2m_0(r)\bigr)}\,\delta m(r) \;+\; \frac{4\pi r^3}{r\bigl(r - 2m_0(r)\bigr)}\,\delta P_{\text{tot}}(r).\;}$$

### 2.2. A $\delta^2 S$ másodrendű hatás és a pozitív definit fizikai kinetikus Hesse-mátrix ($\mathbf{T}_{\text{kin}}$)
A metrikus kényszerek visszahelyettesítése után a másodrendű kvadratikus hatás ($\delta^2 S = \int dt\,(T_{\text{kin}}[\dot{\xi}, \delta\dot{\phi}_L] - U_{\text{pot}}[\xi, \delta\phi_L])$) kinetikus része harmonikus időfüggés ($\xi(r,t) = \xi(r)e^{i\omega t}$, $\delta\phi_L(r,t) = \delta\phi_L(r)e^{i\omega t}$) mellett:
$$\boxed{\;T_{\text{kin}}\bigl[\xi, \delta\phi_L\bigr] \;=\; \frac{\omega^2}{2}\int_0^{R_*} dr\,\mathcal{W}_\xi(r)\,\xi(r)^2 \;+\; \frac{\omega^2}{2}\int_0^\infty dr\,\mathcal{W}_\phi(r)\,\bigl[\delta\phi_L(r)\bigr]^2,\;}$$
ahol a két fizikai tehetetlenségi (kinetikus) súlyfüggvény:
$$\boxed{\;\mathcal{W}_\xi(r) \;=\; 4\pi r^2\,e^{3\lambda_0(r) - \nu_0(r)}\,f_{L,0}(r)\,\bigl(\varepsilon_0(r) + P_0(r)\bigr) \;>\; 0 \quad (0 < r < R_*), \qquad \mathcal{W}_\phi(r) \;=\; \frac{4\pi r^2}{\alpha_Y}\,e^{\lambda_0(r) - \nu_0(r)} \;>\; 0 \quad (r > 0).\;}$$

### 2.3. A csatolt $(\xi, \delta\phi_L)$ és $(\delta h, \delta\phi_L, \delta m)$ sajátérték-rendszer explicit együtthatói (Point 3)
A csatolt radiális pulzációs rendszer két maradéktalanul ekvivalens — másodrendű Lagrange-féle $(\xi, \delta\phi_L)$, illetve elsőrendű/másodrendű kanonikus $(\delta h, \delta\phi_L, \delta m)$ — alakban írható fel:

#### (A) Explicit másodrendű csatolt $(\xi, \delta\phi_L)$ Sturm–Liouville / Klein–Gordon alak:
$$\begin{cases}
  -\dfrac{d}{dr}\!\left(\mathcal{P}_\xi(r)\dfrac{d\xi}{dr}\right) + \mathcal{Q}_{\xi\xi}(r)\,\xi(r) + \mathcal{C}_{\xi\phi}\!\left[\delta\phi_L\right](r) \;=\; \omega^2\,\mathcal{W}_\xi(r)\,\xi(r), \\[8pt]
  -\dfrac{d}{dr}\!\left(\mathcal{P}_\phi(r)\dfrac{d(\delta\phi_L)}{dr}\right) + \mathcal{Q}_{\phi\phi}(r)\,\delta\phi_L(r) + \mathcal{C}_{\phi\xi}\!\left[\xi\right](r) \;=\; \omega^2\,\mathcal{W}_\phi(r)\,\delta\phi_L(r),
\end{cases}$$
ahol a főrész-együtthatók:
$$\mathcal{P}_\xi(r) \;=\; 4\pi r^2\,e^{\nu_0(r) + \lambda_0(r)}\,f_{L,0}(r)\,\Gamma_1(r)\,P_0(r) \;>\; 0, \qquad \mathcal{P}_\phi(r) \;=\; \frac{4\pi r^2}{\alpha_Y}\,e^{\nu_0(r) - \lambda_0(r)} \;>\; 0,$$
a folyadék-folyadék potenciál-együttható ($\mathcal{Q}_{\xi\xi}$, a módosított hidrosztatikai gradienssel $P_0' = -(\varepsilon_0+P_0)(\nu_0' - \phi_{L,0}')$ és az összesített háttér-nyomással/energiával $\rho_{\text{tot},0} = f_{L,0}\varepsilon_0 + \rho_{\phi,0}$, $P_{\text{tot},0} = f_{L,0}P_0 + p_{r,\phi,0}$):
$$\boxed{\;\begin{aligned}
  \mathcal{Q}_{\xi\xi}(r) &\;=\; 4\pi r^2 e^{\nu_0+\lambda_0}f_{L,0}\left\{\frac{4}{r}P_0' + \frac{(P_0')^2}{\varepsilon_0 + P_0} - 8\pi e^{2\lambda_0}(\varepsilon_0 + P_0)\bigl(f_{L,0}P_0 + p_{r,\phi,0}\bigr)\right. \\
  &\qquad\left. -\; \frac{d}{dr}\!\left(\frac{2\Gamma_1 P_0}{r}\right) + (\varepsilon_0 + P_0)\left[\phi_{L,0}'' + \bigl(\nu_0' + \lambda_0' - \phi_{L,0}'\bigr)\phi_{L,0}'\right]\right\},
\end{aligned}\;}$$
a skalár-skalár potenciál-együttható ($\mathcal{Q}_{\phi\phi}$, beleértve a $\delta\nu'$-ből és $\delta B$-ből származó metrikus visszahatást):
$$\boxed{\;\mathcal{Q}_{\phi\phi}(r) \;=\; \frac{4\pi r^2}{\alpha_Y}\,e^{\nu_0+\lambda_0}\left[\frac{1}{L^2} + \alpha_Y\,f_{L,0}(r)\,\varepsilon_0(r) + \frac{8\pi r\,e^{2\lambda_0}}{\alpha_Y L^2}\phi_{L,0}(r)\phi_{L,0}'(r)\right] \;+\; \mathcal{Q}_{\phi\phi}^{\text{metric}}[\delta m],\;} $$
és a kölcsönösen adjungált folyadék–skalár keresztcsatolási operátorok ($\langle \eta, \mathcal{C}_{\xi\phi}[\delta\phi_L]\rangle = \langle \mathcal{C}_{\phi\xi}[\eta], \delta\phi_L\rangle$):
$$\boxed{\;\begin{aligned}
  \mathcal{C}_{\xi\phi}[\delta\phi_L](r) &\;=\; 4\pi r^2 e^{\nu_0+\lambda_0}f_{L,0}(\varepsilon_0+P_0)\left[-\delta\phi_L'(r) + \delta\nu'\big|_{\delta\phi_L, \delta m[\delta\phi_L]}(r)\right], \\[4pt]
  \mathcal{C}_{\phi\xi}[\xi](r) &\;=\; -\frac{d}{dr}\!\Bigl(4\pi r^2 e^{\nu_0+\lambda_0}f_{L,0}(\varepsilon_0+P_0)\,\xi(r)\Bigr) \;+\; 4\pi r^2 e^{\nu_0+\lambda_0}f_{L,0}\frac{\varepsilon_0+P_0}{c_s^2}\,\delta h[\xi](r) \;+\; \mathcal{C}_{\phi\xi}^{\text{metric}}[\delta m[\xi]](r).
\end{aligned}\;}$$

#### (B) Egzakt kanonikus $(\delta h, \delta\phi_L, \delta m)$ rendszer és mozgó felszíni illesztési feltételek (implementálva a `curved_yukawa_adm_and_stability.py` kódban):
A numerikus integráláshoz és a politróp-határok ($r_{\text{trans}}$) egzakt kezeléséhez a $(\delta h, \delta\phi_L, \delta m)$ kanonikus változók a legalkalmasabbak (mivel $\delta h$ és $\delta m$ folytonos a $c_s^2$ ugrási pontjain is!):
$$\boxed{\;\begin{aligned}
  \frac{d(\delta h)}{dr} &\;=\; -\delta\nu'(r) \;+\; \delta\phi_L'(r) \;+\; \omega^2\,e^{2(\lambda_0-\nu_0)}\,\xi(r), \\[4pt]
  \frac{d(\delta m)}{dr} &\;=\; 4\pi r^2\left[f_{L,0}\frac{\varepsilon_0+P_0}{c_s^2}\,\delta h - f_{L,0}\varepsilon_0\,\delta\phi_L + \frac{1}{\alpha_Y}\left(B_0\phi_{L,0}'\,\delta\phi_L' + \frac{\phi_{L,0}}{L^2}\delta\phi_L - \frac{(\phi_{L,0}')^2}{r}\delta m\right)\right], \\[4pt]
  B_0\,\delta\phi_L'' &\;+\; \left(\frac{2B_0}{r} + B_0\nu_0' + \frac{B_0'}{2}\right)\delta\phi_L' - \left(\frac{1}{L^2} + \alpha_Y f_{L,0}\varepsilon_0 - \omega^2 e^{-2\nu_0}\right)\delta\phi_L \\
  &\;=\; -\alpha_Y f_{L,0}\frac{\varepsilon_0+P_0}{c_s^2}\,\delta h + \frac{2\,\delta m}{r}\phi_{L,0}'' - \left[\delta(B\nu') + \frac{\delta B'}{2} - \frac{2\,\delta m}{r^2}\right]\phi_{L,0}'.
\end{aligned}\;}$$
- **Mozgó felszíni perem- és illesztési feltételek ($r = R_*$):**
  1. A szabad felszínen a Lagrange-féle nyomásperturbáció eltűnik: $\Delta P(R_*) = \delta P(R_*) + P_0'(R_*)\xi(R_*) = 0 \implies \delta h(R_*) = -h_0'(R_*)\xi(R_*)$, ami az integrált tömegmegmaradással ($\delta m(R_*^-) = -4\pi R_*^2 f_{L,0}(R_*)(\varepsilon_0+P_0)\big|_{R_*^-}\xi(R_*) = 0$, mivel politróp kéregre $\varepsilon_0(R_*^-) = P_0(R_*^-) = 0$) biztosítja $\delta m(R_*) \in C^0$-t.
  2. A mozgó $r = R_* + \xi(R_*)$ felszínen a skalármező és normálderiváltja folytonos:
     $$\llbracket \delta\phi_L + \xi\phi_{L,0}' \rrbracket_{R_*} \;=\; 0, \qquad \llbracket \delta\phi_L' + \xi\phi_{L,0}'' \rrbracket_{R_*} \;=\; 0.$$
     Mivel $\phi_{L,0} \in C^1([0,\infty))$, $\llbracket\phi_{L,0}'\rrbracket_{R_*} = 0 \implies \llbracket\delta\phi_L\rrbracket_{R_*} = 0$; továbbá $\llbracket\phi_{L,0}''\rrbracket_{R_*} = \frac{\alpha_Y}{B_0(R_*)}f_{L,0}(R_*)\varepsilon_0(R_*^-) = 0$ a politróp kéreg felszínén ($\varepsilon_0(R_*^-) = 0$), így $\llbracket\delta\phi_L'\rrbracket_{R_*} = 0$ is egzaktul fennáll!
  3. A külső görbült vákuumban ($r \in [R_*, R_{\max}]$) $\varepsilon_0 = P_0 = 0$, így a $(\delta\phi_L, \delta m)$ rendszer csatoltan integrálódik egészen $R_{\max} = R_* + 12 L$-ig, ahol az aszimptotikus Yukawa-lecsengési feltétel zárja le a peremérték-feladatot, és az aszimptotikus tömegperturbáció $\delta M_{\text{ADM}}(\omega^2) = \delta m(R_{\max}; \omega^2) = 0$ határozza meg a sajátértéket!

---

## 3. Numerikus ellenőrzés a görbült külső `yukawa_action_el` sorozaton (Point 4)

A `curved_yukawa_adm_and_stability.py` szkript közvetlenül a görbült külső megoldással rendelkező `yukawa_action_el` egyensúlyi sorozatokon oldja meg a fenti csatolt $(\delta h, \delta\phi_L, \delta m)$ sajátérték-feladatot ($r \in [0, R_* + 14L]$, az egzakt bal oldali `scalar_m_sign = -1.0` előjellel). Az eredmények (`curved_yukawa_coupled_stability.csv` és `curved_yukawa_audit_summary.json`):

| Sorozat | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dinamikus $\omega_0^2=0$ | $\Delta_{\text{ADM}}$ dinamikus (dex) | $\Delta_{M_B}$ dinamikus (dex) | Kvázisztatikus $\omega_0^2=0$ | Befagyasztott $\delta\phi_L=0$ | $\Delta_{\text{ADM}}$ befagyasztott (dex) |
|---|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316409` | **`15.316410`** | **`+5.80e-06`** | `+1.33e-06` | `15.316410` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.307169` | **`15.307170`** | **`+4.20e-06`** | `+6.41e-07` | `15.307170` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303082` | **`15.303080`** | **`-6.86e-06`** | `-2.61e-06` | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.293025` | **`15.293025`** | **`+3.57e-06`** | `-1.40e-07` | `15.293025` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227787` | **`15.227783`** | **`-1.34e-06`** | `-4.05e-06` | `15.227783` | `15.251384` | `+2.36e-02` |

*(GL2 Volterra, $N=480$, `scalar_m_sign = -1.0`, dinamikus skalárhullám, véges differenciákkal igazolt külső $\partial\phi_L'(R_*)/\partial m(R_*)=-1{,}68365\times 10^{-4}\;\text{km}^{-2}$ és belső $r=6\;\text{km}$ $\partial\phi_L''/\partial m=-1{,}72833\times 10^{-4}\;\text{km}^{-3}$ deriváltak; $N$- és $R_{\max}$-konvergencia a gyökön: $p_{\text{eff,root}}=4{,}248$, $\log_{10}\rho_{0c,\text{Rich}}(\omega_0^2=0)=15{,}31640981$, lásd `step22`, 4. és 7. szakasz.)*

---

## 4. Ellenőrzések (Checks)

1. **$f_{L,0} \to 1$, $\alpha_Y \to 0$ (GR-határeset):**
   - $\phi_{L,0} \to 0$ és $\delta\phi_L \to 0$ mellett $\mathcal{W}_\xi(r) \to 4\pi r^2 e^{3\lambda_0-\nu_0}(\varepsilon_0 + P_0)$, $\mathcal{P}_\xi(r) \to 4\pi r^2 e^{\nu_0+\lambda_0}\Gamma_1 P_0$ és $\mathcal{Q}_{\xi\xi}(r)$ visszaadja a klasszikus Chandrasekhar (1964) általános relativitáselméleti radiális pulzációs Sturm–Liouville-operátort $\checkmark$.
2. **A csatolt $(\delta h, \delta\phi_L, \delta m)$ sajátérték-zérushely és az ADM-/barionos tömegmaximum összehasonlítása a `yukawa_action_el` ágon:**
   - SLy $L = 10\;\text{km}$, $\alpha_Y = 1$ esetén a befagyasztott-$\phi_L$ ($\delta\phi_L = 0$) közelítés $+1{,}22\times 10^{-2}\;\text{dex}$-szel tolja el a zérushelyet a $dM_{\text{ADM}}/d\rho_c = 0$ fordulóponttól ($15{,}328565$ vs. $15{,}316404$), míg a csatolt $(\delta h, \delta\phi_L, \delta m)$ operátor (`scalar_m_sign = -1.0`) $\log_{10}\rho_{0c}(\omega_0^2=0)=15{,}316410$-et ad ($+5{,}80\times 10^{-6}\;\text{dex}$ az ADM-csúcsra és $+1{,}33\times 10^{-6}\;\text{dex}$ a barionos tömegcsúcsra, $2097\times$-es csökkenés a befagyasztott-$\phi_L$ eltolódáshoz képest; mind az 5 sorozaton $\le 6{,}9\times 10^{-6}\;\text{dex}$ egyezés, $165\times\text{–}17593\times$ csökkenés) $\checkmark$.
3. **A kinetikus súlyok pozitivitása $L \ge 367\text{ m}$ felett:**
   - Mivel $f_{L,0}(r) = e^{-\phi_{L,0}(r)} > 0$ és $\alpha_Y > 0$, $\mathcal{W}_\xi(r) > 0$ és $\mathcal{W}_\phi(r) > 0$ minden $L \in (0, \infty)$ skálán fennáll, előjelváltás nélkül $\checkmark$.
