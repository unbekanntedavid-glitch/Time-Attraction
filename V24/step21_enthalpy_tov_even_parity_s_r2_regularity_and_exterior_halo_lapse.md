# 21. lépés: A 7 szegmenses politróp log-entalpia ($h$) TOV-transzformáció, az $s = r^2$ páros paritású centrális regularizáció (`Audit Fix O1`) és a `gradient` módú $r \in [R_*, R_* + L]$ külső vákuum-halo lapse-integrál ($\Delta\nu_{\text{ext}}$) egzakt levezetése
# Step 21: Exact Derivation of the 7-Segment Piecewise-Polytrope Log-Enthalpy ($h$) TOV Transformation, the $s = r^2$ Even-Parity Central Regularization (`Audit Fix O1`), and the `gradient`-Mode Exterior Vacuum-Halo Lapse Integral ($\Delta\nu_{\text{ext}}$) over $r \in [R_*, R_* + L]$

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step20_eight_coil_elliptic_biot_savart_and_exterior_compactified_quadrature.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p2.tex` (7. fejezet, 102–120. sor; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.3. szakasz és 1. táblázat, 109–156. sor; DOI: `10.5281/zenodo.23639084`)  
- V21/V23 profil- és TOV-integrátor: `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_profiles.py` (14–108. sor) és `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/stellar_tal_v21.py` (33–136. sor; DOI: `10.5281/zenodo.23602501`)

---

## 1. Szó szerinti idézetek a V22 és V23 kéziratokból, valamint a `tal_profiles.py` modulból

A 7., 16., 17. és 18. lépésben felhasználtuk a V21/V23 neutroncsillag-egyensúlyi profilokat (`dense_profile`, SLy, APR4, H4 állapotegyenletek, $1{,}82\times 10^{-11}$ és $1{,}59\times 10^{-11}$ folytonos Bianchi-reziduumok), de még nem vezettük le azt a három analitikus transzformációt a `stellar_tal_v21.py` és a `tal_profiles.py` modulokban, amelyek ezt a pontosságot és a centrális szingularitásmentességet biztosítják:

1. **A 3 állapotegyenlet (SLy, APR4, H4), az oksági $c_s = c$ határ és a két feltételes TOV-mód (`hu_p2.tex`, 7. fejezet, 108–115. sor):**
   $$\begin{equation}
     \textbf{(A) \texttt{mass} mód:}\quad \frac{dm_{\text{eff}}}{dr} = 4\pi r^2 \rho(r)\,f_L(r), \qquad
     \textbf{(B) \texttt{gradient} mód:}\quad \frac{d\Phi_{\text{eff}}}{dr} = f_L(r)\,\frac{G\bigl(m(r) + 4\pi r^3 P(r)/c^2\bigr)}{c^2 r\bigl(r - 2Gm(r)/c^2\bigr)}, \tag{V22 7. fejezet}
   \end{equation}$$
   > *„33 egyensúlyi sorozatot integráltam három nukleáris állapotegyenletre — **SLy** (nukleonikus), **APR4** (variációs nukleonikus; az oksági $c_s = c$ hangsebességhatárt **$\rho_0 \approx 1{,}4458\times 10^{15}\,\text{g/cm}^3$** nyugalmi tömegsűrűségnél éri el) és **H4** (merev hiperonikus) — 11 logaritmikusan elosztott $L \in [1\,\text{m}, 100\,\text{km}]$ simítási skálán. A V21 archív alaptáblázata **182 egyensúlyi modellt** tartalmaz, míg a sűrű diagnosztikai stabilitási rács **összesen 769 csillagmodellből** áll (beleértve mind a **33 egzakt oksági végpontot**, ahol $c_s/c = 1$).”*

2. **Az $s = r^2$ páros paritású centrális leképezés és a `gradient` módú külső lapse-integrál (`tal_profiles.py`, 53–56. és 91–93. sor):**
   - 53–56. sor:
     > *„`# Profile mapping, regular at the centre (audit fix O1, run 2): h and mu = m/r^3 are even, analytic functions of r, so they are interpolated in s = r^2 with their exact centre values h(0)=hc, mu(0)=(4pi/3) eps_c f_m. Then h'(r)=2r dh/ds -> 0, m=r^3 mu, A=r(mu+4pi P)/(1-2 mu r^2) are regular. (Run 1 interpolated h(r), m(r) in r with a single node at r=0 -> wrong slopes near the centre.)`”*
   - 91–93. sor:
     > *„`def exterior_lapse(cur, L, M, R, nq=160): """nu(R) = -int_R^inf f(r) M/(r(r-2M)) dr with f=exp(-L^2 rhobar_L) from the converged profile; f=1 beyond R+L. Returns nu(R) and diagnostic dict."""`”*

3. **A V23 1. táblázatának folytonos Bianchi-megmaradási reziduumai (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 144–147. sor):**
   - `tophat_mass`: $1{,}82\times 10^{-11}$, `tophat_gradient`: $1{,}59\times 10^{-11}$ (az egzakt szimbolikus maradék `0.0`).

---

## 2. A Read et al. (2009) 7 szegmenses politróp log-entalpia ($h$) lezárás és a TOV-rendszer $h$-transzformációja

### 2.1. Termodinamikai konzisztencia, $a_i$ folytonossági rekurzió és a $c_s^2 / c^2 = 1$ oksági küszöb (`stellar_tal_v21.py`, 33–69. sor)
A Read et al. (2009) paraméterezés 4 héj-szegmensből ($\Gamma_0,\dots,\Gamma_3$) és 3 mag-szegmensből ($\Gamma_4,\Gamma_5,\Gamma_6 = g_1,g_2,g_3$, illesztési sűrűségek $\rho_{0,4} = \rho_{\text{match}}$, $\rho_{0,5} = 10^{14{,}7}\;\text{g/cm}^3$, $\rho_{0,6} = 10^{15}\;\text{g/cm}^3$) épül fel. Az $i$-edik szegmensben ($\rho_0 \in [\rho_{0,i-1}, \rho_{0,i}]$) a nyomás (a $c^2$-tel osztott CGS-egységekben) $P(\rho_0) = K_i \rho_0^{\Gamma_i}$. A nyomás folytonossága a szegmenshatárokon:
$$K_{i+1} \;=\; K_i\,\rho_{0,i}^{\Gamma_i - \Gamma_{i+1}}, \qquad \rho_{\text{match}} \;=\; \left(\frac{K_3}{K_4}\right)^{\frac{1}{\Gamma_4 - \Gamma_3}}.$$
A termodinamika első főtétele barionmegmaradás mellett ($d(\varepsilon / \rho_0) = -P\,d(1/\rho_0) = \frac{P}{\rho_0^2}d\rho_0$) integrálva megadja a tömeg-energia sűrűséget (`stellar_tal_v21.py`, 47–49. és 67. sor):
$$\varepsilon(\rho_0) \;=\; (1 + a_i)\,\rho_0 \;+\; \frac{P(\rho_0)}{\Gamma_i - 1}, \qquad a_{i+1} \;=\; a_i \;+\; \frac{K_i}{\Gamma_i - 1}\rho_{0,i}^{\Gamma_i - 1} \;-\; \frac{K_{i+1}}{\Gamma_{i+1} - 1}\rho_{0,i}^{\Gamma_{i+1} - 1} \quad (a_0 = 0).$$
A fajlagos relativisztikus log-entalpia ($dh \equiv \frac{dP}{\varepsilon + P}$) ebből analitikusan integrálható és invertálható (`stellar_tal_v21.py`, 50. és 60–68. sor):
$$\boxed{\;h(\rho_0) \;=\; \ln\!\left(\frac{\varepsilon + P}{\rho_0}\right) \;=\; \ln\!\left(1 + a_i + \frac{\Gamma_i}{\Gamma_i - 1}K_i\,\rho_0^{\Gamma_i - 1}\right) \;\Longleftrightarrow\; \rho_0(h) \;=\; \left[\frac{\bigl(e^h - 1 - a_i\bigr)(\Gamma_i - 1)}{\Gamma_i K_i}\right]^{\frac{1}{\Gamma_i - 1}}.\;}$$
Az adiabatikus hangsebesség négyzete minden szegmensben:
$$\frac{c_s^2(h)}{c^2} \;=\; \left(\frac{\partial P}{\partial\varepsilon}\right)_s \;=\; \frac{\Gamma_i\,P(h)}{\varepsilon(h) + P(h)} \;=\; (\Gamma_i - 1)\,\frac{e^h - 1 - a_i}{e^h}.$$
Mivel $\Gamma_i > 2$ a nagy sűrűségű magban (APR4-re $\Gamma_5 = 3{,}445$, $\Gamma_6 = 3{,}348$), a $c_s^2(h)/c^2$ szigorúan monoton nő $h$-val, és az oksági határt ($c_s^2/c^2 = 1$) a `brentq` gyökkereső $10^{-15}$ pontossággal határozza meg a $h_{\text{causal}}$ entalpiánál (APR4-re $\rho_0(h_{\text{causal}}) \approx 1{,}4458\times 10^{15}\;\text{g/cm}^3$, `hu_p2.tex`, 115. sor).

### 2.2. A TOV-egyenletek átírása $h$ független változóra és a centrális Taylor-indítás (`tal_profiles.py`, 38–50. sor)
Mivel a hidrosztatikai egyensúlyi egyenlet $P'(r) = (\varepsilon + P)h'(r) = -(\varepsilon + P)f_{\text{grad}}(r)A_{\text{GR}}(r)$, a $(\varepsilon + P)$ faktor kiesik:
$$\frac{dh}{dr} \;=\; -f_{\text{grad}}(r)\,A_{\text{GR}}(r) \;=\; -f_{\text{grad}}(r)\,\frac{m(r) + 4\pi r^3 P(h)}{r\bigl(r - 2m(r)\bigr)},$$
ahol $f_{\text{grad}}(r) = f_L(r)$ a `gradient` módban és $1$ a `mass` (valamint GR) módban, míg $f_{\text{mass}}(r) = f_L(r)$ a `mass` módban és $1$ a `gradient` módban. Mivel $h(r)$ a csillag belsejében szigorúan monoton csökken a centrális $h_c$-től a felszíni $h(R_*) = 0$-ig, a $(r(h), m(h))$ függvényekre vett közönséges differenciálegyenlet-rendszer (`tal_profiles.py`, 43–48. sor):
$$\boxed{\;\frac{dr}{dh} \;=\; -\frac{r(r - 2m)}{\bigl(m + 4\pi r^3 P(h)\bigr)\,f_{\text{grad}}(r)}, \qquad \frac{dm}{dh} \;=\; 4\pi r^2 \varepsilon(h)\,f_{\text{mass}}(r)\,\frac{dr}{dh}.\;}$$
**Miért előnyös a $h$-szerinti integrálás?** Ha $r$ a független változó, a csillag $P(R_*) = 0$ felszínét eseménydetektálással (event location) kell megkeresni, ahol $\varepsilon(P) \to 0$ meredeken változik. Ezzel szemben $h$ független változóval a felszín egy előre rögzített végpont ($h_{\text{surf}} = 10^{-10}$), amelynél $R_* = r(h_{\text{surf}})$ és $M = m(h_{\text{surf}})$ közvetlenül, gyökkeresési túllövés nélkül kiolvasódik! A $r = 0$ koordináta-szingularitás elkerülésére $r_0 = 10^{-4}\;\text{km}$-nél a másodrendű centrális Taylor-sor indítja az integrálást (`tal_profiles.py`, 41–42. sor):
$$m(r_0) \;=\; \frac{4\pi}{3}\varepsilon_c f_m\,r_0^3 + \mathcal{O}(r_0^5), \qquad h(r_0) \;=\; h_c - \delta h, \qquad \delta h \;=\; \frac{2\pi}{3}(\varepsilon_c f_m + 3 P_c)\,f_g\,r_0^2 + \mathcal{O}(r_0^4).$$

---

## 3. Az $s = r^2$ páros paritású centrális regularizáció (`Audit Fix O1`, `tal_profiles.py`, 53–79. sor)

### 3.1. Miért okoz $r$-szerinti spline-interpoláció centrális szingularitást, és hogyan javítja ezt az $s = r^2$ változó?
Egy gömbszimmetrikus, az origóban sima ($C^\infty(\mathbb{R}^3)$) skalármező — mint a log-entalpia $h(\mathbf{x})$, az átlagos belső tömegsűrűség $\mu(\mathbf{x}) \equiv \frac{m(r)}{r^3}$, vagy a 3D top-hat kitevő $x(\mathbf{x}) \equiv \gamma_L\bar{\rho}_L(r)$ — csak $|\mathbf{x}|^2 = x^2 + y^2 + z^2 = r^2 \equiv s$ hatványait tartalmazhatja Taylor-sorában, vagyis $r$-ben **páros függvény**:
$$h(r) \;=\; h_s(r^2), \qquad \mu(r) \;\equiv\; \frac{m(r)}{r^3} \;=\; \mu_s(r^2), \qquad x(r) \;=\; x_s(r^2),$$
ahol $h_s(0) = h_c$ és $\mu_s(0) = \frac{4\pi}{3}\varepsilon_c f_m$. Ha $h(r)$-et és $m(r)$-et közvetlenül $r$-ben interpoláljuk PCHIP-spline-nal (`Run 1`), akkor az $r = 0$ csomópontban a harmadfokú Hermite-polinom nem kényszeríti ki sem a $h''(0)$ pontos értékét, sem az $m'(0) = m''(0) = 0$ feltételt, így $\frac{m(r)}{r^3}$ és $\frac{h'(r)}{r}$ az origó közelében numerikus zajt és $\mathcal{O}(1/r)$ hibát kap!

Ezzel szemben a `tal_profiles.py` (64–79. sor, `Audit Fix O1`) a $s = r^2 \in [0, R_*^2]$ rácson építi fel a `PchipInterpolator` függvényeket:
$$h_s(s), \qquad \mu_s(s) \;=\; \frac{m(\sqrt{s})}{s^{3/2}} \quad (\mu_s(0) = \mu_c), \qquad x_s(s) \;=\; x(\sqrt{s}).$$
A láncszabály ($\frac{d}{dr} = 2r\frac{d}{ds}$) alapján ebből minden centrális derivált és görbületi faktor **analitikusan reguláris** az $r \to 0^+$ határesetben (`tal_profiles.py`, 66–79. sor):
$$\boxed{\;\begin{aligned}
  h'(r) &\;=\; 2r\,\frac{dh_s}{ds}\!(r^2) \;=\; \mathcal{O}(r), \qquad &m(r) &\;=\; r^3\,\mu_s(r^2) \;=\; \mathcal{O}(r^3), \\[4pt]
  \frac{d\mu}{dr}(r) &\;=\; 2r\,\frac{d\mu_s}{ds}\!(r^2) \;=\; \mathcal{O}(r), \qquad &\frac{d\ln f_L}{dr}(r) &\;=\; -2r\,\frac{dx_s}{ds}\!(r^2) \;=\; \mathcal{O}(r), \\[4pt]
  A_{\text{GR}}(r) &\;=\; \frac{m(r) + 4\pi r^3 P(r)}{r\bigl(r - 2m(r)\bigr)} \;=\; \frac{r\bigl(\mu_s(r^2) + 4\pi P(r)\bigr)}{1 - 2r^2\,\mu_s(r^2)} \;=\; \mathcal{O}(r).
\end{aligned}\;}$$
Következésképpen a Chandrasekhar-operátorban és a Bianchi-azonosságban fellépő $\frac{h'(r)}{r} = 2\frac{dh_s}{ds}(r^2)$ és $\frac{A_{\text{GR}}(r)}{r} = \frac{\mu_s(r^2) + 4\pi P(r)}{1 - 2r^2\mu_s(r^2)}$ hányadosokból az $r$ faktor **algebrailag kiesik** $r = 0$-ban! Ez a regularizáció csökkenti a V23 1. táblázatában a folytonos Bianchi-megmaradási reziduumot $1{,}82\times 10^{-11}$-re (`tophat_mass`) és $1{,}59\times 10^{-11}$-re (`tophat_gradient`).

---

## 4. A `gradient` módú $r \in [R_*, R_* + L]$ külső vákuum-halo lapse-integrál ($\Delta\nu_{\text{ext}}$) és a belső Bernoulli-integrál (`tal_profiles.py`, 80–106. sor)

### 4.1. Miért nyúlik túl az $f_L(r) < 1$ tartomány a csillag $R_*$ felszínén $R_* + L$-ig?
A csillagon kívül ($r > R_*$) a lokális anyagsűrűség azonosan nulla ($\varepsilon(r) = 0$, $P(r) = 0$), így $m(r) \equiv M$ és $A_{\text{GR}}(r) = \frac{M}{r(r - 2M)}$ minden $r \ge R_*$-re.
- **`mass` módban (és GR-ben):** a külső lapse-gradiens $\nu'(r) = A_{\text{GR}}(r) = \frac{M}{r(r-2M)}$ nem tartalmaz $f_L(r)$ szorzót. Ezért a külső téridő egzaktul Schwarzschild-geometria, és a felszíni lapse (`tal_profiles.py`, 81. sor):
  $$\nu_{\text{mass}}(R_*) \;=\; \nu_{\text{GR}}(R_*) \;=\; -\int_{R_*}^\infty \frac{M}{r(r - 2M)}\,dr \;=\; \frac{1}{2}\ln\!\left(1 - \frac{2M}{R_*}\right).$$
- **`gradient` módban:** a lapse-gradiens $\nu'(r) = f_L(r)\,A_{\text{GR}}(r)$ mindenhol. Bár $\varepsilon(r) = 0$ a csillagon kívül, egy $r \in (R_*, R_* + L)$ külső pont köré írt $L$ sugarú 3D top-hat gömb $B(\mathbf{x}, L)$ **még belemetsz a csillag $[r - L, R_*]$ külső rétegébe** (a 18. lépésben levezetett sapkametszet-képlettel)! Ezért a csillag körül egy $L$ vastagságú **nemlokális vákuum-halo héj** ($r \in [R_*, R_* + L]$) alakul ki, amelyben:
  $$x(r) \;=\; L^2\bar{\rho}_L(r) \;=\; \frac{3}{4r\,L}\int_{r-L}^{R_*} s\,\varepsilon(s)\,\left[L^2 - (r - s)^2\right]\,ds \;>\; 0 \quad\Longrightarrow\quad f_L(r) \;=\; e^{-x(r)} \;<\; 1 \qquad (R_* \le r < R_* + L),$$
  és csak $r \ge R_* + L$-től kezdve válik $x(r) \equiv 0 \implies f_L(r) \equiv 1$-gyé!

### 4.2. A $\Delta\nu_{\text{ext}}$ véges tartományú integrál és a belső $\nu(r)$ profil levezetése (`tal_profiles.py`, 86–105. sor)
Integráljuk a `gradient` módú $\nu'(r) = f_L(r)\,A_{\text{GR}}(r)$ egyenletet az aszimptotikusan lapos végtelentől ($\nu(\infty) = 0$) a csillag $r = R_*$ felszínéig! Mivel $f_L(r) = 1 - (1 - e^{-x(r)})$, és $1 - e^{-x(r)} \equiv 0$ minden $r \ge R_* + L$-re, a végtelen integrál egzaktul szétválik a Schwarzschild-értékre és egy kompakt $[R_*, R_* + L]$ intervallumon vett véges korrekcióra (`tal_profiles.py`, 91–105. sor):
$$\boxed{\;\begin{aligned}
  \nu_{\text{grad}}(R_*) &\;=\; -\int_{R_*}^\infty e^{-x(r)}\,\frac{M}{r(r - 2M)}\,dr \;=\; -\int_{R_*}^\infty \frac{M}{r(r - 2M)}\,dr \;+\; \int_{R_*}^{R_* + L}\bigl(1 - e^{-x(r)}\bigr)\,\frac{M}{r(r - 2M)}\,dr \\[6pt]
  &\;=\; \frac{1}{2}\ln\!\left(1 - \frac{2M}{R_*}\right) \;+\; \Delta\nu_{\text{ext}}, \qquad \Delta\nu_{\text{ext}} \;\equiv\; \int_{R_*}^{R_* + L}\bigl(-\operatorname{expm1}(-x(r))\bigr)\,\frac{M}{r(r - 2M)}\,dr \;>\; 0,
\end{aligned}\;}$$
amelyet az `exterior_lapse()` rutin 200 pontos Gauss–Legendre kvadratúrával értékel ki az $[R_*, R_* + L]$ szakaszon.
Végül a csillag belsejében ($0 \le r \le R_*$) mind a `mass`, mind a `gradient` módban azonosan fennáll a $\nu'(r) = -h'(r)$ relativisztikus Bernoulli-összefüggés, amelyből a belső lapse-profil további numerikus kvadratúra nélkül, egzaktul (`tal_profiles.py`, 86–87. sor):
$$\boxed{\;\nu(r) \;=\; \nu(R_*) \;-\; \bigl[h(r) - h(R_*)\bigr], \qquad \forall\,r \in [0, R_*].\;}$$

---

## 5. Ellenőrzések (Checks)

1. **A $\Delta\nu_{\text{ext}} > 0$ előjelének és határeseteinek fizikai ellenőrzése:**
   - Mivel $x(r) > 0$ az $r \in [R_*, R_* + L)$ tartományban és $r > 2M$, az integráns $(1 - e^{-x(r)})\frac{M}{r(r-2M)} > 0$, így $\Delta\nu_{\text{ext}} > 0$, vagyis a `gradient` mód csillapítja a felszíni gravitációs potenciálmélységet ($|\nu_{\text{grad}}(R_*)| < |\nu_{\text{GR}}(R_*)|$). Amikor $L \to 0^+$, az integrálási tartomány hossza $L \to 0$ és $\Delta\nu_{\text{ext}} \to 0$, folytonosan visszaadva a Schwarzschild-határértéket $\checkmark$.
2. **A belső Bernoulli-integrál ($\nu(r) = \nu(R_*) - [h(r) - h(R_*)]$) deriváltjának ellenőrzése:**
   - Mindkét módban $\frac{d\nu}{dr} = -\frac{dh}{dr} = -2r\frac{dh_s}{ds}(r^2) = f_{\text{grad}}(r)A_{\text{GR}}(r)$, amely $r = 0$-ban egzaktul $\nu'(0) = 0$-t, $r = R_*$-nál pedig $\nu(R_*)$-ot ad $\checkmark$.

---

## 6. Explicit összehasonlítás a V23 lokális Yukawa (`yukawa_action_el`, Branch II-Fluid) ággal (Point 5)

A bírálati egyértelműség érdekében szigorúan el kell különíteni a fent levezetett V21 fenomenologikus 3D top-hat (`tophat_mass`, `tophat_gradient` — **Branch I-A**) tulajdonságait a V23 lokális Yukawa-árnyékolt Einstein–skalár–folyadék (`yukawa_action_el` — **Branch II-Fluid**, lásd `step4`, `step7` és `step22`) ágtól:

| Fizikai tulajdonság | V21 `tophat_mass` / `tophat_gradient` (**Branch I-A**, jelen 21. lépés) | V23 `yukawa_action_el` (**Branch II-Fluid**, 4., 7. és 22. lépés) |
|---|---|---|
| **Belső hidrosztatikai / Bernoulli-egyenlet ($0 \le r \le R_*$)** | $h'(r) = -\nu'(r) \implies \nu(r) = \nu(R_*) - [h(r) - h(R_*)]$ (az $f_L$ közvetlenül $\nu' = f_L A_{\text{GR}}$-be vagy $m' = 4\pi r^2 f_L\varepsilon$-ba lép be) | $h'(r) = -\nu'(r) + \phi_L'(r) \implies \nu(r) - \phi_L(r) = \nu(R_*) - \phi_L(R_*) - [h(r) - h(R_*)]$ (külön skalár ötödik erő $\phi_L'$) |
| **Külső tartomány ($r > R_*$) szerkezete** | Kompakt tartójú 3D top-hat halo $r \in [R_*, R_* + L]$-en; $r \ge R_* + L$-re $x(r) \equiv 0$, $f_L(r) \equiv 1$, egzakt Schwarzschild-vákuum | Végtelen tartományon ($r \in [R_*, \infty)$) exponenciálisan lecsengő, gravitáló skalárfarok: $\rho_\phi = \frac{1}{2\alpha_Y}[B(\phi_L')^2 + \phi_L^2/L^2] > 0$ |
| **Felszíni tömeg ($M_{\text{surf}}$) vs. ADM-tömeg ($M_{\text{ADM}}$)** | $m'(r > R_*) = 0 \implies M_{\text{ADM}} \equiv M_{\text{surf}} = m(R_*)$ | $m'(r > R_*) = 4\pi r^2\rho_\phi(r) > 0 \implies M_{\text{ADM}} = M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}} > M_{\text{surf}}$ (lásd `step22`) |

> **V27 frissítés (görbült `yukawa_action_el` ág):** 4-komponensű külső BVP ($\delta\nu_R=1{,}63\times 10^{-14}$), egzakt külső linearizáció, dinamikus skalárhullám, GL2 $N$-konvergencia ($p_{\text{eff}}=5{,}23$). SLy $L=10$ km, $\alpha_Y=1$: $M_{\text{ADM}}$-fordulópont $15{,}316404$, dinamikus $\omega_0^2=0$ gyök $15{,}316425$ ($2{,}04\times 10^{-5}$ dex). Részletes táblák: `step22`, 4. szakasz.
