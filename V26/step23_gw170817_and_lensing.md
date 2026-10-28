# 23. lépés (1. modul): GW170817 tenzori hullámsebesség ($c_T = 1$), tömeges skalár-diszperzió és gravitációs lencsézés ($\Phi + \Psi$, $\gamma_{\text{PPN}}(r)$, görbült téridőbeli fényelhajlás $\hat{\alpha}(b)$, Shapiro-időkésés és gyenge lencsézési konvergencia $\kappa$) $L = 10\;\text{km}$, $\alpha_Y = 1$ mellett

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Verzió / Version:** V26  
**Dátum / Date:** 2026-10-27  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step23_gw170817_and_lensing.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md` (V24 Zenodo DOI: `10.5281/zenodo.23680715`; V25 Zenodo DOI: `10.5281/zenodo.23699685`)  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források és futtatható kód / Local Sources and Runnable Code:**  
- Futtatható megoldókód: `scalar-tensor-yukawa-fit/gw170817_and_lensing_audit.py` (szinkronizálva: `code/gw170817_and_lensing_audit.py`)  
- Generált kimeneti fájlok: `gw170817_dispersion_and_cT.csv`, `lensing_deflection_and_ppn_slip.csv`, `gw170817_and_lensing_summary.json` (szinkronizálva: `csv/` és `data/v23/`)

---

## 1. Alapvető hatáselv és csatolási szótár ($L = 10\;\text{km}$, $\alpha_Y = 1$)

Az Idővonzástörvény (Time Attraction Law, TAL) kovariáns skalár–tenzor megfogalmazásában (`step1`–`step22`, V25 DOI: `10.5281/zenodo.23699685`) az Einstein-keretbeli gravitációs és tömeges Yukawa-skalár hatás:

$$S \;=\; \int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}\,R \;-\; \frac{c^4}{16\pi G\,\tilde{\alpha}_L}\left(g^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L + 2\,V(\phi_L)\right)\right] \;+\; S_m\!\left[\tilde{g}_{\mu\nu},\,\Psi_m\right],$$

ahol a tömeges Yukawa-árnyékolási potenciál:

$$V(\phi_L) \;=\; \frac{1}{2}\,m_\phi^2\,\phi_L^2, \qquad m_\phi \;=\; \frac{1}{L} \;=\; 0{,}1\;\text{km}^{-1} \;=\; 10^{-4}\;\text{m}^{-1},$$

és a csatolási állandók szótára (`step1`–`step5`, `step14`, `step22`):

$$\alpha_Y \;=\; 4\pi\,\alpha_{\text{std}} \;=\; 8\pi\,\tilde{\alpha}_L \;=\; 8\pi\,\beta^2, \qquad \gamma_L \;=\; \frac{G L^2}{c^2} \;=\; 7{,}42616027\times 10^{-20}\;\text{m}^3\,\text{kg}^{-1}.$$

**Egyetlen hatás, egyetlen mezőváltozó (V26 normálási javítás).** A 23–25. lépés minden állítása a fenti $S[\phi_L]$ hatásból következik, $\tilde{\alpha}_L = \alpha_Y/(8\pi)$ mellett, az $S_m[A^2(\phi_L)g_{\mu\nu},\Psi_m]$ anyagcsatolással, $A(\phi_L) = \exp(-\alpha_Y\phi_L)$, $d\ln A/d\phi_L = -\alpha_Y$ ($\phi_L > 0$ az anyag közelében, $g^{(J)}_{\mu\nu} = \exp(-2\alpha_Y\phi_L)g_{\mu\nu}$). Ha a $\frac{c^4}{16\pi G}[\dots]$ zárójelen belüli kinetikus együtthatót $1/\tilde{\alpha}_L$-ről $2$-re váltjuk (Damour–Esposito-Farèse alak), a mezőt **át kell skálázni**; ugyanaz a $\phi_L$ szimbólum változatlan $A = \exp(-\alpha_Y\phi_L)$-val nem tartható meg:

$$\frac{1}{\tilde{\alpha}_L}(\partial\phi_L)^2 \;=\; 2(\partial\phi_{\text{canon}})^2 \quad\Longleftrightarrow\quad \phi_{\text{canon}} \;=\; \frac{\phi_L}{\sqrt{2\tilde{\alpha}_L}} \;=\; \sqrt{\frac{4\pi}{\alpha_Y}}\,\phi_L, \qquad \phi_L \;=\; \sqrt{\frac{\alpha_Y}{4\pi}}\,\phi_{\text{canon}},$$

$$S_{\text{can}} \;=\; \int d^4x\,\sqrt{-g}\,\frac{c^4}{16\pi G}\Bigl[R - 2\,g^{\mu\nu}\partial_\mu\phi_{\text{canon}}\partial_\nu\phi_{\text{canon}} - 2\,m_\phi^2\phi_{\text{canon}}^2\Bigr] \;+\; S_m\!\left[A^2 g_{\mu\nu},\,\Psi_m\right], \qquad A \;=\; \exp\!\Bigl(-\alpha_Y\sqrt{2\tilde{\alpha}_L}\,\phi_{\text{canon}}\Bigr),$$

$$\frac{d\ln A}{d\phi_{\text{canon}}} \;=\; \sqrt{2\tilde{\alpha}_L}\,\frac{d\ln A}{d\phi_L} \;=\; -\alpha_Y\sqrt{\frac{\alpha_Y}{4\pi}} \;=\; -\frac{1}{\sqrt{4\pi}} \;=\; -0{,}2820947918 \quad (\alpha_Y = 1).$$

Ekvivalensen, a $4\pi$-normált $\hat{\phi}_{\text{canon}} \equiv \phi_{\text{canon}}/\sqrt{4\pi} = \phi_L/\sqrt{\alpha_Y}$ változóval (kinetikus tag: $-\frac{c^4}{2G}(\partial\hat{\phi}_{\text{canon}})^2$) $d\ln A/d\hat{\phi}_{\text{canon}} = -\alpha_Y^{3/2}$. Az ötödik erő erőssége tehát

$$\boxed{\;\beta_{\text{Yuk}} \;=\; \Bigl(\frac{d\ln A}{d\phi_{\text{canon}}}\Bigr)^{2} \;=\; \frac{(d\ln A/d\hat{\phi}_{\text{canon}})^2}{4\pi} \;=\; 2\tilde{\alpha}_L\Bigl(\frac{d\ln A}{d\phi_L}\Bigr)^{2} \;=\; 2\tilde{\alpha}_L\,\alpha_Y^2 \;=\; \frac{\alpha_Y^3}{4\pi} \;=\; \frac{1}{4\pi} \;\approx\; 0{,}0795774715 \quad (\alpha_Y = 1).\;}$$

Az $S[\phi_L]$ variálásából a nemrelativisztikus por- és gyengetér-határban ($\tilde{P} \ll \tilde{\varepsilon} \simeq \rho_m c^2$, $|\phi_L| \ll 1$) a linearizált sztatikus forrásegyenlet és a Jordan-keretbeli newtoni potenciál:

$$\bigl(\nabla^2 - L^{-2}\bigr)\phi_L \;=\; -8\pi\tilde{\alpha}_L\,\alpha_Y\,\frac{G\rho_m}{c^2} \;=\; -\alpha_Y^2\,\frac{G\rho_m}{c^2}, \qquad \tilde{\Psi}(r) \;=\; -\frac{GM}{r}\bigl(1 + \beta_{\text{Yuk}}e^{-r/L}\bigr), \qquad \frac{G_{\text{eff}}(k)}{G} \;=\; 1 + \beta_{\text{Yuk}}\frac{k^2L^2}{1 + k^2L^2}.$$

$\alpha_Y = 1$ mellett a $\beta_{\text{Yuk}} = 1/(4\pi)$ érték változatlan a V25-höz képest; $\alpha_Y \neq 1$ esetén a konform akcióból következő konzisztens skálázás $\beta_{\text{Yuk}} = 2\tilde{\alpha}_L\alpha_Y^2 = \alpha_Y^3/(4\pi)$ (`BETA_YUK = 2*ALPHA_TILDE_L*ALPHA_Y**2`), az effektív Brans–Dicke-paraméter pedig $\omega_{\text{BD}} = \frac{1}{2\beta_{\text{Yuk}}} - \frac{3}{2} = \frac{2\pi}{\alpha_Y^3} - \frac{3}{2} = 4{,}7831853072$ ($\alpha_Y = 1$ mellett; elkülönítve a `step5` fejezet korábbi $4\pi$-konvenciós $\frac{1 - 3\alpha_{\text{std}}}{2\alpha_{\text{std}}} = \frac{4\pi}{\alpha_Y} - \frac{3}{2} = 11{,}0663706144$ képletétől). Itt $\alpha_{\text{DEF}} \equiv d\ln A/d\phi_{\text{canon}} = -1/\sqrt{4\pi} = -0{,}2820947918$ az elsőrendű Damour–Esposito-Farèse anyagcsatolási derivált (második deriváltja $\beta_{\text{DEF}} \equiv d^2\ln A/d\phi_{\text{canon}}^2 = 0$), míg $\beta_{\text{Yuk}} \equiv \alpha_{\text{DEF}}^2 = 1/(4\pi)$ a dimenziótlan ötödikerő-arány.

### 1.1. A `Branch II-Conf` (23. lépés) és a `Branch II-Fluid` (22. lépés) egzakt összevetése: belső csillagegyenletek, keretleképezés és külső vákuum-azonosság

Míg a `Branch II-Fluid` (`step3`, `step4`, `step7`, `step12`, `step22`) és a `Branch II-Conf` (`step1`, `step2`, `step5`, `step14`, `step23`–`step25`) ágak pontosan ugyanazon az Einstein-keretbeli gravitációs és tömeges skalárhatáson ($S_{g+\phi}[g_{\mu\nu},\phi_L]$) osztoznak, egy adott $\tilde{\varepsilon}(h), \tilde{P}(h), \tilde{\rho}_0(h)$ ($dh = d\tilde{P}/(\tilde{\varepsilon} + \tilde{P})$) barotróp állapotegyenletű csillag belsejében két eltérő anyagcsatolást alkalmaznak:

1. **`Branch II-Fluid` (`yukawa_action_el` / a 22. lépés közvetlen folyadéksűrűség-csatolása):**
   - *Anyagi Lagrange-sűrűség:* $\mathcal{L}_m^{\text{(II-Fluid)}} = -e^{-\phi_L}\,\varepsilon(n)$, ahol $P(n) = n\,d\varepsilon/dn - \varepsilon(n)$, az Einstein-keretbeli effektív anyagi energia–impulzus tenzor pedig $T^{\mu\nu}_{\text{matter,(II-Fluid)}} = e^{-\phi_L}[(\varepsilon + P)u^\mu u^\nu + P g^{\mu\nu}]$, azaz $\rho_m^{(E)} = e^{-\phi_L}\varepsilon(h)$ és $P_m^{(E)} = e^{-\phi_L}P(h)$.
   - *Einstein- ($m', \nu'$), skalár- ($\phi_L''$) és Euler–TOV ($h'$) egyenletek ($B(r) = 1 - 2m(r)/r$, $G=c=1$):*
     $$\frac{dm}{dr} = 4\pi r^2\bigl(e^{-\phi_L}\varepsilon + \rho_\phi\bigr), \qquad \frac{d\nu}{dr} = \frac{m + 4\pi r^3\bigl(e^{-\phi_L}P + p_{r,\phi}\bigr)}{r(r - 2m)}, \qquad \frac{dh}{dr} = -\nu'(r) + \phi_L'(r),$$
     $$B\,\phi_L'' + \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L' - \frac{\phi_L}{L^2} \;=\; -\frac{\alpha_Y G}{c^2}\,e^{-\phi_L}\,\varepsilon(r), \qquad \rho_\phi = \frac{B(\phi_L')^2 + \phi_L^2/L^2}{2\alpha_Y}, \quad p_{r,\phi} = \frac{B(\phi_L')^2 - \phi_L^2/L^2}{2\alpha_Y}.$$

2. **`Branch II-Conf` (a 23. lépés univerzális konform Jordan-csatolása):**
   - *Anyagi hatás és keretleképezés:* $S_m^{\text{(II-Conf)}}[\tilde{g}_{\mu\nu},\Psi_m]$, ahol $\tilde{g}_{\mu\nu} = A^2(\phi_L)g_{\mu\nu} = e^{-2\alpha_Y\phi_L}g_{\mu\nu}$ és $\sqrt{-\tilde{g}} = A^4(\phi_L)\sqrt{-g} = e^{-4\alpha_Y\phi_L}\sqrt{-g}$. A Jordan-keretbeli lokális együttmozgó rendszerben mért $\tilde{\varepsilon}(h)$ energiasűrűséggel, $\tilde{P}(h)$ nyomással és $\tilde{u}^\mu = A^{-1}u^\mu$ négyessebességgel jellemzett folyadék Jordan-tenzora $\tilde{T}^{\mu\nu} = (\tilde{\varepsilon} + \tilde{P})\tilde{u}^\mu\tilde{u}^\nu + \tilde{P}\tilde{g}^{\mu\nu}$ (nyoma $\tilde{T} = \tilde{g}_{\mu\nu}\tilde{T}^{\mu\nu} = -\tilde{\varepsilon} + 3\tilde{P}$), Einstein-keretbeli energia–impulzus tenzora pedig $T^{\mu\nu}_{\text{matter,(II-Conf)}} = \frac{2}{\sqrt{-g}}\frac{\delta S_m}{\delta g_{\mu\nu}} = A^6(\phi_L)\tilde{T}^{\mu\nu}$, azaz $T^\mu{}_{\nu,\text{matter,(II-Conf)}} = e^{-4\alpha_Y\phi_L}\,\text{diag}(-\tilde{\varepsilon}, \tilde{P}, \tilde{P}, \tilde{P})$, így $\rho_m^{(E)} = e^{-4\alpha_Y\phi_L}\tilde{\varepsilon}(h)$, $P_m^{(E)} = e^{-4\alpha_Y\phi_L}\tilde{P}(h)$ és $T_{\text{matter}}^{(E)} = -e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$.
   - *Einstein- ($m', \nu'$), skalár- ($\phi_L''$) és Euler–TOV ($h'$) egyenletek:*
     $$\frac{dm}{dr} = 4\pi r^2\bigl(e^{-4\alpha_Y\phi_L}\tilde{\varepsilon} + \rho_\phi\bigr), \qquad \frac{d\nu}{dr} = \frac{m + 4\pi r^3\bigl(e^{-4\alpha_Y\phi_L}\tilde{P} + p_{r,\phi}\bigr)}{r(r - 2m)}, \qquad \frac{dh}{dr} = -\nu'(r) + \alpha_Y\phi_L'(r),$$
     $$B\,\phi_L'' + \left(\frac{2B}{r} + B\nu' + \frac{B'}{2}\right)\phi_L' - \frac{\phi_L}{L^2} \;=\; -\frac{\alpha_Y^2 G}{c^2}\,e^{-4\alpha_Y\phi_L}\bigl(\tilde{\varepsilon}(r) - 3\tilde{P}(r)\bigr).$$

**Szerkezeti következmények a 22–25. lépésre:**
- **Külső vákuum ($r > R_*$, ahol $\tilde{\varepsilon} = \tilde{P} = 0$):** A csillag felszíne felett mindkét ág **pontosan ugyanazt a 4-komponensű görbült külső Einstein–skalár vákuum ODE-rendszert** (`step22`) elégíti ki; eltérés kizárólag a belső integrálásból érkező $(R_*, M_{\text{surf}}, \phi_L(R_*))$ felszíni peremértékekben lép fel.
- **Gyengetérbeli nemrelativisztikus por-határ ($\tilde{P} \ll \tilde{\varepsilon} \simeq \rho_m c^2$, $|\phi_L| \ll 1$):** Az $\alpha_Y = 1$ alapparaméter mellett ($\alpha_Y^2 = \alpha_Y = 1$) a két linearizált forrás egybeesik; $\alpha_Y \neq 1$ esetén azonban a `Branch II-Conf` forrása $\alpha_Y^2$-tel, míg a `Branch II-Fluid` forrása $\alpha_Y$-nal skálázódik.
- **Relativisztikus neutroncsillag-belső ($0 \le r \le R_*$):** Mivel a relativisztikus nyomás csökkenti a nyomot (az SLy $\log_{10}\rho_{0c}=15{,}0$ csillag centrumában $3\tilde{P}_c/\tilde{\varepsilon}_c = 0{,}566065$, így $1 - 3\tilde{P}_c/\tilde{\varepsilon}_c = 0{,}433935$), és a konform térfogati súlyfaktor $e^{-4\alpha_Y\phi_L}$ az $e^{-\phi_L}$ helyett, a két belső BVP **nem azonos**. Ennek megfelelően az alábbi 3.3. szakaszban a V25 22. lépésének hátterén számított fénysugárkövetési táblázatot (**2. táblázat**, $M_{\text{ADM}} = 1{,}50999814\;M_\odot$) explicit módon **a `Branch II-Fluid` háttéren végzett feltételes ágak közötti diagnosztikaként** címkézzük, és közvetlenül mellé tesszük az önkonzisztens **`Branch II-Conf` nyom-forrásolt belső + külső BVP megoldást** (**2B. táblázat**, $M_{\text{ADM}}^{\text{(II-Conf)}} = 1{,}56203491\;M_\odot$).

A rögzített $L = 10\;\text{km}$, $\alpha_Y = 1$ alapparaméterek mellett a csatolási szótár számértékei:

| Mennyiség | Képlet | Számérték ($L = 10\;\text{km}$, $\alpha_Y = 1$) | Mértékegység |
|---|---|---|---|
| Yukawa-árnyékolási hossz ($L$) | $L = 1/m_\phi$ | `10.0` (`1.0e4`) | $\text{km}$ ($\text{m}$) |
| Yukawa-hullámszám ($m_\phi$) | $m_\phi = 1/L$ | `0.1` (`1.0e-4`) | $\text{km}^{-1}$ ($\text{m}^{-1}$) |
| Skalár Compton-tömeg ($m_\phi \hbar c$) | $\hbar c / L$ | `1.9732698e-11` | $\text{eV}$ |
| Skalár levágási körfrekvencia ($\omega_{\text{cut}}$) | $c / L$ | `29979.2458` | $\text{rad}\,\text{s}^{-1}$ |
| Skalár levágási frekvencia ($f_{\text{cut}}$) | $c / (2\pi L)$ | `4771.345159` | $\text{Hz}$ |
| Yukawa forráserősség ($\alpha_Y$) | $4\pi\alpha_{\text{std}} = 8\pi\tilde{\alpha}_L$ | `1.0` | dimenziótlan |
| Standard ötödik erő csatolás ($\alpha_{\text{std}}$) | $\alpha_Y / (4\pi)$ | `0.0795774715` | dimenziótlan |
| Hatás-normálási állandó ($\tilde{\alpha}_L$) | $\alpha_Y / (8\pi)$ | `0.0397887358` | dimenziótlan |
| Kanonikus átskálázás ($\phi_{\text{canon}}/\phi_L$) | $1/\sqrt{2\tilde{\alpha}_L} = \sqrt{4\pi/\alpha_Y}$ | `3.5449077018` | dimenziótlan |
| Elsőrendű DEF-derivált ($\lvert\alpha_{\text{DEF}}\rvert = \lvert d\ln A/d\phi_{\text{canon}}\rvert$; $\beta_{\text{DEF}} = d^2\ln A/d\phi_{\text{canon}}^2 = 0$) | $\alpha_Y\sqrt{2\tilde{\alpha}_L} = 1/\sqrt{4\pi}$ | `0.2820947918` | dimenziótlan |
| Ötödik erő erőssége ($\beta_{\text{Yuk}}$) | $(d\ln A/d\phi_{\text{canon}})^2 = (d\ln A/d\hat{\phi}_{\text{canon}})^2/(4\pi) = 2\tilde{\alpha}_L\alpha_Y^2$ | `0.0795774715` | dimenziótlan |
| Effektív Brans–Dicke paraméter ($\omega_{\text{BD}}$) | $1/(2\beta_{\text{Yuk}}) - 3/2 = 2\pi/\alpha_Y^3 - 3/2$ | `4.7831853072` | dimenziótlan |

---

## 2. A rész: GW170817 / GRB 170817A és a tenzori gravitációs hullámok sebessége ($c_T = 1$)

### 2.1. Kovariáns perturbációs levezetés az Einstein- és a Jordan-keretben (Horndeski-leképezés)

1. **Einstein-keretbeli tenzormódusok ($g_{\mu\nu} = \bar{g}_{\mu\nu} + h_{\mu\nu}$):**
   Az Einstein-keretben a gravitációs kinetikus tag tisztán az Einstein–Hilbert-tag $\frac{c^4}{16\pi G}R$; a skalármező $\phi_L$ nem csatolódik deriváltakkal a Riemann-, Ricci- vagy Einstein-tenzorhoz. A transzverzális-nyommentes (TT) mértékben ($h^\mu{}_\mu = 0$, $\bar{\nabla}^\mu h_{\mu\nu}^{\text{TT}} = 0$) a másodrendű kvadratikus tenzori hatás:

   $$S_T^{(2)} \;=\; \frac{c^4}{64\pi G}\int d^4x\,\sqrt{-\bar{g}}\left[-\bar{g}^{\mu\nu}\partial_\mu h_{ij}^{\text{TT}}\,\partial_\nu h^{ij\,\text{TT}}\right] \quad\Longrightarrow\quad \Box_{\bar{g}}\,h_{ij}^{\text{TT}} \;=\; 0.$$

   A tenzori gravitációs hullámok terjedési sebessége az Einstein-keretben tehát egzaktul:

   $$c_{T,\text{Einstein}}^2 \;=\; c^2 \quad\Longleftrightarrow\quad \frac{c_{T,\text{Einstein}} - c}{c} \;=\; 0.$$

2. **Jordan-keretbeli Horndeski-függvények ($\tilde{g}_{\mu\nu} = A^2(\phi_L) g_{\mu\nu}$):**
   Áttérve a fizikai Jordan-metrikára $\tilde{g}_{\mu\nu} = A^2(\phi_L) g_{\mu\nu}$ (a rögzített V26 konvencióban: $A(\phi_L) = \exp(-\alpha_Y\phi_L)$, $\phi_L > 0$), a skalár–tenzor hatás az általános Horndeski-alakba írható, ahol a görbületi csatolási függvények:

   $$\tilde{G}_4(\phi_L) \;=\; \frac{c^4}{16\pi G}\,A^{-2}(\phi_L), \qquad \tilde{G}_{4,X} \;\equiv\; \frac{\partial\tilde{G}_4}{\partial X} \;=\; 0, \qquad \tilde{G}_5(\phi_L, X) \;=\; 0,$$

   ahol $X \equiv -\frac{1}{2}\tilde{g}^{\mu\nu}\partial_\mu\phi_L\partial_\nu\phi_L$. Mivel $\tilde{G}_4$ kizárólag $\phi_L$ algebrai függvénye és független a kinetikus $X$ invariánstól ($\tilde{G}_{4,X} = 0$, $\tilde{G}_5 = 0$), az általános Horndeski tenzori terjedési sebesség-képlet alapján:

   $$\boxed{\;c_{T,\text{Jordan}}^2 \;=\; \frac{\tilde{G}_4}{\tilde{G}_4 - 2 X \tilde{G}_{4,X}}\,c^2 \;=\; c^2 \quad\Longrightarrow\quad \frac{c_T - c}{c} \;=\; 0 \quad \text{azonosan minden }(L,\,\alpha_Y)\text{ esetén!}\;}$$

   **Érvényességi kör (V26):** a $c_T = c$ eredmény a tisztán konform, minimális Einstein–Hilbert + kanonikus skalár + univerzális Maxwell-csatolású ágon ($G_{4,X} = 0$, $G_5 = 0$) bizonyított; ez nem bizonyítja automatikusan a $c_T = c$ egyenlőséget egy diszformális/refraktív ágon ($D(\phi_L, X) \ne 0$).

3. **Fénykúpok konform invarianciája:**
   A konform transzformáció $\tilde{g}_{\mu\nu} = A^2(\phi_L) g_{\mu\nu}$ azonosan megőrzi a nullkúpokat ($\tilde{g}_{\mu\nu} k^\mu k^\nu = A^2(\phi_L) g_{\mu\nu} k^\mu k^\nu = 0 \Longleftrightarrow g_{\mu\nu} k^\mu k^\nu = 0$). Mivel a négydimenziós vákuumbeli Maxwell-elektrodinamika hatása ($\sqrt{-\tilde{g}}\,\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} = \sqrt{-g}\,g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}$) konform-invariáns, az elektromágneses hullámok (fotonok) és a $h_{+,\times}^{\text{TT}}$ tenzori gravitációs hullámok ugyanazon a nullkúpon terjednek.
   - **Összevetés a GW170817 / GRB 170817A megfigyelési korláttal** ($\Delta t = 1{,}74 \pm 0{,}05\;\text{s}$ időkülönbség $D = 40\;\text{Mpc} = 1{,}2343\times 10^{20}\;\text{km}$ távolságon, $t_{\text{fény}} = 4{,}1171\times 10^{15}\;\text{s} \approx 1{,}3046\times 10^8\;\text{év}$, Abbott et al. 2017):
     $$-3\times 10^{-15} \;\le\; \frac{c_T - c}{c} \;\le\; +7\times 10^{-15}, \qquad \left(\frac{c_T - c}{c}\right)_{\text{TAL}} \;=\; 0.$$
     A szimbolikus (`SymPy`) és numerikus ellenőrzés (`gw170817_and_lensing_summary.json`) egyaránt `0` maradékot ad mindkét keretben.

---

### 2.2. A tömeges skalár „lélegző” (breathing) módus diszperziója és a $4{,}77\;\text{kHz}$-es levágás

A vákuumbeli linearizált skaláregyenlet ($\Box\phi_L - m_\phi^2\phi_L = 0$, $m_\phi = 1/L$) síkhullám-megoldásainak ($\delta\phi_L \propto e^{i(\mathbf{k}\cdot\mathbf{x} - \omega t)}$) diszperziós relációja:

$$\omega^2 \;=\; c^2 k^2 + \omega_{\text{cut}}^2, \qquad \omega_{\text{cut}} \;=\; \frac{c}{L} \;=\; 29\,979{,}2458\;\text{rad}\,\text{s}^{-1}, \qquad f_{\text{cut}} \;=\; \frac{c}{2\pi L} \;=\; 4771{,}345159\;\text{Hz}.$$

1. **Levágás alatti (evaneszcens / Yukawa-árnyékolt) tartomány ($f < f_{\text{cut}} = 4771{,}345\;\text{Hz}$):**
   Minden olyan frekvencián, amely kisebb a $4{,}771\;\text{kHz}$-es Compton-levágásnál — ideértve a kettőspulzárok pályafrekvenciáit ($f \sim 10^{-4}\;\text{Hz}$), a LISA sávot ($10^{-2}\;\text{Hz}$), az egész GW170817 LIGO bespirálozási sávot ($f_{\text{GW}} \in [24, 600]\;\text{Hz}$), sőt a neutroncsillag-összeolvadás utáni alap kvadrupól $f$-módust ($f \sim 1500\text{--}3000\;\text{Hz}$) is — a hullámszám tisztán képzetes:

   $$k \;=\; i\,\kappa_{\text{ev}}(f), \qquad \kappa_{\text{ev}}(f) \;=\; \frac{1}{L}\sqrt{1 - \left(\frac{f}{f_{\text{cut}}}\right)^2}, \qquad \lambda_{\text{ev}}(f) \;=\; \frac{1}{\kappa_{\text{ev}}(f)} \;=\; \frac{L}{\sqrt{1 - (f/f_{\text{cut}})^2}} \;\approx\; 10\;\text{km}.$$

   A skalár lélegző módus nem terjedő hullám, hanem a forrástól $\lambda_{\text{ev}}(f) \approx 10\text{--}12{,}86\;\text{km}$ karakterisztikus hosszon exponenciálisan lecsengő közeli tér. A GW170817 $D = 40\;\text{Mpc} = 1{,}2343\times 10^{20}\;\text{km}$ távolságán az amplitúdó-elnyomás:

   $$\exp\!\left(-\frac{D}{\lambda_{\text{ev}}(f)}\right) \;\sim\; 10^{-5{,}36\times 10^{19}} \;=\; 0.$$

2. **Levágás feletti (terjedő diszperzív) tartomány ($f > f_{\text{cut}} = 4771{,}345\;\text{Hz}$):**
   Ha egy összeolvadás utáni magasfrekvenciás folyamat $f > 4771{,}345\;\text{Hz}$ frekvenciát gerjesztene, annak csoportsebessége és fázissebessége:

   $$v_g(f) \;=\; \frac{\partial\omega}{\partial k} \;=\; c\sqrt{1 - \left(\frac{f_{\text{cut}}}{f}\right)^2} \;<\; c, \qquad v_p(f) \;=\; \frac{\omega}{k} \;=\; \frac{c^2}{v_g(f)} \;>\; c.$$

   $D = 40\;\text{Mpc}$ távolságon a terjedő skalármódus késése a $c_T = c$ tenzorhullámhoz képest:
   - $f = 5000\;\text{Hz}$-en ($v_g/c = 0{,}298949$): $\Delta t_{\text{scalar}} = 9{,}655\times 10^{15}\;\text{s} = 3{,}059\times 10^8\;\text{év}$ ($305{,}9$ millió év),
   - $f = 10\,000\;\text{Hz}$-en ($v_g/c = 0{,}878830$): $\Delta t_{\text{scalar}} = 5{,}676\times 10^{14}\;\text{s} = 1{,}799\times 10^7\;\text{év}$ ($18{,}0$ millió év).
   Így a detektorokhoz $D = 40\;\text{Mpc}$-en kizárólag a $c_T = c$ sebességű $h_+, h_\times$ tenzorpolarizációk jutnak el.

#### 1. táblázat: GW170817 tenzorhullám-sebesség és a tömeges skalármódus diszperziója ($L = 10\;\text{km}$, $\alpha_Y = 1$; `gw170817_dispersion_and_cT.csv`)

| Módus / Sáv (`mode_label`) | $f$ ($\text{Hz}$) | $f / f_{\text{cut}}$ | Skalár rezsim | $c_T / c$ (Einstein & Jordan) | $(c_T - c)/c$ | $k_{\text{real}}$ ($\text{km}^{-1}$) | $\kappa_{\text{ev}}$ ($\text{km}^{-1}$) | $\lambda_{\text{ev}}$ ($\text{km}$) | $\log_{10}\!\left(A_{\text{scalar}}(40\;\text{Mpc})\right)$ | $v_g / c$ | $\Delta t_{\text{scalar}}(40\;\text{Mpc})$ ($\text{év}$) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `PSR_B1913p16_orbital` | `3.583e-5` | `7.509e-9` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.100000` | `10.000000` | `-5.360e+19` | `0.0` | `0.0` |
| `PSR_J0737m3039_orbital` | `1.132e-4` | `2.372e-8` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.100000` | `10.000000` | `-5.360e+19` | `0.0` | `0.0` |
| `LISA_mHz_band` | `1.000e-2` | `2.096e-6` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.100000` | `10.000000` | `-5.360e+19` | `0.0` | `0.0` |
| `GW170817_LIGO_entry_24Hz` | `24.0` | `0.005030` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.099999` | `10.000127` | `-5.360e+19` | `0.0` | `0.0` |
| `GW170817_LIGO_100Hz` | `100.0` | `0.020958` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.099978` | `10.002197` | `-5.359e+19` | `0.0` | `0.0` |
| `GW170817_LIGO_300Hz` | `300.0` | `0.062875` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.099802` | `10.019825` | `-5.350e+19` | `0.0` | `0.0` |
| `GW170817_LIGO_600Hz` | `600.0` | `0.125751` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.099206` | `10.080016` | `-5.318e+19` | `0.0` | `0.0` |
| `NS_post_merger_1500Hz` | `1500.0` | `0.314377` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.094930` | `10.534097` | `-5.089e+19` | `0.0` | `0.0` |
| `NS_f_mode_2500Hz` | `2500.0` | `0.523961` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.085174` | `11.740643` | `-4.566e+19` | `0.0` | `0.0` |
| `NS_high_3000Hz` | `3000.0` | `0.628754` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.077760` | `12.860005` | `-4.168e+19` | `0.0` | `0.0` |
| `Sub_cutoff_4500Hz` | `4500.0` | `0.943130` | `evanescent_Yukawa_screened` | `1.0` | `0.0` | `0.0` | `0.033242` | `30.082119` | `-1.782e+19` | `0.0` | `0.0` |
| `Yukawa_cutoff_f_cut` | `4771.345159` | `1.000000` | `cutoff_threshold` | `1.0` | `0.0` | `0.0` | `0.000000` | `inf` | `0.0` | `0.0` | `inf` |
| `Super_cutoff_4800Hz` | `4800.0` | `1.006006` | `propagating_dispersive` | `1.0` | `0.0` | `0.010976` | `0.0` | `0.0` | `0.0` | `0.109105` | `1.065291e+09` |
| `Super_cutoff_5000Hz` | `5000.0` | `1.047923` | `propagating_dispersive` | `1.0` | `0.0` | `0.031328` | `0.0` | `0.0` | `0.0` | `0.298949` | `3.059412e+08` |
| `Super_cutoff_7500Hz` | `7500.0` | `1.571884` | `propagating_dispersive` | `1.0` | `0.0` | `0.121277` | `0.0` | `0.0` | `0.0` | `0.771541` | `3.863088e+07` |
| `Super_cutoff_10000Hz` | `10000.0` | `2.095845` | `propagating_dispersive` | `1.0` | `0.0` | `0.184189` | `0.0` | `0.0` | `0.0` | `0.878830` | `1.798767e+07` |

---

## 3. B rész: Gravitációs lencsézés ($\Phi + \Psi$, $\gamma_{\text{PPN}}(r)$, görbült téridőbeli $\hat{\alpha}(b)$ és konvergencia $\kappa$) $L = 10\;\text{km}$, $\alpha_Y = 1$ mellett

### 3.1. Gyenge téri potenciálok az Einstein- és a Jordan-keretben: konform kiejtődés a Weyl-potenciálban

Newtoni mértékben az Einstein-keretbeli gyenge téri metrika:

$$ds^2 \;=\; -\left(1 + \frac{2\Psi_E(r)}{c^2}\right)c^2 dt^2 \;+\; \left(1 - \frac{2\Phi_E(r)}{c^2}\right)d\mathbf{x}^2, \qquad \Psi_E(r) \;=\; \Phi_E(r) \;=\; -\frac{G M_{\text{ADM}}}{r}.$$

A fizikai Jordan-keretbeli metrikában ($\tilde{g}_{\mu\nu} = A^2(\phi_L) g_{\mu\nu} = \exp(-2\alpha_Y\phi_L) g_{\mu\nu}$, $A(\phi_L) = \exp(-\alpha_Y\phi_L)$, $\phi_L > 0$) az időbeli ($\tilde{\Psi}$) és térbeli ($\tilde{\Phi}$) metrikus potenciálok első rendben ellentétes előjellel tolódnak el:

$$\tilde{\Psi}(r) \;=\; \Psi_E(r) - \alpha_Y c^2 \phi_L(r), \qquad \tilde{\Phi}(r) \;=\; \Phi_E(r) + \alpha_Y c^2 \phi_L(r).$$

Következésképpen a **nullgeodetikusokat (gravitációs lencsézést és Shapiro-időkésést rögzített $M_{\text{ADM}}$ mellett) kormányzó Weyl-lencsézési potenciálban a konform faktor lineáris rendben egzaktul kiesik**:

$$\boxed{\;\Phi_{\text{lens}}(r) \;\equiv\; \frac{1}{2}\bigl(\tilde{\Psi}(r) + \tilde{\Phi}(r)\bigr) \;=\; \frac{1}{2}\bigl(\Psi_E(r) + \Phi_E(r)\bigr) \;=\; -\frac{G M_{\text{ADM}}}{r}\quad (\text{a Yukawa-energiafelhőn kívül}).\;}$$

---

### 3.2. Skálafüggő PPN-csúszás ($\gamma_{\text{PPN}}(r)$) és a lencsézési/dinamikai tömegarány

Míg a relativisztikus fotonok a $\Phi_{\text{lens}}(r) = \frac{1}{2}(\tilde{\Psi} + \tilde{\Phi}) = -G M_{\text{ADM}}/r$ kombinációra érzékenyek, addig a nemrelativisztikus próbatestek (bolygók, csillagok, gázfelhők) a dinamikai Newton-potenciált érzékelik:

$$\tilde{\Psi}(r) \;=\; -\frac{G M_{\text{ADM}}}{r}\bigl(1 + \beta_{\text{Yuk}}\,e^{-r/L}\bigr), \qquad \tilde{\Phi}(r) \;=\; -\frac{G M_{\text{ADM}}}{r}\bigl(1 - \beta_{\text{Yuk}}\,e^{-r/L}\bigr).$$

Ebből a skálafüggő Eddington–Robertson PPN-paraméter $\gamma_{\text{PPN}}(r) \equiv \tilde{\Phi}(r)/\tilde{\Psi}(r)$ és a lencsézési/dinamikai potenciálarány:

$$\boxed{\;\gamma_{\text{PPN}}(r) \;=\; \frac{1 - \beta_{\text{Yuk}}\,e^{-r/L}}{1 + \beta_{\text{Yuk}}\,e^{-r/L}}, \qquad \frac{\Phi_{\text{lens}}(r)}{\tilde{\Psi}(r)} \;=\; \frac{1 + \gamma_{\text{PPN}}(r)}{2} \;=\; \frac{1}{1 + \beta_{\text{Yuk}}\,e^{-r/L}}.\;}$$

A V26 csomag egyetlen normálást használ (`lensing_deflection_and_ppn_slip.csv`): $\beta_{\text{Yuk}}(\alpha_Y) = 2\tilde{\alpha}_L\alpha_Y^2 = \alpha_Y^3/(4\pi) = 1/(4\pi) \approx 0{,}0795774715$ ($\alpha_Y = 1$; $\tilde{\alpha}_L = \frac{1}{8\pi}$). Párhuzamos $\beta_{\text{Yuk}} = 2{,}0$ oszlop nem szerepel.

A Cassini-szonda 2003-as együttállási mérésénél ($b = 1{,}6\,R_\odot = 1{,}11312\times 10^6\;\text{km}$, $b/L = 1{,}11312\times 10^5$) az exponenciális Yukawa-faktor $e^{-b/L} = 10^{-48\,342{,}19} = 0$:

$$\gamma_{\text{PPN}}(1{,}6\,R_\odot) \;=\; 1, \qquad \bigl|1 - \gamma_{\text{PPN}}(1{,}6\,R_\odot)\bigr| \;<\; 10^{-48\,342} \;\ll\; 2{,}3\times 10^{-5}\;\text{(Bertotti et al. 2003)}.$$

---

### 3.3. Egzakt görbült téridőbeli elhajlási szög ($\hat{\alpha}(b)$) és Shapiro-időkésés: feltételes diagnosztika a V25 22. lépés (`Branch II-Fluid`) hátterén és összevetés az önkonzisztens `Branch II-Conf` nyom-forrásolt megoldással

Amint az 1.1. szakaszban levezettük, a V25 22. lépésének neutroncsillag-háttere (`curved_yukawa_adm_and_stability.py`) a belső tartományt a `Branch II-Fluid` $\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$ csatolással (skalárforrás: $-\alpha_Y e^{-\phi_L}\varepsilon$) integrálja, míg a `Branch II-Conf` ágon a relativisztikus csillag belsejében a skalármezőt a nyom ($-\alpha_Y^2 e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$) forrásozza. A V25 22. lépésbeli metrikáján végzett fénysugárkövetés ezért **a `Branch II-Fluid` háttéren kiértékelt feltételes ágak közötti diagnosztika**, amelyet a **2. táblázatban** mutatunk be, közvetlenül kiegészítve az **önkonzisztens `Branch II-Conf` nyom-forrásolt belső + külső BVP megoldással** a **2B. táblázatban** (`gw170817_and_lensing_summary.json`, `branch_ii_fluid_vs_branch_ii_conf_interior_bvp_comparison` kulcs).

A V25 (`step22`, DOI: `10.5281/zenodo.23699685`, `Branch II-Fluid` háttér) referencia SLy neutroncsillagára ($\log_{10}(\rho_{0c}/[\text{g cm}^{-3}]) = 15{,}0$, $L = 10\;\text{km}$, $\alpha_Y = 1$) az egzakt 4-komponensű külső BVP paraméterei:
- Felszíni sugár: $R_* = 11{,}46572423\;\text{km}$,
- Felszíni koordináta-tömeg: $M_{\text{surf}} = m(R_*) = 1{,}50766122\;M_\odot$ ($2{,}22625027\;\text{km}$),
- Aszimptotikusan lapos ADM-tömeg: $M_{\text{ADM}} = \lim_{r\to\infty} m(r) = 1{,}50999814\;M_\odot$ ($2{,}22970113\;\text{km}$),
- Külső skalárfelhő tömegjáruléka: $\Delta M_{\phi,\text{ext}}^{\text{curved}} = +0{,}00233692\;M_\odot$ ($+0{,}155003\%$ a felszíni tömeghez képest),
- Felszíni skalármező és görbült derivált: $\phi_L(R_*) = 0{,}0054086823$, $\phi_{L,\text{curved}}'(R_*) = -0{,}0012631563\;\text{km}^{-1}$ ($k_{\text{curved}} = 0{,}23354233\;\text{km}^{-1}$),
- Felszíni metrikus potenciál: $\nu_{\text{exact}}(R_*) = -0{,}2463066872$ (szemben a vákuum Schwarzschild-értékekkel: $\nu_{\text{Schw}}(M_{\text{surf}}) = -0{,}2457823129$, $\nu_{\text{Schw}}(M_{\text{ADM}}) = -0{,}2462745914$).

A gömbszimmetrikus külső téridőben ($ds^2 = -e^{2\nu(r)}c^2 dt^2 + B(r)^{-1}dr^2 + r^2 d\Omega^2$, $B(r) = 1 - 2m(r)/r$) a $b = r_0\,e^{-\nu_{\text{eff}}(r_0)}$ impaktparaméterű fénysugár egzakt görbült téridőbeli elhajlási szöge ($u = r_0/r = \sin\theta$, $\theta \in [0, \pi/2]$ helyettesítéssel, kioltásmentes 128-pontos Gauss–Legendre kvadratúrával):

$$\hat{\alpha}(b) \;=\; 2\int_0^{\pi/2}\left[\frac{\cos\theta}{\sqrt{B(r_0/\sin\theta)}\sqrt{\exp\!\bigl(2[\nu_{\text{eff}}(r_0) - \nu_{\text{eff}}(r_0/\sin\theta)]\bigr) - \sin^2\theta}} \;-\; 1\right]d\theta,$$

ahol a konform Jordan-ágon $\nu_{\text{eff}}(r) = \nu(r)$ (mivel $d\tilde{s}^2 = e^{-2\alpha_Y\phi_L}ds^2 = 0 \Longleftrightarrow ds^2 = 0$), míg a diszformális refraktív ágon (`step2`, `step14`, $n_{\text{disf}}(r) = e^{\phi_L(r)}$) $\nu_{\text{eff}}(r) = \nu(r) - \phi_L(r)$.

**Három fizikai következmény a csillagfelszíntől ($b = 14{,}749\;\text{km}$) a galaxishalmaz-skáláig ($b = 1\;\text{Mpc}$):**
1. **A külső Yukawa-energiafelhő Birkhoff-sértő hatása a neutroncsillag közelében ($b \in [14{,}75, 100]\;\text{km}$):**
   - A csillagfelszínt súroló fénysugárnál ($b = 14{,}749032\;\text{km}$, legközelebbi pont $r_0 = 11{,}558475\;\text{km}$) a bezárt tömeg $m(r_0) = 1{,}50772597\;M_\odot < M_{\text{ADM}}$.
   - A konform Yukawa-megoldás elhajlási szöge $\hat{\alpha}_{\text{Yukawa,conf}} = 1{,}29085207\;\text{rad}$ ($73{,}96037^\circ$), ami **$+0{,}36521\%$-kal nagyobb**, mint a felszíni tömeggel számolt vákuum Schwarzschild-érték ($\hat{\alpha}_{\text{Schw}}(M_{\text{surf}}) = 1{,}28615489\;\text{rad}$), viszont **$-0{,}02867\%$-kal kisebb**, mint az aszimptotikus ADM-tömeggel számolt vákuum Schwarzschild-érték ($\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}}) = 1{,}29122231\;\text{rad}$), mivel a fénysugár a külső skalár-energiafelhő belsejében halad át!
   - A diszformális refraktív ágon ($n_{\text{disf}} = e^{\phi_L}$) az elhajlási szög a felszínnél $\hat{\alpha}_{\text{Yukawa,disf}} = 1{,}34618158\;\text{rad}$ ($+4{,}28628\%$-kal nagyobb a konform ágnál), de $b = 100\;\text{km}$-nél az eltérés már csak $+5{,}253\times 10^{-6}$, $b \ge 1000\;\text{km}$-nél pedig $0$.
2. **Aszimptotikus ADM-tömeg-kalibráció ($b \ge 200\;\text{km}$):**
   - $b \ge 200\;\text{km} = 20\,L$ felett a skalármező lecseng ($\phi_L(200\;\text{km}) = 1{,}069\times 10^{-12}$), így $\hat{\alpha}_{\text{Yukawa,conf}}(b)$ gépi pontossággal ($\le 10^{-12}$ relatív eltéréssel, $b \ge R_\odot$-nál $0{,}0$-val) megegyezik az $\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}}, b)$ értékkel, miközben minden távolságon állandó **$+0{,}155003\% = \Delta M_{\phi,\text{ext}}^{\text{curved}}/M_{\text{surf}}$** többletet mutat a naiv $\hat{\alpha}_{\text{Schw}}(M_{\text{surf}}, b)$ görbéhez képest.

#### 2. táblázat: Skálafüggő PPN-csúszás $\gamma_{\text{PPN}}(b)$, görbült téridőbeli fényelhajlás $\hat{\alpha}(b)$ és Shapiro-időkésés $b = 14{,}75\;\text{km}$-től $1\;\text{Mpc}$-ig (`lensing_deflection_and_ppn_slip.csv`)

| Skála (`scale_label`) | $b$ ($\text{km}$) | $r_0$ ($\text{km}$) | $b / L$ | $\log_{10}(e^{-b/L})$ | $m(r_0)$ ($M_\odot$) | $\gamma_{\text{PPN}}^{\text{TAL}}(b)$ | $\lvert 1 - \gamma_{\text{PPN}}^{\text{TAL}}\rvert$ | $\hat{\alpha}_{\text{Yuk,conf}}$ ($\text{rad}$) | $\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})$ ($\text{rad}$) | $\hat{\alpha}_{\text{Schw}}(M_{\text{surf}})$ ($\text{rad}$) | $\frac{\hat{\alpha}_{\text{Yuk,conf}} - \hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})}{\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})}$ | $\frac{\hat{\alpha}_{\text{Yuk,conf}} - \hat{\alpha}_{\text{Schw}}(M_{\text{surf}})}{\hat{\alpha}_{\text{Schw}}(M_{\text{surf}})}$ | $\Delta t_{\text{Shapiro}}^{\text{Yuk,conf}}$ ($\mu\text{s}$) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `NS_surface_grazing` | `14.749032` | `11.558475` | `1.4749` | `-0.6405` | `1.507726` | `0.964236` | `3.576e-2` | `1.29085207` | `1.29122231` | `1.28615489` | `-2.867e-4` | `+3.652e-3` | `420.139048` |
| `NS_near_halo_15km` | `15.000000` | `11.843795` | `1.5000` | `-0.6514` | `1.507913` | `0.965107` | `3.489e-2` | `1.23768742` | `1.23800029` | `1.23332361` | `-2.527e-4` | `+3.538e-3` | `417.501587` |
| `NS_Yukawa_cloud_20km` | `20.000000` | `17.215854` | `2.0000` | `-0.8686` | `1.509517` | `0.978690` | `2.131e-2` | `0.70230838` | `0.70233272` | `0.70052511` | `-3.465e-5` | `+2.546e-3` | `387.157454` |
| `NS_Yukawa_cloud_30km` | `30.000000` | `27.455794` | `3.0000` | `-1.3029` | `1.509957` | `0.992107` | `7.893e-3` | `0.38749878` | `0.38749944` | `0.38670428` | `-1.700e-6` | `+2.055e-3` | `361.964526` |
| `NS_Yukawa_cloud_50km` | `50.000000` | `47.600307` | `5.0000` | `-2.1715` | `1.509998` | `0.998928` | `1.072e-3` | `0.20641763` | `0.20641763` | `0.20604614` | `-1.039e-8` | `+1.803e-3` | `339.066547` |
| `NS_Yukawa_tail_100km` | `100.000000` | `97.690938` | `10.000` | `-4.3429` | `1.509998` | `0.999993` | `7.226e-6` | `0.09556441` | `0.09556441` | `0.09540578` | `-2.099e-12` | `+1.663e-3` | `313.665881` |
| `NS_Yukawa_tail_200km` | `200.000000` | `197.731861` | `20.000` | `-8.6859` | `1.509998` | `1.000000` | `3.280e-10` | `0.04612015` | `0.04612015` | `0.04604631` | `+8.027e-13` | `+1.603e-3` | `290.894079` |
| `NS_far_1000km` | `1000.00000` | `997.762797` | `100.00` | `-43.429` | `1.509998` | `1.000000` | `5.921e-45` | `0.00897785` | `0.00897785` | `0.00896386` | `+2.771e-12` | `+1.560e-3` | `241.379351` |
| `Sirius_B_radius_5844km` | `5844.00000` | `5841.76902` | `584.40` | `-253.80` | `1.509998` | `1.000000` | `2.513e-255` | `0.00152786` | `0.00152786` | `0.00152550` | `+1.062e-12` | `+1.552e-3` | `188.460119` |
| `NS_far_field_b_6p957e5_km` | `6.95700e+5` | `6.95698e+5` | `6.957e+4` | `-3.021e+4` | `1.509998` | `1.000000` | `0.0` | `1.282002e-5` | `1.282002e-5` | `1.280018e-5` | `0.0` | `+1.550e-3` | `102.503515` |
| `NS_far_field_b_1p113e6_km` | `1.11312e+6` | `1.11312e+6` | `1.113e+5` | `-4.834e+4` | `1.509998` | `1.000000` | `0.0` | `8.012485e-6` | `8.012485e-6` | `8.000084e-6` | `0.0` | `+1.550e-3` | `102.503351` |
| `NS_far_field_b_1AU` | `1.49598e+8` | `1.49598e+8` | `1.496e+7` | `-6.497e+6` | `1.509998` | `1.000000` | `0.0` | `5.961853e-8` | `5.961853e-8` | `5.952626e-8` | `0.0` | `+1.550e-3` | `102.503080` |
| `Outer_Solar_10_AU` | `1.49598e+9` | `1.49598e+9` | `1.496e+8` | `-6.497e+7` | `1.509998` | `1.000000` | `0.0` | `5.961852e-9` | `5.961852e-9` | `5.952626e-9` | `0.0` | `+1.550e-3` | `102.503078` |
| `Subparsec_1_pc` | `3.08568e+13` | `3.08568e+13` | `3.086e+12` | `-1.340e+12` | `1.509998` | `1.000000` | `0.0` | `2.890388e-13` | `2.890388e-13` | `2.885914e-13` | `0.0` | `+1.550e-3` | `102.503078` |
| `Galactic_1_kpc` | `3.08568e+16` | `3.08568e+16` | `3.086e+15` | `-1.340e+15` | `1.509998` | `1.000000` | `0.0` | `2.890388e-16` | `2.890388e-16` | `2.885914e-16` | `0.0` | `+1.550e-3` | `102.503078` |
| `NGC3198_30p08_kpc` | `9.28172e+17` | `9.28172e+17` | `9.282e+16` | `-4.031e+16` | `1.509998` | `1.000000` | `0.0` | `9.609001e-18` | `9.609001e-18` | `9.594130e-18` | `0.0` | `+1.550e-3` | `102.503078` |
| `Bullet_Cluster_250_kpc` | `7.71419e+18` | `7.71419e+18` | `7.714e+17` | `-3.350e+17` | `1.509998` | `1.000000` | `0.0` | `1.156155e-18` | `1.156155e-18` | `1.154366e-18` | `0.0` | `+1.550e-3` | `102.503078` |
| `Cluster_virial_1_Mpc` | `3.08568e+19` | `3.08568e+19` | `3.086e+18` | `-1.340e+18` | `1.509998` | `1.000000` | `0.0` | `2.890388e-19` | `2.890388e-19` | `2.885914e-19` | `0.0` | `+1.550e-3` | `102.503078` |

**Megjegyzés (V26):** a 2. táblázat minden sora a 22. lépésbeli $1{,}50999814\;M_\odot$ tömegű `Branch II-Fluid` neutroncsillag-háttéren van kiértékelve (`source_body = NS_SLy_1p51Msun`); az `NS_far_field_b_*` sorok kizárólag távoltéri referencia-impaktparaméterek, nem napadatok. Az alábbi 2B. táblázat ezt a `Branch II-Fluid` diagnosztikai hátteret közvetlenül összeveti az ugyanazon SLy $\log_{10}\rho_{0c} = 15{,}0$, $L = 10\;\text{km}$, $\alpha_Y = 1$ csillagra kiszámított önkonzisztens `Branch II-Conf` nyom-forrásolt ($\tilde{\varepsilon} - 3\tilde{P}$) BVP megoldással, a 3. táblázat pedig a dedikált $M = 1{,}000000\;M_\odot$ Nap–Cassini benchmarkot adja meg.

#### 2B. táblázat: A 22. lépésbeli `Branch II-Fluid` háttér és a 23. lépésbeli önkonzisztens `Branch II-Conf` nyom-forrásolt belső + görbült külső BVP közvetlen összevetése (SLy, $\log_{10}\rho_{0c} = 15{,}0$, $L = 10\;\text{km}$, $\alpha_Y = 1$)

| Mennyiség / Megfigyelhető | `Branch II-Fluid` (22. lépés háttér-diagnosztika) | `Branch II-Conf` (Önkonzisztens nyom-forrásolt BVP) | Az eltérés fizikai eredete |
|---|---|---|---|
| Belső anyagi hatás | $\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$ | $S_m[e^{-2\alpha_Y\phi_L}g_{\mu\nu},\Psi_m]$ ($\sqrt{-\tilde{g}} = e^{-4\alpha_Y\phi_L}\sqrt{-g}$) | Közvetlen energiasűrűség-súlyozás vs. univerzális Jordan-metrika |
| Belső skalárforrás jobb oldala | $-(\alpha_Y G/c^2)\,e^{-\phi_L}\varepsilon(r)$ | $-(\alpha_Y^2 G/c^2)\,e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$ | Relativisztikus nyomás levonása a nyomban ($1 - 3\tilde{P}_c/\tilde{\varepsilon}_c = 0{,}433935$) |
| Belső Euler–TOV egyensúly | $h'(r) = -\nu'(r) + \phi_L'(r)$ | $h'(r) = -\nu'(r) + \alpha_Y\phi_L'(r)$ | $\alpha_Y = 1$-nél azonos alak, a saját $\phi_L(r)$ profilon kiértékelve |
| Külső ($r > R_*$) vákuum ODE-k | 4-komponensű görbült Einstein–skalár BVP | Azonos 4-komponensű görbült Einstein–skalár BVP | A vákuumegyenletek azonosan egybeesnek ($\tilde{\varepsilon} = \tilde{P} = 0$) |
| Csillagsugár $R_*$ ($\text{km}$) | `11.46572423` | `11.66269842` | Az $e^{-4\phi_L}$ súlyfaktor és a nyomforrás eltolja az egyensúlyi sugarat |
| Felszíni tömeg $M_{\text{surf}} = m(R_*)$ ($M_\odot$) | `1.50766122` | `1.56063550` | Belső tömegintegrál az $e^{-4\alpha_Y\phi_L}\tilde{\varepsilon}$ sűrűséggel |
| Aszimptotikus ADM-tömeg $M_{\text{ADM}}$ ($M_\odot$) | `1.50999814` | `1.56203491` | $M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$ |
| Külső skalárfelhő-tömeg $\Delta M_{\phi,\text{ext}}^{\text{curved}}$ ($M_\odot$) | `+0.00233692` (`+0.155003%`) | `+0.00139941` (`+0.089669%`) | Kisebb felszíni skalármező a $-3\tilde{P}$ nyom-csökkentés miatt |
| Centrális / felszíni skalármező $\phi_L(0)$ / $\phi_L(R_*)$ | `0.01777896` / `0.00540868` | `0.01157688` / `0.00414318` | A magban $1 - 3\tilde{P}/\tilde{\varepsilon} < 1$ |
| Görbült felszíni meredekség $\phi_{L,\text{curved}}'(R_*)$ ($\text{km}^{-1}$) / $k_{\text{curved}}$ ($\text{km}^{-1}$) | `-0.00126316` / `0.23354233` | `-0.00096508` / `0.23293146` | Görbült Dirichlet–Neumann külső illesztés |
| Súroló impaktparaméter $b_{\text{graz}}$ ($\text{km}$) / $r_0$ ($\text{km}$) | `14.749032` / `11.558475` | `15.064914` / `11.736304` | $1{,}0001\,R_* e^{-(\nu(R_*) - \phi_L(R_*))}$-nél kiértékelve |
| Súroló elhajlási szög $\hat{\alpha}_{\text{Yuk,conf}}(b_{\text{graz}})$ ($\text{rad}$) | `1.29085207` | `1.33402365` | Konform null-geodetikus a saját felszínt súroló sugáron |
| Azonos impaktparaméterű elhajlás $\hat{\alpha}_{\text{Yuk,conf}}(b = 15{,}064914\;\text{km})$ ($\text{rad}$) | `1.22474419` ($r_0 = 11{,}917144\;\text{km}$) | `1.33402365` ($r_0 = 11{,}736304\;\text{km}$) | Közös külső $b = b_{\text{graz}}^{\text{(II-Conf)}}$-nál kiértékelve ($+8{,}92\%$ a nagyobb $M_{\text{ADM}}$ miatt) |
| Súroló relatív eltérés vs. $\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})$ / $\hat{\alpha}_{\text{Schw}}(M_{\text{surf}})$ | `-2.867e-4` / `+3.652e-3` | `-1.671e-4` / `+2.184e-3` | A fénysugár belemerül a külső skalárfelhőbe ($m(r_0) < M_{\text{ADM}}$) |
| Távoltéri relatív eltérés vs. $\hat{\alpha}_{\text{Schw}}(M_{\text{ADM}})$ ($b \ge 200\;\text{km}$) | `<= 8.03e-13` (`0.0`, ha $b \ge R_\odot$) | `<= 1.79e-13` (`0.0`, ha $b \ge R_\odot$) | Exponenciális leárnyékolás $\exp(-b/L) \to 0$ mindkét ágon |
| Belső–külső felszíni meredekség-eltérés $\lvert\phi'_{L,\text{int}}(R_*) - \phi'_{L,\text{ext}}(R_*)\rvert$ ($\text{km}^{-1}$) | `1.88e-09` (4-lépéses Robin) | Robin: `1.465e-09` (4-lép.) / `3.937e-15` (minden lép.); FD2 ($N=4000$): `4.723e-10` | Kikényszerített Robin-perem vs. független egyoldali FD2 ($\mathcal{O}(h^2)$) |
| Konvergált belső skalár ODE-maradék $\max_{r\in(0,R_*)}\lvert\text{res}_{\text{ODE}}\rvert$ ($\text{km}^{-2}$) | `1.42e-10` (iterált) / `1.95e-14` (sávos) | `1.078e-10` (iter., `tol=2e-8`) / `6.958e-14` (`tol=1e-11`) / `1.480e-14` (sávos) | Diszkrét 2. rendű sávos Helmholtz–Yukawa-rendszer maradéka |

#### 2B.1. A `Branch II-Conf` SLy $\log_{10}(\rho_{0c}) = 15{,}0$ BVP célzott numerikus ellenőrzése: belső egyenletmaradék, felszíni $\phi'_L(R_*)$ illesztés és $N = 2000$ vs. $N = 4000$ hálófinomítás

A `gw170817_and_lensing_audit.py` (`compute_branch_ii_conf_bvp_convergence_and_mesh_refinement` függvény, `gw170817_and_lensing_summary.json` $\to$ `bvp_convergence_and_mesh_refinement` kulcs) közvetlenül naplózza az egyetlen vizsgált $L = 10\;\text{km}$, $\alpha_Y = 1$, SLy $\log_{10}(\rho_{0c}) = 15{,}0$ `Branch II-Conf` nyom-forrásolt referencia-csillag BVP-jének három numerikus ellenőrzési mutatóját:

1. **Konvergált belső diszkrét egyenletmaradék ($\max_{r\in(0,R_*)}\lvert\text{res}_{\text{ODE}}\rvert$, $\text{km}^{-2}$) és külső normalizált kollokációs maradék:**
   - A konvergált $N = 601$ rácson a belső skaláregyenlet $B\phi_L'' + (2B/r + B\nu' + B'/2)\phi_L' - \phi_L/L^2 + \alpha_Y^2 e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P}) = 0$ diszkrét maradéka a közvetlen tridiagonális sávos megoldón (`_solve_curved_scalar_with_k`, amely a saját diszkrét lineáris rendszerének megoldási pontosságát méri az adott rácson) **$\max\lvert\text{res}_{\text{ODE,banded}}\rvert = 1{,}480\times 10^{-14}\;\text{km}^{-2}$**, a relaxált fixpont-iterált profilon (`tol_iter = 2e-8`, 19 iteráció) **$\max\lvert\text{res}_{\text{ODE,iter}}\rvert = 1{,}078\times 10^{-10}\;\text{km}^{-2}$** ($\max\lvert\phi_L - \phi_{L,\text{banded}}\rvert = 3{,}065\times 10^{-9}$), míg a külön `tol_iter = 1e-11` küszöbű, minden lépésben frissített 4-komponensű külső Robin-peremes futáson ($N = 601$, 36 iteráció) **$\max\lvert\text{res}_{\text{ODE,iter}}\rvert = 6{,}958\times 10^{-14}\;\text{km}^{-2}$**. A külső 4-komponensű BVP (`scipy.integrate.solve_bvp`, `tol = 1e-12`) normalizált RMS kollokációs maradéka a külső hálón ($[R_*, R_{\max}]$) **$\max(\text{rms\_residuals}) = 9{,}998\times 10^{-13}$** (dimenziótlan normalizált kollokációs maradék, amely nem azonos a belső $\text{km}^{-2}$ egységű diszkrét ODE-maradékkal).
2. **A kikényszerített Robin-derivált ($\phi'_{L,\text{int,Robin}}(R_*)$) és a külön egyoldali véges differenciás derivált ($\phi'_{L,\text{int,FD2}}(R_*)$) szétválasztása $r = R_*$-nál:**
   - **Kikényszerített Robin-peremderivált ($\phi'_{L,\text{int,Robin}}(R_*) = -k_{\text{int}}\phi_L(R_*)$):** A 4-lépéses Robin-frissítési ütemezésnél (`it in (1, 3, 6, 10)`, `tol_iter = 2e-8`, $N = 601$) a belső Robin-tényező $k_{\text{int}} = 0{,}232931108388\;\text{km}^{-1}$, a külső 4-komponensű görbült Dirichlet–Neumann-tényező $k_{\text{curved}} = 0{,}232931461902\;\text{km}^{-1}$ ($\lvert\Delta k\rvert = 3{,}535\times 10^{-7}\;\text{km}^{-1}$, relatív eltérés $1{,}518\times 10^{-6}$); a kikényszerített Robin-derivált $\phi'_{L,\text{int,Robin}}(R_*) = -0{,}000965074633\;\text{km}^{-1}$ és a külső $\phi'_{L,\text{ext,curved}}(R_*) = -0{,}000965076098\;\text{km}^{-1}$ eltérése **$\lvert\Delta\phi'_{L,\text{Robin}}(R_*)\rvert = 1{,}465\times 10^{-9}\;\text{km}^{-1}$**. Ha a görbült Robin-tényezőt **minden külső iterációban** (`it >= 1`) a 4-komponensű külső BVP-ből frissítjük (`tol_iter = 1e-11`, $N = 601$, 36 iteráció), akkor $k_{\text{int}} = 0{,}232931462480\;\text{km}^{-1}$ és $k_{\text{curved}} = 0{,}232931462481\;\text{km}^{-1}$ (**$\lvert\Delta k\rvert = 9{,}502\times 10^{-13}\;\text{km}^{-1}$**), így a **kikényszerített Robin-derivált** egyezése **$\lvert\phi'_{L,\text{int,Robin}}(R_*) - \phi'_{L,\text{ext,curved}}(R_*)\rvert = 3{,}937\times 10^{-15}\;\text{km}^{-1}$** (a szellempontos centrális sémából visszaszámolt belső deriváltra pedig $\lvert\phi'_{L,\text{int,ghost}}(R_*) - \phi'_{L,\text{ext,curved}}(R_*)\rvert = 1{,}430\times 10^{-14}\;\text{km}^{-1}$), miközben $R_*$ és $M_{\text{ADM}}$ csupán $\lvert\Delta R_*\rvert = 2{,}06\times 10^{-7}\;\text{km}$-rel és $\lvert\Delta M_{\text{ADM}}\rvert = 6{,}07\times 10^{-8}\;M_\odot$-tal mozdul el.
   - **Független, 3-pontos egyoldali másodrendű véges differenciás derivált ($\phi'_{L,\text{int,FD2}}(R_*) = (3\phi_{N-1} - 4\phi_{N-2} + \phi_{N-3})/(2h)$):** A peremfeltételtől elkülönített egyoldali FD2 derivált eltérése a külső $\phi'_{L,\text{ext,curved}}(R_*)$-tól $N = 601$-nél $2{,}1213\times 10^{-8}\;\text{km}^{-1}$, $N = 2000$-nél $1{,}9048\times 10^{-9}\;\text{km}^{-1}$, **$N = 4000$-nél pedig $4{,}7228\times 10^{-10}\;\text{km}^{-1}$** ($\approx 4{,}89\times 10^{-7}$ relatív eltérés a $\lvert\phi'_L(R_*)\rvert = 9{,}6507\times 10^{-4}\;\text{km}^{-1}$ meredekséghez mérve), és a háló sűrítésével másodrendben ($\mathcal{O}(h^2)$) konvergál ($N = 2000 \to 4000$ empirikus javulási arány: `4.0332`, összhangban a másodrendű hálóskálázás elméleti $(h_{2000}/h_{4000})^2 = 4{,}0020$ értékével).
3. **Hálófinomítási próba ($N = 601$ vs. $N = 2000$ vs. $N = 4000$, `run_variant`: `DOP853`, `rtol_ivp = 1e-11`, `atol = 1e-13`, `max_step = 0.02 km`, `tol_iter = 1e-10`, minden lépéses Robin-frissítés):**

| Rács $N$ | $h = \Delta r$ ($\text{km}$) | $R_*$ ($\text{km}$) | $M_{\text{surf}}$ ($M_\odot$) | $M_{\text{ADM}}$ ($M_\odot$) | $\phi_L(0)$ | $\phi_L(R_*)$ | $\phi'_{L,\text{int,Robin}}(R_*)$ ($\text{km}^{-1}$) | $\phi'_{L,\text{ext,curved}}(R_*)$ ($\text{km}^{-1}$) | $\lvert\phi'_{L,\text{int,FD2}} - \phi'_{L,\text{ext}}\rvert$ ($\text{km}^{-1}$) | $\max\lvert\text{res}_{\text{ODE,iter}}\rvert$ ($\text{km}^{-2}$) | $\hat{\alpha}_{\text{Yuk,conf}}(b_{\text{graz}})$ ($\text{rad}$) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `601` | `0.019438` | `11.662698212` | `1.560635444` | `1.562034848` | `0.011576871` | `0.004143168` | `-0.000965074210` | `-0.000965074204` | `2.1213e-08` | `7.568e-13` | `1.334023631` |
| `2000` | `0.005834` | `11.662697490` | `1.560635270` | `1.562034651` | `0.011576835` | `0.004143134` | `-0.000965066327` | `-0.000965066321` | `1.9048e-09` | `7.844e-13` | `1.334023603` |
| `4000` | `0.002916` | `11.662697448` | `1.560635260` | `1.562034640` | `0.011576833` | `0.004143132` | `-0.000965065881` | `-0.000965065874` | `4.7228e-10` | `7.430e-13` | `1.334023601` |
| **$\lvert\Delta\rvert_{2000\to 4000}$** | — | **`4.213e-08`** (`3.61e-9`) | **`9.670e-09`** (`6.20e-9`) | **`1.097e-08`** (`7.03e-9`) | **`2.042e-09`** (`1.76e-7`) | **`1.920e-09`** (`4.63e-7`) | **`4.463e-10`** | **`4.463e-10`** | **Skálázás: `4.0332`** (várt $\mathcal{O}(h^2)$: `4.0020`) | — | **`1.055e-09`** (`7.91e-10`) |

#### 3. táblázat: Dedikált $1{,}000000\;M_\odot$ Nap–Cassini benchmark ($GM_\odot/c^2 = 1{,}4766250385\;\text{km}$, $R_\odot = 695700{,}0\;\text{km}$, $r_E = 1\;\text{AU}$, $r_S = 8{,}43\;\text{AU}$; `source_body = Sun_1p00Msun`)

| Sor (`scale_label`) | $b$ ($\text{km}$) | $b/R_\odot$ | $\hat{\alpha}_{\text{exact}}$ ($''$) | $\hat{\alpha}_{\text{1PN}}$ ($''$) | $\Delta t_{\text{Shapiro}}(10\,r_0)$ ($\mu\text{s}$) | $\Delta t_{\text{Cassini,rt}}$ ($\mu\text{s}$) | $\Delta t_{\text{Cassini,rt}}$ $-1$ geometriai taggal ($\mu\text{s}$) | $\log_{10}\lvert\gamma_{\text{PPN}}-1\rvert$ |
|---|---|---|---|---|---|---|---|---|
| `Solar_1Msun_limb_1_Rsun` | `6.957000e+05` | `1.0000` | `1.751201` | `1.751190` | `67.8831` | `280.9434` | `261.2415` | `-30214.67` |
| `Solar_1Msun_Cassini_1p6_Rsun` | `1.113120e+06` | `1.6000` | `1.094498` | `1.094494` | `67.8830` | `262.4234` | `242.7215` | `-48342.99` |
| `Solar_1Msun_10_Rsun` | `6.957000e+06` | `10.0000` | `0.175119` | `0.175119` | `67.8829` | — | — | `-302139.47` |
| `Solar_1Msun_Earth_1_AU` | `1.495979e+08` | `215.0322` | `0.008144` | `0.008144` | `67.8829` | — | — | `-6496953.77` |

---

### 3.4. Gyenge lencsézési konvergencia ($\kappa(\mathbf{\theta})$) galaktikus és halmazskálákon: átvezetés a 2. modulhoz (Bullet Cluster)

A gyenge gravitációs lencsézés konvergenciája ($\kappa(\mathbf{\theta}) = \frac{1}{2}\nabla_{\mathbf{\theta}}^2\psi(\mathbf{\theta})$) a látóirány mentén integrált Weyl-potenciálból ($\Phi_{\text{lens}} = \frac{1}{2}(\tilde{\Psi} + \tilde{\Phi})$) származik:

$$\kappa(\mathbf{b}_\perp) \;=\; \frac{1}{2\,c^2}\,\frac{D_L D_{LS}}{D_S}\int_{-\infty}^{+\infty}\nabla_\perp^2\bigl(\tilde{\Psi}(\mathbf{b}_\perp, z) + \tilde{\Phi}(\mathbf{b}_\perp, z)\bigr)\,dz \;=\; \frac{\Sigma_{\text{lens}}(\mathbf{b}_\perp)}{\Sigma_{\text{crit}}}, \qquad \Sigma_{\text{crit}} \;=\; \frac{c^2}{4\pi G}\,\frac{D_S}{D_L D_{LS}}.$$

Vizsgáljuk meg egzaktul, mit ad a rögzített $L = 10\;\text{km}$, $\alpha_Y = 1$ paraméterpár galaktikus ($\sim 30\;\text{kpc}$, SPARC NGC 3198) és galaxishalmaz ($\sim 250\;\text{kpc}$, Bullet Cluster 1E 0657-558) skálán (`gw170817_and_lensing_summary.json`):

1. **A lokális, rögzített $L = 10\;\text{km}$-es Yukawa-komponens viselkedése kpc–Mpc skálán:**
   - **Pontforrások (csillagok, kompakt objektumok) külső tere:** $b = 1\;\text{kpc}$ ($b/L = 3{,}086\times 10^{15}$) és $b = 250\;\text{kpc}$ ($b/L = 7{,}714\times 10^{17}$) távolságon az $e^{-b/L} = 10^{-3{,}35\times 10^{17}} = 0$ faktor miatt minden lokális $L = 10\;\text{km}$-es Yukawa-potenciál maradéktalanul leárnyékolódik.
   - **Diffúz barionos közegben (galaktikus korong vagy röntgengáz, $\rho_{\text{ICM}} \sim 10^{-22}\;\text{kg/m}^3$) ébredő lokális egyensúlyi mező:** Homogén háttérben a lokális Helmholtz–Yukawa-egyenlet megoldása $\phi_L^{\text{local}}(L = 10\;\text{km}) = \alpha_Y \gamma_L \rho$:
     - Az NGC 3198 külső korongjában ($\rho_{\text{gal}} \sim 10^{-24}\;\text{kg/m}^3$, $R = 30{,}08\;\text{kpc}$):
       $$\phi_L^{\text{local}}(10\;\text{km}) \;=\; 7{,}426\times 10^{-44}, \qquad \frac{|\Phi_N|}{c^2} \;\approx\; \frac{V_{\text{bar}}^2}{c^2} \;=\; 5{,}499\times 10^{-8} \quad\Longrightarrow\quad \frac{c^2\phi_L^{\text{local}}(10\;\text{km})}{|\Phi_N|} \;=\; 1{,}350\times 10^{-36}.$$
     - A Bullet Cluster forró intracluster plazmájában ($\rho_{\text{ICM}} \sim 10^{-22}\;\text{kg/m}^3$, $b = 250\;\text{kpc}$):
       $$\phi_L^{\text{local}}(10\;\text{km}) \;=\; 7{,}426\times 10^{-42} \quad\Longrightarrow\quad \kappa_{\text{lens}}^{\text{local }L=10\text{ km}}(\mathbf{b}_\perp) \;=\; \kappa_{\text{baryon}}(\mathbf{b}_\perp) \quad \left(\frac{\kappa_{\text{lens}}^{\text{local }L=10\text{ km}}}{\kappa_{\text{baryon}}} = 1{,}000000000000\right).$$

2. **Fizikai következtetés és peremfeltétel a 2. modulhoz (Bullet Cluster 1E 0657-558):**
   - Az $L = 10\;\text{km}$, $\alpha_Y = 1$ lokális Yukawa-tag szigorúan **erős téri / kompakt objektum ($r \lesssim 100\;\text{km}$) ötödik erő**, amely:
     - egzaktul kielégíti a GW170817 $c_T = c$ korlátot ($(c_T - c)/c = 0$),
     - a $4771{,}35\;\text{Hz}$-es levágás révén elnyomja a skalár dipól- és monopólsugárzást a kettőspulzárokban és a LIGO sávban,
     - a Naprendszerben ($b = 1{,}6\,R_\odot$) $\lvert 1 - \gamma_{\text{PPN}}\rvert = 0 < 2{,}3\times 10^{-5}$ pontossággal kielégíti a Cassini-korlátot,
     - ugyanakkor kiloparszekes és megaparszekes skálán ($b/L \sim 10^{16}\text{--}10^{18}$) **nem ad járulékot** sem a galaktikus rotációs görbéhez (ahol az NGC 3198-ban $R = 30{,}08\;\text{kpc}$-nél $C_{\text{req}} = (V_{\text{obs}}/V_{\text{bar}})^2 = 4{,}311$ szükséges $\Upsilon_* = 0{,}5$ mellett; `step6`), sem a Bullet Cluster lencsézési konvergencia-térképéhez (ahol a galaxiskomponens körül $\kappa_{\text{obs}}/\kappa_{\text{baryon,gal}} \approx 0{,}232 / 0{,}0258 \approx 9{,}0$ erősítés látható a $\Sigma_{\text{crit}} \approx 3{,}10\times 10^9\;M_\odot\,\text{kpc}^{-2}$ kritikus felületi sűrűségnél).
   - Ahhoz, hogy egy makroszkopikus óramező ($\phi_{\text{macro}}$, pl. a `step6`, `step11`, `step14` kaméleon/nemlokális $L_{\text{eff}}(\rho) \sim 10\text{--}300\;\text{kpc}$ ága) módosítsa a gyenge lencsézési konvergenciát ($\kappa > \kappa_{\text{baryon}}$), **két feltételnek kell egyszerre fennállnia**:
     1. A csatolás nem lehet tisztán konform ($\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu}$), mert a konform faktor kiesik a $(\tilde{\Psi} + \tilde{\Phi})/2$ Weyl-potenciálból, hanem **diszformális / refraktív csatolásra** ($\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu} + B(\phi)\partial_\mu\phi\partial_\nu\phi$ vagy effektív optikai metrikára $\Phi_{\text{lens}}^{\text{eff}} = \Phi_N - \frac{1}{2}c^2\phi_{\text{macro}}$, `step14`) van szükség; és
     2. A skalármező effektív hatótávolságának sűrűségfüggőnek ($L_{\text{eff}}(\rho) \propto \rho^{-1/2}$) vagy nemlokálisnak kell lennie, mivel a rögzített $L = 10\;\text{km}$-es Yukawa-mag $b \gg 10\;\text{km}$-en azonosan eltűnik.
