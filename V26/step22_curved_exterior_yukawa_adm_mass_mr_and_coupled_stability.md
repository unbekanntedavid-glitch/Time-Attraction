# 22. lépés (V25): A `yukawa_action_el` görbült külső vákuum-peremértékfeladata, az egzakt ADM-tömeg ($M_{\text{ADM}}$), az újraszámolt tömeg–sugár ($M_{\text{ADM}}$–$R_*$) görbék és a csatolt $(\delta h, \delta\phi_L, \delta m)$ radiális pulzációs stabilitás

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-27  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források és futtatható kód / Local Sources and Runnable Code:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.2. és 2. szakasz, 65–153. sor; DOI: `10.5281/zenodo.23639084`)  
- Futtatható megoldókód: `scalar-tensor-yukawa-fit/curved_yukawa_adm_and_stability.py` (és integrálva a `covariant_closure_50aa2db0.py` szkriptbe)  
- Generált kimeneti fájlok: `curved_yukawa_adm_sequences.csv`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_audit_summary.json`, `closure_tov_comparison.csv`, `conservation_Residuals.json`

---

## 1. Feladatkitűzés és a korábbi közelítés feloldása

A V23 (`50aa2db0`) eredeti `yukawa_action_el` számításában (`covariant_closure_50aa2db0.py`) két egyszerűsítő közelítés szerepelt:
1. A csillag $r = R_*$ felszínén a lapos téridőbeli Helmholtz–Yukawa-megoldás ($\phi_L^{\text{flat}}(r) = \phi_L(R_*)\frac{R_*}{r}e^{-(r-R_*)/L}$) deriváltját írta elő peremfeltételként:
   $$\phi_{L,\text{flat}}'(R_*) \;=\; -\left(\frac{1}{R_*} + \frac{1}{L}\right)\phi_L(R_*).$$
2. A tömeg–sugár sorozatok maximumát ($M_{\max}$), a fordulóponti centrális sűrűséget ($\log_{10}\rho_{0c}$) és az $1{,}4\,M_\odot$-hoz tartozó sugarat ($R_{1{,}4}$) a csillag felszínén kiértékelt $M_{\text{surf}} \equiv m(R_*)$ tömegre optimalizálta, míg a külső skalármező energiájából származó $\Delta M_{\phi,\text{ext}}$ járulékot külön oszlopban tüntette fel.

Ebben a 22. lépésben **egzaktul megoldjuk ugyanennek a `yukawa_action_el` ágnak a csatolt görbült külső Einstein–skalár vákuumrendszerét** ($r \in [R_*, \infty)$), meghatározzuk belőle a valódi görbült Dirichlet–Neumann-meredekséget $\phi_{L,\text{curved}}'(R_*)$ és az aszimptotikus **ADM-tömeget** ($M_{\text{ADM}} = \lim_{r\to\infty} m(r)$), majd az ADM-tömeggel (és a konzervált $M_B$ barionos nyugalmi tömeggel) újraszámoljuk az összes tömeg–sugár görbét (SLy, APR4, H4; $L = 1\;\text{km}$ és $L = 10\;\text{km}$; $\alpha_Y = 1$ és $\alpha_Y = 1/3$) és a csatolt $(\delta h, \delta\phi_L, \delta m)$ dinamikus radiális pulzációs stabilitást.

---

## 2. A görbült külső vákuum ($r \ge R_*$) egzakt Einstein–skalár egyenletei

### 2.1. A külső vákuum-azonosság és a másodrendű skaláregyenlet
A csillag felszínén kívül ($r \ge R_*$, ahol $\varepsilon(r) = P(r) = 0$) a `yukawa_action_el` (**Branch II-Fluid**, $G = c = 1$) hatásból származó skalár-energiasűrűség és radiális skalárnyomás:
$$\rho_\phi(r) \;=\; \frac{1}{2\alpha_Y}\left[B(r)\bigl(\phi_L'(r)\bigr)^2 + \frac{\phi_L(r)^2}{L^2}\right], \qquad p_{r,\phi}(r) \;=\; \frac{1}{2\alpha_Y}\left[B(r)\bigl(\phi_L'(r)\bigr)^2 - \frac{\phi_L(r)^2}{L^2}\right], \qquad B(r) \;=\; 1 - \frac{2m(r)}{r}.$$
A $(tt)$ és $(rr)$ Einstein-egyenletek a külső vákuumban:
$$\frac{dm(r)}{dr} \;=\; 4\pi r^2\rho_\phi(r), \qquad \frac{d\nu(r)}{dr} \;=\; \frac{m(r) + 4\pi r^3 p_{r,\phi}(r)}{r\bigl(r - 2m(r)\bigr)}.$$
A görbült téridőbeli Klein–Gordon / Yukawa-egyenletben ($B\phi_L'' + (\frac{2B}{r} + B\nu' + \frac{B'}{2})\phi_L' - \frac{\phi_L}{L^2} = 0$) szereplő metrikus kombinációra a külső vákuumban ($\varepsilon = P = 0$) fennáll az alábbi **egzakt algebrai azonosság**:
$$B(r)\nu'(r) + \frac{B'(r)}{2} \;=\; \frac{m + 4\pi r^3 p_{r,\phi}}{r^2} - \frac{m'}{r} + \frac{m}{r^2} \;=\; \frac{2m(r)}{r^2} + 4\pi r\bigl(p_{r,\phi}(r) - \rho_\phi(r)\bigr) \;=\; \frac{2m(r)}{r^2} - \frac{4\pi r}{\alpha_Y L^2}\phi_L(r)^2.$$
Ezt behelyettesítve a külső görbült vákuum csatolt $(m(r), \phi_L(r))$ egyenletrendszere zárt, szingularitásmentes alakot ölt minden $r \ge R_*$-ra:
$$\boxed{\;\begin{aligned}
  \left(1 - \frac{2m(r)}{r}\right)\phi_L''(r) \;+\; \left[\frac{2}{r}\left(1 - \frac{m(r)}{r}\right) - \frac{4\pi r}{\alpha_Y L^2}\phi_L(r)^2\right]\phi_L'(r) \;-\; \frac{\phi_L(r)}{L^2} &\;=\; 0, \\[6pt]
  \frac{dm(r)}{dr} &\;=\; \frac{2\pi r^2}{\alpha_Y}\left[\left(1 - \frac{2m(r)}{r}\right)\bigl(\phi_L'(r)\bigr)^2 + \frac{\phi_L(r)^2}{L^2}\right].
\end{aligned}\;}$$

### 2.2. Aszimptotikus WKB-lecsengés, $R_{\max}$-függetlenség és a görbült Dirichlet–Neumann-meredekség
Nagy távolságokon ($r \gg R_*, L$, ahol $m(r) \to M_{\text{ADM}}$) a fenti egyenlet aszimptotikus WKB-megoldása:
$$\phi_L(r) \;\underset{r\to\infty}{\sim}\; C\,r^{-\bigl(1 + M_{\text{ADM}}/L\bigr)} e^{-r/L} \quad\Longrightarrow\quad \phi_L'(R_{\max}) \;+\; \left(\frac{1}{L} + \frac{1 + m(R_{\max})/L}{R_{\max}}\right)\phi_L(R_{\max}) \;=\; 0.$$
Adott felszíni állapot $(R_*, M_{\text{surf}} = m(R_*), \phi_L(R_*))$ esetén a $[R_*, R_{\max}]$ ($R_{\max} = R_* + N_{\text{mult}} L$) intervallumon megoldott peremérték-feladat egyértelműen meghatározza a lecsengő ághoz tartozó valódi görbült felszíni deriváltat $\phi_{L,\text{curved}}'(R_*) = -k_{\text{curved}}\phi_L(R_*)$ és az aszimptotikus ADM-tömeget:
$$\boxed{\;M_{\text{ADM}} \;=\; \lim_{r\to\infty} m(r) \;=\; M_{\text{surf}} \;+\; \Delta M_{\phi,\text{ext}}^{\text{curved}}, \qquad \Delta M_{\phi,\text{ext}}^{\text{curved}} \;=\; \int_{R_*}^\infty 4\pi r^2\rho_\phi(r)\,dr.\;}$$

**Numerikus $R_{\max}$-konvergencia és a lapos vs. görbült Robin-meredekség összehasonlítása** (SLy, $\log_{10}\rho_{0c} = 15{,}0$, $L = 10\;\text{km}$, $\alpha_Y = 1$; `curved_yukawa_audit_summary.json`):
- **$R_{\max}$-függetlenség $N_{\text{mult}} \in \{10, 12, 14, 16\}$ között:**
  - A görbült Dirichlet–Neumann-meredekség relatív szórása: $\frac{\max k_{\text{curved}} - \min k_{\text{curved}}}{k_{\text{curved}}} = \mathbf{6{,}17\times 10^{-12}}$,
  - Az ADM-tömeg relatív szórása: $\frac{\max M_{\text{ADM}} - \min M_{\text{ADM}}}{M_{\text{ADM}}} = \mathbf{5{,}54\times 10^{-13}}$.
- **Görbült vs. lapos felszíni meredekség ($\log_{10}\rho_{0c} = 15{,}0$, $R_* = 11{,}4657\;\text{km}$, $B(R_*) = 0{,}6119$):**
  - Valódi görbült meredekség: $k_{\text{curved}} = -\phi_{L,\text{curved}}'(R_*)/\phi_L(R_*) = \mathbf{0{,}233542\;\text{km}^{-1}}$ ($\phi_{L,\text{curved}}'(R_*) = -1{,}263156\times 10^{-3}\;\text{km}^{-1}$),
  - Lapos Robin-közelítés: $k_{\text{flat}} = \frac{1}{R_*} + \frac{1}{L} = \mathbf{0{,}187216\;\text{km}^{-1}}$ ($\phi_{L,\text{flat}}'(R_*) = -1{,}012594\times 10^{-3}\;\text{km}^{-1}$),
  - Relatív eltérés a felszíni deriváltban: **$\mathbf{+24{,}74\%}$** (a téridő-görbület $B(R_*)^{-1/2} \approx 1{,}278$ faktorral meredekebbé teszi a külső Yukawa-lecsengést!).
  - Ennek következtében a valódi görbült külső skalártömeg $\Delta M_{\phi,\text{ext}}^{\text{curved}} = \mathbf{0{,}002337\,M_\odot}$ $\log_{10}\rho_{0c} = 15{,}0$-nál (míg a lapos becslés $0{,}00340\,M_\odot$ volt), és $\Delta M_{\phi,\text{ext}}^{\text{curved}} = \mathbf{0{,}004108\,M_\odot}$ az ADM-tömegmaximumnál (míg a lapos becslés $0{,}00898\,M_\odot$ volt).
- **Felszíni lapse ($\nu(R_*)$) az $r = R_* = 11{,}46572423\;\text{km}$ ($P(R_*)=h(R_*)=0$) határpontban ($\log_{10}\rho_{0c}=15{,}0$, $M_{\text{surf}}=1{,}50766122\,M_\odot$, $M_{\text{ADM}}=1{,}50999814\,M_\odot$):**
  - $\nu_{\text{Schw}}(M_{\text{surf}}) = \tfrac{1}{2}\ln(1 - 2M_{\text{surf}}/R_*) = -0{,}2457823130$,
  - $\nu_{\text{Schw}}(M_{\text{ADM}}) = \tfrac{1}{2}\ln(1 - 2M_{\text{ADM}}/R_*) = -0{,}2462745915$,
  - $\nu_{\text{exact}}(R_*) \equiv \nu_R(15{,}0) = \mathbf{-0{,}2463066872}$ (a 4-komponensű külső BVP $y=[\phi_L,\phi_L',m,\nu]$ megoldásából $\nu(R_{\max}) = \tfrac{1}{2}\ln(1 - 2m(R_{\max})/R_{\max})$ külső peremfeltétellel; egyezés a visszafelé integrált kvadratúrával $\delta\nu_R = 1{,}63\times 10^{-14}$, a külső skalárfarok eltolása $\nu_{\text{Schw}}(M_{\text{surf}})$-hoz képest $-5{,}2437\times 10^{-4}$, $\nu_{\text{Schw}}(M_{\text{ADM}})$-hoz képest $-3{,}2096\times 10^{-5}$).

---

## 3. A tömeg–sugár görbék újraszámolása az ADM-tömeggel ($M_{\text{ADM}}$) és a barionos nyugalmi tömeggel ($M_B$)

Az alábbi táblázat (`curved_yukawa_adm_sequences.csv` és `closure_tov_comparison.csv`) mind a 8 vizsgált `yukawa_action_el` sorozatra megadja az **ADM-tömegre ($M_{\text{ADM}}$) optimalizált** értékeket, együtt a felszíni tömeggel ($M_{\text{surf}}$) és az $s=r^2$ PCHIP + adaptív Gauss–Kronrod kvadratúrával (`epsrel=1e-12`) definiált konzervált barionos nyugalmi tömeg ($M_B^{\text{quad}} = \int_0^{R_*} 4\pi r^2 B(r)^{-1/2}\rho_0(r)\,dr$) maximumával és fordulópontjával:

| EOS | $L$ (km) | $\alpha_Y$ | $M_{\text{ADM},\max}$ ($M_\odot$) | $M_{\text{surf},\max}$ ($M_\odot$) | $M_{B,\max}^{\text{quad}}$ ($M_\odot$) | $\Delta M_{\phi,\text{ext}}^{\text{curved}}(M_{\max})$ ($M_\odot$) | $\log_{10}\rho_{0c}(M_{\text{ADM},\max})$ | $\log_{10}\rho_{0c}(M_{B,\max}^{\text{quad}})$ | $R(M_{\text{ADM},\max})$ (km) | $R(M_{\text{ADM}}=1{,}4\,M_\odot)$ (km) | $R(M_{\text{surf}}=1{,}4\,M_\odot)$ (km) | $M_{\text{ADM}}(15{,}0)$ ($M_\odot$) | $f_L(0; 15{,}0)$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **SLy** | `1.0` | `1.000` ($1/4\pi$) | `2.04684` | `2.04684` | `2.42634` | `0.000002` | `15.303086` | `15.303082` | `9.9689` | `11.6933` | `11.6933` | `1.56171` | `0.999198` |
| **SLy** | `1.0` | `0.333` ($1/12\pi$) | `2.04789` | `2.04789` | `2.42744` | `0.000001` | `15.302656` | `15.302660` | `9.9733` | `11.7018` | `11.7018` | `1.56386` | `0.999733` |
| **SLy** | `10.0` | `1.000` ($1/4\pi$) | **`2.03321`** | `2.02911` | **`2.42339`** | `0.004108` | **`15.316404`** | **`15.316409`** | `9.8853` | **`11.5359`** | `11.5348` | `1.51000` | `0.982378` |
| **SLy** | `10.0` | `0.333` ($1/12\pi$) | `2.04321` | `2.04183` | `2.42630` | `0.001381` | `15.307166` | `15.307169` | `9.9446` | `11.6475` | `11.6471` | `1.54638` | `0.993967` |
| **APR4** | `1.0` | `1.000` ($1/4\pi$) | `2.18693` | `2.18693` | `2.65730` | `0.000003` | `15.282808` | `15.282805` | `9.8947` | `11.3085` | `11.3085` | `1.56316` | `0.999201` |
| **APR4** | `10.0` | `1.000` ($1/4\pi$) | `2.18068` | `2.17624` | `2.66733` | `0.004449` | `15.293022` | `15.293025` | `9.8462` | `11.1830` | `11.1828` | `1.51271` | `0.982098` |
| **H4** | `1.0` | `1.000` ($1/4\pi$) | `2.01139` | `2.01138` | `2.29535` | `0.000001` | `15.209751` | `15.209739` | `11.6917` | `13.9388` | `13.9388` | `1.88120` | `0.999166` |
| **H4** | `10.0` | `1.000` ($1/4\pi$) | `1.97809` | `1.97502` | `2.26306` | `0.003077` | `15.227785` | `15.227787` | `11.5289` | `13.7179` | `13.7169` | `1.82871` | `0.981076` |

**Bariontömeg-definíció és fordulóponti pontosság (`baryonic_mass_and_first_law_audit`):** A `curved_yukawa_adm_sequences.csv` és `closure_tov_comparison.csv` táblákban szereplő `M_baryon_max_Msun` és `log10_rho0c_at_M_baryon_max` oszlopok definíció szerint az $s=r^2$ változóban PCHIP-interpolált belső profilokon végzett adaptív Gauss–Kronrod kvadratúrából ($M_B^{\text{quad}}$, `epsrel=1e-12`) és $\pm 0{,}002\;\text{dex}$-es 7 pontos kvadratikus illesztésből származnak. SLy $L=10\;\text{km}$, $\alpha_Y=1$ esetén $M_{B,\max}^{\text{quad}} = \mathbf{2{,}42338503\,M_\odot}$ ($\log_{10}\rho_{0c}(M_{B,\max}^{\text{quad}}) = 15{,}31640854$, amely $N_{\text{grid}} = 301\to 2401$ között $1{,}2\times 10^{-6}\;\text{dex}$-en belül konvergál), a finomított trapézszabály $N_{\text{grid}}=2401$-nél $M_B^{\text{trapz}} = 2{,}42338667\,M_\odot$ értéket és $15{,}31641814$-es 7 pontos csúcsot ad ($+9{,}6\times 10^{-6}\;\text{dex}$), míg az 5. ODE-állapotváltozóként egyszeri menetben integrált $M_B^{\text{ODE}}$ a kéregszéli ($h\to 0$) $\sqrt{s}\,\rho_0(h)$ meredekség lépésköz-érzékenysége miatt $2{,}42203137\,M_\odot$ értéket és $15{,}31652192$-es 3 pontos csúcsot ($-0{,}001354\,M_\odot$, illetve $+1{,}13\times 10^{-4}\;\text{dex}$ eltolódás) ad.

*(Összehasonlítás a korábbi lapos Robin-peremfeltételű `yukawa_action_el` értékekkel SLy $L = 10\;\text{km}$, $\alpha_Y = 1$ esetén: a lapos Robin-közelítésben $M_{\text{surf},\max}^{\text{flat}} = 2{,}03583\,M_\odot$ és $M_{\text{surf},\max}^{\text{flat}} + \Delta M_{\phi,\text{ext}}^{\text{flat}} = 2{,}04481\,M_\odot$ volt $\log_{10}\rho_{0c} = 15{,}31307$-nél; a valódi görbült külső megoldással a felszíni tömegmaximum $M_{\text{surf},\max}^{\text{curved}} = 2{,}02911\,M_\odot$, az egzakt ADM-tömegmaximum pedig $M_{\text{ADM},\max}^{\text{curved}} = \mathbf{2{,}03321\,M_\odot}$ $\log_{10}\rho_{0c} = \mathbf{15{,}316404}$-nál.)*

---

## 4. A csatolt $(\delta h, \delta\phi_L, \delta m)$ radiális pulzációs stabilitás a görbült `yukawa_action_el` sorozaton (V25)

### V25 re-review javítások: 4-komponensű külső BVP, egzakt külső linearizáció, dinamikus skalárhullám és GL2 N-konvergencia

*(Forrás: `curved_yukawa_adm_and_stability.py`, `curved_yukawa_audit_summary.json`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_adm_sequences.csv`.)*

- **Külső lapse (5. hiba):** a külső $\nu$ negyedik ODE-komponensként együtt integrálódik $(\phi_L,\psi_L,m,\nu)$, `solve_bvp` `tol=1e-12`, végfeltétel $\nu(R_{\max})=\tfrac12\ln(1-2m(R_{\max})/R_{\max})$. SLy, $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $N_{\text{mult}}\in\{10,12,14,16\}$ szórások: $\delta k/k=6{,}17\times 10^{-12}$, $\delta M_{\text{ADM}}/M_{\text{ADM}}=5{,}54\times 10^{-13}$, $\delta\nu_R=1{,}63\times 10^{-14}$.

| $N_{\text{mult}}$ | $R_{\max}$ (km) | $k_{\text{curved}}$ (km$^{-1}$) | $M_{\text{ADM}}$ ($M_\odot$) | $\nu_R$ |
|---|---|---|---|---|
| `10` | `111.4657` | `0.233542329280` | `1.509998138706` | `-0.246306687228` |
| `12` | `131.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `14` | `151.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `16` | `171.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |

- **Külső linearizáció és véges differenciás deriválttesztek (SLy, $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $\alpha_Y=1$):** a `linearized_exterior_response` (593. sor: `a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0`) és a belső skalároperátor bal oldali ($L_u u = R_h\,\delta h$) `a_m_sc` együtthatója (`scalar_m_sign = -1.0`: `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n`, amely a jobb oldali $\partial(\phi_L'')/\partial m = +(2\phi_L''/r + 2\phi_L'/r^2)/B$ derivált bal oldalra viteléből adódik) centrális véges differenciákkal összevetve:
  - Külső válasz $\partial\phi_L'(R_*)/\partial m(R_*)$: linearizált BVP **`-1.68365059e-04`** $\text{km}^{-2}$ vs. 4-komponensű külső BVP centrális véges differenciája **`-1.68365059e-04`** $\text{km}^{-2}$ (relatív eltérés: **$3{,}98\times 10^{-10}$**; $K_\phi = 0{,}23358323\;\text{km}^{-1}$, rel. eltérés $1{,}42\times 10^{-10}$; $\partial M_{\text{ADM}}/\partial m(R_*) = 0{,}99947576$, rel. eltérés $4{,}26\times 10^{-10}$).
  - Belső $r = 6{,}0\;\text{km}$ pontbeli skalár-tömeg derivált ($\partial\phi_L''/\partial m$): analitikus **`-1.72833429e-04`** $\text{km}^{-3}$ vs. belső skalár RHS centrális véges differenciája **`-1.72833428e-04`** $\text{km}^{-3}$ (relatív eltérés: **$6{,}08\times 10^{-9}$**).
- **Dinamikus skalárblokk:** a skalárblokkban szerepel a $-\omega^2 e^{-2\nu}\delta\phi_L$ tag, és a külső választ $\omega$-függő BVP adja (`dynamical_scalar_wave=True`); megoldó: GL2 darabonkénti szegmenses Volterra-diszkretizáció, $N=480$, az egzakt bal oldali `scalar_m_sign = -1.0` előjellel. A kvázisztatikus és a dinamikus gyök a közölt 6 tizedesjegyen egyezik (mivel $\omega^2=0$-nál a két operátor azonos).
- **Fordulópontok és zérushelyek:** `minimize_scalar` + 7 pontos kvadratikus illesztés $\pm0{,}002$ dex-en; a zérushelyek a gyök körüli $\pm0{,}0005$ dex-es 3 pontos kvadratikus illesztésből.

| Sorozat | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dinamikus $\omega_0^2=0$ | $\Delta_{\text{ADM}}$ dinamikus (dex) | $\Delta_{M_B}$ dinamikus (dex) | Kvázisztatikus $\omega_0^2=0$ | Befagyasztott $\delta\phi_L=0$ | $\Delta_{\text{ADM}}$ befagyasztott (dex) |
|---|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316409` | **`15.316410`** | **`+5.80e-06`** | `+1.33e-06` | `15.316410` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.307169` | **`15.307170`** | **`+4.20e-06`** | `+6.41e-07` | `15.307170` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303082` | **`15.303080`** | **`-6.86e-06`** | `-2.61e-06` | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.293025` | **`15.293025`** | **`+3.57e-06`** | `-1.40e-07` | `15.293025` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227787` | **`15.227783`** | **`-1.34e-06`** | `-4.05e-06` | `15.227783` | `15.251384` | `+2.36e-02` |

**Konvergenciavizsgálat magán az $\omega_0^2=0$ zérushelyen ($N$ és $R_{\max}$ szerint)** (SLy, $L=10\;\text{km}$, $\alpha_Y=1$; valamint $\omega_0^2$ konvergenciája $\log_{10}\rho_{0c}=15{,}311404$-nél):
- **Zérushely $N$-konvergenciája ($N_{\text{mult}}=14$):** $p_{\text{eff,root}}=4{,}248$ (a három legfinomabb szintből), Richardson-extrapolált gyök $\log_{10}\rho_{0c,\text{Rich}}(\omega_0^2=0) = `15.31640981`$; $N=480$-nál a diszkretizációs eltérés a Richardson-gyöktől $6{,}72\times 10^{-8}\;\text{dex}$ (és $N=240\to 800$ között mind az 5 sorozaton $\le 2{,}02\times 10^{-5}\;\text{dex}$).
- **Zérushely $R_{\max}$-konvergenciája ($N=480$, $R_{\max}=R_*+N_{\text{mult}}L$, $N_{\text{mult}}\in\{10,12,14,16,20\}$, $R_{\max}\in[109{,}89,\;209{,}89]\;\text{km}$):** a gyök szórása $2{,}42\times 10^{-13}\;\text{dex}$ (mind az 5 sorozaton $\le 2{,}72\times 10^{-13}\;\text{dex}$).
- **Rögzített pontbeli ($\log_{10}\rho_{0c}=15{,}311404$) $\omega_0^2$ konvergencia:** $p_{\text{eff}}=5{,}223$, Richardson-érték $\omega_0^2=6369417{,}611\;\text{s}^{-2}$ ($N=480$-nál a relatív eltérés $1{,}29\times 10^{-5}$).

| $N$ | $\log_{10}\rho_{0c}(\omega_0^2=0)$ | $|\Delta\log_{10}\rho_{0c}|$ a Richardson-gyöktől (dex) | $\Delta_{\text{ADM}}$ (dex) | $\Delta_{M_B}$ (dex) | $\omega_0^2(15{,}311404)$ ($\text{s}^{-2}$) | Rel. eltérés ($\omega_0^2$) |
|---|---|---|---|---|---|---|
| `80` | `15.31645808` | `4.83e-05` | `+5.40e-05` | `+4.95e-05` | `6432300.503336` | `9.87e-03` |
| `120` | `15.31642223` | `1.24e-05` | `+1.82e-05` | `+1.37e-05` | `6385729.603235` | `2.56e-03` |
| `160` | `15.31641426` | `4.46e-06` | `+1.02e-05` | `+5.72e-06` | `6375276.486537` | `9.20e-04` |
| `240` | `15.31641079` | `9.80e-07` | `+6.71e-06` | `+2.24e-06` | `6370703.105233` | `2.02e-04` |
| `360` | `15.31641001` | `2.06e-07` | `+5.94e-06` | `+1.47e-06` | `6369680.451299` | `4.13e-05` |
| `480` | `15.31640987` | `6.72e-08` | `+5.80e-06` | `+1.33e-06` | `6369499.539632` | `1.29e-05` |
| `600` | `15.31640983` | `2.56e-08` | `+5.76e-06` | `+1.29e-06` | `6369446.532399` | `4.54e-06` |
| `800` | `15.31640982` | `7.55e-09` | `+5.74e-06` | `+1.27e-06` | `6369424.048161` | `1.01e-06` |
| `1200` | `15.31640981` | `1.35e-09` | `+5.73e-06` | `+1.26e-06` | `6369418.385568` | `1.22e-07` |

**Értelmezési megjegyzés:** Az egzakt bal oldali előjelű (`scalar_m_sign = -1.0`) GL2-operátor gyöke $N$ és $R_{\max}$ szerint $<10^{-7}\;\text{dex}$ pontossággal konvergál (`15.31640981`), az $\omega_0^2=0$ zérushely mind az 5 sorozaton $\le 6{,}9\times 10^{-6}\;\text{dex}$-en belül egybeesik az ADM-tömegmaximummal (SLy $L=10\;\text{km}$, $\alpha_Y=1$ esetén $15{,}316410$ vs. $15{,}316404$, $\Delta_{\text{ADM}} = +5{,}80\times 10^{-6}\;\text{dex}$, $\Delta_{M_B^{\text{quad}}} = +1{,}33\times 10^{-6}\;\text{dex}$), és a sajátfüggvény belső relatív $L^2$-eltérése a statikus érintőtől SLy $L=10\;\text{km}$, $\alpha_Y=1$-re $\delta h = 5{,}04\times 10^{-4}$, $\delta m = 1{,}77\times 10^{-4}$, $\delta\phi_L = 2{,}27\times 10^{-4}$, $\delta\phi_L' = 3{,}09\times 10^{-4}$ (míg a szabadfelszíni elmozdulások relatív különbsége $\Delta\xi(R_*)/\xi = 3{,}54\%$; lásd a 7.4. szakaszt). Mivel a szabadfelszíni elmozdulás $3{,}54\%$-os eltérése és az ODE vs. kvadratúra bariontömeg-fordulópontok különbsége miatt a mínuszjeles gyök és a rácskonvergencia önmagában nem független bizonyítás a globális radiális alapspektrumra és a szabadfelszín-variációra, **a stabilitási záró állítás feltételesen van megfogalmazva addig, amíg a másodrendű variációs levezetés és a belső–külső peremfeltétel-illesztés független ellenőrzése le nem zárul.**

---

## 5. Ellenőrzések (Checks)

1. **Az ADM-tömeg ($M_{\text{ADM}}$), a barionos tömeg ($M_B$) és a csatolt $\omega_0^2 = 0$ zérushely összehasonlítása:**
   - SLy, $L = 10\;\text{km}$, $\alpha_Y = 1$:
     $$\log_{10}\rho_{0c}\!\left(\frac{dM_{\text{ADM}}}{d\rho_c}=0\right) = 15{,}316404, \quad \log_{10}\rho_{0c}\!\left(\frac{dM_B}{d\rho_c}=0\right) = 15{,}316409, \quad \log_{10}\rho_{0c}\bigl(\omega_0^2=0\bigr) = 15{,}316410,$$
     azaz a csatolt dinamikus stabilitási zérushely $+5{,}80\times 10^{-6}\;\text{dex}$-en belül egybeesik az ADM-tömegmaximummal (és $+1{,}33\times 10^{-6}\;\text{dex}$-en belül a barionos tömegmaximummal), míg a további 4 sorozaton az eltérés $+4{,}20\times 10^{-6}\;\text{dex}$ ($\alpha_Y=1/3$), $-6{,}86\times 10^{-6}\;\text{dex}$ ($L=1\;\text{km}$), $+3{,}57\times 10^{-6}\;\text{dex}$ (APR4) és $-1{,}34\times 10^{-6}\;\text{dex}$ (H4) $\checkmark$.
2. **A skalármező-perturbáció ($\delta\phi_L$) elhagyásának fizikai hatása:**
   - Ha a `yukawa_action_el` háttéren befagyasztjuk a skalármezőt ($\delta\phi_L \equiv 0$), a zérushely SLy $L=10\;\text{km},\;\alpha_Y=1$ esetén $\log_{10}\rho_{0c} = 15{,}328565$-re tolódik el ($+1{,}22\times 10^{-2}\;\text{dex}$ eltérés az ADM-csúcsra), mert a skalármező vonzó önenergiájának kompressziós megváltozása ($\delta\phi_L > 0$) nélkül a csillag merevebbnek mutatkozik. A csatolt $(\delta h, \delta\phi_L, \delta m)$ egyenletrendszer ezt az eltolódást **$2097$-szeresével** ($+5{,}80\times 10^{-6}\;\text{dex}$-re, mind az 5 sorozaton $165\times\text{–}17593\times$ faktorral) csökkenti $\checkmark$.

---

## 6. Reprodukálhatóság és generált fájlok (Reproducibility and Generated Files)

A 22. lépés egyetlen paranccsal reprodukálható:
```bash
python3 curved_yukawa_adm_and_stability.py
```
A futtatás előállítja és frissíti az alábbi fájlokat (a gyökérmappában, a `csv/` és a `data/v23/` almappákban):
- `curved_yukawa_adm_sequences.csv`: az ADM-tömegre ($M_{\text{ADM}}$), felszíni tömegre ($M_{\text{surf}}$) és háromféle barionos nyugalmi tömegre ($M_B^{\text{quad}}$, $M_B^{\text{ODE}}$, $M_B^{\text{trapz}}$) újraszámolt 8 darab tömeg–sugár sorozat összefoglaló táblázata.
- `curved_yukawa_coupled_stability.csv`: az 5 darab `yukawa_action_el` sorozat csatolt $(\delta h, \delta\phi_L, \delta m)$ dinamikus, kvázisztatikus és befagyasztott-$\phi_L$ ($\delta\phi_L = 0$) radiális sajátérték-szkenjei és zérushelyei (`scalar_m_sign = -1.0` fő operátorral).
- `curved_yukawa_audit_summary.json`: az $R_{\max}$-konvergencia, a görbült vs. lapos Robin-meredekség, a véges differenciás deriválttesztek, a 108-pontos $4\times 4$ Jacobi-ellenőrzés, a barionos tömeg / első főtétel audit és a sajátmódus–statikus érintő összevetés gépi olvasható jegyzőkönyve.
- `closure_tov_comparison.csv` és `conservation_Residuals.json`: a `yukawa_action_el_curved_adm` sorokkal és a `curved_exterior_yukawa_adm_and_coupled_stability` blokkal kiegészített V23 referenciafájlok.

---

## 7. V25 diagnosztikai elemzés: az előjeljavítások, a három új diagnosztika és a feltételes stabilitási záró állítás (`curved_yukawa_audit_summary.json`)

**Publikációs státusz és verzióelhatárolás (V24 vs. V25):** A korábbi Zenodo-publikáció (**V24**, DOI: `10.5281/zenodo.23680715`) kizárólag a korábbi V24-es kiadásra mutat. A `curved_yukawa_adm_and_stability.py` előjeljavításai (a `linearized_exterior_response` függvényben `a_m = (2.0 * psi0 / (r * r) + 2.0 * dpsi0 / r) / B0`, valamint a `solve_yukawa_coupled_omega2_gl2` bal oldali $L_u u = R_h\,\delta h$ operátorában az egzakt `scalar_m_sign = -1.0`, azaz `a_m_sc = -(2.0 * ddphi_n / r_nodes + 2.0 * dphi_n / (r_nodes ** 2)) / B_n` alapértelmezetté tétele), a `curved_yukawa_audit_summary.json` három új diagnosztikája (`full_jacobian_multi_eos_and_points`, `baryonic_mass_and_first_law_audit`, `eigenmode_xi_nodes_and_static_tangent_audit`), továbbá az elméleti ellenőrzés lezárása kizárólag az új **V25-ös** verzióhoz tartozik — a korábbi V24-es DOI-ra soha nem szabad úgy hivatkozni, mintha e javítások és lezárások már abban szerepeltek volna.

### 7.1. Véges differenciás deriváltellenőrzés (`derivative_verification_tests_SLy_15p0_L10km`)
A külső `linearized_exterior_response` és a belső skalár-tömeg csatolás véges differenciás ellenőrzése (SLy, $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $\alpha_Y=1$):

| Mennyiség | Linearizált BVP / Analitikus | Centrális véges differencia (FD) | Relatív eltérés |
|---|---|---|---|
| `exterior_dphi_prime_dm_R` ($\partial\phi_L'(R_*)/\partial m(R_*)$, $\text{km}^{-2}$) | `-1.68365059e-04` | `-1.68365059e-04` | `3.98e-10` |
| `exterior_K_phi_R` ($-\partial\phi_L'(R_*)/\partial\phi_L(R_*)$, $\text{km}^{-1}$) | `2.33583231e-01` | `2.33583231e-01` | `1.42e-10` |
| `exterior_M_m_R` ($\partial M_{\text{ADM}}/\partial m(R_*)$) | `9.99475763e-01` | `9.99475764e-01` | `4.26e-10` |
| `exterior_M_phi_R` ($\partial M_{\text{ADM}}/\partial\phi_L(R_*)$, $\text{km}$) | `1.27572889e+00` | `1.27572919e+00` | `2.31e-07` |
| `interior_r6km_a_m_sc` ($\partial\phi_L''/\partial m$ az $r=6\;\text{km}$ pontban, $\text{km}^{-3}$) | `-1.72833429e-04` | `-1.72833428e-04` | `6.08e-09` |

### 7.2. Több-EOS, több-sűrűségű és több-sugárpontú csatolt $4\times 4$ Jacobi-ellenőrzés (`full_jacobian_multi_eos_and_points`)
A `verify_full_jacobian_multi_eos_and_points()` függvény 4 konfiguráción (SLy $\log_{10}\rho_{0c}=15{,}0$ és a $15{,}316404$ fordulópont; APR4 $\log_{10}\rho_{0c}=15{,}0$ és a $15{,}293022$ fordulópont), konfigurációnként 3 belső sugárponton ($r \in \{3{,}0,\;6{,}0,\;8{,}5\}\;\text{km}$) és a felszíni illesztésen ellenőrzi a belső $4\times 4$ Jacobi-mátrixot ($\delta h'$, $\delta\phi_L''$, $\delta m'$, $\delta\nu'$ és az adiabatikus $\xi$ elmozdulás együtthatói) centrális véges differenciákkal szemben. Mind a **108 ellenőrzés** legnagyobb relatív hibája `max_relative_error_across_all_108_checks = 2.74e-05` (`all_checks_pass_tol_1e_4 = true`):

| EOS | $\log_{10}\rho_{0c}$ | $K_m$ rel. hiba | $K_\phi$ rel. hiba | $M_m$ rel. hiba | $M_\phi$ rel. hiba | Belső sugárpontok $r$ (km) | Max. belső rel. hiba |
|---|---|---|---|---|---|---|---|
| SLy | `15.000000` | `3.98e-10` | `1.42e-10` | `4.26e-10` | `2.31e-07` | `3.0, 6.0, 8.5` | `2.74e-05` |
| SLy | `15.316404` | `7.12e-10` | `4.46e-10` | `4.64e-10` | `7.40e-08` | `3.0, 6.0, 8.5` | `1.88e-05` |
| APR4 | `15.000000` | `1.42e-09` | `6.91e-10` | `5.42e-10` | `8.70e-08` | `3.0, 6.0, 8.5` | `1.94e-05` |
| APR4 | `15.293022` | `5.79e-10` | `4.84e-10` | `6.21e-10` | `1.12e-07` | `3.0, 6.0, 8.5` | `1.59e-05` |

### 7.3. Barionos nyugalmi tömeg ($M_B$) és első főtétel-audit (`baryonic_mass_and_first_law_audit`)
A `verify_baryonic_mass_and_first_law()` függvény három független eljárással számítja ki a barionos nyugalmi tömeget ($M_B = \int_0^{R_*} 4\pi r^2 B(r)^{-1/2}\rho_0(r)\,dr$): (i) utólagos trapézszabállyal (`M_B_trapz`), (ii) az $(s=r^2)$ TOV-rendszer 5. állapotváltozójaként (`dM_B/ds = 2*pi*sqrt(s)*rho0/sqrt(B)`, `rtol=1e-12`, `atol=1e-14`, `M_B_ode`), valamint (iii) $s=r^2$ PCHIP-interpolált adaptív Gauss–Kronrod kvadratúrával (`quad`, `epsrel=1e-12`, `M_B_quad`). Rácskonvergencia (SLy, $L=10\;\text{km}$, $\alpha_Y=1$):

| `ngrid` | 7-pontos csúcs $M_B^{\text{trapz}}$ | 7-pontos csúcs $M_B^{\text{quad}}$ | 7-pontos csúcs $M_{\text{ADM}}$ | 3-pontos csúcs $M_B^{\text{ODE}}$ |
|---|---|---|---|---|
| `301` | `15.316389` | `15.316409` | `15.316404` | `15.316522` |
| `601` | `15.316471` | `15.316410` | `15.316414` | `15.316522` |
| `1201` | `15.316433` | `15.316410` | `15.316406` | `15.316523` |
| `2401` | `15.316418` | `15.316410` | `15.316415` | `15.316522` |

Az $s=r^2$ adaptív Gauss–Kronrod kvadratúrával számított $M_B^{\text{quad}}$ és az $M_{\text{ADM}}$ fordulópontja mind a 8 sorozaton $\le 8{,}43\times 10^{-6}\;\text{dex}$ pontossággal egyezik (`delta_dex_ADM_vs_baryon_turning` $-4{,}58\times 10^{-6}$ és $+8{,}43\times 10^{-6}$ dex között). A fordulópont körüli 9 sűrűségponton végzett első főtétel-próba (`dM_ADM/dx - e^{\nu(R_*) - \phi_L(R_*)} dM_B/dx`) számszerűsíti a nem-konform ($\mathcal{L}_m = -e^{-\phi_L}\varepsilon$, `USE_TRACE_SOURCE=False`) csatolásból származó egzakt maradékot:
$$\frac{dM_{\text{ADM}}}{dx} - e^{\nu(R_*) - \phi_L(R_*)}\frac{dM_B}{dx} \;=\; \int_0^{R_*} \mathcal{O}_{\text{nc}}(r)\,\frac{\partial\phi_L}{\partial x}\,dr \;\approx\; +8{,}23\times 10^{-5}\,M_\odot/\text{dex}.$$

### 7.4. Sajátmódus-rekonstrukció ($\xi(r)$, csomópontszám) és összevetés a statikus egyensúlyi érintővel (`eigenmode_xi_nodes_and_static_tangent_audit`)
A `verify_eigenmode_xi_nodes_and_static_tangent()` függvény rekonstruálja a csatolt alapmódus $\xi(r)$ radiális elmozdulásprofilját, ellenőrzi a belső csomópontszámot (`n_nodes = 0`) és a szabad felszíni Lagrange-nyomásfeltételt ($\Delta P(R_*) = -(\varepsilon+P)r\,\delta h(R_*) = 0$), majd összeveti az $\omega_0^2=0$ sajátmódust a statikus egyensúlyi sorozat $\partial_{\log_{10}\rho_{0c}}(h, m, \phi_L, \phi_L')$ érintővektorával mind a fő (`scalar_m_sign = -1.0`, egzakt LHS), mind a pluszjeles (`scalar_m_sign = +1.0`) operátor-előjelkonvencióra:

| Sorozat | Fő zérushely (`scalar_m_sign = -1.0`, egzakt LHS) | $\Delta_{\text{ADM}}$ (dex) | Diagnosztikai zérushely (`scalar_m_sign = +1.0`) | $\Delta_{\text{ADM}}$ (dex) | $\xi(r)$ csomópontok (`-`) | $\xi(R_*)$ szabadfelszín vs. statikus érintő (km) | Felszíni $\xi(R_*)$ rel. eltérés | $\xi/r$ centrális regul. |
|---|---|---|---|---|---|---|---|---|
| `SLy_L10km_alphaY_1.0000` | **`15.316410`** | **`+5.813e-06`** | `15.319129` | `+2.725e-03` | `0` | `-2.8994` vs. `-2.8003` | `3.541e-02` (`3.54%`) | `1.08e-03` |
| `SLy_L10km_alphaY_0.3333` | **`15.307170`** | **`+4.089e-06`** | `15.308114` | `+9.479e-04` | `0` | `-2.9237` vs. `-2.8986` | `8.666e-03` (`0.87%`) | `1.03e-03` |
| `SLy_L1km_alphaY_1.0000` | **`15.303080`** | **`-6.412e-06`** | `15.303088` | `+2.212e-06` | `0` | `-2.9406` vs. `-2.9411` | `1.767e-04` (`0.018%`) | `1.01e-03` |
| `APR4_L10km_alphaY_1.0000` | **`15.293025`** | **`+3.065e-06`** | `15.294175` | `+1.153e-03` | `0` | `-2.1155` vs. `-2.0681` | `2.289e-02` (`2.29%`) | `1.29e-03` |
| `H4_L10km_alphaY_1.0000` | **`15.227783`** | **`-1.658e-06`** | `15.234109` | `+6.324e-03` | `0` | `-5.5394` vs. `-5.5623` | `4.116e-03` (`0.41%`) | `5.30e-04` |

Relatív $L^2$-normabeli eltérés a statikus egyensúlyi érintő ($\partial_{\log_{10}\rho_{0c}}$) és az $\omega_0^2=0$ sajátfüggvény között a csillag térfogatán:

| Sorozat | Operátor-variáns | $\|\delta h - \partial h\|_2$ | $\|\delta m - \partial m\|_2$ | $\|\delta\phi_L - \partial\phi_L\|_2$ | $\|\delta\phi_L' - \partial\phi_L'\|_2$ |
|---|---|---|---|---|---|
| `SLy_L10km_alphaY_1.0000` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, fő eredmény) | `5.043e-04` | `1.767e-04` | `2.273e-04` | `3.093e-04` |
| `SLy_L10km_alphaY_1.0000` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnosztika) | `5.026e-03` | `4.099e-03` | `9.347e-01` | `5.974e-01` |
| `SLy_L10km_alphaY_0.3333` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, fő eredmény) | `1.319e-04` | `1.588e-04` | `5.441e-04` | `3.635e-04` |
| `SLy_L10km_alphaY_0.3333` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnosztika) | `1.525e-03` | `1.341e-03` | `9.397e-01` | `5.967e-01` |
| `SLy_L1km_alphaY_1.0000` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, fő eredmény) | `5.177e-05` | `1.638e-04` | `1.709e-04` | `2.478e-03` |
| `SLy_L1km_alphaY_1.0000` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnosztika) | `4.052e-05` | `1.468e-04` | `2.243e-02` | `2.143e-02` |
| `APR4_L10km_alphaY_1.0000` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, fő eredmény) | `3.048e-04` | `2.108e-04` | `1.239e-03` | `6.328e-04` |
| `APR4_L10km_alphaY_1.0000` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnosztika) | `4.077e-02` | `6.161e-03` | `1.112e+00` | `7.202e-01` |
| `H4_L10km_alphaY_1.0000` | `exact_minus_a_m_sc` (`scalar_m_sign=-1.0`, fő eredmény) | `7.340e-05` | `1.885e-04` | `6.043e-05` | `1.400e-04` |
| `H4_L10km_alphaY_1.0000` | `code_plus_a_m_sc` (`scalar_m_sign=+1.0`, diagnosztika) | `1.059e-01` | `1.964e-03` | `5.713e-01` | `3.420e-01` |

> **Összegzés és az elméleti lezárás határa (`mathematical_origin_of_0p002725_dex_offset`):** A nemlineáris belső skaláregyenletben ($B\phi_L'' + C_0\phi_L' = \phi_L/L^2 - \alpha_Y e^{-\phi_L}\varepsilon$) a jobb oldal (RHS) $m$ szerinti parciális deriváltja rögzített $(h, \phi_L, \phi_L')$ mellett $\partial(\phi_L'')/\partial m = +(2\phi_L''/r + 2\phi_L'/r^2)/B$ (az $r=6\;\text{km}$ pontban $-1{,}728334\times 10^{-4}\;\text{km}^{-3}$). A `solve_yukawa_coupled_omega2_gl2` Volterra-rendszer bal oldalára ($L_u u = R_h\,\delta h$) átrendezve az egzakt bal oldali előjel `scalar_m_sign = -1.0` (`a_m_sc = -\partial(\phi_L'')/\partial m`), amellyel a V25 fő eredményei (`curved_yukawa_coupled_stability.csv`, `curved_yukawa_audit_summary.json`) készültek: mind az 5 sorozaton csomópontmentes ($n=0$) $\xi(r)$ alapmódust és az ADM-tömegmaximummal $\le 6{,}9\times 10^{-6}\;\text{dex}$-en belül egybeeső $\omega_0^2=0$ zérushelyet (`15.316410`) szolgáltat (míg a `scalar_m_sign = +1.0` előjel $+0{,}002725\;\text{dex}$-szel, `15.319129`-re tolja el a gyököt). SLy $L=10\;\text{km}$, $\alpha_Y=1$ esetén a sajátfüggvény és a statikus érintő relatív $L^2$-eltérései $\delta h = 5{,}04\times 10^{-4}$, $\delta m = 1{,}77\times 10^{-4}$, $\delta\phi_L = 2{,}27\times 10^{-4}$, $\delta\phi_L' = 3{,}09\times 10^{-4}$, míg a szabadfelszíni elmozdulások relatív különbsége $3{,}54\%$ (`rel_diff_xi_R_exact_minus_vs_static_tangent = 0.035406`). A jobb, mínuszjeles gyök és a rácskonvergencia komoly előrelépés, de **nem önmagában független bizonyítás a globális radiális alapspektrumra és a szabadfelszín-variációra**, ezért **a stabilitási záró állítás feltételesen van megfogalmazva addig, amíg a variációs és peremfeltételi ellenőrzés függetlenül le nem zárul.**
