# 22. lépés: A lokális Yukawa-ág (`yukawa_action_el`) görbült téridőbeli külső Einstein–skalár vákuummegoldása, az aszimptotikus ADM-tömeg ($M_{\text{ADM}}$), a tömeg–sugár görbék újraszámolása és a csatolt $(\delta h, \delta\phi_L, \delta m)$ radiális stabilitás
# Step 22: Curved-Spacetime Exterior Einstein-Scalar Vacuum Solution of the Local Yukawa Branch (`yukawa_action_el`), Asymptotic ADM Mass ($M_{\text{ADM}}$), Recalculation of Mass-Radius Curves, and Coupled $(\delta h, \delta\phi_L, \delta m)$ Radial Stability

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/step22_curved_exterior_yukawa_adm_mass_mr_and_coupled_stability.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step21_enthalpy_tov_even_parity_s_r2_regularity_and_exterior_halo_lapse.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források és futtatható kód / Local Sources and Runnable Code:**  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1.2. és 2. szakasz, 65–153. sor; DOI: `10.5281/zenodo.23639084`)  
- Futtatható megoldókód: `/workspace/ta-v27-work/scalar-tensor-yukawa-fit/curved_yukawa_adm_and_stability.py` (és integrálva a `covariant_closure_50aa2db0.py` szkriptbe)  
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

---

## 3. A tömeg–sugár görbék újraszámolása az ADM-tömeggel ($M_{\text{ADM}}$) és a barionos nyugalmi tömeggel ($M_B$)

Az alábbi táblázat (`curved_yukawa_adm_sequences.csv` és `closure_tov_comparison.csv`) mind a 8 vizsgált `yukawa_action_el` sorozatra megadja az **ADM-tömegre ($M_{\text{ADM}}$) optimalizált** értékeket, együtt a felszíni tömeggel ($M_{\text{surf}}$) és a konzervált barionos nyugalmi tömeg ($M_B = \int_0^{R_*} 4\pi r^2 B(r)^{-1/2}\rho_0(r)\,dr$) fordulópontjával:

| EOS | $L$ (km) | $\alpha_Y$ | $M_{\text{ADM},\max}$ ($M_\odot$) | $M_{\text{surf},\max}$ ($M_\odot$) | $\Delta M_{\phi,\text{ext}}^{\text{curved}}(M_{\max})$ ($M_\odot$) | $\log_{10}\rho_{0c}(M_{\text{ADM},\max})$ | $\log_{10}\rho_{0c}(M_{B,\max})$ | $R(M_{\text{ADM},\max})$ (km) | $R(M_{\text{ADM}}=1{,}4\,M_\odot)$ (km) | $R(M_{\text{surf}}=1{,}4\,M_\odot)$ (km) | $M_{\text{ADM}}(15{,}0)$ ($M_\odot$) | $f_L(0; 15{,}0)$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **SLy** | `1.0` | `1.000` | **`2.04684`** | `2.04684` | `2.23e-6` | **`15.30309`** | `15.30333` | `9.9689` | **`11.6933`** | `11.6933` | `1.56171` | `0.99920` |
| **SLy** | `1.0` | `0.333` | **`2.04789`** | `2.04789` | `7.43e-7` | **`15.30266`** | `15.30289` | `9.9733` | **`11.7018`** | `11.7018` | `1.56386` | `0.99973` |
| **SLy** | `10.0` | `1.000` | **`2.03321`** | `2.02911` | `4.11e-3` | **`15.31640`** | `15.31639` | `9.8853` | **`11.5359`** | `11.5348` | `1.51000` | `0.98238` |
| **SLy** | `10.0` | `0.333` | **`2.04321`** | `2.04183` | `1.38e-3` | **`15.30717`** | `15.30689` | `9.9446` | **`11.6475`** | `11.6471` | `1.54638` | `0.99397` |
| **APR4** | `1.0` | `1.000` | **`2.18693`** | `2.18693` | `3.04e-6` | **`15.28281`** | `15.28307` | `9.8947` | **`11.3085`** | `11.3085` | `1.56316` | `0.99920` |
| **APR4** | `10.0` | `1.000` | **`2.18068`** | `2.17624` | `4.45e-3` | **`15.29302`** | `15.29299` | `9.8462` | **`11.1830`** | `11.1828` | `1.51271` | `0.98210` |
| **H4** | `1.0` | `1.000` | **`2.01139`** | `2.01138` | `7.71e-7` | **`15.20975`** | `15.20950` | `11.6917` | **`13.9388`** | `13.9388` | `1.88120` | `0.99917` |
| **H4** | `10.0` | `1.000` | **`1.97809`** | `1.97502` | `3.08e-3` | **`15.22778`** | `15.22798` | `11.5289` | **`13.7179`** | `13.7169` | `1.82871` | `0.98108` |

*(Összehasonlítás a korábbi lapos Robin-peremfeltételű `yukawa_action_el` értékekkel SLy $L = 10\;\text{km}$, $\alpha_Y = 1$ esetén: a lapos Robin-közelítésben $M_{\text{surf},\max}^{\text{flat}} = 2{,}03583\,M_\odot$ és $M_{\text{surf},\max}^{\text{flat}} + \Delta M_{\phi,\text{ext}}^{\text{flat}} = 2{,}04481\,M_\odot$ volt $\log_{10}\rho_{0c} = 15{,}31307$-nél; a valódi görbült külső megoldással a felszíni tömegmaximum $M_{\text{surf},\max}^{\text{curved}} = 2{,}02911\,M_\odot$, az egzakt ADM-tömegmaximum pedig $M_{\text{ADM},\max}^{\text{curved}} = \mathbf{2{,}03321\,M_\odot}$ $\log_{10}\rho_{0c} = \mathbf{15{,}316404}$-nál.)*

---

## 4. A csatolt $(\delta h, \delta\phi_L, \delta m)$ radiális pulzációs stabilitás a görbült `yukawa_action_el` sorozaton (V27 újraszámolás)

### V27 re-review javítások: 4-komponensű külső BVP, egzakt külső linearizáció, dinamikus skalárhullám és GL2 N-konvergencia

*(Forrás: `curved_yukawa_adm_and_stability.py`, `curved_yukawa_audit_summary.json`, `curved_yukawa_coupled_stability.csv`, `curved_yukawa_adm_sequences.csv`.)*

- **Külső lapse (5. hiba):** a külső $\nu$ negyedik ODE-komponensként együtt integrálódik $(\phi_L,\psi_L,m,\nu)$, `solve_bvp` `tol=1e-12`, végfeltétel $\nu(R_{\max})=\tfrac12\ln(1-2m(R_{\max})/R_{\max})$. SLy, $\log_{10}\rho_{0c}=15{,}0$, $L=10\;\text{km}$, $N_{\text{mult}}\in\{10,12,14,16\}$ szórások: $\delta k/k=6{,}17\times 10^{-12}$, $\delta M_{\text{ADM}}/M_{\text{ADM}}=5{,}54\times 10^{-13}$, $\delta\nu_R=1{,}63\times 10^{-14}$.

| $N_{\text{mult}}$ | $R_{\max}$ (km) | $k_{\text{curved}}$ (km$^{-1}$) | $M_{\text{ADM}}$ ($M_\odot$) | $\nu_R$ |
|---|---|---|---|---|
| `10` | `111.4657` | `0.233542329280` | `1.509998138706` | `-0.246306687228` |
| `12` | `131.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `14` | `151.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |
| `16` | `171.4657` | `0.233542329278` | `1.509998138707` | `-0.246306687228` |

- **Külső linearizáció (2. hiba):** a `linearized_exterior_response` a $(\delta\phi_L(R_*),\delta m(R_*),\omega^2)$ szerinti egzakt elsőrendű variációs BVP-ből adja $\delta\phi_L'(R_*)$-t, $\delta m(R_*)$-t és $\delta M_{\text{ADM}}=\delta m(\infty)$-t; a sajátérték-feltétel $\delta M_{\text{ADM}}=0$.
- **Dinamikus skalárblokk (1. és 3. hiba):** a skalárblokkban szerepel a $-\omega^2 e^{-2\nu}\delta\phi_L$ tag, és a külső választ $\omega$-függő BVP adja (`dynamical_scalar_wave=True`); megoldó: GL2 darabonkénti szegmenses Volterra-diszkretizáció, $N=480$. A kvázisztatikus és a dinamikus gyök a közölt 6 tizedesjegyen egyezik.
- **Fordulópontok:** `minimize_scalar` + 7 pontos kvadratikus illesztés $\pm0{,}002$ dex-en; a zérushelyek 3 pontos kvadratikus illesztésből.

| Sorozat | $\log_{10}\rho_{0c}(dM_{\text{ADM}}=0)$ | $\log_{10}\rho_{0c}(dM_B=0)$ | Dinamikus $\omega_0^2=0$ | $\Delta$ dinamikus (dex) | Kvázisztatikus $\omega_0^2=0$ | Befagyasztott $\delta\phi_L=0$ | $\Delta$ befagyasztott (dex) |
|---|---|---|---|---|---|---|---|
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.316404` | `15.316388` | **`15.316425`** | **`+2.04e-05`** | `15.316425` | `15.328565` | `+1.22e-02` |
| **SLy, $L=10\;\text{km}$, $\alpha_Y=1/3$** | `15.307166` | `15.306886` | **`15.307172`** | **`+5.84e-06`** | `15.307172` | `15.311340` | `+4.17e-03` |
| **SLy, $L=1\;\text{km}$, $\alpha_Y=1$** | `15.303086` | `15.303330` | **`15.303080`** | **`-6.86e-06`** | `15.303080` | `15.304220` | `+1.13e-03` |
| **APR4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.293022` | `15.292989` | **`15.293042`** | **`+2.07e-05`** | `15.293042` | `15.301564` | `+8.54e-03` |
| **H4, $L=10\;\text{km}$, $\alpha_Y=1$** | `15.227785` | `15.227981` | **`15.227790`** | **`+5.71e-06`** | `15.227790` | `15.251385` | `+2.36e-02` |

**GL2 $N$-konvergencia** (SLy, $L=10\;\text{km}$, $\alpha_Y=1$, $\log_{10}\rho_{0c}=15{,}311404$): $p_{\text{eff}}=5{,}23$ (a három legfinomabb szintből), Richardson-érték $\omega_0^2=6387684{,}936\;\text{s}^{-2}$; $N=480$-nál a relatív hiba $\approx1{,}28\times 10^{-5}$. A sorozat monoton, előjelváltás nélkül konvergál (a korábbi FD-megoldó $N$-függő előjelváltása nem jelentkezik).

| $N$ | $\omega_0^2$ ($\text{s}^{-2}$) | relatív eltérés a Richardson-értéktől |
|---|---|---|
| `80` | `6450570.762819` | `9.84e-03` |
| `120` | `6403998.296759` | `2.55e-03` |
| `160` | `6393544.377742` | `9.17e-04` |
| `240` | `6388970.545464` | `2.01e-04` |
| `360` | `6387947.756886` | `4.11e-05` |
| `480` | `6387766.842088` | `1.28e-05` |
| `600` | `6387713.828704` | `4.52e-06` |
| `800` | `6387691.354875` | `1.00e-06` |
| `1200` | `6387685.705914` | `1.21e-07` |

**Értelmezési korlát:** az $\omega_0^2=0$ és az $M_{\text{ADM}}$-maximum egyezése ugyanazon egyenletrendszer belső konzisztenciavizsgálata; a $\sim10^{-5}$ dex eltérés a 3 pontos gyökillesztés és a 7 pontos fordulópont-illesztés felbontásán belül értendő, nem önálló fizikai teszt.

---

## 5. Ellenőrzések (Checks)

1. **Az ADM-tömeg ($M_{\text{ADM}}$), a barionos tömeg ($M_B$) és a csatolt $\omega_0^2 = 0$ zérushely hármas egyezése:**
   - SLy, $L = 10\;\text{km}$, $\alpha_Y = 1$:
     $$\log_{10}\rho_{0c}\!\left(\frac{dM_{\text{ADM}}}{d\rho_c}=0\right) = 15{,}316404, \quad \log_{10}\rho_{0c}\!\left(\frac{dM_B}{d\rho_c}=0\right) = 15{,}316388, \quad \log_{10}\rho_{0c}\bigl(\omega_0^2=0\bigr) = 15{,}316425,$$
     azaz a csatolt dinamikus stabilitási zérushely **$2{,}04\times 10^{-5}\;\text{dex}$** eltéréssel közelíti az ADM-tömeg maximumával (és $3{,}61\times 10^{-5}\;\text{dex}$ pontossággal a barionos tömeg maximumával) $\checkmark$.
2. **A skalármező-perturbáció ($\delta\phi_L$) elhagyásának fizikai hatása:**
   - Ha a `yukawa_action_el` háttéren befagyasztjuk a skalármezőt ($\delta\phi_L \equiv 0$), a zérushely $\log_{10}\rho_{0c} = 15{,}328565$-re tolódik el ($+1{,}22\times 10^{-2}\;\text{dex}$ hiba), mert a skalármező vonzó önenergiájának kompressziós megváltozása ($\delta\phi_L > 0$) nélkül a csillag látszólag merevebbnek mutatkozik. A csatolt $(\delta h, \delta\phi_L, \delta m)$ egyenletrendszer ezt a $+1{,}22\times 10^{-2}\;\text{dex}$ hibát **$595$-szörösére** csökkenti $\checkmark$.
