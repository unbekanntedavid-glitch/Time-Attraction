# V26 Hiánypótló Számítás: Ágazonosság (`Branch II-Conf` vs. `Branch II-Fluid`) és a SLy $\log_{10}(\rho_{0c}) = 15{,}0$ Csillag-BVP Célzott Numerikus Ellenőrzése

## 1. Ágazonosság és a két 2. ágbeli realizáció (`Branch II-Fluid` vs. `Branch II-Conf`) szétválasztása

A V26 dokumentáció (`step23_gw170817_and_lensing.md` 1.1. és 3.3. szakasz, `step24_bullet_cluster.md` 1.1. szakasz, `step25_cmb_and_cosmology.md` 1.1. szakasz) expliciten rögzíti:
- **`Branch II-Fluid` (22. lépésbeli relativisztikus folyadék-csatolás, $\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$):** a belső skaláregyenlet jobb oldala $-(\alpha_Y G/c^2)\,e^{-\phi_L}\varepsilon(r)$, amely a V25 22. lépésének (`curved_yukawa_adm_and_stability.py`) neutroncsillag-hátterét adja.
- **`Branch II-Conf` (23–25. lépésbeli univerzális konform Jordan-ág, $\tilde{g}_{\mu\nu} = e^{-2\alpha_Y\phi_L}g_{\mu\nu}$):** a $\sqrt{-\tilde{g}} = e^{-4\alpha_Y\phi_L}\sqrt{-g}$ térfogatelemből következően a belső skaláregyenlet forrása a nyom: $-(\alpha_Y^2 G/c^2)\,e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$, miközben a külső vákuumtartományban ($r > R_*$, ahol $\tilde{\varepsilon} = \tilde{P} = 0$) mindkét ág egzaktul ugyanarra a 4-komponensű görbült Einstein–skalár BVP-re redukálódik.

### 1.1. Összehasonlító táblázat (`Branch II-Fluid` vs. `Branch II-Conf`, SLy $\log_{10}\rho_{0c} = 15{,}0$, $L = 10\;\text{km}$, $\alpha_Y = 1$)

| Mennyiség / Megfigyelhető | `Branch II-Fluid` (22. lépés háttér-diagnosztika) | `Branch II-Conf` (Önkonzisztens nyom-forrásolt BVP) | Eltérés fizikai eredete |
|---|---|---|---|
| Belső anyagi hatás | $\mathcal{L}_m = -e^{-\phi_L}\varepsilon(n)$ | $S_m[e^{-2\alpha_Y\phi_L}g_{\mu\nu},\Psi_m]$ | Közvetlen energiasűrűség-súlyozás vs. univerzális Jordan-metrika |
| Belső skalárforrás jobb oldala | $-(\alpha_Y G/c^2)\,e^{-\phi_L}\varepsilon(r)$ | $-(\alpha_Y^2 G/c^2)\,e^{-4\alpha_Y\phi_L}(\tilde{\varepsilon} - 3\tilde{P})$ | Relativisztikus nyomás levonása a nyomban ($1 - 3\tilde{P}_c/\tilde{\varepsilon}_c = 0{,}433935$) |
| Csillagsugár $R_*$ ($\text{km}$) | `11.46572423` | `11.66269842` | Az $e^{-4\phi_L}$ súlyfaktor és a nyomforrás eltolja az egyensúlyi sugarat |
| Felszíni tömeg $M_{\text{surf}} = m(R_*)$ ($M_\odot$) | `1.50766122` | `1.56063550` | Belső tömegintegrál az $e^{-4\alpha_Y\phi_L}\tilde{\varepsilon}$ sűrűséggel |
| Aszimptotikus ADM-tömeg $M_{\text{ADM}}$ ($M_\odot$) | `1.50999814` | `1.56203491` | $M_{\text{surf}} + \Delta M_{\phi,\text{ext}}^{\text{curved}}$ |
| Külső skalárfelhő-tömeg $\Delta M_{\phi,\text{ext}}^{\text{curved}}$ ($M_\odot$) | `+0.00233692` (`+0.155003%`) | `+0.00139941` (`+0.089669%`) | Kisebb felszíni skalármező a $-3\tilde{P}$ nyom-csökkentés miatt |
| Centrális / felszíni skalármező $\phi_L(0)$ / $\phi_L(R_*)$ | `0.01777896` / `0.00540868` | `0.01157688` / `0.00414318` | A magban $1 - 3\tilde{P}/\tilde{\varepsilon} < 1$ |
| Görbült felszíni meredekség $\phi_{L,\text{curved}}'(R_*)$ ($\text{km}^{-1}$) / $k_{\text{curved}}$ ($\text{km}^{-1}$) | `-0.00126316` / `0.23354233` | `-0.00096508` / `0.23293146` | Görbült Dirichlet–Neumann külső illesztés |
| Saját súroló impaktparaméter $b_{\text{graz}}$ ($\text{km}$) / $r_0$ ($\text{km}$) | `14.749032` / `11.558475` | `15.064914` / `11.736304` | $1{,}0001\,R_* e^{-(\nu(R_*) - \phi_L(R_*))}$-nél kiértékelve |
| Saját súroló elhajlási szög $\hat{\alpha}_{\text{Yuk,conf}}(b_{\text{graz}})$ ($\text{rad}$) | `1.29085207` | `1.33402365` | Konform null-geodetikus a saját felszínt súroló sugáron |
| Közös $b = 15{,}064914\;\text{km}$-nél számolt elhajlási szög $\hat{\alpha}_{\text{Yuk,conf}}(b)$ ($\text{rad}$) | `1.22474419` ($r_0 = 11{,}917144\;\text{km}$) | `1.33402365` ($r_0 = 11{,}736304\;\text{km}$) | Azonos impaktparaméternél kiértékelt összehasonlítás |

---

## 2. Nyitott numerikus ellenőrzés lezárása a `Branch II-Conf` SLy $\log_{10}(\rho_{0c}) = 15{,}0$ csillag-BVP-re (`gw170817_and_lensing_audit.py`)

Az alábbi három pont a harmadik ellenőrzésben megjelölt numerikus kérdéseket zárja le közvetlenül a `gw170817_and_lensing_audit.py` (`compute_branch_ii_conf_bvp_convergence_and_mesh_refinement`) futási naplója és a `gw170817_and_lensing_summary.json` alapján.

### 2.1. Konvergált belső diszkrét egyenletmaradék $\max_{r\in(0,R_*)}\lvert\text{res}_{\text{ODE}}\rvert$ ($\text{km}^{-2}$) és külső normalizált kollokációs maradék

A belső tartományban ($r \in (0, R_*)$) a görbült téridőbeli Helmholtz–Yukawa skaláregyenlet diszkrét maradéka:

$$\text{res}_{\text{ODE}}(r_i) \;\equiv\; B(r_i)\,\phi_L''(r_i) + \left(\frac{2B(r_i)}{r_i} + B(r_i)\nu'(r_i) + \frac{1}{2}B'(r_i)\right)\phi_L'(r_i) - \frac{\phi_L(r_i)}{L^2} + \alpha_Y^2\,e^{-4\alpha_Y\phi_L(r_i)}\bigl(\tilde{\varepsilon}(r_i) - 3\tilde{P}(r_i)\bigr) \;=\; 0.$$

A konvergált $N = 601$ rácson ($h = R_*/(N-1) = 0{,}01943783\;\text{km}$) a másodrendű centrális differenciaséma $\text{km}^{-2}$ egységű diszkrét maradéka és a külső 4-komponensű BVP dimenziótlan normalizált kollokációs maradéka:

- **Közvetlen tridiagonális sávos megoldás (`_solve_curved_scalar_with_k`) a konvergált metrikán** (amely a saját diszkrét lineáris rendszerének megoldási pontosságát méri az adott rácson, míg a folytonos BVP diszkretizációs hibáját a 2.3. szakaszbeli $N = 601 \to 2000 \to 4000$ hálófinomítás teszteli):
  $$\max_{1 \le i \le N-2}\lvert\text{res}_{\text{ODE,banded}}(r_i)\rvert \;=\; 1{,}480297\times 10^{-14}\;\text{km}^{-2}, \qquad \text{RMS}\bigl(\text{res}_{\text{ODE,banded}}\bigr) \;=\; 1{,}917229\times 10^{-15}\;\text{km}^{-2}.$$
- **Relaxált fixpont-iterált profil (`it in (1, 3, 6, 10)`, `tol_iter = 2e-8`, 19 iteráció):**
  $$\max_{1 \le i \le N-2}\lvert\text{res}_{\text{ODE,iter}}(r_i)\rvert \;=\; 1{,}078079\times 10^{-10}\;\text{km}^{-2}, \qquad \max_{0 \le i \le N-1}\lvert\phi_L(r_i) - \phi_{L,\text{banded}}(r_i)\rvert \;=\; 3{,}065412\times 10^{-9}.$$
- **Minden iterációban frissített 4-komponensű Robin-perem (`it >= 1`, `tol_iter = 1e-11`, $N = 601$, 36 iteráció):**
  $$\max_{1 \le i \le N-2}\lvert\text{res}_{\text{ODE,iter}}(r_i)\rvert \;=\; 6{,}958172\times 10^{-14}\;\text{km}^{-2}, \qquad \max_{0 \le i \le N-1}\lvert\phi_L(r_i) - \phi_{L,\text{banded}}(r_i)\rvert \;=\; 2{,}174330\times 10^{-12}.$$
- **Külső 4-komponensű görbült Einstein–skalár BVP (`scipy.integrate.solve_bvp`, `tol = 1e-12`, $r \in [R_*, R_{\max}]$):**
  $$\max\,\text{rms\_residuals} \;=\; 9{,}997932\times 10^{-13}$$
  (a külső kollokációs háló dimenziótlan normalizált maradéka, amely nem azonos a belső $\text{km}^{-2}$ egységű diszkrét egyenletmaradékkal).

### 2.2. A belső és külső $\phi'_L(R_*)$ (kikényszerített Robin-derivált vs. külön egyoldali FD2 derivált) számszerű egyezése $r = R_*$-nál és az iterációs konvergenciatörténet

A felszínen ($r = R_*$) élesen meg kell különböztetni a két belső derivált-fogalmat:
1. **Kikényszerített Robin-peremderivált ($\phi'_{L,\text{int,Robin}}(R_*) = -k_{\text{int}}\,\phi_L(R_*)$ és a szellempontos centrális séma $\phi'_{L,\text{int,ghost}}(R_*)$):** azt méri, hogy a belső megoldóban előírt $k_{\text{int}}$ Robin-peremfeltétel milyen pontossággal egyezik a külső 4-komponensű BVP $k_{\text{curved}} \equiv -\phi'_{L,\text{ext,curved}}(R_*)/\phi_L(R_*)$ Dirichlet–Neumann-együtthatójával. A **`3.937e-15 km^-1`** (`1.430e-14 km^-1` szellempontos) egyezés erre a **kikényszerített Robin-deriváltra** igaz.
2. **A peremfeltételtől elkülönített, 3-pontos egyoldali másodrendű véges differenciás derivált ($\phi'_{L,\text{int,FD2}}(R_*) = (3\phi_{N-1} - 4\phi_{N-2} + \phi_{N-3})/(2h)$):** a belső rácspontokból számolt FD2 derivált diszkretizációs eltérése $\lvert\phi'_{L,\text{int,FD2}}(R_*) - \phi'_{L,\text{ext,curved}}(R_*)\rvert = 2{,}1213\times 10^{-8}\;\text{km}^{-1}$ ($N = 601$), $1{,}9048\times 10^{-9}\;\text{km}^{-1}$ ($N = 2000$) és **$4{,}7228\times 10^{-10}\;\text{km}^{-1}$ ($N = 4000$-nél**, $\approx 4{,}89\times 10^{-7}$ relatív eltérés a $\lvert\phi'_L(R_*)\rvert = 9{,}6507\times 10^{-4}\;\text{km}^{-1}$ meredekséghez mérve), és a háló finomításával másodrendben ($\mathcal{O}(h^2)$) konvergál (lásd a 2.3. szakaszt).

| Robin-frissítési ütemezés ($N = 601$) | Iterációk | $k_{\text{int}}$ ($\text{km}^{-1}$) | $k_{\text{curved}}$ ($\text{km}^{-1}$) | $\lvert k_{\text{int}} - k_{\text{curved}}\rvert$ ($\text{km}^{-1}$) | $\phi'_{L,\text{int,Robin}}(R_*)$ ($\text{km}^{-1}$) | $\phi'_{L,\text{ext,curved}}(R_*)$ ($\text{km}^{-1}$) | $\lvert\phi'_{L,\text{int,Robin}} - \phi'_{L,\text{ext,curved}}\rvert$ ($\text{km}^{-1}$) | $\lvert\phi'_{L,\text{int,FD2}} - \phi'_{L,\text{ext,curved}}\rvert$ ($\text{km}^{-1}$) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Korai 4-lépéses (`it in {1,3,6,10}`, `tol_iter=2e-8`) | `19` | `0.232931108388` | `0.232931461902` | `3.535e-07` (`1.52e-6` rel.) | `-0.000965074633` | `-0.000965076098` | `1.465e-09` (`1.52e-6` rel.) | `2.268e-08` |
| Minden lépéses 4-komp. (`it >= 1`, `tol_iter=1e-11`) | `36` | `0.232931462480` | `0.232931462481` | `9.502e-13` (`4.08e-12` rel.) | `-0.000965074204` | `-0.000965074204` | `3.937e-15` (ghost: `1.430e-14`) | `2.121e-08` |
| Minden lépéses 3-komp. (`it >= 1`, `tol_iter=1e-10`, `DOP853`) | `32` | `0.232931463922` | `0.232931462481` | `1.441e-09` (`6.19e-9` rel.) | `-0.000965074210` | `-0.000965074204` | `5.970e-12` (ghost: `5.984e-12`) | `2.121e-08` |

- **A minden lépéses 4-komponensű Robin-frissítés (`tol_iter = 1e-11`, $N = 601$) hatása az egyensúlyi csillagparaméterekre a 4-lépéses ütemezéshez képest:**
  - $R_*$: `11.662698418 km` $\to$ `11.662698212 km` ($\lvert\Delta R_*\rvert = 2{,}06\times 10^{-7}\;\text{km}$, relatív eltérés $1{,}77\times 10^{-8}$),
  - $M_{\text{surf}}$: `1.560635499 M_sun` $\to$ `1.560635444 M_sun` ($\lvert\Delta M_{\text{surf}}\rvert = 5{,}56\times 10^{-8}\;M_\odot$),
  - $M_{\text{ADM}}$: `1.562034909 M_sun` $\to$ `1.562034848 M_sun` ($\lvert\Delta M_{\text{ADM}}\rvert = 6{,}07\times 10^{-8}\;M_\odot$),
  - $\phi_L(0)$: `0.011576878` $\to$ `0.011576871` ($\lvert\Delta\phi_L(0)\rvert = 6{,}80\times 10^{-9}$),
  - $\phi_L(R_*)$: `0.004143176` $\to$ `0.004143168` ($\lvert\Delta\phi_L(R_*)\rvert = 8{,}11\times 10^{-9}$).
  Az összes táblázatban közölt 6–7 értékes jegy változatlan marad, miközben a kikényszerített Robin-peremderivált ($\phi'_{L,\text{int,Robin}}(R_*)$) és a külső $\phi'_{L,\text{ext,curved}}(R_*)$ eltérése $1{,}465\times 10^{-9}\;\text{km}^{-1}$-ről $3{,}937\times 10^{-15}\;\text{km}^{-1}$-re csökken.

#### Külső fixpont-iterációs konvergenciatörténet (`it >= 1`, 4-komponensű külső Robin-frissítés, `tol_iter = 1e-11`, $N = 601$):

| Iteráció (`it`) | $R_*$ ($\text{km}$) | $M_{\text{surf}}$ ($M_\odot$) | $k_{\text{int}}$ ($\text{km}^{-1}$) | $\phi_L(0)$ | $\phi_L(R_*)$ | $\max\lvert\Delta\phi_L\rvert$ | $\lvert\Delta R_*\rvert$ ($\text{km}$) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| `1` | `11.671645529` | `1.454730352` | `0.228661436019` | `0.0056771316` | `0.0020193782` | `5.677e-03` | `1.167e+01` |
| `5` | `11.662535205` | `1.552333804` | `0.232578283400` | `0.0111240559` | `0.0039784369` | `4.545e-04` | `5.441e-05` |
| `10` | `11.662694098` | `1.560260725` | `0.232915538746` | `0.0115564555` | `0.0041357379` | `1.958e-05` | `3.311e-07` |
| `15` | `11.662698043` | `1.560618534` | `0.232930743933` | `0.0115759500` | `0.0041428329` | `8.832e-07` | `2.656e-08` |
| `20` | `11.662698204` | `1.560634681` | `0.232931430060` | `0.0115768296` | `0.0041431530` | `3.985e-08` | `1.224e-09` |
| `25` | `11.662698212` | `1.560635410` | `0.232931461019` | `0.0115768693` | `0.0041431675` | `1.798e-09` | `5.524e-11` |
| `30` | `11.662698212` | `1.560635442` | `0.232931462415` | `0.0115768710` | `0.0041431681` | `8.111e-11` | `2.495e-12` |
| `36` | `11.662698212` | `1.560635444` | `0.232931462480` | `0.0115768711` | `0.0041431681` | `8.742e-12` | `2.665e-13` |

### 2.3. Hálófinomítási próba ($N = 601$ vs. $N = 2000$ vs. $N = 4000$): $R_*$, $M_{\text{ADM}}$, $\phi_L(R_*)$ és az egyoldali FD2 felszíni derivált $\phi'_{L,\text{int,FD2}}(R_*)$ konvergenciája

A `gw170817_and_lensing_audit.py` `run_variant()` függvényében a `mesh_601`, `mesh_2000` és `mesh_4000` hálófinomítási futások beállításai a kód szerint: **`DOP853` integrátor, `rtol_ivp = 1e-11`, `atol = [1e-13, 1e-13, 1e-13]`, `max_step = 0.02 km`, `tol_iter = 1e-10` fixpont-küszöb**, valamint minden lépésben frissített görbült Robin-perem (`update_every = True, use_4comp_robin = False`; a `tol_iter = 1e-11` a 2.2. szakaszbeli külön 4-komponensű $N = 601$ futás küszöbe):

| Mennyiség | $N = 601$ ($h = 0{,}019438\;\text{km}$) | $N = 2000$ ($h = 0{,}005834\;\text{km}$) | $N = 4000$ ($h = 0{,}002916\;\text{km}$) | $\lvert\Delta\rvert_{601\to 2000}$ | $\lvert\Delta\rvert_{2000\to 4000}$ (absz. / rel.) |
|---|---:|---:|---:|---:|---:|
| $R_*$ ($\text{km}$) | `11.662698212` | `11.662697490` | `11.662697448` | `7.224e-07` | **`4.213e-08`** (`3.612e-09`) |
| $M_{\text{surf}}$ ($M_\odot$) | `1.560635444` | `1.560635270` | `1.560635260` | `1.741e-07` | **`9.670e-09`** (`6.196e-09`) |
| $M_{\text{ADM}}$ ($M_\odot$) | `1.562034848` | `1.562034651` | `1.562034640` | `1.976e-07` | **`1.097e-08`** (`7.025e-09`) |
| $\phi_L(0)$ | `0.0115768711` | `0.0115768355` | `0.0115768334` | `3.568e-08` | **`2.042e-09`** (`1.764e-07`) |
| $\phi_L(R_*)$ | `0.0041431681` | `0.0041431342` | `0.0041431323` | `3.394e-08` | **`1.920e-09`** (`4.633e-07`) |
| Kikényszerített $\phi'_{L,\text{int,Robin}}(R_*)$ ($\text{km}^{-1}$) | `-0.000965074210` | `-0.000965066327` | `-0.000965065881` | `7.883e-09` | **`4.463e-10`** (`4.625e-07`) |
| Független egyoldali $\phi'_{L,\text{int,FD2}}(R_*)$ ($\text{km}^{-1}$) | `-0.000965052991` | `-0.000965064416` | `-0.000965065402` | `1.143e-08` | **`9.865e-10`** (`1.022e-06`) |
| Külső $\phi'_{L,\text{ext,curved}}(R_*)$ ($\text{km}^{-1}$) | `-0.000965074204` | `-0.000965066321` | `-0.000965065874` | `7.883e-09` | **`4.463e-10`** (`4.625e-07`) |
| $\lvert\phi'_{L,\text{int,FD2}}(R_*) - \phi'_{L,\text{ext}}(R_*)\rvert$ ($\text{km}^{-1}$) | `2.1213e-08` | `1.9048e-09` | `4.7228e-10` | Skálázás: `11.137` (várt: `11.100`) | **Skálázás: `4.0332`** (várt $\mathcal{O}(h^2)$: `4.0020`) |
| $\max_{r\in(0,R_*)}\lvert\text{res}_{\text{ODE,iter}}\rvert$ ($\text{km}^{-2}$) | `7.5677e-13` | `7.8442e-13` | `7.4301e-13` | — | — |
| $\max_{r\in(0,R_*)}\lvert\text{res}_{\text{ODE,banded}}\rvert$ ($\text{km}^{-2}$) | `1.2719e-14` | `1.4172e-13` | `6.5786e-13` | — | — |
| $\hat{\alpha}_{\text{Yuk,conf}}(b_{\text{graz}})$ ($\text{rad}$) | `1.3340236309` | `1.3340236025` | `1.3340236015` | `2.833e-08` | **`1.055e-09`** (`7.907e-10`) |

**Következtetés és hatókör:**
1. $N = 2000$ és $N = 4000$ között a csillagsugár $\lvert\Delta R_*\rvert = 4{,}213\times 10^{-8}\;\text{km}$ ($3{,}61\times 10^{-9}$ relatív), az ADM-tömeg $\lvert\Delta M_{\text{ADM}}\rvert = 1{,}097\times 10^{-8}\;M_\odot$ ($7{,}03\times 10^{-9}$ relatív), a felszíni skalármező pedig $\lvert\Delta\phi_L(R_*)\rvert = 1{,}920\times 10^{-9}$ ($4{,}63\times 10^{-7}$ relatív) különbséggel konvergált.
2. Míg a $3{,}937\times 10^{-15}\;\text{km}^{-1}$ egyezés (a 2.2. szakasz 4-komponensű `tol_iter = 1e-11` futásában) a **kikényszerített Robin-peremderiváltra** ($\phi'_{L,\text{int,Robin}}(R_*) = -k_{\text{int}}\phi_L(R_*)$) igaz, addig a peremfeltételtől elkülönített **3-pontos egyoldali másodrendű véges differenciás felszíni derivált** ($\phi'_{L,\text{int,FD2}}(R_*) = (3\phi_{N-1} - 4\phi_{N-2} + \phi_{N-3})/(2h)$) és a külső görbült BVP derivált ($\phi'_{L,\text{ext,curved}}(R_*)$) eltérése $2{,}1213\times 10^{-8}\;\text{km}^{-1}$ ($N = 601$) $\to$ $1{,}9048\times 10^{-9}\;\text{km}^{-1}$ ($N = 2000$) $\to$ **$4{,}7228\times 10^{-10}\;\text{km}^{-1}$ ($N = 4000$-nél**, $\approx 4{,}89\times 10^{-7}$ relatív a $9{,}6507\times 10^{-4}\;\text{km}^{-1}$ meredekséghez mérve) szerint csökken, és másodrendben konvergál: az $N = 2000 \to 4000$ empirikus javulási arány **`4.0332`**, ami a másodrendű $\mathcal{O}(h^2)$ hálóskálázás várt $(h_{2000}/h_{4000})^2 = (3999/1999)^2 = 4{,}0020$ értékével összhangban álló empirikus skálázási teszt.
3. **Hatókör:** A fenti célzott numerikus BVP-ellenőrzés az egyetlen vizsgált $L = 10\;\text{km}$, $\alpha_Y = 1$, SLy $\log_{10}\rho_{0c} = 15{,}0$ konform referencia-csillagra és annak statikus, lencsézési diagnosztikájára vonatkozik.

---

## 3. Csomagolás, beágyazott `scalar-tensor-yukawa-fit_audit_package.zip` és frissített SHA-256 hash-ek (`392/392 OK` és `20/20 OK`)

A `/workspace/fix/` csomagolók (`build_step23_25_pdfs.py`, `build_v26_en_audit_pdf.py`, `finalize_manifests_and_zips.py`) újrafuttatásával a nagy `scalar-tensor-yukawa-fit_V26_package.zip` (és `scalar-tensor-yukawa-fit_step1_25_academy_review.zip`) tartalmazza a 3 önálló PDF-et (`step23_25_gw170817_lensing_bullet_cmb_hu.pdf`, `step23_25_gw170817_lensing_bullet_cmb_en.pdf`, `V25_english_audit_package.pdf`) és a visszaállított beágyazott audit-ZIP-et (`scalar-tensor-yukawa-fit/scalar-tensor-yukawa-fit_audit_package.zip`, valamint `V25_english_translations.zip`, összesen 394 bejegyzés a `scalar-tensor-yukawa-fit/` alatt).

### 3.1. A kis modul `SHA256SUMS_step23_25.txt` (`20/20 OK`) hash-ei

```text
c40714093b83206eda4bbec047ba01d217788c595a39d6c73dd43e4beec7b5ce  step23_gw170817_and_lensing.md
3f723eb9f5c06d01f1d3f318c6dd2e2a07df956520b58a253cbb4f1a004a6d02  step23_gw170817_and_lensing_en.md
571036611966f5e1de88828d3d66b1fcd187cf73abe917f9766f98d1a3cbe8f6  step24_bullet_cluster.md
0c57b832dd1be029a891797a17a5bc4f6a96ed6b5b19d202e69abb4435a0e71d  step24_bullet_cluster_en.md
8a58496edc72d40865b94753d0ec688f75836972c73e452b087e9c81c946dda9  step25_cmb_and_cosmology.md
eb4e4cfcfa03890d03505b005c96442b56e92c250327e97ffb5340ab17a4bbc2  step25_cmb_and_cosmology_en.md
c758fadd8baf5ccce0df8f3b047647409463d4cee0b5a65f0e97f7a2198cca29  step23_25_gw170817_lensing_bullet_cmb_hu.pdf
85a1fb3e77b3bd4cd59102ab204346ee9285982e6fdc0423c80378ac5b4fce46  step23_25_gw170817_lensing_bullet_cmb_en.pdf
18435bb1f964c72c08ad333716db5ee7d227a52a668888c6158436ac6c3e5ddd  gw170817_and_lensing_audit.py
8caca2431b9f2cda1847e9d5851ccf4b1b9511327d665ba8101a0b8d58f98914  bullet_cluster_lensing_audit.py
34f98435de79e59d962f34e77fcf52ed91a2e42a57c70ae02fc55b214e04ceb3  cmb_and_cosmology_audit.py
e03014bc93cb6de1d85d82398d737cb31b1da176c90d8ea84b5c3c274901591c  gw170817_dispersion_and_cT.csv
daa67c9a98e57704ea185a1c4682b3a9df47e808494610e739f25e54717e9f46  lensing_deflection_and_ppn_slip.csv
95dda8d26c9c6d8bd8d8c1ca1d15abe7408be41f0f1ea6f1d1e2d7a6efe600d7  gw170817_and_lensing_summary.json
8c28105dda15834a0bac4a05d17436969a13187fe07d04800d15a346785e06d4  bullet_cluster_convergence_profile.csv
1c55be84757dc68f7d4aeecb48430bda1bf07443371fc896783302335f7b758e  bullet_cluster_peak_offsets_summary.csv
5c6453d1fd17c4933839adaf343a40a4475698d471a05d0fe8e740aee77a1841  bullet_cluster_audit_summary.json
0ec03152b731c02b3f6a51f279a5d8a5920e63b48fbc4a95eacc2d60d4836cb9  cosmological_background_models.csv
0f0d86f1139a982f867132d141211d71ee3ee1db78134eaf058580bd55c8c9e2  cmb_geff_and_acoustic_scales.csv
8aeb25850578a519b94eedcd0983be1360eefce7acc71f3f8a012101b8cb7c06  cmb_and_cosmology_summary.json
```

### 3.2. A `SHA256SUMS.txt` (`392/392 OK`) és a csomagfájlok főbb hash-ei

| Fájl | Méret (bájt) | SHA-256 | Megjegyzés |
|---|---:|---|---|
| `SHA256SUMS.txt` | `40 033` | `fdf6692bc65e1ea89e54b41a476f40ea1621f8a75b0d1233ec10cf3d870cf7d2` | `392/392 OK` |
| `SHA256SUMS_step23_25.txt` | `1 989` | `febb1bfe434b22cc54a4f0bdfafb62efa24ad361b96cda05cf9290eafd80ec7c` | `20/20 OK` |
| `MANIFEST.txt` | `45 493` | `841731294c768ed29bc8765510b1109e712edcf0bc5fb3c8fa655eb870106735` | `391/391` sor |
| `MANIFEST_en.txt` | `45 077` | `235300563235319e7dc796611e03819d98942a83a6f75a5c02ca4ba075283dc6` | `388/388` sor |
| `README.md` | `119 011` | `a756186163b28b48654b7d5d28eb9e8fb7fc9cb7dff920a1ff6d8421ca23d498` | 25 HU + 25 EN fejezettábla |
| `README_en.md` | `97 492` | `82d55c6e8f2351451afe7dcd9941f5414e21bbfaea6ebd4d5a7fa0b800b3d82b` | 25 HU + 25 EN fejezettábla |
| `step23_gw170817_and_lensing.md` | `54 254` | `c40714093b83206eda4bbec047ba01d217788c595a39d6c73dd43e4beec7b5ce` | 2B. táblázat + 2B.1. szakasz |
| `step23_gw170817_and_lensing_en.md` | `52 257` | `3f723eb9f5c06d01f1d3f318c6dd2e2a07df956520b58a253cbb4f1a004a6d02` | Table 2B + Section 2B.1 |
| `gw170817_and_lensing_audit.py` | `70 431` | `18435bb1f964c72c08ad333716db5ee7d227a52a668888c6158436ac6c3e5ddd` | `max_step=0.02`, `tol_iter=1e-10` |
| `step23_25_gw170817_lensing_bullet_cmb_hu.pdf` | `890 347` | `c758fadd8baf5ccce0df8f3b047647409463d4cee0b5a65f0e97f7a2198cca29` | 27 oldal |
| `step23_25_gw170817_lensing_bullet_cmb_en.pdf` | `883 156` | `85a1fb3e77b3bd4cd59102ab204346ee9285982e6fdc0423c80378ac5b4fce46` | 27 oldal |
| `V26_english_audit_package.pdf` (`V25_english_audit_package.pdf`) | `1 798 601` | `1fe30c3ebef80be86bef8f7377ab3938ec739c536df981a0be2a81514cd79413` | 167 oldal |
| `scalar-tensor-yukawa-fit_audit_package.zip` | `6 168 154` | `3146431a0e07dc4508c8d370af149bf8645f1dfd027ba1805b0b5bbfbfc215a6` | 393 fájl (beágyazva is a V26 ZIP-ben) |
| `scalar-tensor-yukawa-fit_V26_package.zip` | `12 274 120` | `cd2aa4332a9e18b3d3cb6a1b46bf1325f16da2b9d9b33bad765568ad8c3a88fe` | 394 fájl (3 önálló PDF + 1 beágyazott audit-ZIP) |
| `scalar-tensor-yukawa-fit_step1_25_academy_review.zip` | `12 274 120` | `cd2aa4332a9e18b3d3cb6a1b46bf1325f16da2b9d9b33bad765568ad8c3a88fe` | 394 fájl (bájtonként azonos a V26 ZIP-pel) |
| `step23_25_gw170817_lensing_bullet_cmb_module.zip` | `1 911 725` | `9f8704631aca15db62f6904f6817a2212dc0f63e25c09be895306a0c21167659` | 21 fájl (`20 + SHA256SUMS_step23_25.txt`) |
| `V26_english_translations.zip` (`V25_english_translations.zip`) | `213 335` | `fee352698a64f6f216e577137c464ec6d477ab1ea73a80ebe3ed2c59525b2186` | 27 angol fájl |
