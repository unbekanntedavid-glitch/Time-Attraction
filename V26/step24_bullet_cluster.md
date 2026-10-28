# 24. lépés — Bullet-halmaz (`1E 0657-558`) gravitációs lencsézési konvergencia $\kappa(\mathbf{\theta})$ és gáz–galaxis csúcs-eltolódási (peak offset) analízis $L = 10\;\text{km}$, $\alpha_Y = 1$ mellett (`step24`, V26)

## I. Fizikai és módszertani keret: a Bullet-halmaz (`1E 0657-558`) ütközési tengely menti szintetikus 1D szemléltető tesztje

A `step23_gw170817_and_lensing.md` fejezetben igazoltuk, hogy a kanonikus V24 (`10.5281/zenodo.23680715`) és V25 (`10.5281/zenodo.23699685`) rövid hatótávolságú skalár–tenzor Yukawa-modellben ($L = 10\;\text{km} = 10^4\;\text{m}$, $\alpha_Y = 1$, $\beta_{\text{Yuk}} = 2\tilde{\alpha}_L\alpha_Y^2 = \alpha_Y^3/(4\pi) = 1/(4\pi) \approx 0.0795775$) a gravitációs hullámok tenzori terjedési sebessége a tisztán konformális minimális ágon egzaktul $c_T = c$ ($\Delta c_T/c = 0$), továbbá a konformális Jordan-metrika ($\tilde{g}_{\mu\nu} = A^2(\phi_L)g_{\mu\nu} = e^{-2\alpha_Y\phi_L}g_{\mu\nu}$) alatt a lineáris Weyl-féle lencsézési potenciál $\Phi_{\text{lens}} = \frac{1}{2}(\tilde{\Psi} + \tilde{\Phi}) = \Phi_N$ alakot ölt. Jelen 24. lépésben (`step24`) ezt a keretrendszert alkalmazzuk a **Bullet-halmaz (`1E 0657-558`)** (Markevitch et al. 2004, *ApJ* 606, 819; Clowe et al. 2006, *ApJ* 648, L109; Bradač et al. 2006, *ApJ* 652, 937) ütközési tengely menti **szintetikus 1D szemléltető modelljére**.

### 1.1. Kozmológiai és lencsézési alapadatok (`1E 0657-558`)
- **Lencse vöröseltolódása:** $z_L = 0.296$, anguláris átmérőtávolság $D_L = 910.0\;\text{Mpc}$ ($1' \leftrightarrow 264.7\;\text{kpc}$).
- **Háttérforrások jellemző vöröseltolódása:** $z_S \approx 0.8\text{–}1.0$ (karakterisztikus érték: $z_S = 0.85$).
- **Kritikus felületi tömegsűrűség:**
  $$\Sigma_{\text{crit}} = \frac{c^2 D_S}{4\pi G D_L D_{LS}} = 3.10\times 10^9\;M_\odot\,\text{kpc}^{-2} = 0.6474\;\text{g}\,\text{cm}^{-2}.$$
- **Ütközési sebesség (lökéshullám-front):** $v_{\text{coll}} \approx 4500\;\text{km/s} \approx 0.01501\,c$, amely szuperszonikus torlónyomásos (ram-pressure) leválasztást idéz elő az ütköző röntgengáz (ICM) és az ütközésmentes (collisionless) galaxisok között.

### 1.2. Négykomponensű barionos modell az ütközési tengely mentén ($x \in [-800, +800]\;\text{kpc}$)
A Chandra röntgentérképek és a Hubble/VLT/Magellan optikai megfigyelések karakterisztikus tömeg- és méretarányai alapján a Bullet-halmaz barionos tömegeloszlását az ütközési tengely ($y = 0$, $x \in [-800, +800]\;\text{kpc}$) mentén négy analitikus komponenssel szemléltetjük, rekesz-normált ($R_{\text{ap}} = 650.0\;\text{kpc}$) King / módosított Hubble ($\beta = 1$) profillal $\Sigma_i(R) = \Sigma_{0,i}(1 + R^2/R_{c,i}^2)^{-1}$, ahol $\Sigma_{0,i} = \frac{M_i}{\pi R_{c,i}^2 \ln(1 + R_{\text{ap}}^2/R_{c,i}^2)}$:

| Komponens azonosító | Fizikai leírás | Centroid $x_{c,i}$ ($\text{kpc}$) | Rekesz-tömeg $M_i$ ($M_\odot$) | Magrádiusz $R_{c,i}$ ($\text{kpc}$) | Centrális $\Sigma_{0,i}$ ($M_\odot/\text{kpc}^2$) | Centrális $\kappa_{0,i} = \Sigma_{0,i}/\Sigma_{\text{crit}}$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `gal_main` | Főhalmaz ütközésmentes galaxisok (csillagok + kompakt maradványok) | `-350.0` | `5.50e12` | `120.0` | `3.9329e7` | `0.012687` |
| `gas_main` | Főhalmaz ütköző röntgen ICM plazma (Chandra, torlónyomással lefékezve) | `-170.0` | `2.60e13` | `150.0` | `1.3638e8` | `0.043995` |
| `gas_sub` | Alhalmaz („Bullet”) lökéshullám-tömörített röntgen ICM gáz | `+140.0` | `6.50e12` | `75.0` | `8.4628e7` | `0.027300` |
| `gal_sub` | Alhalmaz („Bullet”) ütközésmentes galaxisok ($160\;\text{kpc}$-cel a lökéshullám előtt) | `+300.0` | `1.80e12` | `65.0` | `2.9296e7` | `0.009450` |

- **Összesített barionos tömegmérleg a rekeszben:**
  $$M_{\text{gas,tot}} = 3.25\times 10^{13}\;M_\odot, \qquad M_{\text{gal,tot}} = 7.30\times 10^{12}\;M_\odot, \qquad M_{\text{bary,tot}} = 3.98\times 10^{13}\;M_\odot, \qquad \frac{M_{\text{gas,tot}}}{M_{\text{gal,tot}}} = 4.4521.$$
- **Megfigyelések által motivált csúcs-helyzetek és szintetikus 1D referencia-görbe (`kappa_synth_target_1D`):**
  A gyenge és erős lencsézési rekonstrukciókban (Clowe et al. 2006; Bradač et al. 2006) a konvergencia-csúcsok nem a domináns tömegű röntgengáz ($x = -170\;\text{kpc}$ és $+140\;\text{kpc}$), hanem az ütközésmentes galaxis-centroidok környezetében helyezkednek el ($x_{\kappa,\text{main}}^{\text{obs}} \approx -350.0\;\text{kpc}$, $\kappa_{\text{peak,main}}^{\text{obs}} \approx 0.320$; $x_{\kappa,\text{sub}}^{\text{obs}} \approx +300.0\;\text{kpc}$, $\kappa_{\text{peak,sub}}^{\text{obs}} \approx 0.220$).
  **Módszertani elhatárolás (V26):** A `bullet_cluster_convergence_profile.csv` táblázatban szereplő `kappa_synth_target_1D` profil **nem közvetlen megfigyelési adatsor**, hanem egy **önkonstruált 1D szintetikus referencia-görbe**:
  $$\kappa_{\text{synth}}(x) = \kappa_{\text{GR,bary}}(x) + \frac{w_{\text{synth,main}}\,\Sigma_{\text{gal,main}}(x) + w_{\text{synth,sub}}\,\Sigma_{\text{gal,sub}}(x)}{\Sigma_{\text{crit}}}, \qquad (w_{\text{synth,main}} = 22.9765,\; w_{\text{synth,sub}} = 21.3095),$$
  amelyet a modell saját négy King-profiljából állítottunk elő úgy, hogy az pontosan reprodukálja az $x = -350\;\text{kpc}$-nél ($\kappa = 0.320$) és $x = +300\;\text{kpc}$-nél ($\kappa = 0.220$) megadott szemléltető csúcsmagasságokat. Ennek megfelelően a `synth_mse_norm_vs_target` oszlopban szereplő
  $$M_{\text{synth}} \equiv \frac{1}{N}\sum_{i=1}^{N}\left[\frac{\kappa(x_i) - \kappa_{\text{synth}}(x_i)}{\sigma_{\text{ref}}}\right]^2 \qquad (\sigma_{\text{ref}} = 0.020,\; N = 321)$$
  mennyiség **normált szintetikus átlagos négyzetes eltérés** a `kappa_synth_target_1D` referencia-görbétől, **nem pedig** a Clowe et al. (2006) vagy Bradač et al. (2006) megfigyelési adatokra végzett valódi statisztikai $\chi^2_\nu$ illesztés. Valódi megfigyelési $\chi^2$ illesztéshez publikált 2D konvergenciatérkép, hibakovariancia-mátrix, előre rögzített térbeli mintavételezés és a becsült szabadsági fokok ($N_{\text{dof}} = N_{\text{bins}} - p$) levonása szükséges.

---

## II. 1. Teszt — Tiszta lokális egy-skálás Yukawa-modell rögzített $L = 10\;\text{km}$, $\alpha_Y = 1$ mellett (Szigorú klaszter-skálájú falszifikáció / No-Go eredmény)

Vizsgáljuk meg egzaktul, milyen 3D-s skalármezőt $\phi_L^{\text{local}}(\mathbf{r})$ és milyen vetített 2D-s lencsézési konvergenciát $\kappa_{\text{local}}(x)$ állít elő a V24/V25 lokális Yukawa-modell rögzített $L = 10\;\text{km} = 3.24077929\times 10^{-16}\;\text{kpc}$ és $\alpha_Y = 1$ paraméterek mellett!

### 2.1. Algebrai kontakt-határérték a klaszter-skálán ($R_c \sim 65\text{–}150\;\text{kpc} \sim (2.0\text{–}4.6)\times 10^{17}\,L$)
A V23–V25 konvenció szerinti statikus 3D-s Helmholtz–Yukawa-egyenlet ($\phi_L > 0$ kompakt források környezetében):
$$\left(\nabla^2 - \frac{1}{L^2}\right)\phi_L(\mathbf{r}) = -\alpha_{\text{std}}\,\frac{4\pi G}{c^2}\,\rho_b(\mathbf{r}), \qquad \alpha_{\text{std}} = \frac{\alpha_Y}{4\pi} \approx 0.0795775.$$
Mivel a halmaz magrádiuszai $R_{c,i} \in [65, 150]\;\text{kpc}$, a dimenziótlan arány:
$$\frac{R_{c,i}}{L} \in \left[2.006\times 10^{17},\; 4.628\times 10^{17}\right] \implies \left(\frac{L}{R_{c,i}}\right)^2 \in \left[4.668\times 10^{-36},\; 2.486\times 10^{-35}\right].$$
A térbeli deriváltak ($|\nabla^2\phi_L| \sim R_c^{-2}|\phi_L|$) több mint 35 nagyságrenddel kisebbek a $-L^{-2}\phi_L$ tömegtagnál, így a 3D-s skalármező a lokális algebrai kontakt-határértékre redukálódik:
$$\left|\phi_L^{\text{local}}(\mathbf{r})\right| \simeq \alpha_{\text{std}}\,L^2\,\frac{4\pi G}{c^2}\,\rho_b(\mathbf{r}) \in \left[10^{-44},\; 3.086\times 10^{-42}\right].$$
Ami a galaxisokat alkotó egyedi csillagokat és kompakt maradványokat illeti: minden egyes csillag $L = 10\;\text{km}$-es Yukawa-burka $e^{-r/(10\;\text{km})}$ szerint cseng le, így $r > 100\;\text{km}$ távolságon ($< 1.5\times 10^{-4}\,R_\odot$, illetve $< 3\times 10^{-18}\;\text{kpc}$) már $e^{-10} \approx 4.5\times 10^{-5}$ alá, halmaz-skálán pedig $\exp(-10^{17}) = 0$ értékre tűnik el.

### 2.2. Konformális vs. diszformális lencsézési konvergencia $L = 10\;\text{km}$ mellett
1. **Konformális Jordan-ág (`step14` Branch II-Conf, `step23`):**
   $\tilde{g}_{\mu\nu} = A^2(\phi_L)g_{\mu\nu} = e^{-2\alpha_Y\phi_L}g_{\mu\nu}$ esetén a lineáris Weyl-potenciálból a skalármező egzaktul kiesik ($\Phi_{\text{lens}}^{\text{conf}} = \Phi_N$), a kompakt objektumok külső skalár-energiaburka ($\Delta M_{\phi,\text{ext}}^{\text{curved}}/M_{\text{surf}} = +1.55\times 10^{-3}$ neutroncsillagokra, és $\lesssim 10^{-16}$ fősorozati csillagokra) pedig már beépül a csillagok aszimptotikus ADM-tömegébe ($M_{\text{gal}}$). Így:
   $$\Delta\kappa_{\text{local }\,L=10\text{ km}}^{\text{conf}}(x) = 0 \quad (\text{azonosan}).$$
2. **Diszformális refraktív óra-ág (`step14` Branch II-Disf):**
   Ha egy külön diszformális optikai ágon a lencsézési potenciálra $\Phi_{\text{lens}}^{\text{disf}} = \Phi_N - \frac{1}{2}c^2\phi_L$ teljesül, a látóirányú ($z \in (-\infty, +\infty)$) integrálás és a 2D transverszális Laplace-operátor alkalmazása után a 2D Fourier-térbeli szűrő $\frac{k_\perp^2}{k_\perp^2 + L^{-2}} \simeq L^2 k_\perp^2 \leftrightarrow -L^2\nabla_\perp^2$, vagyis:
   $$\Delta\kappa_{\text{local }\,L=10\text{ km}}^{\text{disf}}(x) = -\frac{\alpha_{\text{std}}}{2}\,L^2\,\nabla_\perp^2 \kappa_{\text{GR,bary}}(x), \qquad \max_{x}\left|\frac{\Delta\kappa_{\text{local}}^{\text{disf}}}{\kappa_{\text{GR,bary}}}\right| = 2.171\times 10^{-36} \quad (\text{a külső halóban } < 10^{-40}).$$

### 2.3. Transzparens tudományos verdikt az 1. Tesztre
Mindkét csatolási ágon:
$$\kappa_{\text{local }\,L=10\text{ km}}(x) = \kappa_{\text{GR,bary}}(x) = \frac{\Sigma_{\text{gas}}(x) + \Sigma_{\text{gal}}(x)}{\Sigma_{\text{crit}}}.$$
- A konvergencia-csúcsok a **röntgengáz centroidjainál** maradnak: a 4-komponensű szimmetrikus szuperpozícióban $x_{\text{peak,main}}^{\text{local}} = -175.44\;\text{kpc}$ ($\kappa_{\text{peak,main}}^{\text{local}} = 0.04503$) és $x_{\text{peak,sub}}^{\text{local}} = +137.01\;\text{kpc}$ ($\kappa_{\text{peak,sub}}^{\text{local}} = 0.03697$), illetve a lökéshullám-front mentén szeparált / meredekebb külső gázprofilú alhalmaz-párokban $x_{\text{peak,main}}^{\text{local}} = -170.7\;\text{kpc}$ és $x_{\text{peak,sub}}^{\text{local}} = +142.3\;\text{kpc}$ ($\kappa_{\text{peak,main}}^{\text{local}} = 0.0450$).
- A főcsúcs amplitúdója **$7.11$-szer kisebb** a $\kappa = 0.320$ referencia-értéknél, a csúcsok pozíciója pedig **$157.7\text{–}179.3\;\text{kpc}$-cel elcsúszik** a galaxis-centroidoktól ($-350\;\text{kpc}$ és $+300\;\text{kpc}$) az ütköző röntgengáz irányába ($M_{\text{synth}} = 32.089$).
- **Következtetés:** **Rögzített $L = 10\;\text{km}$ és $\alpha_Y = 1$ mellett a rövid hatótávolságú Yukawa-féle ötödik erő tisztán kompakt-objektum / neutroncsillag-léptékű effektus, amely zérus makroszkopikus módosítást ad a Bullet-halmaz lencsézéséhez ($\Delta\kappa/\kappa_{\text{bary}} \le 2.17\times 10^{-36}$; a külső tartományban $< 10^{-40}$).**

---

## III. 2. Teszt — Fenomenológiai makroszkopikus kiterjesztések a szintetikus 1D teszten: mi szükséges ahhoz, hogy a $\kappa(x)$ csúcs a röntgengázról a galaxisokra kerüljön?

Vizsgáljuk meg kvantitatívan a szintetikus 1D King-modellben (`step11`, `step14` Branch II-Disf: $\Phi_{\text{lens}} = \Phi_N - \frac{1}{2}c^2\phi_{\text{macro}}$, szemléltető klaszter-skálájú Helmholtz-hossz $L_{\text{macro}} = 300.0\;\text{kpc}$) a két makroszkopikus forrásolási esetet!

### 3.1. 2A eset (`Macro_unscreened_bulk_baryon`) — Statikus makroszkopikus kernel, amely kizárólag a kisimított összesített barionsűrűségtől ($\bar{\rho}_b = \rho_{\text{gas}} + \rho_{\text{gal}}$) függ
Tegyük fel, hogy egy makroszkopikus skalármező egyetemes $\eta_{\text{bulk}}$ csatolással forrásolódik az összesített barionos tömegsűrűségből:
$$\left(\nabla^2 - L_{\text{macro}}^{-2}\right)\phi_{\text{macro}}(\mathbf{r}) = -\frac{4\pi G}{c^2}\,\eta_{\text{bulk}}\left[\rho_{\text{gas}}(\mathbf{r}) + \rho_{\text{gal}}(\mathbf{r})\right], \qquad \widetilde{\Sigma}_{\phi,\text{bulk}}(\mathbf{k}_\perp) = \eta_{\text{bulk}}\,\frac{k_\perp^2}{k_\perp^2 + L_{\text{macro}}^{-2}}\,\widetilde{\Sigma}_b(\mathbf{k}_\perp).$$
Ha $\eta_{\text{bulk}} = 8.6523$ értéket választunk, a főcsúcs amplitúdója $7.11$-szeresére nő, így eléri a $\kappa_{\text{peak,main}} = 0.3200$ értéket. **Ámde:** mivel $M_{\text{gas,tot}}/M_{\text{gal,tot}} = 4.4521$, a konvergencia-csúcsok **továbbra is a röntgengázhoz maradnak rögzítve**:
$$x_{\text{peak,main}}^{\text{bulk}} = -176.02\;\text{kpc} \quad (\text{eltolódás a galaxisoktól: } +173.98\;\text{kpc}), \qquad x_{\text{peak,sub}}^{\text{bulk}} = +138.35\;\text{kpc} \quad (\text{eltolódás a galaxisoktól: } -161.65\;\text{kpc}),$$
($M_{\text{synth}} = 20.723$). Ez igazolja, hogy **semmilyen olyan statikus módosított gravitációs vagy skalár–tenzor kernel, amelynek forrása egyszerűen arányos a lokális makroszkopikus barionsűrűséggel ($\rho_{\text{gas}} + \rho_{\text{gal}}$), nem képes a konvergencia-csúcsokat a domináns tömegű röntgengázról az ütközésmentes galaxisokra áthelyezni.**

### 3.2. 2B eset (`Macro_collisionless_compact_clock`) — Fenomenológiai ütközésmentes / kompaktság-súlyozott forrásolási hipotézis
Ahhoz, hogy hideg sötét anyag (CDM) részecskék nélkül egy makroszkopikus mező lencsézési csúcsa az ütközésmentes galaxisokat kövesse, a makroszkopikus szektornak az alábbi két fizikai tulajdonság valamelyikével kellene rendelkeznie:
1. **2B-1. hipotézis — Kompaktsági / mély lokális potenciálgödör szerinti súlyozás:** A diffúz ICM plazmában ($n_e \sim 10^{-3}\;\text{cm}^{-3}$) nincsenek csillag-léptékű potenciálgödrök ($\langle (\nabla\nu)^2 \rangle_{\text{micro}} \approx 0$), míg a galaxis-komponensben a barionos tömeg csillagok és kompakt maradványok potenciálgödreibe van bezárva. Ha egy makroszkopikus mező effektív forráserőssége a kis-skálás gradiens-energiával vagy a lokális kompaktsággal súlyozódik, akkor a forrás dominánsan az ütközésmentes ($\sigma_{\text{coll}} \approx 0$) csillag- és kompakt-objektum populációhoz kötődik.
2. **2B-2. hipotézis — Nem-egyensúlyi / retardált skalár-memória:** Ha a mag-áthaladási időskálánál ($t_{\text{cross}} = \Delta x / v_{\text{coll}} \approx 180\;\text{kpc} / (4500\;\text{km/s}) \approx 39.1\;\text{Myr}$) hosszabb a makroszkopikus mező relaxációs ideje ($\tau_{\text{relax}} \gtrsim 150\;\text{Myr} > t_{\text{cross}}$), a skalár-haló tehetetlenül együtt haladhat az ütközésmentes galaxisokkal.

### 3.3. Modellfüggő küszöbérték-számítás az $\eta_{\text{gal}} / \eta_{\text{gas}}$ arányra és módszertani korlátok
Ha a választott 1D King-geometriában ($R_{c,i} \in [65, 150]\;\text{kpc}$) és $L_{\text{macro}} = 300.0\;\text{kpc}$ Helmholtz-szűrés mellett a galaxis- és gázkomponensek fenomenológiai csatolási súlyát $\eta_{\text{gal}}$ és $\eta_{\text{gas}}$ jelöli:
1. **Csúcs-dominancia küszöb ($\kappa(-350\;\text{kpc}) \ge \kappa(-170\;\text{kpc})$ és $\kappa(+300\;\text{kpc}) \ge \kappa(+140\;\text{kpc})$):**
   - Közvetlen felületi sűrűség-arány (Helmholtz-szűrés nélkül):
     $$\frac{\Sigma_{\text{gas,main}}(-170)}{\Sigma_{\text{gal,main}}(-350)} = 3.4592, \qquad \frac{\Sigma_{\text{gas,sub}}(+140)}{\Sigma_{\text{gal,sub}}(+300)} = 2.8895.$$
   - $L_{\text{macro}} = 300\;\text{kpc}$ Helmholtz–Yukawa szűréssel ($\eta_{\text{gas}} = 1$ egységnyi gáz-csatolás mellett):
     $$\left(\frac{\eta_{\text{gal}}}{\eta_{\text{gas}}}\right)_{\text{dom,main}} = 5.1204, \qquad \left(\frac{\eta_{\text{gal}}}{\eta_{\text{gas}}}\right)_{\text{dom,sub}} = 5.0972.$$
2. **$< 5.0\;\text{kpc}$ centroid-eltolódási küszöb ($|x_{\text{peak,main}} - (-350\;\text{kpc})| \le 5.0\;\text{kpc}$):**
   - Ahhoz, hogy ebben az 1D King-geometriában a röntgengáz csóvája $5.0\;\text{kpc}$-nél kisebb mértékben húzza el a főhalmaz lencsézési csúcsát a $-350.0\;\text{kpc}$ galaxis-centroidtól:
     $$\left(\frac{\eta_{\text{gal}}}{\eta_{\text{gas}}}\right)_{|\Delta x| \le 5\;\text{kpc}} \ge 25.6650.$$
3. **Kalibrált 2B szemléltető modell (`Macro_collisionless_compact_clock`, $\eta_{\text{gas}} = 0$, $\eta_{\text{gal,main}} = 31.7448$, $\eta_{\text{gal,sub}} = 24.1556$):**
   - Főhalmaz csúcs: $x_{\text{peak,main}} = -347.86\;\text{kpc}$ (eltolódás a galaxis-centroidtól: **$+2.14\;\text{kpc} < 5\;\text{kpc}$**), $\kappa_{\text{peak,main}} = 0.32011$ ($\kappa(-350\;\text{kpc}) = 0.32000$).
   - Alhalmaz csúcs: $x_{\text{peak,sub}} = +299.45\;\text{kpc}$ (eltolódás a galaxis-centroidtól: **$-0.55\;\text{kpc} < 5\;\text{kpc}$**), $\kappa_{\text{peak,sub}} = 0.22002$ ($\kappa(+300\;\text{kpc}) = 0.22000$).
   - Eltérés a szintetikus 1D referencia-görbétől (`kappa_synth_target_1D`): $\text{RMS}(\Delta\kappa) = 0.02559$, normált szintetikus négyzetes eltérés $M_{\text{synth}} = 1.637$ (szemben a lokális $L = 10\;\text{km}$ modell $M_{\text{synth}} = 32.089$ és a 2A bulk-barion modell $M_{\text{synth}} = 20.723$ értékével).
4. **Szigorú elméleti és módszertani elhatárolás (V26):**
   - Az $\eta_{\text{gal}}/\eta_{\text{gas}} \ge 5.1204$ dominancia-küszöb és az $\eta_{\text{gal}}/\eta_{\text{gas}} \ge 25.6650$ ($<5\;\text{kpc}$) centroid-küszöb **a választott 1D King-geometriára ($R_{c,i} = 65\text{–}150\;\text{kpc}$) és az $L_{\text{macro}} = 300\;\text{kpc}$ Helmholtz-skálára vonatkozó modellfüggő küszöbértékek**, nem pedig modellfüggetlen általános fizikai tételek.
   - Az $\eta_{\text{gas}} = 0$ előírás, valamint a két külön illesztett galaktikus súly ($\eta_{\text{gal,main}} = 31.7448 \neq \eta_{\text{gal,sub}} = 24.1556$) **fenomenológiai, szemléltető feltevések, amelyek nem vezethetők le a kovariáns skalár–tenzor hatásból**.

---

## IV. Szintetikus 1D szemléltető csúcs-eltolódási és konvergencia-táblázatok

### 4.1. A négy modell csúcs-pozícióinak, amplitúdóinak és a szintetikus 1D célprofiltól való eltéréseinek összefoglaló táblázata (`bullet_cluster_peak_offsets_summary.csv` — szintetikus 1D szemléltető teszt)

| Modell (`model_name`) | $x_{\text{peak,main}}$ ($\text{kpc}$) | $\kappa_{\text{peak,main}}$ | $x_{\text{peak,sub}}$ ($\text{kpc}$) | $\kappa_{\text{peak,sub}}$ | $\Delta x_{\text{main,gal}}$ ($\text{kpc}$) | $\Delta x_{\text{main,gas}}$ ($\text{kpc}$) | $\Delta x_{\text{sub,gal}}$ ($\text{kpc}$) | $\Delta x_{\text{sub,gas}}$ ($\text{kpc}$) | RMS $\Delta\kappa_{\text{synth}}$ | $M_{\text{synth}}$ (`synth_mse_norm_vs_target`) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GR_baryon_only` | `-175.44` | `0.04503` | `+137.01` | `0.03697` | `+174.56` | `  -5.44` | `-162.99` | `  -2.99` | `0.11329` | `32.089` |
| `TAL_local_Yukawa_L10km_alphaY1` | `-175.44` | `0.04503` | `+137.01` | `0.03697` | `+174.56` | `  -5.44` | `-162.99` | `  -2.99` | `0.11329` | `32.089` |
| `Macro_unscreened_bulk_baryon` | `-176.02` | `0.32002` | `+138.35` | `0.26531` | `+173.98` | `  -6.02` | `-161.65` | `  -1.65` | `0.09104` | `20.723` |
| `Macro_collisionless_compact_clock` | `-347.86` | `0.32011` | `+299.45` | `0.22002` | `  +2.14` | `-177.86` | `  -0.55` | `+159.45` | `0.02559` | `1.637` |

*(Megjegyzés: $M_{\text{synth}} = \frac{1}{N}\sum_{i=1}^N [(\kappa(x_i) - \kappa_{\text{synth}}(x_i))/\sigma_{\text{ref}}]^2$ az önkonstruált 1D `kappa_synth_target_1D` görbétől mért normált négyzetes eltérés $\sigma_{\text{ref}} = 0.020$ referencia-skálával, nem megfigyelési $\chi^2_\nu$.)*

### 4.2. Vetített felületi sűrűség és $\kappa(x)$ profil az ütközési tengely reprezentatív pontjaiban (`bullet_cluster_convergence_profile.csv` — szintetikus 1D szemléltető teszt)

| $x$ ($\text{kpc}$) | $\Sigma_{\text{gas}}$ ($M_\odot/\text{kpc}^2$) | $\Sigma_{\text{gal}}$ ($M_\odot/\text{kpc}^2$) | $\kappa_{\text{GR,bary}}$ | $\kappa_{\text{local }\,L=10\text{ km}}$ | $\Delta\kappa_{\text{local,disf}}$ | $\kappa_{\text{Macro,2A (bulk)}}$ | $\kappa_{\text{Macro,2B (clock)}}$ | $\kappa_{\text{synth,1D}}$ (`kappa_synth_target_1D`) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `-600.0` | `1.4233e+07` | `6.8239e+06` | `0.006793` | `0.006793` | `-5.740e-40` | `0.010277` | `0.028049` | `0.062017` |
| `-450.0` | `2.8832e+07` | `2.1245e+07` | `0.016154` | `0.016154` | `+2.251e-42` | `0.076874` | `0.165259` | `0.188492` |
| `-350.0` | `5.2453e+07` | `3.5918e+07` | `0.028507` | `0.028507` | `+1.240e-38` | `0.178665` | `0.320000` | `0.320000` |
| `-260.0` | `9.3503e+07` | `2.3192e+07` | `0.037644` | `0.037644` | `+8.951e-39` | `0.254340` | `0.203216` | `0.225535` |
| `-175.0` | `1.2766e+08` | `1.1934e+07` | `0.045031` | `0.045031` | `+2.869e-38` | `0.320000` | `0.104232` | `0.141178` |
| `-170.0` | `1.2794e+08` | `1.1514e+07` | `0.044984` | `0.044984` | `+2.882e-38` | `0.319487` | `0.100405` | `0.137693` |
| ` -80.0` | `9.9460e+07` | `6.7119e+06` | `0.034249` | `0.034249` | `+6.544e-39` | `0.219456` | `0.048876` | `0.087474` |
| `  +0.0` | `7.2882e+07` | `5.0651e+06` | `0.025144` | `0.025144` | `-3.132e-39` | `0.138610` | `0.027722` | `0.064214` |
| ` +70.0` | `7.9995e+07` | `4.8621e+06` | `0.027373` | `0.027373` | `+3.408e-40` | `0.166192` | `0.029235` | `0.063429` |
| `+140.0` | `1.0829e+08` | `6.1784e+06` | `0.036924` | `0.036924` | `+8.015e-38` | `0.265185` | `0.049209` | `0.080507` |
| `+220.0` | `5.5598e+07` | `1.3197e+07` | `0.022192` | `0.022192` | `-2.748e-39` | `0.134507` | `0.088322` | `0.110826` |
| `+300.0` | `2.6688e+07` | `3.0558e+07` | `0.018466` | `0.018466` | `+3.559e-38` | `0.110815` | `0.220000` | `0.220000` |
| `+450.0` | `1.1510e+07` | `5.4294e+06` | `0.005464` | `0.005464` | `-9.488e-40` | `0.006652` | `0.021648` | `0.042188` |
| `+600.0` | `6.7045e+06` | `1.8771e+06` | `0.002768` | `0.002768` | `-1.602e-40` | `-0.004653` | `0.000158` | `0.015925` |

---

## V. Reprodukálhatóság és generált fájlok jegyzéke

- **Numerikus audit szkript:** `bullet_cluster_lensing_audit.py` (szinkronizálva a `code/` könyvtárba).
- **1D ütközési-tengely szintetikus konvergencia-profil (321 rácspont, $x \in [-800, +800]\;\text{kpc}$):** `bullet_cluster_convergence_profile.csv` (szinkronizálva a gyökér-, `code/`, `csv/` és `data/v23/` könyvtárakba).
- **Csúcs-eltolódási és normált szintetikus eltérési ($M_{\text{synth}}$) összefoglaló táblázat:** `bullet_cluster_peak_offsets_summary.csv` (szinkronizálva a gyökér-, `code/`, `csv/` és `data/v23/` könyvtárakba).
- **Strukturált JSON összefoglaló:** `bullet_cluster_audit_summary.json` (szinkronizálva a gyökér-, `code/`, `csv/` és `data/v23/` könyvtárakba).
