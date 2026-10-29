# 26. lépés (V27): Kozmológiai $H(z)$, DESI DR1 BAO, Planck 2018 akusztikus skála, szintetikus SNe Ia alakkontroll, variációs RSD $f\sigma_8$, algebrai BBN-benchmark és inverz távolságlétra no-go audit

**Szerző:** Szabolcsi Dávid (`ORCID: 0009-0005-6300-074X`)  
**Dátum:** 2026-10-28 (a szöveg az Academia Co-scientist bírálata, az azt követő újrabírálat és az F01–F12 második bírálat, valamint a harmadik újraaudit N01, N02, T01 pontjai után átdolgozva: 2026-10-29, lásd 11., 12. és 13. fejezet)  
**Előzmény DOI-k:**

- **V24:** Zenodo DOI `10.5281/zenodo.23680715` (szigorúan V24)
- **V25:** Zenodo DOI `10.5281/zenodo.23699685` (szigorúan V25, `step1`–`step22`)
- **V26:** Zenodo DOI `10.5281/zenodo.23727585` (szigorúan V26, `step1`–`step25`)

**A számok forrása:** a számítási eredményeket a JSON/CSV-forrásokból generáljuk (`tal_hubble_desi_bao_summary.json` és a 26. lépés CSV-táblázatai, a közös `tal_cosmology_common.py` háttérmodullal készítve) a `generate_step26_markdown.py` segítségével. A rögzített bemeneti állandók (adatértékek, priorok, irodalmi benchmarkok, külső bemenetek) a kódban szerepelnek, és a szöveg hivatkozik rájuk.

---

## 0. Hatókör és összefoglaló

A modul egyetlen konkrét ágat vizsgál: a V26 `Step 25` M4 hátteret rögzített $L=10\;\text{km}$, $\alpha_Y=1$, $\beta_{\rm Yuk}=1/(4\pi)=0.079577$ mellett, a skaláregyenlet adiabatikus partikuláris megoldásán, elhanyagolható szabad skaláros oszcillációs energiával. A szűkebb, védhető állítás: **ebben az ágban a $\beta_{\rm Yuk}$ csatolás puszta behelyettesítése $H_0$ szorzójaként nem levezetett mechanizmus a $67.4\to 72.8\;\text{km/s/Mpc}$ eltolásra.** Ez nem a TAL minden lehetséges kozmológiai megoldásának kizárása.

Az átdolgozás a bírálati pontokat így valósítja meg:

- egyetlen hatás, egyetlen skaláregyenlet (SI egységek, kozmikus másodperc), egyetlen attraktorformula és egyetlen $G_{\rm eff}$-normalizálás, a szövegben és a kódban azonosan (1. fejezet);
- explicit ágkiválasztás (adiabatikus attraktor, $C=0$ szabad oszcillációs amplitúdó) a kizárólagosság állítása helyett (1.3);
- DESI DR1 BAO a kerekítetlen Cobaya 12-elemű vektorral és $12\times 12$ kovarianciával, profilozott kontrollokkal (2. fejezet);
- minden CMB-modellhez és az M3-kontrollhoz modellspecifikus, újraintegrált $r_s(z_*)$ és $r_d(z_{\rm drag})$; a korábbi $\theta_*$- és $r_d$-alapú szigma-kizárási címkék törölve (3. és 8. fejezet);
- a szupernóva-rész új neve: **szintetikus távolságmodulus-alakkontroll** (4. fejezet);
- RSD: diagonális benchmark és elsőrendű variációs $\delta f\sigma_8$ (5. fejezet);
- BBN: csak algebrai nagyságrendi benchmark (6. fejezet);
- lokális környezetek: állandó sűrűségű bulk-benchmarkok `expm1`-gyel (7. fejezet);
- az R1–R6 és 5.1–5.3 újrabírálati pontok és a változtatások megfeleltetése a revíziós naplóban található (11. fejezet); a második bírálat F01–F12 pontjait a Revíziós napló 2 (12. fejezet), a harmadik újraaudit N01, N02, T01 pontjait a Revíziós napló 3 (13. fejezet), a P01–P09 nyitott korlátokat a 9.1 fejezet tartalmazza.

---

## 1. Egységes hatás, skaláregyenlet és $G_{\rm eff}$-normalizálás

### 1.1 Hatás és előjelkonvenció

$(-,+,+,+)$ szignatúrával, $g_{\mu\nu}$ Einstein-keretbeli metrikával és $\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}$ Jordan-keretbeli anyagmetrikával, $A(\phi)=e^{+\alpha_Y\phi}$:

$$
S=\int d^4x\,\sqrt{-g}\left[\frac{c^4}{16\pi G}R-\frac{c^4}{4\pi G\beta_{\rm Yuk}}\left(\frac{1}{2}g^{\mu\nu}\partial_\mu\phi\,\partial_\nu\phi+\frac{\phi^2}{2L^2}\right)-\rho_{\rm vac}^{\rm eff}c^2\right]+S_m\!\left[\Psi_m,\,A^2(\phi)g_{\mu\nu}\right].
$$

$T^{\mu\nu}=(2/\sqrt{-g})\,\delta S_m/\delta g_{\mu\nu}$ és $T=-\rho_m c^2+3p_m$ mellett a $\phi$ szerinti variáció (egység: m$^{-2}$):

$$
\Box\phi-\frac{\phi}{L^2}=-\frac{4\pi G}{c^4}\,\alpha_Y\beta_{\rm Yuk}\,T=+\frac{4\pi G}{c^2}\,\alpha_Y\beta_{\rm Yuk}\left(\rho_m-\frac{3p_m}{c^2}\right).
$$

### 1.2 Homogén FLRW-egyenlet (kozmikus másodperc, SI)

Homogén térre $\Box\phi=-c^{-2}(\ddot\phi+3H\dot\phi)$, így mindkét oldal egysége s$^{-2}$:

$$
\ddot{\bar\phi}+3H\dot{\bar\phi}+\frac{c^2}{L^2}\bar\phi=-4\pi G\,\alpha_Y\beta_{\rm Yuk}\left(\rho_m-\frac{3p_m}{c^2}\right).
$$

A skalár Compton-rátája $c/L=29979.2458\;\text{s}^{-1}$, és $H_0L/c=7.286\times 10^{-23}$.

### 1.3 Attraktor és ágkiválasztás

Mivel $4\pi\beta_{\rm Yuk}=1$, az adiabatikus partikuláris megoldás ($\ddot{\bar\phi},\,3H\dot{\bar\phi}\ll (c/L)^2\bar\phi$):

$$
\bar\phi_{\rm attr}(z)=-\alpha_Y\gamma_L\left(\rho_m(z)-\frac{3p_m(z)}{c^2}\right),\qquad \gamma_L=\frac{4\pi\beta_{\rm Yuk}GL^2}{c^2}=\frac{GL^2}{c^2}=7.42616\times 10^{-20}\;\text{m}^3/\text{kg},
$$

így $\bar\phi_{\rm attr}(0)=-1.9949\times 10^{-46}$ és $\bar\phi_{\rm attr}(z_*)\approx-2.589\times 10^{-37}$. Az általános homogén megoldás szabad oszcillációt is tartalmaz:

$$
\delta\phi\simeq C\,a^{-3/2}\cos\!\left[(c/L)\,t+\vartheta\right],\qquad \langle\rho_{\delta\phi}\rangle\propto a^{-3}.
$$

A $H/(c/L)\ll 1$ feltétel az adiabatikus közelítést támasztja alá, de ezt a módust nem tünteti el. **Az elemzés ezért a $C=0$ ($\Omega_{\delta\phi,0}=0$) ágat választja: az adiabatikus partikuláris megoldást elhanyagolható szabad oszcillációs energiával.** A skaláregyenlet kezdetiérték-integrálása nem történik meg, és az alábbi következtetések csak erre az ágra érvényesek.

### 1.4 A skalár energiasűrűsége

A hatásból $\rho_\phi=\dfrac{c^2}{8\pi G\beta_{\rm Yuk}}\left[(\dot\phi/c)^2+(\phi/L)^2\right]$ (tömegsűrűség), így az attraktoron:

$$
\Omega_{\phi,\rm attr}(0)=\frac{(c/L)^2\,\bar\phi_{\rm attr}(0)^2}{3\beta_{\rm Yuk}H_0^2}=3.140\times 10^{-47}.
$$

### 1.5 A $G_{\rm eff}$ normalizálása

Statikus pontszerű tömegre $\phi(r)=-\alpha_Y\beta_{\rm Yuk}\,\dfrac{GM}{c^2r}\,e^{-r/L}<0$, és az $A=e^{\alpha_Y\phi}$ anyagcsatolás $\alpha_Y^2\beta_{\rm Yuk}$ ötödikerő-erősséget ad. Lineáris perturbációkra, $k$ együttmozgó hullámszámmal:

$$
\frac{G_{\rm eff}(k,a)}{G_N}=1+\alpha_Y^2\beta_{\rm Yuk}\,\frac{(kL/a)^2}{1+(kL/a)^2}.
$$

Végig ezt a kanonikus normalizálást használjuk. A $2\alpha_Y^2\beta_{\rm Yuk}$ sor külön, kétszeres erősségű érzékenységi kontroll; nem a publikált V26 Step25 növekedési törvénye (5. fejezet). A lokális sajátidő $d\tau_J=A(\phi)\,d\tau_E=e^{-\alpha_Y|\phi|}d\tau_E$.

**G-szótár (F02 bírálati pont).** Négy jelölést választunk szét:

| Jelölés | Jelentés |
|:-----|:-----------------------------|
| $G_*$ | az 1.1 fejezet hatásában szereplő $G$ (Einstein-keret) |
| $G_T$ | Jordan-kereti tenzorcsatolás: $g_J=A^2g_E$ mellett a Ricci-tag $A^{-2}R/(16\pi G_*)$, így $1/G_T=A^{-2}/G_*$, $G_T=G_*A^2$ ($\phi<0$ esetén $G_T$ csökken) |
| $G_{\rm lab}$ | operatívan mért $G$: $r\ll L$ és árnyékolatlan kis testek esetén $\Phi=-G_*M(1+q\,e^{-r/L})/r$, így $G_{\rm lab}\approx G_*(1+q)$, $q=\alpha_Y^2\beta_{\rm Yuk}=0.0795775$ |
| $G_N$ | az 5. fejezet perturbációs egyenletének jelölése: $G_N:=G_*$ |

$kL/a\ll1$ skálákon ezért, ha a csatolást laboratóriumi érték kalibrálja, a növekedési egyenletben $G_{\rm eff}/G_{\rm lab}=1/(1+q)=0.9262883$ szerepel; az 5. fejezet $10^{-41}$-es állítása formális kontroll ugyanazon csupasz $G_*$ mellett. $f\sigma_8$-kizárási állítást nem teszünk, árnyékolási ágat nem vezetünk be.

### 1.6 Háttérnormalizálás és megfigyelői keret

Minden program egyetlen háttérmodult használ (`tal_cosmology_common.py`). A sima tömeges-neutrínó interpoláció önmagában $E^2(0)=1.00000013546$ értéket ad (eltérés: $1.3546\times 10^{-7}$); a sötétenergia-lezárás ezt elnyeli, így $E^2(0)=1$ pontosan és $H(0)=H_0$. Ez az eltérés ($\Omega_\nu$-ben $1.3546\times 10^{-7}$) az interpoláció $z=0$-beli normalizálási sodródása, nem a közelítés hibája. A tömeges neutrínófaj Fermi–Dirac energiaintegráljához képest (azonos relativisztikus és nemrelativisztikus határesetekkel) a $p=1.83$ interpoláció relatív hibája $z=1000$-nél $+1.149\times 10^{-3}$ a tömeges neutrínók sűrűségében és $+1.795\times 10^{-5}$ az összesített $H$-ban; a $0\le z\le10^4$ rácson a legnagyobb abszolút értékek $3.435\times 10^{-3}$ (tömeges komponens, $z\approx33.5$) és $1.795\times 10^{-5}$ (összesített $H$, $z\approx1000$). Ez a $\Lambda$CDM-szerű háttér közelítési hibája, minden összevetett modellnél közös, és elkülönül az elméleti TAL-korrekciótól ($\sim 10^{-46}$).

| $z$ | $\rho_\nu^{\rm FD}$ [$E^2$-egység] | $\rho_\nu^{p=1.83}$ [$E^2$-egység] | rel. hiba (tömeges) | rel. hiba ($H$) |
|-----:|-------:|-------:|------:|------:|
| 0 | $1.41991\times 10^{-3}$ | $1.41998\times 10^{-3}$ | $+4.449\times 10^{-5}$ | $+3.159\times 10^{-8}$ |
| 1 | $1.13610\times 10^{-2}$ | $1.13626\times 10^{-2}$ | $+1.355\times 10^{-4}$ | $+2.402\times 10^{-7}$ |
| 10 | $1.90137\times 10^{0}$ | $1.90428\times 10^{0}$ | $+1.529\times 10^{-3}$ | $+3.454\times 10^{-6}$ |
| 100 | $2.02529\times 10^{3}$ | $2.02475\times 10^{3}$ | $-2.670\times 10^{-4}$ | $-8.115\times 10^{-7}$ |
| 1000 | $1.27145\times 10^{7}$ | $1.27291\times 10^{7}$ | $+1.149\times 10^{-3}$ | $+1.795\times 10^{-5}$ |
| 3000 | $1.01910\times 10^{9}$ | $1.01945\times 10^{9}$ | $+3.388\times 10^{-4}$ | $+1.083\times 10^{-5}$ |

 A számítások $\Lambda$CDM-szerű $H(z)$ függvényeket és külön órafaktort használnak; a Jordan-keretbe való leképezés ($dt_J=A\,dt_E$, $a_J=A\,a_E$, $H_J=(H_E+\alpha_Y\dot\phi_E)/A$) ennél a $|\phi|$-nél vezetőrendű azonosság, nem numerikusan megoldott, csatolt Friedmann–skalár–anyag rendszer.

### 1.7 A kinetikus normalizálás verzióegyeztetése (R1 újrabírálati pont)

Az alábbi útvonalak a V26 `scalar-tensor-yukawa-fit/` könyvtárfához viszonyítottak (V26, DOI `10.5281/zenodo.23727585`).

- A 26. lépés első vázlata a kinetikus tagot $\frac{c^4}{16\pi G}\frac{1}{2\beta_{\rm Yuk}}(\partial\phi)^2=\frac{c^4}{32\pi G\beta_{\rm Yuk}}(\partial\phi)^2$ alakban írta, az 1.1 hatásában viszont $\frac{c^4}{8\pi G\beta_{\rm Yuk}}(\partial\phi)^2$ áll: $K_{\rm régi}/K_{\rm új}=0.25$.
- A V26 a $\frac{c^4}{16\pi G\tilde\alpha_L}(\partial\phi_L)^2$ tagot használja $\tilde\alpha_L=\alpha_Y/(8\pi)$ mellett, azaz $c^4/(2G\alpha_Y)$-t (`step23_gw170817_and_lensing_en.md:19` és `:29`, `README_en.md:50`, `gw170817_and_lensing_audit.py:104–113`). $\alpha_Y=1$, $\beta_{\rm Yuk}=1/(4\pi)$ mellett ez megegyezik az új együtthatóval ($K_{\rm V26}/K_{\rm új}=1.000000$), és a V26 ötödikerő-erőssége, $2\tilde\alpha_L\alpha_Y^2=0.079577$, egyenlő $\alpha_Y^2\beta_{\rm Yuk}$-kal. A V26 skaláregyenlete, $(\Box-L^{-2})\phi=-(\alpha_Y^2G/c^4)\varepsilon_m$ (`step25_cmb_and_cosmology_en.md:37`), és $\phi_{\rm attr}=\alpha_Y^2\gamma_L\rho$ (`cmb_and_cosmology_audit.py:327`, `:620`) a $\phi_{\rm V26}=-\phi_{\rm V27}$ előjelkonvenciótól eltekintve egyezik az 1.3 fejezettel. **Következtetés:** a V26 az új normalizálást használta; a $c^4/(32\pi G\beta)$ alak az első vázlat átírási inkonzisztenciája volt, a számított értékek már az új normalizálást követték. $\alpha_Y\neq1$ esetén a két verzió eltér, mert a V26 $\beta_{\rm Yuk}=\alpha_Y^3/(4\pi)$ kötést használ, itt pedig $\beta_{\rm Yuk}$ független paraméter.
- A $\phi_{\rm régi}=2\phi_{\rm új}$ térátdefiniálás a régi kinetikus tagot az újba viszi, de $A=e^{\alpha_{\rm régi}\phi_{\rm régi}}=e^{2\alpha_{\rm régi}\phi_{\rm új}}$, ezért $\alpha_{\rm új}=2\alpha_{\rm régi}$ szükséges. Azonos $\alpha_Y=1$ mellett a régi alak $4\alpha_Y^2\beta_{\rm Yuk}=0.318310$ ötödikerő-erősséget adna $\alpha_Y^2\beta_{\rm Yuk}=0.079577$ helyett.
- Változatlanul átvihető: minden háttér-, BAO-, CMB-, SNe-alak- és RSD-szám (a háttér skalárenergia-korrekciója $\Omega_\phi\sim10^{-47}$; az RSD közvetlen Yukawa-perturbációja külön, az 5. fejezetben számszerűsített elhanyagolható tag: $\delta G_{\rm eff}/G_*\approx3.8\times 10^{-41}$, $\delta f\sigma_8\approx1.67\times 10^{-40}$), $\bar\phi_{\rm attr}(0)$, $\Omega_{\phi,\rm attr}(0)$ és a bulk lokális benchmarkok $\alpha_Y=1$ mellett. A V26 háttere emellett egy $(\phi/4\pi)^2$ energiakonvenciót is használt (`cmb_and_cosmology_audit.py:330`, `:623`), amely $\Omega_\phi=1.582\times 10^{-50}$ értéket ad a hatással konzisztens $3.140\times 10^{-47}$ helyett; mindkettő elhanyagolható.
- A kétszeres növekedési konvenció ($2\alpha_Y^2\beta_{\rm Yuk}$) új címkéje **érzékenységi kontroll**, nem egyenértékű konvenció.
- A V21 top-hat belső külső benchmark, $|\phi|=0.0412$ (7. fejezet), a V21 top-hat fixpont-profil ($x=L^2(G/c^2)\bar\rho$) belső értéke (`stellar_tal_v21.py:140–160`, `data/v21/profiles_L10km.csv`: SLy, `gradient` mód, $r=2.0046409680935633$ km, $x=0.041239199439254325$; APR4, `mass` mód, $r=2.75417748649486$ km, $x=0.041233030079658145$), nem felszíni érték, és **nem automatikusan megoldása az új Helmholtz–Yukawa hatásnak**; külső bemenet marad.

---

## 2. DESI DR1 BAO: $12\times 12$ kovariancia és profilozott kontrollok

A 12 elsődleges megfigyelés a Cobaya DESI DR1 vektora és kovarianciája (`desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`, `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`). A külön LRG3 és ELG1 egytracer-sorok csak diagnosztikák, a likelihoodba nem kerülnek; a $z=0.845$ ELG1 sor ($D_V/r_d=18.33\pm0.62$) **igazolt elsődleges forrás nélküli, választott illusztráció**, nem DESI-mérés (nem szerepel a DESI III Table 18-ban és a DESI VI Table 1-ben). Az elsődleges binek $\chi^2$-e a $12\times 12$ kovarianciával, a (dőlt) diagnosztikai soroké a kerekített táblázati hibákkal számolt. A táblázat jóslatai anizotrop binekre $D_M/r_d$; $D_H/r_d$, izotrop binekre $D_V/r_d$.

| Tracer | $z_{\rm eff}$ | DESI DR1 ($X/r_d$) | M4 jóslat | $\chi^2_{\rm bin}$ M4 | $\chi^2_{\rm bin}$ 2P-prof. | $\chi^2_{\rm bin}$ merev 72.8 |
|:-----------|-----:|:------------------------------------------|:--------------------|-----:|-----:|-----:|
| **BGS** | 0.295 | $D_V$: $7.93\pm0.15$ | $8.053$ | 0.724 | 0.090 | 9.69 |
| **LRG1** | 0.510 | $D_M$: $13.62\pm0.25$; $D_H$: $20.98\pm0.61$ ($r=-0.445$) | $13.495$; $22.733$ | 8.963 | 6.010 | 24.43 |
| **LRG2** | 0.706 | $D_M$: $16.85\pm0.32$; $D_H$: $20.08\pm0.60$ ($r=-0.420$) | $17.694$; $20.164$ | 8.967 | 3.097 | 12.92 |
| **LRG3+ELG1** | 0.930 | $D_M$: $21.71\pm0.28$; $D_H$: $17.88\pm0.35$ ($r=-0.389$) | $21.917$; $17.608$ | 0.826 | 2.031 | 75.07 |
| *LRG3 (sep.) (diagn.)* | 0.922 | $D_M$: $21.81\pm0.31$; $D_H$: $17.83\pm0.38$ ($r=-0.393$) | $21.775$; $17.692$ | *0.209* | *2.931* | *69.61* |
| *ELG1 (z=0.955) (diagn.)* | 0.955 | $D_V$: $20.01\pm0.41$ | $20.229$ | *0.285* | *0.004* | *9.77* |
| *ELG1 (választott illusztráció, nem DESI-mérés)* | 0.845 | $D_V$: $18.33\pm0.62$ | $18.666$ | *0.294* | *0.021* | *2.86* |
| **ELG2** | 1.317 | $D_M$: $27.79\pm0.69$; $D_H$: $13.82\pm0.42$ ($r=-0.444$) | $28.017$; $14.094$ | 0.883 | 0.380 | 18.57 |
| **QSO** | 1.491 | $D_V$: $26.07\pm0.67$ | $26.026$ | 0.005 | 0.141 | 8.74 |
| **Ly$\alpha$ QSO** | 2.330 | $D_M$: $39.71\pm0.94$; $D_H$: $8.52\pm0.17$ ($r=-0.477$) | $39.173$; $8.615$ | 0.418 | 0.990 | 44.93 |
| **Összeg, 12 elsődleges megfigyelés (12×12 kov.)** | | | | **20.785** | **12.741** | **194.34** |

- Rögzített M4-pont ($H_0=67.40$, $r_d=147.0922\;\text{Mpc}$): $\chi^2=20.7852$ 12 megfigyelésre, $p=0.0536$ (kerekített Table 1 bemenettel: $20.6996$).
- Merev $H_0=72.8$ változatlan alakkal és $r_d$-vel: $\chi^2=194.342$, $p=4.79\times 10^{-35}$ (kerekített Table 1 bemenettel: $195.178$).

**Profilozott kontrollok** (helyi Fisher-hibák; a $p$-értékek a megadott Gauss-likelihood feltételes gyakorisági értékei, nem a TAL igazságának valószínűségei):

| Kontroll | $\Omega_m$ | $P=H_0r_d$ [km/s] | $\chi^2$ (12×12) | dof | $p$ | Szükséges $r_d$ $H_0=72.8$ mellett [Mpc] |
|:------------------------------------|-----:|------------:|------:|---:|------:|------------:|
| Rögzített M4-pont (nincs szabad BAO-paraméter) | 0.3148 | 9914.02 | 20.7852 | 12 | 0.0536 | 136.18 |
| M4-alak, merev $H_0=72.8$, változatlan $r_d$ | 0.3148 | 10708.32 | 194.3420 | 12 | $4.79\times 10^{-35}$ | — (147.09 rögzített) |
| M4-alak, $P=H_0r_d$ profilozva | 0.3148 | 10029.65 $\pm$ 47.43 | 14.7015 | 11 | 0.1966 | 137.77 $\pm$ 0.65 |
| Sík $\Lambda$CDM-alak, $\Omega_m$ és $P$ profilozva | 0.2939 $\pm$ 0.0145 | 10193.67 $\pm$ 126.47 | 12.7405 | 10 | 0.2385 | 140.02 $\pm$ 1.74 |

- Rögzített $E(z)$ mellett a BAO csak a $P=H_0r_d$ szorzatot korlátozza (rang-1 információ a $(H_0,r_d)$ síkon).
- Kétparaméteres kontroll: $\Omega_m=0.2939\pm0.0145$, $P=10193.67\pm126.47\;\text{km/s}$, $\chi^2_{\min}=12.7405$ (10 dof, $p=0.2385$).
- Szükséges $r_d$ pontosan $H_0=72.8$ mellett: $140.02\pm1.74\;\text{Mpc}$; a $H_0=72.8\pm1.0$ bizonytalanságot is propagálva: $\pm2.59\;\text{Mpc}$. A $72.8\pm1.0$ **feltételezett lokális $H_0$ benchmark-prior** (választott célérték, F05 bírálati pont); nem a Riess és mtsai. 2022-es $73.04\pm1.04$ értéke, és konkrét SH0ES-kiadást nem tulajdonítunk neki.
- A modellfüggő Planck $r_d=147.09\pm0.26\;\text{Mpc}$ priorral: $H_0=69.30\pm0.87\;\text{km/s/Mpc}$.
- Szemléltető közös teszt (BAO + Planck $r_d$ prior + lokális $H_0=72.8\pm1.0$, $\Omega_m,H_0,r_d$ profilozva): legjobb illesztés $\Omega_m=0.2719$, $H_0=70.806$, $r_d=147.025\;\text{Mpc}$, $\chi^2_{\rm közös}=19.713$, felbontva: $\chi^2_{\rm BAO}=15.67431585$, $\chi^2_{r_d\,\rm prior}=0.06234684$, $\chi^2_{H_0\,\rm prior}=3.97660606$; $\Delta\chi^2=6.973$ a csak-BAO minimumhoz képest ($\sqrt{\Delta\chi^2}=2.64$). Ez nem Planck+TAL kizárási szignifikancia.

![DESI DR1 BAO, Step 25 M4 és profilozott kontrollok (`step1_1_desi_dr1_bao_plot.png`)](step1_1_desi_dr1_bao_plot.png){ width=96% }

---

## 3. Planck 2018 akusztikus skála modellspecifikus hanghorizontokkal

Minden modellnél az $r_s(z_*=1089.80)$ és $r_d(z_{\rm drag}=1059.94)$ a modell saját $H(z)$ és $c_s(z)$ függvényével újraintegrált (rögzített $z_*$, $z_{\rm drag}$ kontrollközelítésként). Az $(\omega_b,\omega_{cb},100\theta_*)$ 3D Gauss-prior **szemléltető, feltételezett prior** (Planck 2018 marginális középértékek és hibák, közelítő korrelációs mátrixszal, amelynek láncproveniencája itt nincs dokumentálva); nem a Planck TT,TE,EE+lowE+lensing likelihood. Referencia: $100\theta_*=1.04109\pm0.00030$.

| Modell | $H_0$ | $\Omega_m$ | $\Omega_K$ | $\omega_{cb}$ | $r_s(z_*)$ [Mpc] | $r_d(z_{\rm drag})$ [Mpc] | $100\theta_*$ | $\chi^2_{\rm 3D,szeml.}$ | $\chi^2_{\rm BAO}$ (újraint. $r_d$) |
|:------------------------------|----:|------:|------:|------:|-------:|-------:|-------:|---------:|-------:|
| Step 25 M4 / TAL | 67.40 | 0.3148 | 0.0000 | 0.14237 | 144.443 | 147.092 | 1.04132 | 0.757 | 20.79 |
| Planck 2018 referencia | 67.36 | 0.3152 | 0.0000 | 0.14237 | 144.443 | 147.092 | 1.04120 | 0.173 | 21.19 |
| Minden $z$-n merev $H\to sH$ (kinematikai kontroll) | 72.80 | 0.3148 | 0.0000 | 0.14237 | 133.729 | 136.182 | 1.04132 | 0.757 | 20.79 |
| *$D_M/s$ rögzített korai $r_s,r_d$-vel (inkonzisztens geometriai szemléltetés, nem modell)* | 72.80 | 0.3148 | 0.0000 | 0.14237 | 144.443 | 147.092 | 1.12475 | 101040.4 | 194.34 |
| Sík $H_0=72.8$, azonos $\Omega_m$ | 72.80 | 0.3148 | 0.0000 | 0.16620 | 138.879 | 141.375 | 1.08075 | 24800.4 | 43.73 |
| Sík $H_0=72.8$, $\Omega_m=0.244828$ | 72.80 | 0.2448 | 0.0000 | 0.12911 | 147.903 | 150.651 | 1.04174 | 196.000 | 70.30 |
| Sík $H_0=72.8$, $\theta_*$-hez illesztett $\Omega_m\approx0.2442$ | 72.80 | 0.2442 | 0.0000 | 0.12877 | 147.996 | 150.747 | 1.04132 | 209.592 | 70.84 |
| Nyílt $H_0=72.8$, rögzített $\omega_{cb}$ | 72.80 | 0.2698 | 0.0098 | 0.14237 | 144.442 | 147.092 | 1.04132 | 0.727 | 61.21 |

- Step 25 M4: $r_s=144.4428\;\text{Mpc}$, $r_d=147.0922\;\text{Mpc}$, $100\theta_*=1.041320$ (pull $0.77$), $\chi^2_{\rm 3D}=0.757$.
- Minden $z$-n merev $H\to sH$ átskálázás változatlan $c_s$ mellett: $r_s\to r_s/s=133.7286\;\text{Mpc}$, $r_d=136.1816\;\text{Mpc}$ és $D_M\to D_M/s$, így $100\theta_*=1.041320$ változatlan. Csak $D_M$ skálázása $1.124750$ értéket adna; ezt az inkonzisztens kombinációt már nem használjuk tesztként.
- Sík $H_0=72.8$, $\Omega_m=0.244828$: az újraintegrálás $r_s=147.9028$, $r_d=150.6513\;\text{Mpc}$ és $100\theta_*=1.041745$ értéket ad (rögzített korai $r_s$-sel $1.017374$ lenne). Ez **nem** jelenti, hogy a modell átmegy a Planck-teszten: $\omega_{cb}=0.12911$ eltér a priortól, és a BAO $\chi^2$ az újraintegrált $r_d$-vel $70.30$ (a régi rögzített $r_d$-vel: $29.30$). A korai és késői részt együtt kell újrailleszteni.
- **Két $z_*$-eltolási diagnosztika (F12 bírálati pont).** A Hu–Sugiyama-illesztés $z_*\to z_*+\Delta z_*$ eltolásával (ennél a modellnél $z_*+\Delta z_*=1088.6345$) a CSV két külön nevű értéket tárol: `theta_star_100_rs_only_shifted_fixed_DM` $=100\,r_s(z_*+\Delta z_*)/D_M(1089.80)=1.04248582$ (csak az $r_s$ végpontja eltolva) és `theta_star_100_both_endpoints_shifted` $=100\,r_s(z_*+\Delta z_*)/D_M(z_*+\Delta z_*)=1.04250287$. Egyik sem rekombinációs megoldás.
- **Késői átmenet (R2 újrabírálati pont).** A korábbi „csak késői” sor a $D_M(z_*)/s$ értéket változatlan korai $r_s,r_d$-vel kombinálta ($100\theta_*=1.124750$); nincs átmeneti vöröseltolódása, ezért **inkonzisztens geometriai szemléltetésként sorolódik át** (dőlt sor a táblázatban), és kikerül a konzisztens modellekre vonatkozó állításokból. Definiált kontroll a szakaszos $H(z)=s\,H_{\rm M4}(z)$ ($z<z_t$), illetve $H_{\rm M4}(z)$ ($z\ge z_t$), $s=72.8/67.4=1.080119$, numerikusan integrált $D_M(z_*)$-gal és változatlan $r_s=144.442796\;\text{Mpc}$-vel (M4: $D_M(z_*)=13871.1192\;\text{Mpc}$). $H$ ugrik $z_t$-nél, így ez kinematikai kontroll, nem dinamikai modell; a pull oszlop formális 1D távolság a Planck-értéktől, nem likelihood-szignifikancia. $z_t\ge2.5$ esetén minden DESI-bin az átskálázott tartományba esik, így a BAO $\chi^2$ a merev értékkel egyezik: $194.342$.

| $z_t$ | $D_M(z_*)$ [Mpc] | $100\theta_*$ | pull a $100\theta_*^{\rm Pl}$-hoz (formális) |
|-----:|-----------:|---------:|--------:|
| 2.5 | 13428.2963 | 1.075660 | 115.2 |
| 10 | 13156.6704 | 1.097867 | 189.3 |
| 100 | 12920.1717 | 1.117963 | 256.2 |
| 1000 | 12843.5875 | 1.124630 | 278.5 |
| $D_M/s$ (nincs $z_t$; inkonzisztens szemléltetés) | 12842.2175 | 1.124750 | — |

---

## 4. Szintetikus távolságmodulus-alakkontroll (nem Pantheon+ / Union3 likelihood)

Ez a rész **nem** használ megfigyelt SN-katalógust: 18 vöröseltolódás-csomópont ($z=0.015$–$2.26$) reprezentatív $\sigma_\mu$ hibákkal, és egy $H_0=72.8$, $\Omega_m=0.334$ sík modellből számolt szintetikus referenciagörbe. Nincs SN-azonosító, megfigyelt $\mu$ és statisztikai+szisztematikus kovariancia. Az alábbi $\chi^2$ értékek modell–modell alakmetrikák, nem megfigyelési illeszkedési mértékek.

- Azonos $E(z)$ geometria esetén az eltolás a konstans $5\log_{10}(72.8/67.4)=+0.167357\;\text{mag}$ azonosság lenne (pozitív előjel: M4 mínusz a 72.8-as referencia).
- A táblázat $z=0.015$-beli $0.16782511\;\text{mag}$ különbsége modellkülönbség (a referencia $\Omega_m=0.334$), nem ez az azonosság.
- **Szimmetrikus nullpontkezelés (R5 újrabírálati pont):** mindkét modell ugyanazt a profilozott konstans eltolást kapja (abszolút fényesség marginalizálása), $1/\sigma_\mu^2$ súlyokkal:

| Modell | legjobb eltolás [mag] | $\chi^2$ eltolás nélkül | $\chi^2$ eltolással (17 dof) |
|:--------------------|------:|------:|------:|
| M4 ($H_0=67.4$) | 0.17595029 | 522.00704 | 0.916614 |
| Sík $H_0=72.8$, $\Omega_m=0.2448$ | 0.04159799 | 51.90332 | 22.777527 |

- Az M4 nagy csökkenése csak azt mutatja, hogy a konstans levonása után alakja közel van a referenciához; a csomópontok és a reprezentatív hibák szintetikusak, ezért a csomópontok vagy hibák változtatása nem jelent fizikai javulást.
- Az ábra mindkét modellt profilozott konstans nullpont-eltolással és anélkül is mutatja (az alacsony $\Omega_m$-ű modell eltolással: $\chi^2=22.777527$); $N=18$ csomópont, egy profilozott eltolás után $\nu=17$ szabadsági fok. Az eltolás profilozott zavaró konstans, nem megfigyelési kalibráció (F09 bírálati pont).
- Kalibrálatlan SNe-adatokból önmagában nincs $H_0$-mérés ($H_0$–$M_B$ degeneráció). Valódi katalógusállításhoz a Pantheon+ vagy Union3 likelihoodot és kovarianciát kellene beolvasni; ez a csomag ezt nem teszi.

![Szintetikus SNe Ia távolságmodulus-alakkontroll (`step1_3_pantheon_sne_plot.png`)](step1_3_pantheon_sne_plot.png){ width=96% }

| $z$ csomópont | $\sigma_\mu$ [mag] | $\mu_{\rm ref}$ (72.8, 0.334) | $\mu_{\rm M4}$ (67.4) | $\Delta\mu_{\rm M4-ref}$ | $\Delta\mu_{\rm M4,eltolt-ref}$ | $\Delta\mu_{0.2448-{\rm ref}}$ | $\Delta\mu_{0.2448,{\rm eltolt-ref}}$ |
|------:|------:|----------:|----------:|---------:|---------:|---------:|---------:|
| 0.015 | 0.042 | 33.9780 | 34.1459 | +0.1678 | -0.0081 | +0.0022 | -0.0394 |
| 0.023 | 0.038 | 34.9189 | 35.0869 | +0.1681 | -0.0079 | +0.0033 | -0.0383 |
| 0.035 | 0.035 | 35.8493 | 36.0177 | +0.1684 | -0.0075 | +0.0051 | -0.0365 |
| 0.050 | 0.032 | 36.6467 | 36.8156 | +0.1689 | -0.0070 | +0.0073 | -0.0343 |
| 0.075 | 0.030 | 37.5641 | 37.7338 | +0.1697 | -0.0063 | +0.0109 | -0.0307 |
| 0.100 | 0.028 | 38.2245 | 38.3950 | +0.1704 | -0.0055 | +0.0145 | -0.0271 |
| 0.150 | 0.026 | 39.1726 | 39.3446 | +0.1720 | -0.0040 | +0.0217 | -0.0199 |
| 0.200 | 0.025 | 39.8603 | 40.0337 | +0.1734 | -0.0025 | +0.0287 | -0.0129 |
| 0.250 | 0.025 | 40.4036 | 40.5784 | +0.1749 | -0.0011 | +0.0357 | -0.0059 |
| 0.300 | 0.026 | 40.8543 | 41.0306 | +0.1763 | +0.0003 | +0.0425 | +0.0009 |
| 0.400 | 0.028 | 41.5785 | 41.7575 | +0.1789 | +0.0030 | +0.0556 | +0.0140 |
| 0.500 | 0.030 | 42.1508 | 42.3322 | +0.1814 | +0.0055 | +0.0679 | +0.0263 |
| 0.650 | 0.033 | 42.8339 | 43.0187 | +0.1848 | +0.0088 | +0.0849 | +0.0433 |
| 0.800 | 0.036 | 43.3808 | 43.5685 | +0.1877 | +0.0118 | +0.0999 | +0.0583 |
| 1.000 | 0.042 | 43.9724 | 44.1634 | +0.1910 | +0.0151 | +0.1173 | +0.0757 |
| 1.300 | 0.055 | 44.6699 | 44.8649 | +0.1950 | +0.0191 | +0.1384 | +0.0968 |
| 1.700 | 0.075 | 45.3817 | 45.5806 | +0.1989 | +0.0229 | +0.1593 | +0.1177 |
| 2.260 | 0.110 | 46.1321 | 46.3347 | +0.2026 | +0.0267 | +0.1796 | +0.1380 |

---

## 5. RSD $f\sigma_8(z)$: diagonális benchmark és variációs TAL-korrekció

A $\delta_m''+\left(2+\frac{d\ln H}{d\ln a}\right)\delta_m'-\frac{3}{2}\Omega_{cb}(a)\frac{G_{\rm eff}(k,a)}{G_N}\delta_m=0$ növekedési egyenletet DOP853 integrálja (`rtol=1e-12`, `atol=1e-14`, dense output); $\sigma_{8,0}=0.8111$ rögzített Planck-bemenet, nem TAL-ból levezetett amplitúdó.

- $k=0.1\,h\,\text{Mpc}^{-1}$ mellett: $kL=2.18429\times 10^{-20}$ ($(kL)^2=4.7711\times 10^{-40}$) $z=0$-nál; $\Delta G_{\rm eff}/G_N=3.7967\times 10^{-41}$ (kanonikus $\alpha_Y^2\beta_{\rm Yuk}$) vagy $7.5934\times 10^{-41}$ (kétszeres konvenció).
- float64-ben $1+O(10^{-41})$ pontosan 1, ezért az $f\sigma_8^{\rm TAL}-f\sigma_8^{\rm GR}$ közvetlen különbség lebegőpontos azonosság miatt $0$, nem 40 jegyű mérés eredménye. Az elsőrendű variációs egyenlet a 12 ponton $\max\lvert\delta f\sigma_8\rvert=1.666\times 10^{-40}$ (kanonikus) és $3.333\times 10^{-40}$ (kétszeres érzékenységi kontroll) értéket ad; az újrabírálat független, átskálázott ellenőrzése $1.666\times10^{-40}$-t erősített meg (maximum $1.6663681521\times10^{-40}$, relatív eltérés $\le5.1\times10^{-6}$). A float64-ág $G_{\rm eff}/G_N=1+\Delta G_{\rm eff,kan}/G_N$ alakot használ (5.1 újrabírálati pont). Az ábra jelmagyarázata a kanonikus $1.6664\times 10^{-40}$ maximumot mutatja; a kétszeres érték csak a külön érzékenységi kontroll (F07 bírálati pont). Itt $G_N:=G_*$ (G-szótár, 1.5 fejezet).
- Diagonális $\chi^2$ (bin-kovarianciák, átfedések és AP-kezelés nélkül): M4 $\Omega_{cb}$ csomósodási forrással: $10.1233$; az egész anyaggal: $10.2176$; sík $H_0=72.8$, $\Omega_m=0.2448$: $9.4350$ (12 pont). Ez benchmark, nem kovarianciás likelihood.

![RSD $f\sigma_8(z)$ diagonális benchmark (`step1_4_rsd_fsigma8_plot.png`)](step1_4_rsd_fsigma8_plot.png){ width=96% }

| Felmérés | $z$ | $f\sigma_8$ mért | $f\sigma_8$ M4 | pull M4 | $f\sigma_8$ (72.8, 0.2448) | $\lvert\delta f\sigma_8\rvert$ variációs | $\Delta G_{\rm eff}/G_N$ |
|:------------------------------|-----:|-----------:|-------:|------:|-------:|----------:|----------:|
| 6dFGS (Beutler+ 2012) | 0.067 | $0.423\pm0.055$ | 0.4424 | +0.35 | 0.3887 | $7.57\times 10^{-41}$ | $4.32\times 10^{-41}$ |
| SDSS MGS (Howlett+ 2015) | 0.150 | $0.530\pm0.160$ | 0.4574 | -0.45 | 0.4083 | $8.30\times 10^{-41}$ | $5.02\times 10^{-41}$ |
| BOSS DR12 (Alam+ 2017 z1) | 0.380 | $0.497\pm0.045$ | 0.4750 | -0.49 | 0.4411 | $1.01\times 10^{-40}$ | $7.23\times 10^{-41}$ |
| WiggleZ (Blake+ 2012 z1) | 0.440 | $0.413\pm0.080$ | 0.4750 | +0.77 | 0.4451 | $1.06\times 10^{-40}$ | $7.87\times 10^{-41}$ |
| BOSS DR12 (Alam+ 2017 z2) | 0.510 | $0.459\pm0.038$ | 0.4732 | +0.37 | 0.4478 | $1.11\times 10^{-40}$ | $8.66\times 10^{-41}$ |
| WiggleZ (Blake+ 2012 z2) | 0.600 | $0.390\pm0.063$ | 0.4686 | +1.25 | 0.4487 | $1.17\times 10^{-40}$ | $9.72\times 10^{-41}$ |
| BOSS DR12 (Alam+ 2017 z3) | 0.610 | $0.436\pm0.034$ | 0.4680 | +0.94 | 0.4486 | $1.17\times 10^{-40}$ | $9.84\times 10^{-41}$ |
| eBOSS DR16 LRG (Bautista+21) | 0.700 | $0.473\pm0.041$ | 0.4612 | -0.29 | 0.4467 | $1.23\times 10^{-40}$ | $1.10\times 10^{-40}$ |
| WiggleZ (Blake+ 2012 z3) | 0.730 | $0.437\pm0.072$ | 0.4586 | +0.30 | 0.4457 | $1.25\times 10^{-40}$ | $1.14\times 10^{-40}$ |
| eBOSS DR16 ELG (Tamone+ 20) | 0.850 | $0.315\pm0.095$ | 0.4471 | +1.39 | 0.4395 | $1.32\times 10^{-40}$ | $1.30\times 10^{-40}$ |
| FastSound (Okumura+ 2016) | 1.360 | $0.482\pm0.116$ | 0.3896 | -0.80 | 0.3960 | $1.60\times 10^{-40}$ | $2.11\times 10^{-40}$ |
| eBOSS DR16 QSO (Neveux+ 20) | 1.480 | $0.462\pm0.045$ | 0.3762 | -1.91 | 0.3844 | $1.67\times 10^{-40}$ | $2.34\times 10^{-40}$ |

---

## 6. Algebrai BBN nagyságrendi benchmark

Nukleáris reakcióhálózat (pl. PArthENoPE) nem fut. Az elemgyakoriságok linearizált standard értékek $\eta_{10}=6.1245$ mellett, algebrai TAL-eltolással; a $H(T)$ egyszerű extrapoláció $e^\pm$-termodinamika és változó $g_*$ nélkül, a skalárforrás pedig csak a poranyagot tartalmazza, a termikus plazma trace-ét nem. A pullok csak a megadott megfigyelési hibával osztanak.

**Proxyk és források (R6 újrabírálati pont).** A $|\delta G/G|\approx2|\phi|$ óraproxy és a $0.08$ ($\delta Y_p\approx0.08\,\delta H/H$), illetve $1.6$ ($\delta\ln({\rm D/H})\approx1.6\,\delta H/H$) együtthatók **feltételezett nagyságrendi proxyk**, itt nem levezetettek. Keret és előjel (F01 bírálati pont): $g_J=A^2g_E$ mellett a Jordan-kereti Ricci-tag $A^{-2}/G_*$ szorzót hordoz ($A^{-2}$ a Planck-tömeg négyzetét szorozza), így a tenzorcsatolás $1/G_T=A^{-2}/G_*$, azaz $G_T=G_*A^2=G_*e^{2\alpha_Y\phi}$, és $\phi<0$ esetén $G_T$ **csökken**. A $G_T$ tenzorcsatolás, az operatívan mért $G_{\rm lab}$ (1.5 fejezet) és a csupasz $2|\phi|$ proxy három különböző mennyiség; csak a $2\alpha_Y|\phi|$ nagyságot használjuk, az előjelet és a referenciaállapotot nem propagáljuk. Összevetésként a Kneller–Steigman-illesztések (New J. Phys. 6, 117, 2004) $\delta Y_p\approx0.16\,(S-1)$ és $\delta\ln({\rm D/H})\approx1.6\,(S-1)$ értéket adnak, $S=H'/H$: az 1.6-os együttható egyezik, a 0.08 kétszer kisebb (ez $|\delta H/H|\lesssim10^{-25}$ mellett lényegtelen). Bemeneti értékek: $Y_p^{\rm std}=0.24709$ a PRIMAT központi értéke (Pitrou, Coc, Uzan & Vangioni 2018, Phys. Rep. 754, 1); $10^5\,{\rm D/H}^{\rm std}=2.5147$ **számtani proxy** (P02 bírálati pont): linearizált benchmarkérték igazolt egyedi forrás nélkül, nem forrásolt BBN-jóslat (a PRIMAT 2018 2.459-et ad); $Y_p^{\rm obs}=0.245\pm0.003$ a PDG 2024 BBN-áttekintés ajánlott értéke (Fields, Molaro & Sarkar), összhangban Aver és mtsai. 2021 eredményével (JCAP 03, 027: $0.2453\pm0.0034$); $10^5\,{\rm D/H}^{\rm obs}=2.527\pm0.030$ Cooke, Pettini & Steidel 2018 (ApJ 855, 102) értéke — a `DH_e5_obs_PDG` CSV-mezőnév örökölt címke (a PDG 2024 ajánlása $2.547\pm0.025$).

| Korszak | $z$ | $H/(c/L)$ | $\rho_m$ [kg/m$^3$] | $\lvert\bar\phi_{\rm attr}\rvert$ | $\delta H/H$ (skalárenergia) | $\delta G/G$ (óra) | $\lvert\delta Y_p\rvert$ |
|:------------------------------|------:|------:|------:|------:|------:|------:|------:|
| BBN vége, $T\approx 0.02$ MeV | $1.0\times 10^{8}$ | $6.992\times 10^{-9}$ | $2.686\times 10^{-3}$ | $1.995\times 10^{-22}$ | $1.705\times 10^{-27}$ | $3.990\times 10^{-22}$ | $1.596\times 10^{-23}$ |
| D-szűk keresztmetszet, $T\approx 0.09$ MeV | $4.0\times 10^{8}$ | $1.119\times 10^{-7}$ | $1.719\times 10^{-1}$ | $1.277\times 10^{-20}$ | $2.728\times 10^{-26}$ | $2.553\times 10^{-20}$ | $1.021\times 10^{-21}$ |
| $z=10^9$, $T\approx 0.23$ MeV | $1.0\times 10^{9}$ | $6.992\times 10^{-7}$ | $2.686\times 10^{0}$ | $1.995\times 10^{-19}$ | $1.705\times 10^{-25}$ | $3.990\times 10^{-19}$ | $1.596\times 10^{-20}$ |

- $Y_p=0.24709$, összevetve: $0.245\pm0.003$ (pull $0.70$); $10^5\,{\rm D/H}=2.5147$, összevetve: $2.527\pm0.030$ (pull $-0.41$).
- Maximumok a három korszakon: $H/(c/L)\le6.992\times 10^{-7}$, $|\bar\phi_{\rm attr}|\le1.995\times 10^{-19}$, $\delta H/H\le1.705\times 10^{-25}$, $|\delta Y_p|\le1.596\times 10^{-20}$. A kis adiabatikus korrekció kvalitatív várakozása fennmarad; ez nem általános BBN-korlát.

---

## 7. Lokális Yukawa állandó sűrűségű bulk-benchmarkok (`expm1`)

A sorok állandó sűrűségre vett $|\phi|_{\rm bulk}=\alpha_Y\gamma_L\rho$ szorzatok. Az $R\gg L$ félteres peremérték ($|\phi|_{\rm bulk}/2$) csak ott szerepel, ahol peremgeometria definiált (Nap átlagsűrűség, Föld átlagsűrűség, fehér törpe belseje; státusz: `half-space boundary`). A homogén FLRW-, a kozmikus üreg- és a Tejútrendszer-korong sornak itt nincs peremgeometriája: a peremcellájuk n.a. (CSV/JSON-ban null, `boundary_status` = `n.a. (no boundary geometry)`, a perem-$\Delta H_0$ null); a bulk értéküket nem másoljuk a peremoszlopba, és nem felezzük (N01 bírálati pont). Az értékek nem a Nap vagy a Föld megoldott felszíni terei (ahhoz a valódi forrás és sűrűségprofil kell). Az órafaktorok `expm1`-gyel készülnek, a kerekítési nullák elkerülésére; a látszólagos eltolás $\Delta H_0^{\rm app}=67.40\,{\rm expm1}(\alpha_Y|\phi|)$. A „V21 top-hat belső külső benchmark” sor ($|\phi|=0.0412$) **külső bemenet, itt nem reprodukált**: a V21 top-hat fixpont-profil belső értéke (SLy `gradient`, $r=2.0046$ km, $x=0.041239$; APR4 `mass`, $r=2.7542$ km, $x=0.041233$; 1.7 fejezet), nem felszíni érték (a peremoszlop n.a.), és nem a jelen hatás megoldása.

| Benchmark | $\rho$ [kg/m$^3$] | $\lvert\phi\rvert_{\rm bulk}=\alpha_Y\gamma_L\rho$ | $\lvert\phi\rvert$ ($R\gg L$ perem) | ${\rm expm1}(-\alpha_Y\lvert\phi\rvert)$ | $\Delta H_0^{\rm app}$ [km/s/Mpc] | az 5.4-es rés hányada |
|:------------------------------|------:|------:|------:|------:|------:|------:|
| Homogén FLRW, $z=0$ | $2.686\times 10^{-27}$ | $1.995\times 10^{-46}$ | n.a. | $-1.995\times 10^{-46}$ | $1.345\times 10^{-44}$ | $2.490\times 10^{-45}$ |
| Kozmikus üreg ($0.2\rho_{m,0}$) | $5.373\times 10^{-28}$ | $3.990\times 10^{-47}$ | n.a. | $-3.990\times 10^{-47}$ | $2.689\times 10^{-45}$ | $4.980\times 10^{-46}$ |
| Tejútrendszer-korong (átlag) | $1.000\times 10^{-20}$ | $7.426\times 10^{-40}$ | n.a. | $-7.426\times 10^{-40}$ | $5.005\times 10^{-38}$ | $9.269\times 10^{-39}$ |
| Nap átlagsűrűség (bulk) | $1408.0000$ | $1.046\times 10^{-16}$ | $5.228\times 10^{-17}$ | $-1.046\times 10^{-16}$ | $7.047\times 10^{-15}$ | $1.305\times 10^{-15}$ |
| Föld átlagsűrűség (bulk) | $5514.0000$ | $4.095\times 10^{-16}$ | $2.047\times 10^{-16}$ | $-4.095\times 10^{-16}$ | $2.760\times 10^{-14}$ | $5.111\times 10^{-15}$ |
| Fehér törpe belseje (bulk) | $1.000\times 10^{9}$ | $7.426\times 10^{-11}$ | $3.713\times 10^{-11}$ | $-7.426\times 10^{-11}$ | $5.005\times 10^{-9}$ | $9.269\times 10^{-10}$ |
| V21 top-hat belső külső benchmark | n.a. | $0.0412$ | n.a. | $-0.0404$ | $2.8349$ | $0.5250$ |

- Nap átlagsűrűségű bulk: $\Delta H_0^{\rm app}=7.047\times 10^{-15}\;\text{km/s/Mpc}$ (nem matematikai nulla).
- Csak a külső V21 belső benchmark éri el a $\Delta H_0^{\rm app}=2.835\;\text{km/s/Mpc}$ értéket (az 5.4 km/s/Mpc rés $0.525$ hányada); ez feltételes szemléltetés; egy lokális órafaktor és a kozmológiai $H_0$ azonosítása külön levezetést igényel, definiált megfigyelői konvencióval.

---

## 8. No-go függelék

### 8.1 Az M2 / M3 ágak

- `M2` (sík, kompenzálatlan negatív vákuumsűrűség): $E^2(0)=-0.37020$, így a deklarált, ma pozitív $H_0$-ra vett normalizálás nem fizikai. Az $E^2(z)=0$ gyöke a sugárzást és neutrínót is tartalmazó háttérben $z\approx0.295737$ (csak anyagos közelítés: $0.295878$). A negatív vákuumenergia általában nincs kizárva; ez a konkrét normalizált pont igen.
- `M3` (nyílt, $\Omega_K=1.3702$), az M3 $H(z)$-vel újraintegrált **saját** hanghorizontokkal (F03 bírálati pont): $r_s(z_*)=144.36622618\;\text{Mpc}$, $r_d(z_{\rm drag})=147.01159289\;\text{Mpc}$, $D_M(z_*)=34916.4\;\text{Mpc}$, $100\theta_*=0.41346227$ (formális geometriai pull az 1D $\theta_*$ prioron: $-2092.09$, nem likelihood-szignifikancia) és BAO $\chi^2=2201.29927$: nagyon rossz akusztikus geometriai benchmark. (A felülírt, rögzített M4-horizontos értékek: $0.41368156$ és $2209.46156$.)

### 8.2 A $67.36\,(1+1/4\pi)$ numerikus egybeesés

$67.36\,(1+1/4\pi)=72.7203\;\text{km/s/Mpc}$ és $67.36\,e^{1/4\pi}=72.9394\;\text{km/s/Mpc}$ helyes számtani azonosságok, de dinamikai $H_0$-korrekció nem következik belőlük: az órafaktort a megoldott $|\phi|$ potenciál adja (1. és 7. fejezet), nem a puszta $\beta_{\rm Yuk}$ csatolás.

### 8.3 Inverz távolságlétra profilozott statisztikával

| Forgatókönyv | $H_0$ | $\Omega_m$ | $\Omega_K$ | $r_d$ [Mpc] | $P=H_0r_d$ [km/s] | csak $\chi^2_{\rm BAO}$ (12×12) | $\Delta r_d$ Planckhoz [Mpc] |
|:------------------------------------|-----:|------:|------:|----------:|----------:|-------:|-------:|
| 1. Rögzített Step 25 M4 alapág | 67.40 | 0.3148 | 0.0000 | 147.09 | 9914.02 | 20.785 | +0.00 |
| 2. M4-alak, $P$ profilozva | 72.80 | 0.3148 | 0.0000 | 137.77 $\pm$ 0.65 | 10029.65 $\pm$ 47.43 | 14.701 | -9.32 |
| 3. $\Omega_m$ és $P$ profilozva | 72.80 | 0.2939 | 0.0000 | 140.02 $\pm$ 1.74 | 10193.67 $\pm$ 126.47 | 12.741 | -7.07 |
| 4. Közös BAO + Planck $r_d$ + feltételezett lokális $H_0$-prior | 70.81 | 0.2719 | 0.0000 | 147.03 (prior-$\sigma$ 0.26) | 10410.24 | 15.674 | -0.06 |
| 5. Merev, minden $z$-n $H_0=72.8$, változatlan $r_d$ | 72.80 | 0.3148 | 0.0000 | 147.09 | 10708.32 | 194.342 | +0.00 |
| 6. Rögzített-$P$ szemléltetés | 72.80 | 0.3148 | 0.0000 | 136.18 | 9914.02 | 20.785 | -10.91 |
| 7. Sík 72.8, $\Omega_m=0.2448$, újraint. $r_d$ | 72.80 | 0.2448 | 0.0000 | 150.65 | 10967.41 | 70.304 | +3.56 |
| 8. Nyílt 72.8, rögzített $\omega_{cb}$ | 72.80 | 0.2698 | 0.0098 | 147.09 | 10708.27 | 61.206 | +0.00 |

- A $\chi^2_{\rm BAO}$ oszlop csak a DESI $12\times12$ tagot tartalmazza. A 4. sor a priorokat is hozzáadja: $\chi^2_{\rm közös}=\chi^2_{\rm BAO}+\chi^2_{r_d\,\rm prior}+\chi^2_{H_0\,\rm prior}=15.67431585+0.06234684+3.97660606=19.71326876$; $r_d$-bejegyzése a **bemeneti prior hibáját** ($\pm0.26$ Mpc) közli, nem illesztett bizonytalanságot. Az 5. sor a minden $z$-n merev átskálázás változatlan $r_d$-vel, nem késői átmeneti modell (3. fejezet).
- A rögzített-$P$ szemléltetés (6. sor: $P=9914.02\;\text{km/s}$, $r_d=136.18\;\text{Mpc}$) feltételes számtani példa; a $P=9914.02$ nem a BAO legjobb illesztett szorzata, korábbi szigma-címkéjét visszavontuk.
- A profilozott BAO-igény $H_0=72.8$ mellett $r_d=140.02\pm1.74\;\text{Mpc}$ (kétparaméteres kontroll), a közös BAO + Planck $r_d$ + feltételezett lokális $H_0$ benchmark-prior kontroll pedig $\Delta\chi^2=6.973$ értéket ad.
- A változatlan alakú és $r_d$-jű merev átskálázást a BAO erősen elutasítja ($\chi^2=194.342$, $p=4.79\times 10^{-35}$).

![Inverz távolságlétra a $(H_0,r_d)$ síkon (`step2_3_inverse_ladder_rd_H0_exclusion.png`)](step2_3_inverse_ladder_rd_H0_exclusion.png){ width=96% }

---

## 9. Mi marad meg, és mi nem igazolt

**Megmarad:**

- A rögzített $E(z)$ alak és változatlan $r_d$ mellett minden vöröseltolódáson 8%-kal megemelt $H(z)$ nagyon rosszul illik a BAO-hoz.
- A puszta $\beta_{\rm Yuk}$ nem egyenlő a forrásból megoldott skalárpotenciállal vagy egy órakalibrációs eltolódással.
- A rövid hatótávú, adiabatikus, elhanyagolható oszcillációs energiájú ágban a késői kozmológiai skalárkorrekció rendkívül kicsi.
- A konkrét M2 normalizálás nem fizikai, az M3 görbülete pedig rossz akusztikus geometriai benchmark.
- A profilozott konstans SNe-nullpont-eltolás és a $H_0$–$M_B$ degeneráció geometriai szemléltetése helyes.

**Nem igazolt ebben a csomagban:**

- közös Planck + Pantheon+/Union3 + RSD + BBN adatvalidáció;
- a korábbi szigma-címkék által sugallt általános fizikai kizárási erő;
- a skalár–tenzor elmélet kozmológiai lezárása tetszőleges megengedett kezdeti feltételekkel;
- új, TAL-specifikus empirikus előny a megfelelően illesztett $\Lambda$CDM-hez képest.

**Zárómondat.** A vizsgált L=10 km-es, rögzített csatolású, elhanyagolható szabad skaláros oszcillációs energiájú M4 ágban a lokális Yukawa-faktor nem szolgáltat levezetett mechanizmust a 72.8 km/s/Mpc globális Hubble-állandóhoz. A rögzített korai háttér és rögzített késői alak egyszerű globális átskálázását a BAO erősen elutasítja. Az eredmény nem átfogó, közös Planck–SNe–RSD–BBN likelihood-elemzés, és nem az elmélet összes kozmológiai ágára vonatkozó no-go tétel.

### 9.1 Nyitott korlátok / Open limitations

| Pont | Nyitott korlát | Állapot ebben a csomagban |
|:-----|:-----------------------------|:-----------------------------|
| P01 | A $z=0.845$ ELG1 pont ($18.33\pm0.62$) igazolt elsődleges forrása | igazolt elsődleges forrás nélküli, választott illusztráció; nem DESI-mérés; nincs a 12 pontos likelihoodban |
| P02 | A $10^5\,({\rm D/H})_{\rm std}=2.5147$ egyedi forrása | számtani proxy; empirikus standardmodell-teszthez igazolt bemenet és elméleti hiba kell |
| P03 | Csatolt skalár–Friedmann–anyag kezdetiérték-probléma és szabad homogén amplitúdók | a $C=0$ ág expliciten választott; nincs állítás minden kezdeti feltételről vagy minden kozmológiai TAL-ágról |
| P04 | Planck TT/TE/EE + lensing likelihood és rekombinációs számítás | a 3D prior feltételezett, $z_*$ és $z_{\rm drag}$ rögzített kontrollparaméter; geometriai benchmark, nem CAMB/CLASS-validáció |
| P05 | Valódi Pantheon+/Union3 adatvektor és kovariancia | szintetikus 18 csomópontos metrika, nem megfigyelési modellpreferencia |
| P06 | RSD-kovariancia, átfedések, AP-leképezés, levezetett amplitúdó, $k$-spektrum | $k=0.1\,h$/Mpc variáció rögzített $\sigma_{8,0}$-val; az $1.666\times10^{-40}$ erre az ODE-re, normalizálásra és kezdeti feltételekre vonatkozik |
| P07 | BBN termikus trace, változó $g_*$, $e^\pm$-termodinamika, reakcióhálózat | három mintaepocha feltételezett érzékenységekkel; nem elemgyakoriság-likelihood |
| P08 | Valódi nap- és földi sűrűségprofil, felszíni Yukawa-megoldás, a lokális órák távolságlétrába illesztése | bulk és félteres limit; a $\Delta H_0^{\rm app}$ táblázat számtani szemléltetés |
| P09 | A V26 könyv egésze, minden történeti konvenció és minden csillag-BVP | ez a bírálat a 26. lépésre vonatkozik; a V26 normalizálási azonossága csak célzott ellenőrzés |

---

## 10. Fájlok és reprodukálhatóság

- `desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`: 376 B, SHA-256 `dd2873a0b88459a491af3c0c0307ba059f62df9211d5b976760f310565a1be68`
- `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`: 2,170 B, SHA-256 `bbafa9074b51cf1a45e0d10e4f37db8c0e80a5d1d1788857abb7fc49fb21abcc`
- `tal_cosmology_common.py`: 24,847 B, SHA-256 `54c5636b35631af314f501060e39adf723a46cada78db15526b0591a87c2c109`
- `tal_hubble_desi_bao_summary.json`: 52,790 B, SHA-256 `94ba04e4c09fc8a7d6a94eed04eecf03c9f794444a3a9f646a58f80aa1e3ca53`

- Generátorok: `tal_cosmology_common.py` (közös háttér), `run_step1_1.py`, `run_step1_2_and_1_3.py`, `run_step1_4_to_1_6.py`, `run_step2_nogo_and_report.py`, `run_step26_revision_controls.py` (újrabírálati kontrollok), `hubble_desi_bao_tal_audit.py`; szöveg: `generate_step26_markdown.py` + `step26_md_templates.py`, szöveges regressziós teszt `test_step26_text_regression.py`; PDF: `compile_step26_clean_pdfs.py`; futtató: `run_all.py` (sorrend és adat-hash ellenőrzés, `--skip-block1`), függőségek: `requirements.txt`.
- A csomag önállóan reprodukálható: mindkét DESI bemeneti fájl, az összes CSV- és PNG-kimenet, a szigorú JSON (`allow_nan=False`, a NaN `null`-ként írva) és a `SHA256SUMS.txt` benne van.
- Külső források: DESI VI (arXiv:2404.03002), DESI III (arXiv:2404.03000), Planck 2018 VI (arXiv:1807.06209), Cobaya `bao_data` tároló (DESI 2024 Gauss BAO középérték és kovariancia).

---

## 11. Revíziós napló

| Pont | Változtatás | Hol |
|:-----|:-------|:------|
| R1 | A kinetikus normalizálás egyeztetése a V26-tal (fájl:sor); $\phi_{\rm régi}=2\phi_{\rm új}$ esetén $\alpha_{\rm új}=2\alpha_{\rm régi}$ szükséges; az átvihető eredmények felsorolva; a kompaktcsillag-$\lvert\phi\rvert$ külső V21 belső értékként jelölve; a $2\alpha^2\beta$ új címkéje érzékenységi kontroll | 1.5, 1.7, 7 |
| R2 | Szakaszos $H=sH_{\rm M4}$ ($z<z_t$) kontroll $z_t\in\{2.5,10,100,1000\}$ értékekre; a $D_M/s$ sor inkonzisztens geometriai szemléltetésként átsorolva | 3, `step26_r2_piecewise_transition.csv` |
| R3 | Önállóan reprodukálható csomag: DESI-bemenetek hash-sel, minden CSV/PNG, `requirements.txt`, `run_all.py`, szigorú JSON; a régi `build_step26_pdfs.py` kikerült a csomagból | 10, README |
| R4 | A közös $\chi^2$ felbontása BAO-, $r_d$-prior- és $H_0$-prior-tagra; a táblázat oszlopa csak BAO; az $r_d$ $\pm0.26$ bemeneti prior hibaként jelölve | 2, 8.3 |
| R5 | Szimmetrikus nullpont-eltolás az M4 és az alacsony $\Omega_m$-ű modell számára; a legalsó csomópont értéke elválasztva az $5\log_{10}$ azonosságtól | 4, `step26_r5_sne_symmetric_offsets.csv` |
| R6 | A BBN óraproxy és együtthatók feltételezett proxyként jelölve; keret és előjel tárgyalva; források megadva | 6 |
| 5.1 | A float64-ág $1+\Delta G_{\rm eff,kan}/G_N$ alakot használ; abszolútérték-fejlécek; a független átskálázott ellenőrzés idézve | 5 |
| 5.2 | A $z=0$ sodródás elválasztva a Fermi–Dirac interpolációs hibától, amelyet kiszámoltunk | 1.6 |
| 5.3 | A $\theta_*$-hez illesztett modell átcímkézve: $\Omega_m=0.24419$, $\omega_{cb}=0.12877$; az inverz-létra ábra $\sqrt{\Delta\chi^2}$ megfogalmazást használ; a provenienciaállítás átfogalmazva | 3, 8.3, fejléc |


---

## 12. Revíziós napló 2 / Revision log 2 (második bírálat, F01–F12)

| Pont | Változtatás | Célzott ellenőrzés |
|:-----|:-----------------------------|:-----------------------------|
| F01 | $G_J=GA^{-2}$ helyett $1/G_T=A^{-2}/G_*$, $G_T=G_*A^2$; $\phi<0$ esetén $G_T$ csökken; $G_T$, $G_{\rm lab}$ és a csupasz $2\lvert\phi\rvert$ proxy szétválasztva (6. fejezet, kódkommentek, R6-metaadat) | a szöveg és a JSON `R6_bbn_proxy_provenance.clock_proxy` mezője az új összefüggést tartalmazza |
| F02 | G-szótár az 1.5 fejezet után ($G_*$, $G_T$, $G_{\rm lab}\approx G_*(1+q)$, $G_N:=G_*$) | $1/(1+q)=0.9262883$, számolva a `run_step26_revision_controls.py`-ban |
| F03 | Az M3 saját $r_s$, $r_d$ értékeket használ | $r_s=144.36622618$, $r_d=147.01159289$ Mpc, $100\theta_*=0.41346227$, $\chi^2_{\rm BAO}=2201.29927$ (ellenőrzés a bírálati értékekkel szemben) |
| F04 | A step2 státuszszövegek és a CLI-szöveg átírva: „BAO fit”, „deviation on illustrative prior”, „pull vs fixed target”; nincs Planck-kizárási vagy görbületi-likelihood megfogalmazás | automatikus szövegellenőrzés a `run_step26_revision_controls.py`-ban |
| F05 | A $72.8\pm1.0$ feltételezett lokális $H_0$ benchmark-priorként címkézve; JSON-kulcs `H0_local_benchmark_prior_km_s_Mpc`; ábracímkék módosítva | a címkék már nem tulajdonítják az értéket a SH0ES-nek (CSV/JSON/PNG szöveg) |
| F06 | Az 1.5 fejezet mondata a $2\alpha_Y^2\beta_{\rm Yuk}$ sorról | szövegellenőrzés |
| F07 | Az RSD-jelmagyarázat a kanonikus maximumot mutatja | $1.6664\times 10^{-40}$ a jelmagyarázatban és a JSON-ban |
| F08 | A legalsó csomópontra vonatkozó mondat törölve a 4. fejezetből, az azonosság megmarad | szövegellenőrzés |
| F09 | SNe-ábra: az alacsony $\Omega_m$-ű eltolt görbe hozzáadva, „profilozott konstans nullpont-eltolás”, $N=18$, $\nu=17$ | alacsony $\Omega_m$ eltolással $\chi^2=22.777527$ |
| F10 | A sor új neve „V21 top-hat belső külső benchmark”, a peremoszlop n.a., proveniencia hozzáadva | SLy $r=2.0046409680935633$ km, $x=0.041239199439254325$; APR4 $r=2.75417748649486$ km, $x=0.041233030079658145$ |
| F11 | A PDF-építő és a `run_all.py --pdf` nem nulla kóddal lép ki pandoc/pdflatex vagy újonnan írt PDF hiányában; a `--tex-only` mód külön státuszt ad | teszt pdflatex nélküli `PATH`-szal |
| F12 | `theta_star_100_rs_only_shifted_fixed_DM` (átnevezve) és `theta_star_100_both_endpoints_shifted` (új) | alacsony $\Omega_m$: $1.04248582$ és $1.04250287$ |

## 13. Revíziós napló 3 / Revision log 3 (harmadik újraaudit, N01, N02, T01)

| Pont | Változtatás | Célzott ellenőrzés |
|:-----|:-----------------------------|:-----------------------------|
| N01 | 7. fejezet táblázata: a homogén FLRW-, a kozmikus üreg- és a Tejútrendszer-korong sor „$R\gg L$ perem” cellája n.a. (nincs peremgeometria); a `step1_6_local_yukawa_environments.csv` és a JSON az `abs_phi_surface_R_gg_L` és a `Delta_H0_apparent_surf_km_s_Mpc` mezőben null értéket, valamint új, soronkénti `boundary_status` mezőt tartalmaz; a bulk értékek változatlanok (nincsenek felezve) | ellenőrzés a `run_step26_revision_controls.py`-ban (JSON `review3_N01_N02_T01`); szöveges regressziós teszt |
| N02 | 1.7 fejezet, „Változatlanul átvihető”: a zárójel szétválasztja a háttér skalárenergia-korrekcióját ($\Omega_\phi\sim10^{-47}$) és az RSD közvetlen Yukawa-perturbációját ($\delta G_{\rm eff}/G_*\approx3.8\times 10^{-41}$, $\delta f\sigma_8\approx1.67\times 10^{-40}$, 5. fejezet); a számok a JSON-ból származnak | szöveges regressziós teszt (HU és EN mondat) |
| T01 | Új `test_step26_text_regression.py`, amelyet a `run_all.py` a Markdown-generálás után futtat, hiba esetén nem nulla kilépési kóddal: ellenőrzi a HU/EN F06-mondatot, az F08 azonosság megfogalmazását, a felülírt hibás mondatok hiányát (V26 kétszeres növekedési törvényének tulajdonítás; az azonossági érték a táblázat legalsó csomópontjához rendelve), az N01 táblázatcelláit, az N02 mondatot és az RSD-ábra forrásának kanonikus jelmagyarázat-szövegét (F07); a `run_step26_revision_controls.py` kommentje frissítve | mutációs teszt scratch-másolaton: egy hibás mondat visszaírása nem nulla kilépési kódot ad |
