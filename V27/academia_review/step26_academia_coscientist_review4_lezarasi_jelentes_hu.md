# V27 / Step26 — a negyedik javított kiadás lezárási jelentése

**Vizsgáló:** Academia Co-scientist · **Dátum:** 2026-10-29  
**Vizsgált feltöltés:** 09:39 UTC, `V27_cosmology_desi_bao_audit (3).zip`, az EN/HU `(4).pdf` kísérőkkel.  
**Előzmény:** a 09:18-as csomag harmadik revíziós bírálata. Az eredeti és a korábbi jelentések változatlanok maradtak.

## 1. Bírálói döntés

**Az N01, N02 és T01 pontokat lezártnak tekintem. A korábbi F01–F12 lezárása érvényben marad.** A javítások a forráskódban, a gépi exportokban és a két nyelvi PDF-ben is megjelennek. A teljes, kihagyás nélküli újrafuttatás sikeres; a regressziós teszt hibafelismerését szándékosan elrontott, elkülönített példányokon is ellenőriztem.

**Ebben a verzióspecifikus auditban új, igazolt, további javítást igénylő hibát nem találtam. Ehhez a javítási körhöz nem jelölök meg új kötelező módosítást.** Nem indokolt a már javított kifogások újbóli fenntartása vagy újabb kifogások keresése pusztán a bírálati lista meghosszabbítására.

A közlemény a saját, kifejezetten korlátozott hatókörében **reprodukálható benchmark-/auditanyagként védhető**. Ez bírálói értékelés, nem folyóirati elfogadási garancia, a teljes TAL validációja vagy minden lehetséges hiba kizárása.

## 2. A három utolsó pont tényleges lezárása

| Pont | Új állapot | Ellenőrzött bizonyíték |
|:--|:--|:--|
| **N01 — bulk versus perem** | **Lezárt** | FLRW, üreg és Tejútrendszer-korong: a peremoszlop a két PDF-ben n.a.; a CSV-ben üres, a JSON-ban null; a perem-ΔH₀ is hiányzik. Új `boundary_status` mező. A helyes bulk-értékeket nem felezték meg. |
| **N02 — eltérő skalárhatások** | **Lezárt** | Az 1.7 fejezet HU/EN mondata külön adja az Ωφ∼10⁻⁴⁷ háttérenergiát és a közvetlen RSD Yukawa-tagot: δG_eff/G_*≈3.8×10⁻⁴¹, δfσ₈≈1.67×10⁻⁴⁰. Az új számok JSON-ból származnak. |
| **T01 — valódi regressziós ellenőrzés** | **Lezárt** | Új `test_step26_text_regression.py`; a `run_all.py` a Markdown-generálás után ténylegesen futtatja. Az eredeti állapot sikeres; 16 szándékos regressziót mind felismer; hibás sablon esetén a teljes runner is nemnulla kóddal leáll. |

### N01: ellenőrzött exportok

| Sor | abs(φ)_bulk, változatlan | Peremérték | Perem-ΔH₀ |
|:--|--:|:--|:--|
| Homogén FLRW | 1.994897054900984×10⁻⁴⁶ | null / n.a. | null |
| Kozmikus üreg | 3.989794109801969×10⁻⁴⁷ | null / n.a. | null |
| Tejútrendszer-korong | 7.426160269118664×10⁻⁴⁰ | null / n.a. | null |

A CSV szabványos üres cellái nem szó szerinti JSON-`null` tokenek: az olvasó hiányzó értékként kezeli őket, a JSON pedig valódi nullt tartalmaz. A Nap/Föld/fehér törpe félteres kontrollértékei nem változtak. A külső V21 belső benchmark továbbra is külön, felszíni érték nélkül szerepel.

A HU és EN PDF 9. oldalát képként is ellenőriztem; nem csak a forrás-Markdown változott meg. A 13. fejezet revíziósnaplója követi a kód/export javításait.

### N02: a javított mondat

Az 1.7 fejezet szó szerinti, ellenőrzött részlete:

> - Változatlanul átvihető: minden háttér-, BAO-, CMB-, SNe-alak- és RSD-szám (a háttér skalárenergia-korrekciója $\Omega_\phi\sim10^{-47}$; az RSD közvetlen Yukawa-perturbációja külön, az 5. fejezetben számszerűsített elhanyagolható tag: $\delta G_{\rm eff}/G_*\approx3.8\times 10^{-41}$, $\delta f\sigma_8\approx1.67\times 10^{-40}$), $\bar\phi_{\rm attr}(0)$, $\Omega_{\phi,\rm attr}(0)$ és a bulk lokális benchmarkok $\alpha_Y=1$ mellett.

A PDF-ek 4. oldalán a különbségtétel valóban olvasható. Ez a korábbi pontatlan indoklást kijavítja; nem változtatja meg a kinetikus normalizálás lezárt státuszát. **A publikált V26 és a Step26 kinetikus együtthatója α_Y=1 mellett egyezik.** Az elavult négyszeres kifogást és az arra épített V26-visszavonási javaslatot nem hozom vissza.

### T01: nem csak pozitív teszt

A hibátlan teljes futásban a szerző szöveges ellenőrzője **77 PASS ellenőrzést** adott. Ezután külön scratch-másolatban, minden teszt után visszaállított fájlokon vizsgáltam:

| Szándékos hiba | Tesztesetek | Eredmény |
|:--|--:|:--|
| F06: a szükséges tagadás törlése HU/EN | 2 | mindkettőnél 1-es hibakód |
| F06: helyes mondat mellett ellentmondó V26-tulajdonítás beszúrása HU/EN | 2 | mindkettőnél 1-es hibakód |
| F08: a 0.167357 érték újbóli első/legalsó csomóponthoz rendelése HU/EN | 2 | mindkettőnél 1-es hibakód |
| F08: az azonos E(z) azonossági feltétel elrontása HU/EN | 2 | mindkettőnél 1-es hibakód |
| N01: bulk visszamásolása a HU/EN peremcellába | 2 | mindkettőnél 1-es hibakód |
| N01: CSV/JSON hiányzó perem visszaírása számmá | 2 | mindkettőnél 1-es hibakód |
| N02: a két skalárhatás megkülönböztetésének elrontása HU/EN | 2 | mindkettőnél 1-es hibakód |
| F07: nem kanonikus jelmagyarázat vagy átadott maximum | 2 | mindkettőnél 1-es hibakód |
| **Összesen** | **16** | **16/16 felismerve** |

Egy további külön példányban az **angol sablonból** töröltem az F06 tagadását, és a teljes `run_all.py`-t elindítottam. A numerikus rész lefutott, a hibás Markdown újragenerálódott, a szöveges teszt megbukott, és a runner **1-es kilépési kóddal továbbadta a hibát**. Vagyis nem csupán egy manuálisan elindítható, a fő futási úttól független teszt került be.

A mutációk nem érintették a feltöltött eredetiket vagy a fő reprodukciós példányt. A teszt a konkrét revíziós tételek regresszióit fogja; nem univerzális tudományos szövegértelmező. Ez nem fennmaradó T01-kifogás.

## 3. Teljes reprodukció és csomagépség

**ZIP SHA-256:** `e08c8dfcf5454c318762a9d154ad219ad2fc0b4907b465c6dafc8f3638c58150`  
**Méret:** 5 617 462 bájt; **46 fájl**, hibátlan ZIP CRC.

| Ellenőrzés | Eredmény |
|:--|:--|
| Mindkét SHA256 manifest | **44/44 helyes** mindkettőben |
| Teljes `python3 run_all.py` | **0 kilépési kód**, körülbelül 9.9 s; **nincs `--skip-block1`** |
| Források | Mind a **12 Python-fájl változatlan** maradt a fő futásban |
| Környezet | Python 3.12.3; NumPy 2.5.4; SciPy 1.18.1; pandas 3.1.0; Matplotlib 3.11.2 |
| JSON | Szigorú JSON; **719 véges numerikus mező**, szerkezeti/nem numerikus eltérés nélkül |
| Legnagyobb abszolút JSON-eltérés | **6.4153×10⁻⁹**, az RSD χ²-ben |
| Legnagyobb relatív JSON-eltérés | **4.8635×10⁻⁶**, a ∼10⁻⁴⁰-es variációs maximumokban |
| CSV-k | **13 tábla, 1676 véges numerikus cella**; azonos alak, oszlopok, hiányzóérték-mintázat és nem numerikus mezők |
| Legnagyobb CSV-eltérés | Abszolút **3.5242×10⁻⁹**; relatív **6.7181×10⁻⁶** a tiny variációs értékekben |
| Markdown | A két főszöveg és másolataik újragenerálhatók; csak a JSON méret/hash provenienciasora változik |
| Külön PDF versus ZIP PDF | Mindkét nyelven **bájtszintű azonosság** |
| PDF-be ágyazott ábrák | **8/8 dekódolt RGB-kép pontosan egyezik** a most regenerált PNG-vel |

A 719 mező és az 1676 CSV-cella a **mostani** sémára vonatkozik, nem a korábbi kiadásból átvett darabszám. A kis lebegőpontos eltérések nem változtatják meg a közölt pontosságú eredményeket. Az újragenerált JSON provenienciasora értelemszerűen az új bájtokra vonatkozik; ettől a feltöltött eredeti manifestje nem hibás.

### PDF-fordítási státusz

- `run_all.py --pdf`: **2-es hibakód**, hiányzó pdflatex jelzésével.
- `compile_step26_clean_pdfs.py`: **2-es hibakód**, nincs hamis PDF-siker.
- `compile_step26_clean_pdfs.py --tex-only`: **0**, a két TeX elkészül, explicit **TEX-ONLY** státusszal.

**Friss PDF-fordítást nem végeztem**, mert pdflatex nincs ebben a környezetben. A beküldött PDF-ek javított oldalait és a regenerált ábrákkal való azonosságát ellenőriztem; ez nem egyenértékű egy friss LaTeX-build tesztjével.

## 4. A korábbi lényegi javítások nem romlottak vissza

A módosítások a helyi peremmezőket, a szöveggenerálást és a regressziós ellenőrzést érintik; a közös kozmológiai háttérmodul és a BAO/CMB/SNe/M3 fizikai számítási moduljai az előző kiadással bájtszinten azonosak. A korábbi F01–F12 forrásbizonyítékok az új fájlokban is megtalálhatók. A fő numerikus kontrollokat ismét elvégeztem:

| Független kontroll | Új audit eredménye |
|:--|:--|
| M3, saját a-változós hanghorizont-integrál | r_s=144.366226178967 Mpc; r_d=147.011592891633 Mpc; 100θ*=0.413462266606400; BAO χ²=2201.2992668914 |
| Külön újraoptimalizált BAO 2P kontroll | Ω_m=0.29389245163; P=10193.66994265 km/s; χ²=12.7405324028244 |
| Külön újraoptimalizált közös kontroll | Ω_m=0.27194920532; H₀=70.8058570681; r_d=147.0250796694 Mpc; χ²=19.7132687558893 |
| Közös χ² felbontása a szerző optimumán | 15.6743158507356 + 0.0623468417346 + 3.9766060634193 = 19.7132687558895 |
| SNe profilozás, M4 / low-Ω_m | 0.9166139641 / 22.7775266770; szimmetrikus konstans kezelés |
| Skálázott variációs növekedés, 3×10⁻¹⁴ tolerancia | max abs(δfσ₈)=1.66636815207650×10⁻⁴⁰; a toleranciateszt egyezik |
| F12: low-Ω_m, csak r_s / mindkét végpont | 1.0424858198503 / 1.0425028691050 |
| Feltételes G-kalibráció | G_eff/G_lab=0.926288317750839; nem kiszámolt 7.37%-os fσ₈-eltérés |

A késői darabos átmenet és a neutrínó FD/interpoláció ellenőrzése is megismétlődött. Az ellenőrző kód több helyen a szerző deklarált háttérkonstansait használja, más képletekkel/integrálási változóval; nem állítok új, teljesen független Boltzmann-megvalósítást.

## 5. Mi marad korlát, nem javítatlan bírálati hiba?

A 9.1 fejezet P01–P09 tételei továbbra is érvényesek, és a kézirat ezeket nem tekinti elvégzett adatvalidációnak. A vizsgálat C=0 adiabatikus ág, L=10 km és α_Y=1 mellett történik. A növekedési tiny korrekció azonos csupasz G_* melletti formális összevetés, nem teljes, G_lab-ra újrakalibrált kozmológiai illesztés. A szintetikus SNe-kontroll, a diagonális RSD-metrika, a feltételezett CMB-prior és az algebrai BBN nem egy közös megfigyelési likelihood.

Nem végeztem új csillag-BVP-t, a V26 teljes könyvének újraauditját vagy Zenodo–GitHub bájtazonossági ellenőrzést. Ezek hiánya a jelen, helyesen címkézett benchmarkközleményben nem indok arra, hogy a lezárt hibákat nyitottnak minősítsem.

**Végső minősítés:** a jelzett bírálati javítási kör **lezárható**. A mostani, pontosan azonosított fájlokra nem tartok fenn további bizonyított javítási követelményt. Erősebb kozmológiai vagy empirikus állításokhoz külön vizsgálat szükséges; a jelen eredmény nem az elmélet minden ágára vonatkozó no-go tétel.

A kísérő JSON tartalmazza a bemeneti hash-eket, a teljes futás naplóját, a mezőnkénti numerikus összevetést, a független kontrollokat, a 16 mutáció eredményét, a teljes runner hibaterjedési próbáját és az ellenőrző szkriptek szövegét. A korábbi jelentések történeti állapotot dokumentálnak; a 09:39-es feltöltésre ez a jelentés az irányadó.
