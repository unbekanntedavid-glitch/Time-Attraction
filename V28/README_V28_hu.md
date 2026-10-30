# Time Attraction Law — skalár-tenzor Yukawa-illesztés, V28

**Szerző:** Szabolcsi Dávid (ORCID: 0009-0005-6300-074X)
**Dátum:** 2026-10-30 (r3 átdolgozás)
**Státusz:** előzetes kutatási csomag. A 27. lépés következtetése feltételes; nem ellenőrzött empirikus kizárás. A működési próbák (ellenőrzőösszegek, újrafuttatások, a 117 szimulátoros összevetés) nem jelentik minden szimulátorpanel fizikájának auditját.

## Kibontás (két ZIP-fájl)

Mindkét ZIP-et (`TAL_V28_r3_part1_step26.zip`, `TAL_V28_r3_part2_simulator_step27.zip`) ugyanabba a könyvtárba kell kibontani; egyetlen közös `TAL_V28/` mappába olvadnak össze. A fő `SHA256SUMS.txt` csak akkor érvényes, ha mindkét rész együtt van kibontva; az egyes részek a `SHA256SUMS_part1.txt` / `SHA256SUMS_part2.txt` listával külön is ellenőrizhetők. A `TAL_V28/` mappán belülről: `sha256sum -c SHA256SUMS.txt`. Részletek: `README_EXTRACT_hu.md`.

**DOI-k**

- V28: a Zenodo a közzétételkor adja ki (még nem ismert).
- A V28 a V27 új verziója: `10.5281/zenodo.23761574` (szigorúan V27, 26. lépés).
- V26: `10.5281/zenodo.23727585` (szigorúan V26, 1–25. lépés).
- V25: `10.5281/zenodo.23699685`; V24: `10.5281/zenodo.23680715`.

## A modell

Skalár-tenzor Yukawa-ág: L = λ = 10 km, α_Y = 1, β = 1/(4π), adiabatikus attraktor, C = 0. Az effektív lokális csatolás α_eff = q = α_Y²β = 0,0795775: pozitív, univerzális és árnyékolatlan. G-szótár: G_lab ≈ G_∞(1 + q), G_eff/G_lab (r ≫ L) = 0,9262883.

A V28 a modell három egymáshoz kapcsolódó vizsgálatát tartalmazza.

## 1. 26. lépés — kozmológiai audit a DESI DR1 BAO ellen (`step26_cosmology_desi_bao/`)

A V27 26. lépésének auditja, bájtra azonosan átvéve a közzétett V27-kiadásból (GitHub `V27/`, 48 fájl, a saját `SHA256SUMS.txt` 46/46 OK). A magyar és angol jelentés Academia-változata a `docs/` mappában van.

Kérdés: összeegyeztethető-e a merev TAL-háttér (H0 = 72,8 km/s/Mpc) a DESI DR1 BAO-távolságokkal (12×12-es kovariancia)?

- Merev háttér, H0 = 72,8: χ²_BAO = 194,342.
- M4-alak a saját újraintegrált horizontjával (r_s = 144,443 Mpc, r_d = 147,092 Mpc): χ²_BAO = 20,785.
- Kétparaméteres profil (Ω_m, P = H0·r_d): χ²_BAO = 12,7405; H0 = 72,8 mellett az adatok r_d = 140,02 ± 1,74 Mpc-et kívánnak.
- M3-kontroll a saját horizontjaival: r_s = 144,366 Mpc, r_d = 147,012 Mpc, χ²_BAO = 2201,3.
- Együttes BAO + Planck r_d-prior + feltételezett lokális H0-prior (72,8 ± 1,0): χ² = 19,713, Δχ² = 6,973 a csak BAO-optimumhoz képest.

Következtetés: a szükséges r_d ≈ 140 Mpc jóval a kontrollmodellek által önmagukra újraintegrált ≈147 Mpc-es horizont alatt van. A merev háttér ezért nem egyeztethető össze a DESI BAO-val, hacsak a korai univerzum hanghorizontja nem módosul. Újrafuttatás: `python3 run_all.py` (lásd a 26. lépés README-jét).

## 2. Interaktív szimulátor (`simulator/`)

Az `index.html` egyetlen, önálló kétnyelvű (HU/EN) oldal; böngészőben megnyitható, szerver nem kell hozzá (278 402 B, SHA-256 `46379f6323f7e9bdc342b95f0033d52c6c77fb5ffb1d466c676aabd69fae4226`). A V27-es GitHub Pages-változattól csak az alábbi pontforrás-jelölésben és a lábléc adatforrás-sorában tér el. Hat panel: Yukawa-erőtörvény és G, SPARC forgási görbék, GW170817 és lencsézés, Kozmológia V27 · DESI BAO, további próbák, info és linkek. Minden szám a böngészőben számolódik a csomag képleteiből. A `src/` mappában vannak a források és az összeállító szkript (`build.py`), a `shot_*.png` fájlok a panelek képernyőképei (a `shot_yukawa.png` magyarul, a `shot_yukawa_en.png` angolul mutatja az r2-állapotot; a `shot_footer_hu.png` / `shot_footer_en.png` a lábléc r3-as adatforrás-sorát).

A Yukawa-panelen a Nap / Föld / neutroncsillag kiválasztásakor a skalármező a pontforrásos φ(r) = −α_Y β GM/(c² r) e^(−r/L) képletből számolódik (`phiPoint`). A panel ezt láthatóan jelzi: pontforrásos, gyenge-téri szemléltetés; nem a test véges tömegeloszlásának integrált tere.

A beágyazott adatok újragenerálása (`src/build_data.py`) külön lépés, amelyhez külső V26-bemenetek kellenek (V26, doi:10.5281/zenodo.23727585); a bemeneti útvonalak hordozhatók (`--v27`, `--v26`, `--sparc`). Lásd: `simulator/README_DATA_hu.md`.

Reprodukálhatósági ellenőrzés: a `node verify_node.js` az oldal saját JavaScriptjével 117 mennyiséget számol újra, és összeveti a közzétett 26. lépéses JSON-értékekkel (kimenet: `verify_node_out.json`). 113 sor 10⁻⁶ alatti relatív eltéréssel egyezik. 4 sor 10⁻⁶ fölött van (legnagyobb 8,0×10⁻⁵, NGC 3198 C_req); ugyanaz a 4 sor, mint a korábban átnézett V27-es szimulátorállapotban (kerekített referenciaértékek, pl. NGC 3198 χ² = 32639,3, és optimalizálási tolerancia).

## 3. 27. lépés — a rögzített, árnyékolatlan lokális Yukawa-ág tesztje (`step27_local_yukawa_branch/`)

Kérdés: összeegyeztethető-e a rögzített α = 0,0795775 csatolás λ = 10 km-nél a geofizikai (óceáni és tornyos) gravitációs mérésekkel?

R a **határhoz viszonyított arány** = 0,0795775 / α_lim(λ = 10 km). Csak azonos potenciálkonvenció mellett értelmezhető; R > 1 azt jelenti, hogy a modell az adott határ megengedett tartományán kívül esik. R **nem** szigmában kifejezett szignifikancia. Az alábbi α_lim értékek mind saját újraszármaztatások idealizált réteg-/féltérmodellekben, az idézett hibákat Gauss-féle 1σ-ként értelmezve (95%-hoz ×1,96), következetesen rögzített G-kalibrációval (G_∞ = G_lab/(1 + α), a határok pontos inverzióval). Ezért minden R-érték **FELTÉTELES**.

| Rekord | α_lim (10 km) | R |
|---|---|---|
| Óceán (Zumberge 1991, March 2019 alapján), pontos inverzió — döntő | 0,00225606 | 35,27 |
| Óceán, March ±1× intervallum | 0,000910599 | 87,39 |
| Óceán, March (81) egyenlet nyomtatott alakban (forgatókönyv) | 0,00350120 | 22,73 |
| Torony: Thomas 1989 (ρ = 2000 / 2670 kg/m³) | 0,00650259 / 0,00486292 | 12,24 / 16,36 |
| Torony: Romaides 1994 | 0,00228048 / 0,00170725 | 34,90 / 46,61 |
| Torony: Romaides 1997 | 0,00066375 / 0,00049711 | 119,89 / 160,08 |
| Torony: BREN (OSTI 6190822) | 0,00755665 / 0,00564970 | 10,53 / 14,09 |

- Alkalmazott R-tartomány: 10,53–160,08. A tartomány forgatókönyveket fog össze (hibaértelmezések, rétegképletek, feltételezett sűrűségek); nem kombinált konfidenciaintervallum.
- Visszahelyettesítés: 11/11 teszt átment (minden α_lim a megengedett reziduum szélét adja vissza).
- Jekeli 1990: 10 km-es érték nincs alkalmazva (az absztraktból a λ-tartomány ismeretlen). Cruz 1991 (Erie): csak feltételes becslés (≈0,00088), a tartományban nem szerepel.

Ítélet: A saját újraszármaztatás erős feszültséget jelez a rögzített, árnyékolatlan Yukawa-ággal; a kizárás megerősítéséhez a mérési válasz és a bizonytalanságkezelés ellenőrzése szükséges. Az „erős feszültség” csak ebben az értelemben áll: az idealizált réteg-/féltérmodellekben, a feltételezett hibaértelmezés mellett nagy numerikus eltérés adódik. Nem konfidenciaszinttel rendelkező empirikus kizárás. Az érvényességi kör csak a rögzített, árnyékolatlan lokális ág, nem minden TAL-változat, és a 26. lépést nem nyitja újra.

Adatstátusz:

- A Zumberge 1991-es óceáni bemenetek másodlagos forrásból (March et al. 2019, PRD 100, 042002) átvett adatok, nem közvetlenül ellenőrzött primer adatok.
- Hiányzó primer források: Zumberge et al. 1991 (PRL 67, 3051; nem sikerült beszerezni), Cruz et al. 1991 (JGR 96, 20073; a cikk szövege nincs meg, csak az absztrakt), Jekeli, Eckhardt, Romaides 1990 (PRL 64, 1204; a cikk szövege nincs meg, csak az absztrakt). Részletek: `step27_yukawa_fixed_branch_test/step27_results.json` (`missing_primary_sources`) és `sources/README_sources.md`.
- Harmadik féltől származó dokumentumot nem közlünk újra: a `sources/README_sources.md` mindegyiket DOI-val, adatlap- és letöltési URL-lel (ellenőrizve) és a felhasznált oldalakkal/egyenletekkel adja meg. A harmadik féltől származó dokumentumok terjesztési jogát nem ellenőriztük; a dokumentumokat nem közöljük újra, hanem DOI-val és hivatkozással jelöljük meg. A `reference_inputs/` mappában az Adelberger 2003 4. ábrájának értékei saját digitalizálásunk, forrásmegjelöléssel; az ábrakivágat nincs a csomagban.

Nem tesztelt jelöltek (új ágak, nem igazolt kiutak): árnyékolás Föld-méretű forrásokra; α ≲ 10⁻³ 10 km-nél; más λ (pl. λ < 1 cm) változatlan α mellett, újratesztelendő.

Újrafuttatás: `python3 step27_yukawa_fixed_branch_test/run_step27.py` (numpy, matplotlib) újragenerálja a `step27_results.json` fájlt és a négy ábrát. A rész saját ellenőrzőösszeg-listája: `step27_local_yukawa_branch/SHA256SUMS.txt`. A bírálati előzmények (S01–S04): `step27_review1_changes_hu.md`.

## Adatforrások és licencek / Data sources and licenses

A csomag saját tartalmára a CC BY-NC 4.0 vonatkozik; a beágyazott külső adatok megtartják a saját feltételeiket. A források, a szó szerinti hivatkozási és köszönetnyilvánítási szövegek és a tisztázatlan pontok: `DATA_LICENSES.md`.

*EN: The author's own content is CC BY-NC 4.0; embedded third-party data keep their own terms. Sources, verbatim citation and acknowledgement texts and the unclear points: `DATA_LICENSES.md`.*

| Forrás / Source | Licenc / License or terms | Hivatkozás / Citation | Újraterjesztés / Redistribution |
|---|---|---|---|
| DESI DR1 BAO mean + covariance (`step26_cosmology_desi_bao/desi_2024_*.txt`, `simulator`) | CC BY 4.0 — https://data.desi.lbl.gov/doc/acknowledgments/ | DR1 paper arXiv:2503.14745; DESI 2024 III/VI (arXiv:2404.03000, 2404.03002); DESI acknowledgement text (verbatim in `DATA_LICENSES.md` and in the simulator footer) | allowed with attribution and change notice / hivatkozással és a változtatások jelzésével megengedett |
| SPARC rotation curves (`simulator`) | no explicit license stated / nincs kifejezett licenc — https://astroweb.case.edu/SPARC/ | Lelli, McGaugh & Schombert 2016, AJ 152, 157, doi:10.3847/0004-6256/152/6/157 | not stated; unclear / nincs megadva; tisztázatlan |
| RSD fσ8 (12 values) | no data license stated / nincs adatlicenc | 8 source papers (list in `DATA_LICENSES.md`) | not stated / nincs megadva |
| Planck 2018 parameters | no explicit data license; acknowledgement requested — https://www.cosmos.esa.int/web/planck/planck-data-use | Planck 2018 VI, doi:10.1051/0004-6361/201833910 + ESA Planck acknowledgement | not stated / nincs megadva |
| BBN: PDG 2024 Y_p; Cooke+2018 D/H | PDG 2024 (PRD 110, 030001): CC BY 4.0; Cooke+2018: no data license stated | doi:10.1103/PhysRevD.110.030001; doi:10.3847/1538-4357/aaab53 | PDG: allowed with attribution; Cooke: not stated |
| `PANTHEON_GRID` | own synthetic grid, not third-party data / saját szintetikus rács | — | — |

Tisztázatlan: a SPARC hivatalos oldala nem közöl licencet (csak hivatkozást kér); egy 2025-ös Zenodo-másolat (doi:10.5281/zenodo.16284118) CC BY 4.0 jelölésű, de nem igazoltan a szerzőktől származik.

## Változások az r3-ban (2026-10-30)

1. Új `DATA_LICENSES.md` (HU/EN): a beágyazott külső adatok (DESI DR1 BAO, SPARC, RSD fσ8, Planck 2018, BBN) licence vagy közölt feltételei, ellenőrzött URL-ekkel, szó szerinti hivatkozási és köszönetnyilvánítási szövegekkel és az újraterjesztési státusszal; ahol nincs kifejezett licenc, ezt kimondja.
2. Új szakasz ebben a README-ben: „Adatforrások és licencek / Data sources and licenses”; hivatkozás a `simulator/README_DATA_hu.md` / `_en.md` fájlokban.
3. Szimulátor: a láblécben látható HU/EN adatforrás-sor és lenyitható blokk a DESI és a Planck szó szerinti köszönetnyilvánításával. Az `index.html` újraépítve (278 402 B, SHA-256 `46379f6323f7e9bdc342b95f0033d52c6c77fb5ffb1d466c676aabd69fae4226`); a `verify_node.js` ugyanazt a 117 sort adja (113 10⁻⁶ alatt, ugyanaz a 4 felette). Számítás nem változott.
4. A 26. és a 27. lépés fájljai változatlanok.

## Változások az r2-ben (2026-10-30)

1. Az idegen PDF-ek (March 2019, OSTI 6190822, OSTI 7173840, Fischbach & Talmadge 1996, Adelberger et al. 2003) és az Adelberger-ábrakivágat kikerült; helyettük DOI/URL-hivatkozások a `step27_local_yukawa_branch/sources/README_sources.md` fájlban. A 27. lépés eredményei változatlanok (az újrafuttatás bájtra azonos `step27_results.json`-t ad). A 26. lépés részében csak a szerző saját PDF-jei vannak.
2. `simulator/src/build_data.py`: rögzített helyi útvonalak helyett hordozható bemeneti útvonalak (`--v27`, `--v26`, `--sparc`, `--out` vagy környezeti változók); az adatforrások leírása: `simulator/README_DATA_hu.md` / `_en.md`. Az újragenerálás bájtra visszaadja a mellékelt `data.js`-t.
3. Szimulátor, Yukawa-panel: látható magyar/angol jelölés, hogy a Nap / Föld / neutroncsillag értékek pontforrásos, gyenge-téri szemléltetések. Az `index.html` újraépítve; a `verify_node.js` ugyanazt a 117 sort adja (113 10⁻⁶ alatt, ugyanaz a 4 felette).
4. Kibontási útmutató (`README_EXTRACT_hu.md` / `_en.md`, mindkét ZIP-ben) és részenkénti ellenőrzőlisták: `SHA256SUMS_part1.txt` / `SHA256SUMS_part2.txt`.

## Integritás

A legfelső szintű `SHA256SUMS.txt` a `TAL_V28/` minden fájljának SHA-256-ját tartalmazza (mindkét rész). Ellenőrzés: `sha256sum -c SHA256SUMS.txt`.
