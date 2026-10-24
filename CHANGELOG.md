# Changelog

## V22 Book Edition (Zenodo V22) — 2026-10-24

Zenodo V22 Book Edition: [10.5281/zenodo.23604423](https://doi.org/10.5281/zenodo.23604423)  
Previous version: Zenodo V21 (= GitHub V19/V21, [10.5281/zenodo.23602501](https://doi.org/10.5281/zenodo.23602501)).

- Added [`V22_Book/`](V22_Book/): 15-page, 12-chapter (+ 2-page Round-3 technical reproduction audit appendix) research monograph synthesizing Versions 1–21 of the Time Attraction Law in English (`tal_book_EN_2026-10-24.pdf` & `.tex`, primary, byte-identical to Zenodo record `23604423`), Hungarian (`tal_book_HU_2026-10-24.pdf` & `.tex`, byte-identical to Zenodo record `23604423`), and German (`tal_book_DE_2026-10-24.pdf` & `.tex`), together with `generate_figures.py`, 12 figures in `figures/`, and `SHA256SUMS.txt` (21 verified entries).
- Updated `README.md` and `CITATION.cff` with DOI `10.5281/zenodo.23604423`.

## V18 (GitHub) — 2026-10-21

GitHub V18 = Zenodo V20 (10.5281/zenodo.23540741)

Előző verzió: V19 – 7th model (10.5281/zenodo.23497128).
A Zenodo V20 rekord négy fájlja változatlanul a `V18/zenodo_v20/` mappában van. A `V18/final/` a V20 jelentés-PDF-hez
képest a revision3 review két szövegjavítását (és a review szerinti két pontosítást) tartalmazza; a kód és a CSV/JSON fájlok
bájtra azonosak a V20 csomagéval.

### Változások a V19 óta (a Zenodo V20 jelentésben és a `V18/final`-ban egyaránt)
- **Integrálfolytonosság:** L ≥ 2.6 m-re közvetlen külső különbségi integrál r = L/t helyettesítéssel (32 pontos
  Gauss–Legendre), levágás és radiális interpoláció nélkül; a korábbi 0.01268 J/A²-es mesterséges eltolás megszűnt.
  Folytonossági teszt L = 2.99 / 3.00 / 3.01 m-nél: relatív eltérés 0 / −7.4·10⁻¹⁶ / −1.5·10⁻¹⁵; a 2.6–2.9 m-es
  átfedésben ≤ 2.9·10⁻¹⁵ (`csv/continuity_test.csv`).
- **Háttértér-kereszttag lehatárolása:** a §4 polaritásmondata kiegészítve a 2B_bg·B_coil kereszttaggal; a közölt
  −1.20×10⁻⁴⁸ érték az izolált páros tagra vonatkozik, nem minden lehetséges egyszerű ON/OFF mérés teljes előrejelzése.
- **Biztonsági mondat cseréje:** az Onderdonk-becslés csak idealizált adiabatikus olvadási idő; ebből biztonságos
  impulzushossz nem állapítható meg, a számított nagyáramú esetek nem üzemeltetési javaslatok.
- **291.7 mm mért menetátmérő, mérési jegyzőkönyvvel (1.2 pont):** 2025-11-13, tolómérő (0.02 mm felbontás), mérte
  Szabolcsi Dávid, 5 egymást követő leolvasás, mind 291.7 mm (átlag 291.7 mm, mintaszórás 0); a 300 mm a tekercstest
  pereme (külső méret). A térmodell a mért 291.7 mm-rel számol (B a fő órán 11.63 µT).
- **Konvergencia-közlés:** a szórás definíciója (a 4 futás mintaszórása, n − 1 nevezővel), teljes tartomány és
  mintaszórás táblázatban; a négy MC relatív hiba külön: 2.17 % / 0.99 % / 1.14 % / 4.04 % (L = 0.316 / 0.631 / 1.0 / 2.0 m).
- **„Főóra-egyenérték” megnevezés** a táblázatban (2.20·10⁻⁶ /T; 0.189 /T²).

### A revision3 review pontjai, a `V18/final`-ban alkalmazva (a Zenodo V20 PDF-ben még nem)
- **Távoli „0 szórás” megjegyzés** (konvergenciatáblázat, L = 2.6–100 m): „A négy futás ebben a tartományban azonos nlo
  és nt beállítást használ. A nulla mintaszórás ezért nem távoli hibabecslés; a távoli konvergenciát nlo és nt külön
  változtatásával kell vizsgálni.”
- **Hőkontroll átfogalmazása** (§4, 5. pont): a „Így a be/ki különbség nem tartalmaz hőmérsékleti lépést.” mondat helyett:
  „Az azonos összteljesítményű, kis maradékterű kontroll célja a kapcsolással korrelált termikus változás csökkentése.
  A be- és kikapcsolt állapotok hőmérsékleti azonosságát és a kontroll maradékterét méréssel ellenőrizni kell; a maradék
  eltérést a frekvenciabecslés bizonytalanságába kell foglalni.”
- **L = 0.95 m illesztési hiba elismerése:** a V20 PDF már közölte a teljes energiák 1.4·10⁻³ J/A²-es (0.06 %) eltérését
  és az ebből adódó ugrást; a `V18/final` a pontos értékeket is megadja: belső forma 2.17738947, külső forma
  2.17597607 J/A², illesztési hiba 0.00141340 J/A² (a belső érték 0.0649 %-a).
- **§5 első pontja:** a −1.202·10⁻⁴⁸ az „izolált, polaritásra páros tekercs-B² tag” számított differenciális TAL-jele.

### V20 kézirat (angol, Foundations of Physics)

- Hozzáadva: `V18/manuscript/tal_v20_fop.pdf` és `V18/manuscript/tal_v20_fop_manuscript.zip` (LaTeX-forrás, ábrák, szkript).
- A kézirat tartalmazza a revision3 review két szövegcseréjét (a távoli nulla mintaszórásról szóló megjegyzést és a hőkontroll-mondatot), és hivatkozik a `revision3_checks.json`-ra.
- A Zenodo-feltöltést a szerző maga végzi; ez a commit a Zenodóhoz nem nyúl.

### Ellenőrzések
- `V18/revision3_checks.json`: a szerző pipeline-ja generálta 2026-10-21-én (`V18/make_revision3_checks.py`), nem a bíráló
  eredeti fájlja (a bíráló `check_revision3.py` szkriptje nem áll rendelkezésre). Tartalma: 18/18 SHA-256 egyezés a
  `V18/final/SHA256SUMS.txt` alapján; B a fő és a referenciaórán független köráram-összeggel (8 × 2000 menet) vs. a csomag
  térfüggvénye (relatív eltérés −4.3·10⁻⁸, ill. −7.6·10⁻⁹; 11.6329 µT, ill. 0.43054 µT); célzott főcsúcs-újrafuttatás
  L = 0.630957 m-nél (−1.2017·10⁻⁴⁸, egyezik a JSON-nal); külső integrálok L = 2.6 / 3.0 / 10 / 100 m-nél, nlo és nt
  együttes finomításával (6000 → 12000, 32 → 64; relatív változás −9.1·10⁻⁷ … −3.0·10⁻⁴); az L = 0.95 m-es maradék
  számtana (0.00141340 J/A², 0.0649 %, egyezik a teljes energiák különbségével).

### Repó-szintű változások
- `LICENSE`: a Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) hivatalos jogi szövege.
- `CITATION.cff`: érvényes CFF formátum, licenc `CC-BY-NC-4.0`, DOI 10.5281/zenodo.23540741 (Zenodo V20 / GitHub V18),
  javított GitHub URL (`unbekanntedavid-glitch`).
- Nem része a csomagnak: NIST SP 1065 PDF.
