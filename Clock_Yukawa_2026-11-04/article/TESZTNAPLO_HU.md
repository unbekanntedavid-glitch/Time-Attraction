# Tesztnapló (közlemény és publikációs készlet, 2026-11-04)

Minden időpont Europe/Zurich (CET, UTC+1).

## Környezet
Linux 6.12.94+ x86_64, Python 3.13.5, NumPy 2.5.4, SciPy 1.18.1, matplotlib 3.11.2, Pillow 13.0.0,
pdfTeX 3.141592653-2.6-1.40.29 (TeX Live 2026), poppler-utils (pdftotext, pdfinfo, pdftoppm), 8 CPU.

## 1. Friss modulfuttatások (a közlemény bemenetei)
Mindhárom r4 modularchívumot külön üres mappába csomagoltuk ki, és egyszer lefuttattuk. Naplók: `logs/fresh_run_*.txt`
(az archívumok SHA-256 értékei és a parancsok a `logs/README.md` fájlban).

| modul | parancs | eredmény |
|---|---|---|
| óra (`tal_branch_clock_accel_test_r4.zip`) | `python3 run_check.py --out runs/fresh_article` | `RESULT: PASS`, kilépési kód 0 |
| BREN (`bren_continuation_zenodo_r4.zip`) | `python3 run_check.py` | `all steps passed`, kilépési kód 0 |
| szub-cm (`subcm_yukawa_window_zenodo_r4.zip`) | `python3 run_check.py` | `all steps passed (frozen-digitization mode)`, kilépési kód 0 |

## 2. A közlemény ellenőrzései
- `python3 src/check_claims.py`: 14/14 ellenőrzés sikeres, kilépési kód 0 (295 állítás újraolvasva az archívumokból; a
  generált fájlok bájtra azonosan újraelőállítva; nincs kézzel beírt szám sem az EN, sem a HU szövegben; 7 modulközi
  egyezés; a korlátfájlok és az ábrák bájtra azonosak a modulfájlokkal; csak projekt-DOI; licencláb minden oldalon;
  friss futtatási naplók).
- `python3 src/check_en_hu.py`: 5/5 ellenőrzés sikeres, kilépési kód 0 (374 értékmakró-használat, 28 figyelmeztető jelölő,
  31 címke, 942 számtoken a PDF-szövegben, EN és HU között azonosan).

## 3. Negatív tesztek (az ellenőrzésnek el kell buknia egy elrontott másolaton)
Nyolc szándékosan elrontott másolat; mindegyik a megfelelő ellenőrzés bukását okozta, kilépési kód 1:

| elrontott másolat | bukó ellenőrzés |
|---|---|
| kézzel beírt szám az EN szövegtörzsben | check_claims 4 |
| módosított érték a `generated/values_en.tex` fájlban | check_claims 2 |
| módosított érték az állítástáblában | check_claims 1 és 2 |
| a HU más értékmakrót használ | check_en_hu 1 |
| törölt HU figyelmeztető jelölő | check_en_hu 3 |
| megváltoztatott ábra | check_claims 7 |
| idegen DOI | check_claims 8 |
| módosított HU érték, újrafordított HU PDF | check_en_hu 5 |

## 4. A készlet ellenőrzése friss kicsomagolásból
- `sha256sum -c SHA256SUMS.txt`: minden összeg rendben; nincs a listán kívüli fájl (magát a `SHA256SUMS.txt` fájlt kivéve).
- A kicsomagolt készletben check_claims 14/14 és check_en_hu 5/5; mindkét PDF 10 oldalas.
- Újraépítés egy másolatban (`python3 src/build_article.py ..`): az ellenőrzések ismét sikeresek, és minden generált fájl,
  a `CLAIMS_SOURCES.csv` és a `CLAIMS_SOURCES.md` azonos a mellékelt példánnyal.
- A készlet három modularchívumát újra lefuttattuk, mindegyiket külön üres mappában: óra kilépési kód 0 / `RESULT: PASS`;
  BREN kilépési kód 0 / `all steps passed`; szub-cm kilépési kód 0 / `all steps passed (frozen-digitization mode)`.

## 5. Vizuális ellenőrzés
Mindkét PDF minden oldalát (10 EN + 10 HU) 150 dpi-vel rendereltük (`pdftoppm -r 150`) és átnéztük: nincs túlcsorduló sor,
a táblázatok és ábrák olvashatók, minden oldalon ott a licencet és oldalszámot tartalmazó lábléc.
