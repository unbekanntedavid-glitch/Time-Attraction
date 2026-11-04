# TAL publikációs készlet, 2026-11-04, r2

*Egy rögzített skalár–tenzor ág óra- és gyorsulás-előrejelzései: kalibráció és Yukawa-korlátok* című közlemény a három
futtatható modullal, amelyekből számolva van. Előzetes kutatás. Az ellenőrzések módja: numerikus keresztellenőrzés az
Academia Co-scientist segítségével (MI-alapú). Az ellenőrzések technikai jellegűek (szkriptelt reprodukció, numerikus
összevetés, negatív kontrollok); egy sikeres ellenőrzés nem empirikus megerősítés.

## Tartalom

| elem | mi ez |
|---|---|
| `article/` | a közlemény (angol és magyar PDF), LaTeX-források, állítás–forrás táblázat, ellenőrző szkriptek, friss futtatási naplók, változás- és tesztnapló (lásd `article/README.md`) |
| `tal_branch_clock_accel_test_r4.zip` | óra- és gyorsulásmodul (fő elemzés) |
| `bren_continuation_zenodo_r4.zip` | BREN-torony folytatási modul (korlátmodul, Zenodo-változat) |
| `subcm_yukawa_window_zenodo_r4.zip` | szub-centiméteres Yukawa-ablak modul (korlátmodul, Zenodo-változat) |
| `LICENSE.txt` | licenc, a harmadik fél anyagaira vonatkozó megjegyzéssel |
| `SHA256SUMS.txt` | a készlet minden fájljának SHA-256 összege |

## A modulok kapcsolata

Az óramodul a két korlátmodul kimeneteit használja: az `inputs/bounds/` mappája a BREN-eredményfájl (történeti nevén
`bren_continuation_results_r2.json`) és három szub-cm kimenet bájtazonos másolatát tartalmazza, SHA-256 összegükkel
az `inputs/PROVENANCE_bounds.json` fájlban. A három modul tehát egy lánc, nem független bizonyíték. Minden archívumot
külön mappába kell kicsomagolni; ne csomagolja ki őket egymásra.

## Fő parancsok

A parancsok az angol README-ben vannak (`README_EN.md`, „Main commands”): `sha256sum -c SHA256SUMS.txt`; a három archívum
kicsomagolása három külön mappába; a modulokban `python3 run_check.py`; az `article/` mappában
`python3 src/check_claims.py` és `python3 src/check_en_hu.py`. A modulellenőrzésekhez Python 3, NumPy, SciPy és matplotlib
kell (verziók a modulok `requirements.txt` fájljában; a naplózott futások verziói az `article/logs/` mappában); a
közlemény-ellenőrzésekhez csak Python 3 és poppler (`pdftotext`, `pdfinfo`). A szub-cm futás a befagyasztott digitalizálási
adatokat használja; a BREN DEM újragenerálása a nyers Copernicus-csempékből külön feladat, amelyet a modul nem végez el.

## Azonosítók

| modul | DOI |
|---|---|
| TAL-mű, koncepció-DOI | 10.5281/zenodo.22983681 |
| óra és gyorsulás | V20 10.5281/zenodo.23540741 (forrásgeometria), V26 10.5281/zenodo.23727585 (ág és hatás), V28 10.5281/zenodo.23814053 (kalibrációs konvenció, Step 27) |
| BREN-folytatás | V28 10.5281/zenodo.23814053 (a Step 27 finomítása) |
| szub-cm Yukawa-ablak | V28 10.5281/zenodo.23814053 (a Step 27 jelöltje) |

HU
© 2026 Szabolcsi Dávid (ORCID: 0009-0005-6300-074X)
Licenc: Creative Commons Nevezd meg! – Ne add el! 4.0 Nemzetközi (CC BY-NC 4.0)
https://creativecommons.org/licenses/by-nc/4.0/deed.hu
Valamennyi kereskedelmi, ipari és technológiai hasznosítási jogot kifejezetten a szerző tart fenn.
Kereskedelmi hasznosítás kizárólag a szerző írásos engedélyével. Harmadik fél anyagaira a saját feltételeik vonatkoznak (lásd LICENSE.txt).
