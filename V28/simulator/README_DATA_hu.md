# A szimulátor adatgenerálása (`src/build_data.py`)

Az `index.html` oldal úgy, ahogy van, működik: minden adat be van ágyazva (`src/data.js`, amelyet a `src/build.py` illeszt be). A `src/data.js` újragenerálása az alapbemenetekből **külön lépés, amelyhez külső V26-bemenetek kellenek**; ezek nem részei a V28-csomagnak.

## Bemenetek

| Bemenet | Hely | A V28-ban? |
|---|---|---|
| `tal_cosmology_common.py`, `run_step1_2_and_1_3.py`, `run_step1_4_to_1_6.py`, `tal_hubble_desi_bao_summary.json` (és a `tal_cosmology_common.py` által olvasott DESI átlag-/kovariancia-fájlok) | V27 26. lépéses fa; a V28-ban: `../step26_cosmology_desi_bao/` (alapértelmezés). Továbbá: V27, doi:[10.5281/zenodo.23761574](https://doi.org/10.5281/zenodo.23761574); GitHub `V27/`: https://github.com/unbekanntedavid-glitch/Time-Attraction/tree/main/V27 | igen |
| `gw170817_and_lensing_summary.json`, `lensing_deflection_and_ppn_slip.csv`, `gw170817_dispersion_and_cT.csv` (a V26 23. lépésének kimenetei) | a V26-fa gyökere: V26, doi:[10.5281/zenodo.23727585](https://doi.org/10.5281/zenodo.23727585) (a V26 ZIP `scalar-tensor-yukawa-fit/` mappája), vagy GitHub `V26/`: https://github.com/unbekanntedavid-glitch/Time-Attraction/tree/main/V26 | **nem (külső)** |
| SPARC forgási görbék, `*_rotmod.dat` (175 fájl) | `<V26>/data/Rotmod_LTG/` (ugyanazok a V26-források). Eredeti adat: Lelli, McGaugh & Schombert 2016, AJ 152, 157, doi:[10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157); https://astroweb.case.edu/SPARC/ | **nem (külső)** |

## Parancs

A kibontott `TAL_V28/` mappán belülről:

```
python3 simulator/src/build_data.py --v26 /utvonal/V26 --out /tmp/data.js
cmp /tmp/data.js simulator/src/data.js     # bájtra azonos eredmény várható
python3 simulator/src/build.py              # a src/* fájlokat a simulator/index.html-be illeszti
node simulator/verify_node.js               # 117 újraszámolt mennyiség a közzétett JSON-nal összevetve
```

Kapcsolók (vagy környezeti változók): `--v27` / `TAL_V27_DIR` (alapértelmezés: `../../step26_cosmology_desi_bao` a `src/`-hez képest), `--v26` / `TAL_V26_DIR` (kötelező), `--sparc` / `TAL_SPARC_DIR` (alapértelmezés: `<V26>/data/Rotmod_LTG`), `--out` (alapértelmezés: `src/data.js`). Hiányzó bemenet esetén a szkript a hiányzó fájlok listájával leáll. A bemeneteket csak olvassa, a bemeneti fákba nem ír.

Próba (2026-10-30): egy kibontott `TAL_V28/` fából, a `--v26` kapcsolóval a GitHub `V26/` mappájára mutatva, az újragenerált `data.js` (150 822 B, 175 galaxis) bájtra azonos a mellékelt `src/data.js` fájllal.

Megjegyzés: a `src/data.js` (és így az `index.html`) számszerű bemenetként tartalmazza a SPARC forgásigörbe-táblázatokat és a DESI DR1 BAO átlag-/kovarianciaértékeit; eredeti forrásaik fent és a 26. lépés jelentéseiben szerepelnek.

Adatforrások és licencek (DESI CC BY 4.0, SPARC: nincs kifejezett licenc, RSD, Planck, BBN), szó szerinti hivatkozási és köszönetnyilvánítási szövegek: `../DATA_LICENSES.md`. A szimulátor lábléce ugyanezt röviden, HU/EN nyelven jelzi.
