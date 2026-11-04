# Article / Közlemény

EN. *Clock and acceleration predictions of a fixed scalar–tensor branch: calibration and Yukawa constraints* (EN and HU PDF).

| path | content |
|---|---|
| `TAL_clock_accel_yukawa_article_EN.pdf`, `TAL_clock_accel_yukawa_article_HU.pdf` | the article |
| `tex/` | LaTeX sources (`article_en.tex`, `article_hu.tex`, shared `common.tex`) |
| `generated/` | value macros (`values_en.tex`, `values_hu.tex`) and table bodies, written by `src/make_values.py` |
| `figures/` | module figures, byte-for-byte copies of the module files (checked) |
| `CLAIMS_SOURCES.csv`, `CLAIMS_SOURCES.md` | claim–source table: every number of the article, its archive, file, key and formula, and the sections where it appears |
| `src/` | `make_values.py`, `build_article.py`, `check_claims.py`, `check_en_hu.py` and helpers |
| `logs/` | fresh-run logs of the three modules (see `logs/README.md`) |
| `CHANGELOG_EN.md`, `VALTOZASNAPLO_HU.md` | itemised change list |
| `TEST_LOG_EN.md`, `TESZTNAPLO_HU.md` | test log |

How the numbers get into the text: every number is a macro `\V{id}`. `src/make_values.py` reads each value directly
from the three module archives in the set folder (no unpacking), formats it (HU with decimal comma) and writes the macro
files, the tables and the claim table. The article sources contain no hand-typed values; `src/check_claims.py` enforces this.

Commands (from this folder; Python 3 with the standard library; pdfLaTeX only for rebuilding):

```
python3 src/check_claims.py        # claims against the archives in ..  (exit 0 = all checks pass)
python3 src/check_en_hu.py         # EN/HU agreement in numbers, units, signs and caveats
python3 src/build_article.py ..    # optional: regenerate values and tables and rebuild both PDFs
```

HU. *Egy rögzített skalár–tenzor ág óra- és gyorsulás-előrejelzései: kalibráció és Yukawa-korlátok* (angol és magyar PDF).
A mappák tartalma a fenti táblázatban van. A szövegben minden szám egy `\V{id}` makró; a `src/make_values.py` minden
értéket közvetlenül a készletmappa három modularchívumából olvas ki (kicsomagolás nélkül), formáz (magyarul tizedesvesszővel),
és megírja a makrófájlokat, a táblázatokat és az állítás–forrás táblázatot. A forrásokban nincs kézzel beírt érték; ezt a
`src/check_claims.py` ellenőrzi. Parancsok: lásd fent (`check_claims.py`, `check_en_hu.py`, opcionálisan `build_article.py ..`).

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International · https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission / kereskedelmi felhasználás kizárólag a szerző írásos engedélyével
