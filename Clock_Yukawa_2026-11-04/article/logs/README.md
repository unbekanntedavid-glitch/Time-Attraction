# Fresh-run logs / Friss futtatási naplók

EN. One log per module. Each module archive of this set was unpacked into its own empty folder (three separate
folders, not on top of each other) and its check was run once. The local folder path is replaced by `<unzip folder>`.
Each log starts with the environment (date, Python, platform, NumPy, SciPy, matplotlib, Pillow, CPU count),
then the command, the full output and the exit status.

| log | archive (SHA-256) | command | result |
|---|---|---|---|
| `fresh_run_clock.txt` | `tal_branch_clock_accel_test_r4.zip` (6ba806210d557189fc77777a29d3e9964b0b75f81fd185230a62eedfc0daa105) | `python3 run_check.py --out runs/fresh_article` | `RESULT: PASS`, exit 0 |
| `fresh_run_bren.txt` | `bren_continuation_zenodo_r4.zip` (5f217dc47b6f91756cda9a5b1765930be82e32c5f77a80456825cfe5168edf25) | `python3 run_check.py` | `all steps passed`, exit 0 |
| `fresh_run_subcm.txt` | `subcm_yukawa_window_zenodo_r4.zip` (e3f027587b7462900d607d858315668bc65346eada7489f34f3a06d31096ae46) | `python3 run_check.py` | `all steps passed (frozen-digitization mode)`, exit 0 |

The environment block records what was used for these runs; `requirements.txt` in each module lists versions and is not
by itself a proof of an environment. Figures produced by a run are not expected to be byte-identical; numbers are
compared with the reference results within the tolerance of each module.

HU. Modulonként egy napló. A készlet minden modularchívumát külön üres mappába csomagoltuk ki (három külön mappa,
nem egymásra), és egyszer lefuttattuk az ellenőrzését. A helyi mappaútvonal helyén `<unzip folder>` áll. Minden napló a
környezettel kezdődik (dátum, Python, platform, NumPy, SciPy, matplotlib, Pillow, CPU-szám), ezt követi a parancs, a
kimenet és a kilépési kód. A parancsok és eredmények a fenti táblázatban vannak. A környezeti blokk azt rögzíti, amivel a
futások készültek; a modulok `requirements.txt` fájlja verziókat sorol fel, és önmagában nem bizonyítja a környezetet.
A futás által készített ábrák bájtazonossága nem elvárás; a számokat a modulok saját toleranciájával vetjük össze a
referencia-eredményekkel.
