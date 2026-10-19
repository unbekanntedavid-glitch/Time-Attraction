# Time Attraction Law v17 – SPARC Audit

## Gyors indítás
1. Klónozd a repót.
2. Hozd létre a `data/` mappát.
3. Tedd a `Rotmod_LTG.zip`-et a `data/` mappába (ne csomagold ki).
4. `python3 sparc_full_runner.py`

A futtató ellenőrzi a ZIP hash-t (`0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`),
kicsomagolja, és lefut a hat modellen.

## Hibametrika
σ = max(0.05 · Vobs, 20 km/s)  
χ² = Σ [(Vpred − Vobs)/σ]²  
χ²_red = χ² / (N − 1)

## Licenc
CC BY-NC 4.0
