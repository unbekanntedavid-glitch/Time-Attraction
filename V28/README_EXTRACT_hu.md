# A V28 kibontása (két ZIP-fájl)

A V28 két ZIP-fájlban jelenik meg:

- `TAL_V28_r3_part1_step26.zip` — `TAL_V28/step26_cosmology_desi_bao/` (+ `SHA256SUMS_part1.txt`)
- `TAL_V28_r3_part2_simulator_step27.zip` — `TAL_V28/simulator/`, `TAL_V28/step27_local_yukawa_branch/`, a README-k, a `DATA_LICENSES.md` és a fő `SHA256SUMS.txt` (+ `SHA256SUMS_part2.txt`)

1. **Mindkét** ZIP-et **ugyanabba** a könyvtárba bontsd ki. Mindkettőben a `TAL_V28/` a legfelső mappa; a tartalom egyetlen közös `TAL_V28/` mappába olvad össze (ez a fájl mindkét ZIP-ben azonos tartalommal szerepel, felülírása nem okoz gondot).
2. Ellenőrzés a `TAL_V28/` mappán belülről:

```
cd TAL_V28
sha256sum -c SHA256SUMS.txt          # csak akkor érvényes, ha MINDKÉT rész együtt van kibontva
sha256sum -c SHA256SUMS_part1.txt    # az 1. rész önálló ellenőrzése
sha256sum -c SHA256SUMS_part2.txt    # a 2. rész önálló ellenőrzése
```

A fő `SHA256SUMS.txt` mindkét rész összes fájlját felsorolja; ha csak az egyik rész van kibontva, a másik rész fájljait hiányzónak jelzi. Az egyes részek a saját `SHA256SUMS_partN.txt` listájukkal külön is ellenőrizhetők.
