# Extracting V28 (two ZIP files)

V28 is distributed as two ZIP files:

- `TAL_V28_r3_part1_step26.zip` — `TAL_V28/step26_cosmology_desi_bao/` (+ `SHA256SUMS_part1.txt`)
- `TAL_V28_r3_part2_simulator_step27.zip` — `TAL_V28/simulator/`, `TAL_V28/step27_local_yukawa_branch/`, the READMEs, `DATA_LICENSES.md` and the main `SHA256SUMS.txt` (+ `SHA256SUMS_part2.txt`)

1. Extract **both** ZIP files into the **same** directory. Both contain the top folder `TAL_V28/`; the contents merge into one `TAL_V28/` folder (this file is in both ZIPs with identical content, so overwriting it is harmless).
2. Verify from inside `TAL_V28/`:

```
cd TAL_V28
sha256sum -c SHA256SUMS.txt          # valid only when BOTH parts are extracted together
sha256sum -c SHA256SUMS_part1.txt    # checks part 1 on its own
sha256sum -c SHA256SUMS_part2.txt    # checks part 2 on its own
```

The main `SHA256SUMS.txt` lists every file of both parts; with only one part extracted it reports the files of the other part as missing. Each part can be checked alone with its own `SHA256SUMS_partN.txt`.
