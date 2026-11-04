# Article change log (2026-11-04)

New article, first version in this set. The computational modules were not changed: the three r4 archives are included
byte-for-byte (SHA-256 in `../SHA256SUMS.txt`).

1. Chapter plan fixed before writing: 1 introduction and categories; 2 action, matter coupling, notation; 3 clock and
   acceleration from one source (one potential, laboratory calibration, observables, declared geometry); 4 numerical
   method, checks, reproducibility; 5 results for P1–P3 (with feasibility indicator S); 6 Yukawa bounds (BREN chain,
   sub-cm chain, ratios); 7 conclusion; appendices A BREN operator, B digitization and bound constructions,
   C non-universal addon; data and code availability; references.
2. Claim–source system: `src/claims_spec.py` names, for every number, the archive, file, key or row and the formula;
   `src/make_values.py` writes the EN/HU value macros, the tables and `CLAIMS_SOURCES.csv/.md` (295 claims).
3. Text written with value macros only; tables are generated from the same macros.
4. Figures embedded as PNG copies of the module figures (no text layer, so the EN/HU PDF number comparison is not
   disturbed by figure labels); the clock figure has English labels in both versions, stated in the HU caption.
5. During drafting: hand-typed numbers found by the digit scan were replaced by value macros (the GR-relative
   size of δy, the truncation radii of the convergence configurations, the 1 cm limit of the quadrature grid, now
   written as 10^quadLo m); the sentence “R = 10 does not mean 10σ” was reworded without numbers.
6. Unused entries removed from the claim table (P2/P3 raw columns, two operator self-test entries with no unit in the
   module output, the three PASS flags; the PASS results are now checked directly in the logs by `check_claims.py`, check 10).
7. Points per decade: the quadrature grid and the bound-only grid are stated with their own values.
8. DOIs set in an unbreakable box, so a line break inside a DOI cannot change the PDF number comparison.
9. References limited to sources named in the module packages; DOIs limited to the four TAL records.
10. HU wording of the checking method placed as a separate sentence with the fixed phrase.
