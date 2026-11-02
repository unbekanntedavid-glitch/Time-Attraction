# Time Attraction Law (TAL): Derivation Chain Audit and Step 29 Corrections S1–S3

This record contains the audit of the derivation chain of the Time Attraction Law (TAL), from the first paper (Time Attraction, v1–v7) to the current version V28, together with corrections S1–S3 to the Step 29 analyses and the full verification material.

Main result. The derivation first breaks in Time Attraction v1–v7, §3.1, Equations (9)–(11). Equation (9) correctly gives the Newtonian potential. Equation (10) relabels the ordinary general-relativistic redshift 1 + Φ/c² as a new factor f_L and drops the geometric factor 4π/3. Equation (11) inserts this factor into the action a second time, next to the Einstein–Hilbert term; in the static weak-field, point-particle limit the resulting dynamics count the potential twice (Φ_eff = Φ + c² ln f_L). This double counting produces the universal, gravitational-strength fifth force of the later versions (α = 1/(4π) ≈ 0.0796 at λ = L = 10 km), which laboratory and geophysical fifth-force searches exclude by a factor of about 177. The double counting persists for f_L = exp(Φ/c²) as long as f_L is retained as a separate factor in the matter action (Φ_eff = 2Φ); only if f_L is read as the clock rate already contained in g₀₀ does it coincide with GR, to first order in the weak field, and add no new physics.

Corrections. The audit also found three errors in the author's own Step 29 analyses:
S1 – a spurious long-range allowed window, now removed using the Cassini bound;
S2 – an incorrect "clock-only H1" wording, corrected throughout;
S3 – α_Y = 1/3 values that depended on one coupling convention, now given in both conventions.
None of the corrections changes any conclusion. The Galileo 5/6 reanalysis (report and reproducibility package) was rebuilt after the corrections.

Contents. academy_review_package.zip contains:
- the audit and corrections document in English and Hungarian (PDF and Markdown);
- the Galileo 5/6 H1 reanalysis report (PDF) and its reproducibility package (ZIP);
- the fifth-force analysis and its results;
- the dependency diagram and the machine-readable edit log (fix_log.json);
- the verification/ folder with the scripts, symbolic checks and results referenced in the appendix of the audit;
- SHA256SUMS covering all files.

Reproducibility limit: the verification files of the audit are included, but some scripts require absolute project paths and input data that are not part of this package, so the full analysis cannot be re-run independently from this ZIP alone.

SHA-256 of academy_review_package.zip: 6b0cd3b6ecfa4dd3e6461f475561858b4677b0e365632d4e755a522b8bb640b5

Related records: TAL V28 (10.5281/zenodo.23814053; concept DOI 10.5281/zenodo.22983681); Time Attraction v1–v7 (10.5281/zenodo.23066993).

© 2026 Dávid Szabolcsi. License: CC BY-NC 4.0.
