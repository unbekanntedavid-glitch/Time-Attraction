## GitHub V18 csomag (= Zenodo V20, DOI 10.5281/zenodo.23540741)

**Ez a GitHub V18-as csomag.** GitHub V18 = Zenodo V20: „Time Attraction Law V20 — audited release with independent
third-round review” (2026-10-21), DOI: https://doi.org/10.5281/zenodo.23540741 ·
Licenc: CC BY-NC 4.0 · Szerző: David Szabolcsi (ORCID 0009-0005-6300-074X)

Előző verzió: V19 – 7th model, DOI 10.5281/zenodo.23497128.

### Tartalom (`V18/`)

- [`V18/zenodo_v20/`](V18/zenodo_v20/) — a Zenodo V20 rekord négy fájlja, változatlanul (MD5 a Zenodo API alapján ellenőrizve):
  `final_report (3).zip`, `report_hu (3).pdf`, `revision3_review_HU.md`, `revision3_review_EN.md`.
- [`V18/final/`](V18/final/) — **a V18 végleges jelentéscsomagja**: `final_report.zip` és kicsomagolt tartalma
  (`report_hu.pdf`, `report_hu.md`, `code/final_fixed.py`, `code/energy_integrator.py`, `code/fields.py`,
  `code/tophat_finite_source.py`, `code/build_pdf.py`, `csv/*.csv`, `csv/fixed_case_numbers.json`,
  `plots/fixed_signal_vs_L.png`, `refs/oscilloquartz_8663_ed03-2003.pdf`, `SHA256SUMS.txt`).
  A Zenodo V20 jelentéshez képest a revision3 review két szövegjavítását tartalmazza szó szerint
  (távoli „0 szórás” megjegyzés, hőkontroll átfogalmazása), továbbá a review szerinti pontosításokat
  (L = 0.95 m illesztési hiba pontos értéke, az „izolált, polaritásra páros tekercs-B² tag” megnevezés az 5. fejezetben, Ítélet).
  A kód és a CSV/JSON fájlok bájtra azonosak a Zenodo V20 csomagéval.
- [`V18/revision3_checks.json`](V18/revision3_checks.json) — számított ellenőrzések a `V18/final` csomagon
  (a szerző pipeline-ja generálta 2026-10-21-én a [`V18/make_revision3_checks.py`](V18/make_revision3_checks.py) szkripttel;
  **nem** a bíráló eredeti fájlja, mert a bíráló `check_revision3.py` szkriptje nem áll rendelkezésre).
- A menetátmérő **mért** érték: 291.7 mm (tolómérő, 0.02 mm felbontás, 5 egymást követő leolvasás, mind 291.7 mm;
  2025-11-13, mérte: Szabolcsi Dávid) — lásd a jelentés 1.2 pontját (mérési jegyzőkönyv). A 300 mm a tekercstest pereme (külső méret).
- NIST SP 1065 PDF nem része a csomagnak (sem a Zenodo V20 ZIP-ben, sem a `V18/final`-ban nem szerepel).
- Változások a V19 óta: [`CHANGELOG.md`](CHANGELOG.md).

### Mellékletek és SHA-256

| Fájl | SHA-256 |
|---|---|
| [`V18/zenodo_v20/revision3_review_HU.md`](V18/zenodo_v20/revision3_review_HU.md) | `e2a8f32efaddcf188a7926ef7cd0297db391335bd043beb0b7d1b4396deae8ec` |
| [`V18/zenodo_v20/revision3_review_EN.md`](V18/zenodo_v20/revision3_review_EN.md) | `103e77e72906312f16090822712961494ddf8b502b835d3121bb4277a1598d30` |
| [`V18/revision3_checks.json`](V18/revision3_checks.json) | `836e75b2baa67c9084891f3bb5dba33ae218f6f75025eaa00add8b366114e01f` |
| `V18/make_revision3_checks.py` | `7189e747d1424806e1ee233c3f10164e329df7b5c6a066535a91f85bf7ad0a52` |
| `V18/final/final_report.zip` | `5540ebdb25137fe48c2e0aa05c9e7057d2407c22a9891ca5daa37505b7ad2239` |
| `V18/final/report_hu.pdf` | `c4225b1917b493db03e8733228dd5337a010fa5bbb3d7894e8bd5a7b5fedae81` |
| `V18/zenodo_v20/final_report (3).zip` | `160981d3adb435b4fa719ca1a9d9a04d4e889a27f9d4c7bb0381e5c057ad2ea8` |
| `V18/zenodo_v20/report_hu (3).pdf` | `544567958646b397e77036ade822da75bdf2f79b4dd18235bcff6d5ca14ab20a` |

Ellenőrzés: `cd V18 && sha256sum -c SHA256SUMS_V18.txt`, `cd V18/zenodo_v20 && sha256sum -c ../SHA256SUMS_ZENODO_V20.txt`,
`cd V18/final && sha256sum -c SHA256SUMS.txt`.

---

# Time-Attraction Law v17.0 (2026-10-19) - SPARC Baryonic Baselines Locked Blind Audit

DOI v17: 10.5281/zenodo.23479940 | Parent: 10.5281/zenodo.22983681
DOI V19 – 7th model (independent verification): 10.5281/zenodo.23497128
ORCID: 0009-0005-6300-074X | License: CC BY-NC 4.0
GitHub: https://unbekanntedavid-glitch.github.io/Time-Attraction/

## v17.0 - Locked blind audit, 123 galaxies, 2596 radii

- Blind input: blind_123.txt (hash-verified)
- Rotmod data: Rotmod_data.zip (hash-verified)
- Six model outputs: vbar, vbar_095, vbar_105, no_bulge, no_gas, disk_only
- Runner: sparse_full_runner.py — 123 galaxies, 2596 radius points
- Interactive page: index.html (SPARC Baryonic Baselines)
- simulation.html redirects to index.html (v8.2 Phi_t Protocol retained as archive)
- Full PDF: Time-Attraction_v17 publication zip on Zenodo
- Changelog mapping v1–v17 fixes included

## Strategy

v17 deliberately does NOT claim absence of dark matter. It establishes the baryonic baseline as testable and reproducible. Dark matter implications follow only if independent replication confirms the baseline.

## Previous versions

v16: 10.5281/zenodo.23234825 | v15: 10.5281/zenodo.23220318 | v14: 10.5281/zenodo.23217189

> "The time does not just curve. The time is attracted."
