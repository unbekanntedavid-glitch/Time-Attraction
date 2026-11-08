# Equilibrium binding energy at fixed baryon and ADM mass in a local massive-scalar neutron-star model

YBE-2026-11-08-a. Publication date: 2026-11-08. Author: Szabolcsi Dávid (ORCID [0009-0005-6300-074X](https://orcid.org/0009-0005-6300-074X)).

Zenodo DOI: [10.5281/zenodo.24053942](https://doi.org/10.5281/zenodo.24053942)  
Concept DOI: [10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681) · Zenodo version index 36.

Zenodo creator name: Szabolcsi, David. This record is a distinct study from the time-field mass record [10.5281/zenodo.24022021](https://doi.org/10.5281/zenodo.24022021). The public API also returns `isNewVersionOf` [10.5281/zenodo.23699685](https://doi.org/10.5281/zenodo.23699685), the same concept-series relation as the neighbouring versions.

The deposited manuscripts, data, code and sources are the published bytes. `metadata/provenance.json` inside the archives still has `publication_identifier` null: that file was written before this DOI and stays unchanged. The English and Hungarian manuscripts keep the published Arnowitt–Deser–Misner bibliographic page range (pp. 227–264 in the deposited TeX). A minor page-range discrepancy remains, as stated on the Zenodo record; the manuscript bytes are unchanged, and the discrepancy does not affect the numerical results.

Numerical cross-checks with Academia Co-scientist (AI-assisted); not peer review.

## Scope

Phenomenological local massive-scalar neutron-star equilibria with the SLy, APR4 and H4 equations of state, α_Y = 1, and ranges L = 1 km and 10 km. Binding energy is E_bind = (M_B − M_ADM) c². Baryon mass uses proper spatial volume; ADM mass includes the exterior scalar-field energy.

Fixed baryon mass and fixed ADM mass are two views of one equilibrium-energy result. Fixed baryon mass compares energies. Fixed ADM mass compares baryon capacity. Near-limit comparisons use the causally admissible first-rising branches. At the baryon mass of each EOS’s GR configuration with ADM mass 1.4 solar masses, the 10 km Yukawa−GR binding-energy differences are 10.850, 11.233 and 8.443 × 10^51 erg for SLy, APR4 and H4.

These are model-dependent equilibrium diagnostics. They are not a dynamical-stability result, an observational detection or emitted conversion energy. Numerical reproducibility is not independent physical validation. The Zenodo record reports subsequent AI-assisted reproduction of the 1,391 equilibrium configurations, 114 convergence-control cases and 18 local peak fits. This folder keeps that published report and does not add a separate review.

## Files

| file | role |
|---|---|
| `binding_energy_EN.pdf` | English article, 6 pages. Same bytes as `manuscript/binding_energy_EN.pdf` in the publication set. |
| `binding_energy_HU.pdf` | Hungarian article, 6 pages. Same bytes as `manuscript/binding_energy_HU.pdf` in the publication set. |
| `YBE-2026-11-08-a_publication_set.zip` | Published combined package: both articles, editable sources, data, code, figures, license, provenance and its own `SHA256SUMS.txt`. SHA-256 `a22e22882ff46dc6efa9dc782fed81e0162685539444d76930ee757d2fa7af8b`, 1755699 bytes. |
| `YBE-2026-11-08-a_data_code.zip` | Published data and code archive. Shared scientific members match the publication set. |
| `YBE-2026-11-08-a_editable_sources.zip` | Published editable-source archive. Shared scientific members match the publication set. |
| `LICENSE.txt` | License notice. Same bytes as `LICENSE.txt` in each archive. |
| `SHA256SUMS` | SHA-256 manifest of the imported files in this folder. |

The three archives are the published file set of this one study. Every shared member other than `SHA256SUMS.txt` matches across the archives that contain it. Each archive’s `SHA256SUMS.txt` differs because it lists that archive.

Verification: `cd YukawaBindingEnergy_2026-11-08 && sha256sum -c SHA256SUMS`.

## Safe extraction and reproduction

Reject archive members whose names are absolute or contain `..` before writing them to disk. The published archives have no such members. After a safe extract of `YBE-2026-11-08-a_publication_set.zip`:

```sh
sha256sum -c SHA256SUMS.txt
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 reproduce.py --check --out run/check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 reproduce.py --controls --out run/controls
```

`--check` and `--controls` read the supplied inputs and write only under the chosen `--out` directory. The package records Python 3.13.5 and the versions in `requirements.txt`. Article rebuild, figure regeneration, `reproduce.py --regenerate` and `--archived-peak-validator` write further outputs; they were not used for this import. Input provenance remains in `metadata/provenance.json`, `metadata/source_members.csv` and `metadata/number_provenance.md`.

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International + https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission.

## Magyar

Egyensúlyi kötési energia rögzített barion- és ADM-tömegnél lokális tömeges skalármezős neutroncsillagmodellben. Belső verzió: YBE-2026-11-08-a. Dátum: 2026-11-08. Szerző: Szabolcsi Dávid (ORCID [0009-0005-6300-074X](https://orcid.org/0009-0005-6300-074X)).

Zenodo DOI: [10.5281/zenodo.24053942](https://doi.org/10.5281/zenodo.24053942). Fogalom-DOI: [10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681). Zenodo verzióindex: 36.

Ez a rekord külön tanulmány a [10.5281/zenodo.24022021](https://doi.org/10.5281/zenodo.24022021) időmező-tömeg rekordtól. A kéziratok, az adatok és a forrás a közzétett bájtok. A csomagbeli `metadata/provenance.json` továbbra is `publication_identifier: null`, mert a DOI előtt készült, és így marad. Az ADM-hivatkozás oldalszáma a közzétett kéziratban marad; a számtani eredményeket nem érinti.

SLy, APR4 és H4, α_Y = 1, L = 1 km és 10 km. A közös bariontömeg melletti energiarangsor és a közös ADM-tömeg melletti barionkapacitás ugyanazon egyensúlyi eredmény két nézete, kauzális első emelkedő ágakon. Az eredmény modellfüggő egyensúlyi diagnosztika: nem dinamikai stabilitás, nem észlelés és nem kibocsátott átalakulási energia. A numerikus reprodukció nem független fizikai validáció.

Numerikus keresztellenőrzés az Academia Co-scientist segítségével (MI-alapú); nem szakmai lektorálás.

Ellenőrzés: `cd YukawaBindingEnergy_2026-11-08 && sha256sum -c SHA256SUMS`. Kicsomagolás előtt az abszolút vagy `..` tagú útvonalakat el kell utasítani.

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International + https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission.
