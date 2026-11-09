# Equilibrium binding energy at fixed baryon and ADM mass in a local massive-scalar neutron-star model — revised contribution statement and literature context (revision b)

YBE-2026-11-08-b. Publication date: 2026-11-08. Author: Szabolcsi Dávid (ORCID [0009-0005-6300-074X](https://orcid.org/0009-0005-6300-074X)).

Zenodo DOI: [10.5281/zenodo.24054760](https://doi.org/10.5281/zenodo.24054760)  
Concept DOI: [10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681) · Zenodo version index 37.

This DOI is revision b. Revision a remains [10.5281/zenodo.24053942](https://doi.org/10.5281/zenodo.24053942), folder [`YukawaBindingEnergy_2026-11-08/`](../YukawaBindingEnergy_2026-11-08/), GitHub tag V37, Zenodo version index 36. Zenodo creator name: Szabolcsi, David.

The deposited archives still say that no DOI was assigned to b. `metadata/provenance.json` has `publication_identifier` null, and `DELIVERY_MANIFEST.json` has `doi_assigned` false. Those bytes were written before this DOI and stay unchanged. The DOI is recorded in this repository note.

Numerical cross-checks with Academia Co-scientist (AI-assisted); not peer review.

## Scope

Revision b narrows the contribution statement and the literature context. The contribution is a reproducible equilibrium dataset and a documented comparison protocol for this specified phenomenological local massive-scalar fluid model. It is not a new general scalar–fluid coupling, a new binding-energy method, a physical validation, a priority claim, a stability result, or an observability result. Literature priority is not established.

SLy, APR4 and H4, α_Y = 1, L = 1 km and 10 km. E_bind = (M_B − M_ADM) c². Fixed baryon mass compares equilibrium energies. Fixed ADM mass compares baryon capacity. These are two views of one equilibrium-energy result. Near-limit comparisons stay on the jointly accessible causal portions of the selected branches. At the baryon mass of each EOS’s GR configuration with ADM mass 1.4 solar masses, the 10 km Yukawa−GR differences remain 10.850, 11.233 and 8.443 × 10^51 erg for SLy, APR4 and H4. Those states are not a universal 1.5-solar-mass comparison.

The 16 required scientific inputs listed in `metadata/provenance.json` have the same SHA-256 values as revision a. This import did not change solver, EOS, tolerances, or those tables. Manuscript prose, bibliography, and some figure files differ from a. The results stay model-dependent equilibrium differences. They are not dynamical stability, observational detection, or emitted conversion energy. Numerical reproducibility is not independent physical validation.

The Zenodo record reports that the scientific data, solver, equations of state, and numerical tolerances are unchanged. The author-supplied Academia report of the 1,391 configurations, 114 controls, and 18 peak fits is the same attributed report as in the deposited text. This folder does not add a new simulation or a separate review.

## Files

| file | role |
|---|---|
| `YBE-2026-11-08-b_EN.pdf` | Deposited English article, 7 pages. Same bytes as `manuscript/binding_energy_EN.pdf` in the combined archive. |
| `YBE-2026-11-08-b_HU.pdf` | Deposited Hungarian article, 7 pages. Same bytes as `manuscript/binding_energy_HU.pdf` in the combined archive. |
| `YBE-2026-11-08-b_data_code.zip` | Deposited data and code archive. |
| `YBE-2026-11-08-b_editable_sources.zip` | Deposited editable-source archive. |
| `YBE_literature_comparison_20261108_r1.zip` | Separate literature supplement. It is not a numerical input to `reproduce.py`. |
| `DELIVERY_MANIFEST.json` | Deposited delivery manifest. `doi_assigned` remains false. |
| `DELIVERY_SHA256.txt` | Deposited checksum list for the delivery names. |
| `distribution_validation.json` | Deposited validation record from the delivery checks. Its recorded manuscript-PDF hashes come from a rebuild during that delivery and differ from the deposited article PDFs. |
| `YBE-2026-11-08-b_publication_set.zip` | Combined archive. Its SHA-256 is a line in the deposited `DELIVERY_SHA256.txt`. Zenodo did not store this zip as its own file key. |
| `LICENSE.txt` | License notice. Same bytes as `LICENSE.txt` in the archives. |
| `SHA256SUMS` | SHA-256 manifest of the imported files in this folder. |

Zenodo stored two keys with a ` (1)` suffix: `distribution_validation (1).json` and `YBE_literature_comparison_20261108_r1 (1).zip`. Their bytes match `distribution_validation.json` and `YBE_literature_comparison_20261108_r1.zip` in the deposited checksum list, which are the names used here.

Shared scientific members of the data archive, the editable-source archive, and the combined archive match. Each archive’s own `SHA256SUMS.txt` differs because it lists that archive.

Verification: `cd YukawaBindingEnergy_2026-11-08_b && sha256sum -c SHA256SUMS`.

## Safe extraction and reproduction

Reject archive members whose names are absolute or contain `..` before writing them to disk. The deposited archives have no such members. After a safe extract of `YBE-2026-11-08-b_publication_set.zip`:

```sh
sha256sum -c SHA256SUMS.txt
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 reproduce.py --check --out run/check
```

`--check` reads the supplied inputs and writes only under `--out`. The package records Python 3.13.5 and the versions in `requirements.txt`. The literature zip is not an input to that command. `--controls`, `--regenerate`, `--archived-peak-validator`, figure regeneration, and `build_articles.sh` write further outputs and were not used for this import.

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International + https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission.

## Magyar

Egyensúlyi kötési energia rögzített barion- és ADM-tömegnél lokális tömeges skalármezős neutroncsillagmodellben — pontosított hozzájárulás és irodalmi kontextus (b revízió). Belső verzió: YBE-2026-11-08-b. Dátum: 2026-11-08.

Zenodo DOI: [10.5281/zenodo.24054760](https://doi.org/10.5281/zenodo.24054760). Fogalom-DOI: [10.5281/zenodo.22983681](https://doi.org/10.5281/zenodo.22983681). Zenodo verzióindex: 37.

A [10.5281/zenodo.24053942](https://doi.org/10.5281/zenodo.24053942) azonosító az a revízióé. A csomagok belsejében a b revízióhoz még nincs DOI beírva; ezek a bájtok változatlanok.

A hozzájárulás ennek a megadott modellnek a reprodukálható egyensúlyi adatállománya és dokumentált összevetési eljárása. Nem új általános csatolás, nem új kötésienergia-módszer, nem fizikai validáció, nem elsőbbség, nem stabilitás és nem észlelhetőség. Az irodalmi összevetés külön melléklet, nem numerikus bemenet. SLy, APR4 és H4, α_Y = 1, L = 1 km és 10 km. A közös bariontömeg melletti energia és a közös ADM-tömeg melletti barionkapacitás ugyanazon egyensúlyi eredmény két nézete.

Numerikus keresztellenőrzés az Academia Co-scientist segítségével (MI-alapú); nem szakmai lektorálás.

Ellenőrzés: `cd YukawaBindingEnergy_2026-11-08_b && sha256sum -c SHA256SUMS`.

© 2026 Szabolcsi Dávid (ORCID 0009-0005-6300-074X) · CC BY-NC 4.0 International + https://creativecommons.org/licenses/by-nc/4.0/ · commercial use only with the author's written permission.
