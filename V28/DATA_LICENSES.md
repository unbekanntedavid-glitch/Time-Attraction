# Adatforrások és licencek / Data sources and licenses — TAL V28

Ellenőrizve / Checked: 2026-10-30 (Europe/Zurich). Every URL below was opened on that date (HTTP 200 unless noted). / Minden alábbi URL-t ezen a napon nyitottunk meg (HTTP 200, ha nincs más jelölve).

**HU.** A csomag saját tartalmára (kód, szöveg, ábrák, szimulátor) a CC BY-NC 4.0 International licenc vonatkozik. A beágyazott külső adatok megtartják a saját feltételeiket; a CC BY-NC 4.0 ezekre nem terjed ki és nem szűkíti őket. Ahol a forrás nem közöl kifejezett licencet, ezt alább kimondjuk; ilyenkor csak a forrás által kért hivatkozást adjuk meg, és nem állítunk újraterjesztési engedélyt.

**EN.** The CC BY-NC 4.0 International license covers the author's own content of this package (code, text, figures, simulator). Embedded third-party data keep their own terms; CC BY-NC 4.0 does not extend to them and does not restrict them. Where a source states no explicit license, this is said plainly below; in those cases we only give the citation the source asks for and claim no redistribution permission.

## Összefoglaló táblázat / Summary table

| # | Forrás / Source | Hol van a V28-ban / Where in V28 | Licenc vagy feltételek / License or terms | Kötelező hivatkozás / Required citation | Származtatott számértékek újraterjesztése / Redistribution of derived numeric values |
|---|---|---|---|---|---|
| 1 | DESI DR1 BAO mean + covariance (12-vector, Cobaya format) and DESI 2024 Table values | `step26_cosmology_desi_bao/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`, `..._cov.txt` (byte-identical to the CobayaSampler/bao_data files); `tal_cosmology_common.py` (`DESI_BINS`); `simulator/src/data.js` → `index.html` (`desi`) | **CC BY 4.0** (DESI Data License page) | DR1 paper + acknowledgements text (verbatim below); indicate changes | **Allowed** under CC BY 4.0 with attribution, citation, change notice and acknowledgement text |
| 2 | SPARC rotation curves, 175 `*_rotmod.dat` (Lelli, McGaugh & Schombert 2016) | `simulator/src/data.js` → `index.html` (`sparc`); the raw files are not shipped | **No explicit license stated** on the official SPARC page; citation request only | "please cite the SPARC master paper" (verbatim below) | **Not stated** by the official page (see note on the 2025 Zenodo copy) |
| 3 | RSD fσ8(z), 12 published values (6dFGS, SDSS MGS, BOSS DR12, WiggleZ, eBOSS DR16, FastSound) | `run_step1_4_to_1_6.py` (`RSD_DATA`), `step1_4_rsd_fsigma8_growth.csv`; `data.js` (`rsd`) | **No data license stated**; values quoted from the journal papers | Cite the 12 source papers (list below) | **Not stated**; individual numbers quoted with citation |
| 4 | Planck 2018 parameters and priors (θ*, r_d, r_s, z*, z_drag, H0, ω_b, ω_cb, σ8) | `tal_cosmology_common.py`, `run_step1_4_to_1_6.py`, Step 26 CSV/JSON; simulator core | **No explicit data license stated** on the ESA "Planck Data Use" page; acknowledgement text requested | Planck acknowledgement (verbatim below) + Planck 2018 VI | **Not stated** |
| 5 | BBN reference values: Y_p = 0.245 ± 0.003 (PDG 2024 BBN review); D/H = (2.527 ± 0.030)×10⁻⁵ (Cooke, Pettini & Steidel 2018) | `run_step1_4_to_1_6.py`, `step1_5_bbn_consistency.csv`; simulator probes panel | PDG 2024 (Phys. Rev. D 110, 030001): article license **CC BY 4.0** (Crossref metadata). Cooke et al. 2018 (ApJ 855, 102): IOP copyright page, **no data license stated** | Cite both papers | PDG: allowed under CC BY 4.0 with attribution. Cooke 2018: **not stated** |
| 6 | "Pantheon+-like" SN grid (18 nodes) | `run_step1_2_and_1_3.py` (`PANTHEON_GRID`), `step1_3_*`; `data.js` (`sne_grid`) | **Not third-party data.** A synthetic shape-control grid of our own; it contains no Pantheon+ or Union3 catalogue values (the CSV name `..._pantheon_union3_...` is a legacy label) | — | — |
| 7 | Adelberger et al. 2003, Fig. 4 bound, digitized by us | `step27_local_yukawa_branch/.../reference_inputs/adelberger2003_fig4_digitized.json` | Our own reading of a published figure; the figure itself is **not redistributed** (rights not verified) | Adelberger, Heckel & Nelson 2003 (see sidecar README) | Only our digitized numbers are shipped |
| 8 | Ocean-gravity inputs (Zumberge et al. 1991 via March et al. 2019), Hoskins/Long/LLNL values | Step 27 code and documents | Numbers quoted from papers; **no data license stated**; secondary-source data via March 2019 | See `step27_local_yukawa_branch/sources/README_sources.md` | Individual numbers quoted with citation; source PDFs not redistributed |

## 1. DESI DR1 BAO

- Licenc oldal / License page: https://data.desi.lbl.gov/doc/acknowledgments/ ("Data License and Acknowledgments – DESI Data").
- License text (verbatim): "The Dark Energy Spectroscopic Instrument (DESI) data are licensed under the Creative Commons Attribution 4.0 International License (“CC BY 4.0”, Summary, Full Legal Code). Users are free to share, copy, redistribute, adapt, transform, and build upon the DESI data available through this website for any purpose, including commercially, as long as they 1. Cite the corresponding data release paper. […] 2. Indicate if any changes are made (if re-distributing DESI data). 3. Include the following acknowledgments text in any publications or derived works."
- CC BY 4.0: https://creativecommons.org/licenses/by/4.0/
- Where the likelihood files come from: the official DESI page "DESI DR1 BAO cosmology results" (https://data.desi.lbl.gov/doc/releases/dr1/vac/bao-cosmo-params/) states that "The DESI DR1 BAO likelihoods used to produce these results are publicly available at https://github.com/CobayaSampler/bao_data (in the files named `desi_2024_*`)". Our two files `desi_2024_gaussian_bao_ALL_GCcomb_mean.txt` / `_cov.txt` have the same git blob hashes as the files in that repository (6aff4520…, 7e3a6666…; committed 2024-04-05 "desi 2024 data"). The CobayaSampler/bao_data repository itself has no LICENSE file (GitHub API: license = none); the terms that apply are the DESI terms above.
- Data release paper (required): DESI Collaboration et al., "Data Release 1 of the Dark Energy Spectroscopic Instrument", arXiv:2503.14745, https://arxiv.org/abs/2503.14745.
- BAO papers (cited by the DESI page and by the Cobaya repository): DESI 2024 III, arXiv:2404.03000, doi:10.1088/1475-7516/2025/04/012; DESI 2024 IV, arXiv:2404.03001, doi:10.1088/1475-7516/2025/01/124; DESI 2024 VI, arXiv:2404.03002, doi:10.1088/1475-7516/2025/02/021.
- Changes (CC BY 4.0 notice): the 12-vector and covariance are used unchanged; in the simulator they are transcribed to JSON (`data.js`). `DESI_BINS` holds values from DESI 2024 III Table 18 and DESI 2024 VI Table 1; the entry `4c_ELG1_845` is a chosen illustration, **not** a DESI measurement (see the comment in `tal_cosmology_common.py`).
- Required acknowledgements text (verbatim, from the DESI page):

> This research used data obtained with the Dark Energy Spectroscopic Instrument (DESI). DESI construction and operations is managed by the Lawrence Berkeley National Laboratory. This material is based upon work supported by the U.S. Department of Energy, Office of Science, Office of High-Energy Physics, under Contract No. DE–AC02–05CH11231, and by the National Energy Research Scientific Computing Center, a DOE Office of Science User Facility under the same contract. Additional support for DESI was provided by the U.S. National Science Foundation (NSF), Division of Astronomical Sciences under Contract No. AST-0950945 to the NSF’s National Optical-Infrared Astronomy Research Laboratory; the Science and Technology Facilities Council of the United Kingdom; the Gordon and Betty Moore Foundation; the Heising-Simons Foundation; the French Alternative Energies and Atomic Energy Commission (CEA); the National Council of Humanities, Science and Technology of Mexico (CONAHCYT); the Ministry of Science and Innovation of Spain (MICINN), and by the DESI Member Institutions: www.desi.lbl.gov/collaborating-institutions. The DESI collaboration is honored to be permitted to conduct scientific research on I’oligam Du’ag (Kitt Peak), a mountain with particular significance to the Tohono O’odham Nation. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the U.S. National Science Foundation, the U.S. Department of Energy, or any of the listed funding agencies.

**HU.** A DESI-adatok CC BY 4.0 licenc alatt állnak: újraterjeszthetők és átdolgozhatók, ha megadjuk a DR1-adatközlő cikket, jelezzük a változtatásokat, és a fenti köszönetnyilvánítást szó szerint közöljük. A szöveg a szimulátor láblécében (lenyitható blokk) is szerepel.

## 2. SPARC

- Hivatalos oldal / Official page: https://astroweb.case.edu/SPARC/ (the old address http://astroweb.cwru.edu/SPARC/ returned HTTP 500 on 2026-10-30).
- Stated terms (verbatim, the only usage statement on the page): "If you are using data from this website in your scientific publications, please cite the SPARC master paper. References to the original HI+Halpha rotation curves are given in Table1.mrt. Please consider citing the original rotation-curve sources if relevant to your publication."
- **No explicit license is stated** on the official page, and it says nothing about redistribution.
- Master paper: F. Lelli, S. S. McGaugh, J. M. Schombert, AJ 152, 157 (2016), doi:10.3847/0004-6256/152/6/157; arXiv:1606.09251. The journal page carries "© 2016. The American Astronomical Society. All rights reserved."; this covers the article and states no data license.
- Zenodo copy: record https://zenodo.org/records/16284118 (doi:10.5281/zenodo.16284118, "SPARC Galaxy Rotation Curve", creators Lelli, McGaugh, Schombert) is labeled **CC BY 4.0**. It was created on 2025-07-21 and names "Original source: http://astroweb.cwru.edu/SPARC/". Its `Rotmod_LTG.zip` has the same MD5 (e4c8b927…) as the V26 input. We could not confirm that the SPARC authors made this deposit, so we do not rely on its CC BY 4.0 label as the authors' grant. **Status: unclear.**
- Use in V28: the 175 rotation-curve tables (galaxy distance D and R, V_obs, errV, V_gas, V_disk, V_bul) are embedded as numbers in `simulator/src/data.js` / `index.html`; the raw files are not shipped.

**HU.** A SPARC hivatalos oldala nem közöl kifejezett licencet, csak a fő cikk idézését kéri (és ajánlja az eredeti forgásigörbe-források idézését). Az újraterjesztésről nem nyilatkozik. Egy 2025-ös Zenodo-másolat CC BY 4.0 jelölést visel, de nem tudtuk megerősíteni, hogy azt a szerzők töltötték fel. Ezért a terjesztési státuszt tisztázatlannak tekintjük. Ha szükséges, a szerzőktől írásbeli megerősítés kérhető.

## 3. RSD fσ8 (12 values)

No data license is stated for these published measurements. The numbers are quoted with citation, as transcribed in the V27 code (`RSD_DATA`); in r3 they were not re-checked against the papers.

- 6dFGS: Beutler et al. 2012, MNRAS, doi:10.1111/j.1365-2966.2012.21136.x, arXiv:1204.4725
- SDSS MGS: Howlett et al. 2015, MNRAS, doi:10.1093/mnras/stu2693, arXiv:1409.3238
- BOSS DR12: Alam et al. 2017, MNRAS, doi:10.1093/mnras/stx721, arXiv:1607.03155
- WiggleZ: Blake et al. 2012, MNRAS, doi:10.1111/j.1365-2966.2012.21473.x, arXiv:1204.3674
- eBOSS DR16 LRG: Bautista et al. 2021, MNRAS, doi:10.1093/mnras/staa2800, arXiv:2007.08993
- eBOSS DR16 ELG: Tamone et al. 2020, MNRAS, doi:10.1093/mnras/staa3050, arXiv:2007.09009
- eBOSS DR16 QSO: Neveux et al. 2020, MNRAS, doi:10.1093/mnras/staa2780, arXiv:2007.08999
- FastSound: Okumura et al. 2016, PASJ, doi:10.1093/pasj/psw029, arXiv:1511.08083

## 4. Planck 2018

- Terms page: https://www.cosmos.esa.int/web/planck/planck-data-use ("Use of Planck data"). **No explicit data license is stated**; the page asks: "When publishing a paper using Planck data please include the following text, preferably as a footnote to the title:"

> Based on observations obtained with Planck (http://www.esa.int/Planck), an ESA science mission with instruments and contributions directly funded by ESA Member States, NASA, and Canada.

- Parameter source: Planck Collaboration, "Planck 2018 results. VI. Cosmological parameters", A&A 641, A6 (2020), doi:10.1051/0004-6361/201833910, arXiv:1807.06209 (Table 2, base_plikHM_TTTEEE_lowl_lowE_lensing, as noted in `tal_cosmology_common.py`). We do not use Planck Legacy Archive files; only published parameter values.
- `H0_LOCAL = 72.80 ± 1.00` is a rounded local-H0 prior and **not** a SH0ES data product (see the comment in `tal_cosmology_common.py`).

## 5. BBN

- PDG 2024 BBN review: B. D. Fields, P. Molaro, S. Sarkar, "24. Big Bang Nucleosynthesis", in S. Navas et al. (Particle Data Group), Phys. Rev. D 110, 030001 (2024), doi:10.1103/PhysRevD.110.030001; PDF https://pdg.lbl.gov/2024/reviews/rpp2024-rev-bbang-nucleosynthesis.pdf. Crossref lists the article license as https://creativecommons.org/licenses/by/4.0/.
- D/H: R. J. Cooke, M. Pettini, C. C. Steidel, ApJ 855, 102 (2018), doi:10.3847/1538-4357/aaab53. Crossref lists only the IOP copyright page; **no data license stated**.
- Y_p^std = 0.24709 and 10⁵ D/H^std = 2.5147 in `run_step1_5` are the benchmark's own algebraic reference values, not values taken from a BBN code release.

## 6. Nem külső adat / Not third-party data

`PANTHEON_GRID` (18 nodes, z = 0.015–2.260, σ_μ values) is a synthetic shape-control grid ("Pantheon+-like"). It contains no Pantheon+ (https://github.com/PantheonPlusSH0ES/DataRelease, no LICENSE file in the repository) or Union3 catalogue values, so no third-party data terms apply to it.

**HU.** A `PANTHEON_GRID` szintetikus alakkontroll-rács; Pantheon+ vagy Union3 katalógusértéket nem tartalmaz.

## Nyitott kérdések / Open points

1. SPARC: no explicit license on the official page; the CC BY 4.0 label of Zenodo 16284118 is not confirmed as the authors' own. / SPARC: a hivatalos oldalon nincs kifejezett licenc; a Zenodo 16284118 CC BY 4.0 jelölése nem igazoltan a szerzőktől származik.
2. Planck, RSD papers, Cooke 2018, ocean-gravity papers: no data license stated; we quote individual published numbers with citation. / Planck, RSD-cikkek, Cooke 2018, óceáni gravitációs cikkek: nincs adatlicenc; egyedi közölt számokat idézünk hivatkozással.
3. DESI: terms are explicit (CC BY 4.0). The likelihood files are taken from the CobayaSampler/bao_data repository, which the DESI page names as their public location; that repository has no LICENSE file of its own. / DESI: a feltételek egyértelműek (CC BY 4.0); a fájlok a DESI által megnevezett CobayaSampler/bao_data tárolóból származnak, amelynek nincs saját LICENSE-fájlja.
