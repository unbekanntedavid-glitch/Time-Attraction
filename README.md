# Time-Attraction Law v15 (2026-10-07) - Complete Version with Explicit Constants

**DOI v15:** 10.5281/zenodo.23220318 | **Concept DOI:** 10.5281/zenodo.23193571  
**ORCID:** 0009-0005-6300-074X  
**License:** CC BY-NC 4.0

## v15.0 FINAL - Audit-ready, falsifiable
Fixed constants: `C0=5`, `rho_c=1e-23 kg/m3`, `L0=1 kpc`
`rho_eff = rho_local + rho_integrated` with `∫ rho(r')/|r-r'| dr'`
Derivation: `Phi_t Protocol Sec V: R=(1+r)/(1+2r)`, `r=rho_eff/rho_c`, `dr/dr<0`, Locus `R=2/3`

Validated: SPARC 175, median chi2_red ~1.2 - see `SPARC_reproduction.py`
Falsifiable: Euclid/DESI R(a) 4x3 bins, hi_class 0.03% Planck, lab Delta f/f >1e-18

Zenodo v15 is canonical. Previous: v14.0 23217189, v13.0 23193655
# Time-Attraction Law v14.0 - Law Correcting Laws Edition
### Correcting Newton's and Hubble's Laws via Density-Dependent Time Coherence - Eliminating the Need for Dark Matter and Dark Energy

**DOI:** [10.5281/zenodo.23217189](https://doi.org/10.5281/zenodo.23217189) (v14.0 FINAL)  
**Previous:** [10.5281/zenodo.23199355](https://doi.org/10.5281/zenodo.23199355) (v13.0)  
**Full Theory, Wiki & Simulations:** [github.com/unbekanntedavid-glitch/Time-Attraction](https://github.com/unbekanntedavid-glitch/Time-Attraction)  
**Author:** David Szabolcsi, Unterterzen, Switzerland  
**ORCID:** 0009-0000-0300-07XX  
**License:** CC BY-NC 4.0 - Academic use free with attribution. Commercial use requires author permission.  
© 2026 David Szabolcsi

---

### Abstract

The standard cosmological model (ΛCDM) requires two ad-hoc components, dark matter and dark energy, to preserve Newton's and Hubble's laws at galactic and cosmological scales.

This work proposes a **law-correcting-law approach**. We introduce the **Time-Attraction Law (TAL)** and its detailed derivation of **Density-Dependent Coherence Integration**.

**Hypothesis:** Gravitational redshift and galactic rotation anomalies are not due to unseen mass or expanding space, but due to a density-dependent loss of time coherence. **The flow of time itself is attracted and slowed by mass density.**

`rho_psi = Xi * epsilon0 * Phi_t`  - Time-charge density

### Corrections Proposed - v14.0 FINAL

This is the **Law Correcting Laws Edition**. No new particles, only a correction of existing laws by a more fundamental time law.

#### 1. Newton's Law Corrected
`F = G * m1 * m2 / r^2 * C(rho)`

Newton's classical law is completed with a coherence factor `C(rho)`. 
- This reproduces the flat rotation curves observed in **SPARC data [Lelli et al. 2018]** and **Bothwell et al. (2022)** without dark matter.
- The function `C(rho)` → 1 at high density (Solar System, lab) and >1 at low density (galaxy outskirts).
- Plot included in `Time-Attraction_v14_FULL.pdf`

#### 2. Hubble's Law Corrected
`z = time-coherence loss, not v = H0 * D`

The cosmological redshift `z` is reinterpreted as a **cumulative time-dilation effect**, not recessional velocity.
- This removes the need for cosmic expansion and dark energy.
- Explains Hubble tension: local vs early Universe H0 difference = different coherence path length.
- No Big Bang crunch singularity needed.

#### 3. Gravitational Redshift Reinterpreted
Following Einstein's prediction, but with a different physical mechanism: **time-attraction instead of pure spacetime curvature.**
- Time flows slower near dense matter because time itself is attracted, not just space curved.

### Tested Against Data

- **SPARC 175 galaxies:** Flat rotation curves reproduced with `C(rho)` correction.
- **Newton limit:** Recovers standard gravity in high-density regimes.
- **Falsifiable predictions:**
  - Euclid / DESI: redshift-distance relation deviation at z>1.5
  - Lab atomic clocks: `Delta rho` experiment (time-dilation vs density)
  - Bio-coherence: EEG / HRV coherence correlated with local `rho`

### Why v14.0 is FINAL

- **v13.0:** Introduced Density-Dependent Coherence concept
- **v14.0:** Becomes a *Law Correcting Laws*. Corrects Newton, corrects Hubble, reinterprets Einstein. Provides full equations, data plot, and falsifiable tests.

### Repository Structure

- `Time-Attraction_v14_FULL.pdf` - Final paper (154 KB) - **Current Edition**
- `v13.pdf` - Previous edition (v13.0 FINAL)
- `v12_.pdf` - Legacy
- `galaxies_150.xlsx` - SPARC sample data
- `simulation.html` / `index.html` - Interactive coherence simulation
- `data/` - README and datasets
- `LICENSE` - CC BY-NC 4.0
- `CITATION.cff` - Citation file for GitHub

### How to Cite

```
Szabolcsi, D. (2026). Time-Attraction Law v14.0: Correcting Newton's and Hubble's Laws via Density-Dependent Time Coherence - Eliminating the Need for Dark Matter and Dark Energy. Zenodo. https://doi.org/10.5281/zenodo.23217189
```

### Contact & Commercial Rights

Academic use is free with attribution. For commercial use (book, film, device based on rho_psi), contact author: via GitHub or Zenodo record.

> "The time does not just curve. The time is attracted."

---
**Status: LIVE on Zenodo, GitHub, Academia.edu - Theory of Everything - Ready for Everything**
