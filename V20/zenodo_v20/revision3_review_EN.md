# Third revision — targeted independent review

**Date:** 2026-10-21. **Inputs:** `final_report (3).zip`, `report_hu (3).pdf`. **Requested checks:** the previous integral-matching error near L≈3 m and the treatment of the background magnetic-field cross term. The original files were not modified; no publication platform was modified. The Bothwell-gradient task remains separate.

*English translation of `revision3_review_HU.md`. Numerical results, qualifications, and conclusions are unchanged. Quoted passages from the Hungarian report and proposed replacement wording are translated into English below. No additional calculations were performed for this translation.*

## Verdict

**I accept the two requested corrections within the stated model and the newly specified scope of the report.**

1. The previous large numerical discrepancy near L≈3 m has been eliminated. Checking the new outer integral with a different angular quadrature confirms the correction.
2. The report now explicitly states the background-field cross term 2B_bg·B_coil, the conditions for ideal ±I averaging, and that the reported signal is the isolated coil-B² term that is even under polarity reversal.

The main signal remains **approximately −1.20×10⁻⁴⁸**. This is an acceptable calculated order of magnitude within the chosen ideal current-sheet model, **not a TAL detection and not a measured instrument sensitivity**.

Three minor or experimental-interpretation qualifications remain: “zero standard deviation” in the far-field region is not an error estimate; the explicitly disclosed small matching error at L=0.95 m remains; and a thermal control with the same total power does not automatically guarantee the same temperature field. These are distinguished below from the main issues that have already been corrected.

## 1. What was actually checked in this round?

- **All 18/18 SHA-256 entries match.** The separately uploaded PDF is byte-for-byte identical to the report inside the ZIP.
- I read the revised report, the new `energy_integrator.py` module, and the changes to `final_fixed.py`. The PDF tables were also checked against the rendered page images.
- I performed **targeted reruns** of the archived integrator at the main peak and of the outer A/B integration routes at the relevant L values.
- I checked the field using a separate circular-loop formula, the main signal using Cartesian Sobol cubature, and the far-field energy difference using a different angular and radial quadrature.

I do **not** claim a new complete execution of `final_fixed.py`: the full four-variant convergence analysis, static Monte Carlo calculations, and figure generation were not repeated in this round. The two requested corrections and the main peak were checked through targeted tests.

Archive SHA-256: `160981d3adb435b4fa719ca1a9d9a04d4e889a27f9d4c7bb0381e5c057ad2ea8`.

## 2. Main signal: reproducible and independently checkable

The new measurement record supplied by the author gives a winding diameter of 291.7 mm. The code now actually uses this value rather than deriving the model radius from resistance. It identifies the 300 mm dimension as the flange diameter. This resolves the previous ambiguity about the geometry at the documentation level. I treat the record as author-supplied measurement data, not as an on-site measurement independently certified by me.

| Quantity | Current check |
|---|---:|
| B at the main clock, 14 mA | **11.632900 µT** |
| B at the reference clock | **0.430539 µT** |
| Largest relative discrepancy in the independent circular-loop test, 5 points | **1.04×10⁻¹⁴** |
| Targeted rerun of the archived integrator, L=0.630957 m | **−1.2017108×10⁻⁴⁸** |
| Discrepancy relative to the archived main result | **2.2×10⁻¹⁶ relative**, at rounding level |
| Independent volumetric Sobol integral at the same L | **−1.2008435×10⁻⁴⁸** |
| Difference in independent signal magnitude relative to the package | **−0.0722%** |
| Standard error estimated from Sobol scrambles | **3.04×10⁻⁵²** |

The Sobol check uses 4×2¹⁹ points for the main sphere and 4×2¹² points for the smooth reference-sphere integral. The quoted standard error is not the uncertainty of the full model. It does not include uncertainties in coil geometry, current, the actual winding, or the environment. I do not identify the grid maximum with an accurately determined continuous maximum location.

The resistance calculated from the measured D is **1788.160 Ω**, differing from the stated 1789 Ω by **−0.0470%**: the report’s 0.05% comparison is correct. Five identical readings do indeed have zero sample standard deviation; this does not imply zero total measurement uncertainty in the dimension.

## 3. The correction near L≈3 m is acceptable

The previous **0.01268 J/A²** offset arising from separately estimated total energies no longer enters the outer branch near 3 m. The new route A uses panel integration plus a tail extending to infinity; route B uses the direct substitution r=L/t. Across L=2.6–3.01 m, the targeted rerun finds the two routes consistent to approximately **10⁻¹⁵ in relative terms**.

In addition, I calculated the difference between the outer energies using separate **Gauss–Legendre polar and radial quadrature with uniformly spaced azimuthal angles**. This is not the archived Fibonacci angular average. The infinite domain is likewise mapped using r=L/t; there is no finite truncation or radial interpolation. The field function is shared with the archive, but was checked at five points against a separate circular-loop integral. The coarse/fine grids of the separate integral differ by at most a few 10⁻⁶ in relative terms at the nearby points tabulated here; this is a convergence test, not a proof of an absolute error bound.

| L (m) | New archived integrator ΔW (J/A²) | Separate integral (J/A²) | Relative discrepancy |
|---:|---:|---:|---:|
| 2.6 | 0.00713721108 | **0.00713720352** | +0.000106% |
| 2.818383 | 0.00369084579 | **0.00369083957** | +0.000169% |
| 3.00 | 0.00233131494 | **0.00233130964** | +0.000227% |
| 3.162 | 0.00162105726 | **0.00162105265** | +0.000284% |
| 10 | 2.77763231×10⁻⁶ | **2.77754429×10⁻⁶** | +0.00317% |
| 31.62 | 8.37811251×10⁻⁹ | **8.37717664×10⁻⁹** | +0.0112% |
| 100 | 2.63649938×10⁻¹¹ | **2.63556282×10⁻¹¹** | +0.0355% |

The previous discrepancies of several hundred percent in the transition region and several percent in the far-field region have therefore been eliminated at the checked points. The small residual at 100 m decreases further with angular-quadrature refinement: using 6000/12000/24000 far-field angular points gives **2.6364994 / 2.6357023 / 2.6355308 ×10⁻¹¹ J/A²**, respectively. The radial t-quadrature uses 32 points, then 64 points. I do not claim a certified error bound for every L or for arbitrary distances.

### Two numerical qualifications remain

**A. Far-field “zero standard deviation” does not demonstrate zero error.** All four values in the far-field row of the convergence table are indeed identical. However, the four runs vary the high-resolution angular point count `nhi` and the inner panel spacing; the far-field settings `nlo=6000`, `nt=32` are unchanged in all four. Consequently, this row does not test resolution dependence in that region.

**Correction — high confidence, minor reporting clarification:** add the following note to the “0” entries in the far-field row:

> “The four runs use identical nlo and nt settings in this region. The zero sample standard deviation is therefore not a far-field error estimate; far-field convergence must be tested by varying nlo and nt separately.”

**B. The new, smaller matching error is disclosed and real.** The report now states that the inner/outer switch at L=0.95 m is not exactly continuous. The archived values are **2.17738947 and 2.17597607 J/A²**, differing by **0.00141340 J/A²**, or **0.0649%** of the inner value. This agrees algebraically with the discrepancy between the total-energy estimates; it is not a persistence of the old large error near 3 m. This was checked here using the archived values and the integral identity, not through a new total-energy integration.

**Qualification:** the small error near L=0.95 m must remain identified as a numerical limitation of the curve; the entire curve must not be described as continuous to machine precision. The document essentially already does this, so I do not treat it as a new issue overturning the main conclusion.

## 4. The background-field treatment is now appropriately qualified

**Location:** §4, polarity reversal. The added text, translated:

> The pure coil-B² term is even under polarity reversal. In a constant background field, the change in energy density also contains a 2B_bg·B_coil cross term; ideally, averaging equal-magnitude ±I states relative to OFF cancels this term. The reported −1.20×10⁻⁴⁸ value refers to this isolated even term, not to the full prediction for every possible simple ON/OFF measurement.

**I accept this correction.** Within the linear frequency-shift approximation used for the weak signal, with an unchanged background field and equal-magnitude ±I currents:

\[
\frac{1}{2}\{[|\mathbf B_{bg}+\mathbf b|^2-|\mathbf B_{bg}|^2]
+[|\mathbf B_{bg}-\mathbf b|^2-|\mathbf B_{bg}|^2]\}=|\mathbf b|^2.
\]

The code still calculates the pure coil-B² term. This now agrees with the explicitly specified computational scope of the text. It does not establish an actual background-field measurement, ±I symmetry, or a complete prediction for a simple ON/OFF experiment; those require experimental data. The first item in §5 should also use the wording “isolated even term” so that the condition is not lost in the final conclusion.

## 5. One remaining experimental overstatement: thermal compensation

**Location:** §4, item 5. Introducing a matched-heat-load control is appropriate, but the following conclusion is too strong. Translated passage:

> During the switched-off intervals, the same heat load (2.81 W) must be dissipated in a non-magnetic-field-generating (bifilar) resistor near the coils. Thus the ON/OFF difference contains no temperature step.

Equal total power does not guarantee identical spatial heat deposition or identical conductive and radiative coupling to the clocks. Waiting for settling can reduce time-dependent transients, but the temperatures of the two settled states may still differ. The residual field of the bifilar arrangement also needs to be checked; it is not necessarily exactly zero.

**Correction — high confidence, affects the evidential strength of the experiment:** replace the second sentence of the quotation with:

> “The purpose of the equal-total-power control with a small residual magnetic field is to reduce switching-correlated thermal changes. Temperature matching between the ON and OFF states and the residual field of the control must be verified by measurement; any remaining mismatch must be included in the uncertainty of the frequency estimate.”

This does not reject the thermal control. It replaces an assertion of guaranteed cancellation with a requirement to determine the measurable residual.

## 6. Other corrections and the limits of the measurement claim

**Accepted corrections:**

- The offset is not a noise threshold; the main-clock-equivalent magnetic coefficients are not measured upper bounds.
- The geometry measurement record and the distinction between the flange and winding surface have been added.
- The Sr value is identified as a systematic uncertainty, not a universal noise threshold, and the extraordinary white-FM time is identified as an algebraic extrapolation.
- The comparison between the case temperature and the datasheet specification is now conditional.
- The high-current calculation is no longer used to infer a safe pulse duration: the requested warning has been added.

**No new phase or frequency time series, ON/OFF/±I state log, or switching-correlated thermal measurements have been supplied.** The dimensional measurement record does not replace these. Actual sensitivity therefore remains undetermined. Conversion of the 65-minute cycle, provided it is a true 2π electrical phase cycle, remains **0.256410 mHz** and **2.56410×10⁻¹¹ mean fractional offset**.

Section 4 still ties the switching schedule to the beat period. This choice is not intrinsically invalid, but the offset does not make 65 minutes optimal; the interval must be justified through noise and thermal-settling tests. The next required evidence is an ON/OFF estimate and its uncertainty calculated from timestamped phase and state labels, not another calculated signal-to-beat ratio.

**Final assessment:** the two requested main corrections have been made; the calculated pure even field-energy signal of ≈−1.20×10⁻⁴⁸ is numerically supported. Subject to the remaining numerical limitations and qualifications concerning the control, this is a usable conditional design prediction, not a validated measurement result.

## Reproduction and accompanying data

`revision3_checks.json` contains the 18 integrity tests, the five separate field tests, the targeted main-peak rerun, the two Sobol evaluation points, nine outer-integral checks, the far-field resolution variation, and the arithmetic for the disclosed matching residual at 0.95 m. These are computational checks, not raw clock data.

Local verification code:

```bash
python3 user/workspace/co-scientist-outputs/ocxo-revision3-review/check_revision3.py
```

Python 3, NumPy, and SciPy are required. The code reads the ZIP, copies the two modules needed for execution into scratch space, and regenerates the results JSON. It does not overwrite the original package or manuscript and does not modify publication platforms.
