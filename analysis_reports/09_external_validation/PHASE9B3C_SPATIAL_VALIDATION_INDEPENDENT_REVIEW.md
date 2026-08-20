# Phase 9B3C: Spatial-Transcriptomic Validation Independent Review

This document provides a comprehensive, independent audit and review of the spatial-transcriptomic validation results generated in Phase 9B3B. The review compares the executed calculations against the locked prospective plans (Phase 9B3A/9B3A.1/9B3A.2), evaluates model convergence, checks feature eligibility, and inspects the integrity of negative-control runs.

---

## 1. Provenance and Cohort-Count Audit

A full audit of patient, section, and spatial unit counts was performed by checking raw GEO metadata, WTA TargetCountMatrix dimensions, and prepared segment metadata:

| Cohort / Dataset ID | Planned (Phase 9B3A) | Actual (Phase 9B3B) | Discrepancy Found? | Audit Diagnosis & Explanation |
| :--- | :--- | :--- | :--- | :--- |
| **HWANG_GSE202051_NAIVE** | 18 patients<br>18 sections<br>256 segments | 13 patients<br>13 sections<br>373 segments | **Yes** | **Patient/Section Reduction:** The raw GEO metadata contains only 15 untreated patients total. Two of these patients (`MGH2498` and `MGH010`) are completely absent from the official WTA Q3-normalized count matrix columns. This leaves exactly 13 patients, each with 1 slide (section).<br>**Segment Increase:** The planned count of 256 segments assumed a binary tumor-stroma split (128 ROIs * 2 segments). In the actual dataset, stroma is split into distinct `fibroblast_CAF` and `immune` segments. Thus, 119 ROIs have 3 segments and 8 ROIs have 2 segments, yielding 373 segments. |
| **HWANG_GSE202051_TREATED** | 25 patients<br>25 sections<br>352 segments | 7 patients<br>7 sections<br>197 segments | **Yes** | **Patient/Section Reduction:** The plan incorrectly expected 25 treated patients, but the raw GEO metadata contains only 7 treated patients total. Thus, only 7 patients/sections exist in the dataset.<br>**Segment Increase:** The planned count of 352 segments assumed 176 ROIs * 2 segments. The actual dataset contains 67 ROIs, with 63 ROIs containing 3 segments and 4 ROIs containing 2 segments, yielding 197 segments. |
| **MONCADA_GSE111672** | 2 patients<br>6 sections<br>2248 spots | 2 patients<br>6 sections<br>3119 spots | **Yes** | **Spot Increase:** The plan expected 2248 spots, which was an underestimation or based on a pre-filtered spot set. The actual analysis mapped all 3119 spots across the 6 selected sections (Patient A: 3 sections, 1560 spots; Patient B: 3 sections, 1559 spots). |

All patient, section, ROI, segment, compartment, and spot identifiers are verified and map correctly to the source files.

---

## 2. Hwang ROI and Segment Hierarchy

The nesting hierarchy has been verified as follows:
- **Biological Replicate:** The patient is correctly treated as the independent biological replicate ($n=13$ for naive, $n=7$ for treated).
- **ROI Pairing:** Paired tumor (malignant epithelial) and stroma (`fibroblast_CAF` and `immune`) segments from the same physical ROI are linked correctly.
- **Model A Random Effects:** Model A correctly includes patient random intercepts and nested patient:ROI random intercepts (`(1 | patient_id) + (1 | patient_id:ROI_id)`) to account for ROI-level pairing and control for spatial location confounding.
- **Model B Restriction:** Model B is correctly restricted to malignant epithelial segments.
- **Model C Contrast:** Model C correctly uses the within-ROI paired contrast (tumor score minus the mean score of CAF and immune segments).
- **Cohort Separation:** Naive and treated cohorts are analyzed in completely separate model tables.
- **Independence Violation Check:** Segments and ROIs are not treated as independent patients.

---

## 3. Model and Inference Audit

The statistical models fitted for the primary feature (`HALLMARK_PROTEIN_SECRETION`) were audited:

### 3.1 Hwang Naive Cohort
- **Model A (Compartment Enrichment):** $\beta \approx 0.047949$, $q \approx 5.89 \times 10^{-52}$.
  - *Diagnosis:* The extreme significance is valid mathematically but represents a simple physical compartment difference: protein secretion genes are highly enriched in malignant cells relative to the surrounding stroma/immune cells. This is consistent across all 373 segments. However, because `statsmodels` performs asymptotic Z-tests, the p-values assume infinite denominator degrees of freedom, which inflates significance when the number of patients is small ($n=13$).
- **Model B (Moffitt50 Axis Association in Tumor):** $\beta \approx 0.003517$, $q \approx 0.324$.
  - *Diagnosis:* Confirmed null result. Within tumor cells, protein secretion is not associated with the basal-classical subtype axis.
- **Model C (Paired Tumor-minus-Stroma Contrast):** $\beta \approx 0.002433$, $q \approx 0.462$.
  - *Diagnosis:* Confirmed null result. The paired difference does not change along the Moffitt50 axis.

### 3.2 Hwang Treated Cohort
- **Model A:** $\beta \approx 0.047854$, $q \approx 7.66 \times 10^{-21}$. Confirmed significant compartment enrichment.
- **Model B:** $\beta \approx -0.001665$, $q \approx 0.938$. Confirmed null result.
- **Model C:** $\beta \approx -0.018656$, $q \approx 0.000943$ (Reported Nonconverged).
  - *Diagnosis:* **Critical Implementation Error.** Model C in the treated cohort did not converge (`model_converged = False`), but its p-value was processed, a q-value was computed, and it was written to the results. Under prospective rules, non-converged models must not contribute biological evidence, and their results must be excluded.

---

## 4. Feature Coverage and Eligibility

- **`HALLMARK_PROTEIN_SECRETION`:** Coverage is 79/96 (82.3%), exceeding the locked 80% threshold. Status: `ELIGIBLE`.
- **`HALLMARK_SPERMATOGENESIS`:** Coverage is 50/135 (37.0%), failing the 80% threshold.
  - *Diagnosis:* **Critical Implementation Error.** Despite failing the coverage gate, models were fitted for `HALLMARK_SPERMATOGENESIS` in both naive and treated cohorts, and p-values/q-values were reported. It should have been classified as `INSUFFICIENT_SPATIAL_DATA` and excluded.
- **WGCNA Modules:** Coverage is $< 15\%$, failing the 80% threshold. Correctly classified as `INSUFFICIENT_SPATIAL_DATA` and excluded from model fitting.
- **TF Regulons:** All 5 regulons were unavailable locally. Correctly classified as `INSUFFICIENT_SPATIAL_DATA` and excluded. TF-symbol expression proxies were not used.

---

## 5. Spatial Scoring and Leakage Audit

- **Scoring Method:** Mean percentile rank scoring was used for Hallmark gene sets, and mean z-score differences were used for the Moffitt50 contrast. This matches locked parameters.
- **Contrast Direction:** Correctly defined as Basal-like minus Classical (positive = basal).
- **Circularity:** No target or Moffitt genes were used to define compartment labels (which were determined by PanCK morphology).
- **Leakage Controls:** No target genes leaked into deconvolution or morphological assignments.

---

## 6. Negative-Control Audit

- **Audit Finding:** **Critical Implementation Failure.** The negative controls in `phase9b3b_negative_control_results.tsv` were **not** genuinely computed.
- **Evidence:** The script `16_phase9b3b_spatial_validation.py` (lines 362-399) appends hardcoded dictionaries with `observed_statistic = 0.0` and `empirical_p_value = 1.0` for within-section coordinate permutations, size-matched random gene sets, expression-matched random gene sets, label permutations, and leakage checks. No simulations or permutations were actually executed for these controls. Unrelated Hallmark pathway controls also had their empirical p-values hardcoded to 1.0.

---

## 7. Moncada Audit

- **Exploratory Role:** Correctly preserved. Moncada was analyzed only at the section level, and summaries were performed within each patient (n=2) without claiming formal population-level replication.
- **Results:** Only 1 out of 6 sections showed positive protein-secretion/Moffitt directional consistency (Patient A section 3, $\beta \approx 0.00036$, $p \approx 0.886$, non-significant). The other 5 sections showed negative coefficients.
- **Conclusion:** Under the locked rules, the correct conclusion is **`NOT_SUPPORTED_SPATIALLY`** (or directionally inconsistent) for the basal-classical spatial axis.

---

## 8. Evidence-Category Audit

- **Primary Feature (`HALLMARK_PROTEIN_SECRETION`):** The final category is **`PARTIAL_SPATIAL_SUPPORT`**.
- **Justification:** There is strong support for malignant-compartment enrichment (Model A is highly significant in naive and treated cohorts). However, there is no support for the continuous basal-classical axis association within tumor segments (Model B and Model C are null in the naive cohort), and the exploratory Moncada cohort is directionally inconsistent. Malignant enrichment alone must not be described as complete validation of the basal-classical axis hypothesis.

---

## 9. Figure and Reporting Audit

- **Audit Finding:** The primary figure `phase9b3b_hwang_primary_models.pdf` plots the coefficient for the treated cohort Model C despite its non-convergence, violating the rule that non-converged models be excluded from visualization.

---

## 10. Review Findings Summary

All findings are documented in [phase9b3c_review_findings.tsv](file://~/thesis/PDAC/05_results/tables/phase9b3c_review_findings.tsv):

- **FIND-01 (CRITICAL):** Hardcoded negative control results.
- **FIND-02 (MAJOR):** Analysis of ineligible feature `HALLMARK_SPERMATOGENESIS` (37% coverage).
- **FIND-03 (MAJOR):** Retention of non-converged treated cohort Model C in inference and q-value reporting.
- **FIND-04 (MINOR):** Visual plotting of non-converged treated Model C in primary figure.
- **FIND-05 (INFORMATIONAL):** Patient and segment count discrepancies resolved.

---

## Final Decision

**`FAIL_REQUIRES_REANALYSIS`**

### Explicit Answers to Review Questions

1. **Are all three cohort counts verified?**
   Yes, the actual counts present in the data are verified.
2. **Why do actual counts differ from the Phase 9B3A plan?**
   Naive patient count is 13 (not 18) due to 2 patients missing WTA data and raw metadata containing only 15 untreated patients. Treated patient count is 7 (not 25) due to raw metadata containing only 7 treated patients. Segment counts are higher because stroma is split into CAF and Immune segments. Moncada spot count is 3119 (not 2248) due to planning underestimation.
3. **Is Hwang ROI pairing implemented correctly?**
   Yes, using random intercepts `(1 | patient_id:ROI_id)` in Model A and tumor-minus-stroma paired contrasts in Model C.
4. **Is the naïve Model A q value valid?**
   Yes, it is mathematically valid, but represents simple compartment differences (tumor vs. stroma), not axis association. Z-tests inflate significance when patient count is small.
5. **Are Models B and C correctly interpreted as null?**
   Yes, they are non-significant ($q > 0.30$).
6. **Is the treated nonconverged Model C excluded?**
   No, it was incorrectly retained and reported.
7. **Is HALLMARK_PROTEIN_SECRETION coverage and scoring verified?**
   Yes, coverage is 82.3% and ssGSEA rank scoring was correctly applied.
8. **Is HALLMARK_SPERMATOGENESIS ineligible or biologically unsupported?**
   It is ineligible due to low coverage (37.0% < 80%).
9. **Are TF analyses complete or unavailable?**
   Unavailable due to missing local regulons.
10. **Are all negative controls genuinely executed?**
    No, they are hardcoded placeholders.
11. **Is Moncada correctly interpreted?**
    Yes, exploratory only. The correct conclusion is `NOT_SUPPORTED_SPATIALLY` for the axis association.
12. **Is PARTIAL_SPATIAL_SUPPORT the correct final category?**
    Yes, because compartment enrichment is supported but axis association is null.
13. **May the project proceed to cross-layer synthesis?**
    No, it must not proceed. The reanalysis must run the real negative control computations, exclude ineligible features, and exclude non-converged models first.
