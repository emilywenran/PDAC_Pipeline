# Phase 9B3A.2: Spatial Validation Hierarchy Final Correction

This document acts as a limited amendment to the existing Phase 9B3A/9B3A.1 spatial validation planning documents to lock the final statistical model structures and reclassify cohort roles.

---

## 1. Hwang GeoMx ROI Pairing

For both `HWANG_GSE202051_NAIVE` and `HWANG_GSE202051_TREATED`, the NanoString GeoMx DSP data has a paired design:
- Each patient provides exactly **one section**;
- Each section contains **multiple ROIs**;
- Each ROI is segmented into **paired tumor and stromal segments**.

To model this pairing and control for confounding by intra-tumor heterogeneity and spatial location:
1. **No Independent Segments:** Segments or ROIs must **not** be treated as independent biological replicates.
2. **Metadata Storage:** We store `patient_id`, `section_id`, `ROI_id`, `segment_id`, `compartment` (tumor vs. stroma), and `paired_segment_id` as distinct fields.

We lock three separate Linear Mixed-Effects Models (LMM):

### Model A: Compartment Comparison LMM
Tests the association of the feature score with the compartment (tumor vs. stroma), adjusting for patient subtype contrast and composition covariates, with nested random intercepts to account for ROI-level pairing:
$$\text{feature\_score} \sim \beta_0 + \beta_1 \cdot \text{compartment} + \beta_2 \cdot \text{Moffitt50\_contrast} + \beta_3 \cdot \text{CAF\_fraction} + \beta_4 \cdot \text{myeloid\_fraction} + \beta_5 \cdot \text{lymphoid\_fraction} + (1 \mid \text{patient\_id}) + (1 \mid \text{patient\_id:ROI\_id})$$
*where `(1 | patient_id:ROI_id)` models the ROI-level random intercept, pairing the tumor and stroma segments within each ROI.*

### Model B: Tumor-Segment-Only Axis LMM
Restricts the analysis to tumor/malignant segments only to assess axis association within the malignant cells:
$$\text{protein\_secretion\_score} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast} + \beta_2 \cdot \text{CAF\_fraction} + \beta_3 \cdot \text{myeloid\_fraction} + \beta_4 \cdot \text{lymphoid\_fraction} + (1 \mid \text{patient\_id})$$

### Model C: Paired Tumor–Stroma Contrast LMM
Models the difference between matched segments within each ROI as a function of continuous Moffitt50 contrast:
$$\text{tumor\_score\_minus\_stroma\_score} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast} + (1 \mid \text{patient\_id})$$

---

## 2. Moncada Low-Patient Limitation

`MONCADA_GSE111672` contains only 2 patients with spatial data (6 sections, 2248 spots). Due to this small sample size:
- **Reclassification:** We reclassify `MONCADA_GSE111672` as **`EXPLORATORY_CROSS_PLATFORM_SPATIAL_CONSISTENCY`**.
- **No Formal Replication:** It must **not** be treated as a formal population-level replication, and a two-patient meta-analysis must not be used as formal external replication.
- **Nesting Ban:** Spots or sections must **not** be treated as independent patients.

### 2.1 Locked Exploratory Analysis Protocol
1. **Section-Specific Analysis:** Run spatial association and autocorrelation analyses within each of the 6 tissue sections independently.
2. **Within-Section Permutations:** Run 1,000 spatial coordinate permutations for each section to establish empirical nulls.
3. **Section Summaries:** Compute summary statistics (effect size and direction) for each section.
4. **Direction Consistency:** Verify direction consistency across the sections and both patients (i.e. check if the association is consistently positive).
5. **Descriptive Comparison:** Perform a descriptive comparison of the localized programs with the primary `HWANG_GSE202051_NAIVE` results.

---

## 3. Final Evidence Hierarchy

To maintain scientific rigor:
- **Primary Spatial Inference:** `HWANG_GSE202051_NAIVE` (determines baseline baseline spatial localization of host programs).
- **Treatment-Sensitivity Inference:** `HWANG_GSE202051_TREATED` (determines if neoadjuvant therapy remodels localization).
- **Exploratory Cross-Platform Consistency:** `MONCADA_GSE111672` (cross-platform exploratory consistency only).

### 3.1 Matrix Pooling Ban
Because Hwang uses NanoString GeoMx DSP (targeted ROIs/segments) and Moncada uses first-generation ST (grid of spots), direct merging or pooling of count matrices, expression matrices, or spatial units across platforms is strictly **prohibited**.

---

## 4. Verification Status

All validator scripts have been successfully executed and pass:
- `06_scripts/python/16_validate_phase9b3a_spatial_plan.py`: **`PASS`**
- `06_scripts/python/15_validate_provenance_consistency.py`: **`PASS`**
- `06_scripts/python/00_validate_manifests.py`: **`PASS`**
