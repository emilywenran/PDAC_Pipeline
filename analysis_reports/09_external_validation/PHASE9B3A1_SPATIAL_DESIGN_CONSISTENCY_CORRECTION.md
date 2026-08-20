# Phase 9B3A.1: Limited Spatial-Design Consistency Correction

This document details the prospective corrections to the Phase 9B3 spatial transcriptomic validation plan to resolve inconsistencies in covariates, model formulas, nesting hierarchies, sample counts, and execution structures.

---

## 1. Issue 1 Resolved: Covariate and Model-Formula Inconsistency

The prospective validation framework originally required adjustment for CAF, myeloid, and lymphoid cell-type fractions under Hypothesis 3, but the model formula omitted the lymphoid fraction. We have fully reconciled this:

### 1.1 Corrected Fully Mapped Model
The primary fully adjusted model is defined as:
$$\text{protein\_secretion\_score} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast} + \beta_2 \cdot \text{CAF\_fraction} + \beta_3 \cdot \text{myeloid\_fraction} + \beta_4 \cdot \text{lymphoid\_fraction} + \text{random\_effects}$$
*where the random effects structure depends on the dataset nesting.*

### 1.2 Reduced-Model Hierarchy
To prevent model overfitting, collinearity, or non-convergence, we prospectively lock a reduced-model hierarchy. Covariates must **not** be removed based on their P-values. Instead, a reduced model is permitted only under the following conditions:
1. **Unavailable Covariates:** If a specific cell type compartment is unestimable or absent in a dataset (e.g. lymphoid segments in specific GeoMx plates).
2. **Severe Collinearity:** Variance Inflation Factor (VIF) $> 10$ or condition number $> 30$.
3. **Model Non-Convergence:** The mixed-effects solver fails to converge after maximum iterations (or returns singular fit).

### 1.3 Reduced-Model Decision Hierarchy
- **Level 1 (Full):** $\text{Moffitt50\_contrast} + \text{CAF\_fraction} + \text{myeloid\_fraction} + \text{lymphoid\_fraction}$
- **Level 2 (No Lymphoid):** $\text{Moffitt50\_contrast} + \text{CAF\_fraction} + \text{myeloid\_fraction}$ (applied if lymphoid is unestimable or collinear)
- **Level 3 (Contrast Only):** $\text{Moffitt50\_contrast}$ (applied if multiple covariates cause singular fit or extreme collinearity)

---

## 2. Issue 2 Resolved: Patient, Section, ROI, and Spatial-Unit Nesting

We reject the incorrect term "mixed-effects OLS" and replace it with correct statistical terminology: **linear mixed-effects models (LMM)**. Spots or ROIs are nested within patients; the patient remains the unique biological replicate.

### 2.1 Cohort Nesting Structures

1. **`HWANG_GSE202051_NAIVE`** (NanoString GeoMx DSP)
   - Stated nesting: ROIs/segments nested within sections, which are nested within patients.
   - Stated parameters: 1 section per patient (so patient and section are 1-to-1).
   - Biological replicate: `patient_id` ($n=18$)
   - Inferential observation: `segment_id` (nested within patient)
   - Descriptive observation: Segment-level gene expression.
   - Random effect structure: `(1 | patient_id)`. The section-level random effect is not identifiable because there is exactly one section per patient.
2. **`HWANG_GSE202051_TREATED`** (NanoString GeoMx DSP)
   - Stated nesting: ROIs/segments nested within sections, which are nested within patients.
   - Stated parameters: 1 section per patient.
   - Biological replicate: `patient_id` ($n=25$)
   - Inferential observation: `segment_id`
   - Descriptive observation: Segment-level gene expression.
   - Random effect structure: `(1 | patient_id)`.
3. **`MONCADA_GSE111672`** (Microarray Spatial Transcriptomics ST)
   - Stated nesting: Spots nested within sections, which are nested within patients.
   - Stated parameters: Multiple sections exist per patient (Patient A: 4 sections, Patient B: 2 sections).
   - Biological replicate: `patient` ($n=2$)
   - Inferential observation: `spot_id` (nested within section and patient)
   - Descriptive observation: Spot-level gene expression.
   - Random effect structure: `(1 | patient_id) + (1 | patient_id:section_id)`.
   - *Patient-Aware Alternative:* Given the small patient sample size ($n=2$), fitting a mixed model with nested random intercepts and three covariates may cause non-convergence. In this case, we lock the alternative: perform patient-level regression for each patient independently, followed by fixed-effects meta-analysis of the coefficients.

---

## 3. Issue 3 Resolved: Disambiguating Patient and Section Counts

We have removed all combined "patients/sections" fields. Stated counts have been verified against NCBI GEO metadata:

* **`HWANG_GSE202051_NAIVE`**:
  - `patient_count`: 18
  - `section_count`: 18
  - `ROI_or_spot_count`: 256 segments (derived from 128 ROIs, each segmented into tumor and stroma AOIs)
  - `tumor_region_count`: 256
  - `control_or_adjacent_region_count`: 0
  - `treatment_group`: treatment-naïve
  - `platform`: NanoString GeoMx DSP
* **`HWANG_GSE202051_TREATED`**:
  - `patient_count`: 25
  - `section_count`: 25
  - `ROI_or_spot_count`: 352 segments (derived from 176 ROIs, each segmented into tumor and stroma AOIs)
  - `tumor_region_count`: 352
  - `control_or_adjacent_region_count`: 0
  - `treatment_group`: neoadjuvant-treated
  - `platform`: NanoString GeoMx DSP
* **`MONCADA_GSE111672`**:
  - `patient_count`: 2 (spatial)
  - `section_count`: 6
  - `ROI_or_spot_count`: 2248 spots
  - `tumor_region_count`: 2248
  - `control_or_adjacent_region_count`: 0
  - `treatment_group`: treatment-naïve
  - `platform`: Microarray Spatial Transcriptomics

---

## 4. Issue 4 Resolved: Dataset-Specific Execution Structure

To ensure methodological rigour:
1. **No Pooling:** The three canonical cohorts must not be pooled into a single expression matrix or primary model, as they use different platforms (GeoMx DSP vs. Microarray ST) and treatment groups.
2. **Analysis Roles:**
   - **Primary Spatial Cohort:** `HWANG_GSE202051_NAIVE` (determines baseline spatial localization of host programs).
   - **Treatment-Sensitivity Cohort:** `HWANG_GSE202051_TREATED` (determines if neoadjuvant therapy alters localization).
   - **Secondary Replication Cohort:** `MONCADA_GSE111672` (evaluates independent replication on a different platform).
3. **Cross-Cohort Synthesis Rule:** Synthesis is performed at the *coefficient level* after independent model fitting:
   - Direction-consistency requirement: The beta coefficients for `Moffitt50_contrast` in the primary model must have the same sign (positive) across naïve and replication datasets to support full replication.
   - Eligibility for meta-analysis: Only coefficients from treatment-naïve cohorts (`HWANG_GSE202051_NAIVE` and `MONCADA_GSE111672`) may be combined in a formal patient-level random-effects meta-analysis. Neoadjuvant-treated samples are excluded from baseline synthesis.

---

## 5. Verification Status

All validators have been updated and executed successfully:
- `06_scripts/python/16_validate_phase9b3a_spatial_plan.py`: **`PASS`**
- `06_scripts/python/15_validate_provenance_consistency.py`: **`PASS`**
- `06_scripts/python/00_validate_manifests.py`: **`PASS`**
