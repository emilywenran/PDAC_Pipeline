# Phase 9B3A Prospective Method Lock: Spatial-Transcriptomic Validation Planning

This document prospectively locks the spatial-transcriptomic validation framework and statistical design for Phase 9B3 (Layer 3 Validation). No modifications to thresholds, model formulas, or feature sets are permitted during the Phase 9B3 execution phase.

---

## 1. Locked Inputs & Reconciled Single-Cell Context

This spatial validation plan builds directly upon the authoritative results from Phase 9B2R and the Phase 9B2C2 PASS review:
- **Primary Spatial Feature:** `HALLMARK_PROTEIN_SECRETION` is locked as the primary target. Corrected single-cell analysis classified this feature as `MALIGNANT_CELL_INTRINSIC_SUPPORT` and confirmed its positive association with the Moffitt50 basal–classical contrast at $q \approx 0.03361$.
- **Prespecified Biological Comparator:** `HALLMARK_SPERMATOGENESIS` is locked as a negative-control biological comparator (classified as `CELL_COMPOSITION_EXPLAINED` in single-cell validation).
- **Secondary Localization Features:** Five TF regulons identified with stromal, immune, or partial cellular support are transferred:
  - `ELF1` (Stromal/Immune)
  - `MBD2` (Stromal/Immune)
  - `ZBTB7A` (Stromal/Immune)
  - `ZNF384` (Stromal/Immune)
  - `ZNF740` (Partial Cellular Support)
- **WGCNA Sensitivity Features:** The five WGCNA modules (`MEblack`, `MEblue`, `MEgreen`, `MEtan`, `MEgreenyellow`) had coverages $< 49\%$ in single-cell data and were excluded as `INSUFFICIENT_SINGLE_CELL_DATA`. They are barred from formal spatial inference unless a dataset-specific coverage check independently exceeds the locked 80% threshold.

> [!IMPORTANT]
> **Scientific Boundary:** Spatial host expression profiling evaluates the localization of human transcriptional programs. It does not validate the presence, abundance, or causality of *Ochrobactrum* or other microbes. Host spatial transcripts must not be used to infer the presence or location of microorganisms.

---

## 2. Authoritative Spatial Cohort Set & Qualification

Only spatial datasets qualified and approved in Phase 9A are evaluated. No new cohorts are authorized.

- **Primary Spatial Cohort:** `HWANG_GSE202051_NAIVE` (NanoString GeoMx DSP, 18 patients, 18 sections, 256 segments, treatment-naïve).
- **Exploratory Cross-Platform Spatial Consistency Cohort:** `MONCADA_GSE111672` (Microarray Spatial Transcriptomics, 2 patients, 6 sections, 2248 spots, treatment-naïve). Reclassified to reflect low sample size (n=2 patients) and prevent treatment as formal population-level replication.
- **Treatment-Sensitivity Cohort:** `HWANG_GSE202051_TREATED` (NanoString GeoMx DSP, 25 patients, 25 sections, 352 segments, neoadjuvant-treated).
- **Exploratory / Non-Authorized Cohorts:** `GSE274103` (Wang et al. 2024, 10x Visium) and `GSE272362` (Zhang et al. 2024, 10x Visium) are qualified but **NOT** authorized for execution in this phase (`current_execution_authorized = FALSE`).

---

## 3. Statistical Analysis Unit

The patient is the independent biological replicate ($n = 18$ for naïve primary, $n = 25$ for treated primary, and $n = 2$ for Moncada).
- **Nesting Structure & Design:**
  - **HWANG_GSE202051_NAIVE/TREATED (GeoMx DSP):** One section per patient; multiple ROIs per patient; each ROI contains paired tumor (PanCK+) and stromal (PanCK-) segments. Treating segments or ROIs as independent biological replicates is strictly prohibited.
  - **MONCADA_GSE111672 (ST):** Spots (inferential observation) nested within sections, which are nested within patients (Patient A: 4 sections, Patient B: 2 sections). Spots or sections must NOT be treated as independent patients.
- **Metadata Fields:** We store `patient_id`, `section_id`, `ROI_id`, `segment_id`, `compartment`, and `paired_segment_id` as separate fields to preserve nesting structure.
- **Replication Thresholds:** 
  - Minimum patients per cohort: $n \ge 2$
  - Minimum tissue sections per patient: $\ge 1$
  - Minimum usable spatial units (spots/ROIs) per section: $\ge 20$
  - Minimum malignant-enriched spatial units per patient: $\ge 5$

---

## 4. Spatial Quality Control (QC)

Spatial QC thresholds must be applied project-wide without adjusting for association results:
1. **Minimum Detected Genes:** $\ge 500$ genes per spot/segment (Visium) or $\ge 1000$ genes (GeoMx).
2. **Minimum Library Size:** $\ge 1000$ counts per spot/segment.
3. **Mitochondrial Fraction:** Maximum 15% (Visium). Not applicable for FFPE-based GeoMx WTA.
4. **Tissue Boundaries:** Proactive removal of off-tissue or background spots using histology alignment.
5. **Coordinate Integrity:** Check and fail if duplicate coordinates exist on the same slide.
6. **Exclusion Rules:** Exclude sections where $>80\%$ of spots fail QC; exclude patients if all sections fail.

---

## 5. Cell-Type Localization & Deconvolution

To prevent target-feature leakage and circularity:
- **Circularity Ban:** Do **NOT** use `HALLMARK_PROTEIN_SECRETION` genes or `Moffitt50` signature genes to perform cell-type deconvolution or define histological regions.
- **Independent Panels:** Use independent major cell-type marker panels (e.g. *EPCAM*/*KRT19* for epithelial, *COL1A1*/*ACTA2* for CAFs, *CD68* for myeloid, *CD3D* for T cells).
- **GeoMx DSP Localization:** Utilize PanCK fluorescence-guided morphology segmentation (PanCK+ tumor segments vs. PanCK- stroma segments) to physically isolate compartments.
- **Visium Deconvolution:** Apply RCTD (Robust Cell Type Deconvolution) or MIA (Multimodal Intersection Analysis) using matched single-cell reference datasets.
- **Major Spatial Compartments:** Mapped features must be assigned to one of the following:
  1. Malignant epithelial
  2. Nonmalignant epithelial
  3. Fibroblast or CAF
  4. Endothelial
  5. Myeloid
  6. Lymphoid
  7. Acinar
  8. Other or Ambiguous

---

## 6. Spatial Feature Scoring

Scores will be calculated using the following locked parameters:
- **Algorithms:** ssGSEA via decoupleR for Hallmark pathways; VIPER (using DoRothEA A/B/C regulons) for TF activities; standardized mean rank scores for WGCNA modules.
- **Coverage Gate:** Minimum 80% gene/regulon target coverage is required. WGCNA modules must not enter formal spatial models if their coverage in the dataset falls below 80%.
- **TF Scoring Rule:** TFs must be scored based on their regulon targets. Using the expression level of the TF gene symbol itself is prohibited.

---

## 7. Primary Spatial Hypotheses

### Hypothesis 1 (Compartment Enrichment)
`HALLMARK_PROTEIN_SECRETION` activity is enriched in malignant epithelial compartments relative to nonmalignant, stromal, and immune compartments.
- **Null Hypothesis ($H_{0,1}$):** Mean activity difference is zero.
- **Decision Boundary:** $q < 0.05$ (enrichment in tumor).

### Hypothesis 2 (Malignant-Axis Association)
Within malignant-enriched spatial units, `HALLMARK_PROTEIN_SECRETION` is positively associated with the Moffitt50 basal–classical contrast in the direction observed in Phase 9B2R (higher in basal).
- **Null Hypothesis ($H_{0,2}$):** The coefficient for the basal-classical contrast is $\le 0$.
- **Decision Boundary:** $q < 0.05$ (positive coefficient).

### Hypothesis 3 (Composition Adjustment)
The association between `HALLMARK_PROTEIN_SECRETION` and the Moffitt50 contrast remains significant after adjusting for major non-epithelial cell-type fractions (CAF, myeloid, lymphoid).
- **Null Hypothesis ($H_{0,3}$):** The adjusted coefficient for the contrast is $\le 0$.
- **Decision Boundary:** $q < 0.05$.

---

## 8. Spatial Statistical Models

Linear mixed-effects models (LMM) will be fitted for GeoMx DSP cohorts, treating the patient and ROI as nested random intercepts:

### 8.1 Hwang GeoMx Locked Models (Naive & Treated)

#### Model A: Compartment Comparison LMM
$$\text{feature\_score}_{ij} \sim \beta_0 + \beta_1 \cdot \text{compartment}_{ij} + \beta_2 \cdot \text{Moffitt50\_contrast}_{ij} + \beta_3 \cdot \text{CAF\_fraction}_{ij} + \beta_4 \cdot \text{myeloid\_fraction}_{ij} + \beta_5 \cdot \text{lymphoid\_fraction}_{ij} + (1 \mid \text{patient\_id}_i) + (1 \mid \text{patient\_id}_i:\text{ROI\_id}_j)$$
*where `(1 | patient_id:ROI_id)` models the ROI-level random intercept, pairing the tumor and stroma segments within each ROI to control for spatial location confounding.*

#### Model B: Tumor-Segment-Only Axis LMM
Restricted to tumor/malignant segments only:
$$\text{protein\_secretion\_score}_{ij} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast}_{ij} + \beta_2 \cdot \text{CAF\_fraction}_{ij} + \beta_3 \cdot \text{myeloid\_fraction}_{ij} + \beta_4 \cdot \text{lymphoid\_fraction}_{ij} + (1 \mid \text{patient\_id}_i)$$

#### Model C: Paired Tumor–Stroma Contrast LMM
Models the difference between matched segments within each ROI:
$$\text{tumor\_score\_minus\_stroma\_score}_{ij} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast}_{ij} + (1 \mid \text{patient\_id}_i)$$

### 8.2 Moncada ST Exploratory Protocol
Due to $n=2$ patients, Moncada ST cannot be analyzed using formal population LMM or patient-level meta-analysis. We lock the following exploratory protocol:
1. **Section-Specific Analysis:** Run spatial association models independently within each of the 6 sections:
   $$\text{protein\_secretion\_score}_{jk} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast}_{jk}$$
2. **Within-Section Permutations:** Run 1,000 spatial coordinate permutations per section to generate empirical nulls.
3. **Section Summaries:** Compute effect size and direction summaries within sections for each patient.
4. **Direction Consistency:** Verify direction consistency across all sections and both patients.
5. **Descriptive Comparison:** Perform descriptive comparison with the primary `HWANG_GSE202051_NAIVE` results.

### 8.3 Matrix Pooling & Direct Merge Ban
Because Hwang (targeted GeoMx DSP) and Moncada (grid-based ST) represent different platforms and resolutions, direct merging or pooling of count/expression matrices, coordinates, or spatial units across platforms is strictly **prohibited**.

### 8.4 Reduced-Model Hierarchy
To prevent collinearity, unestimable covariates, or non-convergence in Hwang Model A/B, we lock the following reduced-model hierarchy (covariates must NOT be removed based on P-values):
- **Level 1 (Full):** $\text{Moffitt50\_Contrast} + \text{CAF\_Fraction} + \text{Myeloid\_Fraction} + \text{Lymphoid\_Fraction}$
- **Level 2 (No Lymphoid):** $\text{Moffitt50\_Contrast} + \text{CAF\_Fraction} + \text{Myeloid\_Fraction}$ (used if lymphoid is collinear or unestimable)
- **Level 3 (Contrast Only):** $\text{Moffitt50\_Contrast}$ (used if severe multicollinearity or singular fit occurs)

### Sensitivity Models
- **Section-Level model:** Adjusting for section-to-section variation:
$$\text{Feature\_Score}_{ijk} \sim \beta_1 \cdot \text{Moffitt50\_Contrast}_{ijk} + (1 | \text{Patient\_ID}_i) + (1 | \text{Patient\_ID}_i:\text{Section\_ID}_j)$$
- **Leave-One-Patient-Out (LOPO) diagnostics:** Run to verify that no single patient drives the association.
- **Standard Errors:** Robust sandwich standard errors (HC3-based) at the patient level.

---

## 9. Spatial-Autocorrelation & Co-Localization

- **Autocorrelation:** Moran's I will be calculated for spot-level feature activities within each tissue section. Permutations (1,000 runs) will be used to generate empirical nulls.
- **Co-localization:** Test spatial co-localization of `HALLMARK_PROTEIN_SECRETION` with CAFs (*COL1A1*), myeloid cells (*CD68*), and Moffitt50-high regions.
- **Interpretation Guardrail:** Co-localization must not be used to infer direct cell-to-cell signaling or physical interaction without independent verification.

---

## 10. Negative-Control Plan

To prevent false positives, the following negative controls must be executed:
1. **Size-matched random gene sets (100):** Randomly selected genes scored similarly to WGCNA/Hallmark features; pooled coefficients must cluster around zero.
2. **Unrelated Hallmark pathways (5):** Score `HALLMARK_MYOCARDIUM_DEVELOPMENT`, `HALLMARK_OLFACTORY_TRANSDUCTION`, `HALLMARK_BILE_ACID_METABOLISM`, `HALLMARK_PANCREATIC_BETA_CELLS`, and `HALLMARK_HEME_METABOLISM` which must not show significant malignant enrichment.
3. **Permutations (1,000 runs):** Permute coordinates within-section, patient labels, and deconvolution fractions to establish empirical null distributions.
4. **Seed:** All randomizations are locked to seed `2026`.

---

## 11. Multiple-Testing Framework

Benjamini-Hochberg (BH) FDR corrections will be applied within the following independent families:
- **Family 1 (Primary Hypotheses):** 3 tests ($q < 0.05$ threshold).
- **Family 2 (Secondary TFs):** 5 tests ($q < 0.10$ threshold).
- **Family 3 (Spatial Autocorrelations):** Feature-wise tests per dataset ($q < 0.05$).
- **Family 4 (Co-localizations):** Feature-compartment tests ($q < 0.05$).
- **Family 5 (WGCNA modules):** Gated sensitivity tests ($q < 0.10$).

---

## 12. Spatial Evidence Categories

The final spatial evidence categories are locked:
- **`SPATIAL_MALIGNANT_LOCALIZATION_SUPPORTED`:** Mapped to malignant compartment, significant in primary models.
- **`SPATIAL_STROMAL_OR_IMMUNE_LOCALIZATION_SUPPORTED`:** Mapped to non-malignant/stromal/immune compartments.
- **`SPATIAL_COMPOSITION_CONSISTENT`:** Spatial distribution matches deconvolution fractions but lacks independent axis association.
- **`SPATIAL_AXIS_ASSOCIATION_SUPPORTED`:** Positive association with Moffitt50 contrast in malignant-enriched units.
- **`PARTIAL_SPATIAL_SUPPORT`:** Direction matches but fails significance, or supported in only one dataset.
- **`NOT_SUPPORTED_SPATIALLY`:** Significant opposite direction or no enrichment.
- **`INSUFFICIENT_SPATIAL_DATA`:** Excluded due to dataset-specific coverage $< 80\%$.
- **`TO_VERIFY`:** Flagged for quality, batch, or parser issues.

---

## 13. Resource and Environment Specifications

- **Key Packages:** `decoupleR` (v2.12.0), `StandR` (v1.10.0), `Seurat` (v5.1.0), `statsmodels` (v0.14.6)
- **MacBook Suitability:** Confirmed practical on Apple Silicon (peak RAM $< 16$ GB, run time $< 30$ mins per cohort).
- **HPC Requirement:** Not required for planned processed-matrix runs.

---

*Locked on: 2026-07-03*
*Lock Agent: Antigravity*
