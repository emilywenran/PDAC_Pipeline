# Pancreatic Cancer Spatial Transcriptomic Validation Protocol

This document locks the execution protocol for Layer 3 (Spatial Transcriptomics) validation in Phase 9B3.

---

## Step 1: Data Acquisition & Setup
1. **Download processed data:** Download processed expression matrices and spatial coordinates for:
   - `HWANG_GSE202051` (GeoMx DSP count tables and metadata under GSE199102).
   - `MONCADA_GSE111672` (Spatial Transcriptomics coordinates and spot matrices).
2. **Directory Structure:** Save counts, coordinates, and metadata under `02_data/external/phase9_spatial/`.
3. **Index check:** Verify accession IDs map exactly to the publication record.

---

## Step 2: Quality Control (QC)
1. **Filtering:** Apply the following rules to each section independently:
   - For GeoMx DSP: Filter out ROIs/segments with $< 1000$ detected genes or $< 1000$ library size.
   - For Visium/ST: Filter out spots with $< 500$ detected genes or library size $< 1000$.
2. **Duplicate Coordinates:** Search and remove duplicate coordinate entries per section.
3. **Biological Block Check:** Ensure each section has its matched `patient_id` metadata. Reject sections with missing biological replicate links.

---

## Step 3: Spatial Mapping & Compartment Deconvolution

The mapping step must be executed using the following locked guidelines to prevent target-feature leakage and pseudoreplication:

### 1. Compartment Assignment
Each spatial unit (spot, cell, or ROI segment) must be assigned to one of the following major tissue compartments:
- **Malignant epithelial** (enriched for tumor cell signatures, PanCK+ in GeoMx)
- **Nonmalignant epithelial** (PanCK- normal ductal or acinar cells)
- **Fibroblast / CAF** (enriched for fibroblasts/stellate cells, *COL1A1*/*ACTA2*)
- **Endothelial** (enriched for blood vessel markers, *CD31*/*PECAM1*)
- **Myeloid** (macrophages/monocytes, *CD68*/*CD14*)
- **Lymphoid** (T/B/NK cells, *CD3D*/*CD19*)
- **Acinar** (*PRSS1*/*CEL*)
- **Other or Ambiguous** (mixed signatures or low confidence)

### 2. Leakage Protection (Circularity Ban)
To ensure unbiased validation, target genes must never be used to define histological regions or cell-type deconvolution models:
- Do **NOT** use `HALLMARK_PROTEIN_SECRETION` genes to define malignant spatial regions.
- Do **NOT** use Moffitt50 signature genes to train cell-type deconvolution models or perform reference-based spot decomposition (e.g. RCTD/MIA).
- Deconvolution models must be trained using independent marker panels or cell-type reference libraries.

### 3. WGCNA Module Coverage Gating
The five co-expression modules transferred from Phase 8 are gated:
- Check WGCNA module gene coverage in each dataset independently.
- If gene coverage is $< 80\%$, the WGCNA module must be classified as `INSUFFICIENT_SPATIAL_DATA` and barred from formal spatial inference.
- Only modules with $\ge 80\%$ independent gene coverage are allowed in the spatial sensitivity layer.

### 4. Patient-Aware Replication Constraint
- Spots, cells, segments, or sections must be treated as nested observations, not independent replicates.
- The **patient ID** is the sole biological unit of replication.
- Nesting & Design Rules:
  - **GeoMx DSP (Hwang):** One section per patient; multiple ROIs per patient; each ROI contains paired tumor and stromal segments. Do **not** treat segments or ROIs as independent. We store `patient_id`, `section_id`, `ROI_id`, `segment_id`, `compartment`, and `paired_segment_id` separately.
  - **Spatial Transcriptomics (Moncada):** Spots nested within sections nested within patients (n=2 patients). Do **not** treat spots or sections as independent patients.
- Spot-level or segment-level pseudoreplication is strictly prohibited.

---

## Step 4: Spatial Feature Scoring
1. **Scores Calculation:** Run ssGSEA or Rank score for pathway features. Run VIPER for TFs using regulon-based matrices.
2. **Gene Expression Check:** Do not score TFs using the expression of their own gene symbol.

---

## Step 5: Statistical Model Fitting

### 1. Hwang GeoMx Models (Naive & Treated)
Fit Linear Mixed-Effects Models (LMM) with nested random intercepts to account for pairing:
- **Model A (Compartment comparison LMM):**
  $$\text{feature\_score}_{ij} \sim \beta_0 + \beta_1 \cdot \text{compartment}_{ij} + \beta_2 \cdot \text{Moffitt50\_contrast}_{ij} + \beta_3 \cdot \text{CAF\_fraction}_{ij} + \beta_4 \cdot \text{myeloid\_fraction}_{ij} + \beta_5 \cdot \text{lymphoid\_fraction}_{ij} + (1 \mid \text{patient\_id}_i) + (1 \mid \text{patient\_id}_i:\text{ROI\_id}_j)$$
- **Model B (Tumor-segment-only axis LMM):**
  $$\text{protein\_secretion\_score}_{ij} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast}_{ij} + \beta_2 \cdot \text{CAF\_fraction}_{ij} + \beta_3 \cdot \text{myeloid\_fraction}_{ij} + \beta_4 \cdot \text{lymphoid\_fraction}_{ij} + (1 \mid \text{patient\_id}_i)$$
- **Model C (Paired tumor–stroma contrast LMM):**
  $$\text{tumor\_score\_minus\_stroma\_score}_{ij} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast}_{ij} + (1 \mid \text{patient\_id}_i)$$

### 2. Moncada ST Exploratory Protocol
Due to $n=2$ patients, do **not** use a two-patient meta-analysis as formal external replication. Lock the analysis as:
- **Section-Specific Analysis:** Run spatial models within each section independently:
  $$\text{protein\_secretion\_score}_{jk} \sim \beta_0 + \beta_1 \cdot \text{Moffitt50\_contrast}_{jk}$$
- **Within-Section Permutations:** Run 1,000 coordinate permutations per section to generate empirical nulls.
- **Section Summaries:** Compute summaries within each patient.
- **Direction Consistency:** Check direction consistency across sections and both patients.
- **Descriptive Comparison:** Perform descriptive comparison with Hwang naïve results.

### 3. Matrix Pooling & Direct Merge Ban
Because Hwang (GeoMx) and Moncada (ST) represent different platforms and spatial resolutions, direct merging or pooling of count/expression matrices, coordinates, or spatial units across platforms is strictly **prohibited**.

### 4. Reduced-Model Fallback
If a model fails to converge or exhibits severe collinearity (VIF $> 10$) in Hwang Model A/B, fit a reduced model from the hierarchy:
- *Level 1 (Full):* Contrast + CAF + Myeloid + Lymphoid
- *Level 2 (No Lymphoid):* Contrast + CAF + Myeloid
- *Level 3 (Contrast Only):* Contrast
- *LOPO Diagnostics:* Run Leave-One-Patient-Out to verify robustness.

---

## Step 6: Permutations and Multiple Testing
1. **Permutations:** Run 1,000 within-section coordinate permutations using seed `2026` to generate empirical null distributions.
2. **Correction:** Apply BH FDR correction within the locked testing families.

---

*Locked on: 2026-07-03*
*Lock Agent: Antigravity*
