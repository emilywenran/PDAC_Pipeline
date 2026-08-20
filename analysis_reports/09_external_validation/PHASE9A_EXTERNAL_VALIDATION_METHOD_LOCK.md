# Phase 9A Prospective Method Lock: External Validation Framework

This document prospective locks the external-validation framework for the PDAC host transcriptional-state and host–microbiome findings. No modification of thresholds, models, or feature selections is permitted during the Phase 9B execution phase.

---

## 1. Context & Retained Discovery Findings

As required by the prospective protocol, this validation framework preserves and builds upon the conclusions of the discovery phases without reinterpretation:
- **Phase 3B Subtype Reproduction:** Subtype assignments GSE172356 primary (17 Basal, 23 Hybrid, 22 Classical) exactly reproduced public labels.
- **Phase 4B Discrete Subtype Stability:** Assignment stability conclusion is **`INCONCLUSIVE`**. Discrete clusters are scale- and feature-sensitive; a robust three-cluster (including Hybrid) partition is not stable across all representations.
- **Phase 5B Continuous Subtype Axis:** Continuous axis conclusion is **`INCONCLUSIVE`**. Hybrid samples represent intermediate/heterogeneous states rather than a distinct discrete biology. The continuous basal-classical contrast (`AXIS_MOFFITT50_PRIMARY`) is locked as the primary continuous host phenotype.
- **Phase 7C Independent Association Review:** **`PASS`**. 9 robust associated genera were identified under standard CLR.
- **Phase 8C Independent Host Mechanism Review:** **`PASS`**. The only taxon with robust, verified host-mechanism associations is **`Ochrobactrum`**, with 43 robust evidence rows:
  - **Pathway Layer:** `HALLMARK_PROTEIN_SECRETION` (negative association) and `HALLMARK_SPERMATOGENESIS` (negative association).
  - **Transcription Factor (TF) Layer:** 34 robustly associated TFs (viper activity).
  - **WGCNA Layer:** 7 robustly associated modules: `MEblack` (positive), `MEblue` (positive), `MEgreen` (negative), `MEgreenyellow` (positive), `MEpurple` (negative), `MEred` (negative), and `MEtan` (negative).

### Important Taxonomic Constraint
`Ochrobactrum` is a highly sensitive environmental and reagent contaminant. It must not be described as a confirmed tumor-resident organism. The validation framework is explicitly designed to distinguish biological replication from technical contamination or reagent-associated signal.

---

## 2. Validation Hierarchy & Modality Layers

Validation is organized into four separate, locked hierarchical layers:

### Layer 1: Independent Bulk-Transcriptome Validation
- **Primary Objective:** Validate the continuous basal-classical contrast, the 2 Hallmark pathways, the 34 TF activities, and the 7 WGCNA-derived co-expression signatures.
- **Secondary Objective:** Evaluate clinical outcome (survival) and treatment response associations.
- **Adjustment:** Host features must be tested for sensitivity to tumor purity and immune/stromal cell composition (using ESTIMATE or microdissection annotations) where metadata permit.

### Layer 2: Single-Cell Validation
- **Primary Objective:** Localize the cellular source of the 2 Hallmark pathways, 34 TFs, and 7 WGCNA-module signatures.
- **Key Questions:** Are these programs expressed in malignant epithelial cells, CAFs, endothelial, myeloid, or lymphoid cells? Is the continuous basal-classical score driven by malignant cells intrinsically or by stroma/composition? Do individual malignant cells exhibit basal/classical coactivation, or is the hybrid signature a bulk mixture artifact?

### Layer 3: Spatial-Transcriptomic Validation
- **Primary Objective:** Map the spatial coordinates of the basal-classical axis and localized host program activities.
- **Key Questions:** Do host signatures co-localize with malignant regions, CAFs, immune aggregates, or necrotic/ductal compartments? Does the association arise from physical co-localization of cell types rather than within-cell regulation?

### Layer 4: Independent Microbiome Validation
- **Primary Objective:** Determine if `Ochrobactrum` is detected in independent PDAC tissue sequencing datasets.
- **Key Questions:** Does detection survive host depletion, low-biomass filtering, and negative-control subtraction?
- **Paired Cohorts:** If paired host-microbiome datasets are available, test whether `Ochrobactrum` abundance is associated with the same host programs in the same direction.

---

## 3. Cohort Inclusion and Exclusion Criteria

To ensure validation rigor, candidate datasets must satisfy the following prospective criteria:

### Bulk Transcriptome (Layer 1)
- **Inclusion:** Human primary PDAC tumor specimens; processed sample-level gene expression matrix available; clinical annotations (survival/stage) available; sample size $n \ge 30$ (smaller cohorts are designated as underpowered secondary evidence).
- **Exclusion:** GSE172356 (discovery cohort); any cohort containing overlapping patients/samples; cell-line-only datasets.

### Single-Cell Transcriptomics (Layer 2)
- **Inclusion:** Human primary PDAC specimens; patient-level identifiers; cell-level expression matrix and cell-type annotations; minimum 5 patients and $\ge 500$ malignant epithelial cells.
- **Independence:** Statistical tests must treat the patient, not the individual cell, as the biological replicate (via pseudobulking) to avoid pseudoreplication.

### Spatial Transcriptomics (Layer 3)
- **Inclusion:** Human primary PDAC specimens; spot-level/cell-level expression matrix and corresponding 2D spatial coordinates; histological region annotations.

### Microbiome & Paired (Layer 4)
- **Inclusion:** Primary human PDAC tumor tissue sequencing; raw reads available (FASTQ/SRA).
- **Exclusion:** Stool, oral, pancreatic fluid, or cell-line microbiomes are prohibited as direct tumor tissue replication cohorts (contextual only).
- **Contamination Designation:** If sequenced extraction/reagent negative controls are absent, the dataset must be flagged as *contamination limited*.

---

## 4. Signature Transfer & Harmonization Plan

To prevent post hoc optimization, signatures will be transferred and harmonized across cohorts using the following locked protocols:

1. **Identifier Conversion:** Gene symbols will be converted to HGNC symbols. Outdated symbols will be mapped via NCBI Gene/Ensembl.
2. **Duplicate-Gene Handling:** Duplicate rows will be resolved by keeping the probe/transcript with the highest mean expression across samples.
3. **Missing-Gene Policy:** If a signature gene is missing in an external dataset:
   - If signature coverage is $\ge 80\%$, the score will be calculated on available genes and rescaled.
   - If signature coverage is $< 80\%$, the analysis will be flagged as *failed coverage feasibility* and not reported.
4. **Normalization:** Bulk datasets will be normalized within-cohort (e.g., log2(CPM+1) or log2(FPKM+1) for RNA-seq; median-centered/row-scaled for microarrays). Cross-cohort merging or batch correction is prohibited prior to within-cohort scoring.
5. **Scoring Methods:**
   - **Moffitt50 Axis:** Centroid-anchored reference score using the 50-gene signature.
   - **Hallmark Pathways:** Single-sample GSEA (ssGSEA) via `decoupleR` (minsize = 15).
   - **TF Activity:** VIPER algorithm via `decoupleR` using A/B/C confidence DoRothEA regulons.
   - **WGCNA Modules:** Module-eigengene-like standardized mean rank scores of the available module genes.

---

## 5. Statistical Validation Endpoints

### Continuous Axis and Host Phenotypes
Calculate Spearman rank correlation and OLS regression (with HC3 robust standard errors) between the primary continuous axis score and the sensitivity scores (singscore, PurIST, no-LEMD1). Test association of continuous axis position with survival using Cox proportional hazards models (reporting assumption diagnostics, continuous scores only, no data-driven cutpoint tuning).

### Host Mechanism Validation
For each supported pathway, TF, and module score, run:
$$\text{Feature\_Score} \sim \text{Continuous\_Axis\_Position} + \text{Purity\_Covariate}$$
Evaluate whether the effect direction is consistent with discovery and whether the 95% confidence interval excludes zero.

### Single-Cell Pseudobulking
Single-cell tests must use per-patient pseudobulk (average expression per patient within a cell type) to fit models, using patient as the biological unit:
$$\text{Pseudobulk\_Feature\_Score} \sim \text{Malignant\_Epithelial\_Purity} + (1 | \text{Patient\_ID})$$

### Spatial Autocorrelation
Spatial transcriptomics must adjust for spatial autocorrelation using spatial regression models (SAR/CAR) or spot-level ANOVA with patient blocking:
$$\text{Spot\_Score} \sim \text{Malignant\_Region} + \text{Stromal\_Region} + (1 | \text{Patient\_ID})$$

### Microbiome Validation
Microbiome validation requires:
1. Raw FASTQ quality control (fastp).
2. Host depletion (Bowtie2 mapping to GRCh38).
3. Taxonomic profiling (Kraken2/Bracken) against a standard RefSeq database.
4. Ochrobactrum species/genus abundance extraction and CLR re-normalization.
5. Contamination check: compare tumor abundance against negative-control runs (decontam or prevalence-difference tests).

---

## 6. Validation Evidence Categories

Replication outcomes for each feature must be classified into one of the following categories:

- **`EXTERNALLY_REPLICATED_HOST_FEATURE`:** Feature is available in an independent cohort; effect direction matches discovery; 95% CI excludes zero; replication is consistent across at least two independent bulk cohorts.
- **`PARTIALLY_REPLICATED_HOST_FEATURE`:** Direction is consistent but fails statistical significance, or is replicated in only one underpowered cohort.
- **`CELLULAR_SOURCE_SUPPORTED`:** Single-cell/spatial data localizes feature activity to a biologically plausible compartment (e.g., epithelial intrinsic for basal program).
- **`EXTERNALLY_REPLICATED_MICROBIOME_FEATURE`:** `Ochrobactrum` is detected in tumor tissue with sequencing negative control support, and direction is consistent in paired data.
- **`CONTAMINATION_LIMITED_MICROBIOME_EVIDENCE`:** `Ochrobactrum` is detected but negative controls are absent or show equal signal.
- **`NOT_REPLICATED`:** Adequate independent data show statistically significant opposite direction or null effect.
- **`INSUFFICIENT_EXTERNAL_DATA`:** Accession was not available, or gene/taxon coverage was below feasibility thresholds.

---

## 7. Negative-Control and Falsification Analyses

To ensure that statistical findings are not due to cohort-specific bias or code artifacts, Phase 9B must execute the following negative controls:
1. **Permutation Control:** Permute sample labels 1,000 times to calculate empirical null distributions for host-microbiome and host-subtype associations.
2. **Unrelated Pathway Control:** Score 5 unrelated MSigDB Hallmark pathways (e.g., `HALLMARK_MYOCARDIUM_DEVELOPMENT`, `HALLMARK_OLFACTORY_TRANSDUCTION`) which must not show significant association.
3. **Random Gene Modules:** Generate 100 random gene sets matched for size and mean expression to WGCNA modules; their associations must cluster around zero.
4. **Contaminant Control:** Run associations for 3 known sequencing contaminants (e.g., `Ralstonia`, `Novosphingobium`) that did not show association in discovery; they must remain non-significant.

---

## 8. Phase 9A Readiness Decision

Based on the audit of 12 candidate datasets and the verification of coverage and software resources, the readiness status is:

### **`READY_WITH_MICROBIOME_LIMITATIONS`**

- **Bulk-Host Readiness:** `READY_FOR_EXTERNAL_VALIDATION` (TCGA-PAAD, GSE71729, GSE62452 are fully qualified and locked).
- **Single-Cell Readiness:** `READY_WITH_DATASET_LIMITATIONS` (Updated under Phase 9A.3. Reconciled and locked four datasets for Phase 9B2 single-cell and spatial validation: PENG_CRA001160 as the sole active primary cohort for the Phase 9B2-primary execution; LIN_GSE154778 as secondary Layer 2; MONCADA_GSE111672 as Layer 3 and exploratory Layer 2; HWANG_GSE202051 as Layer 3 and secondary treatment-sensitivity Layer 2. The other single-cell datasets remain planned supplementary analyses, and separate later authorization is required for their execution. This sequencing decision is operational and does not imply that the supplementary datasets are scientifically unsuitable).
- **Spatial Readiness:** `READY_FOR_EXTERNAL_VALIDATION` (GSE202051, GSE274103, GSM3405527 are qualified and locked).
- **Microbiome Readiness:** `READY_WITH_MICROBIOME_LIMITATIONS` (PRJNA542615 and EGAS00001004572 are qualified, but PRJNA542615 lacks negative controls and EGAS00001004572 is controlled access, meaning microbiome replication remains contamination-sensitive).

---

## 9. Prospective Execution Authorization

**Phase 9B execution may proceed immediately under the updated Phase 9A.3 scope.** The primary execution for Phase 9B2 (single-cell) is restricted to PENG_CRA001160 only. Other cohorts (LIN_GSE154778, MONCADA_GSE111672, HWANG_GSE202051) are reserved for separate, subsequent supplementary execution. This operational sequencing preserves resource focus and is not a scientific downgrade.

*Locked on: 2026-07-03 (Phase 9A.3 correction)*
*Agent: Antigravity*

