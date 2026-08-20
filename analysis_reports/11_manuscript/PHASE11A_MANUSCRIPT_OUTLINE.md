# Phase 11A Manuscript Outline

## Abstract
* Background: The tumor microbiome in Pancreatic Ductal Adenocarcinoma (PDAC) has potential interactions with host state.
* Methods: Integrative multi-omics analysis of bulk, single-cell, and spatial transcriptomics data.
* Results: 
    * Global community-level association exists ($R^2 = 0.0534$).
    * 9 robust genera (e.g., Ochrobactrum) associated with host states.
    * Validation across multiple external cohorts shows replication of transcription factor activities (e.g., CTCFL/BORIS) but with composition-sensitive caveats.
    * HALLMARK_PROTEIN_SECRETION identified as a robust malignant-cell intrinsic signature with spatial malignant compartment support, despite lack of continuous basal-classical spatial replication.
* Conclusions: Tumor-associated microbiota exhibit robust but non-causal statistical associations with specific host transcriptional programs.

## Introduction
* Role of the microbiome in PDAC.
* Focus on tumor-intrinsic microbiomes and host transcriptional axes (e.g., basal-classical).
* Rationale for multi-layer validation (bulk, single-cell, spatial).

## Results
### 1. Tumor Microbiome Associates with Host Transcriptional State
* Identification of 9 robust genera.
* Global PERMANOVA community-level association.
* Explicit reporting of transformation sensitivities (rCLR direction reversals).

### 2. Microbial Associations Link to Host Biological Mechanisms
* 43 robust host mechanisms linked to Ochrobactrum.
* Target prioritization according to the Phase 10C2 framework.

### 3. External Bulk Transcriptomic Validation
* Replication of 12 TF activities.
* Identification of partially replicated features (7 pathways/modules and 13 TFs).
* Clear reporting of non-replicated features (e.g., HALLMARK_SPERMATOGENESIS, WGCNA modules MEred, MEpurple).

### 4. Cellular Source and Composition Sensitivity
* Single-cell analysis resolving cellular sources.
* Identification of HALLMARK_PROTEIN_SECRETION as malignant-cell intrinsic.
* Clear penalization of composition-sensitive signatures (e.g., CTCFL/BORIS correctly penalized as composition-explained).

### 5. Spatial Organization and Compartment Support
* Partial spatial support for host features (HALLMARK_PROTEIN_SECRETION enriched in malignant compartment).
* Spatial validation of the basal-classical axis is only PARTIAL_SPATIAL_SUPPORT.
* Exploratory findings in Moncada cohort (inconsistent spatial correlation).

## Discussion
* Strict associative language: the microbiome is associated with, but cannot be claimed to cause, host transcriptional changes.
* Acknowledgment of contamination risks (Herbaspirillum).
* Acknowledgment of lack of direct microbial localization/physical interaction data (no microbial localization is claimed due to lack of sequenced negative controls).
* Importance of rigorous statistical controls (composition sensitivity, multi-layer validation) in computational microbiome analyses.

## Methods
* Detailed cohort inclusion, preprocessing, normalization (CLR/rCLR), modeling (OLS with HC3), multiple testing correction (BH).
* WGCNA, TF activity (DoRothEA/VIPER), Pathway activity (MSigDB/PROGENy).
* Strict negative control iterations and validations.
