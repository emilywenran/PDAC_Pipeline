# PDAC Cross-Layer Synthesis Protocol

## Overview
This protocol outlines the objective methodology for synthesizing evidence across the discovery (host-microbiome) and multi-layer validation (bulk, single-cell, spatial) phases. It aims to prevent narrative bias and ensure that targets promoted for further investigation have rigorous, multi-modal support.

## 1. Evidence Hierarchy

Candidates are categorized into a 10-level hierarchy based on the accumulation of evidence across layers.
1. `MULTI_LAYER_SUPPORTED`
2. `PARTIALLY_REPLICATED`
3. `DISCOVERY_ONLY`
4. `METHOD_SENSITIVE`
5. `COMPOSITION_SENSITIVE`
6. `CONTAMINATION_SENSITIVE`
7. `NOT_EXTERNALLY_SUPPORTED`
8. `INSUFFICIENT_DATA`
9. `EXPLORATORY_ONLY`
10. `NO_SUPPORTED_ASSOCIATION`

## 2. Preserved Constraints

* **Microbiome Causality**: Microbial presence, localization, interaction, and causality were NOT validated by this pipeline. Associations are purely statistical. Association does not imply causation.
* **Ineligible Features**: Features that fail coverage gates (e.g., <80% coverage) must be classified as INSUFFICIENT_DATA and cannot be promoted. This strictly applies to WGCNA modules (black, blue, green, tan, greenyellow).
* **Literature Rescue Prohibition**: Extensive literature backing cannot rescue a computational candidate that fails external validation. 
* **Null Findings**: Null and negative results must be retained and reported to prevent publication bias.

## 3. Specific Validated Features

* **Ochrobactrum**: The only genus retaining a robust host-mechanism association after rigorous sensitivity and contamination checks.
* **HALLMARK_PROTEIN_SECRETION**: Exhibits malignant-cell intrinsic support (single-cell) and malignant-compartment spatial enrichment. However, its association with the continuous basal-classical spatial axis was not replicated. Overall spatial evidence is PARTIAL_SPATIAL_SUPPORT.

## 4. Target Prioritization

Prioritization relies on orthogonal, independent data sources rather than circular reasoning. The framework incorporates:
* Druggability (OpenTargets)
* Genetic Dependency (DepMap)
* Tumor-vs-Normal Selectivity (GTEx, TCGA)
* Pathway Position (Reactome)
* Cell-Type Specificity (scRNA-seq atlases)
* Pharmacological availability (ChEMBL)
