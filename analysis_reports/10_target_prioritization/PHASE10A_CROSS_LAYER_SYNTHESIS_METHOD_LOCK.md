# Phase 10A Cross-Layer Synthesis Method Lock

## 1. Cross-Layer Evidence Hierarchy

The synthesis of multi-layer evidence will classify candidate features into one of the following distinct categories:

1. **MULTI_LAYER_SUPPORTED**: Supported by primary discovery evidence, robust across microbiome sensitivities, replicated in bulk external cohorts, supported by single-cell malignant-cell intrinsic expression, and localized to the spatial malignant compartment.
2. **PARTIALLY_REPLICATED**: Shows robust discovery support and replicates in at least one external layer (bulk or single-cell), but lacks comprehensive multi-layer support.
3. **DISCOVERY_ONLY**: Shows robust host-mechanism support in the discovery cohort but fails to replicate externally.
4. **METHOD_SENSITIVE**: Primary association is sensitive to the choice of modeling or preprocessing method.
5. **COMPOSITION_SENSITIVE**: Primary association is driven primarily by cellular composition rather than true feature expression.
6. **CONTAMINATION_SENSITIVE**: Flagged as a potential contaminant with high risk of false positivity.
7. **NOT_EXTERNALLY_SUPPORTED**: Formally evaluated in bulk or single-cell datasets but failed external validation.
8. **INSUFFICIENT_DATA**: Excluded from formal inference due to coverage below locked thresholds (e.g., <80% coverage for WGCNA modules).
9. **EXPLORATORY_ONLY**: Findings derived exclusively from exploratory small-sample subsets.
10. **NO_SUPPORTED_ASSOCIATION**: Null or negative finding across primary validation analyses.

## 2. Objective Assignment Rules

* **MULTI_LAYER_SUPPORTED**: Primary = ROBUST AND Bulk = REPLICATED AND Single-Cell = MALIGNANT_CELL_INTRINSIC AND Spatial = SPATIAL_SUPPORT.
* **PARTIALLY_REPLICATED**: Primary = ROBUST AND (Bulk = REPLICATED OR Single-Cell = MALIGNANT_CELL_INTRINSIC) AND NOT MULTI_LAYER_SUPPORTED.
* **DISCOVERY_ONLY**: Primary = ROBUST AND Bulk = NOT_REPLICATED AND Single-Cell = NOT_SUPPORTED.
* **METHOD_SENSITIVE**: Assigned in Phase 7B/8B based on covariate/transformation tests.
* **COMPOSITION_SENSITIVE**: Assigned in Phase 9B2 based on pseudobulk fraction adjustment.
* **CONTAMINATION_SENSITIVE**: Assigned in Phase 6/7 based on negative controls.
* **NOT_EXTERNALLY_SUPPORTED**: Primary = ROBUST AND Feature is eligible AND validation yields NOT_REPLICATED.
* **INSUFFICIENT_DATA**: Feature coverage < 0.80 in external datasets.

## 3. Preserved Conclusions

* **Ochrobactrum** is the only taxon with robust host-mechanism associations.
* Microbial presence, localization, physical interaction, and causality were **not validated**.
* **HALLMARK_PROTEIN_SECRETION** has malignant-cell support and malignant-compartment spatial enrichment.
* Its basal–classical spatial-axis association was **not replicated**.
* Spatial evidence is classified as **PARTIAL_SPATIAL_SUPPORT**.
* Ineligible WGCNA modules (MEblack, MEblue, MEgreen, MEtan, MEgreenyellow) must **not be promoted**.
* Null and negative findings must remain in the synthesis.

## 4. Target Prioritization Framework

Targets will be prioritized using independent evidence sources. Literature support alone **cannot** rescue an unsupported computational candidate. Association is **not** causation.

* **Druggability**: Tractability via binding pockets or existing targeted therapies.
* **Genetic Dependency**: Essentiality in PDAC cell lines (e.g., DepMap).
* **Tumor-versus-Normal Selectivity**: Differential expression in tumor vs adjacent normal tissue.
* **Pathway Position**: Position in pathway hierarchy (e.g., Reactome hub or upstream regulator).
* **Cell-Type Specificity**: Expression restricted to specific cell populations.
* **External Replication**: Replicated in independent patient cohorts across validation layers.
* **Safety or Essentiality Concerns**: Potential off-target toxicity or pan-essentiality.
* **Existing Inhibitors or Compounds**: Availability of FDA-approved or clinical-stage compounds.
