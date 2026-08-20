# Phase 9B2R Corrected Single-Cell Cellular-Source Results

**Status:** CORRECTED_PRIMARY_RUN_READY_FOR_FULL_INDEPENDENT_REVIEW

Phase 9B2R corrects the failed Phase 9B2 analysis after Phase 9B2C. The execution scope remains PENG_CRA001160 only: CRA001160, BioProject PRJCA001063, Peng et al. 2019, 24 untreated PDAC tumors, 11 control pancreases, and 57,530 processed cells. No raw FASTQ or BAM files were downloaded. LIN_GSE154778, MONCADA_GSE111672, and HWANG_GSE202051 were not analyzed.

## Phase 9B2C Findings and Corrections
- FIND_01 (CRITICAL): CORRECTED_PHASE9B2R. Required controls are executed where applicable; module random controls are technically inapplicable after coverage exclusion.
- FIND_02 (MAJOR): CORRECTED_PHASE9B2R. All five modules are INSUFFICIENT_SINGLE_CELL_DATA and are not biological support.
- FIND_03 (MINOR): CORRECTED_PHASE9B2R. TF evidence categories are rule-derived with no blanket TO_VERIFY assignment.

## Preserved Versus Invalidated Results
Official provenance, cell/patient counts, metadata alignment, broad cell-type annotations, malignant-cell definitions, patient-aware pseudobulk construction, Moffitt50 scoring, Hallmark scoring, tumor-control context, and composition covariates were rerun or checksum-verified and preserved. Initial module-based biological support, module malignant-axis associations, placeholder negative-control conclusions, and blanket TF negative-control categories are invalidated.

## Module Coverage and Exclusion
- MEblack: coverage=0.2519, eligibility=INELIGIBLE, reason=INSUFFICIENT_SINGLE_CELL_DATA_LOW_COVERAGE_LT_0.80
- MEblue: coverage=0.3442, eligibility=INELIGIBLE, reason=INSUFFICIENT_SINGLE_CELL_DATA_LOW_COVERAGE_LT_0.80
- MEgreen: coverage=0.3864, eligibility=INELIGIBLE, reason=INSUFFICIENT_SINGLE_CELL_DATA_LOW_COVERAGE_LT_0.80
- MEtan: coverage=0.4854, eligibility=INELIGIBLE, reason=INSUFFICIENT_SINGLE_CELL_DATA_LOW_COVERAGE_LT_0.80
- MEgreenyellow: coverage=0.3193, eligibility=INELIGIBLE, reason=INSUFFICIENT_SINGLE_CELL_DATA_LOW_COVERAGE_LT_0.80

All transferred modules are classified as INSUFFICIENT_SINGLE_CELL_DATA. They remain only in descriptive coverage tables and do not enter formal source, malignant-axis, or support claims.

## Negative Controls
- TECHNICALLY_INAPPLICABLE: 10 rows
- EXECUTED: 64 rows

## Corrected TF Classifications
- CELL_COMPOSITION_EXPLAINED: 19
- NOT_SUPPORTED_AT_CELLULAR_LEVEL: 1
- PARTIAL_CELLULAR_SUPPORT: 1
- STROMAL_OR_IMMUNE_SOURCE_SUPPORTED: 4

## Corrected Malignant-Cell Axis Associations
- HALLMARK_PROTEIN_SECRETION: q=0.03361, direction=positive

All tested eligible features, including null results, are reported in `05_results/tables/phase9b2r_malignant_feature_axis_associations.tsv`.

## Corrected Cellular-Source Evidence
- CELL_COMPOSITION_EXPLAINED: 20
- INSUFFICIENT_SINGLE_CELL_DATA: 5
- MALIGNANT_CELL_INTRINSIC_SUPPORT: 1
- NOT_SUPPORTED_AT_CELLULAR_LEVEL: 1
- PARTIAL_CELLULAR_SUPPORT: 1
- STROMAL_OR_IMMUNE_SOURCE_SUPPORTED: 4

## Tumor-Control, Composition, and Null Findings
Control pancreases remain contextual only and do not redefine locked features. Composition-sensitive eligible features are recorded in `phase9b2r_cell_composition_sensitivity.tsv`. Null malignant-axis TF and Hallmark results are retained in the corrected association table.

## Boundary Conditions
This is a single-cohort cellular-source analysis. Ochrobactrum was not tested. No spatial validation, supplementary single-cell cohort analysis, microbiome validation, survival analysis, target prioritization, causal mediation, or manuscript writing was performed.

## Unresolved Items
No Phase 9B2R implementation blocker remains. The q < 0.10 malignant-axis reporting threshold is retained as the Phase 9B2C-reviewed reporting threshold rather than a new post hoc threshold.

## Review Readiness
Phase 9B2R is ready for complete independent review.
