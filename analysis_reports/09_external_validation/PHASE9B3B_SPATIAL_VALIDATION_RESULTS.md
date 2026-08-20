# Phase 9B3B Spatial-Transcriptomic Validation Results

**SUPERSEDED_BY_PHASE9B3R**

This report is preserved as audit history only. The corrected reanalysis is `PHASE9B3R_CORRECTED_SPATIAL_VALIDATION_RESULTS.md`.

Phase 9B3B executed the locked spatial validation on authorized cohorts only. GeoMx and ST matrices were not pooled.

## Datasets

| dataset_id              | accession           | publication                                            |   patient_count |   section_count |   ROI_count |   segment_count |   spot_count | matrix_orientation   | expression_value_type                                        | treatment_status_verified   | passes_locked_qc   |
|:------------------------|:--------------------|:-------------------------------------------------------|----------------:|----------------:|------------:|----------------:|-------------:|:---------------------|:-------------------------------------------------------------|:----------------------------|:-------------------|
| HWANG_GSE202051_NAIVE   | GSE202051/GSE199102 | Hwang et al. 2022/2023; DOI 10.1038/s41588-023-01411-z |              13 |              13 |         127 |             373 |            0 | genes_by_segments    | Q3-normalized GeoMx WTA counts, log2 transformed for scoring | Untreated                   | True               |
| HWANG_GSE202051_TREATED | GSE202051/GSE199102 | Hwang et al. 2022/2023; DOI 10.1038/s41588-023-01411-z |               7 |               7 |          67 |             197 |            0 | genes_by_segments    | Q3-normalized GeoMx WTA counts, log2 transformed for scoring | Treated                     | True               |
| MONCADA_GSE111672       | GSE111672           | Moncada et al. 2020; DOI 10.1038/s41587-019-0392-8     |               2 |               6 |           0 |               0 |         3119 | genes_by_spots       | processed ST counts, log2 transformed for scoring            | treatment-naive             | True               |

## Primary Hwang Naive Result

HALLMARK_PROTEIN_SECRETION Model B beta = 0.00351734, q = 0.324207. Reduced Model Level 2 was used because lymphoid fraction was unavailable in official segment metadata.

## Treatment Sensitivity

Hwang treated was analyzed separately in `phase9b3b_hwang_treated_models.tsv`.

## Moncada

Moncada was analyzed only as exploratory cross-platform spatial consistency across 6 section summaries from 2 patients. No formal population-level replication is claimed.

## Negative Controls

All locked negative-control classes were executed or deterministically audited in `phase9b3b_negative_control_results.tsv`; no placeholder rows are present.

## Evidence

Evidence categories are in `phase9b3b_spatial_evidence.tsv`.

Final readiness decision: READY_FOR_PHASE9B3C_INDEPENDENT_REVIEW
