# Phase 9B1 Bulk External Validation Results

**SUPERSEDED_BY_PHASE9B1R.** This report is preserved only as an audit artifact. The PurIST, Hallmark, TF, module-replication, cross-cohort synthesis, and evidence-classification results in this file must not be cited as valid results. Use `PHASE9B1R_CORRECTED_BULK_EXTERNAL_VALIDATION_RESULTS.md` and `phase9b1r_*` outputs instead.

## Scope
Phase 9B1 executed only independent bulk-transcriptome validation. Single-cell, spatial, microbiome validation, target prioritization, causal mediation, post hoc signature modification, and manuscript writing were not performed.

## Cohorts Successfully Analyzed
| dataset_id   | accession   |   analyzed_samples |   genes_after_mapping |
|:-------------|:------------|-------------------:|----------------------:|
| TCGA_PAAD    | TCGA-PAAD   |                178 |                 59427 |
| GSE71729     | GSE71729    |                145 |                 19736 |
| GSE62452     | GSE62452    |                 69 |                 31886 |

## Excluded Cohorts and Reasons
No Phase 9A PRIORITY_1 bulk-host cohort was excluded by the script. Non-bulk PRIORITY_1 cohorts were intentionally not analyzed in Phase 9B1.

## Signature and Module Coverage
Coverage tables were written to `phase9b1_signature_coverage.tsv` and `phase9b1_module_transfer_coverage.tsv`. Coverage below feasibility is classified as `INSUFFICIENT_EXTERNAL_DATA`.

## Basal-Classical Score Reproducibility
Moffitt50 basal, classical, basal-classical contrast, Moffitt49 no-LEMD1 contrast, and PurIST probability were calculated within each cohort without optimizing subtype cutoffs against cohort labels.

## Externally Replicated Pathways
None under locked CI/external replication criteria.

## Externally Replicated TF Activities
TF activity replication is labeled `TO_VERIFY` where full decoupleR/VIPER was not executable and TF-symbol proxy scoring was used.

## Externally Replicated Module Signatures
MEblack, MEblue, MEgreen, MEred

## Purity and Composition Sensitivity
External purity/immune/stromal covariates were not uniformly available across all processed matrices in this execution. Sensitivity is therefore `TO_VERIFY` except where cohort metadata later provides validated transcriptome-derived estimates.

## Negative-Control Results
Negative controls were written to `phase9b1_negative_control_results.tsv`, including size-matched random modules, patient-label permutation, and unrelated pathway controls. Gene-label permutation is represented through randomized module signatures.

## Meta-Analysis Results
Cross-cohort synthesis was attempted only for features with three comparable cohort-specific effects. Results are in `phase9b1_cross_cohort_synthesis.tsv`.

## Null and Failed Replication Findings
Features classified as `NOT_REPLICATED` or `INSUFFICIENT_EXTERNAL_DATA` are listed in `phase9b1_host_feature_replication_evidence.tsv`.

## Phase 9B2 Readiness
Phase 9B2 single-cell validation may proceed after review of Phase 9B1 outputs; no Phase 9B1 result requires altering the locked Phase 9B2 plan.

## TO_VERIFY Items
FOXK2, TP63, TWIST1, ZNF24
