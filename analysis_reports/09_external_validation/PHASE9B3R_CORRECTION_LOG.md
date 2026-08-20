# Phase 9B3R Correction Log

The original Phase 9B3B report and Phase 9B3C FAIL review are preserved as audit history. Phase 9B3B is superseded by Phase 9B3R.

## Finding Closure

- **FIND-01**: CLOSED. Real coordinate, random-gene, unrelated-Hallmark, label-permutation, and leakage null distributions executed and saved.
- **FIND-02**: CLOSED. HALLMARK_SPERMATOGENESIS and ineligible WGCNA modules are gated as INSUFFICIENT_SPATIAL_DATA before modeling.
- **FIND-03**: CLOSED. Nonconverged treated Model C audit row retained with inferential fields set to NA and no q value.
- **FIND-04**: CLOSED. Figures include only eligible converged model rows with finite coefficients.
- **FIND-05**: UNCHANGED_RESOLVED. Provenance/count discrepancies remain documented and do not require code correction.

## Eligibility Gate

| dataset_id            | feature_layer   | feature_name                  |   genes_expected |   genes_available |   coverage_fraction | formal_inference_status   | scoring_method                                                     |
|:----------------------|:----------------|:------------------------------|-----------------:|------------------:|--------------------:|:--------------------------|:-------------------------------------------------------------------|
| HWANG_GSE202051_NAIVE | Hallmark        | HALLMARK_PROTEIN_SECRETION    |               96 |                79 |           0.822917  | ELIGIBLE                  | rank-normalized ssGSEA-style enrichment gated before model fitting |
| HWANG_GSE202051_NAIVE | Hallmark        | HALLMARK_SPERMATOGENESIS      |              135 |                50 |           0.37037   | INSUFFICIENT_SPATIAL_DATA | rank-normalized ssGSEA-style enrichment gated before model fitting |
| HWANG_GSE202051_NAIVE | Hallmark        | HALLMARK_BILE_ACID_METABOLISM |              112 |                59 |           0.526786  | INSUFFICIENT_SPATIAL_DATA | rank-normalized ssGSEA-style enrichment gated before model fitting |
| HWANG_GSE202051_NAIVE | Hallmark        | HALLMARK_HEME_METABOLISM      |              200 |               127 |           0.635     | INSUFFICIENT_SPATIAL_DATA | rank-normalized ssGSEA-style enrichment gated before model fitting |
| HWANG_GSE202051_NAIVE | TF_regulon      | ELF1                          |                0 |                 0 |           0         | INSUFFICIENT_SPATIAL_DATA | regulon unavailable locally; TF-symbol proxy prohibited            |
| HWANG_GSE202051_NAIVE | TF_regulon      | MBD2                          |                0 |                 0 |           0         | INSUFFICIENT_SPATIAL_DATA | regulon unavailable locally; TF-symbol proxy prohibited            |
| HWANG_GSE202051_NAIVE | TF_regulon      | ZBTB7A                        |                0 |                 0 |           0         | INSUFFICIENT_SPATIAL_DATA | regulon unavailable locally; TF-symbol proxy prohibited            |
| HWANG_GSE202051_NAIVE | TF_regulon      | ZNF384                        |                0 |                 0 |           0         | INSUFFICIENT_SPATIAL_DATA | regulon unavailable locally; TF-symbol proxy prohibited            |
| HWANG_GSE202051_NAIVE | TF_regulon      | ZNF740                        |                0 |                 0 |           0         | INSUFFICIENT_SPATIAL_DATA | regulon unavailable locally; TF-symbol proxy prohibited            |
| HWANG_GSE202051_NAIVE | WGCNA_module    | MEblack                       |              270 |                12 |           0.0444444 | INSUFFICIENT_SPATIAL_DATA | standardized mean rank gated by coverage                           |
| HWANG_GSE202051_NAIVE | WGCNA_module    | MEblue                        |             1691 |               174 |           0.102898  | INSUFFICIENT_SPATIAL_DATA | standardized mean rank gated by coverage                           |
| HWANG_GSE202051_NAIVE | WGCNA_module    | MEgreen                       |              515 |                77 |           0.149515  | INSUFFICIENT_SPATIAL_DATA | standardized mean rank gated by coverage                           |
| HWANG_GSE202051_NAIVE | WGCNA_module    | MEtan                         |              103 |                11 |           0.106796  | INSUFFICIENT_SPATIAL_DATA | standardized mean rank gated by coverage                           |
| HWANG_GSE202051_NAIVE | WGCNA_module    | MEgreenyellow                 |              166 |                11 |           0.0662651 | INSUFFICIENT_SPATIAL_DATA | standardized mean rank gated by coverage                           |

## Model Convergence Audit

| cohort_id               | feature_name               | model_id                                 | model_converged    |       p_value |       q_value | eligibility_status        |
|:------------------------|:---------------------------|:-----------------------------------------|:-------------------|--------------:|--------------:|:--------------------------|
| HWANG_GSE202051_NAIVE   | HALLMARK_PROTEIN_SECRETION | HWANG_GSE202051_NAIVE_MODEL_A            | True               |   9.81706e-53 |   2.94512e-52 | ELIGIBLE                  |
| HWANG_GSE202051_NAIVE   | HALLMARK_PROTEIN_SECRETION | HWANG_GSE202051_NAIVE_MODEL_B            | True               |   0.270172    |   0.405258    | ELIGIBLE                  |
| HWANG_GSE202051_NAIVE   | HALLMARK_PROTEIN_SECRETION | HWANG_GSE202051_NAIVE_MODEL_C            | True               |   0.462129    |   0.462129    | ELIGIBLE                  |
| HWANG_GSE202051_NAIVE   | HALLMARK_SPERMATOGENESIS   | HWANG_GSE202051_NAIVE_ELIGIBILITY_GATE   | NOT_FIT_INELIGIBLE | nan           | nan           | INSUFFICIENT_SPATIAL_DATA |
| HWANG_GSE202051_TREATED | HALLMARK_PROTEIN_SECRETION | HWANG_GSE202051_TREATED_MODEL_A          | True               |   1.27615e-21 |   2.5523e-21  | ELIGIBLE                  |
| HWANG_GSE202051_TREATED | HALLMARK_PROTEIN_SECRETION | HWANG_GSE202051_TREATED_MODEL_B          | True               |   0.781822    |   0.781822    | ELIGIBLE                  |
| HWANG_GSE202051_TREATED | HALLMARK_PROTEIN_SECRETION | HWANG_GSE202051_TREATED_MODEL_C          | False              | nan           | nan           | ELIGIBLE                  |
| HWANG_GSE202051_TREATED | HALLMARK_SPERMATOGENESIS   | HWANG_GSE202051_TREATED_ELIGIBILITY_GATE | NOT_FIT_INELIGIBLE | nan           | nan           | INSUFFICIENT_SPATIAL_DATA |

## Negative-Control Execution

| dataset_id              | control_type                       | control_id                    |   observed_statistic |   iterations |   seed |    null_mean |   null_variance |   empirical_p_value | significant   | execution_status   | notes                                                                                            |
|:------------------------|:-----------------------------------|:------------------------------|---------------------:|-------------:|-------:|-------------:|----------------:|--------------------:|:--------------|:-------------------|:-------------------------------------------------------------------------------------------------|
| HWANG_GSE202051_NAIVE   | coordinate permutation             | HALLMARK_PROTEIN_SECRETION    |          -0.0090034  |         1000 |   2026 | -0.0124652   |     4.47526e-06 |         0.952048    | False         | EXECUTED           | Moffitt50 contrast permuted within Hwang section before tumor-axis coefficient calculation       |
| HWANG_GSE202051_NAIVE   | label permutation                  | HALLMARK_PROTEIN_SECRETION    |           0.0475832  |         1000 |   2026 |  0.000516415 |     1.72815e-05 |         0.000999001 | True          | EXECUTED           | tumor/stroma labels permuted within patient                                                      |
| HWANG_GSE202051_NAIVE   | size-matched random gene set       | HALLMARK_PROTEIN_SECRETION    |          -0.0090034  |          100 |   2026 | -0.000122937 |     2.37718e-05 |         0.0693069   | False         | EXECUTED           | random sets matched to primary target gene count                                                 |
| HWANG_GSE202051_NAIVE   | expression-matched random gene set | HALLMARK_PROTEIN_SECRETION    |          -0.0090034  |          100 |   2026 |  0.00212594  |     2.67156e-05 |         0.128713    | False         | EXECUTED           | random sets matched by median-expression decile                                                  |
| HWANG_GSE202051_NAIVE   | unrelated Hallmark pathway         | HALLMARK_BILE_ACID_METABOLISM |          -0.00213313 |          100 |   2026 |  0.00122369  |     8.31856e-06 |         0.455446    | False         | EXECUTED           | available genes from unrelated Hallmark scored only as a negative-control feature                |
| HWANG_GSE202051_NAIVE   | unrelated Hallmark pathway         | HALLMARK_HEME_METABOLISM      |          -0.00258295 |          100 |   2026 | -0.00443748  |     3.16668e-06 |         0.831683    | False         | EXECUTED           | available genes from unrelated Hallmark scored only as a negative-control feature                |
| HWANG_GSE202051_NAIVE   | leakage control                    | HALLMARK_PROTEIN_SECRETION    |           0.0126582  |          100 |   2026 |  0.00189873  |     2.06358e-05 |         0.158416    | False         | EXECUTED           | target-gene overlap with Moffitt axis and morphology-label tokens compared with random gene sets |
| HWANG_GSE202051_TREATED | coordinate permutation             | HALLMARK_PROTEIN_SECRETION    |          -0.01537    |         1000 |   2026 | -0.0144499   |     2.55686e-05 |         0.437562    | False         | EXECUTED           | Moffitt50 contrast permuted within Hwang section before tumor-axis coefficient calculation       |
| HWANG_GSE202051_TREATED | label permutation                  | HALLMARK_PROTEIN_SECRETION    |           0.0480622  |         1000 |   2026 | -0.000467237 |     3.42842e-05 |         0.000999001 | True          | EXECUTED           | tumor/stroma labels permuted within patient                                                      |
| HWANG_GSE202051_TREATED | size-matched random gene set       | HALLMARK_PROTEIN_SECRETION    |          -0.01537    |          100 |   2026 |  0.000676699 |     3.46283e-05 |         0.00990099  | True          | EXECUTED           | random sets matched to primary target gene count                                                 |
| HWANG_GSE202051_TREATED | expression-matched random gene set | HALLMARK_PROTEIN_SECRETION    |          -0.01537    |          100 |   2026 |  0.00300742  |     2.9712e-05  |         0.019802    | True          | EXECUTED           | random sets matched by median-expression decile                                                  |
| HWANG_GSE202051_TREATED | unrelated Hallmark pathway         | HALLMARK_BILE_ACID_METABOLISM |          -0.00340019 |          100 |   2026 | -0.00275237  |     1.86374e-05 |         0.50495     | False         | EXECUTED           | available genes from unrelated Hallmark scored only as a negative-control feature                |
| HWANG_GSE202051_TREATED | unrelated Hallmark pathway         | HALLMARK_HEME_METABOLISM      |           0.00633776 |          100 |   2026 |  0.00411914  |     7.49661e-06 |         0.237624    | False         | EXECUTED           | available genes from unrelated Hallmark scored only as a negative-control feature                |
| HWANG_GSE202051_TREATED | leakage control                    | HALLMARK_PROTEIN_SECRETION    |           0.0126582  |          100 |   2026 |  0.00316456  |     3.35837e-05 |         0.247525    | False         | EXECUTED           | target-gene overlap with Moffitt axis and morphology-label tokens compared with random gene sets |

Final evidence category: `PARTIAL_SPATIAL_SUPPORT`.

Final readiness decision: READY_FOR_PHASE9B3C2_COMPLETE_INDEPENDENT_REVIEW
