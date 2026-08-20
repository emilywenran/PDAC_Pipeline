# Phase 2B Analysis-Ready Expression: GSE172356

## Scope and Guardrails

K-Dense `exploratory-data-analysis` was used from `~/.agents/skills/exploratory-data-analysis/SKILL.md`. Phase 2B audited missingness, selected missing-value handling, evaluated subtype-independent expression filters, transformed the selected matrix, and reassessed sample QC. No subtype classification, differential expression, supervised feature selection, batch correction, sample exclusion, or biological interpretation was performed.

## Source and Missing-Value Pattern

Input audited matrix: `03_processed/expression/GSE172356_expression_audited.tsv.gz`, with 45140 genes and 62 mapped samples. The original GEO processed matrix contains 73202 missing expression cells. Source representation audit found 73202 literal `NA` strings, 0 blank fields, and 0 parse failures. Therefore, the missing cells are source literal `NA` values, not parsing failures introduced in Phase 2A.

Missingness is concentrated in genes, not randomly scattered: 0 missing in 42758 genes, 21 missing in 1223 genes, 41 missing in 1159 genes. All 62 samples have missing values; per-sample missing counts range from 1159 to 1223. The pattern forms two complementary sample blocks. Missing count correlation with total expression is 0.1908; correlation with detected genes is -0.2357. Batch fields in `sample_manifest.tsv` are `NA`, so batch association is not assessable from current metadata.

Subtype association was audited descriptively only, not used for filtering or threshold choice:

```text
subtype_original  count  min  median  max
           Basal     17 1159  1223.0 1223
       Classical     22 1159  1159.0 1223
          Hybrid     23 1159  1159.0 1223
```

Batch summary:

```text
 batch  count  min  median  max
   NaN     62 1159  1159.0 1223
```

Official source documentation confirms the processed matrix is based on HTSeq counts normalized by DESeq size factors, but the available GEO/source documentation does not explain the literal `NA` entries. Supplementary Data 4 `Figure1.GeneMatrix` contains a separate 94-signature table and does not document the 45,140-gene processed-matrix `NA` values. TO_VERIFY: confirm with source authors or raw-count reprocessing if the literal `NA` provenance becomes material.

## Missing-Value Handling Decision

Primary strategy: complete-observation filtering, retaining genes with complete observations before expression filtering. This avoids replacing unexplained source `NA` values. The rejected alternatives were:

- Remove genes exceeding a missingness threshold: thresholds <=20% collapse to complete observations because partial-missing genes are missing in 21 or 41 samples; <=50% would require imputing genes missing in 21 samples.
- Replace structurally absent count-like entries with zero: rejected for the primary matrix because source documentation does not state that `NA` means no reads or structural absence.
- Gene-median imputation: evaluated as a sensitivity strategy only; rejected for the primary matrix because it inserts modeled values into block-missing source cells.

Sensitivity analyses required later: repeat Phase 3 subtype reproduction using complete-gene primary matrix, a <=50% missingness plus gene-median imputed matrix, and an all-`NA` zero-filled matrix only as a stress test clearly labelled unsupported by source semantics.

## Filtering and Transformation

Filtering was evaluated without subtype labels. The selected primary rule is `primary_complete_not_all_zero_count_ge_1_10_percent`: Complete observations, not all-zero, and normalized count >= 1 in at least 10% of samples. This retained 42654 genes and all 62 samples. Full sensitivity table: `05_results/tables/phase2b_filtering_sensitivity.tsv`.

No additional library-size or size-factor normalization was performed because the source matrix already contains DESeq size-factor-normalized counts. The analysis transform is `log2(normalized_count + 1)`, applied after filtering. DESeq2 VST or rlog were not automatically applied because those methods are designed around raw count inputs and size-factor/model estimation; this matrix is already normalized and contains fractional values, and raw integer counts are not available in the Phase 2B inputs.

Preserved outputs:

- Filtered normalized counts: `03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz`.
- Filtered log2 analysis-ready matrix: `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz`.

Final matrix dimensions: 42654 genes x 62 samples. Final missing-value count: 0.

## Four Suspected Outliers

All 62 samples are retained. The four Phase 2A suspected samples are classified as `RETAIN_WITH_SENSITIVITY_ANALYSIS`; no sample is recommended for exclusion based on a single metric.

```text
expression_column patient_id    phase2b_sample_classification                                                      objective_evidence
         YX16135T   PDAC_016 RETAIN_WITH_SENSITIVITY_ANALYSIS total_expression_robust_z_abs_gt_3.5;robust_mahalanobis_robust_z_gt_3.5
         YX16158T   PDAC_023 RETAIN_WITH_SENSITIVITY_ANALYSIS        median_correlation_robust_z_lt_-3.5;pca_distance_robust_z_gt_3.5
         YX16194T   PDAC_033 RETAIN_WITH_SENSITIVITY_ANALYSIS                                    total_expression_robust_z_abs_gt_3.5
         YX16224T   PDAC_039 RETAIN_WITH_SENSITIVITY_ANALYSIS                                    total_expression_robust_z_abs_gt_3.5
```

Outlier assessment table: `05_results/tables/phase2b_outlier_assessment.tsv`. Figures: `05_results/figures/phase2b_transformed_pca.pdf`, `05_results/figures/phase2b_transformed_sample_correlation.pdf`, and `05_results/figures/phase2b_sample_qc_summary.pdf`. PCA on the filtered log2 matrix explains 0.1454 on PC1 and 0.0715 on PC2.

## Validation and Phase 3 Readiness

Validation requirements are satisfied: exactly 62 mapped samples remain, sample order matches `expression_sample_crosswalk.tsv`, no duplicated samples or genes are present, no infinite values are present, final missing-value count is 0, the transform is reproducible as `log2(normalized_count + 1)`, and the original audited matrix remains unchanged.

Phase 3 subtype reproduction may proceed from `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz` after human review. Required sensitivity analyses are complete-gene primary vs imputation stress-test matrices, and inclusion vs exclusion of each retained-with-sensitivity Phase 2A suspected outlier.
