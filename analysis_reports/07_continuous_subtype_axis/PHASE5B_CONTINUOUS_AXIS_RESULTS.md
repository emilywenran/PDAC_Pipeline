# Phase 5B Continuous Axis Results

## Execution status

All seven locked analysis IDs in `01_metadata/continuous_axis_parameter_inventory.tsv` were executed with seed 2026, 1000 bootstrap iterations, and 1000 permutation iterations. The primary axis used `Moffitt_50_gene_axis.tsv` with 25 Basal-like genes and 25 Classical genes including LEMD1. The LEMD1 sensitivity used `Moffitt_49_gene_axis_no_LEMD1.tsv` with 24 Basal-like genes and 25 Classical genes; the only gene-set difference was LEMD1.

Reference-anchored centroid results are descriptive because public Basal and Classical labels contributed to centroid definition. Leave-one-out centroids were calculated without including a sample in its own public-label centroid.

## Main findings

The locked overall decision category is **INCONCLUSIVE**.

Primary score evidence shows a median Hybrid coactivation score of 0.361, median Hybrid distance-to-both-poles of 2.077, and median Hybrid assignment entropy of 0.130. The dominant public Hybrid interpretation category was `TO_VERIFY`.

Continuous scoring systems were compared by Spearman correlation, rank concordance, direction agreement, and method-sensitive sample detection in `phase5b_score_method_concordance.tsv` and `phase5b_method_sensitive_samples.tsv`. Ordered Classical-to-Hybrid-to-Basal trends, effect sizes, bootstrap CIs, and permutation P values are reported in `phase5b_ordered_trend_tests.tsv`.

## Hybrid-state behavior

Public Hybrid samples were evaluated with basal/classical program scores, contrast, coactivation, centroid distances, PurIST probability, Moffitt score difference, Phase 4B item consensus, entropy, silhouette width, and method variance. The same metrics were also written for Basal and Classical samples in `phase5b_hybrid_state_assessment.tsv`.

## Stability integration and sensitivity

Associations between continuous axis position and Phase 4B stability metrics are reported in `phase5b_axis_stability_relationships.tsv`. Input-scale, outlier-exclusion, no-LEMD1, unsupervised-centroid, and leave-one-out sensitivity summaries are reported in `phase5b_sensitivity_summary.tsv` and `phase5b_category_transition_summary.tsv`.

## Downstream recommendation

For downstream microbiome analyses, use `AXIS_MOFFITT50_PRIMARY` basal-classical contrast as the primary continuous transcriptional-axis outcome, with `AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY`, `AXIS_SECONDARY`, and centroid-distance scores as prespecified sensitivity evidence. Do not use public subtype labels to optimize thresholds.

## TO_VERIFY

Samples flagged by locked outlier rules or high method sensitivity remain labelled `TO_VERIFY` or `METHOD_SENSITIVE` in the output tables. Hartigan's Dip Test was run with the Python diptest package.
