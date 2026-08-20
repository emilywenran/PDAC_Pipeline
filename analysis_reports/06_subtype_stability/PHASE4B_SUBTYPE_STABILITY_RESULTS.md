# Phase 4B Subtype Stability Results

## Execution status

All eight locked analyses were executed for K=2 through K=6 with seed 2026 and 1,000 resampling iterations per analysis: STAB_CSY_PRIMARY, STAB_CSY_LOG2, STAB_UNSUP_HVG, STAB_CSY_OUTLIER_EXCL, STAB_HVG_OUTLIER_EXCL, STAB_CSY_FEAT_RESAMP, STAB_CSY_IMPUTED, STAB_HVG_VAR_FILTER.

Public subtype labels were not used during clustering or K selection. They were used only after clustering for descriptive Hungarian/max-overlap label alignment.

## Evidence by candidate K

| Analysis | Preferred K | K=2 PAC | K=3 PAC | K=3 silhouette | K=3 mean Jaccard | K=3 min Jaccard |
|---|---:|---:|---:|---:|---:|---:|
| STAB_CSY_PRIMARY | 2 | 0.061 | 0.196 | 0.236 | 0.559 | 0.530 |
| STAB_CSY_LOG2 | 2 | 0.178 | 0.256 | 0.276 | 0.505 | 0.373 |
| STAB_UNSUP_HVG | 4 | 0.597 | 0.347 | 0.164 | 0.517 | 0.484 |
| STAB_CSY_OUTLIER_EXCL | 2 | 0.036 | 0.174 | 0.245 | 0.574 | 0.539 |
| STAB_HVG_OUTLIER_EXCL | 3 | 0.574 | 0.236 | 0.155 | 0.542 | 0.526 |
| STAB_CSY_FEAT_RESAMP | 2 | 0.108 | 0.327 | 0.236 | 0.563 | 0.532 |
| STAB_CSY_IMPUTED | 2 | 0.058 | 0.189 | 0.236 | 0.563 | 0.535 |
| STAB_HVG_VAR_FILTER | 2 | 0.246 | 0.281 | 0.167 | 0.549 | 0.532 |

## Primary and independent preferred K

The primary CSY normalized-count analysis preferred K=2. Its K=2 PAC was 0.061, while K=3 PAC was 0.196 with mean silhouette 0.236 and mean bootstrap Jaccard 0.559.

The independent HVG analysis preferred K=4. Its K=3 PAC was 0.347, mean silhouette was 0.164, and mean bootstrap Jaccard was 0.517. Primary-vs-HVG K=3 ARI was 0.332.

## Basal-like, classical, and hybrid stability

| Public subtype | Primary K=3 mean entropy | Primary K=3 item consensus | Primary K=3 silhouette |
|---|---:|---:|---:|
| Basal | 0.005 | 0.999 | 0.374 |
| Classical | 0.239 | 0.827 | 0.201 |
| Hybrid | 0.258 | 0.894 | 0.169 |

Hybrid samples are not supported as a universally stable discrete group across all representations. In the primary K=3 run they show moderate item consensus and positive silhouette, but the log2 and HVG sensitivity results show method dependence. The appropriate locked interpretation is mixed METHOD_SENSITIVE / intermediate-or-heterogeneous evidence rather than a final biological interpretation from one metric.

Interpretation category counts across all K=3 analyses:

| Public subtype | Category | Count |
|---|---|---:|
| Basal | HETEROGENEOUS_OR_UNSTABLE | 17 |
| Basal | STABLE_BASAL | 112 |
| Basal | TO_VERIFY | 7 |
| Classical | HETEROGENEOUS_OR_UNSTABLE | 40 |
| Classical | STABLE_CLASSICAL | 97 |
| Classical | TO_VERIFY | 39 |
| Hybrid | HETEROGENEOUS_OR_UNSTABLE | 32 |
| Hybrid | INTERMEDIATE_STATE | 7 |
| Hybrid | STABLE_HYBRID | 122 |
| Hybrid | TO_VERIFY | 15 |

## Sensitivity evidence

Log2 transformation changed 30 K=3 assignments relative to the primary run (ARI 0.385), increased K=3 PAC by 0.060, and did not change the preferred K from 2 to 2.

Excluding the four Phase 2B outlier candidates changed 3 K=3 primary assignments among common samples (ARI 0.843) and kept preferred K=2.

Alternative missing-value handling changed 0 K=3 assignments (ARI 1.000). Feature resampling changed 0 K=3 assignments but increased K=3 PAC by 0.131. The HVG variance-filter sensitivity changed 8 K=3 assignments (ARI 0.626) and changed the independent preferred K from 4 to 2.

## Recurrently unstable samples

20 samples were recurrently unstable in at least three of eight K=3 analyses: YX16112T, YX16147T, YX16188T, YX16128T, YX16155T, YX16222T, YX15211T, YX16041T, YX16057T, YX16158T, YX16218T, YX16227T, YX16236T, YX15047T, YX16070T, YX16094T, YX16113T, YX16172T, YX16202T, YX16224T.

## Public-label agreement

The primary K=3 clusters aligned exactly to public labels post hoc (ARI 1.000, NMI 1.000, Cohen's kappa 1.000). This is descriptive agreement only and was not used for K selection.

## Decision-rule application

Primary analysis conclusion: K=2 is preferred by the locked multi-metric rank; K=3 has moderate PAC and sub-threshold Jaccard stability.

Independent HVG analysis: K=4 is preferred, and K=3 is not clearly stable by PAC/Jaccard.

Sensitivity-analysis evidence: preprocessing scale and HVG filtering materially affect assignments or preferred K, while imputation and feature-resampling assignments are less disruptive but still alter stability metrics.

Public-label agreement: strong for primary K=3 post hoc, but this is not stability evidence for K selection.

Prespecified overall decision category: `INCONCLUSIVE`.

Proceed to continuous basal-classical axis: Yes. The locked stability results do not provide strong support for three robust discrete clusters, and scale/HVG sensitivity remains material.

## TO_VERIFY

- ConsensusClusterPlus and fpc were not installed in the local R environment, so the locked resampling, consensus, bootstrap Jaccard, and prediction-strength procedures were implemented directly in R using the locked parameters.
- The Phase 4A protocol did not define a formal iteration count for prediction strength separately from the resampling count; this execution used the locked 1,000 iteration count.
- Public-label agreement remains exact for the primary K=3 descriptive comparison, but stability metrics do not independently establish a robust three-cluster structure.

## Runtime and versions

Total analysis runtime across per-analysis records: 59.5 seconds.

| Component | Version |
|---|---|
| R | R version 4.5.3 (2026-03-11) |
| cluster | 2.1.8.2 |
| clue | 0.3.68 |
| mclust | 6.1.2 |
| ggplot2 | 4.0.3 |
