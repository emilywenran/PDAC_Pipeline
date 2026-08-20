# Phase 7A.5 Host TME Covariates

## Status

Phase 7A.5 is complete. No microbiome association tests, differential abundance, PERMANOVA, survival analysis, pathway enrichment, mediation analysis, or target prioritization were performed.

Phase 7B may proceed using Model 0 as primary, Model 1 as the existing technical sensitivity model, and the Phase 7A.5-permitted host TME sensitivity models Model 3P, Model 3I, and Model 3S.

## Method

Scores were calculated with the official MD Anderson ESTIMATE R package installed from the project-recommended R-Forge source repository.

| Component | Version / Source |
|---|---|
| R | R version 4.5.3 (2026-03-11) |
| ESTIMATE | `estimate` 1.0.13, R-Forge official source package |
| R executor | `06_scripts/R/10_phase7a5_host_covariates.R` |
| Python validator | `06_scripts/python/10_validate_phase7a5_host_covariates.py` |

Input expression matrix: `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz`.

Input expression scale: log2 analysis-ready values from prior expression processing. The expression matrix itself was not altered.

Exact sample mapping: the 62 expression matrix columns matched `01_metadata/expression_sample_crosswalk.tsv` in order exactly. Output rows in `01_metadata/host_tme_covariates.tsv` preserve that order and map each `expression_sample_id` to one `patient_id`.

ESTIMATE execution used `filterCommonGenes(id = "GeneSymbol")` followed by `estimateScore(platform = "affymetrix")`. The `platform = "affymetrix"` setting was used so the official package would emit the ESTIMATE cosine-equation tumor-purity estimate. This purity value is transcriptome-derived and should not be interpreted as an independent pathology measurement.

## Gene-Set Coverage

ESTIMATE common-gene coverage was 9,837 / 10,412 genes. The ESTIMATE stromal/immune signature gene coverage was 277 / 282 genes. Missing SI genes were `FYB`, `GPR124`, `LPPR4`, `TXNDC3`, and `ODZ4`.

ESTIMATE package messages recorded:

| Check | Result |
|---|---|
| Common gene merge | 9,837 genes included; 575 mismatched |
| Stromal signature overlap | 137 genes |
| Immune signature overlap | 140 genes |
| Warnings | None |

## Score Completeness

All 62 patients have finite stromal, immune, ESTIMATE, and inferred tumor-purity scores. There were no duplicated patients and no missing or infinite scores. Inferred purity values were bounded from 0.317 to 0.948.

| Score | Mean | SD | Median | Min | Max |
|---|---:|---:|---:|---:|---:|
| Stromal score | 773.912 | 785.010 | 942.815 | -1237.248 | 1877.295 |
| Immune score | 727.974 | 826.222 | 528.584 | -707.631 | 2906.211 |
| ESTIMATE score | 1501.886 | 1434.666 | 1537.403 | -1905.521 | 4382.869 |
| Inferred tumor purity | 0.663 | 0.149 | 0.674 | 0.317 | 0.948 |

No technically extreme values were flagged by the Phase 7A.5 descriptive IQR/z-score screen.

## Distributions And Missingness

Distribution plots were written to `05_results/figures/phase7a5_host_covariate_distributions.pdf`. Missingness details were written to `05_results/tables/phase7a5_missingness_report.tsv`; all required score fields were complete.

## Collinearity Assessment

Pairwise Spearman correlations were descriptive host-covariate QC only; microbiome genera were not inspected.

| Pair | Spearman rho | Warning |
|---|---:|---|
| Stromal vs immune | 0.598 | none |
| Stromal vs ESTIMATE | 0.860 | HIGH_CORRELATION |
| Stromal vs inferred purity | -0.860 | HIGH_CORRELATION |
| Immune vs ESTIMATE | 0.891 | HIGH_CORRELATION |
| Immune vs inferred purity | -0.891 | HIGH_CORRELATION |
| ESTIMATE vs inferred purity | -1.000 | SEVERE_COLLINEARITY |

The combined all-TME screen (`host score + inferred purity + immune score + stromal score + ESTIMATE score`) is not permitted: maximum VIF was infinite and condition number was 1.04e15.

## Permitted Sensitivity Models

Model 0 remains the primary microbiome association model. The following host TME models are permitted only as sensitivity analyses:

| Model | Covariates | Available patients | Maximum VIF | Condition number | Role |
|---|---|---:|---:|---:|---|
| Model 3P | host transcriptional score + inferred tumor purity | 62 | 1.40 | 1.82 | sensitivity |
| Model 3I | host transcriptional score + immune score | 62 | 1.24 | 1.61 | sensitivity |
| Model 3S | host transcriptional score + stromal score | 62 | 1.46 | 1.88 | sensitivity |

The ESTIMATE score itself is not added as a separate sensitivity model because it is the arithmetic sum of stromal and immune scores and is perfectly monotonic with inferred tumor purity under the ESTIMATE purity equation.

## Interpretation Guardrails

1. ESTIMATE-derived purity and cell-composition scores are inferred from the same host transcriptomic matrix as the Moffitt scores.
2. These variables are sensitivity covariates, not independent experimental measurements.
3. Adjustment may remove biological variation genuinely associated with the PDAC transcriptional state.
4. Therefore, Model 0 remains primary and Models 3P/3I/3S assess robustness only.

## Outputs

| Output | Purpose |
|---|---|
| `01_metadata/host_tme_covariates.tsv` | Patient-level ESTIMATE-derived covariates |
| `05_results/tables/phase7a5_host_covariate_qc.tsv` | Descriptive score QC and extreme-value screen |
| `05_results/tables/phase7a5_host_covariate_correlations.tsv` | Pairwise host-covariate Spearman correlations |
| `05_results/tables/phase7a5_covariate_model_feasibility.tsv` | VIF, condition-number, complete-case, and model-permission decisions |
| `05_results/figures/phase7a5_host_covariate_distributions.pdf` | Score distributions |
| `05_results/figures/phase7a5_host_covariate_correlation.pdf` | Host-covariate correlation heatmap |

## Validation

Both scripts ran successfully:

```text
Rscript 06_scripts/R/10_phase7a5_host_covariates.R
python3 06_scripts/python/10_validate_phase7a5_host_covariates.py
```

The validator reported `VALIDATION_PASS`.

## TO_VERIFY

TO_VERIFY: ESTIMATE inferred purity was generated using the package's Affymetrix tumor-purity equation on log2 analysis-ready RNA expression values. This is the official package-supported purity calculation, but it remains an inferred transcriptomic estimate rather than a pathology-derived tumor cellularity measurement.

TO_VERIFY: The ESTIMATE package is no longer available from the current Bioconductor release for R 4.5.3; the official MD Anderson R-Forge installation path was used to install `estimate` 1.0.13.
