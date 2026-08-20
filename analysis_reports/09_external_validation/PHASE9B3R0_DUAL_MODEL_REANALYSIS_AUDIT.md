# Phase 9B3R0: Dual-Model Pre-Reanalysis Audit

This document summarizes the Phase 9B3R0 dual-model audit of the failed Phase 9B3B spatial validation implementation. Two independent reviewers (ChatGPT for implementation, Claude for statistics) evaluated the failed pipeline against the locked Phase 9B3A specifications.

## 1. Audit Methodology

- **ChatGPT Implementation Audit**: Scanned `16_phase9b3b_spatial_validation.py` and `16_validate_phase9b3b_spatial.py` for deviations from intended logic, focusing on code construction, mocked outputs, and validator insufficiencies.
- **Claude Statistical Audit**: Evaluated the statistical properties of the Phase 9B3B outputs, assessing hierarchy correctness, degrees of freedom assumptions, coverage rules, model convergence, and evidence logic.

## 2. Reconciled Findings

The following critical deviations from the locked spatial validation protocol were identified and must be repaired prior to Phase 9B3R reanalysis.

### 2.1 Faked Negative Controls (FIND_9B3_01)
The primary execution script hardcoded `observed_statistic = 0.0` and `empirical_p_value = 1.0` for all negative controls instead of running the mandated permutation loops. The spatial false-positive safeguards were completely bypassed, yet falsely logged as `"EXECUTED"`. The downstream validator logic was insufficient, as it accepted the "EXECUTED" string without verifying numerical variance.

### 2.2 Bypassed Coverage Gates (FIND_9B3_02)
The locked protocol dictates an 80% spatial coverage threshold for feature eligibility. While `coverage_table()` accurately flagged features with low coverage (e.g., `HALLMARK_SPERMATOGENESIS` at 37%), the scoring functions `rank_score()` and `mean_z()` ignored this flag, computing scores as long as at least one gene was present. These invalid scores then entered formal LMMs, polluting the evidence tables.

### 2.3 Silent Retention of Non-Converged Models (FIND_9B3_03)
When `statsmodels.mixedlm` fails to converge via L-BFGS, it retains the last parameter states. The Phase 9B3B script blindly extracted `beta` and `p_value` from `fit.params` without checking `fit.converged`. This resulted in the mathematically invalid Hwang-treated Model C contributing a $q$-value to reporting and coefficients to primary figures.

### 2.4 Small-Sample LMM Inference (FIND_9B3_04)
The locked software (`statsmodels`) defaults to asymptotic Z-tests for fixed effects in mixed models. This implicitly assumes infinite denominator degrees of freedom, significantly inflating Type I error rates for small sample sizes ($n=13$ naive, $n=7$ treated). While this is a statistical flaw (manifesting as extreme $q$-values like $10^{-52}$ for basic compartment differences), it accurately reflects the *locked plan*. The repair specification does not prescribe switching the inference engine post-hoc, but mandates explicit documentation of the anti-conservative Z-test properties.

### 2.5 Hardcoded Evidence Logic (FIND_9B3_05)
The synthesis logic in `16_phase9b3b_spatial_validation.py` hardcoded the `COMPARATOR` category as `NOT_SUPPORTED_SPATIALLY` and `PRIMARY` as `PARTIAL_SPATIAL_SUPPORT` based on partial assumptions, failing to dynamically evaluate the actual fitted betas and $q$-values.

## 3. Recommended Action

The dual-model audit concludes that Phase 9B3B is **computationally invalid**.

A detailed code repair specification has been generated at:
`05_results/tables/phase9b3r0_reconciled_repair_specification.tsv`

**Decision**: READY_FOR_PHASE9B3R_REANALYSIS
The project must proceed to Phase 9B3R, beginning with the modification of `16_phase9b3b_spatial_validation.py` and `16_validate_phase9b3b_spatial.py` strictly according to the repair specification TSV. No new biological data or unapproved statistical methodologies (e.g., ad-hoc KR degrees of freedom) may be introduced.
