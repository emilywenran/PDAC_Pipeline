# Phase 9B1R Corrected Bulk External Validation Results

## Scope
Phase 9B1R reran only independent bulk-transcriptome validation for TCGA_PAAD, GSE71729, and GSE62452. Single-cell validation was not performed.

## Errors Corrected
- FIND_01: Omission of PurIST model intercept (-6.815) in logistic link calculation.
- FIND_02: Violation of missing-gene policy for PurIST in GSE62452.
- FIND_03: Violation of missing-gene policy (80% coverage threshold) for WGCNA modules.
- FIND_04: Use of 15-gene proxy set and rank-mean instead of ssGSEA via decoupleR on the full pathway.
- FIND_05: Proxy single-gene expression used instead of VIPER TF activity, and inappropriate evidence classification.
- FIND_06: Incomplete negative control audit.

## Corrected PurIST
PurIST was recalculated with all available locked gene pairs, intercept beta0 = -6.815, logistic transformation, no cohort-specific refitting, and the locked 0.5 cutoff. Runtime validation: TCGA_PAAD=PASS; GSE71729=PASS; GSE62452=PASS.

## Corrected Hallmark Results
Hallmark scores were recalculated with MSigDB Hallmark 2026.1.Hs and decoupleR ssGSEA on full available pathway gene sets. Previous proxy scores are invalidated.

## Corrected TF Activity Results
DoRothEA A/B/C regulon coverage was evaluated per cohort and VIPER activity scoring was executed using decoupleR. TF evidence categories were derived from the saved cohort replication statistics and matched the Phase 9B1C2 audit counts (12 externally replicated, 13 partially replicated, 9 not replicated, 0 TO_VERIFY). No TF-symbol proxy is used.

## Module Coverage and Replication
The locked 80% external coverage threshold was enforced. Low-coverage cohort-module combinations are excluded from formal replication rather than counted as biological failures.

## Negative Controls
Patient-label permutation, gene-label permutation, size-matched randomized modules, expression-matched randomized modules, and unrelated Hallmark controls were executed where the corresponding feature was technically eligible.

## Cross-Cohort Synthesis and Evidence
Random-effects synthesis was used only where at least three eligible and comparable cohorts existed. Modules with only TCGA_PAAD eligibility are reported as cohort-specific or partial evidence.

## Phase 9B2 Readiness
Phase 9B2 may proceed only after Phase 9B1R validator and manifest validator pass. FIND_05 is now fully corrected and no TF activity remains TO_VERIFY.
