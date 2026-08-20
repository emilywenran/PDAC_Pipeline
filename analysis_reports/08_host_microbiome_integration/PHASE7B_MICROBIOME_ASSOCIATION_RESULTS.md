# Phase 7B Microbiome Association Results

Phase 7B executed the locked continuous association framework between tumor microbiome composition and PDAC transcriptional states. The primary host outcome was the Moffitt50 basal-classical contrast, where higher values indicate the locked Basal direction.

## Primary Global Community Result
The primary Aitchison PERMANOVA used 9,999 permutations with seed 2026 and found R-squared = 0.0534, pseudo-F = 3.3842, P = 0.0001.

## Primary Genus-Level Results
Exactly 122 primary genus tests were run with OLS and HC3 robust standard errors. 33 genera met the locked primary BH FDR threshold q < 0.05.
CLR coefficients are relative compositional associations and are not absolute microbial-load effects.

| Genus | Coefficient | 95% CI | P | q |
|---|---:|---:|---:|---:|
| Sphingobium | -0.4028 | [-0.5468, -0.2589] | 5.721e-07 | 6.979e-05 |
| Erythrobacter | -0.2549 | [-0.3517, -0.158] | 1.991e-06 | 0.0001215 |
| Novosphingobium | -0.2578 | [-0.3637, -0.1519] | 8.508e-06 | 0.0002776 |
| Pandoraea | 0.4328 | [0.2544, 0.6113] | 9.101e-06 | 0.0002776 |
| Candida | -1.124 | [-1.611, -0.6374] | 2.095e-05 | 0.0005112 |

## Supporting-Method Concordance
Spearman, permutation, and bootstrap outputs were generated for all primary genera. MaAsLin2 was not installed locally, so the required MaAsLin2 table records `NOT_RUN_PACKAGE_UNAVAILABLE` with the locked `normalization=NONE` and `transform=NONE` settings; no alternate second normalization was substituted.

## Covariate Sensitivity
Model 0 remains primary. Model 1, Model 3P, Model 3I, and Model 3S were run as separate sensitivity models only. Clinical Model 2 was not generated because age, sex, and stage are unavailable.

## Preprocessing Sensitivity
All locked Phase 6 sensitivity representations were analyzed using their precomputed CLR/rCLR matrices. Contaminant-exclusion analyses used the locked recomputed sensitivity matrices rather than dropping columns from the primary CLR matrix.

## Contamination Sensitivity
Candidate findings were annotated with contamination-risk categories and total-abundance-proxy sensitivity. Flagged genera are reported as potential-risk categories only; no genus is described as confirmed contamination because sequenced negative controls are absent.

## Sample Influence
Cook's distance, DFBETAs, leverage, studentized residuals, leave-one-sample-out ranges, and extreme-sample sensitivity outputs were generated. Influential samples were not automatically removed.

## Secondary Host Outcomes
Coactivation, Moffitt49 no-LEMD1, singscore contrast, PurIST basal probability, and assignment entropy were analyzed in separate BH families. Agreement across correlated host scores is not treated as independent replication.

## Descriptive Public-Subtype Results
Public Basal / Hybrid / Classical labels were used only descriptively for Aitchison PERMANOVA/PERMDISP and genus-level Kruskal-Wallis tests. These outputs are not elevated above the prespecified continuous primary analysis.

## Evidence Classification
Evidence category counts: `{"CONTAMINATION_SENSITIVE": 21, "METHOD_SENSITIVE": 67, "NO_SUPPORTED_ASSOCIATION": 23, "ROBUST_ASSOCIATION": 9, "SUGGESTIVE_ASSOCIATION": 2}`.

## Negative and Null Findings
The null and negative-results table explicitly retains global null results, absence of primary-FDR discoveries where applicable, transformation-dependent findings, contamination-sensitive findings, and covariate-sensitive findings. Nominal P values are not promoted when the primary q-value threshold is not met.

## Limitations
Clinical Model 2 could not be run because age, sex, and stage are unavailable. Sequenced negative controls are absent, so contamination assessments remain sensitivity annotations rather than definitive contamination calls. ESTIMATE-derived purity and immune/stromal scores come from the same host transcriptome and are robustness covariates, not independent measurements.

## Recommendation For Next Host-Mechanism Phase
Proceed to host-mechanism analysis only after carrying forward the Phase 7B evidence categories, prioritizing transcriptional-state interpretation and not treating nominal microbiome associations as discoveries.

## TO_VERIFY
- MaAsLin2 package execution remains `TO_VERIFY` because the package was not installed locally.
- ESTIMATE inferred purity remains `TO_VERIFY` as an inferred transcriptomic estimate rather than pathology-derived cellularity.
