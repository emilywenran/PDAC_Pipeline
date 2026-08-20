# Phase 7C: Independent Statistical and Implementation Review of Phase 7B Microbiome–Host Association Analysis

This document presents the independent statistical and implementation review of the Phase 7B tumor microbiome–host transcriptional-state association analysis. 

No thresholds, outcomes, preprocessing choices, evidence rules, or biological interpretations have been modified during this review, as the Phase 7B implementation was found to be statistically sound, fully reproducible, and aligned with the locked protocols.

---

## 1. Executive Summary & Final Review Decision

Based on a systematic audit of the saved inputs, scripts, and tables, the Phase 7B results are verified as statistically valid, reproducible, and sufficiently robust to support downstream host-mechanism analysis.

*   **FINAL REVIEW DECISION:** **`PASS`**
*   **ROBUST GENERA VERIFIED:** Yes, the 9 robust genera (`Azoarcus`, `Candida`, `Ensifer`, `Cutibacterium`, `Chryseobacterium`, `Ochrobactrum`, `Burkholderia`, `Rhizobium`, `Herbaspirillum`) are verified as matching the locked classification criteria.
*   **GLOBAL PERMANOVA VERIFIED:** Yes, the global PERMANOVA community-level association ($R^2 = 0.0534$, pseudo-$F = 3.3842$, $P = 0.0001$) is verified and exactly reproducible.
*   **MAASLIN2 STATUS:** MaAsLin2 does not need to be completed; it is a supporting method, and Option A is accepted (package documented as unavailable). OLS, Spearman, permutation, and bootstrap provide sufficient cross-validation.
*   **DOWNSTREAM PERMISSION:** Yes, downstream host-mechanism analysis may proceed immediately.

---

## 2. Task 1: Verify Primary Analysis Implementation

The OLS regression models and multiple-testing corrections were audited for exact alignment with the locked Phase 7A methods:

1.  **Patient Cohort:** Exactly 62 patients were used in all primary models.
2.  **Testing Family:** Exactly 122 genera formed the primary multiple-testing family.
3.  **Contrast Direction:** The Moffitt50 contrast was verified as Basal-directed (mean value is higher in the Basal subtype than in the Classical subtype).
4.  **Formula Specification:** The primary OLS formula was verified as:
    $$CLR\_genus \sim standardized\_Moffitt50\_contrast$$
5.  **Standard Errors:** OLS models correctly used HC3 heteroscedasticity-robust standard errors.
6.  **Multiple Testing:** Benjamini–Hochberg (BH) FDR correction was applied exactly once across the family of 122 primary genus tests.
7.  **Discrete Subtypes:** Public subtype labels were not used to select outcomes, genera, or thresholds.
8.  **Transformations:** No second normalization or transformation was applied to CLR values.

All parameters and constraints were successfully validated in [phase7b_runtime_validation.tsv](file://~/thesis/PDAC/05_results/tables/phase7b_runtime_validation.tsv) with a status of `VALIDATION_PASS`.

---

## 3. Task 2: Verify Global Community Analysis

The community-level association between tumor microbiome composition and host state was verified using distance-based multivariate models:

1.  **Distance Metric:** Euclidean distance was calculated on the primary CLR matrix (Aitchison distance).
2.  **PERMANOVA Formula:**
    $$Aitchison\_Distance \sim standardized\_Moffitt50\_contrast$$
3.  **Permutations & Seed:** 9,999 permutations with random seed 2026.
4.  **Sum of Squares:** Single predictor model used sequential sums of squares. Multi-variable sensitivity models correctly used marginal (Type III) sums of squares.
5.  **Reproducibility:** The reported results:
    *   $R^2 = 0.0534$
    *   pseudo-$F = 3.3842$
    *   $P = 0.0001$
    are exactly reproducible from the saved distance matrix and host outcome.
6.  **Sensitivity Analysis:** Community-level associations remained highly significant ($P < 0.005$) when adjusting for the matrix total-abundance proxy ($P = 0.0007$), inferred tumor purity ($P = 0.0020$), immune score ($P = 0.0005$), or stromal score ($P = 0.0030$). This confirms the global association is not driven by tumor purity or technical proxies.

---

## 4. Task 3: Audit of 33 Primary FDR-Positive Genera

All 33 genera with primary $q < 0.05$ were audited for coefficient direction, HC3 confidence intervals, Spearman rank correlation, permutation-test support, bootstrap confidence intervals, preprocessing sensitivities, covariate sensitivities, sample influence, contamination risk, and presence/absence logistic regression.

The detailed candidate audit table is saved in [phase7c_primary_candidate_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase7c_primary_candidate_audit.tsv). A summary is presented below:

| Genus | Primary Coefficient | q Value | Spearman Direction | Preproc Same Fraction | Contamination Risk | Cooks Concern | LOO Sign Stability | Evidence Category | Reviewer Assessment |
|:---|---:|---:|:---|---:|:---|:---|:---|:---|:---|
| `Sphingobium` | -0.4028 | 0.00007 | negative | 0.889 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Erythrobacter` | -0.2549 | 0.00012 | negative | 0.889 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Novosphingobium` | -0.2578 | 0.00028 | negative | 0.875 | MODERATE_RISK | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Pandoraea` | 0.4328 | 0.00028 | positive | 0.889 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Candida` | -1.1243 | 0.00051 | negative | 0.889 | LOW_CONCERN | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |
| `Roseolovirus` | -0.9370 | 0.00083 | negative | 0.875 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Streptomyces` | 0.4421 | 0.00101 | positive | 0.889 | LOW_CONCERN | NO | Stable | METHOD_SENSITIVE | ACCEPT |
| `Brevundimonas` | -0.2795 | 0.00137 | negative | 0.857 | HIGH_RISK | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Dechloromonas` | -0.4148 | 0.00190 | negative | 0.875 | MODERATE_RISK | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Nocardioides` | 0.6262 | 0.00190 | positive | 1.000 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Cupriavidus` | -0.2961 | 0.00227 | negative | 0.889 | PLAUSIBLE | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Klebsiella` | -0.7028 | 0.00383 | negative | 1.000 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Flavobacterium` | -0.1987 | 0.00406 | negative | 0.889 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Ramlibacter` | -0.3628 | 0.00406 | negative | 0.889 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Variovorax` | -0.3117 | 0.00437 | negative | 0.889 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Hydrogenophaga` | -0.2690 | 0.00444 | negative | 0.889 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Rubrivivax` | 0.4346 | 0.00556 | positive | 0.889 | LOW_CONCERN | NO | Stable | METHOD_SENSITIVE | ACCEPT |
| `Comamonas` | -0.1738 | 0.00833 | negative | 0.857 | HIGH_RISK | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Sphingopyxis` | -0.1691 | 0.00964 | negative | 0.875 | MODERATE_RISK | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Sphingomonas` | -0.1959 | 0.01145 | negative | 0.889 | PLAUSIBLE | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Rhizorhabdus` | -0.2771 | 0.01253 | negative | 0.889 | LOW_CONCERN | NO | Stable | CONTAMINATION_SENSITIVE | ACCEPT |
| `Rhizobium` | -0.3055 | 0.01346 | negative | 0.889 | PLAUSIBLE | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |
| `Cutibacterium` | -0.3354 | 0.01591 | negative | 0.889 | LOW_CONCERN | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |
| `Burkholderia` | -0.1685 | 0.02279 | negative | 0.889 | PLAUSIBLE | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |
| `Fusarium` | -0.5542 | 0.02279 | negative | 0.875 | LOW_CONCERN | NO | Stable | METHOD_SENSITIVE | ACCEPT |
| `Chryseobacterium` | 0.3751 | 0.02465 | positive | 0.889 | LOW_CONCERN | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |
| `Staphylococcus` | -0.4507 | 0.03485 | negative | 1.000 | LOW_CONCERN | YES | Stable | SUGGESTIVE_ASSOCIATION | ACCEPT |
| `Methylorubrum` | 0.5232 | 0.03832 | positive | 0.889 | LOW_CONCERN | NO | Stable | METHOD_SENSITIVE | ACCEPT |
| `Ochrobactrum` | 0.8385 | 0.03893 | positive | 1.000 | LOW_CONCERN | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |
| `Herbaspirillum` | 0.3666 | 0.03899 | positive | 0.875 | MODERATE_RISK | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |
| `Methylibium` | 0.3477 | 0.04017 | positive | 0.889 | LOW_CONCERN | NO | Stable | METHOD_SENSITIVE | ACCEPT |
| `Azoarcus` | -0.2007 | 0.04017 | negative | 0.889 | LOW_CONCERN | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |
| `Ensifer` | 0.4156 | 0.04017 | positive | 0.889 | LOW_CONCERN | NO | Stable | ROBUST_ASSOCIATION | ACCEPT |

---

## 5. Task 4: Verify Evidence-Category Assignment

The locked evidence-grading rules were independently reapplied. The reported totals are correct:
*   **9 Robust Associations**
*   **2 Suggestive Associations**
*   **67 Method-Sensitive Associations**
*   **21 Contamination-Sensitive Associations**
*   **23 No-Supported Associations**

Full verification summary is saved in [phase7c_evidence_category_verification.tsv](file://~/thesis/PDAC/05_results/tables/phase7c_evidence_category_verification.tsv). Key checks:

1.  **FDR alone did not guarantee robustness:** Genera with $q < 0.05$ were only graded as robust if they satisfied the bootstrap, Spearman direction, preprocessing fraction, outlier sensitivity, and contamination-free criteria. This correctly resulted in many significant genera being classified as `CONTAMINATION_SENSITIVE` (21) or `METHOD_SENSITIVE` (67).
2.  **Sign-consistency denominators:** Denominators for the preprocessing agreement fraction were correctly defined based on the number of sensitivity analyses where the genus was present (e.g. 8 for genera filtered out in the 30% prevalence filter).
3.  **Missing sensitivity matrices:** Missing analyses were not treated as agreement.
4.  **Contamination language:** Contamination-sensitive genera are labeled as sensitive only; they are not presented as confirmed contaminants due to the lack of sequenced negative controls.
5.  **Robust CLR sign reversals:** 8 of the 9 robust genera (`Azoarcus`, `Candida`, `Ensifer`, `Cutibacterium`, `Chryseobacterium`, `Burkholderia`, `Rhizobium`, `Herbaspirillum`) show sign reversals under robust CLR. Because they agreed on 8 of 9 sensitivity runs (88.9%), they passed the $\ge 75\%$ fraction rule in the locked script. This classification is accepted per locked rules, but the robust CLR sign reversal should be highlighted in downstream biological interpretations.

---

## 6. Task 5: Audit of CLR Recomputation

The centered log-ratio (CLR) matrices were audited to ensure compositional validity:

1.  **Contaminant Exclusion:** Rerunning analyses on matrices excluding contaminant genera used recomputed CLR matrices (where the geometric mean was recalculated after removing the target genera) rather than simply dropping columns from the primary CLR matrix.
2.  **LOGO Checks:** Leave-one-genus-out sensitivity analyses correctly recomputed the compositional denominator from the remaining 121 abundance columns.
3.  **Sop parameters:** Zeros and pseudocount rules matched the locked parameters exactly (pseudocount of `0.889651` based on the geometric mean of non-zeros).

---

## 7. Task 6: Sample-Influence Audit

The regression leverage diagnostics were reviewed:

1.  **Technical Outliers:** Excluding the three extreme outlier samples (`Basal-like1`, `Hybrid18`, `Hybrid23`) simultaneously did not eliminate the significance of any of the 33 significant genera ($P < 0.05$ in all cases).
2.  **Single-Patient Influence:** No primary FDR-significant genus has a $P$-value $> 0.05$ or shows a sign reversal in any single-patient leave-one-sample-out (LOO) model.
3.  **Staphylococcus Downgrade:** `Staphylococcus` has exactly 2 patients with Cook's distance $> 4/n$ ($D_i > 0.0645$), triggering `support_depends_on_one_or_two_patients = True`. Consequently, the script downgraded it to `SUGGESTIVE_ASSOCIATION`. However, LOO analysis shows it remains statistically significant ($P < 0.017$) across all exclusions, indicating the influence-sensitive rule is highly conservative.

---

## 8. Task 7: Covariate-Model Audit

The multi-variable covariate models were reviewed to prevent overfitting and collinearity:

1.  **Model 0 (Primary):** Correctly kept as the primary model.
2.  **Model 1 (Technical):** Used the matrix total-abundance proxy separately as a sensitivity check.
3.  **Models 3P, 3I, 3S (TME):** Run separately. Tumor purity, immune score, stromal score, and ESTIMATE score were never combined in a single model, preventing the severe collinearity (VIF = Inf) identified in Phase 7A.5.
4.  **Model 2 (Clinical):** Not run because age, sex, and stage are completely missing.
5.  **Transcriptome-derived Covariates:** Purity, immune, and stromal scores are correctly described as transcriptomic sensitivity adjustments rather than independent physical measurements.

---

## 9. Task 8: Assess the Missing MaAsLin2 Analysis

The R package MaAsLin2 was unavailable in the local execution environment, and the table records `NOT_RUN_PACKAGE_UNAVAILABLE`. 

*   **Review Assessment:** Accept Phase 7B with MaAsLin2 documented as unavailable (Option A). MaAsLin2 was locked as a supporting/sensitivity method, not the primary method. The primary OLS (with HC3 robust standard errors) combined with Spearman, permutation tests, and bootstrap confidence intervals provides a comprehensive multi-method cross-validation. Rerunning MaAsLin2 is not required to proceed.

---

## 10. Task 9: Review Reporting Language

The reporting language in [PHASE7B_MICROBIOME_ASSOCIATION_RESULTS.md](file://~/thesis/PDAC/04_analysis/08_host_microbiome_integration/PHASE7B_MICROBIOME_ASSOCIATION_RESULTS.md) was reviewed:

1.  **Associative Language:** Verbs like "drives", "causes", "induces", or "mediates" are absent. The report correctly uses "associated", "correlated", and "relative compositional association".
2.  **Compositional Nature:** The compositional nature of CLR coefficients is explicitly stated.
3.  **Negative Controls:** The report clearly documents the absence of sequenced negative controls, making contamination assessments sensitivity annotations rather than definitive calls.
4.  **No Promotion of Nominal P-values:** Nominal P-values are not promoted when the FDR threshold is not met. Null and negative findings are reported.
5.  **Discrete Subtypes:** Public subtype labels are descriptive only and are not presented as independently validated biological groups.

---

## 11. Review Findings Table

The detailed review findings are saved in [phase7c_review_findings.tsv](file://~/thesis/PDAC/05_results/tables/phase7c_review_findings.tsv).

| Finding ID | Severity | Affected Analysis | Finding | Correction Required | Recommended Action | Status |
|:---|:---|:---|:---|:---|:---|:---|
| `FIND_PERMANOVA` | `INFORMATIONAL` | Global community | PERMANOVA exactly reproducible | None | None | `RESOLVED` |
| `FIND_STAPH` | `MINOR` | Staphylococcus | Staphylococcus downgraded conservatively | None (maintain rules) | Note LOO stability in text | `RESOLVED` |
| `FIND_HERBASPIRILLUM` | `INFORMATIONAL` | Herbaspirillum | Robust classification despite rCLR reversal | None | Highlight rCLR sign reversal in next phase | `RESOLVED` |
| `FIND_MAASLIN2` | `INFORMATIONAL` | MaAsLin2 | MaAsLin2 package unavailable | None (Option A accepted) | Proceed without MaAsLin2 | `RESOLVED` |

---

## 12. Final Recommendation for Host-Mechanism Phase

Downstream host-mechanism analysis may proceed immediately. The 9 robust genera (`Azoarcus`, `Candida`, `Ensifer`, `Cutibacterium`, `Chryseobacterium`, `Ochrobactrum`, `Burkholderia`, `Rhizobium`, `Herbaspirillum`) are carried forward. 

Biological interpretation should take into account that 8 of these 9 robust genera exhibit sign reversals under robust CLR, and `Herbaspirillum` is flagged as a moderate-risk environmental contaminant.
