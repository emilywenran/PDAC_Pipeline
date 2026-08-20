# Phase 3B Subtype Reproduction

## Methods executed

The locked primary GSE172356/Chan-Seng-Yue 94-gene hierarchical clustering method was executed on untransformed DESeq2 size-factor normalized counts using row median centering, row scaling, Pearson correlation distance, average linkage, and fixed dendrogram slices of 17 Basal, 23 Hybrid, and 22 Classical samples.

Verified secondary methods executed were the locked 49-active-gene Moffitt hierarchical clustering procedure and the 8-pair/16-gene PurIST classifier with the Phase 3A coefficients, intercept, and 0.5 basal-like probability cutoff.

## Methods not reproducible

Bailey and the full Chan-Seng-Yue 100-gene exploratory framework were not executed as subtype assignment methods because the Phase 3A inventory marks them `TO_VERIFY`/not directly reproducible without a pre-fitted single-sample classifier or exact locked implementation.

## Primary subtype counts

| Subtype | Reproduced n |
|---|---:|
| Basal | 17 |
| Hybrid | 23 |
| Classical | 22 |

## Agreement with public labels

The primary reproduction exactly matched the verified public labels for all 62 patients: exact agreement = 1.000; Cohen's kappa = 1.000. The primary confusion matrix is written to `05_results/tables/phase3b_confusion_matrices.tsv`.

The public labels were not used as model-training inputs; they were used only after the locked assignments were generated to calculate agreement metrics. No supervised model, feature selection, or parameter optimization was performed.

Score direction is documented as follows: higher PurIST basal probability indicates stronger basal-like evidence; higher Moffitt basal-minus-classical score indicates movement toward the Moffitt basal axis; the primary GSE172356 method has no locked continuous confidence score and uses dendrogram order only.

## Hybrid samples

Hybrid samples were preserved for the primary three-class reproduction. For binary secondary frameworks, public Hybrid samples were reported separately rather than counted as automatic errors. PurIST public-Hybrid distribution: {'Classical': 20, 'Basal-like': 3}. Moffitt public-Hybrid distribution: {'Others': 12, 'Classical': 7, 'Basal': 4}.

## Discordant and ambiguous samples

Primary discordant samples: 0. Moffitt `Others` assignments are retained as method-defined non-basal/classical calls and are flagged in `ambiguous_assignment`; PurIST has no locked ambiguous class and reports probability/confidence categories.

## Phase 2B outlier sensitivity

Excluding `YX16135T`, `YX16158T`, `YX16194T`, and `YX16224T` retained exact agreement for the remaining samples under the prespecified 17/19/22 slice sizes. See `05_results/tables/phase3b_sensitivity_summary.tsv` for assignment-change counts and log2 stress-test results.

## Missing signature genes

The primary method used 94 of the original 100 Chan-Seng-Yue-derived genes. The six unavailable genes were `C11orf70`, `C15orf52`, `RP11-400G3.5`, `DPCR1`, `FAM105A`, and `RP11-77K12.7`; they are absent from the source matrix rather than recoverable by imputation. Moffitt used 49 active genes after locked `LEMD1` exclusion. PurIST used all 16 genes.

## Phase 4 readiness

Phase 4 subtype stability analysis may proceed using the reproduced primary labels, with the caveat that exploratory Bailey/full Chan-Seng-Yue frameworks remain unresolved and must not be treated as validated assignment methods.

## Unresolved issues labelled TO_VERIFY

Exploratory TO_VERIFY methods are not executed: Bailey,Chan-Seng-Yue
