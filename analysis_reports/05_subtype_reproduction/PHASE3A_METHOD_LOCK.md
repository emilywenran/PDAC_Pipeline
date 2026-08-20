# Phase 3A Method Lock: PDAC Subtype Reproduction

This document defines the method lock, parameter settings, preprocessing steps, and verification standards for reproducing the PDAC transcriptomic subtypes for the GSE172356 cohort. No subtype assignment or classification is performed in Phase 3A.

---

## 1. Primary Reproduction Method

The primary method is the exact reproduction of the original public subtyping procedure used in Guo et al. (2021).

- **Gene Set**: A 94-gene subset of the Chan-Seng-Yue et al. (2020) signature. Originally derived by selecting the top 25 genes from each of the four signatures (Sig. 1, Sig. 2, Sig. 6, Sig. 10) in Supplementary Table 4 of `41588_2019_566_MOESM3_ESM.xlsx`. Six genes are missing in the GSE172356 dataset (`C11orf70`, `C15orf52`, `RP11-400G3.5`, `DPCR1`, `FAM105A`, `RP11-77K12.7`), leaving exactly 94 genes.
- **Preprocessing**: 
  1. Input scale: Size-factor normalized counts (untransformed DESeq2 normalized counts).
  2. Row median subtraction: Subtract the median value of each gene across all 62 samples.
  3. Row scaling: Center and scale rows within `pheatmap` (which subtracts the mean of the median-centered values and divides by their standard deviation).
- **Algorithm**: Hierarchical clustering using Pearson correlation distance ($1 - r$) and average linkage (hierarchical clustering `hc` via `pheatmap` and `ConsensusClusterPlus`).
- **Assignment Rules**: 
  - Column ordering is defined by the column dendrogram.
  - The 62 samples are partitioned into three groups based on their order in the dendrogram:
    - First 17 samples: **Basal**
    - Next 23 samples: **Hybrid**
    - Last 22 samples: **Classical**

---

## 2. Secondary Comparison Methods

Two established independent subtype classifiers will be evaluated in parallel to assess the robustness of the primary classification.

### A. Moffitt Tumor Basal-like/Classical Framework
- **Gene Set**: 50 genes (25 basal-like, 25 classical) from Moffitt et al. (2015).
  - Outdated or alternative gene symbols are mapped to modern counterparts in GSE172356:
    - `CTSL2` $\rightarrow$ `CTSV`
    - `ANXA8L2` $\rightarrow$ `ANXA8`
    - `ATAD4` $\rightarrow$ `FLAD1`
    - `LOC400573` $\rightarrow$ `TMEM238L`
  - The gene `LEMD1` is excluded for exact reproduction of the Guo et al. validation, resulting in 49 active genes.
- **Preprocessing**: Same as the primary method (untransformed size-factor normalized counts, row median centering, and row scaling in `pheatmap`).
- **Algorithm**: Hierarchical clustering with Pearson correlation distance and average linkage.
- **Assignment Rules**: The column dendrogram order is sliced into three hardcoded groups:
  - First 27 samples: **Classical**
  - Next 17 samples: **Basal**
  - Last 18 samples: **Others**

### B. PurIST (Purity Independent Subtyping of Tumors)
- **Gene Set**: 16 genes representing 8 Top Scoring Pairs (TSPs) from Rashid et al. (2020) (`fitteds_public_2019-02-12.Rdata`).
- **Preprocessing**: None. The method is rank-based and scale-invariant. It does not require row-centering, log-transformation, or scaling.
- **Algorithm**: k-TSP pairwise comparisons with coefficients from a penalized logistic regression model.
- **Assignment Rules**:
  - Calculate the predictor score $S$:
    $$S = -6.815 + 1.994 \cdot I(\text{GPR87} > \text{REG4}) + 2.031 \cdot I(\text{KRT6A} > \text{ANXA10}) + 1.618 \cdot I(\text{BCAR3} > \text{GATA6}) + 0.922 \cdot I(\text{PTGES} > \text{CLDN18}) + 1.059 \cdot I(\text{ITGA3} > \text{LGALS4}) + 0.929 \cdot I(\text{C16orf74} > \text{DDC}) + 2.505 \cdot I(\text{S100A2} > \text{SLC40A1}) + 0.485 \cdot I(\text{KRT5} > \text{CLRN3})$$
  - Calculate the basal-like probability:
    $$p = \frac{1}{1 + e^{-S}}$$
  - If $p > 0.5$, classify as **Basal-like**; if $p \le 0.5$, classify as **Classical**.
  - Confidence categories:
    - $p > 0.9$: **Strong Basal-like**
    - $0.5 < p \le 0.9$: **Lean/Likely Basal-like**
    - $0.1 \le p \le 0.5$: **Lean/Likely Classical**
    - $p < 0.1$: **Strong Classical**

---

## 3. Predefined Sensitivity Analyses

The following sensitivity cohorts are locked and will be analyzed in Phase 3B without optimization:

1. **Full Cohort**: All 62 samples.
2. **Outlier-Excluded Cohort**: Exclude the four Phase 2A suspected outliers: `YX16135T` (PDAC_016), `YX16158T` (PDAC_023), `YX16194T` (PDAC_033), and `YX16224T` (PDAC_039).
3. **Alternative Missingness Matrices**: Stress-test the clustering results using the 50% missingness threshold + gene-median imputed matrix, and the zero-filled matrix.
4. **Gene Scaling**: Evaluate whether performing median subtraction on log2-transformed counts rather than untransformed counts changes dendrogram column order.

---

## 4. Agreement and Comparison Metrics

Agreement between subtyping methods will be evaluated in Phase 3B using:
- **Confusion Matrix**: Visualizing cross-classifications.
- **Exact Agreement**: Percentage of identical labels.
- **Balanced Accuracy**: Accounting for class imbalance.
- **Cohen's Kappa**: Assessing agreement beyond chance.
- **Adjusted Rand Index (ARI)**: Evaluating clustering similarity.
- **Normalized Mutual Information (NMI)**: Information-theoretic overlap.
- **Per-Class Sensitivity**: Characterizing performance for Basal and Classical classes.
- **Assignment Confidence**: Evaluating the distribution of PurIST probabilities.

---

## 5. Exit Conditions for Phase 3B Approval

The transition from Phase 3A to Phase 3B requires:
1. **Dendrogram Verification**: Exact replication of the R script's dendrogram column order for GSE172356.
2. **Nomenclature Match**: Successful mapping of alternative gene symbols for the Moffitt signatures.
3. **Inventory Completion**: Lock of all signature file checksums.

---

## 6. Unresolved Items (TO_VERIFY)

- `TO_VERIFY`: The precise biological or technical reason for the literal `NA` values in the processed matrix remains unverified.
- `TO_VERIFY`: Exploratory frameworks (Bailey, Chan-Seng-Yue 100-gene NMF) are not suitable for direct single-sample reproduction in this project due to the lack of pre-fitted classifiers or exact replication code, and are marked exploratory.
