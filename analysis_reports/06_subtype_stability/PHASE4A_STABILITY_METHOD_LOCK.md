# Phase 4A Method Lock: PDAC Subtype Stability Evaluation

This document defines the statistical framework and parameters for evaluating whether the basal-like, classical, and hybrid pancreatic ductal adenocarcinoma (PDAC) molecular subtypes represent stable, discrete transcriptomic clusters in the GSE172356 cohort. 

In accordance with Phase 4A requirements, no final stability clustering, differential expression, survival analysis, or microbiome association is performed. This document locks the analytical pipeline, metrics, sensitivity runs, and decision rules.

---

## 1. Primary Questions to Address

The stability analysis will systematically answer the following questions:
1. **Discrete Support**: Does the cohort's transcriptomic structure support the existence of three stable, discrete clusters?
2. **Optimal Cluster Number (K)**: Is $K=3$ better supported by stability metrics than other reasonable choices, specifically $K=2$, $K=4$, $K=5$, or $K=6$?
3. **Subtype-Specific Stability**: Are the basal-like and classical samples more stable (exhibiting higher item consensus and lower assignment entropy) than the hybrid samples?
4. **Hybrid Cluster Validity**: Does the "hybrid" subtype form a reproducible, distinct cluster, or does it represent an intermediate transition state or heterogeneous collection of samples?
5. **Sensitivity to Parameters**: How sensitive are the sample assignments and overall stability metrics to changes in expression scale, gene-filtering thresholds, sample exclusion (specifically the four suspected outliers), and resampling?

---

## 2. Required Method Hierarchy

To answer these questions without bias, we establish a structured hierarchy of representation and comparison methods. Unsupervised clustering will be completed in a blind fashion relative to the public labels; public subtype labels may only be compared to the resulting clusters *after* the unsupervised clustering and stability evaluations are complete.

### A. Primary Analysis: Reproduction Framework
- **Objective**: Reproduce the original subtyping framework's unsupervised stability using its exact genes and locked preprocessing settings.
- **Gene Set**: The Chan-Seng-Yue (CSY) 94-gene signature (verified in Phase 3A/3B).
- **Preprocessing**: Untransformed DESeq2 size-factor normalized counts, row median subtraction, and row centering/scaling (mean subtraction followed by standard deviation scaling).
- **Distance Metric**: Pearson correlation distance ($1 - r$).
- **Clustering Algorithm**: Hierarchical clustering with average linkage.
- **Candidate Cluster Numbers**: Evaluated for $K \in \{2, 3, 4, 5, 6\}$.
- **Nomenclature/Label Alignment**: To compare clusters with reference public labels, the Hungarian algorithm (minimum cost bipartite matching) or maximum-overlap assignment will be applied to map unsupervised cluster IDs ($1..K$) to biological labels only *after* clustering is finished.

### B. Secondary Analysis: Independent Unsupervised Representation
- **Objective**: Evaluate whether the natural, unbiased expression structure of the cohort supports three clusters when genes are selected without using public subtype labels.
- **Gene Selection Rule**: The top 1000 High-Variance Genes (HVGs) selected by Median Absolute Deviation (MAD) on the $\log_2(\text{normalized count} + 1)$ matrix. This selection is fixed *a priori* and cannot be tuned.
- **Preprocessing**: $\log_2(\text{normalized count} + 1)$ transformation, followed by row Z-score scaling (mean subtraction, standard deviation scaling).
- **Distance Metric**: Euclidean distance.
- **Clustering Algorithm**: Hierarchical clustering with Ward's linkage (`ward.D2`).
- **Candidate Cluster Numbers**: Evaluated for $K \in \{2, 3, 4, 5, 6\}$.

### C. Independent Biological-Axis Comparison
- **Objective**: Use the independent classifiers, Moffitt and PurIST, strictly to characterize the samples along a continuous basal-classical positioning.
- **Constraint**: No "hybrid" class will be invented or constructed for these binary classification frameworks.
- **Moffitt representation**: Standardized continuous score computed as the mean expression of the 25 basal genes minus the mean expression of the 24 classical genes (using modern gene symbols and row-scaled values).
- **PurIST representation**: Continuous probability of basal-like assignment ($p$) calculated from the locked k-TSP logistic regression formula.

---

## 3. Stability Statistics to Prespecify

We prespecify 12 stability statistics to evaluate the robustness of clusters at the cohort, cluster, and individual sample levels. The resampling protocol is fixed at **1,000 iterations** using a **sample-subsampling fraction of 0.80** (49 out of 62 samples drawn without replacement per run) and a **feature-subsampling fraction of 0.80** where specified. The random seed is fixed at **`2026`** for reproducibility.

1. **Consensus Clustering Matrix ($M$)**: A $N \times N$ matrix (where $N=62$ or $N=58$) tracking the fraction of times sample $i$ and sample $j$ cluster together given that they were both selected in the subsample.
2. **Consensus Cumulative Distribution Function (CDF)**: The empirical CDF of the consensus values in matrix $M$. A highly stable clustering is indicated by a bi-modal distribution (values concentrated near 0 and 1), resulting in a flat CDF curve in the intermediate range $(0.1, 0.9)$.
3. **Proportion of Ambiguous Clustering (PAC)**: Quantified as the fraction of sample pairs with consensus index values falling within the intermediate interval $[0.1, 0.9]$. Lower PAC indicates cleaner cluster boundaries and higher stability.
4. **Within-Cluster Consensus ($I_k$)**: The mean consensus index for all sample pairs assigned to the same cluster $k$. High $I_k$ indicates a tightly defined, cohesive cluster.
5. **Sample-Level Item Consensus ($I_i$)**: The mean consensus index between a specific sample $i$ and all other members of its assigned cluster. Samples with low $I_i$ are poorly integrated into their assigned group.
6. **Silhouette Width ($s_i$)**: Calculated on the full-dataset distance matrix using the final partition, representing how much closer sample $i$ is to its own cluster than to the next closest cluster. We will report the mean silhouette width per cluster and for the overall cohort.
7. **Bootstrap Cluster Jaccard Stability**: Calculated by drawing 1,000 bootstrap samples (sampling with replacement), performing clustering, and calculating the maximum Jaccard coefficient between the bootstrap clusters and the reference full-dataset clusters. A cluster is considered stable if its mean Jaccard coefficient is $\ge 0.75$, and highly stable if $\ge 0.85$.
8. **Prediction Strength**: Evaluated by splitting the dataset into training and test sets, clustering both, and measuring how well training-set centroids predict test-set cluster membership.
9. **Adjusted Rand Index (ARI) across Resampling Runs**: The pairwise ARI calculated between clusterings generated across the 1,000 subsampling iterations to measure overall partition consistency.
10. **Sample Assignment Entropy ($H_i$)**: Calculated from the frequency $p_{ik}$ of sample $i$ being assigned to cluster $k$ across all resampling runs where the sample was present:
    $$H_i = -\sum_{k=1}^K p_{ik} \log_2(p_{ik})$$
    An entropy of 0 indicates absolute stability, while high entropy indicates sample volatility/ambiguity.
11. **Co-Clustering Probability**: The distribution of off-diagonal elements in the consensus matrix $M$ for samples assigned to different clusters, indicating the degree of cluster leakage.
12. **Cluster Size and Minimum Viable Cluster Size**: The number of samples per cluster. The minimum viable cluster size for any candidate $K$ is set at **5 samples ($\ge 8\%$ of the cohort)**. Any partition yielding a cluster with $<5$ samples will be deemed biologically uninterpretable and rejected as an unstable partition.

---

## 4. Required Sensitivity Analyses

To ensure that the results are not artifacts of specific preprocessing choices, the stability evaluation will be executed across the 8 distinct runs defined in the parameter inventory:

1. **Full Cohort Analysis**: Baseline stability using all 62 samples.
2. **Outlier-Excluded Cohort Analysis**: Stability run excluding the four suspected outliers identified in Phase 2A/2B: `YX16135T`, `YX16158T`, `YX16194T`, and `YX16224T`.
3. **Normalized-Count Scale (Original Method)**: Evaluation using untransformed size-factor normalized counts (the primary representation).
4. **Log2-Transformed Scale**: Evaluation using the $\log_2(\text{normalized count} + 1)$ matrix under the primary clustering algorithm.
5. **Gene-wise Centering/Z-Scoring Policy**: Applied strictly when methodologically appropriate (e.g., row-wise median subtraction and scaling for the Pearson primary analysis; row-wise mean subtraction and standard deviation scaling for the Euclidean/Ward secondary analysis). 
6. **Alternative Missing-Value Handling**: Run stability comparisons using the complete-observation matrix versus the alternative gene-median imputed matrix.
7. **Gene-Filtering Threshold Variation**: Unsupervised stability using the secondary representation under alternative filtering thresholds (minimum count $\ge 10$ in $\ge 20\%$ of samples) and a smaller HVG set (top 500 HVGs).
8. **Sample Resampling (Bootstrap vs. Subsampling)**: Comparing the consensus subsampling (fraction 0.80) with bootstrap resampling (with replacement).
9. **Feature Resampling**: Executing simultaneous sample-subsampling (0.80) and feature-subsampling (80% of the 94 CSY signature genes) to evaluate signature robustness.

> [!IMPORTANT]
> Preprocessing and filtering strategies will NOT be selected based on which pipeline best reproduces the public labels. We will lock the primary preprocessing (untransformed, median-centered) as the replication reference, but assess the biological reality of the clusters by their stability across all sensitivity configurations.

---

## 5. Sample-Level Hybrid Assessment

To characterize the "hybrid" samples, we predefine a set of sample-specific metrics and qualitative interpretation categories. No samples are classified or assigned to these categories in Phase 4A.

### A. Metrics for Each Sample
- **Cluster Probability**: The fraction of resampling runs (in which the sample is selected) that it is assigned to cluster $k$ (for $k \in \{1..K\}$).
- **Item Consensus ($I_i$)**: The sample's mean consensus index with its assigned cluster.
- **Bootstrap Assignment Frequency**: The count of bootstrap runs yielding identical assignments.
- **Assignment Entropy ($H_i$)**: The quantitative uncertainty of the sample's cluster assignment.
- **Silhouette Width ($s_i$)**: The individual sample silhouette score.
- **Basal Score**: The continuous signature score (e.g., mean row-scaled expression of Basal signature genes).
- **Classical Score**: The continuous signature score (e.g., mean row-scaled expression of Classical signature genes).
- **Centroid Distance**: The Euclidean distance of the sample to the Basal centroid and the Classical centroid in the primary representation space.

### B. Predefined Interpretation Categories
After the stability metrics are computed, each sample will be assigned to one of the following categories:

*   **`STABLE_BASAL`**: Strong, consistent assignment to the basal-like cluster across resamples (entropy $H_i < 0.2$, item consensus $I_i \ge 0.85$, positive silhouette width).
*   **`STABLE_CLASSICAL`**: Strong, consistent assignment to the classical cluster across resamples (entropy $H_i < 0.2$, item consensus $I_i \ge 0.85$, positive silhouette width).
*   **`STABLE_HYBRID`**: Consistent assignment to a third "hybrid" cluster that is well-separated from the classical and basal poles (entropy $H_i < 0.3$, item consensus $I_i \ge 0.80$, positive silhouette width).
*   **`INTERMEDIATE_STATE`**: Moderate assignment frequency between classical and basal clusters, placing the sample between the two centroids, without forming a stable third cluster (characterized by intermediate PurIST probabilities $0.3 \le p \le 0.7$ and high assignment entropy when $K=3$).
*   **`HETEROGENEOUS_OR_UNSTABLE`**: Highly volatile sample exhibiting high assignment entropy ($H_i \ge 0.5$) across multiple resampling runs, representing noise or poor quality.
*   **`METHOD_SENSITIVE`**: Assignment flips completely depending on the preprocessing pipeline (e.g., untransformed counts vs. log2 counts) but remains stable within each pipeline.
*   **`TO_VERIFY`**: Exceeds outlier threshold, has high missingness, or presents conflicting clinical/taxonomic indicators requiring manual inspection.

---

## 6. Multi-Statistic Decision Rules

To determine whether the transcriptomic structure supports three stable clusters, we define four mutually exclusive qualitative conclusions. These decisions must be based on the intersection of multiple metrics, not on a single statistic (such as PAC or Jaccard stability alone).

### A. Phase 3B Reproduction Baseline Context
We explicitly state how the following Phase 3B results affect the stability design:
1. **Exact Reproduction on Normalized Counts**: The 94-gene primary hierarchical clustering matched the public labels exactly (agreement = 1.000) on the untransformed size-factor normalized count scale. This represents the baseline replication.
2. **Log2 Stress Test Volatility**: Applying the same median-centering and row-scaling procedure on the $\log_2(\text{normalized count} + 1)$ matrix resulted in **26 changed assignments** (agreement drops to **0.581**). 
3. **Implication**: This extreme sensitivity is treated as strong evidence that the primary clustering partition is highly sensitive to the count scale. This log2-scale volatility indicates that a discrete cluster boundaries model may be fragile, and reinforces the necessity of checking stability statistics before concluding that the three subtypes are biologically discrete entities.

### B. Objective Decision Rules

```
                       Consensus Matrix & PAC
                       │
       ┌───────────────┴───────────────┐
       ▼ (PAC < 0.15 for K=3)          ▼ (PAC >= 0.15 for K=3)
  Jaccard >= 0.75?                Evaluate K=2 PAC
       │                               │
 ┌─────┴─────┐                   ┌─────┴─────┐
 ▼ Yes       ▼ No                ▼ PAC < 0.15▼ PAC >= 0.15
[STRONG 3]  [PARTIAL 3]         [POLES/INT] [NO STRUCTURE]
```

1.  **Strong Support for Three Discrete Clusters**
    *   *Required Conditions*:
        *   The Proportion of Ambiguous Clustering (PAC) is minimized at $K=3$, and is $< 0.15$ in both the primary and secondary unsupervised runs.
        *   The bootstrap cluster Jaccard stability is $\ge 0.75$ for all three clusters.
        *   Within-cluster consensus ($I_k$) is $\ge 0.80$ for all three clusters.
        *   Average silhouette width is positive for all three clusters.
        *   Fewer than 10% of samples (6 samples) are classified as `INTERMEDIATE_STATE` or `HETEROGENEOUS_OR_UNSTABLE`.
2.  **Partial Support for Three Clusters**
    *   *Required Conditions*:
        *   The PAC is minimized or low at $K=3$ ($\text{PAC} < 0.25$), but Jaccard stability or within-cluster consensus is low ($\le 0.75$) for at least one cluster (typically the hybrid cluster).
        *   Unsupervised secondary clustering (HVGs) yields moderate agreement (ARI $\ge 0.40$) with the primary CSY classification.
        *   The primary classification is stable on the normalized count scale but exhibits high sample assignment entropy when log2 transformed.
3.  **Stronger Support for Two Major Poles with an Intermediate Group**
    *   *Required Conditions*:
        *   Stability metrics (CDF, PAC, Jaccard) are significantly superior for $K=2$ compared to $K=3$ (e.g., $K=2$ PAC is $< 0.15$, while $K=3$ PAC is $\ge 0.25$).
        *   The hybrid cluster in the $K=3$ run exhibits low Jaccard stability ($< 0.60$) and low within-cluster consensus ($I_k < 0.70$).
        *   Sample-level assignment entropy is concentrated in the hybrid samples (mean $H_{hybrid} > 0.5$, while mean $H_{basal, classical} < 0.2$).
        *   Continuous biological-axis projections (Moffitt and PurIST) reveal a unimodal or flat intermediate density between the two poles, rather than three discrete peaks.
4.  **No Sufficiently Stable Discrete Cluster Structure**
    *   *Required Conditions*:
        *   All candidate $K$ values yield high PAC ($\ge 0.30$) and low bootstrap Jaccard stability (mean Jaccard $< 0.60$).
        *   No stable boundaries can be drawn; consensus matrices display high off-diagonal co-clustering probabilities across all proposed boundaries.
        *   Average cohort silhouette width is near-zero or negative.
        *   Assignments are dominated by technical parameters (missingness, normalization scale) rather than biological features.

---

## 7. Operational Commit and Execution

As required by Phase 4A, this protocol is locked for subsequent automated execution. No scripts will perform clustering or compute stability results until Phase 4A is signed off.
