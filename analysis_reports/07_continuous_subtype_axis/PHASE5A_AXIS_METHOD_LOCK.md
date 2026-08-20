# Phase 5A Method Lock: PDAC Continuous Basal–Classical transcriptional-axis Analysis

This document prespecifies and locks the statistical framework and parameters for analyzing the pancreatic ductal adenocarcinoma (PDAC) transcriptomics in the GSE172356 cohort along a continuous basal–classical transcriptional axis. 

In accordance with Phase 5A requirements, no final scoring, differential expression, survival analysis, microbiome association, pathway enrichment, or target prioritization is performed. This document locks the analytical pipeline, continuous metrics, sensitivity runs, and decision rules.

---

## 1. Context and Baseline Decision

The Phase 4B stability evaluation of the GSE172356 cohort returned an overall decision of **`INCONCLUSIVE`**. 
- The primary Chan-Seng-Yue (CSY) normalized-count analysis preferred $K=2$ (PAC = 0.061 vs $K=3$ PAC = 0.196).
- The independent high-variance gene (HVG) analysis preferred $K=4$ (PAC = 0.347 at $K=3$).
- The primary $K=3$ partition matched the public labels post-clustering (agreement = 1.000), but bootstrap Jaccard stability for the clusters was sub-threshold (mean Jaccard = 0.559), and sample assignments were highly sensitive to count transformations (log2 transformation changed 30 assignments, agreement dropping to 0.385).
- Recurrently unstable samples (20 samples volatile in $\ge 3$ of 8 runs) were identified.

**Locked Policy**: The Phase 4B decision is preserved as **`INCONCLUSIVE`**. The continuous axis analysis will not interpret Phase 4B as proof of either $K=2$ or $K=3$. Instead, it aims to determine whether the publicly labeled "Hybrid" samples occupy an intermediate continuum, a co-activated hybrid state, a heterogeneous/unstable state, or a method-sensitive state in the continuous transcriptomic space.

---

## 2. Input Data and Checksums

The analysis will consume the following verified input matrices, signatures, and prior metrics:

| Input Name | File Path | Expected SHA256 Checksum |
| :--- | :--- | :--- |
| **Filtered Normalized counts** | [GSE172356_expression_filtered_normalized.tsv.gz](file://~/thesis/PDAC/03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz) | `8947ca75c3240177f8daeb8426e4cc9978a94c51ed17b14cb6eaf0146c4d73c1` |
| **Filtered Log2 counts** | [GSE172356_expression_log2_analysis_ready.tsv.gz](file://~/thesis/PDAC/03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz) | `13c16a95c7ef94e59b7d685c85b78f4bc2a2d22b9e6ffaafb929dd2a50c0328a` |
| **Moffitt 50-gene axis signature (Primary)** | [Moffitt_50_gene_axis.tsv](file://~/thesis/PDAC/02_data/reference/PDAC_subtype_signatures/Moffitt_50_gene_axis.tsv) | `3fa1790ff692898d01e2f4f8058d438c1263245c0d5316afab9840c968a2b72f` |
| **Moffitt 49-gene axis signature (Sensitivity)** | [Moffitt_49_gene_axis_no_LEMD1.tsv](file://~/thesis/PDAC/02_data/reference/PDAC_subtype_signatures/Moffitt_49_gene_axis_no_LEMD1.tsv) | `65cadb4c059a4b5be81efe03b8be1b5a6fc88937fd3eadf46f399ee007f1fc61` |
| **PurIST signature & coeffs** | [PurIST_signatures.tsv](file://~/thesis/PDAC/02_data/reference/PDAC_subtype_signatures/PurIST_signatures.tsv) | `b198e583e65c8e4f1da04e2054c24c23d201c7cada330535fe3f3645a11d249f` |
| **Phase 3B Sample-level labels** | [phase3b_primary_subtype_assignments.tsv](file://~/thesis/PDAC/05_results/tables/phase3b_primary_subtype_assignments.tsv) | *Reference label source* |
| **Phase 4B Sample-level stability** | [phase4b_sample_level_stability.tsv](file://~/thesis/PDAC/05_results/tables/phase4b_sample_level_stability.tsv) | *Stability metrics source* |

---

## 3. Seven Continuous Axis Scores

For each sample $i$, we will construct seven independent continuous measures along the basal–classical programs:

### A. Basal Program Score ($S_{basal}$)
Defined as the mean expression of the verified Moffitt basal genes on the row-standardized scale of the $\log_2(\text{normalized count} + 1)$ matrix:
- **Primary Analysis**: $|G_{basal}| = 25$ genes (including `LEMD1`).
- **Sensitivity Analysis (LEMD1-Excluded)**: $|G_{basal}| = 24$ genes (excluding `LEMD1`).
$$S_{basal, i} = \frac{1}{|G_{basal}|} \sum_{g \in G_{basal}} z_{g, i}$$
where $z_{g, i}$ is the row Z-score of gene $g$ in sample $i$.

### B. Classical Program Score ($S_{classical}$)
Defined as the mean expression of the verified Moffitt classical genes on the row-standardized scale of the $\log_2(\text{normalized count} + 1)$ matrix:
- **Both Analyses**: $|G_{classical}| = 25$ genes.
$$S_{classical, i} = \frac{1}{|G_{classical}|} \sum_{g \in G_{classical}} z_{g, i}$$

### C. Basal–Classical Contrast Score ($C_i$)
The primary continuous measure representing the relative position of sample $i$ along the basal-classical axis:
$$C_i = S_{basal, i} - S_{classical, i}$$
A highly positive $C_i$ indicates a strong basal-like bias, while a highly negative $C_i$ indicates a strong classical-like bias.

### D. Basal/Classical Co-activation Score ($A_i$)
Quantifies the simultaneous activation of both basal and classical programs (to test the "co-activated hybrid" hypothesis):
$$A_i = \min(S'_{basal, i}, S'_{classical, i})$$
where $S'_{basal}$ and $S'_{classical}$ are $S_{basal}$ and $S_{classical}$ rescaled to a $[0, 1]$ interval across the cohort:
$$S'_{program, i} = \frac{S_{program, i} - \min(S_{program})}{\max(S_{program}) - \min(S_{program})}$$

### E. Distance to Reference Centroids ($d_{basal, i}, d_{classical, i}$)
The Euclidean distance of sample $i$ to the Basal centroid ($\mathbf{\mu}_{basal}$) and Classical centroid ($\mathbf{\mu}_{classical}$) computed in the log2 Z-score expression space of the active genes (50 genes for Primary, 49 genes for LEMD1-Excluded Sensitivity):
$$d_{basal, i} = \sqrt{\sum_{g \in G_{Moffitt}} (z_{g, i} - \mu_{basal, g})^2}$$
$$d_{classical, i} = \sqrt{\sum_{g \in G_{Moffitt}} (z_{g, i} - \mu_{classical, g})^2}$$

### F. PurIST Basal Probability ($P_{basal, i}$)
The single-sample continuous probability output from the locked PurIST k-TSP logistic regression classifier (retrieved from Phase 3B outputs).

### G. Moffitt Basal–Classical Score Difference ($D_{Moffitt, i}$)
The standardized difference calculated directly from the locked Moffitt classifier (retrieved from Phase 3B outputs).

---

## 4. Scoring Approaches and Centroids

To ensure robustness, the continuous scoring will be conducted using two independent methods, neither of which will be tuned or optimized against public subtype labels.

### Primary Approach: Signature-Mean
- **Methodology**: Row-scaled Z-score means of the locked Moffitt basal (25 genes for Primary; 24 genes for LEMD1-Excluded Sensitivity) and classical (25 genes) sets, as defined in Section 3A-C.
- **Input Matrix**: $\log_2(\text{normalized count} + 1)$ expression.
- **Centering/Scaling**: Genes (rows) centered by median subtraction, scaled by standard deviation.

### Secondary Approach: Singscore (Rank-Based)
- **Methodology**: Single-sample rank-based scoring using the R package `singscore`. It ranks the genes within each sample and computes a normalized score based on the distribution of ranks of the basal and classical signature gene sets (50-gene for Primary; 49-gene for LEMD1-Excluded Sensitivity).
- **Robustness**: Independent of cross-sample normalization scale and robust to batch effects or outlier samples.
- **Outputs**: Independent scores for Basal ($SS_{basal, i}$) and Classical ($SS_{classical, i}$) programs, and a contrast score ($SS_{contrast, i} = SS_{basal, i} - SS_{classical, i}$).

### Reference Centroid Definitions (No Leakage)
To compute centroid distances ($d_{basal, i}$ and $d_{classical, i}$), we establish three distinct centroid definition rules (applied in the 50-gene space for Primary and 49-gene space for LEMD1-Excluded Sensitivity):

1.  **Reference-Anchored Centroids (Descriptive Comparison)**:
    - *Definition*: $\mathbf{\mu}_{basal}$ and $\mathbf{\mu}_{classical}$ are the mean vectors of the active Moffitt genes (50 or 49) across samples publicly labeled as `Basal` and `Classical` respectively (excluding public `Hybrid` samples).
    - *Leakage Label*: Explicitly labeled as **Descriptive Reference-Anchored Analysis** because it utilizes public label metadata.
2.  **Unsupervised $K=2$ Centroids (Zero Leakage)**:
    - *Definition*: Hierarchical clustering ($K=2$) is performed on the active Moffitt genes (50 or 49). The centroids $\mathbf{\mu}_{C1}$ and $\mathbf{\mu}_{C2}$ are computed as the mean vectors of the two resulting clusters. The cluster with the higher mean expression of basal genes is designated the basal centroid.
    - *Leakage Label*: **Zero-Leakage Unsupervised Centroid Analysis**.
3.  **Leave-One-Out (LOO) Centroid Sensitivity**:
    - *Definition*: For each sample $i$, the reference-anchored centroids $\mathbf{\mu}_{basal, -i}$ and $\mathbf{\mu}_{classical, -i}$ are calculated by excluding sample $i$ from the mean calculation. Distances $d_{basal, i}$ and $d_{classical, i}$ are then computed using these sample-specific centroids.
    - *Leakage Label*: **LOO Reference Sensitivity Centroid Analysis**.

---

## 5. Hybrid-State Metrics

To systematically characterize the state of public Hybrid samples, we define the following quantitative metrics:

1.  **Midpoint Proximity ($MP_i$)**:
    Measures how close sample $i$ is to the exact midpoint between the Basal and Classical centroids:
    $$MP_i = 1 - \frac{|d_{basal, i} - d_{classical, i}|}{d_{basal, i} + d_{classical, i}}$$
    A value close to $1.0$ indicates that the sample is equidistant from both poles.
2.  **Simultaneous Activation ($SA_i$)**:
    Equivalent to the co-activation score $A_i$ defined in Section 3D.
3.  **Method Score Discordance ($SD_i$)**:
    The absolute difference in ranks of sample $i$ between the Primary signature-mean contrast score $C_i$ and the Secondary singscore contrast score $SS_{contrast, i}$:
    $$SD_i = |\text{rank}(C_i) - \text{rank}(SS_{contrast, i})|$$
4.  **Phase 4B Resampling Assignment Entropy ($H_i$)**:
    Retrieved from `phase4b_sample_level_stability.tsv`. Higher entropy indicates that the sample volatilely switched clusters under resampling.
5.  **Phase 4B Cluster Item Consensus ($I_i$)**:
    Retrieved from `phase4b_sample_level_stability.tsv`. Measures sample integration into its assigned cluster.
6.  **Distance to Both Poles ($DP_i$)**:
    The sum of Euclidean distances to both centroids, normalized by the distance between the two centroids ($d_{bc}$):
    $$DP_i = \frac{d_{basal, i} + d_{classical, i}}{d_{bc}}$$
    An intermediate sample in a linear continuum will have $DP_i \approx 1.0$. A sample that is poorly described by both signatures will have $DP_i \gg 1.0$.
7.  **Method-to-Method Axis Variance ($MV_i$)**:
    The variance across all continuous scoring metrics (Primary contrast, Singscore contrast, PurIST probability, Moffitt score difference) after each score is standardized (Z-scored) across the cohort.

---

## 6. Predefined Sample Interpretation Categories

For Phase 5B execution, samples will be classified into one of the following mutually exclusive categories based on the locked rules below. No samples are assigned to these categories in Phase 5A.

```
                                  Is sample close to a pole?
                                  (d_pole < 0.4 * d_bc)
                                  │
                  ┌───────────────┴───────────────┐
                  ▼ Yes                           ▼ No
            Assign to POLE                 Evaluate DP_i and SA_i
                                                  │
                 ┌────────────────────────────────┴────────────────────────────────┐
                 ▼ (DP_i <= 1.25)                                                  ▼ (DP_i > 1.25)
         Evaluate Co-activation (SA_i)                                     Evaluate Discordance (MV_i)
                 │                                                                 │
       ┌─────────┴─────────┐                                             ┌─────────┴─────────┐
       ▼ (SA_i < 0.4)      ▼ (SA_i >= 0.4)                               ▼ (MV_i < 0.5)      ▼ (MV_i >= 0.5)
 [INTERMEDIATE_CONTINUUM] [COACTIVATED_HYBRID]                     [HETEROGENEOUS_UNSTABLE] [METHOD_SENSITIVE]
```

*   **`BASAL_POLE`**:
    - *Criteria*: Distance to basal centroid $d_{basal, i} < 0.4 \times d_{bc}$, and contrast score $C_i > 0.5 \times \text{sd}(C)$.
*   **`CLASSICAL_POLE`**:
    - *Criteria*: Distance to classical centroid $d_{classical, i} < 0.4 \times d_{bc}$, and contrast score $C_i < -0.5 \times \text{sd}(C)$.
*   **`INTERMEDIATE_CONTINUUM`**:
    - *Criteria*: Midpoint proximity $MP_i \ge 0.7$, Distance to both poles $DP_i \le 1.25$, and co-activation score $SA_i < 0.4$. (Indicates the sample lies on a linear path between the poles without high co-expression).
*   **`COACTIVATED_HYBRID`**:
    - *Criteria*: Distance to both poles $DP_i \le 1.25$, and co-activation score $SA_i \ge 0.4$. (Indicates simultaneous high expression of both basal and classical marker genes).
*   **`HETEROGENEOUS_OR_UNSTABLE`**:
    - *Criteria*: Distance to both poles $DP_i > 1.25$, and Phase 4B assignment entropy $H_i \ge 0.5$. (Indicates the sample is far from both poles and shows high clustering instability).
*   **`METHOD_SENSITIVE`**:
    - *Criteria*: Method-to-method axis variance $MV_i \ge 0.5$, or rank discordance $SD_i > 15$ positions. (Indicates sample placement depends heavily on the algorithm).
*   **`TO_VERIFY`**:
    - *Criteria*: Sample was flagged as a Phase 2B outlier candidate, or has missing clinical/expression features.

---

## 7. Predefined Statistical Evaluations

To evaluate the mathematical properties of the continuous axis, we predefine 10 statistical checks to be executed in Phase 5B. The random seed is fixed at **`2026`** for all bootstrap and permutation procedures.

1.  **Score Distributions by Public Subtype**:
    - Calculate descriptive statistics (M, SD, median, IQR) of continuous scores ($C_i, SA_i$) grouped by public labels.
    - Evaluate overall differences using Kruskal-Wallis non-parametric ANOVA (alpha = 0.05).
2.  **Ordered Trend Test Across Subtypes**:
    - Test the hypothesis that expression follows an ordered sequence: $\text{Classical} \rightarrow \text{Hybrid} \rightarrow \text{Basal}$.
    - Implement the **Jonckheere-Terpstra (JT) test** (using the `Clinfun` R package or custom permutation) to calculate the standardized JT statistic ($J$) and its p-value.
3.  **Effect Sizes and Confidence Intervals**:
    - Calculate Hedges' $g$ effect size for the difference in contrast scores between the Basal and Classical poles.
    - Compute 95% bootstrap confidence intervals for $g$ using 1,000 bootstrap resamples.
4.  **Spearman Correlations between Scoring Systems**:
    - Compute pairwise Spearman's rank correlation coefficients ($\rho$) between all scoring systems: Primary contrast ($C$), Singscore contrast ($SS_{contrast}$), PurIST probability ($P_{basal}$), and Moffitt difference ($D_{Moffitt}$).
    - Correct p-values for multiple testing using the Benjamini-Hochberg (BH) procedure.
5.  **Concordance with Phase 4B Stability Metrics**:
    - Compute Spearman correlation between continuous hybrid metrics ($MP_i, SA_i, DP_i$) and stability metrics ($H_i, I_i$).
    - Test if intermediate proximity or co-activation correlates with cluster assignment volatility.
6.  **Multimodality Testing**:
    - Test whether the continuous contrast score $C_i$ is distributed unimodally (suggesting a single continuum) or multimodally (suggesting discrete clusters).
    - Implement **Hartigan’s Dip Test** on $C_i$ using the R package `diptest`. A significant p-value ($p < 0.05$) rejects the null hypothesis of unimodality.
7.  **Sensitivity to Phase 2B Outlier Candidates**:
    - Re-calculate centroids and sample scores excluding the four outlier candidates: `YX16135T`, `YX16158T`, `YX16194T`, and `YX16224T`.
    - Compute the Pearson correlation between outlier-included and outlier-excluded scores.
8.  **Sensitivity to Input Scale (Normalized-Count vs. Log2)**:
    - Compare the primary scoring system computed on the log2 counts matrix with the same calculation on untransformed normalized counts.
    - Compute the fraction of samples whose interpretation category changes.
9.  **Bootstrap Confidence Intervals**:
    - Draw 1,000 bootstrap samples (with replacement).
    - Re-estimate the centroids, mean program scores, and Spearman correlation coefficients.
    - Define the 95% CI as the 2.5th and 97.5th percentiles of the bootstrap distribution.
10. **Permutation Tests**:
    - Perform 1,000 random shuffles of sample assignments to estimate the null distribution for the Jonckheere-Terpstra trend test and Spearman correlation tests.

---

## 8. Biological Hypothesis Evaluation

The continuous axis results will be used to test four competing biological models. We predefine the theoretical criteria for each:

*   **Hypothesis A: Two Poles with an Intermediate Continuum**
    - *Supported if*: Hartigan's Dip Test on $C_i$ is non-significant ($p \ge 0.05$, unimodal); midpoint proximity $MP_i$ is high for Hybrid samples; co-activation $SA_i$ is low; and distance to both poles $DP_i \approx 1.0$.
*   **Hypothesis B: Co-activated Hybrid Program**
    - *Supported if*: Hybrid samples exhibit high co-activation scores ($SA_i \ge 0.4$) and distance to both poles $DP_i \le 1.25$.
*   **Hypothesis C: Heterogeneous Collection of Unstable Samples**
    - *Supported if*: Hybrid samples exhibit high distance to both poles ($DP_i \gg 1.0$) and high Phase 4B assignment entropy ($H_i \ge 0.5$).
*   **Hypothesis D: Multiple Distinct Hybrid Mechanisms**
    - *Supported if*: Hybrid samples split into sub-groups, with some exhibiting co-activation ($SA_i \ge 0.4$, $DP_i \le 1.25$) and others exhibiting heterogeneous/unstable profiles ($DP_i > 1.25, H_i \ge 0.5$).

---

## 9. Multi-Metric Decision Rules

To synthesize the statistical findings, we lock five mutually exclusive qualitative conclusions. The final conclusion will be selected according to these exact conditions:

1.  **`TWO_POLES_WITH_INTERMEDIATE_CONTINUUM`**:
    - Hartigan's Dip Test on $C_i$ is non-significant ($p \ge 0.05$).
    - Median co-activation $SA_i$ of public Hybrid samples is $< 0.3$.
    - Median distance to both poles $DP_i$ of public Hybrid samples is $\le 1.15$.
    - More than 60% of public Hybrid samples are classified as `INTERMEDIATE_CONTINUUM`.
2.  **`TWO_POLES_WITH_COACTIVATED_HYBRID`**:
    - Hartigan's Dip Test on $C_i$ is significant ($p < 0.05$) or unimodality is rejected.
    - Median co-activation $SA_i$ of public Hybrid samples is $\ge 0.4$.
    - Median distance to both poles $DP_i$ of public Hybrid samples is $\le 1.15$.
    - More than 60% of public Hybrid samples are classified as `COACTIVATED_HYBRID`.
3.  **`HETEROGENEOUS_HYBRID_STATES`**:
    - Median distance to both poles $DP_i$ of public Hybrid samples is $> 1.25$.
    - Median Phase 4B assignment entropy $H_i$ of public Hybrid samples is $\ge 0.4$.
    - Public Hybrid samples are split between `HETEROGENEOUS_OR_UNSTABLE` and `METHOD_SENSITIVE` categories, with no single dominant continuous state.
4.  **`NO_CLEAR_CONTINUOUS_AXIS`**:
    - Spearman correlation between Primary contrast $C_i$ and Secondary singscore contrast $SS_{contrast, i}$ is weak ($\rho < 0.40$).
    - Ranks of samples are highly volatile and method-dependent (more than 40% of samples classified as `METHOD_SENSITIVE`).
    - The basal and classical program scores do not display inverse correlation (Spearman correlation between $S_{basal}$ and $S_{classical}$ is positive or near zero, $\rho \ge -0.2$).
5.  **`INCONCLUSIVE`**:
    - The statistical metrics violate the thresholds of all rules above, or the sensitivity analyses yield contradictory results (e.g., outlier exclusion or input scale shifts the trend test p-value from highly significant to non-significant).

---

## 10. Operational Status

This method lock is signed off for Phase 5B execution. No scripts will calculate scores or output tables until this lock is verified.
