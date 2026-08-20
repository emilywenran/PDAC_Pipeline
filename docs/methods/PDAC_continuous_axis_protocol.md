# PDAC Continuous Basal–Classical Axis Protocol

## 1. Scope and Purpose

This protocol defines the step-by-step procedures, computational environment, package requirements, script structure, and output schemas for executing the continuous basal–classical transcriptional-axis characterization of the GSE172356 cohort. 

This protocol serves as the execution guide for Phase 5B. In accordance with Phase 5A requirements, no continuous scoring or downstream statistical evaluations are run in this phase.

---

## 2. Input Data and Checksums

The execution phase is restricted to the following verified input matrices, signatures, and prior metrics. No alternative inputs may be introduced.

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

## 3. Computational Environment and Packages

The analysis will be implemented in R (v4.3+) or Python (v3.10+) within the established workspace environment.

### A. R Packages and Versioning
- **`singscore` (v1.20+)**: For secondary rank-based scoring.
- **`diptest` (v0.77+)**: For Hartigan's Dip Test of multimodality.
- **`clinfun` (v1.1+)**: For Jonckheere-Terpstra trend testing.
- **`boot` (v1.3+)**: For bootstrap confidence interval estimation.
- **`ggplot2` (v4.0+)** & **`ggpubr` (v0.6+)**: For density, correlation, and scatter visualizations.

### B. Directory Layout for Outputs
All outputs from the execution script must be written to:
- Reports: [04_analysis/07_continuous_subtype_axis/](file://~/thesis/PDAC/04_analysis/07_continuous_subtype_axis/)
- Data Tables: [05_results/tables/](file://~/thesis/PDAC/05_results/tables/) (prefixed with `phase5b_`)
- Visualizations: [05_results/figures/](file://~/thesis/PDAC/05_results/figures/) (prefixed with `phase5b_`)

---

## 4. Preprocessing and Scoring Workflows

### A. Primary Score Workflow
1. Load `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz`.
2. Subset expression matrix to the active Moffitt signature genes:
   - Primary: 50 genes (25 basal-like, 25 classical, including LEMD1).
   - Sensitivity (LEMD1-Excluded): 49 genes (24 basal-like, 25 classical, excluding LEMD1).
3. Compute row Z-scores for each gene:
   $$z_{g, i} = \frac{x_{g, i} - \mu_g}{\sigma_g}$$
4. Calculate primary scores for each sample $i$:
   - **Basal score**: $S_{basal, i} = \text{mean}(z_{g, i})$ for $g \in G_{basal}$.
     - $|G_{basal}| = 25$ for Primary.
     - $|G_{basal}| = 24$ for LEMD1-Excluded Sensitivity.
   - **Classical score**: $S_{classical, i} = \text{mean}(z_{g, i})$ for $g \in G_{classical}$.
     - $|G_{classical}| = 25$ for both Primary and Sensitivity.
   - **Contrast score**: $C_i = S_{basal, i} - S_{classical, i}$.
5. Calculate co-activation score $A_i$:
   - Normalize $S_{basal}$ and $S_{classical}$ across all samples to a $[0, 1]$ scale ($S'_{basal}, S'_{classical}$).
   - Set $A_i = \min(S'_{basal, i}, S'_{classical, i})$.

### B. Secondary Score Workflow (Singscore)
1. Load `03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz`.
2. Construct singscore GeneSet objects for:
   - Primary: 25 basal and 25 classical genes.
   - Sensitivity (LEMD1-Excluded): 24 basal and 25 classical genes.
3. Compute single-sample rank-based scores:
   - Run `singscore::simpleScore` on the normalized count matrix with the basal and classical GeneSets separately to obtain $SS_{basal}$ and $SS_{classical}$.
   - Calculate singscore contrast: $SS_{contrast, i} = SS_{basal, i} - SS_{classical, i}$.

### C. Centroid Distance Workflow
1. For reference-anchored centroids:
   - Identify samples publicly labeled as `Basal` and `Classical` (from `phase3b_primary_subtype_assignments.tsv`, excluding `Hybrid`).
   - Compute mean Z-score vectors $\mathbf{\mu}_{basal}$ and $\mathbf{\mu}_{classical}$ in the active space (50-gene space for Primary, 49-gene space for Sensitivity).
2. For unsupervised centroids:
   - Perform hierarchical clustering ($K=2$) on the active Z-score matrix of all samples (50-gene for Primary, 49-gene for Sensitivity).
   - Designate the cluster with higher average basal marker Z-scores as the unsupervised basal centroid $\mathbf{\mu}_{C1}$, and the other as $\mathbf{\mu}_{C2}$.
3. For each sample $i$:
   - Calculate Euclidean distance to both centroid vectors: $d_{basal, i}$ and $d_{classical, i}$.
   - Calculate centroid-to-centroid distance $d_{bc} = d(\mathbf{\mu}_{basal}, \mathbf{\mu}_{classical})$.

---

## 5. Statistical Execution Details

### A. Jonckheere-Terpstra Trend Test
1. Set sample groups as `Classical`, `Hybrid`, and `Basal` based on public labels.
2. Order groups as: $\text{Classical} < \text{Hybrid} < \text{Basal}$.
3. Run `clinfun::jonckheere.test` on the continuous scores ($C_i$ and $SS_{contrast, i}$) with the ordered factor.
4. Record the standardized statistic $J$ and two-tailed p-value.

### B. Hartigan's Dip Test of Multimodality
1. Run `diptest::dip.test` on the continuous contrast score vector $C$.
2. Record the dip statistic $D$ and p-value.
3. Repeat for $SS_{contrast}$.

### C. Bootstrap Confidence Intervals
1. Set random seed `2026`.
2. Define a function to draw a bootstrap sample of size $N=62$ with replacement.
3. For each bootstrap iteration:
   - Re-compute row Z-scores and centroids.
   - Re-calculate primary continuous scores.
   - Re-estimate the correlation coefficients between scoring systems and Hedges' $g$ difference between poles.
4. Run for 1,000 iterations. Compute the 95% CI using the percentile method.

---

## 6. Output File Schemas

The Phase 5B script must produce the following files with these exact schemas:

### Table 1: `phase5b_continuous_axis_scores.tsv`
Stores the continuous scoring vectors for all 62 samples:
- `sample_id`: character, GEO accession (e.g., `YX16021T`)
- `public_label`: character, subtype (`Basal`, `Classical`, `Hybrid`)
- `primary_basal_score`: numeric
- `primary_classical_score`: numeric
- `primary_contrast_score`: numeric
- `coactivation_score`: numeric
- `singscore_basal`: numeric
- `singscore_classical`: numeric
- `singscore_contrast`: numeric
- `purist_basal_prob`: numeric
- `moffitt_score_diff`: numeric
- `distance_basal_ref`: numeric
- `distance_classical_ref`: numeric
- `distance_basal_unsup`: numeric
- `distance_classical_unsup`: numeric
- `midpoint_proximity`: numeric
- `distance_to_both_poles`: numeric
- `method_axis_variance`: numeric
- `axis_interpretation_category`: character (`BASAL_POLE`, `CLASSICAL_POLE`, `INTERMEDIATE_CONTINUUM`, `COACTIVATED_HYBRID`, `HETEROGENEOUS_OR_UNSTABLE`, `METHOD_SENSITIVE`, `TO_VERIFY`)

### Table 2: `phase5b_hybrid_metric_summary.tsv`
Summarizes continuous and stability metrics specifically for the public Hybrid samples ($n=23$):
- `metric`: character (e.g., `coactivation_score`, `midpoint_proximity`, `distance_to_both_poles`, `assignment_entropy`, `item_consensus`)
- `median`: numeric
- `mean`: numeric
- `sd`: numeric
- `ci_lower`: numeric
- `ci_upper`: numeric

### Table 3: `phase5b_statistical_evaluations.tsv`
Stores the results of the prespecified statistical checks:
- `test_id`: character (e.g., `JT_trend_primary`, `Dip_test_primary`, `Spearman_corr_1_2`)
- `test_name`: character
- `statistic`: numeric
- `p_value`: numeric
- `effect_size`: numeric
- `ci_lower`: numeric
- `ci_upper`: numeric
- `notes`: character

### Table 4: `phase5b_outlier_sensitivity.tsv`
Logs the effect of excluding the four outlier candidates on scores and centroids:
- `evaluation_metric`: character (e.g., `centroid_shift_basal`, `score_correlation_pearson`)
- `value`: numeric
- `notes`: character

---

## 7. Visualization Requirements

The following four figures must be generated as PDF and saved under `05_results/figures/`:

1.  **`phase5b_basal_classical_scatter.pdf`**:
    - Scatter plot of $S_{basal}$ (y-axis) vs $S_{classical}$ (x-axis) for all samples.
    - Points colored by public subtype (`Classical` = blue, `Hybrid` = green, `Basal` = red).
    - Points styled with different shapes based on their Phase 4B stability interpretation.
2.  **`phase5b_contrast_boxplot.pdf`**:
    - Boxplot of $C_i$ (y-axis) grouped by public subtype (x-axis).
    - Jittered sample points overlaid on the boxplot.
3.  **`phase5b_contrast_density.pdf`**:
    - Density curve of the continuous contrast score $C_i$ across the cohort.
    - Overlay of Hartigan's Dip Test statistic and p-value.
4.  **`phase5b_method_correlation.pdf`**:
    - Heatmap or scatter matrix displaying Spearman correlation coefficients between all continuous scoring systems, with significance stars.
