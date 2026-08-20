# PDAC Subtype Stability Evaluation Protocol

## 1. Scope and Purpose

This protocol defines the step-by-step procedures, computational environment, and metrics calculation guidelines for evaluating the structural stability and reproducibility of PDAC transcriptomic subtypes (Basal-like, Classical, and Hybrid) in the GSE172356 cohort. 

This protocol serves as a guide for executing the Phase 4B stability run. In accordance with Phase 4A requirements, no actual stability calculations are performed in this phase.

---

## 2. Input Data and Checksums

The stability analysis is restricted to the following verified input matrices and signatures. No alternative inputs may be introduced.

| Input Name | File Path | Expected SHA256 Checksum |
| :--- | :--- | :--- |
| **Filtered Normalized counts** | `03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz` | `8947ca75c3240177f8daeb8426e4cc9978a94c51ed17b14cb6eaf0146c4d73c1` |
| **Filtered Log2 counts** | `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz` | `13c16a95c7ef94e59b7d685c85b78f4bc2a2d22b9e6ffaafb929dd2a50c0328a` |
| **Original 94-gene signature** | `02_data/reference/PDAC_subtype_signatures/GSE172356_original_signatures.tsv` | `6288537c06abe45db5251ebcc22dfa9a2df944176652a7a298d391db29664175` |
| **Moffitt 2015 signature** | `02_data/reference/PDAC_subtype_signatures/Moffitt_2015_signatures.tsv` | `29468d394dffe32b08b0959696dda2ac8b9a10f1eef6ff5c1655587638b9afa0` |
| **PurIST signature & coeffs** | `02_data/reference/PDAC_subtype_signatures/PurIST_signatures.tsv` | `b198e583e65c8e4f1da04e2054c24c23d201c7cada330535fe3f3645a11d249f` |

---

## 3. Computational Environment

The analysis will be implemented using R (v4.3+) or Python (v3.10+) within the established workspace environment.

### A. R Packages and Versioning
- **`ConsensusClusterPlus` (v1.64+)**: For consensus clustering matrix generation and CDF metrics.
- **`cluster` (v2.1+)**: For silhouette width calculation.
- **`fpc` (v2.2+)**: For bootstrap Jaccard cluster stability and prediction strength evaluations.
- **`clue` (v0.3+)**: For Hungarian label-alignment algorithm implementation.

### B. Directory Layout for Outputs
Stability results must be written to:
- Reports: `04_analysis/06_subtype_stability/`
- Data Tables: `05_results/tables/` (prefixed with `phase4b_`)
- Visualizations: `05_results/figures/` (prefixed with `phase4b_`)

---

## 4. Preprocessing Workflows

### Workflow A: CSY Primary representation
1. Load `03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz`.
2. Subset matrix to the 94 CSY signature genes.
3. Compute the row median for each gene across the samples.
4. Center the rows by subtracting their respective medians.
5. Scale the centered rows by subtracting the row mean and dividing by the row standard deviation (default behavior of `pheatmap` row scaling).

### Workflow B: Log2 count representation
1. Load `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz`.
2. Subset matrix to the 94 CSY signature genes.
3. Center rows by subtracting row medians, then scale row-wise.

### Workflow C: Unsupervised Representation (HVGs)
1. Load `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz`.
2. Calculate Median Absolute Deviation (MAD) for all 42,654 genes:
   $$\text{MAD}_g = \text{median}(|x_{g,i} - \text{median}(x_g)|)$$
3. Rank genes by MAD and select the top 1000 highest-variance genes.
4. Scale the top 1000 genes by Z-score transformation:
   $$z_{g,i} = \frac{x_{g,i} - \mu_g}{\sigma_g}$$

---

## 5. Unsupervised Clustering Execution

### A. Consensus Clustering (CCP)
- **Sample Sampling Fraction ($p_{sample}$)**: 0.80 (sampling 49 out of 62 samples without replacement).
- **Feature Sampling Fraction ($p_{feature}$)**: 1.00 (except for Run ID `STAB_CSY_FEAT_RESAMP` where $p_{feature} = 0.80$).
- **Iterations ($B$)**: 1,000 runs.
- **Distance Metric**: 
  - Primary runs: Pearson correlation distance ($1 - r$).
  - Secondary runs: Euclidean distance.
- **Clustering Method**:
  - Primary runs: Hierarchical clustering with average linkage.
  - Secondary runs: Hierarchical clustering with Ward's linkage (`ward.D2`).
- **Random Seed**: `2026` initialized before running the CCP loop.

### B. Full-Dataset Reference Clustering
- Run the clustering algorithm on the complete cohort (62 or 58 samples) to obtain the reference partition for each candidate $K \in \{2..6\}$.

### C. Post-Clustering Label Alignment
- To align cluster assignments from arbitrary $K$ partitions to public labels for comparison:
  1. Construct a $K \times K$ confusion matrix between the unsupervised clusters and reference labels (e.g., Basal, Hybrid, Classical).
  2. Implement the **Hungarian Algorithm** to find the bijection that maximizes the trace of the confusion matrix (maximum-overlap matching).
  3. Re-label cluster numbers to the matching reference label names.

---

## 6. Stability Metrics Calculation

### A. Consensus Cumulative Distribution Function (CDF)
For the consensus matrix $M$ where $M(i,j)$ is the co-clustering fraction of samples $i$ and $j$ across resampling iterations:
$$H(x) = \frac{\sum_{i < j} I(M(i,j) \le x)}{\frac{N(N-1)}{2}}$$
Evaluate $H(x)$ over $x \in [0, 1]$ to detect the transition slope.

### B. Proportion of Ambiguous Clustering (PAC)
Calculate the fraction of consensus values in the intermediate interval:
$$\text{PAC} = H(u_2) - H(u_1)$$
where $u_1 = 0.1$ and $u_2 = 0.9$. Lower values of PAC indicate a more stable discrete partition.

### C. Within-Cluster Consensus ($I_k$)
For cluster $k$ containing a subset of samples $C_k$ of size $N_k$:
$$I_k = \frac{\sum_{i, j \in C_k, i < j} M(i,j)}{\frac{N_k(N_k - 1)}{2}}$$

### D. Sample Assignment Entropy ($H_i$)
For each sample $i$, compute the entropy across all resampling runs in which sample $i$ was drawn:
$$H_i = - \sum_{k=1}^K p_{ik} \log_2(p_{ik})$$
where $p_{ik}$ is the probability that sample $i$ is assigned to cluster $k$.

### E. Bootstrap Cluster Jaccard Stability
1. Perform 1,000 bootstrap resamples (sampling with replacement) of size $N$.
2. Run the clustering algorithm on each bootstrap dataset.
3. For each cluster $C_k$ in the full-dataset partition, calculate the maximum Jaccard coefficient with any cluster $C'_j$ in the bootstrap partition:
   $$\text{Jaccard}(C_k, C'_j) = \frac{|C_k \cap C'_j|}{|C_k \cup C'_j|}$$
4. Report the mean Jaccard coefficient across all 1,000 bootstrap runs.

---

## 7. Sample-Level Hybrid Evaluation Metrics

To characterize the intermediate or hybrid state of each sample, compute:

1. **Centroid Distances**:
   - Compute the Basal and Classical cluster centroids in the 94-gene scaled space.
   - Calculate the Euclidean distance of each sample to both centroids: $d(i, \text{Basal})$ and $d(i, \text{Classical})$.
2. **Continuous Signature Scores**:
   - **Moffitt Basal Score**: Mean row-scaled expression of the 25 basal genes.
   - **Moffitt Classical Score**: Mean row-scaled expression of the 24 classical genes.
   - **PurIST Basal Probability**: Single-sample probability $p$ from the TSP logistic model.
3. **Item Consensus ($I_i$)**:
   - The average consensus index between sample $i$ and all other members of its assigned cluster:
     $$I_i(k) = \frac{\sum_{j \in C_k, j \ne i} M(i,j)}{N_k - 1}$$

---

## 8. Output File Schemas

### Table 1: `phase4b_cohort_stability_metrics.tsv`
Stores cohort-level statistics for each candidate K across all 8 configurations:
- `analysis_id`, `candidate_K`, `overall_mean_silhouette`, `PAC`, `mean_Jaccard_stability`, `prediction_strength`, `mean_ARI`, `mean_entropy`.

### Table 2: `phase4b_cluster_specific_stability.tsv`
Stores cluster-level statistics for each cluster under candidate K:
- `analysis_id`, `candidate_K`, `cluster_id`, `cluster_label`, `cluster_size`, `within_cluster_consensus_Ik`, `cluster_mean_silhouette`, `cluster_mean_Jaccard`.

### Table 3: `phase4b_sample_level_stability.tsv`
Stores sample-specific stability and assignment metrics for the $K=3$ run:
- `sample_id`, `assigned_cluster`, `item_consensus`, `assignment_entropy`, `silhouette_width`, `purist_basal_prob`, `moffitt_basal_score`, `moffitt_classical_score`, `dist_to_basal_centroid`, `dist_to_classical_centroid`, `stability_interpretation`.
