# Standard Operating Procedure (SOP): PDAC Tumor Microbiome Preprocessing and Transformation Protocol (Amended)

## Objective
This protocol describes the step-by-step procedures for filtering, zero-replacement, compositional transformation, and contamination-sensitivity analysis of the genus-level tumor microbiome abundance matrix for the PRJNA719915 cohort. 

---

## Step 1: Input Data Verification

1. Load the audited genus-level abundance matrix:
   - File Path: `03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz`
   - Format: Gzipped TSV, rows as genera, columns as samples.
2. Confirm matrix dimensions: **365 genera x 62 samples**.
3. Verify that all values are non-negative ($\ge 0.0$) and that there are no missing (`NaN` or `None`) or infinite values.
4. Verify that the columns map 1-to-1 to the 62 patient identifiers in `01_metadata/sample_manifest.tsv` (`Basal-like1-17`, `Hybrid1-23`, `Classical1-22`).

---

## Step 2: Prevalence Filtering (Primary Rule)

1. Define "detected" as a genus having a normalized count **strictly greater than 0.0** (abundance $> 0.0$).
2. Compute the prevalence fraction for each genus:
   $$\text{Prevalence Fraction} = \frac{\text{Number of samples where Genus } > 0.0}{62}$$
3. Filter out all genera with a prevalence fraction **less than 20%** (i.e. keep only genera detected in $\ge 13$ samples).
4. Programmatically in Python:
   ```python
   import pandas as pd
   df = pd.read_csv("03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz", sep="\t", index_col=0)
   detected_mask = df > 0.0
   prevalence = detected_mask.sum(axis=1) / df.shape[1]
   df_filtered = df.loc[prevalence >= 0.20]
   # Verify: 122 genera retained, 243 removed.
   ```

---

## Step 3: Zero-Replacement (Pseudocount Policy)

1. Retrieve the minimum observed non-zero value across the entire matrix. This value is verified as **`1.77930272379619`**.
2. Calculate the primary pseudocount as half of this minimum observed non-zero value:
   $$\text{Pseudocount} = \frac{1.77930272379619}{2} \approx 0.889651$$
3. Add the pseudocount to all cells of the prevalence-filtered matrix:
   ```python
   df_pseudo = df_filtered + 0.889651361898095
   ```
4. **Cohort-Specificity Warning**: The primary pseudocount of `0.889651` is **strictly specific to this source matrix** and is calculated based on its unique non-zero distribution. It **must not be transferred unchanged to external cohorts or datasets**. For external matrices, the pseudocount must be recalculated as half of their respective minimum observed non-zero value.

---

## Step 4: Centered Log-Ratio (CLR) Transformation

1. Compute the centered log-ratio (CLR) transformation on the pseudocount-adjusted matrix.
2. For each sample column:
   - Compute the log of each feature value.
   - Calculate the arithmetic mean of these log values.
   - Subtract the mean log value from each log value.
3. Programmatic implementation:
   ```python
   import numpy as np
   log_df = np.log(df_pseudo)
   mean_log = log_df.mean(axis=0)
   clr_df = log_df.sub(mean_log, axis=1)
   ```
4. Verify that the column sums of `clr_df` are zero (within floating-point tolerance).

---

## Step 5: Distance Matrix Calculation

1. Calculate the Aitchison distance matrix, which is the pairwise Euclidean distance between samples in the CLR space.
2. Programmatic implementation:
   ```python
   from scipy.spatial.distance import pdist, squareform
   distances = pdist(clr_df.T, metric='euclidean')
   aitchison_dist_matrix = pd.DataFrame(
       squareform(distances),
       index=clr_df.columns,
       columns=clr_df.columns
   )
   ```

---

## Step 6: Contamination-Sensitivity Analysis

Because sequenced negative controls are absent, the primary statistical findings must be validated against potential contamination using three locked sensitivity steps:

1. **Category Flagging**: Apply the evidence-based classifications:
   - **High Risk**: `Elizabethkingia`, `Delftia`, `Brevundimonas`, `Comamonas`, `Caulobacter`, `Ralstonia`.
   - **Moderate Risk**: `Paraburkholderia`, `Mesorhizobium`, `Novosphingobium`, `Dechloromonas`, `Sphingopyxis`, `Herbaspirillum`.
2. **Analysis ID `MICRO_SENS_NO_CONTAMINANTS`**:
   - Exclude all 12 High Risk and Moderate Risk genera from the audited matrix before running the filtering, pseudocount, and transformation steps.
   - Re-run all downstream association tests and confirm that the primary findings remain statistically significant ($p_{adj} < 0.05$).
3. **Leave-One-Genus-Out (LOGO) Check**:
   - For any genus showing significant association, exclude that specific genus and re-compute CLR and Aitchison distances. Re-run associations to ensure that the correlation is not an artifact of that single genus driving the overall composition.
4. **Depth Correlation Control**:
   - Calculate the Spearman correlation between each taxon's CLR abundance and the sample's **matrix total-abundance proxy** (sum of normalized counts before filtering). 
   - *Technical Proxy Warning*: This sum is a **matrix total-abundance proxy**, not a direct biological measurement of absolute microbial load, as raw classified-read depth has not been independently verified.
   - Flag any taxon showing strong correlation ($\rho > 0.5$, $p < 0.01$) for manual inspection.

---

## Step 7: Preprocessing & Transform Sensitivity Audits

To ensure robust results, verify the main findings under these alternative parameters:

### A. Prevalence Thresholds
- **Lower Threshold (10%)**: Run `MICRO_SENS_PREV_10` (retaining 149 genera).
- **Higher Threshold (30%)**: Run `MICRO_SENS_PREV_30` (retaining 99 genera).

### B. Pseudocount Adjustments
- **Fixed 1.0**: Run `MICRO_SENS_PSEUDO_1.0` (adding 1.0 instead of 0.889651).
- **Fixed 0.1**: Run `MICRO_SENS_PSEUDO_0.1` (adding 0.1).
- **Robust CLR (rCLR)**: Run `MICRO_SENS_ROBUST_CLR` (log-transforming non-zero values and centering around their geometric mean, leaving zeros as zero/NaN).

### C. Representation Sensitivity
- **Presence/Absence**: Run `MICRO_SENS_PRESENCE_ABSENCE` using Jaccard distance on binarized data:
  ```python
  binarized_df = (df_filtered > 0.0).astype(int)
  # Calculate Jaccard distance using scipy
  ```

---

## Step 8: Extreme Sample Outlier Validation

1. Identify the three samples flagged for single-metric extreme outliers:
   - `Basal-like1` (extreme matrix total-abundance proxy, $z = 5.01$)
   - `Hybrid23` (extreme matrix total-abundance proxy, $z = 4.58$)
   - `Hybrid18` (extreme richness, $z = 3.23$)
2. Run `MICRO_SENS_EXCLUDE_EXTREME` by removing these three samples from the filtered matrix before CLR transformation.
3. Re-run downstream associations and verify that results are not driven by these extreme technical outliers.

---

## Step 9: Downstream Association Execution Rules

For downstream associations (Phase 7), run the statistical models with the following locked parameters:

### A. Multivariable Linear Models (MaAsLin2)
If MaAsLin2 is used to run associations on CLR-transformed values, you **MUST** configure the run to prevent secondary transformations or normalizations:
```R
library(Maaslin2)
fit_data = Maaslin2(
    input_data = clr_abundance_matrix,
    input_metadata = metadata,
    output = "output_directory",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("covariate1", "covariate2"),
    # ...
)
```
- **Covariate Rule**: The **matrix total-abundance proxy** (log10 total abundance) must **not** be included in every primary model. It must be evaluated only as a technical sensitivity covariate to assess the impact of sequencing depth or normalized abundance scale variability.

### B. Distance Associations (PERMANOVA / PERMDISP)
When evaluating discrete subtype differences on Aitchison distance matrices:
1. Replicate the test using **PERMANOVA** with **9,999 permutations**.
2. Run a corresponding **PERMDISP** (`betadisper` in R or `permdisp` in python/skbio) to check for differences in within-group dispersion.
3. Report both the location shift significance (PERMANOVA p-value) and the dispersion homogeneity significance (PERMDISP p-value) side-by-side to ensure location effects are not confounded by dispersion differences.
