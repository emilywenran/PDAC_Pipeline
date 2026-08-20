# Phase 6B Microbiome Preprocessing, Transformation, and Contamination-Sensitivity Method Lock (Amended)

This document locks the preprocessing pipeline, compositional transformations, contamination-sensitivity protocols, and downstream statistical modeling parameters for the tumor microbiome analysis of PRJNA719915. All specifications are predefined and locked prior to executing any differential-abundance or host-microbiome association tests.

## 1. Abundance-Scale Decision (Task 1)

The audited tumor microbiome matrix (`03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz`) contains:
- **Kraken2/Bracken-derived non-integer normalized counts** as published in Supplementary Data 1 (Genus-level sheet) of the source publication.
- It does **not** contain raw classified sequencing reads, and it does **not** contain relative abundances (percentages or fractions summing to 1).
- The values are positive continuous numbers representing Bracken-reallocated abundance estimates that have been normalized across samples by an unspecified formula (to be resolved under `TO_VERIFY`).

### Method Compatibility Constraints

Given this numerical scale, standard microbiome tools must be applied with care:
- **Valid Methods**:
  - **Centered Log-Ratio (CLR) Transform**: Valid after handling zero values. CLR transforms normalized counts into real space, resolving the compositionality constraint.
  - **Robust CLR (rCLR) Transform**: Valid and preferred as a sensitivity check. It avoids pseudocounts by calculating the geometric mean using only non-zero features.
  - **Spearman or Partial Spearman Correlation**: Valid for non-parametric association analyses on CLR or relative abundance scales.
  - **Multivariable Linear Modeling (OLS)**: Valid when applied to CLR-transformed values, as they occupy an unconstrained Euclidean space.
  - **Presence/Absence (Binarized) Analysis**: Valid for low-abundance or highly sparse features.
- **Invalid Methods**:
  - **ANCOM-BC2 (Standard Mode)**: Invalid on this matrix, as it strictly expects raw, unnormalized integer counts to model library-specific sampling fractions. Applying it to pre-normalized non-integer counts violates its statistical assumptions.
  - **ALDEx2**: Invalid directly on this matrix because it performs Dirichlet-multinomial sampling, which requires raw integer counts to model sampling uncertainty. Rounding these pre-normalized non-integer values is statistically invalid.
  - **DESeq2 / EdgeR**: Invalid because their statistical models assume raw, unnormalized integer counts and apply their own library-size normalization.
  - **Reconstruction of counts from relative abundance**: Prohibited. We must not reconstruct pseudo-counts and describe them as raw sequencing counts.

---

## 2. Prevalence and Abundance Filtering (Task 2)

Based on the Phase 6A audit and calculations, the baseline matrix zero fraction is **75.7%** (overall sparsity). We evaluate multiple filtering rules below.

### Evaluation of Candidate Rules

We evaluated combination rules of detection thresholds and prevalence thresholds. A taxon is "detected" in a sample if its abundance is strictly greater than the detection threshold ($T_{det}$). The candidate prevalence thresholds evaluated are 5%, 10%, 20%, and 30% of samples ($n=62$).

| Analysis ID | Detection Threshold ($T_{det}$) | Prevalence Threshold | Genera Retained | Genera Removed | Matrix Zero Fraction | Median Taxa per Sample | Samples with Low Signal ($<10$ taxa) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **MICRO_SENS_PREV_5** | $>0.0$ | $\ge 5\%$ ($\ge 4$ samples) | 187 | 178 | 54.22% | 80.0 | 0 |
| **MICRO_SENS_PREV_10** | $>0.0$ | $\ge 10\%$ ($\ge 7$ samples) | 149 | 216 | 44.56% | 78.0 | 0 |
| **MICRO_PRIMARY** | $>0.0$ | $\ge 20\%$ ($\ge 13$ samples) | 122 | 243 | 35.56% | 76.0 | 0 |
| **MICRO_SENS_PREV_30** | $>0.0$ | $\ge 30\%$ ($\ge 19$ samples) | 99 | 266 | 26.52% | 73.0 | 0 |
| **MICRO_SENS_DET_10_P20**| $>10.0$ | $\ge 20\%$ ($\ge 13$ samples) | 106 | 259 | 29.70% | 72.0 | 0 |

### Preprocessing Locks
1. **Primary Prevalence Rule**: Detected in **at least 20% of samples** (at abundance $> 0.0$). This retains **122 genera**, removes **243 genera**, reduces the matrix zero fraction to **35.56%**, and maintains a median of **76 detected taxa** per sample.
2. **Abundance Criterion**: For the primary analysis, the detection threshold is set to **abundance $> 0.0$**, as the minimum non-zero value in the matrix is `1.779` (meaning no sub-integer noise values exist in the audited matrix).
3. **Prevalence Sensitivity Rules**:
   - Lower threshold: **at least 10%** prevalence (retains 149 genera).
   - Higher threshold: **at least 30%** prevalence (retains 99 genera).
4. **Abundance Sensitivity Rule**: Detection defined as **abundance $> 10.0$ counts** combined with a **20% prevalence** threshold (retains 106 genera).
5. **Selection Independence**: Filtering thresholds are locked based purely on sparsity reduction and richness maintenance, without examining any associations with PDAC subtypes or continuous transcriptional scores.

---

## 3. Zero and Pseudocount Policy (Task 3)

Compositional log-ratio transformations cannot operate on zero values. We predefine a zero-replacement strategy that respects the scale of the Bracken-normalized abundance matrix.

### Primary Pseudocount Policy
- **Primary Method**: Add a **fixed small pseudocount** of **`0.889651`** to all cells in the filtered matrix before CLR transformation.
- **Rationale**: This value is exactly half of the minimum observed non-zero value in the matrix (`1.77930272379619`).
- **Cohort Specificity Warning**: This pseudocount is **highly specific to this source matrix** and is calculated based on its unique numerical distribution. It **must not be transferred unchanged to external cohorts or matrices** (which will have different sequencing depths, normalization scales, and minimum non-zero values). For any external dataset, the pseudocount must be recalculated using the corresponding matrix-specific minimum non-zero value.

### Pseudocount Sensitivity Policy
To verify that downstream findings are not artifacts of the chosen pseudocount, we lock three sensitivity policies:
1. **MICRO_SENS_PSEUDO_1.0**: A fixed pseudocount of `1.0`.
2. **MICRO_SENS_PSEUDO_0.1**: A very small fixed pseudocount of `0.1` (to assess behavior near zero).
3. **MICRO_SENS_ROBUST_CLR**: A pseudocount-free robust CLR transformation where the geometric mean is calculated using only non-zero values for each sample, and zeros are retained as zeros (or mapped to NaN in distance calculations).

---

## 4. Compositional Representations (Task 4)

We lock the specific mathematical representations and distances to be used in downstream analyses. Ordinary Euclidean distance on untransformed relative abundances is strictly prohibited as a primary distance metric due to compositionality artifacts.

### Primary Representations
- **Primary Feature Matrix**: Centered Log-Ratio (CLR)-transformed abundance matrix, defined as:
  $$CLR(x_i) = \ln\left(\frac{x_i}{g(x)}\right) = \ln(x_i) - \frac{1}{D}\sum_{j=1}^D \ln(x_j)$$
  where $x$ is the sample vector after prevalence filtering and pseudocount addition, and $g(x)$ is the geometric mean of the sample.
- **Primary Distance Metric**: **Aitchison distance**, which is equivalent to the Euclidean distance calculated on the CLR-transformed matrix. This metric is scale-invariant and subcompositionally coherence.

### Secondary or Sensitivity Representations
- **Robust CLR (rCLR)**: Log-ratio transformation calculated by centering log values around the geometric mean of observed (non-zero) features only, without adding a pseudocount.
- **Presence/Absence (Binarized)**: Binarized presence (1 if abundance $> 0.0$) or absence (0). This representation ignores abundance magnitude and captures occurrence patterns, evaluated using **Jaccard distance**.
- **Bray-Curtis Distance**: Restricted to descriptive/comparative visualizations only; it will not be used for primary statistical inference.
- **Count-based models**: Prohibited on this matrix because raw/integer counts are unavailable.

---

## 5. Contamination-Sensitivity Framework (Task 5)

No sequenced negative controls are available for the PRJNA719915 cohort. Thus, automated contaminant removal tools like `decontam` cannot be run, and environmental genera cannot be assumed to be contaminants without evidence. We establish a robust contamination-sensitivity framework.

### Contamination Risk Categories
We classify the 21 flagged environmental/reagent-associated genera from Phase 6A into evidence-based categories:

1. **HIGH_RISK_POTENTIAL_CONTAMINANT**: Common kit/reagent contaminants showing extremely high prevalence or anomalous abundance without biological rationale.
   - *Genera*: `Elizabethkingia` (extremely high abundance, accounting for 90%+ of reads in some samples; discussed in original study as potential database/genus mismatch), `Delftia`, `Brevundimonas`, `Comamonas`, `Caulobacter`, `Ralstonia`.
2. **MODERATE_RISK_ENVIRONMENTAL**: Genera commonly found in soil/water, but not typical human pathogens or symbionts.
   - *Genera*: `Paraburkholderia`, `Mesorhizobium`, `Novosphingobium`, `Dechloromonas`, `Sphingopyxis`, `Herbaspirillum`.
3. **BIOLOGICALLY_PLAUSIBLE_BUT_CONTAMINATION_SENSITIVE**: Flagged environmental/reagent genera that have also been reported in pancreatic tissue, biliary tract microbiota, or clinical infections, representing a mixture of potential biology and background noise.
   - *Genera*: `Pseudomonas` (known opportunistic pathogen, frequently reported in PDAC tumors), `Acinetobacter` (opportunistic pathogen, clinical contaminant), `Burkholderia`, `Stenotrophomonas`, `Sphingomonas`, `Rhizobium`, `Cupriavidus`, `Methylobacterium`, `Bradyrhizobium`.
4. **LOW_CURRENT_CONCERN**: Non-flagged genera with no significant contamination concerns in standard lists.
5. **TO_VERIFY**: Reserved for genera with ambiguous classification.

### Sensitivity Protocol for Contamination
To ensure biological findings are robust to potential contamination:
- **Analyses With and Without Flagged Genera**:
  - The primary analysis (`MICRO_PRIMARY`) will retain all filtered genera (including those flagged).
  - A sensitivity analysis (`MICRO_SENS_NO_CONTAMINANTS`) will exclude all genera classified as `HIGH_RISK_POTENTIAL_CONTAMINANT` and `MODERATE_RISK_ENVIRONMENTAL`.
- **Leave-One-Genus-Out (LOGO) Checks**: Any taxon found to be significantly associated with a host subtype or transcriptional axis must be subjected to a leave-one-genus-out check (i.e. recalculating sample CLR vectors and Aitchison distances without that specific genus to verify the remaining community structure or association is not dependent on that single genus).
- **Matrix Total-Abundance Proxy Correlation**:
  - We predefine an association check between each taxon's CLR abundance and the sample's **matrix total-abundance proxy** (sum of normalized values before filtering).
  - *Clarification on Abundance Proxy*: Because raw classified-read depth has not been independently verified from FASTQ metadata, the sum of normalized values serves strictly as a **matrix total-abundance proxy**, not a direct biological measurement of absolute microbial load.
  - Taxa showing strong positive correlation with this proxy ($\rho > 0.5$, $p < 0.01$) but low biological plausibility will be highlighted as potentially contamination-sensitive.
- **Reference Sources**: Contamination flags are based on literature lists of kit/reagent contaminants (Salter et al., 2014; Eisenhofer et al., 2019) and peer reviews of low-biomass tissue microbiomes.

---

## 6. Extreme-Sample Policy (Task 6)

Using only Phase 6A technical QC, we predefine objective sample exclusion criteria. No sample may be excluded based on its host molecular subtype, transcriptional axis score, or association with clinical outcomes.

### Objective Criteria Definitions
- **RETAIN**: Samples meeting all standard technical criteria:
  - Matrix total-abundance proxy $\ge 10,000$ Bracken-normalized counts.
  - Detected taxonomic richness $\ge 20$ genera.
  - Successfully mapped to a single patient with verified library layout.
- **RETAIN_WITH_SENSITIVITY_ANALYSIS**: Samples that pass basic sequencing quality checks but exhibit extreme values (outliers) on a single technical metric.
  - *Outlier definition*: Metric value $> 3$ standard deviations from the mean/median of log-scaled variables.
  - *Identified samples in this cohort*:
    - `Basal-like1`: Extreme high matrix total-abundance proxy ($1.42 \times 10^7$ counts, $z$-score = 5.01).
    - `Hybrid23`: Extreme high matrix total-abundance proxy ($9.02 \times 10^6$ counts, $z$-score = 4.58).
    - `Hybrid18`: Extreme high detected richness (198 genera, $z$-score = 3.23).
  - *Action*: These samples are **retained** in the primary analysis, but a sensitivity analysis (`MICRO_SENS_EXCLUDE_EXTREME`) must be run excluding them.
- **EXCLUDE_RECOMMENDED**: Excluded from all analyses due to technical failure. Requires violation of multiple technical metrics (e.g. low matrix total-abundance proxy AND low richness, or clear contamination/mapping failure).
  - *Action*: No samples in the current cohort meet this criterion.
- **TO_VERIFY**: Samples with metadata discrepancies or potential library format mismatches.

---

## 7. Downstream Method Compatibility (Task 7)

We prespecify the valid statistical frameworks to be used in Phase 7. These associations are **not** run during Phase 6B.

### A. Continuous Outcomes
Outcomes: `Moffitt50` basal–classical contrast, `coactivation score`, `PurIST probability`, `assignment entropy`.
- **Primary Method**: Multivariable linear regression (OLS) using CLR-transformed abundance as the predictor or outcome.
  - *Technical Sensitivity Covariate*: The **matrix total-abundance proxy** (log10 total abundance) must **not** be forced into every primary association model. Instead, it will be used selectively as a technical sensitivity covariate to assess the robustness of associations against library depth/normalized scale variation. It is explicitly declared as a proxy and not a direct microbial load measurement.
- **Non-parametric Correlation**: Spearman and partial Spearman correlation coefficients to assess monotonic relationships.
- **Multivariable Compositional Tool (MaAsLin2)**:
  - If MaAsLin2 is used to run associations on CLR-transformed values, it must be explicitly configured with:
    ```R
    normalization = "NONE"
    transform = "NONE"
    ```
  - This configuration is locked to **prevent any secondary normalization or logarithmic transformation** of the CLR values, which would invalidate their compositional mathematical properties.
- **Presence/Absence Sensitivity**: Logistic regression on binarized genus data to evaluate whether occurrence (rather than abundance) relates to host transcriptional features.

### B. Descriptive Discrete Outcome
Outcome: Public `Basal-like`, `Hybrid`, and `Classical` labels.
- **Primary Distance Association (PERMANOVA / PERMDISP)**:
  - We require the use of **PERMANOVA** on Aitchison distances with **9,999 permutations** to test for location differences between host subtypes.
  - To ensure that PERMANOVA results are not confounded by differences in within-group dispersion, a corresponding **PERMDISP** (homogeneity of multivariate dispersions / `betadisper`) assessment must be run.
  - Both location (PERMANOVA) and dispersion (PERMDISP) statistical results must be reported side-by-side.
- **Discrete Taxon Differences**: Wilcoxon rank-sum (binary contrasts) or Kruskal-Wallis (three-group comparison) on CLR values.

---

## 8. Multiple-Testing and Covariates (Task 8)

### Multiple-Testing Corrections
- **Procedure**: Benjamini-Hochberg (BH) False Discovery Rate (FDR) control.
- **Significance Threshold**: FDR-adjusted $p$-value $< 0.05$. Nominal $p$-values will be reported for exploratory purposes but flagged as such.
- **Effect-Size Reporting**: Spearman's $\rho$ (correlation) or standardized $\beta$ (regression) with 95% bootstrap confidence intervals.
- **Minimum Prevalence for Testing**: Taxa must be detected in at least 20% of samples to be tested for association, preventing excessive statistical testing on extremely sparse features.
- **Permutations & Bootstrap**: 10,000 permutations for non-parametric significance testing; 1,000 bootstrap resamples for confidence intervals.

### Covariate Hierarchy
To control for confounding, multivariable models will incorporate covariates in a hierarchical manner:
1. **Tier 1 (Technical Control)**:
   - **Matrix Total-Abundance Proxy**: Restructured as a technical sensitivity covariate (not forced in primary models; used to test robustness).
   - **Sequencing Batch**: Currently marked `TO_VERIFY` (unavailable in public SRA/ENA metadata).
2. **Tier 2 (Biological/Microenvironment Control)**:
   - **Tumor Purity**: Marked `TO_VERIFY` (not available in public metadata; must be computed from expression data in a later phase if required).
   - **Immune/Stromal Scores**: Marked `TO_VERIFY`.
3. **Tier 3 (Clinical Control)**:
   - **Age / Sex / Tumor Stage**: Marked `TO_VERIFY` (not present in public metadata).

---

## 9. Decision Rules for Phase 6B

We define the rules to establish overall cohort suitability before starting Phase 7 association analyses:

- **READY_FOR_COMPOSITIONAL_ANALYSIS**:
  - *Condition*: Non-negative abundance matrix with confirmed counts or estimated counts, low sample missingness, and sequenced negative controls available and processed to remove contaminant noise.
  - *Status*: **NOT MET** (due to lack of sequenced negative controls).
- **READY_WITH_CONTAMINATION_LIMITATIONS**:
  - *Condition*: Mapped and validated abundance matrix (62/62 samples), confirmed non-integer normalized counts, zero fraction reducible to $< 40\%$ by locked prevalence filters, but **sequenced negative controls are absent**.
  - *Status*: **MET**. The project is approved to proceed, provided the locked contamination-sensitivity framework (Tier-based risk flagging, leave-one-genus-out, with/without flagged analyses, and matrix total-abundance proxy technical covariates) is strictly applied.
- **INSUFFICIENT_MICROBIAL_SIGNAL**:
  - *Condition*: Post-filtering genus count $< 10$, or median richness per sample $< 10$ taxa.
  - *Status*: **NOT MET**.
- **INCOMPATIBLE_ABUNDANCE_SCALE**:
  - *Condition*: Abundance scale is relative abundance and cannot be converted to counts, or scale is unknown or remains `TO_VERIFY`.
  - *Status*: **NOT MET**.
- **INCONCLUSIVE**:
  - *Condition*: Unreconciled sample counts, duplicate runs, or lack of agreement on baseline QC metrics.
  - *Status*: **NOT MET**.

**Phase 6B Conclusion**: **READY_WITH_CONTAMINATION_LIMITATIONS**.
Preprocessing parameters and statistical methods are officially locked. Proceed to Phase 7 methods.
