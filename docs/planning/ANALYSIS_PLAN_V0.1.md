# Analysis Plan V0.1

This document outlines the step-by-step computational workflow for auditing data, profiling expression and microbiome samples, evaluating subtype stability, integrating host-microbe interactions, and identifying therapeutic targets.

---

## 1. Data Audit & Metadata Mapping

### 1.1 Primary Datasets
*   **Host Transcriptomic Data (NCBI GEO: GSE172356):**
    *   **Platform:** Illumina HiSeq X Ten (GPL20795)
    *   **Organism:** *Homo sapiens*
    *   **Target Sample Size:** $n=62$ resectable PDAC primary tumors.
    *   **Data Types Expected:** Processed expression matrices (counts/FPKM/TPM) as GEO supplementary files, and raw reads in SRA (`TO VERIFY` exact format of GEO supplementary files).
    *   **Clinical Covariates Needed:** Age, sex, tumor stage, survival time, recurrence status, resection margin status (`TO VERIFY` completeness of clinical tables).
*   **Tumor Microbiome Data (NCBI BioProject: PRJNA719915):**
    *   **Platform:** Illumina (`TO VERIFY` exact machine model, e.g., NovaSeq).
    *   **Data Type:** Shotgun metagenomic tumor sequencing data.
    *   **Control Samples:** Negative controls / low-biomass contamination assessment samples (`TO VERIFY` SRA Run accessions for these controls).

### 1.2 External Validation & Prioritization Databases
*   **External Transcriptomic Cohorts:**
    *   **TCGA-PAAD:** Bulk RNA-seq of $n=178$ pancreatic cancer samples.
    *   **ICGC PACA-AU / PACA-CA:** Bulk RNA-seq/microarray datasets.
    *   **Single-Cell RNA-seq:** Peng et al. (2019) or similar public PDAC dataset (`TO VERIFY` accessions).
    *   **Spatial Transcriptomics:** Moncada et al. (2020) or similar public spatial PDAC dataset (`TO VERIFY` accessions).
*   **Target Prioritization Databases:**
    *   **DepMap / Achilles:** CRISPR knockout screens across $n \approx 40$ pancreatic cancer cell lines (Chronos dependency scores).
    *   **PRISM:** Repurposing drug screen dataset across pancreatic cancer cell lines.
    *   **Druggability:** DGIdb (Drug-Gene Interaction Database) and ChEMBL.

---

## 2. Sample Mapping & Verification

### Protocol
1.  Download GEO metadata for GSE172356 and BioProject metadata for PRJNA719915.
2.  Extract sample mapping identifiers (e.g., patient IDs, SRA Run IDs, and GEO GSM IDs).
3.  Construct a mapping master table (`01_metadata/sample_mapping.tsv`) containing:
    *   `Patient_ID` (unique identifier)
    *   `GSM_Expression_ID`
    *   `SRR_Expression_ID` (if raw data downloaded)
    *   `SRR_Microbiome_ID`
    *   `Sample_Type` (Tumor, Adjacent Normal, Kit Blank, PCR Blank, Environmental Control)
4.  **Verification Check:** Assert that $n=62$ matches exist for paired host-microbiome samples. Flag any sample lacking one of the profiles.

---

## 3. Host Expression Quality Control & Deconvolution

### 3.1 Expression QC
1.  **Read Metrics:** If utilizing raw FASTQ files, run FastQC/MultiQC. Align with STAR to GRCh38.
2.  **Count Matrix Audit:** Check for low-expression genes. Filter genes with $< 10$ reads in $> 80\%$ of samples.
3.  **Outlier Detection:** Perform log2(CPM+1) transformation or variance stabilizing transformation (VST) using DESeq2. Run Principal Component Analysis (PCA) and hierarchical clustering to identify transcriptomic outliers.
4.  **Batch Effect Assessment:** Evaluate if sample sequencing batch or run dates correlate with major PCs using guided PCA or the PVCA (Principal Variance Component Analysis) R package. If batch effects are detected, correct using `ComBat` or include batch as a covariate in downstream linear models.

### 3.2 Cell-Type Composition & Tumor Purity Estimation
Because bulk RNA-seq contains stromal, immune, and malignant cells, we must estimate tumor purity and cell-type fractions to prevent stromal and immune cell signals from confounding the subtype and microbiome associations:
1.  **Tumor Purity:** Run the `ESTIMATE` R package to calculate tumor purity (tumoral, stromal, and immune scores).
2.  **Cell-Type Deconvolution:** Apply `ConsensusTME` or `xCell` to estimate fractions of specific stromal cell types (fibroblasts, endothelial cells) and immune cell types (macrophages, T-cells, B-cells).
3.  Save the deconvolution output to `03_processed/expression/tumor_purity_deconvolution.tsv`.

---

## 4. Microbiome Shotgun Metagenomic Workflow

### 4.1 Raw-Read QC & Host-Read Depletion
1.  **Raw-Read QC:** Run FastQC/MultiQC to evaluate base qualities, adapter contamination, and sequence duplicates. Use Trimmomatic or fastp to trim adapters and filter out low-quality reads (e.g., Q < 20, length < 50bp).
2.  **Host-Read Depletion:** Align reads to the human reference genome (GRCh38) using Bowtie2 or BWA-MEM. Retrieve unmapped (non-human microbial) reads for downstream profiling using samtools. Save depletion statistics (percent human reads) to evaluate potential host contamination and batch effects.

### 4.2 Taxonomic & Functional Profiling
1.  **Taxonomic Profiling:** Profile non-human reads using Kraken2 paired with Bracken (for species-level abundance estimation) and/or MetaPhlAn (marker-gene database).
2.  **Sensitivity Analysis across Taxonomic Methods:** Perform sensitivity analysis by comparing the taxonomic results from Kraken2/Bracken vs. MetaPhlAn to check for consistency in identified taxa and overall community composition.
3.  **Microbial Functional Profiling:** Where sequencing depth supports it, perform functional profiling using HUMAnN (or a similar tool) to estimate metabolic pathway abundance and coverage (e.g., UniRef gene families mapped to MetaCyc pathways).

### 4.3 Contamination Assessment & Filtering (Low-Biomass Safeguards)
Since tumor tissue is a low-biomass environment, we must guard against reagent and laboratory contamination:
1.  **Negative-Control and Low-Biomass Contamination Assessment:** Utilize negative controls (extraction blanks, PCR controls, or other available blanks) to identify contaminant features using the `decontam` R package (or equivalent approach) via:
    *   **Prevalence Method:** Compare the prevalence of microbial features in tumor samples vs. negative controls.
    *   **Frequency Method:** If DNA concentration metadata is available for PRJNA719915, test for inverse correlation between feature abundance and measured DNA concentration.
2.  **Prevalence and Abundance Filtering:** Apply filters to retain only robust microbial features:
    *   Remove identified contaminant features.
    *   Filter out low-abundance features (e.g., threshold of $>0.01\%$ relative abundance) and low-prevalence features (e.g., present in $<5\%$ of samples).

### 4.4 Compositional Data Analysis
1.  Shotgun metagenomic profiles are compositional (subject to varying sequencing depth/arbitrary total sum constraints).
2.  Do not use simple rarefaction for downstream statistical modeling, as it discards reads and reduces statistical power.
3.  Apply composition-aware transformations such as **Centered Log-Ratio (CLR)** (using a pseudocount of 1 or multiplicative replacement) for PCA, UMAP, and Euclidean distance-based analyses.
4.  Utilize composition-aware regression and differential abundance tools (e.g., **ANCOM-BC** or **ALDEx2**) for testing differences across subtypes or covariates.

---

## 5. Subtype Reproduction & Stability Analysis

### 5.1 Subtype Reproduction
We will attempt to replicate the three-class clustering (basal-like, classical, hybrid) of Guo et al. (2021):
1.  Obtain the gene expression matrix for GSE172356. Apply log2(CPM+1) or VST normalization.
2.  Retrieve the top 25 ranked genes from the four malignancy-related NMF components (Components 1, 2, 6, and 10 of Moffitt et al., 2015).
    *   *Moffitt Basal-like Components:* Components 6 & 10 (`TO VERIFY` gene lists).
    *   *Moffitt Classical Components:* Components 1 & 6 or 1 & 2 (`TO VERIFY` gene lists).
3.  Extract these signature genes from the GSE172356 matrix.
4.  Run the **ConsensusClusterPlus** R package:
    *   Parameters: `maxK = 6`, `reps = 10000`, `pItem = 0.8`, `pFeature = 1`, `clusterAlg = "hc"` (hierarchical clustering), `distance = "pearson"`.
5.  Assign samples to $k=3$ groups based on consensus clustering. Compare assignments with the classification reported in Guo et al. using a confusion matrix and Rand Index / Cohen's Kappa.

### 5.2 Subtype Stability Analysis
Evaluate if the three-class structure is stable or represents clustering of a continuous gradient:
1.  **Silhouette Width:** Calculate the silhouette width ($s_i$) for each sample in the $k=3$ partition.
    *   Evaluate if the mean silhouette width of the "hybrid" group is significantly lower than the basal-like and classical groups, which would suggest it lacks distinct boundaries.
2.  **Area Under CDF Curve:** Plot the cumulative distribution function (CDF) for consensus matrix values across $k=2$ to $k=6$. Inspect the CDF slope. A lack of flat regions in the CDF suggests the absence of clear natural clusters.
3.  **Bootstrapping Stability:** Perform $1000$ bootstrap iterations. Calculate the Jaccard similarity coefficient for each cluster. Stable clusters should exhibit median Jaccard $> 0.85$.
4.  **Comparison to Random Gene Sets:** Run the same consensus clustering on 1000 random gene sets of size 100. Calculate the proportion of random gene sets that yield similar apparent "clusters" (to evaluate if the 3-class structure is an artifact of forced clustering).

---

## 6. Continuous Subtype scoring

To map PDAC samples along a continuum:
1.  Calculate a continuous score for each sample:
    *   **Option A: ssGSEA Scoring:** Run single-sample GSEA (using the GSVA R package) to calculate independent enrichment scores for the Moffitt classical signature and basal-like signature. Define the continuous subtype axis score as:
        $$\text{Subtype Axis Score} = \text{ssGSEA}_{\text{Basal}} - \text{ssGSEA}_{\text{Classical}}$$
    *   **Option B: NMF Projection:** Project the expression matrix onto the Moffitt NMF components. Calculate the ratio of basal-like components (6 + 10) to classical components (1 + 2).
2.  Plot the distribution of the continuous scores (e.g., density plots, histograms).
3.  Test for multimodality in the continuous score distribution using **Hartigan's Dip Test**. A non-significant dip test ($p > 0.05$) supports a continuous unimodal distribution rather than discrete subtypes.
4.  Perform correlation analysis between the continuous scores and the discrete class assignments. Plot sample silhouette widths against their continuous scores.

---

## 7. Host–Microbiome Integration

Identify which microbial taxa or functional pathways are associated with the host subtype structure.

### 7.1 Taxonomic Associations (Compositional & Confounder-Adjusted)
1.  **Linear Modeling via MaAsLin2:**
    *   *Dependent Variable:* CLR-transformed abundance of microbial genera or species.
    *   *Independent Variable (Predictor):* Host continuous Subtype Axis Score (or discrete subtype labels as contrast).
    *   *Confounders (Covariates):* Tumor purity (ESTIMATE score), immune cell infiltration score, host sequencing library size, sequencing batch (`TO VERIFY` batch variable), age, sex, tumor stage.
    *   *Parameters:* Random effects for clinical batch if applicable; minimum abundance $= 0.001$; minimum prevalence $= 0.05$.
2.  **Composition-Aware Differential Abundance:**
    *   Use **ANCOM-BC** and **ALDEx2** to identify differentially abundant genera between discrete basal-like and classical subtypes, controlling for the same covariates where supported.
3.  **Significance:** Filter results by FDR-adjusted $p$-value ($q < 0.05$).

### 7.2 Host Pathway Associations
1.  Perform differential expression analysis between basal-like and classical tumors (DESeq2).
2.  Run Gene Set Enrichment Analysis (GSEA) against MSigDB Hallmark gene sets.
3.  Correlate the abundance of significantly associated bacterial genera (from Section 7.1) with host pathway enrichment scores (ssGSEA/GSVA) using partial Spearman correlation, adjusting for tumor purity.

---

## 8. External Validation

Validate the host molecular signatures and host-microbial associations in independent data:
1.  **Host Subtype Continuum Validation:**
    *   Download TCGA-PAAD expression data.
    *   Calculate the continuous Subtype Axis Score for TCGA-PAAD samples.
    *   Verify if the continuous score distribution exhibits the same unimodal shape (Hartigan's Dip Test $p > 0.05$).
2.  **Single-Cell/Spatial Validation:**
    *   Analyze public PDAC single-cell RNA-seq datasets (e.g., Peng et al. 2019).
    *   Evaluate if single tumor cells express classical and basal-like genes simultaneously (supporting a single-cell hybrid state) or if classical and basal-like genes are strictly expressed by different subclonal tumor cells (supporting cellular heterogeneity).
    *   Use spatial transcriptomics to verify the spatial co-localization of classical and basal-like programs and their proximity to immune and stromal elements.
3.  **Microbiome Validation (Exploratory):**
    *   TCGA-PAAD should primarily support host transcriptomic validation. Any TCGA-derived microbiome signal must be treated as exploratory because of extensive contamination and batch concerns.
    *   Do not assume that an identical external paired transcriptome-microbiome cohort exists. If attempting exploratory validation, treat any host-microbiome association results in external datasets as tentative and requiring careful verification of batch effects and sequencing platforms.

---

## 9. Therapeutic Target Prioritization

For genes overexpressed in basal-like PDAC (FDR $< 0.05$, log2 Fold Change $> 1.0$), prioritize candidates for therapeutic target discovery:
1.  **Dependency Filter (DepMap):**
    *   Extract CRISPR knockout Chronos scores for candidate genes in PDAC cell lines.
    *   Calculate the median dependency score in basal-like cell lines vs. classical cell lines (`TO VERIFY` classification of DepMap cell lines).
    *   Retain genes with a median Chronos score $< -0.5$ in basal-like cell lines, indicating gene essentiality.
2.  **Druggability Filter (DGIdb / ChEMBL):**
    *   Query prioritized genes in DGIdb to check for known small molecule inhibitors, monoclonal antibodies, or ongoing clinical trials.
    *   Retrieve binding affinities (Ki, IC50) for targets in ChEMBL to identify targets with potent chemical matter.
3.  **Clinical Survival Correlation:**
    *   Test if higher expression of the target gene is associated with significantly worse overall survival in TCGA-PAAD (Cox proportional hazards regression, adjusting for age, sex, stage; threshold hazard ratio $> 1.0$, $p < 0.05$).
4.  **Prioritization Matrix:**
    *   Construct a ranked table containing: Gene Symbol, Basal-Classical log2FC, DepMap Score, DGIdb Druggability (Yes/No), ChEMBL Compounds count, TCGA Survival Hazard Ratio, and Evidence Tier (based on the Evidence Policy).
