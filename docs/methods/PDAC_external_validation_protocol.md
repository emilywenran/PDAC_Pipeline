# PDAC External Validation Protocol (Phase 9B Execution Guide)

This document provides the standard operating protocol for executing Phase 9B (External Validation) of the PDAC transcriptional-state and host-microbiome findings. 

---

## 1. Prerequisites and Environment Setup

All analyses must be run in the project-local R and Python environments. Package management is handled via `renv` (R) and `uv` (Python).

### Software Requirements
- **R Version:** 4.5.3 (as per `renv.lock`)
- **Key R Packages:** `decoupleR`, `GSVA`, `limma`, `edgeR`, `WGCNA`, `tidyverse`, `data.table`, `BiocParallel`
- **Python Version:** >=3.10
- **Key Python Packages:** `pandas`, `numpy`, `statsmodels`, `scipy`

Initialize the environment:
```bash
# Verify R env status
R -e "renv::status()"

# Verify Python stack
python3 -c "import statsmodels; import scipy; print('Python stack ready')"
```

---

## 2. Step 1: Bulk Transcriptome Validation (Layer 1)

This step validates the continuous basal-classical axis, 2 Hallmark pathways, 34 TFs, and 7 WGCNA modules in independent bulk cohorts (TCGA-PAAD, GSE71729, GSE62452).

### Execution Steps
1. **Download Datasets:**
   - For **TCGA-PAAD**, download the processed FPKM-UQ or STAR-count matrix from the GDC portal.
   - For **GSE71729** and **GSE62452**, download the processed series matrix files from NCBI GEO:
     ```bash
     wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71729/matrix/GSE71729_series_matrix.txt.gz -P 02_data/external/GSE71729/
     wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62452/matrix/GSE62452_series_matrix.txt.gz -P 02_data/external/GSE62452/
     ```
2. **Preprocess Matrices:**
   - Map probe/transcript IDs to HGNC symbols.
   - Apply the duplicate-gene resolution policy (keep highest mean probe).
   - Apply the missing-gene policy (rescale scores if coverage is $\ge 80\%$; fail analysis if $< 80\%$).
3. **Calculate Scores:**
   - Compute the Moffitt50 axis reference score using the locked 50-gene contrast centroid model.
   - Run ssGSEA via `decoupleR::run_gsva` to compute Hallmark pathway scores.
   - Run VIPER via `decoupleR::run_viper` to compute TF activity scores.
   - Calculate standardized average rank scores for the 7 WGCNA modules.
4. **Statistical Testing:**
   - Run Spearman correlations and OLS regressions with HC3 robust standard errors:
     $$\text{Feature\_Score} \sim \text{Moffitt50\_Axis} + \text{Tumor\_Purity\_Covariate}$$
   - Apply Benjamini-Hochberg FDR correction separately within each cohort-feature family.
   - Run Cox proportional hazard models for continuous host features against survival:
     $$\text{CoxPH(Survival)} \sim \text{Feature\_Score} + \text{Age} + \text{Stage} + \text{Purity}$$

---

## 3. Step 2: Single-Cell Source Validation (Layer 2)

This step determines the cellular origin of the host programs. 

### Operational Sequencing & Authorization
To manage execution complexity, Phase 9B2 execution is staged:
- **Phase 9B2-primary:** Analyzes **`PENG_CRA001160`** (accession CRA001160) only.
- **Phase 9B2-supplementary:** Planned analyses of `LIN_GSE154778` (accession GSE154778), `MONCADA_GSE111672` (accession GSE111672), and `HWANG_GSE202051` (accession GSE202051) remain valid and planned but require separate execution authorization.
- **Operational Rationale:** This sequencing decision is purely operational to focus initial verification resources and does not imply that the supplementary datasets are scientifically or biologically unsuitable.

### Execution Steps (Phase 9B2-primary)

1. **Prepare Expression Matrices:**
   - Download the processed cell-by-gene matrices and cell-type metadata:
     ```bash
     # Example download accessions:
     # PENG_CRA001160: CNCB GSA CRA001160
     # LIN_GSE154778: GEO GSE154778 (GSE154778_dgeMtx.csv.gz)
     # MONCADA_GSE111672: GEO GSE111672
     # HWANG_GSE202051: GEO GSE202051
     ```
2. **Calculate Cell-Level Scores:**
   - Run decoupleR on the cell-level matrices to generate single-cell pathway, TF, and module scores.
3. **Generate Patient Pseudobulk:**
   - Aggregate expression (or cell-level scores) by patient within each major cell type (malignant epithelium, CAFs, immune cell compartments).
4. **Mixed-Effects Modeling:**
   - Fit patient-level mixed-effects models (treating patient as a random effect) to test cell-type-specific enrichment and patient-level reproducibility:
     $$\text{Cell\_Score} \sim \text{Cell\_Type} + (1 | \text{Patient\_ID})$$
   - Evaluate whether malignant cells exhibit co-activation of basal and classical programs.

---

## 4. Step 3: Spatial Transcriptomics Mapping (Layer 3)

This step maps the spatial localization of host programs using GSE202051, GSE274103, and GSM3405527.

### Execution Steps
1. **Prepare Spatial Coordinates:**
   - Load spot-level expression matrices and corresponding pixel coordinates (X, Y).
2. **Compute Spot-Level Scores:**
   - Calculate Moffitt50, Hallmark pathway, TF activity, and module scores for each spot.
3. **Spatial Autocorrelation and Compartment Contrast:**
   - Fit spatial regression models or spot-level ANOVA with patient blocking:
     $$\text{Spot\_Score} \sim \text{Histological\_Compartment} + (1 | \text{Patient\_ID})$$
   - Test whether pathway and WGCNA activities are enriched in epithelial vs. stromal coordinates.

---

## 5. Step 4: Independent Microbiome Validation (Layer 4)

This step tests whether `Ochrobactrum` is detected and replicates host associations in PRJNA542615 and EGAS00001004572.

### Execution Steps
1. **Download Raw FASTQs:**
   - Download raw sequencing reads from SRA/EGA (do not run during Phase 9A planning):
     ```bash
     # Example SRA download (Phase 9B execution only)
     prefetch PRJNA542615
     fastq-dump --split-files --gzip PRJNA542615/*.sra
     ```
2. **Host Depletion Pipeline:**
   - Run quality control and deplete human host reads by mapping to GRCh38:
     ```bash
     fastp -i read1.fastq.gz -I read2.fastq.gz -o qc_r1.fq.gz -O qc_r2.fq.gz
     bowtie2 -x grch38_index -1 qc_r1.fq.gz -2 qc_r2.fq.gz --un-conc-gz depleted_reads.fq.gz
     ```
3. **Taxonomic Profiling:**
   - Run Kraken2 and Bracken against the Standard RefSeq database:
     ```bash
     kraken2 --db standard_db --threads 16 --paired depleted_reads.1.fq.gz depleted_reads.2.fq.gz --report kraken_report.txt
     bracken -d standard_db -i kraken_report.txt -o bracken_output.txt -r 150 -l G
     ```
4. **Contamination Filtering:**
   - Compare `Ochrobactrum` abundance in tumor tissue against extraction and library preparation negative controls using `decontam` (prevalence or frequency method).
   - If negative controls are absent, tag the dataset as *contamination limited*.
5. **Association Analysis:**
   - For paired cohorts, fit:
     $$\text{Host\_Feature\_Score} \sim \text{CLR\_Ochrobactrum}$$
   - Report whether the coefficient matches the direction and magnitude observed in discovery.

---

## 6. Step 5: Synthesis and Evidence Classification

1. **Meta-Analysis:**
   - If at least three independent bulk cohorts exist, run a random-effects meta-analysis to compute pooled effect sizes and confidence intervals.
   - Calculate heterogeneity statistics ($I^2$, Cochrane's $Q$).
   - Perform leave-one-cohort-out sensitivity analysis.
2. **Evidence Grading:**
   - Map each validated feature into one of the 8 prospectively locked Phase 9A evidence categories.
   - Compile final evidence summaries.
3. **Reporting:**
   - Document all results in `04_analysis/09_external_validation/PHASE9B_EXTERNAL_VALIDATION_RESULTS.md`.
