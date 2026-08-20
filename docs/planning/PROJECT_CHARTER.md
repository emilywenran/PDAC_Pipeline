# Project Charter: PDAC Microbiome–Transcriptome Plasticity

## 1. Scientific Background
Pancreatic ductal adenocarcinoma (PDAC) is one of the most lethal malignancies, characterized by extreme therapeutic resistance, a dense desmoplastic stroma, and a complex tumor microenvironment. Genomic profiling has revealed that despite significant clinical heterogeneity, PDAC exhibits relatively few recurrent driver mutations (e.g., *KRAS*, *TP53*, *CDKN2A*, *SMAD4*). To explain this clinical variation, transcriptomic studies have proposed molecular subtypes. The most widely accepted framework, established by Moffitt et al. (2015), defines two primary tumor-specific subtypes:
*   **Classical:** Associated with ductal differentiation markers, better therapeutic response to fluorouracil-based regimens (e.g., FOLFIRINOX), and improved overall survival.
*   **Basal-like:** Associated with EMT markers, poor differentiation, chemoresistance, and significantly worse clinical outcomes.

Subsequent studies have proposed additional subtypes (e.g., Collisson et al. and Bailey et al.) or intermediate "hybrid" states. However, whether these subtypes represent stable, discrete cellular states or arbitrary partitions of a continuous transcriptomic continuum (plasticity) remains a matter of intense debate. 

Recent research, including Guo et al. (2021) (*Communications Biology*), has suggested that the intratumoral microbiome is not merely a passive bystander but actively shapes the tumor phenotype. Specifically, Guo et al. reported that basal-like tumors harbor distinct microbial communities (e.g., enriched in *Acinetobacter*, *Pseudomonas*, and *Sphingopyxis*) that correlate with inflammatory host programs. 

However, tumor tissue is a low-biomass microbial environment. Metagenomic sequencing of low-biomass tissues is highly susceptible to contamination from reagents (the "kitome"), environment, and clinical handling. Furthermore, bulk transcriptomics is a mixture of malignant, stromal, and immune cells. Without adjusting for tumor purity, immune infiltration, sequencing depth, and technical batch effects, associations between the microbiome and host subtypes may be heavily confounded.

---

## 2. Project Objectives

### Primary Objective
Evaluate whether the basal-like, classical, and hybrid PDAC molecular subtypes represent stable, discrete transcriptomic states or positions along a continuous basal–classical transcriptional gradient. Subsequently, determine which intratumoral microbiome features (taxa or functional modules) remain robustly associated with this subtype structure after accounting for tumor purity, immune/stromal composition, sequencing depth, batch effects, and available clinical covariates.

### Secondary Objectives
1.  **Reproduction:** Reproduce the discrete three-class subtype assignments (basal-like, classical, hybrid) of the 62 primary PDAC tumors from host transcriptomic data (NCBI GEO: GSE172356) using the consensus clustering protocol described by Guo et al. (2021).
2.  **Continuum Comparison:** Compare the discrete three-class classifications against continuous transcriptomic scores derived from the Moffitt classical and basal-like signatures.
3.  **Hybrid Characterization:** Evaluate the stability, heterogeneity, and intermediate nature of the "hybrid" subtype using silhouette analysis, bootstrapping, and expression deconvolution.
4.  **Microbiome Profiling & Integration:** Process raw shotgun metagenomic sequencing data (NCBI BioProject: PRJNA719915) through raw-read QC, host-read depletion, taxonomic profiling (Kraken2/Bracken and/or MetaPhlAn), negative-control contamination assessment, prevalence/abundance filtering, and functional profiling. Identify robust microbiome-host pathway associations using composition-aware statistics while adjusting for confounding variables.
5.  **External Validation:** Validate host transcriptomic signatures and continuous scoring distributions in independent bulk (e.g., TCGA-PAAD, ICGC) and single-cell/spatial transcriptomics datasets. Any validation of microbiome signals in external cohorts (such as TCGA-derived microbiome data) must be treated as exploratory due to batch and contamination concerns, and we do not assume that an identical external paired transcriptome-microbiome cohort exists.
6.  **Therapeutic Prioritization:** Prioritize therapeutic targets for basal-like tumors by integrating subtype-specific expression, dependency scores (DepMap/PRISM), and druggability data (DGIdb/ChEMBL).

---

## 3. Scope and Exclusions

### In Scope
*   **Data Retrieval:** Downloading and auditing expression data (GSE172356), microbiome sequencing data (PRJNA719915), and external validation datasets.
*   **Bioinformatic Pipelines:** Host RNA-seq quality control, normalization, and deconvolution (tumor purity/immune profiling). Microbiome shotgun metagenomic profiling (including raw-read QC, host depletion, taxonomic profiling with Kraken2/Bracken and/or MetaPhlAn, contamination control, and functional profiling).
*   **Statistical Analysis:** Consensus clustering, silhouette analysis, continuous signature scoring, robust regression (e.g., MaAsLin2), and compositional data analysis (e.g., ALDEx2, ANCOM-BC).
*   **Target Prioritization:** Integration of public databases (DepMap, ChEMBL, DGIdb) for vulnerability profiling.
*   **Documentation:** Reproducible code, planning, and decision logs.

### Out of Scope / Exclusions
*   **Wet-Lab Validation:** No generation of new biological samples, no *in vitro* or *in vivo* experiments.
*   **Patient Enrollment:** No collection of new patient cohorts or clinical trials.
*   **Clinical Diagnostic Tool Deployment:** No production of clinical software or diagnostic tests.
*   **Causal Claims:** No assertion of biological causality from cross-sectional association data alone.

---

## 4. Expected Outputs
1.  **Quality Control Reports:** Diagnostic summaries for both host transcriptomics (purity, library size, batch effects) and microbiome sequencing (raw-read QC, host depletion metrics, contamination filters, profiling statistics).
2.  **Subtype Classification Matrix:** A validated classification mapping each sample to its discrete subtype (reproduced) and continuous classical/basal-like score.
3.  **Microbiome Abundance Tables:** Cleaned, contamination-filtered taxonomic and functional abundance tables with composition-aware normalizations.
4.  **Host-Microbiome Association Tables:** Lists of taxa and functional modules significantly associated with host subtypes and continuous programs, including effect sizes, p-values, and FDR values, adjusted for confounders.
5.  **Validation Reports:** Results demonstrating the reproducibility of the subtyping and microbial correlations in independent datasets.
6.  **Therapeutic Target Vulnerability List:** A ranked table of candidate genes prioritized for therapeutic development in basal-like tumors.
7.  **Reproducible Pipeline:** A fully documented, containerized code repository (scripts/notebooks) to reproduce all analyses.

---

## 5. Criteria for Project Success
*   **High Subtype Concordance:** Achieving >95% classification concordance with the original Guo et al. (2021) study when reproducing the discrete three-class labels using the same methodology.
*   **Rigorous Contamination Control:** Successful identification and assessment of kit/environmental contaminants in low-biomass shotgun metagenomic data, leaving only high-confidence biological signals.
*   **Statistically Guarded Associations:** All host-microbiome association testing must utilize composition-aware algorithms, false discovery rate (FDR) control, and adjustment for tumor purity, library size, and batch effects.
*   **Successful External Validation:** Replicating the continuous host transcriptomic gradient and candidate host pathway signatures in at least one independent cohort (e.g., TCGA-PAAD). Note that TCGA-PAAD primarily supports host transcriptomic validation, while any TCGA-derived microbiome signal is strictly exploratory.
*   **Zero Data Leakage:** Maintaining strict isolation between feature selection folds during internal cross-validation and resampling steps.
*   **Fully Documented Decisions:** Clear traceability of all computational decisions and pipeline parameters in the decision log.
