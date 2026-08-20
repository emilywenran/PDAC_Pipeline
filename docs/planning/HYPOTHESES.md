# Project Hypotheses

This document establishes the scientific hypotheses to be formally tested in the project. Each hypothesis is paired with its null hypothesis ($H_0$), required datasets, statistical endpoints, and interpretation boundaries.

---

## Hypothesis 1: Discrete vs. Continuous PDAC Subtyping

### Hypothesis
*   **Primary Hypothesis ($H_1$):** Pancreatic ductal adenocarcinoma (PDAC) tumors exhibit a continuous transcriptomic gradient along a classical-to-basal-like axis. The "hybrid" subtype does not represent a stable, discrete molecular entity, but rather intermediate transcriptional states, cellular mixtures, or tumor-stroma hybrid compositions.
*   **Null Hypothesis ($H_0$):** PDAC tumors partition into three stable, discrete molecular subtypes (classical, basal-like, hybrid) with well-defined, distinct cluster boundaries.

### Required Data
*   **Host Transcriptomics:** Normalized gene expression matrix from GSE172356 ($n=62$).
*   **Signature Gene Lists:** Moffitt et al. (2015) classical and basal-like gene lists (e.g., top genes from components 1, 2, 6, and 10).
*   **Verification Data:** Clinical annotations (to check for associations with clinical endpoints, `TO VERIFY`).

### Statistical Endpoints
*   **Consensus Clustering Metrics:** Consensus cumulative distribution function (CDF) curve, tracking the change in area under CDF (Delta K) across $k=2$ to $k=6$.
*   **Cluster Stability:** Bootstrap-resampled Jaccard similarity index of clusters (cutoff for stable cluster: Jaccard $> 0.85$).
*   **Cluster Separation:** Silhouette widths ($s_i$) calculated for each sample. A mean silhouette width close to 0 or negative values for the hybrid group indicates instability.
*   **Gradient Evaluation:** Fit a principal curve or perform diffusion map projection on the expression space. Evaluate if the distribution of samples along the first trajectory is unimodal (gradient) or multimodal (discrete clusters) using Hartigan's Dip Test ($p < 0.05$ indicating multimodality).

### Interpretation Limits
*   **Bulk Averaging:** GSE172356 consists of bulk RNA-seq data. An intermediate "hybrid" score can arise from either:
    1.  A homogeneous population of tumor cells expressing both classical and basal-like programs (true hybrid state).
    2.  A heterogeneous mixture of classical and basal-like tumor subclones.
    3.  Varying fractions of normal/activated stroma and immune infiltration.
*   This hypothesis cannot distinguish between these sub-cellular configurations using bulk data alone. Resolution requires single-cell or spatial transcriptomic validation.

---

## Hypothesis 2: Intratumoral Microbiome Subtype Association

### Hypothesis
*   **Primary Hypothesis ($H_1$):** Specific intratumoral bacterial taxa or functional modules are associated with the classical-to-basal-like transcriptomic axis, and these associations remain statistically significant after adjusting for tumor purity, immune infiltration, sequencing depth, batch effects, and clinical covariates.
*   **Null Hypothesis ($H_0$):** Intratumoral microbial features show no association with PDAC subtyping or continuous scores once technical and biological confounders are controlled.

### Required Data
*   **Microbiome Profiles:** Shotgun metagenomic profiles from PRJNA719915 ($n=62$).
*   **Negative Controls:** Environmental, kit, or PCR blanks/negative controls from PRJNA719915 (to identify contaminants, `TO VERIFY` availability in repository).
*   **Confounders Matrix:** Estimates of tumor purity (via ESTIMATE or ConsensusTME), stromal/immune infiltration fractions, host library sizes, sequencing depth, and clinical metadata (age, sex, pathological stage).

### Statistical Endpoints
*   **Decontamination:** R `decontam` prevalence/frequency classification (threshold $p < 0.1$ or custom stringency to identify contaminants).
*   **Association Signatures:** 
    *   For discrete subtypes: Differential abundance statistics (FDR-adjusted $p$-value $< 0.05$ using ANCOM-BC or ALDEx2).
    *   For continuous scores: Linear mixed-effect regression models (e.g., via MaAsLin2) tracking centered log-ratio (CLR) transformed taxonomic abundance against continuous host scores, yielding regression coefficients ($\beta$) and Benjamini-Hochberg adjusted $p$-values ($q < 0.05$).

### Interpretation Limits
*   **Taxonomic Resolution:** Shotgun metagenomic sequencing allows species-level and potential strain-level taxonomic profiling, but detection of low-abundance organisms and strain-level variations remains limited by sequencing depth and host read proportion.
*   **Biomass and Viability:** The presence of bacterial DNA does not guarantee the presence of live, metabolically active bacteria inside the tumor cells or microenvironment.
*   **Observational Data:** Association does not equal causation. A high abundance of a taxon in basal-like tumors could simply indicate that the necrotic, hypoxic, or inflammatory microenvironment of basal-like tumors is more permissive to its survival.

---

## Hypothesis 3: Association of Intratumoral Microbiome with Specific Host Pathways

### Hypothesis
*   **Primary Hypothesis ($H_1$):** The abundance of subtype-associated intratumoral taxa correlates with the upregulation of host inflammatory (e.g., IL-6/JAK/STAT3 signaling, TNF-$\alpha$ signaling), epithelial-mesenchymal transition (EMT), and metabolic pathways.
*   **Null Hypothesis ($H_0$):** Host pathway enrichment scores are independent of the abundance of intratumoral microbial taxa.

### Required Data
*   **Host Expression:** Gene expression matrix (GSE172356).
*   **Microbiome Taxonomy:** Composition-corrected taxonomic and functional abundance tables (PRJNA719915).
*   **Pathway Definitions:** Molecular Signatures Database (MSigDB) Hallmark gene sets.

### Statistical Endpoints
*   **Pathway Scoring:** Single-sample GSEA (ssGSEA) or GSVA enrichment scores calculated for each host sample.
*   **Correlation / Regression:** Partial Spearman correlation coefficients ($\rho$) or linear models adjusting for tumor purity, where the pathway score is the dependent variable and the CLR-abundance of the target taxon is the independent variable. Significance threshold: $q$-value $< 0.05$.

### Interpretation Limits
*   **Indirect Signaling:** The pathways enriched in the host may be triggered by stromal/immune cells responding to bacteria, rather than direct tumor cell-bacterial interactions.
*   **Functional Profiling:** Reconstructing microbial functional potential from shotgun metagenomic reads (e.g., via HUMAnN) profiles the presence and abundance of microbial metabolic pathway genes. However, DNA-based functional profiling measures functional potential rather than active microbial transcription (RNA-seq) or actual metabolite production (metabolomics).

---

## Hypothesis 4: Subtype-Specific Vulnerability Prioritization

### Hypothesis
*   **Primary Hypothesis ($H_1$):** Genes that are significantly overexpressed in basal-like tumors relative to classical tumors correspond to high-essentiality vulnerabilities in basal-like PDAC cell lines and represent targets with established chemical/clinical druggability.
*   **Null Hypothesis ($H_0$):** Subtype-specific differentially expressed genes do not show enrichment for cellular dependencies or druggable targets.

### Required Data
*   **Differential Expression Output:** Differentially expressed genes (DEGs) between basal-like and classical tumors (GSE172356).
*   **Dependency Screening:** Cancer Dependency Map (DepMap) Achilles CRISPR knockout screens across human pancreatic cancer cell lines.
*   **Drug Sensitivity:** PRISM repurposing database compound sensitivity profiles.
*   **Druggability Databases:** DGIdb and ChEMBL target annotations.

### Statistical Endpoints
*   **Dependency Contrast:** Distribution of Chronos dependency scores for candidate targets in basal-like vs. classical cell lines. A score $< -0.5$ indicates significant essentiality. Mann-Whitney U test comparison ($p < 0.05$, FDR corrected).
*   **Target Enrichment:** Hypergeometric test ($p < 0.05$) to evaluate if basal-like DEGs are enriched for druggable targets compared to the genomic background.

### Interpretation Limits
*   **In Vitro Disconnect:** DepMap and PRISM profiles are derived from cell lines grown in 2D or 3D monocultures. They lack the stromal desmoplasia, immune components, hypoxia, and intratumoral microbiome that characterize *in vivo* PDAC tumors. Target essentiality might be altered in the actual microenvironment.
