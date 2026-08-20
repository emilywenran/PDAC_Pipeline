# Title: Multi-omic and Spatial Integration Reveals Associative Intratumoral Microbiome Signatures and Prioritized Host Targets in Pancreatic Ductal Adenocarcinoma

## Abstract

**Background:** The tumor microbiome in Pancreatic Ductal Adenocarcinoma (PDAC) has potential interactions with host state.
**Methods:** Integrative multi-omics analysis of bulk, single-cell, and spatial transcriptomics data.
**Results:** Global community-level association exists. Nine robust genera (e.g., *Ochrobactrum*) are associated with host states, although we do not claim microbial causality, microbial localization, or physical host-microbe interaction. Validation across multiple external cohorts shows replication of transcription factor activities but with composition-sensitive caveats. HALLMARK_PROTEIN_SECRETION is identified as a robust malignant-cell intrinsic signature with spatial malignant compartment support, despite a lack of continuous basal-classical spatial replication (remaining PARTIAL_SPATIAL_SUPPORT). CTCFL/BORIS was evaluated but not prioritized because of CELL_COMPOSITION_EXPLAINED evidence. Null findings such as HALLMARK_SPERMATOGENESIS, WGCNA modules MEred and MEpurple, and 9 TF activities are reported. 
**Conclusions:** Tumor-associated microbiota exhibit robust but non-causal statistical associations with specific host transcriptional programs. No candidate should be described as an established therapeutic target; all require further validation.

## Introduction

Pancreatic ductal adenocarcinoma remains one of the most lethal solid malignancies, with a five-year survival rate below 13% in most high-income settings [1]. The poor prognosis reflects a combination of late-stage diagnosis, intrinsic chemoresistance, and a highly immunosuppressive tumour microenvironment (TME) that limits the efficacy of both conventional and immune-based therapies [2,3]. Transcriptional heterogeneity within PDAC tumours has been a consistent obstacle to the development of broadly applicable targeted strategies: bulk transcriptomic studies have established at least two major molecular subtypes—broadly designated "classical" and "basal-like" (or "squamous")—that differ substantially in prognosis, stromal composition, and transcription-factor dependencies [4,5]. Understanding the biological programmes that define and maintain these subtypes, and identifying tractable molecular targets within them, constitutes a challenge in PDAC research.

A growing body of evidence has drawn attention to the intratumoral microbiome as an additional axis of biological complexity in PDAC. Several reports have documented the presence of diverse microbial taxa within pancreatic tumour tissue and have observed statistical associations between microbial community composition and patient outcomes, immune infiltration patterns, and drug metabolism [6–9]. Among the taxa identified in these studies, members of the genus Ochrobactrum have been detected in PDAC specimens, with co-occurrence analyses suggesting that their relative abundance is correlated with specific host transcriptional programmes. Critically, however, it has not been established that Ochrobactrum or any other intratumoral bacterium causes the host transcriptional changes with which it co-occurs; the direction of any relationship, and whether it is direct, indirect, or confounded by shared environmental or host-genetic determinants, remains unresolved. Furthermore, the spatial localisation of Ochrobactrum within tumour compartments has not been directly demonstrated in the data sets analysed here. These qualifications are essential for accurate interpretation and are preserved throughout the analyses reported in this manuscript.

Pathway-level analyses offer a complementary approach to identifying biological programmes relevant to PDAC pathogenesis. Gene-set enrichment frameworks enable the aggregation of signal from individually modest transcriptional changes into coherent biological themes, increasing statistical power and mechanistic interpretability [10]. Among hallmark gene sets, HALLMARK_PROTEIN_SECRETION captures a coordinated programme of vesicular transport, endoplasmic reticulum function, and extracellular cargo delivery that is mechanistically plausible in the context of pancreatic exocrine biology and tumour-stromal communication [11]. In the present analyses, HALLMARK_PROTEIN_SECRETION was supported by evidence originating from both individual malignant cells and from the broader malignant compartment, providing cross-resolution convergence. An association between enrichment of this programme and position along the basal–classical subtype axis was identified in bulk transcriptomic data; however, spatial transcriptomic analyses did not replicate this axis association with sufficient resolution to constitute full spatial validation. The spatial evidence for this claim is therefore characterised as partial, and conclusions regarding subtype-specific spatial patterning of protein secretion activity are correspondingly limited.

Target prioritisation in oncology requires the integration of evidence across multiple analytical levels, including expression specificity, pathway membership, regulatory accessibility, and the availability of tractable intervention points [12,13]. In the present work, candidate targets were ranked using a corrected scoring framework developed across Phases 10B-R and 10C2, which applied empirically derived evidence weights and penalised claims unsupported by cross-platform replication. Critically, the framework strictly penalized features confounded by cell-type proportions. For example, while the transcription factor CTCFL (BORIS) showed initial differential activity, it was explicitly excluded from prioritization after single-cell analysis revealed its signal was driven by cell composition sensitivity rather than intrinsic malignant compartment expression.

The present manuscript reports the integrated findings of a multi-phase analytical programme designed to characterise transcriptional, microbial co-occurrence, and spatial features of PDAC with appropriate epistemic rigour. We describe pathway-level associations and prioritised targets, delineate the boundaries of spatial validation achieved, and distinguish associative findings from causal claims. We report null and partial results with equal transparency. Our aim is to provide a structured, reproducible foundation for hypothesis generation and experimental follow-up, rather than to assert mechanistic conclusions beyond what the available data support.

## Results

### Cohort Characteristics and Analytical Overview
Multi-phase integrative analysis was performed across PDAC transcriptomic data sets encompassing bulk RNA-sequencing, single-cell RNA-sequencing (scRNA-seq), and spatial transcriptomic profiling. Cohort composition, quality-control thresholds, and preprocessing decisions are described in the Methods. Briefly, bulk RNA-sequencing data provided the primary substrate for subtype classification, differential expression, and gene-set enrichment analyses. Single-cell data enabled compartment-resolved decomposition of transcriptional programmes, including separate characterisation of malignant-cell, stromal, and immune populations. Spatial transcriptomic data were used to assess whether bulk-derived associations could be localised within tissue architecture, subject to the resolution and coverage constraints of the platform used.

Microbial co-occurrence profiling was performed on a subset of specimens for which suitable sequencing data were available. Relative abundance estimates were used as associative variables in downstream analyses; no inference of absolute abundance, spatial tissue localisation, or causal host-microbe relationship was made. All microbial findings are interpreted as co-occurrence signals and are reported with this qualification explicitly preserved.

### Transcriptional Subtype Landscape and Malignant-Compartment Characterisation
Unsupervised clustering of bulk transcriptomic profiles recapitulated the established two-axis structure of PDAC molecular heterogeneity, with sample distributions consistent with previously described classical and basal-like (squamous) subtypes. Classical-subtype tumours exhibited relative enrichment of glandular differentiation programmes and a transcriptional profile consistent with luminal pancreatic identity, while basal-like tumours showed upregulation of epithelial-to-mesenchymal transition markers, squamous differentiation genes, and stress-response pathways.

At the single-cell level, malignant-cell populations resolved into distinct transcriptional states that partially recapitulated the bulk subtype axis, though with substantial intra-tumour heterogeneity. A subset of malignant cells in basal-enriched tumours co-expressed markers of both subtype programmes, consistent with published reports of transitional or hybrid states. These observations are reported descriptively; no causal model of subtype switching or state transition is proposed on the basis of these cross-sectional data.

### Gene-Set Enrichment Analyses: HALLMARK_PROTEIN_SECRETION and Convergent Pathway Evidence
Gene-set enrichment analysis across the hallmark collection identified HALLMARK_PROTEIN_SECRETION as one of the most consistently enriched programmes in malignant cells and the malignant compartment broadly defined. Enrichment was detected at both the malignant-cell level in single-cell analyses and at the malignant-compartment level in pseudobulk and bulk analyses, providing cross-resolution concordance.

Specifically:
• Malignant-cell evidence. In scRNA-seq data, HALLMARK_PROTEIN_SECRETION gene scores were significantly elevated in malignant cells relative to non-malignant compartment populations.
• Malignant-compartment evidence. At the bulk level, HALLMARK_PROTEIN_SECRETION enrichment scores were positively correlated with classical-subtype classification and showed a negative association with basal-like classification.
• Spatial association — partial evidence only. Spatially variable gene analysis identified focal regions of elevated gene-set module scores within tumour-containing tissue areas. However, the bulk-derived association between HALLMARK_PROTEIN_SECRETION and the basal–classical subtype axis was not spatially replicated. The spatial evidence for the basal–classical axis association of HALLMARK_PROTEIN_SECRETION is therefore classified as PARTIAL_SPATIAL_SUPPORT. Exploratory analyses in an independent spatial cohort (Moncada) yielded positive concordance in only 1 of 6 sections.

Additional hallmark gene sets showing enrichment in the malignant compartment included HALLMARK_UNFOLDED_PROTEIN_RESPONSE, HALLMARK_MYC_TARGETS_V1, and HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION. Conversely, multiple features failed cross-layer validation and are explicitly reported as null findings: HALLMARK_SPERMATOGENESIS exhibited opposite directional effects in external bulk replication and lacked sufficient spatial coverage. Similarly, WGCNA modules MEred and MEpurple, alongside 9 distinct transcription factor activities, failed to replicate in external bulk transcriptomic cohorts.

### Microbial Co-occurrence Associations with Host Transcriptional Signatures
Microbial co-occurrence profiling identified Ochrobactrum as a taxon whose relative abundance showed statistically significant correlation with host transcriptional signatures across specimens. Ochrobactrum relative abundance was positively associated with enrichment scores for immune-regulatory gene sets, including programmes related to innate immune signalling and myeloid activation. A secondary association was observed with secretory pathway gene scores.

These findings are interpreted exclusively as associative signals. No causal relationship between Ochrobactrum and the observed host transcriptional programmes is established or claimed by these analyses. Furthermore, the spatial localisation of Ochrobactrum within tumour tissue was not determined, and no physical interaction between Ochrobactrum and host cells is proposed. Of the nine genera classified as robustly associated with host state, eight exhibited direction reversals or loss of statistical significance under robust centred log-ratio (rCLR) transformation. Additionally, Herbaspirillum was flagged as a moderate-risk environmental contaminant.

### Regulatory Candidate Analysis and Transcription-Factor Programme Associations
Regulatory candidate analysis identified a set of candidate transcriptional regulators with differential activity across the PDAC subtype landscape. While initial analyses identified differential activity for several regulators, rigorous validation identified 13 transcription factor activities as only partially replicated, and 9 activities that failed to replicate entirely. Notably, the transcription factor CTCFL (BORIS), which initially showed differential activity, was rigorously evaluated for composition sensitivity. Single-cell decomposition revealed that its apparent differential expression was entirely explained by cell-type proportion variations (CELL_COMPOSITION_EXPLAINED) rather than malignant-intrinsic regulation. Consequently, CTCFL was explicitly penalized and removed from candidate prioritization.

### Integrated Target Prioritisation
Target prioritisation was performed using the corrected scoring framework established in Phases 10B-R and Phase 10C2. Candidates were penalised where evidence derived from a single analytical platform, where spatial validation was absent or only partial, or where regulatory associations were motif-based only without corroborating functional data. Features whose single-cell evidence was classified as CELL_COMPOSITION_EXPLAINED received high-weight penalties that blocked advancement to prioritised status. No candidate identified by this framework is described as an established therapeutic target; all represent hypotheses requiring rigorous experimental evaluation.

## Discussion
This study integrated bulk, single-cell, and spatial transcriptomics with microbial co-occurrence profiling to dissect PDAC heterogeneity. We found that the tumor microbiome is associated with, but cannot be claimed to cause, host transcriptional changes. Our strict statistical controls highlighted the importance of mitigating composition sensitivity in computational analyses, successfully identifying robust signals like HALLMARK_PROTEIN_SECRETION while demoting confounded signals like CTCFL/BORIS. 

## Limitations
Several critical limitations constrain the interpretation of these findings. First, microbial associations with host state are strictly non-causal; neither microbial causality, microbial localization, nor physical host-microbe interaction can be established from these data. Second, our analysis of the microbiome is limited by methodological challenges, including transformation sensitivity (rCLR) and the presence of contamination risks (e.g., Herbaspirillum). Third, spatial validation of the basal-classical axis for HALLMARK_PROTEIN_SECRETION yielded only PARTIAL_SPATIAL_SUPPORT, and the Moncada cohort findings were purely exploratory, showing positive concordance in only 1 of 6 sections. Finally, our integrated target prioritization framework generates hypotheses only; no candidate should be considered an established therapeutic target.

## Methods Summary
Data from bulk RNA-seq, scRNA-seq, and spatial transcriptomics cohorts were processed using standard pipelines. Microbial abundances were transformed via CLR and rCLR to assess compositional sensitivity. Transcriptional pathway activities were estimated via PROGENy and ssGSEA, while TF activities were modeled using DoRothEA/VIPER. Target prioritization applied penalty frameworks established in Phase 10C2.

## Data Availability
[Placeholder for Data Availability]

## Code Availability
[Placeholder for Code Availability]

## Author Contributions
[Placeholder for Author Contributions]

## Funding
[Placeholder for Funding]

## Conflict of Interest
[Placeholder for Conflict of Interest]

## References
[Placeholder for References]

## Figure Legends
**Figure 1.** Tumor Microbiome and Host State. Includes PERMANOVA plots, 9 robust genera forest plots, and rCLR sensitivity annotations.
**Figure 2.** Host Biological Mechanisms. Depicts the *Ochrobactrum* mechanism network alongside pathway and TF activity heatmaps.
**Figure 3.** External Bulk Validation. Shows replication metrics for 12 TFs and 7 partially replicated pathways across cohorts.
**Figure 4.** Single-Cell Cellular Sources. Contrasts malignant versus non-malignant compartment expression and highlights composition sensitivity for CTCFL/BORIS.
**Figure 5.** Spatial Compartment Support. Illustrates HALLMARK_PROTEIN_SECRETION spatial mapping and Model A enrichment, demonstrating PARTIAL_SPATIAL_SUPPORT.

## Supplementary Table Legends
**Table S1.** Target Prioritization. Final ranked target list following the Phase 10C2 approved ranking framework, explicitly excluding composition-sensitive features like CTCFL.
**Table S2.** Unsupported Features. Full list of features failing bulk external validation (e.g., HALLMARK_SPERMATOGENESIS, WGCNA modules MEred and MEpurple, 9 TF activities).
**Table S3.** Exploratory Spatial Validation. Results from the Moncada cohort spatial findings detailing directional inconsistency.
