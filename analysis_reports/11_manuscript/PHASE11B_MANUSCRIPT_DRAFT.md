## 1. Introduction

  Pancreatic ductal adenocarcinoma remains one of the most lethal solid malignancies, with a five-year survival rate below 13% in most high-income settings [1]. The poor prognosis reflects a combination
  of late-stage diagnosis, intrinsic chemoresistance, and a highly immunosuppressive tumour microenvironment (TME) that limits the efficacy of both conventional and immune-based therapies [2,3].
  Transcriptional heterogeneity within PDAC tumours has been a consistent obstacle to the development of broadly applicable targeted strategies: bulk transcriptomic studies have established at least two
  major molecular subtypes—broadly designated "classical" and "basal-like" (or "squamous")—that differ substantially in prognosis, stromal composition, and transcription-factor dependencies [4,5].
  Understanding the biological programmes that define and maintain these subtypes, and identifying tractable molecular targets within them, constitutes a challenge in PDAC research.

  A growing body of evidence has drawn attention to the intratumoral microbiome as an additional axis of biological complexity in PDAC. Several reports have documented the presence of diverse microbial
  taxa within pancreatic tumour tissue and have observed statistical associations between microbial community composition and patient outcomes, immune infiltration patterns, and drug metabolism [6–9].
  Among the taxa identified in these studies, members of the genus Ochrobactrum have been detected in PDAC specimens, with co-occurrence analyses suggesting that their relative abundance is correlated
  with specific host transcriptional programmes. Critically, however, it has not been established that Ochrobactrum or any other intratumoral bacterium causes the host transcriptional changes with which
  it co-occurs; the direction of any relationship, and whether it is direct, indirect, or confounded by shared environmental or host-genetic determinants, remains unresolved. Furthermore, the spatial
  localisation of Ochrobactrum within tumour compartments has not been directly demonstrated in the data sets analysed here. These qualifications are essential for accurate interpretation and are
  preserved throughout the analyses reported in this manuscript.

  Pathway-level analyses offer a complementary approach to identifying biological programmes relevant to PDAC pathogenesis. Gene-set enrichment frameworks enable the aggregation of signal from
  individually modest transcriptional changes into coherent biological themes, increasing statistical power and mechanistic interpretability [10]. Among hallmark gene sets, HALLMARK_PROTEIN_SECRETION
  captures a coordinated programme of vesicular transport, endoplasmic reticulum function, and extracellular cargo delivery that is mechanistically plausible in the context of pancreatic exocrine biology
  and tumour-stromal communication [11]. In the present analyses, HALLMARK_PROTEIN_SECRETION was supported by evidence originating from both individual malignant cells and from the broader malignant
  compartment, providing cross-resolution convergence. An association between enrichment of this programme and position along the basal–classical subtype axis was identified in bulk transcriptomic data;
  however, spatial transcriptomic analyses did not replicate this axis association with sufficient resolution to constitute full spatial validation. The spatial evidence for this claim is therefore
  characterised as partial, and conclusions regarding subtype-specific spatial patterning of protein secretion activity are correspondingly limited.

  Target prioritisation in oncology requires the integration of evidence across multiple analytical levels, including expression specificity, pathway membership, regulatory accessibility, and the
  availability of tractable intervention points [12,13]. In the present work, candidate targets were ranked using a corrected scoring framework developed across Phases 10B-R and 10C2, which applied
  empirically derived evidence weights and penalised claims unsupported by cross-platform replication. Critically, the framework strictly penalized features confounded by cell-type proportions. For example,
  while the transcription factor CTCFL (BORIS) showed initial differential activity, it was explicitly excluded from prioritization after single-cell analysis revealed its signal was driven by cell composition
  sensitivity rather than intrinsic malignant compartment expression.

  The present manuscript reports the integrated findings of a multi-phase analytical programme designed to characterise transcriptional, microbial co-occurrence, and spatial features of PDAC with
  appropriate epistemic rigour. We describe pathway-level associations and prioritised targets, delineate the boundaries of spatial validation achieved, and distinguish associative findings from causal
  claims. We report null and partial results with equal transparency. Our aim is to provide a structured, reproducible foundation for hypothesis generation and experimental follow-up, rather than to assert
  mechanistic conclusions beyond what the available data support.
  ──────
  [References to be populated from project bibliography in final submission draft.]
  ──────
  End of Section 1 Draft
  ──────
  │ Compliance note: No microbial causality, localisation, or physical interaction is claimed. Basal–classical spatial validation is explicitly qualified as partial. CTCFL/BORIS is explicitly
  │ excluded due to cell composition sensitivity. Ochrobactrum is described in associative terms throughout. HALLMARK_PROTEIN_SECRETION claims are bounded to malignant-cell and malignant-compartment support.
  │ Target ranking references Phase 10B-R and Phase 10C2 corrected outputs.

 ## 2. Results

  ### 2.1 Cohort Characteristics and Analytical Overview
  Multi-phase integrative analysis was performed across PDAC transcriptomic data sets encompassing bulk RNA-sequencing, single-cell RNA-sequencing (scRNA-seq), and spatial transcriptomic profiling. Cohort
  composition, quality-control thresholds, and preprocessing decisions are described in the Methods. Briefly, bulk RNA-sequencing data provided the primary substrate for subtype classification,
  differential expression, and gene-set enrichment analyses. Single-cell data enabled compartment-resolved decomposition of transcriptional programmes, including separate characterisation of malignant-
  cell, stromal, and immune populations. Spatial transcriptomic data were used to assess whether bulk-derived associations could be localised within tissue architecture, subject to the resolution and
  coverage constraints of the platform used.

  Microbial co-occurrence profiling was performed on a subset of specimens for which suitable sequencing data were available. Relative abundance estimates were used as associative variables in downstream
  analyses; no inference of absolute abundance, spatial tissue localisation, or causal host-microbe relationship was made. All microbial findings are interpreted as co-occurrence signals and are reported
  with this qualification explicitly preserved.
  ──────
  ### 2.2 Transcriptional Subtype Landscape and Malignant-Compartment Characterisation
  Unsupervised clustering of bulk transcriptomic profiles recapitulated the established two-axis structure of PDAC molecular heterogeneity, with sample distributions consistent with previously described
  classical and basal-like (squamous) subtypes. Classical-subtype tumours exhibited relative enrichment of glandular differentiation programmes and a transcriptional profile consistent with luminal
  pancreatic identity, while basal-like tumours showed upregulation of epithelial-to-mesenchymal transition markers, squamous differentiation genes, and stress-response pathways, in line with prior
  classifications [Refs].

  At the single-cell level, malignant-cell populations resolved into distinct transcriptional states that partially recapitulated the bulk subtype axis, though with substantial intra-tumour heterogeneity.
  A subset of malignant cells in basal-enriched tumours co-expressed markers of both subtype programmes, consistent with published reports of transitional or hybrid states. These observations are reported
  descriptively; no causal model of subtype switching or state transition is proposed on the basis of these cross-sectional data.
  Stromal compartment decomposition identified fibroblast, endothelial, and myeloid populations in proportions consistent with the desmoplastic architecture characteristic of PDAC. Immune compartment
  composition, including cytotoxic T-lymphocyte and macrophage fractions, varied across samples and showed nominal associations with subtype classification that are described below in the context of
  pathway analyses.
  ──────
  ### 2.3 Gene-Set Enrichment Analyses: HALLMARK_PROTEIN_SECRETION and Convergent Pathway Evidence
  Gene-set enrichment analysis across the hallmark collection identified HALLMARK_PROTEIN_SECRETION as one of the most consistently enriched programmes in malignant cells and the malignant compartment
  broadly defined. Enrichment was detected at both the malignant-cell level in single-cell analyses and at the malignant-compartment level in pseudobulk and bulk analyses, providing cross-resolution
  concordance that was not observed for all candidate pathways evaluated.

  Specifically:
  • Malignant-cell evidence. In scRNA-seq data, HALLMARK_PROTEIN_SECRETION gene scores were significantly elevated in malignant cells relative to non-malignant compartment populations (fibroblasts,
  endothelial cells, T cells, myeloid cells), with the difference surviving correction for multiple comparisons. This enrichment was not uniform across all malignant cells; a gradient of programme
  activity was observed across the single-cell landscape, with higher scores in classical-aligned transcriptional states.
  • Malignant-compartment evidence. At the bulk level, HALLMARK_PROTEIN_SECRETION enrichment scores were positively correlated with classical-subtype classification and showed a negative association with
  basal-like classification, a pattern consistent with the single-cell findings. This association was robust across multiple enrichment scoring methods evaluated.
  • Spatial association — partial evidence only. Spatial transcriptomic analysis was performed to assess whether the HALLMARK_PROTEIN_SECRETION enrichment gradient observed in bulk and single-cell data
  could be spatially resolved within tissue sections. Spatially variable gene analysis identified focal regions of elevated gene-set module scores within tumour-containing tissue areas; however, the
  spatial distribution of high-scoring spots did not reliably track with morphologically defined tumour architecture in a manner that would permit confident assignment to specific tissue compartments.
  Crucially, the bulk-derived association between HALLMARK_PROTEIN_SECRETION and the basal–classical subtype axis was not spatially replicated: spots with high programme scores did not show consistent
  spatial concordance with subtype-informative marker expression. The spatial evidence for the basal–classical axis association of HALLMARK_PROTEIN_SECRETION is therefore classified as
  PARTIAL_SPATIAL_SUPPORT. Conclusions regarding subtype-specific spatial patterning of this programme must be regarded as provisional pending analysis with higher-resolution spatial platforms, particularly
  as exploratory analyses in an independent spatial cohort (Moncada) yielded positive concordance in only 1 of 6 sections.

  Additional hallmark gene sets showing enrichment in the malignant compartment included HALLMARK_UNFOLDED_PROTEIN_RESPONSE, HALLMARK_MYC_TARGETS_V1, and HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.
  Conversely, multiple features failed cross-layer validation and are reported as null findings: HALLMARK_SPERMATOGENESIS exhibited opposite directional effects in external bulk replication and lacked
  sufficient spatial coverage (<80%) for validation. Similarly, WGCNA modules MEred and MEpurple, alongside 9 distinct transcription factor activities, failed to replicate in external bulk transcriptomic cohorts.
  ──────
  ### 2.4 Microbial Co-occurrence Associations with Host Transcriptional Signatures

  Microbial co-occurrence profiling identified Ochrobactrum as a taxon whose relative abundance showed statistically significant correlation with host transcriptional signatures across specimens for which
  paired data were available. The following findings are reported:
  Ochrobactrum relative abundance was positively associated with enrichment scores for immune-regulatory gene sets, including programmes related to innate immune signalling and myeloid activation. A
  secondary association was observed with secretory pathway gene scores, including a nominal positive correlation with HALLMARK_PROTEIN_SECRETION enrichment. These associations were identified in co-
  occurrence analyses and persisted after adjustment for library size and batch covariates.

  These findings are interpreted exclusively as associative signals. No causal relationship between Ochrobactrum and the observed host transcriptional programmes is established or claimed by these
  analyses. The observed correlations are consistent with multiple non-exclusive explanations, including: (i) shared host or environmental determinants that independently influence both microbial
  community composition and transcriptional state; (ii) host transcriptional programmes that alter the intratumoral microenvironment in ways that affect microbial community structure; (iii) unmeasured
  confounders not accounted for in the available data. The directionality of any relationship cannot be determined from the cross-sectional data analysed here.

  Furthermore, the spatial localisation of Ochrobactrum within tumour tissue was not determined in these analyses. The microbial profiling data used do not permit inference about whether Ochrobactrum is
  present preferentially in malignant, stromal, or other tissue compartments, or about the proximity of detected microbial signal to any specific host cell population. No physical interaction between
  Ochrobactrum and host cells is proposed or supported by these data. These limitations are noted explicitly to prevent over-interpretation of the co-occurrence findings.

  Of the nine genera classified as robustly associated with host state in community-level analyses, eight exhibited direction reversals or loss of statistical significance under robust centred log-ratio
  (rCLR) transformation, indicating that their associations are sensitive to compositional normalisation strategy. Additionally, Herbaspirillum was flagged as a moderate-risk environmental contaminant on
  the basis of its prevalence in reagent and environmental control databases. These sensitivities do not invalidate the co-occurrence findings but impose important interpretive constraints.

  The association between Ochrobactrum relative abundance and host secretory and immune-modulatory programme activity is regarded as a hypothesis-generating finding warranting prospective experimental
  investigation under controlled conditions, including culture-based or gnotobiotic-model approaches capable of establishing directional causality.
  ──────
  ### 2.5 Regulatory Candidate Analysis and Transcription-Factor Programme Associations

  Regulatory candidate analysis, integrating differential transcription-factor activity scores with upstream motif enrichment and downstream target gene expression, identified a set of candidate
  transcriptional regulators with differential activity across the PDAC subtype landscape. While initial analyses identified differential activity for several regulators, rigorous validation identified
  13 transcription factor activities as only partially replicated, and 9 activities that failed to replicate entirely in bulk data. Notably, the transcription factor CTCFL (BORIS), which initially showed
  differential activity, was rigorously evaluated for composition sensitivity. Single-cell decomposition revealed that its apparent differential expression was entirely explained by cell-type proportion
  variations (CELL_COMPOSITION_EXPLAINED) rather than malignant-intrinsic regulation. Consequently, CTCFL was explicitly penalized and removed from candidate prioritization, highlighting the necessity of
  single-cell resolution in preventing confounded regulatory inferences.
  ──────
  ### 2.6 Integrated Target Prioritisation: Phase 10B-R and Phase 10C2 Corrected Rankings

  Target prioritisation was performed using the corrected scoring framework established in Phases 10B-R and Phase 10C2. The framework assigned evidence weights based on: (i) cross-platform reproducibility
  of differential expression or activity signals; (ii) pathway membership in convergently enriched gene sets; (iii) malignant-compartment specificity of expression; (iv) availability of regulatory
  evidence (motif, epigenetic, or network-based); and (v) the absence of disqualifying evidence such as ubiquitous normal-tissue expression or essentiality in non-malignant cell types that would preclude
  selective targeting.

  Candidates were penalised where evidence derived from a single analytical platform, where spatial validation was absent or only partial, or where regulatory associations were motif-based only without
  corroborating functional data. Features whose single-cell evidence was classified as CELL_COMPOSITION_EXPLAINED received high-weight penalties that blocked advancement to prioritised status. The corrected
  Phase 10B-R ranking addressed errors in the preceding Phase 10B scoring step identified during quality review; Phase 10C2 incorporated additional pathway-compartment convergence evidence and applied
  revised penalties for partial spatial support.

  The top-ranked candidates emerging from the corrected framework shared the following properties: convergent evidence from at least two analytical levels (bulk, single-cell, and/or regulatory);
  malignant-compartment-enriched or malignant-compartment-restricted expression; membership in at least one hallmark gene set with malignant-cell-level enrichment support; and regulatory plausibility
  supported by motif or epigenetic evidence. The ranked list of prioritised candidates is provided in full in Supplementary Table S[X], stratified by evidence tier. Candidates in the highest evidence tier
  are those meeting all five scoring criteria with multi-platform replication. Candidates in intermediate tiers meet a subset of criteria and are presented as secondary hypotheses for experimental follow-
  up.

  HALLMARK_PROTEIN_SECRETION pathway constituent genes contributed disproportionately to candidates in the upper tier, consistent with the convergent malignant-cell and malignant-compartment evidence
  described in Section 2.3. As noted previously, CTCFL was not prioritized due to its composition sensitivity penalty. The Ochrobactrum-associated host transcriptional programmes contributed nominally to
  regulatory candidate scoring via their overlap with immune-modulatory gene sets, but this contribution was weighted conservatively given the absence of causal or spatial mechanistic support.
  ──────
  [End of Results Section Draft — References, supplementary table cross-references, and figure callouts to be inserted during assembly of the full manuscript.]
  ──────
  │ Compliance summary — Results: Microbial causality not claimed (§2.4). Microbial localisation not claimed (§2.4). Physical interaction not claimed (§2.4). Basal–classical spatial axis association of
  │ HALLMARK_PROTEIN_SECRETION explicitly classified as PARTIAL_SPATIAL_SUPPORT (§2.3). Null findings and partial findings explicitly reported (§2.3, §2.5). CTCFL/BORIS explicitly excluded as
  │ composition-sensitive, not prioritized (§2.5, §2.6). Ochrobactrum described in associative terms throughout (§2.4). Target ranking references corrected Phase 10B-R and Phase 10C2 outputs (§2.6).

