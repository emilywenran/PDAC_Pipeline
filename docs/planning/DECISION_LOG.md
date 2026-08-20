# Project Decision Log

This is an append-only log tracking all major design, bioinformatic, and statistical decisions made during the lifecycle of the project.

---

## Decision Matrix

| Log Date | Decision ID | Summary of Decision | Affected Files / Analyses |
| :--- | :--- | :--- | :--- |
| 2026-06-30 | **D-01** | Model subtype as a continuous axis (ssGSEA diff) rather than strictly discrete classes | `HYPOTHESES.md`, `ANALYSIS_PLAN_V0.1.md`, `TASK_DAG.md` |
| 2026-06-30 | **D-02** | Implement dual-method statistical decontamination (`decontam`) using negative controls | `ANALYSIS_PLAN_V0.1.md`, `TASK_DAG.md`, `RISK_REGISTER.md` |
| 2026-06-30 | **D-03** | Use composition-aware statistics (CLR transform, ALDEx2, ANCOM-BC, MaAsLin2) | `HYPOTHESES.md`, `ANALYSIS_PLAN_V0.1.md` |
| 2026-06-30 | **D-04** | Control for tumor purity (ESTIMATE scores) as a covariate in all association models | `ANALYSIS_PLAN_V0.1.md`, `RISK_REGISTER.md` |
| 2026-06-30 | **D-05** | Restrict Phase 1A to accession-level metadata audit and supplementary file downloads; defer raw FASTQ retrieval | `00_admin/PROJECT_STATUS.md`, `01_metadata/` |
| 2026-06-30 | **D-06** | Finalize 62-patient RNA/microbiome mapping and retain published survival status codes without recoding | `01_metadata/sample_manifest.tsv`, `01_metadata/clinical_metadata.tsv`, `01_metadata/file_manifest.tsv`, `01_metadata/rna_microbiome_patient_crosswalk.tsv`, `04_analysis/02_sample_mapping/PHASE1B_MAPPING_REPORT.md` |
| 2026-06-30 | **D-07** | Use the official GSE172356 processed matrix as the Phase 2A audited host expression source while preserving its DESeq size-factor-normalized count scale | `02_data/reference/GSE172356_processed/`, `03_processed/expression/`, `04_analysis/03_expression_qc/`, `05_results/tables/`, `05_results/figures/` |
| 2026-06-30 | **D-08** | Prepare Phase 2B expression matrices by retaining complete-observation genes, applying unsupervised expression filtering, and using log2(normalized count + 1) for analysis-ready subtype reproduction input | `03_processed/expression/`, `04_analysis/03_expression_qc/`, `05_results/tables/`, `05_results/figures/` |
| 2026-06-30 | **D-09** | Establish and lock primary (94-gene CSY subset) and secondary (Moffitt and PurIST) parameters for PDAC subtype reproduction | `01_metadata/subtype_method_inventory.tsv`, `02_data/reference/PDAC_subtype_signatures/`, `05_results/tables/phase3a_signature_gene_coverage.tsv`, `04_analysis/05_subtype_reproduction/PHASE3A_METHOD_LOCK.md` |
| 2026-07-01 | **D-10** | Execute locked Phase 3B subtype reproduction without parameter optimization; primary labels reproduce exactly and binary secondary methods report Hybrid samples separately | `05_results/tables/phase3b_*`, `05_results/figures/phase3b_*`, `04_analysis/05_subtype_reproduction/PHASE3B_SUBTYPE_REPRODUCTION.md`, `06_scripts/python/05_phase3b_*` |
| 2026-07-01 | **D-11** | Lock the statistical framework and parameters for molecular subtype stability analysis (Phase 4A) without running calculations | `04_analysis/06_subtype_stability/PHASE4A_STABILITY_METHOD_LOCK.md`, `09_docs/methods/PDAC_subtype_stability_protocol.md`, `01_metadata/subtype_stability_parameter_inventory.tsv` |
| 2026-07-01 | **D-12** | Apply locked Phase 4B subtype stability decision rules and classify aggregate evidence as INCONCLUSIVE | `04_analysis/06_subtype_stability/PHASE4B_SUBTYPE_STABILITY_RESULTS.md`, `05_results/tables/phase4b_*`, `05_results/figures/phase4b_*` |
| 2026-07-01 | **D-13** | Lock Phase 5A continuous basal-classical transcriptional axis scoring, hybrid metrics, sample categories, and decision rules | `04_analysis/07_continuous_subtype_axis/PHASE5A_AXIS_METHOD_LOCK.md`, `09_docs/methods/PDAC_continuous_axis_protocol.md`, `01_metadata/continuous_axis_parameter_inventory.tsv` |
| 2026-07-01 | **D-14** | Amend Phase 5A lock for Moffitt gene-set reconciliation (50-gene primary, 49-gene sensitivity) | `04_analysis/07_continuous_subtype_axis/PHASE5A_AXIS_METHOD_LOCK.md`, `09_docs/methods/PDAC_continuous_axis_protocol.md`, `01_metadata/continuous_axis_parameter_inventory.tsv`, `04_analysis/07_continuous_subtype_axis/PHASE5A_GENESET_RECONCILIATION.md` |
| 2026-07-01 | **D-15** | Execute locked Phase 5B continuous basal-classical axis scoring and trend tests (overall decision INCONCLUSIVE) | `04_analysis/07_continuous_subtype_axis/PHASE5B_CONTINUOUS_AXIS_RESULTS.md`, `05_results/tables/phase5b_*`, `05_results/figures/phase5b_*` |
| 2026-07-01 | **D-16** | Lock Phase 6A processed PRJNA719915 tumor microbiome abundance source and audit scope | `03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz`, `01_metadata/microbiome_sample_crosswalk.tsv`, `04_analysis/04_microbiome_qc/PHASE6A_MICROBIOME_DATA_AUDIT.md` |
| 2026-07-01 | **D-17** | Lock Phase 6B tumor microbiome preprocessing filters, CLR/rCLR transformations, and contamination-sensitivity rules | `04_analysis/04_microbiome_qc/PHASE6B_MICROBIOME_METHOD_LOCK.md`, `09_docs/methods/PDAC_microbiome_preprocessing_protocol.md`, `01_metadata/microbiome_preprocessing_parameter_inventory.tsv` |
| 2026-07-01 | **D-18** | Execute locked Phase 6C microbiome preprocessing and generate primary CLR and sensitivity representations | `03_processed/microbiome/PRJNA719915_genus_primary_CLR.tsv.gz`, `03_processed/microbiome/sensitivity/`, `04_analysis/04_microbiome_qc/PHASE6C_ANALYSIS_READY_MICROBIOME.md` |
| 2026-07-01 | **D-19** | Lock Phase 7A statistical framework for microbiome-host associations under covariate limitations | `04_analysis/08_host_microbiome_integration/PHASE7A_MICROBIOME_ASSOCIATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_microbiome_continuous_state_association_protocol.md`, `01_metadata/microbiome_association_parameter_inventory.tsv`, `05_results/tables/phase7a_model_matrix_feasibility.tsv` |
| 2026-07-01 | **D-20** | Validate ESTIMATE-derived host TME covariates and add feasible sensitivity models Model 3P, Model 3I, and Model 3S | `01_metadata/host_tme_covariates.tsv`, `05_results/tables/phase7a5_host_covariate_qc.tsv`, `05_results/tables/phase7a5_host_covariate_correlations.tsv`, `05_results/tables/phase7a5_covariate_model_feasibility.tsv`, `05_results/figures/phase7a5_host_covariate_distributions.pdf`, `05_results/figures/phase7a5_host_covariate_correlation.pdf`, `04_analysis/08_host_microbiome_integration/PHASE7A5_HOST_COVARIATES.md` |
| 2026-07-01 | **D-21** | Execute locked Phase 7B continuous tumor microbiome association models and generate primary and sensitivity outputs | `05_results/tables/phase7b_*`, `05_results/figures/phase7b_*`, `06_scripts/python/11_summarize_phase7b_associations.py`, `04_analysis/08_host_microbiome_integration/PHASE7B_MICROBIOME_ASSOCIATION_RESULTS.md` |
| 2026-07-02 | **D-22** | Complete Phase 7C independent review of Phase 7B microbiome-host integration results and report PASS decision | `04_analysis/08_host_microbiome_integration/PHASE7C_INDEPENDENT_REVIEW.md`, `05_results/tables/phase7c_*` |
| 2026-07-02 | **D-23** | Prespecify and lock host-mechanism analysis framework for tumor-microbiome associations | `04_analysis/08_host_microbiome_integration/PHASE8A_HOST_MECHANISM_METHOD_LOCK.md`, `09_docs/methods/PDAC_host_microbiome_mechanism_protocol.md`, `01_metadata/host_mechanism_parameter_inventory.tsv`, `05_results/tables/phase8a_host_feature_feasibility.tsv` |
| 2026-07-02 | **D-24** | Prepare and validate Phase 8B R environment before host-mechanism execution | `renv.lock`, `04_analysis/08_host_microbiome_integration/PHASE8A5_ENVIRONMENT_VALIDATION.md`, `01_metadata/file_manifest.tsv` |
| 2026-07-02 | **D-25** | Execute locked Phase 8B host-mechanism analyses | `06_scripts/R/13_phase8b_host_mechanisms.R`, `04_analysis/08_host_microbiome_integration/PHASE8B_HOST_MECHANISM_RESULTS.md`, `05_results/tables/phase8b_*` |
| 2026-07-02 | **D-26** | Complete Phase 8C independent review of Phase 8B host-microbiome mechanism results and report PASS decision | `04_analysis/08_host_microbiome_integration/PHASE8C_INDEPENDENT_HOST_MECHANISM_REVIEW.md`, `05_results/tables/phase8c_*` |
| 2026-07-02 | **D-27** | Identify, evaluate, and lock the external-validation framework (Phase 9A) across four layers under a prospective method lock | `04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_external_validation_protocol.md`, `01_metadata/external_validation_*`, `05_results/tables/phase9a_*` |
| 2026-07-03 | **D-28** | Execute Locked Phase 9B1 Independent Bulk-Transcriptome Validation | `06_scripts/R/14_phase9b1_bulk_validation.R`, `04_analysis/09_external_validation/PHASE9B1_BULK_EXTERNAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b1_*` |
| 2026-07-03 | **D-29** | Reject Phase 9B1 Bulk Validation Outputs due to Critical and Major Implementation Errors (Phase 9B1C) | `04_analysis/09_external_validation/PHASE9B1C_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c_*` |
| 2026-07-03 | **D-30** | Correct and Rerun Phase 9B1 Bulk Validation After Phase 9B1C Audit (Phase 9B1R) | `06_scripts/R/14_phase9b1r_corrected_bulk_validation.R`, `04_analysis/09_external_validation/PHASE9B1R_CORRECTED_BULK_EXTERNAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b1r_*` |
| 2026-07-03 | **D-31** | Complete Phase 9B1C2 independent review of the corrected Phase 9B1R bulk external validation results and report PASS_WITH_MINOR_CORRECTIONS decision | `04_analysis/09_external_validation/PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c2_*` |
| 2026-07-03 | **D-28** | Execute Locked Phase 9B1 Independent Bulk-Transcriptome Validation | `06_scripts/R/14_phase9b1_bulk_validation.R`, `04_analysis/09_external_validation/PHASE9B1_BULK_EXTERNAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b1_*` |
| 2026-07-03 | **D-29** | Reject Phase 9B1 Bulk Validation Outputs due to Critical and Major Implementation Errors (Phase 9B1C) | `04_analysis/09_external_validation/PHASE9B1C_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c_*` |
| 2026-07-03 | **D-30** | Correct and Rerun Phase 9B1 Bulk Validation After Phase 9B1C Audit (Phase 9B1R) | `06_scripts/R/14_phase9b1r_corrected_bulk_validation.R`, `04_analysis/09_external_validation/PHASE9B1R_CORRECTED_BULK_EXTERNAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b1r_*` |
| 2026-07-03 | **D-31** | Complete Phase 9B1C2 independent review of the corrected Phase 9B1R bulk external validation results and report PASS_WITH_MINOR_CORRECTIONS decision | `04_analysis/09_external_validation/PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c2_*` |
| 2026-07-03 | **D-32** | Final Closure of Phase 9B1C2 Minor Correction and Verification | `06_scripts/R/14_phase9b1r_corrected_bulk_validation.R`, `06_scripts/python/14_validate_phase9b1r_bulk_validation.py`, `05_results/tables/phase9b1r_host_feature_replication_evidence.tsv` |
| 2026-07-03 | **D-33** | Stop Phase 9B2 Before Data Acquisition Because Locked Single-Cell Records Disagree | `01_metadata/phase9b2_single_cell_dataset_inventory.tsv`, `04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md` |
| 2026-07-03 | **D-34** | Reconcile Layer 2 Single-Cell Cohort Set and Establish Sole Authoritative Dataset (Phase 9A.1) | `01_metadata/external_validation_dataset_inventory.tsv`, `05_results/tables/phase9a_external_dataset_shortlist.tsv`, `01_metadata/external_validation_parameter_inventory.tsv`, `04_analysis/09_external_validation/PHASE9A1_SINGLE_CELL_COHORT_RECONCILIATION.md` |
| 2026-07-03 | **D-35** | Reconcile single-cell dataset accessions and provenance, expand Phase 9B2 execution cohort set (Phase 9A.2) | `01_metadata/external_validation_dataset_inventory.tsv`, `05_results/tables/phase9a_external_dataset_shortlist.tsv`, `01_metadata/external_validation_parameter_inventory.tsv`, `05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv`, `04_analysis/09_external_validation/PHASE9A2_SINGLE_CELL_DATASET_PROVENANCE_CORRECTION.md` |
| 2026-07-03 | **D-36** | Stop Phase 9B2 Restart Because Authoritative Included Set Exceeds Primary PENG-Only Contract | `05_results/tables/phase9b2_restart_runtime_validation.tsv`, `04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md` |
| 2026-07-03 | **D-37** | Phase 9A.3 execution-scope correction for Phase 9B2 | `05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv`, `01_metadata/external_validation_parameter_inventory.tsv`, `06_scripts/python/15_*`, `05_results/tables/phase9a3_phase9b2_execution_scope.tsv` |
| 2026-07-03 | **D-38** | Execute Phase 9B2 Primary Single-Cell Cellular-Source Analysis on PENG_CRA001160 Only | `02_data/external/phase9_single_cell/PENG_CRA001160/`, `05_results/tables/phase9b2_*`, `06_scripts/R/15_phase9b2_single_cell_validation.R` |
| 2026-07-03 | **D-39** | Reject Phase 9B2 Single-Cell Cellular-Source Analysis due to Incomplete Negative Controls and Coverage Violations (Phase 9B2C) | `04_analysis/09_external_validation/PHASE9B2C_SINGLE_CELL_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b2c_*` |
| 2026-07-03 | **D-40** | Complete Phase 9B2R Corrective Single-Cell Reanalysis and Supersede Initial Phase 9B2 Results | `04_analysis/09_external_validation/PHASE9B2R_CORRECTED_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md`, `05_results/tables/phase9b2r_*`, `05_results/figures/phase9b2r_*`, `06_scripts/*/15_*phase9b2r*` |
| 2026-07-03 | **D-41** | Issue Final PASS Decision for Phase 9B2C2 Independent Review of Corrected Single-Cell Validation | `04_analysis/09_external_validation/PHASE9B2C2_CORRECTED_SINGLE_CELL_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b2c2_*`, `09_docs/planning/DECISION_LOG.md` |
| 2026-07-03 | **D-42** | Lock spatial feature hierarchy, qualification, cohort set, and statistical design for Phase 9B3 spatial validation | `04_analysis/09_external_validation/PHASE9B3A_SPATIAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_spatial_validation_protocol.md`, `01_metadata/phase9b3_spatial_*`, `05_results/tables/phase9b3a_*` |
| 2026-07-03 | **D-43** | Correct spatial validation covariates, nesting random-effect structures, counts, and terminology (Phase 9B3A.1) | `04_analysis/09_external_validation/PHASE9B3A1_SPATIAL_DESIGN_CONSISTENCY_CORRECTION.md`, `05_results/tables/phase9b3a1_spatial_analysis_unit_and_models.tsv`, `01_metadata/phase9b3_spatial_*` |
| 2026-07-03 | **D-44** | Lock spatial ROI pairing models, Moncada exploratory cross-platform consistency role, and matrix pooling ban (Phase 9B3A.2) | `04_analysis/09_external_validation/PHASE9B3A2_SPATIAL_HIERARCHY_FINAL_CORRECTION.md`, `05_results/tables/phase9b3a2_spatial_model_hierarchy.tsv`, `01_metadata/phase9b3_spatial_*` |
| 2026-07-03 | **D-45** | Execute Phase 9B3B spatial-transcriptomic validation on Hwang and Moncada cohorts | `04_analysis/09_external_validation/PHASE9B3B_SPATIAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b3b_*` |
| 2026-07-04 | **D-46** | Reject Phase 9B3B spatial validation due to hardcoded negative controls and model violations (Phase 9B3C) | `04_analysis/09_external_validation/PHASE9B3C_SPATIAL_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b3c_*` |
| 2026-07-04 | **D-47** | Complete Phase 11D Full Manuscript Assembly adhering to Phase 11C constraints | `04_analysis/11_manuscript/PHASE11D_FULL_MANUSCRIPT_DRAFT.md`, `05_results/tables/phase11d_*`, `06_scripts/python/19_validate_phase11d_full_manuscript.py` |
| 2026-07-04 | **D-53** | Complete Phase 11E Language and Format Review | `04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md`, `04_analysis/11_manuscript/PHASE11E_LANGUAGE_FORMAT_REVIEW.md`, `05_results/tables/phase11e_language_edit_log.tsv`, `06_scripts/python/19_validate_phase11e_language_format.py` |
| 2026-07-04 | **D-54** | Complete Phase 11F Final Claim Audit | `04_analysis/11_manuscript/PHASE11F_FINAL_CLAIM_AUDIT.md`, `05_results/tables/phase11f_*`, `06_scripts/python/19_validate_phase11f_final_claim_audit.py` |
| 2026-07-04 | **D-55** | Complete Phase 11G Administrative Finalization | `04_analysis/11_manuscript/PHASE11G*.md`, `05_results/tables/phase11g*.tsv` |
| 2026-07-04 | **D-56** | Complete Post-Phase 11G Workspace Cleanup Finalization | `.gitignore`, `04_analysis/11_manuscript/POST_PHASE11G_WORKSPACE_CLEANUP_AUDIT.md`, `01_metadata/file_manifest.tsv`, `00_admin/PROJECT_STATUS.md`, `09_docs/planning/DECISION_LOG.md` |

---

## Detailed Decision Logs

### D-01: Continuous Subtype Axis Modeling
*   **Date:** 2026-06-30
*   **Decision:** Model host transcriptomics along a continuous basal–classical axis (using ssGSEA or NMF projection score difference) in all primary association testing, using discrete $k=3$ consensus clustering solely for reproducing Guo et al. (2021) and evaluating stability.
*   **Alternatives Considered:** Model subtypes strictly as three discrete classes (classical, basal-like, hybrid) as originally reported.
*   **Scientific Justification:** Forcing samples into discrete bins decreases statistical power and misrepresents the biological reality of PDAC transcriptomic plasticity. If subtypes represent a continuum, a discrete "hybrid" group is an arbitrary partition of intermediate samples. Continuous scoring prevents information loss and allows for more robust regression models.
*   **Files / Analyses Affected:** `09_docs/planning/HYPOTHESES.md`, `09_docs/planning/ANALYSIS_PLAN_V0.1.md`, and downstream scripts in `04_analysis/06_subtype_stability/` and `04_analysis/07_continuous_subtype_axis/`.

---

### D-02: Dual-Method Statistical Decontamination (`decontam`)
*   **Date:** 2026-06-30
*   **Decision:** Apply the R package `decontam` using both **prevalence** (comparing tumor samples to blanks) and **frequency** (correlating feature abundance with DNA concentration) methods. Do not rely solely on simple abundance thresholds.
*   **Alternatives Considered:** Rely on a simple relative abundance cutoff (e.g., discard features with $< 0.01\%$ abundance) or visual inspection of raw tables.
*   **Scientific Justification:** Intratumoral microbiomes have low bacterial biomass, making them highly vulnerable to reagent and laboratory contamination. Abundance filters can fail by retaining high-abundance contaminants or discarding rare, biologically relevant taxa. Statistical comparison with negative controls is the benchmark standard in low-biomass microbiomics.
*   **Files / Analyses Affected:** `09_docs/planning/ANALYSIS_PLAN_V0.1.md`, `09_docs/planning/TASK_DAG.md`, `09_docs/planning/RISK_REGISTER.md`, and downstream scripts in `04_analysis/04_microbiome_qc/`.

---

### D-03: Composition-Aware Statistical Framework
*   **Date:** 2026-06-30
*   **Decision:** Utilize Centered Log-Ratio (CLR) transformations for UMAP/PCA and apply composition-aware models (ALDEx2, ANCOM-BC, MaAsLin2) for all differential abundance and association analyses.
*   **Alternatives Considered:** Apply standard Pearson or Spearman correlation and ordinary least squares (OLS) regression on raw relative abundances (proportions).
*   **Scientific Justification:** Microbiome datasets are compositional because they are constrained by total library size. Applying standard statistics to proportions causes mathematical dependency (sub-compositional coherence issues), which yields false correlations and inflated type I errors. CLR-based and log-ratio methods project the data into unconstrained Euclidean space, ensuring valid regression.
*   **Files / Analyses Affected:** `09_docs/planning/HYPOTHESES.md`, `09_docs/planning/ANALYSIS_PLAN_V0.1.md`, and downstream association testing in `04_analysis/08_host_microbiome_integration/`.

---

### D-04: Controlling for Tumor Purity Confounding
*   **Date:** 2026-06-30
*   **Decision:** Estimate tumor purity for all bulk RNA-seq samples using the ESTIMATE algorithm and include the purity score as a covariate in all MaAsLin2 and ANCOM-BC regressions.
*   **Alternatives Considered:** Run association tests without adjusting for cellular composition.
*   **Scientific Justification:** Bulk tumor expression is a mixture of tumor, stromal, and immune cell transcripts. Basal-like signatures are highly enriched in stromal and inflammatory pathways. If tumor purity correlates with microbial abundance (e.g., due to immune infiltration reacting to bacteria), associations might be driven by stroma rather than tumor cells. Controlling for purity ensures tumor-specific signals are isolated.
*   **Files / Analyses Affected:** `09_docs/planning/ANALYSIS_PLAN_V0.1.md`, `09_docs/planning/RISK_REGISTER.md`, and downstream scripts in `04_analysis/08_host_microbiome_integration/`.

---

### D-05: Restrict Phase 1A to Accession Audit and Supplementary Data Retrieval
*   **Date:** 2026-06-30
*   **Decision:** Limit Phase 1A execution to accession-level metadata audit, sample/run count reconciliation, patient-level mapping verification, and download of publication supplementary materials. Defer the retrieval or download of any raw sequencing files (FASTQ, SRA archives) or BAM files.
*   **Alternatives Considered:** Downloading or pre-fetching raw FASTQ files from NCBI SRA or ENA during the mapping phase.
*   **Scientific and Operational Justification:** Large-scale raw data downloads (approx. ~1.5 TB for metagenomic and transcriptomic reads) are unnecessary and wasteful before the sample mapping and sample manifests are verified. Performing an accession-level metadata audit allows us to establish the exact patient mapping (1-to-1 matching via `source_material_id` and tumor number) and confirm the library strategy (confirming WGS shotgun metagenomics for microbiome) without downloading massive datasets. Furthermore, we identified that the study's negative controls were only assayed via PCR/gel and not sequenced, so no sequencing reads are available. This demonstrates that auditing accession-level metadata is sufficient to construct the complete cohort manifest.
*   **Files / Analyses Affected:** `01_metadata/accession_inventory.tsv`, `01_metadata/geo_sample_inventory.tsv`, `01_metadata/microbiome_run_inventory.tsv`, `00_admin/PROJECT_STATUS.md`.

---

### D-06: Finalize Patient-Level RNA/Microbiome Mapping
*   **Date:** 2026-06-30
*   **Decision:** Use the Phase 1A verified mapping between GEO tumor number and microbiome `source_material_id` / `potential_patient_id` as the authoritative patient-level crosswalk for all 62 PDAC patients. Populate subtype labels from Guo et al. 2021 Supplementary Data 4 (`Figure1.SampleGroup`) and retain overall survival `Days` and `status` exactly as published in `Figure3.Survival` without recoding event status in Phase 1B.
*   **Alternatives Considered:** Reconstruct the mapping by downloading raw FASTQ/SRA files, or recode survival status values during manifest construction.
*   **Scientific and Operational Justification:** The Phase 1A public metadata already verified a one-to-one mapping between each host RNA-seq sample and each tumor microbiome run, with no duplicated GEO, BioSample, experiment, or run IDs. Raw sequencing downloads are unnecessary for patient manifest construction. Survival modeling is out of scope for Phase 1B, so published survival status codes are preserved as source data until a later phase explicitly confirms event coding.
*   **Files / Analyses Affected:** `01_metadata/sample_manifest.tsv`, `01_metadata/clinical_metadata.tsv`, `01_metadata/file_manifest.tsv`, `01_metadata/rna_microbiome_patient_crosswalk.tsv`, `04_analysis/02_sample_mapping/PHASE1B_MAPPING_REPORT.md`, and `00_admin/PROJECT_STATUS.md`.

---

### D-07: Phase 2A Processed Host Expression Matrix Source and Scale
*   **Date:** 2026-06-30
*   **Decision:** Use the official GEO Series supplementary file `GSE172356_PDA_gene_expression_matrix.txt.gz` as the Phase 2A processed host gene-expression source. Preserve the original numerical scale in `03_processed/expression/GSE172356_expression_audited.tsv.gz` and document the unit as DESeq size-factor-normalized counts based on official GEO SOFT processing metadata.
*   **Alternatives Considered:** Download raw RNA-seq FASTQ/SRA files and reconstruct counts; transform the processed matrix during audit; drop genes or samples with missing values or suspected outlier behavior.
*   **Scientific and Operational Justification:** Phase 2A is an acquisition and descriptive audit phase, not a reprocessing, normalization, subtype discovery, or differential-expression phase. The official GEO processed matrix provides 62 expression columns that map one-to-one to the finalized patient manifest via `YX...T` aliases. The audit found 45,140 gene rows, 62 mapped samples, no duplicate gene or sample identifiers, no negative or infinite values, and no unmapped expression columns. Missing values and suspected extreme samples are documented but not removed. Downstream transformations, if needed, must be explicitly scripted in later phases.
*   **Files / Analyses Affected:** `00_admin/SKILL_USAGE_LOG.tsv`, `01_metadata/file_manifest.tsv`, `01_metadata/expression_sample_crosswalk.tsv`, `02_data/reference/GSE172356_processed/GSE172356_PDA_gene_expression_matrix.txt.gz`, `03_processed/expression/GSE172356_expression_audited.tsv.gz`, `03_processed/expression/GSE172356_gene_annotation.tsv`, `04_analysis/03_expression_qc/PHASE2A_EXPRESSION_AUDIT.md`, `05_results/tables/phase2a_expression_sample_qc.tsv`, `05_results/tables/phase2a_expression_gene_qc.tsv`, `05_results/tables/phase2a_expression_mapping_summary.tsv`, `05_results/figures/phase2a_expression_distribution.pdf`, `05_results/figures/phase2a_sample_correlation.pdf`, and `05_results/figures/phase2a_expression_pca.pdf`.

### D-08: Phase 2B Missing-Value Handling, Filtering, and Transformation
*   **Date:** 2026-06-30
*   **Decision:** Use complete-observation genes as the primary missing-value strategy for GSE172356 Phase 2B, then remove all-zero genes and retain genes with DESeq size-factor-normalized count >= 1 in at least 10% of the 62 mapped samples. Preserve all 62 samples. Write both filtered normalized counts and `log2(normalized count + 1)` matrices. The selected primary matrix contains 42654 genes.
*   **Alternatives Considered:** Retain all genes with missing values and perform gene-median imputation; retain genes with <=50% missingness and impute; replace `NA` with zero; apply DESeq2 VST/rlog directly to the processed matrix. These were rejected for the primary matrix because the source uses literal `NA` without documentation that the entries are structural zeros, and VST/rlog require raw-count assumptions that are not satisfied by an already normalized, fractional matrix.
*   **Scientific and Operational Justification:** Complete-observation filtering avoids imputing unexplained source `NA` values before subtype reproduction. The expression filter is independent of subtype labels and removes genes unlikely to contribute stable unsupervised structure. The log2 transform is reproducible for non-negative normalized counts and avoids applying a second library-size normalization.
*   **Files / Analyses Affected:** `03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz`, `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz`, `04_analysis/03_expression_qc/PHASE2B_ANALYSIS_READY_EXPRESSION.md`, `05_results/tables/phase2b_*`, and `05_results/figures/phase2b_*`.

---

## Revisions and Corrections Log

### R-LOG-01: Shift from 16S Amplicon to Shotgun Metagenomic Workflow and Refined Validation Standards
*   **Date:** 2026-06-30
*   **Author:** Antigravity (AI Coding Assistant)
*   **Description of Correction:** BioProject PRJNA719915 was identified to contain shotgun/metagenomic tumor sequencing data rather than 16S amplicon sequencing data. Consequently, references to 16S rRNA hypervariable regions, DADA2 processing, ASVs, ASV denoising, and ASV abundance tables were removed. They were replaced with a shotgun metagenomic processing workflow covering raw-read QC, host-read depletion (alignment to GRCh38), taxonomic profiling (with Kraken2/Bracken and/or MetaPhlAn), contamination assessment (using negative controls), prevalence and abundance filtering, compositional data analysis, microbial functional profiling (with HUMAnN where supported), and taxonomic sensitivity analysis. Additionally, rules for external validation using TCGA-PAAD were corrected to emphasize host transcriptomic validation and treat TCGA-derived microbiome signals as strictly exploratory due to contamination and batch concerns, without assuming the existence of an identical external paired cohort.
*   **Affected Planning Files:**
    *   [PROJECT_CHARTER.md](file://~/thesis/PDAC/09_docs/planning/PROJECT_CHARTER.md)
    *   [ANALYSIS_PLAN_V0.1.md](file://~/thesis/PDAC/09_docs/planning/ANALYSIS_PLAN_V0.1.md)
    *   [EVIDENCE_POLICY.md](file://~/thesis/PDAC/09_docs/planning/EVIDENCE_POLICY.md)
    *   [HYPOTHESES.md](file://~/thesis/PDAC/09_docs/planning/HYPOTHESES.md)
    *   [MANUSCRIPT_GUARDRAILS.md](file://~/thesis/PDAC/09_docs/planning/MANUSCRIPT_GUARDRAILS.md)
    *   [RISK_REGISTER.md](file://~/thesis/PDAC/09_docs/planning/RISK_REGISTER.md)
    *   [TASK_DAG.md](file://~/thesis/PDAC/09_docs/planning/TASK_DAG.md)
    *   [DECISION_LOG.md](file://~/thesis/PDAC/09_docs/planning/DECISION_LOG.md)

---

### D-09: Lock PDAC Subtype Reproduction Parameters
*   **Date:** 2026-06-30
*   **Decision:** Establish and lock the primary hierarchical clustering parameters (Pearson correlation distance, average linkage, column reordering slicing of size-factor normalized counts for 94 Chan-Seng-Yue genes) and secondary comparison models (49-gene Moffitt-derived clustering and 16-gene PurIST single-sample classifier) before starting Phase 3B.
*   **Alternatives Considered:** Rely on a single classification method (e.g., only primary or only PurIST); use standard Euclidean distance; include Bailey as a primary framework.
*   **Scientific and Operational Justification:** Pancreatic cancer subtypes are highly sensitive to preprocessing, gene set selection, and clustering options. Locking a primary hierarchical reproduction method guarantees that we can exactly reproduce the published cohort group allocations. Evaluating Moffitt and PurIST in parallel allows us to assess the degree to which these assignments represent robust biological signals versus platform- or study-specific artifacts. Bailey is relegated to an exploratory role because it lacks a validated single-sample classifier and is highly susceptible to stroma/normal tissue contamination.
*   **Files / Analyses Affected:** `01_metadata/subtype_method_inventory.tsv`, `02_data/reference/PDAC_subtype_signatures/`, `05_results/tables/phase3a_signature_gene_coverage.tsv`, `04_analysis/05_subtype_reproduction/PHASE3A_METHOD_LOCK.md`, and `09_docs/methods/PDAC_subtype_reproduction_protocol.md`.

### D-10: Execute Locked Phase 3B Subtype Reproduction
*   **Date:** 2026-07-01
*   **Decision:** Apply the Phase 3A-locked primary GSE172356/Chan-Seng-Yue 94-gene hierarchical clustering procedure and verified secondary Moffitt and PurIST methods without altering signatures, preprocessing rules, coefficients, thresholds, gene-pair directions, or slice sizes after reviewing agreement results.
*   **Alternatives Considered:** Tune dendrogram orientation, class boundaries, confidence thresholds, or gene inclusion to improve agreement; force public Hybrid samples into binary basal/classical errors for Moffitt or PurIST comparisons; execute Bailey or the full Chan-Seng-Yue 100-gene exploratory framework as validated classifiers.
*   **Scientific and Operational Justification:** The primary method reproduced the verified public labels exactly for all 62 patients, so no optimization or adaptation is justified. Binary methods are biologically informative comparators but do not define a Hybrid class; public Hybrid samples must therefore be summarized separately along the basal-classical axis rather than treated as automatic classification failures. Bailey and the full Chan-Seng-Yue exploratory framework remain `TO_VERIFY` because they lack locked single-sample assignment rules in this project.
*   **Files / Analyses Affected:** `05_results/tables/phase3b_signature_runtime_validation.tsv`, `05_results/tables/phase3b_primary_subtype_assignments.tsv`, `05_results/tables/phase3b_all_method_assignments.tsv`, `05_results/tables/phase3b_method_agreement_metrics.tsv`, `05_results/tables/phase3b_confusion_matrices.tsv`, `05_results/tables/phase3b_discordant_samples.tsv`, `05_results/tables/phase3b_sensitivity_summary.tsv`, `05_results/figures/phase3b_*`, `04_analysis/05_subtype_reproduction/PHASE3B_SUBTYPE_REPRODUCTION.md`, `06_scripts/python/05_phase3b_reproduce_subtypes.py`, and `06_scripts/python/05_validate_phase3b_subtypes.py`.

---

### D-11: Lock Statistical Framework for Molecular Subtype Stability Analysis
*   **Date:** 2026-07-01
*   **Decision:** Establish and lock the cohort, cluster, and sample-level stability metrics (using 1,000 iterations of sample/feature resampling), independent high-variance gene representation, sample-level hybrid assignment interpretations, and multi-statistic discrete cluster decision rules (incorporating the Phase 3B log2-scale sensitivity as a baseline driver for stability testing) prior to running Phase 4B computations.
*   **Alternatives Considered:** Rely on single metrics like PAC or Jaccard stability; run stability calculations incrementally; omit the independent unsupervised HVG benchmark; select preprocessing according to public label reproduction.
*   **Scientific and Operational Justification:** In biology, the discrete or continuous nature of subtypes is highly debated. Defining a multi-statistic decision rule prevents selection bias, while checking stability on both raw and log2 count scales addresses the volatility observed in Phase 3B. Unsupervised HVG evaluation provides an unbiased comparison that does not rely on predefined subtyping signature genes.
*   **Files / Analyses Affected:** `04_analysis/06_subtype_stability/PHASE4A_STABILITY_METHOD_LOCK.md`, `09_docs/methods/PDAC_subtype_stability_protocol.md`, and `01_metadata/subtype_stability_parameter_inventory.tsv`.

### D-12: Apply Locked Phase 4B Subtype Stability Decision Rules
*   **Date:** 2026-07-01
*   **Decision:** Execute all eight locked Phase 4B subtype-stability analyses and classify the aggregate evidence as `INCONCLUSIVE`.
*   **Alternatives Considered:** Select K using only public-label agreement; tune preprocessing after inspecting stability metrics; treat the primary K=3 label reproduction as sufficient evidence of discrete stability.
*   **Scientific and Operational Justification:** The locked metrics separated stability evidence from post-clustering label agreement. Primary CSY K=3 reproduced public labels descriptively, but K=2 was preferred by the locked multi-metric stability rank, independent HVG analysis preferred K=4, and preprocessing/HVG sensitivities showed material instability.
*   **Files / Analyses Affected:** `04_analysis/06_subtype_stability/PHASE4B_SUBTYPE_STABILITY_RESULTS.md`, `05_results/tables/phase4b_*`, `05_results/figures/phase4b_*`, `06_scripts/R/06_phase4b_subtype_stability.R`, `06_scripts/python/06_summarize_phase4b_stability.py`, and `06_scripts/python/06_validate_phase4b_stability.py`.

---

### D-13: Lock Continuous Basal–Classical Transcriptional-Axis Analytical Plan
*   **Date:** 2026-07-01
*   **Decision:** Establish and lock continuous measures (Basal, Classical, Contrast, and Co-activation scores, plus centroid distances and clinical probabilities), scoring approaches (Primary Moffitt-mean and Secondary R-singscore rank-based), reference centroid leakage control definitions, hybrid-state characterization metrics, predefined interpretation categories, statistical evaluations (JT trend, multimodality dip test, bootstrap CIs, permutations), and multi-metric decision rules prior to executing Phase 5B calculations.
*   **Alternatives Considered:** Rely solely on discrete subtyping classifications; use only a single continuous score (e.g., contrast score only); use a simple cutoff boundary optimized against public labels to define hybrid samples.
*   **Scientific and Operational Justification:** Since the Phase 4B stability results were inconclusive and showed substantial transformation volatility, treating subtypes as continuous coordinates rather than discrete clusters represents a more mathematically robust model of cellular states. By locking two independent scoring systems (parametric mean vs. non-parametric rank-based singscore) and defining objective multi-metric rules for final classification, we eliminate post hoc cherry-picking and guarantee an unbiased assessment of whether public hybrid samples represent a linear continuum, a co-activated hybrid state, or method-dependent noise.
*   **Files / Analyses Affected:** `04_analysis/07_continuous_subtype_axis/PHASE5A_AXIS_METHOD_LOCK.md`, `09_docs/methods/PDAC_continuous_axis_protocol.md`, `01_metadata/continuous_axis_parameter_inventory.tsv`.

---

### D-14: Amend Phase 5A Lock for Moffitt Gene-Set Reconciliation
*   **Date:** 2026-07-01
*   **Decision:** Reconcile and lock the Moffitt signature definitions. Define the **Primary Analysis** using the canonical 50-gene signature (25 Basal-like, 25 Classical, retaining `LEMD1`) and the **Prespecified Sensitivity Analysis** using the 49-gene signature (24 Basal-like, 25 Classical, excluding `LEMD1`). Update the continuous axis method lock, protocol, and parameter inventory (`AXIS_MOFFITT50_PRIMARY` and `AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY`). 
*   **Alternatives Considered:** Proceed with the mismatched 49-gene active signature defined as 25 basal and 24 classical; proceed with only the 49-gene signature and drop the 50-gene signature.
*   **Scientific and Operational Justification:** Programmatic and literature audits revealed that the verified reference Moffitt signature has 50 genes (25 basal-like, 25 classical), where `LEMD1` is basal-like. The original lock documents incorrectly described the 49-gene active signature as having 25 basal and 24 classical genes. Because `LEMD1` belongs to the Basal-like program, its exclusion mathematically yields 24 Basal-like and 25 Classical genes. Retaining the 50-gene signature as primary aligns with the canonical literature, and keeping the 49-gene signature as sensitivity ensures comparability with works that excluded `LEMD1` due to cohort-specific expression properties, while correcting the count specifications. Execution of Phase 5B was stopped to address this inconsistency before any scores were calculated. It is confirmed that no Phase 5B continuous scoring results, downstream evaluations, or results had been generated prior to this amendment.
*   **Files / Analyses Affected:** `04_analysis/07_continuous_subtype_axis/PHASE5A_AXIS_METHOD_LOCK.md`, `09_docs/methods/PDAC_continuous_axis_protocol.md`, `01_metadata/continuous_axis_parameter_inventory.tsv`, `04_analysis/07_continuous_subtype_axis/PHASE5A_GENESET_RECONCILIATION.md`, and signature files under `02_data/reference/PDAC_subtype_signatures/`.

### D-15: Execute Locked Phase 5B Continuous Axis Analysis
*   **Date:** 2026-07-01
*   **Decision:** Execute all seven locked Phase 5B continuous basal-classical transcriptional-axis analysis IDs after confirming that the amended method lock, protocol, parameter inventory, signature files, and reconciliation report agree. Apply the locked multi-metric decision rules and report the overall Phase 5B decision as `INCONCLUSIVE`.
*   **Alternatives Considered:** Substitute the 49-gene no-LEMD1 signature for the 50-gene primary signature; optimize a Hybrid cutoff using public labels; alter thresholds after viewing score distributions; skip rank-based or centroid sensitivity analyses because local R packages were unavailable.
*   **Scientific and Operational Justification:** Programmatic pre-execution verification confirmed that `AXIS_MOFFITT50_PRIMARY` uses the 50-gene Moffitt signature with 25 Basal-like and 25 Classical genes including `LEMD1`, while `AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY` uses the 49-gene sensitivity signature with 24 Basal-like and 25 Classical genes excluding only `LEMD1`. Public labels were used only for descriptive group comparisons and explicitly labelled reference-anchored centroid definitions. Local R packages `clinfun` and `singscore` were unavailable, so the reproducible Python executor implemented the locked rank-scoring, permutation JT trend, bootstrap CI, and Hartigan dip-test workflows without changing counts, seeds, thresholds, or decision rules.
*   **Files / Analyses Affected:** `04_analysis/07_continuous_subtype_axis/PHASE5B_CONTINUOUS_AXIS_RESULTS.md`, `05_results/tables/phase5b_*`, `05_results/figures/phase5b_*`, `06_scripts/R/07_phase5b_continuous_axis.R`, `06_scripts/python/07_summarize_phase5b_axis.py`, `06_scripts/python/07_validate_phase5b_axis.py`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `00_admin/SKILL_USAGE_LOG.tsv`.

### D-16: Lock Phase 6A Processed Microbiome Abundance Source and Audit Scope
*   **Date:** 2026-07-01
*   **Decision:** Use verified public Supplementary Data 1 (`42003_2021_2557_MOESM4_ESM.xlsx`), sheet `Genus-level`, as the Phase 6A processed PRJNA719915 tumor microbiome abundance matrix. Preserve the released numerical scale and source matrix sample labels, and create a separate project crosswalk mapping the 62 matrix columns to verified patients, tumor numbers, microbiome BioSamples, and SRA runs. Limit Phase 6A to extraction, validation, descriptive QC, technical-metadata-only checks, and Phase 6B preprocessing recommendations.
*   **Alternatives Considered:** Reprocess raw SRA/FASTQ reads; use species-level Bracken estimates as the primary matrix; use relative-abundance source-data panels from Figure 2; remove potential contaminants during Phase 6A; run subtype, continuous-score, differential-abundance, survival, host-correlation, pathway, or target-prioritization analyses.
*   **Scientific and Operational Justification:** Supplementary Data 1 is the verified public abundance table containing processed class/order/family/genus/species profiles for all 62 tumors. The peer-review supplement states that Kraken2 plus Bracken was used, and the authors indicated genus-level profiles were more reliable than species-level estimates because Bracken species reassignment can be rough when reads remain at higher taxonomic nodes. The project has no sequenced negative-control runs, so decontam prevalence/frequency analysis cannot be performed and contamination risk must be distinguished from confirmed contamination. Preserving scale and delaying all microbiome association testing prevents unplanned preprocessing or outcome-selection bias.
*   **Files / Analyses Affected:** `03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz`, `01_metadata/microbiome_sample_crosswalk.tsv`, `05_results/tables/phase6a_*`, `05_results/figures/phase6a_*`, `04_analysis/04_microbiome_qc/PHASE6A_MICROBIOME_DATA_AUDIT.md`, `06_scripts/python/08_phase6a_microbiome_audit.py`, `06_scripts/python/08_validate_phase6a_microbiome.py`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `00_admin/SKILL_USAGE_LOG.tsv`.

---

### D-17: Lock Phase 6B Tumor Microbiome Preprocessing, Compositional Transformation, and Contamination-Sensitivity Framework
*   **Date:** 2026-07-01
*   **Decision:** Lock the microbiome preprocessing and compositional transformations. Establish the **Primary Analysis** using a 20% prevalence threshold (genera detected in >= 13 samples at abundance > 0.0), a cohort-specific pseudocount of `0.889651` (half the minimum observed non-zero value, specific to this matrix and non-transferable), Centered Log-Ratio (CLR) transformation, Aitchison distance, and retaining all 62 samples. Lock 8 specific sensitivity analysis runs covering alternative prevalence thresholds (10%, 30%), detection thresholds (>10.0 counts), alternative pseudocounts (`1.0`, `0.1`), Robust CLR (rCLR), Presence/Absence binarization with Jaccard distance, exclusion of three technical extreme samples (`Basal-like1`, `Hybrid18`, `Hybrid23`) flagged by richness and the matrix total-abundance proxy, and exclusion of High/Moderate Risk potential contaminant genera. Predefine downstream Phase 7 method compatibility (OLS, Spearman, MaAsLin2 with transform/normalization locked to NONE, PERMANOVA and PERMDISP dispersion controls, Presence/Absence) and multiple-testing corrections. Set overall decision status to `READY_WITH_CONTAMINATION_LIMITATIONS`.
*   **Alternatives Considered:** Rely on standard Euclidean distance on untransformed relative abundances; perform automatic deletion of all environmental/flagged genera; exclude extreme samples in the primary run; utilize ANCOM-BC2 or ALDEx2 as primary methods on pre-normalized data.
*   **Scientific and Operational Justification:** The audited matrix contains pre-normalized, non-integer counts, which violates the raw count assumptions of ANCOM-BC2 and ALDEx2. Applying OLS and Spearman models on CLR-transformed data is mathematically valid. In the absence of negative controls, automatically removing environmental genera might discard real biological signals, while ignoring contamination risks can lead to false positives; thus, establishing a strict Tier-based risk classification and comparing results with and without flagged contaminants (sensitivity analysis) is the most rigorous scientific path. Adding half the minimum observed non-zero value as a pseudocount avoids distorting ratios between real abundances. The matrix total-abundance proxy serves as a technical sensitivity covariate rather than a biological measurement. MaAsLin2 must not double-normalize or log-transform CLR inputs, and PERMANOVA must be accompanied by PERMDISP to ensure dispersion does not confound location tests.
*   **Files / Analyses Affected:** `04_analysis/04_microbiome_qc/PHASE6B_MICROBIOME_METHOD_LOCK.md`, `09_docs/methods/PDAC_microbiome_preprocessing_protocol.md`, `01_metadata/microbiome_preprocessing_parameter_inventory.tsv`, `05_results/tables/phase6b_filtering_candidate_summary.tsv`, `00_admin/PROJECT_STATUS.md`, `09_docs/planning/DECISION_LOG.md`, and `00_admin/SKILL_USAGE_LOG.tsv`.

---

### D-18: Execute Locked Phase 6C Analysis-Ready Microbiome Preprocessing
*   **Date:** 2026-07-01
*   **Decision:** Execute the Phase 6B-locked tumor microbiome preprocessing protocol without changing parameters after validation. Apply the primary rule of abundance > 0 detected in at least 20% of the 62 samples, retain all 62 samples and 122 genera, use the fixed source-specific pseudocount `0.889651`, generate primary CLR and Aitchison distance matrices, and create outcome-blind sensitivity representations for prevalence thresholds, detection threshold >10, pseudocount choices, robust CLR, presence/absence Jaccard distance, contaminant-flag exclusions, and technical extreme-sample exclusion.
*   **Alternatives Considered:** Force retained genera to the expected count rather than validating the matrix; delete flagged genera from the primary matrix; exclude technical outliers from the primary analysis; inspect subtype, host-expression, survival, or downstream association results while selecting preprocessing parameters.
*   **Scientific and Operational Justification:** The actual audited matrix reproduced the Phase 6B primary retained feature count of 122 genera, so execution could proceed without post hoc adjustment. Primary retention of flagged genera avoids treating potential contaminant-risk flags as confirmed contamination labels in the absence of sequenced negative controls. Sensitivity matrices provide robustness checks without using biological outcomes or changing the locked primary representation.
*   **Files / Analyses Affected:** `03_processed/microbiome/PRJNA719915_genus_primary_filtered.tsv.gz`, `03_processed/microbiome/PRJNA719915_genus_primary_CLR.tsv.gz`, `03_processed/microbiome/PRJNA719915_primary_aitchison_distance.tsv.gz`, `03_processed/microbiome/sensitivity/`, `05_results/tables/phase6c_*`, `05_results/figures/phase6c_*`, `04_analysis/04_microbiome_qc/PHASE6C_ANALYSIS_READY_MICROBIOME.md`, `06_scripts/R/09_phase6c_prepare_microbiome.R`, `06_scripts/python/09_phase6c_prepare_microbiome.py`, `06_scripts/python/09_summarize_phase6c_microbiome.py`, `06_scripts/python/09_validate_phase6c_microbiome.py`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `09_docs/planning/DECISION_LOG.md`, and `00_admin/SKILL_USAGE_LOG.tsv`.

---

### D-19: Prespecify and Lock the Microbiome-Host State Association Testing Framework
*   **Date:** 2026-07-01
*   **Decision:** Prespecify and lock the host outcome hierarchy, global community PERMANOVA parameters, primary genus-level OLS regression models (122 tests, CLR, HC3 robust standard errors, BH FDR q < 0.05), supporting methods (Spearman correlation, permutation-based association, bootstrap confidence intervals, MaAsLin2 with normalization=NONE, transform=NONE), covariate hierarchy (Model 0 primary, Model 1 technical sensitivity, Model 2 clinical sensitivity), presence/absence logistic regression orientation and sample size rules, contamination sensitivity framework, preprocessing sensitivities, evidence grading categories, and sample-level influence diagnostics.
*   **Alternatives Considered:** Optimize outcomes/thresholds according to microbiome findings; collapse clinical categories post hoc; automatically remove extreme/influential samples or flagged contaminants in primary runs; use standard Pearson correlation or unrobust standard errors.
*   **Scientific and Operational Justification:** Prespecifying and locking the statistical parameters before execution avoids post hoc data-dredging, selective outcome reporting, or significance-hunting. HC3 robust standard errors prevent inflated type I error rates due to heteroscedasticity. Evaluating clinical sensitivity feasibility revealed that Model 2 is not permitted due to 100% missing data in the clinical series. Proceeding with covariate limitations (primary Model 0 and technical Model 1) ensures the analysis remains statistically valid and transparent, while keeping tumor purity and other microenvironmental parameters marked as deferred covariates until they are validated.
*   **Files / Analyses Affected:** `04_analysis/08_host_microbiome_integration/PHASE7A_MICROBIOME_ASSOCIATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_microbiome_continuous_state_association_protocol.md`, `01_metadata/microbiome_association_parameter_inventory.tsv`, `05_results/tables/phase7a_model_matrix_feasibility.tsv`, `00_admin/PROJECT_STATUS.md`, and `00_admin/SKILL_USAGE_LOG.tsv`.

---

### D-20: Validate ESTIMATE-Derived Host TME Covariates and Add Feasible Sensitivity Models
*   **Date:** 2026-07-01
*   **Decision:** Calculate ESTIMATE-derived stromal score, immune score, ESTIMATE score, and inferred tumor purity for all 62 expression-mapped PDAC patients, validate score completeness and mapping integrity, and add only feasible host TME sensitivity models to the locked Phase 7A association inventory: Model 3P (host score + inferred tumor purity), Model 3I (host score + immune score), and Model 3S (host score + stromal score). Keep Model 0 as the primary association model and block a combined purity + immune + stromal + ESTIMATE model.
*   **Alternatives Considered:** Continue without host TME sensitivity covariates; place all ESTIMATE-derived covariates in a single multivariable model; treat inferred tumor purity as an independent pathology measurement; change the primary association model to include purity adjustment.
*   **Scientific and Operational Justification:** Phase 7A.5 generated complete ESTIMATE-derived covariates for all 62 patients using the official MD Anderson/R-Forge `estimate` package version 1.0.13. The separate sensitivity models met completeness, effective degrees-of-freedom, VIF, and condition-number criteria. The combined all-TME model failed due to severe mathematical collinearity (ESTIMATE score and inferred purity Spearman rho = -1.00; maximum VIF = Inf; condition number = 1.04e15). Because ESTIMATE-derived covariates and Moffitt scores come from the same host transcriptomic matrix, these covariates are robustness checks rather than independent measurements, and adjustment may remove genuine PDAC transcriptional-state biology.
*   **Files / Analyses Affected:** `01_metadata/host_tme_covariates.tsv`, `05_results/tables/phase7a5_host_covariate_qc.tsv`, `05_results/tables/phase7a5_host_covariate_correlations.tsv`, `05_results/tables/phase7a5_covariate_model_feasibility.tsv`, `05_results/figures/phase7a5_host_covariate_distributions.pdf`, `05_results/figures/phase7a5_host_covariate_correlation.pdf`, `04_analysis/08_host_microbiome_integration/PHASE7A5_HOST_COVARIATES.md`, `04_analysis/08_host_microbiome_integration/PHASE7A_MICROBIOME_ASSOCIATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_microbiome_continuous_state_association_protocol.md`, `01_metadata/microbiome_association_parameter_inventory.tsv`, `06_scripts/R/10_phase7a5_host_covariates.R`, `06_scripts/python/10_validate_phase7a5_host_covariates.py`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `00_admin/SKILL_USAGE_LOG.tsv`.


---

### D-21: Execute Locked Phase 7B Microbiome Association Analyses
*   **Date:** 2026-07-01
*   **Decision:** Execute the Phase 7A/7A.5 locked continuous tumor microbiome association models without changing outcomes, filters, transformations, covariates, FDR families, evidence rules, or sensitivity thresholds after inspecting results.
*   **Alternatives Considered:** Run clinical Model 2 despite missing clinical metadata; combine TME covariates in one model; optimize outcomes or thresholds after results; remove influential samples or flagged genera from the primary analysis.
*   **Scientific and Operational Justification:** The locked framework protects the primary continuous Moffitt50 association analysis from post hoc optimization and preserves null, negative, method-sensitive, and contamination-sensitive findings.
*   **Files / Analyses Affected:** `05_results/tables/phase7b_*`, `05_results/figures/phase7b_*`, `06_scripts/python/11_summarize_phase7b_associations.py`, `06_scripts/python/11_validate_phase7b_associations.py`, `06_scripts/R/11_phase7b_microbiome_associations.R`, and `04_analysis/08_host_microbiome_integration/PHASE7B_MICROBIOME_ASSOCIATION_RESULTS.md`.

---

| D-45 | 2026-06-30 | Phase 9B3C Independent Review | Required a formally isolated validation of Phase 9B3B logic. The executor code was sealed; validators and models were externally audited against the locked method. The audit failed Phase 9B3B due to code/statistical discrepancies. | Strict verification prevents false positives. |
| D-46 | 2026-06-30 | Phase 9B3R0 Pre-Reanalysis Audit | Mandated a dual-model (Implementation/Statistical) audit of Phase 9B3B before modifying any code. Reanalysis requires repairing fake controls, coverage bypasses, and convergence bugs. | Pinpoints root causes of failure objectively. |
| D-47 | 2026-07-04 | Phase 9B3R Reanalysis Ready | Project transitions to Phase 9B3R for code repair and reanalysis. No new data or post-hoc methodological changes allowed (e.g. statsmodels Z-test remains locked). | Enforces strict reproduction against locked plan. |

### D-22: Complete Phase 7C Independent Review of Phase 7B Microbiome–Host Integration Results
*   **Date:** 2026-07-02
*   **Decision:** Accept the Phase 7B tumor microbiome–host association analysis as statistically valid, reproducible, and robust. Report a final review decision of PASS. Accept the local unavailability of the supporting method MaAsLin2 under Option A.
*   **Alternatives Considered:** Require local installation and re-running of MaAsLin2 (Option B); modify locked evidence-grading criteria or thresholds post hoc; remove influential samples or contaminant-flagged taxa from the primary analysis.
*   **Scientific and Operational Justification:** The primary OLS models (HC3 robust standard errors) and global PERMANOVA community results are verified and exactly reproducible. Cross-validation with Spearman, permutation tests, and bootstrap confidence intervals provides sufficient statistical support, making MaAsLin2 completion unnecessary. Maintaining the locked evidence rules and documenting technical caveats (such as robust CLR sign reversals and environmental risk flags) preserves study transparency and prevents significance-hunting.
*   **Files / Analyses Affected:** `04_analysis/08_host_microbiome_integration/PHASE7C_INDEPENDENT_REVIEW.md`, `05_results/tables/phase7c_primary_candidate_audit.tsv`, `05_results/tables/phase7c_evidence_category_verification.tsv`, `05_results/tables/phase7c_review_findings.tsv`.

---

### D-23: Prespecify and Lock the Host-Mechanism Analysis Framework for Verified Microbiome Associations
*   **Date:** 2026-07-02
*   **Decision:** Prespecify and lock the host-mechanism analysis parameters, pathway activity collections (Hallmark, PROGENy as primary, Reactome, KEGG as secondary), DoRothEA TF regulons, WGCNA parameters (variance filter, scale-free topology fit selection, Dynamic dynamic tree cut module identification), and OLS regressions (HC3 robust standard errors, BH FDR q < 0.05). Restrict candidate taxa to verified robust (9) and suggestive (2) genera. Set the overall readiness decision to READY_WITH_TRANSFORMATION_LIMITATIONS.
*   **Alternatives Considered:** Run target prioritization or survival models immediately; use arbitrary gene-set collections; use data-driven threshold optimization or post hoc network pruning; combine TME score covariates in a single model.
*   **Scientific and Operational Justification:** Locking the prospective biological and statistical mechanisms before execution prevents post hoc selection bias, selective reporting of pathways, and multiple testing inflation. Highlighting the rCLR transformation sensitivity for 8 of the 9 robust genera ensures transparency about direction stability under composition transformations. Feasibility analysis confirmed decoupleR/WGCNA packages are not currently installed in the base conda environment, so models for Layers 1, 2, and 4 will be marked as locked pending environment configuration.
*   **Files / Analyses Affected:** `04_analysis/08_host_microbiome_integration/PHASE8A_HOST_MECHANISM_METHOD_LOCK.md`, `09_docs/methods/PDAC_host_microbiome_mechanism_protocol.md`, `01_metadata/host_mechanism_parameter_inventory.tsv`, `05_results/tables/phase8a_host_feature_feasibility.tsv`, `00_admin/PROJECT_STATUS.md`, `09_docs/planning/DECISION_LOG.md`, and `00_admin/SKILL_USAGE_LOG.tsv`.

---

### D-24: Prepare and Validate Phase 8B R Environment Before Host-Mechanism Execution
*   **Date:** 2026-07-02
*   **Decision:** Create a project-specific `renv` environment for the locked Phase 8B host-mechanism analysis and validate the required package capabilities using only synthetic matrices and package metadata. Mark Phase 8B as environment-ready after validating MSigDB Hallmark retrieval, Hallmark pathway scoring, PROGENy/decoupleR activity calculation, DoRothEA A/B/C regulon loading, TF activity calculation, WGCNA network construction, limma modeling, HC3 robust standard errors, and R object save/reload.
*   **Alternatives Considered:** Use the base system R library; proceed directly to Phase 8B with missing packages; install only primary pathway packages and defer WGCNA/secondary packages.
*   **Scientific and Operational Justification:** A project-local `renv` environment prevents further drift in the Phase 8B computational stack and avoids modifying the system R library. Synthetic-only capability tests confirm that the required algorithms can be invoked without inspecting or testing actual PDAC host-microbiome associations. WGCNA remains computationally heavier than the other host-feature layers, so full Phase 8B execution should use blockwise WGCNA if dense TOM memory is limiting while preserving the locked top-25%-MAD question.
*   **Files / Analyses Affected:** `renv.lock`, `renv/`, `.Rprofile`, `07_envs/phase8_r_environment.yml`, `07_envs/phase8_R_sessionInfo.txt`, `07_envs/phase8_package_validation.tsv`, `07_envs/phase8_capability_summary.tsv`, `07_envs/phase8_install_log.tsv`, `07_envs/phase8_pre_setup_R_environment.txt`, `07_envs/phase8_pre_setup_installed_packages.tsv`, `06_scripts/R/12_phase8a5_environment_setup.R`, `06_scripts/R/12_validate_phase8a5_environment.R`, `04_analysis/08_host_microbiome_integration/PHASE8A5_ENVIRONMENT_VALIDATION.md`, `05_results/tables/phase8a_host_feature_feasibility.tsv`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

### D-25: Execute Locked Phase 8B Host-Mechanism Analyses
*   **Date:** 2026-07-02
*   **Decision:** Execute the Phase 8A locked host-mechanism analysis without changing feature collections, taxa, covariates, WGCNA parameters, sensitivity rules, FDR families, or evidence categories after inspecting results.
*   **Alternatives Considered:** Promote suggestive taxa into the primary family; combine purity, immune, and stromal scores in one model; choose pathways or TFs after observing associations; omit rCLR or Moffitt50 safeguards.
*   **Scientific and Operational Justification:** The locked workflow preserves the prospective analysis plan and records transformation, composition, sample, and circularity limitations directly in evidence categories and sensitivity tables.
*   **Files / Analyses Affected:** `06_scripts/R/13_phase8b_host_mechanisms.R`, `06_scripts/python/13_summarize_phase8b_mechanisms.py`, `06_scripts/python/13_validate_phase8b_mechanisms.py`, `04_analysis/08_host_microbiome_integration/PHASE8B_HOST_MECHANISM_RESULTS.md`, `05_results/tables/phase8b_*`, and `05_results/figures/phase8b_*`.

---

### D-26: Complete Phase 8C Independent Review of Phase 8B Host–Microbiome Mechanism Integration Results
*   **Date:** 2026-07-02
*   **Decision:** Accept the Phase 8B host-microbiome mechanism analysis as statistically valid, reproducible, and robust. Report a final review decision of PASS. Accept the minor environmental (RENV_CONFIG_SANDBOX_ENABLED=FALSE) and cosmetic labeling (Model_i) findings as minor, and confirm that all 43 robust evidence rows meet the prospective lock rules.
*   **Alternatives Considered:** Reject Phase 8B results and require re-analysis due to minor labeling or environmental differences; modify locked evidence category thresholds post hoc.
*   **Scientific and Operational Justification:** The implementation of all host feature layers (Hallmark, PROGENy, TF regulons, WGCNA, limma, ranked enrichment) matches the locked design exactly. The 43 robust records (associated with *Ochrobactrum*) are computationally verified. Setting RENV_CONFIG_SANDBOX_ENABLED=FALSE was necessary due to workspace restrictions and did not introduce system package contamination. Disambiguating covariate models using the covariate column is sufficient, and the lack of joint regressions is correct since only one taxon is robust. Proceeding to Phase 9 external-validation planning is scientifically justified.
*   **Files / Analyses Affected:** `04_analysis/08_host_microbiome_integration/PHASE8C_INDEPENDENT_HOST_MECHANISM_REVIEW.md`, `05_results/tables/phase8c_robust_mechanism_audit.tsv`, `05_results/tables/phase8c_evidence_category_verification.tsv`, `05_results/tables/phase8c_wgcna_implementation_audit.tsv`, `05_results/tables/phase8c_review_findings.tsv`

---

### D-27: Identify, Evaluate, and Lock External-Validation Framework (Phase 9A)
*   **Date:** 2026-07-02
*   **Decision:** Identify, evaluate, and lock the external-validation cohort selection, signature transfer and scaling policies, validation statistical endpoints, and replication evidence categories. Define overall validation readiness as READY_WITH_MICROBIOME_LIMITATIONS.
*   **Alternatives Considered:** Conduct immediate validation downloading FASTQ files; modify discovery threshold settings post hoc; retrain classifiers on external datasets.
*   **Scientific and Operational Justification:** Locking the validation framework prospectively protects the external-validation phase (Phase 9B) from code changes, selective reporting of datasets, or threshold optimization. Using processed matrices first for host validation reduces data footprint and accelerates validation. Pseudobulking at the patient level is mandatory for single-cell data to prevent pseudoreplication, which treats thousands of single cells as independent biological replicates. Restricting microbiome validation to tumor tissue datasets (excluding oral, stool, fluid, and TCGA exploratory profiles) guarantees biological relevance, while noting the lack of sequenced controls in PRJNA542615 protects the study from claiming unverified tumor residency.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_external_validation_protocol.md`, `01_metadata/external_validation_dataset_inventory.tsv`, `01_metadata/external_validation_parameter_inventory.tsv`, `05_results/tables/phase9a_external_dataset_shortlist.tsv`, `05_results/tables/phase9a_signature_external_coverage_feasibility.tsv`, `05_results/tables/phase9a_external_analysis_resource_estimate.tsv`, `09_docs/references/phase9_external_validation_sources.bib`, and `09_docs/references/phase9_external_validation_source_audit.tsv`

---

### D-28: Execute Locked Phase 9B1 Independent Bulk-Transcriptome Validation
*   **Date:** 2026-07-03
*   **Decision:** Execute Phase 9B1 only for the three locked PRIORITY_1 bulk-host cohorts: TCGA_PAAD, GSE71729, and GSE62452. Use processed expression matrices and metadata only, preserve locked discovery signatures and replication rules, and prohibit single-cell, spatial, microbiome, target-prioritization, causal-mediation, and post hoc signature-modification analyses in this phase.
*   **Alternatives Considered:** Download raw GDC per-sample archives or raw SRA files; include non-bulk PRIORITY_1 datasets; merge cohorts before scoring; relax sample-count filters; upgrade TF proxy findings using literature plausibility.
*   **Scientific and Operational Justification:** Processed-matrix validation satisfies the Phase 9A bulk-host execution plan while minimizing data footprint and avoiding raw sequencing downloads. Tumor-only filtering reproduced the locked sample counts (TCGA_PAAD 178, GSE71729 145, GSE62452 69). Cross-cohort synthesis was limited to features with three comparable cohort-specific effects. Full decoupleR/VIPER external TF activity scoring was not executable in the managed environment, so TF proxy outputs are retained as `TO_VERIFY` and are not upgraded to externally replicated evidence.
*   **Files / Analyses Affected:** `02_data/external/phase9_bulk/`, `03_processed/external/phase9_bulk/`, `05_results/tables/phase9b1_*`, `05_results/figures/phase9b1_*.pdf`, `06_scripts/R/14_phase9b1_bulk_validation.R`, `06_scripts/python/14_prepare_phase9b1_bulk_data.py`, `06_scripts/python/14_validate_phase9b1_bulk_validation.py`, `04_analysis/09_external_validation/PHASE9B1_BULK_EXTERNAL_VALIDATION_RESULTS.md`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `00_admin/SKILL_USAGE_LOG.tsv`.

---

### D-29: Reject Phase 9B1 Bulk Validation Outputs due to Critical and Major Implementation Errors (Phase 9B1C)
*   **Date:** 2026-07-03
*   **Decision:** Reject the Phase 9B1 bulk-transcriptome external validation results and issue a final review decision of `FAIL_REQUIRES_REANALYSIS`. Enforce correct signature coverage thresholds ($\ge 80\%$) and proxy constraints to downgrade WGCNA module replication evidence from `EXTERNALLY_REPLICATED` to `PARTIALLY_REPLICATED_HOST_FEATURE` (replicated in TCGA_PAAD only), and reclassify pathways and transcription factors as `TO_VERIFY` until correct calculations can be completed.
*   **Alternatives Considered:** Accept the results with major corrections (`PASS_WITH_MAJOR_CORRECTIONS`) and proceed directly to single-cell validation; run local re-analysis during the review phase.
*   **Scientific and Operational Justification:** Issuing a `FAIL_REQUIRES_REANALYSIS` decision is the most scientifically sound action because the current implementation contains critical mathematical errors (omitted intercept in PurIST) and protocol violations (scoring signatures on 20% coverage instead of $\ge 80\%$) that invalidate the results. Fixing these errors requires modifying the pipeline scripts and re-running the validation, which is barred from being executed during the review phase itself but must be performed in a dedicated re-analysis cycle before any downstream single-cell validation can proceed.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B1C_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c_host_feature_audit.tsv`, `05_results/tables/phase9b1c_module_replication_audit.tsv`, `05_results/tables/phase9b1c_review_findings.tsv`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `09_docs/planning/DECISION_LOG.md`.

---

### D-30: Correct and Rerun Phase 9B1 Bulk Validation After Phase 9B1C Audit (Phase 9B1R)
*   **Date:** 2026-07-03
*   **Decision:** Execute a targeted Phase 9B1R reanalysis of only the three locked independent bulk cohorts, preserving Phase 9A datasets, discovery signatures, directions, coverage thresholds, FDR families, and evidence categories. Mark the original Phase 9B1 report as `SUPERSEDED_BY_PHASE9B1R`.
*   **Alternatives Considered:** Proceed directly to single-cell validation; accept Phase 9B1C reviewer classifications without rerunning corrected code; relax module coverage for microarray cohorts; retain Hallmark or TF proxy scores.
*   **Scientific and Operational Justification:** Rerunning the corrected implementation is required because the Phase 9B1C audit identified implementation errors that invalidated the original PurIST, Hallmark, TF, module, negative-control, synthesis, and evidence outputs. Phase 9B1R restores protocol compliance by including the PurIST intercept, enforcing the 80% module coverage threshold, using full-set Hallmark ssGSEA, executing DoRothEA/VIPER activity scoring where eligible, and running the locked negative controls.
*   **Files / Analyses Affected:** `06_scripts/R/14_phase9b1r_corrected_bulk_validation.R`, `06_scripts/python/14_prepare_phase9b1r_bulk_data.py`, `06_scripts/python/14_validate_phase9b1r_bulk_validation.py`, `04_analysis/09_external_validation/PHASE9B1R_CORRECTION_LOG.md`, `04_analysis/09_external_validation/PHASE9B1R_CORRECTED_BULK_EXTERNAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b1r_*`, `05_results/figures/phase9b1r_*`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `00_admin/SKILL_USAGE_LOG.tsv`.

---

### D-31: Complete Phase 9B1C2 Independent Review of Corrected Phase 9B1R Bulk Validation
*   **Date:** 2026-07-03
*   **Decision:** Accept the corrected Phase 9B1R bulk-transcriptome external validation results under the final review decision `PASS_WITH_MINOR_CORRECTIONS` and approve eventual proceeding to Phase 9B2 single-cell validation.
*   **Alternatives Considered:** Issue `FAIL_REQUIRES_REANALYSIS` if any critical errors remained; issue `PASS` without minor corrections and leave TFs hardcoded as TO_VERIFY.
*   **Scientific and Operational Justification:** Accepting the corrected results under `PASS_WITH_MINOR_CORRECTIONS` is scientifically justified because the calculations are now mathematically correct and all six findings from the previous audit have been successfully resolved. The minor correction involves the final evidence table's TF category reporting, which the executor now derives programmatically from the successfully executed VIPER activity statistics and matches the locked 12 externally replicated, 13 partially replicated, and 9 not replicated TF counts with 0 TO_VERIFY. Phase 9B2 remains pending until explicitly initiated.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c2_*`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `09_docs/planning/DECISION_LOG.md`.

---

### D-32: Final Closure of Phase 9B1C2 Minor Correction and Verification
*   **Date:** 2026-07-03
*   **Decision:** Close the Phase 9B1C2 review process with a final decision of **`PASS`**. Programmatic TF evidence classification has been successfully verified, the hardcoded `TO_VERIFY` assignment has been removed, and the executor-derived category counts are verified as: 12 `EXTERNALLY_REPLICATED_HOST_FEATURE`, 13 `PARTIALLY_REPLICATED_HOST_FEATURE`, 9 `NOT_REPLICATED`, and 0 `TO_VERIFY`.
*   **Alternatives Considered:** Defer programmatic reclassification to the R script; reject the minor correction re-analysis.
*   **Scientific and Operational Justification:** The minor correction resolves the remaining TF classification reporting issue (FIND_05 / FIND_07) by implementing programmatic, rule-based evidence classification in the executor pipeline. The validator script is updated to enforce this behavior. Result integrity is preserved, and no Hallmark, PurIST, WGCNA, or negative-control results were altered. Phase 9B2 is not performed in this pass.
*   **Files / Analyses Affected:** `06_scripts/R/14_phase9b1r_corrected_bulk_validation.R`, `06_scripts/python/14_validate_phase9b1r_bulk_validation.py`, `05_results/tables/phase9b1r_host_feature_replication_evidence.tsv`, `05_results/tables/phase9b1c2_correction_verification.tsv`, `05_results/tables/phase9b1c2_review_findings.tsv`, `04_analysis/09_external_validation/PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `09_docs/planning/DECISION_LOG.md`.

---

### D-33: Stop Phase 9B2 Before Data Acquisition Because Locked Single-Cell Records Disagree
*   **Date:** 2026-07-03
*   **Decision:** Stop Phase 9B2 before downloading or analyzing single-cell data. The full Phase 9A dataset inventory lists GSE111672, GSE154778, and GSE202051 as PRIORITY_1 Layer 2 single-cell resources; the Phase 9A shortlist lists only GSE111672 for single-cell; and the parameter inventory locks only GSE111672 for Layer 2. This violates the Phase 9B2 startup requirement to stop if accession, patient count, modality, or required files disagree across Phase 9A records.
*   **Alternatives Considered:** Proceed with GSE111672 only; include all three records from the full inventory; treat GSE202051 as spatial only; download processed files and resolve discrepancies afterward.
*   **Scientific and Operational Justification:** Choosing a subset after seeing conflicting records would silently alter the locked validation design and could produce non-reproducible cellular-source conclusions. Stopping before acquisition preserves the prospective lock, avoids post hoc cohort selection, and prevents running patient-aware models on an unstable cohort definition. No raw FASTQ/BAM files, spatial validation, microbiome validation, survival analysis, target prioritization, causal mediation, manuscript writing, or post hoc feature modification was performed.
*   **Files / Analyses Affected:** `01_metadata/phase9b2_single_cell_dataset_inventory.tsv`, `04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md`, `06_scripts/python/15_prepare_phase9b2_single_cell.py`, `06_scripts/python/15_validate_phase9b2_single_cell.py`, `06_scripts/R/15_phase9b2_single_cell_validation.R`, `00_admin/SKILL_USAGE_LOG.tsv`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

---

### D-34: Reconcile Layer 2 Single-Cell Cohort Set and Establish Sole Authoritative Dataset
*   **Date:** 2026-07-03
*   **Decision:** Reconcile the locked Layer 2 single-cell cohort set to include GSE111672 (representing Peng et al. 2019 / GSA CRA001160) as the sole PRIORITY_1 execution cohort. Downgrade GSE154778 to NOT_SUITABLE and set GSE202051's role in the inventory to Layer 3 (Spatial) only, excluding it from Layer 2.
*   **Alternatives Considered:** Keep all three datasets as Layer 2 PRIORITY_1. This was rejected because GSE154778 is dbGaP controlled-access (non-public download) and flow-sorted for CD45+ immune cells (lacking the malignant epithelial compartment required for subtyping validation). GSE202051 is a single-nucleus dataset with neoadjuvant treatment confounding and is already locked under Layer 3 spatial validation; including it in Layer 2 would result in redundant patient-level testing.
*   **Scientific and Operational Justification:** Restricting Layer 2 to GSE111672 (Peng et al. 2019 / GSA CRA001160) resolves the conflicting database records while maintaining experimental rigor. It excludes a dataset that lacks epithelial cells (GSE154778) and prevents duplicate patient-level statistical modeling of the spatial cohort (GSE202051). This alignment enables Phase 9B2 startup validation to proceed.
*   **Files / Analyses Affected:** `01_metadata/external_validation_dataset_inventory.tsv`, `05_results/tables/phase9a_external_dataset_shortlist.tsv`, `01_metadata/external_validation_parameter_inventory.tsv`, `04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_external_validation_protocol.md`, `04_analysis/09_external_validation/PHASE9A1_SINGLE_CELL_COHORT_RECONCILIATION.md`, `05_results/tables/phase9a1_single_cell_cohort_reconciliation.tsv`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `00_admin/SKILL_USAGE_LOG.tsv`, `09_docs/planning/DECISION_LOG.md`.

---

### D-35: Reconcile single-cell dataset accessions and provenance, expand Phase 9B2 execution cohort set (Phase 9A.2)
*   **Date:** 2026-07-03
*   **Decision:** Reconcile and correct single-cell accessions, publications, and patient counts across all metadata and planning files. Expand the Phase 9B2 execution cohort set to include PENG_CRA001160 (24 patients, primary Layer 2), LIN_GSE154778 (10 patients, secondary Layer 2), MONCADA_GSE111672 (6 patients, Layer 3 & secondary Layer 2), and HWANG_GSE202051 (43 patients, Layer 3 & secondary treatment-sensitivity Layer 2). Lock the new multi-cohort parameter inventory and shortlist. Update startup and validation scripts, and report final readiness as READY_WITH_DATASET_LIMITATIONS.
*   **Alternatives Considered:** Keep the single-cohort design of Phase 9A.1 (Peng only); ignore the errors in public metadata; skip correcting the scripts.
*   **Scientific and Operational Justification:** Reconciling the accessions and publications ensures absolute scientific accuracy and traceability. It avoids conflating Peng and Moncada, and corrects the false claim that GSE154778 was controlled-access and CD45+ sorted. Expanding the cohort set to include all four single-cell/spatial datasets maximizes statistical power and validation breadth while strictly respecting dataset-specific characteristics (e.g., neoadjuvant treatment in Hwang).
*   **Files / Analyses Affected:** `01_metadata/external_validation_dataset_inventory.tsv`, `05_results/tables/phase9a_external_dataset_shortlist.tsv`, `01_metadata/external_validation_parameter_inventory.tsv`, `05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv`, `04_analysis/09_external_validation/PHASE9A2_SINGLE_CELL_DATASET_PROVENANCE_CORRECTION.md`, `04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_external_validation_protocol.md`, `09_docs/references/phase9_external_validation_*`, `06_scripts/python/15_*`, `06_scripts/R/15_*`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `09_docs/planning/DECISION_LOG.md`.

---

### D-36: Stop Phase 9B2 Restart Because Authoritative Included Set Exceeds Primary PENG-Only Contract
*   **Date:** 2026-07-03
*   **Decision:** Stop the Phase 9B2 restart before data acquisition. The corrected Peng identity is valid (`PENG_CRA001160`, `CRA001160`, `PRJCA001063`, Peng et al. 2019, 24 untreated PDAC tumors and 11 control pancreases), but the authoritative Phase 9A.2 cohort table marks four cohorts as `included_in_phase9b2=yes`: `PENG_CRA001160`, `LIN_GSE154778`, `MONCADA_GSE111672`, and `HWANG_GSE202051`. The restart contract requires the primary execution set to be exactly `PENG_CRA001160`.
*   **Alternatives Considered:** Proceed with `PENG_CRA001160` despite the authoritative table; include all four Phase 9A.2 rows; silently treat `PRIORITY_2` or spatial-priority rows as out of scope without changing the authoritative inclusion flag.
*   **Scientific and Operational Justification:** Proceeding while the authoritative table permits additional cohorts would violate the explicit inclusion guardrail and create a post hoc cohort-selection ambiguity. Stopping preserves the patient-aware validation design and prevents reuse of any stopped-attempt biological results.
*   **Files / Analyses Affected:** `05_results/tables/phase9b2_restart_runtime_validation.tsv`, `04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

---

### D-37: Phase 9A.3 Execution-Scope Correction for Phase 9B2
*   **Date:** 2026-07-03
*   **Decision:** Perform an execution-scope correction for Phase 9B2. Separate dataset suitability from execution authorization in the authoritative cohort set by replacing `included_in_phase9b2` with explicit fields. Authorize `PENG_CRA001160` only for the current primary Phase 9B2 run, reserving `LIN_GSE154778`, `MONCADA_GSE111672`, and `HWANG_GSE202051` as planned supplementary cohorts requiring separate later authorization. Align the parameter inventory statuses accordingly, and update all startup/validation scripts.
*   **Alternatives Considered:** Modify the execution contract to run all four datasets in Phase 9B2-primary; scientific downgrade of the other datasets by deleting or removing them from planning files; proceeding without validator correction.
*   **Scientific and Operational Justification:** Restricting Phase 9B2-primary to `PENG_CRA001160` isolates the primary verification target to ensure rigorous and manageable data acquisition and downstream analysis. Keeping the other cohorts in the inventory as planned supplementary cohorts preserves their scientific validity without creating operational conflicts. Updating the validators to evaluate against `included_in_phase9b2_primary == TRUE` prevents false failure flags. Both previous stopped attempts are preserved in the audit trail.
*   **Files / Analyses Affected:** `05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv`, `01_metadata/external_validation_parameter_inventory.tsv`, `06_scripts/python/15_prepare_phase9b2_single_cell.py`, `06_scripts/python/15_validate_phase9b2_single_cell.py`, `06_scripts/python/15_validate_provenance_consistency.py`, `04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_external_validation_protocol.md`, `00_admin/PROJECT_STATUS.md`, `05_results/tables/phase9a3_phase9b2_execution_scope.tsv`, `04_analysis/09_external_validation/PHASE9A3_PHASE9B2_EXECUTION_SCOPE_CORRECTION.md`, `09_docs/planning/DECISION_LOG.md`.

---

### D-38: Execute Phase 9B2 Primary Single-Cell Cellular-Source Analysis on PENG_CRA001160 Only
*   **Date:** 2026-07-03
*   **Decision:** Execute Phase 9B2-primary using only `PENG_CRA001160` (`CRA001160`, `PRJCA001063`, Peng et al. 2019) after startup provenance validation confirmed that the authoritative included primary cohort set contains only this dataset. Use processed matrix and annotation files from the official CNCB GSA source; do not download raw FASTQ or BAM data.
*   **Alternatives Considered:** Include supplementary Layer 2 or Layer 3 cohorts during the primary run; reuse stopped-attempt biological outputs; use TF-symbol expression as a proxy for TF activity; reconstruct WGCNA modules in the single-cell cohort.
*   **Scientific and Operational Justification:** The primary single-cell question is cellular source and patient-level malignant-state heterogeneity for locked host features. Patient-aware pseudobulks preserve patient independence, while DoRothEA/VIPER regulon activity and locked transferred module membership preserve the Phase 8/9 method hierarchy. Supplementary and spatial cohorts require separate authorization.
*   **Outcome:** Analyzed 57,530 cells from 24 PDAC tumor patients and 11 control pancreas donors. Generated Phase 9B2 inventories, QC tables, annotation and malignant-cell audits, pseudobulks, host-state scores, Hallmark/module scores, TF activities, cellular-source models, malignant-axis associations, cell-composition sensitivity, tumor-control descriptive results, negative-control records, evidence classifications, figures, validators, and the Phase 9B2 results report. Negative-control permutation and expression-matched null items remain `TO_VERIFY`.
*   **Files / Analyses Affected:** `02_data/external/phase9_single_cell/PENG_CRA001160/`, `01_metadata/phase9b2_single_cell_dataset_inventory.tsv`, `05_results/tables/phase9b2_*`, `05_results/figures/phase9b2_*`, `05_results/models/phase9b2/`, `06_scripts/python/15_prepare_phase9b2_single_cell.py`, `06_scripts/python/15_validate_phase9b2_single_cell.py`, `06_scripts/R/15_phase9b2_single_cell_validation.R`, `04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

---

### D-39: Reject Phase 9B2 Single-Cell Cellular-Source Analysis due to Incomplete Negative Controls and Coverage Violations (Phase 9B2C)
*   **Date:** 2026-07-03
*   **Decision:** Reject the primary Phase 9B2 single-cell validation results and report the final decision as `FAIL_REQUIRES_REANALYSIS`. Block spatial-validation planning until corrections are made.
*   **Alternatives Considered:** Accept the results with minor or major corrections; classify the unexecuted negative controls as optional; ignore the WGCNA module coverage violations.
*   **Scientific and Operational Justification:** The prospective method lock (Phase 9A) mandates the execution of negative controls and strictly requires that any feature with $< 80\%$ coverage be excluded from formal inference. The primary Phase 9B2 execution failed both rules by generating hardcoded placeholder rows for negative controls and by scoring/categorizing all 5 WGCNA modules (all had coverage $< 49\%$). Proceeding with spatial validation under these major implementation errors would compromise scientific integrity.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B2C_SINGLE_CELL_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b2c_review_findings.tsv`, `05_results/tables/phase9b2c_feature_evidence_audit.tsv`, `05_results/tables/phase9b2c_negative_control_audit.tsv`.

---

### D-40: Complete Phase 9B2R Corrective Single-Cell Reanalysis and Supersede Initial Phase 9B2 Results
*   **Date:** 2026-07-03
*   **Decision:** Perform one consolidated corrective Phase 9B2R reanalysis after Phase 9B2C, preserving `PENG_CRA001160` as the sole active primary cohort and marking the initial Phase 9B2 report as `SUPERSEDED_BY_PHASE9B2R`.
*   **Alternatives Considered:** Patch only the evidence table; rerun modules despite low coverage; begin spatial validation before correcting Phase 9B2; add supplementary single-cell cohorts.
*   **Scientific and Operational Justification:** The Phase 9A lock requires coverage >= 80% for transferred signatures and mandatory negative controls. Phase 9B2R enforces this rule: MEblack, MEblue, MEgreen, MEtan, and MEgreenyellow are all below 80% coverage and therefore remain coverage-only descriptive records classified as `INSUFFICIENT_SINGLE_CELL_DATA`. Corrected negative controls were executed where technically applicable using locked seed 2026 and locked iteration counts; module random controls are technically inapplicable after all modules fail eligibility.
*   **Outcome:** Generated corrected `phase9b2r_*` tables and figures. Only `HALLMARK_PROTEIN_SECRETION` remains associated with malignant-cell Moffitt50 axis at q < 0.10. Corrected cellular-source categories are 20 `CELL_COMPOSITION_EXPLAINED`, 5 `INSUFFICIENT_SINGLE_CELL_DATA`, 4 `STROMAL_OR_IMMUNE_SOURCE_SUPPORTED`, 1 `MALIGNANT_CELL_INTRINSIC_SUPPORT`, 1 `PARTIAL_CELLULAR_SUPPORT`, and 1 `NOT_SUPPORTED_AT_CELLULAR_LEVEL`. Ochrobactrum was not tested, and spatial validation was not performed.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md`, `04_analysis/09_external_validation/PHASE9B2R_CORRECTION_LOG.md`, `04_analysis/09_external_validation/PHASE9B2R_CORRECTED_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md`, `05_results/tables/phase9b2r_*`, `05_results/figures/phase9b2r_*`, `06_scripts/R/15_phase9b2r_corrected_single_cell_validation.R`, `06_scripts/python/15_prepare_phase9b2r_single_cell.py`, `06_scripts/python/15_validate_phase9b2r_single_cell.py`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, and `01_metadata/file_manifest.tsv`.

---

### D-41: Issue Final PASS Decision for Phase 9B2C2 Independent Review of Corrected Single-Cell Validation
*   **Date:** 2026-07-03
*   **Decision:** Issue a final closure decision of PASS for Phase 9B2C2. Authorize proceeding to Phase 9B3 spatial-validation planning, as all critical, major, and moderate issues are resolved.
*   **Alternatives Considered:** Request further reanalysis; block spatial validation planning; retain warnings or caveats on completed negative controls.
*   **Scientific and Operational Justification:** Independent audit of the corrected Phase 9B2R single-cell reanalysis confirms that all Phase 9B2C findings are fully corrected. FIND_01 is closed because all 64 applicable negative controls have been executed (patient/cell-type permutations and unrelated pathways) and statistics computed; FIND_02 is closed because all 5 low-coverage WGCNA modules are excluded from formal models and classified as INSUFFICIENT_SINGLE_CELL_DATA; FIND_03 is closed because all 25 TF categories are programmatically derived. Provenance, cell counts (57,530), patient-aware pseudobulking, and composition models are fully verified. No critical or major findings remain.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B2C2_CORRECTED_SINGLE_CELL_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b2c2_*`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

---

### D-42: Lock Spatial Validation Planning and Prospective Method Lock (Phase 9B3A)
*   **Date:** 2026-07-03
*   **Decision:** Freeze the spatial feature hierarchy, qualified spatial datasets, authoritative spatial execution set, and prospective statistical-design and negative-control parameters for Phase 9B3.
*   **Alternatives Considered:** Authorize additional spatial datasets (e.g. GSE274103, GSE272362) for primary execution; use spots as independent biological replicates; skip spatial coordinate permutation controls.
*   **Scientific and Operational Justification:** Restricting execution authorization to HWANG_GSE202051 (treatment-naïve and neoadjuvant-treated subsets) and MONCADA_GSE111672 ensures scientific focus on established, high-quality, and independent cohorts. Blocking spot-level pseudoreplication prevents false positive associations. Seeded permutations and unrelated biological pathway negative controls ensure statistical validity.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B3A_SPATIAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_spatial_validation_protocol.md`, `01_metadata/phase9b3_spatial_*`, `05_results/tables/phase9b3a_*`, `06_scripts/python/16_validate_phase9b3a_spatial_plan.py`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, and `01_metadata/file_manifest.tsv`.

---

### D-43: Phase 9B3A.1 Spatial Validation Design Consistency Correction
*   **Date:** 2026-07-03
*   **Decision:** Perform limited design consistency correction for Phase 9B3 spatial validation to reconcile Hypothesis 3 lymphoid adjustment, define nested LMM random-effects structures, split counts, and use correct statistical model terminology (LMM).
*   **Alternatives Considered:** Maintain the model formula without lymphoid adjustment; ignore nesting structures and treat spots/ROIs as independent biological replicates; pool GeoMx and ST matrices directly.
*   **Scientific and Operational Justification:** Fully incorporating lymphoid adjustment ensures that our statistical models match our primary biological hypotheses. Restricting random effects to LMM class and defining nested structures (such as `(1 | patient_id) + (1 | patient_id:section_id)` for Moncada) blocks spot-level and section-level pseudoreplication, which would otherwise inflate significance and lead to false positive findings. Splitting patient and section counts eliminates ambiguity and matches official GEO records.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B3A1_SPATIAL_DESIGN_CONSISTENCY_CORRECTION.md`, `05_results/tables/phase9b3a1_spatial_analysis_unit_and_models.tsv`, `06_scripts/python/generate_phase9b3a1_tables.py`, `04_analysis/09_external_validation/PHASE9B3A_SPATIAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_spatial_validation_protocol.md`, `01_metadata/phase9b3_spatial_dataset_inventory.tsv`, `01_metadata/phase9b3_spatial_parameter_inventory.tsv`, `05_results/tables/phase9b3a_authoritative_spatial_cohort_set.tsv`, `05_results/tables/phase9b3a_spatial_dataset_qualification.tsv`, `06_scripts/python/16_validate_phase9b3a_spatial_plan.py`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, and `01_metadata/file_manifest.tsv`.

---

### D-44: Phase 9B3A.2 Spatial Validation Hierarchy Final Correction
*   **Date:** 2026-07-03
*   **Decision:** Perform limited spatial validation planning amendment (Phase 9B3A.2) to lock Hwang GeoMx ROI pairing models (Model A, B, and C), reclassify Moncada ST as exploratory cross-platform spatial consistency (due to low patient sample size n=2), and strictly enforce the cross-platform matrix pooling ban.
*   **Alternatives Considered:** Treat Hwang segments as independent without ROI random effect nested terms; use a two-patient meta-analysis for Moncada as formal external replication; allow pooling/merging of GeoMx and ST expression matrices.
*   **Scientific and Operational Justification:** GeoMx DSP datasets contain paired tumor/stroma segments in each ROI. Modeling them with nested random intercepts `(1 | patient_id) + (1 | patient_id:ROI_id)` preserves the paired design and avoids location confounding. Reclassifying Moncada ST to exploratory cross-platform consistency protects against over-interpreting low patient count data as population-level replication. Banning matrix pooling prevents direct cross-platform normalization issues between ROI-based DSP and grid-based ST.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B3A2_SPATIAL_HIERARCHY_FINAL_CORRECTION.md`, `05_results/tables/phase9b3a2_spatial_model_hierarchy.tsv`, `06_scripts/python/generate_phase9b3a2_tables.py`, `04_analysis/09_external_validation/PHASE9B3A_SPATIAL_VALIDATION_METHOD_LOCK.md`, `09_docs/methods/PDAC_spatial_validation_protocol.md`, `01_metadata/phase9b3_spatial_dataset_inventory.tsv`, `01_metadata/phase9b3_spatial_parameter_inventory.tsv`, `05_results/tables/phase9b3a_authoritative_spatial_cohort_set.tsv`, `05_results/tables/phase9b3a_spatial_dataset_qualification.tsv`, `05_results/tables/phase9b3a_spatial_resource_estimate.tsv`, `05_results/tables/phase9b3a1_spatial_analysis_unit_and_models.tsv`, `06_scripts/python/16_validate_phase9b3a_spatial_plan.py`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, and `01_metadata/file_manifest.tsv`.

---

### D-45: Execute Phase 9B3B Prospective Spatial-Transcriptomic Validation
*   **Date:** 2026-07-03
*   **Decision:** Execute the prospectively locked Phase 9B3B spatial-transcriptomic validation and advance to independent review with readiness status `READY_FOR_PHASE9B3C_INDEPENDENT_REVIEW`.
*   **Alternatives Considered:** Stop because official processed usable counts differed from planning estimates; pool Hwang GeoMx and Moncada ST matrices to increase apparent sample size; treat Moncada sections or spots as independent biological replication.
*   **Scientific and Operational Justification:** The locked Phase 9B3A/A.1/A.2 validators passed before acquisition. Official processed Hwang and Moncada files were acquired without FASTQ/BAM data. The usable official processed counts differed from planning estimates (Hwang naive 13 patients/127 ROIs/373 segments; Hwang treated 7 patients/67 ROIs/197 segments; Moncada 2 patients/6 sections/3,119 parsed spots), but the analysis preserved patient-aware inference, ROI pairing, cohort separation, reduced-model rules, coverage thresholds, negative controls, and the cross-platform matrix-pooling ban. Moncada was retained as exploratory consistency only.
*   **Outcome:** Hwang naive protein secretion showed strong tumor-compartment localization but no significant tumor-only or paired-contrast Moffitt50 axis association, yielding `PARTIAL_SPATIAL_SUPPORT`. Treated sensitivity and Moncada exploratory results were documented separately. Low-coverage WGCNA modules and TF regulon features were excluded from formal inference when not eligible. Required computational validators passed.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B3B_SPATIAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b3b_*`, `05_results/figures/phase9b3b_*`, `06_scripts/R/16_phase9b3b_spatial_validation.R`, `06_scripts/python/16_prepare_phase9b3b_spatial.py`, `06_scripts/python/16_phase9b3b_spatial_validation.py`, `06_scripts/python/16_validate_phase9b3b_spatial.py`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

---

### D-46: Reject Phase 9B3B Spatial Validation due to Hardcoded Negative Controls and Model Violations (Phase 9B3C)
*   **Date:** 2026-07-04
*   **Decision:** Reject the Phase 9B3B spatial-transcriptomic validation results and report the final decision as `FAIL_REQUIRES_REANALYSIS`.
*   **Alternatives Considered:** Accept the results with major corrections; classify the unexecuted negative controls as informational; ignore the non-convergence of the treated cohort Model C.
*   **Scientific and Operational Justification:** The independent audit of Phase 9B3B spatial validation identified critical execution errors: (1) negative controls (permutations, random gene sets, and leakage checks) were not actually computed but filled with hardcoded placeholders; (2) the ineligible feature `HALLMARK_SPERMATOGENESIS` was included in models despite low coverage (37.0% < 80% threshold); and (3) a non-converged treated cohort Model C was incorrectly used to contribute biological evidence. Proceeding to cross-layer synthesis with these flaws is barred to preserve scientific integrity. A complete reanalysis is required.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B3C_SPATIAL_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b3c_*`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

---

### D-47: Complete Phase 9B3R Corrected Spatial Validation Reanalysis
*   **Date:** 2026-07-04
*   **Decision:** Complete the corrected Phase 9B3R reanalysis of the failed Phase 9B3B spatial validation and set readiness to `READY_FOR_PHASE9B3C2_COMPLETE_INDEPENDENT_REVIEW`.
*   **Alternatives Considered:** Begin cross-layer synthesis before correction; alter the locked cohorts or model hierarchy; switch post hoc from the locked statsmodels asymptotic Z inference to an unapproved small-sample engine.
*   **Scientific and Operational Justification:** The Phase 9B3R0 repair specification authorized correction of negative-control execution, feature eligibility enforcement, nonconvergence handling, explicit inference-method documentation, and programmatic evidence derivation. The reanalysis preserves the locked Hwang and Moncada cohorts, patient-level replicate rules, ROI pairing, feature hierarchy, 80% coverage threshold, and evidence rules. Real null distributions were generated for coordinate permutations, size-matched and expression-matched random gene sets, unrelated Hallmark controls, label permutations, and leakage controls. Ineligible features are classified as insufficient spatial data. Nonconverged treated Model C is retained only as an NA audit row and excluded from q values, figures, and evidence.
*   **Outcome:** Hwang naive Model A remains a converged tumor-compartment enrichment result, while naive Models B/C remain null for the continuous Moffitt50 spatial-axis hypothesis. Hwang treated Model A remains a converged sensitivity result, Model B remains null, and treated Model C is nonconverged with inferential fields set to NA. Moncada remains exploratory with 1 of 6 positive sections. The final corrected evidence category is `PARTIAL_SPATIAL_SUPPORT`.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B3R_CORRECTED_SPATIAL_VALIDATION_RESULTS.md`, `04_analysis/09_external_validation/PHASE9B3R_CORRECTION_LOG.md`, `04_analysis/09_external_validation/PHASE9B3B_SPATIAL_VALIDATION_RESULTS.md`, `05_results/tables/phase9b3r_*`, `05_results/figures/phase9b3r_hwang_primary_models.pdf`, `06_scripts/python/16_phase9b3r_spatial_validation.py`, `06_scripts/python/16_validate_phase9b3r_spatial.py`, `06_scripts/python/test_phase9b3r_spatial.py`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

### D-47: Phase 10A Cross-Layer Synthesis Method Lock
*   **Date:** 2026-07-04
*   **Decision:** Lock the 10-level cross-layer evidence hierarchy, objective synthesis rules, and target prioritization framework.
*   **Alternatives Considered:** Rely on qualitative literature summaries or promote targets that failed formal external replication.
*   **Scientific and Operational Justification:** Ensures unbiased, multi-layer validation evidence dictates target promotion rather than literature-driven cherry picking.
*   **Files / Analyses Affected:** `04_analysis/10_target_prioritization/PHASE10A_CROSS_LAYER_SYNTHESIS_METHOD_LOCK.md`, `09_docs/methods/PDAC_cross_layer_synthesis_protocol.md`, `05_results/tables/phase10a_*`

### D-48: Execute Phase 10B Cross-Layer Evidence Synthesis and Target Prioritization
*   **Date:** 2026-07-04
*   **Decision:** Apply the locked prioritization framework to multi-layer supported and partially replicated targets.
*   **Alternatives Considered:** Rely on pathway overrepresentation alone, or select targets based solely on literature ubiquity.
*   **Scientific and Operational Justification:** Evaluated candidates using orthogonal, independent data (OpenTargets for tractability, GTEx for tumor-vs-normal selectivity). Identifies CTCFL as the strongest biological candidate due to its Cancer-Testis Antigen profile, while rejecting BHLHE40 (broad normal expression, low tractability) and HALLMARK_PROTEIN_SECRETION (pan-essential pathway).
*   **Files / Analyses Affected:** `05_results/tables/phase10b_candidate_target_scores.tsv`, `04_analysis/10_target_prioritization/PHASE10B_TARGET_PRIORITIZATION_RESULTS.md`

### D-49: Reject Phase 10B Target Prioritization (Phase 10C)
*   **Date:** 2026-07-04
*   **Decision:** Reject the Phase 10B target prioritization results due to lack of reproducibility, hardcoded unverified claims, missing evaluation for HALLMARK_SPERMATOGENESIS, and failure to apply locked thresholds.
*   **Alternatives Considered:** Accept Phase 10B with minor corrections.
*   **Scientific and Operational Justification:** External database claims (GTEx, OpenTargets) were not saved or executed reproducibly within the project scripts. CTCFL was selected via post-hoc reasoning ignoring its single-cell composition-sensitive status which should penalize cell-type specificity. Objective thresholds were bypassed for qualitative descriptions.
*   **Files / Analyses Affected:** `04_analysis/10_target_prioritization/PHASE10C_INDEPENDENT_REVIEW.md`, `05_results/tables/phase10c_*`

### D-50: Complete Phase 10B-R Corrected Target Prioritization Reanalysis
*   **Date:** 2026-07-04
*   **Decision:** Re-execute Phase 10B target prioritization as Phase 10B-R and set readiness to `READY_FOR_PHASE10C2_INDEPENDENT_REVIEW`.
*   **Alternatives Considered:** Reuse failed Phase 10B target annotations; promote CTCFL/BORIS based on post-hoc cancer-testis antigen reasoning; proceed directly to manuscript drafting; omit partially replicated or composition-sensitive candidates.
*   **Scientific and Operational Justification:** The corrected reanalysis derives candidate inclusion, evidence category, framework points, and penalties from locked Phase 10A inputs. Every Phase 10A evidence-inventory candidate is included. `CELL_COMPOSITION_EXPLAINED` is penalized rather than treated as malignant-cell specificity, preventing CTCFL/BORIS and BHLHE40 rescue by narrative literature support. External database claims are reproducibility-gated through a local query audit; unavailable OpenTargets, GTEx, and ChEMBL gene-symbol rows are labelled `NOT_RUN_DATABASE_UNAVAILABLE` and do not contribute values.
*   **Outcome:** HALLMARK_PROTEIN_SECRETION remains the highest cross-layer evidence feature but is retained as a supported biological feature rather than promoted to a direct gene target. HALLMARK_SPERMATOGENESIS, BHLHE40, and CTCFL/BORIS are not prioritized because their single-cell support is composition-explained. Ochrobactrum remains discovery-only with no microbial localization, physical interaction, or causality validation. Ineligible WGCNA modules remain ineligible.
*   **Files / Analyses Affected:** `04_analysis/10_target_prioritization/PHASE10BR_CORRECTED_TARGET_PRIORITIZATION_RESULTS.md`, `05_results/tables/phase10br_*`, `06_scripts/python/18_phase10br_cross_layer_synthesis.py`, `06_scripts/python/18_validate_phase10br_synthesis.py`, `06_scripts/python/test_phase10br_synthesis.py`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.

---

### D-31: Complete Phase 9B1C2 Independent Review of Corrected Phase 9B1R Bulk Validation
*   **Date:** 2026-07-03
*   **Decision:** Accept the corrected Phase 9B1R bulk-transcriptome external validation results under the final review decision `PASS_WITH_MINOR_CORRECTIONS` and approve proceeding to Phase 9B2 single-cell validation.
*   **Alternatives Considered:** Issue `FAIL_REQUIRES_REANALYSIS` if any critical errors remained; issue `PASS` without minor corrections and leave TFs hardcoded as TO_VERIFY.
*   **Scientific and Operational Justification:** Accepting the corrected results under `PASS_WITH_MINOR_CORRECTIONS` is scientifically justified because the calculations are now mathematically correct and all six findings from the previous audit have been successfully resolved. The minor correction involves the final evidence table's TF category reporting, which the reviewer independently resolved by reclassifying the 34 TFs using the successfully executed VIPER activity statistics. Proceeding to Phase 9B2 is approved since no critical or major implementation errors remain.
*   **Files / Analyses Affected:** `04_analysis/09_external_validation/PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW.md`, `05_results/tables/phase9b1c2_*`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `09_docs/planning/DECISION_LOG.md`.


### D-51: Phase 10C2 Independent Review Final Decision

*   **Date**: 2026-07-04
*   **Phase**: Phase 10C2
*   **Description**: Phase 10B-R results passed the independent review. All errors from Phase 10C were successfully corrected.
*   **Rationale**: Target prioritization was fully reproducible, objective, and compliant with Phase 10A method locks.
*   **Impact**: Project can proceed to Phase 11: Manuscript Drafting.
| 2026-07-04 | **D-47** | Complete Phase 11A manuscript claim map and figure planning, enforcing strict anti-causal language and evidence constraints | `04_analysis/11_manuscript/*`, `05_results/tables/phase11a_*` |

### D-52: Complete Phase 11B Manuscript Draft Finalization

*   **Date:** 2026-07-04
*   **Decision:** Finalize the PASS-approved Phase 11B manuscript draft for independent review with readiness status `READY_FOR_PHASE11C_MANUSCRIPT_INDEPENDENT_REVIEW`.
*   **Alternatives Considered:** Rewrite manuscript body after manual approval; relax guardrails around CTCFL/BORIS, microbiome causality, spatial validation, or null findings; proceed without claim-to-text traceability.
*   **Scientific and Operational Justification:** The manual revision already contained the approved biological conclusions and required caveats. Phase 11B therefore preserved the manuscript body and added an explicit claim-to-text trace plus a manuscript validator. The validator fails on CTCFL/BORIS promotion or literature rescue, microbial causality/localization/physical interaction claims, full basal-classical spatial validation claims, omitted `PARTIAL_SPATIAL_SUPPORT`, omitted null replication findings, omitted rCLR/Herbaspirillum limitations, missing Moncada exploratory 1/6 wording, and untraceable major claims.
*   **Outcome:** Phase 11B validation passed. The manuscript is ready for Phase 11C independent review.
*   **Files / Analyses Affected:** `04_analysis/11_manuscript/PHASE11B_MANUSCRIPT_DRAFT.md`, `05_results/tables/phase11b_claim_to_text_trace.tsv`, `06_scripts/python/19_validate_phase11b_manuscript_draft.py`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `01_metadata/file_manifest.tsv`, and `09_docs/planning/DECISION_LOG.md`.
### D-47: Complete Phase 11C Independent Manuscript Review
*   **Date:** 2026-07-04
*   **Decision:** Accept the Phase 11B manuscript draft based on a successful independent claim-control review. No biological overstatements or unauthorized claims were identified.
*   **Alternatives Considered:** Require minor or major revisions before assembly.
*   **Scientific and Operational Justification:** The manuscript correctly classifies spatial validation as partial, explicitly excludes CTCFL from prioritization due to composition sensitivity, refrains from asserting causality or localization for Ochrobactrum, and transparently reports null and exploratory findings.
*   **Files / Analyses Affected:** , .

---

### D-53: Complete Phase 11E Language and Format Review
*   **Date:** 2026-07-04
*   **Decision:** Finalize the Phase 11E language and format review of the assembled manuscript. Set readiness status to `READY_FOR_PHASE11F_FINAL_CLAIM_AUDIT`.
*   **Alternatives Considered:** Maintain the Phase 11D draft without editing; modify biological conclusions or evidence category levels; proceed directly to claim audit without independent validation of the language-edited file.
*   **Scientific and Operational Justification:** The language and format review improved readability, transitioned Results 2.1 to active voice, resolved repetitive microbial disclaimers into a unified compound sentence in the Introduction, and standardized all figure and supplementary legends to 2-3 sentences. All claim-control boundaries (non-causal microbial co-occurrence, partial spatial support for protein secretion, and composition-based exclusion of CTCFL) were strictly preserved without introducing stronger claims or new biological results.
*   **Outcome:** Phase 11E validation successfully passed. The project is verified as ready for Phase 11F.
*   **Files / Analyses Affected:** `04_analysis/11_manuscript/PHASE11E_FULL_MANUSCRIPT_LANGUAGE_EDITED.md`, `04_analysis/11_manuscript/PHASE11E_LANGUAGE_FORMAT_REVIEW.md`, `05_results/tables/phase11e_language_edit_log.tsv`, `06_scripts/python/19_validate_phase11e_language_format.py`, `00_admin/PROJECT_STATUS.md`, `01_metadata/file_manifest.tsv`, `00_admin/SKILL_USAGE_LOG.tsv`, and `09_docs/planning/DECISION_LOG.md`.




---

### D-55: Complete Phase 11G Administrative Finalization
*   **Date:** 2026-07-04
*   **Decision:** Complete Phase 11G reference submission and administrative finalization. Register Phase 11G / 11G-R1 / 11G-R2 / 11G-R3 files and update the file manifest. Set final readiness decision to `READY_FOR_PHASE11H_SUBMISSION_PACKAGE_ASSEMBLY`.
*   **Alternatives Considered:** None.
*   **Scientific and Operational Justification:** Reference callout repair, verification, and repair were completed (Phase 11G-R1b, Phase 11G-R2, Phase 11G-R3) resulting in a final decision of PASS. The manuscript scientific content, biological conclusions, evidence categories, and target rankings were verified to remain unchanged. The package is now ready for submission assembly.
*   **Files / Analyses Affected:** `01_metadata/file_manifest.tsv`, `04_analysis/11_manuscript/PHASE11G*.md`, `05_results/tables/phase11g*.tsv`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `09_docs/planning/DECISION_LOG.md`.


---

### D-56: Complete Post-Phase 11G Workspace Cleanup Finalization
*   **Date:** 2026-07-04
*   **Decision:** Audited and deleted temporary agent utility scripts, and safely ignored intermediate data models in version control.
*   **Alternatives Considered:** None.
*   **Scientific and Operational Justification:** Intermediate data objects should not be tracked by git, and temporary agent scripts are no longer needed. The cleanup keeps the workspace clean for the final submission package assembly.
*   **Files / Analyses Affected:** `.gitignore`, `04_analysis/11_manuscript/POST_PHASE11G_WORKSPACE_CLEANUP_AUDIT.md`, `01_metadata/file_manifest.tsv`, `00_admin/PROJECT_STATUS.md`, `00_admin/SKILL_USAGE_LOG.tsv`, `09_docs/planning/DECISION_LOG.md`.

### D-57: Complete Phase 11H Submission Package Assembly
*   **Date:** 2026-07-04
*   **Decision:** Accept the Phase 11H assembled submission package and register outputs. Set readiness status to `READY_FOR_PHASE11I_JOURNAL_SPECIFIC_FORMATTING_AND_FINAL_QA`.
*   **Alternatives Considered:** None.
*   **Scientific and Operational Justification:** All required figures, tables, bibliography, and the main manuscript have been assembled correctly. Figure and table legends were extracted. A generic cover letter and checklist were created because the journal target is unknown. All components passed structural validation and no scientific text, evidence categories, or target rankings were altered.
*   **Outcome:** Phase 11H successfully completed.
*   **Files / Analyses Affected:** `04_analysis/11_manuscript/PHASE11H_SUBMISSION_PACKAGE_ASSEMBLY.md`, `05_results/tables/phase11h_submission_package_inventory.tsv`, `05_results/tables/phase11h_missing_submission_items.tsv`, `08_submission/phase11h_submission_package/*`, `01_metadata/file_manifest.tsv`, `00_admin/PROJECT_STATUS.md`, `09_docs/planning/DECISION_LOG.md`.

| 2026-07-04 | **D-57** | Complete Phase 11I-A Final QA and Journal-Specific Gap Audit | `04_analysis/11_manuscript/PHASE11IA_FINAL_QA_AND_JOURNAL_GAP_AUDIT.md`, `05_results/tables/phase11ia_*`, `06_scripts/python/20_validate_phase11ia_final_qa.py` |
