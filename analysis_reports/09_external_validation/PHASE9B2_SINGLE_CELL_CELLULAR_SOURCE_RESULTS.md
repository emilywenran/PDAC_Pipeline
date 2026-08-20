# Phase 9B2 Single-Cell Cellular-Source Results

**SUPERSEDED_BY_PHASE9B2R.** This report is preserved only as an audit artifact after the Phase 9B2C `FAIL_REQUIRES_REANALYSIS` decision. The module-based cellular-source claims, module malignant-axis associations, negative-control records, and blanket TF control categories in this file must not be cited as current evidence. Use `PHASE9B2R_CORRECTED_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md` and `phase9b2r_*` outputs instead.

## Status

**PRIMARY_RUN_COMPLETE_WITH_TO_VERIFY_CONTROLS**

Phase 9B2 was restarted and executed on 2026-07-03 as a single-cell cellular-source evaluation of bulk-externally replicated TF activities and partially replicated or discovery-supported host programs. The primary execution cohort was restricted to `PENG_CRA001160`; no LIN, Moncada, Hwang, spatial, microbiome, survival, causal mediation, target-prioritization, or manuscript-writing analyses were performed.

## Authoritative Dataset Identity

The executed cohort is `PENG_CRA001160`, accession `CRA001160`, BioProject `PRJCA001063`, Peng et al. 2019, human PDAC scRNA-seq. The cohort contains 24 untreated PDAC tumors and 11 control pancreases. The official source used for processed data acquisition was CNCB GSA: `https://download.cncb.ac.cn/gsa/CRA001160/`.

Startup validation passed for canonical ID, accession, BioProject, publication, cohort size, and primary-only inclusion. The Peng cohort was not labeled as `GSE111672`; `GSE111672` remains assigned to `MONCADA_GSE111672` and was not analyzed in this primary Layer 2 run.

## Data Acquired

Only processed data were downloaded. No FASTQ, BAM, SRA, or other raw sequencing files were acquired.

| File | Role | Size |
|---|---:|---:|
| `02_data/external/phase9_single_cell/PENG_CRA001160/count-matrix.txt` | processed expression matrix | 2,771,872,913 bytes |
| `02_data/external/phase9_single_cell/PENG_CRA001160/all_celltype.txt` | published cell annotations | 2,101,436 bytes |
| `02_data/external/phase9_single_cell/PENG_CRA001160/md5sum.txt` | official checksum listing | 7,156 bytes |

Checksums and source metadata are recorded in `01_metadata/phase9b2_single_cell_dataset_inventory.tsv` and `01_metadata/file_manifest.tsv`.

## Cohort and QC

The processed matrix is gene-by-cell and was streamed into patient-level and patient-cell-type pseudobulks. The analysis retained 57,530 cells from 35 patients: 24 PDAC tumor patients and 11 control pancreas donors. Patient identifiers were retained, and primary inference used patient-level pseudobulk or patient-level aggregate scores.

Published major cell annotations were reviewed at the broad-class level. The reviewed major categories were malignant epithelial, nonmalignant epithelial, fibroblast/CAF, endothelial, myeloid, T cell, and B cell. No high-resolution subtype reannotation was performed.

Cell counts by reviewed major category were:

| Reviewed category | Cells |
|---|---:|
| nonmalignant epithelial | 12,981 |
| fibroblast/CAF | 12,649 |
| malignant epithelial | 11,315 |
| endothelial | 9,117 |
| myeloid | 5,361 |
| T cell | 3,660 |
| B cell | 2,447 |

Malignant-cell calls used source annotations plus conservative epithelial context. Tumor `Ductal type 2` cells were classified as `MALIGNANT`; tumor `Ductal type 1` cells were retained as `AMBIGUOUS`; acinar/endocrine and control epithelial cells were classified as `NONMALIGNANT_EPITHELIAL`; non-epithelial cells were `NOT_APPLICABLE`. Not every epithelial or ductal cell was classified as malignant.

## Feature Coverage and Scoring

Coverage was calculated before scoring for Moffitt50, Moffitt49 no-LEMD1, PurIST gene pairs, the two locked Hallmark programs, five transferred WGCNA modules, and DoRothEA A/B/C regulons for the 25 selected TFs.

Primary scoring used patient-cell-type pseudobulk log2-CPM profiles. Hallmark scores used the locked `decoupleR::run_gsva` single-sample implementation. TF activities used DoRothEA A/B/C regulons and `decoupleR::run_viper`; no TF-symbol expression proxy was used. WGCNA module transfer used exact Phase 8 discovery gene membership and did not reconstruct single-cell WGCNA modules.

## Cellular Sources

The evidence table contains all 32 locked features: two Hallmark programs, five transferred modules, twelve bulk-externally replicated TF activities, and thirteen partially replicated TF activities.

Final single-cell evidence categories were:

| Category | Features |
|---|---:|
| `STROMAL_OR_IMMUNE_SOURCE_SUPPORTED` | 14 |
| `CELL_COMPOSITION_EXPLAINED` | 10 |
| `PARTIAL_CELLULAR_SUPPORT` | 5 |
| `MALIGNANT_CELL_INTRINSIC_SUPPORT` | 3 |

Features with malignant-intrinsic support in this primary run were `HALLMARK_PROTEIN_SECRETION`, `MEblue`, and `MEgreen`. Stromal or immune localization was common among TF activities and included B-cell, T-cell, fibroblast/CAF, myeloid, and endothelial-dominant signals. Several externally replicated TF activities were composition-sensitive rather than clearly malignant-cell intrinsic.

Full per-feature classifications are in `05_results/tables/phase9b2_cellular_source_evidence.tsv`; model outputs are in `05_results/tables/phase9b2_cellular_source_models.tsv`.

## Malignant-Cell Axis Associations

Within malignant-cell patient pseudobulks, the following features remained associated with malignant-cell Moffitt50 basal-classical contrast at BH q < 0.10:

| Feature | Layer | Direction |
|---|---|---|
| `HALLMARK_PROTEIN_SECRETION` | Hallmark | positive |
| `MEblack` | transferred WGCNA module | positive |
| `MEblue` | transferred WGCNA module | positive |
| `MEgreen` | transferred WGCNA module | negative |
| `MEtan` | transferred WGCNA module | negative |

No TF regulon activity passed the malignant-cell axis association threshold in this run. Continuous malignant-cell state heterogeneity was retained; no new single-cell Hybrid threshold was optimized.

## Cell-Composition Sensitivity

Cell-fraction sensitivity models found composition associations across multiple features. At q < 0.10, significant composition covariates included malignant epithelial fraction, endothelial fraction, lymphoid fraction, fibroblast/CAF fraction, and myeloid fraction. These results support the interpretation that some bulk-supported host signals could reflect cellular composition rather than malignant-cell intrinsic activity.

## Tumor-Control Context

The 11 control pancreases were used only as a contextual comparison. Tumor-versus-control descriptive results were not used to redefine the locked feature list or evidence rules.

## Negative Controls and TO_VERIFY Items

Negative-control records were generated for randomized same-size module controls, expression-matched module controls, unrelated Hallmark pathways, patient-label permutations, and cell-type-label permutations. The computationally intensive expression-matched and permutation controls remain `TO_VERIFY` and should be independently reviewed before using Phase 9B2 for final evidentiary claims beyond this single-cohort cellular-source evaluation.

## Interpretation Boundary

Phase 9B2 evaluated cellular source, malignant-cell specificity, stromal or immune contributions, cell-composition sensitivity, and patient-level malignant-state heterogeneity of host transcriptional programs. It did not test Ochrobactrum presence, microbial abundance, microbial localization, microbial causality, or microbe-cell physical interactions. No result in this report is presented as microbial replication evidence, experimental validation, or causal evidence.

Spatial validation remains recommended to test localization and tissue architecture in an orthogonal modality.

## Key Outputs

Core tables:

- `05_results/tables/phase9b2_restart_runtime_validation.tsv`
- `01_metadata/phase9b2_single_cell_dataset_inventory.tsv`
- `05_results/tables/phase9b2_single_cell_cohort_qc.tsv`
- `05_results/tables/phase9b2_cell_annotation_audit.tsv`
- `05_results/tables/phase9b2_malignant_cell_audit.tsv.gz`
- `05_results/tables/phase9b2_single_cell_feature_coverage.tsv`
- `05_results/tables/phase9b2_pseudobulk_inventory.tsv`
- `05_results/tables/phase9b2_patient_celltype_state_scores.tsv`
- `05_results/tables/phase9b2_patient_celltype_host_program_scores.tsv`
- `05_results/tables/phase9b2_patient_celltype_tf_activity.tsv`
- `05_results/tables/phase9b2_cellular_source_models.tsv`
- `05_results/tables/phase9b2_malignant_feature_axis_associations.tsv`
- `05_results/tables/phase9b2_cell_composition_sensitivity.tsv`
- `05_results/tables/phase9b2_negative_control_results.tsv`
- `05_results/tables/phase9b2_cellular_source_evidence.tsv`

Figures:

- `05_results/figures/phase9b2_cohort_cell_counts.pdf`
- `05_results/figures/phase9b2_cell_annotation_markers.pdf`
- `05_results/figures/phase9b2_malignant_cell_audit.pdf`
- `05_results/figures/phase9b2_moffitt_axis_by_cell_type.pdf`
- `05_results/figures/phase9b2_malignant_axis_by_patient.pdf`
- `05_results/figures/phase9b2_hallmark_cellular_source.pdf`
- `05_results/figures/phase9b2_module_cellular_source.pdf`
- `05_results/figures/phase9b2_tf_activity_cellular_source.pdf`
- `05_results/figures/phase9b2_malignant_feature_axis_heatmap.pdf`
- `05_results/figures/phase9b2_cell_composition_sensitivity.pdf`
- `05_results/figures/phase9b2_tumor_control_descriptive.pdf`
- `05_results/figures/phase9b2_negative_control_summary.pdf`
- `05_results/figures/phase9b2_cellular_source_evidence_summary.pdf`

## Review Readiness

Phase 9B2 is ready for independent review as a primary single-cohort, patient-aware cellular-source analysis of `PENG_CRA001160`. Review should focus on the conservative malignant-cell strategy, the patient-level pseudobulk modeling assumptions, the source-study annotation audit, and the `TO_VERIFY` negative-control items.
