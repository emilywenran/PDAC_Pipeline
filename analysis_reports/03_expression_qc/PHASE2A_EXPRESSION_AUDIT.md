# Phase 2A Expression Audit: GSE172356

## Source and Skill Usage

- Required K-Dense skill loaded: `~/.agents/skills/exploratory-data-analysis/SKILL.md`.
- Applicable skill workflow: tab-delimited scientific data ingestion with row/column counts, type inference, missingness checks, duplicate checks, outlier detection, and correlation/PCA summaries.
- Official processed expression source: `GSE172356_PDA_gene_expression_matrix.txt.gz`.
- Source URL/accession: `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE172356&format=file&file=GSE172356%5FPDA%5Fgene%5Fexpression%5Fmatrix%2Etxt%2Egz`; GEO Series accession `GSE172356`.
- GEO SOFT evidence URL: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE172nnn/GSE172356/soft/GSE172356_family.soft.gz`.
- Download date: `2026-06-30`.
- Source file size: `17364266` bytes.
- Source MD5: `41e6f2e842c9b8cd02eb6665856aabe6`.
- Source SHA256: `861fc9af648e5b77e1bbf6e5f6a0d61711ccb4bc216a714dd327135a8ee45725`.
- Compression format: gzip.

## Matrix Structure

- Final audited matrix: `03_processed/expression/GSE172356_expression_audited.tsv.gz`.
- Matrix dimensions: 45140 gene rows x 62 expression samples, plus one `gene` identifier column.
- Orientation: genes are rows; samples are columns.
- Sample identifier format: `YX...T` tumor aliases, for example `YX15261T`.
- Gene identifier type: gene symbols/gene names from the GEO matrix first column.
- Gene identifiers unique: True.
- Duplicate gene rows: 0; duplicate gene IDs beyond first occurrence: 0.
- Duplicate samples: 0.

## Expression Unit

Expression unit is **DESeq size-factor-normalized counts**. Supporting evidence: official GEO SOFT metadata states reads were aligned to GRCh38 with HISAT2 v2.1.0, raw transcript counts were calculated using HTSeq v0.12.4, and normalization was performed by SizeFactors using DESeq v1.24.0 in R. The matrix contains non-integer fractional values, supporting normalized rather than raw integer counts. It is not FPKM, TPM, CPM, or log-transformed based on available GEO evidence.

## Value Integrity

- Missing values: 73202.
- Non-numeric values in expression cells: 0.
- Infinite values: 0.
- Negative values: 0.
- Zero values: 329457.
- All-zero genes: 15.
- Minimum value: 0.0.
- Maximum value: 4474617.58125138.

## Sample Mapping

- Expression columns evaluated: 62.
- Successfully mapped samples: 62 / 62.
- Mapping output: `01_metadata/expression_sample_crosswalk.tsv`.
- Mapping rule: expression `YX...T` aliases were matched to the 62-patient manifest aliases recorded in `sample_manifest.tsv`; GEO sample IDs and patient IDs were retained from the finalized manifest.

## Descriptive QC

Outputs:

- Per-sample QC: `05_results/tables/phase2a_expression_sample_qc.tsv`.
- Per-gene QC: `05_results/tables/phase2a_expression_gene_qc.tsv`.
- Mapping summary: `05_results/tables/phase2a_expression_mapping_summary.tsv`.
- Distribution figure: `05_results/figures/phase2a_expression_distribution.pdf`.
- Sample correlation figure: `05_results/figures/phase2a_sample_correlation.pdf`.
- PCA figure: `05_results/figures/phase2a_expression_pca.pdf`.

PCA was performed for QC visualization only on `log10(expression + 1)` values after visualization-only gene-median imputation for missing cells and gene-wise scaling. This transformed PCA input was not written as the audited matrix and was not used for subtype discovery, feature selection, differential expression, or biological interpretation. PC1 explained 0.1416 and PC2 explained 0.0689 of variance in this QC space.

## Prespecified Extreme-Sample Criteria

Samples were flagged if any of the following objective criteria were met:

- absolute robust z-score for total expression > 3.5;
- robust z-score for detected genes < -3.5;
- mean sample correlation below Q1 - 3 x IQR;
- absolute robust z-score for PC1 or PC2 > 3.5 in the QC PCA space.

Suspected outliers:

expression_column geo_sample_id patient_id                                                              outlier_flag
         YX16135T    GSM5253102   PDAC_016 total_expression_robust_z_abs_gt_3.5;mean_correlation_below_Q1_minus_3IQR
         YX16158T    GSM5253109   PDAC_023                                            pc1_or_pc2_robust_z_abs_gt_3.5
         YX16194T    GSM5253119   PDAC_033                                      total_expression_robust_z_abs_gt_3.5
         YX16224T    GSM5253125   PDAC_039                                      total_expression_robust_z_abs_gt_3.5

No outliers were removed.

## Transformations and Phase 2B Readiness

The audited matrix preserves the original GEO numerical scale and requires no transformation for archiving or sample mapping. For downstream PCA, clustering, subtype reproduction, and regression diagnostics, variance-stabilizing or log-like transformations should be considered and explicitly scripted in the appropriate later phase because the current normalized-count scale is strongly right-skewed. No normalization, batch correction, sample removal, subtype-driven gene selection, differential expression, or biological interpretation was performed in Phase 2A.

Phase 2B may proceed after human review of the generated files and matrix dimensions. TO_VERIFY: confirm with the project owner that DESeq size-factor-normalized counts are the intended input scale for downstream subtype reproduction, or whether raw-count reprocessing from FASTQ is required in a later phase.
