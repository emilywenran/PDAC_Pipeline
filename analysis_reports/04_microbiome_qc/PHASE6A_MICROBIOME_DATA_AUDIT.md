# Phase 6A Microbiome Data Audit

## Source and Matrix

- Exact source table: Supplementary Data 1 (42003_2021_2557_MOESM4_ESM.xlsx), sheet 'Genus-level'.
- Extracted taxonomic level: Genus.
- Abundance unit: Kraken2/Bracken-derived non-integer normalized counts as released in Supplementary Data 1; not raw classified reads and not relative abundance.
- Numerical scale: preserved exactly from the public supplementary workbook; no filtering, rarefaction, renormalization, log transform, or CLR transform was applied to the stored matrix.
- Final dimensions: 365 genera x 62 tumor samples.
- Source sample identifier format: `Basal-like1-17`, `Hybrid1-23`, and `Classical1-22`.
- Stored matrix: `03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz`.

## Method Provenance

Repository metadata identifies PRJNA719915 as 62 single-end Illumina shotgun metagenomic tumor runs. The peer-review supplement states that read counts for taxa were measured using Kraken2 plus Bracken and then converted to relative abundance for composition displays; Bracken reallocation was specifically discussed for species estimates. Supplementary Data 1 is released as abundance profiles at class, order, family, genus, and species levels. The exact taxonomic database name/version used by Kraken2/Bracken was not found in the verified local public supplementary files and is marked `TO_VERIFY`.

No unclassified, Homo sapiens, human, host, or unmapped categories were present as genus-level feature rows in the extracted sheet. The matrix contains only named genera, including several zero-only rows.

## Mapping and Validation

- Mapping success: 62/62 matrix columns verified to project patients.
- Unique tumor samples: 62.
- One microbiome profile per patient: yes.
- Duplicate samples: none detected.
- Duplicate taxonomic identifiers: none detected.
- Unmatched patients: 0.
- Negative values: none.
- Missing or infinite values: none.
- Feature/sample orientation: rows are genera and columns are tumor microbiome profiles.

## Descriptive QC

- Overall zero fraction: 0.7570.
- Genera detected in at least 50% of samples: 74.
- Genera detected in at most 10% of samples: 216.
- Total abundance range per sample: 36869.109 to 14200330.124.
- Detected-genera range per sample: 43 to 198.
- Suspected extreme samples by descriptive z-score screening: Basal-like1, Hybrid18, Hybrid23.

Sample-sample structure was summarized with Bray-Curtis distances on the released abundance scale and PCA on log10 relative genus abundance with a small display pseudocount. These are descriptive ordination/QC summaries only and are not subtype, continuous-axis, differential-abundance, survival, host-correlation, pathway, or target-prioritization tests.

Technical metadata-only Spearman checks of total abundance were performed against available `bases`, `spots`, and file size fields:

| technical_metadata   |   n |   spearman_rho_total_abundance |   spearman_p_total_abundance | analysis_scope          |
|:---------------------|----:|-------------------------------:|-----------------------------:|:------------------------|
| bases                |  62 |                       0.451537 |                  0.000229902 | technical_metadata_only |
| spots                |  62 |                       0.45607  |                  0.000195091 | technical_metadata_only |
| file_size_numeric    |  62 |                       0.424089 |                  0.00059259  | technical_metadata_only |

## Contamination-Control Limitations

The project has no sequenced negative-control runs. Therefore no decontam prevalence or frequency analysis was performed, no contaminant feature is confirmed by negative-control evidence, and no potential contaminant was automatically deleted. Phase 6A records contamination risk only. Potential reagent/environment-associated genera flagged in the extracted matrix: Brevundimonas, Paraburkholderia, Cupriavidus, Mesorhizobium, Comamonas, Dechloromonas, Caulobacter, Pseudomonas, Methylobacterium, Elizabethkingia, Bradyrhizobium, Sphingopyxis, Ralstonia, Novosphingobium, Burkholderia, Sphingomonas, Rhizobium, Stenotrophomonas, Acinetobacter, Delftia, Herbaspirillum.

Contamination risk is not equivalent to confirmed contamination. Phase 6B should compare results under multiple prevalence and abundance filters and should keep flagged genera visible in sensitivity reports.

## Suitability for Compositional Analysis

The extracted matrix is suitable for exploratory compositional preprocessing evaluation because it is non-negative and contains 62 complete tumor profiles. It is not directly suitable for CLR/Aitchison analysis without an explicit zero-handling policy because zeros are common. CLR must not be applied directly to raw zeros.

## Recommended Phase 6B Preprocessing Candidates

- Prevalence filtering: compare thresholds such as detected in at least 5%, 10%, and 20% of samples, with all thresholds reported.
- Abundance filtering: compare low-total-abundance removal thresholds independent of subtype labels and continuous scores.
- Pseudocount policy: use a documented small pseudocount or multiplicative replacement after filtering; perform sensitivity to the pseudocount choice.
- CLR/Aitchison transformation: apply only after zero handling; use Aitchison distances for compositional sensitivity analyses.
- Count-based methods: only if raw or integer Bracken estimated counts are obtained or regenerated later; do not treat the released non-integer normalized matrix as raw counts.
- Compositional differential-abundance methods: evaluate ALDEx2, ANCOM-BC, or related approaches in a later locked phase, with contamination-risk sensitivity and technical covariates.

## TO_VERIFY

- Exact Kraken2/Bracken database and version used by the original study.
- Exact normalization formula that produced the non-integer Supplementary Data 1 abundance scale.
- Whether zero-only genus rows are intentional retained taxa from the original pipeline or workbook artifacts.
