# Phase 6C Analysis-Ready Tumor Microbiome Matrices

## Scope

Phase 6C executed the locked Phase 6B tumor microbiome preprocessing protocol for PRJNA719915 using only the audited genus abundance matrix, microbiome sample crosswalk, Phase 6A technical QC outputs, and Phase 6B method-lock files. No host-expression scores, continuous-axis results, survival outcomes, subtype comparisons, differential-abundance results, host-microbiome correlations, pathway analyses, or target-prioritization outputs were used to select preprocessing parameters.

## Primary Filtering Result

The primary rule retained genera with abundance strictly greater than 0 in at least 20% of the 62 samples, corresponding to a minimum detection count of 13 samples. The resulting matrix retained **122 genera** across **62 samples**. This matches the Phase 6B expected retained feature count of approximately 122 genera; no forced feature count adjustment was applied.

The zero fraction before zero replacement in the primary filtered matrix was **0.3556**. No sample became degenerate after filtering.

## Primary Retained Genera

Lysobacter, Brevundimonas, Paraburkholderia, Aerococcus, Acidovorax, Cupriavidus, Staphylococcus, Luteimonas, Rhodopseudomonas, Rhodoferax, Blastomonas, Aeromonas, Achromobacter, Azoarcus, Enterobacter, Methyloversatilis, Mycobacterium, Rubrivivax, Mesorhizobium, Nocardioides, Comamonas, Dechloromonas, Chelatococcus, Bordetella, Kocuria, Bosea, Phenylobacterium, Xanthomonas, Azospira, Deinococcus, Hydrogenophaga, Porphyrobacter, Aromatoleum, Microbacterium, Mitsuaria, Pandoraea, Rhodococcus, Streptococcus, Rhizorhabdus, Citromicrobium, Massilia, Candida, Ramlibacter, Caulobacter, Bacteroides, Enterococcus, Laribacter, Variovorax, Pseudomonas, Methylobacterium, Elizabethkingia, Tessaracoccus, Agrobacterium, Azospirillum, Clostridium, Ensifer, Roseomonas, Croceicoccus, Rhodobacter, Micrococcus, Bradyrhizobium, Aquabacterium, Prevotella, Sulfuritalea, Escherichia, Moraxella, Shinella, Cutibacterium, Chryseobacterium, Aminobacter, Sphingopyxis, Thiomonas, Flavobacterium, Verminephrobacter, Citrobacter, Sphingobium, Methylorubrum, Ralstonia, Erythrobacter, Xanthobacter, Paucibacter, Streptomyces, Collimonas, Bacillus, Polaromonas, Ochrobactrum, Fusarium, Serratia, Candidatus Accumulibacter, Methylibium, Novosphingobium, Pseudoxanthomonas, Burkholderia, Rhizobacter, Sphingomonas, Diaphorobacter, Fusobacterium, Sphingobacterium, Alicycliphilus, Devosia, Ottowia, Thauera, Phreatobacter, Dickeya, Rhizobium, Sinorhizobium, Stenotrophomonas, Corynebacterium, Leptothrix, Altererythrobacter, Roseolovirus, Janibacter, Lactobacillus, Acinetobacter, Janthinobacterium, Paracoccus, Delftia, Sphingorhabdus, Klebsiella, Brachybacterium, Herbaspirillum, Melaminivora

## Primary Pseudocount And Justification

The primary CLR representation used the locked fixed pseudocount **0.889651**, derived in Phase 6B as one half of the minimum non-zero value in the audited matrix. This value is source-matrix specific and remains marked as non-transferable to external cohorts.

## CLR And Aitchison Validation

The primary CLR matrix contains no missing or infinite values. Per-sample CLR column sums were within floating-point tolerance, with maximum absolute column sum **1.084e-13**. The primary Aitchison distance matrix is 62 x 62, symmetric, has a zero diagonal, and contains no negative or infinite distances. The primary CLR value range was **-4.5364** to **13.2533**.

## Sensitivity Representations Produced

Sensitivity outputs were generated under `03_processed/microbiome/sensitivity/` for these analysis IDs: MICRO_SENS_PREV_10, MICRO_SENS_PREV_30, MICRO_SENS_DET_10_P20, MICRO_SENS_PSEUDO_0.1, MICRO_SENS_PSEUDO_1.0, MICRO_SENS_ROBUST_CLR, MICRO_SENS_NO_HIGH_RISK, MICRO_SENS_NO_CONTAMINANTS, MICRO_SENS_EXCLUDE_EXTREME, MICRO_SENS_PRESENCE_ABSENCE. These cover the locked 10% and 30% prevalence filters, abundance threshold >10 at 20% prevalence, pseudocounts 0.1 and 1.0, robust CLR, high-risk contaminant exclusion, high- plus moderate-risk contaminant exclusion, presence/absence Jaccard representation, and exclusion of the three Phase 6A technical extreme samples.

## Contamination-Flag Handling

Potential contaminant genera were retained in the primary matrix. Categories in `phase6c_retained_taxa_with_contamination_flags.tsv` are evidence flags, not confirmed contamination labels. Among primary retained genera, **21** had non-low contamination-sensitivity flags. The high- and moderate-risk removal analyses were generated only as sensitivity representations.

## Technical Outlier Handling

The three locked technical extreme samples (`Basal-like1`, `Hybrid18`, `Hybrid23`) were retained in the primary analysis. The `MICRO_SENS_EXCLUDE_EXTREME` sensitivity representation excludes those samples and was used only to evaluate preprocessing robustness.

## Preprocessing Robustness

Outcome-blind sensitivity concordance was assessed using upper-triangle distance-matrix Spearman correlations, Procrustes concordance on PCoA coordinates, sample-order stability, taxon-rank stability, CLR correlations across pseudocount choices, and ordination shifts after contaminant removal. The lowest non-missing distance-matrix correlation across sensitivity comparisons was **0.5000**. Full results are in `05_results/tables/phase6c_preprocessing_sensitivity_concordance.tsv`.

## Remaining Limitations

Sequenced negative controls are unavailable, so contaminant categories remain potential-risk flags rather than confirmed contaminant calls. The abundance scale remains Bracken-normalized non-integer estimates rather than raw integer classified-read counts. Matrix total abundance is therefore a technical proxy, not a direct absolute microbial load measure.

## Proceed Decision

Microbiome association method locking may proceed with contamination limitations, provided Phase 7 continues to use the locked sensitivity framework and does not treat flagged genera as confirmed contaminants.

## TO_VERIFY

- TO_VERIFY: Original publication/source-specific normalization formula for the Bracken-derived abundance estimates.
- TO_VERIFY: Sequencing batch covariates remain unavailable in public metadata.
- TO_VERIFY: Tumor purity, immune/stromal scores, and clinical covariates remain unavailable for Phase 6C preprocessing.
