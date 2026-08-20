# Phase 11A Manuscript Claim Map

This document outlines the strict claim map for the PDAC host-microbiome manuscript, categorized according to the validated evidence levels.

## 1. Supported Findings
* Global community-level association between tumor microbiome composition and host state (PERMANOVA $R^2 = 0.0534$, $P = 0.0001$).
* 9 robust genera (Azoarcus, Candida, Ensifer, Cutibacterium, Chryseobacterium, Ochrobactrum, Burkholderia, Rhizobium, Herbaspirillum) are robustly associated with host state.
* 43 robust biological mechanisms associated with Ochrobactrum.
* 12 TF activities externally replicated in bulk (CTCFL, IRF3, JUNB, KLF13, KLF9, MNT, MXI1, SNAI2, TFAP4, TP63, ZBTB7A, ZNF24).
* HALLMARK_PROTEIN_SECRETION has malignant-cell intrinsic and malignant-compartment spatial support.

## 2. Partially Supported Findings
* Spatial validation of the basal–classical axis is PARTIAL_SPATIAL_SUPPORT.
* 7 features are partially replicated in bulk (HALLMARK_PROTEIN_SECRETION, HALLMARK_SPERMATOGENESIS, MEblack, MEblue, MEgreen, MEtan, MEgreenyellow).
* 13 TFs are partially replicated in bulk (BHLHE40, E2F6, ELF1, GRHL2, KLF1, MBD1, MBD2, OTX2, SIX5, SNAPC4, ZBED1, ZNF384, ZNF740).

## 3. Unsupported or Null Findings
* HALLMARK_SPERMATOGENESIS is not externally replicated in bulk (direction opposite to discovery).
* WGCNA modules MEred and MEpurple are not replicated.
* 9 TFs are not replicated in bulk (GFI1B, STAT1, ZBTB11, ZNF639, TWIST1, FOXK2, KDM5B, MAFF, TEAD4).
* WGCNA modules have insufficient single-cell data (<80% coverage) and thus are not formally supported.
* HALLMARK_SPERMATOGENESIS has insufficient spatial data (<80% coverage).

## 4. Exploratory Findings
* Moncada spatial cohort analysis is exploratory (1/6 sections positive, 5/6 sections non-significant or negative).

## 5. Contamination-Sensitive Findings
* Herbaspirillum is flagged as a moderate-risk environmental contaminant.
* 21 genera are categorized as contamination-sensitive.

## 6. Composition-Sensitive Findings
* 8 of 9 robust genera (Azoarcus, Candida, Ensifer, Cutibacterium, Chryseobacterium, Burkholderia, Rhizobium, Herbaspirillum) show direction reversals or loss of significance under rCLR, making them transformation-sensitive.
* 20 features in single-cell data, including CTCFL/BORIS, are CELL_COMPOSITION_EXPLAINED (composition-sensitive).
* HALLMARK_PROTEIN_SECRETION is composition-sensitive but remains supported as malignant-intrinsic because it is also axis-associated and localized to malignant cells.

## 7. Findings that Must Not Be Stated Causally
* Any association between microbiome and host expression. Do not use causal language (e.g., "drives", "causes", "induces", "mediates").
* Do not claim microbial localization or physical interaction (no sequenced negative controls).
