# Phase 11G-R3: Reference Repair

## 1. Overview
This report documents the repairs applied to references identified as failing during the Phase 11G-R2 verification audit. No scientific content or biological analyses were modified. The changes were strictly limited to correcting metadata inconsistencies and removing a fabricated reference.

## 2. Methodology
- Extracted the four failed references from the 11G-R2 audit ([3], [6], [7], [13]).
- Searched CrossRef and PubMed for the correct PMIDs associated with the DOIs for `[3]`, `[6]`, and `[7]`.
- Replaced the PMIDs in `05_results/tables/phase11gr1_reference_mapping.tsv` with the true PMIDs.
- Addressed reference `[13]` (Boehm 2015), which was found to be fabricated/mismatched. 
- Removed `[13]` from `09_docs/references/phase11g_references.bib` and the mapping table.
- Removed citation `[13]` from the manuscript, modifying `[12,13]` to `[12]`.
- Renumbered citation `[14]` to `[13]` in the manuscript and mapping table to maintain contiguous numbering.

## 3. Results

### 3.1 Metadata Corrections
- **[3] Bear 2020**: The incorrect PMID `32673551` was replaced with the verified PMID `32946773`. The reference is now fully `PASS`.
- **[6] Nejman 2020**: The incorrect PMID `32471947` was replaced with the verified PMID `32467386`. The reference is now fully `PASS`.
- **[7] Geller 2017**: The incorrect PMID `28916614` was replaced with the verified PMID `28912244`. The reference is now fully `PASS`.

### 3.2 Reference Removal
- **[13] Boehm 2015**: Since the reference was fabricated and no verified source was available or required to support the general statement about target prioritization (already supported by `[12]`), the reference was removed. The claim support status for this action is marked as `PROJECT_INTERNAL_METHOD`. The manuscript citation was updated to reflect the removal of `[13]` and the subsequent renumbering of `[14]` to `[13]`.

### 3.3 Manuscript State
- All citation keys resolve to valid, verified references.
- No text changes were made to the manuscript beyond citation number adjustments (`[12,13]` -> `[12]` and `[14]` -> `[13]`).
- Evidence categories (`PARTIAL_SPATIAL_SUPPORT`, `CELL_COMPOSITION_EXPLAINED`) and target rankings remain intact.

## 4. Final Decision
**PASS**. All failing references from 11G-R2 have been resolved, and the manuscript contains no unsupported or unverified citations.
