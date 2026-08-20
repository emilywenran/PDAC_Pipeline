# Phase 10B-R Corrected Target Prioritization Results

## Scope

Phase 10B-R re-executes cross-layer evidence synthesis after the Phase 10C rejection of the first Phase 10B attempt. The reanalysis uses the Phase 10A method lock, the cross-layer synthesis protocol, the Phase 10A parameter inventory, the Phase 10A evidence inventory, the locked target-prioritization framework, and the final PASS reviews from Phases 7C, 8C, 9B1C2, 9B2C2, and 9B3C2 as authoritative inputs.

No manuscript drafting was performed.

## Corrections

- Removed descriptive target overrides from scoring. Candidate rows are generated from `phase10a_cross_layer_evidence_inventory.tsv`.
- Scored every Phase 10A evidence-inventory candidate. Skipped candidates: `NONE`.
- Treated `CELL_COMPOSITION_EXPLAINED` as a high-weight penalty for cell-type specificity, including CTCFL/BORIS.
- Did not use literature support to add points or rescue weak/composition-sensitive candidates.
- Preserved the locked negative and partial evidence for Ochrobactrum, HALLMARK_PROTEIN_SECRETION, spatial support, and ineligible WGCNA modules.

## External Database Reproducibility

OpenTargets, GTEx, and ChEMBL values are not filled from dictionaries. Phase 10B-R generated a local audit table for every candidate-database pair. Gene-symbol database queries without local reproducible result tables are marked `NOT_RUN_DATABASE_UNAVAILABLE`; 6 such rows were recorded.

## Candidate Scores

| Candidate | Derived evidence class | Single-cell evidence | Target score | Decision |
|---|---|---|---:|---|
| HALLMARK_PROTEIN_SECRETION | MULTI_LAYER_SUPPORTED | MALIGNANT_CELL_INTRINSIC_SUPPORT | 5 | RETAINED_AS_SUPPORTED_BIOLOGICAL_FEATURE_NOT_DIRECT_GENE_TARGET |
| BHLHE40 | PARTIALLY_REPLICATED | CELL_COMPOSITION_EXPLAINED | 1 | NOT_PRIORITIZED_COMPOSITION_SENSITIVE_SINGLE_CELL_EVIDENCE |
| CTCFL | PARTIALLY_REPLICATED | CELL_COMPOSITION_EXPLAINED | 1 | NOT_PRIORITIZED_COMPOSITION_SENSITIVE_SINGLE_CELL_EVIDENCE |
| HALLMARK_SPERMATOGENESIS | PARTIALLY_REPLICATED | CELL_COMPOSITION_EXPLAINED | 1 | NOT_PRIORITIZED_COMPOSITION_SENSITIVE_SINGLE_CELL_EVIDENCE |
| Ochrobactrum | DISCOVERY_ONLY | NOT_APPLICABLE | 0 | NOT_PRIORITIZED_DISCOVERY_ONLY_NO_CAUSAL_VALIDATION |
| Lysobacter | METHOD_SENSITIVE | NOT_APPLICABLE | 0 | NOT_PRIORITIZED_LOCKED_EVIDENCE_WEAK_OR_SENSITIVE |
| Brevundimonas | CONTAMINATION_SENSITIVE | NOT_APPLICABLE | 0 | NOT_PRIORITIZED_LOCKED_EVIDENCE_WEAK_OR_SENSITIVE |
| MEpurple | NOT_EXTERNALLY_SUPPORTED | NOT_EVALUATED | 0 | NOT_PRIORITIZED_NOT_EXTERNALLY_SUPPORTED |
| MEred | NOT_EXTERNALLY_SUPPORTED | NOT_EVALUATED | 0 | NOT_PRIORITIZED_NOT_EXTERNALLY_SUPPORTED |
| MEblack | INSUFFICIENT_DATA | INSUFFICIENT_SINGLE_CELL_DATA | 0 | NOT_PRIORITIZED_INELIGIBLE_INSUFFICIENT_DATA |
| MEblue | INSUFFICIENT_DATA | INSUFFICIENT_SINGLE_CELL_DATA | 0 | NOT_PRIORITIZED_INELIGIBLE_INSUFFICIENT_DATA |
| MEgreen | INSUFFICIENT_DATA | INSUFFICIENT_SINGLE_CELL_DATA | 0 | NOT_PRIORITIZED_INELIGIBLE_INSUFFICIENT_DATA |
| MEgreenyellow | INSUFFICIENT_DATA | INSUFFICIENT_SINGLE_CELL_DATA | 0 | NOT_PRIORITIZED_INELIGIBLE_INSUFFICIENT_DATA |
| MEtan | INSUFFICIENT_DATA | INSUFFICIENT_SINGLE_CELL_DATA | 0 | NOT_PRIORITIZED_INELIGIBLE_INSUFFICIENT_DATA |
| Staphylococcus | EXPLORATORY_ONLY | NOT_APPLICABLE | 0 | NOT_PRIORITIZED_LOCKED_EVIDENCE_WEAK_OR_SENSITIVE |
| Paraburkholderia | NO_SUPPORTED_ASSOCIATION | NOT_APPLICABLE | 0 | NOT_PRIORITIZED_LOCKED_EVIDENCE_WEAK_OR_SENSITIVE |

## Key Preserved Evidence

HALLMARK_PROTEIN_SECRETION retains malignant-cell and malignant-compartment support, but the basal-classical spatial-axis association is not replicated and spatial evidence remains `PARTIAL_SPATIAL_SUPPORT`.

CTCFL/BORIS is not promoted using cancer-testis antigen reasoning because reproducible GTEx/OpenTargets/ChEMBL evidence is unavailable in this run and its single-cell evidence is `CELL_COMPOSITION_EXPLAINED`.

Ochrobactrum retains robust host-mechanism association only. Microbial localization, physical interaction, and causality validation remain absent.

The five locked ineligible WGCNA modules remain `INSUFFICIENT_DATA` and are not promoted.

## Audits

- Cross-layer evidence scores: `05_results/tables/phase10br_cross_layer_evidence_scores.tsv`
- Candidate target scores: `05_results/tables/phase10br_candidate_target_scores.tsv`
- External database query audit: `05_results/tables/phase10br_external_database_query_audit.tsv`
- Penalty audit: `05_results/tables/phase10br_penalty_audit.tsv` (3 applied penalties)
- Rank-change audit: `05_results/tables/phase10br_rank_change_audit.tsv`

## Final Readiness Decision

`READY_FOR_PHASE10C2_INDEPENDENT_REVIEW`
