# Project Status

Chronological, phase-by-phase summary of completed work, Phase 1 through Phase 11I-A. Final readiness decision: **`READY_FOR_TARGET_JOURNAL_SELECTION`** (as of 2026-07-04). Next approved task: Target Journal Selection (Phase 11I-B).

## Phase 1 — Data source and provenance

- **Phase 1A (Source audit).** Audited GSE172356 (host expression) and PRJNA719915/PRJNA723013/SRP315394 (tumour microbiome) via NCBI E-utilities and ENA. Reconciled sample/run counts (62 each) and established 1-to-1 patient mapping via `source_material_id`. Confirmed microbiome library layout is single-end WGS shotgun metagenomics. Found that negative controls were only assayed by PCR/gel, not sequenced — no raw control reads exist in SRA.
- **Phase 1B (Patient-level mapping).** Finalized `sample_manifest.tsv` (62 patients), `clinical_metadata.tsv`, and `rna_microbiome_patient_crosswalk.tsv`. Verified subtype counts (Basal 17, Classical 22, Hybrid 23) and survival availability (53/62 patients).

## Phase 2 — Bulk expression QC and preprocessing

- **Phase 2A (Expression audit).** Downloaded the official GEO processed matrix only (no FASTQ/SRA). Verified 45,140 genes × 62 samples; detected 73,202 missing cells (confirmed literal `NA` strings), 0 duplicates, 15 all-zero genes; flagged 4 suspected extreme samples (none removed).
- **Phase 2B (Analysis-ready matrix).** Applied complete-observation filtering, retaining 42,654 genes × 62 samples. Produced filtered normalized-count and log2 matrices.

## Phase 3 — Subtype reproduction

- **Phase 3A (Method lock).** Audited GSE172356-original, Moffitt, PurIST, Bailey, and Chan-Seng-Yue subtyping frameworks; locked signatures, coefficients, and thresholds.
- **Phase 3B (Reproduction).** Reproduced the primary Chan-Seng-Yue 94-gene subtype labels exactly for all 62 patients (Basal 17, Hybrid 23, Classical 22; exact agreement = 1.0 on all metrics). Applied Moffitt and PurIST as verified secondary methods.

## Phase 4 — Subtype stability

- **Phase 4A (Method lock).** Locked 12 stability metrics and 8 primary/sensitivity runs (K=2..6, 1,000 resampling iterations).
- **Phase 4B (Execution).** Overall discrete-subtype stability decision: **`INCONCLUSIVE`** (conflict between CSY-preferred K=2 and independent-HVG-preferred K=4, plus transformation sensitivity).

## Phase 5 — Continuous basal-classical axis

- **Phase 5A (Method lock, amended).** Locked the continuous-axis framework; corrected a Moffitt gene-set inconsistency (LEMD1) by creating separate 50-gene primary and 49-gene no-LEMD1 sensitivity signatures (D-14).
- **Phase 5B (Execution).** Ran all 7 locked continuous-axis analyses. Overall decision: **`INCONCLUSIVE`**.

## Phase 6 — Tumour microbiome QC and preprocessing

- **Phase 6A (Data audit).** Extracted the published genus-level abundance matrix (365 genera × 62 samples) and built the microbiome sample crosswalk.
- **Phase 6B (Method lock, amended).** Locked CLR/rCLR preprocessing, contamination-sensitivity framework, and MaAsLin2 configuration. Status: `READY_WITH_CONTAMINATION_LIMITATIONS`.
- **Phase 6C (Analysis-ready matrices).** Applied the locked prevalence filter (>0 in ≥20% of samples), retaining 122 genera × 62 samples. Generated primary and 9 sensitivity-variant matrices.

## Phase 7 — Host–microbiome association

- **Phase 7A (Method lock).** Locked the OLS/HC3 association framework (122 genus-level tests), PERMANOVA/PERMDISP, and evidence-classification categories. Model 2 blocked (0 complete clinical cases). Status: `READY_WITH_COVARIATE_LIMITATIONS`.
- **Phase 7A.5 (Host covariates).** Computed ESTIMATE-derived purity/immune/stromal scores for all 62 patients; permitted Models 3P/3I/3S as sensitivity-only due to collinearity.
- **Phase 7B (Execution).** Ran the full locked association battery (OLS HC3, PERMANOVA, Spearman, bootstrap, sensitivity, presence/absence models).
- **Phase 7C (Independent review).** Verified all 33 primary FDR-significant candidates and the global PERMANOVA result. Decision: **PASS**.

## Phase 8 — Host-mechanism analysis

- **Phase 8A (Method lock).** Locked 5 host feature layers (Hallmark, PROGENy, DoRothEA/VIPER, WGCNA, genome-wide) and circularity safeguards. Status: `READY_WITH_TRANSFORMATION_LIMITATIONS`.
- **Phase 8A.5 (Environment setup).** Built and validated the project `renv` R environment for all required packages.
- **Phase 8B (Execution).** Ran host-mechanism analyses for the 9 Phase 7C-verified primary taxa.
- **Phase 8C (Independent review).** Verified all 43 robust host-mechanism evidence rows. Decision: **PASS** — authorized Phase 9 external-validation planning.

## Phase 9 — External validation

- **Phase 9A (Planning).** Evaluated 12 candidate external datasets across bulk, single-cell, spatial, and microbiome layers; shortlisted 7 PRIORITY_1 datasets. Status: `READY_WITH_MICROBIOME_LIMITATIONS`.
- **Phase 9B1 → 9B1C → 9B1R → 9B1C2 (External bulk validation).** Initial Phase 9B1 execution (TCGA-PAAD 178, GSE71729 145, GSE62452 69 tumour-only samples) was independently reviewed at Phase 9B1C and returned **`FAIL_REQUIRES_REANALYSIS`** (6 implementation findings, incl. 1 Critical). Phase 9B1R corrected PurIST scoring, Hallmark ssGSEA, and TF activity derivation, and reran negative controls. Phase 9B1C2 verified all corrections and issued a final **PASS** (12 TFs externally replicated, 13 partially replicated, 9 not replicated).
- **Phase 9A.1 / 9A.2 (Single-cell provenance correction).** Two mandatory stops (D-33, D-36) halted single-cell data acquisition before download because Layer-2 cohort records disagreed across the Phase 9A inventory, shortlist, and parameter tables. Reconciled records and locked `PENG_CRA001160` as the sole authorized primary cohort (Phase 9A.3).
- **Phase 9B2 → 9B2C → 9B2R → 9B2C2 (Single-cell validation).** Initial Phase 9B2 execution (57,530 cells, 24 PDAC tumours + 11 control pancreases) was independently reviewed at Phase 9B2C and returned **`FAIL_REQUIRES_REANALYSIS`** (unexecuted negative controls, coverage violations). Phase 9B2R reran the analysis with corrected negative controls and coverage handling. Phase 9B2C2 verified all findings resolved and all 64 applicable negative controls executed. Decision: **PASS** — authorized Phase 9B3 spatial planning.
- **Phase 9B3A / A.1 / A.2 (Spatial planning).** Locked the Hwang GeoMx ROI-pairing models (A/B/C) and reclassified Moncada ST as exploratory-only (n=2 patients); enforced a cross-platform matrix-pooling ban.
- **Phase 9B3B → 9B3C → 9B3R → 9B3C2 (Spatial validation).** Initial Phase 9B3B execution (Hwang naive 13 patients/373 segments; Hwang treated 7 patients/197 segments; Moncada 2 patients/6 sections) was independently reviewed at Phase 9B3C and returned **`FAIL_REQUIRES_REANALYSIS`** (hardcoded negative-control placeholders, non-converged Model C retained, ineligible model fit for a low-coverage feature). Phase 9B3R corrected all three issues and reran with real empirical null distributions. Phase 9B3C2 verified the correction. Result for HALLMARK_PROTEIN_SECRETION: strong tumour-compartment localization (Model A), but no significant basal–classical axis association (Models B/C) — evidence category **`PARTIAL_SPATIAL_SUPPORT`**. Moncada: 1 of 6 sections showed positive concordance (exploratory only, no population-level claim).

## Phase 10 — Cross-layer synthesis and target prioritisation

- **Phase 10A (Planning).** Locked a 10-level cross-layer evidence hierarchy and an independent target-prioritisation framework (druggability, dependency, selectivity, pathway position, cell-type specificity, safety).
- **Phase 10B (Execution).** Scored candidates; identified CTCFL as provisionally the most tissue-selective candidate (GTEx Testis median TPM 7.6 vs <0.1 elsewhere).
- **Phase 10C (Independent review).** Found unreproducible external-database queries, post-hoc CTCFL selection circumventing the composition-sensitivity penalty, qualitative (not quantitative) deprioritisation, and omission of HALLMARK_SPERMATOGENESIS from scoring. Decision: **`FAIL_REQUIRES_REANALYSIS`**.
- **Phase 10B-R (Corrected reanalysis).** Removed all hardcoded overrides; scored every Phase 10A candidate programmatically; applied the locked `CELL_COMPOSITION_EXPLAINED` penalty (CTCFL/BORIS, BHLHE40, HALLMARK_SPERMATOGENESIS not promoted); marked unrunnable external-database rows `NOT_RUN_DATABASE_UNAVAILABLE` rather than filling them in.
- **Phase 10C2 (Independent review).** Verified all corrections, confirmed CTCFL/BORIS correctly penalized and not promoted. Decision: **PASS** — ready for manuscript drafting.

## Phase 11 — Manuscript preparation, audit, and submission package

- **Phase 11A (Claim map and outline).** Locked manuscript structure, claim-control boundaries, and figure/table plan.
- **Phase 11B (Draft finalization).** Produced the manuscript draft and a claim-to-text traceability table mapping all 16 major claims to source evidence and allowed/prohibited wording. AI-assisted drafting used the K-Dense `scientific-writing` skill under the claim-trace constraints locked in Phase 11A.
- **Phase 11C (Independent review).** Independently reviewed the Phase 11B draft against 12 critical domains (causality, localisation, spatial validation, target-prioritisation logic, null findings, limitations); verified the manuscript body was unmodified by the review itself.
- **Phase 11D (Full assembly).** Assembled the complete manuscript (Title, Abstract, Introduction, Results, Discussion, Limitations, Methods summary, legends), preserving all claim-control constraints (non-causal microbiome language, `PARTIAL_SPATIAL_SUPPORT` for protein secretion, CTCFL/BORIS exclusion).
- **Phase 11E (Language and format review).** Improved readability and UK-spelling consistency; expanded the Discussion and figure/table legends; strictly preserved all epistemic constraints.
- **Phase 11F (Final claim audit).** Independently audited 15 required computational-claim constraints: **100% compliance**. Null findings remained visible.
- **Phase 11G / 11G-R1b / R2 / R3 (Reference and callout audit).** Repaired citation and figure/table callout inconsistencies; confirmed scientific content, evidence categories, and target rankings were unchanged.
- **Post-11G (Workspace cleanup).** Removed temporary agent utility scripts and added intermediate model directories to `.gitignore`.
- **Phase 11H (Submission package assembly).** Assembled manuscript, figures, tables, legends, cover letter, and checklist into the submission package. Decision: **PASS**.
- **Phase 11I-A (Final QA and journal-gap audit).** Verified package completeness, callout consistency, and claim-control adherence; logged 11 journal-specific items as `TO_BE_CONFIRMED`. Final readiness decision: **`READY_FOR_TARGET_JOURNAL_SELECTION`**.

## Open items

- `TO_VERIFY`: official documentation does not explain the literal `NA` entries in the GEO expression matrix; raw-count reprocessing remains a possible later check.
- `TO_VERIFY`: Bailey and full Chan-Seng-Yue 100-gene frameworks remain exploratory (no pre-fitted classifier available for direct reproduction).
- `TO_VERIFY`: exact Kraken2/Bracken database version and normalization formula for the non-integer microbiome abundance scale.
- 11 journal-specific formatting/administrative items (word limits, COI, funding statements, author signatures) pending target-journal confirmation.
- No blocking issues.
