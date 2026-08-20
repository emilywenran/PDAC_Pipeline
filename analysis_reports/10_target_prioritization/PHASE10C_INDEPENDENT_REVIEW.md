# Phase 10C: Independent Review of Phase 10B Cross-Layer Synthesis

## Overview
This document represents the independent Phase 10C review of the Phase 10B cross-layer evidence synthesis and target prioritization. The objective is to ensure strict adherence to the locked Phase 10A framework and ensure reproducibility.

## Audit Findings

### 1. Phase 10A Adherence
- **Status**: **FAIL**
- **Details**: Phase 10A was properly locked and committed. However, Phase 10B did not strictly follow the objective thresholds defined in `phase10a_target_prioritization_framework.tsv`. Values were populated manually via hardcoding in `18_phase10b_cross_layer_synthesis.py` rather than programmatically evaluating against the predefined thresholds (e.g. `FDR < 0.05`, `Score < -0.5`). 
- **Required Correction**: The synthesis script must programmatically fetch, score, and evaluate targets against the exact locked thresholds, or the data must be stored and objectively evaluated.

### 2. Selection of CTCFL/BORIS
- **Status**: **FAIL**
- **Details**: CTCFL was selected based on post-hoc biological preference and manual hardcoding. Furthermore, CTCFL's `sc_evidence` is `CELL_COMPOSITION_EXPLAINED`. According to the Phase 10A prioritization framework, `cell_type_specificity` has a `High` weight ("Expression restricted to malignant or specific stromal cells"). Since its single-cell association is composition-sensitive, promoting it without penalizing this criterion violates the framework.
- **Required Correction**: Apply the predefined scoring rules objectively without manual circumvention.

### 3. Verification of External Database Claims
- **Status**: **FAIL**
- **Details**: The claims regarding OpenTargets, GTEx, and ChEMBL were not backed by saved query outputs in the project repository (outputs were temporarily stored in `/tmp/`), and the final execution script `18_phase10b_cross_layer_synthesis.py` contained no reproducible query commands.
- **Required Correction**: External database queries must be reproducible. The script must either invoke the API wrappers or load saved JSON outputs committed to `02_data/external/` or `05_results/`.

### 4. Deprioritization of HALLMARK_PROTEIN_SECRETION and BHLHE40
- **Status**: **FAIL**
- **Details**: These candidates were deprioritized using qualitative statements ("Pan-expressed", "High expression in normal tissues") rather than the objective, quantitative thresholds locked in Phase 10A.
- **Required Correction**: Re-evaluate these targets using the formal Phase 10A framework thresholds.

### 5. Over-Promotion of Partial/Unsupported Evidence
- **Status**: **FAIL**
- **Details**: `HALLMARK_SPERMATOGENESIS` was classified as `PARTIALLY_REPLICATED` but omitted from evaluation (filled with "N/A").
- **Required Correction**: All partially replicated candidates must be formally scored through the prioritization framework.

## Final Decision
**Decision: FAIL_REQUIRES_REANALYSIS**

Phase 10B must be re-executed using reproducible scripts that genuinely query the required databases, save their outputs to the repository, and evaluate all eligible targets objectively against the predefined Phase 10A thresholds.
