# Phase 8A.5 R Environment Validation

No host-mechanism association tests were performed. All capability checks used small synthetic matrices or package metadata only.

## Environment
- R version: R version 4.5.3 (2026-03-11)
- Platform: aarch64-apple-darwin20
- Architecture: aarch64
- Library paths: ~/thesis/PDAC/renv/library/macos/R-4.5/aarch64-apple-darwin20 | /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
- Bioconductor version: 3.22
- Compiler availability: Apple clang 21.0.0 available through Command Line Tools; R reports `clang -arch arm64 -std=gnu2x` and `clang++ -arch arm64 -std=gnu++17`.
- Fortran availability: R reports `/opt/gfortran/bin/gfortran -arch arm64`, but `/opt/gfortran/bin/gfortran` was not present and `gfortran` was not found in shell PATH during the pre-setup audit. The requested Phase 8 package set installed from binaries or source paths that did not require a working local Fortran compiler in this run.
- Package installation errors: initial `renv` install attempt targeted the non-writable system R library and was stopped; `renv` was then installed into `07_envs/R_bootstrap_lib`. Initial network/download attempts timed out or failed for `renv` and `dorothea`, and `msigdbr` initially attempted to cache outside the workspace. Retrying with longer timeout/network access and `R_USER_CACHE_DIR=~/thesis/PDAC/07_envs/R_user_cache` resolved these issues. Final package validation had no failed packages.

## Installed Packages
- Successfully installed/loaded: decoupleR, progeny, dorothea, viper, WGCNA, GSVA, msigdbr, limma, edgeR, matrixStats, dynamicTreeCut, fastcluster, clusterProfiler, ReactomePA, fgsea, tidyverse, data.table, BiocParallel, sandwich, lmtest
- Failed or unavailable: None

## Capability Tests
- msigdbr_hallmark: PASS (50)
- hallmark_scores: PASS (TRUE)
- progeny_activity: PASS (TRUE)
- dorothea_abc: PASS (13223)
- tf_activity: PASS (TRUE)
- wgcna_network: PASS (TRUE)
- limma_models: PASS (TRUE)
- hc3: PASS (TRUE)
- save_reload: PASS (TRUE)

## Runtime and Memory Expectations
- Hallmark and PROGENy scoring for 62 samples: expected seconds to a few minutes; memory <1 GB for the 42,654 x 62 expression matrix plus gene-set objects.
- DoRothEA/VIPER activity estimation: expected minutes; memory generally <2-4 GB depending on regulon expansion and parallel backend.
- WGCNA on top 25% MAD-variable genes: top 25% of 42,654 genes is about 10,664 genes. A dense TOM can require roughly 0.9 GB per numeric matrix and several such matrices during construction; practical peak memory may reach 8-16 GB and runtime may be hours on a MacBook. Blockwise WGCNA is recommended while preserving the locked top-25%-MAD statistical question.
- Exploratory 42,654-gene x 9-taxon regressions: naive per-gene OLS loops are inefficient but feasible; a vectorized matrix regression or limma design per taxon should complete in minutes and preserve the locked model host_gene_expression ~ standardized_CLR_genus with BH correction per taxon.

## Host-Feature Layer Executability
- Layer 1 MSigDB Hallmark: EXECUTABLE
- Layer 1 PROGENy: EXECUTABLE
- Layer 2 DoRothEA/VIPER: EXECUTABLE
- Layer 3 ESTIMATE TME: EXECUTABLE_FROM_PHASE7A5_EXISTING_SCORES
- Layer 4 WGCNA: EXECUTABLE
- Layer 5 host-gene models: EXECUTABLE

## Phase 8B Readiness
- Phase 8B may proceed: YES
- No unresolved package blockers detected by synthetic validation.

## TO_VERIFY
- Confirm full WGCNA memory/runtime on the actual MacBook immediately before Phase 8B execution; use blockwiseModules if a dense all-gene TOM is memory-limited.
- Confirm MSigDB retrieval remains available under local network/cache policy at Phase 8B runtime.
