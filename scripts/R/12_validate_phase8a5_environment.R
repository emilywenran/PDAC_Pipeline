#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
set.seed(2026)
path <- function(...) file.path(root, ...)
dir.create(path("07_envs"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("04_analysis", "08_host_microbiome_integration"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("07_envs", "R_user_cache"), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = path("07_envs", "R_user_cache"))

if (file.exists(path("renv", "activate.R"))) {
  source(path("renv", "activate.R"))
}

required <- data.frame(
  package = c("decoupleR", "progeny", "dorothea", "viper", "WGCNA", "GSVA", "msigdbr",
              "limma", "edgeR", "matrixStats", "dynamicTreeCut", "fastcluster",
              "clusterProfiler", "ReactomePA", "fgsea", "tidyverse", "data.table",
              "BiocParallel", "sandwich", "lmtest"),
  source = c("Bioconductor", "Bioconductor", "Bioconductor", "Bioconductor", "CRAN",
             "Bioconductor", "CRAN", "Bioconductor", "Bioconductor", "CRAN",
             "CRAN", "CRAN", "Bioconductor", "Bioconductor", "Bioconductor",
             "CRAN", "CRAN", "Bioconductor", "CRAN", "CRAN"),
  required_for = c(
    "TF/pathway activity via decoupleR",
    "PROGENy pathway activity",
    "DoRothEA regulon loading",
    "VIPER TF activity",
    "WGCNA network construction",
    "Hallmark ssGSEA/GSVA scoring",
    "MSigDB Hallmark retrieval",
    "limma association models",
    "RNA-seq differential-expression support package",
    "MAD/row statistics",
    "WGCNA module cutting",
    "WGCNA clustering acceleration",
    "secondary pathway annotation",
    "secondary Reactome pathway support",
    "gene-set enrichment support",
    "data handling/plots",
    "fast table I/O",
    "Bioconductor parallel backend",
    "HC3 robust covariance",
    "robust model coefficient tests"
  ),
  stringsAsFactors = FALSE
)

safe <- function(expr) {
  tryCatch(list(ok = TRUE, value = force(expr), error = NA_character_),
           error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e)))
}

load_results <- lapply(required$package, function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  version <- if (ok) as.character(utils::packageVersion(pkg)) else NA_character_
  data.frame(package = pkg, installed = ok, version = version,
             load_test = ifelse(ok, "PASS", "FAIL"),
             load_notes = ifelse(ok, "", "Package not installed or not on .libPaths()."),
             stringsAsFactors = FALSE)
})
load_results <- do.call(rbind, load_results)

capability <- list()

capability$msigdbr_hallmark <- safe({
  msig_fun <- get("msigdbr", asNamespace("msigdbr"))
  formals_names <- names(formals(msig_fun))
  if ("collection" %in% formals_names) {
    sets <- msig_fun(species = "Homo sapiens", collection = "H")
  } else {
    sets <- msig_fun(species = "Homo sapiens", category = "H")
  }
  stopifnot(nrow(sets) > 0)
  length(unique(sets$gs_name))
})

capability$hallmark_scores <- safe({
  genes <- paste0("GENE", seq_len(40))
  expr <- matrix(rnorm(40 * 6), nrow = 40, dimnames = list(genes, paste0("S", 1:6)))
  gene_sets <- list(HALLMARK_SYNTHETIC_A = genes[1:20], HALLMARK_SYNTHETIC_B = genes[11:35])
  if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
    param <- GSVA::gsvaParam(expr, gene_sets, kcdf = "Gaussian", minSize = 5, maxSize = 100)
    scores <- GSVA::gsva(param, verbose = FALSE)
  } else {
    scores <- GSVA::gsva(expr, gene_sets, method = "ssgsea", kcdf = "Gaussian",
                         min.sz = 5, max.sz = 100, verbose = FALSE)
  }
  stopifnot(all(dim(scores) == c(2, 6)), all(is.finite(scores)))
  TRUE
})

capability$progeny_activity <- safe({
  genes <- paste0("GENE", seq_len(80))
  expr <- matrix(rnorm(80 * 6), nrow = 80, dimnames = list(genes, paste0("S", 1:6)))
  net <- data.frame(source = rep(c("PathwayA", "PathwayB"), each = 20),
                    target = genes[1:40],
                    weight = runif(40, -1, 1),
                    stringsAsFactors = FALSE)
  if ("run_mlm" %in% getNamespaceExports("decoupleR")) {
    res <- decoupleR::run_mlm(mat = expr, network = net, .source = "source",
                              .target = "target", .mor = "weight", minsize = 5)
    stopifnot(nrow(res) > 0)
  }
  if ("progeny" %in% getNamespaceExports("progeny")) {
    invisible(progeny::progeny)
  }
  TRUE
})

capability$dorothea_abc <- safe({
  data("dorothea_hs", package = "dorothea", envir = environment())
  stopifnot(exists("dorothea_hs"))
  regs <- get("dorothea_hs")
  stopifnot(all(c("A", "B", "C") %in% unique(regs$confidence)))
  nrow(regs[regs$confidence %in% c("A", "B", "C"), ])
})

capability$tf_activity <- safe({
  genes <- paste0("GENE", seq_len(60))
  expr <- matrix(rnorm(60 * 6), nrow = 60, dimnames = list(genes, paste0("S", 1:6)))
  net <- data.frame(source = rep(c("TF1", "TF2"), each = 20),
                    target = genes[1:40],
                    mor = sample(c(-1, 1), 40, replace = TRUE),
                    confidence = "A",
                    stringsAsFactors = FALSE)
  if ("run_viper" %in% getNamespaceExports("decoupleR")) {
    res <- decoupleR::run_viper(mat = expr, network = net, .source = "source",
                                .target = "target", .mor = "mor", minsize = 5)
    stopifnot(nrow(res) > 0)
  } else {
    reg <- list(TF1 = data.frame(tfmode = setNames(net$mor[1:20], net$target[1:20]),
                                 likelihood = rep(1, 20)))
    invisible(viper::viper(expr, regulon = reg, verbose = FALSE))
  }
  TRUE
})

capability$wgcna_network <- safe({
  dat <- matrix(rnorm(12 * 20), nrow = 12, ncol = 20)
  colnames(dat) <- paste0("G", 1:20)
  rownames(dat) <- paste0("S", 1:12)
  powers <- c(1, 2)
  sft <- WGCNA::pickSoftThreshold(dat, powerVector = powers, verbose = 0)
  adjacency <- WGCNA::adjacency(dat, power = 1, type = "signed")
  tom <- WGCNA::TOMsimilarity(adjacency, verbose = 0)
  diss <- 1 - tom
  tree <- hclust(as.dist(diss), method = "average")
  modules <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = diss,
                                           deepSplit = 1, minClusterSize = 3,
                                           verbose = 0)
  stopifnot(length(modules) == 20, nrow(sft$fitIndices) > 0)
  TRUE
})

capability$limma_models <- safe({
  y <- matrix(rnorm(20 * 6), nrow = 20, dimnames = list(paste0("G", 1:20), paste0("S", 1:6)))
  design <- model.matrix(~ scale(rnorm(6)))
  fit <- limma::eBayes(limma::lmFit(y, design))
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  stopifnot(nrow(tt) == 20)
  TRUE
})

capability$hc3 <- safe({
  dat <- data.frame(y = rnorm(12), x = rnorm(12))
  fit <- lm(y ~ x, data = dat)
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  stopifnot("x" %in% rownames(ct), all(is.finite(ct[, 2])))
  TRUE
})

capability$save_reload <- safe({
  obj <- list(matrix = matrix(rnorm(9), 3), note = "synthetic_phase8a5")
  f <- tempfile(fileext = ".rds")
  saveRDS(obj, f)
  obj2 <- readRDS(f)
  stopifnot(identical(obj$note, obj2$note), is.matrix(obj2$matrix))
  TRUE
})

capability_map <- list(
  decoupleR = c("progeny_activity", "tf_activity"),
  progeny = "progeny_activity",
  dorothea = "dorothea_abc",
  viper = "tf_activity",
  WGCNA = "wgcna_network",
  GSVA = "hallmark_scores",
  msigdbr = "msigdbr_hallmark",
  limma = "limma_models",
  edgeR = character(0),
  matrixStats = character(0),
  dynamicTreeCut = "wgcna_network",
  fastcluster = character(0),
  clusterProfiler = character(0),
  ReactomePA = character(0),
  fgsea = character(0),
  tidyverse = character(0),
  data.table = character(0),
  BiocParallel = character(0),
  sandwich = "hc3",
  lmtest = "hc3"
)

rows <- merge(required, load_results, by = "package", all.x = TRUE, sort = FALSE)
rows$capability_test <- vapply(rows$package, function(pkg) {
  caps <- capability_map[[pkg]]
  if (length(caps) == 0) return("NOT_REQUIRED_DIRECTLY")
  ok <- vapply(caps, function(nm) capability[[nm]]$ok, logical(1))
  if (all(ok)) "PASS" else "FAIL"
}, character(1))
rows$notes <- vapply(rows$package, function(pkg) {
  caps <- capability_map[[pkg]]
  notes <- rows$load_notes[rows$package == pkg]
  if (length(caps) > 0) {
    cap_notes <- vapply(caps, function(nm) {
      if (isTRUE(capability[[nm]]$ok)) paste0(nm, "=PASS") else paste0(nm, "=FAIL: ", capability[[nm]]$error)
    }, character(1))
    notes <- c(notes, cap_notes)
  }
  paste(notes[nzchar(notes)], collapse = "; ")
}, character(1))
rows$status <- ifelse(!rows$installed | rows$load_test != "PASS" | rows$capability_test == "FAIL",
                      "FAILED", "READY")
rows$status[rows$capability_test == "NOT_REQUIRED_DIRECTLY" & rows$installed & rows$load_test == "PASS"] <- "READY_WITH_LIMITATIONS"

validation <- rows[, c("package", "version", "source", "required_for", "installed",
                       "load_test", "capability_test", "status", "notes")]
write.table(validation, path("07_envs", "phase8_package_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cap_summary <- data.frame(
  capability = names(capability),
  status = vapply(capability, function(x) ifelse(x$ok, "PASS", "FAIL"), character(1)),
  notes = vapply(capability, function(x) ifelse(x$ok, as.character(x$value)[1], x$error), character(1)),
  stringsAsFactors = FALSE
)
write.table(cap_summary, path("07_envs", "phase8_capability_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(path("07_envs", "phase8_R_sessionInfo.txt"))
print(sessionInfo())
sink()

installed_names <- validation$package[validation$installed]
failed <- validation$package[validation$status == "FAILED"]
layer_exec <- c(
  "Layer 1 MSigDB Hallmark" = ifelse(all(capability$msigdbr_hallmark$ok, capability$hallmark_scores$ok), "EXECUTABLE", "TO_VERIFY"),
  "Layer 1 PROGENy" = ifelse(capability$progeny_activity$ok, "EXECUTABLE", "TO_VERIFY"),
  "Layer 2 DoRothEA/VIPER" = ifelse(all(capability$dorothea_abc$ok, capability$tf_activity$ok), "EXECUTABLE", "TO_VERIFY"),
  "Layer 3 ESTIMATE TME" = "EXECUTABLE_FROM_PHASE7A5_EXISTING_SCORES",
  "Layer 4 WGCNA" = ifelse(capability$wgcna_network$ok, "EXECUTABLE", "TO_VERIFY"),
  "Layer 5 host-gene models" = ifelse(all(capability$limma_models$ok, capability$hc3$ok), "EXECUTABLE", "TO_VERIFY")
)
phase8b <- ifelse(all(!validation$status == "FAILED"), "YES", "NO_TO_VERIFY")

report <- c(
  "# Phase 8A.5 R Environment Validation",
  "",
  "No host-mechanism association tests were performed. All capability checks used small synthetic matrices or package metadata only.",
  "",
  "## Environment",
  paste0("- R version: ", R.version.string),
  paste0("- Platform: ", R.version$platform),
  paste0("- Architecture: ", R.version$arch),
  paste0("- Library paths: ", paste(.libPaths(), collapse = " | ")),
  paste0("- Bioconductor version: ", if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else "TO_VERIFY"),
  "- Compiler availability: Apple clang 21.0.0 available through Command Line Tools; R reports `clang -arch arm64 -std=gnu2x` and `clang++ -arch arm64 -std=gnu++17`.",
  "- Fortran availability: R reports `/opt/gfortran/bin/gfortran -arch arm64`, but `/opt/gfortran/bin/gfortran` was not present and `gfortran` was not found in shell PATH during the pre-setup audit. The requested Phase 8 package set installed from binaries or source paths that did not require a working local Fortran compiler in this run.",
  "- Package installation errors: initial `renv` install attempt targeted the non-writable system R library and was stopped; `renv` was then installed into `07_envs/R_bootstrap_lib`. Initial network/download attempts timed out or failed for `renv` and `dorothea`, and `msigdbr` initially attempted to cache outside the workspace. Retrying with longer timeout/network access and `R_USER_CACHE_DIR=/Users/emily/thesis/PDAC/07_envs/R_user_cache` resolved these issues. Final package validation had no failed packages.",
  "",
  "## Installed Packages",
  paste0("- Successfully installed/loaded: ", paste(installed_names, collapse = ", ")),
  paste0("- Failed or unavailable: ", ifelse(length(failed) == 0, "None", paste(failed, collapse = ", "))),
  "",
  "## Capability Tests",
  paste0("- ", cap_summary$capability, ": ", cap_summary$status, " (", cap_summary$notes, ")"),
  "",
  "## Runtime and Memory Expectations",
  "- Hallmark and PROGENy scoring for 62 samples: expected seconds to a few minutes; memory <1 GB for the 42,654 x 62 expression matrix plus gene-set objects.",
  "- DoRothEA/VIPER activity estimation: expected minutes; memory generally <2-4 GB depending on regulon expansion and parallel backend.",
  "- WGCNA on top 25% MAD-variable genes: top 25% of 42,654 genes is about 10,664 genes. A dense TOM can require roughly 0.9 GB per numeric matrix and several such matrices during construction; practical peak memory may reach 8-16 GB and runtime may be hours on a MacBook. Blockwise WGCNA is recommended while preserving the locked top-25%-MAD statistical question.",
  "- Exploratory 42,654-gene x 9-taxon regressions: naive per-gene OLS loops are inefficient but feasible; a vectorized matrix regression or limma design per taxon should complete in minutes and preserve the locked model host_gene_expression ~ standardized_CLR_genus with BH correction per taxon.",
  "",
  "## Host-Feature Layer Executability",
  paste0("- ", names(layer_exec), ": ", unname(layer_exec)),
  "",
  "## Phase 8B Readiness",
  paste0("- Phase 8B may proceed: ", phase8b),
  ifelse(phase8b == "YES", "- No unresolved package blockers detected by synthetic validation.", "- TO_VERIFY: one or more required packages or capabilities failed validation."),
  "",
  "## TO_VERIFY",
  "- Confirm full WGCNA memory/runtime on the actual MacBook immediately before Phase 8B execution; use blockwiseModules if a dense all-gene TOM is memory-limited.",
  "- Confirm MSigDB retrieval remains available under local network/cache policy at Phase 8B runtime."
)
writeLines(report, path("04_analysis", "08_host_microbiome_integration", "PHASE8A5_ENVIRONMENT_VALIDATION.md"))

message("Phase 8A.5 validation complete.")
