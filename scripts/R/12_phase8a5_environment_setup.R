#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
set.seed(2026)

path <- function(...) file.path(root, ...)
dir.create(path("07_envs"), recursive = TRUE, showWarnings = FALSE)

bootstrap_lib <- path("07_envs", "R_bootstrap_lib")
dir.create(bootstrap_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(bootstrap_lib, .libPaths()))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 600,
  Ncpus = max(1L, parallel::detectCores(logical = TRUE) - 1L)
)

inspect_environment <- function() {
  compiler_checks <- list(
    cc = tryCatch(system2("R", c("CMD", "config", "CC"), stdout = TRUE, stderr = TRUE), error = conditionMessage),
    cxx = tryCatch(system2("R", c("CMD", "config", "CXX"), stdout = TRUE, stderr = TRUE), error = conditionMessage),
    fc = tryCatch(system2("R", c("CMD", "config", "FC"), stdout = TRUE, stderr = TRUE), error = conditionMessage),
    clang = tryCatch(system2("clang", "--version", stdout = TRUE, stderr = TRUE), error = conditionMessage),
    gfortran_path = Sys.which("gfortran")
  )
  ip <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath", "Priority")])
  ip <- ip[order(tolower(ip$Package)), ]
  write.table(ip, path("07_envs", "phase8_pre_setup_installed_packages.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  sink(path("07_envs", "phase8_pre_setup_R_environment.txt"))
  cat("Phase 8A.5 pre-setup R environment\n")
  cat("timestamp\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
  cat("R_version\t", R.version.string, "\n", sep = "")
  cat("platform\t", R.version$platform, "\n", sep = "")
  cat("arch\t", R.version$arch, "\n", sep = "")
  cat("library_paths\t", paste(.libPaths(), collapse = " | "), "\n", sep = "")
  cat("Bioconductor_version\t",
      if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else "BiocManager_not_installed",
      "\n", sep = "")
  cat("compiler_CC\t", paste(compiler_checks$cc, collapse = " "), "\n", sep = "")
  cat("compiler_CXX\t", paste(compiler_checks$cxx, collapse = " "), "\n", sep = "")
  cat("compiler_FC\t", paste(compiler_checks$fc, collapse = " "), "\n", sep = "")
  cat("shell_gfortran\t", ifelse(nzchar(compiler_checks$gfortran_path), compiler_checks$gfortran_path, "NOT_FOUND_IN_PATH"), "\n", sep = "")
  cat("clang_version\t", paste(compiler_checks$clang, collapse = " | "), "\n", sep = "")
  sink()
}

inspect_environment()

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", lib = bootstrap_lib, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(renv))

if (!file.exists(path("renv", "activate.R"))) {
  renv::init(project = root, bare = TRUE, restart = FALSE)
} else {
  renv::activate(project = root)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  renv::install("BiocManager")
}

cran_packages <- c(
  "renv",
  "BiocManager",
  "WGCNA",
  "msigdbr",
  "matrixStats",
  "dynamicTreeCut",
  "fastcluster",
  "tidyverse",
  "data.table",
  "sandwich",
  "lmtest"
)

bioc_packages <- c(
  "decoupleR",
  "progeny",
  "dorothea",
  "viper",
  "GSVA",
  "limma",
  "edgeR",
  "clusterProfiler",
  "ReactomePA",
  "fgsea",
  "BiocParallel"
)

package_plan <- data.frame(
  package = c(cran_packages, bioc_packages),
  source = c(rep("CRAN", length(cran_packages)), rep("Bioconductor", length(bioc_packages))),
  stringsAsFactors = FALSE
)
write.table(package_plan, path("07_envs", "phase8_required_package_plan.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

install_log <- lapply(seq_len(nrow(package_plan)), function(i) {
  pkg <- package_plan$package[[i]]
  src <- package_plan$source[[i]]
  started <- Sys.time()
  message("Installing/checking ", pkg, " from ", src)
  result <- tryCatch({
    if (src == "Bioconductor") {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      renv::install(pkg)
    }
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  data.frame(
    package = pkg,
    source = src,
    started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    result = result,
    stringsAsFactors = FALSE
  )
})
install_log <- do.call(rbind, install_log)
write.table(install_log, path("07_envs", "phase8_install_log.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

renv::snapshot(project = root, prompt = FALSE)

ip <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")])
project_lib <- normalizePath(.libPaths()[[1]], mustWork = FALSE)
ip$in_project_lib <- normalizePath(ip$LibPath, mustWork = FALSE) == project_lib
ip <- ip[order(ip$Package, !ip$in_project_lib), ]
ip <- ip[!duplicated(ip$Package), ]
inventory <- merge(package_plan, ip[, c("Package", "Version", "LibPath")],
                   by.x = "package", by.y = "Package", all.x = TRUE, sort = FALSE)
inventory$installed <- !is.na(inventory$Version)
names(inventory)[names(inventory) == "Version"] <- "version"
names(inventory)[names(inventory) == "LibPath"] <- "library_path"
write.table(inventory, path("07_envs", "phase8_r_environment.yml"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

sink(path("07_envs", "phase8_R_sessionInfo.txt"))
print(sessionInfo())
sink()

message("Phase 8A.5 environment setup complete. Run 06_scripts/R/12_validate_phase8a5_environment.R next.")
