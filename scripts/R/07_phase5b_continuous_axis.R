#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

# Locked Phase 5B entry point.
# The local R library lacks diptest, clinfun, and singscore, so the
# reproducible executor is implemented in Python and called from here.

args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(args_file)) sub("^--file=", "", args_file[[1]]) else "06_scripts/R/07_phase5b_continuous_axis.R"
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
setwd(root)

pkg_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(utils::packageVersion(pkg))
  } else {
    NA_character_
  }
}

versions <- data.frame(
  r_version = R.version.string,
  ggplot2 = pkg_version("ggplot2"),
  boot = pkg_version("boot"),
  diptest = pkg_version("diptest"),
  clinfun = pkg_version("clinfun"),
  singscore = pkg_version("singscore"),
  random_seed = 2026,
  stringsAsFactors = FALSE
)

dir.create("05_results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  versions,
  file = "05_results/tables/phase5b_r_runtime_versions.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

status <- system2("python3", c("06_scripts/python/07_summarize_phase5b_axis.py"))
if (!identical(status, 0L)) {
  stop("Phase 5B Python executor failed with status ", status)
}

status <- system2("python3", c("06_scripts/python/07_validate_phase5b_axis.py"))
if (!identical(status, 0L)) {
  stop("Phase 5B validation failed with status ", status)
}

cat("Phase 5B execution and validation completed.\n")
