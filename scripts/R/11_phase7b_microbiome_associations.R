#!/usr/bin/env Rscript

# Phase 7B wrapper. The locked execution is implemented in Python because the
# local R environment does not contain vegan or MaAsLin2.

`%||%` <- function(a, b) if (!is.null(a)) a else b
root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(root) != "PDAC") {
  root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "06_scripts/R/11_phase7b_microbiome_associations.R"), "../.."), mustWork = TRUE)
}

tables_dir <- file.path(root, "05_results", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

runtime <- data.frame(
  item = c("R", "vegan_available", "MaAsLin2_available", "wrapper_role"),
  version = c(
    R.version.string,
    as.character(requireNamespace("vegan", quietly = TRUE)),
    as.character(requireNamespace("Maaslin2", quietly = TRUE) || requireNamespace("MaAsLin2", quietly = TRUE)),
    "delegates_to_06_scripts/python/11_summarize_phase7b_associations.py"
  )
)
write.table(runtime, file.path(tables_dir, "phase7b_R_wrapper_runtime.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

script <- file.path(root, "06_scripts", "python", "11_summarize_phase7b_associations.py")
status <- system2("python3", script)
if (!identical(status, 0L)) {
  stop("Phase 7B Python executor failed with status ", status)
}
