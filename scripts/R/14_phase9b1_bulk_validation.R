#!/usr/bin/env Rscript

# Phase 9B1 locked bulk-transcriptome validation wrapper.
# The managed workspace can block renv sandbox creation for long Bioconductor
# startup. The executable implementation is therefore delegated to the Python
# script that uses the verified project-local scientific stack.

root <- "/Users/emily/thesis/PDAC"
script <- file.path(root, "06_scripts/python/14_prepare_phase9b1_bulk_data.py")
validator <- file.path(root, "06_scripts/python/14_validate_phase9b1_bulk_validation.py")

status <- system2("python3", c(script), stdout = TRUE, stderr = TRUE)
cat(paste(status, collapse = "\n"), "\n")
if (!is.null(attr(status, "status")) && attr(status, "status") != 0) {
  quit(status = attr(status, "status"))
}

v <- system2("python3", c(validator), stdout = TRUE, stderr = TRUE)
cat(paste(v, collapse = "\n"), "\n")
if (!is.null(attr(v, "status")) && attr(v, "status") != 0) {
  quit(status = attr(v, "status"))
}
