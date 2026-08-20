#!/usr/bin/env Rscript

Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "false")
root <- "/Users/emily/thesis/PDAC"
status <- system2("python3", file.path(root, "06_scripts/python/16_phase9b3b_spatial_validation.py"), stdout = TRUE, stderr = TRUE)
cat(paste(status, collapse = "\n"), "\n")
exit_status <- attr(status, "status")
if (is.null(exit_status)) exit_status <- 0
quit(status = exit_status)
