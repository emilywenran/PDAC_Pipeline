#!/usr/bin/env Rscript
# Reproducible Phase 6C entry point. The implementation is in Python so the
# same code path generates matrices, QC tables, figures, and validations.

message("Phase 6C microbiome preprocessing is implemented in:")
message("  06_scripts/python/09_phase6c_prepare_microbiome.py")
status <- system2("python3", c("06_scripts/python/09_phase6c_prepare_microbiome.py"))
if (status != 0) {
  stop("Phase 6C Python preprocessing failed", call. = FALSE)
}
