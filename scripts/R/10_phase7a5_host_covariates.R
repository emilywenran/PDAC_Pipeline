#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(estimate)
  library(ggplot2)
})

root <- normalizePath(getwd(), mustWork = TRUE)
set.seed(2026)

path <- function(...) file.path(root, ...)

dir.create(path("05_results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("05_results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("04_analysis", "08_host_microbiome_integration"), recursive = TRUE, showWarnings = FALSE)
dir.create(path("03_processed", "expression", "phase7a5_estimate_work"), recursive = TRUE, showWarnings = FALSE)

expr_path <- path("03_processed", "expression", "GSE172356_expression_log2_analysis_ready.tsv.gz")
annot_path <- path("03_processed", "expression", "GSE172356_gene_annotation.tsv")
crosswalk_path <- path("01_metadata", "expression_sample_crosswalk.tsv")
manifest_path <- path("01_metadata", "sample_manifest.tsv")
axis_path <- path("05_results", "tables", "phase5b_sample_continuous_scores.tsv")

stopifnot(file.exists(expr_path), file.exists(annot_path), file.exists(crosswalk_path),
          file.exists(manifest_path), file.exists(axis_path))

expr <- fread(expr_path)
annot <- fread(annot_path)
crosswalk <- fread(crosswalk_path)
manifest <- fread(manifest_path)
axis <- fread(axis_path)

stopifnot("gene" %in% names(expr), nrow(crosswalk) == 62)
sample_cols <- setdiff(names(expr), "gene")
sample_order_matches_crosswalk <- identical(sample_cols, crosswalk$expression_column)
if (!sample_order_matches_crosswalk) {
  stop("Expression sample order does not match expression_sample_crosswalk.tsv")
}
if (anyDuplicated(crosswalk$patient_id) > 0 || anyDuplicated(crosswalk$expression_column) > 0) {
  stop("Duplicated patient or expression sample in crosswalk")
}
if (anyDuplicated(expr$gene) > 0) {
  stop("Duplicated genes in expression matrix")
}

expr_mat <- as.data.frame(expr)
rownames(expr_mat) <- expr_mat$gene
expr_mat$gene <- NULL

estimate_input <- path("03_processed", "expression", "phase7a5_estimate_work", "phase7a5_estimate_input.tsv")
filtered_gct <- path("03_processed", "expression", "phase7a5_estimate_work", "phase7a5_estimate_common_genes.gct")
scores_gct <- path("03_processed", "expression", "phase7a5_estimate_work", "phase7a5_estimate_scores.gct")

write.table(data.frame(GeneSymbol = rownames(expr_mat), expr_mat, check.names = FALSE),
            estimate_input, sep = "\t", quote = FALSE, row.names = FALSE)

warnings <- character()
capture.output({
  tryCatch(
    filterCommonGenes(input.f = estimate_input, output.f = filtered_gct, id = "GeneSymbol"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
})

capture.output({
  tryCatch(
    estimateScore(input.ds = filtered_gct, output.ds = scores_gct, platform = "affymetrix"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
})

scores <- fread(scores_gct, skip = 2)
score_names <- scores[[1]]
score_matrix <- as.data.frame(scores[, -c(1, 2), with = FALSE])
rownames(score_matrix) <- score_names
score_matrix[] <- lapply(score_matrix, as.numeric)

required_score_rows <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
if (!all(required_score_rows %in% rownames(score_matrix))) {
  stop("ESTIMATE score output lacks required rows: ", paste(setdiff(required_score_rows, rownames(score_matrix)), collapse = ", "))
}

covariates <- data.table(
  patient_id = crosswalk$patient_id,
  expression_sample_id = crosswalk$expression_column,
  stromal_score = as.numeric(score_matrix["StromalScore", crosswalk$expression_column]),
  immune_score = as.numeric(score_matrix["ImmuneScore", crosswalk$expression_column]),
  estimate_score = as.numeric(score_matrix["ESTIMATEScore", crosswalk$expression_column]),
  inferred_tumor_purity = as.numeric(score_matrix["TumorPurity", crosswalk$expression_column])
)

common_total <- nrow(common_genes)
common_overlap <- sum(common_genes$GeneSymbol %in% expr$gene)
common_missing <- setdiff(common_genes$GeneSymbol, expr$gene)
si_values <- unique(unlist(SI_geneset[, -1], use.names = FALSE))
si_values <- si_values[!is.na(si_values) & nzchar(si_values)]
si_overlap <- intersect(si_values, expr$gene)
si_missing <- setdiff(si_values, expr$gene)
gene_coverage_text <- sprintf(
  "ESTIMATE_common_genes=%d/%d; SI_genes=%d/%d; missing_common_genes=%d; missing_SI_genes=%d",
  common_overlap, common_total, length(si_overlap), length(si_values),
  length(common_missing), length(si_missing)
)

validation_notes <- c(
  "OFFICIAL_estimate_R_package_1.0.13_from_MDAnderson_R-Forge",
  "input_expression_scale=log2_analysis_ready_values",
  "sample_order_matches_expression_crosswalk=True",
  "subtype_and_microbiome_not_used_for_score_generation",
  "tumor_purity_from_ESTIMATE_affymetrix_cosine_equation"
)
if (length(warnings) > 0) {
  validation_notes <- c(validation_notes, paste0("warnings=", paste(unique(warnings), collapse = " | ")))
} else {
  validation_notes <- c(validation_notes, "warnings=None")
}
if (length(si_missing) > 0) {
  validation_notes <- c(validation_notes, paste0("missing_SI_genes=", paste(si_missing, collapse = ",")))
} else {
  validation_notes <- c(validation_notes, "missing_SI_genes=None")
}

covariates[, `:=`(
  method = "ESTIMATE_R_package_filterCommonGenes_estimateScore",
  method_version = paste0("estimate_", as.character(packageVersion("estimate")),
                          "; R_", paste(R.version$major, R.version$minor, sep = ".")),
  gene_coverage = gene_coverage_text,
  validation_status = "PASS",
  notes = paste(validation_notes, collapse = "; ")
)]

score_cols <- c("stromal_score", "immune_score", "estimate_score", "inferred_tumor_purity")
finite_scores <- all(is.finite(as.matrix(covariates[, ..score_cols])))
if (!finite_scores || nrow(covariates) != 62 || anyDuplicated(covariates$patient_id) > 0) {
  covariates[, validation_status := "FAIL"]
}

fwrite(covariates, path("01_metadata", "host_tme_covariates.tsv"), sep = "\t")

long_scores <- melt(covariates[, c("patient_id", score_cols), with = FALSE],
                    id.vars = "patient_id", variable.name = "score_name", value.name = "value")

extreme_label <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[[2]] - q[[1]]
  lower <- q[[1]] - 3 * iqr
  upper <- q[[2]] + 3 * iqr
  z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  ifelse(x < lower | x > upper | abs(z) > 3.5, "TECHNICAL_EXTREME", "within_distribution")
}

qc_rows <- rbindlist(lapply(score_cols, function(col) {
  x <- covariates[[col]]
  flags <- extreme_label(x)
  data.table(
    metric = col,
    n = sum(!is.na(x)),
    missing_count = sum(is.na(x)),
    infinite_count = sum(is.infinite(x)),
    mean = mean(x),
    sd = sd(x),
    median = median(x),
    min = min(x),
    q1 = as.numeric(quantile(x, 0.25)),
    q3 = as.numeric(quantile(x, 0.75)),
    max = max(x),
    technically_extreme_count = sum(flags == "TECHNICAL_EXTREME"),
    technically_extreme_patients = paste(covariates$patient_id[flags == "TECHNICAL_EXTREME"], collapse = ","),
    validation_status = ifelse(all(is.finite(x)) && sum(is.na(x)) == 0, "PASS", "FAIL"),
    notes = "Descriptive QC only; no microbiome feature relationships inspected."
  )
}))
fwrite(qc_rows, path("05_results", "tables", "phase7a5_host_covariate_qc.tsv"), sep = "\t")

cor_pairs <- combn(score_cols, 2, simplify = FALSE)
cor_rows <- rbindlist(lapply(cor_pairs, function(pair) {
  x <- covariates[[pair[1]]]
  y <- covariates[[pair[2]]]
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.table(
    variable_1 = pair[1],
    variable_2 = pair[2],
    n = sum(complete.cases(x, y)),
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value,
    correlation_warning = ifelse(abs(unname(ct$estimate)) >= 0.9, "SEVERE_COLLINEARITY",
                                 ifelse(abs(unname(ct$estimate)) >= 0.7, "HIGH_CORRELATION", "none")),
    notes = "Descriptive host-covariate correlation only; not tested against microbiome features."
  )
}))
fwrite(cor_rows, path("05_results", "tables", "phase7a5_host_covariate_correlations.tsv"), sep = "\t")

axis_primary <- axis[analysis_id == "AXIS_MOFFITT50_PRIMARY",
                     .(patient_id, expression_sample_id, host_transcriptional_score = basal_classical_contrast)]
model_data <- merge(axis_primary, covariates, by = c("patient_id", "expression_sample_id"))

standardize <- function(x) as.numeric(scale(as.numeric(x)))
vif_values <- function(df) {
  df <- as.data.frame(lapply(df, as.numeric), check.names = FALSE)
  if (ncol(df) <= 1) return(setNames(1, names(df)))
  out <- numeric(ncol(df))
  names(out) <- names(df)
  for (nm in names(df)) {
    others <- setdiff(names(df), nm)
    fit <- lm(df[[nm]] ~ ., data = as.data.frame(df[, others, drop = FALSE]))
    r2 <- summary(fit)$r.squared
    out[[nm]] <- ifelse(isTRUE(all.equal(r2, 1)), Inf, 1 / (1 - r2))
  }
  out
}

condition_number <- function(df) {
  df <- as.data.frame(lapply(df, as.numeric), check.names = FALSE)
  x <- as.matrix(data.frame(intercept = 1, lapply(df, standardize), check.names = FALSE))
  kappa(x, exact = TRUE)
}

model_specs <- list(
  list(id = "Model_0", covars = c("host_transcriptional_score"), role = "primary"),
  list(id = "Model_3P", covars = c("host_transcriptional_score", "inferred_tumor_purity"), role = "sensitivity"),
  list(id = "Model_3I", covars = c("host_transcriptional_score", "immune_score"), role = "sensitivity"),
  list(id = "Model_3S", covars = c("host_transcriptional_score", "stromal_score"), role = "sensitivity"),
  list(id = "Model_3ALL_TME", covars = c("host_transcriptional_score", "inferred_tumor_purity", "immune_score", "stromal_score", "estimate_score"), role = "not_prespecified_combined_screen")
)

feasibility <- rbindlist(lapply(model_specs, function(spec) {
  covars <- spec$covars
  df <- model_data[, covars, with = FALSE]
  cc <- complete.cases(df)
  df_cc <- df[cc]
  cor_warn <- "none"
  if (length(spec$covars) > 1) {
    cm <- suppressWarnings(cor(df_cc, method = "spearman"))
    max_abs <- max(abs(cm[upper.tri(cm)]))
    cor_warn <- ifelse(max_abs >= 0.9, "SEVERE_COLLINEARITY",
                       ifelse(max_abs >= 0.7, "HIGH_CORRELATION", "none"))
  }
  vifs <- vif_values(as.data.frame(df_cc))
  max_vif <- max(vifs, na.rm = TRUE)
  cond <- condition_number(as.data.frame(df_cc))
  available <- nrow(df_cc)
  edf <- available - (length(spec$covars) + 1)
  permitted <- available >= 50 && edf >= 50 && is.finite(max_vif) && max_vif < 5 &&
    is.finite(cond) && cond < 30 && cor_warn != "SEVERE_COLLINEARITY"
  if (spec$id == "Model_3ALL_TME" && cor_warn != "none") {
    permitted <- FALSE
  }
  reason <- if (permitted) {
    "Meets prespecified completeness, VIF, condition-number, and effective-df criteria."
  } else {
    paste(c(
      if (available < 50) "available_patients_below_50" else NULL,
      if (edf < 50) "effective_df_below_50" else NULL,
      if (!is.finite(max_vif) || max_vif >= 5) "maximum_VIF_ge_5" else NULL,
      if (!is.finite(cond) || cond >= 30) "condition_number_ge_30" else NULL,
      if (cor_warn == "SEVERE_COLLINEARITY") "severe_pairwise_correlation_ge_0.9" else NULL,
      if (spec$id == "Model_3ALL_TME" && cor_warn != "none") "combined_TME_model_blocked_by_correlation_warning" else NULL
    ), collapse = "; ")
  }
  data.table(
    model_id = spec$id,
    covariates = paste(spec$covars, collapse = " + "),
    available_patients = available,
    correlation_warning = cor_warn,
    maximum_VIF = max_vif,
    condition_number = cond,
    model_permitted = ifelse(permitted, "YES", "NO"),
    reason = reason,
    analysis_role = spec$role,
    notes = paste0("effective_df=", edf, "; criteria: complete_cases>=50, effective_df>=50, max_VIF<5, condition_number<30, no severe pairwise collinearity; Model 0 remains primary.")
  )
}))
fwrite(feasibility, path("05_results", "tables", "phase7a5_covariate_model_feasibility.tsv"), sep = "\t")

pdf(path("05_results", "figures", "phase7a5_host_covariate_distributions.pdf"), width = 8.5, height = 6)
print(
  ggplot(long_scores, aes(x = value)) +
    geom_histogram(bins = 18, fill = "#52796f", color = "white") +
    facet_wrap(~ score_name, scales = "free", ncol = 2) +
    theme_bw(base_size = 10) +
    labs(x = "Score", y = "Sample count", title = "Phase 7A.5 Host TME Covariate Distributions")
)
dev.off()

cor_mat <- cor(covariates[, ..score_cols], method = "spearman")
cor_long <- as.data.table(as.table(cor_mat))
names(cor_long) <- c("variable_1", "variable_2", "rho")
pdf(path("05_results", "figures", "phase7a5_host_covariate_correlation.pdf"), width = 7, height = 6)
print(
  ggplot(cor_long, aes(variable_1, variable_2, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    scale_fill_gradient2(low = "#4c78a8", mid = "white", high = "#b05a4a", limits = c(-1, 1)) +
    coord_equal() +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = "Spearman rho", title = "Phase 7A.5 Host Covariate Correlations")
)
dev.off()

runtime <- data.table(
  component = c("R", "estimate", "data.table", "ggplot2"),
  version = c(R.version.string, as.character(packageVersion("estimate")),
              as.character(packageVersion("data.table")), as.character(packageVersion("ggplot2"))),
  source = c("local R runtime", "R-Forge official ESTIMATE source package", "R package", "R package")
)
fwrite(runtime, path("05_results", "tables", "phase7a5_runtime_versions.tsv"), sep = "\t")

missing_report <- data.table(
  category = c("expression_samples", "unique_patients", "duplicated_patients", "missing_scores", "infinite_scores",
               "expression_sample_order_matches_crosswalk", "common_gene_missing_count", "SI_gene_missing_count"),
  value = c(nrow(covariates), uniqueN(covariates$patient_id), anyDuplicated(covariates$patient_id),
            sum(is.na(as.matrix(covariates[, ..score_cols]))), sum(is.infinite(as.matrix(covariates[, ..score_cols]))),
            sample_order_matches_crosswalk, length(common_missing), length(si_missing)),
  notes = c("Expected 62", "Expected 62", "Expected 0", "Expected 0", "Expected 0",
            "Expected TRUE", paste(common_missing, collapse = ","), paste(si_missing, collapse = ","))
)
fwrite(missing_report, path("05_results", "tables", "phase7a5_missingness_report.tsv"), sep = "\t")

cat("Phase 7A.5 host covariate calculation complete.\n")
cat("Permitted sensitivity models:\n")
print(feasibility[model_permitted == "YES" & model_id %in% c("Model_3P", "Model_3I", "Model_3S"),
                  .(model_id, covariates, maximum_VIF, condition_number)])
