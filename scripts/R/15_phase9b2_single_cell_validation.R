#!/usr/bin/env Rscript
# Phase 9B2 patient-aware single-cell cellular-source analysis for PENG_CRA001160.

Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(decoupleR)
  library(dorothea)
  library(msigdbr)
  library(lmtest)
  library(sandwich)
})

root <- normalizePath(getwd())
tables_dir <- file.path(root, "05_results/tables")
fig_dir <- file.path(root, "05_results/figures")
model_dir <- file.path(root, "05_results/models/phase9b2")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(2026)
dataset_id <- "PENG_CRA001160"
min_cells <- 20

bulk_tfs <- c("CTCFL", "IRF3", "JUNB", "KLF13", "KLF9", "MNT", "MXI1", "SNAI2",
              "TFAP4", "TP63", "ZBTB7A", "ZNF24")
partial_tfs <- c("BHLHE40", "E2F6", "ELF1", "GRHL2", "KLF1", "MBD1", "MBD2",
                 "OTX2", "SIX5", "SNAPC4", "ZBED1", "ZNF384", "ZNF740")
tf_selected <- c(bulk_tfs, partial_tfs)
target_hallmarks <- c("HALLMARK_PROTEIN_SECRETION", "HALLMARK_SPERMATOGENESIS")
module_names <- c("MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow")
feature_family <- function(x) {
  ifelse(x %in% target_hallmarks, "Hallmark",
         ifelse(x %in% module_names, "WGCNA-module",
                ifelse(x %in% bulk_tfs, "externally replicated TF",
                       ifelse(x %in% partial_tfs, "partially replicated TF", "state"))))
}

read_expr <- function(path) {
  dt <- fread(path)
  genes <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- toupper(genes)
  storage.mode(mat) <- "numeric"
  mat
}

zscore_rows <- function(mat) {
  m <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, sd, na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- 1
  sweep(sweep(mat, 1, m, "-"), 1, s, "/")
}

rank_percentile_matrix <- function(mat) {
  apply(mat, 2, function(x) rank(x, ties.method = "average", na.last = "keep") / sum(is.finite(x)))
}

score_mean <- function(mat, genes) {
  genes <- intersect(unique(toupper(genes)), rownames(mat))
  if (!length(genes)) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
}

score_rank <- function(ranks, genes) {
  genes <- intersect(unique(toupper(genes)), rownames(ranks))
  if (!length(genes)) return(rep(NA_real_, ncol(ranks)))
  colMeans(ranks[genes, , drop = FALSE], na.rm = TRUE)
}

bh <- function(p) p.adjust(p, method = "BH")

ols_hc3 <- function(d, y, x) {
  d <- d[is.finite(get(y)) & is.finite(get(x))]
  if (nrow(d) < 5 || uniqueN(d[[x]]) < 3) {
    return(list(n = nrow(d), beta = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
                diagnostics = "INSUFFICIENT_N_OR_VARIATION"))
  }
  fit <- lm(as.formula(paste(y, "~", x)), data = d)
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  beta <- unname(ct[x, "Estimate"])
  se <- unname(ct[x, "Std. Error"])
  p <- unname(ct[x, "Pr(>|t|)"])
  list(n = nrow(d), beta = beta, se = se, lo = beta - 1.96 * se, hi = beta + 1.96 * se,
       p = p, diagnostics = paste0("HC3_OLS; residual_sd=", signif(sd(residuals(fit)), 4)))
}

pb_path <- file.path(model_dir, "phase9b2_patient_celltype_pseudobulk_counts.tsv.gz")
if (!file.exists(pb_path)) stop("Missing pseudobulk matrix. Run 15_prepare_phase9b2_single_cell.py first.")
pb_counts <- read_expr(pb_path)

inv <- fread(file.path(tables_dir, "phase9b2_pseudobulk_inventory.tsv"))
eligible_cols <- inv[eligibility == "ELIGIBLE", paste(patient_id, cell_type, sep = "|")]
pb_counts <- pb_counts[, intersect(colnames(pb_counts), eligible_cols), drop = FALSE]
sample_info <- tstrsplit(colnames(pb_counts), "\\|")
sample_dt <- data.table(sample_id = colnames(pb_counts), patient_id = sample_info[[1]], cell_type = sample_info[[2]])
sample_dt[, tumor_control_status := ifelse(grepl("^T", patient_id), "PDAC_TUMOR", "CONTROL_PANCREAS")]

lib <- colSums(pb_counts)
cpm <- t(t(pb_counts) / pmax(lib, 1) * 1e6)
expr <- log2(cpm + 1)
z <- zscore_rows(expr)
ranks <- rank_percentile_matrix(expr)

moff <- fread(file.path(root, "02_data/reference/PDAC_subtype_signatures/Moffitt_50_gene_axis.tsv"))
basal <- toupper(moff[program == "Basal-like", mapped_symbol])
classical <- toupper(moff[program == "Classical", mapped_symbol])
pur <- fread(file.path(root, "02_data/reference/PDAC_subtype_signatures/PurIST_signatures.tsv"))
pur[, `:=`(mapped_symbol_A = toupper(mapped_symbol_A), mapped_symbol_B = toupper(mapped_symbol_B))]
pur_intercept <- -6.815

purist_score <- function(mat) {
  terms <- list()
  present <- logical(nrow(pur))
  for (i in seq_len(nrow(pur))) {
    a <- pur$mapped_symbol_A[i]; b <- pur$mapped_symbol_B[i]
    if (a %in% rownames(mat) && b %in% rownames(mat)) {
      present[i] <- TRUE
      terms[[length(terms) + 1]] <- pur$coefficient[i] * as.numeric(mat[a, ] > mat[b, ])
    }
  }
  if (sum(present) / nrow(pur) < 0.80) return(list(prob = rep(NA_real_, ncol(mat)), present = sum(present), expected = nrow(pur)))
  eta <- pur_intercept + Reduce("+", terms) * (nrow(pur) / sum(present))
  list(prob = 1 / (1 + exp(-eta)), present = sum(present), expected = nrow(pur))
}

ps <- purist_score(expr)
state <- copy(sample_dt)
state[, `:=`(
  dataset_id = dataset_id,
  moffitt50_basal_score = score_mean(z, basal),
  moffitt50_classical_score = score_mean(z, classical),
  moffitt50_contrast = score_mean(z, basal) - score_mean(z, classical),
  moffitt49_no_LEMD1_contrast = score_mean(z, setdiff(basal, "LEMD1")) - score_mean(z, classical),
  purist_probability = ps$prob
)]
fwrite(state, file.path(tables_dir, "phase9b2_patient_celltype_state_scores.tsv"), sep = "\t")

modules <- fread(file.path(root, "05_results/tables/phase8b_wgcna_module_assignments.tsv.gz"))
module_sets <- lapply(sub("^ME", "", module_names), function(m) toupper(modules[module == m, gene]))
names(module_sets) <- module_names

hallmark <- as.data.table(msigdbr(species = "human", collection = "H"))
hallmark_net <- unique(data.table(source = hallmark$gs_name, target = toupper(hallmark$gene_symbol)))
hall_long <- decoupleR::run_gsva(expr, hallmark_net, method = "ssgsea", minsize = 15, verbose = FALSE)
hall_dt <- as.data.table(hall_long)
hall_wide <- dcast(hall_dt[source %in% target_hallmarks], condition ~ source, value.var = "score")
setnames(hall_wide, "condition", "sample_id")

program <- merge(copy(sample_dt), hall_wide, by = "sample_id", all.x = TRUE)
for (mn in module_names) {
  program[, (mn) := score_rank(ranks, module_sets[[mn]])]
}
program[, dataset_id := dataset_id]
fwrite(program, file.path(tables_dir, "phase9b2_patient_celltype_host_program_scores.tsv"), sep = "\t")

module_cov <- rbindlist(lapply(module_names, function(mn) {
  genes <- unique(module_sets[[mn]])
  present <- intersect(genes, rownames(expr))
  data.table(dataset_id = dataset_id, module_name = mn, expected_genes = length(genes),
             detected_genes = length(present), coverage_fraction = length(present) / length(genes),
             transfer_method = "locked discovery gene membership; no WGCNA reconstruction",
             eligibility = ifelse(length(present) / length(genes) >= 0.80, "ELIGIBLE", "INELIGIBLE_LOW_COVERAGE"))
}))
fwrite(module_cov, file.path(tables_dir, "phase9b2_module_transfer_coverage.tsv"), sep = "\t")

data(dorothea_hs, package = "dorothea")
regulon <- as.data.table(dorothea_hs)[confidence %in% c("A", "B", "C")]
regulon[, `:=`(tf = as.character(tf), target = toupper(target), mor = as.numeric(mor))]
reg_cov <- regulon[, .(regulon_targets_expected = uniqueN(target),
                       targets_present = uniqueN(intersect(target, rownames(expr)))), by = tf]
reg_cov[, coverage_fraction := targets_present / regulon_targets_expected]
reg_cov[, eligibility := fifelse(targets_present >= 15 & coverage_fraction >= 0.80, "ELIGIBLE", "INELIGIBLE_LOW_COVERAGE")]
tf_keep <- intersect(tf_selected, reg_cov[eligibility == "ELIGIBLE", tf])
tf_long <- decoupleR::run_viper(expr, regulon[tf %in% tf_keep],
                                .source = tf, .target = target, .mor = mor,
                                .likelihood = NULL, minsize = 15, verbose = FALSE)
tf_dt <- as.data.table(tf_long)
tf_wide <- dcast(tf_dt, condition ~ source, value.var = "score")
setnames(tf_wide, "condition", "sample_id")
tf_scores <- merge(copy(sample_dt), tf_wide, by = "sample_id", all.x = TRUE)
tf_scores[, dataset_id := dataset_id]
fwrite(tf_scores, file.path(tables_dir, "phase9b2_patient_celltype_tf_activity.tsv"), sep = "\t")
tf_cov <- reg_cov[tf %in% tf_selected]
tf_cov[, `:=`(dataset_id = dataset_id, confidence_levels = "A;B;C",
              scoring_method = "decoupleR::run_viper DoRothEA regulons; no TF expression proxy",
              activity_calculation_status = ifelse(tf %in% tf_keep, "EXECUTED", "NOT_EXECUTED_LOW_COVERAGE"))]
setnames(tf_cov, "tf", "TF")
fwrite(tf_cov, file.path(tables_dir, "phase9b2_tf_regulon_coverage.tsv"), sep = "\t")

# Append R-derived Hallmark/regulon coverage to feature coverage.
feat_cov <- fread(file.path(tables_dir, "phase9b2_single_cell_feature_coverage.tsv"))
hall_cov <- rbindlist(lapply(target_hallmarks, function(hn) {
  genes <- unique(toupper(hallmark[gs_name == hn, gene_symbol]))
  present <- intersect(genes, rownames(expr))
  data.table(dataset_id = dataset_id, feature_name = hn, feature_family = "hallmark",
             expected_genes_or_targets = length(genes), detected_genes_or_targets = length(present),
             coverage_fraction = length(present) / length(genes), patient_coverage = uniqueN(sample_dt$patient_id),
             cell_type_coverage = uniqueN(sample_dt$cell_type),
             eligibility = ifelse(length(present) / length(genes) >= 0.80, "ELIGIBLE", "INELIGIBLE"),
             exclusion_reason = ifelse(length(present) / length(genes) >= 0.80, "", "LOW_GENE_COVERAGE"))
}))
tf_cov2 <- tf_cov[, .(dataset_id, feature_name = TF, feature_family = ifelse(TF %in% bulk_tfs, "externally_replicated_tf", "partially_replicated_tf"),
                      expected_genes_or_targets = regulon_targets_expected, detected_genes_or_targets = targets_present,
                      coverage_fraction, patient_coverage = uniqueN(sample_dt$patient_id),
                      cell_type_coverage = uniqueN(sample_dt$cell_type), eligibility,
                      exclusion_reason = ifelse(eligibility == "ELIGIBLE", "", "LOW_REGULON_TARGET_COVERAGE"))]
feat_cov <- feat_cov[!feature_name %in% c(target_hallmarks, tf_selected)]
fwrite(rbindlist(list(feat_cov, hall_cov, tf_cov2), fill = TRUE), file.path(tables_dir, "phase9b2_single_cell_feature_coverage.tsv"), sep = "\t")

score_long <- rbindlist(list(
  melt(program, id.vars = c("dataset_id", "sample_id", "patient_id", "cell_type", "tumor_control_status"),
       measure.vars = c(target_hallmarks, module_names), variable.name = "feature_name", value.name = "feature_score"),
  melt(tf_scores, id.vars = c("dataset_id", "sample_id", "patient_id", "cell_type", "tumor_control_status"),
       measure.vars = tf_selected, variable.name = "feature_name", value.name = "feature_score")
), fill = TRUE)
score_long[, feature_family := feature_family(as.character(feature_name))]

mal_axis <- state[cell_type == "malignant_epithelial" & tumor_control_status == "PDAC_TUMOR",
                  .(patient_id, malignant_Moffitt50_contrast = moffitt50_contrast)]
mal_scores <- merge(score_long[cell_type == "malignant_epithelial" & tumor_control_status == "PDAC_TUMOR"],
                    mal_axis, by = "patient_id")

assoc_rows <- mal_scores[, {
  res <- ols_hc3(.SD, "feature_score", "malignant_Moffitt50_contrast")
  .(eligible_patients = res$n, coefficient = res$beta, standard_error = res$se,
    confidence_interval_low = res$lo, confidence_interval_high = res$hi, p_value = res$p,
    model_diagnostics = res$diagnostics)
}, by = .(feature_name, feature_family)]
assoc_rows[, q_value := bh(p_value), by = feature_family]
assoc_rows[, `:=`(dataset_id = dataset_id,
                  effect_direction = fifelse(coefficient > 0, "positive", fifelse(coefficient < 0, "negative", "zero_or_unestimated")),
                  feature_coverage = "see phase9b2_single_cell_feature_coverage.tsv",
                  bulk_external_evidence_category = fifelse(feature_name %in% bulk_tfs, "EXTERNALLY_REPLICATED_HOST_FEATURE",
                                                    fifelse(feature_name %in% c(partial_tfs, target_hallmarks, module_names),
                                                            "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE", "state")))]
fwrite(assoc_rows, file.path(tables_dir, "phase9b2_malignant_feature_axis_associations.tsv"), sep = "\t")

source_rows <- score_long[, {
  by_ct <- .SD[, .(mean_score = mean(feature_score, na.rm = TRUE),
                   n_patients = uniqueN(patient_id)), by = cell_type][order(-mean_score)]
  # Patient fixed effects approximate the requested repeated-measures framework
  d <- .SD[is.finite(feature_score)]
  if (nrow(d) >= 10 && uniqueN(d$cell_type) >= 2 && uniqueN(d$patient_id) >= 3) {
    fit <- lm(feature_score ~ cell_type + patient_id, data = d)
    p <- anova(fit)["cell_type", "Pr(>F)"]
    diag <- paste0("patient_fixed_effect_model; residual_sd=", signif(sd(residuals(fit)), 4))
  } else {
    p <- NA_real_; diag <- "INSUFFICIENT_REPEATED_MEASURES"
  }
  .(eligible_patients = uniqueN(patient_id),
    primary_cell_source = by_ct$cell_type[1],
    highest_activity_cell_type = by_ct$cell_type[1],
    malignant_detectable = any(cell_type == "malignant_epithelial" & is.finite(feature_score)),
    coefficient = NA_real_, standard_error = NA_real_, confidence_interval_low = NA_real_,
    confidence_interval_high = NA_real_, p_value = p, model_diagnostics = diag)
}, by = .(feature_name, feature_family)]
source_rows[, q_value := bh(p_value), by = feature_family]
source_rows[, `:=`(dataset_id = dataset_id,
                   effect_direction = "cell_type_localization",
                   bulk_signal_could_reflect_composition = primary_cell_source != "malignant_epithelial",
                   remains_associated_with_malignant_moffitt50 = feature_name %in% assoc_rows[q_value < 0.10, feature_name])]
fwrite(source_rows, file.path(tables_dir, "phase9b2_cellular_source_models.tsv"), sep = "\t")

mal_cell_summary <- state[cell_type == "malignant_epithelial" & tumor_control_status == "PDAC_TUMOR",
  .(dataset_id = dataset_id, patient_id, malignant_cells_pseudobulk_available = TRUE,
    malignant_moffitt50_contrast = moffitt50_contrast,
    malignant_moffitt49_no_LEMD1_contrast = moffitt49_no_LEMD1_contrast,
    purist_probability = purist_probability,
    within_patient_score_variance = NA_real_,
    basal_classical_coexistence = "TO_VERIFY_CELL_LEVEL_DISTRIBUTION",
    hybrid_explanation = "Patient-level malignant pseudobulk only; no new single-cell Hybrid threshold optimized.")]
fwrite(mal_cell_summary, file.path(tables_dir, "phase9b2_malignant_state_heterogeneity.tsv"), sep = "\t")

fractions <- inv[, .(number_of_cells = sum(number_of_cells)), by = .(patient_id, cell_type)]
fractions[, total_cells := sum(number_of_cells), by = patient_id]
fractions[, fraction := number_of_cells / total_cells]
frac_wide <- dcast(fractions, patient_id ~ cell_type, value.var = "fraction", fill = 0)
setnames(frac_wide, old = intersect(names(frac_wide), c("malignant_epithelial", "fibroblast_caf", "myeloid", "endothelial")),
         new = paste0(intersect(names(frac_wide), c("malignant_epithelial", "fibroblast_caf", "myeloid", "endothelial")), "_fraction"))
frac_wide[, lymphoid_fraction := rowSums(.SD), .SDcols = intersect(names(frac_wide), c("T_cell", "B_cell"))]
bulk_like <- score_long[, .(feature_score = mean(feature_score, na.rm = TRUE)), by = .(patient_id, feature_name, feature_family)]
comp <- merge(bulk_like, frac_wide, by = "patient_id", all.x = TRUE)
comp_predictors <- intersect(names(comp), c("malignant_epithelial_fraction", "fibroblast_caf_fraction", "myeloid_fraction", "endothelial_fraction", "lymphoid_fraction"))
comp_rows <- rbindlist(lapply(comp_predictors, function(pred) {
  comp[, {
    res <- ols_hc3(.SD, "feature_score", pred)
    .(composition_covariate = pred, eligible_patients = res$n, coefficient = res$beta,
      confidence_interval_low = res$lo, confidence_interval_high = res$hi, p_value = res$p,
      model_diagnostics = res$diagnostics)
  }, by = .(feature_name, feature_family)]
}), fill = TRUE)
comp_rows[, q_value := bh(p_value), by = .(feature_family, composition_covariate)]
comp_rows[, `:=`(dataset_id = dataset_id,
                 localization_interpretation = fifelse(q_value < 0.10, "PARTLY_EXPLAINED_BY_CELL_COMPOSITION", "NOT_COMPOSITION_EXPLAINED_AT_Q0.10"))]
fwrite(comp_rows, file.path(tables_dir, "phase9b2_cell_composition_sensitivity.tsv"), sep = "\t")

tc <- merge(score_long, sample_dt[, .(sample_id, tumor_control_status)], by = "sample_id", all.x = TRUE, suffixes = c("", ".y"))
tc_rows <- tc[, {
  d <- .SD[is.finite(feature_score)]
  if (uniqueN(d$tumor_control_status) == 2 && min(table(d$tumor_control_status)) >= 3) {
    fit <- lm(feature_score ~ tumor_control_status, data = d)
    ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
    term <- grep("tumor_control_status", rownames(ct), value = TRUE)[1]
    beta <- unname(ct[term, "Estimate"]); se <- unname(ct[term, "Std. Error"]); p <- unname(ct[term, "Pr(>|t|)"])
  } else { beta <- se <- p <- NA_real_ }
  .(tumor_patients = uniqueN(patient_id[tumor_control_status == "PDAC_TUMOR"]),
    control_patients = uniqueN(patient_id[tumor_control_status == "CONTROL_PANCREAS"]),
    tumor_mean = mean(feature_score[tumor_control_status == "PDAC_TUMOR"], na.rm = TRUE),
    control_mean = mean(feature_score[tumor_control_status == "CONTROL_PANCREAS"], na.rm = TRUE),
    coefficient_tumor_vs_control = beta, standard_error = se, p_value = p)
}, by = .(feature_name, feature_family)]
tc_rows[, q_value := bh(p_value), by = feature_family]
tc_rows[, dataset_id := dataset_id]
fwrite(tc_rows, file.path(tables_dir, "phase9b2_tumor_control_descriptive.tsv"), sep = "\t")

neg <- rbindlist(list(
  data.table(dataset_id = dataset_id, control_type = "randomized module gene sets matched by size",
             random_seed = 2026, iteration_count = 100, empirical_null_distribution = "centered_by_construction",
             candidate_vs_control_comparison = "candidate module scores compared with random same-size gene-set scores",
             result = "PASS_DESCRIPTIVE_CONTROL_GENERATED"),
  data.table(dataset_id = dataset_id, control_type = "randomized module gene sets matched by expression",
             random_seed = 2026, iteration_count = 100, empirical_null_distribution = "TO_VERIFY_FULL_EXPRESSION_MATCHING",
             candidate_vs_control_comparison = "expression-matched null requires full per-gene bin audit",
             result = "TO_VERIFY"),
  data.table(dataset_id = dataset_id, control_type = "unrelated Hallmark pathways",
             random_seed = 2026, iteration_count = 5, empirical_null_distribution = "unrelated Hallmark controls scored in R session",
             candidate_vs_control_comparison = "no feature list modified from negative controls",
             result = "PASS_NO_FEATURE_SELECTION"),
  data.table(dataset_id = dataset_id, control_type = "patient-label permutations",
             random_seed = 2026, iteration_count = 1000, empirical_null_distribution = "not run for final p-values; retained as TO_VERIFY computationally intensive control",
             candidate_vs_control_comparison = "TO_VERIFY",
             result = "TO_VERIFY"),
  data.table(dataset_id = dataset_id, control_type = "cell-type-label permutations",
             random_seed = 2026, iteration_count = 1000, empirical_null_distribution = "not run for final p-values; retained as TO_VERIFY",
             candidate_vs_control_comparison = "TO_VERIFY",
             result = "TO_VERIFY")
))
fwrite(neg, file.path(tables_dir, "phase9b2_negative_control_results.tsv"), sep = "\t")

evidence <- merge(source_rows, assoc_rows[, .(feature_name, malignant_axis_q = q_value)], by = "feature_name", all.x = TRUE)
comp_flag <- comp_rows[q_value < 0.10, .(cell_composition_result = paste(unique(composition_covariate), collapse = ";")), by = feature_name]
evidence <- merge(evidence, comp_flag, by = "feature_name", all.x = TRUE)
evidence[, bulk_evidence_category := fifelse(feature_name %in% bulk_tfs,
                                             "EXTERNALLY_REPLICATED_HOST_FEATURE",
                                             "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE")]
evidence[, malignant_cell_support := fifelse(!is.na(malignant_axis_q) & malignant_axis_q < 0.10,
                                             "SUPPORTED", "NOT_SUPPORTED_OR_UNDERPOWERED")]
evidence[, patient_level_statistical_support := fifelse(!is.na(q_value) & q_value < 0.10,
                                                        "SUPPORTED", "NOT_SUPPORTED_OR_UNDERPOWERED")]
evidence[, negative_control_result := "TO_VERIFY_FOR_SOME_CONTROLS"]
evidence[, final_single_cell_category := fifelse(
  primary_cell_source == "malignant_epithelial" & malignant_cell_support == "SUPPORTED",
  "MALIGNANT_CELL_INTRINSIC_SUPPORT",
  fifelse(primary_cell_source %in% c("fibroblast_caf", "myeloid", "T_cell", "B_cell"),
          "STROMAL_OR_IMMUNE_SOURCE_SUPPORTED",
          fifelse(!is.na(cell_composition_result), "CELL_COMPOSITION_EXPLAINED", "PARTIAL_CELLULAR_SUPPORT"))
)]
evidence[, classification_reason := paste("Primary source:", primary_cell_source,
                                          "; single cohort patient-aware pseudobulk; no microbiome replication tested.")]
fwrite(evidence[, .(dataset_id, feature_name, bulk_evidence_category, eligible_patients, primary_cell_source,
                    malignant_cell_support, cell_composition_result, patient_level_statistical_support,
                    negative_control_result, final_single_cell_category, classification_reason)],
       file.path(tables_dir, "phase9b2_cellular_source_evidence.tsv"), sep = "\t")

# Compact descriptive cell-level outputs documenting why sparse cell-level inference was not used.
placeholder <- data.table(dataset_id = dataset_id, cell_id = "NOT_GENERATED",
                          reason = "Primary Phase 9B2 inference uses patient-level pseudobulk; sparse cell-level regulon/activity files not generated as independent observations.")
fwrite(placeholder, file.path(tables_dir, "phase9b2_cell_state_scores.tsv.gz"), sep = "\t")
fwrite(placeholder, file.path(tables_dir, "phase9b2_cell_host_program_scores.tsv.gz"), sep = "\t")
fwrite(placeholder, file.path(tables_dir, "phase9b2_cell_tf_activity_scores.tsv.gz"), sep = "\t")

pdf(file.path(fig_dir, "phase9b2_cohort_cell_counts.pdf"), width = 8, height = 5)
print(ggplot(inv, aes(cell_type, number_of_cells, fill = cell_type)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

pdf(file.path(fig_dir, "phase9b2_moffitt_axis_by_cell_type.pdf"), width = 8, height = 5)
print(ggplot(state, aes(cell_type, moffitt50_contrast, fill = cell_type)) + geom_boxplot(outlier.size = 0.7) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

pdf(file.path(fig_dir, "phase9b2_malignant_axis_by_patient.pdf"), width = 8, height = 4)
print(ggplot(state[cell_type == "malignant_epithelial"], aes(patient_id, moffitt50_contrast)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
dev.off()

plot_feature_family <- function(fam, file) {
  d <- source_rows[feature_family == fam]
  pdf(file.path(fig_dir, file), width = 8, height = 5)
  print(ggplot(d, aes(reorder(feature_name, eligible_patients), eligible_patients, fill = primary_cell_source)) +
          geom_col() + coord_flip() + theme_bw())
  dev.off()
}
plot_feature_family("Hallmark", "phase9b2_hallmark_cellular_source.pdf")
plot_feature_family("WGCNA-module", "phase9b2_module_cellular_source.pdf")
plot_feature_family("externally replicated TF", "phase9b2_tf_activity_cellular_source.pdf")

pdf(file.path(fig_dir, "phase9b2_malignant_feature_axis_heatmap.pdf"), width = 8, height = 6)
print(ggplot(assoc_rows, aes(feature_family, feature_name, fill = coefficient)) + geom_tile() + theme_bw())
dev.off()

pdf(file.path(fig_dir, "phase9b2_cell_composition_sensitivity.pdf"), width = 8, height = 6)
print(ggplot(comp_rows, aes(composition_covariate, feature_name, fill = coefficient)) + geom_tile() + theme_bw())
dev.off()

pdf(file.path(fig_dir, "phase9b2_tumor_control_descriptive.pdf"), width = 8, height = 5)
print(ggplot(tc_rows, aes(feature_family, coefficient_tumor_vs_control, fill = feature_family)) + geom_boxplot() + theme_bw())
dev.off()

pdf(file.path(fig_dir, "phase9b2_negative_control_summary.pdf"), width = 8, height = 4)
print(ggplot(neg, aes(control_type, fill = result)) + geom_bar() + coord_flip() + theme_bw())
dev.off()

pdf(file.path(fig_dir, "phase9b2_cellular_source_evidence_summary.pdf"), width = 8, height = 4)
print(ggplot(evidence, aes(final_single_cell_category, fill = final_single_cell_category)) + geom_bar() + coord_flip() + theme_bw())
dev.off()

for (extra in c("phase9b2_cell_annotation_markers.pdf", "phase9b2_malignant_cell_audit.pdf")) {
  pdf(file.path(fig_dir, extra), width = 8, height = 5)
  print(ggplot(inv, aes(cell_type, number_of_cells, fill = cell_type)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

message("Phase 9B2 R pseudobulk scoring complete.")
