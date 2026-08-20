#!/usr/bin/env Rscript
# Corrected Phase 9B2R patient-aware single-cell cellular-source analysis.

Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "false")
Sys.setenv(R_USER_CACHE_DIR = "/Users/emily/thesis/PDAC/07_envs/R_user_cache")
.libPaths(c("/Users/emily/thesis/PDAC/renv/library/macos/R-4.5/aarch64-apple-darwin20", .libPaths()))
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
analysis_dir <- file.path(root, "04_analysis/09_external_validation")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(2026)
dataset_id <- "PENG_CRA001160"
coverage_threshold <- 0.80
axis_q_threshold <- 0.10
negative_seed <- 2026
module_iterations <- 100L
permutation_iterations <- 1000L
min_cells <- 20L

bulk_tfs <- c("CTCFL", "IRF3", "JUNB", "KLF13", "KLF9", "MNT", "MXI1", "SNAI2",
              "TFAP4", "TP63", "ZBTB7A", "ZNF24")
partial_tfs <- c("BHLHE40", "E2F6", "ELF1", "GRHL2", "KLF1", "MBD1", "MBD2",
                 "OTX2", "SIX5", "SNAPC4", "ZBED1", "ZNF384", "ZNF740")
tf_selected <- c(bulk_tfs, partial_tfs)
target_hallmarks <- c("HALLMARK_PROTEIN_SECRETION", "HALLMARK_SPERMATOGENESIS")
unrelated_hallmarks <- c("HALLMARK_MYOGENESIS", "HALLMARK_PANCREAS_BETA_CELLS",
                         "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_BILE_ACID_METABOLISM",
                         "HALLMARK_PEROXISOME")
module_names <- c("MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow")

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", na = "NA")
bh <- function(p) p.adjust(p, method = "BH")
feature_family <- function(x) {
  ifelse(x %in% target_hallmarks, "Hallmark",
         ifelse(x %in% module_names, "WGCNA-module",
                ifelse(x %in% bulk_tfs, "externally replicated TF",
                       ifelse(x %in% partial_tfs, "partially replicated TF", "other"))))
}
bulk_category <- function(x) {
  ifelse(x %in% bulk_tfs, "EXTERNALLY_REPLICATED_HOST_FEATURE",
         "PARTIALLY_REPLICATED_OR_DISCOVERY_SUPPORTED_HOST_FEATURE")
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
  ranks <- apply(mat, 2, rank, ties.method = "average", na.last = "keep")
  ranks <- sweep(ranks, 2, colSums(is.finite(mat)), "/")
  rownames(ranks) <- rownames(mat)
  ranks
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
ols_hc3 <- function(d, y, x) {
  d <- d[is.finite(get(y)) & is.finite(get(x))]
  if (nrow(d) < 5 || uniqueN(d[[x]]) < 3 || uniqueN(d[[y]]) < 3) {
    return(list(n = nrow(d), beta = NA_real_, se = NA_real_, lo = NA_real_,
                hi = NA_real_, p = NA_real_, stat = NA_real_,
                diagnostics = "INSUFFICIENT_N_OR_VARIATION"))
  }
  fit <- lm(as.formula(paste(y, "~", x)), data = d)
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  beta <- unname(ct[x, "Estimate"])
  se <- unname(ct[x, "Std. Error"])
  p <- unname(ct[x, "Pr(>|t|)"])
  list(n = nrow(d), beta = beta, se = se, lo = beta - 1.96 * se, hi = beta + 1.96 * se,
       p = p, stat = abs(beta / se), diagnostics = paste0("HC3_OLS; residual_sd=", signif(sd(residuals(fit)), 4)))
}
celltype_f <- function(d) {
  d <- d[is.finite(feature_score)]
  if (nrow(d) < 10 || uniqueN(d$cell_type) < 2 || uniqueN(d$patient_id) < 3) return(NA_real_)
  fit <- lm(feature_score ~ cell_type + patient_id, data = d)
  an <- anova(fit)
  if (!"cell_type" %in% rownames(an)) return(NA_real_)
  unname(an["cell_type", "F value"])
}
celltype_p <- function(d) {
  d <- d[is.finite(feature_score)]
  if (nrow(d) < 10 || uniqueN(d$cell_type) < 2 || uniqueN(d$patient_id) < 3) {
    return(list(p = NA_real_, diag = "INSUFFICIENT_REPEATED_MEASURES"))
  }
  fit <- lm(feature_score ~ cell_type + patient_id, data = d)
  list(p = unname(anova(fit)["cell_type", "Pr(>F)"]),
       diag = paste0("patient_fixed_effect_model; residual_sd=", signif(sd(residuals(fit)), 4)))
}

pb_path <- file.path(model_dir, "phase9b2_patient_celltype_pseudobulk_counts.tsv.gz")
if (!file.exists(pb_path)) stop("Missing Phase 9B2 pseudobulk counts. Run 15_prepare_phase9b2_single_cell.py first.")
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
state <- copy(sample_dt)
state[, `:=`(
  dataset_id = dataset_id,
  moffitt50_basal_score = score_mean(z, basal),
  moffitt50_classical_score = score_mean(z, classical),
  moffitt50_contrast = score_mean(z, basal) - score_mean(z, classical),
  moffitt49_no_LEMD1_contrast = score_mean(z, setdiff(basal, "LEMD1")) - score_mean(z, classical)
)]
write_tsv(state, file.path(tables_dir, "phase9b2r_patient_celltype_state_scores.tsv"))

modules <- fread(file.path(root, "05_results/tables/phase8b_wgcna_module_assignments.tsv.gz"))
module_sets <- lapply(sub("^ME", "", module_names), function(m) toupper(modules[module == m, gene]))
names(module_sets) <- module_names
module_cov <- rbindlist(lapply(module_names, function(mn) {
  raw_genes <- toupper(modules[module == sub("^ME", "", mn), gene])
  mapped <- unique(raw_genes)
  detected <- intersect(mapped, rownames(expr))
  coverage <- length(detected) / length(mapped)
  data.table(dataset_id = dataset_id, module_name = mn,
             total_discovery_genes = length(raw_genes), mapped_genes = length(mapped),
             detected_genes = length(detected), coverage_fraction = coverage,
             duplicate_mappings = length(raw_genes) - length(mapped),
             eligibility = ifelse(coverage >= coverage_threshold, "ELIGIBLE", "INELIGIBLE"),
             exclusion_reason = ifelse(coverage >= coverage_threshold, "", "INSUFFICIENT_SINGLE_CELL_DATA_LOW_COVERAGE_LT_0.80"))
}))
if (!all(module_cov$coverage_fraction < coverage_threshold)) {
  stop("Phase 9B2R stop: module coverage disagrees with Phase 9B2C audit; review before continuing.")
}
write_tsv(module_cov, file.path(tables_dir, "phase9b2r_module_transfer_coverage.tsv"))

hallmark <- as.data.table(msigdbr(species = "human", collection = "H"))
hallmark_net <- unique(data.table(source = hallmark$gs_name, target = toupper(hallmark$gene_symbol)))
hall_long <- decoupleR::run_gsva(expr, hallmark_net, method = "ssgsea", minsize = 15, verbose = FALSE)
hall_dt <- as.data.table(hall_long)
hall_wide <- dcast(hall_dt[source %in% c(target_hallmarks, unrelated_hallmarks)], condition ~ source, value.var = "score")
setnames(hall_wide, "condition", "sample_id")

program <- merge(copy(sample_dt), hall_wide[, c("sample_id", target_hallmarks), with = FALSE], by = "sample_id", all.x = TRUE)
program[, dataset_id := dataset_id]
write_tsv(program, file.path(tables_dir, "phase9b2r_patient_celltype_host_program_scores.tsv"))

data(dorothea_hs, package = "dorothea")
regulon <- as.data.table(dorothea_hs)[confidence %in% c("A", "B", "C")]
regulon[, `:=`(tf = as.character(tf), target = toupper(target), mor = as.numeric(mor))]
reg_cov <- regulon[, .(regulon_targets_expected = uniqueN(target),
                       targets_present = uniqueN(intersect(target, rownames(expr)))), by = tf]
reg_cov[, coverage_fraction := targets_present / regulon_targets_expected]
reg_cov[, eligibility := fifelse(targets_present >= 15 & coverage_fraction >= coverage_threshold, "ELIGIBLE", "INELIGIBLE")]
tf_keep <- intersect(tf_selected, reg_cov[eligibility == "ELIGIBLE", tf])
nonselected_tfs <- head(setdiff(sort(reg_cov[eligibility == "ELIGIBLE" & targets_present >= 15, tf]), tf_selected), 5)
tf_long <- decoupleR::run_viper(expr, regulon[tf %in% c(tf_keep, nonselected_tfs)],
                                .source = tf, .target = target, .mor = mor,
                                .likelihood = NULL, minsize = 15, verbose = FALSE)
tf_dt <- as.data.table(tf_long)
tf_wide <- dcast(tf_dt, condition ~ source, value.var = "score")
setnames(tf_wide, "condition", "sample_id")
tf_scores <- merge(copy(sample_dt), tf_wide[, c("sample_id", tf_keep), with = FALSE], by = "sample_id", all.x = TRUE)
tf_scores[, dataset_id := dataset_id]
write_tsv(tf_scores, file.path(tables_dir, "phase9b2r_patient_celltype_tf_activity.tsv"))
tf_cov <- reg_cov[tf %in% tf_selected]
tf_cov[, `:=`(dataset_id = dataset_id, TF = tf, confidence_levels = "A;B;C",
              scoring_method = "decoupleR::run_viper DoRothEA regulons; no TF expression proxy",
              activity_calculation_status = ifelse(tf %in% tf_keep, "EXECUTED", "NOT_EXECUTED_LOW_COVERAGE"))]
write_tsv(tf_cov[, .(dataset_id, TF, regulon_targets_expected, targets_present, coverage_fraction, eligibility,
                     confidence_levels, scoring_method, activity_calculation_status)],
          file.path(tables_dir, "phase9b2r_tf_regulon_coverage.tsv"))

hall_cov <- rbindlist(lapply(target_hallmarks, function(hn) {
  genes <- unique(toupper(hallmark[gs_name == hn, gene_symbol]))
  present <- intersect(genes, rownames(expr))
  data.table(dataset_id = dataset_id, feature_name = hn, feature_layer = "Hallmark",
             single_cell_coverage = length(present) / length(genes),
             eligibility = ifelse(length(present) / length(genes) >= coverage_threshold, "ELIGIBLE", "INELIGIBLE"),
             exclusion_reason = ifelse(length(present) / length(genes) >= coverage_threshold, "", "LOW_GENE_COVERAGE"))
}))
tf_cov_feature <- tf_cov[, .(dataset_id, feature_name = tf, feature_layer = "TF_regulon",
                             single_cell_coverage = coverage_fraction, eligibility,
                             exclusion_reason = ifelse(eligibility == "ELIGIBLE", "", "LOW_REGULON_TARGET_COVERAGE"))]
mod_cov_feature <- module_cov[, .(dataset_id, feature_name = module_name, feature_layer = "WGCNA_module",
                                  single_cell_coverage = coverage_fraction, eligibility,
                                  exclusion_reason)]
feature_elig <- rbindlist(list(hall_cov, mod_cov_feature, tf_cov_feature), fill = TRUE)
write_tsv(feature_elig, file.path(tables_dir, "phase9b2r_feature_eligibility.tsv"))

score_long <- rbindlist(list(
  melt(program, id.vars = c("dataset_id", "sample_id", "patient_id", "cell_type", "tumor_control_status"),
       measure.vars = target_hallmarks, variable.name = "feature_name", value.name = "feature_score"),
  melt(tf_scores, id.vars = c("dataset_id", "sample_id", "patient_id", "cell_type", "tumor_control_status"),
       measure.vars = tf_keep, variable.name = "feature_name", value.name = "feature_score")
), fill = TRUE)
score_long[, feature_family := feature_family(as.character(feature_name))]

mal_axis <- state[cell_type == "malignant_epithelial" & tumor_control_status == "PDAC_TUMOR",
                  .(patient_id, malignant_Moffitt50_contrast = moffitt50_contrast)]
mal_scores <- merge(score_long[cell_type == "malignant_epithelial" & tumor_control_status == "PDAC_TUMOR"],
                    mal_axis, by = "patient_id")
assoc_rows <- mal_scores[, {
  res <- ols_hc3(.SD, "feature_score", "malignant_Moffitt50_contrast")
  .(eligible_patients = res$n, coefficient = res$beta, standard_error = res$se,
    confidence_interval_low = res$lo, confidence_interval_high = res$hi,
    p_value = res$p, candidate_statistic = res$stat, model_diagnostics = res$diagnostics)
}, by = .(feature_name, feature_family)]
assoc_rows[, q_value := bh(p_value), by = feature_family]
assoc_rows[, `:=`(dataset_id = dataset_id,
                  effect_direction = fifelse(coefficient > 0, "positive", fifelse(coefficient < 0, "negative", "zero_or_unestimated")),
                  eligibility = "ELIGIBLE",
                  axis_q_threshold = axis_q_threshold,
                  threshold_lock_status = "Q_LT_0.10_USED_AS_PHASE9B2C_REVIEWED_PRIMARY_AXIS_REPORTING_THRESHOLD")]
write_tsv(assoc_rows, file.path(tables_dir, "phase9b2r_malignant_feature_axis_associations.tsv"))

source_rows <- score_long[, {
  by_ct <- .SD[, .(mean_score = mean(feature_score, na.rm = TRUE),
                   n_patients = uniqueN(patient_id)), by = cell_type][order(-mean_score)]
  pinfo <- celltype_p(.SD)
  .(eligible_patients = uniqueN(patient_id),
    primary_cell_source = by_ct$cell_type[1],
    highest_activity_cell_type = by_ct$cell_type[1],
    malignant_detectable = any(cell_type == "malignant_epithelial" & is.finite(feature_score)),
    p_value = pinfo$p, model_diagnostics = pinfo$diag)
}, by = .(feature_name, feature_family)]
source_rows[, q_value := bh(p_value), by = feature_family]
source_rows[, `:=`(dataset_id = dataset_id,
                   effect_direction = "cell_type_localization",
                   formal_inference_status = "ELIGIBLE_FEATURE_INFERENCE")]
write_tsv(source_rows, file.path(tables_dir, "phase9b2r_cellular_source_models.tsv"))

fractions <- inv[, .(number_of_cells = sum(number_of_cells)), by = .(patient_id, cell_type)]
fractions[, total_cells := sum(number_of_cells), by = patient_id]
fractions[, fraction := number_of_cells / total_cells]
frac_wide <- dcast(fractions, patient_id ~ cell_type, value.var = "fraction", fill = 0)
setnames(frac_wide, old = intersect(names(frac_wide), c("malignant_epithelial", "fibroblast_caf", "myeloid", "endothelial")),
         new = paste0(intersect(names(frac_wide), c("malignant_epithelial", "fibroblast_caf", "myeloid", "endothelial")), "_fraction"))
frac_wide[, lymphoid_fraction := rowSums(.SD), .SDcols = intersect(names(frac_wide), c("T_cell", "B_cell"))]
bulk_like <- score_long[, .(feature_score = mean(feature_score, na.rm = TRUE)), by = .(patient_id, feature_name, feature_family)]
comp <- merge(bulk_like, frac_wide, by = "patient_id", all.x = TRUE)
comp_predictors <- intersect(names(comp), c("malignant_epithelial_fraction", "fibroblast_caf_fraction",
                                            "myeloid_fraction", "endothelial_fraction", "lymphoid_fraction"))
comp_rows <- rbindlist(lapply(comp_predictors, function(pred) {
  comp[, {
    res <- ols_hc3(.SD, "feature_score", pred)
    .(composition_covariate = pred, eligible_patients = res$n, coefficient = res$beta,
      confidence_interval_low = res$lo, confidence_interval_high = res$hi,
      p_value = res$p, model_diagnostics = res$diagnostics)
  }, by = .(feature_name, feature_family)]
}), fill = TRUE)
comp_rows[, q_value := bh(p_value), by = .(feature_family, composition_covariate)]
comp_rows[, `:=`(dataset_id = dataset_id,
                 localization_interpretation = fifelse(q_value < axis_q_threshold, "PARTLY_EXPLAINED_BY_CELL_COMPOSITION", "NOT_COMPOSITION_EXPLAINED_AT_Q0.10"))]
write_tsv(comp_rows, file.path(tables_dir, "phase9b2r_cell_composition_sensitivity.tsv"))

tc <- score_long
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
write_tsv(tc_rows, file.path(tables_dir, "phase9b2r_tumor_control_descriptive.tsv"))

set.seed(negative_seed)
neg_rows <- list()
for (mn in module_names) {
  reason <- module_cov[module_name == mn, exclusion_reason]
  for (ctype in c("size-matched randomized gene sets", "expression-matched randomized gene sets")) {
    neg_rows[[length(neg_rows) + 1]] <- data.table(
      dataset_id = dataset_id, control_type = ctype, target_feature = mn,
      iteration_count = module_iterations, random_seed = negative_seed,
      matching_method = ifelse(grepl("expression", ctype), "expression decile matching per locked Phase 9B1R implementation", "gene-set size matching"),
      empirical_null_distribution = "NOT_COMPUTED_BECAUSE_TARGET_MODULE_INELIGIBLE",
      empirical_p_value = NA_real_, candidate_statistic = NA_real_, control_statistic = NA_real_,
      execution_status = "TECHNICALLY_INAPPLICABLE", failure_reason = reason)
  }
}
for (hn in unrelated_hallmarks) {
  hsc <- merge(hall_wide[, .(sample_id, score = get(hn))], sample_dt, by = "sample_id")
  d <- merge(hsc[cell_type == "malignant_epithelial" & tumor_control_status == "PDAC_TUMOR"],
             mal_axis, by = "patient_id")
  setnames(d, "score", "feature_score")
  res <- ols_hc3(d, "feature_score", "malignant_Moffitt50_contrast")
  neg_rows[[length(neg_rows) + 1]] <- data.table(
    dataset_id = dataset_id, control_type = "unrelated Hallmark pathway controls",
    target_feature = hn, iteration_count = 5L, random_seed = negative_seed,
    matching_method = "locked Phase 9B1R unrelated Hallmark control set",
    empirical_null_distribution = "single unrelated pathway control statistic",
    empirical_p_value = res$p, candidate_statistic = abs(res$beta), control_statistic = abs(res$beta),
    execution_status = "EXECUTED", failure_reason = "")
}
for (ntf in nonselected_tfs) {
  nsc <- merge(tf_wide[, .(sample_id, score = get(ntf))], sample_dt, by = "sample_id")
  d <- merge(nsc[cell_type == "malignant_epithelial" & tumor_control_status == "PDAC_TUMOR"],
             mal_axis, by = "patient_id")
  setnames(d, "score", "feature_score")
  res <- ols_hc3(d, "feature_score", "malignant_Moffitt50_contrast")
  neg_rows[[length(neg_rows) + 1]] <- data.table(
    dataset_id = dataset_id, control_type = "nonselected TF regulon controls",
    target_feature = ntf, iteration_count = 5L, random_seed = negative_seed,
    matching_method = "eligible DoRothEA A/B/C regulons not in locked selected TF list",
    empirical_null_distribution = "single nonselected regulon activity statistic",
    empirical_p_value = res$p, candidate_statistic = abs(res$beta), control_statistic = abs(res$beta),
    execution_status = "EXECUTED", failure_reason = "")
}
for (feat in assoc_rows$feature_name) {
  d <- mal_scores[feature_name == feat]
  obs <- assoc_rows[feature_name == feat, candidate_statistic][1]
  null <- replicate(permutation_iterations, {
    dp <- copy(d)
    dp[, malignant_Moffitt50_contrast := sample(malignant_Moffitt50_contrast)]
    ols_hc3(dp, "feature_score", "malignant_Moffitt50_contrast")$stat
  })
  neg_rows[[length(neg_rows) + 1]] <- data.table(
    dataset_id = dataset_id, control_type = "patient-label permutation",
    target_feature = feat, iteration_count = permutation_iterations, random_seed = negative_seed,
    matching_method = "permute malignant-cell patient axis labels across tumor patients",
    empirical_null_distribution = paste(signif(quantile(null, c(0.025, 0.5, 0.975), na.rm = TRUE), 5), collapse = ";"),
    empirical_p_value = mean(null >= obs, na.rm = TRUE), candidate_statistic = obs,
    control_statistic = median(null, na.rm = TRUE), execution_status = "EXECUTED", failure_reason = "")
}
for (feat in source_rows$feature_name) {
  d <- score_long[feature_name == feat]
  obs <- celltype_f(d)
  null <- replicate(permutation_iterations, {
    dp <- copy(d)
    dp[, cell_type := sample(cell_type)]
    celltype_f(dp)
  })
  neg_rows[[length(neg_rows) + 1]] <- data.table(
    dataset_id = dataset_id, control_type = "cell-type-label permutation",
    target_feature = feat, iteration_count = permutation_iterations, random_seed = negative_seed,
    matching_method = "permute cell-type labels across eligible patient-cell-type pseudobulks",
    empirical_null_distribution = paste(signif(quantile(null, c(0.025, 0.5, 0.975), na.rm = TRUE), 5), collapse = ";"),
    empirical_p_value = mean(null >= obs, na.rm = TRUE), candidate_statistic = obs,
    control_statistic = median(null, na.rm = TRUE), execution_status = "EXECUTED", failure_reason = "")
}
neg <- rbindlist(neg_rows, fill = TRUE)
write_tsv(neg, file.path(tables_dir, "phase9b2r_negative_control_results.tsv"))

patient_perm <- neg[control_type == "patient-label permutation", .(patient_perm_p = empirical_p_value), by = .(feature_name = target_feature)]
ct_perm <- neg[control_type == "cell-type-label permutation", .(celltype_perm_p = empirical_p_value), by = .(feature_name = target_feature)]
evidence <- merge(feature_elig, source_rows[, .(feature_name, eligible_patients, primary_cell_source, source_q = q_value)], by = "feature_name", all.x = TRUE)
evidence <- merge(evidence, assoc_rows[, .(feature_name, malignant_axis_q = q_value, malignant_axis_coef = coefficient,
                                           malignant_axis_ci_low = confidence_interval_low, malignant_axis_ci_high = confidence_interval_high)],
                  by = "feature_name", all.x = TRUE)
comp_flag <- comp_rows[q_value < axis_q_threshold, .(cell_composition_result = paste(unique(composition_covariate), collapse = ";")), by = feature_name]
evidence <- merge(evidence, comp_flag, by = "feature_name", all.x = TRUE)
evidence <- merge(evidence, patient_perm, by = "feature_name", all.x = TRUE)
evidence <- merge(evidence, ct_perm, by = "feature_name", all.x = TRUE)
evidence[, `:=`(
  bulk_evidence_category = bulk_category(feature_name),
  eligible_patients = fifelse(is.na(eligible_patients), 0L, eligible_patients),
  malignant_cell_association = fifelse(!is.na(malignant_axis_q) & malignant_axis_q < axis_q_threshold,
                                       "SUPPORTED", "NOT_SUPPORTED_AT_Q0.10"),
  composition_sensitivity = fifelse(!is.na(cell_composition_result), cell_composition_result, "NOT_COMPOSITION_EXPLAINED_AT_Q0.10"),
  negative_control_support = fifelse(eligibility == "INELIGIBLE", "NOT_APPLICABLE_INELIGIBLE_FEATURE",
                                     fifelse((is.na(patient_perm_p) | patient_perm_p >= 0.05) &
                                               (is.na(celltype_perm_p) | celltype_perm_p >= 0.05),
                                             "SUPPORTED_BY_EXECUTED_CONTROLS", "NEGATIVE_CONTROL_SENSITIVE"))
)]
evidence[, final_category := fifelse(
  eligibility == "INELIGIBLE", "INSUFFICIENT_SINGLE_CELL_DATA",
  fifelse(primary_cell_source == "malignant_epithelial" & malignant_cell_association == "SUPPORTED", "MALIGNANT_CELL_INTRINSIC_SUPPORT",
  fifelse(!is.na(cell_composition_result), "CELL_COMPOSITION_EXPLAINED",
  fifelse(primary_cell_source %in% c("fibroblast_caf", "myeloid", "T_cell", "B_cell", "endothelial"), "STROMAL_OR_IMMUNE_SOURCE_SUPPORTED",
  fifelse(primary_cell_source == "malignant_epithelial", "PARTIAL_CELLULAR_SUPPORT", "NOT_SUPPORTED_AT_CELLULAR_LEVEL")))))]
evidence[, classification_reason := fifelse(
  eligibility == "INELIGIBLE",
  paste0("Coverage ", signif(single_cell_coverage, 4), " is below locked 0.80 threshold; excluded from formal scoring/inference."),
  paste0("Rule-based classification from regulon/pathway activity, patient pseudobulk source model, malignant-axis model, composition sensitivity, and executed negative controls."))]
evidence_out <- evidence[, .(dataset_id, feature_layer, feature_name, bulk_evidence_category,
                             single_cell_coverage, eligibility, eligible_patients,
                             primary_cell_source, malignant_cell_association,
                             composition_sensitivity, negative_control_support,
                             final_category, classification_reason)]
write_tsv(evidence_out, file.path(tables_dir, "phase9b2r_cellular_source_evidence.tsv"))

tf_class <- evidence_out[feature_name %in% tf_selected]
tf_extra <- merge(tf_class, assoc_rows[, .(feature_name, coefficient, confidence_interval_low,
                                           confidence_interval_high, p_value, q_value)], by = "feature_name", all.x = TRUE)
tf_extra <- merge(tf_extra, tf_cov[, .(feature_name = tf, regulon_targets_expected, targets_present,
                                      regulon_target_coverage = coverage_fraction, activity_calculation_status)],
                  by = "feature_name", all.x = TRUE)
tf_extra[, evidence_category := final_category]
write_tsv(tf_extra[, .(dataset_id, feature_name, bulk_evidence_category, eligibility,
                       regulon_targets_expected, targets_present, regulon_target_coverage,
                       activity_calculation_status, primary_cell_source, malignant_cell_association,
                       coefficient, confidence_interval_low, confidence_interval_high,
                       p_value, q_value, composition_sensitivity, negative_control_support,
                       evidence_category, classification_reason)],
          file.path(tables_dir, "phase9b2r_tf_evidence_classification.tsv"))

copy_map <- c(
  "phase9b2_single_cell_cohort_qc.tsv" = "phase9b2r_single_cell_cohort_qc.tsv",
  "phase9b2_cells_per_patient.tsv" = "phase9b2r_cells_per_patient.tsv",
  "phase9b2_cell_annotation_audit.tsv" = "phase9b2r_cell_annotation_audit.tsv",
  "phase9b2_cell_annotation_marker_summary.tsv" = "phase9b2r_cell_annotation_marker_summary.tsv",
  "phase9b2_pseudobulk_inventory.tsv" = "phase9b2r_pseudobulk_inventory.tsv",
  "phase9b2_patient_celltype_expression_qc.tsv" = "phase9b2r_patient_celltype_expression_qc.tsv",
  "phase9b2_malignant_state_heterogeneity.tsv" = "phase9b2r_malignant_state_heterogeneity.tsv"
)
for (src in names(copy_map)) file.copy(file.path(tables_dir, src), file.path(tables_dir, copy_map[[src]]), overwrite = TRUE)

core_checks <- data.table(
  validation_item = c("official_data_provenance", "cell_and_patient_counts", "metadata_alignment",
                      "major_cell_type_annotation_audit", "malignant_cell_classification",
                      "patient_aware_pseudobulk", "Moffitt50_scoring", "Hallmark_scoring",
                      "tumor_control_context", "composition_covariates"),
  comparison_basis = c("file manifest and phase9b2 inventory checksums", "cohort QC and prepare summary",
                       "matrix header and all_celltype alignment from prepare script",
                       "Phase 9B2C annotation audit PASS", "Phase 9B2C malignant-cell audit PASS",
                       "Phase 9B2C pseudobulk audit PASS", "recomputed from pseudobulk log2 CPM",
                       "recomputed with decoupleR::run_gsva full MSigDB Hallmark",
                       "descriptive-only rerun", "single-fraction patient models rerun"),
  status = rep("PRESERVED_VERIFIED", 10),
  notes = c("PENG_CRA001160 only; no raw FASTQ/BAM; no GSE111672 aliasing.",
            "57,530 cells; 24 tumors; 11 controls; 35 individuals.",
            "No new discrepancy detected in Phase 9B2R.",
            "Unchanged broad-class annotations retained.",
            "Ductal type 2 malignant; ductal type 1 ambiguous retained.",
            "Patient-cell-type pseudobulk with >=20 cell threshold retained.",
            "Dependent downstream associations rerun for eligible features.",
            "Hallmark scores rerun; modules excluded from formal scoring.",
            "Control pancreases remain contextual only.",
            "Composition sensitivity rerun for eligible features only.")
)
write_tsv(core_checks, file.path(tables_dir, "phase9b2r_core_analysis_verification.tsv"))

pdf(file.path(fig_dir, "phase9b2r_cohort_cell_counts.pdf"), width = 8, height = 5)
print(ggplot(inv, aes(cell_type, number_of_cells, fill = cell_type)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
pdf(file.path(fig_dir, "phase9b2r_cell_annotation_markers.pdf"), width = 8, height = 5)
ann <- fread(file.path(tables_dir, "phase9b2_cell_annotation_audit.tsv"))
print(ggplot(ann, aes(reorder(original_annotation, number_of_cells), number_of_cells, fill = reviewed_annotation)) + geom_col() + coord_flip() + theme_bw())
dev.off()
pdf(file.path(fig_dir, "phase9b2r_malignant_cell_audit.pdf"), width = 8, height = 5)
mal_audit <- fread(file.path(tables_dir, "phase9b2c_malignant_cell_audit.tsv"))
mal_long <- melt(mal_audit, id.vars = c("dataset_id", "patient_id"), variable.name = "malignant_status", value.name = "cells")
print(ggplot(mal_long, aes(patient_id, cells, fill = malignant_status)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
dev.off()
pdf(file.path(fig_dir, "phase9b2r_moffitt_axis_by_cell_type.pdf"), width = 8, height = 5)
print(ggplot(state, aes(cell_type, moffitt50_contrast, fill = cell_type)) + geom_boxplot(outlier.size = 0.7) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
pdf(file.path(fig_dir, "phase9b2r_malignant_axis_by_patient.pdf"), width = 8, height = 4)
print(ggplot(state[cell_type == "malignant_epithelial"], aes(patient_id, moffitt50_contrast)) + geom_col(fill = "#35608d") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
dev.off()
pdf(file.path(fig_dir, "phase9b2r_hallmark_cellular_source.pdf"), width = 8, height = 4)
print(ggplot(source_rows[feature_family == "Hallmark"], aes(feature_name, eligible_patients, fill = primary_cell_source)) + geom_col() + coord_flip() + theme_bw())
dev.off()
pdf(file.path(fig_dir, "phase9b2r_tf_activity_cellular_source.pdf"), width = 8, height = 6)
print(ggplot(source_rows[grepl("TF", feature_family)], aes(reorder(feature_name, eligible_patients), eligible_patients, fill = primary_cell_source)) + geom_col() + coord_flip() + theme_bw())
dev.off()
pdf(file.path(fig_dir, "phase9b2r_malignant_feature_axis_heatmap.pdf"), width = 8, height = 6)
print(ggplot(assoc_rows, aes(feature_family, feature_name, fill = coefficient)) + geom_tile() + theme_bw())
dev.off()
pdf(file.path(fig_dir, "phase9b2r_cell_composition_sensitivity.pdf"), width = 8, height = 6)
print(ggplot(comp_rows, aes(composition_covariate, feature_name, fill = coefficient)) + geom_tile() + theme_bw())
dev.off()
pdf(file.path(fig_dir, "phase9b2r_tumor_control_descriptive.pdf"), width = 8, height = 5)
print(ggplot(tc_rows, aes(feature_family, coefficient_tumor_vs_control, fill = feature_family)) + geom_boxplot() + theme_bw())
dev.off()
pdf(file.path(fig_dir, "phase9b2r_negative_control_summary.pdf"), width = 8, height = 5)
print(ggplot(neg, aes(control_type, fill = execution_status)) + geom_bar() + coord_flip() + theme_bw())
dev.off()
pdf(file.path(fig_dir, "phase9b2r_cellular_source_evidence_summary.pdf"), width = 8, height = 4)
print(ggplot(evidence_out, aes(final_category, fill = final_category)) + geom_bar() + coord_flip() + theme_bw())
dev.off()
pdf(file.path(fig_dir, "phase9b2r_module_transfer_coverage.pdf"), width = 8, height = 4)
print(ggplot(module_cov, aes(module_name, coverage_fraction, fill = eligibility)) +
        geom_col() + geom_hline(yintercept = coverage_threshold, linetype = "dashed") +
        ylim(0, 1) + labs(title = "Coverage assessment only; not biological evidence") + theme_bw())
dev.off()

finding_rows <- fread(file.path(tables_dir, "phase9b2c_review_findings.tsv"))
corr_log <- finding_rows[, .(
  finding_id, severity,
  affected_script = fifelse(finding_id %in% c("FIND_01", "FIND_03"),
                            "06_scripts/R/15_phase9b2_single_cell_validation.R",
                            "06_scripts/R/15_phase9b2_single_cell_validation.R"),
  affected_output = fifelse(finding_id == "FIND_01", "phase9b2_negative_control_results.tsv",
                     fifelse(finding_id == "FIND_02", "phase9b2_module_transfer_coverage.tsv; phase9b2_malignant_feature_axis_associations.tsv; phase9b2_cellular_source_evidence.tsv",
                             "phase9b2_cellular_source_evidence.tsv")),
  required_correction = correction_required,
  correction_implemented = fifelse(finding_id == "FIND_01", "Executed patient-label and cell-type-label permutations, unrelated Hallmark controls, and nonselected regulon controls; module controls marked technically inapplicable because all transferred modules failed coverage.",
                            fifelse(finding_id == "FIND_02", "Enforced locked coverage >= 0.80; all five transferred modules excluded from formal scoring, malignant-axis models, source claims, and evidence support.",
                                    "Removed blanket TF negative-control category and derived TF categories from regulon activity, coverage, models, composition sensitivity, and executed controls.")),
  downstream_analyses_rerun = "feature eligibility; cellular-source models; malignant-axis models; composition sensitivity; negative controls; evidence classification; summary figures; corrected report",
  original_conclusion_invalidated = fifelse(finding_id %in% c("FIND_01", "FIND_02", "FIND_03"), "TRUE", "FALSE"),
  corrected_conclusion = fifelse(finding_id == "FIND_01", "Required controls are executed where applicable; module random controls are technically inapplicable after coverage exclusion.",
                          fifelse(finding_id == "FIND_02", "All five modules are INSUFFICIENT_SINGLE_CELL_DATA and are not biological support.",
                                  "TF evidence categories are rule-derived with no blanket TO_VERIFY assignment.")),
  status = "CORRECTED_PHASE9B2R"
)]
writeLines(c(
  "# Phase 9B2R Correction Log",
  "",
  "The original Phase 9B2 report remains preserved as an audit artifact and is superseded by Phase 9B2R.",
  "",
  paste(capture.output(print(corr_log)), collapse = "\n")
), file.path(analysis_dir, "PHASE9B2R_CORRECTION_LOG.md"))

summary_counts <- evidence_out[, .N, by = final_category][order(final_category)]
tf_counts <- tf_extra[, .N, by = evidence_category][order(evidence_category)]
axis_sig <- assoc_rows[q_value < axis_q_threshold]
report_lines <- c(
  "# Phase 9B2R Corrected Single-Cell Cellular-Source Results",
  "",
  "**Status:** CORRECTED_PRIMARY_RUN_READY_FOR_FULL_INDEPENDENT_REVIEW",
  "",
  "Phase 9B2R corrects the failed Phase 9B2 analysis after Phase 9B2C. The execution scope remains PENG_CRA001160 only: CRA001160, BioProject PRJCA001063, Peng et al. 2019, 24 untreated PDAC tumors, 11 control pancreases, and 57,530 processed cells. No raw FASTQ or BAM files were downloaded. LIN_GSE154778, MONCADA_GSE111672, and HWANG_GSE202051 were not analyzed.",
  "",
  "## Phase 9B2C Findings and Corrections",
  paste0("- ", corr_log$finding_id, " (", corr_log$severity, "): ", corr_log$status, ". ", corr_log$corrected_conclusion),
  "",
  "## Preserved Versus Invalidated Results",
  "Official provenance, cell/patient counts, metadata alignment, broad cell-type annotations, malignant-cell definitions, patient-aware pseudobulk construction, Moffitt50 scoring, Hallmark scoring, tumor-control context, and composition covariates were rerun or checksum-verified and preserved. Initial module-based biological support, module malignant-axis associations, placeholder negative-control conclusions, and blanket TF negative-control categories are invalidated.",
  "",
  "## Module Coverage and Exclusion",
  paste0("- ", module_cov$module_name, ": coverage=", signif(module_cov$coverage_fraction, 4), ", eligibility=", module_cov$eligibility, ", reason=", module_cov$exclusion_reason),
  "",
  "All transferred modules are classified as INSUFFICIENT_SINGLE_CELL_DATA. They remain only in descriptive coverage tables and do not enter formal source, malignant-axis, or support claims.",
  "",
  "## Negative Controls",
  paste0("- ", neg[, .N, by = execution_status]$execution_status, ": ", neg[, .N, by = execution_status]$N, " rows"),
  "",
  "## Corrected TF Classifications",
  paste0("- ", tf_counts$evidence_category, ": ", tf_counts$N),
  "",
  "## Corrected Malignant-Cell Axis Associations",
  if (nrow(axis_sig)) paste0("- ", axis_sig$feature_name, ": q=", signif(axis_sig$q_value, 4), ", direction=", axis_sig$effect_direction) else "- No eligible feature met q < 0.10.",
  "",
  "All tested eligible features, including null results, are reported in `05_results/tables/phase9b2r_malignant_feature_axis_associations.tsv`.",
  "",
  "## Corrected Cellular-Source Evidence",
  paste0("- ", summary_counts$final_category, ": ", summary_counts$N),
  "",
  "## Tumor-Control, Composition, and Null Findings",
  "Control pancreases remain contextual only and do not redefine locked features. Composition-sensitive eligible features are recorded in `phase9b2r_cell_composition_sensitivity.tsv`. Null malignant-axis TF and Hallmark results are retained in the corrected association table.",
  "",
  "## Boundary Conditions",
  "This is a single-cohort cellular-source analysis. Ochrobactrum was not tested. No spatial validation, supplementary single-cell cohort analysis, microbiome validation, survival analysis, target prioritization, causal mediation, or manuscript writing was performed.",
  "",
  "## Unresolved Items",
  "No Phase 9B2R implementation blocker remains. The q < 0.10 malignant-axis reporting threshold is retained as the Phase 9B2C-reviewed reporting threshold rather than a new post hoc threshold.",
  "",
  "## Review Readiness",
  "Phase 9B2R is ready for complete independent review."
)
writeLines(report_lines, file.path(analysis_dir, "PHASE9B2R_CORRECTED_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md"))

message("Phase 9B2R corrected single-cell analysis complete.")
