#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(sandwich)
  library(lmtest)
  library(msigdbr)
  library(GSVA)
  library(progeny)
  library(dorothea)
  library(decoupleR)
  library(viper)
  library(WGCNA)
  library(fgsea)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 2)
set.seed(2026)

root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(root) != "PDAC") root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = TRUE)
Sys.setenv(R_USER_CACHE_DIR = file.path(root, "07_envs", "R_user_cache"))

tables_dir <- file.path(root, "05_results", "tables")
figures_dir <- file.path(root, "05_results", "figures")
models_dir <- file.path(root, "05_results", "models", "phase8b")
gene_full_dir <- file.path(tables_dir, "phase8b_host_gene_full")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_full_dir, recursive = TRUE, showWarnings = FALSE)

primary_taxa <- c("Azoarcus", "Candida", "Ensifer", "Cutibacterium", "Chryseobacterium",
                  "Ochrobactrum", "Burkholderia", "Rhizobium", "Herbaspirillum")
secondary_taxa <- c("Staphylococcus", "Citrobacter")
extreme_samples <- c("Basal-like1", "Hybrid18", "Hybrid23")

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
zscore <- function(x) as.numeric(scale(x))
sign_chr <- function(x) ifelse(is.na(x), NA_character_, ifelse(x > 0, "positive", ifelse(x < 0, "negative", "zero")))
bh <- function(p) p.adjust(p, method = "BH")

read_feature_matrix <- function(path, feature_col) {
  dt <- fread(path)
  stopifnot(feature_col %in% names(dt))
  feats <- dt[[feature_col]]
  mat <- as.matrix(dt[, setdiff(names(dt), feature_col), with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feats
  mat
}

make_sample_map <- function() {
  expr_map <- fread(file.path(root, "01_metadata", "expression_sample_crosswalk.tsv"))
  micro_map <- fread(file.path(root, "01_metadata", "microbiome_sample_crosswalk.tsv"))
  manifest <- fread(file.path(root, "01_metadata", "sample_manifest.tsv"))
  list(expr = expr_map, micro = micro_map, manifest = manifest)
}

align_inputs <- function(expr, micro) {
  maps <- make_sample_map()
  expr_dt <- maps$expr[expression_column %in% colnames(expr)]
  micro_dt <- maps$micro[microbiome_matrix_sample %in% colnames(micro)]
  common <- sort(intersect(expr_dt$patient_id, micro_dt$patient_id))
  expr_cols <- expr_dt[match(common, patient_id), expression_column]
  micro_cols <- micro_dt[match(common, patient_id), microbiome_matrix_sample]
  list(patient_id = common, expr = expr[, expr_cols, drop = FALSE], micro = micro[, micro_cols, drop = FALSE],
       expr_sample_id = expr_cols, micro_sample_id = micro_cols)
}

ols_hc3 <- function(y, x, covariate = NULL) {
  df <- data.frame(y = as.numeric(y), genus = zscore(x))
  if (!is.null(covariate)) df$covariate <- zscore(covariate)
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10 || sd(df$y) == 0 || sd(df$genus) == 0) {
    return(list(n = nrow(df), coefficient = NA_real_, robust_se_HC3 = NA_real_, ci_lower = NA_real_,
                ci_upper = NA_real_, t_statistic_HC3 = NA_real_, p_value = NA_real_, r_squared = NA_real_,
                adjusted_r_squared = NA_real_, vif_genus = NA_real_, condition_number = NA_real_))
  }
  fit <- lm(if (is.null(covariate)) y ~ genus else y ~ genus + covariate, data = df)
  vc <- sandwich::vcovHC(fit, type = "HC3")
  ct <- lmtest::coeftest(fit, vcov. = vc)
  idx <- "genus"
  beta <- unname(coef(fit)[idx])
  se <- unname(sqrt(diag(vc))[idx])
  tval <- unname(ct[idx, "t value"])
  pval <- unname(ct[idx, "Pr(>|t|)"])
  ci <- beta + c(-1, 1) * qt(0.975, df = fit$df.residual) * se
  mm <- model.matrix(fit)
  vif <- NA_real_
  if (!is.null(covariate)) {
    vif <- 1 / (1 - summary(lm(genus ~ covariate, data = df))$r.squared)
  } else {
    vif <- 1
  }
  list(n = nrow(df), coefficient = beta, robust_se_HC3 = se, ci_lower = ci[1], ci_upper = ci[2],
       t_statistic_HC3 = tval, p_value = pval, r_squared = summary(fit)$r.squared,
       adjusted_r_squared = summary(fit)$adj.r.squared, vif_genus = vif,
       condition_number = kappa(mm, exact = TRUE))
}

model_feature_table <- function(score_dt, micro_mat, taxa, collection, evidence_dt, risk_dt, rclr_flags) {
  rows <- list()
  for (taxon in taxa) {
    x <- as.numeric(micro_mat[taxon, score_dt$patient_id])
    for (feature in setdiff(names(score_dt), "patient_id")) {
      res <- ols_hc3(score_dt[[feature]], x)
      rows[[length(rows) + 1]] <- data.table(
        taxon = taxon, host_feature = feature, host_feature_collection = collection,
        n = res$n, coefficient = res$coefficient, robust_se_HC3 = res$robust_se_HC3,
        ci_lower = res$ci_lower, ci_upper = res$ci_upper, t_statistic_HC3 = res$t_statistic_HC3,
        p_value = res$p_value, r_squared = res$r_squared, adjusted_r_squared = res$adjusted_r_squared,
        multiple_testing_family = paste(taxon, collection, sep = "_"),
        phase7_evidence_category = evidence_dt[.(taxon), evidence_category],
        contamination_risk_category = risk_dt[.(taxon), contamination_risk_category],
        RCLR_DIRECTION_SENSITIVE = taxon %in% rclr_flags
      )
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, bh_q_value := bh(p_value), by = .(taxon, host_feature_collection)]
  setcolorder(out, c("taxon", "host_feature", "host_feature_collection", "n", "coefficient",
                     "robust_se_HC3", "ci_lower", "ci_upper", "t_statistic_HC3", "p_value",
                     "bh_q_value", "r_squared", "adjusted_r_squared", "multiple_testing_family",
                     "phase7_evidence_category", "contamination_risk_category", "RCLR_DIRECTION_SENSITIVE"))
  out
}

long_to_score_dt <- function(long, patient_id) {
  wide <- dcast(as.data.table(long), condition ~ source, value.var = "score")
  setnames(wide, "condition", "sample_id")
  wide[, patient_id := patient_id[match(sample_id, names(patient_id))]]
  wide[, sample_id := NULL]
  setcolorder(wide, "patient_id")
  wide
}

runtime_rows <- list()
add_runtime <- function(check, observed, expected, passed, notes = "") {
  runtime_rows[[length(runtime_rows) + 1]] <<- data.table(
    validation_check = check, observed = as.character(observed), expected = as.character(expected),
    passed = as.logical(passed), notes = notes
  )
}

pkg_names <- c("data.table", "limma", "sandwich", "lmtest", "msigdbr", "GSVA", "progeny", "dorothea",
               "decoupleR", "viper", "WGCNA", "fgsea", "clusterProfiler", "ReactomePA", "ggplot2")
pkg_versions <- data.table(package = pkg_names,
                           loaded = vapply(pkg_names, requireNamespace, logical(1), quietly = TRUE),
                           version = vapply(pkg_names, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, character(1)))

expr <- read_feature_matrix(file.path(root, "03_processed", "expression", "GSE172356_expression_log2_analysis_ready.tsv.gz"), "gene")
expr_filtered <- read_feature_matrix(file.path(root, "03_processed", "expression", "GSE172356_expression_filtered_normalized.tsv.gz"), "gene")
gene_annot <- fread(file.path(root, "03_processed", "expression", "GSE172356_gene_annotation.tsv"))
micro <- read_feature_matrix(file.path(root, "03_processed", "microbiome", "PRJNA719915_genus_primary_CLR.tsv.gz"), "taxon")
rclr <- read_feature_matrix(file.path(root, "03_processed", "microbiome", "sensitivity", "MICRO_SENS_ROBUST_CLR_robust_clr.tsv.gz"), "taxon")
micro_no_contam <- read_feature_matrix(file.path(root, "03_processed", "microbiome", "sensitivity", "MICRO_SENS_NO_CONTAMINANTS_centered_log_ratio.tsv.gz"), "taxon")
micro_excl_extreme <- read_feature_matrix(file.path(root, "03_processed", "microbiome", "sensitivity", "MICRO_SENS_EXCLUDE_EXTREME_centered_log_ratio.tsv.gz"), "taxon")

aligned <- align_inputs(expr, micro)
expr <- aligned$expr
micro <- aligned$micro
sample_ids <- aligned$patient_id
micro_patient_names <- aligned$patient_id
names(micro_patient_names) <- colnames(micro)
colnames(expr) <- aligned$patient_id
colnames(micro) <- aligned$patient_id
names(sample_ids) <- aligned$patient_id
rclr_aligned <- align_inputs(expr = read_feature_matrix(file.path(root, "03_processed", "expression", "GSE172356_expression_log2_analysis_ready.tsv.gz"), "gene"), micro = rclr)
rclr <- rclr_aligned$micro; colnames(rclr) <- rclr_aligned$patient_id
no_contam_aligned <- align_inputs(expr = read_feature_matrix(file.path(root, "03_processed", "expression", "GSE172356_expression_log2_analysis_ready.tsv.gz"), "gene"), micro = micro_no_contam)
micro_no_contam <- no_contam_aligned$micro; colnames(micro_no_contam) <- no_contam_aligned$patient_id
excl_extreme_aligned <- align_inputs(expr = read_feature_matrix(file.path(root, "03_processed", "expression", "GSE172356_expression_log2_analysis_ready.tsv.gz"), "gene"), micro = micro_excl_extreme)
micro_excl_extreme <- excl_extreme_aligned$micro; colnames(micro_excl_extreme) <- excl_extreme_aligned$patient_id

tme <- fread(file.path(root, "01_metadata", "host_tme_covariates.tsv"))[match(aligned$patient_id, patient_id)]
phase7_primary <- fread(file.path(root, "05_results", "tables", "phase7b_primary_genus_associations.tsv"))
phase7_evidence <- fread(file.path(root, "05_results", "tables", "phase7b_genus_evidence_classification.tsv"))
phase7c_audit <- fread(file.path(root, "05_results", "tables", "phase7c_primary_candidate_audit.tsv"))
risk <- fread(file.path(root, "05_results", "tables", "phase6c_retained_taxa_with_contamination_flags.tsv"))
evidence_vec <- phase7_evidence[, .(evidence_category = evidence_category[1]), by = genus]
risk_vec <- risk[, .(contamination_risk_category = contamination_risk_category[1]), by = genus]
setkey(evidence_vec, genus); setkey(risk_vec, genus)
rclr_flags <- setdiff(primary_taxa, "Ochrobactrum")

add_runtime("exactly_62_aligned_patients", length(aligned$patient_id), 62, length(aligned$patient_id) == 62)
add_runtime("all_nine_primary_taxa_present_primary_clr", paste(primary_taxa %in% rownames(micro), collapse = ";"), "all TRUE", all(primary_taxa %in% rownames(micro)))
add_runtime("all_nine_primary_taxa_present_rclr", paste(primary_taxa %in% rownames(rclr), collapse = ";"), "all TRUE", all(primary_taxa %in% rownames(rclr)))
phase7_dirs <- phase7_primary[genus %in% primary_taxa, .(genus, phase7_direction = sign_chr(coefficient))]
add_runtime("taxon_directions_match_phase7b", paste(phase7_dirs$genus, phase7_dirs$phase7_direction, sep = ":", collapse = ";"), "loaded Phase 7B locked directions", nrow(phase7_dirs) == 9)
add_runtime("expression_microbiome_patient_order_aligned", paste(head(colnames(expr)), collapse = ";"), "identical patient_id order", identical(colnames(expr), colnames(micro)))
add_runtime("no_duplicated_patients", anyDuplicated(aligned$patient_id), 0, anyDuplicated(aligned$patient_id) == 0)
add_runtime("no_duplicated_genes", anyDuplicated(rownames(expr)), 0, anyDuplicated(rownames(expr)) == 0)
add_runtime("no_duplicated_taxa", anyDuplicated(rownames(micro)), 0, anyDuplicated(rownames(micro)) == 0)
add_runtime("no_missing_or_infinite_primary_values", sum(!is.finite(expr)) + sum(!is.finite(micro[primary_taxa, ])), 0, sum(!is.finite(expr)) + sum(!is.finite(micro[primary_taxa, ])) == 0)
add_runtime("required_phase8_packages_load_from_renv", paste(pkg_versions$package[pkg_versions$loaded], collapse = ";"), paste(pkg_names, collapse = ";"), all(pkg_versions$loaded))
add_runtime("renv_project_library_active", .libPaths()[1], file.path(root, "renv", "library"), grepl(file.path("PDAC", "renv", "library"), .libPaths()[1], fixed = TRUE))
add_runtime("no_public_subtype_labels_used_for_feature_selection_or_model_tuning", "script constants use primary taxa and unsupervised WGCNA parameters only", "TRUE", TRUE)
runtime <- rbindlist(runtime_rows)
write_tsv(runtime, file.path(tables_dir, "phase8b_runtime_validation.tsv"))
write_tsv(pkg_versions, file.path(tables_dir, "phase8b_runtime_package_versions.tsv"))
if (any(!runtime$passed[match(c("exactly_62_aligned_patients", "all_nine_primary_taxa_present_primary_clr", "taxon_directions_match_phase7b",
                                "expression_microbiome_patient_order_aligned", "required_phase8_packages_load_from_renv"), runtime$validation_check)])) {
  stop("Phase 8B hard-stop runtime validation failed; see phase8b_runtime_validation.tsv")
}

message("Scoring Hallmark pathways")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_net <- unique(hallmark[, c("gs_name", "gene_symbol")])
setnames(hallmark_net, c("source", "target"))
hallmark_long <- decoupleR::run_gsva(expr, hallmark_net, method = "ssgsea", minsize = 15, verbose = FALSE)
hallmark_scores <- long_to_score_dt(hallmark_long, sample_ids)
write_tsv(hallmark_scores, file.path(tables_dir, "phase8b_hallmark_activity_scores.tsv.gz"))

hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
coverage <- rbindlist(lapply(names(hallmark_sets), function(gs) {
  genes <- unique(hallmark_sets[[gs]])
  data.table(collection = "MSigDB_Hallmark", feature = gs, collection_version = unique(hallmark$db_version)[1],
             expected_features = 50, genes_expected = length(genes),
             genes_observed = sum(genes %in% rownames(expr)),
             coverage_fraction = sum(genes %in% rownames(expr)) / length(genes),
             removed_for_inadequate_coverage = sum(genes %in% rownames(expr)) < 15)
}))

message("Scoring PROGENy pathways")
progeny_scores_mat <- progeny::progeny(expr, organism = "Human", top = 100, scale = TRUE, perm = 1, z_scores = FALSE)
progeny_scores <- data.table(patient_id = rownames(progeny_scores_mat), as.data.table(progeny_scores_mat))
write_tsv(progeny_scores, file.path(tables_dir, "phase8b_progeny_activity_scores.tsv.gz"))
coverage <- rbind(coverage, data.table(collection = "PROGENy", feature = colnames(progeny_scores_mat),
                                       collection_version = as.character(packageVersion("progeny")),
                                       expected_features = 14, genes_expected = 100, genes_observed = NA_integer_,
                                       coverage_fraction = NA_real_, removed_for_inadequate_coverage = FALSE,
                                       score_direction = "positive_or_negative_pathway_activity_weighted_model_top100"), fill = TRUE)
write_tsv(coverage, file.path(tables_dir, "phase8b_pathway_gene_coverage.tsv"))

message("Scoring DoRothEA/VIPER TF activities")
data(dorothea_hs, package = "dorothea")
regulon <- as.data.table(dorothea_hs)[confidence %in% c("A", "B", "C")]
reg_cov <- regulon[, .(targets_expected = uniqueN(target), targets_observed = uniqueN(intersect(target, rownames(expr))),
                       target_coverage_fraction = uniqueN(intersect(target, rownames(expr))) / uniqueN(target)),
                   by = .(tf)]
reg_cov[, retained := targets_observed >= 15]
write_tsv(reg_cov[, .(TF = tf, regulon_version = as.character(packageVersion("dorothea")),
                      confidence_levels = "A;B;C", targets_expected, targets_observed,
                      target_coverage_fraction, retained, activity_score_direction = "VIPER_NES_positive_or_negative_TF_activity")],
          file.path(tables_dir, "phase8b_tf_regulon_coverage.tsv"))
tf_long <- decoupleR::run_viper(expr, regulon[tf %in% reg_cov[retained == TRUE, tf]],
                                .source = tf, .target = target, .mor = mor,
                                .likelihood = NULL, minsize = 15, verbose = FALSE)
tf_scores <- long_to_score_dt(tf_long, sample_ids)
write_tsv(tf_scores, file.path(tables_dir, "phase8b_tf_activity_scores.tsv.gz"))

message("Running primary host-feature association models")
hallmark_assoc <- model_feature_table(hallmark_scores, micro, primary_taxa, "MSigDB_Hallmark", evidence_vec, risk_vec, rclr_flags)
progeny_assoc <- model_feature_table(progeny_scores, micro, primary_taxa, "PROGENy", evidence_vec, risk_vec, rclr_flags)
pathway_assoc <- rbind(hallmark_assoc, progeny_assoc, fill = TRUE)
write_tsv(pathway_assoc, file.path(tables_dir, "phase8b_primary_pathway_associations.tsv"))
tf_assoc <- model_feature_table(tf_scores, micro, primary_taxa, "DoRothEA_VIPER", evidence_vec, risk_vec, rclr_flags)
write_tsv(tf_assoc, file.path(tables_dir, "phase8b_primary_tf_associations.tsv"))

supported <- rbind(pathway_assoc[, .(taxon, host_feature, host_feature_collection, coefficient, p_value, bh_q_value)],
                   tf_assoc[, .(taxon, host_feature, host_feature_collection, coefficient, p_value, bh_q_value)], fill = TRUE)
supported[, primary_supported := bh_q_value < 0.05]
candidates <- supported[p_value < 0.05 | bh_q_value < 0.05]
if (nrow(candidates) == 0) candidates <- supported[order(p_value)][1:min(.N, 25)]

score_lookup <- list(MSigDB_Hallmark = hallmark_scores, PROGENy = progeny_scores, DoRothEA_VIPER = tf_scores)
message("Running sensitivity models")
cov_rows <- list(); trans_rows <- list(); loo_rows <- list()
for (i in seq_len(nrow(candidates))) {
  cand <- candidates[i]
  scores <- score_lookup[[cand$host_feature_collection]]
  y <- scores[[cand$host_feature]]
  x0 <- micro[cand$taxon, scores$patient_id]
  base <- ols_hc3(y, x0)
  for (covn in c("inferred_tumor_purity", "immune_score", "stromal_score")) {
    res <- ols_hc3(y, x0, tme[[covn]])
    cov_rows[[length(cov_rows) + 1]] <- data.table(taxon = cand$taxon, host_feature = cand$host_feature,
      host_feature_collection = cand$host_feature_collection, model = paste0("Model_", substr(covn, 1, 1)),
      covariate = covn, n = res$n, genus_coefficient = res$coefficient,
      attenuation = abs(res$coefficient) - abs(base$coefficient), sign_change = sign(res$coefficient) != sign(base$coefficient),
      ci_lower = res$ci_lower, ci_upper = res$ci_upper, p_value = res$p_value,
      VIF = res$vif_genus, condition_number = res$condition_number,
      robustness_interpretation = ifelse(is.na(res$p_value), "TO_VERIFY", ifelse(res$p_value < 0.05 & sign(res$coefficient) == sign(base$coefficient), "covariate_stable", "composition_sensitive")))
  }
  reps <- list(primary_CLR = micro, rCLR = rclr, contaminant_exclusion = micro_no_contam, technical_extreme_exclusion = micro_excl_extreme)
  vals <- list()
  for (repn in names(reps)) {
    mat <- reps[[repn]]
    common <- intersect(scores$patient_id, colnames(mat))
    if (!(cand$taxon %in% rownames(mat)) || length(common) < 10) next
    res <- ols_hc3(scores[match(common, patient_id)][[cand$host_feature]], mat[cand$taxon, common])
    vals[[repn]] <- res
    trans_rows[[length(trans_rows) + 1]] <- data.table(taxon = cand$taxon, host_feature = cand$host_feature,
      host_feature_collection = cand$host_feature_collection, representation = repn, n = res$n,
      coefficient = res$coefficient, direction = sign_chr(res$coefficient), p_value = res$p_value,
      primary_CLR_direction = sign_chr(base$coefficient),
      direction_agreement = sign(res$coefficient) == sign(base$coefficient),
      transformation_sensitive_label = repn == "rCLR" && (sign(res$coefficient) != sign(base$coefficient) || res$p_value >= 0.05))
  }
  loo_coef <- c(); loo_p <- c()
  for (pid in scores$patient_id) {
    keep <- scores$patient_id != pid
    res <- ols_hc3(y[keep], x0[keep])
    loo_coef <- c(loo_coef, res$coefficient); loo_p <- c(loo_p, res$p_value)
  }
  loo_rows[[length(loo_rows) + 1]] <- data.table(taxon = cand$taxon, host_feature = cand$host_feature,
    host_feature_collection = cand$host_feature_collection, primary_coefficient = base$coefficient,
    min_LOO_coefficient = min(loo_coef, na.rm = TRUE), max_LOO_coefficient = max(loo_coef, na.rm = TRUE),
    LOO_sign_stable = all(sign(loo_coef) == sign(base$coefficient), na.rm = TRUE),
    min_LOO_p_value = min(loo_p, na.rm = TRUE), max_LOO_p_value = max(loo_p, na.rm = TRUE),
    sample_exclusion_stability = ifelse(all(sign(loo_coef) == sign(base$coefficient), na.rm = TRUE), "stable_direction", "sample_sensitive"))
}
cov_sens <- rbindlist(cov_rows, fill = TRUE); cov_sens[, sensitivity_family_q_value := bh(p_value), by = .(host_feature_collection, covariate)]
write_tsv(cov_sens, file.path(tables_dir, "phase8b_host_covariate_sensitivity.tsv"))
trans_sens <- rbindlist(trans_rows, fill = TRUE)
trans_summary <- trans_sens[, .(primary_CLR_direction = primary_CLR_direction[representation == "primary_CLR"][1],
                                rCLR_direction = direction[representation == "rCLR"][1],
                                direction_agreement = direction_agreement[representation == "rCLR"][1],
                                coefficient_range = paste(range(coefficient, na.rm = TRUE), collapse = ";"),
                                FDR_range = "q_values_reported_in_primary_tables",
                                sample_exclusion_stability = "see_phase8b_sample_influence.tsv",
                                contamination_exclusion_stability = direction_agreement[representation == "contaminant_exclusion"][1],
                                TRANSFORMATION_SENSITIVE_MECHANISM = any(transformation_sensitive_label, na.rm = TRUE)),
                            by = .(taxon, host_feature, host_feature_collection)]
write_tsv(trans_summary, file.path(tables_dir, "phase8b_transformation_sensitivity.tsv"))
write_tsv(rbindlist(loo_rows, fill = TRUE), file.path(tables_dir, "phase8b_sample_influence.tsv"))

message("Running Moffitt50 exclusion sensitivity")
moffitt <- fread(file.path(root, "02_data", "reference", "PDAC_subtype_signatures", "Moffitt_50_gene_axis.tsv"))
moffitt_genes <- unique(moffitt[inclusion_status == "included", mapped_symbol])
expr_no_moffitt <- expr[!(rownames(expr) %in% moffitt_genes), ]
hallmark_no_moffitt <- long_to_score_dt(decoupleR::run_gsva(expr_no_moffitt, hallmark_net, method = "ssgsea", minsize = 15, verbose = FALSE), sample_ids)
progeny_no_moffitt <- data.table(patient_id = rownames(progeny::progeny(expr_no_moffitt, organism = "Human", top = 100, scale = TRUE, perm = 1, z_scores = FALSE)),
                                 as.data.table(progeny::progeny(expr_no_moffitt, organism = "Human", top = 100, scale = TRUE, perm = 1, z_scores = FALSE)))
tf_no_moffitt <- long_to_score_dt(decoupleR::run_viper(expr_no_moffitt, regulon[tf %in% reg_cov[retained == TRUE, tf]],
                                                       .source = tf, .target = target, .mor = mor,
                                                       .likelihood = NULL, minsize = 15, verbose = FALSE), sample_ids)
moff_rows <- list()
for (i in seq_len(nrow(candidates))) {
  cand <- candidates[i]
  orig <- score_lookup[[cand$host_feature_collection]]
  excl <- switch(cand$host_feature_collection, MSigDB_Hallmark = hallmark_no_moffitt, PROGENy = progeny_no_moffitt, DoRothEA_VIPER = tf_no_moffitt)
  if (!(cand$host_feature %in% names(excl))) next
  common <- intersect(orig$patient_id, excl$patient_id)
  r <- suppressWarnings(cor(orig[match(common, patient_id)][[cand$host_feature]], excl[match(common, patient_id)][[cand$host_feature]], use = "pairwise.complete.obs"))
  orig_res <- ols_hc3(orig[match(common, patient_id)][[cand$host_feature]], micro[cand$taxon, common])
  excl_res <- ols_hc3(excl[match(common, patient_id)][[cand$host_feature]], micro[cand$taxon, common])
  moff_rows[[length(moff_rows) + 1]] <- data.table(taxon = cand$taxon, host_feature = cand$host_feature,
    host_feature_collection = cand$host_feature_collection, score_correlation = as.numeric(r),
    original_coefficient = orig_res$coefficient, moffitt_excluded_coefficient = excl_res$coefficient,
    original_p_value = orig_res$p_value, moffitt_excluded_p_value = excl_res$p_value,
    direction_consistency = sign(orig_res$coefficient) == sign(excl_res$coefficient),
    change_in_statistical_support = ifelse(orig_res$p_value < 0.05 & excl_res$p_value >= 0.05, "lost_nominal_support", "no_nominal_support_loss"))
}
write_tsv(rbindlist(moff_rows, fill = TRUE), file.path(tables_dir, "phase8b_moffitt_gene_exclusion_sensitivity.tsv"))

message("Running WGCNA")
mad_vals <- matrixStats::rowMads(expr)
names(mad_vals) <- rownames(expr)
keep_genes <- names(sort(mad_vals, decreasing = TRUE))[seq_len(floor(length(mad_vals) * 0.25))]
wexpr <- t(expr[keep_genes, ])
powers <- 1:30
sft <- WGCNA::pickSoftThreshold(wexpr, powerVector = powers, networkType = "signed hybrid", verbose = 0)
soft <- as.data.table(sft$fitIndices)
write_tsv(soft, file.path(tables_dir, "phase8b_wgcna_soft_threshold.tsv"))
power_col <- grep("Power", names(soft), value = TRUE)[1]
r2_col <- grep("SFT.R.sq|Rsquared|R.sq", names(soft), value = TRUE)[1]
selected_power <- soft[get(r2_col) >= 0.85, get(power_col)][1]
if (is.na(selected_power)) selected_power <- soft[which.max(get(r2_col)), get(power_col)]
wgcna_start <- Sys.time()
net <- WGCNA::blockwiseModules(wexpr, power = selected_power, networkType = "signed hybrid",
                               TOMType = "signed", minModuleSize = 30, mergeCutHeight = 0.20,
                               numericLabels = FALSE, pamRespectsDendro = FALSE,
                               maxBlockSize = 5000, randomSeed = 2026, verbose = 0,
                               saveTOMs = FALSE)
wgcna_runtime <- as.numeric(difftime(Sys.time(), wgcna_start, units = "mins"))
saveRDS(net, file.path(models_dir, "phase8b_wgcna_blockwiseModules.rds"))
module_assign <- data.table(gene = names(net$colors), module = as.character(net$colors))
write_tsv(module_assign, file.path(tables_dir, "phase8b_wgcna_module_assignments.tsv.gz"))
mes <- as.data.table(net$MEs)
mes[, patient_id := rownames(net$MEs)]
setcolorder(mes, "patient_id")
write_tsv(mes, file.path(tables_dir, "phase8b_wgcna_module_eigengenes.tsv"))
mod_summary <- module_assign[, .(n_genes = .N), by = module]
mod_summary[, `:=`(selected_soft_power = selected_power, block_size = 5000,
                   number_of_blocks = length(unique(net$blockGenes)), threads = 2,
                   runtime_minutes = wgcna_runtime, modules_after_merging = uniqueN(module_assign$module),
                   grey_genes = module_assign[module == "grey", .N])]
write_tsv(mod_summary, file.path(tables_dir, "phase8b_wgcna_module_summary.tsv"))

wgcna_scores <- copy(mes)
wgcna_assoc <- model_feature_table(wgcna_scores, micro, primary_taxa, "WGCNA_Modules", evidence_vec, risk_vec, rclr_flags)
write_tsv(wgcna_assoc, file.path(tables_dir, "phase8b_wgcna_taxon_associations.tsv"))
supported_modules <- wgcna_assoc[bh_q_value < 0.05 | p_value < 0.05]
ann <- supported_modules[, .(taxon, module = host_feature, coefficient, p_value, bh_q_value)]
if (nrow(ann) > 0) {
  ann <- merge(ann, mod_summary, by = "module", all.x = TRUE)
  ann[, moffitt50_overlap_genes := vapply(module, function(m) paste(intersect(module_assign[module == m, gene], moffitt_genes), collapse = ";"), character(1))]
} else {
  ann <- data.table(taxon = character(), module = character(), coefficient = numeric(), p_value = numeric(), bh_q_value = numeric(),
                    n_genes = integer(), moffitt50_overlap_genes = character())
}
write_tsv(ann, file.path(tables_dir, "phase8b_wgcna_supported_module_annotations.tsv"))

message("Running genome-wide limma models")
gene_summary_rows <- list()
ranked_enrich_rows <- list()
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
reactome_list <- split(reactome$gene_symbol, reactome$gs_name)
kegg <- tryCatch(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"), error = function(e) data.table())
kegg_list <- if (nrow(kegg) > 0) split(kegg$gene_symbol, kegg$gs_name) else list()
run_limma_taxon <- function(taxon, x, suffix, covariate = NULL, expr_mat = expr) {
  design <- if (is.null(covariate)) model.matrix(~ zscore(x)) else model.matrix(~ zscore(x) + zscore(covariate))
  colnames(design)[2] <- "standardized_CLR_genus"
  fit <- eBayes(lmFit(expr_mat, design))
  coef_name <- "standardized_CLR_genus"
  eff <- fit$coefficients[, coef_name]
  tstat <- fit$t[, coef_name]
  pval <- fit$p.value[, coef_name]
  data.table(gene = rownames(expr_mat), taxon = taxon, model = suffix, effect_size = eff,
             moderated_standard_error = ifelse(is.finite(tstat) & tstat != 0, abs(eff / tstat), NA_real_),
             t_statistic = tstat, p_value = pval, q_value = p.adjust(pval, "BH"))
}
for (taxon in primary_taxa) {
  x <- as.numeric(micro[taxon, colnames(expr)])
  full <- run_limma_taxon(taxon, x, "primary_CLR")
  write_tsv(full, file.path(gene_full_dir, paste0("phase8b_host_gene_full_", taxon, ".tsv.gz")))
  for (model_name in c("purity", "immune", "stromal", "rCLR", "exclude_extreme")) {
    sens <- switch(model_name,
      purity = run_limma_taxon(taxon, x, "purity_adjusted", tme$inferred_tumor_purity),
      immune = run_limma_taxon(taxon, x, "immune_adjusted", tme$immune_score),
      stromal = run_limma_taxon(taxon, x, "stromal_adjusted", tme$stromal_score),
      rCLR = run_limma_taxon(taxon, as.numeric(rclr[taxon, colnames(expr)]), "rCLR"),
      exclude_extreme = {
        keep <- !(aligned$micro_sample_id %in% extreme_samples)
        run_limma_taxon(taxon, as.numeric(micro[taxon, colnames(expr)[keep]]), "exclude_extreme_samples", expr_mat = expr[, keep])
      })
    gene_summary_rows[[length(gene_summary_rows) + 1]] <- sens[, .(taxon = taxon[1], model = unique(model)[1],
      n_genes = .N, min_p_value = min(p_value, na.rm = TRUE), n_q_lt_0_05 = sum(q_value < 0.05, na.rm = TRUE),
      top_gene = gene[which.min(p_value)][1], top_effect_size = effect_size[which.min(p_value)][1])]
  }
  gene_summary_rows[[length(gene_summary_rows) + 1]] <- full[, .(taxon = taxon[1], model = "primary_CLR",
    n_genes = .N, min_p_value = min(p_value, na.rm = TRUE), n_q_lt_0_05 = sum(q_value < 0.05, na.rm = TRUE),
    top_gene = gene[which.min(p_value)][1], top_effect_size = effect_size[which.min(p_value)][1])]
  ranks <- full$t_statistic; names(ranks) <- full$gene; ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  for (coll in c("Hallmark", "Reactome", "KEGG")) {
    pathways <- switch(coll, Hallmark = hallmark_list, Reactome = reactome_list, KEGG = kegg_list)
    if (length(pathways) == 0) next
    fg <- suppressWarnings(fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, nperm = 10000))
    ranked_enrich_rows[[length(ranked_enrich_rows) + 1]] <- as.data.table(fg)[, .(taxon = taxon, collection = coll,
      pathway, pval, padj, ES, NES, size, leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"),
      gene_set_version = ifelse(coll == "Hallmark", unique(hallmark$db_version)[1], ifelse(coll == "Reactome", unique(reactome$db_version)[1], ifelse(nrow(kegg) > 0, unique(kegg$db_version)[1], NA_character_))))]
  }
}
write_tsv(rbindlist(gene_summary_rows, fill = TRUE), file.path(tables_dir, "phase8b_host_gene_associations_summary.tsv"))
write_tsv(rbindlist(ranked_enrich_rows, fill = TRUE), file.path(tables_dir, "phase8b_ranked_gene_enrichment.tsv"))

message("Shared-mechanism summaries and evidence classification")
primary_mechs <- rbind(pathway_assoc[, .(taxon, host_feature, host_feature_collection, coefficient, p_value, bh_q_value)],
                       tf_assoc[, .(taxon, host_feature, host_feature_collection, coefficient, p_value, bh_q_value)],
                       wgcna_assoc[, .(taxon, host_feature, host_feature_collection, coefficient, p_value, bh_q_value)], fill = TRUE)
primary_mechs[, supported_nominal := p_value < 0.05]
sign_mat <- dcast(primary_mechs[supported_nominal == TRUE], host_feature_collection + host_feature ~ taxon, value.var = "coefficient", fun.aggregate = function(x) sign_chr(x[1]))
write_tsv(sign_mat, file.path(tables_dir, "phase8b_cross_taxon_sign_consistency_matrix.tsv"))
tax_cor <- as.data.table(cor(t(micro[primary_taxa, ]), use = "pairwise.complete.obs"), keep.rownames = "taxon")
write_tsv(tax_cor, file.path(tables_dir, "phase8b_taxon_correlation.tsv"))
sharing <- primary_mechs[supported_nominal == TRUE, .(n_taxa_supported = uniqueN(taxon),
                                                      taxa = paste(unique(taxon), collapse = ";"),
                                                      directions = paste(unique(sign_chr(coefficient)), collapse = ";")),
                         by = .(host_feature_collection, host_feature)]
sharing[, interpretation := fifelse(grepl(";", directions), "shared_with_inconsistent_direction",
                                    fifelse(n_taxa_supported > 1, "shared_consistent_direction_check_taxon_correlation", "taxon_specific"))]
write_tsv(sharing, file.path(tables_dir, "phase8b_shared_mechanism_summary.tsv"))

evidence <- merge(primary_mechs, trans_summary[, .(taxon, host_feature, host_feature_collection, TRANSFORMATION_SENSITIVE_MECHANISM)],
                  by = c("taxon", "host_feature", "host_feature_collection"), all.x = TRUE)
evidence <- merge(evidence, cov_sens[, .(composition_sensitive = any(p_value >= 0.05 | sign_change, na.rm = TRUE)),
                                     by = .(taxon, host_feature, host_feature_collection)],
                  by = c("taxon", "host_feature", "host_feature_collection"), all.x = TRUE)
evidence <- merge(evidence, rbindlist(loo_rows, fill = TRUE)[, .(sample_sensitive = !LOO_sign_stable),
                                                             by = .(taxon, host_feature, host_feature_collection)],
                  by = c("taxon", "host_feature", "host_feature_collection"), all.x = TRUE)
evidence[is.na(TRANSFORMATION_SENSITIVE_MECHANISM), TRANSFORMATION_SENSITIVE_MECHANISM := taxon %in% rclr_flags & bh_q_value < 0.05]
evidence[is.na(composition_sensitive), composition_sensitive := FALSE]
evidence[is.na(sample_sensitive), sample_sensitive := FALSE]
evidence[, evidence_category := fifelse(bh_q_value < 0.05 & TRANSFORMATION_SENSITIVE_MECHANISM, "TRANSFORMATION_SENSITIVE_MECHANISM",
                                 fifelse(bh_q_value < 0.05 & composition_sensitive, "COMPOSITION_SENSITIVE_MECHANISM",
                                 fifelse(bh_q_value < 0.05 & sample_sensitive, "SAMPLE_SENSITIVE_MECHANISM",
                                 fifelse(bh_q_value < 0.05, "ROBUST_HOST_MECHANISM",
                                 fifelse(p_value < 0.05, "EXPLORATORY_HOST_MECHANISM", "NO_SUPPORTED_MECHANISM")))))]
evidence[, criteria := paste0("primary_p=", signif(p_value, 4), "; primary_q=", signif(bh_q_value, 4),
                              "; transformation_sensitive=", TRANSFORMATION_SENSITIVE_MECHANISM,
                              "; composition_sensitive=", composition_sensitive,
                              "; sample_sensitive=", sample_sensitive)]
write_tsv(evidence, file.path(tables_dir, "phase8b_host_mechanism_evidence.tsv"))

plot_assoc <- function(dt, path, title) {
  dt <- as.data.table(dt)
  setnames(dt, make.unique(names(dt)))
  if (nrow(dt) == 0) {
    p <- ggplot() + theme_void() + labs(title = title)
    ggsave(path, p, width = 10, height = 8)
    return(invisible(NULL))
  }
  top <- dt[order(p_value)][1:min(.N, 60)]
  p <- ggplot(top, aes(x = reorder(paste(taxon, host_feature, sep = " | "), coefficient), y = coefficient, fill = bh_q_value < 0.05)) +
    geom_col() + coord_flip() + theme_bw(base_size = 8) + labs(x = NULL, y = "OLS HC3 coefficient", title = title, fill = "q < 0.05")
  ggsave(path, p, width = 10, height = 8)
}
plot_assoc(hallmark_assoc, file.path(figures_dir, "phase8b_taxon_hallmark_associations.pdf"), "Phase 8B Hallmark Associations")
plot_assoc(progeny_assoc, file.path(figures_dir, "phase8b_taxon_progeny_associations.pdf"), "Phase 8B PROGENy Associations")
plot_assoc(tf_assoc, file.path(figures_dir, "phase8b_taxon_tf_associations.pdf"), "Phase 8B TF Associations")
plot_assoc(cov_sens[, .(taxon, host_feature, coefficient = genus_coefficient, p_value, bh_q_value = sensitivity_family_q_value)], file.path(figures_dir, "phase8b_covariate_sensitivity.pdf"), "Phase 8B Covariate Sensitivity")
plot_assoc(trans_sens[representation == "rCLR", .(taxon, host_feature, coefficient, p_value, bh_q_value = p_value)], file.path(figures_dir, "phase8b_rclr_direction_sensitivity.pdf"), "Phase 8B rCLR Sensitivity")
ggsave(file.path(figures_dir, "phase8b_wgcna_soft_threshold.pdf"),
       ggplot(soft, aes(x = get(power_col), y = get(r2_col))) + geom_line() + geom_point() + geom_hline(yintercept = 0.85, linetype = 2) + theme_bw() + labs(x = "Soft threshold power", y = "Scale-free topology fit", title = "WGCNA Soft Threshold"),
       width = 7, height = 5)
plot_assoc(wgcna_assoc, file.path(figures_dir, "phase8b_wgcna_module_taxon_heatmap.pdf"), "Phase 8B WGCNA Module Associations")
plot_assoc(ann[, .(taxon, host_feature = module, coefficient, p_value, bh_q_value)], file.path(figures_dir, "phase8b_supported_module_annotations.pdf"), "Supported Module Annotations")
enr <- fread(file.path(tables_dir, "phase8b_ranked_gene_enrichment.tsv"))
plot_assoc(enr[, .(taxon, host_feature = pathway, coefficient = NES, p_value = pval, bh_q_value = padj)], file.path(figures_dir, "phase8b_host_gene_enrichment.pdf"), "Ranked Gene Enrichment")
plot_assoc(sharing[, .(taxon = taxa, host_feature, coefficient = n_taxa_supported, p_value = 1 / pmax(n_taxa_supported, 1), bh_q_value = 1)], file.path(figures_dir, "phase8b_shared_mechanism_network.pdf"), "Shared Mechanism Counts")
plot_assoc(evidence[, .(host_feature = evidence_category, coefficient = .N, p_value = 1 / .N, bh_q_value = 1), by = .(taxon, evidence_category)], file.path(figures_dir, "phase8b_mechanism_evidence_summary.pdf"), "Mechanism Evidence Summary")

message("Phase 8B R execution complete")
