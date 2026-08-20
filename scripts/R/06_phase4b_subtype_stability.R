#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cluster)
  library(clue)
  library(mclust)
  library(ggplot2)
})

ROOT <- normalizePath(getwd())
set.seed(2026)

dir.create("05_results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("05_results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("05_results/models/phase4b", recursive = TRUE, showWarnings = FALSE)

K_VALUES <- 2:6
OUTLIERS <- c("YX16135T", "YX16158T", "YX16194T", "YX16224T")
ITERATIONS <- 1000L
SEED <- 2026L
PAC_LOWER <- 0.1
PAC_UPPER <- 0.9
MIN_CLUSTER_SIZE <- 5L

tbl_path <- function(x) file.path("05_results/tables", x)
fig_path <- function(x) file.path("05_results/figures", x)
model_path <- function(x) file.path("05_results/models/phase4b", x)

read_matrix <- function(path) {
  df <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  genes <- df[[1]]
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  mat
}

row_scale <- function(mat, median_center = FALSE) {
  if (median_center) {
    mat <- mat - apply(mat, 1, median, na.rm = TRUE)
  }
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

gene_median_impute <- function(mat) {
  for (i in seq_len(nrow(mat))) {
    miss <- is.na(mat[i, ])
    if (any(miss)) mat[i, miss] <- median(mat[i, ], na.rm = TRUE)
  }
  mat
}

pearson_dist <- function(mat) {
  cm <- suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = "pearson"))
  cm[is.na(cm)] <- 0
  cm <- pmax(pmin(cm, 1), -1)
  diag(cm) <- 1
  as.dist(1 - cm)
}

euclidean_dist <- function(mat) dist(t(mat), method = "euclidean")

cluster_partition <- function(mat, distance, method, k) {
  d <- if (distance == "pearson") pearson_dist(mat) else euclidean_dist(mat)
  cutree(hclust(d, method = method), k = k)
}

distance_matrix <- function(mat, distance) {
  if (distance == "pearson") as.matrix(pearson_dist(mat)) else as.matrix(euclidean_dist(mat))
}

align_to_reference <- function(ref_labels, run_labels, selected_samples, k) {
  ref <- ref_labels[selected_samples]
  run <- run_labels[selected_samples]
  overlap <- matrix(0, nrow = k, ncol = k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) overlap[i, j] <- sum(run == i & ref == j)
  }
  cost <- max(overlap) - overlap
  assignment <- solve_LSAP(cost)
  aligned <- integer(length(run))
  for (i in seq_len(k)) aligned[run == i] <- assignment[i]
  names(aligned) <- selected_samples
  aligned
}

best_label_crosswalk <- function(clusters, labels) {
  cluster_levels <- sort(unique(clusters))
  label_levels <- sort(unique(labels[!is.na(labels)]))
  n <- max(length(cluster_levels), length(label_levels))
  padded_labels <- c(label_levels, paste0("UNMAPPED_", seq_len(n - length(label_levels))))
  tab <- matrix(0, nrow = n, ncol = n)
  for (i in seq_along(cluster_levels)) {
    for (j in seq_along(label_levels)) tab[i, j] <- sum(clusters == cluster_levels[i] & labels == label_levels[j])
  }
  cost <- max(tab, na.rm = TRUE) - tab
  cost[is.na(cost)] <- max(cost, na.rm = TRUE)
  assignment <- solve_LSAP(cost)
  data.frame(
    cluster_id = cluster_levels,
    aligned_public_label = padded_labels[as.integer(assignment[seq_along(cluster_levels)])],
    overlap_n = tab[cbind(seq_along(cluster_levels), as.integer(assignment[seq_along(cluster_levels)]))],
    stringsAsFactors = FALSE
  )
}

entropy_bits <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log2(p))
}

normalized_mutual_information <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  px <- rowSums(tab) / n
  py <- colSums(tab) / n
  pxy <- tab / n
  mi <- 0
  for (i in seq_len(nrow(pxy))) {
    for (j in seq_len(ncol(pxy))) {
      if (pxy[i, j] > 0) mi <- mi + pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
    }
  }
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

cohens_kappa <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / (n * n)
  if (isTRUE(all.equal(1, pe))) return(NA_real_)
  (po - pe) / (1 - pe)
}

cluster_balance <- function(sizes) min(sizes) / max(sizes)

consensus_cdf_area <- function(vals) {
  x <- sort(vals)
  if (length(x) < 2) return(NA_real_)
  y <- seq_along(x) / length(x)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

prediction_strength_once <- function(mat, distance, method, k, train_idx, test_idx) {
  train_mat <- mat[, train_idx, drop = FALSE]
  test_mat <- mat[, test_idx, drop = FALSE]
  train_labels <- cluster_partition(train_mat, distance, method, k)
  test_labels <- cluster_partition(test_mat, distance, method, k)
  centroids <- sapply(seq_len(k), function(cl) rowMeans(train_mat[, train_labels == cl, drop = FALSE]))
  if (is.null(dim(centroids))) centroids <- matrix(centroids, ncol = 1)
  pred <- apply(test_mat, 2, function(x) which.min(colSums((centroids - x)^2)))
  strengths <- c()
  for (cl in sort(unique(test_labels))) {
    members <- names(test_labels)[test_labels == cl]
    if (length(members) < 2) next
    pairs <- combn(members, 2)
    strengths <- c(strengths, mean(pred[pairs[1, ]] == pred[pairs[2, ]]))
  }
  if (length(strengths) == 0) NA_real_ else min(strengths)
}

make_heatmap_pdf <- function(mat, labels, file, title) {
  ord <- order(labels)
  df <- as.data.frame(as.table(mat[ord, ord]))
  names(df) <- c("sample_i", "sample_j", "consensus")
  pdf(file, width = 7.5, height = 6.8)
  print(
    ggplot(df, aes(sample_i, sample_j, fill = consensus)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0.5, limits = c(0, 1)) +
      labs(title = title, x = NULL, y = NULL, fill = "Consensus") +
      theme_minimal(base_size = 9) +
      theme(axis.text = element_blank(), panel.grid = element_blank())
  )
  dev.off()
}

analysis_specs <- read.delim("01_metadata/subtype_stability_parameter_inventory.tsv", stringsAsFactors = FALSE)
analysis_specs <- analysis_specs[match(c(
  "STAB_CSY_PRIMARY", "STAB_CSY_LOG2", "STAB_UNSUP_HVG", "STAB_CSY_OUTLIER_EXCL",
  "STAB_HVG_OUTLIER_EXCL", "STAB_CSY_FEAT_RESAMP", "STAB_CSY_IMPUTED", "STAB_HVG_VAR_FILTER"
), analysis_specs$analysis_id), ]

primary_assign <- read.delim(tbl_path("phase3b_primary_subtype_assignments.tsv"), stringsAsFactors = FALSE)
public_labels <- setNames(primary_assign$original_public_subtype, primary_assign$expression_sample_id)
patient_ids <- setNames(primary_assign$patient_id, primary_assign$expression_sample_id)

all_methods <- read.delim(tbl_path("phase3b_all_method_assignments.tsv"), stringsAsFactors = FALSE)
moffitt_rows <- all_methods[all_methods$method_name == "Moffitt", ]
purist_rows <- all_methods[all_methods$method_name == "PurIST", ]
patient_to_sample <- setNames(primary_assign$expression_sample_id, primary_assign$patient_id)
purist_prob <- setNames(as.numeric(purist_rows$probability_or_confidence), patient_to_sample[purist_rows$patient_id])
moffitt_basal <- setNames(as.numeric(moffitt_rows$basal_score), patient_to_sample[moffitt_rows$patient_id])
moffitt_classical <- setNames(as.numeric(moffitt_rows$classical_score), patient_to_sample[moffitt_rows$patient_id])

norm_mat <- read_matrix("03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz")
log2_mat <- read_matrix("03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz")

csy_sig <- read.delim("02_data/reference/PDAC_subtype_signatures/GSE172356_original_signatures.tsv", stringsAsFactors = FALSE)
csy_genes <- unique(csy_sig$mapped_symbol[csy_sig$presence_in_GSE172356 == "True" & !is.na(csy_sig$mapped_symbol) & csy_sig$mapped_symbol != "NA"])

if (length(csy_genes) != 94L) stop("CSY signature gene count is not 94 after verified mapping.")

prepare_analysis_matrix <- function(spec) {
  mat <- if (grepl("log2", spec$transformation, fixed = TRUE)) log2_mat else norm_mat
  if (spec$analysis_id %in% c("STAB_CSY_OUTLIER_EXCL", "STAB_HVG_OUTLIER_EXCL")) {
    mat <- mat[, setdiff(colnames(mat), OUTLIERS), drop = FALSE]
  }
  if (spec$analysis_id == "STAB_CSY_IMPUTED") mat <- gene_median_impute(mat)
  if (spec$gene_set == "CSY_94_gene_signature") {
    mat <- mat[intersect(csy_genes, rownames(mat)), , drop = FALSE]
    mat <- row_scale(mat, median_center = TRUE)
  } else {
    if (spec$gene_set == "HVG_500") {
      detected <- rowSums(norm_mat[, colnames(mat), drop = FALSE] >= 10, na.rm = TRUE)
      keep <- names(detected)[detected >= ceiling(0.20 * ncol(mat))]
      mat <- mat[intersect(rownames(mat), keep), , drop = FALSE]
      topn <- 500L
    } else {
      topn <- 1000L
    }
    madv <- apply(mat, 1, mad, na.rm = TRUE)
    genes <- names(sort(madv, decreasing = TRUE))[seq_len(min(topn, length(madv)))]
    mat <- row_scale(mat[genes, , drop = FALSE], median_center = FALSE)
  }
  mat[is.na(mat)] <- 0
  mat
}

metric_rows <- list()
size_rows <- list()
sample_rows <- list()
prob_rows <- list()
comparison_rows <- list()
crosswalk_rows <- list()
all_final <- list()
consensus_store <- list()
runtime_rows <- list()

start_all <- Sys.time()

for (a in seq_len(nrow(analysis_specs))) {
  spec <- analysis_specs[a, ]
  analysis_id <- spec$analysis_id
  message("Running ", analysis_id)
  analysis_start <- Sys.time()
  set.seed(SEED + a)
  mat <- prepare_analysis_matrix(spec)
  samples <- colnames(mat)
  n <- length(samples)
  n_sub <- floor(as.numeric(spec$resampling_fraction) * n)
  distance <- spec$distance
  method <- if (spec$clustering_method == "hierarchical_average") "average" else "ward.D2"
  feat_frac <- if (spec$feature_resampling == "subsample_80_percent") 0.8 else 1.0
  dmat <- distance_matrix(mat, distance)

  ref_labels_by_k <- list()
  sil_by_k <- list()
  for (k in K_VALUES) {
    ref_labels_by_k[[as.character(k)]] <- cluster_partition(mat, distance, method, k)
    silv <- silhouette(ref_labels_by_k[[as.character(k)]], as.dist(dmat))[, "sil_width"]
    names(silv) <- names(ref_labels_by_k[[as.character(k)]])
    sil_by_k[[as.character(k)]] <- silv
  }

  consensus_num <- lapply(K_VALUES, function(k) matrix(0, n, n, dimnames = list(samples, samples)))
  names(consensus_num) <- as.character(K_VALUES)
  consensus_den <- matrix(0, n, n, dimnames = list(samples, samples))
  assign_counts <- lapply(K_VALUES, function(k) matrix(0, n, k, dimnames = list(samples, paste0("cluster_", seq_len(k)))))
  names(assign_counts) <- as.character(K_VALUES)
  present_counts <- setNames(integer(n), samples)
  same_ref_counts <- lapply(K_VALUES, function(k) setNames(integer(n), samples))
  names(same_ref_counts) <- as.character(K_VALUES)
  run_partitions <- lapply(K_VALUES, function(k) vector("list", ITERATIONS))
  names(run_partitions) <- as.character(K_VALUES)

  for (b in seq_len(ITERATIONS)) {
    sel <- sort(sample(samples, n_sub, replace = FALSE))
    present_counts[sel] <- present_counts[sel] + 1L
    consensus_den[sel, sel] <- consensus_den[sel, sel] + 1L
    if (feat_frac < 1) {
      feats <- sort(sample(rownames(mat), max(2L, floor(feat_frac * nrow(mat))), replace = FALSE))
      run_mat <- mat[feats, sel, drop = FALSE]
    } else {
      run_mat <- mat[, sel, drop = FALSE]
    }
    hc <- hclust(if (distance == "pearson") pearson_dist(run_mat) else euclidean_dist(run_mat), method = method)
    for (k in K_VALUES) {
      run_lab <- cutree(hc, k = k)
      names(run_lab) <- sel
      run_partitions[[as.character(k)]][[b]] <- run_lab
      aligned <- align_to_reference(ref_labels_by_k[[as.character(k)]], run_lab, sel, k)
      for (cl in seq_len(k)) assign_counts[[as.character(k)]][sel, cl] <- assign_counts[[as.character(k)]][sel, cl] + as.integer(aligned == cl)
      same_ref_counts[[as.character(k)]][sel] <- same_ref_counts[[as.character(k)]][sel] + as.integer(aligned == ref_labels_by_k[[as.character(k)]][sel])
      for (cl in sort(unique(run_lab))) {
        members <- names(run_lab)[run_lab == cl]
        consensus_num[[as.character(k)]][members, members] <- consensus_num[[as.character(k)]][members, members] + 1L
      }
    }
  }

  consensus_by_k <- lapply(K_VALUES, function(k) {
    m <- consensus_num[[as.character(k)]] / pmax(consensus_den, 1)
    diag(m) <- 1
    m
  })
  names(consensus_by_k) <- as.character(K_VALUES)

  jaccard_by_k <- lapply(K_VALUES, function(k) matrix(NA_real_, ITERATIONS, k))
  names(jaccard_by_k) <- as.character(K_VALUES)
  pred_strength <- lapply(K_VALUES, function(k) numeric(ITERATIONS))
  names(pred_strength) <- as.character(K_VALUES)

  for (b in seq_len(ITERATIONS)) {
    boot_idx <- sample(samples, n, replace = TRUE)
    unique_boot <- unique(boot_idx)
    if (length(unique_boot) < max(K_VALUES)) next
    boot_mat <- mat[, unique_boot, drop = FALSE]
    hc_boot <- hclust(if (distance == "pearson") pearson_dist(boot_mat) else euclidean_dist(boot_mat), method = method)
    split_idx <- sample(samples, n_sub, replace = FALSE)
    test_idx <- setdiff(samples, split_idx)
    for (k in K_VALUES) {
      boot_lab <- cutree(hc_boot, k = k)
      names(boot_lab) <- unique_boot
      ref_lab <- ref_labels_by_k[[as.character(k)]]
      for (cl in seq_len(k)) {
        ref_members <- names(ref_lab)[ref_lab == cl]
        vals <- sapply(seq_len(k), function(j) {
          boot_members <- names(boot_lab)[boot_lab == j]
          length(intersect(ref_members, boot_members)) / length(union(ref_members, boot_members))
        })
        jaccard_by_k[[as.character(k)]][b, cl] <- max(vals)
      }
      pred_strength[[as.character(k)]][b] <- prediction_strength_once(mat, distance, method, k, split_idx, test_idx)
    }
  }

  for (k in K_VALUES) {
    k_chr <- as.character(k)
    ref_lab <- ref_labels_by_k[[k_chr]]
    all_final[[paste(analysis_id, k_chr, sep = "__")]] <- ref_lab
    consensus_store[[paste(analysis_id, k_chr, sep = "__")]] <- consensus_by_k[[k_chr]]
    cm <- consensus_by_k[[k_chr]]
    offdiag <- cm[upper.tri(cm)]
    pac <- mean(offdiag >= PAC_LOWER & offdiag <= PAC_UPPER)
    cdf_area <- consensus_cdf_area(offdiag)
    delta_area <- NA_real_
    if (k > min(K_VALUES)) {
      prev <- metric_rows[[length(metric_rows)]]
      if (!is.null(prev) && prev$analysis_id == analysis_id && prev$candidate_K == k - 1) delta_area <- cdf_area - prev$consensus_cdf_area
    }
    sizes <- as.integer(table(factor(ref_lab, levels = seq_len(k))))
    sil <- sil_by_k[[k_chr]]
    jmeans <- colMeans(jaccard_by_k[[k_chr]], na.rm = TRUE)
    probs <- assign_counts[[k_chr]] / pmax(as.numeric(present_counts), 1)
    ent <- apply(probs, 1, entropy_bits)
    bfreq <- same_ref_counts[[k_chr]] / pmax(as.numeric(present_counts), 1)
    ari_vals <- c()
    set.seed(SEED + a * 100 + k)
    pairs <- replicate(1000L, sample(seq_len(ITERATIONS), 2), simplify = FALSE)
    for (pair in pairs) {
      p1 <- run_partitions[[k_chr]][[pair[1]]]
      p2 <- run_partitions[[k_chr]][[pair[2]]]
      common <- intersect(names(p1), names(p2))
      if (length(common) > k) ari_vals <- c(ari_vals, adjustedRandIndex(p1[common], p2[common]))
    }
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      candidate_K = k,
      consensus_cdf_area = cdf_area,
      delta_area_under_cdf = delta_area,
      PAC = pac,
      overall_mean_silhouette = mean(sil),
      mean_Jaccard_stability = mean(jmeans, na.rm = TRUE),
      min_cluster_Jaccard_stability = min(jmeans, na.rm = TRUE),
      prediction_strength = mean(unlist(pred_strength[[k_chr]]), na.rm = TRUE),
      mean_ARI_across_resampling_runs = mean(ari_vals, na.rm = TRUE),
      min_cluster_size = min(sizes),
      max_cluster_size = max(sizes),
      cluster_size_balance = cluster_balance(sizes),
      empty_cluster_occurrence = 0,
      very_small_cluster_occurrence = as.integer(any(sizes < MIN_CLUSTER_SIZE)),
      mean_assignment_entropy = mean(ent, na.rm = TRUE),
      sample_n = n,
      gene_n = nrow(mat),
      stringsAsFactors = FALSE
    )
    for (cl in seq_len(k)) {
      members <- names(ref_lab)[ref_lab == cl]
      ik <- if (length(members) > 1) mean(cm[members, members][upper.tri(cm[members, members])]) else NA_real_
      size_rows[[length(size_rows) + 1L]] <- data.frame(
        analysis_id = analysis_id,
        candidate_K = k,
        cluster_id = cl,
        cluster_size = length(members),
        within_cluster_consensus = ik,
        cluster_mean_silhouette = mean(sil[members], na.rm = TRUE),
        bootstrap_jaccard_stability = jmeans[cl],
        cluster_size_flag = ifelse(length(members) < MIN_CLUSTER_SIZE, "VERY_SMALL", "OK"),
        stringsAsFactors = FALSE
      )
    }
    alt_dist <- sapply(samples, function(s) {
      centroids <- sapply(seq_len(k), function(cl) rowMeans(mat[, ref_lab == cl, drop = FALSE]))
      if (is.null(dim(centroids))) centroids <- matrix(centroids, ncol = 1)
      d <- colSums((centroids - mat[, s])^2)^0.5
      paste(paste0("cluster_", seq_len(k), "=", sprintf("%.6f", d)), collapse = ";")
    })
    for (s in samples) {
      cl <- ref_lab[s]
      members <- names(ref_lab)[ref_lab == cl]
      item_cons <- if (length(members) > 1) mean(cm[s, setdiff(members, s)], na.rm = TRUE) else NA_real_
      cluster_probs <- probs[s, ]
      centroids <- sapply(seq_len(k), function(x) rowMeans(mat[, ref_lab == x, drop = FALSE]))
      if (is.null(dim(centroids))) centroids <- matrix(centroids, ncol = 1)
      dist_vec <- colSums((centroids - mat[, s])^2)^0.5
      sample_rows[[length(sample_rows) + 1L]] <- data.frame(
        analysis_id = analysis_id,
        candidate_K = k,
        sample_id = s,
        patient_id = patient_ids[s],
        public_subtype = public_labels[s],
        final_cluster_assignment = cl,
        item_consensus = item_cons,
        bootstrap_assignment_frequency = bfreq[s],
        assignment_entropy = ent[s],
        silhouette_width = sil[s],
        co_clustering_probability = mean(cm[s, setdiff(samples, s)], na.rm = TRUE),
        distance_to_assigned_centroid = dist_vec[cl],
        distance_to_alternative_centroids = alt_dist[s],
        assignment_confidence = max(cluster_probs, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      for (cl2 in seq_len(k)) {
        prob_rows[[length(prob_rows) + 1L]] <- data.frame(
          analysis_id = analysis_id,
          candidate_K = k,
          sample_id = s,
          cluster_id = cl2,
          assignment_probability = cluster_probs[cl2],
          stringsAsFactors = FALSE
        )
      }
    }

    cross <- best_label_crosswalk(ref_lab, public_labels[names(ref_lab)])
    cross$analysis_id <- analysis_id
    cross$candidate_K <- k
    crosswalk_rows[[length(crosswalk_rows) + 1L]] <- cross[, c("analysis_id", "candidate_K", "cluster_id", "aligned_public_label", "overlap_n")]
    aligned_labels <- setNames(cross$aligned_public_label[match(ref_lab, cross$cluster_id)], names(ref_lab))
    tab <- table(public_labels[names(ref_lab)], aligned_labels)
    per_class <- sapply(rownames(tab), function(lbl) {
      if (sum(tab[lbl, ]) == 0 || !(lbl %in% colnames(tab))) NA_real_ else tab[lbl, lbl] / sum(tab[lbl, ])
    })
    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      candidate_K = k,
      adjusted_rand_index = adjustedRandIndex(public_labels[names(ref_lab)], aligned_labels),
      normalized_mutual_information = normalized_mutual_information(public_labels[names(ref_lab)], aligned_labels),
      cohens_kappa = cohens_kappa(public_labels[names(ref_lab)], aligned_labels),
      confusion_matrix = paste(capture.output(print(tab)), collapse = " | "),
      per_class_agreement = paste(paste(names(per_class), sprintf("%.6f", per_class), sep = "="), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  saveRDS(list(consensus = consensus_by_k, final_labels = ref_labels_by_k), model_path(paste0(analysis_id, "_consensus_reference.rds")))
  runtime_rows[[length(runtime_rows) + 1L]] <- data.frame(
    analysis_id = analysis_id,
    sample_n = n,
    gene_n = nrow(mat),
    seed = SEED,
    iterations = ITERATIONS,
    runtime_seconds = as.numeric(difftime(Sys.time(), analysis_start, units = "secs")),
    stringsAsFactors = FALSE
  )
}

cluster_metrics <- do.call(rbind, metric_rows)
cluster_sizes <- do.call(rbind, size_rows)
sample_stability <- do.call(rbind, sample_rows)
assignment_probs <- do.call(rbind, prob_rows)
public_compare <- do.call(rbind, comparison_rows)
crosswalk <- do.call(rbind, crosswalk_rows)
runtime <- do.call(rbind, runtime_rows)

preferred <- do.call(rbind, lapply(split(cluster_metrics, cluster_metrics$analysis_id), function(df) {
  viable <- df[df$min_cluster_size >= MIN_CLUSTER_SIZE, ]
  if (nrow(viable) == 0) viable <- df
  viable$rank_pac <- rank(viable$PAC, ties.method = "min")
  viable$rank_sil <- rank(-viable$overall_mean_silhouette, ties.method = "min")
  viable$rank_jac <- rank(-viable$mean_Jaccard_stability, ties.method = "min")
  viable$rank_pred <- rank(-viable$prediction_strength, ties.method = "min")
  viable$multi_metric_rank_sum <- viable$rank_pac + viable$rank_sil + viable$rank_jac + viable$rank_pred
  best <- viable[order(viable$multi_metric_rank_sum, viable$PAC, -viable$overall_mean_silhouette), ][1, ]
  data.frame(
    analysis_id = best$analysis_id,
    preferred_K = best$candidate_K,
    preferred_K_basis = "lowest prespecified multi-metric rank across PAC, silhouette, Jaccard, and prediction strength among viable cluster sizes",
    primary_K3_PAC = df$PAC[df$candidate_K == 3],
    primary_K3_silhouette = df$overall_mean_silhouette[df$candidate_K == 3],
    primary_K3_mean_jaccard = df$mean_Jaccard_stability[df$candidate_K == 3],
    stringsAsFactors = FALSE
  )
}))

write.table(cluster_metrics, tbl_path("phase4b_cluster_stability_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cluster_sizes, tbl_path("phase4b_cluster_size_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(preferred, tbl_path("phase4b_k_selection_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_stability, tbl_path("phase4b_sample_stability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(assignment_probs, tbl_path("phase4b_sample_assignment_probabilities.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(public_compare, tbl_path("phase4b_public_label_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(crosswalk, tbl_path("phase4b_cluster_label_crosswalk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(runtime, tbl_path("phase4b_runtime_versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

k3_sample <- sample_stability[sample_stability$candidate_K == 3, ]
axis_df <- data.frame(
  sample_id = names(public_labels),
  public_subtype = public_labels,
  basal_score = moffitt_basal[names(public_labels)],
  classical_score = moffitt_classical[names(public_labels)],
  moffitt_score_difference = moffitt_basal[names(public_labels)] - moffitt_classical[names(public_labels)],
  purist_basal_probability = purist_prob[names(public_labels)],
  stringsAsFactors = FALSE
)

primary_k3 <- k3_sample[k3_sample$analysis_id == "STAB_CSY_PRIMARY", ]
primary_cross <- crosswalk[crosswalk$analysis_id == "STAB_CSY_PRIMARY" & crosswalk$candidate_K == 3, ]
basal_cluster <- primary_cross$cluster_id[primary_cross$aligned_public_label == "Basal"][1]
classical_cluster <- primary_cross$cluster_id[primary_cross$aligned_public_label == "Classical"][1]
primary_mat <- prepare_analysis_matrix(analysis_specs[analysis_specs$analysis_id == "STAB_CSY_PRIMARY", ])
primary_labels <- all_final[["STAB_CSY_PRIMARY__3"]]
centroids <- sapply(seq_len(3), function(cl) rowMeans(primary_mat[, primary_labels == cl, drop = FALSE]))
dist_basal <- setNames(rep(NA_real_, ncol(primary_mat)), colnames(primary_mat))
dist_classical <- dist_basal
if (!is.na(basal_cluster)) dist_basal <- setNames(colSums((centroids[, basal_cluster] - primary_mat)^2)^0.5, colnames(primary_mat))
if (!is.na(classical_cluster)) dist_classical <- setNames(colSums((centroids[, classical_cluster] - primary_mat)^2)^0.5, colnames(primary_mat))

hybrid_rows <- merge(k3_sample, axis_df, by = c("sample_id", "public_subtype"), all.x = TRUE)
hybrid_rows$distance_to_basal_centroid <- dist_basal[hybrid_rows$sample_id]
hybrid_rows$distance_to_classical_centroid <- dist_classical[hybrid_rows$sample_id]
hybrid_rows$interpretation_category <- "TO_VERIFY"
hybrid_rows$interpretation_category[hybrid_rows$assignment_entropy >= 0.5] <- "HETEROGENEOUS_OR_UNSTABLE"
hybrid_rows$interpretation_category[hybrid_rows$public_subtype == "Hybrid" & hybrid_rows$assignment_entropy < 0.3 & hybrid_rows$item_consensus >= 0.80 & hybrid_rows$silhouette_width > 0] <- "STABLE_HYBRID"
hybrid_rows$interpretation_category[hybrid_rows$public_subtype == "Hybrid" & hybrid_rows$purist_basal_probability >= 0.3 & hybrid_rows$purist_basal_probability <= 0.7 & hybrid_rows$assignment_entropy >= 0.3] <- "INTERMEDIATE_STATE"
hybrid_rows$interpretation_category[hybrid_rows$public_subtype == "Basal" & hybrid_rows$assignment_entropy < 0.2 & hybrid_rows$item_consensus >= 0.85 & hybrid_rows$silhouette_width > 0] <- "STABLE_BASAL"
hybrid_rows$interpretation_category[hybrid_rows$public_subtype == "Classical" & hybrid_rows$assignment_entropy < 0.2 & hybrid_rows$item_consensus >= 0.85 & hybrid_rows$silhouette_width > 0] <- "STABLE_CLASSICAL"
write.table(hybrid_rows, tbl_path("phase4b_hybrid_stability_assessment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sensitivity_pairs <- data.frame(
  comparison_id = c("log2_transformation", "primary_outlier_exclusion", "hvg_outlier_exclusion", "feature_resampling", "imputed_missingness", "hvg_var_filter"),
  reference_analysis = c("STAB_CSY_PRIMARY", "STAB_CSY_PRIMARY", "STAB_UNSUP_HVG", "STAB_CSY_PRIMARY", "STAB_CSY_PRIMARY", "STAB_UNSUP_HVG"),
  sensitivity_analysis = c("STAB_CSY_LOG2", "STAB_CSY_OUTLIER_EXCL", "STAB_HVG_OUTLIER_EXCL", "STAB_CSY_FEAT_RESAMP", "STAB_CSY_IMPUTED", "STAB_HVG_VAR_FILTER"),
  stringsAsFactors = FALSE
)
sens_rows <- list()
for (i in seq_len(nrow(sensitivity_pairs))) {
  p <- sensitivity_pairs[i, ]
  for (k in K_VALUES) {
    rlab <- all_final[[paste(p$reference_analysis, k, sep = "__")]]
    slab <- all_final[[paste(p$sensitivity_analysis, k, sep = "__")]]
    common <- intersect(names(rlab), names(slab))
    rmet <- cluster_metrics[cluster_metrics$analysis_id == p$reference_analysis & cluster_metrics$candidate_K == k, ]
    smet <- cluster_metrics[cluster_metrics$analysis_id == p$sensitivity_analysis & cluster_metrics$candidate_K == k, ]
    sens_rows[[length(sens_rows) + 1L]] <- data.frame(
      comparison_id = p$comparison_id,
      reference_analysis = p$reference_analysis,
      sensitivity_analysis = p$sensitivity_analysis,
      candidate_K = k,
      common_sample_n = length(common),
      assignment_changes = sum(rlab[common] != slab[common]),
      ARI_between_analyses = adjustedRandIndex(rlab[common], slab[common]),
      reference_preferred_K = preferred$preferred_K[preferred$analysis_id == p$reference_analysis],
      sensitivity_preferred_K = preferred$preferred_K[preferred$analysis_id == p$sensitivity_analysis],
      PAC_change = smet$PAC - rmet$PAC,
      silhouette_change = smet$overall_mean_silhouette - rmet$overall_mean_silhouette,
      stringsAsFactors = FALSE
    )
  }
}
sensitivity <- do.call(rbind, sens_rows)
write.table(sensitivity, tbl_path("phase4b_sensitivity_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

unstable <- sample_stability[sample_stability$candidate_K == 3, ]
unstable$unstable_flag <- unstable$assignment_entropy >= 0.5 | unstable$item_consensus < 0.7 | unstable$silhouette_width <= 0
recur <- aggregate(unstable_flag ~ sample_id + patient_id + public_subtype, unstable, sum)
names(recur)[4] <- "unstable_analysis_count"
recur$total_analysis_count <- length(unique(unstable$analysis_id))
recur$recurrently_unstable <- recur$unstable_analysis_count >= 3
write.table(recur[order(-recur$unstable_analysis_count, recur$sample_id), ], tbl_path("phase4b_recurrently_unstable_samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

version_info <- data.frame(
  component = c("R", "cluster", "clue", "mclust", "ggplot2"),
  version = c(R.version.string, as.character(packageVersion("cluster")), as.character(packageVersion("clue")), as.character(packageVersion("mclust")), as.character(packageVersion("ggplot2"))),
  stringsAsFactors = FALSE
)
write.table(version_info, tbl_path("phase4b_package_versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

make_heatmap_pdf(consensus_store[["STAB_CSY_PRIMARY__2"]], all_final[["STAB_CSY_PRIMARY__2"]], fig_path("phase4b_consensus_matrix_primary_K2.pdf"), "Primary CSY consensus matrix, K=2")
make_heatmap_pdf(consensus_store[["STAB_CSY_PRIMARY__3"]], all_final[["STAB_CSY_PRIMARY__3"]], fig_path("phase4b_consensus_matrix_primary_K3.pdf"), "Primary CSY consensus matrix, K=3")
make_heatmap_pdf(consensus_store[["STAB_CSY_PRIMARY__4"]], all_final[["STAB_CSY_PRIMARY__4"]], fig_path("phase4b_consensus_matrix_primary_K4.pdf"), "Primary CSY consensus matrix, K=4")

cdf_df <- do.call(rbind, lapply(names(consensus_store), function(key) {
  parts <- strsplit(key, "__")[[1]]
  vals <- sort(consensus_store[[key]][upper.tri(consensus_store[[key]])])
  data.frame(analysis_id = parts[1], candidate_K = as.integer(parts[2]), consensus = vals, cdf = seq_along(vals) / length(vals))
}))
pdf(fig_path("phase4b_consensus_cdf.pdf"), width = 8, height = 5.5)
print(ggplot(cdf_df, aes(consensus, cdf, color = analysis_id, linetype = factor(candidate_K))) + geom_line(linewidth = 0.45) + theme_minimal(base_size = 9) + labs(x = "Consensus index", y = "CDF", color = "Analysis", linetype = "K"))
dev.off()

plot_metric <- function(y, file, ylab) {
  pdf(fig_path(file), width = 8, height = 5.5)
  print(ggplot(cluster_metrics, aes(candidate_K, .data[[y]], color = analysis_id, group = analysis_id)) + geom_line() + geom_point() + theme_minimal(base_size = 9) + labs(x = "Candidate K", y = ylab, color = "Analysis"))
  dev.off()
}
plot_metric("PAC", "phase4b_pac_by_K.pdf", "PAC")
plot_metric("overall_mean_silhouette", "phase4b_silhouette_by_K.pdf", "Mean silhouette")
plot_metric("mean_Jaccard_stability", "phase4b_jaccard_stability_by_K.pdf", "Mean bootstrap Jaccard")

pdf(fig_path("phase4b_sample_item_consensus.pdf"), width = 8, height = 5.5)
print(ggplot(sample_stability[sample_stability$candidate_K == 3, ], aes(public_subtype, item_consensus, color = analysis_id)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) + theme_minimal(base_size = 9) + labs(x = "Public subtype", y = "Item consensus", color = "Analysis"))
dev.off()

pdf(fig_path("phase4b_sample_assignment_entropy.pdf"), width = 8, height = 5.5)
print(ggplot(sample_stability[sample_stability$candidate_K == 3, ], aes(public_subtype, assignment_entropy, color = analysis_id)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) + theme_minimal(base_size = 9) + labs(x = "Public subtype", y = "Assignment entropy", color = "Analysis"))
dev.off()

primary_labels <- do.call(cbind, lapply(all_final[grepl("__3$", names(all_final))], function(x) x[names(public_labels)]))
colnames(primary_labels) <- sub("__3$", "", names(all_final)[grepl("__3$", names(all_final))])
ari_mat <- matrix(NA_real_, ncol(primary_labels), ncol(primary_labels), dimnames = list(colnames(primary_labels), colnames(primary_labels)))
for (i in seq_len(ncol(primary_labels))) for (j in seq_len(ncol(primary_labels))) {
  ok <- !is.na(primary_labels[, i]) & !is.na(primary_labels[, j])
  ari_mat[i, j] <- adjustedRandIndex(primary_labels[ok, i], primary_labels[ok, j])
}
ari_df <- as.data.frame(as.table(ari_mat))
names(ari_df) <- c("analysis_i", "analysis_j", "ARI")
pdf(fig_path("phase4b_analysis_concordance_heatmap.pdf"), width = 7.2, height = 6)
print(ggplot(ari_df, aes(analysis_i, analysis_j, fill = ARI)) + geom_tile() + geom_text(aes(label = sprintf("%.2f", ARI)), size = 2.5) + scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0.5, limits = c(0, 1)) + theme_minimal(base_size = 8) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + labs(x = NULL, y = NULL, fill = "ARI"))
dev.off()

pdf(fig_path("phase4b_hybrid_stability_summary.pdf"), width = 8, height = 5.5)
print(ggplot(hybrid_rows[hybrid_rows$public_subtype %in% c("Basal", "Classical", "Hybrid"), ], aes(public_subtype, assignment_entropy, fill = interpretation_category)) + geom_boxplot() + facet_wrap(~analysis_id) + theme_minimal(base_size = 8) + labs(x = "Public subtype", y = "Entropy", fill = "Category"))
dev.off()

axis_plot <- merge(axis_df, primary_k3[, c("sample_id", "final_cluster_assignment", "assignment_entropy")], by = "sample_id", all.x = TRUE)
pdf(fig_path("phase4b_basal_classical_axis_with_clusters.pdf"), width = 7.5, height = 5.8)
print(ggplot(axis_plot, aes(moffitt_score_difference, purist_basal_probability, color = factor(final_cluster_assignment), shape = public_subtype, size = assignment_entropy)) + geom_point(alpha = 0.85) + theme_minimal(base_size = 9) + labs(x = "Moffitt basal - classical score", y = "PurIST basal probability", color = "Primary K=3 cluster", shape = "Public subtype", size = "Entropy"))
dev.off()

message("Phase 4B R execution complete in ", round(as.numeric(difftime(Sys.time(), start_all, units = "mins")), 2), " minutes.")
