#!/usr/bin/env Rscript
# 08_build_report.R  (05_models_other — XGBoost + Random Forest soft-voting ENSEMBLE)
#
# Tests whether blending XGBoost with a Random Forest (ranger) improves cell-type classification,
# especially the rare/critical classes that matter under Tumor imbalance. Combination = weighted
# SOFT-VOTING on class probabilities:  P = alpha * P_xgb + (1 - alpha) * P_rf  (argmax to decode).
# alpha = 1 -> XGB-only, alpha = 0 -> RF-only fall out for free, so the same run benchmarks both
# single models and the blend.
#
# STAGED grid (mirrors how 01/02/04 are staged):
#   Stage A — model permutations at the fixed guardrail-best tier split: ~3 XGB x ~3 RF configs,
#             swept over alpha. Each XGB/RF config is trained ONCE (on the fixed split) and the blends
#             are evaluated cheaply, so the alpha sweep is almost free.
#   Stage B — split/knob permutation on the Stage-A winner: a 3x3 tier-split factorial
#             (critical frac {0.20,0.30,0.40} x critical weight {1,40,100}) re-evaluating the chosen
#             ensemble. (Caveat, as in 04: each split yields a DIFFERENT held-out test set, so
#             cross-split metric diffs conflate model quality with which images land in test.)
#
# Selection: guardrail Accuracy_test >= 0.93 & Kappa_test >= 0.85, then rank by critical-class recall
# (same rule as the rest of sys_results). The report answers: does the ensemble beat the best single
# model on critical recall?
#
# Writes (suffixed; leaves 04's tuning_results_split_v2.csv untouched):
#   tuning_results_ensemble.csv  + report_ensemble.rds   (in Panel_2_10/ROut)
#
# RUN:  Rscript 08_build_report.R     (offline; tens of minutes)
# NEW DEPENDENCY: ranger (not used elsewhere in the repo).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(caret)
  library(xgboost)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

if (!requireNamespace("ranger", quietly = TRUE)) {
  stop("This experiment needs the 'ranger' package. Install it with: install.packages('ranger')")
}
message("Using ranger ", as.character(utils::packageVersion("ranger")))

base_path <- "/Users/svitlanamoiseyenko/Repos/genomic/xgboost-imc-panel2-cell-classifier/XGBoost_classifier/Panel_2_10"
out_dir   <- file.path(base_path, "ROut")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
NTHREAD <- max(1, parallel::detectCores() - 1)

# ---------------------------------------------------------------------------
# 1. Rebuild lab_sce from spe.rds (verbatim from the sys_results label pipeline).
# ---------------------------------------------------------------------------
lab_rds <- file.path(out_dir, "lab_sce.rds")
if (file.exists(lab_rds)) {
  message("Loading cached lab_sce ...")
  lab_sce <- readRDS(lab_rds)
} else {
  message("Rebuilding lab_sce from spe.rds ...")
  spe <- readRDS(file.path(out_dir, "spe.rds"))

  gates <- list.files(file.path(base_path, "gates"), pattern = "\\.rds$",
                      full.names = TRUE, recursive = TRUE)
  gate_labels <- purrr::map_dfr(gates, function(path) {
    x <- readRDS(path)
    data.frame(sample_id            = as.character(x$sample_id),
               ObjectNumber         = as.character(x$ObjectNumber),
               cytomapper_CellLabel = as.character(x$cytomapper_CellLabel),
               gate_file            = basename(path), row.names = NULL)
  })
  gate_labels$cytomapper_CellLabel <- sub("_[0-9]+[a-zA-Z]?$", "", gate_labels$cytomapper_CellLabel)
  gate_labels$cytomapper_CellLabel <- sub(" \\([0-9]+\\)$",   "", gate_labels$cytomapper_CellLabel)
  gate_labels$cytomapper_CellLabel[gate_labels$cytomapper_CellLabel == "Cd4"]   <- "CD4"
  gate_labels$cytomapper_CellLabel[gate_labels$cytomapper_CellLabel == "MDSC_"] <- "MDSC"
  gate_labels$cytomapper_CellLabel[gate_labels$cytomapper_CellLabel == "Cd8"]   <- "CD8"
  gate_labels$cell_id <- paste0(gate_labels$sample_id, "_", gate_labels$ObjectNumber)

  matched_gate_labels <- gate_labels[gate_labels$cell_id %in% colnames(spe), ]
  images_28 <- sort(unique(matched_gate_labels$sample_id))

  blank <- "10123883-USL-VAR-TIS-01-IMC-02_005"
  spe_27 <- spe[, spe$sample_id %in% images_28 & spe$sample_id != blank]
  gate_labels <- gate_labels[gate_labels$sample_id != blank, ]

  cur_tab <- unclass(table(gate_labels$cell_id, gate_labels$cytomapper_CellLabel))
  cur_tab_raw <- cur_tab
  cur_tab <- cur_tab_raw; cur_tab[cur_tab > 1] <- 1
  cur_labels <- rep("doublets", nrow(cur_tab)); names(cur_labels) <- rownames(cur_tab)
  single_index <- rowSums(cur_tab) == 1
  cur_labels[single_index] <- colnames(cur_tab)[apply(cur_tab[single_index, , drop = FALSE], 1, which.max)]
  double_tumor_index <- rowSums(cur_tab) == 2 & cur_tab[, "Tumor"] == 1
  if (sum(double_tumor_index) > 0) {
    no_tumor <- cur_tab[, colnames(cur_tab) != "Tumor", drop = FALSE]
    cur_labels[double_tumor_index] <- colnames(no_tumor)[apply(no_tumor[double_tumor_index, , drop = FALSE], 1, which.max)]
  }
  cur_labels <- cur_labels[cur_labels != "doublets"]

  sce_labels <- rep("unlabelled", ncol(spe_27)); names(sce_labels) <- colnames(spe_27)
  matched_ids <- intersect(names(cur_labels), colnames(spe_27))
  sce_labels[matched_ids] <- cur_labels[matched_ids]
  spe_27$cell_labels <- sce_labels

  lab_sce <- spe_27[, !is.na(spe_27$cell_labels) & spe_27$cell_labels != "unlabelled"]
  em <- assay(lab_sce, "exprs"); lab_sce <- lab_sce[, !(lab_sce$cell_labels == "Igg" & em["IgG", ] < 1)]
  em <- assay(lab_sce, "exprs"); lab_sce <- lab_sce[, !(lab_sce$cell_labels == "NK"  & em["CD7", ] < 0.5)]
  em <- assay(lab_sce, "exprs")
  mreg2tum <- lab_sce$cell_labels == "mregDC" & (em["panCK", ] > 0.5 | em["Ecad", ] > 0.8)
  lab_sce$cell_labels[mreg2tum] <- "Tumor"

  saveRDS(lab_sce, lab_rds)
}

marker_keep_global <- !grepl("DNA|Histone", rownames(lab_sce), ignore.case = TRUE)
lab_sce_base <- lab_sce[, !is.na(lab_sce$cell_labels) & lab_sce$cell_labels != "unlabelled"]
n_features <- sum(marker_keep_global)
message("lab_sce_base: ", ncol(lab_sce_base), " cells, ", n_features, " markers")
print(table(lab_sce_base$cell_labels))

# ---------------------------------------------------------------------------
# 2. Harness helpers (verbatim from 08_systematic_tuning_01_svitlana.Rmd).
# ---------------------------------------------------------------------------
cell_priority <- list(
  critical     = c("MDSC", "BnT", "Igg", "mregDC", "PMN_MDSC", "Neutrophil"),
  medium       = c("B", "NK", "DC", "MacCD209", "Treg", "vCAF"),
  non_critical = c("Tumor", "CD4", "CD8", "Mac", "MacCD204", "Fibro",
                   "SMA", "Vasculature", "undefined")
)
priority_of <- function(lab) dplyr::case_when(
  lab %in% cell_priority$critical ~ "critical",
  lab %in% cell_priority$medium   ~ "medium",
  TRUE                            ~ "non_critical")
critical_classes <- cell_priority$critical

split_test_rows_tiered <- function(counts, cell_priority,
                                    tier_test_frac = c(critical = 0.30, medium = 0.25, non_critical = 0.20),
                                    tier_weight    = c(critical = 40,   medium = 10,   non_critical = 1),
                                    min_test_count = c(critical = 8,    medium = 5,    non_critical = 0),
                                    penalty_scale  = 1, n_restarts = 1000, seed = 42) {
  set.seed(seed)
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  labels <- colnames(counts); totals <- colSums(counts); n <- nrow(counts)
  tier <- priority_of(labels)
  target_frac <- as.numeric(tier_test_frac[tier]); names(target_frac) <- labels
  min_tc      <- as.numeric(min_test_count[tier]); names(min_tc)      <- labels
  active <- totals > 0
  w <- as.numeric(tier_weight[tier]) * active
  inv_totals <- ifelse(active, 1 / pmax(totals, 1), 0)
  min_tc <- min_tc * active
  score_of <- function(test_sum) {
    err <- test_sum * inv_totals - target_frac; pen <- pmax(0, min_tc - test_sum)
    sum(w * err * err) + penalty_scale * sum(w * pen * pen)
  }
  row_prob <- apply(counts, 1, function(r) { pres <- r > 0
    if (any(pres)) max(target_frac[pres]) else min(tier_test_frac) })
  best_score <- Inf; best_in_test <- logical(n)
  for (restart in seq_len(n_restarts)) {
    in_test <- runif(n) < row_prob
    test_sum <- colSums(counts[in_test, , drop = FALSE]); cur_score <- score_of(test_sum)
    repeat {
      sgn <- 1 - 2 * in_test
      cand <- sweep(sgn * counts, 2, test_sum, FUN = "+")
      err <- sweep(sweep(cand, 2, inv_totals, "*"), 2, target_frac, "-")
      pen <- pmax(sweep(-cand, 2, min_tc, "+"), 0)
      cand_score <- as.numeric((err * err) %*% w) + penalty_scale * as.numeric((pen * pen) %*% w)
      best_i <- which.min(cand_score)
      if (cand_score[best_i] >= cur_score - 1e-12) break
      test_sum <- test_sum + sgn[best_i] * counts[best_i, ]
      in_test[best_i] <- !in_test[best_i]; cur_score <- cand_score[best_i]
    }
    if (cur_score < best_score) { best_score <- cur_score; best_in_test <- in_test }
  }
  test_sum <- colSums(counts[best_in_test, , drop = FALSE]); train_sum <- totals - test_sum
  achieved <- ifelse(active, test_sum * inv_totals, NA_real_)
  summary_df <- data.frame(label = labels, tier = tier, total = totals,
    test_count = test_sum, train_count = train_sum, target_test_frac = target_frac,
    achieved_test_frac = round(achieved, 3), min_test_count = min_tc,
    meets_min = test_sum >= min_tc, row.names = NULL)
  list(test_idx = which(best_in_test), train_idx = which(!best_in_test),
       test_rownames = rownames(counts)[best_in_test],
       train_rownames = rownames(counts)[!best_in_test],
       summary = summary_df, score = best_score)
}

area_vector <- function(sce) {
  cd <- colData(sce); nm <- grep("^area$", colnames(cd), ignore.case = TRUE, value = TRUE)
  if (length(nm) == 0) return(NULL); as.numeric(cd[[nm[1]]])
}
profile_outlier_cor <- function(sce, marker_keep, min_class_n = 30) {
  mat <- as.matrix(assay(sce, "exprs")[marker_keep, , drop = FALSE])
  labs <- as.character(sce$cell_labels); cors <- rep(NA_real_, ncol(sce))
  for (cl in unique(labs)) {
    idx <- which(labs == cl); if (length(idx) < min_class_n) next
    centroid <- apply(mat[, idx, drop = FALSE], 1, median)
    cors[idx] <- apply(mat[, idx, drop = FALSE], 2, function(v)
      suppressWarnings(stats::cor(v, centroid, method = "spearman")))
  }
  cors
}
apply_cleaning <- function(sce, cfg, marker_keep) {
  if (isTRUE(cfg$drop_undefined)) sce <- sce[, sce$cell_labels != "undefined"]
  if (isTRUE(cfg$area_qc)) {
    area <- area_vector(sce)
    if (!is.null(area)) {
      la <- log(pmax(area, 1))
      keep <- la >= (median(la) - 3 * mad(la)) & la <= (median(la) + 3 * mad(la))
      sce <- sce[, keep]
    }
  }
  if (isTRUE(cfg$outlier_remove)) {
    cors <- profile_outlier_cor(sce, marker_keep)
    keep <- is.na(cors) | cors >= (cfg$outlier_min_cor %||% 0); sce <- sce[, keep]
  }
  sce
}

make_xy <- function(train_sce, test_sce, marker_keep, weight_cap = c(0.2, 5), weight_overrides = NULL) {
  train_sce$cell_id <- paste0(train_sce$sample_id, "_", train_sce$ObjectNumber)
  test_sce$cell_id  <- paste0(test_sce$sample_id,  "_", test_sce$ObjectNumber)
  colnames(train_sce) <- train_sce$cell_id; colnames(test_sce) <- test_sce$cell_id
  X <- t(assay(train_sce, "exprs")[marker_keep, ])
  y_true <- factor(train_sce$cell_labels); label_mapping <- levels(y_true)
  y <- as.integer(y_true) - 1; num_classes <- length(label_mapping)
  test_keep <- !is.na(test_sce$cell_labels) & test_sce$cell_labels != "unlabelled" &
    test_sce$cell_labels %in% label_mapping
  X_test <- t(assay(test_sce, "exprs")[marker_keep, test_keep])
  y_test_true <- factor(test_sce$cell_labels[test_keep], levels = label_mapping)
  y_test <- as.integer(y_test_true) - 1
  stopifnot(identical(colnames(X), colnames(X_test)))
  freq <- as.numeric(table(y_true) / length(y_true)); w_class <- 1 / freq
  names(w_class) <- levels(y_true)
  if (!is.null(weight_overrides)) for (nm in names(weight_overrides))
    if (nm %in% names(w_class)) w_class[nm] <- w_class[nm] * weight_overrides[[nm]]
  w_class <- w_class / mean(w_class)
  w_class <- pmin(pmax(w_class, weight_cap[1]), weight_cap[2])
  w <- as.numeric(w_class[as.character(y_true)]); w <- w / mean(w)
  list(Xy = xgboost::xgb.DMatrix(data = X, label = y, weight = w),
       Xy_test = xgboost::xgb.DMatrix(data = X_test, label = y_test),
       X = X, X_test = X_test, y_true = y_true, y_test_true = y_test_true,
       label_mapping = label_mapping, num_classes = num_classes)
}

make_folds <- function(train_sce, k = 5) {
  tr_idx  <- caret::groupKFold(factor(train_sce$sample_id), k = k)
  all_idx <- seq_len(ncol(train_sce))
  lapply(tr_idx, function(idx) setdiff(all_idx, idx))
}

to_prob_matrix <- function(p, n, num_classes) {
  if (is.list(p)) p <- do.call(rbind, p)
  if (is.null(dim(p))) {
    byrow_m <- matrix(p, ncol = num_classes, byrow = TRUE)
    bycol_m <- matrix(p, ncol = num_classes, byrow = FALSE)
    p <- if (mean(abs(rowSums(byrow_m) - 1)) <= mean(abs(rowSums(bycol_m) - 1))) byrow_m else bycol_m
  }
  as.matrix(p)
}

extract_metrics <- function(pred, truth, label_mapping, critical, suffix) {
  cm <- caret::confusionMatrix(
    data      = factor(as.character(pred),  levels = label_mapping),
    reference = factor(as.character(truth), levels = label_mapping))
  ov <- cm$overall
  bc <- as.data.frame(cm$byClass)
  rownames(bc) <- sub("Class: ", "", rownames(bc))
  crit <- intersect(critical, rownames(bc))
  out <- list(
    Accuracy = unname(ov["Accuracy"]), Kappa = unname(ov["Kappa"]),
    AccuracyNull = unname(ov["AccuracyNull"]),
    BalancedAccuracy_macro    = mean(bc[["Balanced Accuracy"]], na.rm = TRUE),
    F1_macro                  = mean(bc[["F1"]],               na.rm = TRUE),
    Sensitivity_macro         = mean(bc[["Sensitivity"]],      na.rm = TRUE),
    Sensitivity_critical_mean = mean(bc[crit, "Sensitivity"],  na.rm = TRUE),
    Precision_critical_mean   = mean(bc[crit, "Pos Pred Value"], na.rm = TRUE),
    F1_critical_mean          = mean(bc[crit, "F1"],           na.rm = TRUE))
  for (cl in critical)
    out[[paste0("Sens_", cl)]] <- if (cl %in% rownames(bc)) bc[cl, "Sensitivity"] else NA_real_
  stats::setNames(out, paste0(names(out), "_", suffix))
}

# ---------------------------------------------------------------------------
# 2b. NEW — ensemble building blocks.
# ---------------------------------------------------------------------------

# Same per-cell weight formula make_xy() uses, exposed for ranger case.weights.
cell_weights <- function(y_true, weight_cap = c(0.2, 5), weight_overrides = NULL) {
  y_true <- factor(y_true)
  w_class <- 1 / as.numeric(table(y_true) / length(y_true))
  names(w_class) <- levels(y_true)
  if (!is.null(weight_overrides)) for (nm in names(weight_overrides))
    if (nm %in% names(w_class)) w_class[nm] <- w_class[nm] * weight_overrides[[nm]]
  w_class <- w_class / mean(w_class)
  w_class <- pmin(pmax(w_class, weight_cap[1]), weight_cap[2])
  w <- as.numeric(w_class[as.character(y_true)]); w / mean(w)
}

# Prepare one (split x cleaning) partition: features, folds, weights, sizes.
prep <- function(split_cfg, clean_cfg) {
  lab_clean <- apply_cleaning(lab_sce_base, clean_cfg, marker_keep_global)
  dat <- as.data.frame(unclass(table(lab_clean$sample_id, lab_clean$cell_labels)))
  sp  <- do.call(split_test_rows_tiered,
                 c(list(counts = dat, cell_priority = cell_priority), split_cfg))
  train_sce <- lab_clean[, lab_clean$sample_id %in% sp$train_rownames]
  test_sce  <- lab_clean[, lab_clean$sample_id %in% sp$test_rownames]
  d <- make_xy(train_sce, test_sce, marker_keep_global,
               weight_cap = clean_cfg$weight_cap %||% c(0.2, 5),
               weight_overrides = clean_cfg$weight_overrides)
  list(d = d,
       folds = make_folds(train_sce),
       w = cell_weights(d$y_true, clean_cfg$weight_cap %||% c(0.2, 5), clean_cfg$weight_overrides),
       sp = sp, n_train = ncol(train_sce), n_test = ncol(test_sce),
       n_test_images = length(sp$test_rownames))
}

# Train XGB once on a partition: out-of-fold + test probability matrices (cols = label_mapping).
fit_xgb <- function(P, xgb_cfg, nrounds = 500L) {
  d <- P$d
  xparams <- c(list(booster = "gbtree", objective = "multi:softprob", eval_metric = "mlogloss",
                    num_class = d$num_classes, tree_method = "hist", nthread = NTHREAD), xgb_cfg)
  cv <- xgboost::xgb.cv(params = xparams, data = d$Xy, folds = P$folds, nrounds = nrounds,
                        early_stopping_rounds = 20, verbose = FALSE, prediction = TRUE)
  best_it <- which.min(cv$evaluation_log$test_mlogloss_mean)
  oof <- to_prob_matrix(cv$pred %||% cv$cv_predict, length(d$y_true), d$num_classes)
  colnames(oof) <- d$label_mapping
  fm <- xgboost::xgb.train(params = xparams, data = d$Xy, nrounds = best_it, verbose = FALSE)
  test <- to_prob_matrix(predict(fm, d$Xy_test), length(d$y_test_true), d$num_classes)
  colnames(test) <- d$label_mapping
  list(oof = oof, test = test, best_it = best_it)
}

# Train a ranger probability forest: OOF (same folds as XGB) + test prob matrices, aligned to
# label_mapping (a class absent from a fold's training split stays 0).
rf_oof_and_test <- function(P, rf_cfg) {
  d <- P$d; lm <- d$label_mapping; K <- d$num_classes
  Xtr <- as.matrix(d$X)
  df_tr <- data.frame(Xtr, .y = d$y_true, check.names = FALSE)
  rf_args <- function(idx) list(
    dependent.variable.name = ".y", data = df_tr[idx, , drop = FALSE], probability = TRUE,
    num.trees = rf_cfg$num.trees, mtry = rf_cfg$mtry, min.node.size = rf_cfg$min.node.size,
    max.depth = rf_cfg$max.depth %||% 0, case.weights = P$w[idx], num.threads = NTHREAD)

  oof <- matrix(0, nrow(Xtr), K, dimnames = list(NULL, lm))
  for (val_idx in P$folds) {
    tr_idx <- setdiff(seq_len(nrow(Xtr)), val_idx)
    rf <- do.call(ranger::ranger, rf_args(tr_idx))
    pr <- predict(rf, data = df_tr[val_idx, , drop = FALSE])$predictions
    oof[val_idx, colnames(pr)] <- pr
  }
  rf_full <- do.call(ranger::ranger, rf_args(seq_len(nrow(Xtr))))
  df_te <- data.frame(as.matrix(d$X_test), check.names = FALSE)
  prt <- predict(rf_full, data = df_te)$predictions
  test <- matrix(0, nrow(df_te), K, dimnames = list(NULL, lm))
  test[, colnames(prt)] <- prt
  list(oof = oof, test = test)
}

# Soft-voting blend. alpha>=1 -> XGB only; alpha<=0 -> RF only.
blend_probs <- function(xgb_p, rf_p, alpha) {
  if (alpha >= 1) xgb_p else if (alpha <= 0) rf_p else alpha * xgb_p + (1 - alpha) * rf_p
}

# CV + test metrics for one blend.
eval_blend <- function(P, xgb_fit, rf_fit, alpha, critical) {
  d <- P$d
  oof  <- blend_probs(xgb_fit$oof,  rf_fit$oof,  alpha)
  test <- blend_probs(xgb_fit$test, rf_fit$test, alpha)
  cv_pred   <- d$label_mapping[max.col(oof,  ties.method = "first")]
  test_pred <- d$label_mapping[max.col(test, ties.method = "first")]
  c(extract_metrics(cv_pred,   d$y_true,      d$label_mapping, critical, "cv"),
    extract_metrics(test_pred, d$y_test_true, d$label_mapping, critical, "test"))
}

# Held-out confusion matrix for one blend (for the winners).
cm_of <- function(P, xgb_fit, rf_fit, alpha) {
  test <- blend_probs(xgb_fit$test, rf_fit$test, alpha)
  pred <- factor(P$d$label_mapping[max.col(test, ties.method = "first")],
                 levels = P$d$label_mapping)
  caret::confusionMatrix(pred, P$d$y_test_true)
}

# Flatten metadata + metrics into one result row.
make_row <- function(stage, kind, xid, rid, alpha, split_cfg, P, xgb_cfg, rf_cfg,
                     nrounds_best, metrics) {
  scfg <- split_cfg
  meta <- list(
    stage = stage, model_kind = kind,
    model_id = paste0(stage, "_", xid, "_", rid, "_a", alpha),
    xgb_id = xid, rf_id = rid, alpha = alpha,
    tier_test_frac_critical = unname((scfg$tier_test_frac %||% c(critical = 0.30))["critical"]),
    tier_weight_critical    = unname((scfg$tier_weight %||% c(critical = 40))["critical"]),
    min_test_count_critical = unname((scfg$min_test_count %||% c(critical = 8))["critical"]),
    seed = scfg$seed %||% 42,
    n_train_cells = P$n_train, n_test_cells = P$n_test, n_test_images = P$n_test_images,
    xgb_max_depth = xgb_cfg$max_depth %||% NA, xgb_eta = xgb_cfg$eta %||% NA,
    xgb_min_child_weight = xgb_cfg$min_child_weight %||% NA,
    xgb_gamma = xgb_cfg$gamma %||% NA, xgb_colsample = xgb_cfg$colsample_bytree %||% NA,
    xgb_nrounds_best = nrounds_best %||% NA,
    rf_num_trees = rf_cfg$num.trees %||% NA, rf_mtry = rf_cfg$mtry %||% NA,
    rf_min_node_size = rf_cfg$min.node.size %||% NA, rf_max_depth = rf_cfg$max.depth %||% NA)
  as.data.frame(c(meta, metrics), stringsAsFactors = FALSE, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# 3. Configs.
# ---------------------------------------------------------------------------
clean_cfg_baseline <- list(drop_undefined = FALSE, area_qc = FALSE, outlier_remove = FALSE,
  weight_cap = c(0, Inf), weight_overrides = list(undefined = 1/10, Igg = 1/30),
  weight_scheme = "legacy_invfreq")

split_main <- list(tier_test_frac = c(critical = 0.30, medium = 0.25, non_critical = 0.20),
                   tier_weight    = c(critical = 40,   medium = 10,   non_critical = 1),
                   min_test_count = c(critical = 8,    medium = 5,    non_critical = 0),
                   n_restarts = 1000, seed = 42)

xgb_configs <- list(
  conservative = list(max_depth = 3, eta = 0.05, min_child_weight = 10, gamma = 1, lambda = 1,
                      alpha = 0, max_delta_step = 1, subsample = 0.8, colsample_bytree = 0.6),
  flexible     = list(max_depth = 4, eta = 0.10, min_child_weight = 1,  gamma = 0, lambda = 1,
                      alpha = 0, max_delta_step = 0, subsample = 0.8, colsample_bytree = 0.8),
  baseline     = list(max_depth = 3, eta = 0.10, min_child_weight = 5,  gamma = 0, lambda = 1,
                      alpha = 0, max_delta_step = 0, subsample = 0.8, colsample_bytree = 0.8))

rf_configs <- list(
  rf_shallowmtry = list(num.trees = 500,  mtry = floor(sqrt(n_features)), min.node.size = 1,  max.depth = 0),
  rf_widemtry    = list(num.trees = 1000, mtry = floor(n_features / 3),   min.node.size = 5,  max.depth = 0),
  rf_regularised = list(num.trees = 1000, mtry = floor(sqrt(n_features)), min.node.size = 10, max.depth = 0))

alpha_interior <- c(0.25, 0.50, 0.75)   # ensemble blends; singletons (1, 0) emitted separately
NROUNDS <- 500L

# ---------------------------------------------------------------------------
# 4. Stage A — model permutations at the fixed best split (train each model once).
# ---------------------------------------------------------------------------
message("Stage A: preparing fixed split ...")
PA <- prep(split_main, clean_cfg_baseline)

message("Stage A: training ", length(xgb_configs), " XGB configs ...")
xgb_fits <- lapply(names(xgb_configs), function(xid) {
  message("  xgb: ", xid); fit_xgb(PA, xgb_configs[[xid]], NROUNDS)
}); names(xgb_fits) <- names(xgb_configs)

message("Stage A: training ", length(rf_configs), " RF configs ...")
rf_fits <- lapply(names(rf_configs), function(rid) {
  message("  rf: ", rid); rf_oof_and_test(PA, rf_configs[[rid]])
}); names(rf_fits) <- names(rf_configs)

rowsA <- list()
# ensemble blends (interior alpha) over every XGB x RF combo
for (xid in names(xgb_configs)) for (rid in names(rf_configs)) for (a in alpha_interior) {
  m <- eval_blend(PA, xgb_fits[[xid]], rf_fits[[rid]], a, critical_classes)
  rowsA[[length(rowsA) + 1]] <- make_row("A", "ensemble", xid, rid, a, split_main, PA,
                                         xgb_configs[[xid]], rf_configs[[rid]],
                                         xgb_fits[[xid]]$best_it, m)
}
# XGB-only singletons (alpha = 1; one per XGB config)
for (xid in names(xgb_configs)) {
  m <- eval_blend(PA, xgb_fits[[xid]], rf_fits[[1]], 1, critical_classes)
  rowsA[[length(rowsA) + 1]] <- make_row("A", "xgb_only", xid, NA, 1, split_main, PA,
                                         xgb_configs[[xid]], NULL, xgb_fits[[xid]]$best_it, m)
}
# RF-only singletons (alpha = 0; one per RF config)
for (rid in names(rf_configs)) {
  m <- eval_blend(PA, xgb_fits[[1]], rf_fits[[rid]], 0, critical_classes)
  rowsA[[length(rowsA) + 1]] <- make_row("A", "rf_only", NA, rid, 0, split_main, PA,
                                         NULL, rf_configs[[rid]], NA, m)
}
resultsA <- dplyr::bind_rows(rowsA)

# ---------------------------------------------------------------------------
# 5. Stage-A selection (guardrail + critical recall).
# ---------------------------------------------------------------------------
acc_floor <- 0.93; kappa_floor <- 0.85
pick_best <- function(df, kinds) {
  d <- df[df$model_kind %in% kinds, ]
  e <- d[d$Accuracy_test >= acc_floor & d$Kappa_test >= kappa_floor, ]
  if (!nrow(e)) e <- d
  e[order(-e$Sensitivity_critical_mean_test, -e$F1_critical_mean_test, -e$Kappa_test), ][1, ]
}
best_ens <- pick_best(resultsA, "ensemble")
best_xgb <- pick_best(resultsA, "xgb_only")
best_rf  <- pick_best(resultsA, "rf_only")

message(sprintf("Stage A winners | ensemble: %s (crit recall %.3f) | xgb: %s (%.3f) | rf: %s (%.3f)",
                best_ens$model_id, best_ens$Sensitivity_critical_mean_test,
                best_xgb$model_id, best_xgb$Sensitivity_critical_mean_test,
                best_rf$model_id,  best_rf$Sensitivity_critical_mean_test))

cm_store <- list(
  ensemble = cm_of(PA, xgb_fits[[best_ens$xgb_id]], rf_fits[[best_ens$rf_id]], best_ens$alpha),
  xgb_only = cm_of(PA, xgb_fits[[best_xgb$xgb_id]], rf_fits[[1]], 1),
  rf_only  = cm_of(PA, xgb_fits[[1]], rf_fits[[best_rf$rf_id]], 0))

best_xid <- best_ens$xgb_id; best_rid <- best_ens$rf_id; best_alpha <- best_ens$alpha

# ---------------------------------------------------------------------------
# 6. Stage B — split/knob permutation on the winning ensemble (3x3 tier factorial).
# ---------------------------------------------------------------------------
split_grid <- expand.grid(frac_crit = c(0.20, 0.30, 0.40), w_crit = c(1, 40, 100),
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
message("Stage B: ", nrow(split_grid), " split scenarios on the winning ensemble (",
        best_xid, " + ", best_rid, ", alpha=", best_alpha, ") ...")
rowsB <- list()
for (i in seq_len(nrow(split_grid))) {
  g <- split_grid[i, ]
  scfgB <- list(tier_test_frac = c(critical = g$frac_crit, medium = 0.25, non_critical = 0.20),
                tier_weight    = c(critical = g$w_crit,    medium = 10,   non_critical = 1),
                min_test_count = c(critical = 8,           medium = 5,    non_critical = 0),
                n_restarts = 1000, seed = 42)
  message(sprintf("  split %d/%d (frac=%.2f w=%d)", i, nrow(split_grid), g$frac_crit, g$w_crit))
  PB <- prep(scfgB, clean_cfg_baseline)
  xf <- fit_xgb(PB, xgb_configs[[best_xid]], NROUNDS)
  rf <- rf_oof_and_test(PB, rf_configs[[best_rid]])
  m  <- eval_blend(PB, xf, rf, best_alpha, critical_classes)
  rowsB[[i]] <- make_row("B", "ensemble", best_xid, best_rid, best_alpha, scfgB, PB,
                         xgb_configs[[best_xid]], rf_configs[[best_rid]], xf$best_it, m)
}
resultsB <- dplyr::bind_rows(rowsB)

# ---------------------------------------------------------------------------
# 7. Combine, write outputs.
# ---------------------------------------------------------------------------
results <- dplyr::bind_rows(resultsA, resultsB)

config <- list(split_main = split_main, clean_cfg_baseline = clean_cfg_baseline,
               xgb_configs = xgb_configs, rf_configs = rf_configs,
               alpha_interior = alpha_interior, n_features = n_features,
               acc_floor = acc_floor, kappa_floor = kappa_floor,
               best_ens_id = best_ens$model_id, best_xgb_id = best_xgb$model_id,
               best_rf_id = best_rf$model_id,
               winner = list(xgb_id = best_xid, rf_id = best_rid, alpha = best_alpha))

write.csv(results, file.path(out_dir, "tuning_results_ensemble.csv"), row.names = FALSE)
saveRDS(list(results = results, resultsA = resultsA, resultsB = resultsB,
             cm_store = cm_store, config = config,
             best_ens = best_ens, best_xgb = best_xgb, best_rf = best_rf),
        file.path(out_dir, "report_ensemble.rds"))

message("Wrote tuning_results_ensemble.csv (", nrow(results), " rows) and report_ensemble.rds to ",
        out_dir)
message(sprintf("Headline | ensemble crit recall %.3f vs xgb %.3f vs rf %.3f  (Kappa_test: %.3f / %.3f / %.3f)",
                best_ens$Sensitivity_critical_mean_test, best_xgb$Sensitivity_critical_mean_test,
                best_rf$Sensitivity_critical_mean_test,
                best_ens$Kappa_test, best_xgb$Kappa_test, best_rf$Kappa_test))
