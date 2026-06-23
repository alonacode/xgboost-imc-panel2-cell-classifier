#!/usr/bin/env Rscript
# 08_build_report.R  (split-parameter experiment, v2 — factorial)
# Like 02_split_params, the MODEL is fixed at the guardrail-best params (M9) and the tier-aware
# train/test SPLIT is varied. v2 runs a 3x3x3 FACTORIAL over the three deterministic split knobs:
#   critical test-fraction {0.20, 0.30, 0.40} x critical weight {1, 40, 100} x critical floor {0, 8, 15}
# = 27 scenarios (V1-V27). Seed is fixed at 42 because v1 (02_split_params) established the optimiser
# is seed-stable (S1/S9/S10 converged to the same split; residual recall spread ~0.007 was pure
# model-training stochasticity). Each scenario is evaluated end-to-end (image-grouped xgb.cv to pick
# nrounds_best, final fit, held-out test metrics + confusion matrix).
# Writes tuning_results_split_v2.csv + report_cms_split_v2.rds (suffixed; nothing else clobbered).
#
# NOTE: every split scenario produces a DIFFERENT held-out test set, so cross-scenario metric
# differences conflate model quality with which images land in test. The factorial lets us read each
# knob's MARGINAL effect by averaging over the other two (see the report).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(caret)
  library(xgboost)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

base_path <- "/Users/svitlanamoiseyenko/Repos/genomic/xgboost-imc-panel2-cell-classifier/XGBoost_classifier/Panel_2_10"
out_dir   <- file.path(base_path, "ROut")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---------------------------------------------------------------------------
# 1. Rebuild lab_sce from spe.rds (replays 01_read_data.Rmd label pipeline,
#    skipping all image/mask loading which is only used for Shiny/QC).
#    (verbatim from 01_model_params/08_build_report.R)
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

  # consolidate labels (same as 01)
  cur_tab <- unclass(table(gate_labels$cell_id, gate_labels$cytomapper_CellLabel))
  cur_tab_raw <- cur_tab
  same_label_duplicate_ids <- rownames(cur_tab_raw)[rowSums(cur_tab_raw > 1) > 0]
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

  # biological cleanup rules
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
message("lab_sce_base: ", ncol(lab_sce_base), " cells")
print(table(lab_sce_base$cell_labels))

# ---------------------------------------------------------------------------
# 2. Harness helpers (verbatim from 08_systematic_tuning_01_svitlana.Rmd)
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
    AccuracyLower = unname(ov["AccuracyLower"]), AccuracyUpper = unname(ov["AccuracyUpper"]),
    AccuracyNull = unname(ov["AccuracyNull"]), AccuracyPValue = unname(ov["AccuracyPValue"]),
    McnemarPValue = unname(ov["McnemarPValue"]),
    BalancedAccuracy_macro    = mean(bc[["Balanced Accuracy"]], na.rm = TRUE),
    Sensitivity_critical_mean = mean(bc[crit, "Sensitivity"],     na.rm = TRUE),
    Precision_critical_mean   = mean(bc[crit, "Pos Pred Value"],  na.rm = TRUE),
    F1_critical_mean          = mean(bc[crit, "F1"],              na.rm = TRUE))
  for (cl in critical)
    out[[paste0("Sens_", cl)]] <- if (cl %in% rownames(bc)) bc[cl, "Sensitivity"] else NA_real_
  stats::setNames(out, paste0(names(out), "_", suffix))
}

run_experiment <- function(scenario, lab_sce_base, marker_keep, critical) {
  lab_clean <- apply_cleaning(lab_sce_base, scenario$clean_cfg, marker_keep)
  dat <- as.data.frame(unclass(table(lab_clean$sample_id, lab_clean$cell_labels)))
  sp  <- do.call(split_test_rows_tiered,
                 c(list(counts = dat, cell_priority = cell_priority), scenario$split_cfg))
  train_sce <- lab_clean[, lab_clean$sample_id %in% sp$train_rownames]
  test_sce  <- lab_clean[, lab_clean$sample_id %in% sp$test_rownames]
  d <- make_xy(train_sce, test_sce, marker_keep,
               weight_cap = scenario$clean_cfg$weight_cap %||% c(0.2, 5),
               weight_overrides = scenario$clean_cfg$weight_overrides)
  folds  <- make_folds(train_sce)
  params <- c(list(booster = "gbtree", objective = "multi:softprob", eval_metric = "mlogloss",
                   num_class = d$num_classes, tree_method = "hist",
                   nthread = max(1, parallel::detectCores() - 1)), scenario$xgb)
  cv <- xgboost::xgb.cv(params = params, data = d$Xy, folds = folds,
                        nrounds = scenario$nrounds %||% 500L,
                        early_stopping_rounds = 20, verbose = FALSE, prediction = TRUE)
  best_it <- which.min(cv$evaluation_log$test_mlogloss_mean)
  cv_prob <- to_prob_matrix(cv$pred %||% cv$cv_predict, length(d$y_true), d$num_classes)
  cv_pred <- d$label_mapping[max.col(cv_prob, ties.method = "first")]
  final_model <- xgboost::xgb.train(params = params, data = d$Xy, nrounds = best_it, verbose = FALSE)
  test_prob <- to_prob_matrix(predict(final_model, d$Xy_test), length(d$y_test_true), d$num_classes)
  test_pred <- d$label_mapping[max.col(test_prob, ties.method = "first")]
  m_cv   <- extract_metrics(cv_pred,   d$y_true,      d$label_mapping, critical, "cv")
  m_test <- extract_metrics(test_pred, d$y_test_true, d$label_mapping, critical, "test")
  if (exists("cm_store")) {
    cm_obj <- caret::confusionMatrix(factor(test_pred, levels = d$label_mapping), d$y_test_true)
    cm_store[[scenario$model_id]]    <<- cm_obj$table
    split_store[[scenario$model_id]] <<- sp$summary
  }
  scfg <- scenario$split_cfg
  split_params <- list(
    split_id                = scenario$split_id,
    default_test_frac       = unname((scfg$tier_test_frac %||% c(non_critical = 0.20))["non_critical"]),
    tier_test_frac_critical = unname((scfg$tier_test_frac %||% c(critical = 0.30))["critical"]),
    tier_test_frac_medium   = unname((scfg$tier_test_frac %||% c(medium = 0.25))["medium"]),
    tier_test_frac_noncrit  = unname((scfg$tier_test_frac %||% c(non_critical = 0.20))["non_critical"]),
    tier_weight_critical    = unname((scfg$tier_weight %||% c(critical = 40))["critical"]),
    tier_weight_medium      = unname((scfg$tier_weight %||% c(medium = 10))["medium"]),
    tier_weight_noncrit     = unname((scfg$tier_weight %||% c(non_critical = 1))["non_critical"]),
    min_test_count_critical = unname((scfg$min_test_count %||% c(critical = 8))["critical"]),
    n_restarts              = scfg$n_restarts %||% 1000,
    seed                    = scfg$seed %||% 42,
    n_train_cells           = ncol(train_sce),
    n_test_cells            = ncol(test_sce),
    n_test_images           = length(sp$test_rownames),
    min_test_count_observed = min(sp$summary$test_count))
  model_params <- list(
    model_id = scenario$model_id, max_depth = params$max_depth, eta = params$eta,
    subsample = params$subsample, colsample_bytree = params$colsample_bytree,
    min_child_weight = params$min_child_weight, gamma = params$gamma %||% 0,
    lambda = params$lambda %||% 1, alpha = params$alpha %||% 0,
    max_delta_step = params$max_delta_step %||% 0, nrounds_best = best_it,
    drop_undefined = isTRUE(scenario$clean_cfg$drop_undefined),
    weight_scheme = scenario$clean_cfg$weight_scheme %||% "capped_invfreq")
  as.data.frame(c(scenario_label = scenario$label, split_params, model_params, m_cv, m_test),
                stringsAsFactors = FALSE, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# 3. Fixed model (guardrail-best M9) + baseline cleaning. Only the split varies.
# ---------------------------------------------------------------------------
# M9 "conservative" from 01_model_params (best critical recall behind the acc>=.93/kappa>=.85 guardrail)
best_xgb <- list(max_depth = 3, eta = 0.05, min_child_weight = 10, gamma = 1, lambda = 1,
                 alpha = 0, max_delta_step = 1, subsample = 0.8, colsample_bytree = 0.6)

clean_cfg_baseline <- list(drop_undefined = FALSE, area_qc = FALSE, outlier_remove = FALSE,
  weight_cap = c(0, Inf), weight_overrides = list(undefined = 1/10, Igg = 1/30),
  weight_scheme = "legacy_invfreq")

# ---------------------------------------------------------------------------
# 4. Factorial split grid (V1-V27): critical frac x critical weight x critical floor.
#    medium/non_critical fractions (0.25/0.20), weights (10/1) and floors (5/0) held fixed;
#    seed fixed at 42 (optimiser is seed-stable, established in v1).
# ---------------------------------------------------------------------------
split_grid <- expand.grid(
  frac_crit = c(0.20, 0.30, 0.40),
  w_crit    = c(1, 40, 100),
  min_crit  = c(0, 8, 15),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
split_grid <- split_grid[order(split_grid$frac_crit, split_grid$w_crit, split_grid$min_crit), ]
split_grid$seed <- 42
split_grid$note <- with(split_grid,
  sprintf("frac=%.2f w=%d floor=%d", frac_crit, as.integer(w_crit), as.integer(min_crit)))
rownames(split_grid) <- NULL

scenarios <- lapply(seq_len(nrow(split_grid)), function(i) {
  g <- split_grid[i, ]
  list(label = paste0("split_", i, " (", g$note, ")"), split_id = paste0("V", i),
       model_id = paste0("V", i), nrounds = 500L, clean_cfg = clean_cfg_baseline, xgb = best_xgb,
       split_cfg = list(
         tier_test_frac = c(critical = g$frac_crit, medium = 0.25, non_critical = 0.20),
         tier_weight    = c(critical = g$w_crit,    medium = 10,   non_critical = 1),
         min_test_count = c(critical = g$min_crit,  medium = 5,    non_critical = 0),
         n_restarts = 1000, seed = g$seed))
})

# ---------------------------------------------------------------------------
# 5. Run every split scenario end-to-end (CV -> nrounds_best -> final -> test).
# ---------------------------------------------------------------------------
cm_store <- list(); split_store <- list()
results <- purrr::map_dfr(scenarios, function(sc) {
  message("Running ", sc$label)
  run_experiment(sc, lab_sce_base, marker_keep_global, critical_classes)
})

# Guardrail-best split (same rule as the model report; here it selects a split STRATEGY).
acc_floor <- 0.93; kappa_floor <- 0.85
elig <- results[results$Accuracy_test >= acc_floor & results$Kappa_test >= kappa_floor, ]
if (nrow(elig) == 0) elig <- results
elig <- elig[order(-elig$Sensitivity_critical_mean_test, -elig$F1_critical_mean_test,
                   -elig$Kappa_test), ]
best_s_id <- elig$model_id[1]

write.csv(results, file.path(out_dir, "tuning_results_split_v2.csv"), row.names = FALSE)
saveRDS(list(cm_store = cm_store, split_store = split_store, results = results,
             split_grid = split_grid, best_xgb = best_xgb, best_s_id = best_s_id),
        file.path(out_dir, "report_cms_split_v2.rds"))
message("Wrote tuning_results_split_v2.csv (", nrow(results), " rows) and report_cms_split_v2.rds. ",
        "Best split: ", best_s_id)
