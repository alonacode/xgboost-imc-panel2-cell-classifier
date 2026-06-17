#!/usr/bin/env Rscript
# 08_build_report.R
# Recompute per-scenario confusion matrices + split summaries for the report.
# CMs only need each model retrained at its saved `nrounds_best` (from tuning_results.csv);
# the heavy xgb.cv grid search is NOT re-run. lab_sce is rebuilt from spe.rds (image-free).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(caret)
  library(xgboost)
  library(dplyr)
  library(purrr)
})

base_path <- "/Users/svitlanamoiseyenko/Repos/genomic/xgboost-imc-panel2-cell-classifier/XGBoost_classifier/Panel_2_10"
out_dir   <- file.path(base_path, "ROut")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---------------------------------------------------------------------------
# 1. Rebuild lab_sce from spe.rds (replays 01_read_data.Rmd label pipeline,
#    skipping all image/mask loading which is only used for Shiny/QC).
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
to_prob_matrix <- function(p, n, num_classes) {
  if (is.list(p)) p <- do.call(rbind, p)
  if (is.null(dim(p))) {
    byrow_m <- matrix(p, ncol = num_classes, byrow = TRUE)
    bycol_m <- matrix(p, ncol = num_classes, byrow = FALSE)
    p <- if (mean(abs(rowSums(byrow_m) - 1)) <= mean(abs(rowSums(bycol_m) - 1))) byrow_m else bycol_m
  }
  as.matrix(p)
}

# ---------------------------------------------------------------------------
# 3. Scenario grid (same as Part D) + ablations
# ---------------------------------------------------------------------------
split_cfg_main <- list(
  tier_test_frac = c(critical = 0.30, medium = 0.25, non_critical = 0.20),
  tier_weight    = c(critical = 40,   medium = 10,   non_critical = 1),
  min_test_count = c(critical = 8,    medium = 5,    non_critical = 0),
  n_restarts = 1000, seed = 42)
clean_cfg_baseline <- list(drop_undefined = FALSE, area_qc = FALSE, outlier_remove = FALSE,
  weight_cap = c(0, Inf), weight_overrides = list(undefined = 1/10, Igg = 1/30),
  weight_scheme = "legacy_invfreq")

xgb_grid <- tibble::tribble(
  ~max_depth, ~eta,  ~min_child_weight, ~gamma, ~lambda, ~max_delta_step, ~subsample, ~colsample_bytree, ~note,
  3, 0.10, 5,  0, 1, 0, 0.8, 0.8, "baseline (current best)",
  3, 0.10, 5,  0, 1, 1, 0.8, 0.8, "+max_delta_step",
  3, 0.05, 5,  1, 1, 1, 0.8, 0.8, "slower + gamma",
  4, 0.05, 5,  1, 5, 1, 0.8, 0.8, "+L2",
  4, 0.05, 1,  0, 1, 1, 0.7, 0.6, "more sampling",
  3, 0.03, 5,  1, 5, 2, 0.8, 0.8, "strong reg",
  4, 0.10, 1,  0, 1, 0, 0.8, 0.8, "flexible",
  5, 0.05, 5,  1, 5, 1, 0.8, 0.8, "depth probe",
  3, 0.05, 10, 1, 1, 1, 0.8, 0.6, "conservative",
  4, 0.05, 5,  1, 5, 1, 0.7, 0.8, "balanced reg")

scenarios <- lapply(seq_len(nrow(xgb_grid)), function(i) {
  g <- xgb_grid[i, ]
  list(label = paste0("model_", i, " (", g$note, ")"), split_id = "split_main",
       model_id = paste0("M", i), split_cfg = split_cfg_main, clean_cfg = clean_cfg_baseline,
       xgb = list(max_depth = g$max_depth, eta = g$eta, min_child_weight = g$min_child_weight,
                  gamma = g$gamma, lambda = g$lambda, alpha = 0, max_delta_step = g$max_delta_step,
                  subsample = g$subsample, colsample_bytree = g$colsample_bytree))
})

results_csv <- read.csv(file.path(out_dir, "tuning_results.csv"), stringsAsFactors = FALSE,
                        check.names = FALSE)
results_csv <- results_csv[!duplicated(results_csv$model_id), ]
nrounds_of  <- setNames(results_csv$nrounds_best, results_csv$model_id)

best_m_id   <- results_csv$model_id[which.max(results_csv$Sensitivity_critical_mean_test[
                 grepl("^M", results_csv$model_id)])]
best_xgb    <- scenarios[[as.integer(sub("M", "", best_m_id))]]$xgb

cap_clean <- function(extra) modifyList(list(drop_undefined = FALSE, area_qc = FALSE,
  outlier_remove = FALSE, weight_cap = c(0.2, 5), weight_overrides = NULL,
  weight_scheme = "capped_invfreq"), extra)
ablations <- list(
  list(label = "ablation_capped_weights",  model_id = "A1", split_cfg = split_cfg_main,
       xgb = best_xgb, clean_cfg = cap_clean(list())),
  list(label = "ablation_drop_undefined",  model_id = "A2", split_cfg = split_cfg_main,
       xgb = best_xgb, clean_cfg = cap_clean(list(drop_undefined = TRUE))),
  list(label = "ablation_area_qc",         model_id = "A3", split_cfg = split_cfg_main,
       xgb = best_xgb, clean_cfg = cap_clean(list(area_qc = TRUE))),
  list(label = "ablation_outlier_remove",  model_id = "A4", split_cfg = split_cfg_main,
       xgb = best_xgb, clean_cfg = cap_clean(list(outlier_remove = TRUE, outlier_min_cor = 0.1))))

all_scen <- c(scenarios, ablations)

# ---------------------------------------------------------------------------
# 4. Recompute CM + split summary per scenario (train at saved nrounds_best).
#    Cache the main split (shared by M1-M10 + A1).
# ---------------------------------------------------------------------------
cm_store <- list(); split_store <- list(); main_sp <- NULL
nthread <- max(1, parallel::detectCores() - 1)

for (sc in all_scen) {
  message("Recomputing ", sc$model_id, " - ", sc$label)
  lab_clean <- apply_cleaning(lab_sce_base, sc$clean_cfg, marker_keep_global)
  is_main <- identical(sc$clean_cfg$drop_undefined, FALSE) &&
             identical(sc$clean_cfg$area_qc, FALSE) &&
             identical(sc$clean_cfg$outlier_remove, FALSE)
  if (is_main && !is.null(main_sp)) {
    sp <- main_sp
  } else {
    dat <- as.data.frame(unclass(table(lab_clean$sample_id, lab_clean$cell_labels)))
    sp  <- do.call(split_test_rows_tiered,
                   c(list(counts = dat, cell_priority = cell_priority), sc$split_cfg))
    if (is_main) main_sp <- sp
  }
  train_sce <- lab_clean[, lab_clean$sample_id %in% sp$train_rownames]
  test_sce  <- lab_clean[, lab_clean$sample_id %in% sp$test_rownames]
  d <- make_xy(train_sce, test_sce, marker_keep_global,
               weight_cap = sc$clean_cfg$weight_cap %||% c(0.2, 5),
               weight_overrides = sc$clean_cfg$weight_overrides)
  params <- c(list(booster = "gbtree", objective = "multi:softprob", eval_metric = "mlogloss",
                   num_class = d$num_classes, tree_method = "hist", nthread = nthread), sc$xgb)
  nb <- nrounds_of[[sc$model_id]] %||% 100L
  set.seed(42)
  fm <- xgboost::xgb.train(params = params, data = d$Xy, nrounds = nb, verbose = FALSE)
  test_prob <- to_prob_matrix(predict(fm, d$Xy_test), length(d$y_test_true), d$num_classes)
  test_pred <- factor(d$label_mapping[max.col(test_prob, ties.method = "first")],
                      levels = d$label_mapping)
  cm <- caret::confusionMatrix(test_pred, d$y_test_true)
  cm_store[[sc$model_id]]    <- cm$table
  split_store[[sc$model_id]] <- sp$summary
}

saveRDS(list(cm_store = cm_store, split_store = split_store, results = results_csv,
             xgb_grid = xgb_grid, best_m_id = best_m_id),
        file.path(out_dir, "report_cms.rds"))
message("Wrote report_cms.rds with ", length(cm_store), " confusion matrices.")
