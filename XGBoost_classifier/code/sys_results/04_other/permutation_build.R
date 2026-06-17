#!/usr/bin/env Rscript
# permutation_build.R  — Label-permutation significance test for the 0.994 snapshot classifier.
#
# WHY: the held-out test accuracy (~0.994) is inflated by Tumor dominance (NIR ~0.87: a
# "predict Tumor for everything" model already scores ~0.87). To prove the classifier learns real
# marker->cell-type structure (and not just the class imbalance), we shuffle the cell_labels, refit
# the SAME estimator, repeat N times, and build a NULL distribution for each metric. Observed metrics
# are then compared against the null with an empirical p-value.
#
# WHAT THIS REPLICATES: the LEGACY 0.994 estimator from
#   code/results/01_read_data_svitlana_01_0.994.Rmd
#   - split: split_test_rows() (legacy binary rare/common), seed 42;
#   - weights: inverse-frequency with undefined/10, Igg/30;
#   - model: 12-combo grid (max_depth{3,4,5} x eta{0.05,0.1} x min_child_weight{1,5}),
#            multi:softprob / mlogloss, image-grouped 5-fold CV, early stopping, min-mlogloss select.
#
# KEY DESIGN: permutation happens AFTER the split is fixed and does NOT re-derive the split. The
# train/test image assignment, the CV folds, and the tuned params/nrounds are all held fixed; only
# the labels are shuffled. That isolates the marker->label signal as the only thing destroyed.
#
# OUTPUT (written to Panel_2_10/ROut, no-clobber names):
#   permutation_results.csv  — one row per permutation (CV + test metrics) + an "observed" row.
#   permutation_null.rds     — list(null, observed, pvalues, config).
#
# RUN:  Rscript permutation_build.R      (offline; tens of minutes at N_PERM=100)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(caret)
  library(xgboost)
  library(dplyr)
  library(purrr)
})

# ===========================================================================
# 0. Config
# ===========================================================================
N_PERM    <- 100        # number of label permutations (100 -> min p ~ 0.0099; 1000 ~ 10x runtime)
base_seed <- 42         # per-permutation seed = base_seed + i
DO_CV     <- TRUE       # also build the CV (out-of-fold) null; FALSE = test-only, ~2x faster
null_type <- "global"   # "global" (shuffle across all labelled cells) or "within_image" (stricter)

base_path <- "/Users/svitlanamoiseyenko/Repos/genomic/xgboost-imc-panel2-cell-classifier/XGBoost_classifier/Panel_2_10"
out_dir   <- file.path(base_path, "ROut")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
nthread <- max(1, parallel::detectCores() - 1)

message(sprintf("Config: N_PERM=%d, base_seed=%d, DO_CV=%s, null_type=%s",
                N_PERM, base_seed, DO_CV, null_type))

# ===========================================================================
# 1. Rebuild lab_sce from spe.rds (verbatim from 04_other/08_build_report.R label pipeline).
#    Same 3 biological cleanup rules as the 0.994 snapshot (Igg<1, NK CD7<0.5, mregDC->Tumor).
# ===========================================================================
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

lab_sce <- lab_sce[, !is.na(lab_sce$cell_labels) & lab_sce$cell_labels != "unlabelled"]
marker_keep <- !grepl("DNA|Histone", rownames(lab_sce), ignore.case = TRUE)
message("lab_sce: ", ncol(lab_sce), " labelled cells, ", sum(marker_keep), " markers")
print(table(lab_sce$cell_labels))

critical_classes <- c("MDSC", "BnT", "Igg", "mregDC", "PMN_MDSC", "Neutrophil")

# ===========================================================================
# 2. Legacy estimator helpers
# ===========================================================================

# Legacy image-level split optimiser (verbatim from the 0.994 snapshot / original 01).
split_test_rows <- function(counts,
                            default_test_frac = 0.20,
                            rare_test_frac    = 0.40,
                            rare_labels       = c("MDSC"),
                            rare_threshold    = NULL,
                            rare_weight       = 1,
                            n_restarts        = 100,
                            seed              = 42) {
  set.seed(seed)
  if (is.data.frame(counts)) {
    if (!is.numeric(counts[[1]])) {
      rn <- as.character(counts[[1]]); counts <- as.matrix(counts[, -1, drop = FALSE])
      rownames(counts) <- rn
    } else counts <- as.matrix(counts)
  }
  storage.mode(counts) <- "double"
  labels <- colnames(counts); totals <- colSums(counts); n <- nrow(counts)
  if (!is.null(rare_threshold)) rare_labels <- union(rare_labels, labels[totals < rare_threshold])
  rare_labels <- intersect(rare_labels, labels); is_rare <- labels %in% rare_labels
  target_frac <- ifelse(is_rare, rare_test_frac, default_test_frac); names(target_frac) <- labels
  active <- totals > 0
  w <- ifelse(is_rare, rare_weight, 1) * active
  inv_totals <- ifelse(active, 1 / pmax(totals, 1), 0)
  rate <- counts * rep(inv_totals, each = n)
  rate_sq_w <- as.numeric(rate^2 %*% w)
  has_rare <- if (any(is_rare)) rowSums(counts[, is_rare, drop = FALSE]) > 0 else rep(FALSE, n)
  init_prob <- ifelse(has_rare, rare_test_frac, default_test_frac)
  best_score <- Inf; best_in_test <- logical(n)
  for (restart in seq_len(n_restarts)) {
    in_test <- runif(n) < init_prob
    test_sum <- colSums(counts[in_test, , drop = FALSE])
    repeat {
      cur_err <- test_sum * inv_totals - target_frac
      cur_score <- sum(w * cur_err * cur_err)
      sgn <- 1 - 2 * in_test
      linear <- as.numeric(rate %*% (w * cur_err))
      flip_delta <- 2 * sgn * linear + rate_sq_w
      best_i <- which.min(flip_delta)
      if (flip_delta[best_i] >= -1e-12) break
      test_sum <- test_sum + sgn[best_i] * counts[best_i, ]
      in_test[best_i] <- !in_test[best_i]
    }
    if (cur_score < best_score) { best_score <- cur_score; best_in_test <- in_test }
  }
  list(test_rownames  = rownames(counts)[best_in_test],
       train_rownames = rownames(counts)[!best_in_test],
       score = best_score)
}

# Image-grouped CV folds (validation indices), as in the snapshot.
make_folds <- function(train_sce, k = 5) {
  tr_idx  <- caret::groupKFold(factor(train_sce$sample_id), k = k)
  all_idx <- seq_len(ncol(train_sce))
  lapply(tr_idx, function(idx) setdiff(all_idx, idx))
}

# Reshape a softprob vector into a cells x classes matrix.
to_prob_matrix <- function(p, num_classes) {
  if (is.list(p)) p <- do.call(rbind, p)
  if (is.null(dim(p))) {
    byrow_m <- matrix(p, ncol = num_classes, byrow = TRUE)
    bycol_m <- matrix(p, ncol = num_classes, byrow = FALSE)
    p <- if (mean(abs(rowSums(byrow_m) - 1)) <= mean(abs(rowSums(bycol_m) - 1))) byrow_m else bycol_m
  }
  as.matrix(p)
}

# Legacy inverse-frequency weights (undefined/10, Igg/30). Computed over PRESENT classes only
# (factor of the actual labels) so an absent class never produces an Inf weight.
legacy_weights <- function(labels_chr) {
  yl <- factor(labels_chr)
  w_class <- 1 / (table(yl) / length(yl))
  if ("undefined" %in% names(w_class)) w_class["undefined"] <- w_class["undefined"] / 10
  if ("Igg" %in% names(w_class)) w_class["Igg"] <- w_class["Igg"] / 30
  w <- as.numeric(w_class[as.character(yl)])
  w / mean(w)
}

# Build a weighted training DMatrix. Integer labels are coded against the FIXED global mapping.
build_dmatrix <- function(X, labels_chr, label_mapping, weighted = TRUE) {
  y_int <- as.integer(factor(labels_chr, levels = label_mapping)) - 1
  if (weighted) {
    xgboost::xgb.DMatrix(X, label = y_int, weight = legacy_weights(labels_chr))
  } else {
    xgboost::xgb.DMatrix(X, label = y_int)
  }
}

# Metrics: overall + macro (the honest headline) + critical-class recall.
extract_metrics <- function(pred, truth, label_mapping, critical, suffix) {
  cm <- caret::confusionMatrix(
    data      = factor(as.character(pred),  levels = label_mapping),
    reference = factor(as.character(truth), levels = label_mapping))
  ov <- cm$overall; bc <- as.data.frame(cm$byClass)
  rownames(bc) <- sub("Class: ", "", rownames(bc))
  crit <- intersect(critical, rownames(bc))
  out <- list(
    Accuracy                  = unname(ov["Accuracy"]),
    Kappa                     = unname(ov["Kappa"]),
    AccuracyNull              = unname(ov["AccuracyNull"]),
    BalancedAccuracy_macro    = mean(bc[["Balanced Accuracy"]], na.rm = TRUE),
    F1_macro                  = mean(bc[["F1"]],               na.rm = TRUE),
    Sensitivity_macro         = mean(bc[["Sensitivity"]],      na.rm = TRUE),
    Sensitivity_critical_mean = mean(bc[crit, "Sensitivity"],  na.rm = TRUE),
    F1_critical_mean          = mean(bc[crit, "F1"],           na.rm = TRUE))
  stats::setNames(out, paste0(names(out), "_", suffix))
}

# ===========================================================================
# 3. OBSERVED baseline — fix the split, folds, tuned params; record real metrics.
# ===========================================================================
dat <- as.data.frame(unclass(table(lab_sce$sample_id, lab_sce$cell_labels)))

res <- split_test_rows(
  dat,
  default_test_frac = 0.20,
  rare_test_frac    = 0.20,
  rare_labels       = c("B", "BnT", "MDSC", "Igg", "mregDC", "NK"),
  rare_weight       = 20,
  n_restarts        = 1000,
  seed              = 42
)

train_cells <- which(lab_sce$sample_id %in% res$train_rownames)
test_cells   <- which(lab_sce$sample_id %in% res$test_rownames)
train_sce <- lab_sce[, train_cells]
test_sce  <- lab_sce[, test_cells]

# Feature matrices are invariant under label permutation, so build them ONCE.
X_train <- t(assay(train_sce, "exprs")[marker_keep, ])
X_test  <- t(assay(test_sce,  "exprs")[marker_keep, ])

label_mapping <- levels(factor(lab_sce$cell_labels))
num_classes   <- length(label_mapping)

fold_list <- make_folds(train_sce, k = 5)

base_params <- list(booster = "gbtree", objective = "multi:softprob",
                    eval_metric = "mlogloss", num_class = num_classes,
                    tree_method = "hist", nthread = nthread)

# 3a. 12-combo grid, pick params by min CV mlogloss (exactly as the 0.994 snapshot).
param_grid <- expand.grid(max_depth = c(3, 4, 5), eta = c(0.05, 0.1),
                          subsample = 0.8, colsample_bytree = 0.8,
                          min_child_weight = c(1, 5))
real_labels_train <- as.character(train_sce$cell_labels)
Xy_real <- build_dmatrix(X_train, real_labels_train, label_mapping, weighted = TRUE)

message("Tuning on real labels (", nrow(param_grid), " combos) ...")
grid_summary <- purrr::map_dfr(seq_len(nrow(param_grid)), function(i) {
  p <- c(base_params, as.list(param_grid[i, ]))
  cv <- xgboost::xgb.cv(params = p, data = Xy_real, folds = fold_list,
                        nrounds = 300, early_stopping_rounds = 20,
                        verbose = FALSE, prediction = FALSE)
  it <- which.min(cv$evaluation_log$test_mlogloss_mean)
  data.frame(param_grid[i, ], best_iteration = it,
             best_mlogloss = cv$evaluation_log$test_mlogloss_mean[it])
})
best_row    <- grid_summary[which.min(grid_summary$best_mlogloss), ]
best_nrounds <- as.integer(best_row$best_iteration)
params_best <- c(base_params, list(
  max_depth = as.integer(best_row$max_depth), eta = best_row$eta,
  subsample = best_row$subsample, colsample_bytree = best_row$colsample_bytree,
  min_child_weight = as.integer(best_row$min_child_weight)))
message(sprintf("Selected: max_depth=%d eta=%.2f mcw=%d nrounds=%d (mlogloss=%.4f)",
                params_best$max_depth, params_best$eta, params_best$min_child_weight,
                best_nrounds, best_row$best_mlogloss))

# 3b. Observed metrics with the fixed estimator (one row per real fit + test).
real_labels_test <- as.character(test_sce$cell_labels)

# Refit the fixed estimator for one labelling and return CV (optional) + test metrics.
# predict() accepts the plain feature matrix directly, so no test DMatrix is needed.
fit_once <- function(labels_train, labels_test) {
  Xy <- build_dmatrix(X_train, labels_train, label_mapping, weighted = TRUE)
  out <- list()
  if (DO_CV) {
    cv <- xgboost::xgb.cv(params = params_best, data = Xy, folds = fold_list,
                          nrounds = best_nrounds, verbose = FALSE, prediction = TRUE)
    cv_prob <- to_prob_matrix(cv$pred %||% cv$cv_predict, num_classes)
    cv_pred <- label_mapping[max.col(cv_prob, ties.method = "first")]
    out <- c(out, extract_metrics(cv_pred, labels_train, label_mapping, critical_classes, "cv"))
  }
  fm <- xgboost::xgb.train(params = params_best, data = Xy, nrounds = best_nrounds, verbose = FALSE)
  test_prob <- to_prob_matrix(predict(fm, X_test), num_classes)
  test_pred <- label_mapping[max.col(test_prob, ties.method = "first")]
  out <- c(out, extract_metrics(test_pred, labels_test, label_mapping, critical_classes, "test"))
  out
}

message("Fitting observed model ...")
observed <- fit_once(real_labels_train, real_labels_test)
print(as.data.frame(observed))

# ===========================================================================
# 4. Permutation loop — shuffle labels, refit fixed estimator, record metrics.
# ===========================================================================
all_labels <- as.character(lab_sce$cell_labels)
image_id   <- as.character(lab_sce$sample_id)

shuffle_labels <- function(seed) {
  set.seed(seed)
  if (null_type == "within_image") {
    # Stricter null: shuffle only within each image (preserves per-image label composition).
    out <- all_labels
    for (img in unique(image_id)) {
      idx <- which(image_id == img)
      out[idx] <- sample(all_labels[idx])
    }
    out
  } else {
    sample(all_labels)  # global shuffle: preserves overall class frequencies only
  }
}

message("Running ", N_PERM, " permutations ...")
null_rows <- purrr::map_dfr(seq_len(N_PERM), function(i) {
  shuffled <- shuffle_labels(base_seed + i)
  m <- fit_once(shuffled[train_cells], shuffled[test_cells])
  if (i %% 10 == 0) message("  perm ", i, "/", N_PERM,
                            sprintf("  (test acc=%.3f kappa=%.3f f1macro=%.3f)",
                                    m$Accuracy_test, m$Kappa_test, m$F1_macro_test))
  data.frame(perm = i, m, check.names = FALSE)
})

# ===========================================================================
# 5. Empirical one-sided p-values (higher is better).
# ===========================================================================
metric_cols <- setdiff(names(observed), character(0))
pvalues <- purrr::map_dfr(metric_cols, function(mc) {
  obs <- observed[[mc]]
  nul <- null_rows[[mc]]
  data.frame(metric    = mc,
             observed  = round(obs, 4),
             null_mean = round(mean(nul, na.rm = TRUE), 4),
             null_sd   = round(sd(nul, na.rm = TRUE), 4),
             p_value   = (1 + sum(nul >= obs, na.rm = TRUE)) / (N_PERM + 1))
})
print(pvalues)

# ===========================================================================
# 6. Persist (no-clobber names).
# ===========================================================================
results_long <- dplyr::bind_rows(
  data.frame(perm = 0L, as.data.frame(observed), check.names = FALSE),  # observed row
  null_rows
)
results_long$is_observed <- results_long$perm == 0L

config <- list(N_PERM = N_PERM, base_seed = base_seed, DO_CV = DO_CV, null_type = null_type,
               params_best = params_best, best_nrounds = best_nrounds,
               n_train_cells = ncol(train_sce), n_test_cells = ncol(test_sce),
               num_classes = num_classes, label_mapping = label_mapping,
               critical_classes = critical_classes,
               train_images = res$train_rownames, test_images = res$test_rownames)

write.csv(results_long, file.path(out_dir, "permutation_results.csv"), row.names = FALSE)
saveRDS(list(null = null_rows, observed = observed, pvalues = pvalues, config = config),
        file.path(out_dir, "permutation_null.rds"))

message("Wrote permutation_results.csv (", nrow(results_long), " rows) and permutation_null.rds to ",
        out_dir)
