# =============================================================================
# 09_lightgbm_comparison_build.R
# Trains XGBoost + Random Forest + LightGBM on the SAME image-level train/test
# split and saves all metrics/artifacts for the comparison report.
# Faithful re-implementation of "01_read_data kos.Rmd" (LightGBM fork).
# =============================================================================

suppressMessages({
  library(SingleCellExperiment)
  library(xgboost)
  library(randomForest)
  library(lightgbm)
  library(caret)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

base_path  <- normalizePath(file.path(dirname(dirname(dirname(getwd()))), "Panel_2_10"),
                            mustWork = FALSE)
# Robust fallback to the known repo location.
if (!dir.exists(base_path)) {
  base_path <- "/Users/svitlanamoiseyenko/Repos/genomic/xgboost-imc-panel2-cell-classifier/XGBoost_classifier/Panel_2_10"
}
out_dir    <- file.path(base_path, "Rout")
report_dir <- "/Users/svitlanamoiseyenko/Repos/genomic/xgboost-imc-panel2-cell-classifier/XGBoost_classifier/code/results/lightgbt"
message("base_path: ", base_path)

# -----------------------------------------------------------------------------
# 1. Load cleaned labelled cells, re-apply biological cleanup rules (idempotent)
# -----------------------------------------------------------------------------
lab_sce <- readRDS(file.path(out_dir, "lab_sce.rds"))

exprs_mat <- assay(lab_sce, "exprs")
low_igg     <- lab_sce$cell_labels == "Igg"    & exprs_mat["IgG", ] < 1
lab_sce     <- lab_sce[, !low_igg]
exprs_mat   <- assay(lab_sce, "exprs")
low_cd7_nk  <- lab_sce$cell_labels == "NK"     & exprs_mat["CD7", ] < 0.5
lab_sce     <- lab_sce[, !low_cd7_nk]
exprs_mat   <- assay(lab_sce, "exprs")
mreg_to_tum <- lab_sce$cell_labels == "mregDC" &
  (exprs_mat["panCK", ] > 0.5 | exprs_mat["Ecad", ] > 0.8)
lab_sce$cell_labels[mreg_to_tum] <- "Tumor"

# -----------------------------------------------------------------------------
# 2. Reproduce the image-level train/test split (split_test_rows, kos.Rmd)
# -----------------------------------------------------------------------------
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
      rn <- as.character(counts[[1]]); counts <- as.matrix(counts[, -1, drop = FALSE]); rownames(counts) <- rn
    } else counts <- as.matrix(counts)
  }
  storage.mode(counts) <- "double"
  labels <- colnames(counts); totals <- colSums(counts); n <- nrow(counts)
  if (!is.null(rare_threshold)) rare_labels <- union(rare_labels, labels[totals < rare_threshold])
  rare_labels <- intersect(rare_labels, labels); is_rare <- labels %in% rare_labels
  target_frac <- ifelse(is_rare, rare_test_frac, default_test_frac); names(target_frac) <- labels
  active <- totals > 0; w <- ifelse(is_rare, rare_weight, 1) * active
  inv_totals <- ifelse(active, 1 / pmax(totals, 1), 0)
  rate <- counts * rep(inv_totals, each = n); rate_sq_w <- as.numeric(rate^2 %*% w)
  has_rare <- if (any(is_rare)) rowSums(counts[, is_rare, drop = FALSE]) > 0 else rep(FALSE, n)
  init_prob <- ifelse(has_rare, rare_test_frac, default_test_frac)
  best_score <- Inf; best_in_test <- logical(n)
  for (restart in seq_len(n_restarts)) {
    in_test <- runif(n) < init_prob; test_sum <- colSums(counts[in_test, , drop = FALSE])
    repeat {
      cur_err <- test_sum * inv_totals - target_frac; cur_score <- sum(w * cur_err * cur_err)
      sgn <- 1 - 2 * in_test; linear <- as.numeric(rate %*% (w * cur_err))
      flip_delta <- 2 * sgn * linear + rate_sq_w; best_i <- which.min(flip_delta)
      if (flip_delta[best_i] >= -1e-12) break
      test_sum <- test_sum + sgn[best_i] * counts[best_i, ]; in_test[best_i] <- !in_test[best_i]
    }
    if (cur_score < best_score) { best_score <- cur_score; best_in_test <- in_test }
  }
  test_sum <- colSums(counts[best_in_test, , drop = FALSE])
  achieved <- ifelse(active, test_sum * inv_totals, NA_real_)
  summary_df <- data.frame(label = labels, total = totals, test_count = test_sum,
                           train_count = totals - test_sum, target_test_frac = target_frac,
                           achieved_test_frac = round(achieved, 3), is_rare = is_rare, row.names = NULL)
  list(test_idx = which(best_in_test), train_idx = which(!best_in_test),
       test_rownames = rownames(counts)[best_in_test], train_rownames = rownames(counts)[!best_in_test],
       summary = summary_df, score = best_score)
}

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
split_summary <- res$summary

test_sce  <- lab_sce[, lab_sce$sample_id %in% res$test_rownames]
train_sce <- lab_sce[, lab_sce$sample_id %in% res$train_rownames]

train_sce$cell_id <- paste0(train_sce$sample_id, "_", train_sce$ObjectNumber)
test_sce$cell_id  <- paste0(test_sce$sample_id, "_", test_sce$ObjectNumber)
colnames(train_sce) <- train_sce$cell_id
colnames(test_sce)  <- test_sce$cell_id

# -----------------------------------------------------------------------------
# 3. Feature matrices, labels, class weights (identical to kos.Rmd)
# -----------------------------------------------------------------------------
marker_keep <- !grepl("DNA|Histone", rownames(train_sce))

X      <- t(assay(train_sce, "exprs")[marker_keep, ])
y_true <- factor(train_sce$cell_labels)
label_mapping <- levels(y_true)
y      <- as.integer(y_true) - 1
num_classes <- length(label_mapping)

test_eval_keep <- !is.na(test_sce$cell_labels) &
  test_sce$cell_labels != "unlabelled" &
  test_sce$cell_labels %in% label_mapping
X_test      <- t(assay(test_sce, "exprs")[marker_keep, test_eval_keep])
y_test_true <- factor(test_sce$cell_labels[test_eval_keep], levels = label_mapping)
y_test      <- as.integer(y_test_true) - 1
stopifnot(identical(colnames(X), colnames(X_test)))

w_class <- 1 / (table(y_true) / length(y_true))
if ("undefined" %in% names(w_class)) w_class["undefined"] <- w_class["undefined"] / 10
if ("Igg"       %in% names(w_class)) w_class["Igg"]       <- w_class["Igg"]       / 30
w <- as.numeric(w_class[y_true]); w <- w / mean(w)

split_stats <- data.frame(
  cell_label = label_mapping,
  train_n = as.integer(table(factor(as.character(train_sce$cell_labels), levels = label_mapping))),
  test_n  = as.integer(table(y_test_true))
)

# -----------------------------------------------------------------------------
# 4. XGBoost final model (documented best params, nrounds = 108)
# -----------------------------------------------------------------------------
message("Training XGBoost ...")
Xy      <- xgb.DMatrix(data = X, label = y, weight = w)
Xy_test <- xgb.DMatrix(data = X_test, label = y_test)

params_best <- list(
  booster = "gbtree", objective = "multi:softprob", eval_metric = "mlogloss",
  num_class = num_classes, max_depth = 3, eta = 0.1, subsample = 0.8,
  colsample_bytree = 0.8, min_child_weight = 5, tree_method = "hist",
  nthread = max(1, parallel::detectCores() - 1)
)
final_model <- xgb.train(params = params_best, data = Xy, nrounds = 108, verbose = FALSE)

xgb_raw  <- predict(final_model, Xy_test)
xgb_prob <- matrix(xgb_raw, ncol = num_classes, byrow = TRUE)
if (mean(abs(rowSums(xgb_prob) - 1)) > 1e-6)
  xgb_prob <- matrix(xgb_raw, ncol = num_classes, byrow = FALSE)
colnames(xgb_prob) <- label_mapping
xgb_pred <- factor(label_mapping[max.col(xgb_prob)], levels = label_mapping)
cm_test  <- caret::confusionMatrix(data = xgb_pred, reference = y_test_true)

# -----------------------------------------------------------------------------
# 5. Random Forest
# -----------------------------------------------------------------------------
message("Training Random Forest ...")
rf_model <- randomForest(x = X, y = y_true, ntree = 300, mtry = 8,
                         classwt = w_class, importance = TRUE)
rf_pred  <- predict(rf_model, newdata = X_test)
cm_rf    <- caret::confusionMatrix(data = rf_pred, reference = y_test_true)

# -----------------------------------------------------------------------------
# 6. LightGBM
# -----------------------------------------------------------------------------
message("Training LightGBM ...")
dtrain <- lgb.Dataset(data = X, label = as.integer(y_true) - 1,
                      weight = as.numeric(w_class[as.character(y_true)]))
dtest  <- lgb.Dataset(data = X_test, label = y_test, reference = dtrain)

lgb_params <- list(
  objective = "multiclass", num_class = num_classes, num_leaves = 31, max_depth = 4,
  learning_rate = 0.1, feature_fraction = 0.8, bagging_fraction = 0.8, bagging_freq = 5,
  min_child_weight = 5, lambda_l2 = 1.0, lambda_l1 = 0.0, seed = 42, num_threads = 0
)
lgb_model <- lgb.train(params = lgb_params, data = dtrain, nrounds = 500,
                       valids = list(test = dtest), early_stopping_rounds = 20, verbose = -1)
lgb_best_iter  <- lgb_model$best_iter
lgb_best_score <- lgb_model$best_score   # multiclass logloss on test (early-stopping metric)

lgb_pred_out <- predict(lgb_model, X_test, num_iteration = lgb_best_iter)
lgb_prob <- if (is.matrix(lgb_pred_out) && ncol(lgb_pred_out) == num_classes) {
  lgb_pred_out
} else {
  matrix(lgb_pred_out, ncol = num_classes, byrow = TRUE)
}
colnames(lgb_prob) <- label_mapping
rownames(lgb_prob) <- colnames(test_sce)[test_eval_keep]
lgb_pred <- factor(label_mapping[max.col(lgb_prob)], levels = label_mapping)
cm_lgb   <- caret::confusionMatrix(data = lgb_pred, reference = y_test_true)

# Test-set multiclass logloss for XGB & RF (comparable to LightGBM best_score)
clip <- function(p) pmin(pmax(p, 1e-15), 1 - 1e-15)
true_idx <- as.integer(y_test_true)
mlogloss <- function(prob) -mean(log(clip(prob[cbind(seq_len(nrow(prob)), true_idx)])))
rf_prob_raw <- predict(rf_model, newdata = X_test, type = "prob")
rf_prob <- matrix(0, nrow = nrow(X_test), ncol = num_classes, dimnames = list(NULL, label_mapping))
rf_prob[, colnames(rf_prob_raw)] <- rf_prob_raw

# -----------------------------------------------------------------------------
# 7. Metrics table, consensus, entropy
# -----------------------------------------------------------------------------
metrics <- data.frame(
  Model        = c("XGBoost", "Random Forest", "LightGBM"),
  Accuracy     = c(cm_test$overall["Accuracy"], cm_rf$overall["Accuracy"], cm_lgb$overall["Accuracy"]),
  Kappa        = c(cm_test$overall["Kappa"],    cm_rf$overall["Kappa"],    cm_lgb$overall["Kappa"]),
  Test_mlogloss = c(mlogloss(xgb_prob), mlogloss(rf_prob), as.numeric(lgb_best_score)),
  row.names = NULL
)

get_majority <- function(a, b, c) {
  tbl <- table(c(as.character(a), as.character(b), as.character(c)))
  if (max(tbl) < 2) return("undefined")
  names(which.max(tbl))
}
consensus3 <- data.frame(
  cell_id = colnames(test_sce)[test_eval_keep],
  XGBoost = xgb_pred, RandomForest = rf_pred, LightGBM = lgb_pred,
  True_Label = y_test_true, stringsAsFactors = FALSE
)
consensus3$Majority_Label <- mapply(get_majority, consensus3$XGBoost, consensus3$RandomForest, consensus3$LightGBM)
consensus3$full_agreement <- consensus3$XGBoost == consensus3$RandomForest &
                             consensus3$RandomForest == consensus3$LightGBM

lgb_entropy <- -rowSums(lgb_prob * log(lgb_prob + 1e-10))
ent_thr <- quantile(lgb_entropy, 0.90)
consensus3$high_entropy <- lgb_entropy > ent_thr

consensus_stats <- list(
  full_agreement_pct  = mean(consensus3$full_agreement) * 100,
  majority_agree_pct  = mean(consensus3$Majority_Label != "undefined") * 100,
  majority_accuracy   = mean(consensus3$Majority_Label == as.character(consensus3$True_Label)),
  n_high_entropy      = sum(consensus3$high_entropy),
  high_entropy_full_agree_pct = mean(consensus3$full_agreement[consensus3$high_entropy]) * 100,
  low_entropy_full_agree_pct  = mean(consensus3$full_agreement[!consensus3$high_entropy]) * 100
)

# -----------------------------------------------------------------------------
# 8. Variable importance (normalized) across the three models
# -----------------------------------------------------------------------------
xgb_imp <- xgb.importance(feature_names = colnames(X), model = final_model) |>
  transmute(Feature, XGBoost = Gain / sum(Gain))
lgb_imp <- lgb.importance(lgb_model, percentage = TRUE) |>
  transmute(Feature, LightGBM = Gain / sum(Gain))
rf_imp_raw <- importance(rf_model)[, "MeanDecreaseGini"]
rf_imp <- data.frame(Feature = names(rf_imp_raw), RandomForest = rf_imp_raw / sum(rf_imp_raw),
                     row.names = NULL)

varimp <- xgb_imp |>
  full_join(lgb_imp, by = "Feature") |>
  full_join(rf_imp,  by = "Feature") |>
  mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0))) |>
  mutate(mean_importance = (XGBoost + LightGBM + RandomForest) / 3) |>
  arrange(desc(mean_importance))

# Per-class sensitivity comparison
byclass <- data.frame(
  class = sub("Class: ", "", rownames(cm_test$byClass)),
  XGBoost       = cm_test$byClass[, "Sensitivity"],
  RandomForest  = cm_rf$byClass[,  "Sensitivity"],
  LightGBM      = cm_lgb$byClass[, "Sensitivity"],
  row.names = NULL
)

# -----------------------------------------------------------------------------
# 9. Save artifacts
# -----------------------------------------------------------------------------
results <- list(
  metrics = metrics, split_summary = split_summary, split_stats = split_stats,
  cm_test = cm_test, cm_rf = cm_rf, cm_lgb = cm_lgb,
  byclass = byclass, varimp = varimp, consensus3 = consensus3,
  consensus_stats = consensus_stats, lgb_params = lgb_params, params_best = params_best,
  lgb_best_iter = lgb_best_iter, label_mapping = label_mapping,
  n_train = ncol(train_sce), n_test = sum(test_eval_keep),
  n_train_images = length(res$train_rownames), n_test_images = length(res$test_rownames)
)
saveRDS(results, file.path(out_dir, "lgb_comparison.rds"))
saveRDS(results, file.path(report_dir, "lgb_comparison.rds"))
write.csv(metrics, file.path(report_dir, "lgb_comparison_metrics.csv"), row.names = FALSE)

cat("\n==== METRICS ====\n"); print(metrics)
cat("\nLightGBM best iteration:", lgb_best_iter, "\n")
cat("Full 3-model agreement: ", round(consensus_stats$full_agreement_pct, 2), "%\n")
cat("Saved -> ", file.path(report_dir, "lgb_comparison.rds"), "\n")
message("BUILD_DONE")
