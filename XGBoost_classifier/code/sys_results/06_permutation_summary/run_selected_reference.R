# =============================================================================
# run_selected_reference.R
# Reproduces the selected reference version: adaptive
# label-weight split + 12-model CV grid, evaluated with the SAME metric set as
# the permutation CSVs, so it can be used as the reference baseline in the
# permutation comparison report.
# =============================================================================
suppressMessages({
  library(SingleCellExperiment); library(xgboost); library(caret); library(dplyr)
})
set.seed(42)

base_path <- "/Users/svitlanamoiseyenko/Repos/genomic/xgboost-imc-panel2-cell-classifier/XGBoost_classifier/Panel_2_10"
out_dir   <- file.path(base_path, "Rout")
report_dir <- "/Users/svitlanamoiseyenko/Repos/genomic/xgboost-imc-panel2-cell-classifier/XGBoost_classifier/code/sys_results/06_permutation_summary"
critical <- c("MDSC", "BnT", "Igg", "mregDC", "PMN_MDSC", "Neutrophil")

# ---- split_test_rows with label_weights support (selected version) -----------
split_test_rows <- function(counts, default_test_frac = 0.20, rare_test_frac = 0.40,
                            rare_labels = c("MDSC"), rare_threshold = NULL,
                            rare_weight = 1, label_weights = NULL,
                            n_restarts = 100, seed = 42) {
  set.seed(seed)
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  labels <- colnames(counts); totals <- colSums(counts); n <- nrow(counts)
  if (!is.null(rare_threshold)) rare_labels <- union(rare_labels, labels[totals < rare_threshold])
  rare_labels <- intersect(rare_labels, labels); is_rare <- labels %in% rare_labels
  target_frac <- ifelse(is_rare, rare_test_frac, default_test_frac); names(target_frac) <- labels
  active <- totals > 0
  if (!is.null(label_weights)) {
    w <- label_weights[labels]; unnamed <- is.na(w)
    w[unnamed] <- ifelse(is_rare[unnamed], rare_weight, 1); w <- w * active
  } else w <- ifelse(is_rare, rare_weight, 1) * active
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
  list(test_rownames = rownames(counts)[best_in_test],
       train_rownames = rownames(counts)[!best_in_test],
       summary = data.frame(label = labels, total = totals, test_count = test_sum,
                            achieved_test_frac = round(achieved, 3), is_rare = is_rare, row.names = NULL),
       score = best_score)
}

# ---- data + cleanup ---------------------------------------------------------
lab_sce <- readRDS(file.path(out_dir, "lab_sce.rds"))
em <- assay(lab_sce, "exprs")
lab_sce <- lab_sce[, !(lab_sce$cell_labels == "Igg" & em["IgG", ] < 1)]
em <- assay(lab_sce, "exprs")
lab_sce <- lab_sce[, !(lab_sce$cell_labels == "NK" & em["CD7", ] < 0.5)]
em <- assay(lab_sce, "exprs")
mreg <- lab_sce$cell_labels == "mregDC" & (em["panCK", ] > 0.5 | em["Ecad", ] > 0.8)
lab_sce$cell_labels[mreg] <- "Tumor"

dat <- as.matrix(unclass(table(lab_sce$sample_id, lab_sce$cell_labels)))

# ---- Selected adaptive split: auto-tune label_weights, then seed search --------
all_cls <- setdiff(unique(as.character(lab_sce$cell_labels)), c(NA, "unlabelled"))
lw <- setNames(rep(1, length(all_cls)), all_cls)
for (iter in 1:8) {
  s <- split_test_rows(dat, 0.20, 0.20, rare_threshold = 150, label_weights = lw,
                       n_restarts = 100, seed = 42)$summary
  s <- s[s$label %in% names(lw), ]
  zero_cls <- s$label[!is.na(s$achieved_test_frac) & s$achieved_test_frac == 0]
  over_cls <- s$label[!is.na(s$achieved_test_frac) & s$achieved_test_frac > 0.50]
  if (length(zero_cls) == 0) break
  lw[zero_cls] <- lw[zero_cls] * 4; lw[over_cls] <- pmax(lw[over_cls] / 2, 1)
}
best_seed <- 1; best_score <- Inf
for (s in 1:200) {
  r <- split_test_rows(dat, 0.20, 0.20, rare_threshold = 150, label_weights = lw,
                       n_restarts = 50, seed = s)
  if (r$score < best_score) { best_score <- r$score; best_seed <- s }
}
res <- split_test_rows(dat, 0.20, 0.20, rare_threshold = 150, label_weights = lw,
                       n_restarts = 500, seed = best_seed)
cat("Selected best split seed:", best_seed, "\n")

test_sce  <- lab_sce[, lab_sce$sample_id %in% res$test_rownames]
train_sce <- lab_sce[, lab_sce$sample_id %in% res$train_rownames]
train_sce$cell_id <- paste0(train_sce$sample_id, "_", train_sce$ObjectNumber)
test_sce$cell_id  <- paste0(test_sce$sample_id, "_", test_sce$ObjectNumber)
colnames(train_sce) <- train_sce$cell_id; colnames(test_sce) <- test_sce$cell_id

# ---- features / labels / weights -------------------------------------------
mk <- !grepl("DNA|Histone", rownames(train_sce))
X <- t(assay(train_sce, "exprs")[mk, ]); y_true <- factor(train_sce$cell_labels)
label_mapping <- levels(y_true); y <- as.integer(y_true) - 1; num_classes <- length(label_mapping)
keep <- !is.na(test_sce$cell_labels) & test_sce$cell_labels != "unlabelled" & test_sce$cell_labels %in% label_mapping
X_test <- t(assay(test_sce, "exprs")[mk, keep])
y_test_true <- factor(test_sce$cell_labels[keep], levels = label_mapping); y_test <- as.integer(y_test_true) - 1
w_class <- 1 / (table(y_true) / length(y_true))
w_class["undefined"] <- w_class["undefined"] / 10; w_class["Igg"] <- w_class["Igg"] / 30
w <- as.numeric(w_class[y_true]); w <- w / mean(w)
Xy <- xgb.DMatrix(X, label = y, weight = w); Xy_test <- xgb.DMatrix(X_test, label = y_test)

# ---- Selected 12-model CV grid ----------------------------------------------
cv_train_index <- caret::groupKFold(factor(train_sce$sample_id), k = 5)
all_idx <- seq_len(ncol(train_sce))
folds <- lapply(cv_train_index, function(tr) setdiff(all_idx, tr))
grid <- expand.grid(max_depth = c(3, 4, 5), eta = c(0.05, 0.1),
                    subsample = 0.8, colsample_bytree = 0.8, min_child_weight = c(1, 5))
cvres <- lapply(seq_len(nrow(grid)), function(i) {
  p <- list(booster = "gbtree", objective = "multi:softprob", eval_metric = "mlogloss",
            num_class = num_classes, max_depth = grid$max_depth[i], eta = grid$eta[i],
            subsample = grid$subsample[i], colsample_bytree = grid$colsample_bytree[i],
            min_child_weight = grid$min_child_weight[i], tree_method = "hist",
            nthread = max(1, parallel::detectCores() - 1))
  cv <- xgb.cv(p, Xy, folds = folds, nrounds = 300, early_stopping_rounds = 20,
               verbose = FALSE, prediction = FALSE)
  bi <- which.min(cv$evaluation_log$test_mlogloss_mean)
  data.frame(grid[i, ], best_iteration = bi, best_mlogloss = cv$evaluation_log$test_mlogloss_mean[bi])
})
cv_summary <- bind_rows(cvres) |> arrange(best_mlogloss)
bp <- cv_summary[1, ]
cat("Selected best config: depth", bp$max_depth, "eta", bp$eta, "mcw", bp$min_child_weight,
    "nrounds", bp$best_iteration, "\n")

# ---- final model + evaluation ----------------------------------------------
params_best <- list(booster = "gbtree", objective = "multi:softprob", eval_metric = "mlogloss",
                    num_class = num_classes, max_depth = as.integer(bp$max_depth), eta = bp$eta,
                    subsample = bp$subsample, colsample_bytree = bp$colsample_bytree,
                    min_child_weight = bp$min_child_weight, tree_method = "hist",
                    nthread = max(1, parallel::detectCores() - 1))
fm <- xgb.train(params_best, Xy, nrounds = as.integer(bp$best_iteration), verbose = FALSE)
prob <- matrix(predict(fm, Xy_test), ncol = num_classes, byrow = TRUE)
if (mean(abs(rowSums(prob) - 1)) > 1e-6) prob <- matrix(predict(fm, Xy_test), ncol = num_classes, byrow = FALSE)
pred <- factor(label_mapping[max.col(prob)], levels = label_mapping)
cm <- caret::confusionMatrix(pred, y_test_true)

bc <- cm$byClass; rownames(bc) <- sub("Class: ", "", rownames(bc))
crit_sens <- bc[critical, "Sensitivity"]
alona <- list(
  source = "Selected reference version",
  split = list(default_test_frac = 0.20, rare_test_frac = 0.20, rare_threshold = 150,
               scheme = "auto-tuned label_weights + 200-seed search", best_seed = best_seed,
               n_train_cells = ncol(train_sce), n_test_cells = sum(keep),
               n_test_images = length(res$test_rownames),
               min_test_count_observed = min(res$summary$test_count[res$summary$total > 0])),
  config = list(grid = "depth{3,4,5} x eta{.05,.1} x mcw{1,5}, sub/col=0.8 (no gamma/lambda/alpha/mds)",
                max_depth = bp$max_depth, eta = bp$eta, min_child_weight = bp$min_child_weight,
                nrounds_best = bp$best_iteration),
  metrics = data.frame(
    Accuracy_test = cm$overall["Accuracy"], Kappa_test = cm$overall["Kappa"],
    BalancedAccuracy_macro_test = mean(bc[, "Balanced Accuracy"], na.rm = TRUE),
    Sensitivity_critical_mean_test = mean(crit_sens, na.rm = TRUE), row.names = NULL),
  crit_sens = crit_sens, cv_summary = cv_summary, split_summary = res$summary
)
saveRDS(alona, file.path(report_dir, "selected_reference.rds"))
cat("\n==== SELECTED REFERENCE ====\n"); print(alona$metrics)
cat("critical sens:\n"); print(round(crit_sens, 3))
message("SELECTED_DONE")
