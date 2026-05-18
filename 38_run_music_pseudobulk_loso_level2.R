# scripts/38_run_music_pseudobulk_loso_level2.R

source("scripts/00_setup_paths.R")

required_pkgs <- c("MuSiC", "SingleCellExperiment", "Seurat", "SeuratObject", "data.table")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(MuSiC)
  library(SingleCellExperiment)
  library(Seurat)
  library(SeuratObject)
  library(data.table)
})

set.seed(1)

val_dir    <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso_level2")
seu_path   <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_level2.rds")
bulk_path  <- file.path(val_dir, "pseudobulk_counts_by_sample_level2.rds")
truth_path <- file.path(val_dir, "ground_truth_proportions_wide_level2.csv")
folds_path <- file.path(val_dir, "loso_folds_level2.csv")

stopifnot(dir.exists(val_dir))
stopifnot(file.exists(seu_path))
stopifnot(file.exists(bulk_path))
stopifnot(file.exists(truth_path))
stopifnot(file.exists(folds_path))

log_msg("Loading level2 scRNA object:", seu_path)
seu <- readRDS(seu_path)

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

stopifnot("cellcompartment_lv2" %in% colnames(seu@meta.data))
stopifnot("Sample" %in% colnames(seu@meta.data))

pb_counts <- readRDS(bulk_path)
truth <- data.table::fread(truth_path)
folds <- data.table::fread(folds_path)

celltypes <- c(
  "CD4_T",
  "NK_cells",
  "B_cells",
  "Plasma_cells",
  "Macrophage",
  "DC",
  "Epithelial",
  "Stromal",
  "Endothelial",
  "Other"
)

stopifnot(all(celltypes %in% colnames(truth)))
stopifnot("Sample" %in% colnames(truth))
stopifnot(all(c("test_sample", "train_samples", "n_train_samples") %in% colnames(folds)))

truth <- truth[match(colnames(pb_counts), truth$Sample)]
stopifnot(identical(colnames(pb_counts), truth$Sample))

meta <- seu@meta.data
meta$Sample <- as.character(meta$Sample)
meta$cellcompartment_lv2 <- factor(as.character(meta$cellcompartment_lv2), levels = celltypes)

if (any(is.na(meta$cellcompartment_lv2))) {
  stop("Missing values found in cellcompartment_lv2.")
}

pred_list <- list()
metrics_list <- list()

for (i in seq_len(nrow(folds))) {
  test_sample <- folds$test_sample[i]
  train_samples <- unlist(strsplit(folds$train_samples[i], ";", fixed = TRUE))
  
  log_msg("Running MuSiC level2 LOSO fold | test =", test_sample,
          "| train =", paste(train_samples, collapse = " | "))
  
  # -----------------------------
  # Subset scRNA training cells
  # -----------------------------
  keep_cells <- rownames(meta)[meta$Sample %in% train_samples]
  if (length(keep_cells) == 0) {
    stop("No training cells found for fold with test sample: ", test_sample)
  }
  
  seu_train <- subset(seu, cells = keep_cells)
  
  seu_train$Sample <- droplevels(factor(as.character(seu_train$Sample)))
  seu_train$cellcompartment_lv2 <- droplevels(
    factor(as.character(seu_train$cellcompartment_lv2), levels = celltypes)
  )
  
  train_ct_tab <- table(seu_train$cellcompartment_lv2)
  log_msg("Training cells per level2 compartment:",
          paste(names(train_ct_tab), as.integer(train_ct_tab), sep = "=", collapse = " | "))
  
  # -----------------------------
  # Build test bulk matrix
  # MuSiC bug workaround for single sample
  # -----------------------------
  stopifnot(test_sample %in% colnames(pb_counts))
  bulk_test <- pb_counts[, test_sample, drop = FALSE]
  bulk_test <- as.matrix(bulk_test)
  stopifnot(ncol(bulk_test) == 1)
  
  dup_name <- paste0(test_sample, "_dup")
  bulk_test_music <- cbind(bulk_test, bulk_test)
  colnames(bulk_test_music) <- c(test_sample, dup_name)
  
  # -----------------------------
  # Convert to SingleCellExperiment
  # -----------------------------
  sce_train <- as.SingleCellExperiment(seu_train)
  
  # -----------------------------
  # Run MuSiC
  # -----------------------------
  res <- music_prop(
    bulk.mtx = bulk_test_music,
    sc.sce = sce_train,
    clusters = "cellcompartment_lv2",
    samples = "Sample",
    verbose = FALSE
  )
  
  est <- as.data.frame(res$Est.prop.weighted)
  est$sample <- rownames(est)
  est$sample <- as.character(est$sample)
  
  # keep only the real held-out sample
  est <- est[est$sample == test_sample, , drop = FALSE]
  
  for (ct in celltypes) {
    if (!ct %in% colnames(est)) est[[ct]] <- 0
  }
  est <- est[, c("sample", celltypes), drop = FALSE]
  
  # -----------------------------
  # Truth for held-out sample
  # -----------------------------
  truth_row <- truth[truth$Sample == test_sample]
  if (nrow(truth_row) != 1) {
    stop("Expected exactly one truth row for sample: ", test_sample)
  }
  
  # -----------------------------
  # Save fold-level predictions
  # -----------------------------
  pred_long <- data.table::melt(
    as.data.table(est),
    id.vars = "sample",
    variable.name = "celltype",
    value.name = "estimated_proportion"
  )
  
  truth_long <- data.table::melt(
    truth_row[, c("Sample", celltypes), with = FALSE],
    id.vars = "Sample",
    variable.name = "celltype",
    value.name = "true_proportion"
  )
  setnames(truth_long, "Sample", "sample")
  
  fold_pred <- merge(pred_long, truth_long, by = c("sample", "celltype"), all.x = TRUE)
  fold_pred[, method := "MuSiC"]
  fold_pred[, fold_test_sample := test_sample]
  fold_pred[, abs_error := abs(estimated_proportion - true_proportion)]
  fold_pred[, sq_error := (estimated_proportion - true_proportion)^2]
  
  pred_list[[test_sample]] <- fold_pred
  
  fold_mae <- mean(fold_pred$abs_error, na.rm = TRUE)
  fold_rmse <- sqrt(mean(fold_pred$sq_error, na.rm = TRUE))
  
  metrics_list[[test_sample]] <- data.table(
    method = "MuSiC",
    test_sample = test_sample,
    n_train_samples = length(train_samples),
    n_train_cells = ncol(seu_train),
    fold_mae = fold_mae,
    fold_rmse = fold_rmse
  )
  
  out_fold_dir <- file.path(val_dir, "music", "folds")
  dir.create(out_fold_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, file.path(out_fold_dir, paste0("music_level2_result_", test_sample, ".rds")))
  
  log_msg("Completed MuSiC level2 fold | test =", test_sample,
          "| MAE =", round(fold_mae, 4),
          "| RMSE =", round(fold_rmse, 4))
}

# ----------------------------------------
# Combine all folds
# ----------------------------------------
pred_all <- rbindlist(pred_list, use.names = TRUE, fill = TRUE)
metrics_all <- rbindlist(metrics_list, use.names = TRUE, fill = TRUE)

pred_wide <- dcast(
  pred_all[, .(sample, celltype, estimated_proportion)],
  sample ~ celltype,
  value.var = "estimated_proportion"
)

pred_wide <- pred_wide[match(truth$Sample, pred_wide$sample)]
stopifnot(identical(pred_wide$sample, truth$Sample))

# ----------------------------------------
# Overall metrics
# ----------------------------------------
overall_mae <- pred_all[, mean(abs_error, na.rm = TRUE)]
overall_rmse <- sqrt(pred_all[, mean(sq_error, na.rm = TRUE)])

per_celltype <- pred_all[, .(
  mae = mean(abs_error, na.rm = TRUE),
  rmse = sqrt(mean(sq_error, na.rm = TRUE))
), by = .(method, celltype)]

overall_pearson <- suppressWarnings(cor(
  pred_all$true_proportion,
  pred_all$estimated_proportion,
  method = "pearson",
  use = "complete.obs"
))

overall_spearman <- suppressWarnings(cor(
  pred_all$true_proportion,
  pred_all$estimated_proportion,
  method = "spearman",
  use = "complete.obs"
))

overall_summary <- data.table(
  method = "MuSiC",
  overall_mae = overall_mae,
  overall_rmse = overall_rmse,
  overall_pearson = overall_pearson,
  overall_spearman = overall_spearman
)

# ----------------------------------------
# Save outputs
# ----------------------------------------
out_dir <- file.path(val_dir, "music")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data.table::fwrite(pred_all, file.path(out_dir, "music_level2_loso_predictions_long.csv"))
data.table::fwrite(pred_wide, file.path(out_dir, "music_level2_loso_predictions_wide.csv"))
data.table::fwrite(metrics_all, file.path(out_dir, "music_level2_loso_fold_metrics.csv"))
data.table::fwrite(per_celltype, file.path(out_dir, "music_level2_loso_metrics_by_celltype.csv"))
data.table::fwrite(overall_summary, file.path(out_dir, "music_level2_loso_overall_metrics.csv"))

log_msg("Saved:", file.path(out_dir, "music_level2_loso_predictions_long.csv"))
log_msg("Saved:", file.path(out_dir, "music_level2_loso_predictions_wide.csv"))
log_msg("Saved:", file.path(out_dir, "music_level2_loso_fold_metrics.csv"))
log_msg("Saved:", file.path(out_dir, "music_level2_loso_metrics_by_celltype.csv"))
log_msg("Saved:", file.path(out_dir, "music_level2_loso_overall_metrics.csv"))

log_msg("MuSiC level2 LOSO overall MAE =", round(overall_mae, 4),
        "| overall RMSE =", round(overall_rmse, 4),
        "| Pearson =", round(overall_pearson, 4),
        "| Spearman =", round(overall_spearman, 4))

print(overall_summary)
print(metrics_all)
print(per_celltype)