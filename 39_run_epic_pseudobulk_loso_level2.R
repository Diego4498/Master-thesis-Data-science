# scripts/39_run_epic_pseudobulk_loso_level2.R

source("scripts/00_setup_paths.R")

required_pkgs <- c("EPIC", "Seurat", "SeuratObject", "Matrix", "data.table")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(EPIC)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

set.seed(1)

val_dir    <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso_level2")
seu_path   <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_level2.rds")
bulk_path  <- file.path(val_dir, "pseudobulk_cpm_by_sample_level2.rds")
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

pb_cpm <- readRDS(bulk_path)
truth  <- data.table::fread(truth_path)
folds  <- data.table::fread(folds_path)

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

truth <- truth[match(colnames(pb_cpm), truth$Sample)]
stopifnot(identical(colnames(pb_cpm), truth$Sample))

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
  
  log_msg("Running EPIC level2 LOSO fold | test =", test_sample,
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
  
  missing_ct <- setdiff(celltypes, names(train_ct_tab)[train_ct_tab > 0])
  if (length(missing_ct) > 0) {
    log_msg("WARNING: missing level2 compartments in training fold:",
            paste(missing_ct, collapse = ", "))
  }
  
  # -----------------------------
  # Build EPIC reference from training cells
  # -----------------------------
  counts <- GetAssayData(seu_train, layer = "counts")
  groups <- factor(seu_train$cellcompartment_lv2, levels = celltypes)
  
  # CPM per cell
  cpm <- t(t(counts) / Matrix::colSums(counts)) * 1e6
  
  # Average CPM per compartment
  ref_epic <- sapply(levels(groups), function(g) {
    idx <- which(groups == g)
    if (length(idx) == 0) {
      rep(NA_real_, nrow(cpm))
    } else {
      Matrix::rowMeans(cpm[, idx, drop = FALSE])
    }
  })
  
  ref_epic <- as.matrix(ref_epic)
  rownames(ref_epic) <- rownames(cpm)
  
  # Gene-wise variance per compartment
  ref_var <- sapply(levels(groups), function(g) {
    idx <- which(groups == g)
    if (length(idx) > 1) {
      apply(cpm[, idx, drop = FALSE], 1, var)
    } else if (length(idx) == 1) {
      rep(0, nrow(cpm))
    } else {
      rep(NA_real_, nrow(cpm))
    }
  })
  
  ref_var <- as.matrix(ref_var)
  rownames(ref_var) <- rownames(cpm)
  
  # remove compartments absent in training fold
  keep_comp <- colnames(ref_epic)[colSums(!is.na(ref_epic)) == nrow(ref_epic)]
  ref_epic <- ref_epic[, keep_comp, drop = FALSE]
  ref_var  <- ref_var[, keep_comp, drop = FALSE]
  
  # -----------------------------
  # Build test bulk matrix
  # -----------------------------
  stopifnot(test_sample %in% colnames(pb_cpm))
  bulk_test <- pb_cpm[, test_sample, drop = FALSE]
  bulk_test <- as.matrix(bulk_test)
  stopifnot(ncol(bulk_test) == 1)
  
  # align genes
  common <- intersect(rownames(bulk_test), rownames(ref_epic))
  if (length(common) == 0) {
    stop("No overlapping genes between pseudobulk and EPIC reference for fold: ", test_sample)
  }
  
  bulk2 <- bulk_test[common, , drop = FALSE]
  ref2  <- ref_epic[common, , drop = FALSE]
  var2  <- ref_var[common, colnames(ref2), drop = FALSE]
  
  log_msg("Common genes:", length(common))
  log_msg("Bulk dims (genes x samples):", nrow(bulk2), "x", ncol(bulk2))
  log_msg("Ref dims  (genes x compartments):", nrow(ref2), "x", ncol(ref2))
  
  ref_list <- list(
    refProfiles = ref2,
    refProfiles.var = var2,
    sigGenes = rownames(ref2)
  )
  
  # -----------------------------
  # Run EPIC
  # -----------------------------
  res <- EPIC::EPIC(
    bulk = bulk2,
    reference = ref_list
  )
  
  est <- as.data.frame(res$cellFractions)
  est$sample <- rownames(est)
  est$sample <- as.character(est$sample)
  
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
  fold_pred[, method := "EPIC"]
  fold_pred[, fold_test_sample := test_sample]
  fold_pred[, abs_error := abs(estimated_proportion - true_proportion)]
  fold_pred[, sq_error := (estimated_proportion - true_proportion)^2]
  
  pred_list[[test_sample]] <- fold_pred
  
  fold_mae <- mean(fold_pred$abs_error, na.rm = TRUE)
  fold_rmse <- sqrt(mean(fold_pred$sq_error, na.rm = TRUE))
  
  metrics_list[[test_sample]] <- data.table(
    method = "EPIC",
    test_sample = test_sample,
    n_train_samples = length(train_samples),
    n_train_cells = ncol(seu_train),
    n_common_genes = length(common),
    fold_mae = fold_mae,
    fold_rmse = fold_rmse
  )
  
  # save raw EPIC result per fold
  out_fold_dir <- file.path(val_dir, "epic", "folds")
  dir.create(out_fold_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, file.path(out_fold_dir, paste0("epic_level2_result_", test_sample, ".rds")))
  
  if (!is.null(res$fit.gof)) {
    fg <- as.data.frame(res$fit.gof)
    fg$sample <- rownames(fg)
    data.table::fwrite(
      as.data.table(fg),
      file.path(out_fold_dir, paste0("epic_level2_fit_gof_", test_sample, ".csv"))
    )
  }
  
  log_msg("Completed EPIC level2 fold | test =", test_sample,
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
  method = "EPIC",
  overall_mae = overall_mae,
  overall_rmse = overall_rmse,
  overall_pearson = overall_pearson,
  overall_spearman = overall_spearman
)

# ----------------------------------------
# Save outputs
# ----------------------------------------
out_dir <- file.path(val_dir, "epic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data.table::fwrite(pred_all, file.path(out_dir, "epic_level2_loso_predictions_long.csv"))
data.table::fwrite(pred_wide, file.path(out_dir, "epic_level2_loso_predictions_wide.csv"))
data.table::fwrite(metrics_all, file.path(out_dir, "epic_level2_loso_fold_metrics.csv"))
data.table::fwrite(per_celltype, file.path(out_dir, "epic_level2_loso_metrics_by_celltype.csv"))
data.table::fwrite(overall_summary, file.path(out_dir, "epic_level2_loso_overall_metrics.csv"))

log_msg("Saved:", file.path(out_dir, "epic_level2_loso_predictions_long.csv"))
log_msg("Saved:", file.path(out_dir, "epic_level2_loso_predictions_wide.csv"))
log_msg("Saved:", file.path(out_dir, "epic_level2_loso_fold_metrics.csv"))
log_msg("Saved:", file.path(out_dir, "epic_level2_loso_metrics_by_celltype.csv"))
log_msg("Saved:", file.path(out_dir, "epic_level2_loso_overall_metrics.csv"))

log_msg("EPIC level2 LOSO overall MAE =", round(overall_mae, 4),
        "| overall RMSE =", round(overall_rmse, 4),
        "| Pearson =", round(overall_pearson, 4),
        "| Spearman =", round(overall_spearman, 4))

print(overall_summary)
print(metrics_all)
print(per_celltype)