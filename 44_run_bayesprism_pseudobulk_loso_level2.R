# scripts/44_run_bayesprism_pseudobulk_loso_level2.R

source("scripts/00_setup_paths.R")

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "data.table")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))

if (!requireNamespace("BayesPrism", quietly = TRUE)) {
  stop("Package 'BayesPrism' is not installed.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(BayesPrism)
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

truth <- truth[match(colnames(pb_counts), truth$Sample)]
stopifnot(identical(colnames(pb_counts), truth$Sample))

meta <- seu@meta.data
meta$Sample <- as.character(meta$Sample)
meta$cellcompartment_lv2 <- factor(as.character(meta$cellcompartment_lv2), levels = celltypes)

pred_list <- list()
metrics_list <- list()

out_dir <- file.path(val_dir, "bayesprism")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(folds))) {
  test_sample <- folds$test_sample[i]
  train_samples <- unlist(strsplit(folds$train_samples[i], ";", fixed = TRUE))
  
  log_msg("Running BayesPrism level2 LOSO fold | test =", test_sample,
          "| train =", paste(train_samples, collapse = " | "))
  
  # -----------------------------
  # Training scRNA only
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
  # Counts
  # -----------------------------
  sc_counts <- GetAssayData(seu_train, layer = "counts")
  sc_counts <- as.matrix(sc_counts)
  
  bulk_test <- pb_counts[, test_sample, drop = FALSE]
  bulk_test <- as.matrix(bulk_test)
  
  common <- intersect(rownames(sc_counts), rownames(bulk_test))
  if (length(common) == 0) {
    stop("No overlapping genes for fold: ", test_sample)
  }
  
  sc_counts <- sc_counts[common, , drop = FALSE]
  bulk_test <- bulk_test[common, , drop = FALSE]
  
  # BayesPrism expects cells x genes for reference and samples x genes for mixture
  reference_mat <- t(sc_counts)
  mixture_mat <- t(bulk_test)
  
  cell_type_labels <- as.character(seu_train$cellcompartment_lv2)
  cell_state_labels <- as.character(seu_train$cellcompartment_lv2)
  
  if (!"Epithelial" %in% unique(cell_type_labels)) {
    stop("BayesPrism key 'Epithelial' is not present in cell_type_labels for fold: ", test_sample)
  }
  
  if (!"Epithelial" %in% unique(cell_state_labels)) {
    stop("BayesPrism key 'Epithelial' is not present in cell_state_labels for fold: ", test_sample)
  }
  
  log_msg("Common genes:", length(common))
  log_msg("Reference dims (cells x genes):", nrow(reference_mat), "x", ncol(reference_mat))
  log_msg("Mixture dims (samples x genes):", nrow(mixture_mat), "x", ncol(mixture_mat))
  
  # -----------------------------
  # Build BayesPrism object
  # -----------------------------
  prism_obj <- BayesPrism::new.prism(
    reference = reference_mat,
    mixture = mixture_mat,
    input.type = "count.matrix",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_state_labels,
    key = "Epithelial",
    outlier.cut = 0.01,
    outlier.fraction = 0.1
  )
  
  # -----------------------------
  # Run BayesPrism
  # -----------------------------
  bp_res <- BayesPrism::run.prism(prism_obj, n.cores = 1)
  
  # -----------------------------
  # Extract proportions
  # -----------------------------
  theta <- NULL
  
  if ("posterior.theta_f" %in% slotNames(bp_res)) {
    theta <- bp_res@posterior.theta_f
  } else if ("posterior.initial.cellTypeInfo" %in% slotNames(bp_res)) {
    if (!is.null(bp_res@posterior.initial.cellTypeInfo$theta_f)) {
      theta <- bp_res@posterior.initial.cellTypeInfo$theta_f
    }
  }
  
  if (is.null(theta)) {
    stop("Could not extract BayesPrism theta object for fold: ", test_sample)
  }
  
  if (inherits(theta, "thetaPost")) {
    if ("theta" %in% slotNames(theta)) {
      est_mat <- theta@theta
    } else if ("theta_f" %in% slotNames(theta)) {
      est_mat <- theta@theta_f
    } else {
      stop("thetaPost object found, but no usable theta/theta_f slot for fold: ", test_sample)
    }
  } else if (is.matrix(theta) || is.data.frame(theta)) {
    est_mat <- theta
  } else {
    stop("Unsupported theta object class for fold ", test_sample, ": ", paste(class(theta), collapse = ", "))
  }
  
  est_mat <- as.matrix(est_mat)
  est <- as.data.frame(est_mat)
  est$sample <- rownames(est)
  est$sample <- as.character(est$sample)
  
  for (ct in celltypes) {
    if (!ct %in% colnames(est)) est[[ct]] <- 0
  }
  
  est <- est[, c("sample", celltypes), drop = FALSE]
  est <- est[est$sample == test_sample, , drop = FALSE]
  
  truth_row <- truth[truth$Sample == test_sample]
  if (nrow(truth_row) != 1) {
    stop("Expected exactly one truth row for sample: ", test_sample)
  }
  
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
  fold_pred[, method := "BayesPrism"]
  fold_pred[, fold_test_sample := test_sample]
  fold_pred[, abs_error := abs(estimated_proportion - true_proportion)]
  fold_pred[, sq_error := (estimated_proportion - true_proportion)^2]
  
  pred_list[[test_sample]] <- fold_pred
  
  fold_mae <- mean(fold_pred$abs_error, na.rm = TRUE)
  fold_rmse <- sqrt(mean(fold_pred$sq_error, na.rm = TRUE))
  
  metrics_list[[test_sample]] <- data.table(
    method = "BayesPrism",
    test_sample = test_sample,
    n_train_samples = length(train_samples),
    n_train_cells = ncol(seu_train),
    n_common_genes = length(common),
    fold_mae = fold_mae,
    fold_rmse = fold_rmse
  )
  
  fold_dir <- file.path(out_dir, "folds")
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(bp_res, file.path(fold_dir, paste0("bayesprism_level2_result_", test_sample, ".rds")))
  
  log_msg("Completed BayesPrism level2 fold | test =", test_sample,
          "| MAE =", round(fold_mae, 4),
          "| RMSE =", round(fold_rmse, 4))
}

pred_all <- rbindlist(pred_list, use.names = TRUE, fill = TRUE)
metrics_all <- rbindlist(metrics_list, use.names = TRUE, fill = TRUE)

pred_wide <- dcast(
  pred_all[, .(sample, celltype, estimated_proportion)],
  sample ~ celltype,
  value.var = "estimated_proportion"
)

pred_wide <- pred_wide[match(truth$Sample, pred_wide$sample)]

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
  method = "BayesPrism",
  overall_mae = overall_mae,
  overall_rmse = overall_rmse,
  overall_pearson = overall_pearson,
  overall_spearman = overall_spearman
)

fwrite(pred_all, file.path(out_dir, "bayesprism_level2_loso_predictions_long.csv"))
fwrite(pred_wide, file.path(out_dir, "bayesprism_level2_loso_predictions_wide.csv"))
fwrite(metrics_all, file.path(out_dir, "bayesprism_level2_loso_fold_metrics.csv"))
fwrite(per_celltype, file.path(out_dir, "bayesprism_level2_loso_metrics_by_celltype.csv"))
fwrite(overall_summary, file.path(out_dir, "bayesprism_level2_loso_overall_metrics.csv"))

log_msg("Saved BayesPrism level2 outputs in:", out_dir)
print(overall_summary)
print(metrics_all)
print(per_celltype)