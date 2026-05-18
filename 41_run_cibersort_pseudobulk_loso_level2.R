# scripts/41_run_cibersort_pseudobulk_loso_level2.R

source("scripts/00_setup_paths.R")

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "data.table", "CIBERSORT")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(CIBERSORT)
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

out_dir <- file.path(val_dir, "cibersort")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(folds))) {
  test_sample <- folds$test_sample[i]
  train_samples <- unlist(strsplit(folds$train_samples[i], ";", fixed = TRUE))
  
  log_msg("Running CIBERSORT level2 LOSO fold | test =", test_sample,
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
  # Build CIBERSORT signature from training scRNA
  # -----------------------------
  counts <- GetAssayData(seu_train, layer = "counts")
  logdata <- GetAssayData(seu_train, layer = "data")
  
  counts <- as.matrix(counts)
  logdata <- as.matrix(logdata)
  
  groups <- as.character(seu_train$cellcompartment_lv2)
  cells_by_type <- split(colnames(seu_train), groups)
  
  present_celltypes <- intersect(celltypes, names(cells_by_type))
  if (length(present_celltypes) < 2) {
    stop("Too few level2 cell types present in training fold for CIBERSORT: ", test_sample)
  }
  
  sig_full <- sapply(present_celltypes, function(ct) {
    cells <- cells_by_type[[ct]]
    Matrix::rowMeans(counts[, cells, drop = FALSE])
  })
  sig_full <- as.matrix(sig_full)
  
  marker_list <- lapply(present_celltypes, function(ct) {
    cells_in <- cells_by_type[[ct]]
    cells_out <- setdiff(colnames(seu_train), cells_in)
    
    m_in <- Matrix::rowMeans(logdata[, cells_in, drop = FALSE])
    
    if (length(cells_out) > 0) {
      m_out <- Matrix::rowMeans(logdata[, cells_out, drop = FALSE])
    } else {
      m_out <- rep(0, nrow(logdata))
      names(m_out) <- rownames(logdata)
    }
    
    score <- m_in - m_out
    names(score) <- rownames(logdata)
    
    head(names(sort(score, decreasing = TRUE)), 150)
  })
  
  marker_genes <- unique(unlist(marker_list))
  marker_genes <- intersect(marker_genes, rownames(sig_full))
  
  sig <- sig_full[marker_genes, , drop = FALSE]
  
  keep <- apply(sig, 1, function(x) sd(x) > 0)
  sig <- sig[keep, , drop = FALSE]
  
  if (nrow(sig) == 0) {
    stop("No signature genes left after filtering for fold: ", test_sample)
  }
  
  log_msg("Signature dims:", nrow(sig), "genes x", ncol(sig), "celltypes")
  
  # -----------------------------
  # Use ALL pseudobulk samples together as mixture input
  # -----------------------------
  bulk_all <- as.matrix(pb_counts)
  
  common <- intersect(rownames(bulk_all), rownames(sig))
  if (length(common) == 0) {
    stop("No overlapping genes between pseudobulk mixture and CIBERSORT signature for fold: ", test_sample)
  }
  
  bulk2 <- bulk_all[common, , drop = FALSE]
  sig2  <- sig[common, , drop = FALSE]
  
  bulk2 <- as.matrix(bulk2)
  sig2  <- as.matrix(sig2)
  
  log_msg("Common genes:", length(common))
  log_msg("Mixture dims:", nrow(bulk2), "genes x", ncol(bulk2), "samples")
  log_msg("Signature dims after alignment:", nrow(sig2), "genes x", ncol(sig2), "celltypes")
  
  # -----------------------------
  # Write temporary files
  # -----------------------------
  fold_dir <- file.path(out_dir, "folds")
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  
  sig_file <- file.path(fold_dir, paste0("cibersort_level2_signature_", test_sample, ".txt"))
  mix_file <- file.path(fold_dir, paste0("cibersort_level2_mixture_", test_sample, ".txt"))
  
  data.table::fwrite(
    data.table(Gene = rownames(sig2), sig2),
    sig_file,
    sep = "\t"
  )
  
  data.table::fwrite(
    data.table(Gene = rownames(bulk2), bulk2),
    mix_file,
    sep = "\t"
  )
  
  # -----------------------------
  # Run CIBERSORT
  # -----------------------------
  res <- cibersort(
    sig_file,
    mix_file,
    perm = 0,
    QN = FALSE
  )
  
  # -----------------------------
  # Standardize output robustly
  # -----------------------------
  res_df <- as.data.frame(res, check.names = FALSE)
  
  if (!is.null(rownames(res_df)) &&
      length(rownames(res_df)) == nrow(res_df) &&
      any(rownames(res_df) %in% colnames(pb_counts))) {
    res_df$sample <- rownames(res_df)
  } else {
    first_col <- colnames(res_df)[1]
    res_df$sample <- as.character(res_df[[first_col]])
    res_df[[first_col]] <- NULL
  }
  
  res_dt <- as.data.table(res_df)
  setcolorder(res_dt, c("sample", setdiff(colnames(res_dt), "sample")))
  res_dt$sample <- as.character(res_dt$sample)
  
  log_msg("CIBERSORT returned samples:", paste(res_dt$sample, collapse = " | "))
  
  for (ct in celltypes) {
    if (!ct %in% colnames(res_dt)) res_dt[[ct]] <- 0
  }
  
  extra_cols <- intersect(c("P.value", "Correlation", "RMSE"), colnames(res_dt))
  res_keep <- res_dt[, c("sample", celltypes, extra_cols), with = FALSE]
  
  # Keep only held-out sample for evaluation
  res_test <- res_keep[res_keep$sample == test_sample]
  if (nrow(res_test) != 1) {
    stop("Expected exactly one CIBERSORT prediction row for held-out sample: ", test_sample,
         " | found: ", nrow(res_test),
         " | available samples: ", paste(res_keep$sample, collapse = ", "))
  }
  
  # -----------------------------
  # Truth for held-out sample
  # -----------------------------
  truth_row <- truth[truth$Sample == test_sample]
  if (nrow(truth_row) != 1) {
    stop("Expected exactly one truth row for sample: ", test_sample)
  }
  
  pred_long <- data.table::melt(
    res_test[, c("sample", celltypes), with = FALSE],
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
  fold_pred[, method := "CIBERSORT"]
  fold_pred[, fold_test_sample := test_sample]
  fold_pred[, abs_error := abs(estimated_proportion - true_proportion)]
  fold_pred[, sq_error := (estimated_proportion - true_proportion)^2]
  
  if ("P.value" %in% colnames(res_test)) fold_pred[, P.value := res_test$P.value[1]]
  if ("Correlation" %in% colnames(res_test)) fold_pred[, Correlation := res_test$Correlation[1]]
  if ("RMSE" %in% colnames(res_test)) fold_pred[, CIBERSORT_RMSE := res_test$RMSE[1]]
  
  pred_list[[test_sample]] <- fold_pred
  
  fold_mae <- mean(fold_pred$abs_error, na.rm = TRUE)
  fold_rmse <- sqrt(mean(fold_pred$sq_error, na.rm = TRUE))
  
  metrics_row <- data.table(
    method = "CIBERSORT",
    test_sample = test_sample,
    n_train_samples = length(train_samples),
    n_train_cells = ncol(seu_train),
    n_mixture_samples_used = ncol(bulk2),
    n_signature_genes = nrow(sig2),
    fold_mae = fold_mae,
    fold_rmse = fold_rmse
  )
  
  if ("P.value" %in% colnames(res_test)) metrics_row[, P.value := res_test$P.value[1]]
  if ("Correlation" %in% colnames(res_test)) metrics_row[, Correlation := res_test$Correlation[1]]
  if ("RMSE" %in% colnames(res_test)) metrics_row[, CIBERSORT_RMSE := res_test$RMSE[1]]
  
  metrics_list[[test_sample]] <- metrics_row
  
  fwrite(res_keep, file.path(fold_dir, paste0("cibersort_level2_allmixture_result_", test_sample, ".csv")))
  
  log_msg("Completed CIBERSORT level2 fold | test =", test_sample,
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
  method = "CIBERSORT",
  overall_mae = overall_mae,
  overall_rmse = overall_rmse,
  overall_pearson = overall_pearson,
  overall_spearman = overall_spearman
)

# ----------------------------------------
# Save outputs
# ----------------------------------------
data.table::fwrite(pred_all, file.path(out_dir, "cibersort_level2_loso_predictions_long.csv"))
data.table::fwrite(pred_wide, file.path(out_dir, "cibersort_level2_loso_predictions_wide.csv"))
data.table::fwrite(metrics_all, file.path(out_dir, "cibersort_level2_loso_fold_metrics.csv"))
data.table::fwrite(per_celltype, file.path(out_dir, "cibersort_level2_loso_metrics_by_celltype.csv"))
data.table::fwrite(overall_summary, file.path(out_dir, "cibersort_level2_loso_overall_metrics.csv"))

log_msg("Saved:", file.path(out_dir, "cibersort_level2_loso_predictions_long.csv"))
log_msg("Saved:", file.path(out_dir, "cibersort_level2_loso_predictions_wide.csv"))
log_msg("Saved:", file.path(out_dir, "cibersort_level2_loso_fold_metrics.csv"))
log_msg("Saved:", file.path(out_dir, "cibersort_level2_loso_metrics_by_celltype.csv"))
log_msg("Saved:", file.path(out_dir, "cibersort_level2_loso_overall_metrics.csv"))

log_msg("CIBERSORT level2 LOSO overall MAE =", round(overall_mae, 4),
        "| overall RMSE =", round(overall_rmse, 4),
        "| Pearson =", round(overall_pearson, 4),
        "| Spearman =", round(overall_spearman, 4))

print(overall_summary)
print(metrics_all)
print(per_celltype)