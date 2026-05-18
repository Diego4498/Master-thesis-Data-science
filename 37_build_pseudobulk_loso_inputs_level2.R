# scripts/37_build_pseudobulk_loso_inputs_level2.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

set.seed(1)

in_rds <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_level2.rds")
stopifnot(file.exists(in_rds))

log_msg("Loading level2 Seurat object:", in_rds)
seu <- readRDS(in_rds)

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

meta <- seu@meta.data

required_cols <- c("Sample", "cellcompartment_lv2")
missing_cols <- required_cols[!required_cols %in% colnames(meta)]
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

counts <- GetAssayData(seu, layer = "counts")
stopifnot(ncol(counts) == nrow(meta))
stopifnot(identical(colnames(counts), rownames(meta)))

samples <- as.character(meta$Sample)

celltypes_lv2 <- c(
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

meta$cellcompartment_lv2 <- factor(as.character(meta$cellcompartment_lv2), levels = celltypes_lv2)

if (any(is.na(meta$cellcompartment_lv2))) {
  bad_n <- sum(is.na(meta$cellcompartment_lv2))
  stop("Found ", bad_n, " cells with missing cellcompartment_lv2 after factor conversion.")
}

sample_levels <- sort(unique(samples))
if (length(sample_levels) < 2) {
  stop("Need at least 2 samples for LOSO validation. Found: ", length(sample_levels))
}

log_msg("Samples found:", paste(sample_levels, collapse = " | "))
log_msg("Level2 compartments:", paste(celltypes_lv2, collapse = " | "))

# ----------------------------------------
# 1) True proportions per sample
# ----------------------------------------
truth_tab <- table(meta$Sample, meta$cellcompartment_lv2)
truth_prop <- prop.table(truth_tab, margin = 1)

truth_dt <- as.data.table(as.table(truth_prop))
setnames(truth_dt, c("Sample", "cellcompartment_lv2", "true_proportion"))

truth_wide <- dcast(
  truth_dt,
  Sample ~ cellcompartment_lv2,
  value.var = "true_proportion",
  fill = 0
)

truth_wide <- truth_wide[, c("Sample", celltypes_lv2), with = FALSE]

# ----------------------------------------
# 2) True cell counts per sample
# ----------------------------------------
truth_counts_dt <- as.data.table(as.table(truth_tab))
setnames(truth_counts_dt, c("Sample", "cellcompartment_lv2", "n_cells"))

truth_counts_wide <- dcast(
  truth_counts_dt,
  Sample ~ cellcompartment_lv2,
  value.var = "n_cells",
  fill = 0
)

truth_counts_wide <- truth_counts_wide[, c("Sample", celltypes_lv2), with = FALSE]

# ----------------------------------------
# 3) Pseudobulk counts by sample
# ----------------------------------------
log_msg("Building level2 pseudobulk count matrix by sample")

pb_mat <- sapply(sample_levels, function(s) {
  idx <- which(samples == s)
  Matrix::rowSums(counts[, idx, drop = FALSE])
})

pb_mat <- as.matrix(pb_mat)
colnames(pb_mat) <- sample_levels

stopifnot(nrow(pb_mat) == nrow(counts))
stopifnot(ncol(pb_mat) == length(sample_levels))
stopifnot(identical(rownames(pb_mat), rownames(counts)))

# ----------------------------------------
# 4) Pseudobulk CPM
# ----------------------------------------
log_msg("Computing level2 pseudobulk CPM matrix")

lib_sizes <- colSums(pb_mat)
if (any(lib_sizes == 0)) {
  stop("At least one pseudobulk sample has library size 0.")
}

pb_cpm <- t(t(pb_mat) / lib_sizes) * 1e6

# ----------------------------------------
# 5) LOSO fold table
# ----------------------------------------
folds_dt <- rbindlist(lapply(sample_levels, function(test_sample) {
  train_samples <- setdiff(sample_levels, test_sample)
  data.table(
    test_sample = test_sample,
    train_samples = paste(train_samples, collapse = ";"),
    n_train_samples = length(train_samples)
  )
}))

folds_long_dt <- rbindlist(lapply(sample_levels, function(test_sample) {
  train_samples <- setdiff(sample_levels, test_sample)
  data.table(
    test_sample = test_sample,
    train_sample = train_samples
  )
}))

# ----------------------------------------
# 6) Save outputs
# ----------------------------------------
out_dir <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso_level2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(pb_mat, file.path(out_dir, "pseudobulk_counts_by_sample_level2.rds"))
saveRDS(pb_cpm, file.path(out_dir, "pseudobulk_cpm_by_sample_level2.rds"))

fwrite(
  data.table(gene = rownames(pb_mat), pb_mat),
  file.path(out_dir, "pseudobulk_counts_by_sample_level2.csv")
)

fwrite(
  data.table(gene = rownames(pb_cpm), pb_cpm),
  file.path(out_dir, "pseudobulk_cpm_by_sample_level2.csv")
)

fwrite(truth_dt, file.path(out_dir, "ground_truth_proportions_long_level2.csv"))
fwrite(truth_wide, file.path(out_dir, "ground_truth_proportions_wide_level2.csv"))

fwrite(truth_counts_dt, file.path(out_dir, "ground_truth_cellcounts_long_level2.csv"))
fwrite(truth_counts_wide, file.path(out_dir, "ground_truth_cellcounts_wide_level2.csv"))

fwrite(folds_dt, file.path(out_dir, "loso_folds_level2.csv"))
fwrite(folds_long_dt, file.path(out_dir, "loso_folds_long_level2.csv"))

# ----------------------------------------
# 7) Console summary
# ----------------------------------------
log_msg("Saved: pseudobulk_counts_by_sample_level2.rds")
log_msg("Saved: pseudobulk_cpm_by_sample_level2.rds")
log_msg("Saved: ground_truth_proportions_long_level2.csv")
log_msg("Saved: ground_truth_proportions_wide_level2.csv")
log_msg("Saved: ground_truth_cellcounts_long_level2.csv")
log_msg("Saved: ground_truth_cellcounts_wide_level2.csv")
log_msg("Saved: loso_folds_level2.csv")
log_msg("Saved: loso_folds_long_level2.csv")

log_msg("Level2 pseudobulk counts dims (genes x samples):", nrow(pb_mat), "x", ncol(pb_mat))
log_msg("Truth table rows:", nrow(truth_dt))
log_msg("LOSO folds:", nrow(folds_dt))

print(truth_wide)
print(folds_dt)