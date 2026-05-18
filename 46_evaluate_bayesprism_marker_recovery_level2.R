# scripts/46_evaluate_bayesprism_marker_recovery_level2.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

set.seed(1)

val_dir  <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso_level2")
seu_path <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_level2.rds")
bp_dir   <- file.path(val_dir, "bayesprism", "folds")
folds_path <- file.path(val_dir, "loso_folds_level2.csv")

stopifnot(file.exists(seu_path))
stopifnot(dir.exists(bp_dir))
stopifnot(file.exists(folds_path))

seu <- readRDS(seu_path)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

meta <- seu@meta.data
counts <- GetAssayData(seu, layer = "counts")

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

# marker panel simple y defendible
marker_panel <- list(
  CD4_T = c("IL7R", "LTB", "MALAT1"),
  NK_cells = c("NKG7", "GNLY", "CCL5"),
  B_cells = c("MS4A1", "CD79A", "CD74"),
  Plasma_cells = c("MZB1", "JCHAIN", "SDC1"),
  Macrophage = c("C1QA", "C1QB", "APOE"),
  DC = c("FCER1A", "CD1C", "CLEC10A"),
  Epithelial = c("EPCAM", "KRT8", "KRT18"),
  Stromal = c("COL1A1", "COL1A2", "DCN"),
  Endothelial = c("PECAM1", "VWF", "KDR")
)

# ----------------------------------------
# helper: build true held-out pseudobulk expression by cell type
# rows = genes, cols = celltypes
# ----------------------------------------
build_truth_expr_for_sample <- function(seu, counts, sample_id, celltypes) {
  md <- seu@meta.data
  keep_cells <- rownames(md)[as.character(md$Sample) == sample_id]
  
  if (length(keep_cells) == 0) {
    stop("No cells found for sample: ", sample_id)
  }
  
  ct <- as.character(md[keep_cells, "cellcompartment_lv2", drop = TRUE])
  
  out <- lapply(celltypes, function(lbl) {
    cells_lbl <- keep_cells[ct == lbl]
    if (length(cells_lbl) == 0) {
      return(rep(0, nrow(counts)))
    } else {
      return(Matrix::rowSums(counts[, cells_lbl, drop = FALSE]))
    }
  })
  
  out <- do.call(cbind, out)
  rownames(out) <- rownames(counts)
  colnames(out) <- celltypes
  out
}

# ----------------------------------------
# helper: extract BayesPrism deconvolved expression
# returns matrix with rows = genes, cols = celltypes
# ----------------------------------------
extract_bayesprism_expr <- function(bp_res, celltypes) {
  if (!"reference.update" %in% slotNames(bp_res)) {
    stop("BayesPrism result does not contain 'reference.update' slot.")
  }
  
  ref_upd <- bp_res@reference.update
  
  if (!all(c("psi_mal", "psi_env", "key") %in% slotNames(ref_upd))) {
    stop("reference.update does not contain expected slots psi_mal / psi_env / key.")
  }
  
  psi_mal <- ref_upd@psi_mal   # matrix: rows = malignant sample(s), cols = genes
  psi_env <- ref_upd@psi_env   # matrix: rows = env celltypes, cols = genes
  key_lbl <- ref_upd@key       # should be "Epithelial"
  
  psi_mal <- as.matrix(psi_mal)
  psi_env <- as.matrix(psi_env)
  
  # Convert to genes x celltypes
  env_mat <- t(psi_env)
  mal_vec <- as.numeric(psi_mal[1, ])
  names(mal_vec) <- colnames(psi_mal)
  mal_mat <- matrix(mal_vec, ncol = 1)
  rownames(mal_mat) <- names(mal_vec)
  colnames(mal_mat) <- key_lbl
  
  expr_mat <- cbind(env_mat, mal_mat)
  
  # Ensure all expected celltypes exist
  for (ct in celltypes) {
    if (!ct %in% colnames(expr_mat)) {
      expr_mat <- cbind(expr_mat, setNames(data.frame(rep(0, nrow(expr_mat))), ct))
      expr_mat <- as.matrix(expr_mat)
      rownames(expr_mat) <- rownames(env_mat)
    }
  }
  
  expr_mat <- expr_mat[, celltypes, drop = FALSE]
  expr_mat
}

folds <- fread(folds_path)

marker_results <- list()

for (i in seq_len(nrow(folds))) {
  test_sample <- folds$test_sample[i]
  
  message("Evaluating marker recovery for ", test_sample)
  
  bp_rds <- file.path(bp_dir, paste0("bayesprism_level2_result_", test_sample, ".rds"))
  if (!file.exists(bp_rds)) {
    stop("Missing BayesPrism result file: ", bp_rds)
  }
  
  bp_res <- readRDS(bp_rds)
  
  # true held-out expression
  truth_expr <- build_truth_expr_for_sample(
    seu = seu,
    counts = counts,
    sample_id = test_sample,
    celltypes = celltypes
  )
  
  # BayesPrism deconvolved expression
  pred_expr <- extract_bayesprism_expr(bp_res, celltypes)
  
  common_genes <- intersect(rownames(truth_expr), rownames(pred_expr))
  common_ct <- intersect(colnames(truth_expr), colnames(pred_expr))
  
  truth_expr <- truth_expr[common_genes, common_ct, drop = FALSE]
  pred_expr  <- pred_expr[common_genes, common_ct, drop = FALSE]
  
  for (ct in names(marker_panel)) {
    if (!ct %in% common_ct) next
    
    genes <- intersect(marker_panel[[ct]], common_genes)
    if (length(genes) == 0) next
    
    truth_vals <- truth_expr[genes, ct, drop = TRUE]
    pred_vals  <- pred_expr[genes, ct, drop = TRUE]
    
    res_dt <- data.table(
      sample = test_sample,
      celltype = ct,
      gene = genes,
      truth_expr = as.numeric(truth_vals),
      pred_expr = as.numeric(pred_vals)
    )
    
    marker_results[[paste(test_sample, ct, sep = "__")]] <- res_dt
  }
}

marker_dt <- rbindlist(marker_results, use.names = TRUE, fill = TRUE)

if (nrow(marker_dt) == 0) {
  stop("No marker results produced. Check marker genes and BayesPrism expression extraction.")
}

marker_dt[, abs_error := abs(pred_expr - truth_expr)]

summary_dt <- marker_dt[, .(
  mean_abs_error = mean(abs_error, na.rm = TRUE),
  pearson = suppressWarnings(cor(truth_expr, pred_expr, method = "pearson", use = "complete.obs")),
  spearman = suppressWarnings(cor(truth_expr, pred_expr, method = "spearman", use = "complete.obs"))
), by = .(celltype)]

out_dir <- file.path(val_dir, "bayesprism_marker_recovery")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(marker_dt, file.path(out_dir, "bayesprism_marker_recovery_long.csv"))
fwrite(summary_dt, file.path(out_dir, "bayesprism_marker_recovery_summary.csv"))

p <- ggplot(marker_dt, aes(x = truth_expr, y = pred_expr)) +
  geom_point(size = 2) +
  facet_wrap(~ celltype, scales = "free", ncol = 3) +
  theme_bw() +
  labs(
    title = "BayesPrism marker recovery by cell type",
    x = "True held-out expression",
    y = "BayesPrism deconvolved expression"
  )

ggsave(
  filename = file.path(out_dir, "bayesprism_marker_recovery_scatter.png"),
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)

print(summary_dt)