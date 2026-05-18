# scripts/47_run_bayesprism_gsea_level2.R

source("scripts/00_setup_paths.R")

required_pkgs <- c("data.table", "fgsea", "msigdbr", "ggplot2")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

set.seed(1)

val_dir <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso_level2")
bp_dir  <- file.path(val_dir, "bayesprism", "folds")
folds_path <- file.path(val_dir, "loso_folds_level2.csv")

stopifnot(dir.exists(bp_dir))
stopifnot(file.exists(folds_path))

folds <- fread(folds_path)

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

# ----------------------------------------
# Hallmark gene sets
# ----------------------------------------
msig_h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")

pathways <- split(msig_h$gene_symbol, msig_h$gs_name)

# Optional: focus on a smaller subset first
hallmark_keep <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

pathways <- pathways[names(pathways) %in% hallmark_keep]

# ----------------------------------------
# Extract BayesPrism deconvolved expression
# returns genes x celltypes
# ----------------------------------------
extract_bayesprism_expr <- function(bp_res, celltypes) {
  if (!"reference.update" %in% slotNames(bp_res)) {
    stop("BayesPrism result does not contain 'reference.update' slot.")
  }
  
  ref_upd <- bp_res@reference.update
  psi_mal <- as.matrix(ref_upd@psi_mal)  # rows = malignant sample(s), cols = genes
  psi_env <- as.matrix(ref_upd@psi_env)  # rows = env celltypes, cols = genes
  key_lbl <- ref_upd@key
  
  env_mat <- t(psi_env) # genes x env_celltypes
  
  mal_vec <- as.numeric(psi_mal[1, ])
  names(mal_vec) <- colnames(psi_mal)
  mal_mat <- matrix(mal_vec, ncol = 1)
  rownames(mal_mat) <- names(mal_vec)
  colnames(mal_mat) <- key_lbl
  
  expr_mat <- cbind(env_mat, mal_mat)
  
  for (ct in celltypes) {
    if (!ct %in% colnames(expr_mat)) {
      zero_col <- matrix(0, nrow = nrow(expr_mat), ncol = 1)
      rownames(zero_col) <- rownames(expr_mat)
      colnames(zero_col) <- ct
      expr_mat <- cbind(expr_mat, zero_col)
    }
  }
  
  expr_mat <- expr_mat[, celltypes, drop = FALSE]
  expr_mat
}

# ----------------------------------------
# Build gene ranking for one cell type
# target cell type vs mean of all others
# ----------------------------------------
build_rank_stat <- function(expr_mat, target_ct) {
  stopifnot(target_ct %in% colnames(expr_mat))
  
  others <- setdiff(colnames(expr_mat), target_ct)
  if (length(others) == 0) {
    stop("No other cell types available for comparison.")
  }
  
  target <- expr_mat[, target_ct]
  background <- rowMeans(expr_mat[, others, drop = FALSE])
  
  stat <- target - background
  stat <- sort(stat, decreasing = TRUE)
  stat <- stat[!is.na(stat)]
  stat
}

gsea_results <- list()

for (i in seq_len(nrow(folds))) {
  test_sample <- folds$test_sample[i]
  message("Running BayesPrism GSEA for fold ", test_sample)
  
  bp_rds <- file.path(bp_dir, paste0("bayesprism_level2_result_", test_sample, ".rds"))
  if (!file.exists(bp_rds)) {
    stop("Missing BayesPrism result file: ", bp_rds)
  }
  
  bp_res <- readRDS(bp_rds)
  expr_mat <- extract_bayesprism_expr(bp_res, celltypes)
  
  # restrict to genes present in at least one pathway
  pathway_genes <- unique(unlist(pathways))
  common_genes <- intersect(rownames(expr_mat), pathway_genes)
  expr_mat <- expr_mat[common_genes, , drop = FALSE]
  
  for (ct in setdiff(celltypes, "Other")) {
    stats <- build_rank_stat(expr_mat, ct)
    
    fg <- fgsea::fgsea(
      pathways = pathways,
      stats = stats,
      minSize = 5,
      maxSize = 500,
      nproc = 1
    )
    
    fg <- as.data.table(fg)
    fg[, sample := test_sample]
    fg[, celltype := ct]
    
    gsea_results[[paste(test_sample, ct, sep = "__")]] <- fg
  }
}

gsea_dt <- rbindlist(gsea_results, use.names = TRUE, fill = TRUE)

setcolorder(gsea_dt, c(
  "sample", "celltype", "pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge"
))

out_dir <- file.path(val_dir, "bayesprism_gsea")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(gsea_dt, file.path(out_dir, "bayesprism_level2_gsea_hallmark_long.csv"))

# ----------------------------------------
# Summaries
# ----------------------------------------
sig_dt <- gsea_dt[padj < 0.05]

summary_dt <- sig_dt[, .(
  n_sig_pathways = .N,
  top_pathway = pathway[which.max(abs(NES))],
  top_NES = NES[which.max(abs(NES))]
), by = .(celltype, sample)]

fwrite(summary_dt, file.path(out_dir, "bayesprism_level2_gsea_hallmark_summary.csv"))

mean_nes_dt <- gsea_dt[, .(
  mean_NES = mean(NES, na.rm = TRUE)
), by = .(celltype, pathway)]

fwrite(mean_nes_dt, file.path(out_dir, "bayesprism_level2_gsea_hallmark_meanNES.csv"))

heat_dt <- dcast(mean_nes_dt, pathway ~ celltype, value.var = "mean_NES", fill = 0)
fwrite(heat_dt, file.path(out_dir, "bayesprism_level2_gsea_hallmark_meanNES_wide.csv"))

if (nrow(sig_dt) > 0) {
  top_plot_dt <- gsea_dt[
    pathway %in% unique(sig_dt$pathway)
  ]
  
  p <- ggplot(top_plot_dt, aes(x = celltype, y = NES, fill = sample)) +
    geom_col(position = "dodge") +
    facet_wrap(~ pathway, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "BayesPrism level 2 GSEA (Hallmark)",
      x = NULL,
      y = "NES"
    )
  
  ggsave(
    filename = file.path(out_dir, "bayesprism_level2_gsea_hallmark_sig_pathways.png"),
    plot = p,
    width = 14,
    height = 9,
    dpi = 300
  )
}

print(summary_dt)