# scripts/22_build_cibersort_signature.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

set.seed(1)

in_rds <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_compartments_singleR.rds")
stopifnot(file.exists(in_rds))

seu <- readRDS(in_rds)

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

meta <- seu@meta.data
stopifnot("cellcompartment" %in% colnames(meta))

counts <- GetAssayData(seu, layer = "counts")
logdata <- GetAssayData(seu, layer = "data")

celltypes <- levels(meta$cellcompartment)
cells_by_type <- split(colnames(seu), meta$cellcompartment)

log_msg("Cell types:", paste(celltypes, collapse = " | "))

# -----------------------------
# Build full signature by mean counts
# -----------------------------
sig_full <- sapply(celltypes, function(ct) {
  cells <- cells_by_type[[ct]]
  Matrix::rowMeans(counts[, cells, drop = FALSE])
})

sig_full <- as.matrix(sig_full)

# -----------------------------
# Select informative genes
# top markers per compartment based on log-expression difference
# -----------------------------
marker_list <- lapply(celltypes, function(ct) {
  cells_in  <- cells_by_type[[ct]]
  cells_out <- setdiff(colnames(seu), cells_in)
  
  m_in  <- Matrix::rowMeans(logdata[, cells_in, drop = FALSE])
  m_out <- Matrix::rowMeans(logdata[, cells_out, drop = FALSE])
  
  score <- m_in - m_out
  names(score) <- rownames(logdata)
  
  names(sort(score, decreasing = TRUE))[1:150]
})

marker_genes <- unique(unlist(marker_list))
marker_genes <- intersect(marker_genes, rownames(sig_full))

sig <- sig_full[marker_genes, , drop = FALSE]

# remove zero-variance genes across compartments
keep <- apply(sig, 1, function(x) sd(x) > 0)
sig <- sig[keep, , drop = FALSE]

log_msg("Reduced signature dims:", nrow(sig), "genes x", ncol(sig), "celltypes")

out_dir <- file.path(DIR_RESULTS, "deconv", "cibersort")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sig_file <- file.path(out_dir, "cibersort_signature_matrix.txt")

data.table::fwrite(
  data.table(Gene = rownames(sig), sig),
  sig_file,
  sep = "\t"
)

log_msg("Saved:", sig_file)