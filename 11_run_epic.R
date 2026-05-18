# scripts/11_run_epic.R

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

bulk_path <- file.path(DIR_RESULTS, "bulk", "bulk_tpm_matrix_M0M1.rds")
seu_path  <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_compartments_singleR.rds")

stopifnot(file.exists(bulk_path))
stopifnot(file.exists(seu_path))

bulk <- readRDS(bulk_path)  # genes x samples (TPM)
seu  <- readRDS(seu_path)

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

stopifnot("cellcompartment" %in% colnames(seu@meta.data))

log_msg("Building EPIC reference (CPM) from scRNA counts")

counts <- GetAssayData(seu, layer = "counts")
groups <- factor(seu$cellcompartment)

# CPM per cell
cpm <- t( t(counts) / Matrix::colSums(counts) ) * 1e6

# Average CPM per compartment
ref_epic <- sapply(levels(groups), function(g) {
  idx <- which(groups == g)
  Matrix::rowMeans(cpm[, idx, drop = FALSE])
})

# Gene-wise variance per compartment (for EPIC weighting)
ref_var <- sapply(levels(groups), function(g) {
  idx <- which(groups == g)
  if (length(idx) > 1) {
    apply(cpm[, idx, drop = FALSE], 1, var)
  } else {
    rep(0, nrow(cpm))
  }
})

rownames(ref_var) <- rownames(cpm)

common <- intersect(rownames(bulk), rownames(ref_epic))
stopifnot(length(common) > 0)

bulk2 <- bulk[common, , drop = FALSE]
ref2  <- ref_epic[common, , drop = FALSE]

log_msg("Common genes:", length(common))
log_msg("Bulk dims (genes x samples):", nrow(bulk2), "x", ncol(bulk2))
log_msg("Ref dims  (genes x compartments):", nrow(ref2), "x", ncol(ref2))

# EPIC reference format
ref_list <- list(
  refProfiles = ref2,
  refProfiles.var = ref_var[rownames(ref2), colnames(ref2), drop = FALSE],
  sigGenes = rownames(ref2)
)

log_msg("Running EPIC")
res <- EPIC::EPIC(bulk = bulk2, reference = ref_list)

out_dir <- file.path(DIR_RESULTS, "deconv", "epic")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(res, file.path(out_dir, "epic_result.rds"))

props <- as.data.frame(res$cellFractions)
props$sample <- rownames(props)
data.table::fwrite(as.data.table(props), file.path(out_dir, "epic_proportions.csv"))

if (!is.null(res$fit.gof)) {
  fg <- as.data.frame(res$fit.gof)
  fg$sample <- rownames(fg)
  data.table::fwrite(as.data.table(fg), file.path(out_dir, "epic_fit_gof.csv"))
}

log_msg("Saved: epic_result.rds")
log_msg("Saved: epic_proportions.csv")
log_msg("Saved: epic_fit_gof.csv")