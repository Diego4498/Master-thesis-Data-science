# scripts/08_build_reference_from_singleR.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

in_rds <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_annot_singleR.rds")
stopifnot(file.exists(in_rds))

seu <- readRDS(in_rds)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

# Use cluster-level SingleR label (stable)
lab <- as.character(seu$SingleR_cluster_label)
stopifnot(length(lab) == ncol(seu))

# Collapse SingleR labels to broad compartments for joint modeling
lab <- as.character(seu$SingleR_cluster_label)
stopifnot(length(lab) == ncol(seu))

comp <- rep("Other", length(lab))

# Immune
comp[lab %in% c("T_cells")] <- "T_cells"
comp[lab %in% c("B_cell")] <- "B_cells"
comp[lab %in% c("Macrophage", "CMP")] <- "Myeloid"

# Major non-immune
comp[lab %in% c("Epithelial_cells")] <- "Epithelial"
comp[lab %in% c("Endothelial_cells")] <- "Endothelial"
comp[lab %in% c("Smooth_muscle_cells")] <- "Stromal"

# Everything else stays "Other" (e.g., Tissue_stem_cells, Neurons, etc.)
seu$cellcompartment <- factor(
  comp,
  levels = c("T_cells", "B_cells", "Myeloid", "Epithelial", "Stromal", "Endothelial", "Other")
)

tab <- sort(table(seu$cellcompartment), decreasing = TRUE)
log_msg("Cells per compartment:", paste(names(tab), as.integer(tab), sep="=", collapse=" | "))

# Build reference: average log-normalized expression per compartment
data_mat <- GetAssayData(seu, layer = "data")
stopifnot(ncol(data_mat) == ncol(seu))

groups <- seu$cellcompartment
ref <- sapply(levels(groups), function(g) {
  cells_g <- which(groups == g)
  Matrix::rowMeans(data_mat[, cells_g, drop = FALSE])
})

out_dir <- file.path(DIR_RESULTS, "references")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(ref, file.path(out_dir, "reference_logexpr_compartments_singleR.rds"))
data.table::fwrite(
  data.table(gene = rownames(ref), ref),
  file.path(out_dir, "reference_logexpr_compartments_singleR.csv")
)

# Save Seurat with compartments
saveRDS(seu, file.path(DIR_RESULTS, "scrna", "scrna_seurat_compartments_singleR.rds"))

log_msg("Saved: reference_logexpr_compartments_singleR.rds")
log_msg("Saved: reference_logexpr_compartments_singleR.csv")
log_msg("Saved: scrna_seurat_compartments_singleR.rds")
log_msg("Reference dims (genes x compartments):", nrow(ref), "x", ncol(ref))