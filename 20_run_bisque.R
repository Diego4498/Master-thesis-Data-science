# scripts/20_run_bisque.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(BisqueRNA)
  library(Biobase)
  library(data.table)
})

set.seed(1)

log_msg("Running Bisque deconvolution")

# -----------------------------
# Load aligned matrices
# -----------------------------

bulk <- readRDS(file.path(DIR_RESULTS, "aligned", "bulk_tpm_common_genes.rds"))

sce  <- readRDS(file.path(DIR_RESULTS, "scrna", "scrna_seurat_annot_singleR.rds"))

log_msg("Bulk dims (genes x samples):", nrow(bulk), "x", ncol(bulk))
log_msg("scRNA dims (genes x cells):", nrow(sce), "x", ncol(sce))

# -----------------------------
# Extract counts from Seurat
# -----------------------------

DefaultAssay(sce) <- "RNA"
sce <- JoinLayers(sce)

DefaultAssay(sce) <- "RNA"
sce <- JoinLayers(sce)

counts <- GetAssayData(sce, layer = "counts")

meta <- sce@meta.data

meta$cellcompartment <- as.character(meta$SingleR_cluster_label)

# Immune
meta$cellcompartment[meta$cellcompartment == "B_cell"] <- "B_cells"
meta$cellcompartment[meta$cellcompartment == "Macrophage"] <- "Myeloid"
meta$cellcompartment[meta$cellcompartment == "CMP"] <- "Myeloid"

# Major non-immune
meta$cellcompartment[meta$cellcompartment == "Epithelial_cells"] <- "Epithelial"
meta$cellcompartment[meta$cellcompartment == "Endothelial_cells"] <- "Endothelial"
meta$cellcompartment[meta$cellcompartment == "Smooth_muscle_cells"] <- "Stromal"

# Everything else
meta$cellcompartment[meta$cellcompartment %in% c("Tissue_stem_cells", "Neurons")] <- "Other"

meta$cellcompartment <- factor(
  meta$cellcompartment,
  levels = c(
    "T_cells",
    "B_cells",
    "Myeloid",
    "Epithelial",
    "Stromal",
    "Endothelial",
    "Other"
  )
)

sce$cellcompartment <- meta$cellcompartment

print(table(sce$cellcompartment))

stopifnot("cellcompartment" %in% colnames(sce@meta.data))
stopifnot("Sample" %in% colnames(sce@meta.data))

# -----------------------------
# Build ExpressionSet objects
# -----------------------------

log_msg("Building ExpressionSet objects")

bulk_eset <- ExpressionSet(
  assayData = as.matrix(bulk)
)

sc_eset <- ExpressionSet(
  assayData = as.matrix(counts),
  phenoData = AnnotatedDataFrame(
    data.frame(
      cellcompartment = meta$cellcompartment,
      Sample = meta$Sample,
      row.names = colnames(counts)
    )
  )
)

log_msg("bulk_eset dims:", nrow(exprs(bulk_eset)), "x", ncol(exprs(bulk_eset)))
log_msg("sc_eset dims:", nrow(exprs(sc_eset)), "x", ncol(exprs(sc_eset)))

# -----------------------------
# Run Bisque
# -----------------------------

log_msg("Running ReferenceBasedDecomposition")

res <- BisqueRNA::ReferenceBasedDecomposition(
  bulk.eset = bulk_eset,
  sc.eset = sc_eset,
  cell.types = "cellcompartment",
  subject.names = "Sample",
  use.overlap = FALSE,
  verbose = TRUE
)

# -----------------------------
# Extract proportions
# -----------------------------

props <- as.data.frame(t(res$bulk.props))

props$sample <- rownames(props)

# reorder columns (sample first)
props <- props[, c("sample", setdiff(colnames(props), "sample"))]

log_msg("Output rows (samples):", nrow(props))
log_msg("Output cols:", ncol(props))

# -----------------------------
# Save results
# -----------------------------

out_dir <- file.path(DIR_RESULTS, "deconv", "bisque")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(res, file.path(out_dir, "bisque_result.rds"))

data.table::fwrite(
  props,
  file.path(out_dir, "bisque_proportions.csv")
)

log_msg("Saved: bisque_result.rds")
log_msg("Saved: bisque_proportions.csv")

# -----------------------------
# Quick summary
# -----------------------------

log_msg("Summary of proportions")

print(summary(props[, setdiff(colnames(props), "sample")]))

bisque <- read.csv("results/deconv/bisque/bisque_proportions.csv",
                   row.names = 1,
                   check.names = FALSE)

bisque_mat <- as.matrix(bisque)

summary(t(bisque_mat))