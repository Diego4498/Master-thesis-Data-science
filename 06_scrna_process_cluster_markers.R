# scripts/06_scrna_process_cluster_markers.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

in_qc <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_qc.rds")
stopifnot(file.exists(in_qc))

seu <- readRDS(in_qc)

log_msg("Processing scRNA QC object")
log_msg("Input dims (genes x cells):", nrow(seu), "x", ncol(seu))

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)

out_dir <- file.path(DIR_RESULTS, "scrna")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_proc <- file.path(out_dir, "scrna_seurat_proc.rds")
saveRDS(seu, out_proc)
log_msg("Saved:", out_proc)

cl_tab <- table(Idents(seu))
log_msg("Clusters:", paste(names(cl_tab), as.integer(cl_tab), sep = "=", collapse = " | "))

seu <- JoinLayers(seu)
log_msg("Joined assay layers")

log_msg("Computing cluster markers (positive only)")
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

m_path <- file.path(out_dir, "cluster_markers_res0.8.csv")
data.table::fwrite(as.data.table(markers), m_path)
log_msg("Saved:", m_path)