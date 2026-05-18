# scripts/07_scrna_annotate_singleR.R

source("scripts/00_setup_paths.R")

required_pkgs <- c("Seurat", "SeuratObject", "SingleR", "celldex", "SummarizedExperiment", "data.table")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleR)
  library(celldex)
  library(SummarizedExperiment)
  library(data.table)
})

in_proc <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_proc.rds")
stopifnot(file.exists(in_proc))

seu <- readRDS(in_proc)
DefaultAssay(seu) <- "RNA"

# Seurat v5: ensure layers are joined so GetAssayData works for DE/exports
seu <- JoinLayers(seu)

# Use log-normalized data for SingleR
logcounts <- GetAssayData(seu, layer = "data")
stopifnot(ncol(logcounts) == ncol(seu))

# Reference: Human Primary Cell Atlas (broad cell types)
ref <- celldex::HumanPrimaryCellAtlasData()

log_msg("Running SingleR (cell-level)")
pred <- SingleR(
  test = logcounts,
  ref = ref,
  labels = ref$label.main,
  assay.type.test = "logcounts"
)

# Store in Seurat metadata
seu$SingleR_label <- pred$labels
seu$SingleR_pruned <- pred$pruned.labels

# Cluster-level majority label (using pruned if available)
cl <- as.character(Idents(seu))
lab_use <- ifelse(is.na(seu$SingleR_pruned), seu$SingleR_label, seu$SingleR_pruned)

cluster_labels <- tapply(lab_use, cl, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
})

seu$SingleR_cluster_label <- cluster_labels[cl]

out_dir <- file.path(DIR_RESULTS, "scrna")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Save predictions
cell_out <- data.table(
  cell = colnames(seu),
  cluster = cl,
  SingleR_label = seu$SingleR_label,
  SingleR_pruned = seu$SingleR_pruned,
  SingleR_cluster_label = seu$SingleR_cluster_label,
  Sample = if ("Sample" %in% colnames(seu@meta.data)) as.character(seu$Sample) else NA_character_
)

data.table::fwrite(cell_out, file.path(out_dir, "singleR_cell_labels.csv"))

cl_out <- data.table(
  cluster = names(cluster_labels),
  SingleR_cluster_label = as.character(cluster_labels),
  n_cells = as.integer(table(cl)[names(cluster_labels)])
)

data.table::fwrite(cl_out, file.path(out_dir, "singleR_cluster_labels.csv"))

# Save updated Seurat object
out_rds <- file.path(out_dir, "scrna_seurat_annot_singleR.rds")
saveRDS(seu, out_rds)

log_msg("Saved: singleR_cell_labels.csv")
log_msg("Saved: singleR_cluster_labels.csv")
log_msg("Saved:", out_rds)

# Quick summary
tab <- sort(table(seu$SingleR_cluster_label), decreasing = TRUE)
log_msg("Cluster-label summary:", paste(names(tab), as.integer(tab), sep = "=", collapse = " | "))
