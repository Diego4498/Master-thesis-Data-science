# scripts/35_scrna_annotate_singleR_fine.R

source("scripts/00_setup_paths.R")

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "SingleR",
  "celldex",
  "SummarizedExperiment",
  "data.table"
)

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

set.seed(1)

in_proc <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_proc.rds")
stopifnot(file.exists(in_proc))

log_msg("Loading processed Seurat object:", in_proc)
seu <- readRDS(in_proc)

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

# Use log-normalized data for SingleR
logcounts <- GetAssayData(seu, layer = "data")
stopifnot(ncol(logcounts) == ncol(seu))

# -----------------------------
# Reference
# -----------------------------
# We reuse HPCA because it is already consistent with your previous pipeline.
# Here we try to leverage a finer label level if available.
ref <- celldex::HumanPrimaryCellAtlasData()

available_label_cols <- intersect(
  c("label.fine", "label.main"),
  colnames(SummarizedExperiment::colData(ref))
)

if (length(available_label_cols) == 0) {
  stop("No suitable label columns found in HPCA reference.")
}

label_col <- if ("label.fine" %in% available_label_cols) "label.fine" else "label.main"
labels_use <- SummarizedExperiment::colData(ref)[[label_col]]

log_msg("Running SingleR with label column:", label_col)

pred <- SingleR(
  test = logcounts,
  ref = ref,
  labels = labels_use,
  assay.type.test = "logcounts"
)

# Store fine labels in Seurat metadata
seu$SingleR_label_fine <- pred$labels
seu$SingleR_pruned_fine <- pred$pruned.labels

# Cluster-level majority fine label
cl <- as.character(Idents(seu))
lab_use <- ifelse(is.na(seu$SingleR_pruned_fine), seu$SingleR_label_fine, seu$SingleR_pruned_fine)

cluster_labels_fine <- tapply(lab_use, cl, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
})

seu$SingleR_cluster_label_fine <- cluster_labels_fine[cl]

# Confidence summaries by cluster
cluster_summary <- data.table(
  cluster = names(cluster_labels_fine),
  SingleR_cluster_label_fine = as.character(cluster_labels_fine),
  n_cells = as.integer(table(cl)[names(cluster_labels_fine)])
)

# add top-3 fine labels per cluster
cluster_top3 <- rbindlist(lapply(names(cluster_labels_fine), function(k) {
  idx <- which(cl == k)
  x <- lab_use[idx]
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(data.table(
      cluster = k,
      label_rank = integer(),
      label = character(),
      n = integer(),
      prop = numeric()
    ))
  }
  tab <- sort(table(x), decreasing = TRUE)
  tab <- head(tab, 3)
  data.table(
    cluster = k,
    label_rank = seq_along(tab),
    label = names(tab),
    n = as.integer(tab),
    prop = as.numeric(tab) / sum(as.integer(table(x)))
  )
}))

# per-cell export
cell_out <- data.table(
  cell = colnames(seu),
  cluster = cl,
  SingleR_label_fine = seu$SingleR_label_fine,
  SingleR_pruned_fine = seu$SingleR_pruned_fine,
  SingleR_cluster_label_fine = seu$SingleR_cluster_label_fine,
  Sample = if ("Sample" %in% colnames(seu@meta.data)) as.character(seu$Sample) else NA_character_
)

# -----------------------------
# Save outputs
# -----------------------------
out_dir <- file.path(DIR_RESULTS, "scrna", "level2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(cell_out, file.path(out_dir, "singleR_fine_cell_labels.csv"))
fwrite(cluster_summary, file.path(out_dir, "singleR_fine_cluster_labels.csv"))
fwrite(cluster_top3, file.path(out_dir, "singleR_fine_cluster_top3_labels.csv"))

# Save updated Seurat object
out_rds <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_annot_singleR_fine.rds")
saveRDS(seu, out_rds)

log_msg("Saved: singleR_fine_cell_labels.csv")
log_msg("Saved: singleR_fine_cluster_labels.csv")
log_msg("Saved: singleR_fine_cluster_top3_labels.csv")
log_msg("Saved:", out_rds)

# Quick summaries
tab_fine <- sort(table(seu$SingleR_cluster_label_fine), decreasing = TRUE)
log_msg(
  "Fine cluster-label summary:",
  paste(names(tab_fine), as.integer(tab_fine), sep = "=", collapse = " | ")
)

print(cluster_summary)
print(cluster_top3)