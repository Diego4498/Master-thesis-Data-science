# scripts/36_build_level2_from_fine_labels.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
})

set.seed(1)

in_rds <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_annot_singleR_fine.rds")
stopifnot(file.exists(in_rds))

log_msg("Loading fine-annotated Seurat object:", in_rds)
seu <- readRDS(in_rds)

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

meta <- seu@meta.data

required_cols <- c(
  "seurat_clusters",
  "SingleR_label_fine",
  "SingleR_pruned_fine",
  "SingleR_cluster_label_fine"
)

missing_cols <- required_cols[!required_cols %in% colnames(meta)]
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------------------
# Helper: map fine SingleR labels to conservative level 2 labels
# ----------------------------------------
map_singleR_fine_to_lv2 <- function(lbl) {
  if (is.na(lbl) || !nzchar(lbl)) return("Other")
  
  x <- tolower(trimws(lbl))
  
  # Explicit noisy / irrelevant labels -> Other
  if (grepl("ips", x) ||
      grepl("foreskin", x) ||
      grepl("crl", x) ||
      grepl("embry", x) ||
      grepl("hesc", x)) {
    return("Other")
  }
  
  # Neural / rare unexpected labels -> Other
  if (grepl("schwann", x) || grepl("neuron", x) || grepl("glia", x)) {
    return("Other")
  }
  
  # T / NK
  if (grepl("nk", x)) return("NK_cells")
  
  # Keep only CD4-like T as stable T class here
  if (grepl("cd4", x)) return("CD4_T")
  if (grepl("t_cell", x) && grepl("central_memory|helper|regulatory|naive", x)) return("CD4_T")
  
  # CD8-like labels exist, but if they do not form robust clusters we collapse to Other
  if (grepl("cd8", x)) return("Other")
  if (grepl("t_cell", x) && grepl("cytotoxic|effector", x)) return("Other")
  
  # B lineage
  if (grepl("plasma", x)) return("Plasma_cells")
  if (grepl("b_cell", x) || grepl("\\bb cell\\b", x)) return("B_cells")
  
  # Myeloid
  if (grepl("dendritic", x) || grepl("\\bdc\\b", x)) return("DC")
  if (grepl("macrophage", x)) return("Macrophage")
  
  # Monocyte-like labels are kept conservative here unless clearly present
  if (grepl("monocyte", x)) return("Other")
  if (grepl("\\bcmp\\b", x)) return("Other")
  
  # Non-immune
  if (grepl("endothelial", x)) return("Endothelial")
  
  if (grepl("epithelial", x) ||
      grepl("bronchial", x) ||
      grepl("keratinocyte", x)) {
    return("Epithelial")
  }
  
  # Stromal:
  # allow clearly mesenchymal / smooth muscle / fibro / BM_MSC labels
  if (grepl("smooth_muscle", x) ||
      grepl("smooth muscle", x) ||
      grepl("fibro", x) ||
      grepl("mesenchymal", x) ||
      grepl("\\bbm_msc\\b", x) ||
      grepl("bm_msc", x)) {
    return("Stromal")
  }
  
  return("Other")
}

lv2_levels <- c(
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
# Build cell-level annotation
# Prefer pruned fine labels if available
# ----------------------------------------
meta_dt <- data.table(
  cell = rownames(meta),
  cluster = as.character(meta$seurat_clusters),
  Sample = if ("Sample" %in% colnames(meta)) as.character(meta$Sample) else NA_character_,
  SingleR_label_fine = as.character(meta$SingleR_label_fine),
  SingleR_pruned_fine = as.character(meta$SingleR_pruned_fine),
  SingleR_cluster_label_fine = as.character(meta$SingleR_cluster_label_fine)
)

meta_dt[, fine_label_used := ifelse(
  !is.na(SingleR_pruned_fine) & nzchar(SingleR_pruned_fine),
  SingleR_pruned_fine,
  SingleR_label_fine
)]

meta_dt[, cellcompartment_lv2_auto := vapply(fine_label_used, map_singleR_fine_to_lv2, character(1))]
meta_dt[, cellcompartment_lv2_auto := factor(cellcompartment_lv2_auto, levels = lv2_levels)]

# ----------------------------------------
# Cluster-majority stable level 2 labels
# ----------------------------------------
cluster_lv2_majority <- meta_dt[, {
  x <- as.character(cellcompartment_lv2_auto)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    .(cellcompartment_lv2_cluster = NA_character_)
  } else {
    tab <- sort(table(x), decreasing = TRUE)
    .(cellcompartment_lv2_cluster = names(tab)[1])
  }
}, by = cluster]

# ----------------------------------------
# Manual review overrides for clearly interpretable clusters
# Based on your current review:
# - cluster 15 = Smooth_muscle -> Stromal
# - cluster 19 = NK_cell -> NK_cells
# - clusters 23,24,25,28,35 = BM_MSC-like -> Stromal
# - clusters 10,17,34 remain Other
# ----------------------------------------
cluster_lv2_majority[cluster == "15", cellcompartment_lv2_cluster := "Stromal"]
cluster_lv2_majority[cluster == "19", cellcompartment_lv2_cluster := "NK_cells"]
cluster_lv2_majority[cluster %in% c("23", "24", "25", "28", "35"), cellcompartment_lv2_cluster := "Stromal"]
cluster_lv2_majority[cluster %in% c("10", "17", "34"), cellcompartment_lv2_cluster := "Other"]

# Top-3 lv2 labels per cluster
cluster_lv2_top3 <- meta_dt[, {
  x <- as.character(cellcompartment_lv2_auto)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    .(rank = integer(), label = character(), n = integer(), prop = numeric())
  } else {
    tab <- sort(table(x), decreasing = TRUE)
    tab <- head(tab, 3)
    .(
      rank = seq_along(tab),
      label = names(tab),
      n = as.integer(tab),
      prop = as.numeric(tab) / sum(as.integer(table(x)))
    )
  }
}, by = cluster]

# Merge back
meta_dt <- merge(meta_dt, cluster_lv2_majority, by = "cluster", all.x = TRUE)

# Final stable label for downstream reference building
meta_dt[, cellcompartment_lv2 := factor(cellcompartment_lv2_cluster, levels = lv2_levels)]

# Write back to Seurat
seu$cellcompartment_lv2_auto <- meta_dt$cellcompartment_lv2_auto[match(colnames(seu), meta_dt$cell)]
seu$cellcompartment_lv2 <- meta_dt$cellcompartment_lv2[match(colnames(seu), meta_dt$cell)]

# ----------------------------------------
# Summaries
# ----------------------------------------
cluster_summary <- meta_dt[, .(
  n_cells = .N,
  n_samples = uniqueN(Sample),
  SingleR_cluster_label_fine_majority = names(sort(table(SingleR_cluster_label_fine), decreasing = TRUE))[1],
  fine_label_used_majority = names(sort(table(fine_label_used), decreasing = TRUE))[1],
  cellcompartment_lv2_auto_majority = names(sort(table(cellcompartment_lv2_auto), decreasing = TRUE))[1],
  cellcompartment_lv2_majority = names(sort(table(cellcompartment_lv2), decreasing = TRUE))[1]
), by = cluster]

sample_summary <- meta_dt[, .N, by = .(Sample, cellcompartment_lv2)]
sample_summary_wide <- dcast(
  sample_summary,
  Sample ~ cellcompartment_lv2,
  value.var = "N",
  fill = 0
)

global_summary <- meta_dt[, .N, by = cellcompartment_lv2][order(-N)]

review_dt <- meta_dt[, .(
  n_cells = .N,
  fine_labels_top = paste(head(names(sort(table(fine_label_used), decreasing = TRUE)), 5), collapse = ";"),
  lv2_auto_top = paste(head(names(sort(table(cellcompartment_lv2_auto), decreasing = TRUE)), 5), collapse = ";")
), by = cluster]

review_dt <- merge(review_dt, cluster_summary, by = "cluster", all = TRUE)

# ----------------------------------------
# Save outputs
# ----------------------------------------
out_dir <- file.path(DIR_RESULTS, "scrna", "level2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(seu, file.path(DIR_RESULTS, "scrna", "scrna_seurat_level2.rds"))

fwrite(meta_dt, file.path(out_dir, "level2_cell_annotations.csv"))
fwrite(cluster_summary, file.path(out_dir, "level2_cluster_summary.csv"))
fwrite(cluster_lv2_top3, file.path(out_dir, "level2_cluster_top3.csv"))
fwrite(review_dt, file.path(out_dir, "level2_cluster_review.csv"))
fwrite(sample_summary, file.path(out_dir, "level2_sample_summary_long.csv"))
fwrite(sample_summary_wide, file.path(out_dir, "level2_sample_summary_wide.csv"))
fwrite(global_summary, file.path(out_dir, "level2_global_summary.csv"))

log_msg("Saved: scrna_seurat_level2.rds")
log_msg("Saved: level2_cell_annotations.csv")
log_msg("Saved: level2_cluster_summary.csv")
log_msg("Saved: level2_cluster_top3.csv")
log_msg("Saved: level2_cluster_review.csv")
log_msg("Saved: level2_sample_summary_long.csv")
log_msg("Saved: level2_sample_summary_wide.csv")
log_msg("Saved: level2_global_summary.csv")

print(cluster_summary[order(as.numeric(cluster))])
print(global_summary)