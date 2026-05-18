# scripts/25_reference_comparison.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(ggplot2)
  library(reshape2)
})

# -----------------------------
# Load scRNA reference
# -----------------------------
seu_path <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_annot_singleR.rds")
stopifnot(file.exists(seu_path))

seu <- readRDS(seu_path)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

meta <- seu@meta.data
stopifnot("SingleR_cluster_label" %in% colnames(meta))

# Rebuild same compartments used in the deconvolution pipeline
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
  levels = c("T_cells", "B_cells", "Myeloid", "Epithelial", "Stromal", "Endothelial", "Other")
)

seu$cellcompartment <- meta$cellcompartment

celltypes <- c("T_cells", "B_cells", "Myeloid", "Epithelial", "Stromal", "Endothelial", "Other")

# Reference proportions from scRNA
ref_tab <- prop.table(table(seu$cellcompartment))
ref_vec <- as.numeric(ref_tab[celltypes])
names(ref_vec) <- celltypes

# -----------------------------
# Load deconvolution outputs
# -----------------------------
music_path  <- file.path(DIR_RESULTS, "deconv", "music", "music_proportions.csv")
epic_path   <- file.path(DIR_RESULTS, "deconv", "epic", "epic_proportions.csv")
bisque_path <- file.path(DIR_RESULTS, "deconv", "bisque", "bisque_proportions.csv")
ciber_path  <- file.path(DIR_RESULTS, "deconv", "cibersort", "cibersort_proportions.csv")

stopifnot(file.exists(music_path))
stopifnot(file.exists(epic_path))
stopifnot(file.exists(bisque_path))
stopifnot(file.exists(ciber_path))

music  <- fread(music_path)
epic   <- fread(epic_path)
bisque <- fread(bisque_path)
ciber  <- fread(ciber_path)

# CIBERSORT may include fit metrics
drop_cols <- c("P.value", "Correlation", "RMSE")
ciber <- ciber[, !colnames(ciber) %in% drop_cols, with = FALSE]

# -----------------------------
# Mean proportions by method
# -----------------------------
stopifnot(all(c("sample", celltypes) %in% colnames(music)))
stopifnot(all(c("sample", celltypes) %in% colnames(epic)))
stopifnot(all(c("sample", celltypes) %in% colnames(bisque)))
stopifnot(all(c("sample", celltypes) %in% colnames(ciber)))

music_mean  <- colMeans(music[, ..celltypes])
epic_mean   <- colMeans(epic[, ..celltypes])
bisque_mean <- colMeans(bisque[, ..celltypes])
ciber_mean  <- colMeans(ciber[, ..celltypes])

comparison <- data.table(
  celltype = celltypes,
  Reference_scRNA = as.numeric(ref_vec[celltypes]),
  MuSiC = as.numeric(music_mean[celltypes]),
  EPIC = as.numeric(epic_mean[celltypes]),
  Bisque = as.numeric(bisque_mean[celltypes]),
  CIBERSORT = as.numeric(ciber_mean[celltypes])
)

# -----------------------------
# Save table
# -----------------------------
out_dir <- file.path(DIR_RESULTS, "deconv", "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_table <- file.path(out_dir, "reference_vs_methods_mean_proportions.csv")
fwrite(comparison, out_table)

# -----------------------------
# Plot
# -----------------------------
plot_df <- melt(
  comparison,
  id.vars = "celltype",
  variable.name = "source",
  value.name = "proportion"
)

plot_df$celltype <- factor(plot_df$celltype, levels = celltypes)
plot_df$source <- factor(
  plot_df$source,
  levels = c("Reference_scRNA", "MuSiC", "EPIC", "Bisque", "CIBERSORT")
)

p <- ggplot(plot_df, aes(x = celltype, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(
    title = "Reference scRNA proportions vs deconvolution methods",
    x = "Cell type",
    y = "Mean proportion",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_file <- file.path(out_dir, "reference_vs_methods_barplot.png")
ggsave(plot_file, p, width = 10, height = 6, dpi = 300)

# -----------------------------
# Console output
# -----------------------------
log_msg("Saved:", out_table)
log_msg("Saved:", plot_file)
print(comparison)
