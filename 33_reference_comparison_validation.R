# scripts/33_reference_comparison_validation.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
})

val_dir <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso")

truth_path  <- file.path(val_dir, "ground_truth_proportions_wide.csv")
music_path  <- file.path(val_dir, "music",     "music_loso_predictions_wide.csv")
epic_path   <- file.path(val_dir, "epic",      "epic_loso_predictions_wide.csv")
bisque_path <- file.path(val_dir, "bisque",    "bisque_loso_predictions_wide.csv")
ciber_path  <- file.path(val_dir, "cibersort", "cibersort_loso_predictions_wide.csv")

stopifnot(file.exists(truth_path))
stopifnot(file.exists(music_path))
stopifnot(file.exists(epic_path))
stopifnot(file.exists(bisque_path))
stopifnot(file.exists(ciber_path))

truth  <- fread(truth_path)
music  <- fread(music_path)
epic   <- fread(epic_path)
bisque <- fread(bisque_path)
ciber  <- fread(ciber_path)

celltypes <- c("T_cells", "B_cells", "Myeloid", "Epithelial", "Stromal", "Endothelial", "Other")
setnames(truth, "Sample", "sample")

stopifnot(all(c("sample", celltypes) %in% colnames(truth)))
stopifnot(all(c("sample", celltypes) %in% colnames(music)))
stopifnot(all(c("sample", celltypes) %in% colnames(epic)))
stopifnot(all(c("sample", celltypes) %in% colnames(bisque)))
stopifnot(all(c("sample", celltypes) %in% colnames(ciber)))

truth_mean  <- colMeans(truth[, ..celltypes])
music_mean  <- colMeans(music[, ..celltypes])
epic_mean   <- colMeans(epic[, ..celltypes])
bisque_mean <- colMeans(bisque[, ..celltypes])
ciber_mean  <- colMeans(ciber[, ..celltypes])

comparison <- data.table(
  celltype = celltypes,
  Truth_pseudobulk = as.numeric(truth_mean[celltypes]),
  MuSiC = as.numeric(music_mean[celltypes]),
  EPIC = as.numeric(epic_mean[celltypes]),
  Bisque = as.numeric(bisque_mean[celltypes]),
  CIBERSORT = as.numeric(ciber_mean[celltypes])
)

out_dir <- file.path(val_dir, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_table <- file.path(out_dir, "truth_vs_methods_mean_proportions.csv")
fwrite(comparison, out_table)

plot_df <- melt(
  comparison,
  id.vars = "celltype",
  variable.name = "source",
  value.name = "proportion"
)

plot_df$celltype <- factor(plot_df$celltype, levels = celltypes)
plot_df$source <- factor(
  plot_df$source,
  levels = c("Truth_pseudobulk", "MuSiC", "EPIC", "Bisque", "CIBERSORT")
)

p <- ggplot(plot_df, aes(x = celltype, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(
    title = "Pseudobulk truth proportions vs deconvolution methods",
    x = "Cell type",
    y = "Mean proportion",
    fill = NULL
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_file <- file.path(out_dir, "truth_vs_methods_barplot.png")
ggsave(plot_file, p, width = 10, height = 6, dpi = 300)

log_msg("Saved:", out_table)
log_msg("Saved:", plot_file)
print(comparison)