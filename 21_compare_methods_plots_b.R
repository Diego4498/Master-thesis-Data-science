# scripts/21_compare_methods_plots_b.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(grid)
})

# -----------------------------
# Load deconvolution outputs
# -----------------------------
music_path  <- file.path(DIR_RESULTS, "deconv", "music",  "music_proportions.csv")
epic_path   <- file.path(DIR_RESULTS, "deconv", "epic",   "epic_proportions.csv")
bisque_path <- file.path(DIR_RESULTS, "deconv", "bisque", "bisque_proportions.csv")

stopifnot(file.exists(music_path))
stopifnot(file.exists(epic_path))
stopifnot(file.exists(bisque_path))

music  <- data.table::fread(music_path)
epic   <- data.table::fread(epic_path)
bisque <- data.table::fread(bisque_path)

# -----------------------------
# Standardize formats
# -----------------------------
# MuSiC already: rows = samples, cols = cell types + sample
stopifnot("sample" %in% colnames(music))

# EPIC already: rows = samples, cols = cell types + sample
stopifnot("sample" %in% colnames(epic))

# Bisque should already be fixed, but if not, repair it
if (!"sample" %in% colnames(bisque)) {
  # try to repair orientation if cell types are rownames in first column
  bisque_raw <- read.csv(bisque_path, row.names = 1, check.names = FALSE)
  bisque_fix <- as.data.frame(t(as.matrix(bisque_raw)))
  bisque_fix$sample <- rownames(bisque_fix)
  bisque <- as.data.table(bisque_fix)
}

# -----------------------------
# Common cell types
# -----------------------------
celltypes <- c("T_cells", "B_cells", "Myeloid", "Epithelial", "Stromal", "Endothelial", "Other")

stopifnot(all(c("sample", celltypes) %in% colnames(music)))
stopifnot(all(c("sample", celltypes) %in% colnames(epic)))
stopifnot(all(c("sample", celltypes) %in% colnames(bisque)))

music_df  <- as.data.frame(music[, c("sample", celltypes), with = FALSE])
epic_df   <- as.data.frame(epic[,  c("sample", celltypes), with = FALSE])
bisque_df <- as.data.frame(bisque[, c("sample", celltypes), with = FALSE])

# align order by sample
music_df  <- music_df[order(music_df$sample), ]
epic_df   <- epic_df[order(epic_df$sample), ]
bisque_df <- bisque_df[order(bisque_df$sample), ]

stopifnot(identical(music_df$sample, epic_df$sample))
stopifnot(identical(music_df$sample, bisque_df$sample))

samples <- music_df$sample

# -----------------------------
# Output dir
# -----------------------------
out_dir <- file.path(DIR_RESULTS, "deconv", "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Save mean/median tables
# -----------------------------
means_df <- data.table(
  celltype = celltypes,
  MuSiC  = colMeans(music_df[,  celltypes, drop = FALSE]),
  EPIC   = colMeans(epic_df[,   celltypes, drop = FALSE]),
  Bisque = colMeans(bisque_df[, celltypes, drop = FALSE])
)

medians_df <- data.table(
  celltype = celltypes,
  MuSiC  = apply(music_df[,  celltypes, drop = FALSE], 2, median),
  EPIC   = apply(epic_df[,   celltypes, drop = FALSE], 2, median),
  Bisque = apply(bisque_df[, celltypes, drop = FALSE], 2, median)
)

data.table::fwrite(means_df,   file.path(out_dir, "method_mean_proportions.csv"))
data.table::fwrite(medians_df, file.path(out_dir, "method_median_proportions.csv"))

# -----------------------------
# Pairwise correlations by cell type
# -----------------------------
pairwise_corr <- data.table(
  celltype = celltypes,
  MuSiC_vs_EPIC   = sapply(celltypes, function(ct) cor(music_df[[ct]],  epic_df[[ct]],   method = "spearman")),
  MuSiC_vs_Bisque = sapply(celltypes, function(ct) cor(music_df[[ct]],  bisque_df[[ct]], method = "spearman")),
  EPIC_vs_Bisque  = sapply(celltypes, function(ct) cor(epic_df[[ct]],   bisque_df[[ct]], method = "spearman"))
)

data.table::fwrite(pairwise_corr, file.path(out_dir, "method_pairwise_correlations.csv"))

# -----------------------------
# Long format for plotting
# -----------------------------
music_long <- melt(music_df, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
music_long$method <- "MuSiC"

epic_long <- melt(epic_df, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
epic_long$method <- "EPIC"

bisque_long <- melt(bisque_df, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
bisque_long$method <- "Bisque"

plot_df <- rbind(music_long, epic_long, bisque_long)
plot_df$celltype <- factor(plot_df$celltype, levels = celltypes)
plot_df$method <- factor(plot_df$method, levels = c("MuSiC", "EPIC", "Bisque"))

# -----------------------------
# 1) Boxplot by method and cell type
# -----------------------------
p_box <- ggplot(plot_df, aes(x = celltype, y = proportion, fill = method)) +
  geom_boxplot(outlier.size = 0.4) +
  theme_bw() +
  labs(x = "Cell type", y = "Estimated proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(out_dir, "boxplot_methods_by_celltype.png"),
  plot = p_box,
  width = 10,
  height = 6,
  dpi = 300
)

# -----------------------------
# 2) Faceted boxplot, one panel per cell type
# -----------------------------
p_box_facet <- ggplot(plot_df, aes(x = method, y = proportion, fill = method)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  theme_bw() +
  labs(x = NULL, y = "Estimated proportion") +
  theme(legend.position = "none")

ggsave(
  filename = file.path(out_dir, "boxplot_faceted_by_celltype.png"),
  plot = p_box_facet,
  width = 11,
  height = 7,
  dpi = 300
)

# -----------------------------
# 3) Stacked barplots per method
# -----------------------------
make_stacked_barplot <- function(df, method_name, file_name) {
  long_df <- melt(df, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
  long_df$celltype <- factor(long_df$celltype, levels = celltypes)
  
  p <- ggplot(long_df, aes(x = sample, y = proportion, fill = celltype)) +
    geom_bar(stat = "identity", width = 1) +
    theme_bw() +
    labs(x = "Samples", y = "Estimated proportion", title = method_name) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot = p,
    width = 12,
    height = 5,
    dpi = 300
  )
}

make_stacked_barplot(music_df,  "MuSiC",  "stackedbar_music.png")
make_stacked_barplot(epic_df,   "EPIC",   "stackedbar_epic.png")
make_stacked_barplot(bisque_df, "Bisque", "stackedbar_bisque.png")

# -----------------------------
# 4) Heatmap of mean proportions
# -----------------------------
mean_mat <- as.matrix(data.frame(
  MuSiC  = means_df$MuSiC,
  EPIC   = means_df$EPIC,
  Bisque = means_df$Bisque,
  row.names = means_df$celltype
))

png(file.path(out_dir, "heatmap_mean_proportions.png"), width = 1600, height = 1200, res = 200)
pheatmap(
  mean_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  fontsize = 11,
  main = "Mean estimated proportions by method"
)
dev.off()

# -----------------------------
# 5) Heatmap of pairwise correlations
# -----------------------------
corr_mat <- as.matrix(data.frame(
  MuSiC_vs_EPIC   = pairwise_corr$MuSiC_vs_EPIC,
  MuSiC_vs_Bisque = pairwise_corr$MuSiC_vs_Bisque,
  EPIC_vs_Bisque  = pairwise_corr$EPIC_vs_Bisque,
  row.names = pairwise_corr$celltype
))

png(file.path(out_dir, "heatmap_pairwise_correlations.png"), width = 1600, height = 1200, res = 200)
pheatmap(
  corr_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 11,
  main = "Spearman correlation between methods"
)
dev.off()

# -----------------------------
# 6) Scatter plots for key compartments
# -----------------------------
make_scatter <- function(x, y, xlab, ylab, title, file_name) {
  df <- data.frame(x = x, y = y)
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    labs(x = xlab, y = ylab, title = title)
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot = p,
    width = 5,
    height = 4,
    dpi = 300
  )
}

key_ct <- c("T_cells", "Myeloid", "Epithelial", "Stromal")

for (ct in key_ct) {
  make_scatter(
    music_df[[ct]], epic_df[[ct]],
    paste("MuSiC", ct), paste("EPIC", ct),
    paste(ct, ": MuSiC vs EPIC"),
    paste0("scatter_music_vs_epic_", ct, ".png")
  )
  
  make_scatter(
    music_df[[ct]], bisque_df[[ct]],
    paste("MuSiC", ct), paste("Bisque", ct),
    paste(ct, ": MuSiC vs Bisque"),
    paste0("scatter_music_vs_bisque_", ct, ".png")
  )
  
  make_scatter(
    epic_df[[ct]], bisque_df[[ct]],
    paste("EPIC", ct), paste("Bisque", ct),
    paste(ct, ": EPIC vs Bisque"),
    paste0("scatter_epic_vs_bisque_", ct, ".png")
  )
}

# -----------------------------
# 7) PCA on compositions per method
# -----------------------------
make_pca_plot <- function(df, method_name, file_name) {
  mat <- as.matrix(df[, celltypes, drop = FALSE])
  pca <- prcomp(mat, center = TRUE, scale. = TRUE)
  
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    sample = df$sample
  )
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.8) +
    theme_bw() +
    labs(title = paste("PCA of estimated compositions -", method_name))
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot = p,
    width = 5,
    height = 4,
    dpi = 300
  )
}

make_pca_plot(music_df,  "MuSiC",  "pca_music.png")
make_pca_plot(epic_df,   "EPIC",   "pca_epic.png")
make_pca_plot(bisque_df, "Bisque", "pca_bisque.png")

# -----------------------------
# 8) Rank correlation of samples by epithelial fraction
# -----------------------------
rank_df <- data.table(
  sample = samples,
  MuSiC_Epithelial  = rank(-music_df$Epithelial, ties.method = "average"),
  EPIC_Epithelial   = rank(-epic_df$Epithelial, ties.method = "average"),
  Bisque_Epithelial = rank(-bisque_df$Epithelial, ties.method = "average")
)

data.table::fwrite(rank_df, file.path(out_dir, "epithelial_ranks_by_method.csv"))

# -----------------------------
# 9) Save standardized tables
# -----------------------------
data.table::fwrite(as.data.table(music_df),  file.path(out_dir, "music_standardized.csv"))
data.table::fwrite(as.data.table(epic_df),   file.path(out_dir, "epic_standardized.csv"))
data.table::fwrite(as.data.table(bisque_df), file.path(out_dir, "bisque_standardized.csv"))

# -----------------------------
# 10) Console output
# -----------------------------
log_msg("Saved:", file.path(out_dir, "boxplot_methods_by_celltype.png"))
log_msg("Saved:", file.path(out_dir, "boxplot_faceted_by_celltype.png"))
log_msg("Saved:", file.path(out_dir, "stackedbar_music.png"))
log_msg("Saved:", file.path(out_dir, "stackedbar_epic.png"))
log_msg("Saved:", file.path(out_dir, "stackedbar_bisque.png"))
log_msg("Saved:", file.path(out_dir, "heatmap_mean_proportions.png"))
log_msg("Saved:", file.path(out_dir, "heatmap_pairwise_correlations.png"))
log_msg("Saved:", file.path(out_dir, "method_mean_proportions.csv"))
log_msg("Saved:", file.path(out_dir, "method_pairwise_correlations.csv"))

print(means_df)
print(pairwise_corr)
