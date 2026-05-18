#scripts/21_compare_methods_plots.R
source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
})

music  <- fread(file.path(DIR_RESULTS, "deconv", "music",  "music_proportions.csv"))
epic   <- fread(file.path(DIR_RESULTS, "deconv", "epic",   "epic_proportions.csv"))
bisque <- fread(file.path(DIR_RESULTS, "deconv", "bisque", "bisque_proportions.csv"))

celltypes <- c("T_cells", "B_cells", "Myeloid", "Epithelial", "Stromal", "Endothelial", "Other")

# Keep common columns and a sample column with the same name
music_df <- as.data.frame(music[, c("sample", celltypes), with = FALSE])
epic_df  <- as.data.frame(epic[,  c("sample", celltypes), with = FALSE])
bisque_df <- as.data.frame(bisque[, c("sample", celltypes), with = FALSE])

# Align sample order
music_df  <- music_df[order(music_df$sample), ]
epic_df   <- epic_df[order(epic_df$sample), ]
bisque_df <- bisque_df[order(bisque_df$sample), ]

stopifnot(identical(music_df$sample, epic_df$sample))
stopifnot(identical(music_df$sample, bisque_df$sample))

out_dir <- file.path(DIR_RESULTS, "deconv", "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Boxplot by method and cell type
# -----------------------------
music_long <- melt(music_df, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
music_long$method <- "MuSiC"

epic_long <- melt(epic_df, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
epic_long$method <- "EPIC"

bisque_long <- melt(bisque_df, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
bisque_long$method <- "Bisque"

plot_df <- rbind(music_long, epic_long, bisque_long)

p_box <- ggplot(plot_df, aes(x = celltype, y = proportion, fill = method)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  labs(x = "Cell type", y = "Estimated proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(out_dir, "method_comparison_boxplot.png"),
  plot = p_box,
  width = 10,
  height = 6,
  dpi = 300
)

# -----------------------------
# 2) Mean proportions table
# -----------------------------
music_means  <- colMeans(music_df[, celltypes])
epic_means   <- colMeans(epic_df[, celltypes])
bisque_means <- colMeans(bisque_df[, celltypes])

means_df <- data.table(
  celltype = celltypes,
  MuSiC = as.numeric(music_means[celltypes]),
  EPIC = as.numeric(epic_means[celltypes]),
  Bisque = as.numeric(bisque_means[celltypes])
)

fwrite(means_df, file.path(out_dir, "method_mean_proportions.csv"))

# -----------------------------
# 3) Heatmap of mean proportions
# -----------------------------
heat_df <- melt(means_df, id.vars = "celltype", variable.name = "method", value.name = "mean_prop")

p_heat <- ggplot(heat_df, aes(x = method, y = celltype, fill = mean_prop)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", mean_prop)), size = 3) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Mean proportion")

ggsave(
  filename = file.path(out_dir, "method_mean_heatmap.png"),
  plot = p_heat,
  width = 6,
  height = 5,
  dpi = 300
)

# -----------------------------
# 4) Pairwise correlations between methods
# -----------------------------
pairwise_corr <- data.table(
  celltype = celltypes,
  MuSiC_vs_EPIC   = sapply(celltypes, function(ct) cor(music_df[[ct]], epic_df[[ct]])),
  MuSiC_vs_Bisque = sapply(celltypes, function(ct) cor(music_df[[ct]], bisque_df[[ct]])),
  EPIC_vs_Bisque  = sapply(celltypes, function(ct) cor(epic_df[[ct]], bisque_df[[ct]]))
)

fwrite(pairwise_corr, file.path(out_dir, "method_pairwise_correlations.csv"))

# -----------------------------
# 5) Scatter plots for key compartments
# -----------------------------
make_scatter <- function(x, y, xlab, ylab, title, file) {
  df <- data.frame(x = x, y = y)
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.7) +
    theme_bw() +
    labs(x = xlab, y = ylab, title = title)
  ggsave(file.path(out_dir, file), p, width = 5, height = 4, dpi = 300)
}

key_ct <- c("T_cells", "Myeloid", "Epithelial")

for (ct in key_ct) {
  make_scatter(
    music_df[[ct]], epic_df[[ct]],
    paste("MuSiC", ct), paste("EPIC", ct),
    paste(ct, ": MuSiC vs EPIC"),
    paste0("scatter_MuSiC_vs_EPIC_", ct, ".png")
  )
  
  make_scatter(
    music_df[[ct]], bisque_df[[ct]],
    paste("MuSiC", ct), paste("Bisque", ct),
    paste(ct, ": MuSiC vs Bisque"),
    paste0("scatter_MuSiC_vs_Bisque_", ct, ".png")
  )
  
  make_scatter(
    epic_df[[ct]], bisque_df[[ct]],
    paste("EPIC", ct), paste("Bisque", ct),
    paste(ct, ": EPIC vs Bisque"),
    paste0("scatter_EPIC_vs_Bisque_", ct, ".png")
  )
}

log_msg("Saved:", file.path(out_dir, "method_comparison_boxplot.png"))
log_msg("Saved:", file.path(out_dir, "method_mean_heatmap.png"))
log_msg("Saved:", file.path(out_dir, "method_mean_proportions.csv"))
log_msg("Saved:", file.path(out_dir, "method_pairwise_correlations.csv"))

print(means_df)
print(pairwise_corr)