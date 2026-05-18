# scripts/42_compare_all_methods_validation_level2.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
})

val_dir <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso_level2")

music_wide_path   <- file.path(val_dir, "music",      "music_level2_loso_predictions_wide.csv")
epic_wide_path    <- file.path(val_dir, "epic",       "epic_level2_loso_predictions_wide.csv")
bisque_wide_path  <- file.path(val_dir, "bisque",     "bisque_level2_loso_predictions_wide.csv")
ciber_wide_path   <- file.path(val_dir, "cibersort",  "cibersort_level2_loso_predictions_wide.csv")
bp_wide_path      <- file.path(val_dir, "bayesprism", "bayesprism_level2_loso_predictions_wide.csv")
truth_wide_path   <- file.path(val_dir, "ground_truth_proportions_wide_level2.csv")

music_overall_path  <- file.path(val_dir, "music",      "music_level2_loso_overall_metrics.csv")
epic_overall_path   <- file.path(val_dir, "epic",       "epic_level2_loso_overall_metrics.csv")
bisque_overall_path <- file.path(val_dir, "bisque",     "bisque_level2_loso_overall_metrics.csv")
ciber_overall_path  <- file.path(val_dir, "cibersort",  "cibersort_level2_loso_overall_metrics.csv")
bp_overall_path     <- file.path(val_dir, "bayesprism", "bayesprism_level2_loso_overall_metrics.csv")

music_cell_path  <- file.path(val_dir, "music",      "music_level2_loso_metrics_by_celltype.csv")
epic_cell_path   <- file.path(val_dir, "epic",       "epic_level2_loso_metrics_by_celltype.csv")
bisque_cell_path <- file.path(val_dir, "bisque",     "bisque_level2_loso_metrics_by_celltype.csv")
ciber_cell_path  <- file.path(val_dir, "cibersort",  "cibersort_level2_loso_metrics_by_celltype.csv")
bp_cell_path     <- file.path(val_dir, "bayesprism", "bayesprism_level2_loso_metrics_by_celltype.csv")

stopifnot(file.exists(music_wide_path))
stopifnot(file.exists(epic_wide_path))
stopifnot(file.exists(bisque_wide_path))
stopifnot(file.exists(ciber_wide_path))
stopifnot(file.exists(bp_wide_path))
stopifnot(file.exists(truth_wide_path))

stopifnot(file.exists(music_overall_path))
stopifnot(file.exists(epic_overall_path))
stopifnot(file.exists(bisque_overall_path))
stopifnot(file.exists(ciber_overall_path))
stopifnot(file.exists(bp_overall_path))

stopifnot(file.exists(music_cell_path))
stopifnot(file.exists(epic_cell_path))
stopifnot(file.exists(bisque_cell_path))
stopifnot(file.exists(ciber_cell_path))
stopifnot(file.exists(bp_cell_path))

celltypes <- c(
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

# -----------------------------
# Load prediction tables
# -----------------------------
music  <- fread(music_wide_path)
epic   <- fread(epic_wide_path)
bisque <- fread(bisque_wide_path)
ciber  <- fread(ciber_wide_path)
bp     <- fread(bp_wide_path)
truth  <- fread(truth_wide_path)

setnames(truth, "Sample", "sample")

# -----------------------------
# Standardize helper
# -----------------------------
standardize_pred_table <- function(dt, celltypes) {
  dt <- as.data.table(dt)
  
  # standardize sample column
  if (!"sample" %in% colnames(dt)) {
    if ("Sample" %in% colnames(dt)) {
      setnames(dt, "Sample", "sample")
    } else {
      setnames(dt, colnames(dt)[1], "sample")
    }
  }
  
  # add any missing celltype columns as zero
  for (ct in celltypes) {
    if (!ct %in% colnames(dt)) {
      dt[, (ct) := 0]
    }
  }
  
  # keep only needed columns and order them
  dt <- dt[, c("sample", celltypes), with = FALSE]
  
  return(dt)
}

music  <- standardize_pred_table(music, celltypes)
epic   <- standardize_pred_table(epic, celltypes)
bisque <- standardize_pred_table(bisque, celltypes)
ciber  <- standardize_pred_table(ciber, celltypes)
bp     <- standardize_pred_table(bp, celltypes)
truth  <- truth[, c("sample", celltypes), with = FALSE]

# align order
music  <- music[order(sample)]
epic   <- epic[order(sample)]
bisque <- bisque[order(sample)]
ciber  <- ciber[order(sample)]
bp     <- bp[order(sample)]
truth  <- truth[order(sample)]

# optional debug prints
cat("\n=== Column names after standardization ===\n")
print(colnames(music))
print(colnames(epic))
print(colnames(bisque))
print(colnames(ciber))
print(colnames(bp))
print(colnames(truth))

cat("\n=== Samples after standardization ===\n")
print(music$sample)
print(epic$sample)
print(bisque$sample)
print(ciber$sample)
print(bp$sample)
print(truth$sample)

stopifnot(identical(music$sample, truth$sample))
stopifnot(identical(epic$sample, truth$sample))
stopifnot(identical(bisque$sample, truth$sample))
stopifnot(identical(ciber$sample, truth$sample))
stopifnot(identical(bp$sample, truth$sample))

# align order
music  <- music[order(sample)]
epic   <- epic[order(sample)]
bisque <- bisque[order(sample)]
ciber  <- ciber[order(sample)]
bp     <- bp[order(sample)]
truth  <- truth[order(sample)]

stopifnot(identical(music$sample, truth$sample))
stopifnot(identical(epic$sample, truth$sample))
stopifnot(identical(bisque$sample, truth$sample))
stopifnot(identical(ciber$sample, truth$sample))
stopifnot(identical(bp$sample, truth$sample))

out_dir <- file.path(val_dir, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Mean proportions table
# -----------------------------
means <- data.table(
  celltype = celltypes,
  Truth = colMeans(truth[, ..celltypes]),
  MuSiC = colMeans(music[, ..celltypes]),
  EPIC = colMeans(epic[, ..celltypes]),
  Bisque = colMeans(bisque[, ..celltypes]),
  CIBERSORT = colMeans(ciber[, ..celltypes]),
  BayesPrism = colMeans(bp[, ..celltypes])
)

fwrite(means, file.path(out_dir, "level2_mean_proportions.csv"))

mean_mat <- as.matrix(means[, -1])
rownames(mean_mat) <- means$celltype

png(file.path(out_dir, "level2_mean_proportions_heatmap.png"),
    width = 2000, height = 1400, res = 200)
pheatmap(
  mean_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  main = "Level 2 mean cell type proportions"
)
dev.off()

# -----------------------------
# 2) Long table for boxplots
# -----------------------------
truth_long <- melt(truth, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
truth_long$method <- "Truth"

music_long <- melt(music, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
music_long$method <- "MuSiC"

epic_long <- melt(epic, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
epic_long$method <- "EPIC"

bisque_long <- melt(bisque, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
bisque_long$method <- "Bisque"

ciber_long <- melt(ciber, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
ciber_long$method <- "CIBERSORT"

bp_long <- melt(bp, id.vars = "sample", variable.name = "celltype", value.name = "proportion")
bp_long$method <- "BayesPrism"

plot_dt <- rbindlist(
  list(truth_long, music_long, epic_long, bisque_long, ciber_long, bp_long),
  use.names = TRUE
)

plot_dt$celltype <- factor(plot_dt$celltype, levels = celltypes)
plot_dt$method <- factor(
  plot_dt$method,
  levels = c("Truth", "MuSiC", "EPIC", "Bisque", "CIBERSORT", "BayesPrism")
)

p_box <- ggplot(plot_dt, aes(x = celltype, y = proportion, fill = method)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  labs(
    title = "Level 2 pseudobulk validation: all methods",
    x = "Cell type",
    y = "Proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(out_dir, "level2_boxplot_all_methods.png"),
  plot = p_box,
  width = 12,
  height = 6,
  dpi = 300
)

plot_methods_only <- plot_dt[plot_dt$method != "Truth", ]

p_box_facet <- ggplot(plot_methods_only, aes(x = method, y = proportion, fill = method)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  theme_bw() +
  labs(
    title = "Level 2 method distributions by cell type",
    x = NULL,
    y = "Estimated proportion"
  ) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(out_dir, "level2_boxplot_faceted_by_celltype.png"),
  plot = p_box_facet,
  width = 13,
  height = 8,
  dpi = 300
)

# -----------------------------
# 3) Overall metrics summary
# -----------------------------
overall <- rbindlist(list(
  fread(music_overall_path),
  fread(epic_overall_path),
  fread(bisque_overall_path),
  fread(ciber_overall_path),
  fread(bp_overall_path)
), use.names = TRUE, fill = TRUE)

setorder(overall, overall_mae)
fwrite(overall, file.path(out_dir, "level2_overall_metrics_all_methods.csv"))

# -----------------------------
# 4) Cell-type metrics summary
# -----------------------------
per_celltype <- rbindlist(list(
  fread(music_cell_path),
  fread(epic_cell_path),
  fread(bisque_cell_path),
  fread(ciber_cell_path),
  fread(bp_cell_path)
), use.names = TRUE, fill = TRUE)

fwrite(per_celltype, file.path(out_dir, "level2_metrics_by_celltype_all_methods.csv"))

mae_dt <- dcast(per_celltype, celltype ~ method, value.var = "mae")
mae_mat <- as.matrix(mae_dt[, -1])
rownames(mae_mat) <- mae_dt$celltype

png(file.path(out_dir, "level2_mae_heatmap.png"),
    width = 2000, height = 1400, res = 200)
pheatmap(
  mae_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  main = "Level 2 MAE by method and cell type"
)
dev.off()

rmse_dt <- dcast(per_celltype, celltype ~ method, value.var = "rmse")
rmse_mat <- as.matrix(rmse_dt[, -1])
rownames(rmse_mat) <- rmse_dt$celltype

png(file.path(out_dir, "level2_rmse_heatmap.png"),
    width = 2000, height = 1400, res = 200)
pheatmap(
  rmse_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  main = "Level 2 RMSE by method and cell type"
)
dev.off()

# -----------------------------
# 5) Truth-vs-method correlations by cell type
# -----------------------------
corr_dt <- data.table(
  celltype = celltypes,
  Truth_vs_MuSiC = sapply(celltypes, function(ct) cor(truth[[ct]], music[[ct]], method = "pearson")),
  Truth_vs_EPIC = sapply(celltypes, function(ct) cor(truth[[ct]], epic[[ct]], method = "pearson")),
  Truth_vs_Bisque = sapply(celltypes, function(ct) cor(truth[[ct]], bisque[[ct]], method = "pearson")),
  Truth_vs_CIBERSORT = sapply(celltypes, function(ct) cor(truth[[ct]], ciber[[ct]], method = "pearson")),
  Truth_vs_BayesPrism = sapply(celltypes, function(ct) cor(truth[[ct]], bp[[ct]], method = "pearson"))
)

fwrite(corr_dt, file.path(out_dir, "level2_truth_method_correlations_by_celltype.csv"))

corr_mat <- as.matrix(corr_dt[, -1])
rownames(corr_mat) <- corr_dt$celltype

png(file.path(out_dir, "level2_truth_method_correlation_heatmap.png"),
    width = 2000, height = 1400, res = 200)
pheatmap(
  corr_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  main = "Level 2 Pearson correlation with truth by cell type"
)
dev.off()

# -----------------------------
# 6) Scatter plots
# -----------------------------
make_scatter <- function(df, x_col, y_col, title, file_name) {
  x <- df[[x_col]]
  y <- df[[y_col]]
  
  pear <- cor.test(x, y, method = "pearson")
  spear <- cor.test(x, y, method = "spearman")
  
  label_txt <- paste0(
    "Pearson r = ", round(pear$estimate, 3),
    "\nPearson p = ", signif(pear$p.value, 3),
    "\nSpearman rho = ", round(spear$estimate, 3),
    "\nSpearman p = ", signif(spear$p.value, 3)
  )
  
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(size = 3) +
    geom_text(aes(label = celltype), vjust = -0.7, size = 3) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    labs(title = title, x = x_col, y = y_col) +
    annotate(
      "text",
      x = min(x, na.rm = TRUE),
      y = max(y, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      label = label_txt,
      size = 3.5
    )
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
  
  data.table(
    comparison = title,
    pearson_r = unname(pear$estimate),
    pearson_p = pear$p.value,
    spearman_rho = unname(spear$estimate),
    spearman_p = spear$p.value
  )
}

scatter_stats <- rbindlist(list(
  make_scatter(means, "Truth", "MuSiC", "Level 2 Truth vs MuSiC", "level2_scatter_truth_vs_music.png"),
  make_scatter(means, "Truth", "EPIC", "Level 2 Truth vs EPIC", "level2_scatter_truth_vs_epic.png"),
  make_scatter(means, "Truth", "Bisque", "Level 2 Truth vs Bisque", "level2_scatter_truth_vs_bisque.png"),
  make_scatter(means, "Truth", "CIBERSORT", "Level 2 Truth vs CIBERSORT", "level2_scatter_truth_vs_cibersort.png"),
  make_scatter(means, "Truth", "BayesPrism", "Level 2 Truth vs BayesPrism", "level2_scatter_truth_vs_bayesprism.png")
))

fwrite(scatter_stats, file.path(out_dir, "level2_truth_method_scatter_stats.csv"))

# -----------------------------
# 7) Console output
# -----------------------------
log_msg("Saved level2 comparison outputs in:", out_dir)
print(means)
print(overall)
print(per_celltype)
print(corr_dt)
print(scatter_stats)