# scripts/43_compare_coarse_vs_level2.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# ----------------------------------------
# Paths
# ----------------------------------------
coarse_dir <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso", "comparison")
level2_dir <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso_level2", "comparison")

coarse_overall_path <- file.path(coarse_dir, "validation_overall_metrics_all_methods.csv")
level2_overall_path <- file.path(level2_dir, "level2_overall_metrics_all_methods.csv")

coarse_scatter_path <- file.path(coarse_dir, "truth_method_correlations.csv")
level2_scatter_path <- file.path(level2_dir, "level2_truth_method_scatter_stats.csv")

coarse_cell_path <- file.path(coarse_dir, "validation_truth_method_correlations_by_celltype.csv")
level2_cell_path <- file.path(level2_dir, "level2_truth_method_correlations_by_celltype.csv")

coarse_mae_path <- file.path(coarse_dir, "validation_overall_metrics_all_methods.csv")
level2_mae_path <- file.path(level2_dir, "level2_overall_metrics_all_methods.csv")

stopifnot(file.exists(coarse_overall_path))
stopifnot(file.exists(level2_overall_path))
stopifnot(file.exists(coarse_scatter_path))
stopifnot(file.exists(level2_scatter_path))
stopifnot(file.exists(coarse_cell_path))
stopifnot(file.exists(level2_cell_path))

out_dir <- file.path(DIR_RESULTS, "validation", "comparison_coarse_vs_level2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------
# 1) Load overall metrics
# ----------------------------------------
coarse_overall <- fread(coarse_overall_path)
level2_overall <- fread(level2_overall_path)

coarse_overall[, resolution := "coarse"]
level2_overall[, resolution := "level2"]

overall_all <- rbindlist(list(coarse_overall, level2_overall), use.names = TRUE, fill = TRUE)
fwrite(overall_all, file.path(out_dir, "coarse_vs_level2_overall_metrics_long.csv"))

# wide comparison
overall_cmp <- merge(
  coarse_overall[, .(method, coarse_mae = overall_mae, coarse_rmse = overall_rmse,
                     coarse_pearson = overall_pearson, coarse_spearman = overall_spearman)],
  level2_overall[, .(method, level2_mae = overall_mae, level2_rmse = overall_rmse,
                     level2_pearson = overall_pearson, level2_spearman = overall_spearman)],
  by = "method",
  all = TRUE
)

overall_cmp[, delta_mae := level2_mae - coarse_mae]
overall_cmp[, delta_rmse := level2_rmse - coarse_rmse]
overall_cmp[, delta_pearson := level2_pearson - coarse_pearson]
overall_cmp[, delta_spearman := level2_spearman - coarse_spearman]

fwrite(overall_cmp, file.path(out_dir, "coarse_vs_level2_overall_metrics_wide.csv"))

# ----------------------------------------
# 2) Load scatter stats (mean truth vs mean method)
# ----------------------------------------
coarse_scatter <- fread(coarse_scatter_path)
level2_scatter <- fread(level2_scatter_path)

# standardize method names from comparison strings
extract_method <- function(x) {
  x <- trimws(x)
  x <- sub("^Reference vs ", "", x)
  x <- sub("^Truth vs ", "", x)
  x <- sub("^Level 2 Truth vs ", "", x)
  x
}

coarse_scatter[, method := vapply(comparison, extract_method, character(1))]
level2_scatter[, method := vapply(comparison, extract_method, character(1))]

coarse_scatter[, resolution := "coarse"]
level2_scatter[, resolution := "level2"]

scatter_all <- rbindlist(list(coarse_scatter, level2_scatter), use.names = TRUE, fill = TRUE)
fwrite(scatter_all, file.path(out_dir, "coarse_vs_level2_scatter_stats_long.csv"))

scatter_cmp <- merge(
  coarse_scatter[, .(method,
                     coarse_scatter_pearson = pearson_r,
                     coarse_scatter_spearman = spearman_rho)],
  level2_scatter[, .(method,
                     level2_scatter_pearson = pearson_r,
                     level2_scatter_spearman = spearman_rho)],
  by = "method",
  all = TRUE
)

scatter_cmp[, delta_scatter_pearson := level2_scatter_pearson - coarse_scatter_pearson]
scatter_cmp[, delta_scatter_spearman := level2_scatter_spearman - coarse_scatter_spearman]

fwrite(scatter_cmp, file.path(out_dir, "coarse_vs_level2_scatter_stats_wide.csv"))

# ----------------------------------------
# 3) Barplots overall metrics by resolution
# ----------------------------------------
plot_overall <- melt(
  overall_all,
  id.vars = c("method", "resolution"),
  measure.vars = c("overall_mae", "overall_rmse", "overall_pearson", "overall_spearman"),
  variable.name = "metric",
  value.name = "value"
)

plot_overall$metric <- factor(
  plot_overall$metric,
  levels = c("overall_mae", "overall_rmse", "overall_pearson", "overall_spearman")
)

p_overall <- ggplot(plot_overall, aes(x = method, y = value, fill = resolution)) +
  geom_col(position = "dodge") +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(
    title = "Coarse vs level 2: overall method performance",
    x = NULL,
    y = "Value"
  )

ggsave(
  filename = file.path(out_dir, "coarse_vs_level2_overall_metrics_barplot.png"),
  plot = p_overall,
  width = 10,
  height = 7,
  dpi = 300
)

# ----------------------------------------
# 4) Heatmap of deltas
# ----------------------------------------
delta_mat <- as.matrix(overall_cmp[, .(
  delta_mae,
  delta_rmse,
  delta_pearson,
  delta_spearman
)])
rownames(delta_mat) <- overall_cmp$method

png(file.path(out_dir, "coarse_vs_level2_overall_delta_heatmap.png"),
    width = 1600, height = 1200, res = 200)
pheatmap(
  delta_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  main = "Level 2 minus coarse: overall metric deltas"
)
dev.off()

# ----------------------------------------
# 5) Robustness ranking
# smaller MAE/RMSE increases are better
# smaller Pearson/Spearman drops are better
# ----------------------------------------
robustness <- copy(overall_cmp)

robustness[, rank_delta_mae := frank(delta_mae, ties.method = "average")]
robustness[, rank_delta_rmse := frank(delta_rmse, ties.method = "average")]
robustness[, rank_delta_pearson := frank(-delta_pearson, ties.method = "average")]
robustness[, rank_delta_spearman := frank(-delta_spearman, ties.method = "average")]

robustness[, robustness_score := rowMeans(.SD), .SDcols = c(
  "rank_delta_mae", "rank_delta_rmse", "rank_delta_pearson", "rank_delta_spearman"
)]

setorder(robustness, robustness_score)
fwrite(robustness, file.path(out_dir, "coarse_vs_level2_robustness_ranking.csv"))

# ----------------------------------------
# 6) Cell-type truth correlation comparison
# note: coarse and level2 cell-type universes differ, so we do not merge them directly
# we only save both tables side-by-side for interpretation
# ----------------------------------------
coarse_cell <- fread(coarse_cell_path)
level2_cell <- fread(level2_cell_path)

fwrite(coarse_cell, file.path(out_dir, "coarse_truth_method_correlations_by_celltype.csv"))
fwrite(level2_cell, file.path(out_dir, "level2_truth_method_correlations_by_celltype.csv"))

# ----------------------------------------
# 7) Console output
# ----------------------------------------
log_msg("Saved coarse vs level2 comparison outputs in:", out_dir)

cat("\n=== Overall comparison ===\n")
print(overall_cmp)

cat("\n=== Scatter comparison ===\n")
print(scatter_cmp)

cat("\n=== Robustness ranking ===\n")
print(robustness[, .(
  method,
  delta_mae,
  delta_rmse,
  delta_pearson,
  delta_spearman,
  robustness_score
)])