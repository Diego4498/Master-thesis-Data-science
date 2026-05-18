# scripts/26_scatter_correlations.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

comp_path <- file.path(DIR_RESULTS, "deconv", "comparison", "reference_vs_methods_mean_proportions.csv")
stopifnot(file.exists(comp_path))

comp <- fread(comp_path)

out_dir <- file.path(DIR_RESULTS, "deconv", "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper function (modern ggplot)
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
    labs(
      title = title,
      x = x_col,
      y = y_col
    ) +
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

# -----------------------------
# Run plots
# -----------------------------
stats_list <- list()

stats_list[[1]] <- make_scatter(comp, "Reference_scRNA", "MuSiC",
                                "Reference vs MuSiC",
                                "scatter_reference_vs_music.png")

stats_list[[2]] <- make_scatter(comp, "Reference_scRNA", "EPIC",
                                "Reference vs EPIC",
                                "scatter_reference_vs_epic.png")

stats_list[[3]] <- make_scatter(comp, "Reference_scRNA", "Bisque",
                                "Reference vs Bisque",
                                "scatter_reference_vs_bisque.png")

stats_list[[4]] <- make_scatter(comp, "Reference_scRNA", "CIBERSORT",
                                "Reference vs CIBERSORT",
                                "scatter_reference_vs_cibersort.png")

stats_df <- rbindlist(stats_list)

stats_file <- file.path(out_dir, "reference_method_correlations.csv")
fwrite(stats_df, stats_file)

log_msg("Saved scatter plots + stats")
print(stats_df)
