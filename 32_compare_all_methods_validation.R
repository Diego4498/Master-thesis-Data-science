# scripts/32_compare_all_methods_validation.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
})

val_dir <- file.path(DIR_RESULTS, "validation", "pseudobulk_loso")

music_path  <- file.path(val_dir, "music",     "music_loso_predictions_wide.csv")
epic_path   <- file.path(val_dir, "epic",      "epic_loso_predictions_wide.csv")
bisque_path <- file.path(val_dir, "bisque",    "bisque_loso_predictions_wide.csv")
ciber_path  <- file.path(val_dir, "cibersort", "cibersort_loso_predictions_wide.csv")
truth_path  <- file.path(val_dir, "ground_truth_proportions_wide.csv")

stopifnot(file.exists(music_path))
stopifnot(file.exists(epic_path))
stopifnot(file.exists(bisque_path))
stopifnot(file.exists(ciber_path))
stopifnot(file.exists(truth_path))

music  <- fread(music_path)
epic   <- fread(epic_path)
bisque <- fread(bisque_path)
ciber  <- fread(ciber_path)
truth  <- fread(truth_path)

celltypes <- c(
  "T_cells",
  "B_cells",
  "Myeloid",
  "Epithelial",
  "Stromal",
  "Endothelial",
  "Other"
)

# Standardize sample column name
setnames(truth, "Sample", "sample")

stopifnot(all(c("sample", celltypes) %in% colnames(music)))
stopifnot(all(c("sample", celltypes) %in% colnames(epic)))
stopifnot(all(c("sample", celltypes) %in% colnames(bisque)))
stopifnot(all(c("sample", celltypes) %in% colnames(ciber)))
stopifnot(all(c("sample", celltypes) %in% colnames(truth)))

# Align sample order
music  <- music[order(sample)]
epic   <- epic[order(sample)]
bisque <- bisque[order(sample)]
ciber  <- ciber[order(sample)]
truth  <- truth[order(sample)]

stopifnot(identical(music$sample, truth$sample))
stopifnot(identical(epic$sample, truth$sample))
stopifnot(identical(bisque$sample, truth$sample))
stopifnot(identical(ciber$sample, truth$sample))

out_dir <- file.path(val_dir, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Means per method + truth
# -----------------------------
truth_mean  <- colMeans(truth[, ..celltypes])
music_mean  <- colMeans(music[, ..celltypes])
epic_mean   <- colMeans(epic[, ..celltypes])
bisque_mean <- colMeans(bisque[, ..celltypes])
ciber_mean  <- colMeans(ciber[, ..celltypes])

means <- data.table(
  celltype = celltypes,
  Truth = truth_mean,
  MuSiC = music_mean,
  EPIC = epic_mean,
  Bisque = bisque_mean,
  CIBERSORT = ciber_mean
)

fwrite(means, file.path(out_dir, "validation_mean_proportions.csv"))

# -----------------------------
# Heatmap of mean composition
# -----------------------------
mat <- as.matrix(means[, -1])
rownames(mat) <- means$celltype

png(file.path(out_dir, "validation_mean_composition_heatmap.png"),
    width = 1600, height = 1200, res = 200)
pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Mean cell type proportions - pseudobulk validation"
)
dev.off()

# -----------------------------
# Prepare long format
# -----------------------------
truth$method <- "Truth"
music$method <- "MuSiC"
epic$method <- "EPIC"
bisque$method <- "Bisque"
ciber$method <- "CIBERSORT"

all_dt <- rbind(
  truth[,  c("sample", celltypes, "method"), with = FALSE],
  music[,  c("sample", celltypes, "method"), with = FALSE],
  epic[,   c("sample", celltypes, "method"), with = FALSE],
  bisque[, c("sample", celltypes, "method"), with = FALSE],
  ciber[,  c("sample", celltypes, "method"), with = FALSE]
)

long <- melt(
  all_dt,
  id.vars = c("sample", "method"),
  variable.name = "celltype",
  value.name = "proportion"
)

# -----------------------------
# Boxplot comparison
# -----------------------------
p_box <- ggplot(long, aes(celltype, proportion, fill = method)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Pseudobulk validation: deconvolution comparison",
    x = "Cell type",
    y = "Estimated proportion"
  )

ggsave(
  file.path(out_dir, "validation_boxplot_all_methods.png"),
  p_box, width = 10, height = 6, dpi = 300
)

# -----------------------------
# Method correlations vs truth by cell type
# -----------------------------
cor_truth_music <- sapply(celltypes, function(ct) cor(truth[[ct]], music[[ct]], method = "pearson"))
cor_truth_epic <- sapply(celltypes, function(ct) cor(truth[[ct]], epic[[ct]], method = "pearson"))
cor_truth_bisque <- sapply(celltypes, function(ct) cor(truth[[ct]], bisque[[ct]], method = "pearson"))
cor_truth_ciber <- sapply(celltypes, function(ct) cor(truth[[ct]], ciber[[ct]], method = "pearson"))

corr_dt <- data.table(
  celltype = celltypes,
  Truth_vs_MuSiC = cor_truth_music,
  Truth_vs_EPIC = cor_truth_epic,
  Truth_vs_Bisque = cor_truth_bisque,
  Truth_vs_CIBERSORT = cor_truth_ciber
)

fwrite(corr_dt, file.path(out_dir, "validation_truth_method_correlations_by_celltype.csv"))

# -----------------------------
# Overall metrics summary
# -----------------------------
music_overall  <- fread(file.path(val_dir, "music", "music_loso_overall_metrics.csv"))
epic_overall   <- fread(file.path(val_dir, "epic", "epic_loso_overall_metrics.csv"))
bisque_overall <- fread(file.path(val_dir, "bisque", "bisque_loso_overall_metrics.csv"))
ciber_overall  <- fread(file.path(val_dir, "cibersort", "cibersort_loso_overall_metrics.csv"))

overall <- rbindlist(list(music_overall, epic_overall, bisque_overall, ciber_overall), use.names = TRUE, fill = TRUE)
fwrite(overall, file.path(out_dir, "validation_overall_metrics_all_methods.csv"))

log_msg("Saved validation comparison outputs in:", out_dir)
print(means)
print(corr_dt)
print(overall)