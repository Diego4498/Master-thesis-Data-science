# scripts/24_compare_all_methods.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
})

# -----------------------------
# Load results
# -----------------------------

music  <- fread(file.path(DIR_RESULTS,"deconv","music","music_proportions.csv"))
epic   <- fread(file.path(DIR_RESULTS,"deconv","epic","epic_proportions.csv"))
bisque <- fread(file.path(DIR_RESULTS,"deconv","bisque","bisque_proportions.csv"))
ciber  <- fread(file.path(DIR_RESULTS,"deconv","cibersort","cibersort_proportions.csv"))

# -----------------------------
# Remove CIBERSORT metrics
# -----------------------------

drop_cols <- c("P.value","Correlation","RMSE")

ciber <- ciber[, !colnames(ciber) %in% drop_cols, with=FALSE]

# -----------------------------
# Cell types
# -----------------------------

celltypes <- c(
  "T_cells",
  "B_cells",
  "Myeloid",
  "Epithelial",
  "Stromal",
  "Endothelial",
  "Other"
)

# -----------------------------
# Means per method
# -----------------------------

music_mean  <- colMeans(music[, ..celltypes])
epic_mean   <- colMeans(epic[, ..celltypes])
bisque_mean <- colMeans(bisque[, ..celltypes])
ciber_mean  <- colMeans(ciber[, ..celltypes])

means <- data.table(
  celltype = celltypes,
  MuSiC = music_mean,
  EPIC = epic_mean,
  Bisque = bisque_mean,
  CIBERSORT = ciber_mean
)

print(means)

# -----------------------------
# Heatmap of mean composition
# -----------------------------

mat <- as.matrix(means[, -1])
rownames(mat) <- means$celltype

pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Mean cell type proportions"
)

# -----------------------------
# Prepare long format
# -----------------------------

music$method <- "MuSiC"
epic$method <- "EPIC"
bisque$method <- "Bisque"
ciber$method <- "CIBERSORT"

all <- rbind(
  music[, c("sample",celltypes,"method"), with=FALSE],
  epic[, c("sample",celltypes,"method"), with=FALSE],
  bisque[, c("sample",celltypes,"method"), with=FALSE],
  ciber[, c("sample",celltypes,"method"), with=FALSE]
)

long <- melt(
  all,
  id.vars = c("sample","method"),
  variable.name = "celltype",
  value.name = "proportion"
)

# -----------------------------
# Boxplot comparison
# -----------------------------

ggplot(long, aes(celltype, proportion, fill=method)) +
  geom_boxplot(outlier.size=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(
    title="Deconvolution comparison",
    x="Cell type",
    y="Estimated proportion"
  )

# -----------------------------
# Method correlations
# -----------------------------

cor_music_epic <- cor(
  music[, ..celltypes],
  epic[, ..celltypes]
)

cor_music_bisque <- cor(
  music[, ..celltypes],
  bisque[, ..celltypes]
)

cor_music_ciber <- cor(
  music[, ..celltypes],
  ciber[, ..celltypes]
)

log_msg("Mean proportions:")
print(means)

log_msg("Correlation MuSiC vs EPIC:")
print(diag(cor_music_epic))

log_msg("Correlation MuSiC vs Bisque:")
print(diag(cor_music_bisque))

log_msg("Correlation MuSiC vs CIBERSORT:")
print(diag(cor_music_ciber))