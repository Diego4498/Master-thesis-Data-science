# scripts/03_bulk_make_m_labels.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
})

meta_path <- file.path(DIR_RESULTS, "bulk", "bulk_metadata_clean.csv")
mat_path  <- file.path(DIR_RESULTS, "bulk", "bulk_tpm_matrix.rds")

stopifnot(file.exists(meta_path))
stopifnot(file.exists(mat_path))

meta <- data.table::fread(meta_path)
bulk_mat <- readRDS(mat_path)

stopifnot("Clinical_stage_M" %in% names(meta))

y_raw <- meta[["Clinical_stage_M"]]
y_raw <- trimws(as.character(y_raw))

keep <- y_raw %in% c("M0", "M1")
if (!all(keep)) {
  log_msg("Dropping samples with missing/other Clinical_stage_M values:",
          sum(!keep), "out of", length(y_raw))
}

bulk_mat <- bulk_mat[, keep, drop = FALSE]
meta <- meta[keep]
y <- factor(y_raw[keep], levels = c("M0", "M1"))

stopifnot(ncol(bulk_mat) == length(y))
stopifnot(nrow(meta) == length(y))

# Summary
tab <- table(y)
log_msg("Label counts:", paste(names(tab), as.integer(tab), sep = "=", collapse = " | "))

# Save
out_dir <- file.path(DIR_RESULTS, "bulk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(bulk_mat, file.path(out_dir, "bulk_tpm_matrix_M0M1.rds"))
data.table::fwrite(meta, file.path(out_dir, "bulk_metadata_M0M1.csv"))
saveRDS(y, file.path(out_dir, "bulk_labels_M0M1.rds"))

log_msg("Saved: bulk_tpm_matrix_M0M1.rds")
log_msg("Saved: bulk_metadata_M0M1.csv")
log_msg("Saved: bulk_labels_M0M1.rds")
log_msg("Final bulk dims (genes x samples):", nrow(bulk_mat), "x", ncol(bulk_mat))