# scripts/01_bulk_load_prepare.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
})

in_expr <- file.path(DIR_DATA, "bulk", "GSE297742_YUHS_TPM_EXPR_Biopsy_processed.csv.gz")
in_meta <- file.path(DIR_DATA, "bulk", "GSE297742_YUHS_TPM_EXPR_Biopsy_clinical_metadata.csv.gz")

stopifnot(file.exists(in_expr))
stopifnot(file.exists(in_meta))

log_msg("Reading bulk expression:", in_expr)
expr <- data.table::fread(in_expr)

log_msg("Reading bulk metadata:", in_meta)
meta <- data.table::fread(in_meta)

# Expectation: first column is gene, remaining columns are samples
gene_col <- names(expr)[1]
genes <- expr[[gene_col]]

# Collapse duplicated gene names by mean TPM (required: unique rownames)
if (anyDuplicated(genes)) {
  log_msg("Duplicated gene names detected. Collapsing by mean TPM.")
  expr <- expr[, lapply(.SD, mean), by = gene_col, .SDcols = setdiff(names(expr), gene_col)]
}
genes <- expr[[gene_col]]
stopifnot(!anyDuplicated(genes))

sample_cols <- setdiff(names(expr), gene_col)
stopifnot(length(sample_cols) > 0)

# Build matrix: genes x samples
bulk_mat <- as.matrix(expr[, ..sample_cols])
rownames(bulk_mat) <- genes

# Basic sanity checks
stopifnot(is.numeric(bulk_mat))
stopifnot(all(is.finite(bulk_mat)))
stopifnot(all(bulk_mat >= 0))

# Ensure sample IDs are unique
stopifnot(!anyDuplicated(colnames(bulk_mat)))

# Standardize metadata key: try common sample id columns
meta_cols <- names(meta)
candidate_keys <- c("sample", "sample_id", "Sample", "SampleID", "geo_accession", "GEO_Accession")
key <- candidate_keys[candidate_keys %in% meta_cols][1]

if (is.na(key) || !nzchar(key)) {
  stop("No sample ID column found in metadata. Columns are: ", paste(meta_cols, collapse = ", "))
}

meta_key <- meta[[key]]
if (anyDuplicated(meta_key)) {
  stop("Metadata sample key column has duplicates: ", key)
}

# Keep only samples present in both
common_samples <- intersect(colnames(bulk_mat), meta_key)
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between expression matrix and metadata using key: ", key)
}

bulk_mat <- bulk_mat[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta_key)]

# Save outputs
out_dir <- file.path(DIR_RESULTS, "bulk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_rds <- file.path(out_dir, "bulk_tpm_matrix.rds")
out_meta <- file.path(out_dir, "bulk_metadata_clean.csv")

saveRDS(bulk_mat, out_rds)
data.table::fwrite(meta, out_meta)

log_msg("Saved:", out_rds)
log_msg("Saved:", out_meta)
log_msg("Bulk matrix dims (genes x samples):", nrow(bulk_mat), "x", ncol(bulk_mat))