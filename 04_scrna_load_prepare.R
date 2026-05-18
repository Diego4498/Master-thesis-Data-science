# scripts/04_scrna_load_prepare.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
})

in_expr  <- file.path(DIR_DATA, "scrna", "GSE297652_expression_matrix_data.csv.gz")
in_cells <- file.path(DIR_DATA, "scrna", "GSE297652_cell_metadata.csv.gz")
in_genes <- file.path(DIR_DATA, "scrna", "GSE297652_gene_annotation.csv.gz")

stopifnot(file.exists(in_expr))
stopifnot(file.exists(in_cells))
stopifnot(file.exists(in_genes))

log_msg("Reading scRNA expression:", in_expr)
expr <- data.table::fread(in_expr)

log_msg("Reading scRNA cell metadata:", in_cells)
cell_meta <- data.table::fread(in_cells)

log_msg("Reading scRNA gene annotation:", in_genes)
gene_anno <- data.table::fread(in_genes)

log_msg("Expression table dims:", nrow(expr), "x", ncol(expr))
log_msg("Cell metadata rows:", nrow(cell_meta))
log_msg("Gene annotation rows:", nrow(gene_anno))

# Save quick summaries for the next step
out_dir <- file.path(DIR_RESULTS, "scrna")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(names(cell_meta), file.path(out_dir, "cell_metadata_columns.txt"))
writeLines(names(gene_anno), file.path(out_dir, "gene_annotation_columns.txt"))

log_msg("Saved: cell_metadata_columns.txt")
log_msg("Saved: gene_annotation_columns.txt")

# Identify likely sample identifier columns in cell metadata
nm <- names(cell_meta)
cand <- nm[grepl("sample|orig|patient|tissue|gsm|donor|biospec", nm, ignore.case = TRUE)]
log_msg("Candidate sample-id columns:", if (length(cand) == 0) "NONE" else paste(cand, collapse = ", "))

# Print example values for first few candidates
top <- head(cand, 6)
for (col in top) {
  u <- unique(cell_meta[[col]])
  u <- u[!is.na(u)]
  u <- u[seq_len(min(length(u), 15))]
  log_msg("Column:", col, "| example values:", paste(u, collapse = " ; "))
}