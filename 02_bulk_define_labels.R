# scripts/02_bulk_define_labels.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
})

meta_path <- file.path(DIR_RESULTS, "bulk", "bulk_metadata_clean.csv")
stopifnot(file.exists(meta_path))

meta <- data.table::fread(meta_path)

# Find a metastasis column (robust search by column name)
nm <- names(meta)
hits <- nm[grepl("met", nm, ignore.case = TRUE) |
             grepl("\\bm\\b", nm, ignore.case = TRUE) |
             grepl("stage", nm, ignore.case = TRUE) |
             grepl("tnm", nm, ignore.case = TRUE)]

if (length(hits) == 0) {
  stop("No obvious metastasis-related column found in metadata. Columns are: ",
       paste(nm, collapse = ", "))
}

log_msg("Candidate label columns:", paste(hits, collapse = ", "))

# Print unique values for the top candidates
top <- head(hits, 6)
for (col in top) {
  u <- unique(meta[[col]])
  u <- u[!is.na(u)]
  u <- u[seq_len(min(length(u), 20))]
  log_msg("Column:", col, "| example values:", paste(u, collapse = " ; "))
}

out_dir <- file.path(DIR_RESULTS, "bulk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Save the candidates list for traceability
cand_path <- file.path(out_dir, "bulk_label_candidates.txt")
writeLines(c(hits), cand_path)
log_msg("Saved:", cand_path)