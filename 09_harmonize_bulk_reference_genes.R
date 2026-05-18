# scripts/09_harmonize_bulk_reference_genes.R

source("scripts/00_setup_paths.R")

bulk_path <- file.path(DIR_RESULTS, "bulk", "bulk_tpm_matrix_M0M1.rds")
ref_path  <- file.path(DIR_RESULTS, "references", "reference_logexpr_compartments_singleR.rds")

stopifnot(file.exists(bulk_path))
stopifnot(file.exists(ref_path))

bulk <- readRDS(bulk_path)   # genes x samples (TPM)
ref  <- readRDS(ref_path)    # genes x compartments (log-normalized)

bulk_genes <- rownames(bulk)
ref_genes  <- rownames(ref)

common <- intersect(bulk_genes, ref_genes)
stopifnot(length(common) > 0)

bulk2 <- bulk[common, , drop = FALSE]
ref2  <- ref[common, , drop = FALSE]

out_dir <- file.path(DIR_RESULTS, "aligned")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(bulk2, file.path(out_dir, "bulk_tpm_common_genes.rds"))
saveRDS(ref2,  file.path(out_dir, "ref_logexpr_common_genes.rds"))

log_msg("Common genes:", length(common))
log_msg("Aligned dims bulk (genes x samples):", nrow(bulk2), "x", ncol(bulk2))
log_msg("Aligned dims ref  (genes x compartments):", nrow(ref2), "x", ncol(ref2))
log_msg("Saved: bulk_tpm_common_genes.rds")
log_msg("Saved: ref_logexpr_common_genes.rds")