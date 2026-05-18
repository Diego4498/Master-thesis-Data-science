# scripts/23_run_cibersort.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(CIBERSORT)
})

set.seed(1)

sig_file  <- file.path(DIR_RESULTS, "deconv", "cibersort", "cibersort_signature_matrix.txt")
bulk_file <- file.path(DIR_RESULTS, "aligned", "bulk_tpm_common_genes.rds")

stopifnot(file.exists(sig_file))
stopifnot(file.exists(bulk_file))

bulk <- readRDS(bulk_file)
bulk <- as.matrix(bulk)

# align bulk genes to signature genes
sig_dt <- fread(sig_file, sep = "\t")
sig_genes <- sig_dt$Gene

common <- intersect(rownames(bulk), sig_genes)
stopifnot(length(common) > 0)

bulk2 <- bulk[common, , drop = FALSE]

mix_file <- file.path(DIR_RESULTS, "deconv", "cibersort", "bulk_mixture.txt")

fwrite(
  data.table(Gene = rownames(bulk2), bulk2),
  mix_file,
  sep = "\t"
)

log_msg("Running CIBERSORT")
log_msg("Mixture dims:", nrow(bulk2), "genes x", ncol(bulk2), "samples")

res <- cibersort(
  sig_file,
  mix_file,
  perm = 0,
  QN = TRUE
)

res <- as.data.table(res)
res$sample <- rownames(res)
setcolorder(res, c("sample", setdiff(colnames(res), "sample")))

out_dir <- file.path(DIR_RESULTS, "deconv", "cibersort")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "cibersort_proportions.csv")
fwrite(res, out_file)

log_msg("Saved:", out_file)
log_msg("Output dims:", nrow(res), "samples x", ncol(res), "columns")

cibersort <- data.table::fread("results/deconv/cibersort/cibersort_proportions.csv")

colnames(cibersort)

summary(cibersort)