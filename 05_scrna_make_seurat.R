# scripts/05_scrna_make_seurat.R

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(Seurat)
})

DIR_10X <- file.path(DIR_DATA, "scrna", "10x")
stopifnot(dir.exists(DIR_10X))

sample_dirs <- list.dirs(DIR_10X, full.names = TRUE, recursive = FALSE)
stopifnot(length(sample_dirs) > 0)

sample_dirs <- sample_dirs[grepl("^mPCa_M[1-5]_filtered_feature_bc_matrix$", basename(sample_dirs))]
stopifnot(length(sample_dirs) > 0)

sample_dirs <- sample_dirs[order(sample_dirs)]

objs <- list()

for (d in sample_dirs) {
  sample_id <- sub("_filtered_feature_bc_matrix$", "", basename(d))  # mPCa_M1 etc.
  
  log_msg("Reading 10x:", d)
  counts <- Read10X(data.dir = d)
  
  seu <- CreateSeuratObject(counts = counts, project = "GSE297652_10x")
  seu$Sample <- sample_id
  
  mt_genes <- grep("^MT-", rownames(seu), value = TRUE)
  if (length(mt_genes) > 0) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    log_msg("Computed percent.mt for", sample_id)
  } else {
    log_msg("No MT- genes found for", sample_id)
  }
  
  objs[[sample_id]] <- seu
  log_msg("Loaded", sample_id, "| cells:", ncol(seu), "| genes:", nrow(seu))
}

log_msg("Merging samples")
seu_merged <- Reduce(function(x, y) merge(x, y), objs)

out_dir <- file.path(DIR_RESULTS, "scrna")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_rds <- file.path(out_dir, "scrna_seurat_raw.rds")
saveRDS(seu_merged, out_rds)

tab <- table(seu_merged$Sample)
log_msg("Cells per Sample:", paste(names(tab), as.integer(tab), sep = "=", collapse = " | "))
log_msg("Merged dims (genes x cells):", nrow(seu_merged), "x", ncol(seu_merged))
log_msg("Saved:", out_rds)

# QC filtering
log_msg("QC filtering")

# keep only cells with at least 200 detected genes
seu_qc <- subset(seu_merged, subset = nFeature_RNA >= 200)

# remove extreme high-feature cells (likely doublets); threshold by 99th percentile
nf <- seu_qc$nFeature_RNA
upper_nf <- as.numeric(stats::quantile(nf, probs = 0.99, na.rm = TRUE))
seu_qc <- subset(seu_qc, subset = nFeature_RNA <= upper_nf)

# remove high mitochondrial cells; fixed threshold
seu_qc <- subset(seu_qc, subset = percent.mt <= 20)

log_msg("QC thresholds:",
        "min_nFeature=200",
        "| max_nFeature(p99)=", round(upper_nf),
        "| max_percent_mt=20")

tab_qc <- table(seu_qc$Sample)
log_msg("Cells per Sample after QC:", paste(names(tab_qc), as.integer(tab_qc), sep="=", collapse=" | "))
log_msg("QC dims (genes x cells):", nrow(seu_qc), "x", ncol(seu_qc))

out_qc <- file.path(out_dir, "scrna_seurat_qc.rds")
saveRDS(seu_qc, out_qc)
log_msg("Saved:", out_qc)