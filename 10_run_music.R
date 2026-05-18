# scripts/10_run_music.R

source("scripts/00_setup_paths.R")

required_pkgs <- c("MuSiC", "SingleCellExperiment", "Seurat", "SeuratObject", "data.table")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(MuSiC)
  library(SingleCellExperiment)
  library(Seurat)
  library(SeuratObject)
  library(data.table)
})

bulk_path <- file.path(DIR_RESULTS, "aligned", "bulk_tpm_common_genes.rds")
seu_path  <- file.path(DIR_RESULTS, "scrna", "scrna_seurat_compartments_singleR.rds")

stopifnot(file.exists(bulk_path))
stopifnot(file.exists(seu_path))

bulk <- readRDS(bulk_path)  # genes x samples (TPM), already aligned to reference genes
seu <- readRDS(seu_path)

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

stopifnot("cellcompartment" %in% colnames(seu@meta.data))
stopifnot("Sample" %in% colnames(seu@meta.data))

# Convert to SingleCellExperiment for MuSiC
sce <- as.SingleCellExperiment(seu)

# MuSiC expects:
# - bulk: ExpressionSet-like matrix with genes as rownames and samples as colnames
# - sc: single-cell expression with celltype and sample identifiers
celltypes <- as.character(colData(sce)$cellcompartment)
samples   <- as.character(colData(sce)$Sample)

stopifnot(length(celltypes) == ncol(sce))
stopifnot(length(samples) == ncol(sce))

log_msg("Running MuSiC")
res <- music_prop(
  bulk.mtx = bulk,
  sc.sce = sce,
  clusters = "cellcompartment",
  samples = "Sample",
  verbose = FALSE
)

out_dir <- file.path(DIR_RESULTS, "deconv", "music")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(res, file.path(out_dir, "music_result.rds"))

# Extract proportions
props <- as.data.frame(res$Est.prop.weighted)
props$sample <- rownames(props)

data.table::fwrite(as.data.table(props), file.path(out_dir, "music_proportions.csv"))

log_msg("Saved: music_result.rds")
log_msg("Saved: music_proportions.csv")
log_msg("Output rows (samples):", nrow(props))
log_msg("Output cols (compartments + sample):", ncol(props))