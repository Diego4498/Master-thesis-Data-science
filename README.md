# Master Thesis — Immune Deconvolution Benchmark

This repository contains the code and analysis workflow for a master’s thesis project focused on benchmarking immune deconvolution methods in metastatic prostate cancer.

The project evaluates how different deconvolution frameworks behave when applied to **bulk RNA-seq** data using an external **single-cell RNA-seq (scRNA-seq)** reference in a **cross-cohort, non-patient-matched setting**. The main methods explored are **MuSiC** and **EPIC**, with the workflow structured to support comparative extension to additional approaches. The study focuses on how **preprocessing, annotation, gene harmonization, and reference construction** influence inferred cell-type proportions. :contentReference[oaicite:0]{index=0} :contentReference[oaicite:1]{index=1}

## Project Background

Bulk RNA-seq is widely used to profile tumor transcriptomes, but it captures averaged expression across heterogeneous cell populations. This makes direct characterization of the tumor microenvironment difficult. Deconvolution methods aim to address this by estimating cell-type composition from bulk expression profiles using reference signatures derived from purified or single-cell data. :contentReference[oaicite:2]{index=2}

In this project, bulk RNA-seq data from **metastatic prostate adenocarcinoma** are combined with a scRNA-seq reference derived from **metastatic prostate lesions and PBMC-related data**. Because the bulk and single-cell datasets are not patient-matched, the repository reflects a realistic translational scenario in which public single-cell resources are used as external references for independent tumor cohorts. :contentReference[oaicite:3]{index=3} :contentReference[oaicite:4]{index=4}

## Objectives

- Build a reproducible immune deconvolution workflow in R/R Markdown
- Preprocess and harmonize bulk RNA-seq and scRNA-seq datasets
- Annotate single-cell data and construct reference expression profiles
- Apply deconvolution methods to bulk tumor samples
- Compare inferred immune composition across methods
- Assess the impact of preprocessing and reference design choices on robustness and reproducibility :contentReference[oaicite:5]{index=5} :contentReference[oaicite:6]{index=6}

## Datasets

### Bulk RNA-seq
- **GSE297742**
- Biopsy-derived bulk RNA-seq expression data from metastatic prostate adenocarcinoma samples

### Single-cell RNA-seq
- **GSE297652**
- scRNA-seq data from metastatic prostate cancer samples used for reference construction

## Current Workflow

The analysis is implemented as a modular R Markdown pipeline. The current repository includes the following stages:

### 1. Path setup and project structure
Initial project root detection, output folder creation, and logging utilities.

### 2. Bulk RNA-seq loading and preparation
- Load bulk expression matrix
- Load clinical metadata
- Remove duplicated genes
- Match expression samples with metadata
- Save cleaned bulk matrix and metadata

### 3. Bulk label definition
- Inspect metadata for metastasis-related labels
- Build an M0/M1 subset based on `Clinical_stage_M`

### 4. scRNA-seq loading and preparation
- Load expression matrix, cell metadata, and gene annotations
- Inspect metadata columns and candidate identifiers

### 5. Seurat object construction and QC
- Read 10x matrices
- Merge samples
- Compute mitochondrial content
- Filter low-quality cells
- Save raw and QC-filtered Seurat objects

### 6. scRNA-seq processing and clustering
- Normalization
- Variable feature selection
- Scaling
- PCA
- Nearest-neighbor graph construction
- Clustering
- Marker detection

### 7. Single-cell annotation
- Cell-type annotation with **SingleR**
- Cluster-level label assignment
- Export cell- and cluster-level labels

### 8. Reference construction
- Collapse SingleR labels into broader compartments
- Build average reference expression profiles per compartment

### 9. Gene harmonization
- Intersect genes between bulk RNA-seq and scRNA-seq reference
- Save aligned matrices for downstream deconvolution

### 10. Deconvolution
Currently implemented:
- **MuSiC**
- **EPIC**

### 11. Validation-oriented extensions
The workflow also includes additional PBMC-oriented steps for annotation, reference construction, and pseudobulk generation, used as technical support for validation and method checking. :contentReference[oaicite:7]{index=7}

## Repository Structure

```bash
.
├── data/
│   ├── bulk/
│   └── scrna/
│       └── 10x/
├── scripts/
├── results/
│   ├── bulk/
│   ├── scrna/
│   ├── references/
│   ├── aligned/
│   ├── deconv/
│   ├── pbmc/
│   └── logs/
├── figures/
├── logs/
└── Master thesis Deconvolution.Rmd
