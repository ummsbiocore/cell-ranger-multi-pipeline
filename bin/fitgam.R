#!/usr/bin/env Rscript

library(tradeSeq)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(BiocParallel)
library(Seurat)
library(fastglm)
library(SingleCellExperiment)
library(slingshot)

args<-commandArgs(trailingOnly=TRUE)

obj_path <- args[[1]]
sce_path <- args[[2]]
num_genes <- args[[4]]

so <- readRDS(obj_path)
sce <- readRDS(sce_path)
# ---------------------------
# starting variables (speed + outputs)
# ---------------------------
n_knots      <- 6
top_n        <- 250
padj_method  <- "BH"
min_wt       <- 0.0
heat_scale   <- "log1p"

# ---- fitGAM speed controls (from fitGAM docs) ----
use_parallel <- TRUE
n_workers    <- as.numeric(args[[3]])  # set to your CPU core count
genes_subset <- NULL               # e.g. VariableFeatures(so) for faster run, or NULL = all genes
gam_nthreads <- 1                  # mgcv threads per worker
gam_maxit    <- 50                 # fewer iterations = faster (keep reasonable)

# Parallel backend (Windows-safe)
if (.Platform$OS.type == "windows") {
  BPPARAM <- BiocParallel::SnowParam(workers = n_workers, type = "SOCK")
} else {
  # on Linux/macOS: use MulticoreParam for lower overhead
  # but only if fork-safe for your code/packages
  # BPPARAM <- BiocParallel::SnowParam(workers = n_workers, type = "SOCK")
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_workers)
}
register(BPPARAM)

genes_subset <- VariableFeatures(so)
if (num_genes != ""){
  genes_subset <- genes_subset[1:as.numeric(num_genes)]
}

print(length(genes_subset))
# asdfasdf

# Gene subset index (fitGAM expects indices if genes= used)
if (is.null(genes_subset)) {
  genes_use_idx <- seq_len(nrow(assays(sce)$counts))
} else {
  genes_use_idx <- which(rownames(assays(sce)$counts) %in% genes_subset)
  if (length(genes_use_idx) == 0) {
    stop("genes_subset did not match any genes in sce.")
  }
}

# mgcv control to speed fitting
gam_ctrl <- mgcv::gam.control(
  nthreads = gam_nthreads,
  maxit    = gam_maxit
)

# variable_genes <- VariableFeatures(so)
# genes_use_idx <- variable_genes

# ---------------------------
# Fit NB-GAMs (FASTER: fastglm + parallel BPPARAM + optional gene subset)
# ---------------------------
sce <- fitGAM(
  counts   = assays(sce)$counts,
  sds      = SlingshotDataSet(sce),
  genes    = genes_use_idx,
  nknots   = n_knots,
  verbose  = TRUE,
  parallel = use_parallel,
  BPPARAM  = BPPARAM,
  control  = gam_ctrl,
  sce      = TRUE,
  family   = "nb"
)

saveRDS(sce, "sce_fitgam.rds")