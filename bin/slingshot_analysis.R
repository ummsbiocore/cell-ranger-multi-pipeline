#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

args<-commandArgs(trailingOnly=TRUE)

obj_path <- args[[1]]
reduc <- args[[2]]
cluster_var <- args[[3]]

so <- readRDS(obj_path)
rd <- Embeddings(so, reduction = reduc)[, 1:2, drop = FALSE]
cl <- factor(so[[cluster_var]][, 1])

sce <- as.SingleCellExperiment(so)

reducedDims(sce) <- SimpleList()
reducedDims(sce)[[reduc]] <- rd
colData(sce)$cluster <- cl

set.seed(777)
sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim = reduc
)

saveRDS(sce, "sce.rds")
