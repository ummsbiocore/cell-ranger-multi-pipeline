#!/usr/bin/env Rscript

library(CellChat)

args<-commandArgs(trailingOnly=TRUE)

seurat_obj1 <- args[[1]]
seurat_obj2 <- args[[2]]
sample_name_1 <- args[[3]]
sample_name_2 <- args[[4]]
script <- args[[5]]

cellchat.NL <- readRDS(seurat_obj1)
cellchat.LS <- readRDS(seurat_obj2)

object.list <- list(sample_name_1 = cellchat.NL, sample_name_2 = cellchat.LS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

rmarkdown::run(script, 
               params = list(file1 = file1, file2 = file2, merged_file = "merged_data.rds"), 
               shiny_args = list(port = 8789, host = '0.0.0.0')
)