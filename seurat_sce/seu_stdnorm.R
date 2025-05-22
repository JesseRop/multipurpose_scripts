#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

##load libraries
library(Seurat)
library(BPCells)
library(findPC)
library(ggplot2)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

SEURAT_OBJ_P <- args[1]
HVG <- as.numeric(args[2])
OUT_PATH <- args[3]
UMAP_NM <- args[4]
RESOLN <- as.numeric(args[5])
CL_NM <- args[6]
UMAP_REDCTN_NM <- args[7]
REDCTN <- args[8]
INTEGRTD_PC30 <- args[9]
ELBO_NM <- args[10]
FINDPC_MTHD <- as.character(args[11])
PCA_REDCTN_NM <- args[12]
MIN_DIST <- as.numeric(args[13])

# Handle NULL
my_function <- function(x) {
  if (identical(x, "NULL")) {
    x <- NULL
  }
  return(x)
}

# Call the function with "NULL" as a character
CL_NM <- my_function(CL_NM)

## Read seurat object
Seurat_object <- readRDS(SEURAT_OBJ_P)

## Runs seurat standard normalisation
seurat_norm_stand = function(Seurat_object, hvg = HVG){
  Seurat_object = NormalizeData(Seurat_object, verbose = F)
  Seurat_object = FindVariableFeatures(Seurat_object, nfeatures = hvg, verbose = F)
  Seurat_object = ScaleData(Seurat_object, verbose = F)
  Seurat_object = RunPCA(Seurat_object, verbose = F, reduction.name = PCA_REDCTN_NM)
  return(Seurat_object)
}

Seurat_object <- seurat_norm_stand(Seurat_object, hvg = HVG)

sdev <- prcomp(t(Seurat_object@reductions$pca@cell.embeddings),scale. = T)$sdev[1:30]
sdev <- sort(sdev, decreasing = TRUE)
# npc <- findPC(sdev = sdev, number = 30, figure = T,method = 'all',aggregate = 'voting')
npc <- findPC(sdev = sdev, number = 30, figure = T, method = FINDPC_MTHD)

saveRDS(last_plot(), ELBO_NM)

npc <- ifelse(as.numeric(npc) < 5, 5, as.numeric(npc))
npc <- ifelse(INTEGRTD_PC30 == "yes", 30, npc)

print(npc)


Seurat_object <- RunUMAP(Seurat_object, reduction = REDCTN, dims = 1:npc, n.components = 3L, verbose=F, seed.use = 949, reduction.key = UMAP_NM, reduction.name = UMAP_REDCTN_NM, min.dist = MIN_DIST)
Seurat_object <- FindNeighbors(Seurat_object, reduction = REDCTN, dims = 1:npc, verbose = FALSE)
Seurat_object <- FindClusters(Seurat_object, resolution = RESOLN, verbose = FALSE, cluster.name = CL_NM)

saveRDS(Seurat_object, OUT_PATH)
