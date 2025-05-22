#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

##load libraries
library(Seurat)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

SEURAT_OBJ_P <- args[1]
MTHD <- args[2]
OUT_PATH <- args[3]
# UMAP_NM <- args[4]
RESOLN <- as.numeric(args[4])
CL_NM <- args[5]
UMAP_REDCTN_NM <- args[6]
NEW_REDCTN <- args[7]
INTGRTN_PC <- as.numeric(args[8])
K_WT <- as.numeric(args[9])
NORM_MTHD <- args[10]

print(NORM_MTHD)
print(K_WT)
## Read seurat object
Seurat_object <- readRDS(SEURAT_OBJ_P)

## Integration
Seurat_object <- IntegrateLayers(
  object = Seurat_object, method = MTHD, normalization.method = NORM_MTHD,
  orig.reduction = "pca", new.reduction = NEW_REDCTN,
  verbose = T,k.weight = K_WT
)

## Clustering
Seurat_object <- FindNeighbors(Seurat_object, reduction = NEW_REDCTN, dims = 1:INTGRTN_PC, verbose = FALSE)
Seurat_object <- FindClusters(Seurat_object, resolution = RESOLN, cluster.name = CL_NM, verbose = FALSE)
Seurat_object <- RunUMAP(Seurat_object, reduction = NEW_REDCTN, dims = 1:INTGRTN_PC, seed.use = 949, reduction.name = UMAP_REDCTN_NM, n.components = 3L)

saveRDS(Seurat_object, OUT_PATH)
