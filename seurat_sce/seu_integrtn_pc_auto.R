#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

##load libraries
library(Seurat)
library(findPC)
library(ggplot2)

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
INTEGRTD_PC_MAN <- args[9]
ELBO_NM <- args[10]
FINDPC_MTHD <- as.character(args[11])
MIN_DIST <- as.numeric(args[12])
NNEIGHBS <- as.numeric(args[13])

## Read seurat object
Seurat_object <- readRDS(SEURAT_OBJ_P)

## Integration
Seurat_object <- IntegrateLayers(
  object = Seurat_object, method = MTHD,
  orig.reduction = "pca", new.reduction = NEW_REDCTN,
  verbose = T
)

print(FINDPC_MTHD)
sdev <- prcomp(t(Seurat_object@reductions$pca@cell.embeddings),scale. = T)$sdev[1:30]
sdev <- sort(sdev, decreasing = TRUE)
# npc <- findPC(sdev = sdev, number = 30, figure = T,method = 'all',aggregate = 'voting')
npc <- findPC(sdev = sdev, number = 30, figure = T,method = FINDPC_MTHD)


saveRDS(last_plot(), ELBO_NM)

npc <- ifelse(INTEGRTD_PC_MAN == "yes", INTGRTN_PC, npc)

print(paste0("Manually supplied PCs: ", INTGRTN_PC))
print(paste0("Do we use manually supplied PCs: ", INTEGRTD_PC_MAN))
print(paste0("Number of PCs used: ", npc))


## Clustering
Seurat_object <- FindNeighbors(Seurat_object, reduction = NEW_REDCTN, dims = 1:npc, verbose = FALSE)
Seurat_object <- FindClusters(Seurat_object, resolution = RESOLN, cluster.name = CL_NM, verbose = FALSE)
Seurat_object <- RunUMAP(Seurat_object, reduction = NEW_REDCTN, dims = 1:npc, seed.use = 949, reduction.name = UMAP_REDCTN_NM, n.components = 3L, min.dist = MIN_DIST, n.neighbors = NNEIGHBS)

saveRDS(Seurat_object, OUT_PATH)

