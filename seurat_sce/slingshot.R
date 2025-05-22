#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

##load libraries
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(stringr)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

SEURAT_OBJ_P <- args[1]
O_FILE_NAME <- args[2]
UMAP_NM <- args[3]
CLSTR <- args[4]
S_CLSTR <- args[5]
# E_CLSTR <- args[6]
E_CLSTR <- eval(parse(text = args[6]))
REDCTN <- args[7]
PCA_NM <- args[8]


# Handle NULL
my_function <- function(x) {
  if (identical(x, "NULL")) {
    x <- NULL
  }
  return(x)
}

# Call the function with "NULL" as a character
E_CLSTR <- my_function(E_CLSTR)

print(S_CLSTR)
print(E_CLSTR)

## Read seurat object
Seurat_object <- readRDS(SEURAT_OBJ_P)

## Convert seurat object to SCE
sce_obj <- SingleCellExperiment(assays = list(counts = as(Seurat_object[["RNA"]]$counts, Class = "dgCMatrix")), colData = Seurat_object@meta.data)
    
# print(sce_obj)
logcounts(sce_obj) <- log2((counts(sce_obj) + 1))
normcounts(sce_obj) <- as(Seurat_object[["RNA"]]$data, Class = "dgCMatrix")
rowData(sce_obj)$feature_symbol <- str_replace(rownames(sce_obj), "-", "_")
reducedDims(sce_obj) <- SimpleList(PCA = Seurat_object@reductions[[PCA_NM]]@cell.embeddings, UMAP = Seurat_object@reductions[[UMAP_NM]]@cell.embeddings)
# reducedDims(sce_obj) <- SimpleList(REDCTN = Seurat_object@reductions[[REDCTN]]@cell.embeddings)

# print(sce_obj)
print(colnames(colData(sce_obj)))

## Convert seurat object to SCE
sce_obj <- slingshot(sce_obj, 
            clusterLabels = CLSTR, 
            reducedDim = REDCTN, 
            start.clus = S_CLSTR,
            end.clus = E_CLSTR,
            shrink =1, 
            stretch = 1,
            allow.break = F)

saveRDS(sce_obj, O_FILE_NAME)
