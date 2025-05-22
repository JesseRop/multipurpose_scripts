#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

##load libraries
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(SingleR)
library(scuttle)
library(purrr)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

REF_OBJ_P <- args[1]
MARKER_OBJ_P <- args[2]
QUERY_OBJ_P <- args[3]
O_SCE <- args[4]
REF_STAGE <- args[5]
O_PRED <- args[6]


# Handle NULL
my_function <- function(x) {
  if (identical(x, "NULL")) {
    x <- NULL
  }
  return(x)
}

# Call the function with "NULL" as a character
# E_CLSTR <- my_function(E_CLSTR)

# print(E_CLSTR)


# Read reference list and marker list
ref_list_obj <- readRDS(REF_OBJ_P)
marker_list_obj <- readRDS(MARKER_OBJ_P)

# Read and convert seurat query object to SCE
query_obj <- readRDS(QUERY_OBJ_P)

DefaultAssay(query_obj) <- "RNA"


## Convert seurat object to SCE
obj_sce <- SingleCellExperiment(assays = list(counts = as(query_obj[["RNA"]]$counts, Class = "dgCMatrix")), colData = query_obj@meta.data)

obj_sce <- logNormCounts(obj_sce) 
print(obj_sce)

# normcounts(obj_sce) <- as(query_obj[["RNA"]]$data, Class = "dgCMatrix")
# rowData(obj_sce)$feature_symbol <- str_replace(rownames(obj_sce), "-", "_")
# reducedDims(obj_sce) <- SimpleList(PCA = query_obj@reductions[[PCA_NM]]@cell.embeddings, UMAP = query_obj@reductions[[UMAP_NM]]@cell.embeddings)


# Save SCE object
saveRDS(obj_sce, O_SCE)

# Performing predictions
pref_pred_hm <- SingleR(
  test = obj_sce, 
  assay.type.test="logcounts", 
  ref = ref_list_obj, 
  assay.type.ref = "logcounts", 
  labels = map(ref_list_obj, ~as.character(.x@colData[[REF_STAGE]])),
  genes = rep(list(marker_list_obj), length(ref_list_obj))
  )

# Performing predictions
saveRDS(pref_pred_hm, O_PRED)

