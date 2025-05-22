##Commonly used seurat R functions

## Seurat normal standardization
seurat_norm_stand = function(Seurat_object, hvg){
  Seurat_object = NormalizeData(Seurat_object, verbose = F)
  Seurat_object = FindVariableFeatures(Seurat_object, nfeatures = hvg, verbose = F)
  Seurat_object = ScaleData(Seurat_object, verbose = F)
  Seurat_object = RunPCA(Seurat_object, verbose = F)
  return(Seurat_object)
}

## Seurat SCT standardization
seurat_SCT = function(Seurat_object) {
  
  Seurat_object %>%
    SCTransform(., vst.flavor = "v2", verbose = FALSE, variable.features.n) %>%
    RunPCA(npcs = 30, verbose = F)
  
}

#### Average expression of gene sets 
makeSetScore <- function(seur, gene.set, assay='RNA') {
  # Get mean expression of genes of interest per cell
  mean.exp <- colMeans(x = as.matrix(GetAssayData(seur, slot = 'data', assay = assay))[gene.set, ], na.rm = TRUE)
  
  # Add mean expression values in 'object@meta.data$gene.set.score'
  if (all(names(x = mean.exp) == rownames(x = seur@meta.data))) {
    return(mean.exp)
  }
}

