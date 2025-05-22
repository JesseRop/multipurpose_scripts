
#!/usr/bin/env

##load libraries
library(Seurat)
library(BPCells)
library(SoupX)
library(stringr)
library(ggplot2)
library(cowplot)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

INPUT_SEU <- args[1]
INPUT_RAW_MTX <- args[2]
O_FILE_NAME <- args[3]
HVG <- as.numeric(args[4]) ## For function seurat_norm_stand

## source common function
source("/nfs/users/nfs_j/jr35/multipurpose_scripts/seurat_sce/seu_sce_common_fns.R")

## Soupx function
get_SoupX <- function(data, path){
  
  # data <- FindNeighbors(data, k.param = 20)
  # data <- FindClusters(data, resolution = 0.5)
  
  # print(DimPlot(data, group.by <- 'seurat_clusters'))
  
  # toc <- data@assays[["RNA"]]@counts
  # toc <- data[["RNA"]]$counts
  toc <- as(object = data[["RNA"]]$counts, Class = "dgCMatrix")
  
  # tod <- Seurat::Read10X(paste0(path)) switch to h5 for nextflow which needs path to one file - reading many files complicates issues
  tod <- Seurat::Read10X_h5(paste0(path), use.names = F)
  
  # Convert underscore to dashes since createseuratobject does this for the filtered TOC (table of cells) above
  rownames(tod) <- str_replace(rownames(tod), "_", "-")
  print(head(rownames(tod)))
  print(length(rownames(tod)))
  print(head(rownames(toc)))
  print(length(rownames(tod)))
  
  sc <- SoupChannel(tod, toc)
  sc <- setClusters(sc, data$seurat_clusters)
  
  png("rho_plt.png")
  sc <- autoEstCont(sc, doPlot = TRUE)
  dev.off()

  sc_mdata <- sc$metaData
  sc_profile <- sc$soupProfile
  
  # rho_plot <- last_plot()
  
  # saveRDS(plot_object, file = paste0("rho_", O_FILE_NAME)
  # save_plot("rho_plt.png", last_plot())

  adjusted <- adjustCounts(sc)
  data <- CreateSeuratObject(adjusted)
  # data <- seurat_SCT(data)
  data <- seurat_norm_stand(Seurat_object = data, hvg = HVG)
  # return(data)
  return(list(data, sc_mdata, sc_profile))
}


## read seurat
Seu_obj <- readRDS(INPUT_SEU)

##convert to assay3
# Seu_obj[["RNA3"]] <- as(object = Seu_obj[["RNA"]], Class = "Assay")

## Apply function
all_msc_nosct_supx_4_dblt <- get_SoupX(Seu_obj, INPUT_RAW_MTX)

## Apply function
saveRDS(all_msc_nosct_supx_4_dblt, O_FILE_NAME)

