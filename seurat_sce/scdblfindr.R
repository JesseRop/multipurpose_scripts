#!/usr/bin/env

library(Seurat)
library(dplyr) 
library(stringr)
library(tibble)
library(scDblFinder)
library(SingleCellExperiment)
library(purrr)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

FILT <- args[1]
SUPX_FILT <- args[2]
OUT_PATH <- args[3]


## scDblFinder Doublet removal

## Calculate doublets using scDblFinder for raw counts

### Seurat to SCE
dblt_filt <- readRDS(FILT)
DefaultAssay(dblt_filt) <- "RNA"
dblt_sce <- SingleCellExperiment(assays = list(counts = as(object = dblt_filt[["RNA"]]$counts, Class = "dgCMatrix")), colData = dblt_filt@meta.data)

### random - if cluster is only one use random
dblt_rand <- scDblFinder(dblt_sce, returnType = "table", nfeatures = 500)

### cluster - if cluster is only one use random
# dblt_clst <- scDblFinder(dblt_sce, clusters = "seurat_clusters", returnType = "table", nfeatures = 500)
dblt_clst <- if(length(unique(dblt_sce$seurat_clusters)) > 1) { scDblFinder(dblt_sce, clusters = "seurat_clusters", returnType = "table", nfeatures = 500) } else { dblt_rand }

###  stages
# dblt_stg <- scDblFinder(dblt_sce, clusters = "StagePrelimC", returnType = "table", nfeatures = 500)
dblt_stg <- if(length(unique(dblt_sce$StagePrelimC)) > 1) { scDblFinder(dblt_sce, clusters = "StagePrelimC", returnType = "table", nfeatures = 500)} else { dblt_rand}

## Calculate doublets using scDblFinder for soupx adjusted counts - soup correction seems to improve doublet detection based on Elias MCA analysis

## scdblfinder better and sct minimally improves perfomance https://github.com/plger/scDblFinder/issues/67 

### Seurat to SCE
corrected_dblt_filt <- readRDS(SUPX_FILT)
DefaultAssay(corrected_dblt_filt) <- "RNA"
corrected_dblt_sce <- SingleCellExperiment(assays = list(counts = as(object = corrected_dblt_filt[["RNA"]]$counts, Class = "dgCMatrix")), colData = corrected_dblt_filt@meta.data)

### random
corrected_dblt_rand <- scDblFinder(corrected_dblt_sce, returnType = "table", nfeatures = 500)

### cluster
# corrected_dblt_clst <- scDblFinder(corrected_dblt_sce, clusters = "seurat_clusters", returnType = "table", nfeatures = 500)
corrected_dblt_clst <- if(length(unique(corrected_dblt_sce$seurat_clusters)) > 1) { scDblFinder(corrected_dblt_sce, clusters = "seurat_clusters", returnType = "table", nfeatures = 500)} else { dblt_rand}

###  stages
# corrected_dblt_stg <- scDblFinder(corrected_dblt_sce, clusters = "StagePrelimC", returnType = "table", nfeatures = 500)
corrected_dblt_stg <- if(length(unique(corrected_dblt_sce$StagePrelimC)) > 1) { scDblFinder(corrected_dblt_sce, clusters = "StagePrelimC", returnType = "table", nfeatures = 500)} else { dblt_rand}


## Put the different iterations together
scf_outputs <- pmap(list(list(corrected_dblt_clst, corrected_dblt_stg, corrected_dblt_rand, dblt_clst, dblt_stg, dblt_rand), 
                  list("_supx_norm", "_supx_stg", "_supx_rand", "_norm", "_stg", "_rand")), ~ {
      tbl <- data.frame(..1[,c("score", "class", "type")]) %>% 
        filter(type == "real") %>%
        select(score, class) %>%
        mutate(across(class, str_to_sentence)) 
      
      names(tbl) <- paste0("scf_", names(tbl), ..2)
      tbl <- tbl %>% 
        rownames_to_column("bcode")
       }) %>%
      purrr::reduce(., .f = full_join, by = "bcode")


saveRDS(scf_outputs, OUT_PATH)

