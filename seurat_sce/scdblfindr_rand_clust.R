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
OUT_PATH <- args[2]
HVG <- as.numeric(args[3])


## scDblFinder Doublet removal

## Calculate doublets using scDblFinder for raw counts

### Seurat to SCE
dblt_filt <- readRDS(FILT)
DefaultAssay(dblt_filt) <- "RNA"
dblt_sce <- SingleCellExperiment(assays = list(counts = as(object = dblt_filt[["RNA"]]$counts, Class = "dgCMatrix")), colData = dblt_filt@meta.data)

### random - if cluster is only one use random
dblt_rand <- scDblFinder(dblt_sce, returnType = "table", nfeatures = HVG)

### cluster - if cluster is only one use random
# dblt_clst <- scDblFinder(dblt_sce, clusters = "seurat_clusters", returnType = "table", nfeatures = HVG)
dblt_clst <- if(length(unique(dblt_sce$seurat_clusters)) > 1) { scDblFinder(dblt_sce, clusters = "seurat_clusters", returnType = "table", nfeatures = HVG) } else { dblt_rand }

###  stages
# dblt_stg <- scDblFinder(dblt_sce, clusters = "StagePrelimC", returnType = "table", nfeatures = HVG)
# dblt_stg <- if(length(unique(dblt_sce$StagePrelimC)) > 1) { scDblFinder(dblt_sce, clusters = "StagePrelimC", returnType = "table", nfeatures = HVG)} else { dblt_rand}

## Calculate doublets using scDblFinder for soupx adjusted counts - soup correction seems to improve doublet detection based on Elias MCA analysis


## Put the different iterations together
# scf_outputs <- pmap(list(list(dblt_clst, dblt_stg, dblt_rand), 
#                   list("_norm", "_stg", "_rand")), ~ {

scf_outputs <- pmap(list(list(dblt_clst, dblt_rand), 
                  list("_norm", "_rand")), ~ {
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

