#!/usr/bin/env

## Assigning stages

##load libraries
library(Seurat)
library(SingleCellExperiment)
library(scmap)
library(dplyr)
library(stringr)
library(tibble)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

SEU_OBJ <- args[1]
SCE_REF <- args[2]
W <- args[3]
COS_SIM <- args[4]
O_FILE <- args[5]
STG_LABL <- args[6]

##Read MCA V3 + field blood stage reference dataset
v3labmscs.ref.sce <- readRDS(SCE_REF)

## Add the column to use as stage labels
v3labmscs.ref.sce@colData[,"stageHL_ref"] <- v3labmscs.ref.sce@colData[,STG_LABL]


##scmap feature selection
v3labmscs.ref.sce <- selectFeatures(v3labmscs.ref.sce, suppress_plot = T)

##Cell indexing
set.seed(155)
v3labmscs.ref.sce <- indexCell(v3labmscs.ref.sce)


## Read in seurat object
all_msc <- readRDS(SEU_OBJ)

### Preprocessing MSC query datasets

## Convert seurat object into single cell experiment filling in the slots appropriately

##This command gives error below when running the entire notebook or this whole chunk so run the command independently - coulb be and issue with rmd
##Error in x$.self$finalize() : attempt to apply non-function - need blank line after

msc_sce <- SingleCellExperiment(assays = list(counts = as(object = all_msc[["RNA"]]$counts, Class = "dgCMatrix")), colData = all_msc@meta.data)

logcounts(msc_sce) <- log2((SingleCellExperiment::counts(msc_sce) + 1))

rowData(msc_sce)$feature_symbol <- rownames(msc_sce)

msc_sce <- msc_sce[!duplicated(rownames(msc_sce)), ]

### Stage label assignment (projection)
# transferring labels from  reference datasets

##Seurat standard transformation

msc_jc <-  scmapCell(msc_sce, list(
  Stage_MSC = metadata(v3labmscs.ref.sce)$scmap_cell_index
)
)

##Projection
# Assigning stage labels only when the top 3 most closest reference cells to a query cell all have a similar stage label and a cosine similarity greater than 0.4
msc_lbls <- scmapCell2Cluster(msc_jc, list(
  # as.character(SummarizedExperiment::colData(msc_combi.sce)$stageHL_ref)
  as.character(SummarizedExperiment::colData(v3labmscs.ref.sce)$stageHL_ref)
), w = as.numeric(W), threshold = as.numeric(COS_SIM))


# Obtaining the reference cell that is closest to each query cell in neighbourhood space regardless of cosine similarity. For those where the top 10 most closest cells constitute several stages, we also obtain the second most similar stage.
## Get the top V3 cell and top 3 similarity scores that is most similar to each MSC14 cell 

sc_topv3_cells_sims <- msc_jc$Stage_MSC$similarities %>% 
  .[1:3,] %>% 
  t() %>% 
  data.frame() %>% 
  rename_all(., ~str_replace(., 'X', 'sim_lab')) %>% 
  rownames_to_column('bcode') %>%
  left_join(., msc_jc$Stage_MSC$cells %>% 
              .[1:3,] %>% t() %>% 
              data.frame() %>% 
              rename_all(., ~str_replace(., 'X', 'msc_ref_cells')) %>% 
              rownames_to_column('bcode'),
            by = 'bcode') %>% 
  # mutate(across(starts_with('msc_ref_cells'), ~colnames(msc_combi.sce)[.]))
  mutate(across(starts_with('msc_ref_cells'), ~colnames(v3labmscs.ref.sce)[.]))



# Adding annotations to MSC query datasets

## Adding annotations to query cells
stg_anots <- v3labmscs.ref.sce@colData[,c('stageHL_ref'), drop=F] %>% data.frame()

## scmap merge the top labels and similarity scores
top_stg_sim <- sc_topv3_cells_sims %>% 
  left_join( stg_anots %>% rownames_to_column('msc_ref_cells1') %>% rename('stage_lab1' = stageHL_ref), by = 'msc_ref_cells1') %>%
  left_join(stg_anots %>% rownames_to_column('msc_ref_cells2') %>% rename('stage_lab2' = stageHL_ref), by = 'msc_ref_cells2') %>%
  left_join(stg_anots %>% rownames_to_column('msc_ref_cells3') %>% rename('stage_lab3' = stageHL_ref), by = 'msc_ref_cells3') %>% 
  column_to_rownames('bcode')


# Put the annotations together and write out
msc_anots <- msc_lbls$scmap_cluster_labs %>% 
  data.frame()  %>% 
  mutate(cell_bc = colnames(msc_jc$Stage_MSC$cells))%>%
  full_join(., top_stg_sim %>% rownames_to_column("cell_bc"), by = "cell_bc")


# saveRDS(msc_anots, paste0(O_FILE,"/msc_anots.RDS"))
saveRDS(msc_anots, O_FILE)



