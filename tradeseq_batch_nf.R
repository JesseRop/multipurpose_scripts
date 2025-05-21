#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

##load libraries
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(BiocParallel)
library(purrr)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

SCE_OBJ <- args[1]
PTIME_OBJ <- args[2]
CW_OBJ <- args[3]
U_DON_OBJ <- args[4]
CONDTN <- args[5]
# COMPS <- args[5]
KNOTS <- as.numeric(args[6])
GNS_CUT <- as.numeric(args[7])
TS_OUT <- args[8]
COND_DE_OUT <- args[9]


SCE <- readRDS(SCE_OBJ)
PTIME <-readRDS(PTIME_OBJ)
CW <- readRDS(CW_OBJ)
U_DON <- readRDS(U_DON_OBJ)

print(CONDTN)
# print(COMPS)
print(KNOTS)

## Subset the relevant cells for each run
SCE = SCE[,!is.na(SCE@colData[,CONDTN])]
PTIME = PTIME[!is.na(SCE@colData[,CONDTN])]
CW = CW[!is.na(SCE@colData[,CONDTN])]

## Tradeseq script
ts_sce_ub_de  <- fitGAM(counts = counts(SCE), 
                    pseudotime = PTIME, 
                    cellWeights = CW, 
                    nknots =KNOTS,
                    genes = rownames(SCE)[rowSums(counts(SCE) != 0) > GNS_CUT],
                    conditions = factor(SummarizedExperiment::colData(SCE)[, CONDTN]), 
                    verbose = T, 
                    parallel = T,
                    U = U_DON
                    )
            

## Save tradeseq object
saveRDS(ts_sce_ub_de, TS_OUT)

##Running condition test
condRes_ub_de <- conditionTest(ts_sce_ub_de, pairwise = T)
        
## Save conidition test object    
saveRDS(condRes_ub_de, COND_DE_OUT)
