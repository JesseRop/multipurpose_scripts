#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

##load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(Matrix)
library(scuttle)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

RAW_P <- args[1]
LIMT <- as.numeric(args[2])
OUT_BCRANK <- args[3]
OUT_BCRANK_KNEE <- args[4]
OUT_BCRANK_INFLCTN <- args[5]
OUT_BCRANK_UNIQ <- args[6]
 

## Read sce object and run emptydrops
sce <- read10xCounts(RAW_P, type = c("HDF5"))

##!!! NOTE The hdf5 matrix makes the running of emptydrops very slow so need to convert to dgCMatrix - https://github.com/MarioniLab/DropletUtils/issues/55
sce_counts <- as(counts(sce), "dgCMatrix")

## check if the input is a SingleCellExperiment object
print(class(sce))

set.seed(100)

## Calculate barcode rank and save 
bcrank <- barcodeRanks(sce_counts, lower=LIMT)
bcrank[, "Barcode"] <- sce$Barcode
saveRDS(bcrank, OUT_BCRANK)

## Also save the estimated inflection and knee points
saveRDS(metadata(bcrank)$knee, OUT_BCRANK_KNEE)
saveRDS(metadata(bcrank)$inflection, OUT_BCRANK_INFLCTN)

## save the barcode rank uniq
bcrank_uniq <- bcrank[!duplicated(bcrank$rank), ]
saveRDS(bcrank_uniq, OUT_BCRANK_UNIQ)


