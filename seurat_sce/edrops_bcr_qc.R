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
OUT_QC_STATS <- args[2]
 

## Read sce object and run emptydrops
sce <- read10xCounts(RAW_P, type = c("HDF5"))

##!!! NOTE The hdf5 matrix makes the running of emptydrops very slow so need to convert to dgCMatrix - https://github.com/MarioniLab/DropletUtils/issues/55
sce_counts <- as(counts(sce), "dgCMatrix")

## check if the input is a SingleCellExperiment object
print(class(sce))

## Get nCountRNA and nFeatureRNA as well as mitochondrial/apicoplast percentage
mito <- which(grepl("MIT", rowData(sce)$ID))
api <- which(grepl("API", rowData(sce)$ID))
stats <- perCellQCMetrics(sce, subsets=list("Mt"=mito, "Api"=api))
stats[, "Barcode"] <- sce$Barcode
saveRDS(stats, OUT_QC_STATS)


