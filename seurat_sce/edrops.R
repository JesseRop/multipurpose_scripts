#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

##load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(Matrix)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

RAW_P <- args[1]
LIMT <- as.numeric(args[2])
RETAIN <- args[3]
FDR_L <- as.numeric(args[4])
OUT_NM <- args[5]
OUT_NM_SLIM <- args[6]
 

# Handle NULL
my_function <- function(x) {
  if (identical(x, "NULL")) {
    x <- NULL
  } else if (identical(x, "Inf")) {
    x <- Inf
  }
  return(x)
}

# Call the function with "NULL" as a character
RETAIN <- my_function(RETAIN)

## Read sce object and run emptydrops
sce <- read10xCounts(RAW_P, type = c("HDF5"))

## check if the input is a SingleCellExperiment object
print(class(sce))

##!!! NOTE The hdf5 matrix makes the running of emptydrops very slow so need to convert to dgCMatrix - https://github.com/MarioniLab/DropletUtils/issues/55
##!!! reading the folder with the matrices does not have this issue

sce_counts <- as(counts(sce), "dgCMatrix")
print(sce_counts[1:5, 1:5])

set.seed(100)
# e.out <- emptyDrops(sce_counts[1:500, 1:50000], lower=LIMT, test.ambient=TRUE, retain = RETAIN)

e.out <- emptyDrops(sce_counts, lower=LIMT, test.ambient=TRUE, retain = RETAIN)
print(class(e.out))
## add new column with whether cell or empty based on FDR
e.out["is_cell_ed"] <- ifelse(e.out$FDR <= FDR_L, "CellEd", "Non-Cell")

## save full emptydrops output and that with only identified cells
saveRDS(e.out, OUT_NM)
saveRDS(e.out[e.out[,"is_cell_ed"] == "CellEd" & !is.na(e.out[,"is_cell_ed"]), ], OUT_NM_SLIM)





