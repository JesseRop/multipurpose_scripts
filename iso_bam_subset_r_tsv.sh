#!/bin/bash

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

if [ $# -ne 5 ]; then
    echo ""
    echo "Run souporcell on Plasmodium cellranger output"
    echo ""
    echo "Usage: `basename $0` barcodes_to_subset_file_name folder_to_output_results raw_bam"
    echo ""
    exit 1
fi

##Activate necessary modules/environments
module load samtools/1.14--hb421002_0

##Initializes the positional variables
set_name=${1}
b_codes=${2}
processed_dir=${3}
BAM_FILE=${4}
dedup_codes=${5}

##Subsetting

mkdir -p $processed_dir

##Extract isoseq tag from the code csv file indicating which 10X cell barcode corresponds to which iso-seq tag https://stackoverflow.com/questions/24516141/processing-2-files-with-different-field-separators-using-awk and https://unix.stackexchange.com/questions/331851/grep-line-with-specific-word-from-file-in-a-specific-column
awk 'NR==FNR{a[$1]; next} ($5 in a){print $1}' FS='-' $b_codes  FS='\t' $dedup_codes > $processed_dir/${set_name}_filter.txt

##Faster way to subset the bam files 
samtools view -N $processed_dir/${set_name}_filter.txt -o $processed_dir/${set_name}_filtered.bam $BAM_FILE

##sort
samtools sort $processed_dir/${set_name}_filtered.bam > $processed_dir/${set_name}_filtered_sorted.bam

##index
samtools index $processed_dir/${set_name}_filtered_sorted.bam 

##Remove unnecessary files
rm $processed_dir/${set_name}_filtered.bam
