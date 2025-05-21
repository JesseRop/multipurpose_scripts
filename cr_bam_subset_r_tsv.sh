#!/bin/bash

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

if [ $# -ne 4 ]; then
    echo ""
    echo "Run souporcell on Plasmodium cellranger output"
    echo ""
    echo "Usage: `basename $0` barcodes_to_subset_file_name folder_to_output_results raw_bam"
    echo ""
    exit 1
fi

##Activate necessary modules/environments
module load samtools/1.14--hb421002_0
##older samtools because of disk quota exceed error
#module load samtools/1.9--h91753b0_8 

##Initializes the positional variables
set_name=${1}
b_codes=${2}
processed_dir=${3}
BAM_FILE=${4}

##Subsetting

#subset_job=$(bsub -J $jobID -o $jobID.o -e $jobID.e -q normal -n 8 -M 50000 -R "span[hosts=1] select[mem>50000] rusage[mem=50000]" singularity exec $souporcell_path souporcell_pipeline.py -i $bam -b $barcodes -f $reference -t 8 -o $out -k 2 -p 1 --known_genotypes $vcf) 
mkdir -p $processed_dir

##copy barcodes 
#cp $b_codes $processed_dir
#gunzip $processed_dir/*.gz

##Assigning barcodes to text file
#sed -r 's/$/-1/g' $processed_dir/*.tsv | sed -r 's/^/CB:Z:/g' > $processed_dir/filter.txt
sed -r 's/^/CB:Z:/g' $b_codes > $processed_dir/${set_name}_filter.txt

# Save the header lines
samtools view -H $BAM_FILE > $processed_dir/${set_name}_SAM_header

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $BAM_FILE | LC_ALL=C grep -F -f $processed_dir/${set_name}_filter.txt > $processed_dir/${set_name}_filtered_SAM_body

# Combine header and body
cat $processed_dir/${set_name}_SAM_header $processed_dir/${set_name}_filtered_SAM_body > $processed_dir/${set_name}_filtered.sam

##remove unnecessary files
rm $processed_dir/${set_name}_SAM_header
rm $processed_dir/${set_name}_filtered_SAM_body

# Convert filtered.sam to BAM format
samtools view -b $processed_dir/${set_name}_filtered.sam > $processed_dir/${set_name}_filtered.bam

##remove unnecessary files
rm $processed_dir/${set_name}_filtered.sam

##sort
samtools sort -o $processed_dir/${set_name}_filtered_sorted.bam $processed_dir/${set_name}_filtered.bam
#samtools sort -T $processed_dir -o $processed_dir/${set_name}_filtered_sorted.bam $processed_dir/${set_name}_filtered.bam

##remove unnecessary files
rm $processed_dir/${set_name}_filtered.bam

##index
samtools index $processed_dir/${set_name}_filtered_sorted.bam 

