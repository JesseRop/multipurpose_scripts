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
# module load samtools/1.14--hb421002_0
# module load common-apps/samtools/1.17
module load samtools/1.20--h50ea8bc_0
# module load samtools/1.17--hd87286a_2
# module load samtools/1.19.2--h50ea8bc_1
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

##Run subseting using subset_bam tool copied from sunil but found here
/lustre/scratch126/tol/teams/lawniczak/users/jr35/sware/subset-bam_linux  -b $BAM_FILE -c $b_codes -o $processed_dir/${set_name}_sset.bam

###sort
samtools sort -o $processed_dir/${set_name}_sset_sorted.bam $processed_dir/${set_name}_sset.bam
##samtools sort -T $processed_dir -o $processed_dir/${set_name}_sset_sorted.bam $processed_dir/${set_name}_sset.bam
#
###remove unnecessary files
rm $processed_dir/${set_name}_sset.bam

##index
samtools index $processed_dir/${set_name}_sset_sorted.bam 

