#!/bin/bash

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

if [ $# -ne 7 ]; then
    echo ""
    echo "JesserRopScrpt Split bam files and generate counts for dexseq analysis"
    echo ""
    echo "Usage: `basename $0` input_bam bam_replicate1 bam_replicate2 bam_sorted_replicate1 bam_sorted_replicate2 replicate1_count_txt replicate2_count_txt"
    echo ""
    exit 1
fi

##!!! NOTE - activate the conda environment on the command line as it doesnt work within the script. Use pip to install latest versions of numpy, pysam and HTSeq in the environment - DO NOT USE CONDA
##Activate necessary modules/environments
# source /software/team222/jr35/miniconda3/bin/activate
# conda activate dexseq_py

module load common-apps/samtools/1.17

##Declare variables pointing to directories
sware_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts

##Prepare gtf for the dexseq mapping DEU analysis - Done only once
#python $sware_dir/dexseq_prepare_annotation.py /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7_ManFromGFF_exons.gtf /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7_ManFromGFF_exons_dexseq.gff

dexseq_gff=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7_ManFromGFF_exons_dexseq.gff

##Initializes the positional variables
in_bam=${1}
bam1=${2}
bam2=${3}
bam_sorted1=${4}
bam_sorted2=${5}
rep1_count=${6}
rep2_count=${7}


##Split bam files into 2 replicates
python $sware_dir/split_bam_test.py -b $in_bam -o1 $bam1 -o2 $bam2 -d . -p 0.5

##sort
samtools sort $bam1 > $bam_sorted1
samtools sort $bam2 > $bam_sorted2

##index
samtools index $bam_sorted1
samtools index $bam_sorted2

## Getting counts for the bam files using dexseq's dexseq_count script
python $sware_dir/dexseq_count.py -f bam $dexseq_gff $bam_sorted1 $rep1_count
python $sware_dir/dexseq_count.py -f bam $dexseq_gff $bam_sorted2 $rep2_count



