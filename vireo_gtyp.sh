#!/bin/bash

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

if [ $# -ne 6 ]; then
    echo ""
    echo "Run vireo on Plasmodium cellranger output"
    echo ""
    echo "Usage: `basename $0` compute_cores bamfile barcodes_to_consider folder_to_output_results chromosomes number_of_donors(k)"
    echo ""
    exit 1
fi

##Activate necessary modules/environments
source /lustre/scratch126/tol/teams/lawniczak/users/jr35/sware/miniconda3/bin/activate
conda activate /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/envs_py/gtypers_aligners

ncores=${1}
bam=${2}
bcodes=${3}
out_dir=${4}
chroms=${5}
n_donors=${6}

cellsnp-lite -s $bam -b $bcodes -O $out_dir/csnp -p $ncores --minMAF 0.1 --minCOUNT 100 --gzip --chrom $chroms
vireo -c $out_dir/csnp -N $n_donors -o $out_dir/gtypes

