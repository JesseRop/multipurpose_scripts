#!/usr/bin/env bash 

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

if [ $# -ne 6 ]; then
    echo ""
    echo "Run souporcell on Plasmodium cellranger output"
    echo ""
    echo "Usage: `basename $0` bam_input_file barcodes_of_cells number_of_expected_clusters_k folder_to_output_soupc_matrices number_of_cores memory  jobID (%J-%I)"
    echo ""
    exit 1
fi

##Activate necessary modules/environments

#source "/software/team222/jr35/miniconda3/bin/activate"
#eval "$(conda shell.bash hook)"
#conda activate souporcell
export PATH="$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell"
export PATH="$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell/souporcell/target/release"
export PATH="$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell/troublet/target/release"
export PATH="$PATH:/software/team222/jr35/vartrix/vartrix-1.1.22/target/release"

##Initializes the positional variables

#sw_dir=/software/team222/jr35/souporcell_hard_inst/souporcell
alt_mtx=${1}
ref_mtx=${2}
cell_barcodes=${3}
xpctd_clusters=${4}
ncores=${5}
var_vcf=${6}

## Souporcell on vcf
souporcell -a $alt_mtx -r $ref_mtx -b $cell_barcodes -k $xpctd_clusters -t $ncores > clusters_tmp.tsv 2> clusters.err

## Doublet detction
troublet -a $alt_mtx -r $ref_mtx --clusters clusters_tmp.tsv > clusters.tsv 2> doublets.err

## Genotype and ambient RNA - tips on saving output https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file
consensus.py -c clusters.tsv -a $alt_mtx -r $ref_mtx --soup_out soup.txt -v $var_vcf -p 1 --vcf_out cluster_genotypes.vcf --output_dir . 2> consensus.err

## Put the log likelihood for knee plot in separate file for easy retrieval in R
grep -H 'best total log probability' clusters.err > clusters_log_likelihoods_${xpctd_clusters}.txt
