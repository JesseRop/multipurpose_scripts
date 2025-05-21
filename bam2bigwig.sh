#!/bin/bash

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

if [ $# -ne 6 ]; then
    echo ""
    echo "JesserRopScrpt Convert bam to bigwig"
    echo ""
    echo "Usage: `basename $0` genome_size input_bam plus_strand_bedgraph_output_location plus_strand_bigwig_output_location minus_strand_bedgraph_output_location minus_strand_bigwig_output_location"
    echo ""
    exit 1
fi

##Activate necessary modules/environments
source /lustre/scratch126/tol/teams/lawniczak/users/jr35/sware/miniconda3/bin/activate
conda activate /lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/envs_py/gtypers_aligners

module load bedtools/2.29.0--hc088bd4_3


##Initializes the positional variables
genome_size=${1}
in_bam=${2}
plus_bedgraph=${3}
plus_bgwg=${4}
minus_bedgraph=${5}
minus_bgwg=${6}


##Get genome size - already done hence commenting out but keeping code for documentation
##cut -f1,2 /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7_Genome.fasta.fai > /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7_Genome.size

##convert positive strands of bam to bedgraph then to bigwig
bedtools genomecov -split -ibam $in_bam -bg -strand '+' | LC_COLLATE=C sort -k1,1 -k2,2n > $plus_bedgraph
bedGraphToBigWig $plus_bedgraph $genome_size $plus_bgwg

##convert minus strands of bam to bedgraph then to bigwig
bedtools genomecov -split -ibam $in_bam -bg -strand '-' | LC_COLLATE=C sort -k1,1 -k2,2n > $minus_bedgraph
bedGraphToBigWig $minus_bedgraph $genome_size $minus_bgwg


