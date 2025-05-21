#!/bin/bash

##Activate necessary modules/environments
# module load ISG/singularity/3.9.0
module load ISG/singularity/3.11.4
export PATH="$PATH:/lustre/scratch126/tol/teams/lawniczak/users/jr35/sware/hisat2"

##Initializes the positional variables
# Download latest souporcell v2.5 from https://drive.google.com/file/d/1_KIevXI1MvkoXtuiMFv8amlWlC0hTP7p 
# Mount directory has to be upstream of all input and output folders. souporcell script doesn't have to be within tree

# souporcell_mount_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/
souporcell_mount_dir=/lustre/scratch126/tol/teams/lawniczak/
souporcell_exec=/software/team222/jr35/souporcell_singu/souporcell2.5.sif 
bam=${1}
bcodes=${2}
# ref=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/PlasmoDB-52_Pfalciparum3D7_Genome.fasta
xpctd_clusters=${3}
algnr=${4}
ref=${5}
o_dir=${6}
ncores=${7}
tag=${8}



##Running souporcell
singularity exec -B $souporcell_mount_dir $souporcell_exec souporcell_pipeline.py -i $bam -b $bcodes -f $ref -t $ncores  -o $o_dir -k $xpctd_clusters -p 1 --aligner $algnr --umi_tag $tag
