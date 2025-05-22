#!/usr/bin/env bash 

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

if [ $# -ne 4 ]; then
    echo ""
    echo "Remap cellranger bam with minimap for barcode subsets of stage strain groups"
    echo ""
    echo "Usage: `basename $0` bam_input_file barcodes_of_cells folder_to_output_soupc_matrices number_of_cores"
    echo ""
    exit 1
fi

##Activate necessary modules/environments
#module load ISG/singularity/3.10.0

#source "/software/team222/jr35/miniconda3/bin/activate"
#eval "$(conda shell.bash hook)"
## Provide path to souporcell which has the renamer.py and retag.py
export PATH="$PATH:/software/team222/jr35/hisat2/hisat2"

module load minimap2/2.16=h84994c4_1-c1
module load samtools/1.9
module load freebayes/1.3.6--h346b5cb_1

##Initializes the positional variables

#sw_dir=/software/team222/jr35/souporcell_hard_inst/souporcell
bam=${1}
cell_barcodes=${2}
reference=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/hisat_refs/Pfalciparum.genome
reference_tran=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/hisat_pb63_ref/Pf3D7_DB63_tran
o_file=${3}
ncores=${4}

cd "$(dirname "$o_file"})"

##1. Remapping
##Create fastqs with cell barcodes (CB) and UMIs encoded in readnames
python /software/team222/jr35/souporcell_hard_inst/souporcell/renamer.py --bam $bam --barcodes $cell_barcodes --out ${o_file}/fq.fq

## Remap reads to reference using minimap2
# minimap2 -ax splice -t $ncores -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no $reference ${o_file}_fq.fq > ${o_file}_minimap.sam
# hisat2 -p $ncores -q ${o_file}_fq.fq -x $reference -S ${o_file}_minimap.sam 
hisat2 -p $ncores -q ${o_file}/fq.fq -x $reference_tran -S ${o_file}/tran_hst.sam 

## Retag reads with CBs and UMIs
# python /software/team222/jr35/souporcell_hard_inst/souporcell/retag.py --sam ${o_file}_minimap.sam --out ${o_file}_minitagged.bam
python /software/team222/jr35/souporcell_hard_inst/souporcell/retag.py --sam ${o_file}/tran_hst.sam --out ${o_file}/tran_hsttagged.bam

## Sort and index bam
# samtools sort ${o_file}_minitagged.bam > ${o_file}_minitagged_sorted.bam 
samtools sort ${o_file}/tran_hsttagged.bam > ${o_file}/tran_hsttagged_sorted.bam 

# samtools index ${o_file}_minitagged_sorted.bam
samtools index ${o_file}/tran_hsttagged_sorted.bam

##remove fastq and sam to save space
# rm ${o_file}_fq.fq
# rm ${o_file}_minimap.sam

rm ${o_file}/tran_hst.sam 
