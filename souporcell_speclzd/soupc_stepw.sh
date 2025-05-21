#!/usr/bin/env bash 

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

if [ $# -ne 5 ]; then
    echo ""
    echo "Run souporcell on Plasmodium cellranger output"
    echo ""
    echo "Usage: `basename $0` bam_input_file barcodes_of_cells number_of_expected_clusters_k folder_to_output_soupc_matrices number_of_cores memory  jobID (%J-%I)"
    echo ""
    exit 1
fi

##Activate necessary modules/environments
#module load ISG/singularity/3.10.0

#source "/software/team222/jr35/miniconda3/bin/activate"
#eval "$(conda shell.bash hook)"
#conda activate souporcell
export PATH="$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell"
export PATH="$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell/souporcell/target/release"
export PATH="$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell/troublet/target/release"
export PATH="$PATH:/software/team222/jr35/vartrix/vartrix-1.1.22/target/release"

module load minimap2/2.16=h84994c4_1-c1
module load samtools/1.14--hb421002_0
module load freebayes/1.3.6--h346b5cb_1

##Initializes the positional variables

#sw_dir=/software/team222/jr35/souporcell_hard_inst/souporcell
souporcell_mount_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/
souporcell_exec=/lustre/scratch126/tol/teams/lawniczak/users/jr35/sware/souporcell/souporcell_latest.sif
bam=${1}
cell_barcodes=${2}
reference=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pfalciparum.genome.fasta
xpctd_clusters=${3}
o_dir=${4}
ncores=${5}
#mem=${6}
#jobID=${7}

cd $o_dir

##1. Remapping
##Create fastqs with cell barcodes (CB) and UMIs encoded in readnames
python renamer.py --bam $bam --barcodes $cell_barcodes --out fq.fq

## Remap reads to reference using minimap2
minimap2 -ax splice -t $ncores -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no $reference $o_dir/fq.fq > $o_dir/minimap.sam

## Retag reads with CBs and UMIs
python retag.py --sam $o_dir/minimap.sam --out $o_dir/minitagged.bam

## Sort and index bam
samtools sort $o_dir/minitagged.bam > $o_dir/minitagged_sorted.bam


##remove fastq and sam to save space
rm $o_dir/fq.fq
rm $o_dir/minimap.sam

## Call variants
#freebayes -f $reference -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --max-coverage 100000 $o_dir/minitagged_sorted.bam > $o_dir/var.vcf
freebayes -f $reference -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 $o_dir/minitagged_sorted.bam > $o_dir/var.vcf

## Assign variants to cells
vartrix --umi --mapq 30 -b $o_dir/minitagged_sorted.bam -c $cell_barcodes --scoring-method coverage --threads $ncores --ref-matrix $o_dir/ref.mtx --out-matrix $o_dir/alt.mtx -v $o_dir/var.vcf --fasta $reference

## Souporcell on vcf
souporcell -a $o_dir/alt.mtx -r $o_dir/ref.mtx -b $cell_barcodes -k $xpctd_clusters -t $ncores > $o_dir/clusters_tmp.tsv

## Doublet detction
troublet -a $o_dir/alt.mtx -r $o_dir/ref.mtx --clusters $o_dir/clusters_tmp.tsv > $o_dir/clusters.tsv

## Genotype and ambient RNA
consensus.py -c $o_dir/clusters.tsv -a $o_dir/alt.mtx -r $o_dir/ref.mtx --soup_out $o_dir/soup.txt -v $o_dir/var.vcf -p 1 --vcf_out $o_dir/cluster_genotypes.vcf --output_dir $o_dir/


#soupc_job=$(bsub -J $jobID -o  $o_dir/logs/$jobID.o -e $o_dir/logs/$jobID.e -q normal -n $ncores -M $mem -R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" souporcell_pipeline.py -i $bam -b $cell_barcodes -f $reference -t $ncores  -o $o_dir -k $xpctd_clusters -p 1 --ignore True)

#soupc_job=$(bsub -J $jobID -o  $o_dir/logs/$jobID.o -e $o_dir/logs/$jobID.e -q normal -n $ncores -M $mem -R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" singularity exec -B $souporcell_mount_dir $souporcell_exec souporcell_pipeline.py -i $bam -b $cell_barcodes -f $reference -t $ncores  -o $o_dir -k $xpctd_clusters -p 1)

#soupc_ID=$(echo $soupc_job | sed "s/[^0-9]*//g")

#echo "Job $jobID submitted with ID $soupc_ID"


