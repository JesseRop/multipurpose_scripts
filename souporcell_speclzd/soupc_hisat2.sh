#!/usr/bin/env bash 

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

# if [ $# -ne 7 ]; then
#     echo ""
#     echo "Run souporcell on Plasmodium cellranger output"
#     echo ""
#     echo "Usage: `basename $0` bam_input_file barcodes_of_cells number_of_expected_clusters_k folder_to_output_soupc_matrices number_of_cores memory  jobID (%J-%I)"
#     echo ""
#     exit 1
# fi

##Activate necessary modules/environments
#module load ISG/singularity/3.10.0

#source "/software/team222/jr35/miniconda3/bin/activate"
#eval "$(conda shell.bash hook)"
#conda activate souporcell
export PATH="$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell"
export PATH="$PATH:/software/team222/jr35/vartrix/vartrix-1.1.22/target/release"

module load minimap2/2.16=h84994c4_1-c1
module load samtools/1.9
module load freebayes/1.3.6--h346b5cb_1
module load hisat2/2.1.0--py36pl5.22.0_0
module load bedtools/2.29.0--hc088bd4_3

##Initializes the positional variables

souporcell_mount_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/
souporcell_exec=/lustre/scratch126/tol/teams/lawniczak/users/jr35/sware/souporcell/souporcell_latest.sif
bam=${1}
cell_barcodes=${2}
# reference=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pfalciparum.genome.fasta
# hisat_reference=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/hisat_transcript_ref/Pfalciparum.genome_tran.fasta
hisat_reference=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/hisat_refs/Pfalciparum.genome.fasta
xpctd_clusters=${3}
o_dir=${4}
ncores=${5}
# mem=${6}
# jobID=${7}


##Running souporcell
souporcell_pipeline.py -i $bam -b $cell_barcodes -f $hisat_reference -t $ncores  -o $o_dir -k $xpctd_clusters -p 1 --aligner HISAT2 --common_variants /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/pf7_data/mali/Pf3D7_v3_g_pass_bisplit_snp_001.vcf


# soupc_job=$(bsub -J $jobID -o  $o_dir/logs/$jobID.o -e $o_dir/logs/$jobID.e -q normal -n $ncores -M $mem -R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" souporcell_pipeline.py -i $bam -b $cell_barcodes -f $reference -t $ncores  -o $o_dir -k $xpctd_clusters -p 1 --aligner HISAT2 --common_variants /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/pf7_data/mali/pf7_mali_snp_pass.vcf)

#soupc_job=$(bsub -J $jobID -o  $o_dir/logs/$jobID.o -e $o_dir/logs/$jobID.e -q normal -n $ncores -M $mem -R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" singularity exec -B $souporcell_mount_dir $souporcell_exec souporcell_pipeline.py -i $bam -b $cell_barcodes -f $reference -t $ncores  -o $o_dir -k $xpctd_clusters -p 1)

soupc_ID=$(echo $soupc_job | sed "s/[^0-9]*//g")

echo "Job $jobID submitted with ID $soupc_ID"


