#!/bin/bash

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

# if [ $# -ne 7 ]; then
#     echo ""
#     echo "Run velocyto on cellranger original"
#     echo ""
#     echo "Usage: `basename $0` barcodes_of_filtered_cells folder_to_output_vcyto_matrices bam_input_file gtf number_of_cores memory  jobID (%J-%I)"
#     echo ""
#     exit 1
# fi

##Activate necessary modules/environments
#source /software/team222/jr35/miniconda3/etc/profile.d/conda.sh
# conda activate vcyto_tools
# module load common-apps/samtools/1.17
module load cellgen/samtools/1.17

##Initializes the positional variables

#index_name=${1}
bam_in=${1}
gtf=${2}
filt_cell_bc=${3}
count_o_dir=${4}
ncores=${5}
mem=${6}
jobID=${7}

##Make log folder
#mkdir $count_o_dir/logs

##Running velocyto
vcyto_job=$(bsub -J $jobID -o $count_o_dir/logs/$jobID.o -e $count_o_dir/logs/$jobID.e -q normal -n $ncores -M $mem -R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" velocyto run -b $filt_cell_bc -o $count_o_dir $bam_in $gtf)

vcyto_ID=$(echo $vcyto_job | sed "s/[^0-9]*//g")

echo "Job $jobID submitted with ID $vcyto_ID"


