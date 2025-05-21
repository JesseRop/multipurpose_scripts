#!/bin/bash

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

# Load the necessary modules
module load python-3.9.18/perl-5.38.0
module load cuda-12.1.1

# activate python env
source /software/team222/jr35/cellbender_gpu/cellbender_gpu_env/bin/activate 

##Initializes the positional variables
INPUT_RAW_MTX=${1}
O_DIR=${2}
O_FILE=${3}
ncores=${4}
mem=${5}
jobID=${6}


##Cellbender run
## GPU job launch instructions here - https://ssg-confluence.internal.sanger.ac.uk/display/FARM/How+to+submit+a+job+to+the+GPU+queues+in+LSF
cbendr_job=$(bsub -J $jobID -o $O_DIR/logs/$jobID.o -e $O_DIR/logs/$jobID.e -n $ncores -M $mem -R "span[hosts=1] select[mem>$mem] rusage[mem=$mem]"  -q gpu-normal -gpu 'num=1:j_exclusive=yes' cellbender remove-background --cuda --input $INPUT_RAW_MTX --output $O_FILE)

## Asking for 2 gpus
# -gpu 'num=2:j_exclusive=yes:glink=yes'

cbendr_ID=$(echo $cbendr_job | sed "s/[^0-9]*//g")

echo "Job $jobID submitted with ID $cbendr_ID"
