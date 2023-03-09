#!/bin/bash
#SBATCH -J P005_FIG11 # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one Node
#SBATCH -p shared # Partition
#SBATCH --mem 2000 # Memory request (4Gb)
#SBATCH -t 0-10:00 # Maximum execution time (D-HH:MM)
#SBATCH -o "/n/home02/ndeporzio/projects/P005_FuzzyCdmBias/results/FIG11_%a.out" # Standard output
#SBATCH -e "/n/home02/ndeporzio/projects/P005_FuzzyCdmBias/results/FIG11_%a.err" # Standard error
#SBATCH --account=dvorkin_lab
#SBATCH --mail-type=ALL        ## BEGIN, END, FAIL, ALL                         
#SBATCH --mail-user=nicholasdeporzio@g.harvard.edu 

## SETUP STORAGE ##                                                             
PROJECT_DIR="n/home02/ndeporzio/projects/P005_FuzzyCdmBias/"                          
STORAGE_DIR="/n/home02/ndeporzio/projects/P005_FuzzyCdmBias/results/"           
RELAXIFAST_DIR="${SCRATCH}/dvorkin_lab/ndeporzio/RelAxiFast_${SLURM_ARRAY_TASK_ID}/"        

mkdir -p $STORAGE_DIR

if [ -d "$RELAXIFAST_DIR" ]; then
  echo "RelAxiFast already installed. Deleting..."
  rm -rf $RELAXIFAST_DIR
fi

echo "Installing RelAxiFast..."
cd "${SCRATCH}/dvorkin_lab/ndeporzio"
git clone git@github.com:ndeporzio/RelAxiFast.git $RELAXIFAST_DIR
cd $RELAXIFAST_DIR
git checkout cluster
module load gcc/10.2.0-fasrc01
python scripts/compile_RelAxiFast.py
echo $RELAXIFAST_DIR


echo "PROJECT_DIR: $PROJECT_DIR"                                                
echo "STORAGE_DIR: $STORAGE_DIR"                                                
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"                                
echo "SLURM_CLUSTER_NAME: $SLURM_CLUSTER_NAME"                                  
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"                                
echo "SLURM_JOB_ACCOUNT: $SLURM_JOB_ACCOUNT"                                    
echo "SLURM_JOB_ID: $SLURM_JOB_ID"                                              
echo "SLURM_JOB_NAME: $SLURM_JOB_NAME"                                          
echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES"                                
echo "SLURM_STEP_NUM_TASKS: $SLURM_STEP_NUM_TASKS"                              
echo "SLURM_SUBMIT_DIR : $SLURM_SUBMIT_DIR" 
echo "FORECAST_INDEX : ${SLURM_ARRAY_TASK_ID}"

export PROJECT_DIR                                                              
export STORAGE_DIR                                                              
export RELAXIFAST_DIR                                                                
 

## RUN PROGRAM ##                          
mkdir -p "$STORAGE_DIR/FIG11_$SLURM_ARRAY_TASK_ID"                                     
python $RELAXIFAST_DIR/scripts/Figure_11.py >> "$STORAGE_DIR/FIG11_${SLURM_ARRAY_TASK_ID}/FIG11_${SLURM_ARRAY_TASK_ID}.scriptout" 
cp $RELAXIFAST_DIR/output/result-0/* $STORAGE_DIR/FIG11_$SLURM_ARRAY_TASK_ID
cp $RELAXIFAST_DIR/Figure_11_b1*.txt $STORAGE_DIR/FIG11_$SLURM_ARRAY_TASK_ID
mv "$STORAGE_DIR/FIG11_${SLURM_ARRAY_TASK_ID}.out" "$STORAGE_DIR/FIG11_$SLURM_ARRAY_TASK_ID"
mv "$STORAGE_DIR/FIG11_${SLURM_ARRAY_TASK_ID}.err" "$STORAGE_DIR/FIG11_$SLURM_ARRAY_TASK_ID"  

 
