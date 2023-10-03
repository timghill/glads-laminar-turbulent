#!/bin/bash
#SBATCH --job-name="ISSM"
#SBATCH --time=00-06:00
#SBATCH --mem=4G
#SBATCH --account=def-gflowers
#SBATCH --mail-user=tha111@sfu.ca
#SBATCH --mail-type=FAIL,END

module load matlab/2022b
matlab -nodesktop -nosplash -r "run_job; exit"
