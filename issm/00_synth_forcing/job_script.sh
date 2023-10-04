#!/bin/bash
#SBATCH --job-name="ISSM"
#SBATCH --time=00-03:00
#SBATCH --mem=2G
#SBATCH --account=def-gflowers
#SBATCH --mail-user=tha111@sfu.ca
#SBATCH --mail-type=FAIL,END

module load matlab/2022b
matlab -nodesktop -nosplash -r "run_job(0.05, 0.5, 5./4., 3./2., 1e-4, 0); exit"


