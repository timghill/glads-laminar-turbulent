#!/bin/bash
#SBATCH --job-name="Convert"
#SBATCH --time=00-01:00
#SBATCH --mem=2G
#SBATCH --account=def-gflowers
#SBATCH --mail-user=tha111@sfu.ca
#SBATCH --mail-type=FAIL,END

module load matlab/2022b
matlab -nodesktop -nosplash -r "do_convert; exit"
