#!/bin/bash
#SBATCH --job-name="Seas farm"
#SBATCH --time=00-36:0:00
#SBATCH --mem=4G
#SBATCH --account=def-gflowers
#SBATCH --mail-user=tha111@sfu.ca
#SBATCH --mail-type=FAIL,END



# Don't change this line:
task.run
