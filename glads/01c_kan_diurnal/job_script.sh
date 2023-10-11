#!/bin/bash
#SBATCH --job-name="Diurnal"
#SBATCH --time=03-0:0:00
#SBATCH --mem=16G
#SBATCH --account=def-gflowers
#SBATCH --mail-user=tha111@sfu.ca
#SBATCH --mail-type=FAIL,END



# Don't change this line:
task.run
