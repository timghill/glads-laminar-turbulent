#!/bin/bash
#SBATCH --job-name="glads-alpha"
#SBATCH --time=00-24:0:00
#SBATCH --mem=8G
#SBATCH --account=def-gflowers
#SBATCH --mail-user=tha111@sfu.ca
#SBATCH --mail-type=FAIL,END

run_job(0.2, 10, 4,   2,   0,      6)

