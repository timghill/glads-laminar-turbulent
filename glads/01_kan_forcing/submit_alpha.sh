#!/bin/bash
#SBATCH --job-name="glads-alpha"
#SBATCH --time=00-48:0:00
#SBATCH --mem=8G
#SBATCH --account=def-gflowers
#SBATCH --mail-user=tha111@sfu.ca
#SBATCH --mail-type=FAIL,END

module load matlab/2022a; matlab -nodisplay -r 'try; run_job_alpha(0.2, 10, 4, 2, 0, 6); quit; catch ME; throw(ME); quit(1); end;'

