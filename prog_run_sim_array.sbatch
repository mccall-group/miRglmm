#!/bin/bash
#SBATCH -J prog_run_sim_array
#SBATCH -e prog_run_sim_array_%j_err.log
#SBATCH -o prog_run_sim_array_%j.log
#SBATCH --mem=20G
#SBATCH -t 1-00:00:00
#SBATCH --array=1-100
#SBATCH -p standard
#SBATCH -n 1 -c 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrea_baran@urmc.rochester.edu
module load r/4.1.1/b1
Rscript prog_run_sim_array.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK

