#!/bin/bash

# Runs iSNPcaller_MT on files in IMPB directory

#SBATCH --time 10-00:00:00
#SBATCH --job-name=iSNPcaller_MT
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --partition=<partition> 
#SBATCH --mem=180GB
#SBATCH --mail-type ALL
#SBATCH	-A <user_account>
#SBATCH --mail-type ALL
#SBATCH --mail-user <user@email.com>

echo "SLURM_NODELIST: "$SLURM_NODELIST
echo "PWD :" $PWD

perl iSNPcaller_MT.pl <project_name> run
