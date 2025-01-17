#!/bin/bash

# Run ARGweaver MCMC

#SBATCH --time 14-00:00:00
#SBATCH --job-name=ARGweaver
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 32
#SBATCH --partition=normal
#SBATCH --mem=180GB
#SBATCH --mail-type ALL
#SBATCH	-A coa_farman_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST
echo "PWD :" $PWD


## ARGUMENTS

sitesfile=$1

mutrate=$2

recombrate=$3

region=$4

outprefix=${sitesfile/\.sites/}_${mutrate}_${recombrate}_${region}/${sitesfile/\.sites/}_${region}

## RUN

source /project/farman_uksr/miniconda/etc/profile.d/conda.sh

conda activate bioinfo

arg-sample -s $sitesfile --ntimes 100  -n 200 --sample-step 1 -m $mutrate -r $recombrate -o $outprefix --region $region --overwrite

conda deactivate
