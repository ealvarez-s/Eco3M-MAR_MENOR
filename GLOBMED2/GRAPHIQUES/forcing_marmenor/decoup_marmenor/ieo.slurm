#!/bin/bash
#SBATCH --job-name=R_meanYY
#SBATCH -N 1
#SBATCH --nodelist=c02
#SBATCH --partition=hpc
#SBATCH --time=06:00:00
####SBATCH --time=2-01:00:00
####SBATCH --time=01:00:00
####SBATCH -t 1-00:00:00
####SBATCH --time=00:05:00

### cd $SLURM_SUBMIT_DIR
conda activate seamless-bb-r
export LD_LIBRARY_PATH=/home/alvarez/miniconda3/envs/seamless-bb-r/lib:$LD_LIBRARY_PATH
export SLURM_CONF=/home/alvarez/slurm.conf

date
Rscript make_mean_year.R 2015 /home/alvarez/EVA/Eco3M-MAR_MENOR/GLOBMED2/GRAPHIQUES/forcing_marmenor/mean_year/
date
