#!/bin/bash -l
#SBATCH --job-name=get_WF
#SBATCH --account=rrg-bojana_cpu
#SBATCH --time=8:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=2000G

module load matlab
matlab -nodisplay -r  "FILE NAME"
