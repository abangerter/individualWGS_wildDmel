#! /bin/bash 

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=62000
#SBATCH --time=48:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab

module load gcc/7.1.0
module load R/3.5.1

Rscript /scratch/ab5dr/wildDmel2016/scripts/modifyToMissing.R
