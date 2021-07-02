#! /bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8000
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab

module load gcc/7.1.0
module load R/3.5.1


echo "load in array text"
## switch between which is the input argument file depending on if doing Full Sample vs population simulations
argFile="/scratch/ab5dr/wildDmel2016/simulations/scripts/fs_array.txt"
#argFile="/scratch/ab5dr/wildDmel2016/simulations/scripts/pop_array.txt"

echo "define line of array to use as argument"
lineNo=$SLURM_ARRAY_TASK_ID

echo "pick out the sample file ID"
myArgs=$(sed -n ${lineNo}p ${argFile})

echo "starting simulation"
Rscript /scratch/ab5dr/wildDmel2016/simulations/scripts/simulations.R ${myArgs}