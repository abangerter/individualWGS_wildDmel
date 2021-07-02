#! /bin/bash 

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=62000
#SBATCH --time=1:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab

module load gcc
module load R/3.5.1

# define file name & argument from array
echo "load in array text"
argFile="/scratch/ab5dr/wildDmel2016/simulations/scripts/slide_array.txt"

echo "define line of array to use as argument"
lineNo=$SLURM_ARRAY_TASK_ID

echo "pick out the args"
myArgs=$(sed -n ${lineNo}p ${argFile})

echo "starting the analysis"
Rscript /scratch/ab5dr/wildDmel2016/simulations/scripts/slidingWindowFis.R ${myArgs}
