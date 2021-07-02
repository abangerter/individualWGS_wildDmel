#! /bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8000
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab

# define file name & argument from array
echo "load in array text"
argFile="/scratch/ab5dr/wildDmel2016/scripts/allSampList.txt"

echo "define line of array to use as argument"
lineNo=$SLURM_ARRAY_TASK_ID

echo "pick out the sample file ID"
sampID=$(sed -n ${lineNo}p ${argFile})

echo "start the mapping"
./map_reads2016.sh ${sampID}