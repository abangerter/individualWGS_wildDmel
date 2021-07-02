#! /bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=40000
#SBATCH --time=48:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab

module load manta/1.4.0
module load samtools

lineNo=${SLURM_ARRAY_TASK_ID}

sample=$( sed -n ${lineNo}p /scratch/ab5dr/wildDmel2016/redoAnalysis/keepList_v4.txt )

### DEFINING DIRECTORIES
	inputDir=/scratch/ab5dr/wildDmel2016/mappedBamNoDup
	listDir=/scratch/ab5dr/wildDmel2016/fileList
	interDir=/scratch/ab5dr/wildDmel2016/interCNVbam
	outDir=/scratch/ab5dr/wildDmel2016/cnvOut/v4manta
	finalDir=/scratch/ab5dr/wildDmel2016/cnvOut/v4outManta
	ref=/scratch/ab5dr/wildDmel2016/misc/combined.edit2.fa

### MERGE BAM FILES
echo "merging files"
ls ${inputDir}/*${sample}*.bam > ${listDir}/${sample}.txt
samtools merge -b ${listDir}/${sample}.txt ${interDir}/${sample}.merged.bam
samtools index ${interDir}/${sample}.merged.bam

### MANTA SV CONFIG
echo "manta prep"
mkdir ${outDir}/${sample}
python $EBROOTMANTA/bin/configManta.py \
--bam ${interDir}/${sample}.merged.bam \
--referenceFasta ${ref} \
--runDir ${outDir}/${sample}

### MANTA SV CALLING
echo "manta run"
${outDir}/${sample}/runWorkflow.py -m local -j 1

### KEEP ONLY THE IMPORTANT INFOß
mkdir ${finalDir}/${sample}
mv ${outDir}/${sample}/results/variants/* ${finalDir}/${sample}
rm -r ${outDir}/${sample}/*
