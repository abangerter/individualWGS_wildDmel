#! /bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=62000
#SBATCH --time=70:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab


### DEFINING DIRECTORIES 
	echo 'defining variables'
	ref="/scratch/ab5dr/wildDmel2016/misc/combined.edit2.fa"
	inputDir="/scratch/ab5dr/wildDmel2016/gvcfNoDup"
	miscDir="/scratch/ab5dr/wildDmel2016/misc"
	GATKDir="/home/ab5dr/GATK/gatk-4.1.3.0"
	outDir="/scratch/ab5dr/wildDmel2016/vcf"
	
### MAKING THE LIST
ls ${inputDir}/*.noDup.g.vcf > ${miscDir}/gvcfbatch_goodSamp_revised_v4.list

### COMBINING GVCFs FOR NEXT STEP
	echo "combine GVCFs"
	${GATKDir}/gatk CombineGVCFs \
		-R ${ref} \
		-L ${miscDir}/wild2016chr_revised.intervals \
		--variant ${miscDir}/gvcfbatch_goodSamp_revised_v4.list \
		-O ${outDir}/batch_revised_v4.g.vcf

### GENOTYPE GVCF
	echo "genotype GVCFs"
	${GATKDir}/gatk GenotypeGVCFs \
		-R ${ref} \
		-L ${miscDir}/wild2016chr_revised.intervals \
		--heterozygosity 0.01 \
		--indel-heterozygosity 0.001 \
		-V ${outDir}/batch_revised_v4.g.vcf \
		-O ${outDir}/wild2016raw_revised_v4.vcf
