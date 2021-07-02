#! /bin/bash


echo "${1}"


### define variables
	echo "defining variables"
	ref="/scratch/ab5dr/wildDmel2016/misc/combined.edit2.fa"
	inputDir="/scratch/ab5dr/wildDmel2016/mappedBamNoDup"
	interDir="/scratch/ab5dr/wildDmel2016/interBamNoDup"
	gvcfOutDir="/scratch/ab5dr/wildDmel2016/gvcfNoDup"
	miscDir="/scratch/ab5dr/wildDmel2016/misc"
	GATKDir="/home/ab5dr/GATK/gatk-4.1.3.0"
	fileDir="/scratch/ab5dr/wildDmel2016/fileList"

### BASE RECALIBRATION
	fileType=($( ls ${inputDir}/${1}*.bam | cut -d "." -f4 | sort | uniq ))
	for type in "${fileType[@]}"; do
		echo "${type}"
			### base recalibration
			echo 'BaseRecalibrator'
			${GATKDir}/gatk BaseRecalibrator \
				-L ${miscDir}/wild2016chr_revised.intervals \
				-R ${ref} \
				-I ${inputDir}/${1}.${type}.merge.noDup.bam \
				--known-sites ${miscDir}/dgrp2.v6lift.regexFix.vcf \
				-O ${interDir}/${1}_noDup_recal_data.${type}.table

			### print reads step - merges bam file with recalibration table
			echo 'apply base recalibration'
			${GATKDir}/gatk ApplyBQSR \
	   			-L ${miscDir}/wild2016chr_revised.intervals \
				-R ${ref} \
				-I ${inputDir}/${1}.${type}.merge.noDup.bam \
				-bqsr ${interDir}/${1}_noDup_recal_data.${type}.table \
				-O ${interDir}/${1}.${type}.noDupRecal.bam
	done


### FILE PREP
	ls ${interDir}/${1}.*.noDupRecal.bam > ${fileDir}/${1}.noDup.list


### GVCF GENERATION
	echo 'HaplotypeCaller'
	${GATKDir}/gatk HaplotypeCaller \
		-R ${ref} \
		-I ${fileDir}/${1}.noDup.list \
  		-ERC GVCF \
		-L ${miscDir}/wild2016chr_revised.intervals \
		--heterozygosity 0.01 \
		--indel-heterozygosity 0.001 \
		-O ${gvcfOutDir}/${1}.noDup.g.vcf
