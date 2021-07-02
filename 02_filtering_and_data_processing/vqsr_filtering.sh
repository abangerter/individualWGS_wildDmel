#! /bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=62000
#SBATCH --time=48:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab


### LOADING MODULES
	echo "loading modules"
	module load gcc/7.1.0
	module load vcftools/0.1.15
	module load bedtools/2.26.0
	module load bcftools/1.9


### DEFINE VARIABLES
	echo "Defining some variables"
	GATKDir="/home/ab5dr/GATK/gatk-4.1.3.0"
	ref="/scratch/ab5dr/wildDmel2016/misc/combined.edit2.fa"
	vcfDir="/scratch/ab5dr/wildDmel2016/vcf"
	miscDir="/scratch/ab5dr/wildDmel2016/misc"
	outDir="/scratch/ab5dr/wildDmel2016/redoAnalysis"


### MAKE TRUTHINESS SET FROM DGRP2-RAW INTERSECT
	echo "Making Truthiness Set"
	bedtools intersect -header \
	-a ${vcfDir}/wild2016raw_revised_v4.vcf \
	-b ${miscDir}/dgrp2.v6lift.regexFix.noIndelRegions.noRep.randomSubset.vcf > \
	${vcfDir}/wild2016raw_revised_v4.truthiness.vcf


### INDEX TRUTHINESS VCF
	${GATKDir}/gatk IndexFeatureFile -F ${vcfDir}/wild2016raw_revised_v4.truthiness.vcf


### VQSR
	## Variant Recalibrator
		echo "Variant Recalibrator"
		${GATKDir}/gatk VariantRecalibrator \
		-R ${ref} \
		-V ${vcfDir}/wild2016raw_revised_v4.vcf \
		--resource:dgrp,known=false,training=true,truth=true,prior=15.0 ${vcfDir}/wild2016raw_revised_v4.truthiness.vcf \
		-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
		-mode SNP \
		--tranches-file ${outDir}/recalibrate_SNP2016subset_revised_v4.tranches \
		-O ${outDir}/recalibrate_SNP2016subset_revised_v4.recal \
	## Apply Recalibration
		echo "Apply Recalibration"
		${GATKDir}/gatk ApplyVQSR \
		-R ${ref} \
		-V ${vcfDir}/wild2016raw_revised_v4.vcf \
		-mode SNP \
		--truth-sensitivity-filter-level 99.0 \
		--tranches-file ${outDir}/recalibrate_SNP2016subset_revised_v4.tranches \
		--recal-file ${outDir}/recalibrate_SNP2016subset_revised_v4.recal \
		-O ${vcfDir}/wild2016raw_revised_v4.vqsr.vcf


### FILTER BY PASS & REMOVING INDELS
	echo "filtering out INDEL"
	vcftools \
	--vcf ${vcfDir}/wild2016raw_revised_v4.vqsr.vcf \
	--remove-indels \
	--recode \
	--recode-INFO-all \
	--out ${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel

	echo "filtering based on PASS"
	vcftools \
	--vcf ${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.recode.vcf \
	--remove-filtered-all \
	--recode \
	--recode-INFO-all \
	--out ${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.pass


### FILTER OUT REGIONS AROUND INDELS
	echo "removing regions around indels"
		## make indel only vcf
			vcftools --vcf ${vcfDir}/wild2016raw_revised_v4.vcf --keep-only-indels --recode --recode-INFO-all --out ${vcfDir}/wild2016raw_revised_v4.indelsOnly
		## pull out INDELs and make into BED
			cat ${vcfDir}/wild2016raw_revised_v4.indelsOnly.recode.vcf | sed -n -e '/^#CHROM/,$p' | cut -f1,2 | awk -v pad=50 '{print $1"\t"$2-pad"\t"$2+pad}' > ${outDir}/wild2016_revised_v4_indel_50bp.bed
		##filter them out
			bedtools intersect -sorted -v -header \
			-a ${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.pass.recode.vcf \
			-b ${outDir}/wild2016_revised_v4_indel_50bp.bed > \
			${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.vcf


### FILTER OUT REPETITIVE REGIONS
	echo "filtering out repetitive regions"
	bedtools intersect -v -header \
	-a ${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.vcf \
	-b ${miscDir}/Dmelv6.repmasker.bed > \
	${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.vcf


### ADD IN SNP ID ANNOTATION 
	bcftools annotate --set-id '%CHROM\_%POS' ${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.vcf --output ${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.vcf --output-type v


### EXTRACT READ DEPTH INFORMATION FOR FILTERING 
	vcftools --vcf ${vcfDir}/wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.vcf  --site-depth --out /scratch/ab5dr/wildDmel2016/redoAnalysis/preRDQF_RD_v4



