#! /bin/bash


### LOAD MODULES
	echo "loading modules"
	module load gcc/7.1.0
	module load vcftools/0.1.15

### FILTER
vcftools \
	--vcf wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.vcf \
	--remove /scratch/ab5dr/wildDmel2016/redoAnalysis/rmList_v4.txt \
	--recode \
	--recode-INFO-all \
	--out wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt

vcftools \
	--vcf wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.recode.vcf \
	--exclude-positions /scratch/ab5dr/wildDmel2016/redoAnalysis/RD_siteFilterList_v3.txt \
	--recode \
	--recode-INFO-all \
	--out wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q
