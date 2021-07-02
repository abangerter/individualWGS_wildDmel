#! /bin/bash


### LOAD MODULES
	echo "loading modules"
	module load gcc/7.1.0
	module load vcftools/0.1.15

### REMOVE SNPS FROM rd03_filtering.R
vcftools \
	--vcf wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.recode.vcf \
	--exclude /scratch/ab5dr/wildDmel2016/redoAnalysis/snpsForRemoval_indRDfilt_v3.txt \
	--recode \
	--recode-INFO-all \
	--out wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR

### REMOVE NON-BIALLELIC SITES
vcftools \
	--vcf wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.recode.vcf \
	--min-alleles 2 \
	--max-alleles 2 \
	--recode \
	--recode-INFO-all \
	--out wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic

### SUBMIT SCRIPT TO MODIFY DATA TO MISSING DATA & FIX THE HEADER
sbatch ./submitRscript.slurm
grep "^##" wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.recode.vcf >> wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.MIedit.vcf
cat wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.MIedit_noheader.vcf >> wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.MIedit.vcf
