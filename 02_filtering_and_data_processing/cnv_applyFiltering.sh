#! /bin/bash


### LOAD MODULES
	echo "loading modules"
	module load gcc/7.1.0
	module load vcftools/0.1.15


vcftools \
	--vcf wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.MIedit.vcf \
	--positions /scratch/ab5dr/wildDmel2016/redoAnalysis/cnv_manta_v4_valrMade.txt \
	--recode \
	--recode-INFO-all \
	--out wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.MIedit.CNVfilt
