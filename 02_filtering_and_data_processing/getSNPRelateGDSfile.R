### GENERATE SNPRELATE GDS FILE ###


### libraries 
library(gdsfmt)
library(SNPRelate)


### make a SNPRelate GDS file 
  vcf.fn <- "wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.MIedit.CNVfilt.recode.vcf"
  snpgdsVCF2GDS(vcf.fn, "/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_snpr.gds", method=c("copy.num.of.ref"), snpfirstdim = FALSE)
