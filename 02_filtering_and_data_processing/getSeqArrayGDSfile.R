### GENERATE GDS FILE FOR USE IN ANALYSES  


### libraries
library(gdsfmt)
library(SeqArray)


### make a GDS file from a less filtered file for filtering and other things
  vcf.fn <- "wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.vcf"
  seqVCF2GDS(vcf.fn, "/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_noFilt_seqA.gds", info.import="ID")


### make the main and final GDS file 
  vcf.fn <- "wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.MIedit.CNVfilt.recode.vcf"
  seqVCF2GDS(vcf.fn, "/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_seqA.gds", info.import="ID")
