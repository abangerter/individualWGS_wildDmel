#! /usr/bin/env Rscript
### turning problematic site-sample combos into missing data


print("start")
### libraries
library(data.table)
library(foreach)


### read in relevant data
  print("reading in data")
  vcf <- fread("/scratch/ab5dr/wildDmel2016/vcf/wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.recode.vcf", skip="#CHROM", header=T)
  sites <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/snpsToBeMissData_indRDfilt_v4.txt", header=T)
  setnames(sites, old="snp.id",new="ID")
  sites[,replace:=TRUE]
  samps <- scan("/scratch/ab5dr/wildDmel2016/redoAnalysis/keepList_v4.txt", character(), quote="")
  
    
### loop through to make edits
  print("making the edits")
  newvcf <- foreach(i=samps, .combine="cbind")%do%{
    print(i)
    # subset vcf to just the snp ID column and the sample 
      vcf.temp <- vcf[,c("ID", i), with=FALSE]
      setkey(vcf.temp, ID)
    # subset filter sites to sample
      sites.temp <- sites[sampID==i]
      setkey(sites.temp, ID)
    # merge vcf subset with the sites to replace
      temp <- merge(vcf.temp, sites.temp, all.x=TRUE)
      setnames(temp, old=i, new="genotype_before")
    # replace GT with missing data
      temp[replace==TRUE, genotype_after := gsub("^[0-9.]/[0-9.]:", "./.:", genotype_before)]
      temp[is.na(replace), genotype_after := genotype_before]
    # restructure data table
      temp[,c("chr","pos"):=tstrsplit(ID,"_")]
      temp[,pos:=as.numeric(pos)]
      temp[,chr:=factor(chr, levels=c("X","2L","2R","3L","3R","4","Y"))]
      setkey(temp, chr, pos)
      temp <- temp[,c("genotype_after")]
      setnames(temp, old="genotype_after", new=i)
    # output in some useable format 
      temp
  }
  saveRDS(newvcf, file="vcf_almost")
  
#  newvcf <- cbindlist(newvcf)
  newvcf.noheader <- cbind(vcf[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")], newvcf)


### write out the new vcf file
  print("writing out new vcf without the header")
  write.table(newvcf.noheader, "/scratch/ab5dr/wildDmel2016/vcf/wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.PFsnpR.biAllelic.MIedit_noheader.vcf", col.names=T, row.names=F, sep="\t", quote=FALSE)
  
  
  
  