### EXTRACT INFORMATION ABOUT HETEROZYGOSITY LEVELS ###


### libraries 
library(data.table)
library(gdsfmt)
library(SeqArray)


### read in genofile 
  genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_seqA.gds")
  
### extract heterozygosity information 
  ## pull out dosage 
    gt.dt <- t(seqGetData(genofile, "$dosage"))
    gt.dt <- as.data.table(gt.dt)
    samps <- seqGetData(genofile, "sample.id")
    setnames(gt.dt, samps)
    gt.dt[,snp.id:=seqGetData(genofile, "annotation/id")]
    gt.dt[,c("chrom","pos"):=tstrsplit(snp.id,"_")]
    gt.dt <- melt(gt.dt, measure.vars=samps, value.name="dosage",variable.name="sampID")
  ## aggregate to count number of heterozygotes per sample per chrom
    gt.ag <- gt.dt[,list(Nhets=sum(dosage==1, na.rm=T)), by=list(sampID, chrom)]

### write out a heterozygosity table
  write.table(gt.ag, "/scratch/ab5dr/wildDmel2016/redoAnalysis/nHets_byChr_v4.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)
