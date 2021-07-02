### EXTRACT READ DEPTH INFORMATION ###


### libraries
library(data.table)
library(gdsfmt)
library(SeqArray)
library(foreach)


### pull out the read depth info
  genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_noFilt_seqA.gds") 
  samps <- seqGetData(genofile, "sample.id")
  rd.ag.all <- foreach(i=c("2L","2R","3L","3R","X"))%do%{
    seqSetFilterChrom(genofile, include=i)
    ## pull out read depth
      temp.rd <- seqGetData(genofile, "annotation/format/DP")
      RD.dt <- foreach(j=1:dim(temp.rd$data)[1], .combine="cbind") %do% {data.table(temp.rd$data[j,])}
      colnames(RD.dt) <- seqGetData(genofile, "sample.id")
      RD.dt[,snp.id:=seqGetData(genofile, "annotation/id")]
      RD.dt[,c("chrom","pos"):=tstrsplit(snp.id,"_")]
      RD.dt <- melt(RD.dt, measure.vars=samps, value.name="RD",variable.name="sampID")
      RD.dt[,RD:=as.numeric(RD)]
    ## get RD summary info
      rd.ag <- RD.dt[,list(meanRD=mean(RD, na.rm=T),
                           medianRD=median(RD, na.rm=T),
                           varianceRD=var(RD, na.rm=T),
                           quantile25=quantile(RD, 0.25, na.rm=T),
                           quantile75=quantile(RD, 0.75, na.rm=T),
                           maxRD=max(RD, na.rm=T)),
                     by=list(sampID, chrom)]
      seqResetFilter(genofile)
      rd.ag
  }
  rd.ag.all <- rbindlist(rd.ag.all)
  
### write out the table
  write.table(rd.ag.all, "/scratch/ab5dr/wildDmel2016/redoAnalysis/v4RDinfo.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)
