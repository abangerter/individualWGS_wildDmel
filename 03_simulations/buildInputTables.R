### BUILD INPUT ALLELE FREQUENCY AND READ DEPTH TABLES FOR SIMULATIONS ###


### libraries
library(data.table)
library(foreach)
library(gdsfmt)
library(SeqArray)


### read in important information
  array.list <- readRDS("/scratch/ab5dr/wildDmel2016/simulations/misc/simArray.rds")
  all.snp <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/LDpruned_snpset.txt", header=T)
  codes <- names(array.list)
  
### read in the genofile 
  genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_seqA.gds")

### loop through and build the tables
  for(code in codes) {
    for(chrom in c("2L","2R","3L","3R")){
      ## subset to chromosome & samples
        sampList <- array.list[[code]]
        seqSetFilterPos(genofile, all.snp[chromosome==chrom]$chromosome, all.snp[chromosome==chrom]$position)
        seqSetFilter(genofile, sample.id=sampList, action="intersect")
      ## build the AF info table
        # pull out samps
          samps <- seqGetData(genofile, "sample.id")
        # pull out dosages 
          gt.dt <- t(seqGetData(genofile, "$dosage"))
          gt.dt <- as.data.table(gt.dt)
          setnames(gt.dt, samps)
          gt.dt[,snp.id:=seqGetData(genofile, "annotation/id")]
        # melt 
          gt.dt <- melt(gt.dt, measure.vars=samps, value.name="RefDosage",variable.name="sampID")
        # aggregate to a form where allele frequency can be calculated 
          gt.ag <- gt.dt[,list(refDos=sum(RefDosage, na.rm=T),
                               totalChroms=2*sum(RefDosage<3, na.rm=T)),
                         by=list(snp.id)]
        # calculate allele frequency
          gt.ag[,refAF:=refDos/totalChroms]
          gt.ag[,pos:=tstrsplit(snp.id,"_")[2]]
        # write the AF info table out 
          write.table(gt.ag[,c("snp.id","pos","refAF")], paste("/scratch/ab5dr/wildDmel2016/simulations/rafTables/raf_",code,"_",chrom,".txt",sep=""),
                      col.names=T, row.names=F, sep="\t", quote=FALSE)
      ## build the RD info table
        # pull out the info
          temp.rd <- seqGetData(genofile, "annotation/format/DP")
          RD.dt <- foreach(j=1:dim(temp.rd$data)[1], .combine="cbind") %do% {data.table(temp.rd$data[j,])}
          colnames(RD.dt) <- seqGetData(genofile, "sample.id")
          RD.dt[,snp.id:=seqGetData(genofile, "annotation/id")]
        # melt
          RD.dt <- melt(RD.dt, measure.vars=samps, value.name="RD",variable.name="sampID")
          RD.dt[,pos:=tstrsplit(snp.id,"_")[2]]
        # write the RD info table out 
          write.table(RD.dt, paste("/scratch/ab5dr/wildDmel2016/simulations/RDinfo/RDinfo_",code,"_",chrom,".txt",sep=""), 
                      col.names=T, row.names=F, sep="\t", quote=FALSE)
      seqResetFilter(genofile)
    }
  }

