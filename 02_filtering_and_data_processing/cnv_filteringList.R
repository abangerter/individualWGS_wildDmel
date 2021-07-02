### GENERATE BEDFILE OF DE NOVO STRUCTURAL VARIANTS TO REMOVE ###


### libraries 
library(data.table)
library(gdsfmt)
library(SeqArray)
library(tidyr)
library(valr)


### read in the data
  cnv.dt <- fread("/scratch/ab5dr/wildDmel2016/cnvOut/processMantaV4/all_SVscored_subset.txt", header=F)
  setnames(cnv.dt, c("sampID","chrom","start","SVtype","ref","alt","qual","pass","info"))
  genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_mostlyFilt_seqA.gds")

### prep snp info
  snps.dt <- data.table(chrom=seqGetData(genofile, "chromosome"),
                        start=seqGetData(genofile, "position"))
  snps.dt[,end:=start]
  snps.tbl <- as.tbl_interval(snps.dt)
  
### filter out non-passing individuals & chroms we don't want
    cnv.dt <- cnv.dt[chrom %in% c("X","2L","2R","3L","3R")]
    cnv.dt <- cnv.dt[pass=="PASS"]
    
### manipulate the data into a more workable form
  ## get the SV type
    cnv.dt[,SVtype:=tstrsplit(SVtype, ":")[1]]
    cnv.dt[,SVtype:=tstrsplit(SVtype,"Manta")[2]]
    cnv.dt <- cnv.dt[SVtype %in% c("DEL","INS","DUP")]
  ## get the start, stop, and length of the SV 
    cnv.dt[,stop:=tstrsplit(info,";")[1]]
    cnv.dt[,length:=tstrsplit(info,";")[3]]
    cnv.dt[,c("stopType","stop"):=tstrsplit(stop,"=")]
    cnv.dt[,c("lengthType","length"):=tstrsplit(length,"=")]
    cnv.dt <- cnv.dt[lengthType=="SVLEN"]
    cnv.dt[,stop:=as.numeric(stop)]
    cnv.dt[,length:=as.numeric(length)]
  ## aggregate 
    cnv.ag <- cnv.dt[,.N, by=list(chrom, start, stop, SVtype, length)]
  ## pull out the snps that will be kept after filtering
    # pull out CNVs to filter with
      cnv.ag.tmp1 <- cnv.ag[abs(length) <= 10000]
      cnv.ag.tmp2 <- cnv.ag[abs(length) > 10000][N>=15]
      cnv.ag.tmp <- rbind(cnv.ag.tmp1, cnv.ag.tmp2)
    # get intervals around CNVs 
      cnv.ag.tmp[,start:=as.numeric(start)]
      cnv.ag.tmp[,begin:=start-50]
      cnv.ag.tmp[,end:=stop+50]
      cnv.ag.sub <- cnv.ag.tmp[,c("chrom","begin","end")]
      setnames(cnv.ag.sub, old="begin", new="start")
    # turn into valr compatible format 
      cnv.tbl <- as.tbl_interval(cnv.ag.sub)
    # collapse overlapping CNVs
      cnv.tbl.merge <- bed_merge(cnv.tbl)
    # inverse intersect with snps.tbl
      snps.mat <- as.data.table(as.matrix(bed_intersect(snps.tbl, cnv.tbl.merge, invert=T)))
      snps.mat.ag <- snps.mat[,list(nSNPs=length(start)), chrom]
      
### write output file of good snps to keep
  write.table(snps.mat[,c("chrom","start")], "/scratch/ab5dr/wildDmel2016/redoAnalysis/cnv_manta_v4_valrMade.txt", col.names=F, row.names=F, sep="\t", quote=FALSE)
