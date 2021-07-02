### READ DEPTH BASED FILTERING, PART 1 ###


### libraries
library(gdsfmt)
library(SeqArray)
library(data.table)
library(foreach)



###############################################
### STEP 1: INDIVIDUAL READ DEPTH FILTERING ###
###############################################

### generate GDS from V4 VCF - aka rm dup one 
  vcf.fn <- "wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.vcf"
  seqVCF2GDS(vcf.fn, "/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_noFilt_seqA.gds", info.import="ID")

### pull out the basic info
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
  write.table(rd.ag.all, "/scratch/ab5dr/wildDmel2016/redoAnalysis/v4RDinfo.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)

### who will be kept?
  ## read in the table just made 
    rd.dt <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/v4RDinfo.txt", header=T)
  ## samples to remove
    intersect(rd.dt[chrom=="2L"][medianRD<4][,sampID], intersect(rd.dt[chrom=="2R"][medianRD<4][,sampID], intersect(rd.dt[chrom=="3L"][medianRD<4][,sampID], rd.dt[chrom=="3R"][medianRD<4][,sampID])))
      #107 samples to remove
    rm.vec <- intersect(rd.dt[chrom=="2L"][medianRD<4][,sampID], intersect(rd.dt[chrom=="2R"][medianRD<4][,sampID], intersect(rd.dt[chrom=="3L"][medianRD<4][,sampID], rd.dt[chrom=="3R"][medianRD<4][,sampID])))
    rm.dt <- data.table(sampID=rm.vec)
    write.table(rm.dt, "/scratch/ab5dr/wildDmel2016/redoAnalysis/rmList_v4.txt", col.names=F, row.names=F, sep="\t", quote=FALSE)
  ## what samples will be kept? 
    temp <- data.table(sampID=unique(rd.dt[,sampID]))
    temp <- temp[! sampID %in% rm.vec]
    write.table(temp, "/scratch/ab5dr/wildDmel2016/redoAnalysis/keepList_v4.txt", col.names=F, row.names=F, sep="\t", quote=FALSE)


#################################################
### STEP 2: OVERALL SITE READ DEPTH FILTERING ###
#################################################

### pull out RD data 
  totalRD <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/preRDQF_RD_v4.ldepth", header=T)
    
### look at histogram of per-site total read depth
  ggplot(totalRD[! CHROM %in% c("X","4","Y")], aes(x=SUM_DEPTH)) + geom_histogram(bins=100)
  ggplot(totalRD[! CHROM %in% c("4","Y")], aes(x=SUM_DEPTH)) + geom_histogram(bins=100) + facet_grid(CHROM~.)
  totalRD[,chrType:=ifelse(CHROM %in% c("4","Y"), "other", ifelse(CHROM=="X", "sex", "autosome"))]
  ggplot(totalRD[! CHROM %in% c("4","Y")], aes(x=SUM_DEPTH)) + geom_histogram(bins=100) + facet_grid(chrType~.)

### toy around with what the quantile limits should be 
  ggplot(totalRD[! CHROM %in% c("X","4","Y")], aes(x=SUM_DEPTH)) + geom_histogram(bins=100) + geom_vline(xintercept=quantile(totalRD[! CHROM %in% c("X","4","Y")]$SUM_DEPTH, 0.985)) + geom_vline(xintercept=quantile(totalRD[! CHROM %in% c("X","4","Y")]$SUM_DEPTH, 0.015))
  ggplot(totalRD[CHROM=="X"], aes(x=SUM_DEPTH)) + geom_histogram(bins=100) + geom_vline(xintercept=quantile(totalRD[CHROM=="X"]$SUM_DEPTH, 0.985)) + geom_vline(xintercept=quantile(totalRD[CHROM=="X"]$SUM_DEPTH, 0.015))

### quantile values
  quantile(totalRD[! CHROM %in% c("X","4","Y")]$SUM_DEPTH, 0.985)   # 6274
  quantile(totalRD[! CHROM %in% c("X","4","Y")]$SUM_DEPTH, 0.015)   # 282
  quantile(totalRD[CHROM=="X"]$SUM_DEPTH, 0.985)  # 3462
  quantile(totalRD[CHROM=="X"]$SUM_DEPTH, 0.015)  # 125

### pull out the list of sites
  upper015 <- totalRD[! CHROM %in% c("X","4","Y")][SUM_DEPTH>6274]
  lower015 <- totalRD[! CHROM %in% c("X","4","Y")][SUM_DEPTH<282]
  upper015X <- totalRD[CHROM=="X"][SUM_DEPTH>3462]
  lower015X <- totalRD[CHROM=="X"][SUM_DEPTH<125]
  exclude.dt <- rbindlist(list(upper015, lower015, upper015X, lower015X))
  exclude.dt <- exclude.dt[,c("CHROM","POS")]
  exclude.dt[,list(total=sum(POS>0)), by=CHROM]
  write.table(exclude.dt, "/scratch/ab5dr/wildDmel2016/redoAnalysis/RD_siteFilterList_v4.txt", col.names=F, row.names=F, sep="\t", quote=FALSE)



