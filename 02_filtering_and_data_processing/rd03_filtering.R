### READ DEPTH BASED FILTERING, PART 3 ###


### libraries
library(gdsfmt)
library(SeqArray)
library(data.table)
library(foreach)


############################################################################
### BY-INDIVIDUAL, BY-SNP FILTERING & PREP FOR MISSING DATA MODIFICATION ###
############################################################################


### generate GDS from partially filtered V4 VCF - aka rm dup one 
  vcf.fn <- "wild2016raw_revised_v4.vqsr.noIndel.pass.noIndelReg50bp.noRep.snpID.RDindFilt.UL015q.recode.vcf"
  seqVCF2GDS(vcf.fn, "/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_partFilt_seqA.gds", info.import="ID")


### read in info
  genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_partFilt_seqA.gds")
  samps <- seqGetData(genofile, "sample.id")
  ## pull out the RD data
    temp.rd <- seqGetData(genofile, "annotation/format/DP")
    RD.dt <- foreach(j=1:dim(temp.rd$data)[1], .combine="cbind") %do% {data.table(temp.rd$data[j,])}
    colnames(RD.dt) <- seqGetData(genofile, "sample.id")
    RD.dt[,snp.id:=seqGetData(genofile, "annotation/id")]
    RD.dt[,c("chrom","pos"):=tstrsplit(snp.id,"_")]
  ## melt 
    rd.melt <- melt(RD.dt, measure.vars=samps, value.name="RD",variable.name="sampID")

### aggregate medians and quantiles - turn into pass/fail measure
  ## aggregate to get median and quantile
    rd.ag <- rd.melt[RD>0][,list(medRD=median(as.numeric(RD), na.rm=T),
                                 lower01=quantile(RD, 0.01, na.rm=T),
                                 upper01=quantile(RD, 0.99, na.rm=T)), 
                           by=list(chrom, sampID)]
  ## figure out how to code the pass/fail based on quantile
    rd.ag <- rd.ag[,c("chrom","sampID","lower01","upper01")]
    setkey(rd.melt, sampID, chrom)
    setkey(rd.ag, sampID, chrom)
    rd.out <- merge(rd.melt, rd.ag)
    rd.out[,filter:= ifelse(RD<=lower01, "FAIL", ifelse(RD>=upper01, "FAIL", "PASS"))]
  ## aggregate to counts of failure rates
    rd.out.ag <- rd.out[,list(countFail=sum(filter=="FAIL")), by=snp.id]
    write.table(rd.out.ag, "/scratch/ab5dr/wildDmel2016/redoAnalysis/readDepthSiteSampleFilt_agg_v4.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)
  ## visualize the aggregate of failure rates 
    hist(rd.out.ag$countFail, breaks=50)
  ## where to put the lower threshold of one-off samples that need to turn into missing data? 
    quantile(rd.out.ag$countFail, 0.90, na.rm=T)
      # 90% - 33
      # 85% - 23
    rd.out.ag[countFail>20]
      # results in 245379 SNPs being removed 
    rd.out.ag[,c("chrom","pos"):=tstrsplit(snp.id,"_")]
    rd.out.ag[countFail>20][,list(count=sum(countFail<152)), by=list(chrom)]
      #    chrom count
      # 1:    2L 45660
      # 2:    2R 31019
      # 3:    3L 46760
      # 4:    3R 47323
      # 5:     4  2397
      # 6:     X 71652
      # 7:     Y   374
    write.table(rd.out.ag[countFail>20][,c("snp.id")], "/scratch/ab5dr/wildDmel2016/redoAnalysis/snpsForRemoval_indRDfilt_v4.txt", col.names=F, row.names=F, sep="\t", quote=FALSE)

    
##################################################
### GET LIST OF SNPS TO TURN INTO MISSING DATA ###
##################################################
    
### continue to use the tables made in prior code block, but for a different filtering step
### list of snps to turn into missing data
  ## merge the rd ag and the rd table
    setkey(rd.out.ag, snp.id)
    setkey(rd.out, snp.id)
    rd.again <- merge(rd.out, rd.out.ag[,c("snp.id","countFail")])
  ## get that list
    write.table(rd.again[countFail<21][filter=="FAIL"][,c("snp.id","sampID")], "/scratch/ab5dr/wildDmel2016/redoAnalysis/snpsToBeMissData_indRDfilt_v4.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)

