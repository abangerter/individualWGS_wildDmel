### INVERSION TYPING & FIGURE 2 ###


### libraries
library(data.table)
library(gdsfmt)
library(SeqArray)
library(LEA)


####################################
### PROBALISTIC INVERSION TYPING ###
####################################

### Load in data
  ## Open up the file 
    genofile <- snpgdsOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_snpr.gds")
  ## Get sample info into a data table
    snp.dt <- data.table(snp.id=read.gdsn(index.gdsn(genofile, "snp.id")),
                         chr=read.gdsn(index.gdsn(genofile, "snp.chromosome")),
                         pos=read.gdsn(index.gdsn(genofile, "snp.position")),
                         alleles=read.gdsn(index.gdsn(genofile, "snp.allele")))
  ## load in inversion markers
    fixedInv <- read.delim("/scratch/ab5dr/wildDmel2016/redoAnalysis/v6liftOver_inv_fixed_coord.txt", header=F, sep=" ", as.is=T)
    fixedInv <- as.data.frame(t(fixedInv))
    fixedInv$chr <- substr(fixedInv[,1], 4, 5)
    fixedInv$pos <- as.numeric(as.character(fixedInv[,2]))
    fixedInv <- as.data.table(fixedInv)
  ## merge the tables
    setkey(snp.dt, chr, pos)
    setkey(fixedInv, chr, pos)
    snp.inv <- merge(snp.dt, fixedInv, all.y=T)

    
### Parse genotypes
    ## get allelic status
    genomat.orig <- snpgdsGetGeno(genofile, snp.id=na.omit(snp.inv$snp.id), with.id=T)
    genomat <- data.table(genotype=expand.grid(as.matrix(genomat.orig$genotype))$Var1,
                          sample.id=rep(genomat.orig$sample.id, length(genomat.orig$snp.id)),
                          snp.id=rep(genomat.orig$snp.id, each=length(genomat.orig$sample.id)))
    setkey(genomat, snp.id)
  ## merge back with inversion names
    setkey(snp.inv, snp.id)
    genomat <- merge(genomat, snp.inv)
  ## split alt and ref
    genomat[,c("ref","alt"):=tstrsplit(alleles,"/")]
  ## throw out remaining indels
    genomat <- genomat[nchar(alt)%in%c(1)] 
  ## is "inversion" allele in alt?
    genomat[nchar(alt)==1, invCheck := V3==alt]
    genomat <- genomat[invCheck==T]
  ## deal with monomorphic sites from fixedInv
    setkey(genomat, chr, pos)
    genomat.sites <- genomat[!duplicated(genomat)][,c("chr", "pos"), with=F]
    genomat.sites[,poly:=T]
    setkey(genomat.sites, chr, pos)
    setkey(fixedInv, chr, pos)
    gf <- merge(genomat.sites, fixedInv, all.x=T, all.y=T)
    gf[is.na(poly), poly:=F]  
  ## gather information into table
    mono <- data.table(snp.id=NA,
                       genotype=2,
                       sample.id=rep(unique(genomat$sample.id), each=sum(!gf$poly)),
                       chr=rep(gf[poly==F]$chr, length(unique(genomat$sample.id))),
                       pos=rep(gf[poly==F]$pos, length(unique(genomat$sample.id))),
                       alleles=NA,
                       V1=rep(gf[poly==F]$V1, length(unique(genomat$sample.id))),
                       V2=rep(gf[poly==F]$V2, length(unique(genomat$sample.id))),
                       V3=rep(gf[poly==F]$V3, length(unique(genomat$sample.id))),
                       ref=NA,
                       alt=NA,
                       invCheck=T)
    genomat <- rbind(genomat, mono)

    
### Make inversion status calls
    ## tabulate calls
    genomat.ag <- genomat[,list(n.std=sum(genotype==2, na.rm=T),
                                n.het=sum(genotype==1, na.rm=T),
                                n.inv=sum(genotype==0, na.rm=T)),
                          list(sample.id, inversion=V1)]
    genomat.ag <- as.data.table(gather(genomat.ag, kary, n, n.std:n.inv))
  ## make most-likely call
    kary <- genomat.ag[,list(kary=kary[which.max(n)],
                             evidence=max(n)/sum(n),
                             n=sum(n),
                             prStd=(2*n[kary=="n.std"] + n[kary=="n.het"])/(2*sum(n))),
                       list(sample.id, inversion)]
  ## make a plot
    ggplot(data=kary, aes(x=inversion, y=prStd, color=inversion)) + geom_boxplot() + geom_hline(yintercept=.9)
  ## check double inversion status on one chromosome
    kary[,chr:=gsub("In\\(([2,3,X,L,R]{1,})\\)[A-Za-z]{1,}", "\\1", inversion)]						
    kary.ag <- kary[,list(n=sum(prStd<.9)), list(chr, sample.id)]
    table(kary.ag$n)
  ## MAKE THE FINAL CALL
    karyCall <- function(prStds, inv, cutVal=.75) {
      if(all(prStds>=cutVal)) {
        o <- "std"
      } else if (sum(prStds <= (1-cutVal))==1) {
        o <- as.character(inv[prStds<=(1-cutVal)])
      } else if (sum(prStds>=(1-cutVal) & prStds<=cutVal)==2) {
        o <- paste(inv[prStds>=(1-cutVal) & prStds<=cutVal], collapse=";")
      } else {
        o <- paste("std", as.character(inv[prStds>=(1-cutVal) & prStds<=cutVal]), sep=";")
      }
      return(o)
    }
    kary.ag <- kary[!is.nan(prStd), list(chr2L=karyCall(prStds=prStd[chr=="2L"], inv=inversion[chr=="2L"]),
                                         chr2R=karyCall(prStds=prStd[chr=="2R"], inv=inversion[chr=="2R"]),
                                         chr3L=karyCall(prStds=prStd[chr=="3L"], inv=inversion[chr=="3L"]),
                                         chr3R=karyCall(prStds=prStd[chr=="3R"], inv=inversion[chr=="3R"])),
                    list(sample.id)]
  ## write out the file
    write.table(kary.ag, "/scratch/ab5dr/wildDmel2016/redoAnalysis/wild2016_v6_invStatus.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)

    
##########################################
### VALIDATE INVERSION TYPING WITH PCA ###
##########################################
    
### read in the data 
  genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_seqA.gds")
  all.snp <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/LDpruned_snpset.txt", header=T)
  all.snp <- all.snp[chromosome=="2L"]
  
### filter down to SNPset   
  gt.ag <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/nHets_byChr_v4.txt", header=T)
  gt.ag <- gt.ag[! sampID %in% c("CM.162.0902","CM.024.0708","CM.039.0722")]
  good.samps <- gt.ag[chrom=="X"][Nhets<1000][,c("sampID")]
  good.samps[,pop:=tstrsplit(sampID,"[.]")[3]]

### get input information into better format 
  ## subset to LD pruned SNPS & good samples
    seqSetFilterPos(genofile, all.snp$chromosome, all.snp$position)
    seqSetFilter(genofile, sample.id=good.samps$sampID, action="intersect")
  ## get out genotype information per sample 
    gt.dt <- t(seqGetData(genofile, "$dosage"))
    gt.dt <- as.data.table(gt.dt)
    colnames(gt.dt) <- seqGetData(genofile, "sample.id")
    gt.dt[,snp.id:=seqGetData(genofile, "annotation/id")]
    gt.dt <- melt(gt.dt, measure.vars=seqGetData(genofile, "sample.id"), value.name="RefDosage",variable.name="sampID")
  
### run PCA on 2L SNPs 
  ## prep 
    samp.ids <- unique(gt.dt$sampID)
  ## do PCA
    # reformat the data 
      gt.dt <- gt.dt[,c("sampID","snp.id","RefDosage")]
      gt.dt <- dcast(gt.dt, snp.id ~ sampID, value.var="RefDosage")
      gt.dt[is.na(gt.dt)] <- 9
      gt.dt <- gt.dt[, ..samp.ids]
      gt.dt <- as.matrix(gt.dt)
      gt.dt <- t(gt.dt)
    # get LEA data format 
      write.lfmm(gt.dt, "/scratch/ab5dr/wildDmel2016/redoAnalysis/pca/invTyping_2L_lea.lfmm")
    # run PCA 
      pc = pca("/scratch/ab5dr/wildDmel2016/redoAnalysis/pca/invTyping_2L_lea.lfmm")
    # get pcs into output data format 
      pc.dt <- as.data.table(pc$projections)
      setnames(pc.dt, paste("PC",1:ncol(pc.dt), sep=""))
      pc.dt[,rn:= "pca"]
    # get PVE 
      pve.dt <- as.data.table(summary(pc))
      pve.dt[,rn:=c("stdev","pve","cumulativePVE")]
      pve.dt <- pve.dt[rn=="pve"]
    # get final output table 
      pc.dt <- rbindlist(list(pc.dt, pve.dt))
      pc.dt[rn=="pca",sampID:=samp.ids]

### write out PC loading values and PVE
  write.table(pc.dt, "/scratch/ab5dr/wildDmel2016/redoAnalysis/pca/invTyping_2L_lea_PCA.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)
    
      
######################
### MAKE FIGURE 02 ###
######################
  
### merge info from PCA and inversion probability-based calls 
  ## read in the files generated above
    pc.out <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/pca/invTyping_2L_lea_PCA.txt", header=T)
    inv.dt <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/wild2016_v6_invStatus.txt", header=T)
  ## adjust inversion typing for 2L
    inv.dt[,chr2L:=ifelse(chr2L=="std;In(2L)t", "heterozygous",
                          ifelse(chr2L=="std","homozygous standard","homozygous inverted"))]
    inv.dt <- inv.dt[,c("sample.id","chr2L")]
    setnames(inv.dt, old="sample.id",new="sampID")
  ## adjust format and subset down pca table
    pc.out <- pc.out[,c("PC1","PC2","PC3","PC4","sampID","rn")]
  ## merge together into a plot-able format
    setkey(inv.dt, sampID)
    setkey(pc.out, sampID)
    inv.dt <- merge(inv.dt, pc.out)
    
    
### FIGURE 02
  cb_palette <- c("#440154FF","#20A387FF","#95D840FF") # teal purple green
  ggplot(inv.dt, aes(x=PC1, y=PC2, color=chr2L)) + geom_point() + labs(x="PC1 (pve=9.24%)", y="PC2 (pve=2.22%)", color="Inv(2L)t status") + 
    theme_classic() + scale_color_manual(values=cb_palette)
    
    