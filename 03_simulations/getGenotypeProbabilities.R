### GET GENOTYPE PROBABILITIES ###


### libraries 
library(data.table)
library(foreach)
library(gdsfmt)
library(SeqArray)


### open genofile
  genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_seqA.gds")

### filter out snps from X, 4 and Y
  snp.id <- seqGetData(genofile, "annotation/id")
  snps.dt <- as.data.table(snp.id)
  snps.dt[,c("chr","pos"):=tstrsplit(snp.id,"_")]
  snps.dt <- snps.dt[chr %in% c("2L","2R","3L","3R")]
  seqSetFilterPos(genofile, chr=snps.dt$chr, pos=snps.dt$pos)

### calculate the probabilities based on allele depths
  ## get allele depth info
    geno <- expand.grid(seqGetData(genofile, var.name="$dosage", .useraw=FALSE))[,1]
    sample.id <- rep(seqGetData(genofile, var.name="sample.id", .useraw=FALSE), length(snps.dt$pos))									
    ad <- seqGetData(genofile, var.name="annotation/format/AD", .useraw=FALSE)
    ad.ref = expand.grid(as.matrix(ad$data[,c(1:(2*length(snps.dt$pos)))%%2==1], ncol=length(snps.dt$pos)))[,1]
    ad.alt = expand.grid(as.matrix(ad$data[,c(1:(2*length(snps.dt$pos)))%%2==0], ncol=length(snps.dt$pos)))[,1]
    snp.id <- rep(snps.dt$snp.id, each=length(seqGetData(gdsfile=genofile, var.name="sample.id", .useraw=FALSE)))
    gt.dt <- data.table(snp.id=snp.id, sample.id=sample.id, GT=geno, ref.rd=ad.ref, alt.rd=ad.alt)
  ## get the probabilities
    out.dt <- gt.dt[,list(nRR=sum(GT==2, na.rm=T),
                          nRA=sum(GT==1,na.rm=T),
                          nAA=sum(GT==0,na.rm=T)),
                    by=list(ref.rd,alt.rd)]
    setkey(out.dt, ref.rd, alt.rd)
    out.dt[,nT:=nRR+nRA+nAA]
    out.dt[,GTprob2:=nRR/nT]
    out.dt[,GTprob1:=nRA/nT]
    out.dt[,GTprob0:=nAA/nT]

### write out the table
  write.table(out.dt, "/scratch/ab5dr/wildDmel2016/simulations/misc/gtProbabilities.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)
