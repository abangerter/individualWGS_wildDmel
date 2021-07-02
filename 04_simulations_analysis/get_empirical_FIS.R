### GET EMPIRICAL FIS PER SNP ###


### libraries
library(data.table)
library(foreach)
library(GWASExactHW)
library(gdsfmt)
library(SeqArray)


### do the full sample HWE & FIS calculations 
  ## read in data
    genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_seqA.gds")
    samps <- seqGetData(genofile, "sample.id")
    sam.pop <- data.table(sampID=samps)
    sam.pop[,pop:= tstrsplit(sampID, "[.]")[3]]
    pops <- as.vector(unique(sam.pop$pop))
  ## grab the list of samps to include
    gt.ag <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/nHets_byChr_v4.txt", header=T)
    gt.ag <- gt.ag[! sampID %in% c("CM.162.0902","CM.024.0708","CM.039.0722")]
    lowXhet.samps <- gt.ag[chrom=="X"][Nhets<1000][,c("sampID")]
    lowXhet.samps[,pop:=tstrsplit(sampID,"[.]")[3]]
  ## subset genofile to lowXhet samps
    seqSetFilter(genofile, sample.id=lowXhet.samps$sampID)
  ## pull out dosage 
    gt.dt <- t(seqGetData(genofile, "$dosage"))
    gt.dt <- as.data.table(gt.dt)
    setnames(gt.dt, lowXhet.samps$sampID)
    # continue
      gt.dt[,snp.id:=seqGetData(genofile, "annotation/id")]
      gt.dt <- melt(gt.dt, measure.vars=lowXhet.samps$sampID, value.name="RefDosage",variable.name="sampID")
  ## aggregate into GT counts for HWE and FIS 
    gt.fs <- gt.dt[,list(nRR=sum(RefDosage==2, na.rm=T),
                         nRA=sum(RefDosage==1, na.rm=T),
                         nAA=sum(RefDosage==0, na.rm=T)),
                   by=list(snp.id)]
    gt.fs[,nT:=nRR+nRA+nAA]
  ## pull out only the autosomes
    gt.fs[,chrom:=tstrsplit(snp.id,"_")[1]]
    gt.fs <- gt.fs[chrom %in% c("2L","2R","3L","3R")]
  ## run HWE test
    temp.dt <- gt.fs[,c("nRR","nRA","nAA")]
    setnames(temp.dt, old=c("nRR","nRA","nAA"), new=c("nAA","nAa","naa"))
    HWE.fs <- data.table(snp.id=gt.fs$snp.id, HWEpval=HWExact(temp.dt))
  ## calculate F
    gt.fs[,p:=2*(nRR/(2*nT)) + (nRA/(2*nT))]
    gt.fs[,q:=2*(nAA/(2*nT)) + (nRA/(2*nT))]
    gt.fs[,Ho:=(nRA/nT)]
    gt.fs[,He:=(nT/(nT-1)) * ((2*p*q) - (Ho/(2*nT)))]
    gt.fs[,neiFis:=(1-(Ho/He))]
  ## output table
    setkey(gt.fs, snp.id)
    setkey(HWE.fs, snp.id)
    gt.fs <- merge(gt.fs, HWE.fs)
    write.table(gt.fs, "/scratch/ab5dr/wildDmel2016/redoAnalysis/hweAndF/emp_fs_f_HWE_lowXhet_neiFis.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)


### do the per time-point HWE & FIS calculations 
  ## use data table generated in above code section
  ## aggregate into GT counts for HWE and FIS 
    gt.dt[,pop:=tstrsplit(sampID, "[.]")[3]]
    gt.pop <- gt.dt[,list(nRR=sum(RefDosage==2, na.rm=T),
                          nRA=sum(RefDosage==1, na.rm=T),
                          nAA=sum(RefDosage==0, na.rm=T)),
                    by=list(snp.id, pop)]
    gt.pop[,nT:=nRR+nRA+nAA]
  ## pull out only the autosomes
    gt.pop[,chrom:=tstrsplit(snp.id,"_")[1]]
    gt.pop <- gt.pop[chrom %in% c("2L","2R","3L","3R")]
  ## run HWE test
    HWE.pop <- foreach(i=pops)%do%{
      temp.dt <- gt.pop[pop==i][,c("nRR","nRA","nAA")]
      setnames(temp.dt, old=c("nRR","nRA","nAA"), new=c("nAA","nAa","naa"))
      HWE.temp <- data.table(snp.id=gt.pop[pop==i]$snp.id, pop=i, HWEpval=HWExact(temp.dt))
      HWE.temp
    }
    HWE.pop <- rbindlist(HWE.pop)
  ## calculate F
    gt.pop[,p:=2*(nRR/(2*nT)) + (nRA/(2*nT))]
    gt.pop[,q:=2*(nAA/(2*nT)) + (nRA/(2*nT))]
    gt.pop[,Ho:=(nRA/nT)]
    gt.pop[,He:=(nT/(nT-1)) * ((2*p*q) - (Ho/(2*nT)))]
    gt.pop[,neiFis:=(1-(Ho/He))]
  ## output table
    setkey(gt.pop, snp.id, pop)
    setkey(HWE.pop, snp.id, pop)
    gt.pop <- merge(gt.pop, HWE.pop)
    write.table(gt.pop, "/scratch/ab5dr/wildDmel2016/redoAnalysis/hweAndF/emp_byPop_f_HWE_lowXhet_neiFis.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)
