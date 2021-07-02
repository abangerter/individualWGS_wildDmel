### RECALCULATE FIS ON POPULATIONS AFTER REMOVAL OF Inv(2L)t HOMOZYGOTES ###


### libraries
library(data.table)
library(foreach)
library(GWASExactHW)
library(gdsfmt)
library(SeqArray)


#######################
### FOR SIMULATIONS ###
#######################

### per time-point and full sample FIS after removal of 3 Inv(2L)t homozygotes for SIMULATONS
  for(i in c("lowXhet")){
    for(j in c("0722","0819","fs")) {
      if (j=="fs") { dir <- "fs" } else { dir <- "byPop" }
      for(k in c("2L","2R","3L","3R")) {
        print(paste("filenamesStep_outerLoop",i,j,k, sep="_"))
        filenames=list.files(path=paste("/project/berglandlab/alyssa/simulations/",dir,sep=""), pattern=paste("raw",j,i,k, sep="_"))
        sim.out <- foreach(m=filenames, .combine="rbind")%do%{
          ## read in the raw
            print("readTables")
            sim.temp <- readRDS(paste("/project/berglandlab/alyssa/simulations",dir,m, sep="/"))
            sim.temp <- sim.temp[,c("snp.id", "sampID","simNo","ref.rd", "alt.rd","TsimGT", "Ocall")]
            sim.temp <- sim.temp[! sampID %in% c("CM.025.0722", "CM.017.0819","CM.029.0819")]
          ## do HWE calcs
            print("HWEandF")
            # aggregate by snp
              temp.ag <- sim.temp[,list(nRR=sum(Ocall=="RR", na.rm=T),
                                        nRA=sum(Ocall=="RA", na.rm=T),
                                        nAA=sum(Ocall=="AA", na.rm=T)),
                                  by=list(snp.id, simNo)]
              temp.ag[,nT:=(nRR+nRA+nAA)]
            # subset info to feed into HWE package 
              hwe.temp <- temp.ag[,c("nRR","nRA","nAA")]
              setnames(hwe.temp, old=c("nRR","nRA","nAA"), new=c("nAA","nAa","naa"))
            # test for HWE deviations
              HWEout.dt <- data.table(snp.id=temp.ag$snp.id, pval=HWExact(hwe.temp), simNo=temp.ag$simNo)
          ## do F calcs
            # get p and q
              temp.ag[,p:=2*(nRR/(2*nT)) + (nRA/(2*nT))]
              temp.ag[,q:=2*(nAA/(2*nT)) + (nRA/(2*nT))]
            # get f
              temp.ag[,Ho:=(nRA/nT)]
              temp.ag[,He:=(nT/(nT-1)) * ((2*p*q) - (Ho/(2*nT)))]
              temp.ag[,neiFis:=(1-(Ho/He))]
          ## get final table 
            temp.ag <- temp.ag[,c("snp.id","simNo","nT","p","q","neiFis")]
            setkey(temp.ag, snp.id, simNo)
            setkey(HWEout.dt, snp.id, simNo)
            out.dt <- merge(temp.ag, HWEout.dt)
            out.dt[,class:=paste(j,i,sep="_")]
          ## on to the next 
          out.dt
        }
        print("write output table")
        write.table(sim.out, paste("/scratch/ab5dr/wildDmel2016/simulations/simSumm/neiFis_no2LtHoms_",j,"_",i,"_",k,".txt", sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
      }
    }
  }


### Get average FIS for the full sample and per time-point for SIMULATIONS
  for(i in c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs")) {
    print(i)
    temp.dt <- foreach(j=c("2L","2R","3L","3R"))%do%{
      print(j)
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/simSumm/neiFis_no2LtHoms_",i,"_lowXhet_",j,".txt",sep=""), header=T)
      temp
    }
    temp.dt <- rbindlist(temp.dt)
    print("analysis bit")
    sim.f <- temp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T), stdev=sd(neiFis, na.rm=T)), by=list(class, simNo)]
    temp.dt[,chrom:=tstrsplit(snp.id, "_")[1]]
    sim.f.chrom <- temp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T), stdev=sd(neiFis, na.rm=T)), by=list(class, simNo, chrom)]
    sim.f[,chrom:="WG"]
    sim.f <- sim.f[,c("class","simNo","chrom","fhatmu","fhatmed","stdev")]
    sim.f <- rbindlist(list(sim.f, sim.f.chrom))
    rm(temp.dt)
    write.table(sim.f, paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_plusStdev_",i,"_lowXhet_no2LtHoms.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
  }


##########################
### FOR EMPIRICAL DATA ###
##########################

### full sample FIS after removal of 3 Inv(2L)t homozygotes for EMPIRICAL DATA
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
      gt.dt <- gt.dt[! sampID %in% c("CM.025.0722", "CM.017.0819","CM.029.0819")]
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
    write.table(gt.fs, "/scratch/ab5dr/wildDmel2016/simulations/analysis/emp_fs_f_HWE_lowXhet_no2LtHoms_neiFis.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)

    
### per time-point FIS after removal of 3 Inv(2L)t homozygotes for EMPIRICAL DATA
  ## uses gt.dt generated in above code block
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
    write.table(gt.pop, "/scratch/ab5dr/wildDmel2016/redoAnalysis/hweAndF/emp_byPop_f_HWE_lowXhet_no2LtHoms_neiFis.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)

    
### get average FIS estimates
  ## read in the data 
    emp.fs.dt <- fread("/scratch/ab5dr/wildDmel2016/simulations/analysis/emp_fs_f_HWE_lowXhet_no2LtHoms_neiFis.txt", header=T)
    emp.fs.dt[,pop:="fs"]
    emp.fs.dt <- emp.fs.dt[,c("snp.id","pop","nRR","nRA","nAA","nT","chrom","p","q","neiFis","HWEpval")]
    emp.pop.dt <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/hweAndF/emp_byPop_f_HWE_lowXhet_no2LtHoms_neiFis.txt", header=T,
                        colClasses=c("character","character","numeric","numeric","numeric","numeric","character","numeric","numeric","numeric","numeric","numeric","numeric"))
    emp.pop.dt <- emp.pop.dt[,c("snp.id","pop","nRR","nRA","nAA","nT","chrom","p","q","neiFis","HWEpval")]
    emp.dt <- rbindlist(list(emp.fs.dt, emp.pop.dt))
  ## subset empirical data to the LD-pruned snpset 
    all.snp <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/LDpruned_snpset.txt", header=T)
    all.snp[,snp.id:=paste(chromosome, position, sep="_")]
    all.snp <- all.snp[,c("snp.id")]
    setkey(all.snp, snp.id)
    setkey(emp.dt, snp.id)
    emp.dt <- merge(all.snp, emp.dt)
  ## generate emp FIS averages
    emp.f <- emp.dt[,list(fhatmu=mean(neiFis, na.rm=T),
                          fhatmed=median(neiFis, na.rm=T),
                          stdev=sd(neiFis, na.rm=T)), 
                    by=list(pop)]
    emp.f.chrom <- emp.dt[,list(fhatmu=mean(neiFis, na.rm=T), 
                                fhatmed=median(neiFis, na.rm=T), 
                                stdev=sd(neiFis, na.rm=T)), 
                          by=list(pop, chrom)]
    emp.f[,chrom:="WG"]
    emp.f <- emp.f[,c("pop","chrom","fhatmu","fhatmed","stdev")]
    emp.f <- rbindlist(list(emp.f, emp.f.chrom))
    emp.f[,simNo:="NA"]
    emp.f[,type:="empirical"]
    emp.f <- emp.f[,c("pop","chrom","type","simNo","fhatmu","fhatmed","stdev")]
    write.table(emp.f, "/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_emp_lowXhet_no2LtHoms.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)

