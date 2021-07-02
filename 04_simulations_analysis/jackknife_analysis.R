### PERFORM JACKKNIFE ANALYSIS OF FIS ESTIMATES


### libraries
library(data.table)
library(foreach)
library(GWASExactHW)
library(gdsfmt)
library(SeqArray)


### get a jackknifed per-SNP FIS estimates on subset of simulations
  ## some prep
    gt.ag <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/nHets_byChr_v4.txt", header=T)
    gt.ag <- gt.ag[! sampID %in% c("CM.162.0902","CM.024.0708","CM.039.0722")]
    lowXhet.samps <- gt.ag[chrom=="X"][Nhets<1000][,c("sampID")]
    lowXhet.samps[,pop:=tstrsplit(sampID,"[.]")[3]]
    # selecting the simulation sample now so it's the same across everything
      sim.samp <- sample(1:1000, 10)    
        ### 842 911 195 437 663 245 859 960 285 891  ### simulations used in this manuscript
  ## loop
    for(i in c("lowXhet")){
      for(j in c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203")) {
        dir <- "byPop"
        sim.out <- foreach(k=c("2L","2R","3L","3R")) %do% {
          print(paste("filenamesStep_outerLoop",i,j,k, sep="_"))
          filenames=list.files(path=paste("/scratch/ab5dr/wildDmel2016/simulations/rawOut/",dir,sep=""), pattern=paste("raw",j,i,k, sep="_"))
          foreach(m=filenames, .combine="rbind")%do%{
            ## read in the raw
              print("next file")
              sim.temp <- readRDS(paste("/scratch/ab5dr/wildDmel2016/simulations/rawOut",dir,m, sep="/"))
              sim.temp <- sim.temp[,c("snp.id", "sampID","simNo","ref.rd", "alt.rd","TsimGT", "Ocall")]
            ## subset to 5 simulations total and get samps picked 
              sim.temp <- sim.temp[simNo %in% sim.samp]
              samps <- lowXhet.samps[pop==j]$sampID
            ## jackknife loop
              foreach(o=samps, .combine="rbind")%do%{
                print(o)
                  # subset 
                    temp.gt <- sim.temp[sampID!=o]
                  # do HWE calcs
                    # aggregate by snp
                      temp.ag <- temp.gt[,list(nRR=sum(Ocall=="RR", na.rm=T),
                                               nRA=sum(Ocall=="RA", na.rm=T),
                                               nAA=sum(Ocall=="AA", na.rm=T)),
                                         by=list(snp.id, simNo)]
                      temp.ag[,nT:=(nRR+nRA+nAA)]
                    # subset info to feed into HWE package 
                      hwe.temp <- temp.ag[,c("nRR","nRA","nAA")]
                      setnames(hwe.temp, old=c("nRR","nRA","nAA"), new=c("nAA","nAa","naa"))
                    # test for HWE deviations
                      HWEout.dt <- data.table(snp.id=temp.ag$snp.id, pval=HWExact(hwe.temp), simNo=temp.ag$simNo)
                  # do F calcs
                    # get p and q
                      temp.ag[,p:=2*(nRR/(2*nT)) + (nRA/(2*nT))]
                      temp.ag[,q:=2*(nAA/(2*nT)) + (nRA/(2*nT))]
                    # get f
                      temp.ag[,Ho:=(nRA/nT)]
                      temp.ag[,He:=(nT/(nT-1)) * ((2*p*q) - (Ho/(2*nT)))]
                      temp.ag[,neiFis:=(1-(Ho/He))]
                  # get final table 
                    temp.ag <- temp.ag[,c("snp.id","simNo","nT","p","q","neiFis")]
                    setkey(temp.ag, snp.id, simNo)
                    setkey(HWEout.dt, snp.id, simNo)
                    out.dt <- merge(temp.ag, HWEout.dt)
                    out.dt[,class:=paste(j,i,sep="_")]
                    out.dt[,jackDrop:=o]
                  # on to the next 
                    out.dt
              }
          }
        } 
        sim.out <- rbindlist(sim.out)
        print("write output table")
        write.table(sim.out, paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/simSubset_byPop_",j,"_f_HWE_lowXhet_jackknife.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
      }
    }


### get a jackknifed per-SNP FIS estimates on empirical data
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
  ## subset genofile to lowXhet samps and the LD pruned snps 
    all.snp <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/LDpruned_snpset.txt", header=T)
    seqSetFilterPos(genofile, all.snp$chromosome, all.snp$position)
    seqSetFilter(genofile, sample.id=lowXhet.samps$sampID, action="intersect")
  ## pull out dosage 
    gt.dt <- t(seqGetData(genofile, "$dosage"))
    gt.dt <- as.data.table(gt.dt)
    setnames(gt.dt, lowXhet.samps$sampID)
    # continue
      gt.dt[,snp.id:=seqGetData(genofile, "annotation/id")]
      gt.dt <- melt(gt.dt, measure.vars=lowXhet.samps$sampID, value.name="RefDosage",variable.name="sampID")
      gt.dt[,pop:=tstrsplit(sampID, "[.]")[3]]
  ## loop 
    for(i in pops) {
      print(i)
      samps <- lowXhet.samps[pop==i]$sampID
      gt.jack <- foreach(j=samps, .combine="rbind")%do%{
        print(j)
        # subset 
          temp.gt <- gt.dt[pop==i][sampID!=j]
        # get genotype counts for HWE
          gt.pop <- temp.gt[,list(nRR=sum(RefDosage==2, na.rm=T),
                                  nRA=sum(RefDosage==1, na.rm=T),
                                  nAA=sum(RefDosage==0, na.rm=T)),
                            by=list(snp.id, pop)]
          gt.pop[,nT:=nRR+nRA+nAA]
        # run HWE test
          temp.dt <- gt.pop[pop==i][,c("nRR","nRA","nAA")]
          setnames(temp.dt, old=c("nRR","nRA","nAA"), new=c("nAA","nAa","naa"))
          HWE.temp <- data.table(snp.id=gt.pop[pop==i]$snp.id, pop=i, HWEpval=HWExact(temp.dt))
        # calculate F
          gt.pop[,p:=2*(nRR/(2*nT)) + (nRA/(2*nT))]
          gt.pop[,q:=2*(nAA/(2*nT)) + (nRA/(2*nT))]
          gt.pop[,Ho:=(nRA/nT)]
          gt.pop[,He:=(nT/(nT-1)) * ((2*p*q) - (Ho/(2*nT)))]
          gt.pop[,neiFis:=(1-(Ho/He))]
        # output table
          setkey(gt.pop, snp.id, pop)
          setkey(HWE.temp, snp.id, pop)
          gt.pop <- merge(gt.pop, HWE.temp)
          gt.pop[,jackDrop:=j]
          gt.pop
      }
      write.table(gt.jack, paste("/scratch/ab5dr/wildDmel2016/redoAnalysis/hweAndF/emp_byPop_",i,"_neiFis_HWE_lowXhet_jackknife.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
    }


### get jackknife FIS averages 
  ## simulated
    for(i in c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203")) {
      print(i)
      temp.dt <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/simSubset_byPop_",i,"_NeiFis_HWE_lowXhet_jackknife.txt",sep=""), header=T)
      print("analysis bit")
      sim.f <- temp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T)), by=list(class, jackDrop, simNo)]
      rm(temp.dt)
      write.table(sim.f, paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_sim_jackknife_",i,"_lowXhet.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
    }
  ## empirical
    for(i in c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203")) {
      print(i)
      temp.dt <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/emp_byPop_",i,"_neiFis_HWE_lowXhet_jackknife.txt",sep=""), header=T,
                       colClasses=c("character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character"))
      print("analysis bit")
      sim.f <- temp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T)), by=list(pop, jackDrop)]
      rm(temp.dt)
      write.table(sim.f, paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_emp_jackknife_",i,"_lowXhet.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
    }

      
      
      