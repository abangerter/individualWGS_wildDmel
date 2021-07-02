### GATHER SIMULATION OUTPUT & CALCULATE FIS ###

### libraries
library(foreach)
library(data.table)
library(GWASExactHW)


### gather information and calculate per-SNP FIS with Nei's estimator for the OBSERVED genotypes of the full sample and for each time-point
  for(i in c("lowXhet")){
    for(j in c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs")) {
      if (j=="fs") { dir <- "fs" } else { dir <- "byPop" }
      for(k in c("2L","2R","3L","3R")) {
        print(paste("filenamesStep_outerLoop",i,j,k, sep="_"))
        filenames=list.files(path=paste("/scratch/ab5dr/wildDmel2016/simulations/rawOut/",dir,sep=""), pattern=paste("raw",j,i,k, sep="_"))
        sim.out <- foreach(m=filenames, .combine="rbind")%do%{
          ## read in the raw
            print("readTables")
            sim.temp <- readRDS(paste("/scratch/ab5dr/wildDmel2016/simulations/rawOut",dir,m, sep="/"))
            sim.temp <- sim.temp[,c("snp.id", "sampID","simNo","ref.rd", "alt.rd","TsimGT", "Ocall")]
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
        write.table(sim.out, paste("/scratch/ab5dr/wildDmel2016/simulations/simSumm/neiFis_",j,"_",i,"_",k,".txt", sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
      }
    }
  }


### gather information and calculate per-SNP FIS with Nei's estimator for the TRUE genotypes of the full sample and for each time-point
  for(i in c("lowXhet")){
    for(j in c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs")) {
      if (j=="fs") { dir <- "fs" } else { dir <- "byPop" }
      for(k in c("2L","2R","3L","3R")) {
        print(paste("filenamesStep_outerLoop",i,j,k, sep="_"))
        filenames=list.files(path=paste("/scratch/ab5dr/wildDmel2016/simulations/rawOut/",dir,sep=""), pattern=paste("raw",j,i,k, sep="_"))
        sim.out <- foreach(m=filenames, .combine="rbind")%do%{
          ## read in the raw
            print("readTables")
            sim.temp <- readRDS(paste("/scratch/ab5dr/wildDmel2016/simulations/rawOut",dir,m, sep="/"))
            sim.temp <- sim.temp[,c("snp.id", "sampID","simNo","TsimGT")]
          ## do HWE calcs
            print("HWEandF")
            # aggregate by snp
              temp.ag <- sim.temp[,list(nRR=sum(TsimGT==2, na.rm=T),
                                        nRA=sum(TsimGT==1, na.rm=T),
                                        nAA=sum(TsimGT==0, na.rm=T)),
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
        write.table(sim.out, paste("/scratch/ab5dr/wildDmel2016/simulations/trueSumm/neiFis_TsimGT_",j,"_",i,"_",k,".txt", sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
      }
    }
  }

