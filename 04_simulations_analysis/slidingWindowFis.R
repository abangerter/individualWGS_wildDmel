#! /usr/bin/env Rscript


### libraries 
library(methods)
library(data.table)
library(foreach)


### bring in args 
  args <- commandArgs(trailingOnly=TRUE)
  i <- args[1]
  j <- args[2]
  print(paste(i,j,sep="_"))


### read in data and windows
  ## read in the  data
    temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/simSumm/neiFis_",i,"_lowXhet_",j,".txt",sep=""))
    temp[,c("chrom","pos"):=tstrsplit(snp.id,"_")]
    temp[,pos:=as.numeric(pos)]
    setkey(temp, pos, simNo)
  ## read in position array
    window.dt <- fread(paste(j,"posArray.txt",sep="_"), header=T)
    windows <- nrow(window.dt)
    
    
### loop
print("start loop")
  slide.dt <- foreach(k=1:windows, .combine="rbind") %do%{
    start <- window.dt[k]$start
    stop <-  window.dt[k]$stop
    # subset down to window for the sliding window
      temp.sub <- temp[J(start:stop)]
    # get fhat
      temp.out <- temp.sub[,list(fhatmu=mean(neiFis, na.rm=T),
                                 fhatmed=median(neiFis, na.rm=T),
                                 NposF=sum(neiFis>0, na.rm=T),
                                 NnegF=sum(neiFis<0, na.rm=T),
                                 NzeroF=sum(neiFis==0, na.rm=T),
                                 nT=sum(neiFis<5, na.rm=T),
                                 nWindow=.N,
                                 minPos=min(pos),
                                 maxPos=max(pos)),
                           by=simNo]
    # prep output table
      temp.out[,chrom:=j]
      temp.out[,pop:=i]
    # out
      temp.out
  }

  
### write output 
write.table(slide.dt, paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/slidingWindow_v2_NeiFis_",i,"_",j,"_lowXhet_1000.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)

