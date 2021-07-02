#! /usr/bin/env Rscript
### SIMULATIONS TO GENERATE A NULL DISTRIBUTION ###
### RUNS 1000 SIMULATIONS PER BATCH OF SNPS AND SAMPLES ###


print("start")
### libraries 
  library(methods)
  library(data.table)
  library(foreach)


### bring in the args
  args <- commandArgs(trailingOnly=TRUE)
  # chrom
    m <- args[1]
  # code
    p <- args[2]
  # start & stop 
    start <- args[3]
    start <- as.numeric(as.character(start))
    stop <- args[4]
    stop <- as.numeric(as.character(stop))
    print(paste(m,p,start,stop, sep="_"))
  
  
### read in pre-fabricated tables of empirical data 
  print("reading in tables")
  ## allele frequencies
    af.dt <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/rafTables/raf_",p,"_",m,".txt",sep=""), header=T)
    af.dt <- af.dt[pos >= start][pos <= stop]
    af.dt <- af.dt[refAF!=1][refAF!=0] # remove monomorphic sites
    af.dt <- af.dt[,c("snp.id", "refAF")]
  ## site read depths
    rd.dt <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/RDinfo/RDinfo_",p,"_",m,".txt",sep=""), header=T)
    rd.dt <- rd.dt[pos >= start][pos <= stop]
    rd.dt <- rd.dt[,c("snp.id","sampID","RD")]
  ## gatk probability table 
    gatkProbs.dt <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/gtProbabilities.txt", header=T)
    gatkProbs.dt <- gatkProbs.dt[,c("ref.rd","alt.rd","GTprob2","GTprob1","GTprob0")]
    gatkProbs.dt <- na.omit(gatkProbs.dt)
    setkey(gatkProbs.dt, ref.rd, alt.rd)
  

### prep simulation dt
  setkey(af.dt, snp.id)
  setkey(rd.dt, snp.id)
  sim.dt <- merge(af.dt, rd.dt)
    

### simulation loop
  sim.OUT <- foreach(j=1:1000, .combine="rbind")%do%{
    print(paste("simulation run",j,sep=" "))
    sim.dt[,idx := 1:.N]
    ## simulate a true genotype
      sim.dt[,TsimGT:=rbinom(1,2,refAF), by=idx]
    ## simulate reference allele read depth
      sim.dt[TsimGT==1, ref.rd:=rbinom(1,RD,0.5), by=idx]
      sim.dt[TsimGT==2, ref.rd:=RD]
      sim.dt[TsimGT==0, ref.rd:=0]
    ## calculate simulated alt allele read depth
      sim.dt[,alt.rd:= RD-ref.rd]
    ## simulate the observed genotype
      setkey(sim.dt, ref.rd, alt.rd)
      setkey(gatkProbs.dt, ref.rd, alt.rd)
      sim.temp <- merge(sim.dt, gatkProbs.dt)
      sim.temp[,Ocall:= sample(c("RR","RA","AA"), prob=c(GTprob2, GTprob1, GTprob0), size=1), by=idx]
      sim.temp[ref.rd==0 & alt.rd==0, Ocall:= NA]
    ## prep output
      sim.temp[,simNo:=j]
      sim.temp <- sim.temp[,c("snp.id", "sampID","simNo","ref.rd", "alt.rd","TsimGT", "Ocall")]
    ## print and on to the next round
      sim.temp
  }

  
### write simulations out 
  print("writing raw rds")
  saveRDS(sim.OUT, file = paste("/scratch/ab5dr/wildDmel2016/simulations/rawOut/raw_",p,"_",m,"_",start,"_",stop,".rds",sep=""))

