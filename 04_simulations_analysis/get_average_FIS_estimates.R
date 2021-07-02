### CALCULATE FIS AVERAGES ###


### libraries
library(data.table)
library(foreach)



### full sample and per time-point on OBSERVED genotype FIS estimates
  for(i in c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs")) {
    print(i)
    temp.dt <- foreach(j=c("2L","2R","3L","3R"))%do%{
      print(j)
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/simSumm/neiFis_",i,"_lowXhet_",j,".txt",sep=""), header=T)
      temp
    }
    temp.dt <- rbindlist(temp.dt)
    print("analysis bit")
    sim.f <- temp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T)), by=list(class, simNo)]
    temp.dt[,chrom:=tstrsplit(snp.id, "_")[1]]
    sim.f.chrom <- temp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T)), by=list(class, simNo, chrom)]
    sim.f[,chrom:="WG"]
    sim.f <- sim.f[,c("class","simNo","chrom","fhatmu","fhatmed")]
    sim.f <- rbindlist(list(sim.f, sim.f.chrom))
    rm(temp.dt)
    write.table(sim.f, paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_",i,"_all.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
  }


### full sample and per time-point on TRUE genotype FIS estimates
  for(i in c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs")) {
    print(i)
    temp.dt <- foreach(j=c("2L","2R","3L","3R"))%do%{
      print(j)
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/trueSumm/neiFis_TsimGT_",i,"_lowXhet_",j,".txt",sep=""), header=T)
      temp
    }
    temp.dt <- rbindlist(temp.dt)
    print("analysis bit")
    sim.f <- temp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T)), by=list(class, simNo)]
    temp.dt[,chrom:=tstrsplit(snp.id, "_")[1]]
    sim.f.chrom <- temp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T)), by=list(class, simNo, chrom)]
    sim.f[,chrom:="WG"]
    sim.f <- sim.f[,c("class","simNo","chrom","fhatmu","fhatmed")]
    sim.f <- rbindlist(list(sim.f, sim.f.chrom))
    rm(temp.dt)
    write.table(sim.f, paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_TsimGT_",i,"_lowXhet.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
  }
 
 
### average FIS on actual data, not simulations  
  emp.fs.dt <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/hweAndF/emp_fs_f_HWE_lowXhet_neiFis.txt", header=T)
  emp.fs.dt[,pop:="fs"]
  emp.fs.dt <- emp.fs.dt[,c("snp.id","pop","nRR","nRA","nAA","nT","chrom","p","q","neiFis","HWEpval")]
  emp.pop.dt <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/hweAndF/emp_byPop_f_HWE_lowXhet_neiFis.txt", header=T,
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
  ## generate emp fhats 
    emp.f <- emp.dt[,list(fhatmu=mean(neiFis, na.rm=T),
                          fhatmed=median(neiFis, na.rm=T),
                          stdev=sd(neiFis, na.rm=T)), 
                    by=list(pop)]
    emp.f.chrom <- emp.dt[,list(fhatmu=mean(neiFis, na.rm=T), fhatmed=median(neiFis, na.rm=T), stdev=sd(neiFis, na.rm=T)), by=list(pop, chrom)]
    emp.f[,chrom:="WG"]
    emp.f <- emp.f[,c("pop","chrom","fhatmu","fhatmed","stdev")]
    emp.f <- rbindlist(list(emp.f, emp.f.chrom))
    emp.f[,simNo:="NA"]
    emp.f[,type:="empirical"]
    emp.f <- emp.f[,c("pop","chrom","type","simNo","fhatmu","fhatmed","stdev")]
    write.table(emp.f, "/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_emp_lowXhet.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)


