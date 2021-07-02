### FIGURE 7 ###


### libraries
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)


### read in the data for after Inv(2L)t homozygote removal
  ## simulations 
    sim.f <- foreach(i=c("0722","0819","fs"), .combine="rbind") %do% {
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_plusStdev_",i,"_lowXhet_no2LtHoms.txt",sep=""))
      temp
    }
    sim.f[,pop:=tstrsplit(class,"_")[1]]
    setnames(sim.f, old="fhatmu", new="simFhatmu")
    sim.f <- sim.f[,c("pop","simNo","chrom","simFhatmu")]
  ## empirical
    emp.f <- fread("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_emp_lowXhet_no2LtHoms.txt", header=T,
                   colClasses=c("character","character","character","numeric","numeric","numeric","numeric"))
    setnames(emp.f, old="fhatmu", new="empFhatmu")
    emp.f <- emp.f[,c("pop","chrom","empFhatmu")]
  ## bind together
    setkey(sim.f, pop, chrom)
    setkey(emp.f, pop, chrom)
    no2Lt.f.dt <- merge(emp.f, sim.f)
  ## calculate the proportions
    no2Lt.f.dt[,difference:= empFhatmu-simFhatmu]
    no2Lt.f.dt[,ratio:= empFhatmu/simFhatmu]
    no2Lt.f.dt[,normRatio:= difference/simFhatmu]
  ## prep
    pop.dt <- data.table(pop=c("0722", "0819", "fs"),
                         newPop=c("July \n 22", "August \n 19", "all \n individuals"))
    setkey(pop.dt, pop)
    setkey(no2Lt.f.dt, pop)
    no2Lt.f.dt <- merge(pop.dt, no2Lt.f.dt)
    no2Lt.f.dt[,newPop:= factor(newPop, levels=c("July \n 22", "August \n 19", "all \n individuals"))]
 
    
### read in the data before Inv(2L)t homozygote removal
  ## simulations
    sim.f.dt <- foreach(i=c("0722","0819","fs"), .combine="rbind")%do%{
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_",i,"_all.txt",sep=""), header=T)
      temp
    }
    sim.f.dt[,pop:=tstrsplit(class,"_")[1]]
    sim.f.dt[,type:="simulation"]
    sim.f.dt <- sim.f.dt[,c("pop","type","simNo","chrom","fhatmu","fhatmed")]
  ## empirical
    emp.f <- fread("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_emp_lowXhet.txt", header=T,
                   colClasses=c("character","character","character","numeric","numeric","numeric","numeric"))
    emp.f <- emp.f[,c("pop","type","simNo","chrom","fhatmu","fhatmed")]
    emp.f <- emp.f[pop %in% c("0722","0819","fs")]
  ## rearrange the data
    setnames(sim.f.dt, old="fhatmu", new="simFhatmu")
    sim.f.dt <- sim.f.dt[,c("pop","simNo","chrom","simFhatmu")]
    setnames(emp.f, old="fhatmu", new="empFhatmu")
    emp.f <- emp.f[,c("pop","chrom","empFhatmu")]
  ## bind together
    setkey(sim.f.dt, pop, chrom)
    setkey(emp.f, pop, chrom)
    f.norm.dt <- merge(emp.f, sim.f.dt)
  ## get normalized ratio
    f.norm.dt[,difference:= empFhatmu-simFhatmu]
    f.norm.dt[,normRatio:= difference/simFhatmu]
  ## change from number for month-date labels for pops
    pop.dt <- data.table(pop=c("0722", "0819", "fs"),
                         newPop=c("July \n 22", "August \n 19", "all \n individuals"))
    setkey(pop.dt, pop)
    setkey(f.norm.dt, pop)
    f.norm.dt <- merge(pop.dt, f.norm.dt)
    f.norm.dt[,newPop:= factor(newPop, levels=c("July \n 22", "August \n 19", "all \n individuals"))]

    
### FIGURE 7
  a <- ggplot(f.norm.dt[chrom=="WG"], aes(x=newPop, y=normRatio)) + geom_violin(fill="gray") + theme_half_open() + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylim(0.2,1.3) +
    labs(x=NULL, y=expression("relative difference of mean F"["IS"])) + facet_grid(.~newPop, scales="free_x") +panel_border()
  b <- ggplot(no2Lt.f.dt[chrom=="WG"], aes(x=newPop, y=normRatio)) + geom_violin(fill="gray") + theme_half_open() + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylim(0.2,1.3) +
    labs(x=NULL, y=NULL) + facet_grid(.~newPop, scales="free_x") +panel_border()
  plot_grid(a,b, nrow=1, labels="AUTO")
  

  
  