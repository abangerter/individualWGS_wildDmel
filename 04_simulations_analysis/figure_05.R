### FIGURE 5 ###


### libraries
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
library(scales)


### read in data and combine for A 
  ## simulated observed 
    sim.f.dt <- foreach(i=c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs"), .combine="rbind")%do%{
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_",i,"_all.txt",sep=""), header=T)
      temp
    }
    sim.f.dt[,pop:=tstrsplit(class,"_")[1]]
    sim.f.dt[,type:="simulation"]
    sim.f.dt <- sim.f.dt[,c("pop","type","simNo","chrom","fhatmu","fhatmed")]
  ## simulated True
    sim.T.f <- foreach(i=c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs"), .combine="rbind")%do%{
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_TsimGT_",i,"_lowXhet.txt",sep=""), header=T)
      temp
    }
    sim.T.f[,pop:=tstrsplit(class,"_")[1]]
    sim.T.f[,type:="TrueSim"]
    sim.T.f <- sim.T.f[,c("pop","type","simNo","chrom","fhatmu","fhatmed")]
  ## empirical
    emp.f <- fread("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_emp_lowXhet.txt", header=T,
                   colClasses=c("character","character","character","numeric","numeric","numeric","numeric"))
    emp.f <- emp.f[,c("pop","type","simNo","chrom","fhatmu","fhatmed")]
  ## merge
    f.dt <- rbindlist(list(emp.f, sim.f.dt, sim.T.f))
  ## prep for plotting 
    f.dt[,emp:=ifelse(type=="empirical",T, F)]
    alpha <- ifelse(f.dt[chrom=="WG"]$emp, 1, 0.3)
    f.dt[,type:=ifelse(type=="empirical","Empirical", ifelse(type=="TrueSim","Simulated True Genotype", "Simulated Observed Genotype"))]
  ## change from number for month-date labels for pops
    pop.dt <- data.table(pop=c("0624", "0708", "0722", "0819", "0902", "0916", "1003", "1014", "1028", "1111", "1203", "fs"),
                         newPop=c("June \n 24", "July \n 8", "July \n 22", "August \n 19", "September \n 2", "September \n 16", "October \n 3", 
                                  "October \n 14", "October \n 28", "November \n 11", "December \n 3", "all \n individuals"))
    setkey(pop.dt, pop)
    setkey(f.dt, pop)
    f.dt <- merge(pop.dt, f.dt)
    f.dt[,newPop:= factor(newPop, levels=c("June \n 24", "July \n 8", "July \n 22", "August \n 19", "September \n 2", "September \n 16", "October \n 3", 
                                           "October \n 14", "October \n 28", "November \n 11", "December \n 3", "all \n individuals"))]

    
### read in data and combine for B
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


### FIGURE 5
  ## A
    a <- ggplot() + geom_violin(dat=f.dt[chrom=="WG"][type=="Simulated True Genotype"], aes(x=newPop, y=fhatmu), color="#95D840FF", fill="#95D840FF") + 
      geom_violin(dat=f.dt[chrom=="WG"][type=="Simulated Observed Genotype"], aes(x=newPop, y=fhatmu), color="#20A387FF", fill="#20A387FF") + 
      geom_point(dat=f.dt[chrom=="WG"][type=="Empirical"], aes(x=newPop, y=fhatmu), color="#440154FF") + 
      theme_half_open() + guides(alpha=FALSE) + labs(x=NULL, y=expression("mean F"["IS"])) + 
      theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
      facet_grid(.~newPop, scales="free_x") + panel_border()
  ## B
    b <- ggplot(f.norm.dt[chrom=="WG"], aes(x=pop, y=normRatio)) + geom_violin(fill="gray") + theme_half_open() +
      labs(x=NULL, y=expression("relative difference of mean F"["IS"])) + scale_y_continuous(labels=scales::number_format(accuracy=0.01)) + 
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background=element_blank(), strip.text.x=element_blank()) + 
      facet_grid(.~pop, scales="free_x") + panel_border()
  ## full figure
    plot_grid(a,b, ncol=1, labels="AUTO")


