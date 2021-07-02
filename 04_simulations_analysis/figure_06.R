### FIGURE 6 ###


### libraries 
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)


### read in the data for A
  ## simulations 
    sim.jack <- foreach(i=c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203"), .combine="rbind") %do% {
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_sim_jackknife_",i,"_lowXhet.txt",sep=""), header=T)
      temp
    }
    sim.jack[,pop:=tstrsplit(class,"_")[1]]
    sim.jack <- sim.jack[,c("pop","jackDrop","simNo","fhatmu","fhatmed")]
  ## empirical 
    emp.jack <- foreach(i=c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203"), .combine="rbind") %do% {
      temp <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/neiFis_fhat_emp_jackknife_",i,"_lowXhet.txt",sep=""), header=T,
                    colClasses=c("character","character","numeric","numeric"))
      temp
    }
    emp.jack[,simNo:="emp"]
    emp.jack <- emp.jack[,c("pop","jackDrop","simNo","fhatmu","fhatmed")]
  ## combine
    jack.dt <- rbindlist(list(emp.jack, sim.jack))
    jack.dt[,type:=ifelse(simNo=="emp","empirical","simulation")]
  ## prep for plotting  
    sampOrder <- unique(jack.dt[type=="empirical"][order(pop, fhatmu)][,jackDrop])
    jack.dt[,jackDrop:= factor(jackDrop, levels = sampOrder)]
  ## change from number for month-date labels for pops
    pop.dt <- data.table(pop=c("0624", "0708", "0722", "0819", "0902", "0916", "1003", "1014", "1028", "1111", "1203"),
                         newPop=c("June \n 24", "July \n 8", "July \n 22", "August \n 19", "September \n 2", "September \n 16", "October \n 3", 
                                  "October \n 14", "October \n 28", "November \n 11", "December \n 3"))
    setkey(pop.dt, pop)
    setkey(jack.dt, pop)
    jack.dt <- merge(pop.dt, jack.dt)
    jack.dt[,newPop:= factor(newPop, levels=c("June \n 24", "July \n 8", "July \n 22", "August \n 19", "September \n 2", "September \n 16", "October \n 3", 
                                              "October \n 14", "October \n 28", "November \n 11", "December \n 3"))]
    

### prep for B
  ## rearrange
    setnames(sim.jack, old=c("fhatmu","fhatmed"), new=c("simFhatmu","simFhatmed"))
    sim.jack <- sim.jack[,c("pop","jackDrop","simNo","simFhatmu","simFhatmed")]
    setnames(emp.jack, old=c("fhatmu","fhatmed"), new=c("empFhatmu","empFhatmed"))
    emp.jack <- emp.jack[,c("pop","jackDrop","empFhatmu","empFhatmed")]
  ## merge
    setkey(sim.jack, pop, jackDrop)
    setkey(emp.jack, pop, jackDrop)
    jack.norm.dt <- merge(sim.jack, emp.jack)
  ## get normalized ratio
    jack.norm.dt[,difference:= empFhatmu-simFhatmu]
    jack.norm.dt[,normRatio:= difference/simFhatmu]
  ## prep to plot
    jack.norm.dt[,simNo:=factor(simNo)]
    jack.norm.dt[,jackDrop:= factor(jackDrop, levels = sampOrder)]
    

### plot figure 6
  ## A
    a <- ggplot() + geom_point(data=jack.dt[type=="simulation"], aes(x=jackDrop, y=fhatmu, alpha=0.3), color="#20A387FF") + 
      geom_point(data=jack.dt[type=="empirical"], aes(x=jackDrop, y=fhatmu), color="#440154FF") + 
      facet_grid(.~newPop, scales="free_x") + theme_half_open() + 
      theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
      labs(x=NULL, y=expression("mean F"["IS"]), color="Data Type") + panel_border()
  ## B
    b <- ggplot(jack.norm.dt, aes(x=jackDrop, y=normRatio)) + geom_point(alpha=0.5) + theme_half_open() + facet_grid(.~pop, scales="free_x") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background=element_blank(), strip.text.x=element_blank()) + 
      labs(x=NULL, y=expression("relative difference of mean F"["IS"])) + panel_border()
  ## FIGURE 6
    plot_grid(a,b, ncol=1, labels="AUTO")

    