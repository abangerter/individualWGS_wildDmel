### FIGURE 9 ### 


### notes to self, requires the following files:
  # slidingWindow_v2_NeiFis_POP_CHROM_lowXhet_1000.txt
  # slidingWindow_v2_NeiFis_emp_lowXhet_1000.txt


### libraries 
library(gdsfmt)
library(SeqArray)
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)


### Read in the sliding window already completed where the windows match perfectly
  ## read in the simulations and empirical SWs
    slide.dt <- foreach(i=c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs"), .combine="rbind") %do% {
      foreach(j=c("2L","2R","3L","3R"), .combine="rbind") %do% {
        print(paste(i,j,sep="_"))
        temp.dt <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/analysis/slidingWindow_v2_NeiFis_",i,"_",j,"_lowXhet_1000.txt",sep=""),
                         colClasses=c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character"))
        temp.dt
      }
    }
    slide.emp <- fread("/scratch/ab5dr/wildDmel2016/simulations/analysis/slidingWindow_v2_NeiFis_emp_lowXhet_1000.txt", header=T)
  ## prep to combine 
    slide.dt[,type:="simulation"]
    slide.emp[,simNo:="NA"]
    slide.emp[,type:="empirical"]
    slide.emp <- slide.emp[,c("simNo","fhatmu","fhatmed","NposF","NnegF","NzeroF","nT","nWindow","minPos","maxPos","chrom","pop","type")]
  ## combine to plot SW simply
    slide.rbind <- rbindlist(list(slide.dt,slide.emp))
    slide.rbind[,avgPos:= (minPos+maxPos)/2]
  ## combine to compare 
    setnames(slide.dt, old=c("fhatmu","fhatmed","NposF","NnegF","NzeroF","nT","nWindow"), new=c("fhatmuSim","fhatmedSim","NposFsim","NnegFsim","NzeroFsim","nTsim","nWindowSim"))
    setnames(slide.emp, old=c("fhatmu","fhatmed","NposF","NnegF","NzeroF","nT","nWindow"), new=c("fhatmuEmp","fhatmedEmp","NposFemp","NnegFemp","NzeroFemp","nTemp","nWindowEmp"))
    setkey(slide.dt, pop, chrom, minPos, maxPos)
    setkey(slide.emp, pop, chrom, minPos, maxPos)
    slide.merge <- merge(slide.dt[simNo<1001], slide.emp)


### get summary of SW
  ## aggregate to get sums of when observed > empirical per simulation -- get counts per 1000 comparisons 
    slide.ag <- slide.merge[,list(propOgtS=sum(fhatmuEmp>fhatmuSim),
                                  totalComps=sum(fhatmuSim<1000)),
                            by=list(pop, chrom, minPos, maxPos)]
  ## get proportion of times that empirical >> simulations FIS 
    slide.ag[,actPropOgtS:=propOgtS/totalComps]
  ## summarize those proportions
    slide.summ <- slide.ag[,list(meanActPropOgtS=mean(actPropOgtS),
                                 sdActPropOgtS=sd(actPropOgtS)),
                           by=list(pop,chrom)]

    
### Reassign pop
  pop.dt <- data.table(pop=c("0624", "0708", "0722", "0819", "0902", "0916", "1003", "1014", "1028", "1111", "1203", "fs"),
                       newPop1=c("June \n24", "July \n8", "July \n22", "Aug. \n19", "Sept. \n2", "Sept. \n16", "Oct. \n3", 
                                "Oct. \n14", "Oct. \n28", "Nov. \n11", "Dec. \n3", "all \nindividuals"),
                       newPop2=c("June 24", "July 8", "July 22", "Aug.  19", "Sept. 2", "Sept. 16", "Oct. 3", 
                                 "Oct. 14", "Oct. 28", "Nov. 11", "Dec. 3", "all individuals"))
  setkey(pop.dt, pop)
  setkey(slide.summ, pop)
  setkey(slide.rbind, pop)
  slide.summ <- merge(pop.dt, slide.summ)
  slide.rbind <- merge(pop.dt, slide.rbind)
  slide.summ[,newPop1:= factor(newPop1, levels=c("June \n24", "July \n8", "July \n22", "Aug. \n19", "Sept. \n2", "Sept. \n16", "Oct. \n3", 
                                                "Oct. \n14", "Oct. \n28", "Nov. \n11", "Dec. \n3", "all \nindividuals"))]
  slide.rbind[,newPop2:= factor(newPop2, levels=c("June 24", "July 8", "July 22", "Aug.  19", "Sept. 2", "Sept. 16", "Oct. 3", 
                                                  "Oct. 14", "Oct. 28", "Nov. 11", "Dec. 3", "all individuals"))]
    

### plot 
  ## panel A - SW summary 
    a <- ggplot(slide.summ) + geom_point(aes(x=newPop1, y=meanActPropOgtS)) + 
      geom_linerange(aes(x=newPop1, ymin=meanActPropOgtS-sdActPropOgtS, ymax=meanActPropOgtS+sdActPropOgtS)) + 
      theme_half_open() + facet_grid(chrom~.) + labs(x=NULL, y="Proportion of Emp mean FIS > Sim FIS") +
      geom_hline(yintercept=0.5, linetype="dashed")
  ## panel B - FS example 
    cb_palette <- c("#440154FF","#95D840FF") 
    b <- ggplot(slide.rbind[newPop2=="all individuals"][chrom=="2L"], aes(x=avgPos, y=fhatmu, group=simNo, color=type)) + 
      geom_line(data=slide.rbind[newPop2=="all individuals"][chrom=="2L"][type=="simulation"], alpha=0.3) + 
      geom_line(data=slide.rbind[newPop2=="all individuals"][chrom=="2L"][type=="empirical"], alpha=1) + theme_half_open() + theme(legend.position = "none") +
      labs(x=NULL, y="") + facet_grid(chrom~newPop2) + scale_color_manual(values=cb_palette)
  ## panel C - 1028 example
    c <- ggplot(slide.rbind[newPop2=="Oct. 28"][chrom=="3L"], aes(x=avgPos, y=fhatmu, group=simNo, color=type)) + 
      geom_line(data=slide.rbind[newPop2=="Oct. 28"][chrom=="3L"][type=="simulation"], alpha=0.3) + 
      geom_line(data=slide.rbind[newPop2=="Oct. 28"][chrom=="3L"][type=="empirical"], alpha=1) + theme_half_open() + theme(legend.position = "none") +
      labs(x=NULL, y=expression("window mean F"["IS"])) + facet_grid(chrom~newPop2) + scale_color_manual(values=cb_palette)
  ## panel D - 0624 example
    d <- ggplot(slide.rbind[newPop2=="June 24"][chrom=="2R"], aes(x=avgPos, y=fhatmu, group=simNo, color=type)) + 
      geom_line(data=slide.rbind[newPop2=="June 24"][chrom=="2R"][type=="simulation"], alpha=0.3) + 
      geom_line(data=slide.rbind[newPop2=="June 24"][chrom=="2R"][type=="empirical"], alpha=1) + theme_half_open() + theme(legend.position = "none") +
      labs(x=NULL, y="") + facet_grid(chrom~newPop2) + scale_color_manual(values=cb_palette)
  ## plot 
    plot_grid(a,b,c,d, ncol=1, rel_heights=c(3,1,1,1), labels=c("A","B","",""))
    
