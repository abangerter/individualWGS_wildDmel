### KINSHIP ESTIMATION/ANALYSIS & FIGURE 8 ###


### libraries
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(cowplot)
library(foreach)
library(ggsignif)


### run SNPRelate kinship
  ## open genofile
    genofile <- snpgdsOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_snpr.gds")
  ## filter down to only 119 individuals 
    gt.ag <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/nHets_byChr_v4.txt", header=T)
    gt.ag <- gt.ag[! sampID %in% c("CM.162.0902","CM.024.0708","CM.039.0722")]
    lowXhet.samps <- gt.ag[chrom=="X"][Nhets<1000]$sampID
  ## king robust 
    ibd.king <- snpgdsIBDKING(genofile, sample.id=lowXhet.samps)
    king.dat <- snpgdsIBDSelection(ibd.king)
    king.dt <- as.data.table(king.dat)
  ## get within vs between
    king.dt[,pop1:=tstrsplit(ID1,"[.]")[3]]
    king.dt[,pop2:=tstrsplit(ID2,"[.]")[3]]
    king.dt[,category:=ifelse(pop1==pop2, "within","between")]

    
### get the table for the per-population within vs between 
### the output table produced by snpgdsIBDKING is only half of the matrix, 
### so the way it's arranged plotting by population will miss a portion of the "between" comparisons
  ## rearrange data to get all betweens represented in a plot-able format
    king.fold <- king.dt[category=="between"]
    setnames(king.fold, old=c("pop1","pop2"), new=c("pop2","pop1"))
    king.fold <- king.fold[,c("ID1","ID2", "IBS0","kinship","pop1","pop2","category")]
    king.fold <- rbindlist(list(king.fold, king.dt))
  ## prep population labels for plotting
    pop.dt <- data.table(pop1=c("0624", "0708", "0722", "0819", "0902", "0916", "1003", "1014", "1028", "1111", "1203"),
                         newPop=c("June \n 24", "July \n 8", "July \n 22", "August \n 19", "September \n 2", "September \n 16", "October \n 3", 
                                  "October \n 14", "October \n 28", "November \n 11", "December \n 3"))
    setkey(pop.dt, pop1)
    setkey(king.fold, pop1)
    king.fold <- merge(pop.dt, king.fold)
    king.fold[,newPop:= factor(newPop, levels=c("June \n 24", "July \n 8", "July \n 22", "August \n 19", "September \n 2", "September \n 16", "October \n 3", 
                                                "October \n 14", "October \n 28", "November \n 11", "December \n 3", "all \n individuals"))]
    
    
### test similarity of distributions of within vs between estimates 
  ## do ttests
    test.all <- t.test(king.dt[category=="within"]$kinship, king.dt[category=="between"]$kinship)
    test.bypop <- foreach(i=c("0624", "0708", "0722", "0819", "0902", "0916", "1003", "1014", "1028", "1111", "1203"), .combine="rbind") %do% {
      bypop.ttest <- t.test(king.fold[category=="within"][pop1==i]$kinship, king.fold[category=="between"][pop1==i]$kinship)
      ttest.out <- data.table(newPop=i, start="between", end="within", pval=bypop.ttest$p.value)
      ttest.out
    }
    test.bypop[,pvalLab:=ifelse(pval<0.001, "**", ifelse(pval<0.05, "*", "NS"))]

  ## merge per population ttests with table for plotting 
    #setkey(test.bypop, pop1)
    #setkey(king.fold, pop1)
    #king.fold <- merge(king.fold, test.bypop)

 
### FIGURE 8
  ## panel A
    a <- ggplot(king.dt, aes(x=IBS0, y=kinship)) + geom_point(alpha=0.3) + theme_classic() + ylab("KING kinship") + xlab("IBS 0")
      # + xlim(0,0.08)
  ## panel B
    cb_palette <- c("#440154FF","#95D840FF") 
    b <- ggplot(king.dt, aes(x=category, y=kinship)) + geom_boxplot(aes(fill=category, alpha=0.2)) + theme_classic() + ylab("KING kinship") + 
      theme(legend.position="none") + xlab("Time Point Comparison") + scale_fill_manual(values=cb_palette) +
      geom_signif(comparisons = list(c("within", "between")), annotation=c("NS"))
  ## panel C --- STILL TROUBLESHOOTING
    c <- ggplot(king.fold, aes(x=pop1, y=kinship)) + geom_boxplot(aes(fill=category, alpha=0.2)) + 
      theme_half_open() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.position="none") + 
      ylab("KING kinship") + facet_grid(.~newPop, scales="free_x") + scale_fill_manual(values=cb_palette)
      
    ggplot(king.fold, aes(x=pop1, y=kinship)) + geom_boxplot(aes(fill=category, alpha=0.2)) + 
      theme_half_open() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.position="none") + 
      ylab("KING kinship") + scale_fill_manual(values=cb_palette) + 
      geom_signif(data=test.bypop, aes(xmin=start, xmax=end, annotations=pvalLab, y=0.1), textsize=3, vjust=-0.2, manual=TRUE) + 
      facet_grid(.~newPop, scales="free_x") 
      
    
    ggplot(king.fold[pop1 %in% c("0624","0722")], aes(x=newPop, y=kinship)) + geom_boxplot(aes(fill=category, alpha=0.2)) + 
      theme_half_open() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.position="none") + 
      ylab("KING kinship") + scale_fill_manual(values=cb_palette) + 
      geom_signif(data=as.data.frame(test.bypop), aes(xmin=start, xmax=end, annotations=pvalLab, y=0), textsize=3, vjust=-0.2, manual=TRUE) + 
      facet_grid(.~newPop, scales="free_x") 
    
    ggplot(king.fold[pop1 %in% c("0624","0722")], aes(x=newPop, y=kinship)) + geom_boxplot(aes(fill=category, alpha=0.2)) + 
      theme_half_open() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.position="none") + 
      ylab("KING kinship") + scale_fill_manual(values=cb_palette) + 
      geom_signif(data=as.data.frame(test.bypop), aes(comparisons=list(c("within","between")), annotation=pvalLab, y=0), manual=TRUE) +
      facet_grid(.~newPop, scales="free_x") 
    
    
    
    
