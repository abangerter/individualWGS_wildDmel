### FIGURE 1 ###


### libraries 
library(gdsfmt)
library(SeqArray)
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)


### read in and manipulate the proper data 
  ## read in genofile 
    genofile <- seqOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_seqA.gds")
  ## read in heterozygosity information 
    gt.ag <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/nHets_byChr_v4.txt", header=T)
  ## read in RD info & merge with heterozygosity info
    rd.dt <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/filterFiles/v4RDinfo.txt", header=T)
    setkey(rd.dt, sampID, chrom)
    setkey(gt.ag, sampID, chrom)
    het.dt <- merge(gt.ag, rd.dt)
  ## rearrange het.dt to compare X to autosomes 
    het.auto <- het.dt[chrom!="X"][,c("sampID","chrom","medianRD")]
    setnames(het.auto, old=c("chrom","medianRD"), new=c("autosome","autosomeRD"))
    het.x <- het.dt[chrom=="X"][,c("sampID","chrom","medianRD")]
    setnames(het.x, old="medianRD", new="XRD")
    het.x <- het.x[,c("sampID","XRD")]
    setkey(het.auto, sampID)
    setkey(het.x, sampID)
    het.vs <- merge(het.auto, het.x)


### panel A
  a <- ggplot(het.vs, aes(XRD, autosomeRD, color=autosome)) + geom_jitter(shape=21) + geom_abline(slope=1, color="red") + geom_abline(slope=2, color="blue") + 
    theme_classic() + xlab("X-chromosome Read Depth") + ylab("Autosomal Read Depth") + labs(color="Autosome") 

### panel B
  het.dt[,hetCat:=ifelse(sampID=="CM.039.0722", "highlyHet", "avgHet")]
  b <- ggplot(het.dt[chrom!="X"], aes(x=log2(medianRD), y=Nhets, shape=chrom, color=hetCat)) + geom_point() + theme_classic() + 
    scale_y_continuous(expand=c(0,0), limits=c(0,90000)) + xlab("log(Median Read Depth)") + ylab("N Heterozygous Sites") + 
    guides(color=FALSE) + labs(shape="Autosome") 

### panel C
  c <- ggplot(het.dt[chrom=="X"], aes(x=medianRD, y=Nhets)) + geom_point(alpha=0.5) + geom_hline(yintercept=1000, linetype="dashed") + theme_classic() + 
    xlab("Median Read Depth") + ylab("N Heterozygous Sites") + scale_y_log10() + coord_trans(y="log10")
  
### FIGURE 1
  plot_grid(a,b,c, labels="AUTO")



