### BUILD ARRAY FOR SIMULATION POPULATIONS ###


### libraries 
library(data.table)
library(foreach)


### get the subset of individuals without the high X het individuals 
  ## read in relevant data and manipulate
    gt.ag <- fread("/scratch/ab5dr/wildDmel2016/redoAnalysis/nHets_byChr_v4.txt", header=T)
    gt.ag <- gt.ag[! sampID %in% c("CM.162.0902","CM.024.0708","CM.039.0722")]
  ## get the list of samples to keep 
    lowXhet.samps <- gt.ag[chrom=="X"][Nhets<1000][,c("sampID")]
    lowXhet.samps[,pop:=tstrsplit(sampID,"[.]")[3]]

### add full sample to array
  array.list <- list(fs_lowXhet=lowXhet.samps$sampID)

### add in each population to array
  pops <- unique(lowXhet.samps$pop)
  for(i in pops) {
    array.list[[paste(i,"lowXhet",sep="_")]] <- lowXhet.samps[pop==i]$sampID
  }

### write out array
  saveRDS(array.list, file="/scratch/ab5dr/wildDmel2016/simulations/misc/simArray.rds")
