### BUILD ARRAYS FOR SIMULATION SUBMISSION ON COMPUTING CLUSTER 


### libraries
library(data.table)


### read in the needed info
  array.list <- readRDS("/scratch/ab5dr/wildDmel2016/simulations/misc/simArray.rds")
  codes <- names(array.list)
  all.snp <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/LDpruned_snpset.txt", header=T)
  
  
### build array for Full Sample - every 100 SNPs
  ## build out the start and stop and chromosome stuff
    window <- 100
    all.snp[,idx:=ceiling((1/window)*1:.N), by=chromosome]
    snp.ag.100 <- all.snp[,list(start=min(position),
                                stop=max(position)),
                          by=list(chromosome,idx)]
  ## add in the simulation code 
    snp.out.100 <- data.table(chrom=snp.ag.100$chromosome, 
                              code=codes[1],
                              start=snp.ag.100$start,
                              stop=snp.ag.100$stop)
  ## write it out
    write.table(snp.out.100, "/scratch/ab5dr/wildDmel2016/simulations/scripts/fs_array.txt", col.names=F, row.names=F, sep="\t", quote=FALSE)

    
### build array for each population - every 900 SNPs
  ## build out the start and stop and chromosome stuff
    window <- 900
    all.snp[,idx:=ceiling((1/window)*1:.N), by=chromosome]
    snp.ag.pop <- all.snp[,list(start=min(position),
                                stop=max(position)),
                          by=list(chromosome,idx)]
  ## add in the simulation code 
    snp.out.pop <- data.table(chrom=rep(snp.ag.pop$chromosome, 11), 
                              code=rep(codes[2:12], each=length(snp.ag.pop$chromosome)),
                              start=rep(snp.ag.pop$start, 11),
                              stop=rep(snp.ag.pop$stop, 11))
  ## write it out
    write.table(snp.out.pop, "/scratch/ab5dr/wildDmel2016/simulations/scripts/pop_array.txt", col.names=F, row.names=F, sep="\t", quote=FALSE)
