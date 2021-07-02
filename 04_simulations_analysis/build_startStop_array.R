
### BUILD POSITION ARRAY ###


### libraries
library(data.table)


### build the input start and stop array by chromosome 
  ## read in the data 
    emp.fs.dt <- fread("/scratch/ab5dr/wildDmel2016/simulations/analysis/emp_fs_f_HWE_lowXhet_neiFis.txt", header=T)
    emp.fs.dt[,pop:="fs"]
    emp.fs.dt <- emp.fs.dt[,c("snp.id","pop","nRR","nRA","nAA","nT","chrom","p","q","neiFis","HWEpval")]
    emp.pop.dt <- fread("/scratch/ab5dr/wildDmel2016/simulations/analysis/emp_byPop_f_HWE_lowXhet_neiFis.txt", header=T,
                        colClasses=c("character","character","numeric","numeric","numeric","numeric","character","numeric","numeric","numeric","numeric","numeric","numeric"))
    emp.pop.dt <- emp.pop.dt[,c("snp.id","pop","nRR","nRA","nAA","nT","chrom","p","q","neiFis","HWEpval")]
    emp.dt <- rbindlist(list(emp.fs.dt, emp.pop.dt))
    emp.dt[,pos:=tstrsplit(snp.id,"_")[2]]
    emp.dt[,pos:=as.numeric(pos)]
    setkey(emp.dt, pos, pop)
  ## subset to LD pruned snpset
    all.snp <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/LDpruned_snpset.txt", header=T)
    all.snp[,snp.id:=paste(chromosome, position, sep="_")]
    all.snp <- all.snp[,c("snp.id")]
    setkey(all.snp, snp.id)
    setkey(emp.dt, snp.id)
    emp.dt <- merge(all.snp, emp.dt)
  ## make the start and stop window
    # prep
      window_size <- 1000
      window_step <- 500
    # loop
      for(j in c("2L","2R","3L","3R")) {
        # subset 
          temp <- emp.dt[pop=="fs"][chrom==j]
          setkey(temp, pos)
        # figure out the steps
          nSNPs <- nrow(temp)
          starts <- seq(1,(nSNPs-window_step),window_step)
        # build in an index for picking windows
          temp[,N:=1:nSNPs]
          setkey(temp, N)
        # grab windows 
          stops <- starts+window_size
          temp.out <- data.table(start=temp[N %in% starts]$pos, stop=temp[N %in% stops]$pos)
          temp.out[nrow(temp.out),stop:=max(temp$pos)]
        # out
          write.table(temp.out, paste("/scratch/ab5dr/wildDmel2016/simulations/scripts/",j,"_posArray.txt",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE)
      }
