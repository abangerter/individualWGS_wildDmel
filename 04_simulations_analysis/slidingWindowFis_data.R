### DO SLIDING WINDOW ON EMPIRICAL DATA ###


### libraries
library(data.table)
library(foreach)


### read in the data 
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

  
### subset to LD pruned snpset
  all.snp <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/LDpruned_snpset.txt", header=T)
  all.snp[,snp.id:=paste(chromosome, position, sep="_")]
  all.snp <- all.snp[,c("snp.id")]
  setkey(all.snp, snp.id)
  setkey(emp.dt, snp.id)
  emp.dt <- merge(all.snp, emp.dt)
  emp.dt[,c("chrom","pos"):=tstrsplit(snp.id,"_")]
  emp.dt[,pos:=as.numeric(pos)]

  
### do the sliding window
  ## loop
    slide.emp <- foreach(i=c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs")) %do%{
      foreach(j=c("2L","2R","3L","3R"), .combine="rbind") %do%{
        window.dt <- fread(paste("/scratch/ab5dr/wildDmel2016/simulations/scripts/",j,"_posArray.txt",sep=""), header=T)
        windows <- nrow(window.dt)
        print(paste(i,j,sep="_"))
        foreach(k=1:windows, .combine="rbind") %do%{
          start <- window.dt[k]$start
          stop <- window.dt[k]$stop
          # subset down to window for the sliding window
            temp.sub <- emp.dt[pop==i][chrom==j][pos>=start][pos<=stop]
          # get fhat
            temp.out <- temp.sub[,list(fhatmu=mean(neiFis, na.rm=T),
                                       fhatmed=median(neiFis, na.rm=T),
                                       NposF=sum(neiFis>0, na.rm=T),
                                       NnegF=sum(neiFis<0, na.rm=T),
                                       NzeroF=sum(neiFis==0, na.rm=T),
                                       nT=sum(neiFis<5, na.rm=T),
                                       nWindow=.N,
                                       minPos=min(pos),
                                       maxPos=max(pos))]
          # prep output table
            temp.out[,chrom:=j]
            temp.out[,pop:=i]
          # out
            temp.out
        }
      }
    }
    slide.emp <- rbindlist(slide.emp)
  ## write out the table 
    write.table(slide.emp, "/scratch/ab5dr/wildDmel2016/simulations/analysis/slidingWindow_v2_NeiFis_emp_lowXhet_1000.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)
