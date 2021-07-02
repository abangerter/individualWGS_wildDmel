### BUILD ARRAY FOR COMPUTING CLUSTER SLIDING WINDOW ANALYSIS ###


### libraries
library(data.table)


### build array 
  pops <- c("0624","0708","0722","0819","0902","0916","1003","1014","1028","1111","1203","fs")
  chroms <- c("2L","2R","3L","3R")
  array <- data.table(boot=rep(pops,each=4), chrom=rep(chroms, 12))  ### length of 48
  write.table(array, "/scratch/ab5dr/wildDmel2016/simulations/scripts/slide_array.txt", col.names=F, row.names=F, sep="\t", quote=FALSE)
