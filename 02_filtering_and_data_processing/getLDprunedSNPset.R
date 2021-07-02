### GET LIST OF LD PRUNED SNPS ###


### libraries
library(data.table)
library(gdsfmt)
library(SNPRelate)


### read in the data 
  genofile_snpr <- snpgdsOpen("/scratch/ab5dr/wildDmel2016/vcf/wild2016_v4_preFinal_snpr.gds")

### set filtering parameters 
  missing.rate <- 0.2
  threads <- 10
  bp <- 1000
  ld.thresh <- 0.3
  maf <- 0.05

### do the pruning
  set.seed(1000)
  snpset <- snpgdsLDpruning(genofile_snpr, ld.threshold=ld.thresh, autosome.only=TRUE, maf=maf, missing.rate=missing.rate, slide.max.bp=bp)

### get the list of snps 
  snpset.id <- unlist(unname(snpset))
  all.snp <- as.data.table(snpgdsSNPList(genofile_snpr))
  all.snp <- all.snp[snp.id %in% snpset.id]
  all.snp <- all.snp[,c("chromosome", "position")][chromosome %in% c("2L","2R","3L","3R")]

### write out the snpset 
  write.table(all.snp, "/scratch/ab5dr/wildDmel2016/simulations/misc/LDpruned_snpset.txt", col.names=T, row.names=F, sep="\t", quote=FALSE)
