### SUPPLEMENTAL FIGURE Y ###


### libraries 
library(data.table)
library(viridisLite)
library(viridis)
library(ggplot2)
library(cowplot)


### read in the data 
  prob.dt <- fread("/scratch/ab5dr/wildDmel2016/simulations/misc/gtProbabilities.txt", header=T)

  
### Supplemental figure Y
  ## panel A - probability called RR homozygote
    a <- ggplot(prob.dt[! is.na(GTprob2)], aes(x=ref.rd, y=alt.rd, fill=GTprob2)) + geom_tile() + scale_fill_viridis() + 
      labs(x="reference allele read depth", y="alternate allele read depth", title="Probability of being called an RR homozygote") + 
      theme_classic() + theme(legend.position="none")
  ## panel B - probability called RA heterozygote
    b <- ggplot(prob.dt[! is.na(GTprob1)], aes(x=ref.rd, y=alt.rd, fill=GTprob1)) + geom_tile() + scale_fill_viridis() +
      labs(x="reference allele read depth", y="alternate allele read depth", title="Probability of being called a heterozygote") +
      theme_classic() + theme(legend.position="none")
  ## panel C - probability called AA homozygote
    c <- ggplot(prob.dt[! is.na(GTprob0)], aes(x=ref.rd, y=alt.rd, fill=GTprob0)) + geom_tile() + scale_fill_viridis() + 
      labs(x="reference allele read depth", y="alternate allele read depth", title="Probability of being called an AA homozygote", fill="probability") +
      theme_classic()
  ## make the final figure
    ## still not quite final because the widths are messed up and I can't figure it out
    plot_grid(a,b,c, labels="AUTO")
    